# VoP (variation-of-parameters) transport for elliptic 2×2 sectors.
#
# Packages the per-graph Python/mpmath quadrature (~15 min/point) as a single fast
# q-series call (~ms/point) on top of the existing iterated.jl LogQSeries machinery,
# with a certified-quadrature fallback for kernels known only as functions of t.
#
# Physics (docs/conventions.md §4, §6):  for a 2×2 elliptic top sector with periods
# (Ψ₁, Ψ₂) and Wronskian W,  dτ = (W/Ψ₁²) dt  so τ(t) is the period ratio, and the
# canonical solution of the inhomogeneous DE is
#
#     g_i(t) = pref(t) · ( c₁ᵢ + c₂ᵢ·τ(t) + Eichler[τ; Sᵢ] ),
#     Eichler[τ; S] = ∫_{i∞}^{τ} S(τ') dτ'  (tangential-basepoint regularised),
#
# with S a (quasi-)modular q-series.  The cusp regularisation is exactly the
# LogQSeries one (iterated.jl): the constant Fourier term a₀ of S integrates to
# a₀·τ (folded into c₂), the q^{n≥1} terms to a finite q-series (folded into the
# Eichler value), and the integration constant is c₁.  The (c₁, c₂) per source are
# fixed by ≥2 oracle anchors via `vop_fit_constants`.

# ---- sector descriptor --------------------------------------------------------

"""
    EllipticSector

Descriptor for one elliptic 2×2 sector.  Carries:

  * `prec`, `level` (the Γ₁(N) level of the q-series basis), `ctx` (the shared
    `SeriesContext` for kernel q-series),
  * `cusp` — the kinematic t at the q → 0 cusp (the tangential basepoint), and
  * `tau_of_t :: t::Acb → (q, L, τ)` on the cusp-connected branch, with L = log q
    = 2πi τ supplied explicitly (so the t+i0 prescription is honoured).

Construct via one of the three convenience constructors below; for an arbitrary
graph supply `tau_of_t` directly.
"""
struct EllipticSector
    prec::Int
    level::Int
    ctx::SeriesContext
    cusp::Acb
    tau_of_t::Function     # t::Acb -> (q::Acb, L::Acb, τ::Acb)
end

"""Sector backed by the Γ₁(6) layer (Hauptmodul inversion, certified)."""
EllipticSector(G::Gamma16; cusp = Acb(0; prec = G.ctx.prec)) =
    EllipticSector(G.ctx.prec, 6, G.ctx, cusp,
        t -> nome_and_log_from_t(G, Acb(t; prec = G.ctx.prec)))

"""Sector from a raw period map `periods(t) = (Ψ₁, Ψ₂, W)`.  τ = Ψ₂/Ψ₁; the nome
branch is the principal q = exp(2πiτ) (caller is responsible for ensuring the
cusp-connected sheet)."""
function EllipticSector(prec::Int, Nq::Int, level::Int, periods;
                        cusp, x0::Real = 0.25)
    ctx = SeriesContext(prec, Nq; x0 = x0)
    τmap = t -> begin
        ψ1, ψ2, _ = periods(Acb(t; prec))
        τ = ψ2 / ψ1
        Arblib.is_positive(imag(τ)) || error("EllipticSector: Im τ ≤ 0 at t = $t")
        L = 2 * Acb(0, 1; prec) * Acb(Arb(π; prec); prec) * τ
        (nome(τ), L, τ)
    end
    EllipticSector(prec, level, ctx, Acb(cusp; prec), τmap)
end

"""Sector from a parametric `QuarticCurve` family `curve_at(t)::QuarticCurve` (uses
`curve_periods` for (Ψ₁, Ψ₂); W is not needed for the τ-map)."""
function EllipticSector(prec::Int, Nq::Int, level::Int; curve_at, cusp, x0::Real = 0.25)
    periods = t -> begin
        C = curve_at(t)
        ψ1, ψ2 = curve_periods(C)
        (ψ1, ψ2, Acb(0; prec))
    end
    EllipticSector(prec, Nq, level, periods; cusp = cusp, x0 = x0)
end

"""(τ, q, log q) at kinematic point `t` (cusp-connected branch)."""
function tau_map(sec::EllipticSector, t, prec::Int = sec.prec)
    q, L, τ = sec.tau_of_t(Acb(t; prec))
    return (τ, q, L)
end

# ---- Eichler primitive & cusp regularisation ---------------------------------

"""Period-polynomial constants of a `LogQSeries` primitive: the constant Fourier
term of each log-level, `pp[j+1] = F_j(0)`.  For depth 1 this is `[0, a₀]` where
a₀ is the kernel's constant term (the τ-linear piece that folds into c₂)."""
period_polynomial(V::LogQSeries) = [qcoeff(lev, 0) for lev in V.levels]

"""Tangentially-regularised value of a `LogQSeries` at the cusp q = 0: only the
pure q-power-series level `F₀(0)` survives (the log^j q pieces are absorbed into
the period polynomial / boundary constants)."""
function evaluate_cusp(V::LogQSeries)
    prec = V.levels[1].ctx.prec
    return Acb(qcoeff(V.levels[1], 0); prec)
end

"""
    eichler_transport(sec, S, t_start, t_end, prec; cusp_reg=true)

Certified Eichler transport ∫_{τ(t_start)}^{τ(t_end)} S(τ) dτ of the q-series
kernel `S` (on `sec.ctx`).  Returns `(value, period_poly)`:

  * `value` — the transported integral (Acb ball);
  * `period_poly` — the period-polynomial constants `[F_j(0)]_{j≥0}` of the
    primitive (independent of the endpoints; returned so the caller can fold
    `F_1(0) = a₀` into c₂ explicitly when desired).

If `cusp_reg` and `t_start == sec.cusp`, the lower endpoint is taken as the
tangential basepoint (q → 0, regularised); otherwise both endpoints are mapped
through `tau_map` and the primitive is evaluated at each.

Fast path: builds the `LogQSeries` primitive once via `iterated_integral([S])`
and evaluates it — no quadrature.
"""
function eichler_transport(sec::EllipticSector, S::QSeries, t_start, t_end,
                           prec::Int = sec.prec; cusp_reg::Bool = true)
    @assert S.ctx === sec.ctx "kernel must live on the sector's SeriesContext"
    V = iterated_integral([S])
    pp = period_polynomial(V)
    _, q_end, L_end = tau_map(sec, t_end, prec)
    val_end = evaluate(V, q_end, L_end)
    if cusp_reg && Arblib.overlaps(Acb(t_start; prec), sec.cusp)
        return (val_end - evaluate_cusp(V), pp)
    end
    _, q_start, L_start = tau_map(sec, t_start, prec)
    val_start = evaluate(V, q_start, L_start)
    return (val_end - val_start, pp)
end

# ---- one-call VoP assembly ---------------------------------------------------

"""
    vop_assemble(sec, pref, sources, c1, c2, t, prec=sec.prec)

One-call replacement for the per-graph Python VoP:  for each source kernel
`sources[i]` returns

    g_i(t) = pref(t) · ( c1[i] + c2[i]·τ(t) + Eichler_reg[τ(t); sources[i]] ),

with the Eichler integral taken from the sector's cusp (tangential basepoint).
`pref(t)::Acb` is the algebraic prefactor (typically 1/Ψ₁ or Ψ₁ depending on the
chosen UT normalisation).  All values certified Acb.
"""
function vop_assemble(sec::EllipticSector, pref, sources::Vector{QSeries},
                      c1::Vector, c2::Vector, t, prec::Int = sec.prec)
    @assert length(sources) == length(c1) == length(c2)
    τ, q, L = tau_map(sec, t, prec)
    P = pref(Acb(t; prec))
    out = Vector{Acb}(undef, length(sources))
    for i in eachindex(sources)
        @assert sources[i].ctx === sec.ctx
        V = iterated_integral([sources[i]])
        eich = evaluate(V, q, L) - evaluate_cusp(V)
        out[i] = P * (Acb(c1[i]; prec) + Acb(c2[i]; prec) * τ + eich)
    end
    return out
end

"""
    vop_fit_constants(sec, pref, sources, anchors, prec=sec.prec)
        → (c1::Vector{Acb}, c2::Vector{Acb}, residuals)

Fit the two boundary constants per source from oracle anchors
`anchors = [(t_a, g_a::Vector), …]` (each `g_a` of length `length(sources)`).
The first two anchors fix (c₁, c₂) by the 2×2 linear system

    g_a / pref(t_a)  −  Eichler_reg[τ(t_a); S_i]  =  c₁ᵢ + c₂ᵢ · τ(t_a),

solved by Cramer (certified Acb arithmetic; errors if the τ's coincide).  Any
remaining anchors are returned as `residuals[k][i]` = predicted − supplied
(should contain 0 for held-out CV).
"""
function vop_fit_constants(sec::EllipticSector, pref, sources::Vector{QSeries},
                           anchors::Vector, prec::Int = sec.prec)
    length(anchors) >= 2 || error("vop_fit_constants: need ≥2 anchors")
    ns = length(sources)
    # Precompute (τ, eich_i, rhs_i) at every anchor.
    prims = [iterated_integral([S]) for S in sources]
    data = map(anchors) do (ta, ga)
        τ, q, L = tau_map(sec, ta, prec)
        P = pref(Acb(ta; prec))
        eich = [evaluate(prims[i], q, L) - evaluate_cusp(prims[i]) for i in 1:ns]
        rhs = [Acb(ga[i]; prec) / P - eich[i] for i in 1:ns]
        (τ, rhs)
    end
    τ1, r1 = data[1]; τ2, r2 = data[2]
    det = τ2 - τ1
    Arblib.contains_zero(det) && error("vop_fit_constants: anchor τ's coincide")
    c1 = Vector{Acb}(undef, ns); c2 = Vector{Acb}(undef, ns)
    for i in 1:ns
        c2[i] = (r2[i] - r1[i]) / det
        c1[i] = r1[i] - c2[i] * τ1
    end
    residuals = Vector{Vector{Acb}}()
    for k in 3:length(anchors)
        τk, rk = data[k]
        push!(residuals, [(c1[i] + c2[i] * τk) - rk[i] for i in 1:ns])
    end
    return c1, c2, residuals
end

# ---- quadrature fallback (slow path, validation only) ------------------------

"""Gauss–Legendre nodes and weights on [-1,1] at `prec` bits, order `n`, computed by
Newton on the three-term recurrence (non-certified; for cross-check only)."""
function _gl_nodes(n::Int, prec::Int)
    setprecision(BigFloat, prec)
    x = Vector{Arb}(undef, n); w = Vector{Arb}(undef, n)
    for k in 1:n
        # Tricomi initial guess
        ξ = (1 - (n - 1) / (8n^3)) * cos(big(π) * (4k - 1) / (4n + 2))
        for _ in 1:64
            p0 = big(1.0); p1 = ξ
            for j in 2:n
                p0, p1 = p1, ((2j - 1) * ξ * p1 - (j - 1) * p0) / j
            end
            dp = n * (ξ * p1 - p0) / (ξ^2 - 1)
            δ = p1 / dp; ξ -= δ
            abs(δ) < big(2.0)^(-prec) && break
        end
        # weight: 2 / ((1-ξ²) P_n'(ξ)²)
        p0 = big(1.0); p1 = ξ
        for j in 2:n
            p0, p1 = p1, ((2j - 1) * ξ * p1 - (j - 1) * p0) / j
        end
        dp = n * (ξ * p1 - p0) / (ξ^2 - 1)
        x[k] = Arb(ξ; prec); w[k] = Arb(2 / ((1 - ξ^2) * dp^2); prec)
    end
    return x, w
end

"""
    vop_quadrature(sec, R_of_t, t_start, t_end, prec; method=:certified, n=64)

Slow-path transport ∫_{t_start}^{t_end} R(t) dt for a kernel given as a function
of t (not a q-series).  Methods:

  * `:certified` — `Arblib.integrate` with the analytic-box flag (rigorous Acb
    enclosure; `R_of_t` must accept `(t; analytic::Bool)` and return an
    indeterminate ball when analyticity on the box cannot be certified, cf.
    `curve.jl`);
  * `:gauss_legendre` — fixed-order GL on `n` nodes (NOT certified; matches the
    Python mpmath path for regression).

Kept for validation only; the q-series path in `eichler_transport` is the fast one.
"""
function vop_quadrature(sec::EllipticSector, R_of_t, t_start, t_end,
                        prec::Int = sec.prec; method::Symbol = :certified, n::Int = 64)
    a = Acb(t_start; prec); b = Acb(t_end; prec)
    if method === :certified
        return Arblib.integrate(R_of_t, a, b; check_analytic = true, prec = prec,
                                rtol = exp10(-0.301 * prec))
    elseif method === :gauss_legendre
        xs, ws = _gl_nodes(n, prec)
        half = (b - a) / 2; mid = (a + b) / 2
        s = Acb(0; prec)
        for k in 1:n
            tk = mid + half * Acb(xs[k]; prec)
            s += Acb(ws[k]; prec) * R_of_t(tk)
        end
        return half * s
    else
        error("vop_quadrature: unknown method $method")
    end
end

# ==== multi-curve (coupled-elliptic) VoP transport ============================
#
# Next-frontier sectors (gg→tt̄ 3-curve top sector, soaring kite) couple SEVERAL
# elliptic curves over a common kinematic base t: the canonical DE mixes periods
# of curve i with sources living on curve j.  In the τ-picture this is
#
#     dg_i/dτ_i  =  S_{ii}(τ_i) g_i  +  Σ_{j≠i} S_{ij}(τ_i) · K_{ij}(t) · g_j ,
#
# i.e. the off-diagonal kernel factorises into a (quasi-)modular piece in τ_i
# and a piece that depends on t only through the OTHER curves' data (periods,
# Hauptmoduli, …).  Fully general iterated integrals over multiple modular
# parameters are Brown's multiple modular values / elliptic MZVs and are out of
# scope here.  What the physical 3-curve sectors actually need — and what is
# implemented — is the **shared-base** scheme:
#
#   (i)  DIAGONAL  S_{ii}: per-curve `eichler_transport` exactly as in the
#        single-sector code (fast q-series path).
#   (ii) OFF-DIAGONAL S_{ij}, i≠j: change variables back to the common base t,
#             ∫ S_{ij}(q_i(t)) · K_{ij}(t) · (dL_i/dt) dt ,   L_i = log q_i = 2πi τ_i,
#        and evaluate by the certified-quadrature fallback `vop_quadrature`
#        (so the integration variable matches the diagonal `iterated_integral`
#        convention dq/q = dL).  When K is t-independent this reduces
#        analytically to K_{ij} · eichler_transport(curve_i, S_{ij}, …) and is
#        short-circuited.
#
# Scope/limitations: depth-1 only (one Eichler integral per matrix entry — no
# nested ∫∫ across different τ's); off-diagonal endpoints must be interior
# (cusp-regularised lower limit is supported only on the diagonal); the Jacobian
# dτ_i/dt must be supplied per curve (no default for the raw-period sector).

"""
    MultiEllipticSector

Descriptor for n coupled elliptic 2×2 sectors sharing one kinematic base t.

Fields:
  * `curves`   — the n single-curve `EllipticSector`s (each with its own τ-map
    and q-series context);
  * `coupling` — `t::Acb → Matrix{Acb}` giving the off-diagonal kernel factor
    K_{ij}(t) (size n×n; diagonal entries are ignored).  Pass `nothing` for the
    fully decoupled case;
  * `dlogq_dt` — per-curve Jacobian `t::Acb → dL_i/dt::Acb` with L_i = log q_i
    = 2πi τ_i (length-n vector; entries may be `nothing` if that curve never
    appears as the row index of an off-diagonal quadrature).  This matches the
    integration variable of `iterated_integral`/`eichler_transport`.

Construct via `MultiEllipticSector(curves; coupling, dtau_dt)`.
"""
struct MultiEllipticSector
    curves::Vector{EllipticSector}
    coupling::Union{Nothing,Function}
    dlogq_dt::Vector{Any}
    prec::Int
end

function MultiEllipticSector(curves::Vector{EllipticSector};
                             coupling = nothing,
                             dlogq_dt = fill(nothing, length(curves)))
    isempty(curves) && error("MultiEllipticSector: need ≥1 curve")
    prec = curves[1].prec
    MultiEllipticSector(curves, coupling, collect(Any, dlogq_dt), prec)
end

"""(τᵢ, qᵢ, log qᵢ) for every curve at the common kinematic point `t`."""
multi_tau_map(msec::MultiEllipticSector, t, prec::Int = msec.prec) =
    [tau_map(c, t, prec) for c in msec.curves]

"""d(log q)/dt = 2πi·dτ/dt for the Γ₁(6) Hauptmodul sector.  Since
dt/dL = q·dt/dq (L = log q), this is `(evaluate(qdq(G.t), q))⁻¹` on the
cusp-connected branch."""
function dlogq_dt_gamma16(G::Gamma16)
    tder = qdq(G.t)
    return t -> begin
        q, _, _ = nome_and_log_from_t(G, Acb(t; prec = G.ctx.prec))
        inv(evaluate(tder, q))
    end
end

# For the `:certified` quadrature path we need an integrand that honours the
# `analytic` keyword.  The τ-map (`nome_from_t`) does interval Newton on a wide
# box and may legitimately fail to certify; in that case return an indeterminate
# ball so `Arblib.integrate` bisects.
function _safe_analytic(f, t::Acb, analytic::Bool, prec::Int)
    if analytic
        try
            return f(t)
        catch
            ind = Acb(prec = prec); Arblib.indeterminate!(ind); return ind
        end
    end
    return f(t)
end

"""
    multi_eichler_transport(msec, S, t_start, t_end, prec=msec.prec;
                            method=:certified, n=64)
        → Matrix{Tuple{Acb,Vector}}

Matrix-valued Eichler transport for a `MultiEllipticSector`.  `S` is an n×n
matrix whose (i,j) entry is a `QSeries` on `msec.curves[i].ctx` (or `nothing`
to skip).  Returns an n×n matrix of `(value, period_poly)` pairs:

  * diagonal (i,i): `eichler_transport(curves[i], S[i,i], t_start, t_end)`
    (fast q-series; cusp-regularised if `t_start == curves[i].cusp`);
  * off-diagonal (i,j), i≠j: with K = `msec.coupling`,
      - K ≡ nothing  →  `(0, [])` (decoupled);
      - otherwise  →  ∫_{t_start}^{t_end} S_{ij}(q_i(t)) · K(t)_{ij} · dL_i/dt dt
        via `vop_quadrature` (slow path; `method`/`n` forwarded), L_i = log q_i.
        `period_poly` is returned as `Acb[]` (no clean period polynomial across
        curves).

If `K` is detected to be t-independent (probed at `t_start`, `t_end`, midpoint)
the off-diagonal short-circuits to K_{ij}·eichler_transport(curves[i], S_{ij}, …),
which is exact and avoids the Jacobian.
"""
function multi_eichler_transport(msec::MultiEllipticSector, S::AbstractMatrix,
                                 t_start, t_end, prec::Int = msec.prec;
                                 method::Symbol = :certified, n::Int = 64)
    nc = length(msec.curves)
    size(S) == (nc, nc) || error("multi_eichler_transport: S must be $(nc)×$(nc)")
    out = Matrix{Tuple{Acb,Vector}}(undef, nc, nc)
    # constant-coupling probe
    Kconst = nothing
    if msec.coupling !== nothing
        ta = Acb(t_start; prec); tb = Acb(t_end; prec); tm = (ta + tb)/2
        Ka, Kb, Km = msec.coupling(ta), msec.coupling(tb), msec.coupling(tm)
        if all(Arblib.overlaps(Ka[i,j], Kb[i,j]) && Arblib.overlaps(Ka[i,j], Km[i,j])
               for i in 1:nc, j in 1:nc)
            Kconst = Ka
        end
    end
    for i in 1:nc, j in 1:nc
        Sij = S[i, j]
        if Sij === nothing
            out[i, j] = (Acb(0; prec), Acb[]); continue
        end
        Sij isa QSeries || error("multi_eichler_transport: S[$i,$j] must be a QSeries or nothing")
        Sij.ctx === msec.curves[i].ctx ||
            error("multi_eichler_transport: S[$i,$j] must live on curve $i's context")
        if i == j
            out[i, j] = eichler_transport(msec.curves[i], Sij, t_start, t_end, prec)
        elseif msec.coupling === nothing
            out[i, j] = (Acb(0; prec), Acb[])
        elseif Kconst !== nothing
            v, pp = eichler_transport(msec.curves[i], Sij, t_start, t_end, prec;
                                      cusp_reg = false)
            out[i, j] = (Kconst[i, j] * v, [Kconst[i, j] * c for c in pp])
        else
            J = msec.dlogq_dt[i]
            J === nothing && error("multi_eichler_transport: off-diagonal ($i,$j) with " *
                "t-dependent coupling needs dlogq_dt[$i]")
            ci = msec.curves[i]; K = msec.coupling
            R = (t; analytic::Bool = false) -> _safe_analytic(t, analytic, prec) do tt
                _, q, _ = tau_map(ci, tt, prec)
                evaluate(Sij, q) * K(tt)[i, j] * J(tt)
            end
            v = vop_quadrature(ci, R, t_start, t_end, prec; method = method, n = n)
            out[i, j] = (v, Acb[])
        end
    end
    return out
end

"""
    multi_vop_assemble(msec, pref, sources, c1, c2, t, prec=msec.prec)

Vector-valued analogue of `vop_assemble` for a `MultiEllipticSector`: for each
curve i,

    g_i(t) = pref[i](t) · ( c1[i] + c2[i]·τ_i(t)
                            + Σ_j Eich_{ij}[t; sources[i,j]] ),

with the diagonal Eichler integral cusp-regularised on curve i and the
off-diagonal pieces taken from `multi_eichler_transport` between the curve-i
cusp and t (constant-coupling short-circuit applies).  `pref`, `c1`, `c2` are
length-n vectors; `sources` is the n×n kernel matrix as in
`multi_eichler_transport`.
"""
function multi_vop_assemble(msec::MultiEllipticSector, pref::Vector,
                            sources::AbstractMatrix, c1::Vector, c2::Vector,
                            t, prec::Int = msec.prec; method::Symbol = :certified, n::Int = 64)
    nc = length(msec.curves)
    length(pref) == length(c1) == length(c2) == nc ||
        error("multi_vop_assemble: pref/c1/c2 must have length $nc")
    τs = multi_tau_map(msec, t, prec)
    g = Vector{Acb}(undef, nc)
    for i in 1:nc
        τ, q, L = τs[i]
        P = pref[i](Acb(t; prec))
        acc = Acb(c1[i]; prec) + Acb(c2[i]; prec) * τ
        # diagonal Eichler (cusp-regularised on curve i)
        if sources[i, i] !== nothing
            V = iterated_integral([sources[i, i]])
            acc += evaluate(V, q, L) - evaluate_cusp(V)
        end
        # off-diagonal contributions from curve-i cusp to t
        for j in 1:nc
            (j == i || sources[i, j] === nothing) && continue
            row = multi_eichler_transport(msec,
                [k==i && l==j ? sources[i,j] : nothing for k in 1:nc, l in 1:nc],
                msec.curves[i].cusp, t, prec; method = method, n = n)
            acc += row[i, j][1]
        end
        g[i] = P * acc
    end
    return g
end

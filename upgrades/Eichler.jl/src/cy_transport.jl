# Calabi–Yau / K3 period + transport layer.
#
# Extends the elliptic VoP (vop_transport.jl) up the CY ladder n = 1 (elliptic),
# n = 2 (K3), n = 3 (CY3), …  For a CY n-fold the Picard–Fuchs operator L has
# order r = n+1 and a MUM (maximal-unipotent-monodromy) point at z = 0: the
# indicial polynomial is (ρ − α)^r for an integer α, and the Frobenius basis is
#
#     ϖ_k(z) = z^α · Σ_{j=0}^{k} (log^j z / j!) · h_{k−j}(z),   k = 0,…,r−1,
#     h_0(0) = 1,  h_{j>0}(0) = 0,
#
# the standard log tower (1, log, …, log^{r−1}).  The h_j are the ρ-jet
# coefficients of the deformed holomorphic solution b_m(ρ): h_j(z) = Σ_m
# (∂_ρ^j b_m / j!)|_{ρ=0} z^m, computed in exact ℚ from the θ-form recursion.
#
# Transport of the period vector Π = (ϖ_0,…,ϖ_{r−1}) and the inhomogeneous
# (variation-of-parameters) solution of L·g = S(z) reuse the certified Taylor
# engine in `frobenius.jl` (`transport_step` / `transport`).
#
# RIGOR.  Within the MUM disk |z| < R (R = distance to the nearest nonzero
# singularity) the Frobenius series are summed with an EMPIRICAL geometric tail
# bound (last-window ratio, ×4 safety) — exact rational coefficients, Acb sum,
# but the tail is *not* a proved envelope.  Outside the disk the period vector
# is obtained by certified `transport` from an interior anchor; that step IS
# rigorous (frobenius.jl envelope).  Analytic continuation past a finite
# singularity requires the caller to supply a path through `period_vector(…;
# path=…)` that stays a certified distance away (`transport_step` errors if a
# leg crosses one).  This module does NOT choose monodromy paths for you.
#
# References: Bönisch–Fischbach–Klemm–Nega–Safari 2008.10574 (banana periods),
# Pögel–Wang–Weinzierl 2212.08908 (l-loop banana CY ladder).

# ---- exact ρ-jet helpers (truncated ℚ[[ρ]] of length J) ---------------------

const QJ = Vector{Rational{BigInt}}

qjet(J::Int, v::Rational{BigInt}) = (c = zeros(Rational{BigInt}, J); c[1] = v; c)
qjet(J::Int, v) = qjet(J, Rational{BigInt}(v))
qjet_var(J::Int, v) = (c = zeros(Rational{BigInt}, J); c[1] = Rational{BigInt}(v);
                       J ≥ 2 && (c[2] = Rational{BigInt}(1)); c)   # v + ρ

function qjet_mul(a::QJ, b::QJ)
    J = length(a); c = zeros(Rational{BigInt}, J)
    @inbounds for i in 1:J, j in 1:(J - i + 1)
        c[i+j-1] += a[i] * b[j]
    end
    c
end
qjet_add(a::QJ, b::QJ) = a .+ b
qjet_sub(a::QJ, b::QJ) = a .- b
qjet_scal(s, a::QJ) = Rational{BigInt}(s) .* a

function qjet_inv(a::QJ)
    J = length(a)
    a[1] != 0 || error("qjet_inv: zero constant term")
    c = zeros(Rational{BigInt}, J); c[1] = inv(a[1])
    for k in 2:J
        s = zero(Rational{BigInt})
        for j in 2:k
            s += a[j] * c[k-j+1]
        end
        c[k] = -c[1] * s
    end
    c
end

"""Evaluate polynomial Σ c_k x^k (c::Vector{Rational}) on a ρ-jet `J` (Horner)."""
function qjet_polyeval(coeffs::Vector{Rational{BigInt}}, x::QJ)
    Jl = length(x)
    r = qjet(Jl, 0)
    for c in Iterators.reverse(coeffs)
        r = qjet_add(qjet_mul(r, x), qjet(Jl, c))
    end
    r
end

# ---- Picard–Fuchs operator --------------------------------------------------

"""
    PFOperator

Monic-in-∂ order-`r` linear differential operator  L = Σ_{i=0}^{r} p_i(z) ∂^i  with
exact rational polynomial coefficients (Fuchsian at z = 0, i.e. p_i(z) = O(z^i),
so z = 0 is a regular singularity).  Carries:

  * `order`, `prec`, `coeffs[i+1][j+1] = [z^j] p_i(z) :: Rational{BigInt}`;
  * `sing` — certified Acb enclosures of the finite singular points (roots of
    p_r(z), plus 0).

Construct via `PFOperator(prec, coeffs)` (singularities are isolated from p_r),
or via `banana_pf`.
"""
struct PFOperator
    order::Int
    prec::Int
    coeffs::Vector{Vector{Rational{BigInt}}}
    sing::Vector{Acb}
end

"""Squarefree part of an exact ℚ polynomial (coeff vector, c[k+1] = [z^k])."""
function _squarefree(c::Vector{Rational{BigInt}})
    Qx, x = polynomial_ring(Nemo.QQ, "x")
    p = sum(Nemo.QQ(c[k+1]) * x^k for k in 0:length(c)-1)
    dp = Nemo.derivative(p)
    g = gcd(p, dp)
    sf = divexact(p, g)
    [Rational{BigInt}(Nemo.coeff(sf, k)) for k in 0:Nemo.degree(sf)]
end

function PFOperator(prec::Int, coeffs::Vector{<:Vector})
    r = length(coeffs) - 1
    cf = [Rational{BigInt}.(p) for p in coeffs]
    # singularities: nonzero roots of p_r (certified) ∪ {0}
    pr = copy(cf[r+1])
    k0 = 0
    while !isempty(pr) && pr[1] == 0
        popfirst!(pr); k0 += 1
    end
    isempty(pr) && error("PFOperator: leading coefficient is identically zero")
    # squarefree part — acb_poly_find_roots cannot isolate repeated roots, and
    # we only need the SET of singular points (e.g. PM K3.0 has (1−z)², CY3′
    # has (1−2¹⁰z)³).
    pr = _squarefree(pr)
    p = AcbPoly(prec = prec)
    for (j, c) in enumerate(pr)
        Arblib.set_coeff!(p, j - 1, Acb(c; prec = prec))
    end
    d = Arblib.degree(p)
    sing = Acb[Acb(0; prec)]
    if d > 0
        roots = AcbVector(d; prec = prec)
        n = ccall((:acb_poly_find_roots, Arblib.libflint), Clong,
            (Ptr{Arblib.acb_struct}, Ref{Arblib.acb_poly_struct}, Ptr{Cvoid}, Clong, Clong),
            roots, p, C_NULL, 0, prec)
        n == d || error("PFOperator: could not certify nonzero roots of leading coefficient")
        for i in 1:d
            ri = Acb(roots[i]; prec = prec)
            any(Arblib.overlaps(ri, s) for s in sing) || push!(sing, ri)
        end
    end
    PFOperator(r, prec, cf, sing)
end

"""Wrap as the certified-transport `ODEOperator` (Jeps = 1)."""
to_ode(L::PFOperator) =
    ODEOperator(L.order, 1, L.prec,
        [[[Acb(c; prec = L.prec)] for c in p] for p in L.coeffs], L.sing)

"""Evaluate the i-th coefficient polynomial p_i at an Acb point."""
function _pcoeff(L::PFOperator, i::Int, z)
    zA = Acb(z; prec = L.prec)
    s = Acb(0; prec = L.prec); zp = Acb(1; prec = L.prec)
    for c in L.coeffs[i+1]
        s += c * zp; zp *= zA
    end
    s
end

"""θ-form at z = 0:  z^{−s_min}·L = Σ_{s≥0} z^s R_s(θ).  Returns the `R_s` as
polynomials in θ (Vector of Vector{Rational{BigInt}} coefficients, R[s+1]),
where the global shift s_min = min_{p_{ij}≠0}(j−i) has been removed so that
R[1] = R_0 is the indicial polynomial.  Uses z^j ∂^i = z^{j−i}·θ(θ−1)⋯(θ−i+1)."""
function theta_form(L::PFOperator)
    r = L.order
    pairs = [(i, j - 1) for i in 0:r for j in 1:length(L.coeffs[i+1]) if L.coeffs[i+1][j] != 0]
    isempty(pairs) && error("theta_form: zero operator")
    smin = minimum(j - i for (i, j) in pairs)
    smax = maximum(j - i for (i, j) in pairs)
    R = [zeros(Rational{BigInt}, r + 1) for _ in 0:(smax - smin)]
    # falling factorial θ(θ−1)…(θ−i+1) as a polynomial in θ
    ff = Vector{Vector{Rational{BigInt}}}(undef, r + 1)
    ff[1] = Rational{BigInt}[1]                                  # i = 0
    for i in 1:r
        prev = ff[i]
        cur = zeros(Rational{BigInt}, i + 1)
        for (k, c) in enumerate(prev)                            # multiply by (θ − (i−1))
            cur[k+1] += c
            cur[k]   -= (i - 1) * c
        end
        ff[i+1] = cur
    end
    for (i, j) in pairs
        cij = L.coeffs[i+1][j+1]
        s = j - i - smin
        for (k, fk) in enumerate(ff[i+1])
            R[s+1][k] += cij * fk
        end
    end
    return R
end

# ---- MUM Frobenius basis -----------------------------------------------------

"""
    MUMFrobenius

Frobenius log tower at the MUM point z = 0.  Fields:

  * `order`  — r = n+1 (number of solutions);
  * `alpha`  — the (integer) indicial root, multiplicity r;
  * `ctx`    — `SeriesContext` with reference radius `x0 < R` (R = nearest
    nonzero singularity);
  * `h`      — `h[j+1] :: QSeries`, the j-th ρ-derivative tower h_j(z),
    j = 0,…,r−1 (h_0 = ϖ_0/z^α is the holomorphic period, h_0(0) = 1);
  * `R`      — certified lower bound on the radius of convergence.
"""
struct MUMFrobenius
    order::Int
    alpha::Int
    ctx::SeriesContext
    h::Vector{QSeries}
    R::Arb
end

"""
    mum_frobenius_basis(L, prec; nterms, x0frac = 1//2) → MUMFrobenius

Exact-ℚ Frobenius recursion at the MUM point z = 0.  Asserts the indicial
polynomial is (θ − α)^r; builds the r-jet b_m(ρ) via

    R_0(m+α+ρ)·b_m = − Σ_{s≥1} R_s(m−s+α+ρ)·b_{m−s},   b_0 = 1,

and stores h_j(z) = Σ_m [ρ^j] b_m · z^m as `QSeries` on a context with
x0 = x0frac · R.  The tail bound on each h_j is the empirical last-window
geometric envelope (max ratio over the last r+S terms, ×4 safety) — NOT a
proved bound; see module docstring.
"""
function mum_frobenius_basis(L::PFOperator, prec::Int = L.prec;
                             nterms::Int = 128, x0frac::Rational = 1//2)
    r = L.order
    Rθ = theta_form(L)
    S = length(Rθ) - 1
    # indicial polynomial = R_0(θ); must be c·(θ − α)^r
    ind = Rθ[1]
    while length(ind) > 1 && ind[end] == 0; pop!(ind); end
    length(ind) == r + 1 || error("mum_frobenius_basis: indicial degree ≠ order ($(length(ind)-1) ≠ $r)")
    lc = ind[end]
    α_rat = -ind[end-1] / (r * lc)
    denominator(α_rat) == 1 || error("mum_frobenius_basis: non-integer indicial root")
    α = Int(numerator(α_rat))
    # verify ind = lc·(θ−α)^r: jet at θ = α+ρ must be lc·ρ^r
    test = qjet_polyeval(ind, qjet_var(r + 1, α))
    (all(iszero, test[1:r]) && test[r+1] == lc) ||
        error("mum_frobenius_basis: indicial poly is not (θ−α)^$r — not MUM")
    # ρ-jet recursion
    b = Vector{QJ}(undef, nterms + 1)
    b[1] = qjet(r, 1)
    for m in 1:nterms
        acc = qjet(r, 0)
        for s in 1:min(S, m)
            x = qjet_var(r, m - s + α)                         # (m−s+α) + ρ
            acc = qjet_add(acc, qjet_mul(qjet_polyeval(Rθ[s+1], x), b[m-s+1]))
        end
        # R_0(m+α+ρ) = lc·(m+ρ)^r
        den = qjet(r, 1); mr = qjet_var(r, m)
        for _ in 1:r; den = qjet_mul(den, mr); end
        b[m+1] = qjet_scal(-1 // lc, qjet_mul(qjet_inv(den), acc))
    end
    # radius of convergence: nearest nonzero singularity
    Rconv = Arb(prec = prec); first = true
    for s in L.sing
        Arblib.contains_zero(s) && continue
        d = abs(s)
        if first; Rconv = real(d); first = false; else; Rconv = min(Rconv, real(d)); end
    end
    first && (Rconv = Arb(1; prec))   # no finite nonzero sing — unit disk fallback
    x0 = Mag(Float64(Arblib.lbound(Arb, Rconv)) * Float64(x0frac))
    Arblib.cmp(x0, Mag(1)) < 0 || (x0 = Mag(0.9))
    ctx = SeriesContext(prec, nterms, x0)
    # build QSeries h_j with empirical geometric tail
    h = Vector{QSeries}(undef, r)
    win = max(r + S, 8)
    for j in 0:r-1
        p = AcbPoly(prec = prec)
        for m in 0:nterms
            Arblib.set_coeff!(p, m, Acb(b[m+1][j+1]; prec))
        end
        # empirical tail: max_{last win} |c_{m+1}/c_m|, then geometric sum
        ratio = 0.0
        for m in max(1, nterms - win):nterms-1
            cm, cm1 = b[m+1][j+1], b[m+2][j+1]
            cm == 0 && continue
            ratio = max(ratio, Float64(abs(cm1 / cm)))
        end
        x0f = Float64(x0)
        θ = ratio * x0f
        last = abs(Float64(b[nterms+1][j+1])) * x0f^nterms
        tail = (θ < 0.95 && last > 0) ? Mag(4 * last * θ / (1 - θ)) : Mag(0)
        h[j+1] = QSeries(ctx, p, tail)
    end
    MUMFrobenius(r, α, ctx, h, Rconv)
end

"""ϖ_k as a `LogQSeries` in z (with the z^α prefactor folded in via `qshift`):
levels[j+1] = z^α · h_{k−j} / j!."""
function frobenius_logseries(frob::MUMFrobenius, k::Int)
    0 <= k < frob.order || error("frobenius_logseries: k out of range")
    levels = QSeries[]
    for j in 0:k
        fac = Acb(1 // factorial(big(j)); prec = frob.ctx.prec)
        push!(levels, fac * qshift(frob.h[k-j+1], frob.alpha))
    end
    LogQSeries(levels, Tuple{Mag,Int}[])
end

# ---- period vector / matrix --------------------------------------------------

"""∂^d (z^p · log^j z) at the point z (closed Leibniz on the log tower).
Returns an Acb.  p may be any integer."""
function _dlogmono(z::Acb, p::Int, j::Int, d::Int, prec::Int)
    # c[l+1] = coefficient of log^l z in ∂^i (z^p log^j z) / z^{p−i}
    c = zeros(Acb, j + 1); c[j+1] = Acb(1; prec)
    for i in 0:d-1
        c2 = zeros(Acb, j + 1)
        for l in 0:j
            iszero(c[l+1]) && continue
            c2[l+1] += (p - i) * c[l+1]
            l >= 1 && (c2[l] += l * c[l+1])
        end
        c = c2
    end
    L = log(z); s = Acb(0; prec); Lp = Acb(1; prec)
    for l in 0:j
        s += c[l+1] * Lp; Lp *= L
    end
    return s * z^(p - d)
end

"""Direct (in-disk) evaluation of ∂^d ϖ_k at z, by summing the truncated h-series
termwise through `_dlogmono`.  O(nterms · k) work."""
function _frob_eval(frob::MUMFrobenius, k::Int, d::Int, z)
    prec = frob.ctx.prec; α = frob.alpha; N = frob.ctx.N
    z = Acb(z; prec)
    s = Acb(0; prec)
    for j in 0:k
        hj = frob.h[k-j+1]
        invjf = Acb(1 // factorial(big(j)); prec)
        for m in 0:min(Arblib.degree(hj.coeffs), N)
            cm = Acb(Arblib.ref(hj.coeffs, m); prec)
            Arblib.is_zero(cm) && continue
            s += invjf * cm * _dlogmono(z, m + α, j, d, prec)
        end
    end
    # crude tail (same empirical envelope as the QSeries; derivative gains ≲ N^d/R^d)
    err = Mag(0)
    zr = Mag(); Arblib.get!(zr, z); ratio = Mag(); Arblib.div!(ratio, zr, frob.ctx.x0)
    rpow = Mag(); Arblib.pow!(rpow, ratio, UInt(N + 1))
    fac = Mag(Float64(big(N + 1)^d))
    for j in 0:k
        e = Mag(frob.h[k-j+1].tailF); Arblib.mul!(e, e, rpow); Arblib.mul!(e, e, fac)
        Arblib.add!(err, err, e)
    end
    Arblib.add_error!(s, err)
    s
end

"""
    period_vector(L, z, prec=L.prec; frob, anchor, path) → Vector{Acb}

(ϖ_0,…,ϖ_{r−1}) at z.  If `path === nothing` and |z| is within the Frobenius
context radius, evaluate the MUM series directly.  Otherwise transport the full
Wronskian from `anchor` along the given `path` (default path `[anchor, z]`).
`anchor` defaults to `path[1]` when a path is given, else to `x0/4`; the seed is
evaluated at `anchor` and the transport starts at `path[1]`, so an explicit
`anchor` that does not overlap `path[1]` is a hard error (seed misattribution).
Returns the 0-th-derivative row of the transported Wronskian.
"""
function period_vector(L::PFOperator, z, prec::Int = L.prec;
                       frob::Union{Nothing,MUMFrobenius} = nothing,
                       anchor = nothing, path = nothing)
    period_matrix(L, z, prec; frob = frob, anchor = anchor, path = path)[1, :]
end

"""
    period_matrix(L, z, prec=L.prec; frob, anchor, path) → Matrix{Acb} (r×r)

Wronskian W with W[i+1, k+1] = ϖ_k^{(i)}(z), i,k = 0,…,r−1.  Same evaluation /
transport logic as `period_vector`.
"""
function period_matrix(L::PFOperator, z, prec::Int = L.prec;
                       frob::Union{Nothing,MUMFrobenius} = nothing,
                       anchor = nothing, path = nothing)
    r = L.order
    frob === nothing && (frob = mum_frobenius_basis(L, prec))
    zA = Acb(z; prec)
    zr = Mag(); Arblib.get!(zr, zA)
    in_disk = Arblib.cmp(zr, frob.ctx.x0) < 0
    if path === nothing && in_disk
        W = Matrix{Acb}(undef, r, r)
        for k in 0:r-1, i in 0:r-1
            W[i+1, k+1] = _frob_eval(frob, k, i, zA)
        end
        return W
    end
    # transport each Frobenius column from `anchor` to z along `path`.
    # The Frobenius seed is EVALUATED at `anchor` but the transport STARTS at
    # `path[1]` — these must be the same point.  Pitfall (found in a
    # production run): with an explicit path and the old unconditional default
    # anchor = x0/4, the seed was silently misattributed; the continuation
    # comes back a constant basis conjugation off (first digit wrong) yet
    # passes every internal monodromy diagnostic.  Now: when a path is given
    # and no anchor, the anchor IS path[1]; an explicit anchor that does not
    # overlap path[1] is a hard error.
    pth = path === nothing ? nothing : Acb[Acb(p; prec) for p in path]
    aA = anchor === nothing ?
         (pth === nothing ? Acb(Float64(frob.ctx.x0) / 4; prec) : pth[1]) :
         Acb(anchor; prec)
    pth === nothing && (pth = Acb[aA, zA])
    Arblib.overlaps(aA, pth[1]) ||
        error("period_matrix: anchor ≠ path[1] — the Frobenius seed would be " *
              "silently misattributed (constant basis conjugation; wrong values " *
              "that pass internal monodromy diagnostics). Pass anchor = path[1], " *
              "or omit `anchor` to default it to path[1].")
    W = Matrix{Acb}(undef, r, r)
    for k in 0:r-1
        Y0d = [_frob_eval(frob, k, i, aA) for i in 0:r-1]
        Yθ = _deriv_to_theta(aA, Y0d)
        Ytθ, _ = _theta_transport(L, pth, Yθ)
        Ytd = _theta_to_deriv(zA, Ytθ)
        for i in 0:r-1
            W[i+1, k+1] = Ytd[i+1]
        end
    end
    return W
end

"""Certified-transport variant (frobenius.jl `transport`).  Requires the path to
stay at |z| ≳ R/2 (away from the MUM zero of p_r) — see the module note on the
envelope bracket near regular singularities."""
function period_matrix_certified(L::PFOperator, z, prec::Int = L.prec;
                                 frob::Union{Nothing,MUMFrobenius} = nothing,
                                 anchor, path)
    r = L.order
    frob === nothing && (frob = mum_frobenius_basis(L, prec))
    op = to_ode(L)
    aA = Acb(anchor; prec)
    pth = Acb[Acb(p; prec) for p in path]
    Arblib.overlaps(aA, pth[1]) ||
        error("period_matrix_certified: anchor ≠ path[1] — the Frobenius seed " *
              "would be silently misattributed. Pass anchor = path[1].")
    W = Matrix{Acb}(undef, r, r)
    for k in 0:r-1
        Y0 = [[_frob_eval(frob, k, i, aA)] for i in 0:r-1]
        Yt, _ = transport(op, pth, Y0)
        for i in 0:r-1
            W[i+1, k+1] = Yt[i+1][1]
        end
    end
    return W
end

# ---- θ-companion Taylor transport (non-certified analytic continuation) -----
#
# The certified `transport` engine in frobenius.jl normalises by p_r(x₀)⁻¹ and
# verifies a sufficient envelope bracket ≤ 1.  Near the MUM point p_r(z) ~ z^r
# so p_r'(z)/p_r(z) ~ r/z and the (r, j≥1) bracket terms sum to ≈ (1+θ)^r − 1,
# which exceeds 1 for any θ > 2^{1/r} − 1 once r ≥ 4 — the trial-θ floor 0.18 is
# too high.  The Taylor series itself is fine; only the *bound* is pessimistic.
#
# For analytic continuation we therefore step the θ-companion system
#     θ·Y = B(z)·Y,   Y = (g, θg, …, θ^{r−1}g),   B = companion of −R_s/R_r,
# directly in z by high-order Taylor (Cauchy product on the Taylor expansion of
# B(z)/z about each step point), exactly as in period_transport.py.  This is
# ball arithmetic throughout but carries NO rigorous tail; convergence is
# checked by step-doubling.  Use `method = :certified` in `period_matrix` to
# force the frobenius.jl path (requires the path to stay at |z| ≳ R/2).

"""High-order Taylor step of d/dz Y = (B(z)/z)·Y from z₀ to z₀+h.  B(z) is the
θ-companion matrix with entries rational in z (polynomials over the common
leading θ-coefficient).  `ord` Taylor terms; ball arithmetic, no certified tail."""
function _theta_taylor_step(L::PFOperator, Rθ, z0::Acb, h::Acb, Y::Vector{Acb}; ord::Int = 40)
    r = L.order; prec = L.prec
    S = length(Rθ) - 1
    # B(z)_{ij}: rows 1..r-1 are sub-diagonal shift; last row = −R_s(θ)→ −Σ_s z^s R_s[j]/R_lead.
    # We need Taylor coefficients of A(z) = B(z)/z about z₀ to order `ord`.
    # Build numerator/denominator polynomials in z for each last-row entry; the
    # leading θ^r coefficient is q(z) = Σ_s R_s[r+1] z^s.
    qpoly = AcbPoly(prec = prec)
    for s in 0:S; Arblib.set_coeff!(qpoly, s, Acb(Rθ[s+1][r+1]; prec)); end
    # Taylor of 1/(z·q(z)) about z₀: shift then invert as series.
    function taylor_ratinv(p::AcbPoly, n::Int)      # Taylor of 1/p(z₀+u) to u^n
        sh = AcbPoly(prec = prec); Arblib.taylor_shift!(sh, p, z0)
        out = AcbPoly(prec = prec); Arblib.inv_series!(out, sh, n + 1)
        return out
    end
    zpoly = AcbPoly([Acb(0; prec), Acb(1; prec)]; prec = prec)
    zq = AcbPoly(prec = prec); Arblib.mul!(zq, zpoly, qpoly)
    invzq = taylor_ratinv(zq, ord)
    # last-row numerators: N_j(z) = Σ_s R_s[j+1] z^s, j = 0..r-1
    Nj = Vector{AcbPoly}(undef, r)
    for j in 0:r-1
        p = AcbPoly(prec = prec)
        for s in 0:S; Arblib.set_coeff!(p, s, Acb(Rθ[s+1][j+1]; prec)); end
        sh = AcbPoly(prec = prec); Arblib.taylor_shift!(sh, p, z0)
        Nj[j+1] = AcbPoly(prec = prec); Arblib.mullow!(Nj[j+1], sh, invzq, ord + 1)
    end
    # 1/z Taylor (for the shift rows)
    invz = taylor_ratinv(zpoly, ord)
    # Taylor recursion: a[0]=Y, a[n+1] = (1/(n+1)) Σ_{k=0}^{n} A_k · a[n−k]
    a = Vector{Vector{Acb}}(undef, ord + 1); a[1] = copy(Y)
    for n in 0:ord-1
        nxt = [Acb(0; prec) for _ in 1:r]
        for k in 0:n
            ak = a[n-k+1]
            invzk = Acb(Arblib.ref(invz, k); prec)
            for i in 1:r-1
                nxt[i] += invzk * ak[i+1]
            end
            for j in 0:r-1
                nxt[r] -= Acb(Arblib.ref(Nj[j+1], k); prec) * ak[j+1]
            end
        end
        a[n+2] = [x / (n + 1) for x in nxt]
    end
    out = [Acb(0; prec) for _ in 1:r]; hp = Acb(1; prec)
    for n in 0:ord
        for i in 1:r; out[i] += a[n+1][i] * hp; end
        hp *= h
    end
    out
end

"""Non-certified θ-companion transport along `path`.  Step size is `frac × dist`
where dist = min(|z|, |z − s| for s ∈ sing\\{0}); with `ord` Taylor terms the
truncation error per step is ≲ frac^ord."""
function _theta_transport(L::PFOperator, path::Vector{Acb}, Y::Vector{Acb};
                          ord::Int = 0, frac::Float64 = 0.25)
    prec = L.prec
    ord == 0 && (ord = ceil(Int, prec * log(2.0) / log(1 / frac)) + 8)
    Rθ = theta_form(L)
    z = path[1]
    sing = [s for s in L.sing if !Arblib.contains_zero(s)]
    nstep = 0
    for target in path[2:end]
        while true
            δ = target - z
            d = Float64(abs(δ))
            dist = Float64(Arblib.lbound(Arb, abs(z)))
            for s in sing
                dist = min(dist, Float64(Arblib.lbound(Arb, abs(s - z))))
            end
            nstep += 1
            nstep < 100000 || error("_theta_transport: step limit exceeded")
            if d <= frac * dist
                Y = _theta_taylor_step(L, Rθ, z, δ, Y; ord = ord)
                z = target; break
            end
            h = δ * Acb(frac * dist / d; prec)
            Y = _theta_taylor_step(L, Rθ, z, h, Y; ord = ord)
            z = z + h
        end
    end
    return Y, z
end

_stirling2(n, k) = (k == 0 ? big(n == 0 ? 1 : 0) :
                    n == 0 ? big(0) : k * _stirling2(n - 1, k) + _stirling2(n - 1, k - 1))
_stirling1s(n, k) = (k == 0 ? big(n == 0 ? 1 : 0) :       # signed s(n,k)
                     n == 0 ? big(0) : _stirling1s(n - 1, k - 1) - (n - 1) * _stirling1s(n - 1, k))

"""Convert (g, g', …, g^{(r−1)}) ↔ (g, θg, …, θ^{r−1}g) at z (Stirling-2nd-kind)."""
function _deriv_to_theta(z::Acb, Y::Vector{Acb})
    r = length(Y); prec = Arblib.precision(z)
    out = [Acb(0; prec) for _ in 1:r]; out[1] = Y[1]
    for k in 1:r-1
        s = Acb(0; prec)
        for d in 1:k
            s += Acb(_stirling2(k, d); prec) * z^d * Y[d+1]
        end
        out[k+1] = s
    end
    out
end
function _theta_to_deriv(z::Acb, Yθ::Vector{Acb})
    r = length(Yθ); prec = Arblib.precision(z)
    out = [Acb(0; prec) for _ in 1:r]; out[1] = Yθ[1]
    for k in 1:r-1
        s = Acb(0; prec)
        for j in 1:k
            s += Acb(_stirling1s(k, j); prec) * Yθ[j+1]
        end
        out[k+1] = s / z^k
    end
    out
end

# ---- mirror map --------------------------------------------------------------

"""
    mirror_map(L_or_frob, z, prec) → (q, logq)

CY mirror coordinate q = exp(ϖ_1/ϖ_0) and its log at z (MUM-disk evaluation).
For n = 1 this is the elliptic nome up to a multiplicative constant fixed by the
leading Hauptmodul coefficient (q_mirror/z → 1 by construction); for the K3
banana it equals the underlying sunrise mirror map (Sym² identity).
"""
mirror_map(L::PFOperator, z, prec::Int = L.prec) =
    mirror_map(mum_frobenius_basis(L, prec), z, prec)

function mirror_map(frob::MUMFrobenius, z, prec::Int = frob.ctx.prec)
    zA = Acb(z; prec)
    w0 = _frob_eval(frob, 0, 0, zA)
    w1 = _frob_eval(frob, 1, 0, zA)
    L = w1 / w0
    return (exp(L), L)
end

# ---- CY sector descriptor (parallel to EllipticSector) -----------------------

"""
    CYSector

Descriptor for one CY (order-r) sector.  Carries the `PFOperator`, its
`MUMFrobenius` tower, an interior anchor point with cached Wronskian, and the
wrapped certified `ODEOperator`.  Construct via `CYSector(L; anchor, nterms)`.
"""
struct CYSector
    L::PFOperator
    frob::MUMFrobenius
    op::ODEOperator
    anchor::Acb
    W_anchor::Matrix{Acb}
end

function CYSector(L::PFOperator; anchor = nothing, nterms::Int = 160)
    frob = mum_frobenius_basis(L, L.prec; nterms = nterms)
    a = anchor === nothing ? Acb(Float64(frob.ctx.x0) / 4; prec = L.prec) :
                             Acb(anchor; prec = L.prec)
    W = period_matrix(L, a, L.prec; frob = frob)
    CYSector(L, frob, to_ode(L), a, W)
end

period_matrix(sec::CYSector, z; path = nothing) =
    period_matrix(sec.L, z, sec.L.prec; frob = sec.frob, anchor = sec.anchor, path = path)
period_vector(sec::CYSector, z; path = nothing) =
    period_matrix(sec, z; path = path)[1, :]
mirror_map(sec::CYSector, z) = mirror_map(sec.frob, z, sec.L.prec)

# ---- variation-of-parameters transport --------------------------------------

# Small certified Acb-matrix solve (Gaussian elimination with partial pivoting on
# midpoints; ball arithmetic throughout).
function _acb_solve(A::Matrix{Acb}, b::Vector{Acb})
    n = length(b); prec = Arblib.precision(b[1])
    M = [Acb(A[i,j]; prec) for i in 1:n, j in 1:n]
    v = [Acb(b[i]; prec) for i in 1:n]
    mg = Mag()
    for k in 1:n
        piv = k; best = -1.0
        for i in k:n
            Arblib.get!(mg, M[i, k]); m = Float64(mg)
            m > best && (best = m; piv = i)
        end
        if piv != k
            for j in 1:n
                M[k,j], M[piv,j] = M[piv,j], M[k,j]
            end
            v[k], v[piv] = v[piv], v[k]
        end
        Arblib.contains_zero(M[k,k]) && error("_acb_solve: singular (ball) pivot")
        inv_p = inv(M[k,k])
        for i in 1:n
            i == k && continue
            f = M[i,k] * inv_p
            for j in k:n; M[i,j] -= f * M[k,j]; end
            v[i] -= f * v[k]
        end
    end
    [v[k] / M[k,k] for k in 1:n]
end

"""
    cy_transport(sec, source, z_start, z_end, prec=sec.L.prec;
                 cusp_reg=true, method=:certified, n=96)
        → (g::Vector{Acb}, I::Vector{Acb})

CY analogue of `eichler_transport`: the particular solution of L·g = source(z)
with zero initial data at `z_start`, evaluated at `z_end`, via variation of
parameters

    g_p^{(i)}(z) = Σ_k W_{ik}(z) · I_k,
    I_k = ∫_{z_start}^{z_end} [W^{-1}(s)]_{k,r} · source(s) / p_r(s) ds.

Returns `(g = (g_p, g_p', …, g_p^{(r−1)}), I = (I_0,…,I_{r−1}))`.  The integrals
are evaluated by `Arblib.integrate` (`method=:certified`, rigorous) or by
Gauss–Legendre (`method=:gauss_legendre`, non-certified cross-check).  The
period matrix W(s) at each quadrature node is the in-disk Frobenius evaluation
when |s| < x0, else certified transport from `sec.anchor`.

Reduction to n = 1.  For order 2, det W = Wronskian and [W⁻¹]_{1,2} = −ϖ_1/det W,
[W⁻¹]_{2,2} = ϖ_0/det W, so I = (−∫ ϖ_1 f, ∫ ϖ_0 f)/det W and g_p = ϖ_0 I_0 +
ϖ_1 I_1 — the textbook order-2 VoP that the elliptic `eichler_transport` packages
in the τ variable (after the change dτ = (det W/ϖ_0²) dz).  The two routes are
checked to agree in `test_cy_transport.jl`.

LIMITATION (`cusp_reg`).  The MUM-regularised lower limit z_start = 0 is handled
by replacing 0 → ε with ε = 2^{−prec/r}: the VoP integrand near 0 behaves like
z^{α + r·something − r}·log^{r−1} z (integrable), and the omitted [0, ε] piece is
exponentially small in `prec`; this is the practical regularisation, NOT a proved
tangential-basepoint formula.  For a rigorous cusp constant use the Frobenius
series of g_p directly.
"""
function cy_transport(sec::CYSector, source::Function, z_start, z_end,
                      prec::Int = sec.L.prec; cusp_reg::Bool = true,
                      method::Symbol = :certified, n::Int = 96)
    r = sec.L.order
    za = Acb(z_start; prec); zb = Acb(z_end; prec)
    if cusp_reg && Arblib.contains_zero(za)
        ε = Acb(2; prec)^(-(prec ÷ max(r, 2)))
        za = ε * (Arblib.contains_zero(zb) ? Acb(1; prec) : zb / abs(zb))
    end
    # integrand: k-th component of W^{-1}(s) · e_r · source(s)/p_r(s)
    function col(s::Acb)
        W = period_matrix(sec, s)
        rhs = [i == r ? source(s) / _pcoeff(sec.L, r, s) : Acb(0; prec) for i in 1:r]
        _acb_solve(W, rhs)
    end
    I = Vector{Acb}(undef, r)
    if method === :certified
        for k in 1:r
            f = (s; analytic::Bool = false) -> begin
                if analytic
                    try
                        return col(s)[k]
                    catch
                        ind = Acb(prec = prec); Arblib.indeterminate!(ind); return ind
                    end
                end
                col(s)[k]
            end
            I[k] = Arblib.integrate(f, za, zb; check_analytic = true, prec = prec,
                                    rtol = exp10(-0.301 * prec))
        end
    elseif method === :gauss_legendre
        xs, ws = _gl_nodes(n, prec)
        half = (zb - za) / 2; mid = (za + zb) / 2
        for k in 1:r; I[k] = Acb(0; prec); end
        for j in 1:n
            s = mid + half * Acb(xs[j]; prec)
            c = col(s)
            for k in 1:r; I[k] += Acb(ws[j]; prec) * c[k]; end
        end
        for k in 1:r; I[k] *= half; end
    else
        error("cy_transport: unknown method $method")
    end
    Wend = period_matrix(sec, zb)
    g = [sum(Wend[i, k] * I[k] for k in 1:r) for i in 1:r]
    return (g, I)
end

# Convenience overload mirroring `eichler_transport(L, source, …)`.
cy_transport(L::PFOperator, source::Function, z_start, z_end, prec::Int = L.prec;
             kwargs...) = cy_transport(CYSector(L), source, z_start, z_end, prec; kwargs...)

# ---- equal-mass banana operators --------------------------------------------

"""
    banana_pf(nloops; prec=256, masses=:equal) → PFOperator

The equal-mass l-loop banana Picard–Fuchs operator (l = nloops) in its standard
literature variable, MUM at z = 0:

  * l = 2 (sunrise, elliptic, r = 2): Laporta–Remiddi t = p²/m², singularities
    {0, 1, 9}.  L = t(t−1)(t−9)∂² + (3t²−20t+9)∂ + (t−3); indicial θ².
    ϖ_0 = 1 + t/3 + 5t²/27 + 31t³/243 + …
  * l = 3 (K3, r = 3): Bönisch et al. variable z (sing {0, 1/16, 1/4}); ϖ_0 = Domb
    series Σ A002895(k) z^k.  Coefficients from k3_banana_geometry.py STEP 2A.
  * l = 4 (CY3, r = 4): PWW 2212.08908 eq. 26 / AESZ #34, variable y
    (sing {0, −1/25, −1/9, −1}); indicial (θ−1)⁴, ψ_0 = y(1 − 5y + 45y² − …).
    Polynomial coefficients are p_i = y⁴(1+y)(1+9y)(1+25y)·r_i.

Unequal masses (`masses ≠ :equal`) are not encoded — use `PFOperator(prec, coeffs)`
directly with the operator from `build_pf_uneq.py`.
"""
function banana_pf(nloops::Int; prec::Int = 256, masses = :equal)
    masses === :equal || error("banana_pf: only :equal masses are encoded; " *
        "build the unequal-mass operator externally and call PFOperator(prec, coeffs)")
    if nloops == 2
        return PFOperator(prec, [[-3, 1], [9, -20, 3], [0, 9, -10, 1]])
    elseif nloops == 3
        return PFOperator(prec, [[0, -4, 64], [0, 1, -68, 448],
                                 [0, 0, 3, -90, 384], [0, 0, 0, 1, -20, 64]])
    elseif nloops == 4
        # p_k = y^4 (1+y)(1+9y)(1+25y) · r_k(y),  r_k from PWW eq. 26.
        # Expanded exactly (degree-7 polynomials in y):
        p4 = [0, 0, 0, 0, 1, 35, 259, 225]
        p3 = [0, 0, 0, 2, 140, 1554, 1800, 0]
        p2 = [0, 0, 1, 98, 1839, 3150, 0, 0]
        p1 = [0, -1, 0, 285, 900, 0, 0, 0]
        p0 = [1, 5, 0, 0, 0, 0, 0, 0]
        return PFOperator(prec, [p0, p1, p2, p3, p4])
    else
        error("banana_pf: nloops = $nloops not encoded (have 2,3,4)")
    end
end

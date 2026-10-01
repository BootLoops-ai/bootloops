# Core G(a₁,…,aₙ; z) evaluator.
#
# Algorithm (Vollinga–Weinzierl hep-ph/0410259, plus a rigorous tail bound):
#   0. all aᵢ = 0                → logⁿz/n!.
#   1. trailing zeros            → shuffle-regularise to words ending in aₙ≠0.
#   2. ∃ aᵢ on the open segment  → analytic_continue: Chen path-concatenation
#      (0,z)                       through a complex midpoint zₘ (VW §5 /
#                                   Frellesvig–Tommasini–Wever 1601.02649 §4.2).
#                                   The ±iπ monodromy pieces emerge automatically
#                                   from the principal logs along the deformed
#                                   path; the side is fixed by `prescription`.
#   3. depth-1 / Li closed form  → Arblib log / polylog!.
#   4. |z| ≤ ρ_thresh · r        → direct power series at 0, certified tail.
#   5. otherwise                 → rescale G(a;z)=G(a/z;1), Hölder-convolute at
#                                   p=2 to G(⋯;1/2) pieces, recurse on each.
#
# Step 2 handles the on-path case natively — a real index w ∈ (0,z) yields
# the regularised G (matching GiNaC's 2-arg G default), never the silent
# Cauchy principal value that naive real-axis integration produces there.
#
# Series-core tail bound (rigorous): for any word w with no trailing zero and
# r := min_{aᵢ≠0}|aᵢ|, the Taylor coefficient of G(w;t) at 0 satisfies
# |c_m| ≤ r⁻ᵐ (proof by induction on the O(N) recurrence below).  Hence
#   |G(w;z) − Σ_{m≤N} c_m zᵐ| ≤ (|z|/r)^{N+1}/(1−|z|/r).
#
# A small per-call cache avoids recomputing shared sub-words / Hölder pieces.

const _ACBVEC = Vector{Acb}

"Memoisation cache keyed on a string fingerprint of (a⃗, z) at fixed prec."
struct GPLCache
    prec::Int
    table::Dict{String,Acb}
end
GPLCache(prec::Int) = GPLCache(prec, Dict{String,Acb}())

_fp(x::Arb) = string(x; digits = 32)
_cache_key(a::_ACBVEC, z::Acb, presc::Int) =
    join((_fp(real(ai)) * "|" * _fp(imag(ai)) for ai in a), ";") * "@" *
    _fp(real(z)) * "|" * _fp(imag(z)) * "#" * string(presc)

"Float64 upper bound of an Arb."
_ubf(x::Arb) = Arblib.get_d(Arblib.ubound(Arf, x), RoundUp)
"Float64 lower bound of an Arb."
_lbf(x::Arb) = Arblib.get_d(Arblib.lbound(Arf, x), RoundDown)

# ---------------------------------------------------------------------------
# Closed forms / primitives
# ---------------------------------------------------------------------------

"Liₙ(z) via Arblib (certified ball)."
function li(n::Integer, z; prec::Int)
    w = Acb(prec = prec)
    Arblib.polylog!(w, Int(n), Acb(z, prec = prec))
    return w
end

_ac(x, prec) = Acb(x, prec = prec)

# G(0ᵏ; z) = logᵏz/k!
_gpl_allzero(k::Int, z::Acb, prec::Int) =
    k == 0 ? _ac(1, prec) : Arblib.log!(Acb(prec = prec), z)^k / factorial(big(k))

# ---------------------------------------------------------------------------
# Trailing-zero shuffle regularisation
# ---------------------------------------------------------------------------
# G(0;z)·G(w;z) = Σ_{σ∈ 0⧢w} G(σ;z); the term where 0 is appended LAST is the
# trailing-zero word.  Iterating gives a polynomial in log z with coefficients
# that are G's of words with one fewer trailing zero. Recursion terminates.

_trailing_zeros(a::_ACBVEC) = begin
    k = 0
    @inbounds while k < length(a) && Arblib.contains_zero(a[end-k])
        iszero(a[end-k]) || error("gpl: letter $(a[end-k]) overlaps 0; supply exact 0 or move it off the origin")
        k += 1
    end
    k
end

# Reduce one trailing zero. From  G(0;z)·G(w;z) = Σ_{σ ∈ (0)⊔w} G(σ;z)  with
# w = (a₁,…,aₘ,0^{k−1}) = w0[1:end-1], the σ's that equal w0 = (a,0^k) occur
# with multiplicity k (insert before any of the k−1 trailing zeros, or append),
# hence
#     k · G(w0;z) = log(z)·G(w;z) − Σ_{pos=1}^{m} G(a₁..a_{pos−1},0,a_pos..aₘ,0^{k−1}; z),
# every RHS word having ≤ k−1 trailing zeros (recursion terminates in k steps).
function _reduce_one_trailing_zero(w0::_ACBVEC, z::Acb, prec::Int, cache, depth::Int,
                                   presc::Int)
    k = _trailing_zeros(w0)
    w = w0[1:end-1]                      # drop ONE trailing 0 (k−1 remain)
    m = length(w0) - k                   # length of the nonzero prefix a₁..aₘ
    L = Arblib.log!(Acb(prec = prec), z)
    s = L * _gpl(w, z, prec, cache, depth, presc)
    z0 = Acb(0, prec = prec)
    for pos in 1:m
        v = vcat(w[1:pos-1], [z0], w[pos:end])
        s -= _gpl(v, z, prec, cache, depth, presc)
    end
    return s / k
end

# ---------------------------------------------------------------------------
# Certified power series at 0 (no trailing zero, |z| < r)
# ---------------------------------------------------------------------------
# Recurrence: with G(a, w; t) = Σ gₘ tᵐ, G(w; t) = Σ cₘ tᵐ,
#   a = 0 :  gₘ = cₘ/m,
#   a ≠ 0 :  gₘ = ((m−1) g_{m−1} − c_{m−1})/(a m),  g₀ = 0.
# Innermost (aₙ ≠ 0):  cₘ = −1/(m aₙᵐ).

function _series_coeffs(a::_ACBVEC, N::Int, prec::Int)
    n = length(a)
    @assert n ≥ 1 && !iszero(a[end])
    c = Vector{Acb}(undef, N + 1)
    c[1] = _ac(0, prec)
    inva = inv(a[end]); p = _ac(1, prec)
    for m in 1:N
        p *= inva                                  # = aₙ^{-m}
        c[m+1] = -p / m
    end
    for k in (n-1):-1:1
        ak = a[k]
        if iszero(ak)
            for m in 1:N
                c[m+1] = c[m+1] / m
            end
        else
            g = Vector{Acb}(undef, N + 1)
            g[1] = _ac(0, prec)
            invak = inv(ak)
            for m in 1:N
                g[m+1] = ((m - 1) * g[m] - c[m]) * invak / m
            end
            c = g
        end
    end
    return c
end

"Minimum |aᵢ| over nonzero letters (Arb lower bound)."
function _rmin(a::_ACBVEC, prec::Int)
    r = Arb(Inf, prec = prec)
    for ai in a
        iszero(ai) && continue
        r = min(r, abs(ai))
    end
    return r
end

function _gpl_series(a::_ACBVEC, z::Acb, prec::Int)
    r   = _rmin(a, prec)
    az  = abs(z)
    rho = Arb(Arblib.ubound(Arf, az / r), prec = 64)
    rho < 1 || error("series: |z|/r = $rho ≥ 1")
    # number of terms for prec bits + safety
    digits = ceil(Int, prec * log10(2.0))
    N = ceil(Int, (digits + 20) * log(10.0) / log(1.0 / Float64(rho))) + 10
    N = clamp(N, 8, 200_000)
    c = _series_coeffs(a, N, prec)
    # Horner
    s = _ac(0, prec)
    for m in N:-1:0
        Arblib.mul!(s, s, z)
        Arblib.add!(s, s, c[m+1])
    end
    # rigorous tail: |c_m| ≤ r^{-m}  ⇒  tail ≤ (|z|/r)^{N+1}/(1−|z|/r)
    tail = rho^(N + 1) / (1 - rho)
    Arblib.add_error!(s, tail)
    return s
end

# ---------------------------------------------------------------------------
# Hölder convolution (Vollinga–Weinzierl eq. 34)
# ---------------------------------------------------------------------------
# For zₙ ≠ 0 (no trailing zero) and z₁ ≠ 1 (convergent at upper endpoint):
#   G(z₁,…,zₙ; 1) = Σ_{j=0}^{n} (−1)ʲ G(1−z_j,…,1−z₁; 1−1/p) · G(z_{j+1},…,zₙ; 1/p).
# Each factor is a length-≤n GPL at a smaller argument; recurse via _gpl.
# p is chosen so that 1/p and 1−1/p stay well clear of every zᵢ (otherwise an
# individual factor would be endpoint-divergent even though the sum is finite).

const _HOELDER_PS = (1//2, 1//3, 2//3, 1//4, 3//4, 2//5, 3//5, 1//5, 4//5,
                     3//7, 4//7, 5//11, 6//11, 7//13, 8//17)

"Pick 1/p ∈ (0,1) maximising min_i |zᵢ − 1/p| (so no Hölder factor is endpoint-singular)."
function _pick_hoelder_q(b::_ACBVEC, prec::Int)
    best_q, best_d = _ac(1, prec) / 2, -1.0
    for q in _HOELDER_PS
        qc = _ac(q, prec)
        d = minimum((_lbf(abs(bi - qc)) for bi in b if !iszero(bi)); init = Inf)
        if q == 1//2 && d > 0.15
            return qc                         # prefer p=2 when safely clear
        end
        if d > best_d
            best_q, best_d = qc, d
        end
    end
    return best_q
end

const _HOELDER_MAXDEPTH = 12

function _hoelder(b::_ACBVEC, prec::Int, cache, depth::Int)
    n = length(b)
    one = _ac(1, prec)
    # Endpoint divergence at z=1 ↔ b₁=1: not (yet) regularised.
    if n ≥ 1 && Arblib.contains_zero(b[1] - one)
        error("gpl: leading letter equals z (a₁=z); G is log-divergent at the upper endpoint — supply the shuffle-regularised form")
    end
    q  = _pick_hoelder_q(b, prec)             # = 1/p
    q1 = one - q                              # = 1 − 1/p
    s = _ac(0, prec)
    for j in 0:n
        left  = Acb[one - b[k] for k in j:-1:1]      # reversed, complemented
        right = b[j+1:end]
        gl = isempty(left)  ? _ac(1, prec) : _gpl(left,  q1, prec, cache, depth + 1)
        gr = isempty(right) ? _ac(1, prec) : _gpl(right, q,  prec, cache, depth + 1)
        s += (isodd(j) ? -1 : 1) * gl * gr
    end
    return s
end

# ---------------------------------------------------------------------------
# Analytic continuation across on-path letters (VW §5; FTW 1601.02649 §4.2)
# ---------------------------------------------------------------------------
# When some aₖ lies on the open straight-line segment (0,z), the defining
# integral has a 1/(t−aₖ) pole on the contour. The value of G is fixed by an
# i0-prescription on the letter: aₖ → aₖ + i·s·0 with s = ±1. We realise this
# by deforming the contour to the OPPOSITE side (Chen path-concatenation
# through a complex midpoint zₘ), so the pole is avoided and every sub-G is an
# ordinary convergent GPL handled by steps 3–5 above. The ±iπ pieces are then
# exact (no finite shift, no precision loss).
#
# Prescription: `:ginac` / `:plus_i0` (DEFAULT, s=+1) matches GiNaC's two-arg
# G(lst{a},z) — e.g. G(1/2;1) = log(−1) = +iπ. `:feynman` / `:minus_i0` (s=−1)
# is the Feynman aₖ−i0 side, giving −iπ for the same example.

_presc_sign(p::Symbol) =
    p === :ginac    || p === :plus_i0  ?  1 :
    p === :feynman  || p === :minus_i0 ? -1 :
    error("gpl: unknown prescription :$p (use :ginac/:plus_i0 or :feynman/:minus_i0)")

"""
    on_path_indices(a, z) -> Vector{Int}

Indices `i` for which `a[i]` lies on the open segment (0,z), i.e. `a[i]/z` is
real (its imaginary-part ball contains zero) with real part certainly in (0,1).
Exact zeros are skipped (they are the regularised dlog at the base point).
"""
function on_path_indices(a::_ACBVEC, z::Acb)
    idx = Int[]
    Arblib.contains_zero(z) && return idx
    invz = inv(z)
    for (i, ai) in enumerate(a)
        iszero(ai) && continue
        r  = ai * invz
        rr = real(r)
        if Arblib.contains_zero(imag(r)) && _lbf(rr) > 0 && _ubf(rr) < 1
            push!(idx, i)
        end
    end
    idx
end

# Euclidean distance from point w to the closed segment [p,q] (Arb).
function _dist_to_seg(w::Acb, p::Acb, q::Acb)
    d  = q - p
    L2 = real(d * conj(d))
    t  = real((w - p) * conj(d)) / L2
    pr = precision(real(t))
    tc = min(max(t, Arb(0, prec = pr)), Arb(1, prec = pr))
    abs(w - (p + tc * d))
end

const _MIDPOINT_HS = (1//2, 1//3, 2//3, 1//4, 3//4, 2//5, 3//5,
                      5//11, 6//11, 7//13, 3//17, 8//19)

# Pick zₘ = 1/2 − i·s·h such that no rescaled letter lies on [0,zₘ] or [zₘ,1].
function _path_midpoint(b::_ACBVEC, prec::Int, presc::Int)
    one1 = _ac(1, prec); zero1 = _ac(0, prec)
    best_zm, best_clr = one1 / 2, -1.0
    for h in _MIDPOINT_HS
        zm = Acb(1//2, -presc * h, prec = prec)
        clr = Inf
        for bi in b
            iszero(bi) && continue
            clr = min(clr, _lbf(_dist_to_seg(bi, zero1, zm)),
                           _lbf(_dist_to_seg(bi, zm, one1)))
        end
        # prefer h=1/2 when safely clear (keeps Hölder pieces well-conditioned)
        if h == 1//2 && clr > 0.08
            return zm
        end
        if clr > best_clr
            best_zm, best_clr = zm, clr
        end
    end
    best_clr > 0 || error("analytic_continue: no clear deformed midpoint found " *
                          "(letters tile both half-planes near (0,z)); pass " *
                          "letters with explicit small imaginary parts instead")
    best_zm
end

"""
    analytic_continue(a, z, prec; prescription = :ginac, cache = nothing) -> Acb

G(a₁,…,aₙ; z) when one or more aₖ lie on the open segment (0,z). Implements the
Vollinga–Weinzierl §5 reduction via Chen path-concatenation through a complex
midpoint: with b = a/z and zₘ = ½ − i·s·h (s = prescription sign),

    G(b; 1) = Σ_{k=0}^{n} G(b₁−zₘ,…,b_k−zₘ; 1−zₘ) · G(b_{k+1},…,b_n; zₘ),

each factor being an ordinary off-path GPL. The ±iπ monodromy is exact.

`method = :hybrid` instead applies the outermost-IBP / innermost-closed-form
trick (log kernels in place of pole kernels) and integrates the depth-(n−2)
middle with `Arblib.integrate` along the same deformed contour — a rigorous,
algorithmically independent fallback (also certified).
"""
function analytic_continue(a::_ACBVEC, z::Acb, prec::Int;
                           prescription::Symbol = :ginac,
                           method::Symbol = :path,
                           cache::Union{Nothing,GPLCache} = nothing)
    s = _presc_sign(prescription)
    method === :path   && return _analytic_continue(a, z, prec, cache, 0, s)
    method === :hybrid && return _ac_hybrid(a, z, prec, cache, 0, s)
    error("analytic_continue: unknown method :$method (use :path or :hybrid)")
end
analytic_continue(a::AbstractVector, z; prec::Int, kw...) =
    analytic_continue(Acb[iszero(x) ? Acb(0, prec = prec) : Acb(x, prec = prec) for x in a],
                      Acb(z, prec = prec), prec; kw...)

function _analytic_continue(a::_ACBVEC, z::Acb, prec::Int, cache, depth::Int,
                            presc::Int)
    n = length(a)
    invz = inv(z)
    b = Acb[iszero(ai) ? _ac(0, prec) : ai * invz for ai in a]
    one1 = _ac(1, prec)
    zm = _path_midpoint(b, prec, presc)
    dz = one1 - zm
    s = _ac(0, prec)
    for k in 0:n
        left  = Acb[b[i] - zm for i in 1:k]       # γ₂: zₘ→1, letters shifted (none = 0)
        right = b[k+1:end]                        # γ₁: 0→zₘ
        gl = k == 0 ? _ac(1, prec) : _gpl(left,  dz, prec, cache, depth + 1, presc)
        gr = k == n ? _ac(1, prec) : _gpl(right, zm, prec, cache, depth + 1, presc)
        s += gl * gr
    end
    s
end

# ---------------------------------------------------------------------------
# Hybrid analytic+numeric continuation (the log-kernel trick)
# ---------------------------------------------------------------------------
# Do the OUTERMOST integral by parts and the INNERMOST in closed form, then
# integrate the depth-(n−2) middle numerically with Arblib's certified
# quadrature along the SAME deformed contour γ = [0,zₘ]∪[zₘ,z].  IBP turns the
# 1/(t−a₁) pole into a log kernel L₁(t)=log(1−t/a₁), so the integrand
#     L₁(t) · F'(t) ,   F(t) = G(a₂..aₙ;t),   F'(t) = G(a₃..aₙ;t)/(t−a₂),
# is HOLOMORPHIC on γ (the on-path poles are off the deformed contour, and
# log-type singularities at the real aᵢ are off γ too).  This is independent of
# the algebraic Chen decomposition above and serves as a rigorous cross-check
# and a robust fallback when the algebraic split gets combinatorially heavy.

function _ac_hybrid(a::_ACBVEC, z::Acb, prec::Int, cache, depth::Int, presc::Int)
    n = length(a)
    n == 0 && return _ac(1, prec)
    if n == 1                              # closed form on the prescribed side
        return iszero(a[1]) ? Arblib.log!(Acb(prec = prec), z) :
               _analytic_continue(a, z, prec, cache, depth, presc)
    end
    a1, a2 = a[1], a[2]
    (iszero(a1) || iszero(a2) || iszero(a[end])) &&
        error("analytic_continue(method=:hybrid): requires a₁,a₂,aₙ ≠ 0 " *
              "(shuffle out leading/trailing zeros first, or use method=:path)")
    n > 3 &&
        error("analytic_continue(method=:hybrid): depth ≤ 3 only — at depth ≥ 4 " *
              "the inner G(a₃..aₙ;t) must be enclosed on a complex BALL t and " *
              "the series/Hölder dispatcher is not ball-input efficient. " *
              "Use method=:path (the algebraic Chen split handles all depths).")
    invz = inv(z)
    b  = Acb[iszero(ai) ? _ac(0, prec) : ai * invz for ai in a]
    zm = z * _path_midpoint(b, prec, presc)       # midpoint in original t-plane
    one1 = _ac(1, prec)
    # boundary term  L₁(z) · G(a₂..aₙ; z)  — both via the algebraic path so the
    # ±iπ side of L₁(z) and any on-path letters in a₂..aₙ are handled.
    L1z = _gpl(Acb[a1], z, prec, cache, depth + 1, presc)
    G2z = _gpl(a[2:end], z, prec, cache, depth + 1, presc)
    # integrand  L₁(t) · G(a₃..aₙ; t) / (t − a₂),   t on γ = [0,zₘ]∪[zₘ,z].
    # L₁(0)=0 ⇒ integrand vanishes at the lower endpoint; on the interior of γ
    # every factor is holomorphic (poles and log branch points are off γ).
    rest = a[3:end]
    inner(t) = isempty(rest) ? _ac(1, prec) :
               _gpl(rest, Acb(t, prec = prec), prec, nothing, depth + 1, presc)
    log1(t)  = Arblib.log!(Acb(prec = prec), one1 - t / a1)
    function integrand(t; analytic::Bool = false)
        # `analytic=true` ⇒ certify holomorphy on the BALL t; the integrand has
        # a log branch point at every aᵢ and a pole at a₂, so signal
        # indeterminate when the ball overlaps any of them (forces bisection).
        # The recursive inner G is evaluated through the full dispatcher, which
        # is tuned for tight inputs — on a WIDE ball loose ball-arithmetic can
        # make a rescaled letter spuriously overlap 1; treat that as
        # indeterminate too (still rigorous: bisection narrows the ball).
        for ai in a
            Arblib.contains_zero(t - ai) &&
                return Arblib.indeterminate!(Acb(prec = prec))
        end
        try
            return log1(t) * inner(t) / (t - a2)
        catch e
            e isa ErrorException || rethrow()
            return Arblib.indeterminate!(Acb(prec = prec))
        end
    end
    seg(p, q) = Arblib.integrate(integrand, p, q;
                                 check_analytic = true, prec = prec)
    I = seg(_ac(0, prec), zm) + seg(zm, z)
    return L1z * G2z - I
end

# ---------------------------------------------------------------------------
# Master dispatcher (recursive, cached)
# ---------------------------------------------------------------------------

const _RHO_THRESH = 0.35   # series accepted when |z|/r ≤ this (≈ 220 terms for 200d)

function _gpl(a::_ACBVEC, z::Acb, prec::Int, cache::Union{Nothing,GPLCache},
              depth::Int = 0, presc::Int = 1)
    n = length(a)
    n == 0 && return _ac(1, prec)
    # cache lookup
    if cache !== nothing
        k = _cache_key(a, z, presc)
        haskey(cache.table, k) && return cache.table[k]
    end
    val = _gpl_dispatch(a, z, prec, cache, depth, presc)
    cache !== nothing && (cache.table[_cache_key(a, z, presc)] = val)
    return val
end

function _gpl_dispatch(a::_ACBVEC, z::Acb, prec::Int, cache, depth::Int, presc::Int)
    n = length(a)
    # 0. all-zero
    if all(iszero, a)
        return _gpl_allzero(n, z, prec)
    end
    # 1. trailing zeros → shuffle out
    if iszero(a[end])
        return _reduce_one_trailing_zero(a, z, prec, cache, depth, presc)
    end
    # 2. on-path letters → analytic continuation (must precede closed forms so
    #    that log/polylog never see their branch cut from an undefined side)
    if !isempty(on_path_indices(a, z))
        depth ≥ _HOELDER_MAXDEPTH &&
            error("gpl: analytic-continuation recursion did not terminate " *
                  "(letter ball straddles the deformed path at every midpoint)")
        return _analytic_continue(a, z, prec, cache, depth, presc)
    end
    # 3. closed forms
    if n == 1
        return Arblib.log!(Acb(prec = prec), _ac(1, prec) - z / a[1])
    end
    if all(iszero, @view a[1:end-1])           # G(0ᵐ⁻¹, aₙ; z) = −Liₘ(z/aₙ)
        return -li(n, z / a[end]; prec = prec)
    end
    # 4. direct series if well inside radius
    r = _rmin(a, prec)
    rho = _ubf(abs(z) / r)
    if rho ≤ _RHO_THRESH
        return _gpl_series(a, z, prec)
    end
    # 5. rescale + Hölder
    if rho < 0.7 || (depth ≥ _HOELDER_MAXDEPTH && rho < 1)
        # series still converges (slowly); use it — avoids Hölder term blow-up
        # when letters cluster near |a|=|z|, and is the depth-cap fallback
        return _gpl_series(a, z, prec)
    end
    if depth ≥ _HOELDER_MAXDEPTH
        error("gpl: Hölder recursion exhausted at depth $_HOELDER_MAXDEPTH " *
              "without reaching the series radius (a letter ball clusters " *
              "against every Hölder split point)")
    end
    invz = inv(z)
    b = Acb[iszero(ai) ? _ac(0, prec) : ai * invz for ai in a]
    return _hoelder(b, prec, cache, depth)
end

# ---------------------------------------------------------------------------
# Public entry
# ---------------------------------------------------------------------------

"""
    gpl(a, z; prec, cache=nothing, prescription=:ginac) -> Acb

Certified ball for G(a₁,…,aₙ; z). Letters may be any `Number` (promoted to
`Acb` at `prec` bits). Exact zeros must be passed as `0` (not a tiny ball).
A `GPLCache(prec)` may be supplied to memoise across calls (recommended for
basis evaluation).

If a nonzero letter lies on the open straight-line segment (0,z), the value is
fixed by `prescription`:
  * `:ginac` / `:plus_i0` (default) — aₖ → aₖ + i0; matches GiNaC's two-argument
    `G(lst{a},z)` (e.g. `G(1/2;1) = +iπ`).
  * `:feynman` / `:minus_i0` — aₖ → aₖ − i0 (e.g. `G(1/2;1) = −iπ`).
Letters with an explicit nonzero imaginary part are never on-path and are
unaffected by `prescription`.

`method` selects the evaluation strategy:
  * `:auto` (default) — series → on-path Chen split → Hölder, as in the header.
  * `:hybrid` — analytic outermost (IBP) + innermost (closed log), certified
    `Arblib.integrate` for the middle along a deformed contour; pole kernels
    become log kernels, so on-path indices are handled natively. Independent of
    the algebraic Chen split (used as a cross-check). Currently depth ≤ 3 with
    `a₁,a₂,aₙ ≠ 0` — at depth ≥ 4 the inner G must be enclosed on a complex
    ball and the dispatcher is not yet ball-input efficient.
  * `:ginac` — delegate to `gpl_ginac` (returns an Acb wrapping the BigFloat).
"""
function gpl(a::AbstractVector, z; prec::Int, cache::Union{Nothing,GPLCache} = nothing,
             prescription::Symbol = :ginac, method::Symbol = :auto)
    cache !== nothing && cache.prec != prec &&
        error("GPLCache precision $(cache.prec) ≠ requested $prec")
    av = Acb[iszero(x) ? Acb(0, prec = prec) : Acb(x, prec = prec) for x in a]
    zc = Acb(z, prec = prec)
    s  = _presc_sign(prescription)
    method === :auto   && return _gpl(av, zc, prec, cache, 0, s)
    method === :hybrid && return _ac_hybrid(av, zc, prec, cache, 0, s)
    method === :ginac  && return Acb(gpl_ginac(a, z; digits = ceil(Int, prec * log10(2.0)) + 8),
                                     prec = prec)
    error("gpl: unknown method :$method (use :auto, :hybrid, or :ginac)")
end

"""
    gpl_hybrid(a, z; prec, prescription=:ginac, cache=nothing) -> Acb

The analytic-first+last hybrid: outermost IBP turns the 1/(t−a₁) pole into a
log kernel, innermost ∫dtₙ/(tₙ−aₙ) is closed-form, and the middle is integrated
rigorously by `Arblib.integrate` along a contour deformed off any on-path
index. Robust to on-path letters; an algorithmically independent check of
`method=:auto`. Currently depth ≤ 3 with `a₁,a₂,aₙ ≠ 0`.
"""
gpl_hybrid(a::AbstractVector, z; prec::Int, prescription::Symbol = :ginac,
           cache::Union{Nothing,GPLCache} = nothing) =
    gpl(a, z; prec = prec, cache = cache, prescription = prescription, method = :hybrid)

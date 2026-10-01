# Genus-2 / hyperelliptic + Siegel-modular layer.
#
# This is the geometry class beyond the elliptic layer: a hyperelliptic curve y² = f(x) of degree
# 5 or 6 has genus 2, a 2-dimensional space of holomorphic differentials
# (dx/y, x·dx/y), and a 2×2 period matrix τ ∈ ℍ₂ (the Siegel upper half-space).
# The transport kernels on ℍ₂ are Siegel modular forms — genus-2 theta functions
# θ[a;b](z,τ) and Igusa's (ψ₄,ψ₆,χ₁₀,χ₁₂).
#
# Feynman targets that need this layer (all currently OUT of bootstrap reach):
#   • 2-loop non-planar double box, exact-m_t (Duhr–Porkert–Stawinski 2412.02300:
#     ε-factorised DE entries are Siegel modular forms)
#   • H+jet 2-loop nonplanar top sector (genus-2 hyperelliptic in Baikov,
#     Marzucca et al 2307.11497)
#   • 3-loop non-planar massless crossed box, genus-3→2 after genus drop
#     (Marzucca–McLeod–Page–Pögel–Weinzierl 2307.11497)
#
# RIGOR.  Period integrals are certified via `Arblib.integrate` with the
# `check_analytic` flag (acb_calc rigorous quadrature, same discipline as
# curve.jl).  Theta functions and Siegel modular forms are FLINT's `acb_theta`
# (Kieffer 2023, certified ball arithmetic with proved error bounds — quasi-
# linear in precision); we wrap, we do NOT reimplement.
#
# CONVENTION.  For real-ordered branch points e₁ < … < e_n (n ∈ {5,6}), the
# branch of y is the analytic continuation from y > 0 on (e_n, ∞) through the
# upper half-plane: on gap i = (e_i, e_{i+1}), 1/y = (-i)^{n-i}/√|f|.  The
# symplectic homology basis is the standard hyperelliptic one (Mumford, Tata
# Lectures II §IIIa): a_j = CCW loop around the cut [e_{2j-1}, e_{2j}], j=1,2;
# b_j = path from cut j to cut g+1 across the even gaps,
#
#     A_{kj} = ∮_{a_j} ω_k = 2·(-i)^{n-2j+1}·I^{(k)}_{2j-1},
#     B_{kj} = ∮_{b_j} ω_k = 2·Σ_{l=j}^{g} (-i)^{n-2l}·I^{(k)}_{2l},
#
# where I^{(k)}_i = ∫_{e_i}^{e_{i+1}} x^{k-1}/√|f| dx and ω_k = x^{k-1}dx/y.
# τ = A⁻¹B is then symmetric with Im τ > 0 (Riemann bilinear relations); both
# properties are ASSERTED at construction time.  This basis matches Rosenhain's
# formula in the FLINT theta-characteristic encoding (see `rosenhain_from_tau`
# and the test gauntlet).  Complex branch-point configurations are assembled
# along a simple chain through the branch points — see the chain-assembly
# section before `hyper_big_period_matrix`.
#
# FRONTIER STUBS.  The full route-E transport — the genus-2 analogue of
# `eichler_transport` (iterated integrals of Siegel modular forms over ℍ₂, or
# equivalently the Schottky–Kronecker hyperelliptic-polylog kernels) — is NOT
# in the literature in evaluable form; `siegel_transport` raises
# `SiegelNotImplemented` with the reference list.

# ---- curve descriptor --------------------------------------------------------

"""
    HyperellipticCurve

Genus-2 hyperelliptic curve y² = f(x) = lead·∏ᵢ(x - eᵢ), degree 5 or 6, with
certified branch-point enclosures.  For degree 5 the sixth branch point is ∞.

  * `prec`, `degree`, `lead :: Acb`, `roots :: Vector{Acb}` (length = degree),
  * `real_ordered :: Bool` — true iff all roots are certified real and the
    stored ordering is e₁ < e₂ < … (the period-matrix path).

Construct from exact rational coefficients [c₀,…,c_d] of f = Σ cᵢ xⁱ via
`HyperellipticCurve(coeffs, prec)`, or directly from a root vector via
`HyperellipticCurve(lead, roots, prec)`.
"""
struct HyperellipticCurve
    prec::Int
    degree::Int
    lead::Acb
    roots::Vector{Acb}
    real_ordered::Bool
end

function HyperellipticCurve(lead, roots::Vector, prec::Int)
    d = length(roots)
    d == 5 || d == 6 || error("HyperellipticCurve: need degree 5 or 6 (genus 2)")
    rts = [Acb(r; prec) for r in roots]
    realr = all(r -> Arblib.contains_zero(imag(r)), rts)
    if realr
        # tighten + certify ordering
        rts = [Acb(real(r), Arb(0; prec); prec) for r in rts]
        idx = sortperm(Float64.(real.(rts)))
        rts = rts[idx]
        for i in 1:d-1
            real(rts[i]) < real(rts[i+1]) ||
                error("HyperellipticCurve: cannot certify strict root ordering")
        end
    end
    HyperellipticCurve(prec, d, Acb(lead; prec), rts, realr)
end

function HyperellipticCurve(coeffs::Vector, prec::Int)
    d = length(coeffs) - 1
    d == 5 || d == 6 || error("HyperellipticCurve: need degree 5 or 6")
    p = AcbPoly(prec = prec)
    for (i, c) in enumerate(coeffs)
        Arblib.set_coeff!(p, i - 1, Acb(c; prec))
    end
    roots = AcbVector(d; prec)
    n = ccall((:acb_poly_find_roots, Arblib.libflint), Clong,
        (Ptr{Arblib.acb_struct}, Ref{Arblib.acb_poly_struct}, Ptr{Cvoid}, Clong, Clong),
        roots, p, C_NULL, 0, prec)
    n == d || error("HyperellipticCurve: could not certify $d isolated roots")
    rts = [Acb(roots[i]; prec) for i in 1:d]
    # rigorous realness tightening (as in QuarticCurve): real-coeff poly ⇒ conj
    # of an isolated root is a root; if it overlaps no other ball, it's its own
    # conjugate.
    if all(c -> Arblib.contains_zero(imag(Acb(c; prec))), coeffs)
        for i in 1:d
            Arblib.contains_zero(imag(rts[i])) || continue
            ci = conj(rts[i])
            if !any(j != i && Arblib.overlaps(ci, rts[j]) for j in 1:d)
                rts[i] = Acb(real(rts[i]), Arb(0; prec); prec)
            end
        end
    end
    HyperellipticCurve(Acb(Arblib.ref(p, d); prec), rts, prec)
end

genus(C::HyperellipticCurve) = 2

# ---- certified gap integrals -------------------------------------------------

"""
    hyper_gap_integral(C, i, k) :: Acb

Certified ∫_{eᵢ}^{e_{i+1}} x^k / √|f(x)| dx for the real-ordered curve C, via
the trig desingularisation x = m + δ·sinθ, θ ∈ [-π/2, π/2] (both 1/√ endpoints
absorbed into dθ).  The integrand √|f|⁻¹ is built as the principal square root
of the manifestly-positive product |lead|·∏_{j<i}(x-e_j)·∏_{j>i+1}(e_j-x),
so analyticity on every quadrature box is checked rigorously.
"""
function hyper_gap_integral(C::HyperellipticCurve, i::Int, k::Int; rtol = 0.0)
    C.real_ordered ||
        error("hyper_gap_integral: real-ordered branch points required")
    prec = C.prec; e = C.roots; n = C.degree
    1 ≤ i ≤ n - 1 || error("hyper_gap_integral: gap index out of range")
    a, b = e[i], e[i+1]
    m = (a + b) / 2; δ = (b - a) / 2
    al = abs(C.lead)
    integrand = let m = m, δ = δ, e = e, n = n, i = i, al = al, k = k, prec = prec
        (θ; analytic::Bool = false) -> begin
            x = m + δ * sin(θ)
            pr = Acb(al; prec)
            for j in 1:n
                (j == i || j == i + 1) && continue
                pr *= (j < i ? (x - e[j]) : (e[j] - x))
            end
            if analytic && Arblib.contains_zero(imag(pr)) && !Arblib.is_positive(real(pr))
                ind = Acb(prec = prec); Arblib.indeterminate!(ind); return ind
            end
            x^k / sqrt(pr)
        end
    end
    Arblib.integrate(integrand, Acb(-Arb(π; prec) / 2; prec), Acb(Arb(π; prec) / 2; prec);
                     check_analytic = true, prec = prec,
                     rtol = rtol == 0.0 ? exp10(-0.28 * prec) : rtol)
end

# ---- complex branch-point configurations: chain assembly ---------------------
#
# For branch points off the real line the same Mumford construction applies
# along any SIMPLE polygonal chain e_{σ(1)} → … → e_{σ(n)} through the branch
# points (no self-crossings, no other branch point on a segment): such a chain
# is ambient-isotopic to the real-ordered configuration, so the cuts are the
# odd segments [e_{σ(2j-1)}, e_{σ(2j)}], the a/b combinatorics are unchanged,
# and only the per-segment sheet factors differ.  On each segment the certified
# primitive is `hyper_chain_integral` (the trig-desingularised integral against
# a DETERMINISTIC smooth branch of √f built from per-factor principal square
# roots); the true 1/y on the segment is ±1 times that branch (both square to
# 1/f and are continuous, so the ratio is a constant sign).  The 2^{2g-1} sheet
# assignments (one global flip quotiented out) are screened by the Riemann
# bilinear relations — τ = A⁻¹B symmetric with Im τ ≻ 0 — which certifiably
# exclude wrong assignments at moderate precision; `period_matrix` additionally
# cross-checks the survivor against the curve's Igusa invariants through
# FLINT's Thomae inverse.  Exactly one survivor is required; anything else is a
# loud error, never a guess.

"""
    hyper_default_chain(C) → Vector{Int}

Heuristic non-crossing chain through the branch points of `C`: the permutation
of `1:degree` minimising the total Euclidean length of the polygonal chain
through the root midpoints (shortest open paths in the plane are non-crossing).
Brute force over all orderings (n ≤ 6), reversal-symmetric duplicates removed.
The result is a HEURISTIC ordering only — every certified statement about the
periods is re-derived from it by `hyper_chain_integral` and the assembly
screens, so a bad chain can only produce a loud failure, not a wrong number.
"""
function hyper_default_chain(C::HyperellipticCurve)
    n = C.degree
    pts = [ComplexF64(Float64(real(r)), Float64(imag(r))) for r in C.roots]
    best = collect(1:n); bestlen = Inf
    for p in _all_orderings(collect(1:n))
        p[1] < p[end] || continue           # quotient out reversal
        len = sum(abs(pts[p[i+1]] - pts[p[i]]) for i in 1:n-1)
        if len < bestlen
            bestlen = len; best = p
        end
    end
    return best
end

# All orderings of a small index vector (n ≤ 6 here; stdlib only).
function _all_orderings(v::Vector{Int})
    length(v) <= 1 && return [v]
    out = Vector{Vector{Int}}()
    for (i, x) in enumerate(v)
        rest = [v[1:i-1]; v[i+1:end]]
        for p in _all_orderings(rest)
            push!(out, [x; p])
        end
    end
    return out
end

# A fixed square root of a CONSTANT ball, well-defined even when the ball
# straddles the principal cut (tiny root-enclosure radii put e.g. an exactly
# negative real constant a hair off the axis, where acb sqrt returns the
# useless both-branches ball).  Which root is returned is irrelevant — every
# caller determines the overall sheet sign separately — but it must be ONE
# root, tightly.
function _const_sqrt(q::Acb, prec::Int)
    if Arblib.contains_zero(imag(q)) && !Arblib.is_positive(real(q))
        Arblib.is_negative(real(q)) ||
            error("_const_sqrt: constant not certified away from zero")
        return Acb(0, 1; prec) * sqrt(-q)
    end
    return sqrt(q)
end

# Per-factor cut-safe orientation for a root at scaled position u_j relative to
# the segment (u ∈ [-1,1]): σ = +1 uses the factor (u - u_j), σ = -1 uses
# (u_j - u); the choice keeps the factor certifiably off the principal-branch
# cut (-∞, 0] for all u ∈ [-1,1].  Errors when the root cannot be certified off
# the closed segment (the chain crosses a branch point).
function _segment_factor_sign(uj::Acb)
    if Arblib.is_negative(imag(uj))
        return 1
    elseif Arblib.is_positive(imag(uj))
        return -1
    elseif Arblib.is_negative(real(uj) + 1)
        return 1
    elseif Arblib.is_positive(real(uj) - 1)
        return -1
    end
    error("hyper_chain_integral: a branch point cannot be certified off the " *
          "chain segment — supply a different chain ordering")
end

"""
    hyper_chain_integral(C, chain, i, k; rtol=0.0) :: Acb

Certified ∫ x^k dx/ŷ over the straight segment e_{chain[i]} → e_{chain[i+1]},
where ŷ is the DETERMINISTIC smooth branch of √f(x) on the segment built as
δ·cosθ·√q·∏_j √(σ_j(u-u_j)) under the trig desingularisation x = m + δ·sinθ
(u = sinθ, u_j the scaled other roots, σ_j the cut-safe orientations,
q = -lead·δ^{n-2}·∏σ_j).  The true dx/y on the segment for any global branch
of y equals ±1 times the integrand (both square to the same function), so this
is the sheet-sign-free primitive the chain assembly combines.  Analyticity of
every per-factor principal square root is checked rigorously box-by-box
(`check_analytic`), exactly as in `hyper_gap_integral`.
"""
function hyper_chain_integral(C::HyperellipticCurve, chain::Vector{Int}, i::Int, k::Int;
                              rtol = 0.0)
    prec = C.prec; e = C.roots; n = C.degree
    (length(chain) == n && isperm(chain)) ||
        error("hyper_chain_integral: chain must be a permutation of 1:$n")
    1 ≤ i ≤ n - 1 || error("hyper_chain_integral: segment index out of range")
    a, b = e[chain[i]], e[chain[i+1]]
    m = (a + b) / 2; δ = (b - a) / 2
    others = [e[chain[j]] for j in 1:n if j != i && j != i + 1]
    us = [(ej - m) / δ for ej in others]
    σs = [_segment_factor_sign(uj) for uj in us]
    q = -C.lead * δ^(n - 2) * prod(σs)
    sq = _const_sqrt(q, prec)          # constant: any fixed square root works
    integrand = let m = m, δ = δ, us = us, σs = σs, sq = sq, k = k, prec = prec
        (θ; analytic::Bool = false) -> begin
            u = sin(θ)
            x = m + δ * u
            w = Acb(1; prec)
            for (uj, σj) in zip(us, σs)
                φ = σj > 0 ? (u - uj) : (uj - u)
                if analytic && Arblib.contains_zero(imag(φ)) && !Arblib.is_positive(real(φ))
                    ind = Acb(prec = prec); Arblib.indeterminate!(ind); return ind
                end
                w *= sqrt(φ)
            end
            x^k / (sq * w)
        end
    end
    Arblib.integrate(integrand, Acb(-Arb(π; prec) / 2; prec), Acb(Arb(π; prec) / 2; prec);
                     check_analytic = true, prec = prec,
                     rtol = rtol == 0.0 ? exp10(-0.28 * prec) : rtol)
end

# Float64 midpoint pre-check that the chain is a simple (non-self-crossing)
# polyline.  Heuristic screen only (certification happens downstream); a
# crossing chain is refused loudly here because the homology combinatorics
# assume a simple arc.
function _chain_is_simple(C::HyperellipticCurve, chain::Vector{Int})
    pts = [ComplexF64(Float64(real(r)), Float64(imag(r))) for r in C.roots]
    p = [pts[c] for c in chain]
    orient(a, b, c) = sign(real(conj(b - a) * im * (c - a)))  # sign of cross product
    n = length(p)
    for i in 1:n-1, j in i+2:n-1
        # segments (i, i+1) and (j, j+1), non-adjacent (j ≥ i+2)
        o1 = orient(p[i], p[i+1], p[j]); o2 = orient(p[i], p[i+1], p[j+1])
        o3 = orient(p[j], p[j+1], p[i]); o4 = orient(p[j], p[j+1], p[i+1])
        if o1 != o2 && o3 != o4 && o1 != 0 && o2 != 0 && o3 != 0 && o4 != 0
            return false
        end
    end
    return true
end

function _hyper_chain_big_periods(C::HyperellipticCurve, chain::Vector{Int}; rtol = 0.0)
    prec = C.prec; n = C.degree; g = 2
    (length(chain) == n && isperm(chain)) ||
        error("hyper_big_period_matrix: chain must be a permutation of 1:$n")
    _chain_is_simple(C, chain) ||
        error("hyper_big_period_matrix: chain self-crosses — supply a " *
              "non-crossing ordering (see hyper_default_chain)")
    Ĩ = Dict((i, k) => hyper_chain_integral(C, chain, i, k; rtol = rtol)
             for i in 1:2g, k in 0:g-1)
    # sheet-sign search over the 2g used segments, global flip quotiented out
    # (s and -s give the same τ); Riemann bilinear relations screen the rest.
    survivors = Tuple{NTuple{4,Int},AcbMatrix,AcbMatrix}[]
    for s2 in (1, -1), s3 in (1, -1), s4 in (1, -1)
        s = (1, s2, s3, s4)
        A = AcbMatrix(g, g; prec); B = AcbMatrix(g, g; prec)
        for k in 1:g, j in 1:g
            A[k, j] = 2 * s[2j-1] * Ĩ[(2j - 1, k - 1)]
            acc = Acb(0; prec)
            for l in j:g
                acc += s[2l] * Ĩ[(2l, k - 1)]
            end
            B[k, j] = 2 * acc
        end
        Ainv = AcbMatrix(g, g; prec); Arblib.inv!(Ainv, A)
        τ = Ainv * B
        Arblib.overlaps(τ[1, 2], τ[2, 1]) || continue
        detIm = imag(τ[1, 1]) * imag(τ[2, 2]) - imag(τ[1, 2])^2
        (Arblib.is_positive(imag(τ[1, 1])) && Arblib.is_positive(detIm)) || continue
        push!(survivors, (s, A, B))
    end
    length(survivors) == 1 ||
        error("hyper_big_period_matrix: $(length(survivors)) sheet assignments " *
              "pass the Riemann bilinear screen (need exactly 1) — increase " *
              "prec or supply a different chain")
    return survivors[1][2], survivors[1][3]
end

"""
    hyper_big_period_matrix(C; rtol=0.0, chain=nothing) → (A, B)

Certified 2×2 a- and b-period matrices A_{kj}=∮_{a_j}ω_k, B_{kj}=∮_{b_j}ω_k of
the holomorphic basis ω_k = x^{k-1}dx/y on the symplectic homology basis
described in the module docstring.  Real-ordered branch points take the direct
gap-integral path; complex configurations are assembled along a simple chain
through the branch points (`chain` — a permutation of root indices — overrides
`hyper_default_chain`), with the per-segment sheet signs resolved by the
Riemann-bilinear screen (see the chain-assembly notes above).  Passing an
explicit `chain` forces the chain path even for real-ordered curves.
"""
function hyper_big_period_matrix(C::HyperellipticCurve; rtol = 0.0, chain = nothing)
    if !(C.real_ordered && chain === nothing)
        ch = chain === nothing ? hyper_default_chain(C) : chain
        return _hyper_chain_big_periods(C, ch; rtol = rtol)
    end
    prec = C.prec; n = C.degree; g = 2
    I = Dict((i, k) => hyper_gap_integral(C, i, k; rtol = rtol) for i in 1:2g, k in 0:g-1)
    mi = Acb(0, -1; prec)
    J = Dict((i, k) => mi^(n - i) * I[(i, k)] for i in 1:2g, k in 0:g-1)
    A = AcbMatrix(g, g; prec); B = AcbMatrix(g, g; prec)
    for k in 1:g, j in 1:g
        A[k, j] = 2 * J[(2j - 1, k - 1)]
        s = Acb(0; prec)
        for l in j:g
            s += J[(2l, k - 1)]
        end
        B[k, j] = 2 * s
    end
    return A, B
end

"""
    period_matrix(C::HyperellipticCurve; rtol=0.0, chain=nothing) → τ :: AcbMatrix (2×2)

Small period matrix τ = A⁻¹B ∈ ℍ₂.  Asserts (with certified ball checks) that
τ is symmetric and Im τ is positive-definite — both are theorems (Riemann
bilinear relations) so a failure here is a bug, not a numerical accident.
Complex branch-point configurations go through the chain assembly (`chain`
passed to `hyper_big_period_matrix`) and are additionally gated by the
Igusa-invariant round-trip against FLINT's Thomae sextic (Torelli/Igusa: the
absolute invariants agree iff τ is a period matrix of C up to Sp(4,ℤ)); the
gate is skipped on the locus where the Clebsch degree-10 invariant vanishes
(the `igusa_absolute` ratios are undefined there, e.g. the Bolza curve).
Multiple-dispatch overload alongside `period_matrix(::PFOperator, …)`.
"""
function period_matrix(C::HyperellipticCurve; rtol = 0.0, chain = nothing)
    prec = C.prec; g = 2
    chained = !(C.real_ordered && chain === nothing)
    A, B = hyper_big_period_matrix(C; rtol = rtol, chain = chain)
    Ainv = AcbMatrix(g, g; prec); Arblib.inv!(Ainv, A)
    τ = Ainv * B
    # Riemann bilinear: τ symmetric, Im τ ≻ 0
    Arblib.overlaps(τ[1, 2], τ[2, 1]) ||
        error("period_matrix(HyperellipticCurve): τ not symmetric — basis bug")
    # symmetrise to tighten (component-wise ball intersection with the transpose)
    sre = Arb(prec = prec); sim = Arb(prec = prec)
    Arblib.intersection!(sre, real(τ[1, 2]), real(τ[2, 1]))
    Arblib.intersection!(sim, imag(τ[1, 2]), imag(τ[2, 1]))
    s12 = Acb(sre, sim; prec)
    τ[1, 2] = s12; τ[2, 1] = s12
    detIm = imag(τ[1, 1]) * imag(τ[2, 2]) - imag(τ[1, 2])^2
    (Arblib.is_positive(imag(τ[1, 1])) && Arblib.is_positive(detIm)) ||
        error("period_matrix(HyperellipticCurve): Im τ not certified ≻ 0")
    chained && _igusa_roundtrip_gate(C, τ)
    return τ
end

# Cross-check of a chain-assembled τ: the absolute Igusa invariants of the
# input sextic must overlap those of the Thomae sextic recovered from τ
# (Torelli + Igusa: they agree exactly when τ is Sp(4,ℤ)-equivalent to a
# period matrix of C).  Skipped when either Clebsch degree-10 invariant
# contains zero — the `igusa_absolute` ratios are undefined on that locus
# (e.g. the maximal-automorphism Bolza curve); there the Riemann-bilinear
# assertions above remain the certificate.
function _igusa_roundtrip_gate(C::HyperellipticCurve, τ::AcbMatrix)
    prec = C.prec
    ic_in = igusa_clebsch(C)
    f̃, _ = thomae_sextic(τ, prec)
    cf(i) = (a = Acb(prec = prec); Arblib.get_coeff!(a, f̃, i); a)
    ic_rec = igusa_clebsch([cf(i) for i in 0:6], prec)
    (Arblib.contains_zero(ic_in[4]) || Arblib.contains_zero(ic_rec[4])) &&
        return nothing
    j_in = igusa_absolute(ic_in); j_rec = igusa_absolute(ic_rec)
    all(Arblib.overlaps(j_in[k], j_rec[k]) for k in 1:3) ||
        error("period_matrix: chain-assembled τ fails the Igusa-invariant " *
              "round-trip against the Thomae sextic — sheet/homology assembly " *
              "inconsistency or insufficient precision")
    return nothing
end

# ---- Siegel theta layer (FLINT acb_theta wrap) -------------------------------

"FLINT genus-g characteristic encoding: 2g-bit integer (a₁…a_g b₁…b_g), MSB
first; for g=2, char (a₁,a₂;b₁,b₂)∈{0,1}⁴ ↦ 8a₁+4a₂+2b₁+b₂."
siegel_char(a1, a2, b1, b2) = UInt(8a1 + 4a2 + 2b1 + b2)
siegel_char(ab::NTuple{4,Int}) = siegel_char(ab...)

"Parity of a genus-g half-integer characteristic ab (a·b mod 2)."
function siegel_char_is_even(ab::Unsigned, g::Int)
    s = 0
    for j in 1:g
        s += ((ab >> (2g - j)) & 1) * ((ab >> (g - j)) & 1)
    end
    return iseven(s)
end

"The 10 even genus-2 characteristics in FLINT encoding."
const SIEGEL_G2_EVEN_CHARS = UInt[ab for ab in 0:15 if siegel_char_is_even(UInt(ab), 2)]

"""
    siegel_theta(ab, z, τ, prec) :: Acb

Certified genus-g Riemann/Siegel theta function θ[a;b](z,τ) for half-integer
characteristic `ab` (FLINT encoding, see `siegel_char`), z ∈ ℂ^g, τ ∈ ℍ_g.
Thin wrapper over FLINT `acb_theta_one` (Kieffer 2023): error bounds are
proved, complexity is quasi-linear in `prec`.

`ab` may be a `UInt` (FLINT encoding) or an `NTuple{2g,Int}`.
"""
function siegel_theta(ab::Unsigned, z::AcbVector, τ::AcbMatrix, prec::Int)
    th = Acb(prec = prec)
    Arblib.theta_one!(th, z, τ, ab, prec)
    return th
end
siegel_theta(ab::NTuple{4,Int}, z, τ, prec) = siegel_theta(siegel_char(ab), z, τ, prec)
function siegel_theta(ab, z::Vector{Acb}, τ::AcbMatrix, prec::Int)
    zv = AcbVector(length(z); prec)
    for (i, zi) in enumerate(z); zv[i] = zi; end
    siegel_theta(ab, zv, τ, prec)
end

"""
    siegel_theta_all(z, τ, prec; sqr=false) → Vector{Acb} (length 2^{2g})

All 2^{2g} theta values θ[ab](z,τ) (or their squares if `sqr=true`), indexed by
the FLINT characteristic 0…2^{2g}-1.  Wraps FLINT `acb_theta_all`.
"""
function siegel_theta_all(z::AcbVector, τ::AcbMatrix, prec::Int; sqr::Bool = false)
    g = size(τ, 1)
    th = AcbVector(1 << (2g); prec)
    Arblib.theta_all!(th, z, τ, sqr ? 1 : 0, prec)
    return [Acb(th[i]; prec) for i in 1:(1 << (2g))]
end
function siegel_theta_all(z::Vector, τ::AcbMatrix, prec::Int; sqr::Bool = false)
    g = size(τ, 1)
    zv = AcbVector(g; prec)
    for (i, zi) in enumerate(z); zv[i] = Acb(zi; prec); end
    siegel_theta_all(zv, τ, prec; sqr = sqr)
end

"""
    igusa_invariants(τ, prec) → (ψ₄, ψ₆, χ₁₀, χ₁₂) :: NTuple{4,Acb}

The four scalar generators of the ring of even-weight genus-2 Siegel modular
forms, normalised as in FLINT/Streng (ψ₄ = 1 + 240(q₁+q₃)+…, χ₁₀ =
2⁻¹²·∏_{even}θ², …).  Wraps FLINT `acb_theta_g2_even_weight`.

These are the genus-2 analogue of (E₄,E₆,Δ); under the covariants↔modular-forms
correspondence (CFG2017) the absolute Igusa invariants of the underlying sextic
are recovered by `igusa_clebsch` / `igusa_absolute` (see those for the
algebraic side and the round-trip test).
"""
function igusa_invariants(τ::AcbMatrix, prec::Int)
    g = size(τ, 1)
    g == 2 || error("igusa_invariants: genus-2 only")
    z = AcbVector(2; prec)
    th2 = AcbVector(16; prec)
    Arblib.theta_all!(th2, z, τ, 1, prec)
    ψ4 = Acb(prec = prec); ψ6 = Acb(prec = prec)
    χ10 = Acb(prec = prec); χ12 = Acb(prec = prec)
    Arblib.theta_g2_even_weight!(ψ4, ψ6, χ10, χ12, th2, prec)
    return (ψ4, ψ6, χ10, χ12)
end

"""
    igusa_clebsch(C::HyperellipticCurve) → (C₂₀, C₄₀, C₆₀, C₁₀₀) :: NTuple{4,Acb}
    igusa_clebsch(coeffs, prec)          → idem

Algebraic scalar covariants (degrees 2,4,6,10) of the binary sextic f, computed
via the FLINT transvectant chain (`acb_theta_g2_covariants`).  Intermediate
scalings are dropped (homogeneous, so cancel in `igusa_absolute`).  For a
degree-5 `HyperellipticCurve` the branch point at ∞ is brought to a finite
point by x ↦ c + 1/u first (c chosen off the roots), so the binary sextic is
non-degenerate.

`igusa_absolute(C)` returns the three weight-0 ratios (C₂₀⁵/C₁₀₀, C₂₀³C₄₀/C₁₀₀,
C₂₀²C₆₀/C₁₀₀) — isomorphism invariants of the curve, hence equal for the input
sextic and the FLINT-recovered Thomae sextic `thomae_sextic(τ)` whenever
τ = `period_matrix(C)`.  CAVEAT (test-helper, not bootstrap data): C₁₀₀ is the
Clebsch degree-10 invariant (NOT the discriminant); it vanishes on a
codimension-1 locus including the Bolza curve y² = x⁵−x, where the ratios are
undefined — use `igusa_invariants(τ)` (χ₁₀ ≠ 0 off the diagonal ℍ₁²) instead.
"""
function igusa_clebsch(coeffs::Vector, prec::Int)
    f = AcbPoly(prec = prec)
    for (i, c) in enumerate(coeffs); Arblib.set_coeff!(f, i - 1, Acb(c; prec)); end
    # Transvectant chain (FLINT g2_covariants).  Intermediate scalings (75,30…)
    # are omitted: they contribute homogeneous powers that cancel in the
    # absolute invariants `igusa_absolute`, and scale the bare (C₂₀,C₄₀,C₆₀,
    # C₁₀₀) by fixed rationals — irrelevant for any comparison made here.
    C24 = AcbPoly(prec = prec); Arblib.theta_g2_transvectant!(C24, f, f, 6, 6, 4, 0, prec)
    C32 = AcbPoly(prec = prec); Arblib.theta_g2_transvectant!(C32, f, C24, 6, 4, 4, 0, prec)
    C20 = AcbPoly(prec = prec); Arblib.theta_g2_transvectant!(C20, f, f, 6, 6, 6, 0, prec)
    C40 = AcbPoly(prec = prec); Arblib.theta_g2_transvectant!(C40, C24, C24, 4, 4, 4, 0, prec)
    C60 = AcbPoly(prec = prec); Arblib.theta_g2_transvectant!(C60, C32, C32, 2, 2, 2, 0, prec)
    C32_3 = C32 * C32 * C32
    C100 = AcbPoly(prec = prec); Arblib.theta_g2_transvectant!(C100, f, C32_3, 6, 6, 6, 0, prec)
    cf(p) = (c = Acb(prec = prec); Arblib.get_coeff!(c, p, 0); c)
    return (cf(C20), cf(C40), cf(C60), cf(C100))
end
function igusa_clebsch(C::HyperellipticCurve)
    prec = C.prec
    rts = C.roots; lead = C.lead
    if C.degree == 5
        # Branch point at ∞ ⇒ the binary sextic is degenerate (C₁₀₀ = 0).  Bring
        # ∞ to a finite point by x = c + 1/u with c ∉ {roots}: new roots are
        # 1/(eᵢ-c) and 0, new lead = lead·∏(eᵢ-c).
        c = Acb(0; prec)
        while any(Arblib.overlaps(c, r) for r in rts)
            c += 1
        end
        lead = lead * prod(r - c for r in rts)
        rts = [Acb(0; prec); [inv(r - c) for r in rts]]
    end
    p = AcbPoly([lead]; prec)
    for r in rts
        p = p * AcbPoly([Acb(-r; prec), Acb(1; prec)]; prec)
    end
    cf(i) = (a = Acb(prec = prec); Arblib.get_coeff!(a, p, i); a)
    igusa_clebsch([cf(i) for i in 0:6], prec)
end

function igusa_absolute(ic::NTuple{4,Acb})
    C20, C40, C60, C100 = ic
    return (C20^5 / C100, C20^3 * C40 / C100, C20^2 * C60 / C100)
end
igusa_absolute(C::HyperellipticCurve) = igusa_absolute(igusa_clebsch(C))

"""
    rosenhain_from_tau(τ, prec) → (λ₁, λ₂, λ₃) :: NTuple{3,Acb}

Rosenhain invariants (three moduli of the genus-2 curve in the normal form
y² = x(x-1)(x-λ₁)(x-λ₂)(x-λ₃)) from theta-nulls, in the homology-basis
convention of `period_matrix(::HyperellipticCurve)`:

    λ₁ = e₃ = θ₀²θ₄² / (θ₈²θ₁₂²),
    λ₂ = e₄ = θ₁²θ₄² / (θ₉²θ₁₂²),
    λ₃ = e₅ = θ₀²θ₁² / (θ₈²θ₉²),

with FLINT characteristic indices (0=[00;00], 1=[00;01], 4=[01;00], 8=[10;00],
9=[10;01], 12=[11;00]).  This is Thomae's formula specialised to e₁=0, e₂=1
(degree 5, e₆=∞).
"""
function rosenhain_from_tau(τ::AcbMatrix, prec::Int)
    th = siegel_theta_all(AcbVector(2; prec), τ, prec)
    T2(i) = th[i + 1]^2
    return (T2(0)*T2(4)/(T2(8)*T2(12)),
            T2(1)*T2(4)/(T2(9)*T2(12)),
            T2(0)*T2(1)/(T2(8)*T2(9)))
end

"""
    thomae_sextic(τ, prec) → (f::AcbPoly, χ₅::Acb)

FLINT's Thomae/Rosenhain inverse: returns the binary sextic χ_{-2,6}(τ) (the
covariant corresponding to the curve itself) such that y² = f(x) has Jacobian
with period matrix Sp(4,ℤ)-equivalent to τ.  Wraps `acb_theta_g2_sextic_chi5`.
"""
function thomae_sextic(τ::AcbMatrix, prec::Int)
    f = AcbPoly(prec = prec); χ5 = Acb(prec = prec)
    Arblib.theta_g2_sextic_chi5!(f, χ5, τ, prec)
    return f, χ5
end

# ---- Abel–Jacobi -------------------------------------------------------------

"""
    abel_jacobi(C, xp, prec=C.prec; base = C.degree, rtol=0.0) → Vector{Acb} (length 2)

Abel–Jacobi image  u = A⁻¹·∫_{e_base}^{x_p} (dx/y, x·dx/y)ᵀ  ∈ ℂ²/(ℤ² + τℤ²),
for `xp` with Re(xp) strictly to the right of `e_base` (inside the basepoint's
gap when base < n).  `base` is the index of the basepoint branch point
(default: the rightmost finite branch point e_n).  The representative is the
one produced by the documented path; the caller reduces mod the period lattice.

Real xp is integrated on the +i0 branch of the gap-integral convention.
Complex xp (the Feynman-deformed marked points of the module conventions) is
reached by the two-leg path e_base → Re(xp) → xp: the real leg on the +i0
branch as before, then a straight vertical leg along which y is continued
analytically via a deterministic per-factor branch whose sheet sign is
CERTIFIED against the +i0 convention at the anchor Re(xp) (the certification
fails loudly, never silently, if the precision cannot separate the sheets).
Im(xp) of either sign is accepted; for Im(xp) < 0 on a cut gap the value is
the continuation of the +i0 branch through the cut (second sheet), which is
the convention's meaning of approaching from above.

Restricted to real-ordered curves (for complex branch-point configurations the
chain period matrix exists; the Abel–Jacobi chain path does not yet).
"""
function abel_jacobi(C::HyperellipticCurve, xp; base::Int = C.degree, rtol = 0.0)
    C.real_ordered ||
        error("abel_jacobi: real-ordered branch points required (complex " *
              "branch-point configurations: chain period matrix available, " *
              "Abel–Jacobi chain path not yet)")
    prec = C.prec; e = C.roots; n = C.degree; g = 2
    1 ≤ base ≤ n || error("abel_jacobi: base index out of range")
    xpA = Acb(xp; prec)
    deformed = !Arblib.contains_zero(imag(xpA))
    xbA = deformed ? Acb(real(xpA), Arb(0; prec); prec) : xpA   # real anchor
    real(xbA) > real(e[base]) ||
        error("abel_jacobi: need Re(xp) > e_base on this branch")
    base < n && (real(xbA) < real(e[base + 1]) ||
        error("abel_jacobi: Re(xp) must lie in the basepoint's gap"))
    # ∫_{e_base}^{x_b} x^k/√|f| via x = e_base + u²·(x_b - e_base): single 1/√
    # endpoint at e_base, integrand smooth at u=1.
    Δ = xbA - e[base]
    al = abs(C.lead)
    # Desingularisation:  √|f| = √((x-e_base)·∏_{j≠base}|x-e_j|),
    #   x-e_base = u²Δ ⇒ dx/√(x-e_base) = 2√Δ du.
    incomplete = k -> begin
        ig2 = (u; analytic::Bool = false) -> begin
            x = e[base] + u * u * Δ
            pr = Acb(al; prec)
            for j in 1:n
                j == base && continue
                pr *= (j < base ? (x - e[j]) : (e[j] - x))
            end
            if analytic && Arblib.contains_zero(imag(pr)) && !Arblib.is_positive(real(pr))
                ind = Acb(prec = prec); Arblib.indeterminate!(ind); return ind
            end
            2 * sqrt(Δ) * x^k / sqrt(pr)
        end
        Arblib.integrate(ig2, Acb(0; prec), Acb(1; prec);
                         check_analytic = true, prec = prec,
                         rtol = rtol == 0.0 ? exp10(-0.28 * prec) : rtol)
    end
    mi = Acb(0, -1; prec)
    # 1/y = (-i)^{n-base}/√|f| on (e_base, e_{base+1}) (or (e_n,∞) where the
    # exponent is 0); the same UHP-continuation phase as the period gaps.
    phase = mi^(n - base)
    raw = [phase * incomplete(k) for k in 0:g-1]
    if deformed
        # Vertical continuation leg x_b → xp.  y along the leg is the smooth
        # branch  sgn·√q̂·∏_j √(σ_j(x − e_j))  with σ_j = +1 for j ≤ base and
        # −1 for j > base (every factor positive at the real anchor, hence off
        # the principal cut on a neighbourhood of the leg — the roots are real,
        # the open leg is not), q̂ = al·(−1)^{n−base}, and sgn ∈ {±1} fixed by
        # a certified comparison with the +i0 convention value at x_b.
        ih = Acb(0, 1; prec) * imag(xpA)
        σv = [j <= base ? 1 : -1 for j in 1:n]
        sqv = _const_sqrt(Acb(al * (-1)^(n - base); prec), prec)
        φb = [σv[j] > 0 ? (xbA - e[j]) : (e[j] - xbA) for j in 1:n]
        ybranch = sqv * prod(sqrt.(φb))
        yplus = Acb(0, 1; prec)^(n - base) * sqrt(al * prod(real.(φb)))
        r = ybranch / yplus
        o1 = Arblib.overlaps(r, Acb(1; prec))
        o2 = Arblib.overlaps(r, Acb(-1; prec))
        xor(o1, o2) ||
            error("abel_jacobi: cannot certify the sheet sign at the " *
                  "deformation anchor — increase prec")
        sgn = o1 ? 1 : -1
        vert = k -> begin
            igv = let xbA = xbA, ih = ih, e = e, σv = σv, n = n, k = k,
                      sqv = sqv, sgn = sgn, prec = prec
                (v; analytic::Bool = false) -> begin
                    x = xbA + v * ih
                    w = Acb(1; prec)
                    for j in 1:n
                        φ = σv[j] > 0 ? (x - e[j]) : (e[j] - x)
                        if analytic && Arblib.contains_zero(imag(φ)) &&
                           !Arblib.is_positive(real(φ))
                            ind = Acb(prec = prec)
                            Arblib.indeterminate!(ind)
                            return ind
                        end
                        w *= sqrt(φ)
                    end
                    ih * x^k / (sgn * sqv * w)
                end
            end
            Arblib.integrate(igv, Acb(0; prec), Acb(1; prec);
                             check_analytic = true, prec = prec,
                             rtol = rtol == 0.0 ? exp10(-0.28 * prec) : rtol)
        end
        raw = [raw[k + 1] + vert(k) for k in 0:g-1]
    end
    A, _ = hyper_big_period_matrix(C; rtol = rtol)
    Ainv = AcbMatrix(g, g; prec); Arblib.inv!(Ainv, A)
    u = AcbMatrix(g, 1; prec)
    for k in 1:g; u[k, 1] = raw[k]; end
    v = Ainv * u
    return [Acb(v[k, 1]; prec) for k in 1:g]
end

# ---- Gauss–Manin connection (Rauch variational formula) ----------------------
#
# The genus-2 analogue of the elliptic  dτ = (W/Ψ₁²) dt.  For the single scalar
# modulus τ of an elliptic curve the period ratio's t-derivative is the Wronskian
# over Ψ₁²; for the 2×2 Siegel τ the correct object is the Gauss–Manin connection
# ∂τ/∂t, a symmetric 2×2 matrix.  It is given EXACTLY (no finite differencing) by
# Rauch's variational formula (Rauch 1959; Fay, Theta Functions §III; Mumford,
# Tata Lectures II): under a deformation of the branch points e_m,
#
#     ∂τ_{jk}/∂e_m  =  4πi · ĉ_j(e_m) · ĉ_k(e_m) / f'(e_m),
#
# where ω̂_j = Σ_k (A⁻¹)_{jk} ω_k = (ĉ_{j,0} + ĉ_{j,1} x) dx/y is the NORMALISED
# holomorphic differential (∮_{a_i} ω̂_j = δ_{ij}), ĉ_j(x) = (A⁻¹)_{j1}+(A⁻¹)_{j2}x
# is its numerator polynomial evaluated at the branch point, and f'(e_m) =
# lead·∏_{n≠m}(e_m−e_n) is the derivative of the sextic at that branch point.
# The constant 4πi was pinned numerically to 48 digits against a certified finite
# difference and matches the classical normalisation (the residue of dx/y at a
# branch point is the source of the 2π, the i is the Riemann-bilinear period
# phase).  This formula is EXACT and certified — the only enclosures are the
# already-certified periods (via A⁻¹) and the algebraic f'(e_m).
#
# For a kinematic family the chain rule gives
#     ∂τ_{jk}/∂t = Σ_m (∂τ_{jk}/∂e_m) · (de_m/dt),
# with the branch-point velocities de_m/dt supplied exactly when the curve is a
# coefficient family (implicit differentiation of f(e_m,t)=0:
# de_m/dt = −(∂_t f)(e_m)/f'(e_m)), or by a certified finite difference on the
# root vector otherwise.

"""
    siegel_rauch_de(C::HyperellipticCurve, m::Int) → AcbMatrix (2×2)

The Rauch variational derivative  ∂τ/∂e_m  of the genus-2 period matrix with
respect to the m-th branch point, as a symmetric 2×2 `AcbMatrix`:

    ∂τ_{jk}/∂e_m = 4πi · ĉ_j(e_m) ĉ_k(e_m) / f'(e_m),

with ω̂ = A⁻¹ω the normalised holomorphic differentials.  Real-ordered curves
only (uses `hyper_big_period_matrix`).  EXACT — no finite differencing.
"""
function siegel_rauch_de(C::HyperellipticCurve, m::Int; rtol = 0.0)
    C.real_ordered ||
        error("siegel_rauch_de: real-ordered branch points required")
    prec = C.prec; g = 2; e = C.roots; n = C.degree
    1 ≤ m ≤ n || error("siegel_rauch_de: branch index out of range")
    A, _ = hyper_big_period_matrix(C; rtol = rtol)
    Ainv = AcbMatrix(g, g; prec); Arblib.inv!(Ainv, A)
    em = e[m]
    ĉ = [Ainv[j, 1] + Ainv[j, 2] * em for j in 1:g]   # ω̂_j numerator at e_m
    # f'(e_m) = lead · ∏_{n≠m}(e_m − e_n)   (degree-5: ∞ branch point contributes
    # no finite-derivative factor; the moving-finite-point formula is unchanged).
    fp = Acb(C.lead; prec)
    for nn in 1:n
        nn == m && continue
        fp *= (em - e[nn])
    end
    fourpii = Acb(0, 4 * Arb(π; prec); prec)
    D = AcbMatrix(g, g; prec)
    for j in 1:g, k in 1:g
        D[j, k] = fourpii * ĉ[j] * ĉ[k] / fp
    end
    return D
end

"""
    branch_velocities_fd(curve_at, t; h, prec) → Vector{Acb} (length n)

Certified branch-point velocities de_m/dt of a curve family `curve_at(t)` by a
symmetric finite difference on the ordered root vector.  Truncation error is
NOT enclosed by this helper (the FD is the non-exact link); prefer the exact
coefficient-family path `branch_velocities_coeffs` when f is polynomial in t.
"""
function branch_velocities_fd(curve_at, t; h = nothing, prec::Int)
    hh = h === nothing ? Arb(2; prec)^(-div(prec, 3)) : Arb(h; prec)
    Cp = curve_at(Acb(t; prec) + hh)
    Cm = curve_at(Acb(t; prec) - hh)
    n = Cp.degree
    n == Cm.degree || error("branch_velocities_fd: degree changed across step")
    return [(Cp.roots[m] - Cm.roots[m]) / (2hh) for m in 1:n]
end

"""
    branch_velocities_coeffs(coeffs_at, roots, t; prec) → Vector{Acb}

EXACT branch-point velocities for a coefficient family f(x,t)=Σ cᵢ(t) xⁱ via
implicit differentiation of f(e_m(t), t)=0:

    de_m/dt = −(∂_t f)(e_m) / f'(e_m),   f'(x)=∂_x f.

`coeffs_at(t)` returns the coefficient vector [c₀,…,c_d]; `roots` are the (already
certified) branch points at this t (e.g. `C.roots`).  `∂_t f` is taken by a
certified central difference of the COEFFICIENTS (rationals/algebraic → exact to
working precision; the x-polynomial evaluation is exact), so the only inexactness
is in ∂_t cᵢ when the cᵢ are themselves transcendental in t.
"""
function branch_velocities_coeffs(coeffs_at, roots::Vector, t; prec::Int,
                                  h = nothing)
    hh = h === nothing ? Arb(2; prec)^(-div(prec, 3)) : Arb(h; prec)
    cp = coeffs_at(Acb(t; prec) + hh); cm = coeffs_at(Acb(t; prec) - hh)
    c0 = coeffs_at(Acb(t; prec))
    dtc = [(Acb(cp[i]; prec) - Acb(cm[i]; prec)) / (2hh) for i in eachindex(c0)]
    out = Vector{Acb}(undef, length(roots))
    for (m, em) in enumerate(roots)
        emA = Acb(em; prec)
        # ∂_t f at e_m  = Σ_i (dc_i/dt) e_m^i
        ftm = Acb(0; prec)
        for i in eachindex(dtc); ftm += dtc[i] * emA^(i - 1); end
        # f'(e_m) = Σ_i i c_i e_m^{i-1}
        fpm = Acb(0; prec)
        for i in 2:length(c0); fpm += (i - 1) * Acb(c0[i]; prec) * emA^(i - 2); end
        out[m] = -ftm / fpm
    end
    return out
end

"""
    siegel_gauss_manin(C, de_dt; rtol=0.0) → AcbMatrix (2×2)

Gauss–Manin connection ∂τ/∂t (symmetric 2×2) of the genus-2 period matrix for a
deformation with branch-point velocities `de_dt[m] = de_m/dt`, via Rauch:

    ∂τ/∂t = Σ_m (∂τ/∂e_m) · de_dt[m].

This is the route-E analogue of the elliptic `dτ/dt = W/Ψ₁²` — the kinematic
derivative of the modulus that any `siegel_transport` must integrate against.
EXACT given exact `de_dt` (e.g. from `branch_velocities_coeffs`).
"""
function siegel_gauss_manin(C::HyperellipticCurve, de_dt::Vector; rtol = 0.0)
    prec = C.prec; g = 2
    length(de_dt) == C.degree ||
        error("siegel_gauss_manin: need one velocity per finite branch point")
    D = AcbMatrix(g, g; prec)
    for m in 1:C.degree
        Dm = siegel_rauch_de(C, m; rtol = rtol)
        v = Acb(de_dt[m]; prec)
        for j in 1:g, k in 1:g
            D[j, k] += Dm[j, k] * v
        end
    end
    return D
end

# ---- sector descriptor (parallel to EllipticSector / CYSector) ---------------

"""
    Genus2Sector

Descriptor for one genus-2 hyperelliptic top sector (route E).  Carries:

  * `prec`, `curve_at :: t → HyperellipticCurve` (the kinematic family),
  * `tau_of_t :: t → AcbMatrix` (cached period-matrix map),
  * `base :: Int` — Abel–Jacobi basepoint index.

Construct via `Genus2Sector(prec; curve_at, base)`.  `siegel_tau_map(sec, t)`
returns τ(t); `igusa_invariants(sec, t)` returns (ψ₄,ψ₆,χ₁₀,χ₁₂)(τ(t)).

The transport step `siegel_transport` is the FRONTIER STUB.
"""
struct Genus2Sector
    prec::Int
    curve_at::Function
    tau_of_t::Function
    base::Int
end

function Genus2Sector(prec::Int; curve_at, base::Int = 0)
    τmap = t -> period_matrix(curve_at(Acb(t; prec)))
    b = base == 0 ? curve_at(Acb(1; prec)).degree : base
    Genus2Sector(prec, curve_at, τmap, b)
end

siegel_tau_map(sec::Genus2Sector, t) = sec.tau_of_t(t)
igusa_invariants(sec::Genus2Sector, t) = igusa_invariants(sec.tau_of_t(t), sec.prec)
abel_jacobi(sec::Genus2Sector, t, xp) =
    abel_jacobi(sec.curve_at(Acb(t; prec = sec.prec)), xp; base = sec.base)

"""
    siegel_gauss_manin(sec::Genus2Sector, t; coeffs_at=nothing) → AcbMatrix (2×2)

Gauss–Manin connection ∂τ/∂t at kinematic point `t` for the sector's family.
If `coeffs_at(t)::Vector` (the t-dependent sextic/quintic coefficients) is
supplied, the branch-point velocities are taken EXACTLY via implicit
differentiation (`branch_velocities_coeffs`); otherwise they fall back to the
certified-FD-on-roots path (`branch_velocities_fd`).  This is the kinematic
derivative the (still-stubbed) `siegel_transport` integrates against.
"""
function siegel_gauss_manin(sec::Genus2Sector, t; coeffs_at = nothing)
    prec = sec.prec
    C = sec.curve_at(Acb(t; prec))
    de = coeffs_at === nothing ?
        branch_velocities_fd(sec.curve_at, t; prec = prec) :
        branch_velocities_coeffs(coeffs_at, C.roots, t; prec = prec)
    return siegel_gauss_manin(C, de)
end

# ---- frontier stub -----------------------------------------------------------

struct SiegelNotImplemented <: Exception
    msg::String
end
Base.showerror(io::IO, e::SiegelNotImplemented) =
    print(io, "SiegelNotImplemented: ", e.msg)

"""
    siegel_transport(sec::Genus2Sector, source, t_start, t_end, prec; …)

Genus-2 analogue of `eichler_transport` / `cy_transport`: would transport the
inhomogeneous solution of the ε-factorised DE along the kinematic path, with
kernels = Siegel modular forms (entries of the connection are weight-k forms in
the (ψ₄,ψ₆,χ₁₀,χ₁₂,χ₃₅) ring, cf. Duhr–Porkert–Stawinski 2412.02300).

NOT IMPLEMENTED.  Of the three ingredients an `eichler_transport` needs, two are
now in place and one is the wall:

  (1) the MODULUS map  t ↦ τ(t) ∈ ℍ₂  — `siegel_tau_map` (certified, ≥60 d);
  (2) the KINEMATIC DERIVATIVE  ∂τ/∂t  (the genus-2 analogue of dτ=(W/Ψ₁²)dt) —
      `siegel_gauss_manin` via Rauch's variational formula (EXACT, certified;
      this is NEW, see that function); together with the modular-form generators
      ψ₄,ψ₆,χ₁₀,χ₁₂,θ[ab] (`igusa_invariants`, `siegel_theta`), so the integrand
      of any Siegel–Eichler integral  ∫ S(τ(t)) (∂τ/∂t) dt  is fully evaluable;
  (3) the ITERATED-INTEGRAL FUNCTION SPACE — a closed, PSLQ-ready basis of
      iterated integrals of Siegel modular forms over ℍ₂ (the genus-2 analogue
      of the Γ₁(N) Eichler/eMPL basis with its cusp/tangential-basepoint
      regularisation and period-polynomial bookkeeping).  THIS DOES NOT YET
      EXIST in usable form — it is the route-E wall.

So the gap is now SHARP: it is NOT "we cannot evaluate the kernels" (we can, to
hundreds of digits) and NOT "we cannot differentiate τ(t)" (Rauch gives ∂τ/∂t
exactly).  It is the absence of a *basis*: there is no genus-2 analogue of
`iterated_integral`/`LogQSeries` that (a) closes under the connection, (b) has a
regularised value at the ℍ₂ cusps (the Fourier–Jacobi / Siegel-Φ boundary), and
(c) supplies period-polynomial constants for the `vop_fit_constants` step.  The
two visible closure paths, both at the "worked examples, not a library" stage:
  • Schottky–Kronecker hyperelliptic polylogarithms (Broedel et al,
    refs.bib `SchottkyKronecker`) — a g=2 Kronecker–Eisenstein generating series;
  • Brown's multiple modular values / iterated Siegel integrals.
The MINIMAL unblock is a depth-1, single-form primitive ∫_{i∞}^{τ} ψ_k dτ over
ℍ₂ with a certified Fourier–Jacobi cusp regularisation — the genus-2 counterpart
of `eichler_transport`'s depth-1 q-series path; everything in this file is
arranged to feed exactly that primitive once its function space is fixed.

References: arXiv:2307.11497 (genus drop), 2412.02300 (genus-2 ε-form),
SchottkyKronecker (refs.bib), and Paper/sections/scope-outlook.tex Rung 4.
"""
function siegel_transport(sec::Genus2Sector, source, t_start, t_end, prec::Int; kwargs...)
    throw(SiegelNotImplemented(
        "genus-2 iterated Siegel-modular integral (the route-E `eichler_transport`): " *
        "kernels (ψ₄,ψ₆,χ₁₀,χ₁₂,θ[ab]) are available via `igusa_invariants` / " *
        "`siegel_theta`, but a PSLQ-ready hyperelliptic-polylog / Siegel-Eichler " *
        "function space is not — see 2307.11497, 2412.02300, SchottkyKronecker, " *
        "and Paper/sections/scope-outlook.tex Rung 4."))
end

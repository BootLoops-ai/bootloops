# Post-Minkowskian period rings — the four CY/K3 Picard–Fuchs operators of
# 2401.07899 (Klemm–Nega–Sauer–Plefka) that enter 4PM and 5PM black-hole
# scattering.  Thin wrapper over the generic MUM-Frobenius engine in
# `cy_transport.jl`; this file only encodes the four operators (in θ-form),
# their closed-form holomorphic periods where known, and the small θ→∂
# conversion needed to build a `PFOperator` from a θ-polynomial.
#
# References (eq. numbers are 2401.07899v3):
#   K3.0  Legendre Sym²   eq.(8)   — 4PM-1SF, all five 5PM-1SF K3.i
#   CY3   ₄F₃ hypergeom   eq.(10)  — 5PM-1SF top.3 (dissipative)
#   K3′   Apéry / Beukers eq.(14)  — 5PM-2SF top.40
#   CY3′  Hadamard χ=80   eq.(17)  — 5PM-2SF top.37
#
# All four are MUM at z = 0 with integer indicial exponent α = 0 *after* the
# z^{γ₁} prefactor is stripped (for K3.0 the paper's γ₁ = ½, so we conjugate
# eq.(8) by z^{1/2} once and record `alpha_shift = 1//2`; the other three have
# γ₁ = 0 already).  The engine in `mum_frobenius_basis` rejects non-integer α,
# so this conjugation is done analytically here, not in the engine.

# ---- θ-form → ∂-form (PFOperator) -------------------------------------------

"""
    pf_from_theta(prec, R) → PFOperator

Build a `PFOperator` from a θ-form  L = Σ_{s=0}^{S} z^s R_s(θ),  θ = z d/dz,
given as `R[s+1][j+1] = [θ^j] R_s :: Rational`.  Uses θ^j = Σ_i S₂(j,i) z^i ∂^i
(Stirling-2nd-kind), so `p_i(z) = Σ_{s,j} R[s+1][j+1] S₂(j,i) z^{s+i}`.

Round-trips: `theta_form(pf_from_theta(prec, R)) == R` (up to trailing-zero
trimming), since `theta_form` removes the global `s_min` shift and all four PM
operators have R_0 ≠ 0 (s_min = 0).
"""
function pf_from_theta(prec::Int, R::Vector{<:Vector})
    S = length(R) - 1
    r = maximum(length(Rs) for Rs in R) - 1
    Rq = [vcat(Rational{BigInt}.(Rs), zeros(Rational{BigInt}, r + 1 - length(Rs))) for Rs in R]
    coeffs = [zeros(Rational{BigInt}, S + r + 1) for _ in 0:r]
    for s in 0:S, j in 0:r
        c = Rq[s+1][j+1]
        c == 0 && continue
        for i in 0:j
            coeffs[i+1][s+i+1] += c * Rational{BigInt}(_stirling2(j, i))
        end
    end
    for p in coeffs
        while length(p) > 1 && p[end] == 0; pop!(p); end
    end
    PFOperator(prec, coeffs)
end

"""Exact ℚ holomorphic-series coefficients c_0..c_N of ϖ_0 = Σ c_m z^m for the
θ-form operator `R` (assumed MUM with α = 0, so R_0(θ) = lc·θ^r)."""
function holo_series_coeffs(R::Vector{Vector{Rational{BigInt}}}, N::Int)
    S = length(R) - 1
    pv(p, x) = sum(p[k+1] * x^k for k in 0:length(p)-1; init = zero(Rational{BigInt}))
    lc = R[1][end]
    c = Vector{Rational{BigInt}}(undef, N + 1); c[1] = Rational{BigInt}(1)
    for m in 1:N
        acc = zero(Rational{BigInt})
        for s in 1:min(S, m)
            acc += pv(R[s+1], Rational{BigInt}(m - s)) * c[m-s+1]
        end
        c[m+1] = -acc / pv(R[1], Rational{BigInt}(m))
    end
    c
end

# ---- PMPeriod descriptor -----------------------------------------------------

"""
    PMPeriod

One PM Picard–Fuchs sector.  Carries:
  * `name`, `theta_coeffs` (the θ-form `R[s+1][j+1]`), `L::PFOperator`;
  * `frob::MUMFrobenius` (the full log-tower at z = 0, reused on every call);
  * `alpha_shift` — the rational z-exponent stripped from the *paper's* operator
    to land at α = 0 (so the physical period is `z^alpha_shift · ϖ_0`);
  * `holo_coeffs` — exact ℚ coefficients of ϖ_0;
  * `holo_closed` — independent closed-form ϖ_0(z) (Arb hypergeom / direct sum),
    or `nothing` if no elementary closed form is known.
"""
struct PMPeriod
    name::String
    theta_coeffs::Vector{Vector{Rational{BigInt}}}
    L::PFOperator
    frob::MUMFrobenius
    alpha_shift::Rational{Int}
    holo_coeffs::Vector{Rational{BigInt}}
    holo_closed::Union{Nothing,Function}
end

"""Period vector (ϖ_0,…,ϖ_{r−1}) at `z`, to (at most) `prec` bits.  `prec` must
not exceed the precision the sector was built at."""
frobenius_basis(pm::PMPeriod, z; prec::Int = pm.L.prec, path = nothing) =
    period_vector(pm.L, z, prec; frob = pm.frob, path = path)

"""Full Wronskian W[i+1,k+1] = ϖ_k^{(i)}(z)."""
frobenius_wronskian(pm::PMPeriod, z; prec::Int = pm.L.prec, path = nothing) =
    period_matrix(pm.L, z, prec; frob = pm.frob, path = path)

"""Residual L·ϖ_0(z) = Σ_i p_i(z) ∂^i ϖ_0(z) (should enclose 0; magnitude
reports the series-truncation floor)."""
function pf_residual(pm::PMPeriod, z)
    r = pm.L.order
    W = frobenius_wronskian(pm, z)
    s = Acb(0; prec = pm.L.prec)
    for i in 0:r
        s += _pcoeff(pm.L, i, z) * (i < r ? W[i+1, 1] :
             # ∂^r ϖ_0 from the ODE itself — but that's circular; instead bump
             # the Frobenius eval to one extra derivative:
             _frob_eval(pm.frob, 0, r, Acb(z; prec = pm.L.prec)))
    end
    s
end

"""Closed-form ϖ_0(z) via the independent route stored on the sector."""
holo_period_closed(pm::PMPeriod, z; prec::Int = pm.L.prec) =
    pm.holo_closed === nothing ? nothing : pm.holo_closed(Acb(z; prec = prec), prec)

# small Arb pFq wrapper (a,b are tuples/vectors of rationals)
function _pfq(a, b, z::Acb, prec::Int)
    av = AcbVector(length(a); prec); bv = AcbVector(length(b); prec)
    for (i, x) in enumerate(a); av[i] = Acb(x; prec); end
    for (i, x) in enumerate(b); bv[i] = Acb(x; prec); end
    res = Acb(prec = prec)
    Arblib.hypgeom_pfq!(res, av, length(a), bv, length(b), z, 0, prec)
    res
end

# ---- the four operators ------------------------------------------------------

"""
    pm_k3_legendre(prec; nterms = 300) → PMPeriod

3-loop / 4PM-1SF K3 (Sym² Legendre, Γ₀(4)).  Variable z = x².  Paper form
eq.(8) is (2θ−1)³ + z²(2θ+1)³ − 4zθ(4θ²+1) with γ₁ = ½; conjugating by z^{1/2}
(θ → θ+½) and dividing by 8 gives the α = 0 operator stored here:

    L = θ³  −  z (2θ³ + 3θ² + 2θ + ½)  +  z² (θ+1)³,

with ϖ_0(z) = ₂F₁(½,½;1;z)² = (2/π)² K(z)² (Legendre-K parameter convention,
`Arblib.elliptic_k!`).  Singularities {0, 1, ∞}; Frobenius radius R = 1.
"""
function pm_k3_legendre(prec::Int; nterms::Int = 300)
    R = Vector{Rational{BigInt}}[[0, 0, 0, 1],
                                 [-1//2, -2, -3, -2],
                                 [1, 3, 3, 1]]
    L = pf_from_theta(prec, R)
    frob = mum_frobenius_basis(L, prec; nterms = nterms)
    holo = holo_series_coeffs(R, nterms)
    closed = (z::Acb, p::Int) -> begin
        K = Acb(prec = p); Arblib.elliptic_k!(K, z)
        (Acb(2; prec = p) / Acb(π; prec = p))^2 * K^2
    end
    PMPeriod("K3.0 (Legendre Sym²)", R, L, frob, 1//2, holo, closed)
end

"""
    pm_cy3_4F3(prec; nterms = 300) → PMPeriod

4-loop / 5PM-1SF CY3 (eq.10): L = θ⁴ − 2⁸ z (θ+½)⁴ in the variable z = x⁴/2⁸.
ϖ_0(z) = ₄F₃(½,½,½,½; 1,1,1; 2⁸ z); since (½)_n/n! = C(2n,n)/4ⁿ the series
coefficients are c_n = C(2n,n)⁴.  Singularities {0, 2⁻⁸, ∞}; R = 2⁻⁸.
"""
function pm_cy3_4F3(prec::Int; nterms::Int = 300)
    R = Vector{Rational{BigInt}}[[0, 0, 0, 0, 1],
                                 [-16, -128, -384, -512, -256]]
    L = pf_from_theta(prec, R)
    frob = mum_frobenius_basis(L, prec; nterms = nterms)
    holo = holo_series_coeffs(R, nterms)
    closed = (z::Acb, p::Int) ->
        _pfq((1//2, 1//2, 1//2, 1//2), (1, 1, 1), Acb(256; prec = p) * z, p)
    PMPeriod("CY3 (₄F₃ hypergeometric)", R, L, frob, 0//1, holo, closed)
end

"""Exact Apéry numbers A005259: a_0..a_N via the 3-term recurrence
n³ a_n = (2n−1)(17n²−17n+5) a_{n−1} − (n−1)³ a_{n−2}."""
function apery_numbers(N::Int)
    a = Vector{BigInt}(undef, N + 1); a[1] = big(1); N ≥ 1 && (a[2] = big(5))
    for n in 2:N
        a[n+1] = div((2n - 1) * (17 * big(n)^2 - 17n + 5) * a[n] - big(n - 1)^3 * a[n-1],
                     big(n)^3)
    end
    a
end

"""
    pm_k3_apery(prec; nterms = 300) → PMPeriod

4-loop / 5PM-2SF top.40 K3′ — the Beukers–Peters Apéry operator (eq.14):

    L = θ³  +  z² (θ+1)³  −  z (2θ+1)(17θ²+17θ+5),

variable z = x².  ϖ_0 = Σ_n a_n zⁿ with a_n = Σ_k C(n,k)² C(n+k,k)² (A005259,
the ζ(3)-irrationality Apéry numbers).  Singularities {0, (3∓2√2)², ∞}; the
nearest is z₊ = (3−2√2)² ≈ 0.02944 — **inside** the unit disk and inside the
physical region 0 < x < 1 (paper's "γ = 3, PHYSICAL" point), so the MUM
Frobenius series alone covers only |z| < z₊.  Use `period_matrix(…; path=…)`
(θ-companion transport) to continue past it.
"""
function pm_k3_apery(prec::Int; nterms::Int = 300)
    R = Vector{Rational{BigInt}}[[0, 0, 0, 1],
                                 [-5, -27, -51, -34],
                                 [1, 3, 3, 1]]
    L = pf_from_theta(prec, R)
    frob = mum_frobenius_basis(L, prec; nterms = nterms)
    holo = holo_series_coeffs(R, nterms)
    A = apery_numbers(nterms)
    closed = (z::Acb, p::Int) -> begin
        # direct Apéry-recurrence sum (independent of holo_series_coeffs);
        # tail: |a_{N+1} z^{N+1}|/(1 − |z|/z₊) since a_{n+1}/a_n → 1/z₊.
        N = length(A) - 1
        s = Acb(0; prec = p); zp = Acb(1; prec = p)
        for n in 0:N
            s += A[n+1] * zp; zp *= z
        end
        zplus = (Arb(3; prec = p) - 2 * sqrt(Arb(2; prec = p)))^2
        θ = abs(z) / zplus
        last = Arb(abs(A[N+1]); prec = p) * abs(zp)   # |a_N z^{N+1}| ~ next term scale
        Arblib.add_error!(s, last * (1 / zplus) / (1 - θ))   # crude geometric envelope
        s
    end
    PMPeriod("K3′ (Apéry / Beukers–Peters)", R, L, frob, 0//1, holo, closed)
end

"""
    pm_cy3_hadamard(prec; nterms = 300) → PMPeriod

4-loop / 5PM-2SF top.37 CY3′ (eq.17), variable z = −t²/(2¹²(1+t)), t = x²−1:

    L = θ⁴ − 2⁴ z (192θ⁴+128θ³+112θ²+48θ+7)
            + 2¹⁴ z² (192θ⁴+256θ³+208θ²+64θ+7)
            − 2³⁰ z³ (θ+½)⁴.

Singularities {0, 2⁻¹⁰ (triple), ∞}; R = 2⁻¹⁰.  ϖ_0 series: 1, 112, 47376,
27846400, … (verified against the paper's Baikov T³ residue eq.(A6) via the
t-reexpansion 1 − 7/2⁸·t² + 7/2⁸·t³ − 25711/2²⁰·t⁴ + …).

No elementary closed form.  The independent oracle stored in `holo_closed` is
the **Hadamard identity** (paper, App. A; van Straten case 2.33):
    ϖ_0(z) = (1 + 2¹⁰ z)^{−1/2} · ϖ_0^{Had}( z / (1+2¹⁰z)² ),
with L_Had = θ⁴ + 2¹⁶ w² (4θ+1)(4θ+3)(4θ+5)(4θ+7) − 2⁴ w (4θ+1)(4θ+3)(32θ²+32θ+13).
This is a *different* PF operator solved by the same MUM engine — a non-trivial
cross-operator identity, not a self-check.
"""
function pm_cy3_hadamard(prec::Int; nterms::Int = 300)
    R = Vector{Rational{BigInt}}[[0, 0, 0, 0, 1],
                                 Rational{BigInt}.(-16 .* [7, 48, 112, 128, 192]),
                                 Rational{BigInt}.(2^14 .* [7, 64, 208, 256, 192]),
                                 Rational{BigInt}.(-(big(2)^30) .* [1//16, 1//2, 3//2, 2, 1])]
    L = pf_from_theta(prec, R)
    frob = mum_frobenius_basis(L, prec; nterms = nterms)
    holo = holo_series_coeffs(R, nterms)
    # Hadamard companion operator (App. A)
    Rhad = Vector{Rational{BigInt}}[[0, 0, 0, 0, 1],
                                    Rational{BigInt}.(-16 .* [39, 304, 816, 1024, 512]),
                                    Rational{BigInt}.(2^16 .* [105, 704, 1376, 1024, 256])]
    Lhad = pf_from_theta(prec, Rhad)
    frob_had = mum_frobenius_basis(Lhad, prec; nterms = nterms)
    closed = (z::Acb, p::Int) -> begin
        u = Acb(1; prec = p) + Acb(1024; prec = p) * z
        w = z / u^2
        ϖhad = period_vector(Lhad, w, p; frob = frob_had)[1]
        ϖhad / sqrt(u)
    end
    PMPeriod("CY3′ (Hadamard χ=80)", R, L, frob, 0//1, holo, closed)
end

# convenience: dispatch all four
pm_all_periods(prec::Int; nterms::Int = 300) =
    (pm_k3_legendre(prec; nterms), pm_cy3_4F3(prec; nterms),
     pm_k3_apery(prec; nterms), pm_cy3_hadamard(prec; nterms))

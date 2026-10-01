# Constant dictionaries for PSLQ. All entries certified balls.

# ---- Dirichlet characters (Kronecker symbol) ---------------------------------

"""Kronecker symbol (a|n) for arbitrary integers a, n."""
function kronecker_symbol(a::Integer, n::Integer)
    a = Int(a); n = Int(n)
    n == 0 && return (abs(a) == 1 ? 1 : 0)
    res = 1
    if n < 0
        n = -n
        a < 0 && (res = -res)
    end
    while iseven(n)
        n >>= 1
        iseven(a) && return 0
        r = mod(a, 8)
        (r == 3 || r == 5) && (res = -res)
    end
    a = mod(a, n)
    while a != 0
        while iseven(a)
            a >>= 1
            r = mod(n, 8)
            (r == 3 || r == 5) && (res = -res)
        end
        a, n = n, a
        (mod(a, 4) == 3 && mod(n, 4) == 3) && (res = -res)
        a = mod(a, n)
    end
    return n == 1 ? res : 0
end

"""Squarefree test (trial division — only used to validate small conductors)."""
function _is_squarefree(n::Int)
    n = abs(n)
    p = 2
    while p * p <= n
        n % (p * p) == 0 && return false
        p += 1
    end
    return true
end

"""Is D a fundamental discriminant?"""
function is_fundamental_discriminant(D::Int)
    D == 1 && return false
    if mod(D, 4) == 1
        return _is_squarefree(D)
    elseif mod(D, 4) == 0
        m = div(D, 4)
        return (mod(m, 4) == 2 || mod(m, 4) == 3) && _is_squarefree(m)
    end
    return false
end

"""Unique primitive odd real character mod N (Kronecker (−N | ·)), as the value
vector [χ(1), …, χ(N)]. Requires −N to be a fundamental discriminant; supported
in particular for N ∈ {3,4,7,8,11,15,20,24}."""
function chi_minus(N::Int)
    is_fundamental_discriminant(-N) ||
        error("chi_minus: -$N is not a fundamental discriminant (no primitive odd real χ mod $N)")
    return Int[kronecker_symbol(-N, a) for a in 1:N]
end

# ---- Dirichlet L-function (certified via Hurwitz zeta / digamma) -------------

"""L(χ, s) for a primitive real character χ given as the value vector [χ(1)…χ(N)]
and integer s ≥ 1, certified to `prec` bits.

  s ≥ 2 : L(χ,s) = N^{-s} Σ_{a=1}^N χ(a) ζ(s, a/N)   (Arblib.hurwitz_zeta!)
  s = 1 : L(χ,1) = -(1/N) Σ_{a=1}^N χ(a) ψ(a/N)       (digamma; pole cancels for Σχ=0)

Arblib does not currently wrap `acb_dirichlet_l` (needs the dirichlet_group/char
C structs), so we use the Hurwitz route, which is ball-certified."""
function dirichlet_L(chi::Vector{Int}, s::Int, prec::Int)
    N = length(chi)
    if s == 1
        sum(chi) == 0 || error("dirichlet_L: s=1 requires a non-principal character")
        acc = Arb(0; prec)
        for a in 1:N
            chi[a] == 0 && continue
            d = Acb(prec = prec); Arblib.digamma!(d, Acb(a // N; prec))
            acc -= chi[a] * real(d)
        end
        return acc / N
    end
    sA = Acb(s; prec)
    acc = Acb(0; prec)
    for a in 1:N
        chi[a] == 0 && continue
        z = Acb(prec = prec); Arblib.hurwitz_zeta!(z, sA, Acb(a // N; prec))
        acc += chi[a] * z
    end
    return real(Acb(N; prec)^(-s) * acc)
end

"""L(χ₋₃, s) — kept for backward compatibility (now a thin wrapper)."""
dirichlet_L_chim3(s::Int, prec::Int) = dirichlet_L(chi_minus(3), s, prec)

# ---- polylog / Clausen --------------------------------------------------------

"""Li_n at the primitive `root_order`-th root of unity e^{2πi/root_order} (certified)."""
function polylog_at_root(n::Int, root_order::Int, prec::Int)
    z = exp(2 * Acb(0, π; prec) / root_order)
    res = Acb(prec = prec)
    Arblib.polylog!(res, Acb(n; prec), z)
    return res
end

"""Clausen function Cl_w(θ) = Im Li_w(e^{iθ}) (certified)."""
function clausen(w::Int, theta::Arb, prec::Int)
    z = exp(Acb(0, theta; prec))
    res = Acb(prec = prec); Arblib.polylog!(res, Acb(w; prec), z)
    return imag(res)
end

# ---- cusp dictionary ----------------------------------------------------------

"""Cusp dictionary: Dirichlet-L / root-of-unity polylog ring generators up to
`maxweight`, for the given primitive odd-character `conductors`. Returns
Vector{Pair{String,Acb}}. Default conductors are [3,4] (Γ₁(6) ∪ Γ₀(4))."""
function cusp_dictionary(prec::Int; maxweight::Int = 4, conductors::Vector{Int} = [3, 4])
    π_ = Arb(π; prec)
    out = Pair{String,Acb}[]
    push!(out, "1" => Acb(1; prec))
    for w in 1:maxweight
        push!(out, "pi^$w" => Acb(π_^w; prec))
    end
    push!(out, "log2" => Acb(log(Arb(2; prec)); prec))
    push!(out, "log3" => Acb(log(Arb(3; prec)); prec))
    for N in conductors
        N in (2, 3) && continue
        push!(out, "log$N" => Acb(log(Arb(N; prec)); prec))
    end
    for w in 2:maxweight
        z = Arb(prec = prec); Arblib.zeta!(z, Arb(w; prec))
        push!(out, "zeta$w" => Acb(z; prec))
    end
    for N in conductors
        chi = chi_minus(N)
        for w in 2:maxweight
            push!(out, "L(chi-$N,$w)" => Acb(dirichlet_L(chi, w, prec); prec))
        end
    end
    if 4 in conductors
        push!(out, "Catalan" => Acb(dirichlet_L(chi_minus(4), 2, prec); prec))
    end
    for w in 2:maxweight
        li = polylog_at_root(w, 3, prec)
        push!(out, "ImLi$w(r3)" => Acb(imag(li); prec))
        push!(out, "ReLi$w(r3)" => Acb(real(li); prec))
        li6 = polylog_at_root(w, 6, prec)
        push!(out, "ImLi$w(r6)" => Acb(imag(li6); prec))
    end
    for N in conductors, w in 2:maxweight
        push!(out, "Cl$w(pi/$N)" => Acb(clausen(w, π_ / N, prec); prec))
    end
    return out
end

# ---- weight-2/3 Eisenstein cusp constants for Γ₀(N) --------------------------

# Cusp data: cusp => (width, constant Fourier coefficient of E₂^{(N)} = E₂(τ)-N·E₂(Nτ)
# in the local q_h expansion after slashing). Derived from
#   E₂^{(N)} |₂ γ_c (τ)  =  E₂(τ) - (N/c²)·E₂((Nτ+const)/c²)       (γ_c: i∞ ↦ a/c, c|N)
# whose constant term as τ→i∞ is 1 - N/c². See Diamond–Shurman §4.6, Cohen–Strömberg §5.
const _E2N_CUSPDATA = Dict(
    3 => Dict(1//0 => (1, -2//1),   0//1 => (3,  2//3)),
    4 => Dict(1//0 => (1, -3//1),   0//1 => (4,  3//4),   1//2 => (1,  0//1)),
    6 => Dict(1//0 => (1, -5//1),   0//1 => (6,  5//6),   1//2 => (3,  1//3),  1//3 => (2, -1//2)),
)

# ---- general cusp-constant via slash + certified eta evaluation --------------

"""SL₂(ℤ) matrix (a, b; c, d) with first column (a, c), i.e. γ·∞ = a/c, for a
cusp given as a `Rational` (use `1//0` for i∞)."""
function _cusp_matrix(cusp::Rational{Int})
    a, c = numerator(cusp), denominator(cusp)
    c == 0 && return (1, 0, 0, 1)
    g, d, mb = gcdx(a, c)          # d·a + mb·c = g
    g == 1 || error("_cusp_matrix: cusp $cusp not in lowest terms")
    return (a, -mb, c, d)
end

"""
    cusp_constant_E_k_chi(k, rs, cusp, prec; width, y_mult = 1)

Constant Fourier coefficient a₀ of the weight-`k` modular form given by the
eta quotient `rs = [(d₁,r₁),…]` (so f(τ) = ∏ η(dᵢτ)^{rᵢ}) at the cusp `cusp`
of width `width`, as a certified `Acb` ball.

Algorithm (Diamond–Shurman Prop. 4.6.5 / Miyake Thm. 7.1.3, realised
numerically): pick γ ∈ SL₂(ℤ) with γ·∞ = cusp and evaluate
    (f|ₖγ)(τ) = (cτ+d)^{-k} · f(γτ)
at τ = iy with y large.  Since f|ₖγ ∈ Mₖ(Γ(N)) is holomorphic on |q_h| < 1
(q_h = e^{2πiτ/width}), it equals a₀ + O(q_h).  We take
    y₂ = width · prec · log 2 / π · y_mult,   y₁ = y₂/2,
so |q_h(y₂)| ≤ 2^{-2·prec·y_mult}; the returned ball is v(y₂) inflated by the
observed |v(y₁) - v(y₂)| (which dominates the residual tail at y₂ for any form
holomorphic at the cusp; for the eta-quotient Eisenstein series tabulated below
this has been verified against the exact closed forms).

The eta evaluations themselves are ball-certified for arbitrary τ ∈ ℍ via
`Arblib.modular_eta!` (which reduces to the fundamental domain internally), so
no slow q-series summation at small Im γτ is needed.
"""
function cusp_constant_E_k_chi(k::Int, rs::Vector{Tuple{Int,Int}}, cusp::Rational{Int},
                               prec::Int; width::Int, y_mult::Real = 1)
    a, b, c, d = _cusp_matrix(cusp)
    a*d - b*c == 1 || error("cusp_constant_E_k_chi: bad cusp matrix")
    yA = Arb(width; prec) * prec * log(Arb(2; prec)) / Arb(π; prec) * y_mult
    vals = map((yA / 2, yA)) do y
        τ = Acb(Arb(0; prec), y; prec)
        z = (a*τ + b) / (c*τ + d)
        (c*τ + d)^(-k) * eta_quotient_point(z, rs)
    end
    err = Mag(); Arblib.get!(err, vals[1] - vals[2])
    out = Acb(vals[2]); Arblib.add_error!(out, err)
    return out
end

# Eta-quotient bases for the weight-3 Eisenstein space on Γ₀(f), character χ_{-f}
# (f ∈ {3,4}).  :A is the a₀(∞) = 1 element, :B the a₀(∞) = 0 element; the
# Eisenstein series in the L(χ,1-k) normalisation are
#     E₃(χ,𝟏) = L(χ,-2)·A   (a₀(∞) = -B_{3,χ}/3),     E₃(𝟏,χ) = B   (a₀(∞) = 0).
# B for f=3 coincides with `eisenstein_e3_qseries`.
const _E3chi_ETAREP = Dict(
    3 => Dict(:A => [(1,9),(3,-3)],          :B => [(3,9),(1,-3)]),
    4 => Dict(:A => [(1,4),(2,6),(4,-4)],    :B => [(1,-4),(2,6),(4,4)]),
)

# Hard-tabled cusp constants of (A, B) at every cusp of Γ₀(N), N ∈ {3,4,6}
# (χ = χ_{-3} for N ∈ {3,6}, χ = χ_{-4} for N = 4).  Values verified against
# `cusp_constant_E_k_chi` to > 200 digits.  Stored as `cusp => (width, a₀(A),
# a₀(B))` with each a₀ a closure `prec → Acb` so that irrational constants
# (i/3^{9/2} etc.) are materialised at the requested precision.
const _E3chi_CUSPDATA = Dict(
    3 => (3, Dict(
        1//0 => (1, p -> Acb(1; prec=p),  p -> Acb(0; prec=p)),
        0//1 => (3, p -> Acb(0; prec=p),  p -> Acb(0, 1; prec=p) / (81 * sqrt(Arb(3; prec=p)))),
    )),
    4 => (4, Dict(
        1//0 => (1, p -> Acb(1; prec=p),  p -> Acb(0; prec=p)),
        0//1 => (4, p -> Acb(0; prec=p),  p -> Acb(0, 1//128; prec=p)),
        1//2 => (1, p -> Acb(0; prec=p),  p -> Acb(0; prec=p)),
    )),
    6 => (3, Dict(
        1//0 => (1, p -> Acb(1; prec=p),  p -> Acb(0; prec=p)),
        0//1 => (6, p -> Acb(0; prec=p),  p -> Acb(0, 1; prec=p) / (81 * sqrt(Arb(3; prec=p)))),
        1//2 => (3, p -> Acb(0; prec=p),  p -> Acb(0, -1; prec=p) / (81 * sqrt(Arb(3; prec=p)))),
        1//3 => (2, p -> Acb(1; prec=p),  p -> Acb(0; prec=p)),
    )),
)

"""Bernoulli numbers B₀ … B_k (convention B₁ = -1/2) as exact `Rational{BigInt}`,
via the defining recursion Σ_{j=0}^{m} C(m+1,j) B_j = 0."""
function _bernoulli_numbers(k::Int)
    k >= 0 || error("_bernoulli_numbers: k must be ≥ 0")
    B = Vector{Rational{BigInt}}(undef, k + 1)
    B[1] = 1
    for m in 1:k
        s = zero(Rational{BigInt})
        for j in 0:m-1
            s += binomial(BigInt(m + 1), BigInt(j)) * B[j+1]
        end
        B[m+1] = -s // (m + 1)
    end
    return B
end

"""Coefficients [c₀, …, c_k] of the Bernoulli polynomial B_k(x) = Σ_j c_j x^{k-j},
c_j = C(k,j) B_j, as exact rationals (so B₁(x) = x - 1/2, B₂(x) = x² - x + 1/6, …)."""
function _bernoulli_poly_coeffs(k::Int)
    B = _bernoulli_numbers(k)
    return [binomial(BigInt(k), BigInt(j)) * B[j+1] for j in 0:k]
end

"""Generalised Bernoulli number B_{k,χ} = f^{k-1} Σ_{a=1}^f χ(a) B_k(a/f), as an
exact Arb (rational) ball, for any k ≥ 1.  The Bernoulli polynomial is built
from the recursion B_k(x) = Σ_j C(k,j) B_j x^{k-j} with the B_j exact rationals
(convention B₁ = -1/2), evaluated by Horner in ball arithmetic — every input is
rational, so the only outward rounding is the final Arb conversion at `prec`.
The k = 2, 3, 4 closed forms are pinned against this recursion in the tests."""
function generalized_bernoulli(k::Int, chi::Vector{Int}, prec::Int)
    k >= 1 || error("generalized_bernoulli: k must be ≥ 1")
    f = length(chi)
    co = _bernoulli_poly_coeffs(k)          # B_k(x) = Σ_j co[j+1] x^{k-j}
    Bk = x -> begin
        acc = Arb(co[1]; prec)
        for j in 1:k
            acc = acc * x + Arb(co[j+1]; prec)
        end
        acc
    end
    s = Arb(0; prec)
    for a in 1:f
        chi[a] == 0 && continue
        s += chi[a] * Bk(Arb(a // f; prec))
    end
    return s * Arb(f; prec)^(k-1)
end

"""Regularised constant Fourier coefficients of the weight-2/3 Eisenstein series
on Γ₀(N) at every cusp, N ∈ {3,4,6}.  Returns a `Dict` keyed as follows:

  * `(2, cusp)::Arb` — constant of the holomorphic E₂^{(N)} = E₂ - N·E₂(N·);
    rational, hard-tabled in `_E2N_CUSPDATA`.  Satisfies the residue relation
    Σ_cusps width·a₀ = 0.

  * `(3, cusp)::Arb` — backward-compatible scalar entry: at the cusps where the
    eta-quotient basis element :A (a₀(∞)=1) is non-vanishing this is
    a₀(E₃(χ,𝟏)) = L(χ,-2)·a₀(A); at the remaining cusps it is the imaginary
    part of a₀(B) (the i is the weight-3 odd-character automorphy phase).

  * `(3, :A, cusp)::Acb`, `(3, :B, cusp)::Acb` — full complex constants of the
    two eta-quotient basis Eisenstein series (`_E3chi_ETAREP`) at every cusp,
    including the intermediate cusps 1/2, 1/3 of Γ₀(6).  These are the
    well-defined, normalisation-free data and are what new code should use.

For N = 6 the character is χ_{-3} lifted mod 6 (χ_{-6} is not primitive)."""
function eisenstein_constants(N::Int, prec::Int; maxweight::Int = 3)
    haskey(_E2N_CUSPDATA, N) ||
        error("eisenstein_constants: only N ∈ {3,4,6} are hard-tabled (got N=$N)")
    out = Dict{Any,Any}()
    # weight 2
    for (cusp, (w, a0)) in _E2N_CUSPDATA[N]
        out[(2, cusp)] = Arb(a0; prec)
    end
    maxweight >= 3 || return out
    haskey(_E3chi_CUSPDATA, N) || return out
    cond, tbl = _E3chi_CUSPDATA[N]
    chi = chi_minus(cond)
    Lchi = -generalized_bernoulli(3, chi, prec) / 3      # L(χ, -2)
    for (cusp, (w, fA, fB)) in tbl
        aA = fA(prec); aB = fB(prec)
        out[(3, :A, cusp)] = aA
        out[(3, :B, cusp)] = aB
        # backward-compatible scalar key (see docstring)
        out[(3, cusp)] = Arblib.contains_zero(aA) ? imag(aB) : Lchi * real(aA)
    end
    # legacy keys (kept exactly for N ∈ {3,4}; (3,1//0) is covered above)
    return out
end

# ---- tangential-basepoint regularization tail (period polynomial) ------------

"""Period-polynomial L-value content of the tangential-basepoint regularization
of ∫_{i∞}^{cusp} f(τ)(τ-X)^{k-2} dτ for f a weight-k Eisenstein series on Γ₀(N)
with character χ₋ₙ. Returns the vector of critical L-values
[L(χ₋ₙ, k-1), L(χ₋ₙ, k-2), …, L(χ₋ₙ, 1)], which (up to explicit (2πi)^j/j! and
cusp-dependent rational cocycle factors that the caller supplies) are exactly
the coefficients of the degree-(k-2) period polynomial (Eichler–Shimura for
Eisenstein series; Kohnen–Zagier 1984, Brown 1707.01230 §7).

`cusp` is currently informational: for Eisenstein series the transcendental
content is cusp-independent — only the rational cocycle changes with the cusp
representative."""
function cusp_reg_tail(N::Int, cusp, weight::Int, prec::Int)
    weight >= 2 || error("cusp_reg_tail: weight must be ≥ 2")
    chi = chi_minus(N)
    return Arb[dirichlet_L(chi, j, prec) for j in (weight-1):-1:1]
end

# ---- interior dictionary (unchanged) -----------------------------------------

"""Interior-point dictionary at kinematic point s (ACKM parent-curve variables,
m_t² = 1): periods, quasi-periods (via s-derivatives), the third-kind period G at
(s,t), and π factors. Returns Vector{Pair{String,Acb}}."""
function interior_dictionary(s::Arb, t::Arb, prec::Int)
    sA = Acb(s; prec = prec); tA = Acb(t; prec = prec)
    ψ0, ψ1 = acm_periods(s)
    _, ψ0p = _acm_psi0_with_derivative(s)   # (ψ0, dψ0/ds) — take the derivative slot
    G = third_kind_G(sA, tA)
    # marked point sits on the arcsin cut for Euclidean (s,t): use the documented
    # exact-dyadic Feynman deformation (entry name records the convention)
    Z = abel_image_deformed(sA, tA; delta_exp = div(prec, 2))
    out = Pair{String,Acb}[]
    push!(out, "1" => Acb(1; prec))
    push!(out, "pi" => Acb(Arb(π; prec); prec))
    push!(out, "psi0" => ψ0)
    push!(out, "psi1" => ψ1)
    push!(out, "dpsi0" => ψ0p)
    push!(out, "G" => G)
    push!(out, "Z_xp(t+i*2^-$(div(prec,2)))" => Z)
    return out
end

# ---- serialization helpers ---------------------------------------------------

"""Decimal string of the midpoint of an Arb ball (PSLQ-ready)."""
function decimal_midpoint(v::Arb, digits::Int)
    old = precision(BigFloat)
    s = try
        setprecision(BigFloat, ceil(Int, digits * 3.33) + 64)
        bf = BigFloat(v)
        Base.MPFR.string_mpfr(bf, "%." * string(digits) * "Re")
    finally
        setprecision(BigFloat, old)
    end
    return s
end

"""Serialize a dictionary to a text file: name, midpoint (decimal), radius."""
function write_dictionary(path::AbstractString, dict::Vector{Pair{String,Acb}}; digits::Int = 0)
    open(path, "w") do io
        for (name, v) in dict
            d = digits == 0 ? floor(Int, Arblib.precision(v) * 0.301) : digits
            println(io, name)
            println(io, "  re_mid = ", decimal_midpoint(real(v), d))
            println(io, "  im_mid = ", decimal_midpoint(imag(v), d))
            println(io, "  rad    = ", Arblib.radius(Arb, real(v)), " ", Arblib.radius(Arb, imag(v)))
        end
    end
end

"""Approximate count of certified decimal digits of a ball (for tests/reporting)."""
function certified_digits_test(v::Arb)
    r = Float64(Arblib.radius(Arb, v))
    m = Float64(v)
    r == 0 && return 9999
    round(Int, -log10(r / max(abs(m), 1e-300)))
end

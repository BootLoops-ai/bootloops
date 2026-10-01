# cuspvals.jl — EXACT constant terms ("cusp values") of the Γ₁(6) Eisenstein
# bases of src/modular.jl at all four cusps of Γ₁(6), as exact elements of the
# number field K = Q(√−3) (s3 = √−3 = i√3, s3² = −3).
#
# Included by the LandauAlphabet module; also loadable standalone after the
# package:
#
#     using LandauAlphabet
#     include(joinpath(pkgroot, "src", "cuspvals.jl"))
#
# It references only Nemo and LandauAlphabet names (exported: B2K,
# eisenstein_E, kronecker_char, level6_basis, sunrise_data, solve_in_basis,
# qs_substitute_negq; qualified: LandauAlphabet.gen_bernoulli,
# LandauAlphabet.char_conductor).
#
# ───────────────────────────────────────────────────────────────────────────
# CONVENTIONS
#
# Cusps of Γ₁(6): representatives a/c with γ = [a b; c d] ∈ SL₂(Z), γ(∞)=a/c:
#     (a,c) ∈ {(1,0)=∞, (0,1), (1,2), (1,3)},   widths {1, 6, 3, 2}.
# Cusp value of a weight-k form f at (a,c) := n=0 Fourier coefficient of
#     (f|_k γ)(τ) = (cτ+d)^{−k} f(γτ),
# i.e. lim_{Im τ→∞} (f|_k γ)(τ).  It is independent of the completion (b,d):
# replacing (b,d) → (b+ja, d+jc) translates τ, which fixes the n=0 mode.
# All q-expansions are in the Γ₁(6) nome q = q_C = e^{2πiτ} of modular.jl.
#
# ───────────────────────────────────────────────────────────────────────────
# 1) B_{2,K} (weight 2, quasi-modular E₂)  — Hermite-decomposition method
#
# With E2std(τ) = 1 − 24 Σ σ₁(n) qⁿ  (so qseries.jl's E₂ = −E2std/24):
#     E2std(γτ) = (cτ+d)² E2std(τ) − (6i/π)·c·(cτ+d),     γ ∈ SL₂(Z).     (∗)
# For the dilated piece write
#     [K·a  K·b; c  d] = γ₁ · [A B; 0 D],   γ₁ = [a₁ b₁; c₁ d₁] ∈ SL₂(Z),
# with A = gcd(K·a, c) > 0 (the gcd of the first column, since the first
# column of γ₁ is primitive), D = K/A, and any B making the factorization
# integral.  Note gcd(K·a, c) = gcd(K, c) because gcd(a,c) = 1.  Comparing
# entries: c₁A = c and c₁B + d₁D = d, hence the key identity
#     j_{γ₁}((Aτ+B)/D) := c₁(Aτ+B)/D + d₁ = (cτ+d)/D.
# Then K·γτ = γ₁((Aτ+B)/D) and (∗) gives
#     (cτ+d)^{−2}·K·E2std(K·γτ) = (A/D)·E2std((Aτ+B)/D) − (6i/π)·c/(cτ+d),
#     (cτ+d)^{−2}·E2std(γτ)     =        E2std(τ)        − (6i/π)·c/(cτ+d),
# so in B_{2,K} = −(1/24)[E2std(τ) − K·E2std(Kτ)] the anomalies cancel:
#     (B_{2,K}|₂γ)(τ) = −(1/24)[E2std(τ) − (A/D)·E2std((Aτ+B)/D)].
# E2std((Aτ+B)/D) = 1 + Σ_{n≥1} c_n e^{2πin(Aτ+B)/D} has n=0 coefficient 1
# (the fractional powers nA/D > 0 never hit 0), hence EXACTLY
#     a₀(B_{2,K} at (a,c)) = (A − D)/(24 D),   A = gcd(K,c), D = K/A.
#
# ───────────────────────────────────────────────────────────────────────────
# 2) E_k(K·τ; φ, ψ), k ≥ 3 — absolutely convergent lattice sums
#
# Characters: φ = Kronecker(Dφ/·) mod L, ψ = Kronecker(Dψ/·) mod M, both real
# (D ∈ {1, −3}), primitive, parity φψ(−1) = (−1)^k.  Conventions: trivial
# character ≡ 1 on ALL integers including 0; χ₋₃(0) = 0.
#
# Twisted Lipschitz formula (ψ primitive mod M, Im z > 0):
#     Σ_{n∈Z} ψ(n)/(z+n)^k = C_{ψ,k} Σ_{r≥1} ψ̄(r) r^{k−1} e^{2πi r z/M},
#     C_{ψ,k} = (−2πi)^k τ(ψ)/(M^k (k−1)!),  Gauss sums τ(1)=1, τ(χ₋₃)=s3.
# For Im z < 0 substitute n → −n: Σ_n ψ(n)/(z+n)^k = ψ(−1)(−1)^k × (same with
# z → −z).  Build the lattice sum (absolutely convergent for k ≥ 3)
#     S(τ) := Σ'_{(m,n)∈Z²} φ(m) ψ(n) / (m·M·τ + n)^k        (origin omitted).
# Row m = 0 (present iff φ(0)=1, i.e. cond φ = 1):
#     Σ_{n≠0} ψ(n)/n^k = (1 + ψ(−1)(−1)^k) L(k,ψ) = 2 L(k,ψ)   [parity].
# Rows m ≠ 0, pairing ±m and using the parity condition:
#     Σ_{m≥1}[φ(m) + (−1)^k φ(−m)ψ(−1)]·C_{ψ,k} Σ_r ψ̄(r) r^{k−1} e^{2πi r m τ}
#       = 2 C_{ψ,k} Σ_{N≥1} (Σ_{d|N} ψ̄(d) φ(N/d) d^{k−1}) q^N.
# Since ψ is real (ψ̄ = ψ) this is Stein's coefficient law, so with
#     c_norm := 2 C_{ψ,k}:           E_k(τ; φ, ψ) = S(τ)/c_norm,
# and the constant term δ_{cond φ=1}·2L(k,ψ)/c_norm = −B_{k,ψ}/(2k) requires
# the classical evaluation (χ primitive mod F, χ(−1) = (−1)^k)
#     L(k, χ) = −C_{χ,k} · B_{k,χ}/(2k),                                  (L)
# which for χ = 1, k = 4 reads ζ(4) = (8π⁴/3)·(1/30)/8 = π⁴/90 ✓ and which is
# verified numerically by the acceptance tests (horocycle averages, test E).
#
# Cusp value: with f(τ) = E_k(K·τ; φ, ψ) and γ = [a b; c d],
#     (f|_k γ)(τ) = (1/c_norm) Σ' φ(m)ψ(n) / ((mMKa + nc)τ + (mMKb + nd))^k.
# As Im τ → ∞ every term with mMKa + nc ≠ 0 vanishes (dominated convergence,
# k ≥ 3), so the constant term is the sub-sum over the lattice line
#     (MKa)·m + c·n = 0  ⇔  (m,n) = t·(m₀, n₀),  t ∈ Z\{0},
#     (m₀, n₀) = (c, −MKa)/g,  g = gcd(MKa, c) > 0.
# On the line the second linear form is t·β with
#     β = m₀MKb + n₀d = (MK/g)(cb − ad) = −MK/g  (an exact integer, b,d-free).
# Complete multiplicativity of Kronecker symbols gives
# φ(m₀t)ψ(n₀t) = φ(m₀)ψ(n₀)·(φψ)(t), so with the product character φψ
# (primitive of conductor F; here 1·1 = 1, 1·χ₋₃ = χ₋₃·1 = χ₋₃, χ₋₃² = 1):
#     a₀ = φ(m₀)ψ(n₀)/(c_norm β^k) · Σ_{t≠0} (φψ)(t)/t^k
#        = 2 φ(m₀)ψ(n₀) L(k, φψ)/(c_norm β^k)                  [parity again]
#  (L)  = −φ(m₀)ψ(n₀) · B_{k,φψ}/(2k β^k) · C_{φψ,k}/C_{ψ,k},
#     C_{φψ,k}/C_{ψ,k} = τ(φψ)·M^k / (τ(ψ)·F^k)   ∈ K   ((−2πi)^k cancels).
# All transcendental factors cancel; the result is an exact element of K.
#
# ───────────────────────────────────────────────────────────────────────────
# 3) Level-1 dilations g(K·τ) (used for the independent E₄ cross-check)
#
# If g is SL₂(Z)-modular of weight k with constant term a₀(g), then with the
# same Hermite data as in 1): g(K·γτ) = ((cτ+d)/D)^k · g((Aτ+B)/D), so
#     (g(K·)|_k γ)(τ) = D^{−k} g((Aτ+B)/D)  ⇒  cusp value = a₀(g)/D^k,
# D = K/gcd(K·a, c).  (Lattice route: β = −K/g = −D, a₀ = −B_k/(2k β^k) ✓.)
#
# ───────────────────────────────────────────────────────────────────────────
# 4) η-quotient cusp orders (vanishing of 6.4.a.a at every cusp)
#
# ord of η(mτ) at the cusp (a,c) of Γ₁(6), measured in the LOCAL variable
# q_h = e^{2πiτ'/width} (Ligozat):  gcd(c,m)²·width/(24m).

# Loadable both inside the LandauAlphabet module (the package includes this
# file) and standalone after `using LandauAlphabet`; in the in-module case
# the names are already present.
if !isdefined(@__MODULE__, :B2K)
    using Nemo
    using LandauAlphabet
end

# ---------- the coefficient field K = Q(√−3) ----------

# The number field must be constructed lazily AT RUNTIME: FLINT/antic field
# objects do not survive Julia precompilation serialization (segfaults in
# nf_elem arithmetic when baked into the module image as top-level consts).
const _CUSPVALS_FIELD = Ref{Any}(nothing)

"the exact field K = Q(√−3) used for cusp values, and s3 = √−3"
function cuspvals_field()
    if _CUSPVALS_FIELD[] === nothing
        x = polynomial_ring(QQ, "x")[2]
        _CUSPVALS_FIELD[] = number_field(x^2 + 3, "s3")
    end
    _CUSPVALS_FIELD[]
end

# ---------- cusp bookkeeping ----------

"Γ₁(6) cusp representatives (a, c), aligned with GAMMA16_WIDTHS"
const GAMMA16_CUSP_TUPLES = [(1, 0), (0, 1), (1, 2), (1, 3)]
"Γ₁(6) cusp widths, aligned with GAMMA16_CUSP_TUPLES"
const GAMMA16_WIDTHS = [1, 6, 3, 2]

"""
    cusp_completion(cusp) -> (b, d)

Completion of the primitive vector (a, c) to γ = [a b; c d] ∈ SL₂(Z) with
γ(∞) = a/c.  Cusp values are independent of this choice (see header).
"""
function cusp_completion(cusp::Tuple{Int,Int})
    a, c = cusp
    g, d, b = gcdx(a, -c)          # a·d + (−c)·b = 1  ⇒  ad − bc = 1
    g == 1 || error("cusp $(cusp) is not primitive (gcd(a,c) ≠ 1)")
    @assert a * d - b * c == 1
    (b, d)
end

# ---------- 1) weight 2: B_{2,K} ----------

"""
    cusp_value_B2K(Kdil, cusp) -> QQFieldElem

Exact constant term of (B_{2,Kdil}|₂γ) at the cusp (a,c), derivation in the
file header:  (A − D)/(24D) with A = gcd(Kdil·a, c) (= gcd(Kdil, c) for a
primitive cusp; = Kdil at ∞ where c = 0) and D = Kdil/A.
"""
function cusp_value_B2K(Kdil::Int, cusp::Tuple{Int,Int})::QQFieldElem
    Kdil >= 1 || error("Kdil must be a positive dilation")
    a, c = cusp
    gcd(a, c) == 1 || error("cusp $(cusp) is not primitive")
    A = gcd(Kdil * a, c)
    A > 0 || error("invalid cusp (0,0)")
    Kdil % A == 0 || error("internal: A = gcd(Kdil·a,c) must divide Kdil")
    D = div(Kdil, A)
    QQ(A - D, 24 * D)
end

# ---------- 2) generalized Eisenstein series, lattice-sum method ----------

"""
    cusp_value_E(k, Dphi, Dpsi, Kdil, cusp) -> elem of K = Q(√−3)

Exact constant term of (E_k(Kdil·τ; φ, ψ)|_k γ) at the cusp (a,c), for the
Stein-normalized series of `eisenstein_E` with φ = Kronecker(Dphi/·),
ψ = Kronecker(Dpsi/·), Dphi, Dpsi ∈ {1, −3}.  Derivation in the file header:

    a₀ = −φ(m₀) ψ(n₀) · B_{k,φψ}/(2k β^k) · τ(φψ)·M^k/(τ(ψ)·F^k),

with M = cond(ψ), F = cond(φψ), g = gcd(M·Kdil·a, c), (m₀,n₀) = (c, −MKa)/g,
β = −M·Kdil/g.  Rigorous for k ≥ 3 (absolute convergence); k ≤ 2 would need
Hecke regularization and is rejected.
"""
function cusp_value_E(k::Int, Dphi::Int, Dpsi::Int, Kdil::Int, cusp::Tuple{Int,Int})
    k >= 3 || error("lattice-sum cusp values are rigorous only for k ≥ 3 " *
                    "(k = $k needs Hecke regularization); use cusp_value_B2K for weight 2")
    Kdil >= 1 || error("Kdil must be a positive dilation")
    a, c = cusp
    gcd(a, c) == 1 || error("cusp $(cusp) is not primitive")
    b, d = cusp_completion(cusp)
    phi = kronecker_char(Dphi)
    psi = kronecker_char(Dpsi)
    # parity condition φψ(−1) = (−1)^k — outside it the Eisenstein series
    # vanishes identically and silent zeros would be misleading
    parity = (Dphi == -3 ? -1 : 1) * (Dpsi == -3 ? -1 : 1)
    parity == (-1)^k || error("parity violation: φψ(−1) = $parity ≠ (−1)^$k")
    M = LandauAlphabet.char_conductor(Dpsi)
    # lattice line (M·Kdil·a)m + c·n = 0 through the dilated, γ-pulled lattice
    P = M * Kdil * a
    g = gcd(P, c)
    m0 = div(c, g)
    n0 = -div(P, g)
    beta = m0 * M * Kdil * b + n0 * d
    @assert beta * g == -M * Kdil * (a * d - b * c)   # β = −MK/g exactly
    # product character φψ. For Dphi = Dpsi = −3 the pointwise product is the
    # IMPRIMITIVE principal character mod 3 (χ₋₃², zero on 3Z), whose line sum
    # Σ_{t≠0} χ₋₃²(t)/t^k = 2(1 − 3^{−k})ζ(k) carries an Euler factor relative
    # to the primitive (trivial) lift.
    Dprod = (Dphi == Dpsi) ? 1 : -3
    F = LandauAlphabet.char_conductor(Dprod)
    Bk = LandauAlphabet.gen_bernoulli(k, Dprod)
    rat = -QQ(phi(m0) * psi(n0)) * Bk // (QQ(2 * k) * QQ(beta)^k)
    if Dphi == -3 && Dpsi == -3
        rat *= 1 - QQ(1, 3)^k
    end
    KK, S3 = cuspvals_field()
    tau_prod = Dprod == 1 ? one(KK) : S3   # Gauss sum τ(φψ)
    tau_psi = Dpsi == 1 ? one(KK) : S3     # Gauss sum τ(ψ)
    ratio = tau_prod * KK(QQ(M)^k) * inv(tau_psi * KK(QQ(F)^k))
    KK(rat) * ratio
end

# ---------- 3) Hermite-decomposition cross-check for level-1 dilations ----------

"""
    cusp_value_dilated_level1(a0, k, Kdil, cusp) -> QQFieldElem

Constant term at the cusp (a,c) of g(Kdil·τ)|_k γ for g SL₂(Z)-modular of
weight k with constant term a0.  Derivation in the file header (item 3):
value = a0 / D^k with D = Kdil/gcd(Kdil·a, c).
"""
function cusp_value_dilated_level1(a0::QQFieldElem, k::Int, Kdil::Int,
                                   cusp::Tuple{Int,Int})::QQFieldElem
    a, c = cusp
    gcd(a, c) == 1 || error("cusp $(cusp) is not primitive")
    A = gcd(Kdil * a, c)
    Kdil % A == 0 || error("internal: A must divide Kdil")
    D = div(Kdil, A)
    a0 // QQ(D)^k
end

"E₄(Kdil·τ) cusp value via the Hermite method (independent of the lattice sum)"
cusp_value_E4_hermite(Kdil::Int, cusp::Tuple{Int,Int})::QQFieldElem =
    cusp_value_dilated_level1(QQ(1, 240), 4, Kdil, cusp)   # a₀(E₄ Stein) = −B₄/8 = 1/240

# ---------- 4) η-quotient cusp orders ----------

"""
    eta_order_at_cusp(m, cusp, width) -> QQFieldElem

Vanishing order of η(mτ) at the cusp (a,c) of Γ₁(6), in units of the local
variable q_h = e^{2πiτ'/width} (Ligozat): gcd(c,m)²·width/(24m).
"""
eta_order_at_cusp(m::Int, cusp::Tuple{Int,Int}, width::Int)::QQFieldElem =
    QQ(gcd(cusp[2], m)^2 * width, 24 * m)

"""
    cusp4aa_order_at_cusp(cusp, width) -> QQFieldElem

Vanishing order of the cusp form 6.4.a.a = η(τ)²η(2τ)²η(3τ)²η(6τ)² at a
Γ₁(6) cusp in local-variable units; equals 1 at every cusp (so it vanishes
at all of them and contributes an exactly zero cusp-value column).
"""
cusp4aa_order_at_cusp(cusp::Tuple{Int,Int}, width::Int)::QQFieldElem =
    sum(QQ(2) * eta_order_at_cusp(m, cusp, width) for m in (1, 2, 3, 6))

# ---------- assembled cusp-value matrices ----------

"""
    cusp_values_level6(k) -> 4×dim(M_k) Matrix over K = Q(√−3)

Rows = cusps in the order [(1,0)=∞, (0,1), (1,2), (1,3)] (widths [1,6,3,2]);
columns in `level6_basis` order:
  k = 2: [B_{2,2}, B_{2,3}, B_{2,6}]
  k = 3: [E₃(τ;χ₋₃,1), E₃(2τ;χ₋₃,1), E₃(τ;1,χ₋₃), E₃(2τ;1,χ₋₃)]
  k = 4: [E₄(τ), E₄(2τ), E₄(3τ), E₄(6τ), 6.4.a.a]   (cusp-form column ≡ 0,
         justified by cusp4aa_order_at_cusp > 0 at every cusp).
"""
function cusp_values_level6(k::Int)
    KK = cuspvals_field()[1]
    cusps = GAMMA16_CUSP_TUPLES
    if k == 2
        Ks = (2, 3, 6)
        return [KK(cusp_value_B2K(Ks[j], cusps[i])) for i in 1:4, j in 1:3]
    elseif k == 3
        specs = ((-3, 1, 1), (-3, 1, 2), (1, -3, 1), (1, -3, 2))   # (Dphi, Dpsi, Kdil)
        return [cusp_value_E(3, specs[j]..., cusps[i]) for i in 1:4, j in 1:4]
    elseif k == 4
        out = Matrix{elem_type(KK)}(undef, 4, 5)
        for i in 1:4
            for (j, Kd) in enumerate((1, 2, 3, 6))
                out[i, j] = cusp_value_E(4, 1, 1, Kd, cusps[i])
            end
            ord = cusp4aa_order_at_cusp(cusps[i], GAMMA16_WIDTHS[i])
            ord > 0 || error("cusp form 6.4.a.a fails to vanish at cusp $(cusps[i])")
            out[i, 5] = zero(KK)
        end
        return out
    end
    error("weight $k not implemented (have k ∈ {2,3,4})")
end

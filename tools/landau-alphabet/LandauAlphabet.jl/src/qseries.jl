# Exact q-series over Q (Nemo QQFieldElem vectors), and the modular-form
# expansions needed for the Γ₁(6)/Γ₁(12) sunrise–kite alphabet of 1704.08895:
#   * Dedekind-eta quotients (integer leading power tracked separately),
#   * generalized Eisenstein series E_k(τ; φ, ψ) in Stein's normalization,
#   * B_{2,K}(τ) = E₂(τ) − K E₂(Kτ).
#
# All series are coefficient vectors c[1:N] for q⁰…q^{N−1} in a fixed nome.

export QS, qs_mul, qs_add, qs_scal, qs_inv, qs_div, qs_qdq, qs_sqrt_unit,
       eta_quotient_q2, eisenstein_E, B2K, E2series, kronecker_char,
       qs_shift_power, qs_substitute_negq

const QS = Vector{QQFieldElem}

qzero(N) = QQFieldElem[QQ(0) for _ in 1:N]

function qs_mul(a::QS, b::QS)
    N = min(length(a), length(b))
    c = qzero(N)
    for i in 1:N
        iszero(a[i]) && continue
        for j in 1:(N - i + 1)
            iszero(b[j]) && continue
            c[i+j-1] += a[i] * b[j]
        end
    end
    c
end

qs_add(a::QS, b::QS) = [a[i] + b[i] for i in 1:min(length(a), length(b))]
qs_scal(s, a::QS) = [QQ(s) * x for x in a]

function qs_inv(a::QS)
    iszero(a[1]) && error("series not invertible")
    N = length(a)
    b = qzero(N)
    b[1] = inv(a[1])
    for n in 2:N
        s = QQ(0)
        for k in 2:n
            s += a[k] * b[n-k+1]
        end
        b[n] = -s / a[1]
    end
    b
end

qs_div(a::QS, b::QS) = qs_mul(a, qs_inv(b))

"q d/dq"
qs_qdq(a::QS) = [QQ(n - 1) * a[n] for n in 1:length(a)]

"sqrt of a series with a[1] a square of a rational (Newton); returns s with s²=a"
function qs_sqrt_unit(a::QS)
    r0 = sqrt(a[1])   # exact in QQ or throws
    N = length(a)
    b = qzero(N)
    b[1] = r0
    for n in 2:N
        s = QQ(0)
        for k in 2:n-1
            s += b[k] * b[n-k+1]
        end
        b[n] = (a[n] - s) / (2 * r0)
    end
    b
end

"multiply by q^m (shift); drops overflow"
function qs_shift_power(a::QS, m::Int)
    N = length(a)
    c = qzero(N)
    for i in 1:N-m
        c[i+m] = a[i]
    end
    c
end

"a(q) → a(−q)"
qs_substitute_negq(a::QS) = [isodd(n - 1) ? -a[n] : a[n] for n in 1:length(a)]

# ---------- eta quotients in the q₂ = exp(iπτ) variable ----------

"∏_{n≥1} (1−q^{k n}) as a series in q (Euler series via pentagonal numbers)"
function euler_product(k::Int, N::Int)
    e = qzero(N)
    e[1] = QQ(1)
    # generalized pentagonal numbers: m(3m−1)/2, m ∈ Z\{0}
    out = copy(e)
    # use recursive multiplication by sparse pentagonal series of (1−q^k)(1−q^{2k})…
    # via Euler's pentagonal number theorem: ∏(1−x^n) = Σ (−1)^m x^{m(3m−1)/2}
    p = qzero(N)
    p[1] = QQ(1)
    m = 1
    while true
        g1 = div(m * (3m - 1), 2) * k
        g2 = div(m * (3m + 1), 2) * k
        g1 >= N && g2 >= N && break
        if g1 < N
            p[g1+1] += QQ((-1)^m)
        end
        if g2 < N
            p[g2+1] += QQ((-1)^m)
        end
        m += 1
    end
    p
end

"""
    eta_quotient_q2(factors, N) -> (leadpow::QQFieldElem, series::QS)

η-quotient ∏ η(k·τ/2)^{r_k} expanded in q₂ = exp(iπτ):
η(kτ/2) = q₂^{k/24} ∏ (1−q₂^{k n}).  `factors` = vector of (k, r_k).
Returns the total leading power Σ k·r_k/24 (must be handled by caller) and the
integer-power series ∏(Euler_k)^{r_k}.
"""
function eta_quotient_q2(factors::Vector{Tuple{Int,Int}}, N::Int)
    lead = QQ(0)
    s = qzero(N)
    s[1] = QQ(1)
    for (k, r) in factors
        lead += QQ(k * r, 24)
        e = euler_product(k, N)
        if r > 0
            for _ in 1:r
                s = qs_mul(s, e)
            end
        elseif r < 0
            ei = qs_inv(e)
            for _ in 1:(-r)
                s = qs_mul(s, ei)
            end
        end
    end
    (lead, s)
end

# ---------- Dirichlet characters and Eisenstein series ----------

"""
    kronecker_char(D) -> function n -> Int

Real character n ↦ Kronecker symbol (D/n). Implemented for the cases needed:
D = 1 (trivial) and D = −3.
"""
function kronecker_char(D::Int)
    if D == 1
        return n -> 1
    elseif D == -3
        return n -> (r = mod(n, 3); r == 0 ? 0 : (r == 1 ? 1 : -1))
    else
        error("character D=$D not implemented")
    end
end

"conductor of the implemented characters"
char_conductor(D::Int) = D == 1 ? 1 : (D == -3 ? 3 : error("unknown"))

"generalized Bernoulli number B_{k,ψ} for ψ = Kronecker(D/·)"
function gen_bernoulli(k::Int, D::Int)
    M = char_conductor(D)
    chi = kronecker_char(D)
    # Σ_{m=1}^{M} ψ(m) x e^{mx}/(e^{Mx}−1) = Σ B_{k,ψ} x^k/k!  — extract via
    # exact power series in QQ.
    Npow = k + 2
    # e^{mx} series, (e^{Mx}−1)/x series
    function expser(c)  # e^{cx} truncated
        v = [QQ(c)^(j) / factorial(big(j)) for j in 0:Npow]
        QQFieldElem[QQ(x) for x in v]
    end
    den = [QQ(M)^(j + 1) / factorial(big(j + 1)) for j in 0:Npow]  # (e^{Mx}−1)/x
    num = qzero(Npow + 1)
    for m in 1:M
        c = chi(m)
        c == 0 && continue
        num = qs_add(num, qs_scal(c, expser(m)))
    end
    ser = qs_div(num, QQFieldElem[QQ(d) for d in den])   # = Σ B_{k,ψ} x^k / k!
    ser[k+1] * factorial(big(k))
end

"""
    eisenstein_E(k, Dphi, Dpsi, K, N) -> QS

E_k(K·τ; φ, ψ) in the nome q = e^{2πiτ}, Stein normalization:
E_k(τ;φ,ψ) = a₀ + Σ_{m≥1} (Σ_{d|m} ψ(d) φ(m/d) d^{k−1}) q^m,
a₀ = −B_{k,ψ}/(2k) if cond(φ) = 1 else 0.   (Integer q-powers; K rescales m→Km.)
"""
function eisenstein_E(k::Int, Dphi::Int, Dpsi::Int, K::Int, N::Int)
    phi = kronecker_char(Dphi)
    psi = kronecker_char(Dpsi)
    s = qzero(N)
    if char_conductor(Dphi) == 1
        s[1] = -gen_bernoulli(k, Dpsi) / (2k)
    end
    for m in 1:fld(N - 1, K)
        c = QQ(0)
        for d in 1:m
            m % d == 0 || continue
            c += QQ(psi(d) * phi(div(m, d))) * QQ(d)^(k - 1)
        end
        s[K*m+1] += c
    end
    s
end

"E₂(τ) = −1/24 + Σ σ₁(m) qᵐ"
function E2series(K::Int, N::Int)
    s = qzero(N)
    s[1] = QQ(-1, 24)
    for m in 1:fld(N - 1, K)
        c = QQ(0)
        for d in 1:m
            m % d == 0 && (c += d)
        end
        s[K*m+1] += c
    end
    s
end

"B_{2,K}(τ) = E₂(τ) − K E₂(Kτ)"
B2K(K::Int, N::Int) = qs_add(E2series(1, N), qs_scal(-K, E2series(K, N)))

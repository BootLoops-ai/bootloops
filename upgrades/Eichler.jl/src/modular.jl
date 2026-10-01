# Certified building blocks: Dedekind eta (pointwise + q-expansions of eta quotients
# with rigorous majorant tails) and generalised Eisenstein series with exact coefficients.

# ---- pointwise eta -----------------------------------------------------------

"""Certified η(τ) for Im τ > 0."""
function eta_point(tau::Acb)
    res = Acb(prec = Arblib.precision(tau))
    Arblib.modular_eta!(res, tau)
    return res
end

"""Certified eta quotient ∏_d η(d·τ)^{r[d]} ; rs = vector of (d, r_d)."""
function eta_quotient_point(tau::Acb, rs::Vector{Tuple{Int,Int}})
    prec = Arblib.precision(tau)
    val = Acb(1; prec)
    for (d, r) in rs
        e = eta_point(Acb(d * tau; prec))
        val *= e^r
    end
    return val
end

# ---- Euler factor ∏(1-x^n) at certified real points (for majorants) ----------

"""Certified ∏_{n≥1}(1 - xⁿ) for real 0 < x < 1, via η: ∏(1-xⁿ) = η(τ_x)·x^{-1/24},
x = e^{2πiτ_x}, τ_x = i·ln(1/x)/(2π)."""
function euler_factor(x::Arb)
    prec = Arblib.precision(x)
    @assert Arblib.is_positive(x) && x < 1
    lt = log(inv(x))
    tau = Acb(0, lt / (2 * Arb(π; prec)); prec)
    e = eta_point(tau)
    return real(e) * exp(lt / 24)
end

# ---- eta-quotient q-series with exact integer coefficients --------------------

"""Exact expansion of ∏_d ∏_{n≥1} (1 - q^{dn})^{r_d} to order N as fmpz polynomial."""
function _euler_product_exact(rs::Vector{Tuple{Int,Int}}, N::Int)
    R, q = polynomial_ring(Nemo.ZZ, "q")
    # numerator and denominator separately, then series inversion
    num = one(R); den = one(R)
    for (d, r) in rs
        r == 0 && continue
        f = one(R)
        for n in 1:div(N, d)
            f = mullow(f, (1 - q^(d * n))^abs(r), N + 1)
        end
        if r > 0
            num = mullow(num, f, N + 1)
        else
            den = mullow(den, f, N + 1)
        end
    end
    # invert den mod q^{N+1}
    S, qq = power_series_ring(Nemo.QQ, N + 1, "q")
    nums = S(map(QQFieldElem, collect(coefficients_padded(num, N))), N + 1, N + 1, 0)
    dens = S(map(QQFieldElem, collect(coefficients_padded(den, N))), N + 1, N + 1, 0)
    ser = divexact(nums, dens)
    return [Nemo.coeff(ser, n) for n in 0:N]
end

coefficients_padded(p, N) = [Nemo.coeff(p, n) for n in 0:N]

"""QSeries of the eta quotient ∏_d η(dτ)^{r_d} = q^w ∏∏(1-q^{dn})^{r_d}, w = Σ d·r_d/24
(must be a nonnegative integer). Rigorous tail via the positive-coefficient majorant
q^w ∏∏(1-q^{dn})^{-|r_d|}, whose value at x0 is certified through `euler_factor`."""
function eta_quotient_qseries(ctx::SeriesContext, rs::Vector{Tuple{Int,Int}})
    w24 = sum(d * r for (d, r) in rs)
    @assert w24 % 24 == 0 && w24 >= 0
    w = div(w24, 24)
    c = _euler_product_exact(rs, ctx.N - w)
    # Majorant mass of the UNSHIFTED Euler product ∏∏(1-q^{dn})^{r_d}:
    # F = ∏_d E(x0^d)^{-|r_d|}  (no x0^w factor here — qshift applies it once).
    prec = ctx.prec
    x0a = Arb(prec = prec); Arblib.set!(x0a, ctx.x0)
    Fm = one(x0a)
    for (d, r) in rs
        Fm *= euler_factor(x0a^d)^(-abs(r))
    end
    tail = Mag(); Arblib.get!(tail, Arblib.ubound(Fm))
    f = qseries_exact(ctx, c, tail)
    return qshift(f, w)
end

# ---- Dirichlet characters and Eisenstein series -------------------------------

"""χ₋₃(n) = Kronecker symbol (-3/n): period 3, values 1, -1, 0 for n ≡ 1, 2, 0 mod 3."""
chi_m3(n::Integer) = (r = mod(n, 3); r == 1 ? 1 : (r == 2 ? -1 : 0))

"""Rigorous bound for Σ_{n>N} n^k x0ⁿ (k ≥ 0), via n^k x0ⁿ ≤ (N+1)^k x0^{N+1} θ^{n-N-1}
with θ = x0·((N+2)/(N+1))^k, valid when θ < 1 (each step multiplies by ≤ θ)."""
function polygeom_tail(prec::Int, k::Int, N::Int, x0::Mag)
    @assert k >= 0   # the per-step ratio bound needs (1+1/n)^k non-increasing
    x0a = Arb(prec = prec); Arblib.set!(x0a, x0)
    θ = x0a * (Arb(N + 2; prec) / (N + 1))^k
    @assert θ < 1
    bound = (Arb(N + 1; prec))^k * x0a^(N + 1) / (1 - θ)
    t = Mag(); Arblib.get!(t, Arblib.ubound(bound))
    return t
end

"""E1(τ; χ̄0, χ̄1) = 1/6 + Σ_m (Σ_{d|m} χ₋₃(d)) q^m  — exact coefficients.
|c_m| ≤ σ0(m) ≤ m+1 ≤ 2m for m≥1, so tail ≤ 2·polygeom_tail(1)."""
function eisenstein_e1_qseries(ctx::SeriesContext)
    c = zeros(Int, ctx.N + 1)
    for d in 1:ctx.N
        χ = chi_m3(d)
        χ == 0 && continue
        for m in d:d:ctx.N
            c[m+1] += χ
        end
    end
    cq = [n == 0 ? QQFieldElem(1, 6) : QQFieldElem(c[n+1]) for n in 0:ctx.N]
    t = polygeom_tail(ctx.prec, 1, ctx.N, ctx.x0)
    Arblib.mul!(t, t, Mag(2))
    qseries_exact(ctx, cq, t)
end

"""E2(τ) = -1/24 + Σ σ1(m) q^m; σ1(m) ≤ m·σ0(m) ≤ 2m², tail ≤ 2·polygeom_tail(2)."""
function eisenstein_e2_qseries(ctx::SeriesContext)
    c = zeros(BigInt, ctx.N + 1)
    for d in 1:ctx.N, m in d:d:ctx.N
        c[m+1] += d
    end
    cq = [n == 0 ? QQFieldElem(-1, 24) : QQFieldElem(c[n+1]) for n in 0:ctx.N]
    t = polygeom_tail(ctx.prec, 2, ctx.N, ctx.x0)
    Arblib.mul!(t, t, Mag(2))
    qseries_exact(ctx, cq, t)
end

"""B_{2,K}(τ) = E2(τ) - K·E2(Kτ)."""
function b2K_qseries(ctx::SeriesContext, K::Int)
    e2 = eisenstein_e2_qseries(ctx)
    e2K = qsubst_pow(e2, K)
    return e2 - K * Acb(1; prec = ctx.prec) * e2K
end

"""E3(τ; χ̄1, χ̄0) = Σ_m (Σ_{d|m} χ₋₃(m/d) d²) q^m (a0 = 0); |c_m| ≤ m²σ0(m) ≤ 2m³."""
function eisenstein_e3_qseries(ctx::SeriesContext)
    c = zeros(BigInt, ctx.N + 1)
    for d in 1:ctx.N, m in d:d:ctx.N
        χ = chi_m3(div(m, d))
        χ == 0 && continue
        c[m+1] += χ * d^2
    end
    cq = [n == 0 ? QQFieldElem(0) : QQFieldElem(c[n+1]) for n in 0:ctx.N]
    t = polygeom_tail(ctx.prec, 3, ctx.N, ctx.x0)
    Arblib.mul!(t, t, Mag(2))
    qseries_exact(ctx, cq, t)
end

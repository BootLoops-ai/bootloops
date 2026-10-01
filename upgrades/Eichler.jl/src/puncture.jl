# Puncture layer: Abel image of the marked point and the third-kind period G
# (docs/conventions.md §7–8). Units m_t² = 1 unless stated.

"""Incomplete elliptic integral F(φ|m) (Mathematica m-convention) via Carlson RF:
F(φ|m) = sinφ · RF(cos²φ, 1 - m sin²φ, 1)."""
function ellF(φ::Acb, m::Acb)
    prec = Arblib.precision(φ)
    sφ, cφ = sin(φ), cos(φ)
    rf = Acb(prec = prec)
    Arblib.elliptic_rf!(rf, cφ^2, 1 - m * sφ^2, Acb(1; prec); flags = 0, prec = prec)
    return sφ * rf
end

"""Abel image Z_{x_p}(s,t) of the marked point x_p = -s-t on the ACKM parent curve
(conventions §7, printed closed form):
Z = F(arcsin(u)|k²)/(2K(k²)), u = √(1 + (t-8)√s/(t√(s+16)))/√2,
k² = 2√s√(s+16)/(s + √s√(s+16) + 8).
Branches: principal, with the Feynman t → t+i0 prescription resolving the arcsin cut
when u is real with u > 1: asin(u + i0⁺) = π/2 + i·acosh(u)."""
function abel_image(s::Acb, t::Acb)
    prec = Arblib.precision(s)
    ss = sqrt(s)
    s16 = sqrt(s + 16)
    k2 = 2 * ss * s16 / (s + ss * s16 + 8)
    u = sqrt(1 + (t - 8) * ss / (t * s16)) / sqrt(Acb(2; prec))
    if Arblib.contains_zero(imag(u)) && !Arblib.is_negative(abs(real(u)) - 1)
        # u on (or beyond) the asin cut: the +i0 limit is NOT representable as an
        # honest ball here. The caller must opt in to the documented
        # deformation explicitly.
        error("abel_image: marked point lies on the arcsin branch cut at this (s,t); " *
              "use abel_image_deformed(s, t; delta_exp) — its result encloses " *
              "Z(s, t + i·2^-delta_exp), with the +i0 distance documented, not enclosed")
    end
    φ = asin(u)
    K = Acb(prec = prec); Arblib.elliptic_k!(K, k2)
    return ellF(φ, k2) / (2 * K)
end

"""Abel image at real Euclidean t through the explicit Feynman deformation
t → t + iδ with EXACT dyadic δ = 2^(-delta_exp) (default prec÷2). The returned ball
honestly encloses Z(s, t+iδ); the distance to the +i0 limit is O(δ·|∂_t Z|) and is
documented, not enclosed — choose delta_exp ≫ target digits·3.33 accordingly."""
function abel_image_deformed(s::Acb, t::Acb; delta_exp::Int = 0)
    prec = Arblib.precision(s)
    de = delta_exp == 0 ? div(prec, 2) : delta_exp
    δ = Arb(2; prec = prec)^(-de)
    td = t + Acb(0, 1; prec) * Acb(δ; prec)
    ss = sqrt(s)
    s16 = sqrt(s + 16)
    k2 = 2 * ss * s16 / (s + ss * s16 + 8)
    u = sqrt(1 + (td - 8) * ss / (td * s16)) / sqrt(Acb(2; prec))
    (Arblib.contains_zero(imag(u)) && !Arblib.is_negative(abs(real(u)) - 1)) &&
        error("abel_image_deformed: u ball still touches the asin cut; increase precision")
    φ = asin(u)
    K = Acb(prec = prec); Arblib.elliptic_k!(K, k2)
    return ellF(φ, k2) / (2 * K)
end

"""Holomorphic Frobenius period ϖ0(s, x) of the BCNTW curve with mass² = x
(conventions §8): normalized so ϖ0 = (1/√(sx))(1 - s/(64x) + ...) for small s/x.
Implemented through the vanishing-cycle complete elliptic integral:
the curve P4 has roots X = x, x+s, x + (s ± √(s² + 16xs))/2; the cycle that vanishes
as s→0 surrounds {x, x+s}. For s, x > 0 the four roots are real ordered
x < x+s < ... (the quadratic roots are x + (s±√(s²+16xs))/2, one below x? — the minus
root: (s - √(s²+16xs))/2 < 0, so it lies below x... and the plus root above x+s).
Ordering for s,x>0: r₋ = x + (s-√)/2 < x < x+s < r₊. Byrd–Friedman with
(a,b,c,d) = (r₊, x+s, x, r₋): ∮_{(x,x+s)} dX/√(∏) = 4K(k²)/√((a-c)(b-d)).
The quartic leading coefficient of P4 is +1 in X… with P4 = (x-X)(x+s-X)(X²-...):
expanding, leading term +X⁴? (−X)(−X)(X²) = +X⁴ ✓. On the cycle the product is
positive-oriented; the normalization and phase are fixed by the certified match to the
printed series (test suite)."""
function varpi0(s::Acb, x::Acb)
    prec = Arblib.precision(s)
    disc = sqrt(s * s + 16 * x * s)
    rm = x + (s - disc) / 2
    rp = x + (s + disc) / 2
    a, b, c, d = rp, x + s, x, rm
    k2 = ((b - c) * (a - d)) / ((a - c) * (b - d))
    K = Acb(prec = prec); Arblib.elliptic_k!(K, k2)
    pref = 4 / sqrt((a - c) * (b - d))
    # ∮ = pref·K; ϖ0 = ∮/(2π)·(phase): vanishing cycle ∮ → 2π/√(4xs)·(phase) as s→0;
    # series says ϖ0 → 1/√(xs): ϖ0 = ∮/π · (1/2)·... fixed: ϖ0 = pref·K/π
    return pref * K / Acb(π; prec = prec)
end

"""∂ϖ0(s,x)/∂x (derivative in the mass² argument), certified, via the dK/dm chain on
the Byrd–Friedman representation used in `varpi0`."""
function varpi0_prime(s::Acb, x::Acb)
    prec = Arblib.precision(s)
    disc = sqrt(s * s + 16 * x * s)
    ddisc = 8 * s / disc
    rm = x + (s - disc) / 2;  drm = 1 - ddisc / 2
    rp = x + (s + disc) / 2;  drp = 1 + ddisc / 2
    a, da = rp, drp
    b, db = x + s, Acb(1; prec)
    c, dc = x, Acb(1; prec)
    d, dd = rm, drm
    ac, bd = a - c, b - d
    dac, dbd = da - dc, db - dd
    k2 = ((b - c) * (a - d)) / (ac * bd)
    dk2 = ((db - dc) * (a - d) + (b - c) * (da - dd)) / (ac * bd) -
          k2 * (dac / ac + dbd / bd)
    K, E_, dK, _ = _dKE(k2)
    pref = 4 / sqrt(ac * bd)
    dpref = -2 * (dac * bd + ac * dbd) / (ac * bd)^(Acb(3//2; prec))
    return (dpref * K + pref * dK * dk2) / Acb(π; prec)
end

"""Third-kind period G(s,t) at m²=1 (Euclidean region s>0, -s < t < 0, with
t(s+t) < 4s certified). Definition: the Zenodo-2502.00118 differential relation
  ∂_t G = -s²(16+s)·ω0·r9/(t(s+t)(t(s+t)-4s)²) + s(16+s)·r9·ω0′/(2t(s+t)(t(s+t)-4s)),
  r9 = √(t(s+t)(t(s+t)-4s))  [≡ √(P4(m²-t)), positive on the stated region],
anchored at G(s,0) = 0, integrated along t' = -v² (the endpoint 1/√|t| singularity is
removed by the substitution). Certified quadrature; branch hypotheses enforced via the
analytic flag. Cross-checked against the printed double series of arXiv:2502.00118."""
function third_kind_G(s::Acb, t::Acb; rtol = 0.0)
    prec = Arblib.precision(s)
    # entry hypotheses: certified s > 0, -s < t < 0, t(s+t) < 4s; outside this
    # region the r9 branch identification has not been validated
    (Arblib.is_positive(real(s)) && Arblib.is_negative(real(t)) &&
     Arblib.is_positive(real(s + t)) && Arblib.is_positive(real(4 * s - t * (s + t)))) ||
        error("third_kind_G: hypotheses s>0, -s<t<0, t(s+t)<4s not certified")
    vmax = sqrt(-t)
    w0 = varpi0(s, Acb(1; prec))
    w0p = varpi0_prime(s, Acb(1; prec))
    # Desingularized integrand: with t' = -v², u = t'(s+t') = -v²(s+t'),
    # R := r9/v = √(-(s+t')·den) (positive real on the stated region),
    # (∂_t G)·dt = [ s(16+s)·R/((s+t')·den) · ( -2s·w0/den + w0p ) ] dv — analytic in v
    # on the whole closed path, including the endpoint v = 0.
    integrand = let s = s, w0 = w0, w0p = w0p
        (v; analytic::Bool = false) -> begin
            tt = -v * v
            spt = s + tt
            den = tt * spt - 4 * s
            arg = -spt * den
            if analytic
                bad = Arblib.contains_zero(imag(arg)) && !Arblib.is_positive(real(arg))
                bad |= Arblib.contains_zero(den) || Arblib.contains_zero(spt)
                if bad
                    ind = Acb(prec = prec); Arblib.indeterminate!(ind)
                    return ind
                end
            end
            R = sqrt(arg)
            s * (16 + s) * R / (spt * den) * (-2 * s * w0 / den + w0p)
        end
    end
    Arblib.integrate(integrand, Acb(0; prec = prec), vmax;
                     check_analytic = true, prec = prec,
                     rtol = rtol == 0.0 ? exp10(-0.301 * prec) : rtol)
end

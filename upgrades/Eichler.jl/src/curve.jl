# Curve layer: certified periods/quasi-periods of y² = quartic.
#
# Primary rigorous definition: a period is a certified contour integral of x^k dx/y over
# a circle enclosing exactly two branch points, with the square-root branch fixed by the
# midpoint form u_{ab}(x) = (x-m)·sqrt(1 - (δ/(x-m))²), m = (a+b)/2, δ = (b-a)/2, which
# is analytic in x off the segment [a,b] (the argument of sqrt crosses ℝ≤0 exactly when
# δ/(x-m) is real with modulus ≥ 1). Analyticity on integration boxes is enforced through
# the acb_calc `analytic` flag, so the enclosure is rigorous end to end.
#
# Fast path: complete-elliptic (K/E ≡ AGM) formulas for verified real-ordered
# configurations (Byrd–Friedman 256.00/256.12), cross-validated against the quadrature
# in the gauntlet.

struct QuarticCurve
    prec::Int
    lead::Acb              # leading coefficient c of y² = c·∏(x - r_i)
    roots::Vector{Acb}     # certified root enclosures (length 4)
end

"""Quartic curve from exact rational coefficients [c₀, c₁, c₂, c₃, c₄] of
y² = Σ cᵢ xⁱ, with certified roots via Arb's polynomial root isolation."""
function QuarticCurve(coeffs::Vector{QQFieldElem}, prec::Int)
    @assert length(coeffs) == 5
    p = AcbPoly(prec = prec)
    for (i, c) in enumerate(coeffs)
        num = Arb(BigInt(numerator(c)); prec = prec)
        den = Arb(BigInt(denominator(c)); prec = prec)
        Arblib.set_coeff!(p, i - 1, Acb(num / den; prec = prec))
    end
    roots = AcbVector(4; prec = prec)
    n = ccall((:acb_poly_find_roots, Arblib.libflint), Clong,
        (Ptr{Arblib.acb_struct}, Ref{Arblib.acb_poly_struct}, Ptr{Cvoid}, Clong, Clong),
        roots, p, C_NULL, 0, prec)
    n == 4 || error("QuarticCurve: could not certify 4 isolated roots")
    rts = [Acb(roots[i]; prec = prec) for i in 1:4]
    # Rigorous realness tightening: coefficients are real, so conj(root) is a root.
    # If conj(ball_i) intersects no other isolated enclosure, the root in ball_i is
    # its own conjugate, i.e. exactly real.
    for i in 1:4
        Arblib.contains_zero(imag(rts[i])) || continue
        ci = conj(rts[i])
        if !any(j != i && Arblib.overlaps(ci, rts[j]) for j in 1:4)
            rts[i] = Acb(real(rts[i]), Arb(0; prec = prec); prec = prec)
        end
    end
    QuarticCurve(prec, Acb(Arblib.ref(p, 4); prec = prec), rts)
end

"""u_{ab}(x): branch of sqrt((x-a)(x-b)) analytic off segment [a,b]. `analytic=true`
returns an indeterminate ball if analyticity on the box cannot be certified."""
function _pair_sqrt(x::Acb, a::Acb, b::Acb, analytic::Bool)
    prec = Arblib.precision(x)
    m = (a + b) / 2
    δ = (b - a) / 2
    xm = x - m
    w = δ / xm
    v = 1 - w * w
    if analytic
        # certified: v must avoid ℝ≤0; also xm must avoid 0
        bad = Arblib.contains_zero(imag(v)) && !Arblib.is_positive(real(v))
        bad |= Arblib.contains_zero(xm)
        if bad
            ind = Acb(prec = prec)
            Arblib.indeterminate!(ind)
            return ind
        end
    end
    return xm * sqrt(v)
end

"""Branch of sqrt((x-e)(x-f)) analytic OFF the horizontal rays {e + t, t≥0} and
{f - t, t≥0}: v(x) = i·√(e-x)·√(x-f) (principal square roots). Used for the root pair
NOT enclosed by the contour (cut pushed through infinity). `e` should be the root with
the larger real part."""
function _pair_sqrt_outside(x::Acb, e::Acb, f::Acb, analytic::Bool)
    prec = Arblib.precision(x)
    w1 = e - x
    w2 = x - f
    if analytic
        bad = (Arblib.contains_zero(imag(w1)) && !Arblib.is_positive(real(w1))) ||
              (Arblib.contains_zero(imag(w2)) && !Arblib.is_positive(real(w2)))
        if bad
            ind = Acb(prec = prec)
            Arblib.indeterminate!(ind)
            return ind
        end
    end
    return Acb(0, 1; prec = prec) * sqrt(w1) * sqrt(w2)
end

"""Distance (certified lower bound, Arb) from point z to the horizontal ray
{w0 + t·dir, t ≥ 0} with dir = ±1."""
function _dist_ray(z::Acb, w0::Acb, dir::Int)
    dz = (z - w0) * dir
    re, im = real(dz), imag(dz)
    if Arblib.is_negative(re)
        return Arblib.lbound(Arb, abs(dz))
    end
    return Arblib.lbound(Arb, abs(im))
end

"""Certified ∮ x^k dx / y over the circle |x - m_{ab}| = R enclosing branch points a,b
and excluding the other two roots. Branch: y = sqrt(lead)·u_{ab}(x)·v_{ef}(x), with the
(a,b) cut on the segment [a,b] inside the contour and the (e,f) cut pushed to infinity
along horizontal rays. R chosen automatically; analyticity on every quadrature box is
enforced through the acb_calc analytic flag (rigorous regardless of the R heuristic).
Orientation: counterclockwise in x."""
function period_quadrature(C::QuarticCurve, pair::Tuple{Int,Int}, k::Int = 0;
                           rtol = 0.0)
    prec = C.prec
    (ia, ib) = pair
    rest = setdiff(1:4, [ia, ib])
    a, b = C.roots[ia], C.roots[ib]
    r1, r2 = C.roots[rest[1]], C.roots[rest[2]]
    # e = larger real part (cut rightward), f = smaller (cut leftward)
    e, f = Float64(real(r1)) >= Float64(real(r2)) ? (r1, r2) : (r2, r1)
    m = (a + b) / 2
    δm = abs((b - a) / 2)
    Rmax = min(_dist_ray(m, e, +1), _dist_ray(m, f, -1))
    Rmin_f = Float64(Arblib.ubound(Arb, δm))
    Rmax_f = Float64(Rmax)
    Rmax_f > Rmin_f * 1.02 || error("period_quadrature: no admissible contour radius")
    R = Arb(sqrt(max(Rmin_f, 1e-300) * Rmax_f); prec = prec)
    sc = sqrt(C.lead)
    integrand = let m = m, R = R, a = a, b = b, e = e, f = f, sc = sc, k = k
        (θ; analytic::Bool = false) -> begin
            ph = exp(Acb(0, 1; prec = prec) * θ)
            x = m + R * ph
            uab = _pair_sqrt(x, a, b, analytic)
            vef = _pair_sqrt_outside(x, e, f, analytic)
            dx = Acb(0, 1; prec = prec) * R * ph
            x^k * dx / (sc * uab * vef)
        end
    end
    val = Arblib.integrate(integrand, Acb(0; prec = prec), Acb(2 * Arb(π; prec); prec = prec);
                           check_analytic = true, prec = prec,
                           rtol = rtol == 0.0 ? exp10(-0.301 * prec) : rtol)
    return val
end

"""Complete-elliptic fast path for real-ordered roots a > b > c > d (certified order
check): cycle around (b,c):  ∮ = 2·∫_c^b dx/√((a-x)(b-x)(x-c)(x-d)) = 4 K(k²)/√((a-c)(b-d)),
k² = (b-c)(a-d)/((a-c)(b-d)) (Byrd–Friedman 256.00), all for lead coefficient +1 quartic
opened up as -(x-a)(x-b)(x-c)(x-d) > 0 on (c,b)... — returned value is for
y² = lead·∏(x-rᵢ) with the same branch only up to a documented phase; gauntlet (d)
fixes and verifies the phase against `period_quadrature`."""
function curve_periods(C::QuarticCurve)
    prec = C.prec
    r = [real(x) for x in C.roots]
    all(Arblib.is_real, C.roots) || error("curve_periods fast path: roots not all real; use period_quadrature")
    idx = sortperm([Float64(x) for x in r], rev = true)  # a > b > c > d heuristically
    a, b, c, d = (r[i] for i in idx)
    (a > b && b > c && c > d) || error("curve_periods: could not certify root ordering")
    k2 = ((b - c) * (a - d)) / ((a - c) * (b - d))
    k2p = ((a - b) * (c - d)) / ((a - c) * (b - d))
    pref = 4 / sqrt((a - c) * (b - d))
    K1 = Acb(prec = prec); Arblib.elliptic_k!(K1, Acb(k2; prec = prec))
    K2 = Acb(prec = prec); Arblib.elliptic_k!(K2, Acb(k2p; prec = prec))
    ψ1 = Acb(pref; prec = prec) * K1
    ψ2 = Acb(0, 1; prec = prec) * Acb(pref; prec = prec) * K2
    # principal sqrt of the leading coefficient (audit: |lead| would drop the phase
    # for negative/complex leads; gauntlet (d) ties this branch to period_quadrature)
    sl = sqrt(C.lead)
    return ψ1 / sl, ψ2 / sl
end

# ---- ACKM parent-curve periods and Picard–Fuchs gauntlet (conventions §7) ------

"""ψ0(s), ψ1(s) of arXiv:2402.07311 (printed eqs.), with K(m), E(m) m-convention.
Valid for s > 0 (balls)."""
function acm_periods(s::Arb)
    prec = Arblib.precision(s)
    sA = Acb(s; prec = prec)
    E1 = Acb(prec = prec); Arblib.elliptic_e!(E1, Acb(-s / 16; prec = prec))
    ψ0 = 32 * E1 / (Acb(π; prec = prec) * sA^(Acb(3//2; prec = prec)) * (sA + 16))
    m2 = Acb(s / 16 + 1; prec = prec)
    E2 = Acb(prec = prec); Arblib.elliptic_e!(E2, m2)
    K2 = Acb(prec = prec); Arblib.elliptic_k!(K2, m2)
    ψ1 = 32 * (E2 - K2) / (sA^(Acb(3//2; prec = prec)) * (sA + 16))
    return ψ0, ψ1
end

"""dK/dm = (E - (1-m)K)/(2m(1-m)), dE/dm = (E - K)/(2m)."""
function _dKE(m::Acb)
    prec = Arblib.precision(m)
    K = Acb(prec = prec); Arblib.elliptic_k!(K, m)
    E = Acb(prec = prec); Arblib.elliptic_e!(E, m)
    dK = (E - (1 - m) * K) / (2 * m * (1 - m))
    dE = (E - K) / (2 * m)
    return K, E, dK, dE
end

"""Certified residual of the printed PF operator L0 = ∂² + (4/s + 2/(16+s))∂ +
6(6+s)/(s²(16+s)) applied to ψ0(s) (gauntlet (c)); returns a ball that must contain 0."""
function picard_fuchs_residual_acm(s::Arb)
    prec = Arblib.precision(s)
    sA = Acb(s; prec = prec)
    m = -sA / 16
    K, E, dK, dE = _dKE(m)
    dm = Acb(-1//16; prec = prec)
    π_ = Acb(π; prec = prec)
    den = π_ * sA^(Acb(3//2; prec)) * (sA + 16)
    f = 32 * E / den
    # f' = 32 [E' dm / den - E d(den)/ds / den²]
    dden = π_ * ((Acb(3//2; prec)) * sA^(Acb(1//2; prec)) * (sA + 16) + sA^(Acb(3//2; prec)))
    fp = 32 * (dE * dm / den - E * dden / den^2)
    # f'' : differentiate again
    d2E = _d2E(m) * dm^2
    d2den = π_ * ((Acb(3//4; prec)) * sA^(Acb(-1//2; prec)) * (sA + 16) +
                  2 * (Acb(3//2; prec)) * sA^(Acb(1//2; prec)))
    fpp = 32 * (d2E / den - 2 * dE * dm * dden / den^2 -
                E * d2den / den^2 + 2 * E * dden^2 / den^3)
    L0 = fpp + (4 / sA + 2 / (16 + sA)) * fp + 6 * (6 + sA) / (sA^2 * (16 + sA)) * f
    return L0
end

"""d²E/dm² = -(E - (1+m... ) — derived from dE/dm = (E-K)/(2m):
d²E/dm² = (dE - dK)/(2m) - (E-K)/(2m²)."""
function _d2E(m::Acb)
    K, E, dK, dE = _dKE(m)
    return (dE - dK) / (2 * m) - (E - K) / (2 * m^2)
end

"""ODEOperator for L0 (cleared denominators): s²(16+s)y'' + (4s(16+s)+2s²)y' + 6(6+s)y."""
function acm_pf_operator(prec::Int)
    P(coeffs...) = [epoly_const(prec, 1, c) for c in coeffs]
    p2 = P(0, 0, 16, 1)              # s²(16+s) = 16s² + s³
    p1 = P(0, 64, 6)                 # 4s(16+s) + 2s² = 64s + 6s²
    p0 = P(36, 6)                    # 6(6+s)
    sing = [Acb(0; prec), Acb(-16; prec)]
    ODEOperator(2, 1, prec, [p0, p1, p2], sing)
end

"""Gauntlet (c) transport check: take certified (ψ0, ψ0') at s0 from the closed form,
transport with the printed L0 to s1, compare against the closed form at s1.
Returns (transported, direct) — overlapping balls certify that L0 annihilates ψ0
along the path."""
function frobenius_transport_acm(s0::Arb, s1::Arb; M::Int = 64)
    prec = Arblib.precision(s0)
    op = acm_pf_operator(prec)
    ψ0, dψ0 = _acm_psi0_with_derivative(s0)
    Y = [[ψ0], [dψ0]]
    Yt, _ = transport(op, [Acb(s0; prec = prec), Acb(s1; prec = prec)], Y; M = M)
    ψ0b, _ = _acm_psi0_with_derivative(s1)
    return Yt[1][1], ψ0b
end

function _acm_psi0_with_derivative(s::Arb)
    prec = Arblib.precision(s)
    sA = Acb(s; prec = prec)
    m = -sA / 16
    K, E, dK, dE = _dKE(m)
    π_ = Acb(π; prec = prec)
    den = π_ * sA^(Acb(3//2; prec)) * (sA + 16)
    dden = π_ * ((Acb(3//2; prec)) * sA^(Acb(1//2; prec)) * (sA + 16) + sA^(Acb(3//2; prec)))
    ψ0 = 32 * E / den
    dψ0 = 32 * (dE * Acb(-1//16; prec) / den - E * dden / den^2)
    return ψ0, dψ0
end

"""Certified incomplete integral ∫_{r_b}^{x_p} dx/y from the branch point r_b = roots[ib]
to x_p along the straight path, with endpoint desingularization x = b + u²·(x_p-b):
∫ = ∫₀¹ 2Δ du / (sc·σ·√(Δ(x(u)-a))·v_ef(x(u))), Δ = x_p - b. The sign σ ∈ {±1} ties the
desingularized square root to the SAME y-branch used by `period_quadrature` (certified by
a midpoint comparison against u_ab). For x_p on a branch cut, pass x_p with an explicit
+iδ deformation (exact dyadic), as in `abel_image_deformed`."""
function incomplete_quadrature(C::QuarticCurve, pair::Tuple{Int,Int}, xp::Acb; rtol = 0.0)
    prec = C.prec
    (ia, ib) = pair
    rest = setdiff(1:4, [ia, ib])
    a, b = C.roots[ia], C.roots[ib]
    r1, r2 = C.roots[rest[1]], C.roots[rest[2]]
    e, f = Float64(real(r1)) >= Float64(real(r2)) ? (r1, r2) : (r2, r1)
    Δ = xp - b
    sc = sqrt(C.lead)
    # σ: match w(u) = σ·√(Δ(x-a)) against u_ab(x)/u at u = 1/2
    uhalf = Acb(1//2; prec)
    xh = b + uhalf^2 * Δ
    w_ref = _pair_sqrt(xh, a, b, false) / uhalf
    w_can = sqrt(Δ * (xh - a))
    ratio = w_ref / w_can
    σ = Acb(0; prec)
    if Arblib.contains(ratio, Acb(1; prec)) && !Arblib.contains(ratio, Acb(-1; prec))
        σ = Acb(1; prec)
    elseif Arblib.contains(ratio, Acb(-1; prec)) && !Arblib.contains(ratio, Acb(1; prec))
        σ = Acb(-1; prec)
    else
        error("incomplete_quadrature: could not certify branch sign σ")
    end
    integrand = let b = b, Δ = Δ, a = a, e = e, f = f, sc = sc, σ = σ
        (u; analytic::Bool = false) -> begin
            x = b + u * u * Δ
            arg = Δ * (x - a)
            if analytic
                bad = Arblib.contains_zero(imag(arg)) && !Arblib.is_positive(real(arg))
                if bad
                    ind = Acb(prec = prec); Arblib.indeterminate!(ind)
                    return ind
                end
            end
            w = σ * sqrt(arg)
            vef = _pair_sqrt_outside(x, e, f, analytic)
            2 * Δ / (sc * w * vef)
        end
    end
    Arblib.integrate(integrand, Acb(0; prec = prec), Acb(1; prec = prec);
                     check_analytic = true, prec = prec,
                     rtol = rtol == 0.0 ? exp10(-0.301 * prec) : rtol)
end

"""Abel image of a marked point on an arbitrary quartic curve: Z = ∫_{r_b}^{x_p}(dx/y) / ∮(dx/y),
both with the same branch conventions (cycle = `pair`). Defined modulo the period lattice;
the representative is the one produced by the straight desingularized path."""
function abel_image_quartic(C::QuarticCurve, pair::Tuple{Int,Int}, xp::Acb)
    I = incomplete_quadrature(C, pair, xp)
    ψ = period_quadrature(C, pair)
    return I / ψ
end

"""Incomplete integral ∫_{r_b}^{x_p} dx/y for x_p REAL and strictly inside the open
segment (a,b) (both roots real), taken on the +i0 edge of the [a,b] cut. On the edge,
u_ab continues to u_edge(x) = σ_e·i·u·√(-Δ(x-a)) after the desingularization
x = b + u²Δ (using (x-m)²(w²-1) = -(x-a)(x-b)); -Δ(x-a) > 0 on the whole path, so the
integrand is analytic on [0,1]. The edge sign σ_e is certified by comparing the two
analytic-in-η functions u_ab(x₀+iη) and σ·i·(x₀+iη-m)·√(w(x₀+iη)²-1) at one point
η > 0 (they agree identically on the band once separated at one point, and both are
continuous at η = 0⁺)."""
function incomplete_edge_quadrature(C::QuarticCurve, pair::Tuple{Int,Int}, xp::Acb; rtol = 0.0)
    prec = C.prec
    (ia, ib) = pair
    rest = setdiff(1:4, [ia, ib])
    a, b = C.roots[ia], C.roots[ib]
    # Contract (audit): the +i0-edge orientation is defined for a < b (path starts at
    # the larger root b and runs left to xp ∈ (a,b)); the reversed order returns the
    # negative of the documented integral — reject it.
    (real(a) < real(b)) || error("incomplete_edge_quadrature: pass pair = (smaller, larger) root indices")
    r1, r2 = C.roots[rest[1]], C.roots[rest[2]]
    e, f = Float64(real(r1)) >= Float64(real(r2)) ? (r1, r2) : (r2, r1)
    m = (a + b) / 2
    δab = (b - a) / 2
    Δ = xp - b
    sc = sqrt(C.lead)
    # σ_e certification at x₀ = midpoint of (b, xp), η = 2^-10
    η = Acb(0, 1; prec) * Acb(Arb(2; prec = prec)^(-10); prec)
    x0 = (b + xp) / 2 + η
    w0 = δab / (x0 - m)
    cand = Acb(0, 1; prec) * (x0 - m) * sqrt(w0 * w0 - 1)
    refv = _pair_sqrt(x0, a, b, false)
    ratio = refv / cand
    local σe
    if Arblib.contains(ratio, Acb(1; prec)) && !Arblib.contains(ratio, Acb(-1; prec))
        σe = Acb(1; prec)
    elseif Arblib.contains(ratio, Acb(-1; prec)) && !Arblib.contains(ratio, Acb(1; prec))
        σe = Acb(-1; prec)
    else
        error("incomplete_edge_quadrature: could not certify edge sign")
    end
    integrand = let b = b, Δ = Δ, a = a, e = e, f = f, sc = sc, σe = σe
        (u; analytic::Bool = false) -> begin
            x = b + u * u * Δ
            arg = -Δ * (x - a)
            if analytic
                bad = Arblib.contains_zero(imag(arg)) && !Arblib.is_positive(real(arg))
                if bad
                    ind = Acb(prec = prec); Arblib.indeterminate!(ind)
                    return ind
                end
            end
            w = σe * Acb(0, 1; prec) * sqrt(arg)
            vef = _pair_sqrt_outside(x, e, f, analytic)
            2 * Δ / (sc * w * vef)
        end
    end
    Arblib.integrate(integrand, Acb(0; prec = prec), Acb(1; prec = prec);
                     check_analytic = true, prec = prec,
                     rtol = rtol == 0.0 ? exp10(-0.301 * prec) : rtol)
end

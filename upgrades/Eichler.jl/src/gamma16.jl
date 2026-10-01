# The Γ₁(6) modular layer (conventions: docs/conventions.md §3, from arXiv:1704.08895
# sect. "The elliptic curve related to the maximal cut"). Everything is a function of
# τ = τ_C with nome q = q_C = e^{2πiτ_C}. Units: μ² = m² = 1 internally.

"""All Γ₁(6) data as certified q_C-series on a common SeriesContext."""
struct Gamma16
    ctx::SeriesContext
    e1::QSeries        # E1(τ; χ̄0, χ̄1)
    e2::QSeries        # E1(2τ; χ̄0, χ̄1)
    t::QSeries         # Hauptmodul t(q_C) = 9 η(6τ)⁸η(τ)⁴/(η(2τ)⁸η(3τ)⁴)
    psi1::QSeries      # ψ₁/π = 2√3 (e1+e2)
    f1::QSeries        # 3√2 e1
    f2::QSeries        # -6(e1² + 6e1e2 - 4e2²)
    f3::QSeries        # 36√3 (e1³ - e1²e2 - 4e1e2² + 4e2³)
    f4::QSeries        # f1⁴
    g20::QSeries       # -12(e1² - 4e2²)
    g21::QSeries       # -18(e1² + e1e2 - 2e2²)
    g29::QSeries       # 6(e1² - 3e1e2 + 2e2²)
    g31::QSeries       # -108√3 (e1³ - 3e1e2² + 2e2³)
    sqrt2::Acb
    sqrt3::Acb
end

function Gamma16(prec::Int, N::Int)
    ctx = SeriesContext(prec, N)
    e1 = eisenstein_e1_qseries(ctx)
    e2 = qsubst_pow(e1, 2)
    s2 = sqrt(Acb(2; prec))
    s3 = sqrt(Acb(3; prec))
    t = eta_quotient_qseries(ctx, [(6, 8), (1, 4), (2, -8), (3, -4)])
    t = 9 * t
    e1e1 = e1 * e1
    e1e2 = e1 * e2
    e2e2 = e2 * e2
    psi1 = (2 * s3) * (e1 + e2)
    f1 = (3 * s2) * e1
    f2 = Acb(-6; prec) * (e1e1 + 6 * e1e2 - 4 * e2e2)
    f3 = (36 * s3) * (e1e1 * e1 - e1e1 * e2 - 4 * (e1e2 * e2) + 4 * (e2e2 * e2))
    f4 = f1^4
    g20 = Acb(-12; prec) * (e1e1 - 4 * e2e2)
    g21 = Acb(-18; prec) * (e1e1 + e1e2 - 2 * e2e2)
    g29 = Acb(6; prec) * (e1e1 - 3 * e1e2 + 2 * e2e2)
    g31 = (-108 * s3) * (e1e1 * e1 - 3 * (e1e2 * e2) + 2 * (e2e2 * e2))
    Gamma16(ctx, e1, e2, t, psi1, f1, f2, f3, f4, g20, g21, g29, g31, s2, s3)
end

kernel(G::Gamma16, name::Symbol) = getfield(G, name)

# ---- Hauptmodul, pointwise --------------------------------------------------

"""Certified t(τ) = 9 η(6τ)⁸η(τ)⁴/(η(2τ)⁸η(3τ)⁴) (t in units of m²)."""
hauptmodul_point(tau::Acb) =
    9 * eta_quotient_point(tau, [(6, 8), (1, 4), (2, -8), (3, -4)])

# ---- τ_C(t) via complete elliptic integrals (AW maximal-cut conventions) ------

"""K with AW modulus convention: K_AW(k) = acb_elliptic_k(k²)."""
function ellK(m::Acb)
    res = Acb(prec = Arblib.precision(m))
    Arblib.elliptic_k!(res, m)
    return res
end
function ellE(m::Acb)
    res = Acb(prec = Arblib.precision(m))
    Arblib.elliptic_e!(res, m)
    return res
end

"""(ψ_{1,C}/π, τ_C) from the AW maximal-cut K-formulas with principal branches.
WARNING: the principal-branch K-ratio lands in the correct Γ₁(6) orbit only in a
neighbourhood of the cusp t = 0 (audit: it silently leaves it for t ≲ −0.47). The
result is therefore SELF-CHECKED here against the certified Hauptmodul and errors out
on mismatch. For general kinematic points use `nome_and_log_from_t` (certified Newton
on the cusp-connected branch), which is what `sunrise2` uses."""
function tau_from_t(t::Acb; check::Bool = true)
    prec = Arblib.precision(t)
    st = sqrt_upper(t)
    z1 = t - 4
    z2 = -1 - 2 * st
    z3 = -1 + 2 * st
    z4 = t
    kC2 = ((z3 - z2) * (z4 - z1)) / ((z3 - z1) * (z4 - z2))
    kC2p = ((z2 - z1) * (z4 - z3)) / ((z3 - z1) * (z4 - z2))
    pref = (1 + st)^(Acb(3; prec) / 2) * sqrt(Acb(3; prec) - st)
    psi1C = 4 * ellK(kC2) / (pref * Acb(π; prec))
    psi2C = 4 * Acb(0, 1; prec) * ellK(kC2p) / (pref * Acb(π; prec))
    tauC = psi2C / psi1C
    Arblib.is_positive(imag(tauC)) || error("tau_from_t: Im τ ≤ 0 — outside validated region")
    if check
        Arblib.overlaps(hauptmodul_point(tauC), t) ||
            error("tau_from_t: principal-branch τ is not in the cusp-connected orbit at this t; use nome_and_log_from_t")
    end
    return psi1C, tauC
end

"""√t with the t → t+i0 prescription: for Re t < 0, Im t = 0 take the upper branch."""
function sqrt_upper(t::Acb)
    prec = Arblib.precision(t)
    if Arblib.is_real(t) && !Arblib.is_positive(real(t))
        return Acb(0, 1; prec) * sqrt(-t)
    end
    return sqrt(t)
end

"""Nome q_C = e^{2πiτ}."""
nome(tauC::Acb) = exp(2 * Acb(0, π; prec = Arblib.precision(tauC)) * tauC)

# ---- certified Newton inverse of the Hauptmodul --------------------------------

"""Midpoint of an Acb ball as an exact Acb."""
function _midpoint(q::Acb, prec::Int)
    re = Arb(prec = prec); Arblib.set!(re, Arblib.midref(Arblib.realref(q)))
    im = Arb(prec = prec); Arblib.set!(im, Arblib.midref(Arblib.imagref(q)))
    Acb(re, im; prec = prec)
end

"""Float seed for the nome at kinematic point t: Newton on the (midpoint) series
from q₀ = t/9 — the branch continuously connected to the cusp q(0) = 0."""
function _seed_nome(G::Gamma16, t::Acb)
    nf = min(Arblib.degree(G.t.coeffs), 80)
    cf = [ComplexF64(Float64(real(Arblib.ref(G.t.coeffs, n))), Float64(imag(Arblib.ref(G.t.coeffs, n)))) for n in 0:nf]
    tt = ComplexF64(Float64(real(t)), Float64(imag(t)))
    q = tt / 9
    for _ in 1:60
        f = evalpoly(q, cf) - tt
        fp = sum(n * cf[n+1] * q^(n - 1) for n in 1:nf)
        q -= f / fp
    end
    return q
end

"""Certified inverse: given target t (ball) and the Gamma16 series, find a ball q with
t(exact q) = t certified by the univariate complex interval-Newton contraction test
on the certified series enclosure. The returned branch is the one continuously
connected to the cusp q(t=0) = 0 (seeded by `_seed_nome`). Returns q_C ball."""
function nome_from_t(G::Gamma16, t::Acb; seed::Union{Nothing,Acb} = nothing)
    prec = G.ctx.prec
    q = seed
    if q === nothing
        qf = _seed_nome(G, t)
        # keep the nome exactly real for exactly-real t (series has real coefficients)
        q = Arblib.is_real(t) ? Acb(real(qf), 0; prec = prec) :
            Acb(real(qf), imag(qf); prec = prec)
    end
    q = Acb(q; prec = prec)
    ts = G.t
    tder = qdq(ts)   # q·dt/dq on halved reference radius
    # Newton iterations with point center, then a final certified contraction check.
    for it in 1:ceil(Int, log2(prec)) + 4
        qm = _midpoint(q, prec)
        fval = evaluate(ts, qm) - t
        fder = evaluate(tder, qm) / qm     # dt/dq = (q dt/dq)/q
        q = qm - fval / fder
    end
    # certified step: q_ball := midpoint ± inflated radius; check Newton image ⊂ q_ball
    qm = _midpoint(q, prec)
    fval = evaluate(ts, qm) - t
    fder0 = evaluate(tder, qm) / qm
    step = fval / fder0
    r = Mag(); Arblib.get!(r, step)
    Arblib.mul!(r, r, Mag(16))   # inflation
    qball = Acb(qm); Arblib.add_error!(qball, r)
    fderB = evaluate(qdq(ts), qball) / qball
    img = qm - fval / fderB              # N(B) = m - f(m)/f'(B), f'(B) enclosure
    (Arblib.contains(qball, img) && !Arblib.contains_zero(fderB)) ||
        error("nome_from_t: certified Newton contraction failed")
    return img
end

"""Nome q_C, logarithm L = 2πiτ_C and τ_C at kinematic point t, on the branch
continuously connected to the cusp t=0 with the Feynman t+i0 prescription
(docs/conventions.md §3–4). For real t < 0 the nome is on the negative real axis and
L = log(-q) + iπ exactly; otherwise the principal log is used (the ball must avoid
the cut). The result is certified: t(q) ∋ t by interval Newton."""
function nome_and_log_from_t(G::Gamma16, t::Acb)
    prec = G.ctx.prec
    q = nome_from_t(G, t)
    iπ = Acb(0, 1; prec) * Acb(Arb(π; prec); prec)
    local L
    if Arblib.is_real(t) && Arblib.is_negative(real(t))
        # The true nome is exactly real: the series has real coefficients, so conj(q*)
        # is also a root of t(q) = t inside the conjugated (= same, since the Newton
        # ball is centred on the real axis) ball; uniqueness from the interval-Newton
        # contraction forces q* = conj(q*). Tighten to the real part.
        Arblib.contains_zero(imag(q)) ||
            error("nome_and_log_from_t: nome ball does not touch the real axis")
        q = Acb(real(q), Arb(0; prec); prec = prec)
        Arblib.is_negative(real(q)) || error("nome_and_log_from_t: expected negative nome for t<0")
        L = Acb(log(-real(q)); prec) + iπ
    else
        # principal log; require the ball to avoid ℝ≤0
        (Arblib.contains_zero(imag(q)) && !Arblib.is_positive(real(q))) &&
            error("nome_and_log_from_t: nome ball touches the log cut; use higher precision")
        L = log(q)
    end
    tauC = L / (2 * iπ)
    return q, L, tauC
end

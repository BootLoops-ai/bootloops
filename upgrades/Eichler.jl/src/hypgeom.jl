# Certified ε-expansion of ₂F₁(-2ε, -ε; 1-ε; x) at x = r₃ = e^{2πi/3} (and conj),
# via a certified series start at x = 1/8 and Frobenius transport along 1/8 → r₃.
# Used for the sunrise boundary constants (docs/conventions.md §5).

"""₂F₁(-2ε,-ε;1-ε;x) and its x-derivative at a point x with |x| ≤ 1/4, as ε-polys
(length Jeps), from the defining series with a rigorous ℓ₁ tail bound.
Term recurrence: t_0 = 1, t_{n+1} = t_n (n-2ε)(n-ε)/((n+1-ε)(n+1)) x."""
function hyp2f1_eps_series_smallx(prec::Int, Jeps::Int, x::Acb, M::Int = 200)
    f = epoly_zero(prec, Jeps)
    fp = epoly_zero(prec, Jeps)             # derivative d/dx
    t = epoly_const(prec, Jeps, 1)          # t_n · xⁿ accumulated separately
    xa = Acb(x; prec = prec)
    xn = Acb(1; prec)                        # xⁿ
    for n in 0:M
        f = epoly_add(f, epoly_scal(xn, t))
        if n >= 1
            fp = epoly_add(fp, epoly_scal(Acb(n; prec) * xn / xa, t))
        end
        # t_{n+1} = t_n (n-2ε)(n-ε) / ((n+1-ε)(n+1))
        a = epoly_const(prec, Jeps, n); Jeps >= 2 && (a[2] = Acb(-2; prec))
        b = epoly_const(prec, Jeps, n); Jeps >= 2 && (b[2] = Acb(-1; prec))
        cden = epoly_const(prec, Jeps, n + 1); Jeps >= 2 && (cden[2] = Acb(-1; prec))
        t = epoly_mul(t, epoly_mul(a, b))
        t = epoly_mul(t, epoly_inv(cden))
        t = epoly_scal(Acb(1; prec) / (n + 1), t)
        xn *= xa
    end
    # Tail. After the loop t = t_{M+1} (as ε-poly) and xn = x^{M+1}, so the first
    # neglected term has ℓ₁ norm a_{M+1} := ‖t_{M+1}‖·|x|^{M+1}. Successive ratio:
    # ‖t_{n+1}x^{n+1}/(t_n xⁿ)‖ ≤ |x|·(n+2)(n+1)/n² ≤ |x|·(M+3)(M+2)/(M+1)² =: θ
    # using ‖(n-2ε)(n-ε)‖ ≤ (n+2)(n+1) and ‖(n+1-ε)⁻¹‖ ≤ 1/n (geometric series in 1/(n+1)),
    # so ratio ≤ |x|(n+2)(n+1)/(n(n+1)) ≤ |x|(n+2)/n, decreasing; θ at n = M+1.
    absx = abs(xa)                       # Arb ball
    θa = absx * Arb((M + 3) * (M + 2); prec) / Arb((M + 1)^2; prec)
    θa < 1 || error("hyp2f1 small-x series: tail ratio ≥ 1")
    tnorm = Arb(prec = prec); Arblib.set!(tnorm, epoly_norm(t))
    aM1 = tnorm * abs(xn)
    tailf = aM1 / (1 - θa)
    tm = Mag(); Arblib.get!(tm, Arblib.ubound(tailf))
    for k in 1:Jeps
        Arblib.add_error!(f[k], tm)
    end
    # Derivative tail: Σ_{n>M} n‖t_n‖|x|ⁿ⁻¹; term ratio gains (n+1)/n ≤ (M+2)/(M+1):
    θd = θa * Arb(M + 2; prec) / (M + 1)
    θd < 1 || error("hyp2f1 small-x series: derivative tail ratio ≥ 1")
    taild = Arb(M + 1; prec) * aM1 / absx / (1 - θd)
    tdm = Mag(); Arblib.get!(tdm, Arblib.ubound(taild))
    for k in 1:Jeps
        Arblib.add_error!(fp[k], tdm)
    end
    return f, fp
end

"""ODE operator for ₂F₁(a,b;c;x): x(1-x)y'' + (c-(a+b+1)x)y' - ab·y, with
a = -2ε, b = -ε, c = 1-ε as ε-polys."""
function hyp2f1_operator(prec::Int, Jeps::Int)
    a = epoly_zero(prec, Jeps); Jeps >= 2 && (a[2] = Acb(-2; prec))
    b = epoly_zero(prec, Jeps); Jeps >= 2 && (b[2] = Acb(-1; prec))
    c = epoly_const(prec, Jeps, 1); Jeps >= 2 && (c[2] = Acb(-1; prec))
    one_ = epoly_const(prec, Jeps, 1)
    ab1 = epoly_add(epoly_add(a, b), one_)
    p2 = [epoly_zero(prec, Jeps), epoly_const(prec, Jeps, 1), epoly_const(prec, Jeps, -1)]   # x - x²
    p1 = [c, epoly_scal(Acb(-1; prec), ab1)]                                                  # c - (a+b+1)x
    p0 = [epoly_scal(Acb(-1; prec), epoly_mul(a, b))]                                          # -ab
    sing = [Acb(0; prec), Acb(1; prec)]
    ODEOperator(2, Jeps, prec, [p0, p1, p2], sing)
end

"""Certified ε-expansion (length Jeps) of ₂F₁(-2ε,-ε;1-ε; r₃), r₃ = e^{2πi/3}."""
function hyp2f1_at_r3(prec::Int, Jeps::Int; M::Int = 0)
    # start-series tail ~ 8^{-(M+1)} = 2^{-3(M+1)}: take 3M ≥ prec + 80 so the series
    # start does not cap the certified output precision (audit finding)
    M == 0 && (M = max(80, ceil(Int, (prec + 80) / 3)))
    x_start = Acb(1//8; prec)
    f, fp = hyp2f1_eps_series_smallx(prec, Jeps, x_start, M)
    op = hyp2f1_operator(prec, Jeps)
    r3 = exp(2 * Acb(0, π; prec) / 3)
    # path avoiding both singularities: 1/8 → i/2 → (r3 + i·something)... straight legs
    path = [x_start, Acb(0, 1//2; prec) + Acb(1//8; prec), r3]
    Y, _ = transport(op, path, [f, fp]; M = M)
    return Y[1]
end

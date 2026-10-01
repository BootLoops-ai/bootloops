@testset "CY transport: K3 banana (n=2)" begin
    prec = 400
    L3 = banana_pf(3; prec = prec)
    @test L3.order == 3
    frob = mum_frobenius_basis(L3, prec; nterms = 120)
    @test frob.alpha == 0

    # ϖ_0 = Domb series Σ A002895(k) z^k.
    domb = [1, 4, 28, 256, 2716, 31504, 387136, 4951552, 65218204, 878536624]
    for (k, d) in enumerate(domb)
        @test Arblib.overlaps(E.qcoeff(frob.h[1], k - 1), Acb(d; prec))
    end
    # log-tower coefficients (cross-checked against an independent Frobenius-period computation)
    g1 = [0, 6, 57, 584, 13081//2]
    for (k, g) in enumerate(g1)
        @test Arblib.overlaps(E.qcoeff(frob.h[2], k - 1), Acb(g; prec))
    end

    # ϖ_0(0.01) against the certified Python reference (48 digits).
    z = Acb(1//100; prec)
    Π = period_vector(L3, z, prec; frob = frob)
    ref = Arb("1.043086754591857283080419623179049381651482240127"; prec)
    @test Float64(abs(Π[1] - ref)) < 1e-45

    # frobenius_logseries packages ϖ_k as a LogQSeries; evaluate matches _frob_eval.
    V0 = frobenius_logseries(frob, 0)
    V1 = frobenius_logseries(frob, 1)
    @test Arblib.overlaps(E.evaluate(V0, z, log(z)), Π[1])
    @test Arblib.overlaps(E.evaluate(V1, z, log(z)), Π[2])
end

@testset "CY transport: 4-loop CY3 banana (n=3)" begin
    prec = 300
    L4 = banana_pf(4; prec = prec)
    @test L4.order == 4
    frob = mum_frobenius_basis(L4, prec; nterms = 100)
    @test frob.alpha == 1                          # PWW indicial (θ−1)⁴

    # holomorphic series ψ_0 = y(1 − 5y + 45y² − …)  (PWW eq.32 / step01_operator.json)
    holo = [1, -5, 45, -545, 7885, -127905, 2241225, -41467725, 798562125, -15855173825]
    for (k, c) in enumerate(holo)
        @test Arblib.overlaps(E.qcoeff(frob.h[1], k - 1), Acb(c; prec))
    end
    g1 = [0, -8, 100, -4148//3, 64198//3]
    for (k, c) in enumerate(g1)
        @test Arblib.overlaps(E.qcoeff(frob.h[2], k - 1), Acb(c; prec))
    end

    # ψ_0(0.01) and ψ_1/(2πi)(0.01) ≥ 30 digits (independent ref: step23_periods_Y.json)
    y = Acb(1//100; prec)
    Π = period_vector(L4, y, prec; frob = frob)
    twoπi = 2 * Acb(0, 1; prec) * Acb(π; prec)
    ref0 = Arb("0.00954022952720891421234112408863"; prec)
    ref1 = Acb(Arb(0; prec), Arb("0.00710569104267874076467076602053"; prec); prec)
    @test Float64(abs(Π[1] - ref0)) < 1e-30
    @test Float64(abs(Π[2] / twoπi - ref1)) < 1e-29

    # transport: period_matrix at y=1/10 (outside the MUM disk R=1/25) via the
    # θ-companion stepper from two independent in-disk anchors must agree to ≥70d.
    sec = CYSector(L4; anchor = 1//100, nterms = 140)
    @test Float64(abs(sec.W_anchor[1, 1] - ref0)) < 1e-30
    yt = Acb(1//10; prec)
    Wt  = period_matrix(sec, yt; path = [Acb(1//100; prec), yt])
    Wt2 = period_matrix(L4, yt, prec; frob = sec.frob, anchor = Acb(1//200; prec),
                        path = [Acb(1//200; prec), yt])
    @test maximum(Float64(abs(Wt[i, k] - Wt2[i, k])) for i in 1:4, k in 1:4) < 1e-70
    @test Arblib.rel_accuracy_bits(Wt[1, 1]) > 200
    @info "CY3 ψ_0(1/10) via θ-transport" val = Wt[1, 1]
end

@testset "CY transport: n=1 reduces to elliptic / Γ₁(6)" begin
    prec = 400
    L2 = banana_pf(2; prec = prec)
    @test L2.order == 2
    frob = mum_frobenius_basis(L2, prec; nterms = 200)
    @test frob.alpha == 0
    for (k, c) in enumerate([1, 1//3, 5//27, 31//243])
        @test Arblib.overlaps(E.qcoeff(frob.h[1], k - 1), Acb(c; prec))
    end

    G = Gamma16(prec, 200)
    t1, t2 = Acb(1//7; prec), Acb(1//3; prec)
    q1, Lq1, _ = nome_and_log_from_t(G, t1)
    q2, Lq2, _ = nome_and_log_from_t(G, t2)

    # (a) ϖ_0(t) / ψ₁(q_Γ(t)) is t-independent (= √3/2); 100+ digits.
    s3 = sqrt(Acb(3; prec))
    r1 = E._frob_eval(frob, 0, 0, t1) / E.evaluate(G.psi1, q1)
    r2 = E._frob_eval(frob, 0, 0, t2) / E.evaluate(G.psi1, q2)
    @test Arblib.overlaps(r1, r2)
    @test Arblib.overlaps(r1, s3 / 2)
    @test -log10(Float64(Arblib.radref(real(r1))[])) > 100

    # (b) mirror_map: log q_mirror(t) − log q_Γ(t) is the constant log 9
    #     (Hauptmodul leading coefficient); 100+ digits.
    _, lq1 = mirror_map(frob, t1, prec)
    _, lq2 = mirror_map(frob, t2, prec)
    @test Arblib.overlaps(lq1 - Lq1, lq2 - Lq2)
    @test Arblib.overlaps(lq1 - Lq1, log(Acb(9; prec)))
    @test -log10(Float64(Arblib.radref(real(lq1 - Lq1))[])) > 100

    # (c) period_vector at t=1/10 via direct Frobenius ≡ via CERTIFIED transport
    #     from t=1/100 (frobenius.jl engine) ≡ via θ-companion transport.
    Wd = period_matrix(L2, Acb(1//10; prec), prec; frob = frob)
    Wc = period_matrix_certified(L2, Acb(1//10; prec), prec; frob = frob,
            anchor = Acb(1//100; prec), path = [Acb(1//100; prec), Acb(1//10; prec)])
    Wθ = period_matrix(L2, Acb(1//10; prec), prec; frob = frob,
            anchor = Acb(1//100; prec), path = [Acb(1//100; prec), Acb(1//10; prec)])
    @test all(Arblib.overlaps(Wd[i,k], Wc[i,k]) for i in 1:2, k in 1:2)
    @test maximum(Float64(abs(Wθ[i,k] - Wd[i,k])) for i in 1:2, k in 1:2) < 1e-100

    # (d) cy_transport (order-2 VoP) ≡ elliptic-layer integral.
    #     Abel ⇒ det W = 9/(t(t−1)(t−9)); with source S = p₂ one has f = S/p₂ = 1
    #     and I_2 = ∫ ϖ_0/detW dt = (1/9)·∫ ϖ_0 · t(t−1)(t−9) dt = (√3/18)·∫ ψ₁ · p₂ dt,
    #     evaluated independently through the Γ₁(6) layer.  100+ digit agreement.
    sec = CYSector(L2; anchor = 1//7, nterms = 220)
    src = t -> E._pcoeff(L2, 2, t)                 # S(t) = p₂(t) = t(t−1)(t−9)
    g, I = cy_transport(sec, src, t1, t2; method = :gauss_legendre, n = 160)
    secE = EllipticSector(G)
    Rell = (t; analytic::Bool = false) -> begin
        qn, = nome_and_log_from_t(G, t)
        E.evaluate(G.psi1, qn) * t * (t - 1) * (t - 9)
    end
    Iell = (s3 / 18) * vop_quadrature(secE, Rell, t1, t2; method = :gauss_legendre, n = 160)
    @test Float64(abs(I[2] - Iell)) < 1e-100
    # I_1 = −∫ ϖ_1/detW dt; check against direct quadrature of the Frobenius ϖ_1.
    I1ref = Arblib.integrate(
        s -> -E._frob_eval(frob, 1, 0, s) * s * (s - 1) * (s - 9) / 9, t1, t2; prec = prec)
    @test Float64(abs(I[1] - I1ref)) < 1e-100
    # g_p(z₂) = ϖ_0 I_1 + ϖ_1 I_2
    Wend = period_matrix(sec, t2)
    @test Float64(abs(g[1] - (Wend[1,1]*I[1] + Wend[1,2]*I[2]))) < 1e-100
end

@testset "CY transport: K3 inhomogeneous VoP self-consistency" begin
    prec = 256
    L3 = banana_pf(3; prec = prec)
    sec = CYSector(L3; anchor = 1//100, nterms = 140)
    z1, z2, zm = Acb(1//200; prec), Acb(1//40; prec), Acb(1//80; prec)
    src = z -> Acb(1; prec)
    _, I_a = cy_transport(sec, src, z1, z2; method = :gauss_legendre, n = 96)
    _, I_b = cy_transport(sec, src, z1, z2; method = :gauss_legendre, n = 160)
    for k in 1:3
        @test Float64(abs(I_a[k] - I_b[k])) < 1e-40
    end
    # additivity of the VoP integrals
    _, Ia = cy_transport(sec, src, z1, zm; method = :gauss_legendre, n = 96)
    _, Ib = cy_transport(sec, src, zm, z2; method = :gauss_legendre, n = 96)
    for k in 1:3
        @test Float64(abs(Ia[k] + Ib[k] - I_a[k])) < 1e-40
    end
end

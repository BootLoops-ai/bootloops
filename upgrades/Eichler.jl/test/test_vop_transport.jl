@testset "VoP transport (Γ₁(6) sunrise)" begin
    prec = 256
    G = Gamma16(prec, 120)
    sec = EllipticSector(G)

    # 1. tau_map round-trips through the Hauptmodul.
    for tval in (Acb(-1; prec), Acb(-1//3; prec), Acb(1//2; prec))
        τ, q, L = tau_map(sec, tval)
        @test Arblib.overlaps(hauptmodul_point(τ), tval)
        @test Arblib.overlaps(nome(τ), q)
    end

    # 2. eichler_transport between two interior points = direct iterated_integral_value
    #    difference (the LogQSeries is the SAME object).
    f3 = E.kernel(G, :f3)
    t1, t2 = Acb(-1//3; prec), Acb(-1; prec)
    val, pp = eichler_transport(sec, f3, t1, t2)
    _, q1, L1 = tau_map(sec, t1); _, q2, L2 = tau_map(sec, t2)
    direct = iterated_integral_value([f3], q2, L2) - iterated_integral_value([f3], q1, L1)
    @test Arblib.overlaps(val, direct)
    # f3 has zero constant term ⇒ period polynomial = [0, 0].
    @test all(Arblib.contains_zero, pp)

    # 3. Cusp-regularised transport from t=0: f2 has constant term -1/2, so the
    #    primitive's log-level constant is -1/2 (folds into c₂); the value itself
    #    must be finite.
    f2 = E.kernel(G, :f2)
    valc, ppc = eichler_transport(sec, f2, sec.cusp, t2)
    @test Arblib.overlaps(ppc[2], Acb(-1//2; prec))
    @test Arblib.is_finite(valc)

    # 4. vop_fit_constants recovers synthetic (c1, c2) from anchors and the
    #    held-out residual contains 0.
    pref = t -> Acb(1; prec)
    c1_true = [Acb(7//3; prec)]; c2_true = [Acb(-5//4; prec)]
    mkanchor = t -> (t, vop_assemble(sec, pref, [f3], c1_true, c2_true, t))
    anchors = [mkanchor(Acb(-1//5; prec)), mkanchor(Acb(-1; prec)),
               mkanchor(Acb(-2; prec))]
    c1, c2, res = vop_fit_constants(sec, pref, [f3], anchors)
    @test Arblib.overlaps(c1[1], c1_true[1])
    @test Arblib.overlaps(c2[1], c2_true[1])
    @test length(res) == 1 && Arblib.contains_zero(res[1][1])

    # 5. Timing: q-series transport at one point is sub-second (vs ~15 min for a
    #    direct Python quadrature of the unequal-mass kite).
    t0 = time()
    for _ in 1:10
        eichler_transport(sec, f3, sec.cusp, Acb(-1; prec))
    end
    el = (time() - t0) / 10
    @test el < 1.0
    @info "eichler_transport: $(round(el*1000; digits=2)) ms/point @ prec=$prec, N=120"

    # 6. Quadrature fallback cross-check: ∫_{t1}^{t2} 1/(1-t) dt = log((1-t1)/(1-t2)).
    R = (t; analytic::Bool = false) -> begin
        d = 1 - t
        if analytic && Arblib.contains_zero(d)
            ind = Acb(prec = prec); Arblib.indeterminate!(ind); return ind
        end
        inv(d)
    end
    Iq = vop_quadrature(sec, R, t1, t2; method = :certified)
    Ig = vop_quadrature(sec, R, t1, t2; method = :gauss_legendre, n = 48)
    exact = log((1 - t1) / (1 - t2))
    @test Arblib.overlaps(Iq, exact)
    @test Float64(Arblib.ubound(Arb, abs(Ig - exact))) < 1e-40
end

@testset "Multi-curve VoP transport" begin
    prec = 256
    G = Gamma16(prec, 120)
    secA = EllipticSector(G)
    secB = EllipticSector(G)               # same Γ₁(6) curve, second copy
    f3 = E.kernel(G, :f3)
    f2 = E.kernel(G, :f2)
    t1, t2 = Acb(-1//3; prec), Acb(-1; prec)

    # (a) DIAGONAL: zero coupling → two independent single-curve transports.
    msec0 = MultiEllipticSector([secA, secB])
    @test msec0.coupling === nothing
    S = Matrix{Any}(nothing, 2, 2); S[1,1] = f3; S[2,2] = f2
    M0 = multi_eichler_transport(msec0, S, t1, t2)
    vA, ppA = eichler_transport(secA, f3, t1, t2)
    vB, ppB = eichler_transport(secB, f2, t1, t2)
    @test Arblib.overlaps(M0[1,1][1], vA) && all(Arblib.overlaps.(M0[1,1][2], ppA))
    @test Arblib.overlaps(M0[2,2][1], vB) && all(Arblib.overlaps.(M0[2,2][2], ppB))
    @test Arblib.contains_zero(M0[1,2][1]) && Arblib.contains_zero(M0[2,1][1])
    # multi_tau_map returns one (τ,q,L) per curve
    τs = multi_tau_map(msec0, t2)
    @test length(τs) == 2 && Arblib.overlaps(τs[1][1], τs[2][1])

    # (b) Synthetic CONSTANT coupling: curve A = curve B, S_AB = f3, K_AB = c.
    #     Off-diagonal must equal c · ∫ f3 dτ exactly (short-circuit path).
    c = Acb(7//3; prec)
    Kc = t -> [Acb(0;prec) c; Acb(0;prec) Acb(0;prec)]
    msecK = MultiEllipticSector([secA, secB]; coupling = Kc)
    S2 = Matrix{Any}(nothing, 2, 2); S2[1,1] = f3; S2[2,2] = f3; S2[1,2] = f3
    Mk = multi_eichler_transport(msecK, S2, t1, t2)
    @test Arblib.overlaps(Mk[1,2][1], c * vA)
    @test Arblib.overlaps(Mk[1,1][1], vA)

    # (c) t-DEPENDENT coupling via quadrature: K_AB(t) = t  ⇒  off-diagonal
    #     ∫ f3(q(t))·t·(dτ/dt) dt = ∫ f3(τ)·t(τ) dτ, and t(τ) is the Hauptmodul
    #     q-series G.t, so the closed form is eichler_transport(sec, f3·G.t,…).
    #     K(t1)=-1/3 ≠ K(t2)=-1 so the constant-coupling probe is NOT triggered
    #     and the genuine quadrature path runs.
    J = dlogq_dt_gamma16(G)
    Kv = t -> [Acb(0;prec) Acb(t;prec); Acb(0;prec) Acb(0;prec)]
    msecQ = MultiEllipticSector([secA, secB]; coupling = Kv, dlogq_dt = Any[J, J])
    Mq = multi_eichler_transport(msecQ, S2, t1, t2; method = :gauss_legendre, n = 64)
    f3t = f3 * G.t
    vclosed, _ = eichler_transport(secA, f3t, t1, t2)
    @test Float64(Arblib.ubound(Arb, abs(Mq[1,2][1] - vclosed))) < 1e-30
    # diagonal entries unchanged by coupling
    @test Arblib.overlaps(Mq[1,1][1], vA)

    # (d) multi_vop_assemble with zero coupling reduces to two independent
    #     vop_assemble calls.
    pref = [t -> Acb(1; prec), t -> Acb(1; prec)]
    c1 = [Acb(2; prec), Acb(-1; prec)]; c2 = [Acb(0; prec), Acb(3; prec)]
    g = multi_vop_assemble(msec0, pref, S, c1, c2, t2)
    g1 = vop_assemble(secA, pref[1], [f3], [c1[1]], [c2[1]], t2)[1]
    g2 = vop_assemble(secB, pref[2], [f2], [c1[2]], [c2[2]], t2)[1]
    @test Arblib.overlaps(g[1], g1)
    @test Arblib.overlaps(g[2], g2)

    # (e) multi_vop_assemble with NONZERO t-dependent coupling: g[1] must carry
    #     the off-diagonal Eichler term.  Closed form as in (c): with
    #     K_AB(t) = t the cusp-to-t2 contribution is
    #     eichler_transport(secA, f3·G.t, cusp, t2), so
    #     g[1] = vop_assemble(secA,…)[1] + that value.  Guards the one-hot
    #     kernel-matrix orientation in multi_vop_assemble (a transposed one-hot
    #     silently drops the whole coupling term — (d) alone cannot see that).
    gq = multi_vop_assemble(msecQ, pref, S2, c1, c2, t2;
                            method = :gauss_legendre, n = 64)
    base1 = vop_assemble(secA, pref[1], [f3], [c1[1]], [c2[1]], t2)[1]
    voff, _ = eichler_transport(secA, f3t, secA.cusp, t2)
    @test !Arblib.contains_zero(voff)          # coupling term genuinely nonzero
    @test Float64(Arblib.ubound(Arb, abs(gq[1] - (base1 + voff)))) < 1e-40
    @test !Arblib.overlaps(gq[1], base1)       # …and it actually landed in g[1]
    # curve 2 has no off-diagonal source (S2[2,1] === nothing)
    @test Arblib.overlaps(gq[2], vop_assemble(secB, pref[2], [f3], [c1[2]], [c2[2]], t2)[1])
end

using Test, Arblib, Nemo, Eichler
const E = Eichler

@testset "Eichler.jl" begin

@testset "QSeries core" begin
    ctx = SeriesContext(128, 40)
    f = E.qseries_const(ctx, 2)
    g = E.qseries_const(ctx, 3)
    @test Float64(real(E.qcoeff(f * g, 0))) == 6.0
    # eta quotient: Hauptmodul expansion = 9(q - 4q² + 10q³ - 20q⁴ + 39q⁵ - 76q⁶ + 140q⁷)
    G = Gamma16(128, 40)
    expect = [0, 9, -36, 90, -180, 351, -684, 1260]
    for (n, c) in enumerate(expect)
        @test Arblib.overlaps(E.qcoeff(G.t, n - 1), Acb(c; prec = 128))
    end
end

@testset "Gamma16 kernels vs eta quotients (pointwise)" begin
    prec = 256
    G = Gamma16(prec, 60)
    tau = Acb(1//10, 1; prec)
    q = nome(tau)
    s3 = sqrt(Acb(3; prec))
    checks = [
        (G.t,   9 * E.eta_quotient_point(tau, [(6,8),(1,4),(2,-8),(3,-4)])),
        (G.f3,  -3*s3 * E.eta_quotient_point(tau, [(1,5),(3,1),(6,4),(2,-4)])),
        (G.g31, -54*s3 * E.eta_quotient_point(tau, [(6,9),(2,-3)])),
        (G.psi1, 2/s3 * E.eta_quotient_point(tau, [(3,1),(2,6),(1,-3),(6,-2)])),
        (G.g20, E.eta_quotient_point(tau, [(1,4),(3,4),(2,-2),(6,-2)])),
        (G.g21, -9*E.eta_quotient_point(tau, [(3,3),(6,3),(1,-1),(2,-1)])),
        (G.g29, -E.eta_quotient_point(tau, [(1,7),(6,7),(2,-5),(3,-5)])),
    ]
    for (ser, pt) in checks
        @test Arblib.overlaps(E.evaluate(ser, q), pt)
    end
    # cusp values
    @test Arblib.overlaps(E.qcoeff(G.f2, 0), Acb(-1//2; prec))
    @test Arblib.overlaps(E.qcoeff(G.g20, 0), Acb(1; prec))
    @test Arblib.contains_zero(E.qcoeff(G.f3, 0))
end

@testset "Certified nome inversion" begin
    prec = 256
    G = Gamma16(prec, 80)
    for tval in [Acb(-1; prec), Acb(-1//10; prec), Acb(1//2; prec)]
        q, L, tau = nome_and_log_from_t(G, tval)
        @test Arblib.overlaps(hauptmodul_point(tau), tval)
    end
end

@testset "2F1 transport vs Arb polylog" begin
    prec = 256
    F = E.hyp2f1_at_r3(prec, 4)
    r3 = exp(2 * Acb(0, π; prec) / 3)
    li2 = Acb(prec = prec); Arblib.polylog!(li2, Acb(2; prec), r3)
    @test Arblib.overlaps(F[1], Acb(1; prec))
    @test Arblib.contains_zero(F[2])
    @test Arblib.overlaps(F[3], 2 * li2)
end

@testset "Gauntlet (a): t=0 boundary = (3√3/2)L(χ₋₃,2)" begin
    prec = 384
    B = sunrise_boundary_constants(prec, 3)
    L2 = E.dirichlet_L_chim3(2, prec)
    @test Arblib.overlaps(B[1], Acb(3 * sqrt(Arb(3; prec = prec)) / 2 * L2; prec = prec))
    @test E.certified_digits_test(real(B[1])) > 70
end

@testset "Gauntlet (e): AMFlow anchor at t=-1" begin
    prec = 768
    G = Gamma16(prec, 700)
    S4 = sunrise4(G, Acb(-1; prec); J = 10)
    path = get(ENV, "EICHLER_AMFLOW_ANCHOR",
               joinpath(@__DIR__, "data", "anchor_sunrise_dps130_M240.json"))
    if isfile(path)
        txt = read(path, String)
        for k in -2:7
            v = real(-E.eps_coeff(S4, k))
            m = match(Regex("\"$(k)\":\\s*\\{\\s*\"re\":\\s*\"([0-9.e+-]+)\""), txt)
            m === nothing && continue
            d = Float64(Arblib.ubound(Arb, abs(v - Arb(m.captures[1]; prec = prec))))
            @test d < 1e-126
        end
    end
end

@testset "Gauntlet (c): Picard-Fuchs" begin
    prec = 512
    r = picard_fuchs_residual_acm(Arb(3; prec = prec))
    @test Arblib.contains_zero(real(r)) && Arblib.contains_zero(imag(r))
    tr, dir = frobenius_transport_acm(Arb(2; prec = prec), Arb(5; prec = prec))
    @test Arblib.overlaps(tr[1], dir)
end

@testset "Gauntlet (d): periods quad vs K/AGM" begin
    prec = 512
    s = Nemo.QQ(7, 3)
    C = QuarticCurve([Nemo.QQ(0), -4*s^2, s^2 - 4*s, 2*s, Nemo.QQ(1)], prec)
    ord = sortperm([Float64(real(r)) for r in C.roots])
    Pq = period_quadrature(C, (ord[2], ord[3]))
    ψ1, _ = curve_periods(C)
    @test Arblib.overlaps(Pq, ψ1)
end

@testset "Puncture layer" begin
    prec = 256
    # varpi0 vs printed Frobenius series
    s = Acb(1//1000; prec); x = Acb(1; prec)
    w = E.varpi0(s, x)
    ser = 1/sqrt(s*x)*(1 - s/(64*x) + 9*s^2/(16384*x^2) - 25*s^3/(1048576*x^3))
    @test Float64(abs(real(w - ser)) / abs(real(ser))) < 1e-11
    # Abel image: deformed evaluation returns a finite certified ball
    Z = E.abel_image_deformed(Acb(2; prec), Acb(-1//2; prec))
    @test Arblib.is_finite(Z)
end

@testset "Dictionaries" begin
    prec = 256
    L2 = E.dirichlet_L_chim3(2, prec)
    # L(χ₋₃,2) = (ψ'(1/3) - ψ'(2/3))/9 — known value 0.78130...
    @test abs(Float64(L2) - 0.7813024128964864) < 1e-12
    d = cusp_dictionary(prec; maxweight = 4)
    @test length(d) > 10
end

@testset "Kronecker symbol & chi_minus" begin
    @test [E.kronecker_symbol(-3, n) for n in 1:6] == [1, -1, 0, 1, -1, 0]
    @test [E.kronecker_symbol(-4, n) for n in 1:8] == [1, 0, -1, 0, 1, 0, -1, 0]
    @test [E.kronecker_symbol(-7, n) for n in 1:7] == [1, 1, -1, 1, -1, -1, 0]
    @test [E.kronecker_symbol(-8, n) for n in 1:8] == [1, 0, 1, 0, -1, 0, -1, 0]
    @test chi_minus(3) == [1, -1, 0]
    @test chi_minus(4) == [1, 0, -1, 0]
    @test sum(chi_minus(15)) == 0
    @test_throws ErrorException chi_minus(5)   # -5 not a fundamental discriminant
end

@testset "General Dirichlet L: Catalan and cross-checks" begin
    prec = 512
    # L(χ₋₄,2) = Catalan, 60-digit reference value
    cat_ref = Arb("0.915965594177219015054603514932384110774149374281672134266498"; prec)
    cat = dirichlet_L(chi_minus(4), 2, prec)
    @test Float64(Arblib.ubound(Arb, abs(cat - cat_ref))) < 1e-55
    @test E.certified_digits_test(cat) > 130
    # L(χ₋₃,2) must agree with the legacy dirichlet_L_chim3
    @test Arblib.overlaps(dirichlet_L(chi_minus(3), 2, prec), E.dirichlet_L_chim3(2, prec))
    # L(χ₋₄,3) = π³/32 (closed form) and matches the JSON reference
    pi_ = Arb(π; prec)
    L43 = dirichlet_L(chi_minus(4), 3, prec)
    @test Arblib.overlaps(L43, pi_^3 / 32)
    L43_ref = Arb("0.968946146259369380483634845846918600069540267683909615442017"; prec)
    @test Float64(Arblib.ubound(Arb, abs(L43 - L43_ref))) < 1e-55
    # s=1 closed forms: L(χ₋₄,1)=π/4, L(χ₋₃,1)=π/(3√3)
    @test Arblib.overlaps(dirichlet_L(chi_minus(4), 1, prec), pi_ / 4)
    @test Arblib.overlaps(dirichlet_L(chi_minus(3), 1, prec), pi_ / (3 * sqrt(Arb(3; prec))))
end

@testset "Dirichlet L vs polylog/Clausen (independent route)" begin
    prec = 512
    # L(χ,2) = (1/√N) Σ_a χ(a) Im Li₂(e^{2πia/N}) for primitive odd real χ (g(χ)=i√N)
    for N in (7, 8)
        chi = chi_minus(N)
        Lhz = dirichlet_L(chi, 2, prec)
        acc = Arb(0; prec)
        for a in 1:N
            chi[a] == 0 && continue
            z = exp(2 * Acb(0, π; prec) * a / N)
            li = Acb(prec = prec); Arblib.polylog!(li, Acb(2; prec), z)
            acc += chi[a] * imag(li)
        end
        Lpl = acc / sqrt(Arb(N; prec))
        @test Arblib.overlaps(Lhz, Lpl)
        @test E.certified_digits_test(Lhz) > 130
    end
end

@testset "Extended cusp_dictionary keys" begin
    prec = 256
    d = cusp_dictionary(prec; maxweight = 4, conductors = [3, 4, 7, 8])
    keys = first.(d)
    for N in (3, 4, 7, 8), w in 2:4
        @test "L(chi-$N,$w)" in keys
    end
    @test "Catalan" in keys
    @test "log4" in keys && "log7" in keys && "log8" in keys
    @test "Cl2(pi/4)" in keys && "Cl3(pi/7)" in keys
    dd = Dict(d)
    @test Arblib.overlaps(dd["Catalan"], dd["L(chi-4,2)"])
end

@testset "Eisenstein cusp constants & reg tail" begin
    prec = 256
    ec4 = eisenstein_constants(4, prec; maxweight = 3)
    @test Arblib.overlaps(ec4[(2, 1//0)], Arb(-3; prec))
    @test Arblib.overlaps(ec4[(2, 0//1)], Arb(3//4; prec))
    @test Arblib.overlaps(ec4[(2, 1//2)], Arb(0; prec))
    # b2K_qseries normalization check: a₀ = (-1/24)(1-4) = 1/8, i.e. -3/(-24)
    ctx = SeriesContext(prec, 60)
    b24 = E.b2K_qseries(ctx, 4)
    @test Arblib.overlaps(E.qcoeff(b24, 0), Acb(1//8; prec))
    ec6 = eisenstein_constants(6, prec)
    @test Arblib.overlaps(ec6[(2, 1//3)], Arb(-1//2; prec))
    # weight-3: a₀(E₃(χ₋₄,𝟏)) = L(χ₋₄,-2) = -B_{3,χ}/3; B_{3,χ₋₄} = 3/2 ⇒ -1/2
    @test Arblib.overlaps(ec4[(3, 1//0)], Arb(-1//2; prec))
    # cusp_reg_tail returns critical L-values
    tail = cusp_reg_tail(4, 0//1, 3, prec)
    @test length(tail) == 2
    @test Arblib.overlaps(tail[1], dirichlet_L(chi_minus(4), 2, prec))   # L(χ₋₄,2)=Catalan
    @test Arblib.overlaps(tail[2], Arb(π; prec) / 4)                      # L(χ₋₄,1)=π/4
end

@testset "Weight-3 Eisenstein cusp constants (Γ₀(6) intermediate cusps)" begin
    prec = 256
    # 1. Weight-2 modular-symbol / residue relation: Σ_cusps width·a₀(E₂^{(N)}) = 0.
    for N in (3, 4, 6)
        s = sum(h * a0 for (cusp, (h, a0)) in E._E2N_CUSPDATA[N])
        @test s == 0
    end
    # 2. Hard-tabled (A,B) cusp constants agree with the slash/eta evaluation
    #    `cusp_constant_E_k_chi` at every cusp (including the new 1/2, 1/3).
    for N in (3, 4, 6)
        cond, tbl = E._E3chi_CUSPDATA[N]
        rsA = E._E3chi_ETAREP[cond][:A]
        rsB = E._E3chi_ETAREP[cond][:B]
        for (cusp, (h, fA, fB)) in tbl
            vA = cusp_constant_E_k_chi(3, rsA, cusp, prec; width = h)
            vB = cusp_constant_E_k_chi(3, rsB, cusp, prec; width = h)
            @test Arblib.overlaps(vA, fA(prec))
            @test Arblib.overlaps(vB, fB(prec))
        end
    end
    # 3. {∞,0} cross-check against the generalised-Bernoulli formula:
    #    E₃(χ,𝟏) = L(χ,-2)·A  ⇒  a₀(E₃(χ,𝟏); ∞) = -B_{3,χ}/3 · a₀(A;∞).
    for (N, cond) in ((3,3), (4,4), (6,3))
        ec = eisenstein_constants(N, prec)
        B3 = generalized_bernoulli(3, chi_minus(cond), prec)
        @test Arblib.overlaps(ec[(3, 1//0)], -B3/3 * real(ec[(3, :A, 1//0)]))
        # at cusp 0, A vanishes and B is purely imaginary (odd-character phase)
        @test Arblib.contains_zero(ec[(3, :A, 0//1)])
        @test Arblib.contains_zero(real(ec[(3, :B, 0//1)]))
    end
    # 4. Γ₀(6) intermediate cusps explicitly: a₀(B;1/2) = -a₀(B;0) = -i/3^{9/2},
    #    a₀(A;1/3) = a₀(A;∞) = 1, a₀(B;1/3) = 0.
    ec6 = eisenstein_constants(6, prec)
    inv392 = Acb(0, 1; prec) / (81 * sqrt(Arb(3; prec)))
    @test Arblib.overlaps(ec6[(3, :B, 0//1)],  inv392)
    @test Arblib.overlaps(ec6[(3, :B, 1//2)], -inv392)
    @test Arblib.overlaps(ec6[(3, :A, 1//3)], Acb(1; prec))
    @test Arblib.contains_zero(ec6[(3, :B, 1//3)])
    @test Arblib.contains_zero(ec6[(3, :A, 1//2)])
    # 5. The :B eta quotient for χ₋₃ IS the package's eisenstein_e3_qseries.
    ctx = SeriesContext(prec, 30)
    e3 = E.eisenstein_e3_qseries(ctx)
    bs = E.eta_quotient_qseries(ctx, E._E3chi_ETAREP[3][:B])
    @test all(Arblib.overlaps(E.qcoeff(e3, n), E.qcoeff(bs, n)) for n in 0:20)
end

@testset "Generalised Bernoulli: general-k recursion" begin
    prec = 256
    # evaluate B_k at a rational point through the recursion coefficients
    ev(k, x) = begin
        co = E._bernoulli_poly_coeffs(k)
        acc = Arb(co[1]; prec)
        for j in 1:k; acc = acc * x + Arb(co[j+1]; prec); end
        acc
    end
    x = Arb(3//7; prec)
    # 1. matches the k = 2, 3, 4 closed forms
    @test Arblib.overlaps(ev(2, x), x^2 - x + Arb(1//6; prec))
    @test Arblib.overlaps(ev(3, x), x^3 - 3*x^2/2 + x/2)
    @test Arblib.overlaps(ev(4, x), x^4 - 2*x^3 + x^2 - Arb(1//30; prec))
    # 2. defining identities at higher k (theorems, so exact containment):
    #    translation B_k(x+1) - B_k(x) = k·x^{k-1}, reflection B_k(1-x) = (-1)^k B_k(x)
    for k in (5, 8, 11)
        @test Arblib.contains_zero(ev(k, x + 1) - ev(k, x) - k * x^(k - 1))
        @test Arblib.contains_zero(ev(k, 1 - x) - (-1)^k * ev(k, x))
    end
    # 3. character sums: for the odd characters χ₋₃, χ₋₄ the generalised
    #    Bernoulli numbers vanish at even k ≥ 2 (parity), stay nonzero at odd k,
    #    and B_{1,χ} = (1/f)·Σ_a a·χ(a) (Σχ = 0 kills the -1/2 shift).
    for cond in (3, 4)
        chi = chi_minus(cond)
        for k in (2, 4, 6)
            @test Arblib.contains_zero(generalized_bernoulli(k, chi, prec))
        end
        for k in (3, 5, 7)
            @test !Arblib.contains_zero(generalized_bernoulli(k, chi, prec))
        end
        b1 = Arb(sum(a * chi[a] for a in 1:cond) // cond; prec)
        @test Arblib.overlaps(generalized_bernoulli(1, chi, prec), b1)
    end
    # 4. the documented anchor value B_{3,χ₋₄} = 3/2 (see the weight-3 testset)
    @test Arblib.overlaps(generalized_bernoulli(3, chi_minus(4), prec), Arb(3//2; prec))
    # 5. honest domain fence
    @test_throws ErrorException generalized_bernoulli(0, chi_minus(4), prec)
end

end

@testset "Abel map: general quartic vs ACKM closed form (2-torsion shift)" begin
    prec = 768
    s = Nemo.QQ(2)
    C = QuarticCurve([Nemo.QQ(0), -4*s^2, s^2-4*s, 2*s, Nemo.QQ(1)], prec)
    ord = sortperm([Float64(real(r)) for r in C.roots])
    pair = (ord[2], ord[3])
    Z1 = E.incomplete_edge_quadrature(C, pair, Acb(-3//2; prec)) / period_quadrature(C, pair)
    Z2 = E.abel_image_deformed(Acb(2; prec), Acb(-1//2; prec); delta_exp = 700)
    ψ1, ψ2 = curve_periods(C)
    d = Z1 + Z2 - (1 + ψ2/ψ1)/2
    @test Arblib.is_finite(d)
    @test Float64(Arblib.ubound(Arb, abs(d))) < 2.0^-690
end

include("test_vop_transport.jl")
include("test_cy_transport.jl")
include("test_pm_periods.jl")
include("test_siegel.jl")

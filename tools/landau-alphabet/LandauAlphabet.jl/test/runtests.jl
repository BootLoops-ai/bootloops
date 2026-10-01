using Test
using LandauAlphabet
using Nemo
using Arblib

@testset "LandauAlphabet" begin

@testset "GPL evaluator vs classical polylogs" begin
    prec = 350
    ctx = GPLContext([Acb(0, prec = prec), Acb(1, prec = prec), Acb(-1, prec = prec)];
                     nterms = 600, prec = prec)
    x = Acb(1 // 2, prec = prec)
    tol = Arb(BigFloat(10)^(-95))
    li(s, z) = Arblib.polylog!(Acb(prec = prec), Acb(s, prec = prec), z)
    @test abs(gpl(ctx, [1, 2], x) + li(2, x)) < tol            # G(0,1) = -Li2
    @test abs(gpl(ctx, [2, 2], x) - log(1 - x)^2 / 2) < tol      # G(1,1)
    @test abs(gpl(ctx, [1, 1, 2], x) + li(3, x)) < tol          # G(0,0,1) = -Li3
    @test abs(gpl(ctx, [1, 3], x) + li(2, -x)) < tol            # G(0,-1) = -Li2(-x)
    @test abs(gpl(ctx, [1, 1], x) - log(x)^2 / 2) < tol          # trailing-zero reg
    # shuffle: G(1)G(0) = G(1,0)+G(0,1)
    @test abs(gpl(ctx, [2, 1], x) -
              (gpl(ctx, [2], x) * gpl(ctx, [1], x) - gpl(ctx, [1, 2], x))) < tol
end

@testset "basis tower counts" begin
    @test length(basis_tower(3, 2; variant = :full)) == 17      # 3 letters, weight <= 2
    @test length(basis_tower(3, 3; variant = :full)) == 58      # 3 letters, weight <= 3
end

@testset "exact modular identities (1704.08895)" begin
    sd = sunrise_data(140)
    N = sd.N
    eq(a, b) = (m = min(length(a), length(b)); a[1:m] == b[1:m])
    # t-series
    @test sd.t[1:8] == QQFieldElem[0, -9, -36, -90, -180, -351, -684, -1260]
    # curve side == Eisenstein side (level 12 in q2)
    eis_f2 = qs_add(qs_add(qs_scal(14, B2K(2, N)), qs_scal(-4, B2K(3, N))),
                    qs_add(qs_add(qs_scal(-8, B2K(4, N)), qs_scal(10, B2K(6, N))),
                           qs_scal(-4, B2K(12, N))))
    @test eq(sd.f2, eis_f2)
    @test eq(sd.f2, qs_add(qs_add(sd.g21, sd.g29), qs_scal(QQ(-1, 2), sd.g20)))
    eis_f3 = qs_add(qs_add(eisenstein_E(3, -3, 1, 1, N), qs_scal(2, eisenstein_E(3, -3, 1, 2, N))),
                    qs_scal(-8, eisenstein_E(3, -3, 1, 4, N)))
    @test eq(sd.f3hat, eis_f3)
    tp3 = copy(sd.t); tp3[1] += 3
    @test eq(qs_mul(tp3, sd.u),
             qs_scal(-18, qs_add(eisenstein_E(1, 1, -3, 1, N),
                                 qs_scal(-2, eisenstein_E(1, 1, -3, 4, N)))))
    @test eq(sd.g31hat, qs_scal(-81, eisenstein_E(3, -3, 1, 2, N)))
    f1q = qs_mul(tp3, sd.u)
    @test eq(sd.f4, qs_scal(QQ(1, 324), qs_mul(qs_mul(f1q, f1q), qs_mul(f1q, f1q))))
end

@testset "oracle (iii) weight 2: Landau dlogs == M2(Gamma1(6)), exact" begin
    sd = sunrise_data(140)
    n2, b2 = level6_basis(sd, 2)
    c2n, c2 = landau_weight2_candidates(sd)
    sols = [solve_in_basis(b2, c; nsolve = 20, ncheck = min(130, length(c))) for c in c2]
    @test all(s -> s !== nothing, sols)
    @test sols[1] == QQFieldElem[-8, -4, 8]
    @test sols[2] == QQFieldElem[-9, -3, 3]
    @test sols[3] == QQFieldElem[-5, 5, -1]
    # and they span M2: each basis element solves in the dlogs
    @test all(b -> solve_in_basis(c2, b; nsolve = 20, ncheck = 130) !== nothing, b2)
end

@testset "oracle (iii)/(iv) structure: u^k t^d families == M3, M4" begin
    sd = sunrise_data(140)
    _, b3 = level6_basis(sd, 3)
    _, c3 = landau_weight3_candidates(sd)
    @test length(intersect_with_space(c3, b3)) == 4   # candidates ⊆ M3 (containment)
    # spanning direction (containment count alone does not establish
    # "exactly M3/M4" — solve each basis element in the candidates)
    @test all(b -> solve_in_basis(c3, b; nsolve = 20, ncheck = 130) !== nothing, b3)
    _, b4 = level6_basis(sd, 4)
    _, c4 = landau_weight4_candidates(sd)
    @test length(intersect_with_space(c4, b4)) == 5   # candidates ⊆ M4 (containment)
    @test all(b -> solve_in_basis(c4, b; nsolve = 20, ncheck = 130) !== nothing, b4)
    # f3hat and f4 exact level-6 coordinates
    @test solve_in_basis(b3, qs_substitute_negq(sd.f3hat); nsolve = 20, ncheck = 130) ==
          QQFieldElem[-1, 8, 0, 0]
    @test solve_in_basis(b4, qs_substitute_negq(sd.f4); nsolve = 20, ncheck = 130) ==
          QQFieldElem[6, 0, 54, 0, 0]   # ZERO cusp-form component
    # weight-3 selection rule: the pole-free density (p = t(t-1)(t-9)) is f3hat
    # candidates are u^3 t^d / 9, d = 0..3; f3hat = -(cand2 + ... ) checked above
end

@testset "integrable symbols: known alphabet counts" begin
    ladder = symbol_alphabet(["z", "zb"], ["z", "zb", "1-z", "1-zb", "z-zb"],
        [("zzb", [1, 1, 0, 0, 0]), ("(1-z)(1-zb)", [0, 0, 1, 1, 0]),
         ("z-zb", [0, 0, 0, 0, 1]), ("zb/z", [-1, 1, 0, 0, 0]),
         ("(1-z)/(1-zb)", [0, 0, 1, -1, 0])])
    @test ncols(integrable_symbols(ladder, 1)[2]) == 5
    @test ncols(integrable_symbols(ladder, 2)[2]) == 19   # 5 at weight 1 + 19 at weight 2
    seven = symbol_alphabet(["x", "z"], ["x", "z", "1+x", "1-x", "1-z", "x+z", "1+x*z"],
        [(n, [i == j ? 1 : 0 for i in 1:7]) for (j, n) in
         enumerate(["x", "z", "1+x", "1-x", "1-z", "x+z", "1+xz"])])
    @test ncols(integrable_symbols(seven, 2)[2]) == 37    # full w2 basis: 37+7+1+8 = 53
    @test ncols(integrable_symbols(seven, 3)[2]) == 175   # full w3 basis: 221+53 = 274
end

@testset "constraint assembly" begin
    cs = ConstraintSystem(3)
    add_constraint!(cs, QQFieldElem[QQ(1), QQ(1), QQ(0)], QQ(2))
    add_constraint!(cs, QQFieldElem[QQ(0), QQ(1), QQ(-1)], QQ(0))
    part, nulls = reduce_space(cs)
    @test length(nulls) == 1
    @test part[1] + part[2] == 2 && part[2] == part[3]
end

@testset "lattice fit recovers planted relation" begin
    prec = 256
    # f = 3*B1 - (5/2)*B2 on fake "function values"
    p = 6
    B = Matrix{Arb}(undef, 2, p)
    f = Vector{Arb}(undef, p)
    for j in 1:p
        B[1, j] = sin(Arb(j, prec = prec))
        B[2, j] = exp(Arb(j, prec = prec) / 7)
        f[j] = 3 * B[1, j] - Arb(5 // 2, prec = prec) * B[2, j]
    end
    fr = lattice_fit(f, B; digits = 25)
    @test fr.ok
    @test fr.coeffs == [QQ(3), QQ(-5, 2)]
    @test Float64(fit_residual(fr, f, B)) < 1e-20
end

@testset "integer relation (PSLQ-equivalent)" begin
    prec = 512
    P = Arb(π, prec = prec)
    z3 = Arblib.zeta!(Arb(prec = prec), Arb(3, prec = prec))
    v = [P^2 / 6 + z3, P^2, z3]          # v1 = v2/6 + v3
    ok, rel, res = integer_relation(v; digits = 60)
    @test ok
    # v1 = v2/6 + v3  ⇒  rel ∝ (6, -1, -6)
    @test 6 * rel[2] == -rel[1] && rel[3] == -rel[1]
end

@testset "boundary constants via Cauchy extraction" begin
    # S^(0)(2,0) = (sqrt3/i)[Li2(r3)-Li2(r3inv)] = 2*sqrt3*Im Li2(r3) = 2*sqrt3*Cl2(2pi/3)
    prec = 1024
    cs = sunrise_boundary_constants(2; prec = prec, M = 96, r = 1 // 16)
    r3 = exp(Acb(0, 2, prec = prec) * Arb(π, prec = prec) / 3)
    li2 = Arblib.polylog!(Acb(prec = prec), Acb(2, prec = prec), r3)
    target = 2 * sqrt(Arb(3, prec = prec)) * imag(li2)
    @test abs(real(cs[1]) - target) < Arb(BigFloat(10)^(-50))
    @test abs(imag(cs[1])) < Arb(BigFloat(10)^(-50))
end

end

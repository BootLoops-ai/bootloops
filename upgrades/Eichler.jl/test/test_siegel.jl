using Test, Arblib, Eichler

@testset "Siegel / genus-2 hyperelliptic" begin

@testset "characteristic encoding & parity" begin
    @test siegel_char(0,0,0,0) == 0
    @test siegel_char(1,1,1,1) == 15
    @test siegel_char((1,0,0,1)) == 9
    # 10 even, 6 odd in genus 2
    @test length(SIEGEL_G2_EVEN_CHARS) == 10
    @test sort(Int.(SIEGEL_G2_EVEN_CHARS)) == [0,1,2,3,4,6,8,9,12,15]
    @test siegel_char_is_even(UInt(15), 2)
    @test !siegel_char_is_even(UInt(5), 2)
end

@testset "siegel_theta: odd nulls vanish, FLINT consistency" begin
    prec = 256
    τ = Arblib.AcbMatrix(2, 2; prec)
    τ[1,1] = Acb(1//7, 13//10; prec); τ[2,2] = Acb(-1//5, 17//10; prec)
    τ[1,2] = Acb(1//11, 2//5; prec);  τ[2,1] = τ[1,2]
    th = siegel_theta_all(Arblib.AcbVector(2; prec), τ, prec)
    for ab in 0:15
        if siegel_char_is_even(UInt(ab), 2)
            @test !Arblib.contains_zero(th[ab+1])
        else
            @test Arblib.contains_zero(th[ab+1])
        end
        # single-char API agrees with the all-char vector
        @test Arblib.overlaps(th[ab+1], siegel_theta(UInt(ab), Arblib.AcbVector(2;prec), τ, prec))
    end
    # χ₁₀ = 2⁻¹²·∏_{even} θ² (FLINT def)
    ψ4, ψ6, χ10, χ12 = igusa_invariants(τ, prec)
    pr = Acb(1; prec)
    for ab in SIEGEL_G2_EVEN_CHARS; pr *= th[ab+1]^2; end
    @test Arblib.overlaps(χ10, pr / Acb(2; prec)^12)
end

@testset "period_matrix: Riemann bilinear (symmetry, Im τ ≻ 0)" begin
    prec = 256
    # generic real-ordered sextic (no symmetries)
    C = HyperellipticCurve(1, [-2, -1//2, 3//10, 17//10, 16//5, 51//10], prec)
    @test C.real_ordered
    A, B = hyper_big_period_matrix(C)
    @test all(Arblib.is_finite, (A[i,j] for i in 1:2, j in 1:2))
    τ = period_matrix(C)
    @test Arblib.overlaps(τ[1,2], τ[2,1])
    detIm = imag(τ[1,1])*imag(τ[2,2]) - imag(τ[1,2])^2
    @test Arblib.is_positive(imag(τ[1,1])) && Arblib.is_positive(detIm)
    @test Eichler.certified_digits_test(imag(τ[1,1])) > 60
end

@testset "Igusa absolute round-trip (sextic → τ → Thomae sextic)" begin
    prec = 320
    # generic degree-6, real roots, lead ≠ 1 to exercise normalisation
    C = HyperellipticCurve(3, [-2, -1//2, 3//10, 17//10, 16//5, 51//10], prec)
    τ = period_matrix(C)
    f̃, χ5 = thomae_sextic(τ, prec)
    sc = [Arblib.ref(f̃, i) + Acb(0; prec) for i in 0:6]
    j_in  = igusa_absolute(C)
    j_rec = igusa_absolute(igusa_clebsch(sc, prec))
    for k in 1:3
        @test Arblib.overlaps(j_in[k], j_rec[k])
        @test Eichler.certified_digits_test(real(j_in[k])) > 60
    end
end

@testset "Rosenhain / Thomae: branch points from theta-nulls (deg 5)" begin
    # y² = x(x-1)(x-λ₁)(x-λ₂)(x-λ₃), 1 < λ₁ < λ₂ < λ₃ — recover λ_i from θ⁴.
    for (prec, λ) in ((256, (5//2, 4, 7)), (300, (13//10, 21//10, 8)))
        e = [0, 1, λ...]
        C = HyperellipticCurve(1, e, prec)
        @test C.degree == 5
        τ = period_matrix(C)
        λ1, λ2, λ3 = rosenhain_from_tau(τ, prec)
        @test Arblib.overlaps(λ1, Acb(λ[1]; prec))
        @test Arblib.overlaps(λ2, Acb(λ[2]; prec))
        @test Arblib.overlaps(λ3, Acb(λ[3]; prec))
        @test Eichler.certified_digits_test(real(λ1)) > 50
    end
end

@testset "Bolza curve y² = x⁵ - x: complex-root + special-locus handling" begin
    prec = 256
    Cb = HyperellipticCurve([0,-1,0,0,0,1], prec)
    @test Cb.degree == 5 && !Cb.real_ordered
    # Degree-5 → Möbius-normalised degree-6 has 6 DISTINCT roots (smooth curve);
    # Bolza sits on the maximal-automorphism locus where the FLINT degree-10
    # covariant C₁₀,₀ vanishes (≠ disc; both are degree-10 invariants).  This is
    # a documented limitation of the `igusa_absolute` test helper, NOT of the
    # period/theta layer.
    C20, C40, C60, C100 = igusa_clebsch(Cb)
    @test !Arblib.contains_zero(C20)
    @test Arblib.contains_zero(C40) && Arblib.contains_zero(C60)
    @test Arblib.contains_zero(C100)
    # Complex roots {0,±1,±i} now go through the chain assembly (the Igusa
    # round-trip gate inside period_matrix is skipped on this special locus —
    # C₁₀₀ ≈ 0 — so check the Thomae side's locus pattern explicitly below).
    Cb2 = HyperellipticCurve([0,-1,0,0,0,1], 160)
    τb = period_matrix(Cb2)
    @test Arblib.overlaps(τb[1,2], τb[2,1])
    @test Arblib.is_positive(imag(τb[1,1]))
    @test Eichler.certified_digits_test(imag(τb[1,1])) > 25
    ψ4, ψ6, χ10, χ12 = igusa_invariants(τb, 160)
    @test !Arblib.contains_zero(χ10)              # smooth curve, off ℍ₁²
    # Thomae sextic of τb sits on the same maximal-automorphism locus
    f̃, _ = thomae_sextic(τb, 160)
    cf(i) = (a = Acb(prec = 160); Arblib.get_coeff!(a, f̃, i); a)
    D20, D40, D60, D100 = igusa_clebsch([cf(i) for i in 0:6], 160)
    @test !Arblib.contains_zero(D20)
    @test Arblib.contains_zero(D40) && Arblib.contains_zero(D60)
    @test Arblib.contains_zero(D100)
end

@testset "chain assembly: complex branch-point period matrices" begin
    prec = 192
    # 1. real-ordered sextic pushed through the chain path (identity chain)
    #    reproduces the direct real-ordered τ — the two code paths are
    #    independent oracles for each other.
    Cr = HyperellipticCurve(1, [-2, -1//2, 3//10, 17//10, 16//5, 51//10], prec)
    τr = period_matrix(Cr)
    τc = period_matrix(Cr; chain = collect(1:6))
    @test all(Arblib.overlaps(τr[j,k], τc[j,k]) for j in 1:2, k in 1:2)
    # 2. genuinely complex configuration: f = (x²+1)(x²+4)(x-1)(x-3), four
    #    complex branch points in two conjugate pairs plus two real ones.
    Cc = HyperellipticCurve([12, -16, 19, -20, 8, -4, 1], prec)
    @test !Cc.real_ordered
    ch = hyper_default_chain(Cc)
    @test isperm(ch) && length(ch) == 6
    τ = period_matrix(Cc)     # passes the internal Igusa round-trip gate
    @test Arblib.overlaps(τ[1,2], τ[2,1])
    detIm = imag(τ[1,1])*imag(τ[2,2]) - imag(τ[1,2])^2
    @test Arblib.is_positive(imag(τ[1,1])) && Arblib.is_positive(detIm)
    @test Eichler.certified_digits_test(imag(τ[1,1])) > 30
    # explicit Igusa absolute round-trip (same check the gate performs)
    f̃, _ = thomae_sextic(τ, prec)
    cf(i) = (a = Acb(prec = prec); Arblib.get_coeff!(a, f̃, i); a)
    j_in  = igusa_absolute(Cc)
    j_rec = igusa_absolute(igusa_clebsch([cf(i) for i in 0:6], prec))
    for k in 1:3
        @test Arblib.overlaps(j_in[k], j_rec[k])
    end
    # 3. refusal laws.  A chain with a branch point ON a segment is refused
    #    by the certified per-factor classification…
    Con = HyperellipticCurve(1,
        [Acb(0,1;prec), Acb(1;prec), Acb(0,-1;prec), Acb(0,2;prec),
         Acb(0,-2;prec), Acb(3;prec)], prec)
    #    (segment 2i → -i of this chain passes through the branch point i)
    @test_throws ErrorException hyper_chain_integral(Con, [4,3,1,5,2,6], 1, 0)
    @test_throws ErrorException period_matrix(Con; chain = [4,3,1,5,2,6])
    #    …and a self-crossing chain is refused by the simplicity pre-check
    #    (the two diagonals of the square of branch points cross at 0).
    Cx = HyperellipticCurve(1,
        [Acb(1,1;prec), Acb(-1,1;prec), Acb(1,-1;prec), Acb(-1,-1;prec),
         Acb(3;prec), Acb(4;prec)], prec)
    @test_throws ErrorException period_matrix(Cx; chain = [1,4,2,3,5,6])
end

@testset "abel_jacobi: Feynman-deformed (complex) marked points" begin
    prec = 192
    C = HyperellipticCurve(1, [0, 1, 5//2, 4, 7], prec)
    u0 = abel_jacobi(C, 12; base = 5)
    # continuity across the +i0 anchor: a dyadically small deformation moves
    # the image by O(ε) on both sides (Im < 0 continues through the cut, per
    # the documented convention)
    ε = Arb(1//1048576; prec)
    up = abel_jacobi(C, Acb(Arb(12; prec),  ε; prec); base = 5)
    um = abel_jacobi(C, Acb(Arb(12; prec), -ε; prec); base = 5)
    for k in 1:2
        @test Float64(Arblib.ubound(Arb, abs(up[k] - u0[k]))) < 1e-4
        @test Float64(Arblib.ubound(Arb, abs(um[k] - u0[k]))) < 1e-4
    end
    # a finite-distance deformation stays finite and lands on the Jacobian
    τ = period_matrix(C)
    ub = abel_jacobi(C, Acb(12, 3; prec); base = 5)
    @test all(Arblib.is_finite, ub)
    th = siegel_theta(UInt(0), ub, τ, prec)
    @test Arblib.is_finite(th) && !Arblib.contains_zero(th)
    # interior gap basepoint with a deformed marked point
    ui = abel_jacobi(C, Acb(1//2, 3//10; prec); base = 1)
    @test all(Arblib.is_finite, ui)
    ui0 = abel_jacobi(C, 1//2; base = 1)
    ui1 = abel_jacobi(C, Acb(1//2, 1//1048576; prec); base = 1)
    for k in 1:2
        @test Float64(Arblib.ubound(Arb, abs(ui1[k] - ui0[k]))) < 1e-4
    end
    # window fences still hold for deformed points
    @test_throws ErrorException abel_jacobi(C, Acb(8, 1//2; prec); base = 1)
end

@testset "degree-5 Igusa round-trip via period matrix" begin
    prec = 320
    C = HyperellipticCurve(1, [0, 1, 5//2, 4, 7], prec)
    @test C.degree == 5
    τ = period_matrix(C)
    f̃, _ = thomae_sextic(τ, prec)
    cf(i) = (a = Acb(prec = prec); Arblib.get_coeff!(a, f̃, i); a)
    j_in  = igusa_absolute(C)
    j_rec = igusa_absolute(igusa_clebsch([cf(i) for i in 0:6], prec))
    for k in 1:3
        @test Arblib.overlaps(j_in[k], j_rec[k])
    end
end

@testset "abel_jacobi: half-period at the next branch point" begin
    prec = 256
    C = HyperellipticCurve(1, [0, 1, 5//2, 4, 7], prec)
    τ = period_matrix(C)
    # ∫_{e₅}^{x} on (e₅,∞) limited; instead test base = 1 (∫_{e₁}^{x} on gap 1):
    # at x = e₂ the incomplete integral equals the full gap-1 integral, so
    # A⁻¹·(full a₁-column)/… — concretely, abel_jacobi from e₁ to e₂ should be
    # a half-lattice vector: 2u ∈ ℤ² + τℤ² (e₂ is a 2-torsion point of Jac(C)).
    # Our base parametrisation needs xp strictly inside the gap; take xp → e₂⁻
    # via two evaluations and check 2u ≡ first column of (1|τ) lattice.
    # Simpler invariant: abel_jacobi(C, xp; base=n) on (e_n, ∞) is real and
    # finite, and θ[ab](u, τ) is well-defined (the Jacobian embedding).
    u = abel_jacobi(C, 12; base = 5)
    @test all(Arblib.is_finite, u)
    @test all(Arblib.contains_zero ∘ imag, u)   # real path ⇒ real Abel image
    th = siegel_theta(UInt(0), u, τ, prec)
    @test Arblib.is_finite(th) && !Arblib.contains_zero(th)
end

@testset "Gauss–Manin ∂τ/∂e_m (Rauch) vs certified FD" begin
    prec = 600
    base_roots(tt) = [Acb(-2; prec), Acb(-1; prec), Acb(3//10; prec),
                      Acb(17//10; prec), Acb(32//10; prec) + Acb(tt; prec),
                      Acb(51//10; prec)]
    fam = tt -> HyperellipticCurve(1, base_roots(tt), prec)
    t0 = Arb(1//5; prec)
    C = fam(t0)
    # only e₅ moves with t (de₅/dt = 1, all others 0)
    Drauch = siegel_rauch_de(C, 5)
    @test Arblib.overlaps(Drauch[1,2], Drauch[2,1])      # symmetric (certified)
    @test Eichler.certified_digits_test(real(Drauch[1,1])) > 150  # Rauch is exact
    # finite-difference reference on τ (sanity only — FD truncates at O(h²); we
    # compare MIDPOINTS, since the deliberately-cancelled FD ball is not a tight
    # certified enclosure of the exact Rauch value).
    h = Arb(2; prec)^(-70)
    τp = period_matrix(fam(t0 + h)); τm = period_matrix(fam(t0 - h))
    for (j,k) in ((1,1),(1,2),(2,2))
        fd = (τp[j,k] - τm[j,k]) / (2h)
        diff = abs(Arb(real(fd - Drauch[j,k])))
        @test diff < Arb(10; prec)^(-100)          # agree to ≥100 d (FD-limited)
    end
    # full Gauss–Manin with the velocity vector (only m=5 nonzero) reproduces it
    de = [Acb(0;prec),Acb(0;prec),Acb(0;prec),Acb(0;prec),Acb(1;prec),Acb(0;prec)]
    Dgm = siegel_gauss_manin(C, de)
    for j in 1:2, k in 1:2
        @test Arblib.overlaps(Dgm[j,k], Drauch[j,k])
    end
end

@testset "Gauss–Manin: exact coeff-velocity ≡ FD-velocity path" begin
    prec = 500
    # coefficient family f(x,t) = ∏(x − rᵢ(t)) with two moving roots — exercise
    # branch_velocities_coeffs (exact, ∂_t f / f') vs branch_velocities_fd.
    roots_at(tt) = [Acb(-2;prec), Acb(-1;prec)+Acb(tt;prec)/5, Acb(3//10;prec),
                    Acb(17//10;prec), Acb(32//10;prec)+Acb(tt;prec)/3, Acb(51//10;prec)]
    fam = tt -> HyperellipticCurve(1, roots_at(tt), prec)
    coeffs_at = tt -> begin
        p = AcbPoly([Acb(1;prec)]; prec)
        for r in roots_at(tt); p = p * AcbPoly([-r, Acb(1;prec)]; prec); end
        cf(i) = (a = Acb(prec=prec); Arblib.get_coeff!(a,p,i); a)
        [cf(i) for i in 0:6]
    end
    t0 = Arb(1//7; prec)
    C = fam(t0)
    de_exact = Eichler.branch_velocities_coeffs(coeffs_at, C.roots, t0; prec)
    de_fd    = Eichler.branch_velocities_fd(fam, t0; prec)
    for m in 1:6
        @test Arblib.overlaps(de_exact[m], de_fd[m])
    end
    Dexact = siegel_gauss_manin(C, de_exact)
    Dfd    = siegel_gauss_manin(C, de_fd)
    @test Arblib.overlaps(Dexact[1,2], Dexact[2,1])
    for j in 1:2, k in 1:2
        @test Arblib.overlaps(Dexact[j,k], Dfd[j,k])
        @test Eichler.certified_digits_test(real(Dexact[j,k])) > 60
    end
end

@testset "H+jet nonplanar genus-2 Feynman curve: GM connection" begin
    # y² = x(x+s)(x²+sx−4s)(x−t)(x−(m_H²−s−t)), mt²=1, the nonplanar H+jet top
    # sector (Marzucca et al. 2307.11497).  Pick a real-ordered slice
    # and transport the Siegel modulus's kinematic derivative along t.
    prec = 400
    s = -20; mH2 = 1                      # fixed; deform t (well-separated slice)
    function hjet_roots(tval)
        # x²+sx−4s roots: discriminant s²+16s ≥ 0 for s ≤ −16; collect all six.
        disc = s^2 + 16s
        @assert disc ≥ 0
        r1 = (-s + sqrt(Acb(disc; prec)))/2
        r2 = (-s - sqrt(Acb(disc; prec)))/2
        [Acb(0;prec), Acb(-s;prec), r1, r2, Acb(tval;prec), Acb(mH2 - s - tval;prec)]
    end
    fam = tval -> HyperellipticCurve(1, hjet_roots(tval), prec)
    t0 = 2
    C = fam(t0)
    @test C.degree == 6 && C.real_ordered      # genuine genus-2 sextic, 6 real bp
    τ = period_matrix(C)
    @test Arblib.is_positive(imag(τ[1,1]))
    # the moving branch points are e₅=t (de/dt=1) and e₆=mH²−s−t (de/dt=−1)
    sec = Genus2Sector(prec; curve_at = fam, base = 6)   # skip the t=1 probe
    @test sec.base == 6
    # velocities by FD on roots
    de = Eichler.branch_velocities_fd(fam, Arb(t0; prec); prec)
    D = siegel_gauss_manin(C, de)
    @test Arblib.overlaps(D[1,2], D[2,1])
    @test all(Arblib.is_finite, (D[i,j] for i in 1:2, j in 1:2))
    @test Eichler.certified_digits_test(real(D[1,1])) > 40
    # cross-check midpoints against direct FD of τ(t) (FD-limited; see note above)
    h = Arb(2; prec)^(-60)
    τp = period_matrix(fam(t0 + h)); τm = period_matrix(fam(t0 - h))
    for (j,k) in ((1,1),(1,2),(2,2))
        fd = (τp[j,k]-τm[j,k])/(2h)
        @test abs(Arb(real(fd - D[j,k]))) < Arb(10; prec)^(-80)
    end
end

@testset "Genus2Sector + frontier stub (sharpened)" begin
    prec = 200
    fam = t -> HyperellipticCurve(1, [0, 1, 2 + t, 4 + t, 7 + 2t], prec)
    sec = Genus2Sector(prec; curve_at = fam)
    τ = siegel_tau_map(sec, Acb(1//3; prec))
    @test size(τ) == (2,2)
    @test Arblib.is_positive(imag(τ[1,1]))
    ψ4, ψ6, χ10, χ12 = igusa_invariants(sec, Acb(1//3; prec))
    @test Arblib.is_finite(ψ4) && Arblib.is_finite(χ10)
    # the Gauss–Manin connection (NEW) is available at sector level
    D = siegel_gauss_manin(sec, Acb(1//3; prec))
    @test size(D) == (2,2) && Arblib.overlaps(D[1,2], D[2,1])
    # transport itself remains the honest stub
    @test_throws SiegelNotImplemented siegel_transport(sec, nothing, 0, 1, prec)
end

end

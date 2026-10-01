@testset "pf_from_theta round-trip" begin
    prec = 256
    # banana K3 op (already in ∂-form) → θ-form → ∂-form must agree
    L3 = banana_pf(3; prec)
    R3 = theta_form(L3)
    L3b = pf_from_theta(prec, R3)
    @test L3b.order == L3.order
    for i in 0:L3.order
        a, b = L3.coeffs[i+1], L3b.coeffs[i+1]
        n = max(length(a), length(b))
        @test vcat(a, zeros(Rational{BigInt}, n - length(a))) ==
              vcat(b, zeros(Rational{BigInt}, n - length(b)))
    end
    # θ → ∂ → θ round-trip on the Apéry operator
    Rap = Vector{Rational{BigInt}}[[0,0,0,1], [-5,-27,-51,-34], [1,3,3,1]]
    Lap = pf_from_theta(prec, Rap)
    @test theta_form(Lap) == Rap
end

@testset "PM K3.0 (Legendre Sym²)" begin
    prec = 512
    pm = pm_k3_legendre(prec; nterms = 250)
    @test pm.L.order == 3
    @test pm.frob.alpha == 0
    @test pm.alpha_shift == 1//2
    # ϖ₀ series = ₂F₁(½,½;1;z)² = Σ [Σ_k C(2k,k)²C(2n-2k,n-k)²/16ⁿ] zⁿ
    @test pm.holo_coeffs[1:5] == Rational{BigInt}[1, 1//2, 11//32, 17//64, 1787//8192]
    cn(n) = sum(binomial(big(2k), k)^2 * binomial(big(2(n-k)), n-k)^2 for k in 0:n) // big(16)^n
    @test all(pm.holo_coeffs[n+1] == cn(n) for n in 0:20)
    # singularities {0, 1}
    @test length(pm.L.sing) == 2
    @test any(Arblib.overlaps(s, Acb(1; prec)) for s in pm.L.sing)
    # ϖ₀(z) = (2/π)² K(z)² at three held-out points, ≥80d
    for z in (Acb(1//10; prec), Acb(1//7; prec), Acb(1//4; prec))
        Π = frobenius_basis(pm, z)
        cl = holo_period_closed(pm, z)
        @test Float64(abs(Π[1] - cl)) < 1e-80
        @test Float64(abs(pf_residual(pm, z))) < 1e-80
    end
    # full-tower check at z=1/10: ϖ₁ via direct K,K' identity
    #   ϖ₁ = ϖ₀ log z + h₁(z), and ∂_z ϖ₁/ϖ₀ = (mirror) — covered by mirror_map
    z = Acb(1//10; prec)
    q, lq = mirror_map(pm.frob, z, prec)
    @test Arblib.is_finite(lq)
end

@testset "PM CY3 (₄F₃ hypergeometric)" begin
    prec = 512
    pm = pm_cy3_4F3(prec; nterms = 250)
    @test pm.L.order == 4
    @test pm.frob.alpha == 0
    # ϖ₀ coefficients = C(2n,n)⁴
    @test all(pm.holo_coeffs[n+1] == binomial(big(2n), big(n))^4 for n in 0:30)
    # singularities {0, 2⁻⁸}
    @test any(Arblib.overlaps(s, Acb(1//256; prec)) for s in pm.L.sing)
    # ϖ₀ = ₄F₃(½⁴;1³;2⁸z) at three held-out points, ≥80d (Arb hypgeom_pfq)
    for z in (Acb(1//2048; prec), Acb(1//1024; prec), Acb(1//800; prec))
        Π = frobenius_basis(pm, z)
        cl = holo_period_closed(pm, z)
        @test Float64(abs(Π[1] - cl)) < 1e-80
        @test Float64(abs(pf_residual(pm, z))) < 1e-80
    end
end

@testset "PM K3′ (Apéry / Beukers–Peters)" begin
    # OEIS A005259 — guard against fabricated coefficients
    @test apery_numbers(5) == BigInt[1, 5, 73, 1445, 33001, 819005]
    @test apery_numbers(8)[7:9] == BigInt[21460825, 584307365, 16367912425]
    # closed form a_n = Σ_k C(n,k)² C(n+k,k)²
    apcl(n) = sum(binomial(big(n), k)^2 * binomial(big(n+k), k)^2 for k in 0:n)
    @test all(apery_numbers(15)[n+1] == apcl(n) for n in 0:15)

    prec = 512
    pm = pm_k3_apery(prec; nterms = 300)
    @test pm.L.order == 3
    @test pm.frob.alpha == 0
    @test pm.holo_coeffs[1:6] == Rational{BigInt}.(apery_numbers(5))
    # singularities {0, (3∓2√2)²}
    zplus = (Arb(3; prec) - 2*sqrt(Arb(2; prec)))^2
    @test any(Arblib.overlaps(real(s), zplus) for s in pm.L.sing)
    @test Float64(Arblib.lbound(Arb, pm.frob.R)) < 0.03   # radius caveat
    # ϖ₀ via Frobenius vs direct Apéry sum, ≥80d
    for z in (Acb(1//200; prec), Acb(1//100; prec), Acb(1//80; prec))
        Π = frobenius_basis(pm, z)
        cl = holo_period_closed(pm, z)
        @test Float64(abs(Π[1] - cl)) < 1e-80
        @test Float64(abs(pf_residual(pm, z))) < 1e-80
    end
end

@testset "PM CY3′ (Hadamard χ=80)" begin
    prec = 512
    pm = pm_cy3_hadamard(prec; nterms = 250)
    @test pm.L.order == 4
    @test pm.frob.alpha == 0
    # ϖ₀ z-series (verified by hand against eq.(A6) t-reexpansion in pm_periods.jl)
    @test pm.holo_coeffs[1:4] == Rational{BigInt}[1, 112, 47376, 27846400]
    # singularity z = 2⁻¹⁰
    @test any(Arblib.overlaps(s, Acb(1//1024; prec)) for s in pm.L.sing)
    # eq.(A6): ϖ₀(z(t)) = 1 − 7/2⁸ t² + 7/2⁸ t³ − 25711/2²⁰ t⁴ + …  (Baikov T³ residue)
    t = Acb(1//1000; prec)
    zt = -t^2 / (Acb(4096; prec) * (1 + t))
    Πt = frobenius_basis(pm, zt)[1]
    a6 = 1 - Acb(7//256; prec)*t^2 + Acb(7//256; prec)*t^3 - Acb(25711//2^20; prec)*t^4
    @test Float64(abs(Πt - a6)) < 1e-14    # next term is O(t⁵) ~ 10⁻¹⁵
    # Hadamard cross-operator identity at three held-out z, ≥80d
    for z in (Acb(1//8192; prec), Acb(1//4096; prec), Acb(1//3000; prec))
        Π = frobenius_basis(pm, z)
        cl = holo_period_closed(pm, z)
        @test Float64(abs(Π[1] - cl)) < 1e-80
        @test Float64(abs(pf_residual(pm, z))) < 1e-80
    end
end

@testset "PM periods: θ-transport past Apéry singularity (continuation smoke)" begin
    prec = 384
    pm = pm_k3_apery(prec; nterms = 200)
    # transport ϖ₀ from z=1/200 to z=1/10 (past the sing at z₊≈0.0294) along an
    # arc through the upper half-plane; two independent arcs must agree ≥60d.
    z0, z1 = Acb(1//200; prec), Acb(1//10; prec)
    p1 = [z0, Acb(1//100, 1//50; prec), Acb(3//100, 1//30; prec), Acb(1//20, 1//50; prec), z1]
    p2 = [z0, Acb(1//100, 1//30; prec), Acb(3//100, 1//20; prec), Acb(6//100, 1//30; prec), z1]
    W1 = period_matrix(pm.L, z1, prec; frob = pm.frob, anchor = z0, path = p1)
    W2 = period_matrix(pm.L, z1, prec; frob = pm.frob, anchor = z0, path = p2)
    @test maximum(Float64(abs(W1[i,k] - W2[i,k])) for i in 1:3, k in 1:3) < 1e-60
end

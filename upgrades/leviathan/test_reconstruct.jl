#!/usr/bin/env julia
# test_reconstruct.jl — tests for the finite-field reconstruction route.
# Run: julia +1.10 --project=<dir with Oscar+JSON, e.g. this one> test_reconstruct.jl
include(joinpath(@__DIR__, "src", "Leviathan.jl"))
using .Leviathan, Random, Test

@testset "ratrec dense-ansatz" begin
    R, (a, b) = Leviathan.Oscar.polynomial_ring(Leviathan.Nemo.QQ, ["a","b"])
    # f = (2a-b)/(a^2-3b^2): sample mod p, reconstruct, check factors
    p = P31[1]; Fp = Leviathan.Nemo.GF(p)
    samp = [(Int[k1, k2], Int(Leviathan.Nemo.lift(Leviathan.Nemo.ZZ,
            (Fp(2k1)-Fp(k2)) * inv(Fp(k1)^2 - 3*Fp(k2)^2)))) for k1 in 3:12 for k2 in 3:12]
    Rp,_ = Leviathan.Oscar.polynomial_ring(Fp, ["a","b"])
    r = Leviathan.ratrec_modp(samp, Rp, p; deg_cap=4)
    @test r !== nothing
    @test r[4] == 2   # deg D
end

@testset "bubble reconstruct == exact" begin
    fam = family("x1 + x2 + mm*(x1+x2)^2 - s*x1*x2", ["s","mm"], ["x1","x2"])
    rr = landau_reconstruct(fam; nprimes=3, npoints=20, deg_cap=6, lanes=(:nu0,:regulated))
    @test Set(rr.genuine) == Set(["mm","s","s - 4*mm"])
end

# Provenance: the family below is the sunrise with equal internal masses (mm1) and external
# invariant MM1; the expected factor set is its standard singular set — normal threshold
# MM1 = 9*mm1, pseudo-threshold MM1 = mm1, and the vanishing loci mm1 = 0, MM1 = 0 —
# independently reproduced here by both routes (exact and reconstruction); see PATCHES.md.
@testset "sunrise (masses mm1, invariant MM1): reconstruct == exact == thresholds" begin
    fam = family("(x1*x2+x1*x3+x2*x3) + (x1*x2+x1*x3+x2*x3)*mm1*(x1+x2+x3) - MM1*x1*x2*x3",
                 ["mm1","MM1"], ["x1","x2","x3"])
    rr = landau_reconstruct(fam; nprimes=3, npoints=40, deg_cap=8, lanes=(:nu0,:regulated))
    re = landau(fam; rng=MersenneTwister(20260630))
    @test Set(rr.genuine) == Set(re.genuine)
    @test Set(rr.genuine) == Set(["MM1","mm1","mm1 - 1//9*MM1","mm1 - MM1"])
end
println("ALL RECONSTRUCT TESTS PASSED")

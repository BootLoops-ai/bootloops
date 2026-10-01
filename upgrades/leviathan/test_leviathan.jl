# test_leviathan.jl — unit + gate tests for Leviathan.
# Run: julia +1.10 --project=<dir with Oscar+JSON, e.g. this one> test_leviathan.jl [--full]
include(joinpath(@__DIR__, "src", "Leviathan.jl"))
using .Leviathan, Random, Test

@testset "staircase" begin
    # unit ideal
    @test Leviathan.staircase_count([[0, 0]], 2) == 0
    # <x^2, y^3> -> 6 standard monomials
    @test Leviathan.staircase_count([[2, 0], [0, 3]], 2) == 6
    # <x^2, xy, y^3> -> {1,x,y,y^2} = 4
    @test Leviathan.staircase_count([[2, 0], [1, 1], [0, 3]], 2) == 4
    # positive-dimensional: no pure power in y
    @test Leviathan.staircase_count([[2, 0], [1, 1]], 2) === nothing
    # DRL comparator: deg first, then last-nonzero-negative wins
    @test Leviathan._drl_less([1, 0], [0, 2])
    @test Leviathan._drl_less([0, 1], [1, 0])   # s=(1,0) > mm=(0,1) in DRL
    # sliced counter: exact for boxes far beyond direct point-walk range
    @test Leviathan.staircase_count([[10^6, 0], [0, 10^6]], 2) == 10^12
    @test Leviathan.staircase_count([[10^9, 0], [1, 1], [0, 10^9]], 2) == 2 * 10^9 - 1
    w = Leviathan.staircase_count([[10^7, 0, 0], [0, 10^7, 0], [0, 0, 10^7]], 3)
    @test w == Int128(10)^21 && w isa Int128    # counts wider than Int come back Int128
end

@testset "parse+family" begin
    fam = family("x1 + x2 + mm*(x1+x2)^2 - s*x1*x2", ["s", "mm"], ["x1", "x2"])
    @test nedges(fam) == 2
    @test admissible_sectors(fam) == [[1], [2], [1, 2]]
    @test iszero(sector_G(fam, Int[]))
end

@testset "bubble gate (exact)" begin
    rng = MersenneTwister(20260630)
    fam = family("x1 + x2 + mm*(x1+x2)^2 - s*x1*x2", ["s", "mm"], ["x1", "x2"])
    r = landau(fam; rng=rng)
    @test r.chi_generic == 3
    @test r.chi_sectors[[1, 2]] == 1 && r.chi_sectors[[1]] == 1 && r.chi_sectors[[2]] == 1
    @test Set(r.genuine) == Set(["s", "mm", "s - 4*mm"])
    @test isempty(r.spurious)
    d = diagnose(fam; rng=rng, candidates=r.candidates)
    # type1.2 correctly FAILS on bubble: s=0 is a nu=0-degenerate (second-type) locus the
    # nu=0 lane provably misses (previously false-PASSed when run on the regulated ideal).
    @test all(v[1] for (k, v) in d if k != "type1.2_nodegeneracy") && !d["type1.2_nodegeneracy"][1]
    # negative control: a fake locus must be flagged spurious (no chi drop)
    pass, ev = diag_type22(fam, ["s - 17*mm"]; rng=rng)
    @test !pass && ev.spurious == ["s - 17*mm"]
end

if "--full" in ARGS
    @testset "sunrise gate" begin
        rng = MersenneTwister(20260630)
        fam = family("x1*x2 + x1*x3 + x2*x3 + s*x1*x2*x3 - mm*(x1+x2+x3)*(x1*x2+x1*x3+x2*x3)",
                     ["s", "mm"], ["x1", "x2", "x3"])
        r = landau(fam; rng=rng)
        @test r.chi_sectors[[1, 2, 3]] == 4
        @test Set(r.genuine) == Set(["s", "mm", "s - mm", "s - 9*mm"])
        @test isempty(r.spurious)
        @test diag_type21(fam; rng=rng)[1]
    end
end

println("ALL TESTS PASSED")

# test_staircase.jl — battery for the interval-sliced staircase counter.
# Engine-free: includes ONLY src/staircase.jl (pure Julia; Oscar is touched by
# lead_exps/chi_of_ideal at call time only, never here), so this leg runs on a bare
# Julia with no Oscar install. Run: julia test_staircase.jl
# Oracle: an independent brute-force box walker, defined here and never used by the
# shipped counter (only viable for small boxes).
include(joinpath(@__DIR__, "src", "staircase.jl"))
using Random, Test

"Brute-force oracle: walk the pure-power bounding box, count non-divisible points."
function brute_count(les::Vector{Vector{Int}}, nvars::Int)
    isempty(les) && return nothing
    any(e -> all(==(0), e), les) && return 0
    bounds = fill(-1, nvars)
    for e in les
        supp = findall(!=(0), e)
        if length(supp) == 1
            j = supp[1]
            bounds[j] = bounds[j] < 0 ? e[j] : min(bounds[j], e[j])
        end
    end
    any(<(0), bounds) && return nothing
    count = 0
    for idx in Iterators.product((0:b-1 for b in bounds)...)
        v = collect(idx)
        any(e -> all(v .>= e), les) || (count += 1)
    end
    count
end

t0 = time()

@testset "laws (unit / positive-dimensional / 1-var)" begin
    @test staircase_count(Vector{Vector{Int}}(), 2) === nothing
    @test staircase_count([[0, 0]], 2) == 0
    @test staircase_count([[0, 0], [3, 1]], 2) == 0
    @test staircase_count([[2, 0], [1, 1]], 2) === nothing   # no pure power in y
    @test staircase_count([[10^9, 0, 0], [0, 10^9, 0], [10, 10, 10]], 3) === nothing
    @test staircase_count([[5]], 1) == 5
    @test staircase_count([[5], [3], [7]], 1) == 3
end

@testset "known small staircases" begin
    @test staircase_count([[2, 0], [0, 3]], 2) == 6
    @test staircase_count([[2, 0], [1, 1], [0, 3]], 2) == 4
    @test staircase_count([[3, 0, 0], [0, 3, 0], [0, 0, 3]], 3) == 27
    @test staircase_count([[2, 0, 0], [0, 2, 0], [0, 0, 2], [1, 1, 1]], 3) == 7
end

@testset "randomized cross-check vs brute-force oracle" begin
    rng = MersenneTwister(20260901)
    for trial in 1:300
        nv = rand(rng, 1:4)
        les = Vector{Vector{Int}}()
        for j in 1:nv                        # pure powers (drop one 10% of the time)
            rand(rng) < 0.1 && nv > 1 && continue
            e = zeros(Int, nv); e[j] = rand(rng, 1:6)
            push!(les, e)
        end
        for _ in 1:rand(rng, 0:6)            # mixed generators
            e = [rand(rng, 0:5) for _ in 1:nv]
            all(==(0), e) && continue
            push!(les, e)
        end
        isempty(les) && continue
        shuffle!(rng, les)                   # order invariance comes free
        @test staircase_count(les, nv) == brute_count(les, nv)
    end
end

@testset "huge boxes, exact closed forms" begin
    a = 10^6
    @test staircase_count([[a, 0], [0, a]], 2) == 10^12                    # box 1e12
    @test staircase_count([[10^9, 0], [1, 1], [0, 10^9]], 2) == 2 * 10^9 - 1
    c = 10^3                                                               # corner cut
    @test staircase_count([[a, 0], [c, c], [0, a]], 2) == a^2 - (a - c)^2
    # three thresholds per axis, checked against the same shape at brute-force scale
    big = [[40_000, 0], [30_000, 10_000], [10_000, 20_000], [0, 30_000]]
    sml = [[4, 0], [3, 1], [1, 2], [0, 3]]
    @test staircase_count(big, 2) == 10_000^2 * brute_count(sml, 2)
    # count wider than Int64 comes back exact as Int128
    w = staircase_count([[10^7, 0, 0], [0, 10^7, 0], [0, 0, 10^7]], 3)
    @test w == Int128(10)^21 && w isa Int128
end

println("STAIRCASE PASS in ", round(time() - t0, digits=1), " s")
flush(stdout)

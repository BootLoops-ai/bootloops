# Tests for src/symmetry.jl (port of upstream Symmetry / SymmetryQuotient,
# scr.m:1107 / scr.m:1224).
#
# Self-contained: `src/SOFIA.jl` includes `symmetry.jl`; the guard below
# also injects the file into the module when this test is run against a
# SOFIA build without it, so the tests run unmodified either way.
#
# Run with:
#   julia --project=. -e 'using SOFIA, Test; include("test/symmetry_tests.jl")'

using Test
using SOFIA
import Nemo

if !isdefined(SOFIA, :symmetry_map)
    Base.include(SOFIA, joinpath(@__DIR__, "..", "src", "symmetry.jl"))
end

const smap = SOFIA.symmetry_map
const squot = SOFIA.symmetry_quotient

"Generator of `parent(p)` with name `nm` (for checking rules by name)."
function genbyname(p, nm)
    R = Nemo.parent(p)
    i = findfirst(==(Symbol(nm)), collect(Nemo.symbols(R)))
    i === nothing && error("no generator named $nm in $(Nemo.symbols(R))")
    return Nemo.gen(R, i)
end

# shared fixtures -----------------------------------------------------------

# one-loop triangle, massless internal lines, massive legs M1 M2 M3
triA = Diagram([((1, 2), big(0)), ((2, 3), big(0)), ((1, 3), big(0))],
               [(1, :M1), (2, :M2), (3, :M3)])
# same triangle with external masses M1 <-> M2 swapped on the legs
triB = Diagram([((1, 2), big(0)), ((2, 3), big(0)), ((1, 3), big(0))],
               [(1, :M2), (2, :M1), (3, :M3)])
# one-loop bubbles with massless legs, distinct internal masses
bubA = Diagram([((1, 2), :m1), ((1, 2), :m2)], [(1, big(0)), (2, big(0))])
bubB = Diagram([((1, 2), :m1), ((1, 2), :m3)], [(1, big(0)), (2, big(0))])

@testset "symmetry_map: triangles related by leg relabeling" begin
    rules = smap(triA, triB)
    @test rules !== nothing
    # 3-point kinematics has no Mandelstam invariants, only masses
    @test sort(collect(keys(rules))) == ["MM1", "MM2", "MM3"]
    @test rules["MM1"] == genbyname(rules["MM1"], "MM2")
    @test rules["MM2"] == genbyname(rules["MM2"], "MM1")
    @test rules["MM3"] == genbyname(rules["MM3"], "MM3")
end

@testset "symmetry_map: diagram vs itself is the identity" begin
    for d in (triA, bubA)
        rules = smap(d, d)
        @test rules !== nothing
        for (k, v) in rules
            @test v == genbyname(v, k)
        end
    end
end

@testset "symmetry_map: bubbles differing in one internal mass" begin
    rules = smap(bubA, bubB)
    @test rules !== nothing
    @test rules["mm1"] == genbyname(rules["mm1"], "mm1")
    @test rules["mm2"] == genbyname(rules["mm2"], "mm3")
    # 2-point kinematics: single invariant s maps to itself
    @test rules["s"] == genbyname(rules["s"], "s")
end

@testset "symmetry_map: different loop counts -> nothing" begin
    # 2-loop sunrise (3 internal lines) vs 1-loop bubble: edge counts differ
    sun = Diagram([((1, 2), :m1), ((1, 2), :m2), ((1, 2), :m3)],
                  [(1, big(0)), (2, big(0))])
    @test smap(bubA, sun) === nothing
    @test smap(sun, bubA) === nothing
    # same edge and node counts, different loop number (1 vs 0)
    chain = Diagram([((1, 2), :m1), ((2, 3), :m2)], [(1, big(0)), (3, big(0))])
    @test smap(bubA, chain) === nothing
end

@testset "symmetry_map: box vs cyclically rotated box" begin
    boxedges = [((1, 2), big(0)), ((2, 3), big(0)),
                ((3, 4), big(0)), ((1, 4), big(0))]
    boxA = Diagram(boxedges, [(1, big(0)), (2, big(0)), (3, big(0)), (4, big(0))])
    # legs rotated one step around the loop: leg j attaches at vertex j+1
    boxB = Diagram(boxedges, [(2, big(0)), (3, big(0)), (4, big(0)), (1, big(0))])
    rules = smap(boxA, boxB)
    @test rules !== nothing
    @test sort(collect(keys(rules))) == ["s12", "s23"]
    # the leg map sends (1,2) of A to (4,1) of B, so s12 -> (p4+p1)^2 = s23
    @test rules["s12"] == genbyname(rules["s12"], "s23")
    @test rules["s23"] == genbyname(rules["s23"], "s12")
end

@testset "symmetry_map: cross-bundle label consistency backtracking" begin
    # bundle (1,2) of A could pair m1 -> m3 locally, but bundle (2,3) forces
    # the global map m1 -> m4; the search must revise the first choice
    A = Diagram([((1, 2), :m1), ((1, 2), :m2), ((2, 3), :m1), ((1, 3), big(0))],
                [(1, big(0)), (2, big(0)), (3, big(0))])
    B = Diagram([((1, 2), :m3), ((1, 2), :m4), ((2, 3), :m4), ((1, 3), big(0))],
                [(1, big(0)), (2, big(0)), (3, big(0))])
    rules = smap(A, B)
    @test rules !== nothing
    @test rules["mm1"] == genbyname(rules["mm1"], "mm4")
    @test rules["mm2"] == genbyname(rules["mm2"], "mm3")
    # and a genuinely inconsistent labeling has no map (all quick rejects
    # pass): bundle (1,2) forces f(m1) = m3 = f(m2) while bundle (2,3)
    # forces f(m1) = m4
    C = Diagram([((1, 2), :m3), ((1, 2), :m3), ((2, 3), :m4), ((1, 3), big(0))],
                [(1, big(0)), (2, big(0)), (3, big(0))])
    @test smap(A, C) === nothing
end

@testset "symmetry_map transports singularities A -> B" begin
    # the defining property: substituting the rules (by generator name) into
    # the singularity polynomials of A reproduces those of B
    for (A, B) in ((bubA, bubB), (triA, triB))
        rules = smap(A, B)
        @test rules !== nothing
        polsA, RA = sofia_singularities(A; include_subtopologies=false)
        polsB, RB = sofia_singularities(B; include_subtopologies=false)
        Rmap = Nemo.parent(first(values(rules)))
        gd = Dict(string(s) => Nemo.gen(Rmap, i)
                  for (i, s) in enumerate(Nemo.symbols(Rmap)))
        mappedA = [Nemo.evaluate(p, Any[rules[string(s)] for s in Nemo.symbols(RA)])
                   for p in polsA]
        inB = [Nemo.evaluate(p, Any[gd[string(s)] for s in Nemo.symbols(RB)])
               for p in polsB]
        @test length(mappedA) == length(inB)
        @test all(any(isproportional(a, b) for b in inB) for a in mappedA)
    end
end

@testset "symmetry_quotient: greedy classes" begin
    tops = [triA, triB, bubA]
    reps, classes = squot(tops)
    @test length(reps) == 2
    @test reps[1] == triA
    @test reps[2] == bubA
    @test length(classes) == 2
    @test length(classes[1]) == 2     # identity for triA + map triA -> triB
    @test length(classes[2]) == 1     # identity for bubA
    # the first entry of each class is the representative's identity map
    for cls in classes
        for (k, v) in cls[1]
            @test v == genbyname(v, k)
        end
    end
    # the second entry of class 1 is the triA -> triB map
    @test classes[1][2]["MM1"] == genbyname(classes[1][2]["MM1"], "MM2")
end

# PLD bridge tests (src/pld.jl).
#
# The live backend (PLD.jl with its Oscar stack) is optional and heavy, so
# the bridge plumbing — Δ assembly, name-based polynomial transport into a
# coefficient ring + Laurent ring pair, the getSpecializedPAD call with
# upstream's option set, and the harvest back into a Nemo kinematic ring —
# is exercised against a mock backend built on the same AbstractAlgebra
# generic-ring interface Oscar implements. The live leg at the bottom runs
# whenever PLD.jl is installed, and otherwise records an explicit skip
# naming the missing backend.

module PLDMockBackend
import Nemo

# The Oscar surface the bridge touches, duck-typed on Nemo/AbstractAlgebra
# (deliberately the older constructor names, matching the Oscar versions
# the upstream bridge drives).
module MockRings
import Nemo
const QQ = Nemo.QQ
PolynomialRing(F, names::Vector{String}) = Nemo.polynomial_ring(F, names)
LaurentPolynomialRing(R, names::Vector{String}) = Nemo.polynomial_ring(R, names)
end

const Oscar = MockRings
const CALLS = Vector{Any}()

"""Mock solver: records the call, then returns the coefficients of Δ with
respect to the active/α variables — elements of the coefficient ring
QQ[kinematics], the same ring the real backend's specialized discriminants
live in. For a known Δ the output is exactly predictable, which is what
makes the round trip assertable."""
function getSpecializedPAD(delta, rgens, sgens; kwargs...)
    push!(CALLS, (delta=delta, rgens=rgens, sgens=sgens,
                  kwargs=Dict{Symbol,Any}(kwargs)))
    return [c for c in Nemo.coefficients(delta)]
end
end # module PLDMockBackend

@testset "PLD bridge: transplant" begin
    A, (p_, q_) = Nemo.polynomial_ring(Nemo.QQ, ["p", "q"])
    f = 3 // 2 * p_^3 * q_ - q_^2 + 7
    B, gb = Nemo.polynomial_ring(Nemo.QQ, ["w", "q", "p"])
    fB = SOFIA.transplant(f, Dict("p" => gb[3], "q" => gb[2]), one(B))
    @test !iszero(fB)
    # round trip: unused target variables need no image
    back = SOFIA.transplant(fB, Dict("p" => p_, "q" => q_), one(A))
    @test back == f
    # a variable that is USED but unmapped fails loudly
    @test_throws ErrorException SOFIA.transplant(f, Dict{String,Any}("p" => gb[3]), one(B))
    # zero maps to zero
    @test iszero(SOFIA.transplant(zero(A), Dict{String,Any}(), one(B)))
end

@testset "PLD bridge: delta assembly" begin
    R0, (x_, s_, t_, u_) = Nemo.polynomial_ring(Nemo.QQ, ["x", "s", "t", "u"])
    sys = SOFIA.pld_delta([x_ - u_, x_^2 - s_ * x_ + t_], [1])
    @test sys.kin == ["s", "t", "u"]
    @test sys.active == ["x"]
    @test sys.alphas == ["α0", "α1"]      # ALL α's, matching upstream allActive
    @test sys.ngrams == 2
    nm = Dict(string(sy) => Nemo.gen(sys.ring, k)
              for (k, sy) in enumerate(Nemo.symbols(sys.ring)))
    # Δ = α0·(x - u) + (x² - s·x + t): the LAST α is set to 1 (scrj.m:1306)
    expected = nm["α0"] * (nm["x"] - nm["u"]) +
               nm["x"]^2 - nm["s"] * nm["x"] + nm["t"]
    @test sys.delta == expected
    @test !SOFIA.free_of(sys.delta, findfirst(==("α0"), map(string, Nemo.symbols(sys.ring))))
    # no α1 in Δ, but α1 IS a ring variable
    @test SOFIA.free_of(sys.delta, findfirst(==("α1"), map(string, Nemo.symbols(sys.ring))))
end

@testset "PLD bridge: mock backend round trip" begin
    empty!(PLDMockBackend.CALLS)
    R0, (x_, s_, t_, u_) = Nemo.polynomial_ring(Nemo.QQ, ["x", "s", "t", "u"])
    sys = SOFIA.pld_delta([x_ - u_, x_^2 - s_ * x_ + t_], [1])
    discs, Rk = pld_discriminants(PLDMockBackend, sys)
    # Δ = x²·1 + x·α0·1 + x·(-s) + α0·(-u) + t over QQ[s,t,u]:
    # non-constant coefficients are exactly {-s, -u, t}
    @test map(string, Nemo.symbols(Rk)) == ["s", "t", "u"]
    @test length(discs) == 3
    nm = Dict(string(sy) => Nemo.gen(Rk, k)
              for (k, sy) in enumerate(Nemo.symbols(Rk)))
    for target in (nm["s"], nm["t"], nm["u"])
        @test any(p -> isproportional(p, target), discs)
    end
    # the backend saw the two-ring structure and the upstream option defaults
    @test length(PLDMockBackend.CALLS) == 1
    call = PLDMockBackend.CALLS[end]
    @test length(call.rgens) == 3          # QQ[s,t,u]
    @test length(call.sgens) == 3          # (x, α0, α1) — unused α1 included
    @test call.kwargs == Dict{Symbol,Any}(:method => :sym, :homogeneous => true,
                                          :high_prec => false, :codim_start => -1,
                                          :face_start => 1, :single_face => false)
end

@testset "PLD bridge: diagram lane (mock backend)" begin
    empty!(PLDMockBackend.CALLS)
    mz = BigInt(0)
    bub = Diagram([((1, 2), :m1), ((1, 2), :m2)], [(1, mz), (2, mz)])
    sys = pld_system(bub)
    @test all(a -> occursin(r"^x\d+$", a), sys.active)
    for nmreq in ("s", "mm1", "mm2")
        @test nmreq in sys.kin
    end
    sing, Rk = pld_singularities(PLDMockBackend, bub; method=:num, face_start=2)
    @test !isempty(sing)
    @test all(p -> !SOFIA.is_constant(p), sing)
    @test map(string, Nemo.symbols(Rk)) == sys.kin
    # explicit solver options pass through; the rest keep upstream defaults
    call = PLDMockBackend.CALLS[end]
    @test call.kwargs[:method] == :num
    @test call.kwargs[:face_start] == 2
    @test call.kwargs[:homogeneous] == true
end

@testset "PLD bridge: live backend" begin
    if Base.find_package("PLD") === nothing
        @info "SKIP PLD live leg: optional backend PLD.jl is not installed (mock legs above ran)"
        @test true
    else
        @eval Main import PLD
        mz = BigInt(0)
        bub = Diagram([((1, 2), :m1), ((1, 2), :m2)], [(1, mz), (2, mz)])
        sing, Rk = pld_singularities(Main.PLD, bub)
        @test !isempty(sing)
        @test all(p -> !SOFIA.is_constant(p), sing)
    end
end

using Test
using SOFIA
import Nemo

const REFDIR = joinpath(@__DIR__, "..", "reference")

@testset "WL parser" begin
    e = wlparse("m[1]^2 - 2*m[1]*m[2] + m[2]^2")
    @test e isa SOFIA.WLCall
    vars = SOFIA.wlvariables(e)
    @test length(vars) == 2

    # numeric edge cases
    @test wlparse("3") == BigInt(3)
    @test wlparse("-3") == BigInt(-3)
    @test wlparse("1/2") == 1 // big(2)
    # precedence: 2*s^2 + t
    e2 = wlparse("2*s^2 + t")
    @test e2.head == :Plus

    # right-associativity of ^
    e3 = wlparse("s^2^3")  # s^(2^3)
    R, (s,) = Nemo.polynomial_ring(Nemo.QQ, ["s"])
    @test SOFIA.tonemo(e3, Dict(:s => s)) == s^8
end

@testset "PLD_database loads (all files)" begin
    dbdir = joinpath(REFDIR, "PLD_database")
    files = filter(f -> endswith(f, ".m"), readdir(dbdir; join=true))
    @test length(files) >= 90
    nploys = 0
    for f in files
        polys, R, keys_ = loadsingularities(f)
        @test !isempty(polys)
        @test all(!iszero, polys)
        nploys += length(polys)
    end
    @info "parsed $(length(files)) database files, $nploys singularity polynomials"
end

@testset "proportionality and factors" begin
    R, (x, y) = Nemo.polynomial_ring(Nemo.QQ, ["x", "y"])
    @test isproportional(2x + 4y, x + 2y)
    @test !isproportional(x + y, x - y)
    @test !isproportional(x + y, x + y + 1)
    d = dedup_proportional([x + y, 3x + 3y, x - y])
    @test length(d) == 2
    fl = factor_list_unique([x^2 - y^2, (x + y)^2])
    @test length(fl) == 2   # x+y and x-y, up to constants
end

@testset "fubini: one-loop bubble (analytic check)" begin
    # Maximal-cut Baikov polynomial of the massive bubble in one variable x:
    # the Cayley/Gram polynomial  P(x) = x^2 - 2(m1+m2-s)... — we use the
    # standard form: P(x) = ((s - m1 - m2)^2 - 4 m1 m2) ... Simplest honest
    # test: eliminate x from the quadratic  x^2 + b x + c  with symbolic
    # coefficients chosen as the bubble's: b = -(s - m1 - m2), c = m1*m2.
    # Discriminant: (s - m1 - m2)^2 - 4 m1 m2 = (s - (sqrt m1 + sqrt m2)^2)(s - (sqrt m1 - sqrt m2)^2)
    # in squared-mass variables: the classic Landau curve (Kallen function).
    R, (x, s, m1, m2) = Nemo.polynomial_ring(Nemo.QQ, ["x", "s", "m1", "m2"])
    P = x^2 - (s - m1 - m2) * x + m1 * m2
    sing = fubini([P], [1])
    # expect exactly the Kallen polynomial lambda(s,m1,m2), irreducible over Q
    kallen = s^2 + m1^2 + m2^2 - 2s * m1 - 2s * m2 - 2m1 * m2
    @test any(p -> isproportional(p, kallen), sing)
    @test all(p -> SOFIA.free_of(p, 1), sing)
end

@testset "fastfubini: two-variable system" begin
    # Eliminate x,y from {x*y - t, x + y - s}: resultant chain gives the
    # discriminant of z^2 - s z + t -> s^2 - 4 t (and leading coeffs).
    R, (x, y, s, t) = Nemo.polynomial_ring(Nemo.QQ, ["x", "y", "s", "t"])
    polys = [x * y - t, x + y - s]
    sing = fastfubini(polys, [1, 2])
    target = s^2 - 4t
    @test any(p -> isproportional(p, target), sing)
    @test all(p -> SOFIA.free_of(p, 1) && SOFIA.free_of(p, 2), sing)
end

@testset "kinematics" begin
    for n in 2:7
        k = generate_kinematics(n)
        # momentum conservation: sum_j p_i·p_j == 0 for every i
        for i in 1:n
            @test iszero(sum(SOFIA.dot(k, i, j) for j in 1:n))
        end
        # diagonal entries are the squared masses
        if n >= 3
            for i in 1:n
                @test SOFIA.dot(k, i, i) == k.masses[i]
            end
        end
    end
    # 4-point sanity: (p1+p2)^2 == s12 and (p3+p4)^2 == s12 (Mandelstam s)
    k4 = generate_kinematics(4)
    R = k4.ring
    s12 = Nemo.gens(R)[findfirst(==("s12"), k4.basis)]
    lhs = SOFIA.dot(k4, 1, 1) + 2 * SOFIA.dot(k4, 1, 2) + SOFIA.dot(k4, 2, 2)
    rhs = SOFIA.dot(k4, 3, 3) + 2 * SOFIA.dot(k4, 3, 4) + SOFIA.dot(k4, 4, 4)
    @test lhs == s12
    @test rhs == s12
    # basis count = n(n-3)/2 + n masses
    for n in 4:7
        k = generate_kinematics(n)
        @test length(k.basis) == div(n * (n - 3), 2)
    end
end

@testset "loop edges" begin
    mz = BigInt(0)
    box = Diagram([((1, 2), mz), ((2, 3), mz), ((3, 4), mz), ((4, 1), mz)],
                  [(1, mz), (2, mz), (3, mz), (4, mz)])
    le = fix_loop_edges(box)
    @test length(le) == 1
    banana = Diagram([((1, 2), mz), ((1, 2), mz), ((1, 2), mz)],
                     [(1, mz), (2, mz)])
    @test length(fix_loop_edges(banana)) == 2
    # double box: 7 internal edges, 2 loops
    dbox = Diagram([((1, 2), mz), ((2, 3), mz), ((3, 6), mz), ((6, 5), mz),
                    ((5, 4), mz), ((4, 1), mz), ((2, 5), mz)],
                   [(1, mz), (3, mz), (4, mz), (6, mz)])
    @test nloops(dbox) == 2
    @test length(fix_loop_edges(dbox)) == 2
end

@testset "end-to-end: massive bubble" begin
    m(i) = SOFIA.WLCall(:m, Any[BigInt(i)])
    bub = Diagram([((1, 2), m(1)), ((1, 2), m(2))],
                  [(1, BigInt(0)), (2, BigInt(0))])
    sing, R = sofia_singularities(bub; include_subtopologies=false)
    nm = Dict(string(s) => Nemo.gen(R, k) for (k, s) in enumerate(Nemo.symbols(R)))
    s, mm1, mm2 = nm["s"], nm["mm1"], nm["mm2"]
    kallen = s^2 + mm1^2 + mm2^2 - 2s * mm1 - 2s * mm2 - 2mm1 * mm2
    @test any(p -> isproportional(p, kallen), sing)
end

@testset "end-to-end: massless box vs intuition" begin
    mz = BigInt(0)
    box = Diagram([((1, 2), mz), ((2, 3), mz), ((3, 4), mz), ((4, 1), mz)],
                  [(1, mz), (2, mz), (3, mz), (4, mz)])
    sing, R = sofia_singularities(box)
    nms = Dict(string(s) => Nemo.gen(R, k) for (k, s) in enumerate(Nemo.symbols(R)))
    s12, s23 = nms["s12"], nms["s23"]
    # the massless box's Landau singularities: s, t, and s+t (u-channel)
    @test any(p -> isproportional(p, s12), sing)
    @test any(p -> isproportional(p, s23), sing)
    @test any(p -> isproportional(p, s12 + s23), sing)
end

@testset "effortless odd letters: bubble" begin
    R, (s, mm1, mm2) = Nemo.polynomial_ring(Nemo.QQ, ["s", "mm1", "mm2"])
    kallen = s^2 + mm1^2 + mm2^2 - 2s * mm1 - 2s * mm2 - 2mm1 * mm2
    alphabet = [s, mm1, mm2, kallen]
    odd = find_odd_letters(kallen, alphabet)
    # the three classic bubble odd letters: P in {s-mm1-mm2, s+mm1-mm2, s-mm1+mm2}
    for P in (s - mm1 - mm2, s + mm1 - mm2, s - mm1 + mm2)
        @test any(l -> isproportional(l.P, P), odd)
    end
    @test length(odd) == 3
    # their dlogs satisfy one relation with the even letters: exactly 2 survive
    indep = independent_letters(odd, alphabet)
    @test length(indep) == 2
    # full pipeline finds letters for the kallen root
    all_ = effortless_odd_letters(alphabet)
    @test count(l -> isproportional(l.Q, kallen), all_) == 2
end

include("symmetry_tests.jl")
include("pld_tests.jl")

@testset "diagrams and subtopologies" begin
    mz = BigInt(0)
    m(i) = SOFIA.WLCall(:m, Any[BigInt(i)])
    M(i) = SOFIA.WLCall(:M, Any[BigInt(i)])
    # one-loop box: vertices 1..4, four external legs
    box = Diagram([((1, 2), mz), ((2, 3), mz), ((3, 4), mz), ((4, 1), mz)],
                  [(1, M(1)), (2, M(2)), (3, M(3)), (4, M(4))])
    @test nloops(box) == 1
    subs = subtopologies(box)
    # box -> triangles (4) -> bubbles (3 pairings appear as contractions) ...
    @test box in subs
    @test length(subs) > 4
    @test all(one_vertex_irreducible, subs)
    # sunset/banana: two vertices, three edges
    banana = Diagram([((1, 2), m(1)), ((1, 2), m(2)), ((1, 2), m(3))],
                     [(1, M(1)), (2, M(1))])
    @test nloops(banana) == 2
    @test one_vertex_irreducible(banana)
    # tadpole deletion
    tad = Diagram([((1, 1), m(1)), ((1, 2), m(2))], [(1, M(1)), (2, M(2))])
    @test nloops(delete_tadpoles(tad)) == 0
end

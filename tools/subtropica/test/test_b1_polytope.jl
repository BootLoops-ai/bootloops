# test_b1_polytope.jl — C4 polytope-layer tests (unit-P; DESIGN_B1 §2).
#
# Standalone (DESIGN_B1 §3): DEFAULT global Julia env (has Oscar/Nemo/JSON),
# include-by-path — run as
#   ulimit -v 32505856; julia tools/subtropica/test/test_b1_polytope.jl
# Under the shared suite (runtests.jl auto-include) the subtropica project env
# has NO Oscar until integration: the guard below then records a LOUD skip
# (@warn + @test_skip), never a silent pass (DESIGN_B1 §3 consequence rule).
#
# Spec: SPEC_B1.md §2 (C4), §6.2, gates G1/G7 (+ G8 vertex-set/incidence
# round-trip). Fixtures: SPEC_B1 §8.1/§8.3 hand tables; paper eq 3.19 toy
# control (paper:1017-1036) and the external ray-census fixtures
# ({raycheck_results.json, lbl3se_uf.json}; not shipped — the census testsets
# skip loudly without them).
#
# Structure note: sources are include()d and helpers defined in guarded
# TOP-LEVEL expressions, and each @testset is its own guarded top-level
# expression — so every testset runs in a world that already contains the
# included methods (no define-and-call inside one top-level expression).

using Test
using Nemo, JSON

const _SUBTROPICA = joinpath(@__DIR__, "..")
# Raycheck census fixtures are NOT vendored. Env unset ⇒ the dependent
# cross-check testsets record LOUD skips (same pattern as the Oscar guard
# below — @warn + @test_skip, never a silent pass); env SET but files
# missing ⇒ still a hard FAIL (a misconfigured fixture copy must not skip).
const _RAYCHECK = get(ENV, "SUBTROPICA_RAYCHECK_DIR", "")  # fixture root (not vendored); set env to a fixture copy
const _B1P_RAYCHECK_OK = !isempty(_RAYCHECK)
if !_B1P_RAYCHECK_OK
    @warn "SUBTROPICA_RAYCHECK_DIR unset — raycheck census fixtures not " *
          "vendored; the raycheck cross-check testsets are SKIPPED LOUDLY " *
          "(set the env var to a fixture copy to arm them)"
    @testset "raycheck census cross-checks — FIXTURES MISSING" begin
        @test_skip "SUBTROPICA_RAYCHECK_DIR unset — raycheck-dependent testsets skipped (loud skip)"
    end
end

# ---- LOUD Oscar guard (never a silent pass) --------------------------------
const _B1P_OSCAR_OK = try
    # load Oscar into THIS module — the shared suite
    # runs each test file in its own sandbox module (runtests.jl), where
    # `@eval Main using Oscar` would leave newton_polytope & co. invisible
    # here and trip _c4_require_oscar. Standalone: @__MODULE__ === Main.
    Core.eval(@__MODULE__, :(using Oscar))
    true
catch err
    @warn "test_b1_polytope: OSCAR NOT AVAILABLE in this env — C4 polytope " *
          "tests SKIPPED LOUDLY (expected under the Phase-A subtropica project " *
          "env; Oscar enters at integration, DESIGN_B1 §4 step 1)" err
    false
end

if !_B1P_OSCAR_OK
    @testset "C4 polytope layer (unit-P) — OSCAR MISSING" begin
        @test_skip "Oscar unavailable — every C4 test skipped (loud skip)"
    end
end

# ---- sources + test-local helpers (guarded; definitions only, no calls) ----
if _B1P_OSCAR_OK
    include(joinpath(_SUBTROPICA, "src", "types.jl"))
    include(joinpath(_SUBTROPICA, "src", "b1_types.jl"))
    include(joinpath(_SUBTROPICA, "src", "polytope.jl"))

    # test-local exact TropI (classification cross-checks only):
    # TropI(ρ) = Σ ν_i ρ_i + Σ_j e_j · max_{m ∈ supp P_j}(m·ρ) (SPEC_B1 §3.1,
    # paper eq 3.12/3.14; frozen EpsExp arithmetic from b1_types.jl). unit-T
    # owns the production version (src/tropical.jl); this local copy exists
    # only so C4 ray output can be checked against the raycheck artifacts.
    _supp(p) = [BigInt.(ev) for ev in Nemo.exponent_vectors(p)]
    _tropmax(sup, rho) =
        maximum(sum(m[i] * rho[i] for i in eachindex(rho)) for m in sup)
    function _tropI(nu, pes, rho)
        acc = EpsExp(0 // 1, 0 // 1)
        for (nui, ri) in zip(nu, rho)
            acc += nui * ri
        end
        for (sup, e) in pes
            acc += e * _tropmax(sup, rho)
        end
        return acc
    end
    _rat(s::AbstractString) = occursin("/", s) ?
        (x -> parse(BigInt, x[1]) // parse(BigInt, x[2]))(split(s, "/")) :
        parse(BigInt, s) // 1

    # supporting-hyperplane consistency: vertex_sets[i] must equal
    # {v : rays[i]·v == facet_rhs[i]} — ties argmax (wl:10259-10264) to the
    # (independent) facet offsets. Part of gate G8 round-trips.
    function _incidence_consistent(td)
        for (i, r) in enumerate(td.rays)
            onface = sort!(findall(v -> sum(v[k] * r[k] for k in eachindex(r)) ==
                                        td.facet_rhs[i], td.vertices))
            onface == td.vertex_sets[i] || return false
        end
        return true
    end

    # coeff-free poly from an integer support list (JSON or literal)
    function _from_support(R, sup)
        B = MPolyBuildCtx(R)
        for m in sup
            push_term!(B, one(Nemo.QQ), Int.(m))
        end
        finish(B)
    end
end

# ---------------------------------------------------------------------------
_B1P_OSCAR_OK && @testset "eq 3.19 toy control (paper:1017-1036; raycheck)" begin
    R, (a1, a2) = polynomial_ring(Nemo.QQ, [:a1, :a2])
    P = a1^2 + a2 + a1 * a2                       # eq 3.19 polynomial factor
    clear_tropical_data_cache!()
    td = tropical_data([P], [:a1, :a2])

    @test td.ambient_dim == 2
    @test isempty(td.equations)                    # full-dimensional triangle
    # outer normals, primitive, lex-sorted (paper:1030: ρ1=(−1,−2), ρ2=(1,1),
    # ρ3=(0,1); Oscar rows are ALREADY outer — no flip, SPEC_B1 §2.2 pin):
    @test td.rays == [[-1, -2], [0, 1], [1, 1]]
    @test td.vertices == [[0, 1], [1, 1], [2, 0]]  # lattice triangle, lex
    # per-ray argmax vertex sets (hand): (−1,−2)→{(0,1),(2,0)},
    # (0,1)→{(0,1),(1,1)}, (1,1)→{(1,1),(2,0)}:
    @test td.vertex_sets == [[1, 3], [1, 2], [2, 3]]
    @test td.facet_rhs == [-2 // 1, 1 // 1, 2 // 1]
    @test _incidence_consistent(td)

    # full ε-exact TropI on every ray (eq 3.20: I = a1^{2+ε} a2^{1+ε} P^{−(2+ε)}):
    nu = [EpsExp(2 // 1, 1 // 1), EpsExp(1 // 1, 1 // 1)]
    pes = [(_supp(P), EpsExp(-2 // 1, -1 // 1))]
    tI = Dict(r => _tropI(nu, pes, r) for r in td.rays)
    @test tI[[-1, -2]] == EpsExp(0 // 1, -1 // 1)  # TropI(ρ1) = −ε: log-div
                                                   # (paper:1031 prints +ε —
                                                   # display sign; a = 0 and
                                                   # |b| = 1 either way)
    @test tI[[0, 1]] == EpsExp(-1 // 1, 0 // 1)    # convergent
    @test tI[[1, 1]] == EpsExp(-1 // 1, 0 // 1)    # convergent
    @test classify_ray(tI[[-1, -2]]) == log_divergent
    @test classify_ray(tI[[0, 1]]) == convergent
    @test classify_ray(tI[[1, 1]]) == convergent

    # raycheck artifact cross-check (ray set + a-part + class; the artifact's
    # b_eps deliberately tracked only the POLY exponents' ε — its nu was the
    # a-part (2,1), as the fixture records it — so b is compared via `tI` above):
    if _B1P_RAYCHECK_OK
        rc = JSON.parsefile(joinpath(_RAYCHECK, "raycheck_results.json"))["control"]
        rc_rays = sort([[parse(BigInt, x) for x in r["rho"]] for r in rc["rays"]])
        @test rc_rays == td.rays
        for r in rc["rays"]
            rho = [parse(BigInt, x) for x in r["rho"]]
            @test tI[rho].a == _rat(r["a"])
            @test (r["class"] == "log") == (classify_ray(tI[rho]) == log_divergent)
        end
    else
        @test_skip "SUBTROPICA_RAYCHECK_DIR unset — control-census artifact cross-check skipped (loud skip; hand-worked eq 3.19 assertions above still ran)"
    end
end

# ---------------------------------------------------------------------------
_B1P_OSCAR_OK && @testset "F1/F3 fixture geometry (SPEC_B1 §8.1/§8.3 — G1)" begin
    R, (x1, x2) = polynomial_ring(Nemo.QQ, [:x1, :x2])

    # F1: 1 + x1 + x2 — triangle (0,0),(1,0),(0,1) (§8.1 table)
    td1 = tropical_data([1 + x1 + x2], [:x1, :x2])
    @test td1.rays == [[-1, 0], [0, -1], [1, 1]]
    @test td1.vertices == [[0, 0], [0, 1], [1, 0]]
    @test td1.vertex_sets == [[1, 2], [1, 3], [2, 3]]  # §8.1: ρ1↦{(0,0),(0,1)},
                                                       # ρ2↦{(0,0),(1,0)}, ρ3↦rest
    @test td1.facet_rhs == [0 // 1, 0 // 1, 1 // 1]
    @test isempty(td1.equations)
    @test _incidence_consistent(td1)

    # F3: 1 + x2 + x1·x2 — triangle (0,0),(0,1),(1,1) (§8.3 table; rays
    # A=(−1,0) log, B=(0,1) conv, C=(1,−1) log there)
    td3 = tropical_data([1 + x2 + x1 * x2], [:x1, :x2])
    @test td3.rays == [[-1, 0], [0, 1], [1, -1]]
    @test td3.vertices == [[0, 0], [0, 1], [1, 1]]
    @test td3.vertex_sets == [[1, 2], [2, 3], [1, 3]]  # §8.3 vertex_set column
    @test td3.facet_rhs == [0 // 1, 1 // 1, 0 // 1]
    @test _incidence_consistent(td3)
end

# ---------------------------------------------------------------------------
_B1P_OSCAR_OK && @testset "Minkowski multi-factor == expanded product (§2.1)" begin
    R, (x1, x2) = polynomial_ring(Nemo.QQ, [:x1, :x2])
    f1 = 1 + x1 + x2
    f2 = 1 + x1
    prod12 = f1 * f2            # = 1 + 2x1 + x2 + x1^2 + x1x2 (coeff 2 — the
                                #   coeff-free canonicalization must kill it)
    tdF = tropical_data([f1, f2], [:x1, :x2]; use_cache = false)
    tdP = tropical_data([prod12], [:x1, :x2]; use_cache = false)
    @test tdF == tdP            # frozen fieldwise == (b1_types.jl)
    @test length(tdF.rays) == 4 # quadrilateral: hull{(0,0),(2,0),(0,1),(1,1)}
    @test tdF.rays == [[-1, 0], [0, -1], [0, 1], [1, 1]]
    @test tdF.vertices == [[0, 0], [0, 1], [1, 1], [2, 0]]
    @test tdF.vertex_sets == [[1, 2], [1, 4], [2, 3], [3, 4]]
    @test tdF.facet_rhs == [0 // 1, 0 // 1, 1 // 1, 2 // 1]
    @test _incidence_consistent(tdF)
    # refinement ⇒ same cache key (G7; keyed on the coeff-free product,
    # wl:10419 + SPEC_B1 §2.4 pin):
    @test tropical_data_cache_key([f1, f2], [:x1, :x2]) ==
          tropical_data_cache_key([prod12], [:x1, :x2])
end

# ---------------------------------------------------------------------------
_B1P_OSCAR_OK && @testset "1-var analytic branch (wl:9697-9709; §2.3a — G1)" begin
    R1, (x,) = polynomial_ring(Nemo.QQ, [:x])

    # segment poly 1 + x + x^3 (G1's named case)
    td = tropical_data([1 + x + x^3], [:x])
    @test td.rays == [[-1], [1]]                 # hard-coded, wl:9701
    @test td.vertices == [[0], [1], [3]]         # ALL support points — the
                                                 # interior 1 kept (wl:9704-9705)
    @test td.vertex_sets == [[1], [3]]           # argmax = endpoints
    @test td.facet_rhs == [0 // 1, 3 // 1]       # canonical (−min, max); the
                                                 # upstream {min,+1} offset
                                                 # quirk is D8 — not replicated
    @test isempty(td.equations)
    @test td.ambient_dim == 1
    @test _incidence_consistent(td)

    # multi-factor 1-var: product path (1+x)(1+x^2) → support {0,1,2,3}
    tdm = tropical_data([1 + x, 1 + x^2], [:x])
    @test tdm.vertices == [[0], [1], [2], [3]]
    @test tdm.vertex_sets == [[1], [4]]
    @test tdm.facet_rhs == [0 // 1, 3 // 1]

    # single-monomial 1-var poly: degenerate (scaleless) ⇒ nonempty equations
    # (wl:9703) ⇒ :DegeneratePolytope refusal (SPEC_B1 §2.3a/§6.2)
    tds = tropical_data([x^2], [:x])
    @test tds.rays == [[-1], [1]]                # parity: upstream still emits
                                                 # the ray pair when degenerate
    @test tds.equations == [([1 // 1], 2 // 1)]  # canonical a·x = b (x = 2)
    @test_throws B1Refusal refuse_if_degenerate(tds)
    err = try; refuse_if_degenerate(tds); nothing; catch e; e; end
    @test err isa B1Refusal && err.kind == :DegeneratePolytope
    @test !isempty(err.detail)                   # detail ALWAYS populated (§6)
    @test err.detail["dim"] == 0 && err.detail["ambient_dim"] == 1

    # constant-only factor list: DeleteCases[...,1] (wl:10415) leaves nothing
    # ⇒ point at 0 ⇒ degenerate, same refusal
    tdc = tropical_data([R1(7)], [:x])
    @test !isempty(tdc.equations)
    @test_throws B1Refusal refuse_if_degenerate(tdc)
end

# ---------------------------------------------------------------------------
_B1P_OSCAR_OK && @testset "degenerate 2d segment refusal (§6.2; G6(iv) geom)" begin
    R, (x1, x2) = polynomial_ring(Nemo.QQ, [:x1, :x2])
    td = tropical_data([1 + x1 * x2], [:x1, :x2])   # segment (0,0)–(1,1) in 2d
    @test !isempty(td.equations)                    # dim 1 < n = 2
    @test length(td.equations) == 1
    a, b = td.equations[1]
    for v in td.vertices                            # every vertex satisfies a·x = b
        @test sum(a[i] * v[i] for i in 1:2) == b    # (hull is x1 − x2 = 0 up
    end                                             #  to scale/sign)
    err = try; refuse_if_degenerate(td); nothing; catch e; e; end
    @test err isa B1Refusal && err.kind == :DegeneratePolytope
    @test err.detail["dim"] == 1 && err.detail["ambient_dim"] == 2
end

# ---------------------------------------------------------------------------
_B1P_OSCAR_OK && _B1P_RAYCHECK_OK && @testset "lbl3se census cross-check (20 facets, all conv)" begin
    uf = JSON.parsefile(joinpath(_RAYCHECK, "lbl3se_uf.json"))
    rc = JSON.parsefile(joinpath(_RAYCHECK, "raycheck_results.json"))["lbl3se"]
    vars7 = Symbol.(uf["vars"])                  # x0..x6 (Cheng–Wu x7 = 1)
    R7 = polynomial_ring(Nemo.QQ, vars7)[1]

    U = _from_support(R7, uf["U_support"])       # 32 monomials
    F = _from_support(R7, uf["F_support_phys"])  # 127 monomials, physical point
    @test length(U) == 32 && length(F) == 127

    td = tropical_data([U, F], vars7)
    @test td.ambient_dim == 7
    @test isempty(td.equations)
    @test length(td.rays) == 20                  # census headline: 20 facets
    @test length(td.vertices) == rc["n_vertices"]          # 68
    rc_rays = sort([[parse(BigInt, x) for x in r["rho"]] for r in rc["rays"]])
    @test td.rays == rc_rays                     # exact primitive-ray match
    @test _incidence_consistent(td)

    # classification parity: ν_i = 1 (dlog form), U^{4ε}, F^{−2−3ε} —
    # all 20 rays convergent with the artifact's exact (a, b_eps):
    nu = fill(EpsExp(1 // 1, 0 // 1), 7)
    pes = [(_supp(U), EpsExp(0 // 1, 4 // 1)),
           (_supp(F), EpsExp(-2 // 1, -3 // 1))]
    byray = Dict([parse(BigInt, x) for x in r["rho"]] => r for r in rc["rays"])
    for rho in td.rays
        t = _tropI(nu, pes, rho)
        @test t.a == _rat(byray[rho]["a"])
        @test t.b == _rat(byray[rho]["b_eps"])
        @test classify_ray(t) == convergent
    end
end

# ---------------------------------------------------------------------------
_B1P_OSCAR_OK && _B1P_RAYCHECK_OK && @testset "eq 4.3 eikonal-web census (3-factor Minkowski, 5v)" begin
    # The B1 divergent-path benchmark: 17 facets =
    # 13 convergent + 4 log, ν = 1, factors Q_N^{+1} U^{−1+3ε} G^{−2−2ε}
    # (supports as in the external ray-census fixture; gauge x5 = 1 → vars x1,x2,x3,x4,x6).
    rc = JSON.parsefile(joinpath(_RAYCHECK, "raycheck_results.json"))["eq4.3"]
    vars5 = [:x1, :x2, :x3, :x4, :x6]
    R5 = polynomial_ring(Nemo.QQ, vars5)[1]
    drop5(m) = [m[1], m[2], m[3], m[4], m[6]]
    QN6 = [[0,1,0,1,0,0],[0,0,1,1,0,0],[0,1,0,0,0,1],[0,0,1,0,1,0]]
    U6  = [[1,1,0,0,0,0],[1,0,1,0,0,0],[0,1,1,0,0,0]]
    G6  = [[1,1,0,1,0,0],[1,0,1,1,0,0],[0,1,1,1,0,0],
           [1,1,0,0,1,0],[1,0,1,0,1,0],[0,1,1,0,1,0],
           [0,1,0,2,0,0],[0,0,1,2,0,0],[0,1,0,1,0,1],
           [0,0,1,1,1,0],[1,0,0,0,2,0],[0,0,1,0,2,0],
           [1,0,0,0,1,1]]
    QN = _from_support(R5, unique([drop5(v) for v in QN6]))
    U  = _from_support(R5, unique([drop5(v) for v in U6]))
    G  = _from_support(R5, unique([drop5(v) for v in G6]))

    td = tropical_data([QN, U, G], vars5)
    @test length(td.rays) == 17
    @test length(td.vertices) == rc["n_vertices"]          # 37
    @test isempty(td.equations)
    rc_rays = sort([[parse(BigInt, x) for x in r["rho"]] for r in rc["rays"]])
    @test td.rays == rc_rays
    @test _incidence_consistent(td)

    nu = fill(EpsExp(1 // 1, 0 // 1), 5)
    pes = [(_supp(QN), EpsExp(1 // 1, 0 // 1)),
           (_supp(U), EpsExp(-1 // 1, 3 // 1)),
           (_supp(G), EpsExp(-2 // 1, -2 // 1))]
    byray = Dict([parse(BigInt, x) for x in r["rho"]] => r for r in rc["rays"])
    nlog = 0
    for rho in td.rays
        t = _tropI(nu, pes, rho)
        @test t.a == _rat(byray[rho]["a"])
        @test t.b == _rat(byray[rho]["b_eps"])
        cls = classify_ray(t)
        @test (cls == log_divergent) == (byray[rho]["class"] == "log")
        cls == log_divergent && (nlog += 1)
    end
    @test nlog == 4                              # the ε^{-4} tower's 4 log rays
end

# ---------------------------------------------------------------------------
_B1P_OSCAR_OK && @testset "vertex_sets argmax == Oscar incidence (G8)" begin
    # SPEC_B1 §2.2: the argmax rule (wl:10259-10264) MUST reproduce Oscar's
    # own facet-vertex incidence. Rebuild the control + quadrilateral fresh
    # through Oscar's first-class API and compare, mapping Oscar's emission
    # order onto our sorted arrays by (primitive ray, vertex) content.
    R, (a1, a2) = polynomial_ring(Nemo.QQ, [:a1, :a2])
    for polys in ([a1^2 + a2 + a1 * a2], [1 + a1 + a2, 1 + a1])
        td = tropical_data(polys, [:a1, :a2]; use_cache = false)
        P = reduce(Oscar.minkowski_sum,
                   (Oscar.newton_polytope(forget_coeffs(f)) for f in polys))
        A, b = Oscar.halfspace_matrix_pair(Oscar.facets(P))
        inc = Oscar.facets(Oscar.IncidenceMatrix, P)
        Vosc = collect(Oscar.vertices(P))
        # vertex map: Oscar order → our sorted order
        vmap = [findfirst(==([BigInt(numerator(x)) for x in v]), td.vertices)
                for v in Vosc]
        @test all(!isnothing, vmap)
        for i in 1:Nemo.nrows(A)
            r, _ = primitive_integer([A[i, j] for j in 1:2])
            k = findfirst(==(r), td.rays)
            @test !isnothing(k)
            oscar_set = sort!([vmap[j] for j in Oscar.row(inc, i)])
            @test oscar_set == td.vertex_sets[k]
        end
    end
end

# ---------------------------------------------------------------------------
_B1P_OSCAR_OK && @testset "cache keying (SPEC_B1 §2.4 — G7)" begin
    R, (x1, x2) = polynomial_ring(Nemo.QQ, [:x1, :x2])
    clear_tropical_data_cache!()

    p_cf   = 1 + x1 + x2
    p_coef = Nemo.QQ(5) + 3 * x1 + Nemo.QQ(7 // 2) * x2   # same support,
                                                          # different coeffs
    k_cf   = tropical_data_cache_key([p_cf], [:x1, :x2])
    k_coef = tropical_data_cache_key([p_coef], [:x1, :x2])
    @test k_cf == k_coef                        # coeff values NEVER in the key

    td1 = tropical_data([p_cf], [:x1, :x2])
    @test haskey(TROPICAL_DATA_CACHE, k_coef)
    td2 = tropical_data([p_coef], [:x1, :x2])
    @test td2 === td1                           # HIT: identical stored object

    # different support ⇒ MISS
    k_other = tropical_data_cache_key([1 + x1 + x1 * x2], [:x1, :x2])
    @test k_other != k_cf

    # permuted vars ⇒ MISS (same abstract poly, different coordinate frame)
    Rp, (y2, y1) = polynomial_ring(Nemo.QQ, [:x2, :x1])
    k_perm = tropical_data_cache_key([1 + y1 + y2], [:x2, :x1])
    @test k_perm != k_cf

    # factor refinement ⇒ HIT (also asserted in the Minkowski testset)
    @test tropical_data_cache_key([p_cf, 1 + x1], [:x1, :x2]) ==
          tropical_data_cache_key([p_cf * (1 + x1)], [:x1, :x2])

    # constant factors are dropped from the key (wl:10415)
    @test tropical_data_cache_key([p_cf, R(42)], [:x1, :x2]) == k_cf

    # input-contract refusals (SPEC_B1 §2.1)
    @test_throws ArgumentError tropical_data([p_cf], [:x1])       # wrong vars
    @test_throws ArgumentError tropical_data([p_cf], [:x2, :x1])  # wrong order
    @test_throws ArgumentError tropical_data([zero(R)], [:x1, :x2])
    @test_throws ArgumentError tropical_data(QQMPolyRingElem[], [:x1, :x2])
end

# ---------------------------------------------------------------------------
_B1P_OSCAR_OK && @testset "normal_fan_diagnostic: flag-gated (§2.2, D7)" begin
    R, (x1, x2) = polynomial_ring(Nemo.QQ, [:x1, :x2])
    @test_throws ErrorException normal_fan_diagnostic([1 + x1 + x2], [:x1, :x2])
    nf = normal_fan_diagnostic([1 + x1 + x2], [:x1, :x2]; enable = true)
    # diagnostic-only sanity: one maximal cone per vertex of the triangle
    @test Oscar.n_maximal_cones(nf) == 3
end

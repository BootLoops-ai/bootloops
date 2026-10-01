# test_b1_tropical.jl — OWNER: unit-T (DESIGN_B1 §2). Tests for
# src/tropical.jl (C3+C5): TropI evaluation, ray classification, Σ_div,
# w-search, u-functions, initial-form restriction, §6 refusal ladder.
#
# Self-contained by construction: everything lives in its own module
# (TestB1Tropical), with types.jl / b1_types.jl / tropical.jl included BY
# PATH — so it runs BOTH
#   * standalone under the DEFAULT global Julia env (no Oscar needed —
#     Nemo + JSON only):
#       ulimit -v 32505856; julia --startup-file=no \
#         tools/subtropica/test/test_b1_tropical.jl
#   * under the shared suite (test/runtests.jl auto-include): the module
#     wall prevents any clash with `using SubTropica` exports.
#
# Coverage (task spec): SPEC_B1 §8 micro-fixture ray/GP data (hand-worked);
# eq 4.3 eikonal census cross-check vs the external ray-census fixtures (not shipped)
# (13 conv / 4 log / 0 power, GP holds ×4); lbl3se convergent census;
# paper control eq 3.19/3.20; synthetic GP-violation (paper eq 3.37) that
# must throw typed; all §6 refusal kinds owned by tropical.jl (1/3/4/5/6);
# ordering pins (G2); restriction order-independence + w→u→w round-trip
# (G8, C5 slice); mutation controls so the census comparator provably bites
# (synthetic-truth controls are mandatory).

module TestB1Tropical

using Test
using Nemo
import JSON

const _SRC = normpath(joinpath(@__DIR__, "..", "src"))
include(joinpath(_SRC, "types.jl"))
include(joinpath(_SRC, "b1_types.jl"))
include(joinpath(_SRC, "tropical.jl"))

# Raycheck census fixtures are NOT vendored. Env unset ⇒ the census
# cross-check testsets record LOUD skips (@warn + @test_skip — the same
# consequence rule as the suite's missing-Oscar/engine guards, never a
# silent pass); env SET but files missing ⇒ still a hard FAIL (the
# cross-checks are gate feed on an armed fixture copy).
const RAYCHECK_DIR = get(ENV, "SUBTROPICA_RAYCHECK_DIR", "")  # fixture root (not vendored); set env to a fixture copy
const RAYCHECK_OK = !isempty(RAYCHECK_DIR)
RAYCHECK_OK || @warn "SUBTROPICA_RAYCHECK_DIR unset — raycheck census " *
    "fixtures not vendored; the census cross-check testsets are SKIPPED " *
    "LOUDLY (set the env var to a fixture copy to arm them)"
const NOEQ = Tuple{Vector{Rational{BigInt}},Rational{BigInt}}[]

# exact-rational parse of raycheck JSON number strings ("-2", "3/2", …)
parseq(s::AbstractString) = begin
    if occursin('/', s)
        p = split(s, '/')
        parse(BigInt, p[1]) // parse(BigInt, p[2])
    else
        parse(BigInt, s) // 1
    end
end
parseq(x::Integer) = BigInt(x) // 1

bigvec(v) = BigInt[BigInt(x) for x in v]

mkE(nu, polys, vars, kin, ring) = EulerIntegrand(Prefactor(), nu, polys, vars, kin, ring)

mktd(rays, verts, vsets, rhs, eqs, n) = TropicalData(
    [bigvec(r) for r in rays], [bigvec(v) for v in verts],
    [Int.(s) for s in vsets], [Rational{BigInt}(q) for q in rhs], eqs, n)

# w → u → w round-trip: the .wl step-4 W-vector recovery (wl:10770-10775):
# A = numerator(1−u), B = denominator(1−u) − A, w = exp(B) − exp(A).
# The production recovery is C7's (unit-S, round-trip unit test per SPEC_B1
# §4 step 4); this test-local copy pins the u-convention from THIS side.
function recover_w(xs, w)
    A, den = one_minus_u_pair(xs, w)
    B = den - A
    @assert length(A) == 1 && length(B) == 1   # monomials by construction
    ea = exponent_vector(A, 1)
    eb = exponent_vector(B, 1)
    return [eb[i] - ea[i] for i in eachindex(xs)]
end

# shared 2-var ring + the F1/F2 triangle geometry (SPEC_B1 §8.1/§8.2:
# Newt = triangle (0,0),(1,0),(0,1); rays and argmax vertex sets hand-worked)
R2, R2x = polynomial_ring(QQ, [:x1, :x2])
x1, x2 = R2x
td_F1 = mktd([[-1,0],[0,-1],[1,1]], [[0,0],[1,0],[0,1]],
             [[1,3],[1,2],[2,3]], [0,0,1], NOEQ, 2)

@testset "unit-T tropical.jl (C3+C5)" begin

# ---------------------------------------------------------------------------
@testset "paper control eq 3.19/3.20 (STtropicalize∘STEvalRay semantics)" begin
    Rc, (a1, a2) = polynomial_ring(QQ, [:a1, :a2])
    Ec = mkE([EpsExp(2,1), EpsExp(1,1)], [(a1^2 + a2 + a1*a2, EpsExp(-2,-1))],
             [:a1, :a2], Symbol[], Rc)
    # Paper (below eq 3.21): rho1 = (−1,−2) log-divergent, rho2 = (1,1) and
    # rho3 = (0,1) convergent. NOTE the paper TEXT prints Trop I(rho1) = +ε;
    # plugging rho1 into the paper's own eq 3.20 (and into the .wl formula)
    # gives −ε — the a-part (classification) is identical; code wins
    # (SPEC_B1 preamble). This test pins the code value.
    t1 = trop_on_ray(Ec, [-1, -2])
    @test t1 == EpsExp(0, -1)
    @test classify_ray(t1) == log_divergent
    @test trop_on_ray(Ec, [1, 1]) == EpsExp(-1, 0)
    @test trop_on_ray(Ec, [0, 1]) == EpsExp(-1, 0)
    # a-part parity vs the raycheck control census. Only `a` is comparable:
    # raycheck.jl:105 passed the ε-FREE part of ν (classification-only
    # convention), so its control b_eps (+2 on rho1) differs from the true
    # integrand's b (−1) by the dropped ν ε-parts.
    if RAYCHECK_OK
        ctrl = JSON.parsefile(joinpath(RAYCHECK_DIR, "raycheck_results.json"))["control"]
        for r in ctrl["rays"]
            rho = [parse(BigInt, x) for x in r["rho"]]
            @test trop_on_ray(Ec, rho).a == parseq(r["a"])
        end
    else
        @test_skip "SUBTROPICA_RAYCHECK_DIR unset — control-census a-part parity skipped (loud skip; paper-pinned trop values above still ran)"
    end
end

# ---------------------------------------------------------------------------
@testset "F1 fixture (SPEC_B1 §8.1 — paper worked example, post-NP)" begin
    E = mkE([EpsExp(0,1), EpsExp(0,1)], [(1 + x1 + x2, EpsExp(-1,-3))],
            [:x1, :x2], Symbol[], R2)
    trv = trop_values(E, td_F1)
    @test trv == [EpsExp(0,-1), EpsExp(0,-1), EpsExp(-1,-1)]      # §8.1 table
    @test [classify_ray(t) for t in trv] ==
          [log_divergent, log_divergent, convergent]
    dd = divergence_data(E, td_F1)
    @test dd.div_facets == [1, 2]
    @test [f.facets for f in dd.sigma_div] == [Int[], [1], [2], [1,2]]
    @test dd.sigma_div[1].common_vertices == [1,2,3]   # empty-face convention
    @test dd.sigma_div[2].common_vertices == [1,3]
    @test dd.sigma_div[3].common_vertices == [1,2]
    @test dd.sigma_div[4].common_vertices == [1]       # common vertex v1 (§8.1)
    @test dd.w[1] == [1, 0] && dd.w[2] == [0, 1]       # §8.1 w-search, w·rho = −1
    @test u_pair([x1,x2], dd.w[1]) == (x1, 1 + x1)     # u1 = x1/(1+x1)
    @test u_pair([x1,x2], dd.w[2]) == (x2, 1 + x2)
    @test dd.facets ==
          [DivergentFacet(1, log_divergent, EpsExp(0,-1), bigvec([1,0])),
           DivergentFacet(2, log_divergent, EpsExp(0,-1), bigvec([0,1]))]
    # restrictions (§8.1): along rho1: 1+x1+x2 → 1+x2; rho2: → 1+x1; both: → 1
    P = 1 + x1 + x2
    @test restrict_poly(P, [-1, 0], 2) == 1 + x2
    @test restrict_poly(P, [0, -1], 2) == 1 + x1
    @test restrict_poly(P, [[-1,0], [0,-1]], 2) == one(R2)
    @test restrict_poly(P, [[0,-1], [-1,0]], 2) == one(R2)  # order-independent (G8)
    # ORIGINAL coefficients kept (wl:10510)
    @test restrict_poly(2 + 3*x1 + 5*x2, [-1, 0], 2) == 2 + 5*x2
    # w → u → w round-trip (G8 slice; production recovery is C7's step 4)
    @test recover_w([x1,x2], dd.w[1]) == [1, 0]
    @test recover_w([x1,x2], dd.w[2]) == [0, 1]
end

# ---------------------------------------------------------------------------
@testset "F2 fixture (SPEC_B1 §8.2 — single facet, simplest-w branch)" begin
    E = mkE([EpsExp(1,1), EpsExp(0,1)], [(1 + x1 + x2, EpsExp(-2,-3))],
            [:x1, :x2], Symbol[], R2)
    trv = trop_values(E, td_F1)          # same triangle geometry as F1
    @test trv == [EpsExp(-1,-1), EpsExp(0,-1), EpsExp(-1,-1)]     # §8.2 table
    dd = divergence_data(E, td_F1)
    @test dd.div_facets == [2]
    @test [f.facets for f in dd.sigma_div] == [Int[], [2]]
    # compDiv(2) = {} ⇒ simplest-w (wl:10603-10607): firstNonZero of (0,−1)
    # is slot 2 ⇒ w = (0,1), w·rho = −1 (D2 pin satisfied naturally)
    @test dd.w[2] == [0, 1]
    @test u_pair([x1,x2], dd.w[2]) == (x2, 1 + x2)
end

# ---------------------------------------------------------------------------
@testset "F3 fixture (SPEC_B1 §8.3 — non-axis w, twist geometry)" begin
    td = mktd([[-1,0], [0,1], [1,-1]], [[0,0], [0,1], [1,1]],
              [[1,2], [2,3], [1,3]], [0, 1, 0], NOEQ, 2)
    E = mkE([EpsExp(0,1), EpsExp(0,2)], [(1 + x2 + x1*x2, EpsExp(-1,-5))],
            [:x1, :x2], Symbol[], R2)
    trv = trop_values(E, td)
    @test trv == [EpsExp(0,-1), EpsExp(-1,-3), EpsExp(0,-1)]      # §8.3 (A,B,C)
    dd = divergence_data(E, td)
    @test dd.div_facets == [1, 3]
    @test [f.facets for f in dd.sigma_div] == [Int[], [1], [3], [1,3]]
    @test dd.sigma_div[4].common_vertices == [1]        # common vertex v1
    @test dd.w[1] == [1, 1]   # NON-AXIS w_A: nullspace of rho_C = (1,−1); w·rho_A = −1
    @test dd.w[3] == [0, 1]   # w_C: nullspace of rho_A = (−1,0); w·rho_C = −1
    @test u_pair([x1,x2], dd.w[1]) == (x1*x2, 1 + x1*x2)    # u_A = x1x2/(1+x1x2)
    @test u_pair([x1,x2], dd.w[3]) == (x2, 1 + x2)          # u_C = x2/(1+x2)
    @test one_minus_u_pair([x1,x2], dd.w[1]) == (one(R2), 1 + x1*x2)
    # restrictions (§8.3): along A: 1+x2+x1x2 → 1+x2; along C: → 1+x1x2; both: 1
    P = 1 + x2 + x1*x2
    @test restrict_poly(P, [-1, 0], 2) == 1 + x2
    @test restrict_poly(P, [1, -1], 2) == 1 + x1*x2
    @test restrict_poly(P, [[-1,0], [1,-1]], 2) == one(R2)
    @test restrict_poly(P, [[1,-1], [-1,0]], 2) == one(R2)  # order-independent (G8)
    @test recover_w([x1,x2], dd.w[1]) == [1, 1]
end

# ---------------------------------------------------------------------------
@testset "ordering pins (G2) + u-function Factor form" begin
    # Mathematica Subsets order (wl:10267-10268) — hand-listed Subsets[Range[4],{2}]
    @test lex_combinations(collect(1:4), 2) ==
          [[1,2], [1,3], [1,4], [2,3], [2,4], [3,4]]
    @test lex_combinations([9,10,12,15], 3) ==
          [[9,10,12], [9,10,15], [9,12,15], [10,12,15]]
    @test lex_combinations(Int[], 0) == [Int[]]
    @test isempty(lex_combinations([1, 2], 3))
    @test isempty(sigma_div(td_F1, Int[]))     # wl:10256 parity (empty sel ⇒ {})
    @test_throws ArgumentError sigma_div(td_F1, [2, 1])   # sel must ascend
    # negative-w Factor form (wl:10636): w = (−1,1) ⇒ u = x2/(x1+x2)
    @test u_pair([x1,x2], [-1, 1]) == (x2, x1 + x2)
    @test one_minus_u_pair([x1,x2], [-1, 1]) == (x1, x1 + x2)
    @test recover_w([x1,x2], [-1, 1]) == [-1, 1]
end

# ---------------------------------------------------------------------------
@testset "kinematic variables are tropically inert (SPEC_B1 §2.1/§3.1)" begin
    Rk, (kx1, kx2, ks) = polynomial_ring(QQ, [:x1, :x2, :s])
    Pk = ks + kx1 + ks^2*kx2           # kinvar-dressed 1 + x1 + x2
    @test trop_poly(Pk, [-1, 0], 2) == 0
    @test restrict_poly(Pk, [-1, 0], 2) == ks + ks^2*kx2  # x-part selection only,
                                                          # kinvar coefficients kept
    Ek = mkE([EpsExp(0,1), EpsExp(0,1)], [(Pk, EpsExp(-1,-3))],
             [:x1, :x2], [:s], Rk)
    @test trop_values(Ek, td_F1) == [EpsExp(0,-1), EpsExp(0,-1), EpsExp(-1,-1)]
end

# ---------------------------------------------------------------------------
@testset "typed refusals (SPEC_B1 §6 — unit-T kinds 1/3/4/5/6 + ladder order)" begin
    # ---- §6.1 TropIllDefined --------------------------------------------------
    @test_throws B1Refusal classify_ray(EpsExp(0//1, 0//1))
    E0 = mkE([EpsExp(0), EpsExp(0,1)], [(1 + x1 + x2, EpsExp(-1,-1))],
             [:x1, :x2], Symbol[], R2)
    err = try divergence_data(E0, td_F1); nothing catch e; e end
    @test err isa B1Refusal && err.kind == :TropIllDefined
    @test err.detail["ray_index"] == 1
    # ladder order: ill-defined (ray 1) fires BEFORE power (ray 3)
    Epow0 = mkE([EpsExp(0), EpsExp(2)], [(1 + x1 + x2, EpsExp(0,-1))],
                [:x1, :x2], Symbol[], R2)
    err2 = try divergence_data(Epow0, td_F1); nothing catch e; e end
    @test err2 isa B1Refusal && err2.kind == :TropIllDefined

    # ---- §6.2 DegeneratePolytope (authoritative construction test = unit-P
    #      G6-iv; here: the C5 entry refuses on nonempty equations) -------------
    tdseg = mktd([[-1,-1], [1,1]], [[0,0], [1,1]], [[1], [2]], [0, 2],
                 [(Rational{BigInt}[1//1, -1//1], 0//1)], 2)   # 1 + x1·x2 segment
    Eseg = mkE([EpsExp(1), EpsExp(1)], [(1 + x1*x2, EpsExp(-1,-1))],
               [:x1, :x2], Symbol[], R2)
    err3 = try divergence_data(Eseg, tdseg); nothing catch e; e end
    @test err3 isa B1Refusal && err3.kind == :DegeneratePolytope
    @test err3.detail["dim"] == 1 && err3.detail["ambient_dim"] == 2

    # ---- §6.3 PowerDivergentRay (G6-ii): x1^{1+ε} x2^ε (1+x1+x2)^{−ε} ---------
    Epow = mkE([EpsExp(1,1), EpsExp(0,1)], [(1 + x1 + x2, EpsExp(0,-1))],
               [:x1, :x2], Symbol[], R2)
    errp = try divergence_data(Epow, td_F1); nothing catch e; e end
    @test errp isa PowerDivergentRefusal
    @test errp.ray == 3 && errp.tropI == EpsExp(1, 1)
    @test errp.detail["rays"] == [3]
    @test errp.detail["orders"] == [1//1]
    @test errp.detail["facets"][1].w === missing    # transient record, B2 payload
    @test occursin("Nilsson-Passare", sprint(showerror, errp))  # message law
    @test occursin("allow_continuation", sprint(showerror, errp))  # B2 route named

    # ---- §6.4 GeometricPropertyViolated: paper eq 3.37 RAW integrand (G6-i) ---
    # x1^ε x2^ε (1+x1+x2)^{−3ε}: ALL THREE rays log-divergent, every pair
    # compatible, no valid w exists (paper:1315-1325: all three v_rho NotFound).
    E37 = mkE([EpsExp(0,1), EpsExp(0,1)], [(1 + x1 + x2, EpsExp(0,-3))],
              [:x1, :x2], Symbol[], R2)
    trv37 = trop_values(E37, td_F1)
    @test trv37 == [EpsExp(0,-1), EpsExp(0,-1), EpsExp(0,-1)]
    sig37 = sigma_div(td_F1, [1, 2, 3])
    # size-lex order; the triple has NO common vertex ⇒ absent
    @test [f.facets for f in sig37] ==
          [Int[], [1], [2], [3], [1,2], [1,3], [2,3]]
    errg = try produce_ws(td_F1, sig37); nothing catch e; e end
    @test errg isa GeometricPropertyViolated
    @test errg.facet == 1 && errg.compatible == [2, 3]
    # trivial nullspace case: rank{rho2,rho3} = 2 = n ⇒ empty candidate list
    @test haskey(errg.detail, "nullspace") && isempty(errg.detail["nullspace"])
    @test haskey(errg.detail, "hint")
    @test occursin("Nilsson-Passare", sprint(showerror, errg))  # message law
    @test occursin("allow_continuation", sprint(showerror, errg))  # B2 route named
    @test_throws GeometricPropertyViolated divergence_data(E37, td_F1)

    # ---- §6.5 WNotUnitNormalized (D2 pin), simplest-w branch ------------------
    # primitive ray (2,−1), no compatible partner ⇒ w = −2·e1, w·rho = −4 ≠ −1
    # (upstream with normalizeQ=False, wl:10758, would silently proceed)
    tdw = mktd([[2,-1]], [[0,0]], [[1]], [0], NOEQ, 2)
    sigw = [SigmaFace(Int[], [1]), SigmaFace([1], [1])]
    errw = try produce_ws(tdw, sigw); nothing catch e; e end
    @test errw isa WNotUnitNormalized
    @test errw.facet == 1 && errw.w == [-2, 0] && errw.rho == [2, -1]
    @test errw.product == -4

    # ---- §6.5 WNotUnitNormalized, nullspace branch ----------------------------
    # rho_f = (1,−2), compatible ray (1,1): kernel candidate (−1,1), dot = −3
    tdw2 = mktd([[1,-2], [1,1]], [[0,0], [1,1]], [[1], [1]], [0, 2], NOEQ, 2)
    sigw2 = [SigmaFace(Int[], [1]), SigmaFace([1], [1]), SigmaFace([2], [1]),
             SigmaFace([1,2], [1])]
    errw2 = try produce_ws(tdw2, sigw2); nothing catch e; e end
    @test errw2 isa WNotUnitNormalized
    @test errw2.facet == 1 && errw2.w == [-1, 1] && errw2.product == -3

    # ---- §6.6 SigmaDivTooLarge -------------------------------------------------
    nbig = 25
    tdbig = mktd([[1, 0] for _ in 1:nbig], [[0, 0]], [[1] for _ in 1:nbig],
                 fill(0, nbig), NOEQ, 2)
    errs = try sigma_div(tdbig, collect(1:nbig)); nothing catch e; e end
    @test errs isa B1Refusal && errs.kind == :SigmaDivTooLarge
    @test errs.detail["count"] == 25
end

# ---------------------------------------------------------------------------
RAYCHECK_OK || @testset "eq 4.3 eikonal census cross-check — FIXTURES MISSING" begin
    @test_skip "SUBTROPICA_RAYCHECK_DIR unset — census cross-check skipped (loud skip)"
end
RAYCHECK_OK && @testset "eq 4.3 eikonal census cross-check (RAYCHECK parity — B1 gate target)" begin
    rcpath = joinpath(RAYCHECK_DIR, "raycheck_results.json")
    @test isfile(rcpath)    # hard fail, not a skip: with the env SET this
                            # cross-check is gate feed (missing file = broken
                            # fixture copy, never a skip)
    RC = JSON.parsefile(rcpath)
    e43 = RC["eq4.3"]
    @test e43["n_facets"] == 17 && e43["n_vertices"] == 37

    # Integrand transliterated from raycheck.jl:124-147 (Gardi–Zhu 2-loop
    # eikonal timelike web, paper eq 4.3; gauge x5 = 1, vars (x1,x2,x3,x4,x6)):
    #   ∏ x_i · Q_N^{+1} · U^{−1+3ε} · G^{−(2+2ε)}   (dlog measure)
    # Generic positive kinematics c_a = 5/4, y = 3 — support-exact (all G
    # coefficients positive; Q_N's mixed sign cannot cancel: distinct monomials).
    R5, xs5 = polynomial_ring(QQ, [:x1, :x2, :x3, :x4, :x6])
    X1, X2, X3, X4, X6 = xs5
    ca = QQ(5, 4)
    yy = QQ(3)
    QN = ca*(X2*X4 + X3*X4) - X2*X6 + X3
    UU = X1*X2 + X1*X3 + X2*X3
    GG = (X4 + 1)*UU + QQ(1, 4)*(X2*X4^2 + X3*X4^2 + 2*yy*X2*X4*X6 +
                                 2*ca*X3*X4 + X1 + X3 + 2*X1*X6)
    # supports must equal raycheck.jl:133-142 lists (drop5 applied) EXACTLY
    @test Set(support_x(QN, 5)) ==
          Set([[0,1,0,1,0], [0,0,1,1,0], [0,1,0,0,1], [0,0,1,0,0]])
    @test Set(support_x(UU, 5)) ==
          Set([[1,1,0,0,0], [1,0,1,0,0], [0,1,1,0,0]])
    @test Set(support_x(GG, 5)) ==
          Set([[1,1,0,1,0], [1,0,1,1,0], [0,1,1,1,0], [1,1,0,0,0], [1,0,1,0,0],
               [0,1,1,0,0], [0,1,0,2,0], [0,0,1,2,0], [0,1,0,1,1], [0,0,1,1,0],
               [1,0,0,0,0], [0,0,1,0,0], [1,0,0,0,1]])
    E43 = mkE(fill(EpsExp(1), 5),
              [(QN, EpsExp(1, 0)), (UU, EpsExp(-1, 3)), (GG, EpsExp(-2, -2))],
              [:x1, :x2, :x3, :x4, :x6], Symbol[], R5)

    rays43 = [bigvec([parse(BigInt, x) for x in r["rho"]]) for r in e43["rays"]]
    @test length(rays43) == 17
    trv = EpsExp[trop_on_ray(E43, rho) for rho in rays43]
    for (i, r) in enumerate(e43["rays"])       # exact per-ray parity (a AND b)
        @test trv[i].a == parseq(r["a"])
        @test trv[i].b == parseq(r["b_eps"])
    end
    cls = RayClass[classify_ray(t) for t in trv]
    @test count(==(convergent), cls) == 13     # RAYCHECK headline: 13/4/0
    @test count(==(log_divergent), cls) == 4
    @test count(==(power_divergent), cls) == 0
    dsel = div_facets(trv)
    gpf = sort(Int[g["facet"] for g in e43["gp"]])
    @test dsel == gpf
    for g in e43["gp"]                         # JSON internal consistency
        @test bigvec([parse(BigInt, x) for x in g["rho"]]) == rays43[g["facet"]]
    end
    # hand-pinned divergent table — independent of JSON parsing,
    # so the comparator cannot be vacuously true through a parse bug:
    for (rho, b) in [([-1,-1,-1,0,0], -4), ([1,1,1,0,1], 2),
                     ([-1,0,-1,-1,0], -1), ([0,0,1,0,1], 1)]
        t = trop_on_ray(E43, rho)
        @test t.a == 0 && t.b == b
    end
    # mutation control (must bite): ν₅ 1→2 pushes ray (0,0,1,0,1) to a = 1 > 0
    E43mut = mkE([fill(EpsExp(1), 4); EpsExp(2)], E43.polys, E43.vars, Symbol[], R5)
    tmut = trop_on_ray(E43mut, bigvec([0, 0, 1, 0, 1]))
    @test tmut.a > 0
    @test classify_ray(tmut) == power_divergent

    # Polytope data WITHOUT Oscar (C4 is unit-P's): Minkowski-sum support of
    # the coeff-free product + exact vertex rank test against the facet rows
    # from the raycheck JSON (v is a vertex ⇔ its tight facet normals span R^5).
    msum(A, B) = unique([a .+ b for a in A for b in B])
    sup = msum(msum(support_x(QN, 5), support_x(UU, 5)), support_x(GG, 5))
    rhs43 = [maximum(_xdot(v, rho, 5) for v in sup) for rho in rays43]
    verts = Vector{Int}[]
    for v in sup
        act = [i for i in 1:length(rays43) if _xdot(v, rays43[i], 5) == rhs43[i]]
        length(act) < 5 && continue
        M = matrix(QQ, length(act), 5, [QQ(rays43[i][j]) for i in act for j in 1:5])
        rank(M) == 5 && push!(verts, v)
    end
    @test length(verts) == 37                  # == Oscar's count (JSON n_vertices)
    vsets43 = [sort(findall(v -> _xdot(v, rays43[i], 5) == rhs43[i], verts))
               for i in 1:length(rays43)]
    td43 = TropicalData(rays43, [bigvec(v) for v in verts], vsets43,
                        [Rational{BigInt}(q) for q in rhs43], NOEQ, 5)

    dd = divergence_data(E43, td43)            # must NOT refuse: log-only + GP holds
    @test dd.div_facets == gpf
    # GP holds ×4, AND the D2 pin is satisfied: every selected w has w·rho == −1
    for f in dd.div_facets
        @test _dotbig(dd.w[f], rays43[f]) == -1
    end
    # mutual compatibility: all 6 pairs are faces; per-facet compatible count
    # == 3 == JSON n_compat_div
    prs = [fc.facets for fc in dd.sigma_div if length(fc.facets) == 2]
    @test length(prs) == 6
    for g in e43["gp"]
        f = g["facet"]
        @test count(p -> f in p, prs) == g["n_compat_div"] == 3
    end
    # the full 4-set is a face (common vertex ×4) ⇒ Vol carries 1/ε⁴ —
    # consistent with the published leading pole −1/(4ε⁴y) of eq 4.3
    # (RAYCHECK headline: "maximal compatible set size 4 ⇒ pole order ε⁻⁴")
    @test any(fc -> fc.facets == gpf, dd.sigma_div)
    @test maximum(length(fc.facets) for fc in dd.sigma_div) == 4
end

# ---------------------------------------------------------------------------
RAYCHECK_OK || @testset "lbl3se census cross-check — FIXTURES MISSING" begin
    @test_skip "SUBTROPICA_RAYCHECK_DIR unset — census cross-check skipped (loud skip)"
end
RAYCHECK_OK && @testset "lbl3se census cross-check (convergent path — 20/0/0)" begin
    ufpath = joinpath(RAYCHECK_DIR, "lbl3se_uf.json")
    rcpath = joinpath(RAYCHECK_DIR, "raycheck_results.json")
    @test isfile(ufpath) && isfile(rcpath)   # hard fail with env SET (see header)
    uf = JSON.parsefile(ufpath)
    se = JSON.parsefile(rcpath)["lbl3se"]
    n7 = length(uf["vars"])
    @test n7 == 7
    R7, xs7 = polynomial_ring(QQ, [Symbol("t", i) for i in 1:n7])
    mkpoly(sup) = sum(prod(xs7[j]^Int(m[j]) for j in 1:n7) for m in sup)
    Upoly = mkpoly(uf["U_support"])      # coeff-free: support is what matters,
    Fpoly = mkpoly(uf["F_support_phys"]) # and it is tight at the physical point
                                         # (all-positive, no loss)
    @test length(Upoly) == 32 && length(Fpoly) == 127   # ray-census fixture counts
    E_se = mkE(fill(EpsExp(1), n7),
               [(Upoly, EpsExp(0, 4)), (Fpoly, EpsExp(-2, -3))],  # U^{4ε} F^{−2−3ε}
               [Symbol("t", i) for i in 1:n7], Symbol[], R7)
    @test se["n_facets"] == 20 && se["convergent"] == 20
    trv = EpsExp[]
    for r in se["rays"]
        rho = bigvec([parse(BigInt, x) for x in r["rho"]])
        t = trop_on_ray(E_se, rho)
        @test t.a == parseq(r["a"])
        @test t.b == parseq(r["b_eps"])
        @test classify_ray(t) == convergent
        push!(trv, t)
    end
    # zero divergent rays: the tropical census certifies ε⁰-finiteness and the
    # subtraction machinery never engages (locally-finite route; C7 step 0)
    @test div_facets(trv) == Int[]
end

end # top testset

end # module TestB1Tropical

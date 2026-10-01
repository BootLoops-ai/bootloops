# SPDX-License-Identifier: GPL-3.0-or-later
#
# Derived from HyperInt (Erik Panzer, "Algorithms for the symbolic integration of
# hyperlogarithms with applications to Feynman integrals", Comput. Phys. Commun. 188
# (2015) 148, arXiv:1403.3385; https://bitbucket.org/PanzerErik/hyperint),
# Copyright (C) 2014 Erik Panzer, licensed under the GNU General Public License
# v3.0 or later. This Julia transliteration of HyperInt's compatibility-graph
# reduction procedures (with F. Brown, arXiv:0910.0114, as the mathematical source)
# is Copyright (c) 2026 Anthropic, PBC; created by Matthew D. Schwartz, code written
# by Claude (Anthropic) under his supervision, and is distributed under the same
# GNU General Public License v3.0 or later. The rest of this package is MIT-licensed;
# see ../NOTICE.
#
# test_lr_refine.jl — tests for the compatibility-graph LR refinement,
# src/lr_refine.jl.
#
# Self-contained BOTH ways (DESIGN_B1.md §3): everything is loaded inside a
# private module, so this file runs standalone in the DEFAULT global Julia
# env (Nemo + JSON, no Oscar needed) AND under the shared suite without
# namespace clashes:
#   ulimit -v 32505856; julia tools/subtropica/test/test_lr_refine.jl
#
# Fixtures:
#  (0) unit tests: canon / irreducibles / resultant / discriminant / parser.
#  (1) PUBLISHED-OUTPUT control — HyperInt paper (Panzer arXiv:1403.3385)
#      Example 4.2, S = {(1+x)^2+y, y+z^2}: reproduce L[{y}][1] =
#      {x+1, x+1+z, x+1-z} and L[{x,y}][1] = {1+z, z-1} EXACTLY, and the
#      non-computability of L[{x}] (S_∅ nonlinear in x).
#  (2) SPURIOUS-LETTER fixture (constructed from Brown arXiv:0910.0114
#      eq (71): resultants with four distinct grandparents are spurious):
#      naive Fubini keeps W = Res_y(Res_x(f1,f2), Res_x(f3,f4)); the
#      compatibility-graph reduction certifies W spurious.  Positive AND
#      negative controls included.
#  (3) recorded adjudication: blocking letter P = rho^2*tau^2-5*rho^2+3*tau^2+1
#      from the external probe fixtures, probe set fiber_constant_P (a recorded
#      NOLR certificate). Verdict JSON written to a temporary directory
#      (or SUBTROPICA_VERDICT_DIR when set).
#  (4) J2L cross-validation: the recorded analytic NOLR certificate
#      (probe set lbl3m_j2l of the external fixtures) must SURVIVE refinement —
#      forced x6 first step, then blocked in every remaining variable.
#  (5) typed refusals (algebraic letters etc.) and multi-group init.

module LRRefineTestEnv
using Nemo, JSON
include(joinpath(@__DIR__, "..", "src", "types.jl"))
include(joinpath(@__DIR__, "..", "src", "b1_types.jl"))
include(joinpath(@__DIR__, "..", "src", "lr_refine.jl"))
end # module

using Test, Nemo, JSON
const LR = LRRefineTestEnv

# Probe-bank fixtures are NOT vendored. Env unset ⇒ the two probe-bank
# adjudication testsets record LOUD skips (@warn + @test_skip, matching the
# suite's missing-Oscar/engine consequence rule — never a silent pass);
# env SET but files missing ⇒ still a hard FAIL (broken fixture copy).
const PROBES = get(ENV, "SUBTROPICA_PROBES_DIR", "")  # fixture root (not vendored); set env to a fixture copy
const PROBES_OK = !isempty(PROBES)
PROBES_OK || @warn "SUBTROPICA_PROBES_DIR unset — probe-bank fixtures not " *
    "vendored; the probe-fixture testsets are SKIPPED LOUDLY (set the env " *
    "var to a fixture copy to arm them)"
const VERDICT_DIR = get(ENV, "SUBTROPICA_VERDICT_DIR", mktempdir())

@testset "lr_refine" begin

# ---------------------------------------------------------------------------
@testset "unit: canon / irreducibles / resultant / discriminant" begin
    R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
    # canon: monic, constants/monomials -> nothing
    @test LR.canon_letter(QQ(5, 2) * (x + 1)) == x + 1
    @test LR.canon_letter(R(7)) === nothing
    @test LR.canon_letter(3 * x * y^2) === nothing     # monomial dropped
    @test LR.canon_letter(zero(R)) === nothing
    # irreducibles: factor + drop monomials (HyperInt.mpl:2614-2628)
    S = LR._irreducibles([x^2 * y - y])                # y*(x-1)*(x+1)
    @test S == Set([x - 1, x + 1])
    # bounded variant: only confirms members of `bound`
    Sb = LR._irreducibles([x^2 * y - y]; bound = Set([x - 1]))
    @test Sb == Set([x - 1])
    # resultant of two linears in x: [f,g]_x ~ g_x f_0 - g_0 f_x (Brown eq 66)
    r = LR._resultant_in(x + y, x + 2y + 1, 1)
    @test LR.canon_letter(r) == y + 1
    # resultant with an x-free polynomial: g^deg(f)
    @test LR._resultant_in(x + y, y + 1, 1) == y + 1
    # discriminant (Brown eq 65): D_y(y^2 + 2z y + 1) = 4z^2 - 4
    d = LR._discriminant_in(y^2 + 2 * z * y + 1, 2)
    @test LR._irreducibles([d]) == Set([z - 1, z + 1])
    # parser: **, rationals, refusals
    @test LR.parse_poly_qq("x**2 + (5/2)*y - 1", R) == x^2 + QQ(5, 2) * y - 1
    @test_throws LR.LRRefineRefusal LR.parse_poly_qq("sqrt(x)+1", R)
    @test_throws LR.LRRefineRefusal LR.parse_poly_qq("x + w", R)
    err = try LR.parse_poly_qq("sqrt(x)", R) catch e; e end
    @test err.kind === :AlgebraicLetter
end

# ---------------------------------------------------------------------------
# (1) HyperInt paper Example 4.2 (arXiv:1403.3385 §4.2) — published output:
#     S:={x^2+2*x+1+y, y+z^2}; L[{}]:=[S,{S}]; cgReduction(L):
#       L[{y}][1]   = {x+1, x+1+z, x+1-z}
#       L[{x,y}][1] = {1+z, z-1}
#       L[{x}]      not computed (S_∅ is not linear in x).
# ---------------------------------------------------------------------------
@testset "published control: HyperInt Example 4.2" begin
    R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
    polys = [(1 + x)^2 + y, y + z^2]
    vars = [:x, :y]                                    # z = kinematic

    cert, removed, log = LR.refine_letters(polys, vars, [:y, :x])
    @test log["order_walkable"] === true
    stages = log["stages"]
    @test stages[2]["computed"] === true
    @test Set(stages[2]["letters_refined"]) ==
          Set(string.([x + 1, x + 1 + z, x + 1 - z]))
    @test Set(stages[3]["letters_refined"]) ==
          Set(string.([z + 1, z - 1]))
    # no pruning happens in this example: naive == refined at every stage
    @test isempty(removed)

    # L[{x}] is NOT computable: order (x, y) is not walkable at max_degree 1
    _c2, _r2, log2 = LR.refine_letters(polys, vars, [:x, :y])
    @test log2["order_walkable"] === false
    @test log2["stages"][2]["computed"] === false

    # reachable stage sets: {}, {y}, {x,y} — exactly as in the paper
    c = LR.blocking_certificate(polys, vars, x + 1 + z)
    @test c["stages_reachable_refined"] == [String[], ["y"], ["x", "y"]]
    @test c["real"] === true                            # genuine letter survives
end

# ---------------------------------------------------------------------------
# (2) Spurious-letter fixture (constructed from Brown arXiv:0910.0114 eq (71):
#     [[f1,f2]_x, [f3,f4]_x]_y has FOUR distinct grandparents => spurious).
#
#     f_i = x + b_i(y,z), all monic-quadratic b_i in y  =>  x-first is the
#     unique admissible start (every f_i is quadratic in y), and ALL pairwise
#     x-resultants b_j - b_i are LINEAR in y, so the y step is admissible for
#     both bounds.  Hand-computed descendants (b_i coefficients checked below):
#       b1 = y(y+1),  b2 = (y+2)(y+z),  b3 = y(y+z),  b4 = (y+3)(y+5)
#       A  := Res_x(f1,f2) ~ b2-b1 = (z+1)y + 2z         (parents {f1,f2})
#       B  := Res_x(f3,f4) ~ b4-b3 = (8-z)y + 15         (parents {f3,f4})
#       W  := Res_y(A,B) = (8-z)(2z) - 15(z+1) = -2z^2 + z - 15
#     W is irreducible (disc 1-120 < 0), appears in the NAIVE stage-{x,y}
#     bound, and is EXCLUDED by the compatibility-graph rules (A,B share no
#     parent and no triangle) => certified spurious.
# ---------------------------------------------------------------------------
@testset "fixture: naive keeps spurious W, cg removes it" begin
    R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
    f1 = x + y^2 + y
    f2 = x + (y + 2) * (y + z)
    f3 = x + y^2 + z * y
    f4 = x + (y + 3) * (y + 5)
    polys = [f1, f2, f3, f4]
    vars = [:x, :y]
    W = LR.canon_letter(-2 * z^2 + z - 15)     # = z^2 - z/2 + 15/2

    cert, removed, log = LR.refine_letters(polys, vars, [:x, :y])
    @test log["order_walkable"] === true
    # forced start: y is inadmissible first (all f_i quadratic in y)
    @test log["stages"][1]["admissible_next_vars"] == ["x"]

    # the spurious letter: kept by naive, removed by cg
    @test W in removed
    @test !(W in cert)
    # positive controls: genuine letters survive refinement
    for ctrl in [z - 1, z + 1, z - 2, z - 3, z - 5,
                 LR.canon_letter((z + 1) * y + 2 * z),   # A itself (stage 1)
                 LR.canon_letter((8 - z) * y + 15)]      # B itself (stage 1)
        @test ctrl in cert
    end
    # refined ⊆ naive at every computed stage (upper-bound sanity)
    for stg in log["stages"]
        stg["computed"] || continue
        @test isempty(setdiff(Set(stg["letters_refined"]),
                              Set(stg["letters_naive"])))
    end

    # certificate + Bool API, francis (HyperInt default) and clique variants
    cW = LR.blocking_certificate(polys, vars, -2 * z^2 + z - 15)
    @test cW["status"] == "SPURIOUS_REMOVED"
    @test ["x", "y"] in cW["stages_naive"]      # negative control: W IS in
    @test isempty(cW["stages_refined"])          #   the naive bound (pruning
    @test LR.is_blocking_real(polys, vars, -2 * z^2 + z - 15) === false
    @test LR.is_blocking_real(polys, vars, -2 * z^2 + z - 15;
                              algorithm = :clique) === false
    @test LR.is_blocking_real(polys, vars, z - 1) === true
    @test LR.is_blocking_real(polys, vars, z - 1; algorithm = :clique) === true

    # string-letter entry point
    @test LR.is_blocking_real(polys, vars, "-2*z^2 + z - 15") === false

    # multi-group init (SubTropica.wl:15688: pairs per term): with groups
    # [[f1,f2],[f3,f4]], the CROSS-group resultant letter r24 ~ b4-b2 =
    # (6-z)y + (15-2z) is never formed — differential vs the flat group.
    r24 = LR.canon_letter((6 - z) * y + (15 - 2 * z))
    @test LR.is_blocking_real(polys, vars, (6 - z) * y + (15 - 2 * z)) === true
    cG = LR.blocking_certificate([[f1, f2], [f3, f4]], vars,
                                 (6 - z) * y + (15 - 2 * z))
    @test cG["real"] === false
    @test cG["status"] == "SPURIOUS_REMOVED"    # naive baseline has no groups
end

# ---------------------------------------------------------------------------
# (3) recorded probe: adjudicate P = rho^2*tau^2 - 5*rho^2 + 3*tau^2 + 1
#     Inputs = the recorded engine request (lr_req_sim.json: 9 polynomials,
#     xvars tau,rho,x3, no kinematic symbols; SUBSTITUTION_BANK.json pins
#     P's provenance as lc_x3(Q) at the FORCED first elimination x3).
# ---------------------------------------------------------------------------
PROBES_OK || @testset "recorded probe: blocking letter P — FIXTURES MISSING" begin
    @test_skip "SUBTROPICA_PROBES_DIR unset — probe adjudication skipped (loud skip)"
end
PROBES_OK && @testset "recorded probe: blocking letter P (fiber_constant_P)" begin
    reqf = joinpath(PROBES, "fiber_constant_P", "lr_req_sim.json")
    bankf = joinpath(PROBES, "fiber_constant_P", "SUBSTITUTION_BANK.json")
    @test isfile(reqf) && isfile(bankf)   # hard fail with env SET (see header)
    req = JSON.parsefile(reqf)
    xv = Symbol.(req["xvars"])                     # [:tau, :rho, :x3]
    R = polynomial_ring(QQ, xv)[1]
    polys = [LR.parse_poly_qq(s, R) for s in req["polys"]]
    Pstr = "rho^2*tau^2 - 5*rho^2 + 3*tau^2 + 1"
    P = LR.parse_poly_qq(Pstr, R)

    # provenance pin: P == lc_x3(Q), Q = the rationalized denominator letter
    Q = polys[2]
    i_x3 = findfirst(==(:x3), xv)
    @test LR.canon_letter(LR._coeff_in(Q, i_x3, 1)) == LR.canon_letter(P)

    certL, removed, log = LR.refine_letters(polys, xv, [:x3, :rho, :tau])
    # forced first step: x3 is the ONLY admissible pivot (bank: "FORCED
    # first elimination x3")
    @test log["stages"][1]["admissible_next_vars"] == ["x3"]
    # P survives the refined stage-{x3} bound (singleton lc descendant along
    # the unique projection: cg refinement CANNOT prune it)
    @test string(LR.canon_letter(P)) in log["stages"][2]["letters_refined"]
    @test !(LR.canon_letter(P) in removed)
    # and the reduction stays BLOCKED after x3 (P and 3*rho^2+1 quadratic)
    @test log["order_walkable"] === false
    @test log["stages"][2]["admissible_next_vars"] == String[]
    @test string(LR.canon_letter(P)) in log["stages"][2]["blocking_letters_next"]
    # bank cross-pins: Q(x3=0) factors (rho±1, 3tau^2+1) and Q(x3=1) factors
    # (tau±1) are in the refined bound
    tau, rho = gens(R)[1], gens(R)[2]
    for p in [rho - 1, rho + 1, tau - 1, tau + 1,
              LR.canon_letter(3 * tau^2 + 1), LR.canon_letter(3 * rho^2 + 1)]
        @test string(p) in log["stages"][2]["letters_refined"]
    end

    cP = LR.blocking_certificate(polys, xv, P)
    @test cP["status"] == "SURVIVES"
    @test ["x3"] in cP["stages_refined"]
    @test LR.is_blocking_real(polys, xv, P) === true

    # ---- deliverable: the verdict JSON ---------------------------------
    verdict = Dict{String,Any}(
        "target" => "fiber-constant integral (probe set fiber_constant_P)",
        "question" => "is the NOLR blocking letter P = rho^2*tau^2-5*rho^2+3*tau^2+1 real or spurious under compatibility-graph refinement?",
        "verdict" => "REAL (NOT certifiable spurious): P SURVIVES the compatibility-graph refinement — the recorded NOLR verdict STANDS.",
        "why" => string(
            "P is born as the LEADING COEFFICIENT lc_x3(Q) (a singleton descendant, ",
            "Brown 0910.0114 Def 73: singleton operations are never gated by ",
            "compatibility edges) along the FORCED first elimination x3 (the only ",
            "variable all nine letters are linear in), and the Fubini intersection ",
            "(Def 74) has a single parent there — no second reduction path exists ",
            "to intersect P away.  After x3, P is degree 2 in BOTH tau and rho, ",
            "so no further linear projection is admissible: blocked, refined and naive alike."),
        "caveat" => "survival in the refined bound is upper-bound survival, not a positive proof that P is a singularity of the partial integral; residual pruning routes OUTSIDE this tool's scope: (i) the bank's domain-positivity note (P > 0 on the open simplex — no interior pinch), (ii) rational changes of variables.",
        "certificate" => cP,
        "inputs" => Dict("request" => reqf, "bank" => bankf,
                         "polys" => req["polys"], "xvars" => req["xvars"]),
        "tool" => Dict(
            "file" => "tools/subtropica/src/lr_refine.jl",
            "algorithm" => "cgSingleReductionFrancis transliteration (HyperInt.mpl:2913-2986, default per HyperInt.mpl:79) + cgReduction subset DP (HyperInt.mpl:2989-3086); law: Brown arXiv:0910.0114 Def 68/73/74 + eq (71); Panzer arXiv:1403.3385 §4.2",
            "max_degree" => 1))
    open(joinpath(VERDICT_DIR, "P_VERDICT.json"), "w") do io
        JSON.print(io, verdict, 2)
    end
    @test isfile(joinpath(VERDICT_DIR, "P_VERDICT.json"))
end

# ---------------------------------------------------------------------------
# (4) J2L: the recorded analytic NOLR certificate must survive refinement.
#     Reconstruct the UNGAUGED 7-var U, F from the recorded pivot certificates
#     (F = F_A*x6 + F_B, U = U_A*x6 + U_B — cert_obstruction_structure.json),
#     tie them to the recorded gauge-x6 strings (step1_uf_strings.json), then
#     refine in the Cheng-Wu gauge x0 = 1:
#       * step 1 forced through x6 (cert_first_pivot.json: 6/42 admissible,
#         ALL with pivot x6),
#       * after x6: NO admissible second pivot (cert_second_pivot.json: 0/30)
#         — blocked under refinement too => NOLR upheld (cross-validation).
# ---------------------------------------------------------------------------
PROBES_OK || @testset "recorded probe: J2L — FIXTURES MISSING" begin
    @test_skip "SUBTROPICA_PROBES_DIR unset — probe cross-validation skipped (loud skip)"
end
PROBES_OK && @testset "recorded probe: J2L stays blocked under refinement" begin
    obsf = joinpath(PROBES, "lbl3m_j2l", "cert_obstruction_structure.json")
    uff  = joinpath(PROBES, "lbl3m_j2l", "step1_uf_strings.json")
    @test isfile(obsf) && isfile(uff)     # hard fail with env SET (see header)
    obs = JSON.parsefile(obsf)
    uf  = JSON.parsefile(uff)

    syms7 = [:x0, :x1, :x2, :x3, :x4, :x5, :x6]
    R = polynomial_ring(QQ, syms7)[1]
    x6 = gen(R, 7)
    FA = LR.parse_poly_qq(obs["F_A"]["factors"][1][1], R)
    FB = LR.parse_poly_qq(obs["F_B"]["factors"][1][1], R)
    UA = LR.parse_poly_qq(obs["U_A"]["factors"][1][1], R)
    UB = LR.parse_poly_qq(obs["U_B"]["factors"][1][1], R) *
         LR.parse_poly_qq(obs["U_B"]["factors"][2][1], R)
    F7 = FA * x6 + FB          # content 1/2 dropped: letters live up to units
    U7 = UA * x6 + UB

    # provenance tie to the recorded gauge-x6=1 strings (EXACT, up to units)
    Fg6 = evaluate(F7, [7], [one(R)])
    Ug6 = evaluate(U7, [7], [one(R)])
    @test LR.canon_letter(Fg6) == LR.canon_letter(LR.parse_poly_qq(uf["F_phys"], R))
    @test LR.canon_letter(Ug6) == LR.canon_letter(LR.parse_poly_qq(uf["U"], R))

    # Cheng-Wu gauge x0 = 1 (a gauge that KEEPS the unique linear variable)
    Fg = evaluate(F7, [1], [one(R)])
    Ug = evaluate(U7, [1], [one(R)])
    vars6 = [:x1, :x2, :x3, :x4, :x5, :x6]
    polys = [Ug, Fg]

    certL, removed, log = LR.refine_letters(polys, vars6, [:x6, :x1])
    # forced first step x6 (cert_first_pivot.json, gauge x0 rows)
    @test log["stages"][1]["admissible_next_vars"] == ["x6"]
    # after x6: refined bound contains the two recorded obstruction letters
    # (checked factor-wise in case the x0=1 gauge splits them)
    s1 = Set(log["stages"][2]["letters_refined"])
    for q in LR._irreducibles([evaluate(FA, [1], [one(R)])])
        @test string(q) in s1
    end
    for q in LR._irreducibles([evaluate(FB, [1], [one(R)])])
        @test string(q) in s1
    end
    # ... and NO second pivot is admissible (cert_second_pivot.json: 0/30):
    # refinement does NOT spuriously unblock J2L.
    @test log["stages"][2]["admissible_next_vars"] == String[]
    @test log["order_walkable"] === false

    # the F_A obstruction is not certifiable-spurious (blocked-for-real);
    # adjudicate on an irreducible factor in case gauging split it
    fac = first(sort!(collect(LR._irreducibles([evaluate(FA, [1], [one(R)])]));
                      by = LR._letter_sortkey, rev = true))
    cFA = LR.blocking_certificate(polys, vars6, fac)
    @test cFA["status"] == "SURVIVES"
    @test LR.is_blocking_real(polys, vars6, fac) === true
    # reachable refined stages are exactly {} and {x6} — every deeper subset
    # is degree-blocked (matches the recorded depth<=2 death certificate)
    @test cFA["stages_reachable_refined"] == [String[], ["x6"]]
end

# ---------------------------------------------------------------------------
# (5) typed refusals + guards
# ---------------------------------------------------------------------------
@testset "refusals" begin
    R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
    polys = [x + y + 1, x + z * y + z]

    # algebraic letters: typed refusal (scope pin — rational letters only)
    alg = LR.AlgLetter(1, "u^2+3*u+1", :u, "1", "-3", "1", "5")
    @test_throws LR.LRRefineRefusal LR.is_blocking_real(polys, [:x, :y], alg)
    err = try LR.is_blocking_real(polys, [:x, :y], alg) catch e; e end
    @test err.kind === :AlgebraicLetter
    err = try LR.is_blocking_real(polys, [:x, :y], "sqrt(z)-1") catch e; e end
    @test err.kind === :AlgebraicLetter

    # unknown variable / bad order
    @test_throws LR.LRRefineRefusal LR.refine_letters(polys, [:x, :q], [:x])
    @test_throws LR.LRRefineRefusal LR.refine_letters(polys, [:x, :y], [:z])
    err = try LR.refine_letters(polys, [:x, :y], Symbol[]) catch e; e end
    @test err.kind === :EmptyInput

    # constant / monomial letter is not adjudicable
    err = try LR.blocking_certificate(polys, [:x, :y], x * y) catch e; e end
    @test err.kind === :NotAPolynomial

    # mixed rings (4-var ring: Nemo caches same-signature rings, so a
    # same-symbols ring would be === R and NOT a mixed-ring case)
    R2 = polynomial_ring(QQ, [:x, :y, :z, :w])[1]
    err = try LR.blocking_certificate(polys, [:x, :y], gen(R2, 3) + 1) catch e; e end
    @test err.kind === :MixedRings

    # 2^n cap (refused BEFORE any DP work)
    Rbig = polynomial_ring(QQ, [Symbol("v$i") for i in 1:15])[1]
    pbig = sum(gens(Rbig)) + 1
    err = try LR.is_blocking_real([pbig], [Symbol("v$i") for i in 1:15],
                                  gen(Rbig, 1) + 1) catch e; e end
    @test err.kind === :TooManyVars
end

end # top testset

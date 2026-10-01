# test_b2_continue.jl — standalone checks for src/continue.jl (B2 core, C6).
# Suite-integrated;
# the suite's runtests.jl auto-includes it in its own sandbox module.
#
# Run standalone, DEFAULT global Julia env (DESIGN_B1 §3 pattern):
#   ulimit -v 32505856; julia --project=@. tools/subtropica/test/test_b2_continue.jl
#
# Coverage (task spec (a)-(e)):
#  (a) paper eq 3.19-family toy with ONE power-divergent ray, continued by
#      hand — continued integrand + exact prefactor verified field-by-field;
#      plus the paper's own positive control eq 3.37 -> eq 3.38 (prefactor 3).
#  (b) GP-violation synthetic = the paper's eq 3.37 itself (every divergent
#      ray in the span of its compatible partners; paper:1315-1325 confirms
#      all v_rho NotFound) — find_first_np_continuation returns the minimal
#      subset (size 1; deterministically [1] in Subsets order; the paper's
#      pedagogical choice {1,1} = ray 3 is asserted viable too), and the
#      continued output feeds subtract.jl's subtraction_terms (9 CounterTerms).
#  (c) EXACTNESS invariant, 1-d closed forms (Beta/Gamma algebra):
#      original == sum of continued integrands, >=30 floored digits via
#      mpmath shell-out, WITH the mandatory non-rational synthetic-truth
#      mutation control (known pitfall: rational perturbations alias).
#  (d) regulator bookkeeping: the eq 4.10 two-regulator integrand structure
#      reproduces the recorded ray-census Trop values (lambda^q, lambda^-eps,
#      lambda^{-eps+q} <=> TropI = -q, eps, eps-q) on all three divergent
#      rays; q -> q*eps promotion recorded in provenance and matching the
#      paper's promoted trops {eps, -eps(-1+q), -eps q}.
#  (e) refusal battery: >2 regulators, double promotion, oversized subset
#      search (size + budget guards), continuation-depth cap, ill-defined
#      ray, non-integer divergence order.

using Test, Nemo, JSON

include(joinpath(@__DIR__, "..", "src", "types.jl"))
include(joinpath(@__DIR__, "..", "src", "b1_types.jl"))
include(joinpath(@__DIR__, "..", "src", "tropical.jl"))
include(joinpath(@__DIR__, "..", "src", "subtract.jl"))
include(joinpath(@__DIR__, "..", "src", "continue.jl"))

const EE = EpsExp

# ---------------------------------------------------------------------------
# Hand-built fixtures (SPEC_B1 §8 conventions: dlog measure, hand-pinned
# TropicalData — no Oscar dependence; ray order is fixture-local).
# ---------------------------------------------------------------------------

# 2-var ring (no kinvars)
R2, (x1, x2) = polynomial_ring(Nemo.QQ, ["x1", "x2"])

# F1 geometry (SPEC_B1 §8.1 = paper eq 3.37 triangle (0,0),(1,0),(0,1)):
td_F1 = TropicalData(
    [BigInt[-1, 0], BigInt[0, -1], BigInt[1, 1]],
    [BigInt[0, 0], BigInt[1, 0], BigInt[0, 1]],
    [[1, 3], [1, 2], [2, 3]],
    Rational{BigInt}[0, 0, 1],
    Tuple{Vector{Rational{BigInt}},Rational{BigInt}}[], 2)

# eq 3.37 integrand (dlog): x1^eps x2^eps (1+x1+x2)^{-3eps}
E_337 = EulerIntegrand(Prefactor(), [EE(0, 1//1), EE(0, 1//1)],
    [(1 + x1 + x2, EE(0, -3//1))], [:x1, :x2], Symbol[], R2)

# eq 3.19-family toy geometry: P = x1^2 + x2 + x1*x2, Newt = triangle
# (2,0),(0,1),(1,1); outer normals rho1=(-1,-2), rho2=(1,1), rho3=(0,1)
# (paper Fig. 3 / eq 3.19-3.21).
P_toy = x1^2 + x2 + x1 * x2
td_toy = TropicalData(
    [BigInt[-1, -2], BigInt[1, 1], BigInt[0, 1]],
    [BigInt[2, 0], BigInt[0, 1], BigInt[1, 1]],
    [[1, 2], [1, 3], [2, 3]],
    Rational{BigInt}[-2, 2, 1],
    Tuple{Vector{Rational{BigInt}},Rational{BigInt}}[], 2)

# The toy itself: nu = (1+eps, 1+eps) makes rho1 POWER-divergent, order 1:
# TropI(rho1) = -(1+eps) - 2(1+eps) + (2+eps)*2 = 1 - eps  (a = 1 > 0);
# TropI(rho2) = -2, TropI(rho3) = -1 (convergent). One power ray, GP holds
# for the (empty-removal) surviving set via the simplest-w branch.
E_toy = EulerIntegrand(Prefactor(), [EE(1, 1//1), EE(1, 1//1)],
    [(P_toy, EE(-2, -1//1))], [:x1, :x2], Symbol[], R2)

# 1-var ring + segment geometry for the exactness cases (SPEC_B1 §2.3(a)
# convention: rays [[-1],[1]]).
R1, (xx,) = polynomial_ring(Nemo.QQ, ["x"])
td_1d = TropicalData(
    [BigInt[-1], BigInt[1]],
    [BigInt[0], BigInt[1]],
    [[1], [2]],
    Rational{BigInt}[0, 1],
    Tuple{Vector{Rational{BigInt}},Rational{BigInt}}[], 1)

# ---------------------------------------------------------------------------
@testset "(a) eq 3.19-family toy: one power ray, hand-verified step" begin
    tv = trop_values(E_toy, td_toy)
    @test tv == [EE(1, -1//1), EE(-2, 0//1), EE(-1, 0//1)]   # eq 3.20 by hand

    # single manual step along rho1 = (-1,-2):
    # Q = x1^2 + x2 + lambda*x1*x2 (t = -2; x1*x2 has m.rho = -3), so
    # D = dQ/dlambda|_1 = x1*x2 — a pure monomial: folded into nu; prefactor
    # e_P / TropI(rho1) = (-(2+eps))/(1-eps)  [SIGN PIN — see continue.jl header]
    cts = continue_ray(E_toy, td_toy, 1)
    @test length(cts) == 1
    ct = cts[1]
    @test ct.pref.c == 1//1
    @test ct.pref.num == [EE(-2, -1//1)]                     # e_P = -(2+eps)
    @test ct.pref.den == [EE(1, -1//1)]                      # TropI(rho1) = 1-eps
    @test b2pref_eval(ct.pref, big(1)//5) == -11//4          # exact rational
    @test ct.integrand.nu == [EE(2, 1//1), EE(2, 1//1)]      # x^(1,1) folded in
    @test length(ct.integrand.polys) == 1
    @test ct.integrand.polys[1][1] == P_toy                  # SAME poly support
    @test ct.integrand.polys[1][2] == EE(-3, -1//1)          # exponent shift -1
    @test ct.provenance["steps"][1]["deriv_poly_kept"] == false

    # continued integrand is now log-divergent at worst on rho1 (shaved by 1):
    tv2 = trop_values(ct.integrand, td_toy)
    @test tv2[1] == EE(0, -1//1)                             # 1-eps -> -eps
    @test tv2[2].a < 0 && tv2[3].a < 0

    # orchestrator: GP already holds at zero removal (simplest-w), so the
    # power ray gets exactly `order a = 1` step and nothing else:
    outs = expand_integral_b2(E_toy, td_toy)
    @test length(outs) == 1
    @test outs[1].pref == ct.pref
    @test outs[1].integrand.nu == ct.integrand.nu
    @test outs[1].integrand.polys == ct.integrand.polys
    @test outs[1].provenance["np_search"]["rays_removed"] == Int[]
    @test outs[1].provenance["np_search"]["np_ray_sequence"] == [1]
    @test outs[1].provenance["trdata_reusable"] == true

    # POSITIVE CONTROL (paper eq 3.37 -> 3.38, paper:1329-1331): continuing
    # eq 3.37 along ray 3 = (1,1) gives {{3, x1^eps x2^eps (1+x1+x2)^{-1-3eps}}}
    pc = continue_ray(E_337, td_F1, 3)
    @test length(pc) == 1
    sp = b2pref_simplify(pc[1].pref)
    @test sp.c == 3//1 && isempty(sp.num) && isempty(sp.den) # the literal 3
    @test pc[1].integrand.nu == E_337.nu                     # D = const: no shift
    @test pc[1].integrand.polys == [(1 + x1 + x2, EE(-1, -3//1))]
end

# ---------------------------------------------------------------------------
@testset "(b) GP-violation synthetic (eq 3.37): minimal subset + B1 handoff" begin
    tv = trop_values(E_337, td_F1)
    @test tv == [EE(0, -1//1), EE(0, -1//1), EE(0, -1//1)]   # all log-divergent
    sel = div_facets(tv)
    @test sel == [1, 2, 3]

    # raw input violates GP (paper:1315-1325: all three v_rho NotFound):
    @test_throws GeometricPropertyViolated produce_ws(td_F1, sigma_div(td_F1, sel))

    # minimal-removal search: size-1 removal suffices; Subsets order makes
    # the FIRST size-1 subset [1] the deterministic winner (scan: k=0 fails
    # GP, then S=[1] succeeds => 2 subsets scanned).
    fr = find_first_np_continuation(td_F1, sel)
    @test fr.rays == [1]
    @test fr.scanned == 2
    @test Set(keys(fr.ws)) == Set([2, 3])
    @test fr.ws[2] == BigInt[-1, 1] && fr.ws[3] == BigInt[-1, 0]

    # the paper's pedagogical choice (remove ray 3 = {1,1}) is also viable:
    @test produce_ws(td_F1, sigma_div(td_F1, [1, 2])) isa Dict

    # full orchestration: ray 1 is log (order 0) => exactly ONE step; D = x1
    # monomial => nu shift; prefactor (-3eps)/(-eps) = 3.
    outs = expand_integral_b2(E_337, td_F1)
    @test length(outs) == 1
    o = outs[1]
    @test b2pref_simplify(o.pref).c == 3//1
    @test o.integrand.nu == [EE(1, 1//1), EE(0, 1//1)]       # x1^{1+eps} x2^eps
    @test o.integrand.polys == [(1 + x1 + x2, EE(0, -3//1) - 1)]
    @test o.provenance["np_search"]["rays_removed"] == [1]
    @test o.provenance["trdata_reusable"] == true            # same poly support

    # B1 handoff: the continued integrand passes the FULL C5 ladder and
    # subtract.jl's C7 assembly (task item 4: consume subtraction_terms).
    dd = divergence_data(o.integrand, td_F1)
    @test dd.div_facets == [2, 3]
    @test [f.facets for f in dd.sigma_div] == [Int[], [2], [3], [2, 3]]
    cts = subtraction_terms(o.integrand, td_F1, dd.sigma_div, dd.w)
    @test length(cts) == 9            # 4 (∅) + 2 ([2]) + 2 ([3]) + 1 ([2,3])
    @test [c.sign for c in cts if c.face == Int[]] == [1, -1, -1, 1]  # WARNING-1
    @test all(c.sign == warning1_sign(c.subface) for c in cts)
end

# ---------------------------------------------------------------------------
@testset "(c) exactness invariant: 1-d closed forms, >=30d + mutation control" begin
    # Case A (order-1 power ray, TWO poly entries sharing the same base —
    # exercises the per-poly product-rule split):
    #   I_A = ∫ dlog x · x^{-1+eps} (1+x)^{-1-eps} (1+x)^{-2eps}
    #       = B(-1+eps, 2+2eps)
    EA = EulerIntegrand(Prefactor(), [EE(-1, 1//1)],
        [(1 + xx, EE(-1, -1//1)), (1 + xx, EE(0, -2//1))], [:x], Symbol[], R1)
    @test trop_values(EA, td_1d) == [EE(1, -1//1), EE(-2, -2//1)]
    outsA = expand_integral_b2(EA, td_1d)
    @test length(outsA) == 2                                 # one per poly entry
    # both terms: D = x (monomial), nu -> eps; shifted exponent on its own poly
    @test outsA[1].pref.num == [EE(-1, -1//1)]               # e_1 = -1-eps
    @test outsA[2].pref.num == [EE(0, -2//1)]                # e_2 = -2eps
    @test outsA[1].pref.den == [EE(1, -1//1)] == outsA[2].pref.den
    @test outsA[1].integrand.nu == [EE(0, 1//1)] == outsA[2].integrand.nu
    @test outsA[1].integrand.polys ==
          [(1 + xx, EE(-2, -1//1)), (1 + xx, EE(0, -2//1))]
    @test outsA[2].integrand.polys ==
          [(1 + xx, EE(-1, -1//1)), (1 + xx, EE(-1, -2//1))]
    # each term integrates to B(eps, 2+2eps); exact coefficient sum:
    cA1 = b2pref_eval(outsA[1].pref, big(1)//7)
    cA2 = b2pref_eval(outsA[2].pref, big(1)//7)
    @test cA1 + cA2 == -(1 + 3 * (big(1)//7)) / (1 - big(1)//7)  # -(1+3e)/(1-e)

    # Case B (order-2 power ray => TWO chained steps through continue_rays):
    #   I_B = ∫ dlog x · x^{-2+eps} (1+x)^{-1-3eps} = B(-2+eps, 3+2eps)
    EB = EulerIntegrand(Prefactor(), [EE(-2, 1//1)],
        [(1 + xx, EE(-1, -3//1))], [:x], Symbol[], R1)
    outsB = expand_integral_b2(EB, td_1d)
    @test length(outsB) == 1
    @test outsB[1].pref.num == [EE(-1, -3//1), EE(-2, -3//1)]
    @test outsB[1].pref.den == [EE(2, -1//1), EE(1, -1//1)]
    @test outsB[1].integrand.nu == [EE(0, 1//1)]
    @test outsB[1].integrand.polys == [(1 + xx, EE(-3, -3//1))]
    cB = b2pref_eval(outsB[1].pref, big(1)//7)
    @test cB == 85//39                                       # hand value at eps=1/7

    # Case F1 (2-d pipeline sum): A(eps) = Gamma(eps)^3/Gamma(3eps) vs
    # 3 * Gamma(1+eps)Gamma(eps)^2/Gamma(1+3eps)  (the (b) continuation).
    cF1 = b2pref_eval(b2pref_simplify(
        expand_integral_b2(E_337, td_F1)[1].pref), big(1)//7)
    @test cF1 == 3//1

    # mpmath oracle (dps 80, eps = 1/7): >=30 floored digits on each case,
    # AND the non-rational (1 + 1e-12*pi) mutation MUST fail by a
    # strictly-greater-than-slack margin (a control that only prints is not a
    # gate; rational perturbations alias — hence the pi).
    rs(x::Rational) = "Fraction($(numerator(x)),$(denominator(x)))"
    py = """
import json
from fractions import Fraction
from mpmath import mp, mpf, gamma, beta, log10, fabs, pi
mp.dps = 80
def R(fr): return mpf(fr.numerator)/mpf(fr.denominator)
e = mpf(1)/7
def digits(lhs, rhs):
    d = fabs(lhs - rhs)/fabs(lhs)
    return 999 if d == 0 else int(-log10(d))
out = {}
mut = 1 + mpf(10)**-12 * pi
# case A
lhs = beta(-1+e, 2+2*e)
rhs = (R($(rs(cA1))) + R($(rs(cA2)))) * beta(e, 2+2*e)
out["A"] = [digits(lhs, rhs), digits(lhs, rhs*mut)]
# case B
lhs = beta(-2+e, 3+2*e)
rhs = R($(rs(cB))) * beta(e, 3+2*e)
out["B"] = [digits(lhs, rhs), digits(lhs, rhs*mut)]
# case F1
lhs = gamma(e)**3/gamma(3*e)
rhs = R($(rs(cF1))) * gamma(1+e)*gamma(e)*gamma(e)/gamma(1+3*e)
out["F1"] = [digits(lhs, rhs), digits(lhs, rhs*mut)]
print(json.dumps(out))
"""
    dir = mktempdir()
    script = joinpath(dir, "b2_exactness_oracle.py")
    write(script, py)
    res = JSON.parse(read(`bash -c "ulimit -v 32505856; python3 $script"`, String))
    for k in ("A", "B", "F1")
        dok, dmut = res[k]
        @test dok >= 30                       # the >=30d exactness bar
        @test dmut < 30                       # mutation control must FAIL
        @test dok - dmut > 20                 # strictly-greater-than-slack
    end
end

# ---------------------------------------------------------------------------
@testset "(d) two-regulator bookkeeping: eq 4.10 census + q->q*eps promotion" begin
    # eq 4.10 structure (gauge x3=1, dlog, 2 vars,
    # kinvars X11,X12,X22):  x1^{1-q} x2 * P1^{-1-eps} * P2^{2eps}
    Rq, (y1, y2, X11, X12, X22) =
        polynomial_ring(Nemo.QQ, ["x1", "x2", "X11", "X12", "X22"])
    P1 = y1 * y2 + y1 + y2
    P2 = y2 * X11 + X11 + 2 * y2 * X12 + y1 * X22 + y2 * X22
    E410 = Euler2Integrand(Prefactor(),
        [Eps2(1//1, 0//1, -1//1), Eps2(1//1, 0//1, 0//1)],
        [(P1, Eps2(-1//1, -1//1, 0//1)), (P2, Eps2(0//1, 2//1, 0//1))],
        [:x1, :x2], [:X11, :X12, :X22], Rq, [:eps, :q])

    # recorded ray-census Trop values (lambda^{-Trop} exponents lambda^-eps,
    # lambda^{-eps+q}, lambda^q reproduced by the paper's own STFactor):
    @test trop2_on_ray(E410, [0, 1])  == Eps2(0//1, 1//1,  0//1)   #  eps
    @test trop2_on_ray(E410, [1, 0])  == Eps2(0//1, 1//1, -1//1)   #  eps - q
    @test trop2_on_ray(E410, [1, 1])  == Eps2(0//1, 0//1, -1//1)   # -q
    # paper footnote control: ray (-1,-1) scales as lambda^{1-eps-q}:
    @test trop2_on_ray(E410, [-1, -1]) == Eps2(-1//1, 1//1, 1//1)

    # naive eps-only view is BLIND to the (1,1) divergence (paper: the
    # computation aborts because it is not regulated by eps):
    t11 = trop2_on_ray(E410, [1, 1])
    @test unregulated_by_eps(t11)
    @test is_log(trop2_on_ray(E410, [0, 1])) && !unregulated_by_eps(trop2_on_ray(E410, [0, 1]))

    # q -> q*eps promotion: same stored triples, reinterpreted a + (b+c*q)eps
    # (paper's promoted trops {eps, -eps(-1+q), -eps q}); provenance recorded.
    E410p = promote_second_regulator(E410)
    @test is_promoted(E410p)
    @test E410p.provenance["promotion"] == "q -> q*eps"
    tp = [trop2_on_ray(E410p, r) for r in ([0, 1], [1, 0], [1, 1])]
    @test [(t.a, t.b, t.c) for t in tp] ==
          [(0//1, 1//1, 0//1),      # eps
           (0//1, 1//1, -1//1),     # eps(1-q) = -eps(-1+q)
           (0//1, 0//1, -1//1)]     # -eps*q
    @test all(is_log, tp)           # promoted view: ALL three regulated

    # continuation along the q-only ray (1,1) — the paper's
    # STTropicalContinuation[..., {{1,1}}] call shape: one term per poly.
    cq = continue_ray(E410, [1, 1])
    @test length(cq) == 2
    @test cq[1].pref.num == [Eps2(-1//1, -1//1, 0//1)]        # e_1 = -1-eps
    @test cq[1].pref.den == [Eps2(0//1, 0//1, -1//1)]         # TropI = -q
    @test b2pref_eval(cq[1].pref, big(1)//3, big(1)//5) == 20//3
    @test cq[1].integrand.polys[1] == (P1, Eps2(-2//1, -1//1, 0//1))
    @test cq[1].integrand.polys[end] == (y1 + y2, Eps2(1//1, 0//1, 0//1))
    # P2's derivative along (1,1) is the kinvar-only X11 — stays a poly factor:
    @test cq[2].integrand.polys[end] == (X11, Eps2(1//1, 0//1, 0//1))
end

# ---------------------------------------------------------------------------
@testset "(d2) two-reg driver glue (§C6b): Eps2 chaining + forced NP search" begin
    Rq, (y1, y2, X11, X12, X22) =
        polynomial_ring(Nemo.QQ, ["x1", "x2", "X11", "X12", "X22"])
    P1 = y1 * y2 + y1 + y2
    P2 = (y2 + 1) * X11 + 2 * y2 * X12 + (y1 + y2) * X22
    E410 = Euler2Integrand(Prefactor(),
        [Eps2(1//1, 0//1, -1//1), Eps2(1//1, 0//1, 0//1)],
        [(P1, Eps2(-1//1, -1//1, 0//1)), (P2, Eps2(0//1, 2//1, 0//1))],
        [:x1, :x2], [:X11, :X12, :X22], Rq, [:eps, :q])

    # Eps2 continue_rays: ONE (1,1) step through the chaining overload must
    # reproduce the two direct continue_ray outputs (prefactors multiplied
    # against the identity, provenance threaded through the integrand).
    direct = continue_ray(E410, [1, 1])
    chained = continue_rays(E410, [[1, 1]])
    @test length(chained) == length(direct) == 2
    for k in 1:2
        @test chained[k].pref == direct[k].pref
        @test chained[k].integrand.nu == direct[k].integrand.nu
        @test chained[k].integrand.polys == direct[k].integrand.polys
        @test length(chained[k].provenance["steps"]) == 1
    end
    # two-step chaining multiplies prefactors (identity: continue term 2
    # again along (1,1) — its shaved TropI is -1-q, a=-1: still steppable)
    two = continue_rays([chained[2]], [[1, 1]])
    @test all(length(t.pref.den) == 2 for t in two)
    @test all(t.pref.den[1] == Eps2(0//1, 0//1, -1//1) for t in two)

    # forced NP search on the eq 4.10 hexagon census (divergent set {1,2,3},
    # fixture order (0,1),(1,0),(1,1)): the GP infeasibility (1,1)=(1,0)+(0,1)
    # clears when ANY one ray is removed, so the UNFORCED lex-minimal removal
    # is {1} — NOT the paper's continued ray. The paper continues (1,1)
    # because it is the q-only-regulated ray: exactly the FORCED-base law the
    # two-reg driver applies (force = unregulated_by_eps set).
    td_410 = TropicalData(
        [BigInt[0, 1], BigInt[1, 0], BigInt[1, 1],
         BigInt[-1, 0], BigInt[0, -1], BigInt[-1, -1]],
        [BigInt[0, 1], BigInt[0, 2], BigInt[1, 2],
         BigInt[2, 1], BigInt[2, 0], BigInt[1, 0]],
        [[2, 3], [4, 5], [3, 4], [1, 2], [5, 6], [1, 6]],
        Rational{BigInt}[2, 2, 3, 0, 0, -1],
        Tuple{Vector{Rational{BigInt}},Rational{BigInt}}[], 2)
    fr0 = find_first_np_continuation(td_410, [1, 2, 3])
    @test fr0.rays == [1]                      # unforced lex-minimal
    frf = find_np_continuation_forced(td_410, [1, 2, 3], [3])
    @test frf.rays == [3]                      # forcing (1,1) wins — the paper's step
    fr1 = find_np_continuation_forced(td_410, [1, 2, 3], [1])
    @test fr1.rays == [1]                      # forcing the lex-minimal is a no-op
    # force ⊄ div refuses; guards propagate
    @test_throws ArgumentError find_np_continuation_forced(td_410, [1, 2, 3], [4])
    err = try find_np_continuation_forced(td_410, [1, 2, 3], [3]; subset_budget=0)
        nothing catch e; e end
    @test err isa B2Refusal && err.kind == :ContinuationSearchTooLarge
end

# ---------------------------------------------------------------------------
@testset "(e) refusal battery" begin
    Rq, (y1, y2, X11, X12, X22) =
        polynomial_ring(Nemo.QQ, ["x1", "x2", "X11", "X12", "X22"])
    P1 = y1 * y2 + y1 + y2

    # >2 regulators: typed refusal at construction.
    err = try
        Euler2Integrand(Prefactor(), [Eps2(1), Eps2(1)],
            [(P1, Eps2(0//1, 1//1, 0//1))], [:x1, :x2],
            [:X11, :X12, :X22], Rq, [:eps, :q, :r])
        nothing
    catch e; e end
    @test err isa B2Refusal && err.kind == :TooManyRegulators
    @test err.detail["count"] == 3

    # double promotion refuses (eps^2 has no exponent slot):
    E2 = Euler2Integrand(Prefactor(), [Eps2(1//1, 0//1, -1//1), Eps2(1)],
        [(P1, Eps2(-1//1, -1//1, 0//1))], [:x1, :x2],
        [:X11, :X12, :X22], Rq, [:eps, :q])
    err = try promote_second_regulator(promote_second_regulator(E2)); nothing
    catch e; e end
    @test err isa B2Refusal && err.kind == :TooManyRegulators

    # oversized subset search — size guard (17 > max_div_facets = 16),
    # BEFORE any geometry is touched:
    td_big = TropicalData(
        [BigInt[i, 1] for i in 1:17], [BigInt[0, 0]], [[1] for _ in 1:17],
        zeros(Rational{BigInt}, 17),
        Tuple{Vector{Rational{BigInt}},Rational{BigInt}}[], 2)
    err = try find_first_np_continuation(td_big, collect(1:17)); nothing
    catch e; e end
    @test err isa B2Refusal && err.kind == :ContinuationSearchTooLarge
    @test err.detail["count"] == 17 && err.detail["limit"] == 16

    # oversized subset search — scan budget guard:
    err = try find_first_np_continuation(td_F1, [1, 2, 3]; subset_budget = 1)
        nothing
    catch e; e end
    @test err isa B2Refusal && err.kind == :ContinuationSearchTooLarge
    @test err.detail["subset_budget"] == 1

    # continuation-depth cap (`order`): the toy needs 1 step, cap at 0:
    err = try expand_integral_b2(E_toy, td_toy; order = 0); nothing
    catch e; e end
    @test err isa B2Refusal && err.kind == :ContinuationDepthExceeded
    @test err.detail["steps_required"] == 1

    # ill-defined ray (TropI == 0): 1/TropI undefined:
    E0 = EulerIntegrand(Prefactor(), [EE(0, 0//1)],
        Tuple{QQMPolyRingElem,EpsExp}[], [:x], Symbol[], R1)
    err = try continue_ray(E0, [1]); nothing
    catch e; e end
    @test err isa B2Refusal && err.kind == :ContinuationIllDefined

    # non-integer positive divergence order (repeat count ill-defined):
    Eh = EulerIntegrand(Prefactor(), [EE(-3//2, 1//1)],
        [(1 + xx, EE(-1, -3//1))], [:x], Symbol[], R1)
    err = try expand_integral_b2(Eh, td_1d); nothing
    catch e; e end
    @test err isa B2Refusal && err.kind == :NonIntegerDivergenceOrder
    @test err.detail["order"] == 3//2

    # sanity: refusals print with the B2 banner (loud, typed):
    @test occursin("subtropica B2: ContinuationDepthExceeded",
                   sprint(showerror, B2Refusal(:ContinuationDepthExceeded,
                                               Dict{String,Any}())))
end

println("b2_continue_check: all testsets finished")

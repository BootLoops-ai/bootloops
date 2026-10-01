# test_b1_subtract.jl — OWNER: unit-S. Tests for src/subtract.jl (C7 Möbius
# counterterm assembly). Self-contained: fixtures F1/F2/F3 are the hand-worked
# tables of SPEC_B1 §8 (TropicalData/SigmaFace/ws built literally from the
# spec — NO dependence on unit-P/unit-T code). Closed-form Gamma/Beta expected
# values evaluated with Nemo arb at exact rational ε points.
#
# Standalone run (DEFAULT global Julia env — DESIGN_B1.md §3; NO subtropica
# project env before integration):
#   ulimit -v 32505856; julia tools/subtropica/test/test_b1_subtract.jl
#
# Under the shared suite (post-integration) runtests.jl auto-includes this
# file after `using SubTropica`; the @isdefined guards below then skip the
# by-path includes and the tests run against the module's definitions.

using Test
using Nemo
using JSON   # unused here, but the DESIGN_B1 §3 standalone-header pattern

if !@isdefined(EpsExp)
    include(joinpath(@__DIR__, "..", "src", "types.jl"))
end
if !@isdefined(RayClass)
    include(joinpath(@__DIR__, "..", "src", "b1_types.jl"))
end
if !@isdefined(trop_on_ray)   # standalone-run support:
    # subtract.jl now DELEGATES _c7_trop_poly/_c7_trop_on_ray/
    # _c7_restrict_poly to unit-T's tropical.jl exports — standalone runs
    # need tropical.jl in scope first (parity test: test_b1_parity.jl).
    include(joinpath(@__DIR__, "..", "src", "tropical.jl"))
end
if !@isdefined(_c7_J_of)   # PRIVATE sentinel: `subtraction_terms` may get
    # exported at integration while the _c7_* helpers (exercised below for
    # the WARNING-2 negative control) stay module-private — keying the guard
    # on a private name guarantees the by-path include runs whenever those
    # helpers are not already in scope.
    include(joinpath(@__DIR__, "..", "src", "subtract.jl"))
end

# ---------------------------------------------------------------------------
# helpers (as_ prefix — runtests.jl includes every unit's tests into Main)
# ---------------------------------------------------------------------------

const AS_RR = ArbField(384)                       # ~115 decimal digits

as_eps(t::EpsExp, e) = parent(e)(QQ(t.a)) + parent(e)(QQ(t.b)) * e
as_gamma(x) = gamma(x)
as_B(p, q)  = gamma(p) * gamma(q) / gamma(p + q)
as_pow(c, q) = exp(q * log(c))                    # c^q, c > 0 arb

# exponent-vector support of a (gauge-fixed) Nemo poly, projected onto the
# CounterTerm's J slots; asserts no dependence outside J (gauge-fix worked)
function as_support(P, J::Vector{Int})
    R = parent(P)
    sup = Tuple{Rational{BigInt},Vector{Int}}[]
    for k in 1:length(P)
        c  = coeff(P, k)
        ev = exponent_vector(P, k)
        for i in 1:nvars(R)
            (i in J) || @assert ev[i] == 0 "poly depends on gauged-out/kin var $i"
        end
        push!(sup, (BigInt(numerator(c)) // BigInt(denominator(c)),
                    Int[ev[j] for j in J]))
    end
    return sup
end

# Closed-form (Mellin/Beta-table) evaluation of ONE CounterTerm at numeric
# arb ε. Iterated 1-var rule (∫₀^∞ x^{p−1}(A + B·x)^{e} dx =
# A^{e+p}·B^{−p}·Β(p, −e−p), B a monomial in the later vars) — covers every
# fixture shape (each poly is multilinear; exactly one x_v-dependent factor
# at each stage). Throws (loudly) on scaleless/unsupported shapes.
function as_ct_value(ct::CounterTerm, e)
    RR  = parent(e)
    val = RR(ct.sign) * RR(ct.vol_det)
    for t in ct.vol_trops                         # Vol = vol_det/∏(−TropI)
        val = val / as_eps(t, e)
    end
    k    = length(ct.J)
    mono = [as_eps(m, e) for m in ct.mono]
    facs = Any[Any[as_support(P, ct.J), as_eps(ex, e)] for (P, ex) in ct.polys]
    for v in 1:k
        p   = mono[v] + 1                         # flat measure: ∫x^{mono}dx
        dep = [i for i in eachindex(facs) if any(t[2][v] != 0 for t in facs[i][1])]
        @assert length(dep) == 1 "as_ct_value: var $v has $(length(dep)) dependent factors (scaleless or unsupported)"
        sup, ex = facs[dep[1]]
        @assert all(t[2][v] in (0, 1) for t in sup) "as_ct_value: not multilinear in var $v"
        A  = [t for t in sup if t[2][v] == 0]
        Bs = [t for t in sup if t[2][v] == 1]
        @assert length(Bs) == 1 && !isempty(A) "as_ct_value: unsupported A+Bx shape"
        cB, evB = Bs[1]
        @assert all(evB[w] == 0 for w in 1:v-1) "B-monomial touches integrated var"
        val = val * as_B(p, -ex - p)
        val = val * as_pow(RR(QQ(cB)), -p)
        for w in v+1:k
            mono[w] = mono[w] - p * evB[w]
        end
        facs[dep[1]] = Any[A, ex + p]
    end
    for (sup, ex) in facs                         # leftovers must be constants
        @assert length(sup) == 1 && all(==(0), sup[1][2]) "leftover non-constant factor"
        val = val * as_pow(RR(QQ(sup[1][1])), ex)
    end
    return val
end

as_face_value(cts, face, e) =
    sum(as_ct_value(ct, e) for ct in cts if ct.face == face)

as_close(a, b; tol = AS_RR(10)^-40) = abs(a - b) < tol

# EpsExp shorthands
as_ee(a, b) = EpsExp(Rational{BigInt}(a), Rational{BigInt}(b))
const AS_EPS  = as_ee(0, 1)                        # ε
const AS_ZERO = as_ee(0, 0)

# ---------------------------------------------------------------------------
# Shared geometry: Newt(1+x1+x2) = triangle (0,0),(1,0),(0,1) (F1 = F2;
# SPEC_B1 §8.1 table, matches paper:1303-1312). Ray order fixture-local.
# ---------------------------------------------------------------------------
as_td_tri() = TropicalData(
    [BigInt[-1, 0], BigInt[0, -1], BigInt[1, 1]],          # outer normals
    [BigInt[0, 0], BigInt[1, 0], BigInt[0, 1]],            # v1, v2, v3
    [[1, 3], [1, 2], [2, 3]],                              # argmax vertex sets
    Rational{BigInt}[0, 0, 1],
    Tuple{Vector{Rational{BigInt}},Rational{BigInt}}[],
    2)

# F3 geometry: Newt(1+x2+x1x2) = triangle (0,0),(0,1),(1,1) (SPEC_B1 §8.3).
# Ray indices: 1 = A = (−1,0), 2 = B = (0,1) convergent, 3 = C = (1,−1).
as_td_f3() = TropicalData(
    [BigInt[-1, 0], BigInt[0, 1], BigInt[1, -1]],
    [BigInt[0, 0], BigInt[0, 1], BigInt[1, 1]],
    [[1, 2], [2, 3], [1, 3]],
    Rational{BigInt}[0, 1, 0],
    Tuple{Vector{Rational{BigInt}},Rational{BigInt}}[],
    2)

@testset "C7 subtract.jl (unit-S)" begin

    R, (x1, x2) = polynomial_ring(QQ, ["x1", "x2"])

    # F1: I = x1^ε x2^ε (1+x1+x2)^{−1−3ε}, dlog measure (SPEC_B1 §8.1;
    # prefactor 1 — the paper's NP factor 3 lives in Prefactor, D9).
    E1  = EulerIntegrand(Prefactor(), [AS_EPS, AS_EPS],
                         [(1 + x1 + x2, as_ee(-1, -3))],
                         [:x1, :x2], Symbol[], R)
    td1 = as_td_tri()
    sd1 = [SigmaFace(Int[], [1, 2, 3]), SigmaFace([1], [1, 3]),
           SigmaFace([2], [1, 2]), SigmaFace([1, 2], [1])]
    ws1 = Dict{Int,Vector{BigInt}}(1 => BigInt[1, 0], 2 => BigInt[0, 1])

    # F2: I = x1^{1+ε} x2^ε (1+x1+x2)^{−2−3ε} (SPEC_B1 §8.2) — WARNING-2
    # discriminator (I = [2], J = [1] for face {2}).
    E2  = EulerIntegrand(Prefactor(), [as_ee(1, 1), AS_EPS],
                         [(1 + x1 + x2, as_ee(-2, -3))],
                         [:x1, :x2], Symbol[], R)
    td2 = as_td_tri()
    sd2 = [SigmaFace(Int[], [1, 2, 3]), SigmaFace([2], [1, 2])]
    ws2 = Dict{Int,Vector{BigInt}}(2 => BigInt[0, 1])

    # F3: I = x1^ε x2^{2ε} (1+x2+x1x2)^{−1−5ε} (SPEC_B1 §8.3) — nontrivial
    # x^{u_σ} twist, non-axis w_A = (1,1), double pole.
    E3  = EulerIntegrand(Prefactor(), [AS_EPS, as_ee(0, 2)],
                         [(1 + x2 + x1 * x2, as_ee(-1, -5))],
                         [:x1, :x2], Symbol[], R)
    td3 = as_td_f3()
    sd3 = [SigmaFace(Int[], [1, 2, 3]), SigmaFace([1], [1, 2]),
           SigmaFace([3], [1, 3]), SigmaFace([1, 3], [1])]
    ws3 = Dict{Int,Vector{BigInt}}(1 => BigInt[1, 1], 3 => BigInt[0, 1])

    cts1 = subtraction_terms(E1, td1, sd1, ws1)
    cts2 = subtraction_terms(E2, td2, sd2, ws2)
    cts3 = subtraction_terms(E3, td3, sd3, ws3)

    # -----------------------------------------------------------------------
    @testset "WARNING-1 sign law (wl:10823/10833)" begin
        # net sign = −(−1)^(|s|+1) = (−1)^|s|; the wrong v1 reading
        # ((−1)^(|s|+1)) ADDS single-facet counterterms.
        @test warning1_sign(Int[])   == +1
        @test warning1_sign([7])     == -1
        @test warning1_sign([2, 5])  == +1
        @test all(ct.sign == warning1_sign(ct.subface) for ct in cts1)
        @test all(ct.sign == warning1_sign(ct.subface) for ct in cts3)
    end

    @testset "helpers: lex subsets (G2), det, jacobian factor" begin
        # Mathematica Subsets[Range[4],{2}] hand list (gate G2 pin)
        @test _c7_subsets_lex(4, 2) ==
              [[1,2],[1,3],[1,4],[2,3],[2,4],[3,4]]
        @test _c7_subsets_lex(2, 1) == [[1],[2]]
        @test _c7_det([Rational{BigInt}[0,-1], Rational{BigInt}[1,0]]) == 1
        @test _c7_det([Rational{BigInt}[0,-1], Rational{BigInt}[0,1]]) == 0
        # mixed-sign w: (1+x^{(1,−1)}) → poly x2 + x1, compensation w⁻=(0,1)
        jp, wminus = _c7_jac_factor(R, 2, BigInt[1, -1])
        @test jp == x1 + x2 && wminus == BigInt[0, 1]
        # all-nonnegative w (fixture case): poly = 1 + x^w, compensation 0
        jp0, wm0 = _c7_jac_factor(R, 2, BigInt[0, 1])
        @test jp0 == 1 + x2 && all(iszero, wm0)
    end

    @testset "restriction: initial forms + order-independence (G8)" begin
        rA = BigInt[-1, 0]; rB = BigInt[0, -1]
        @test _c7_restrict_poly(1 + x1 + x2, rA, 2) == 1 + x2   # wl:10506-10511
        @test _c7_restrict_poly(1 + x1 + x2, rB, 2) == 1 + x1
        # sequential face restriction commutes (SPEC_B1 §4 step 9 note)
        @test _c7_restrict_poly(_c7_restrict_poly(1 + x1 + x2, rA, 2), rB, 2) ==
              _c7_restrict_poly(_c7_restrict_poly(1 + x1 + x2, rB, 2), rA, 2)
        P3 = 1 + x2 + x1 * x2
        rC = BigInt[1, -1]
        @test _c7_restrict_poly(P3, rA, 2) == 1 + x2            # F3 table
        @test _c7_restrict_poly(P3, rC, 2) == 1 + x1 * x2
        @test _c7_restrict_poly(_c7_restrict_poly(P3, rA, 2), rC, 2) ==
              _c7_restrict_poly(_c7_restrict_poly(P3, rC, 2), rA, 2)
    end

    # -----------------------------------------------------------------------
    @testset "F1 structure — exact CounterTerm lists (SPEC_B1 §8.1)" begin
        @test length(cts1) == 4 + 2 + 2 + 1
        @test [ct.face for ct in cts1] ==
              [Int[], Int[], Int[], Int[], [1], [1], [2], [2], [1, 2]]
        @test [ct.subface for ct in cts1] ==
              [Int[], [1], [2], [1, 2], Int[], [2], Int[], [1], Int[]]
        @test [ct.sign for ct in cts1] == [1, -1, -1, 1, 1, -1, 1, -1, 1]

        em1 = as_ee(-1, 1)   # ε − 1
        zz  = [AS_ZERO, AS_ZERO]

        # ∅ face, four terms (WARNING-1 discriminator: the [1,2] subface
        # enters with +1 = (−1)²; wrong reading flips terms 2-4):
        @test cts1[1] == CounterTerm(Int[], Int[], 1, big(1), EpsExp[],
                  [1, 2], [:x1, :x2], [em1, em1], zz,
                  [(1 + x1 + x2, as_ee(-1, -3))])
        @test cts1[2] == CounterTerm(Int[], [1], -1, big(1), EpsExp[],
                  [1, 2], [:x1, :x2], [em1, em1], zz,
                  [(1 + x2, as_ee(-1, -3)), (1 + x1, as_ee(-1, -1))])
        @test cts1[3] == CounterTerm(Int[], [2], -1, big(1), EpsExp[],
                  [1, 2], [:x1, :x2], [em1, em1], zz,
                  [(1 + x1, as_ee(-1, -3)), (1 + x2, as_ee(-1, -1))])
        @test cts1[4] == CounterTerm(Int[], [1, 2], 1, big(1), EpsExp[],
                  [1, 2], [:x1, :x2], [em1, em1], zz,
                  [(1 + x1, as_ee(-1, -1)), (1 + x2, as_ee(-1, -1))])

        # face [1] — matches the paper's printed subtraction[[2]] VERBATIM
        # (paper:1362-1365: {x2^{−1+ε}(1+x2)^{−1−3ε}/ε,
        #                    −x2^{−1+ε}((1+x2)^{-1})^{1+ε}/ε}, vars {x2}):
        @test cts1[5] == CounterTerm([1], Int[], 1, big(1), [AS_EPS],
                  [2], [:x2], [em1], zz, [(1 + x2, as_ee(-1, -3))])
        @test cts1[6] == CounterTerm([1], [2], -1, big(1), [AS_EPS],
                  [2], [:x2], [em1], zz, [(1 + x2, as_ee(-1, -1))])

        # face [2] — mirror; NOTE its I-scan already exercises I ≠ [1]:
        # I = [1] fails (ρ2[1] = 0) ⇒ I = [2], J = [1] (WARNING-2 in vivo)
        @test cts1[7] == CounterTerm([2], Int[], 1, big(1), [AS_EPS],
                  [1], [:x1], [em1], zz, [(1 + x1, as_ee(-1, -3))])
        @test cts1[8] == CounterTerm([2], [1], -1, big(1), [AS_EPS],
                  [1], [:x1], [em1], zz, [(1 + x1, as_ee(-1, -1))])

        # face [1,2]: pure number, Vol = 1/ε² via vol_trops = [ε, ε]
        @test cts1[9] == CounterTerm([1, 2], Int[], 1, big(1),
                  [AS_EPS, AS_EPS], Int[], Symbol[], EpsExp[], zz,
                  Tuple{QQMPolyRingElem,EpsExp}[])
    end

    @testset "F1 values — closed Beta/Γ per face + Γ-identity total (G3)" begin
        for eq in (1//64, 3//128)                  # two exact rational points
            e  = AS_RR(QQ(eq))
            tot_exp   = as_gamma(e)^2 * as_gamma(1 + e) / as_gamma(1 + 3e)
            f1_exp    = (as_B(e, 1 + 2e) - as_B(e, AS_RR(1))) / e
            f12_exp   = 1 / e^2
            # ∅ face independently: t1 = Γ-form raw, t2 = t3 = −B(ε,1)B(ε,1+2ε),
            # t4 = +B(ε,1)² (all four factorize into 1-d Betas)
            fem_exp   = tot_exp - 2 * as_B(e, AS_RR(1)) * as_B(e, 1 + 2e) +
                        as_B(e, AS_RR(1))^2
            @test as_close(as_face_value(cts1, Int[],   e), fem_exp)
            @test as_close(as_face_value(cts1, [1],     e), f1_exp)
            @test as_close(as_face_value(cts1, [2],     e), f1_exp)
            @test as_close(as_face_value(cts1, [1, 2],  e), f12_exp)
            # Γ-function identity: assembled total == raw integral (G3(b))
            tot_eng = sum(as_ct_value(ct, e) for ct in cts1)
            @test as_close(tot_eng, tot_exp)
        end
    end

    @testset "F1 pole cancellation — ct pole parts kill the raw divergence" begin
        # raw ~ Γ(ε)² ~ 1/ε²; the ∅-face counterterms must cancel it so that
        # face ∅ stays BOUNDED (→ ζ₂) as ε → 0, poles living ONLY in Vol.
        e   = AS_RR(QQ(1//1024))
        raw = as_ct_value(cts1[1], e)                       # identity term
        cti = sum(as_ct_value(ct, e) for ct in cts1[2:4])   # its counterterms
        @test abs(raw) > AS_RR(10)^6                        # divergent alone
        @test abs(raw + cti) < AS_RR(5)                     # cancelled sum
        @test abs(raw + cti - AS_RR(QQ(1//6)) * const_pi(AS_RR)^2) < AS_RR(1//20)
        # Laurent ledger (SPEC_B1 §8.1): ε²·total → 1, ε·(total − 1/ε²) → 0
        tot = sum(as_ct_value(ct, e) for ct in cts1)
        @test abs(e^2 * tot - 1) < AS_RR(1//10000)
        @test abs(e * (tot - 1/e^2)) < AS_RR(1//100)
    end

    @testset "F1 WARNING-1 NEGATIVE control (G3(d): gate must be able to fail)" begin
        # flip the sign of ONE single-facet counterterm at the eval layer
        bad = copy(cts1)
        c   = cts1[2]                                # ∅-face, subface [1]
        bad[2] = CounterTerm(c.face, c.subface, +1,  # WRONG: −1 → +1
                             c.vol_det, c.vol_trops, c.J, c.face_vars,
                             c.mono, c.u_sigma, c.polys)
        e_id = AS_RR(QQ(1//64))
        tot_exp = as_gamma(e_id)^2 * as_gamma(1 + e_id) / as_gamma(1 + 3e_id)
        good_err = abs(sum(as_ct_value(ct, e_id) for ct in cts1) - tot_exp)
        bad_err  = abs(sum(as_ct_value(ct, e_id) for ct in bad)  - tot_exp)
        @test good_err < AS_RR(10)^-40
        @test bad_err  > AS_RR(1//2)                 # strictly > slack — LOUD
        @test bad_err  > AS_RR(10)^30 * good_err
        # and the local-finiteness probe blows up: ∅ face no longer bounded
        e_sm = AS_RR(QQ(1//1024))
        @test abs(as_face_value(bad,  Int[], e_sm)) > AS_RR(10)^4
        @test abs(as_face_value(cts1, Int[], e_sm)) < AS_RR(5)
    end

    # -----------------------------------------------------------------------
    @testset "F2 — WARNING-2 I/J discriminator (G4) + values" begin
        # face {2}: I-scan in lex order: I=[1] fails (det = ρ2[1] = 0) ⇒
        # I = [2], J = [1] — the completion/division/gauge subset is the
        # COMPLEMENT (WARNING-2; wl:10777-10782)
        @test _c7_first_I([2], td2.rays, 2) == [2]
        @test _c7_J_of([2], td2.rays, 2) == [1]
        @test _c7_vol_det([2], [1], td2.rays, 2) == 1
        # NEGATIVE control — the deliberately WRONG reading J = I = [2]:
        # completion row e_2 ⇒ det((0,−1),(0,1)) = 0 — wrongness DETECTED
        @test _c7_vol_det([2], [2], td2.rays, 2) == 0
        # ... and the wrong gauge (vars ∉ J_wrong = {x1} → 1) leaves the
        # restricted poly 1+x1 CONSTANT ⇒ scaleless ∫x2^{ε−1}dx2 leftover
        @test _c7_gauge_fix(1 + x1, Set([2]), 2) == R(2)
        @test is_constant(_c7_gauge_fix(1 + x1, Set([2]), 2))

        # emitted structure (SPEC_B1 §8.2 expected FaceTerms)
        @test length(cts2) == 2 + 1
        @test [ct.face for ct in cts2] == [Int[], Int[], [2]]
        f2ct = cts2[3]
        @test f2ct.J == [1] && f2ct.face_vars == [:x1]
        @test f2ct.mono == [AS_EPS]                       # x1^ε flat
        @test f2ct.polys == [(1 + x1, as_ee(-2, -3))]
        @test f2ct.vol_trops == [AS_EPS]
        @test f2ct.u_sigma == [AS_ZERO, AS_ZERO]          # π_J kills x2^{−ε}
        @test cts2[2].polys == [(1 + x1, as_ee(-2, -3)), (1 + x2, as_ee(-1, -1))]

        for eq in (1//64, 3//128)
            e = AS_RR(QQ(eq))
            tot_exp = as_gamma(e) * as_gamma(1 + e)^2 / as_gamma(2 + 3e)
            f2_exp  = as_B(1 + e, 1 + 2e) / e
            @test as_close(as_face_value(cts2, [2], e), f2_exp)
            @test as_close(as_face_value(cts2, Int[], e), tot_exp - f2_exp)
            @test as_close(sum(as_ct_value(ct, e) for ct in cts2), tot_exp)
        end
        # face ∅ = −ζ₂ε + O(ε²) (SPEC_B1 §8.2): bounded and small
        e  = AS_RR(QQ(1//1024))
        fe = as_face_value(cts2, Int[], e)
        @test abs(fe) < AS_RR(1//100)
        @test abs(fe / e + const_pi(AS_RR)^2 / 6) < AS_RR(1//20)
    end

    # -----------------------------------------------------------------------
    @testset "F3 — x^{u_σ} twist (eq 3.36), non-axis w, double-entry" begin
        @test length(cts3) == 4 + 2 + 2 + 1
        @test [ct.face for ct in cts3] ==
              [Int[], Int[], Int[], Int[], [1], [1], [3], [3], [1, 3]]

        # u_σ: faces [A]=[1] and [C]=[3] both carry x2^{−ε} (π_J-projected);
        # this is the ONLY fixture with a nontrivial twist — it dies if the
        # −TropI·w sign or the π_J projection is wrong (SPEC_B1 §8.3)
        for ct in cts3
            if ct.face == [1] || ct.face == [3]
                @test ct.u_sigma == [AS_ZERO, as_ee(0, -1)]
                @test ct.J == [2] && ct.face_vars == [:x2]
            elseif ct.face == [1, 3]
                @test ct.u_sigma == [AS_ZERO, AS_ZERO]    # all vars gauged
            end
        end

        em1 = as_ee(-1, 1)
        # face [A]: [ +x2^{ε−1}(1+x2)^{−1−5ε} ; −x2^{ε−1}(1+x2)^{−1−ε} ]
        @test cts3[5].mono == [em1] && cts3[5].polys == [(1 + x2, as_ee(-1, -5))]
        @test cts3[6].sign == -1 && cts3[6].polys == [(1 + x2, as_ee(-1, -1))]
        # double-entry check: face [C] reaches the NUMERICALLY IDENTICAL list
        # through a DIFFERENT restriction/jacobian path (1+x1x2 gauged x1→1)
        keyof(ct) = (ct.sign, ct.mono, ct.polys, ct.J, ct.face_vars,
                     ct.vol_det, ct.vol_trops, ct.u_sigma)
        @test [keyof(ct) for ct in cts3 if ct.face == [1]] ==
              [keyof(ct) for ct in cts3 if ct.face == [3]]
        # the non-axis w_A = (1,1) jacobian inside face [3], subface [1]:
        # (1 + x1x2)^{−1−ε} gauge-fixed to (1 + x2)^{−1−ε}
        @test cts3[8].subface == [1] && cts3[8].polys == [(1 + x2, as_ee(-1, -1))]

        # ∅-face term list (SPEC_B1 §8.3 display)
        @test cts3[2].polys == [(1 + x2, as_ee(-1, -5)), (1 + x1 * x2, as_ee(-1, -1))]
        @test cts3[3].polys == [(1 + x1 * x2, as_ee(-1, -5)), (1 + x2, as_ee(-1, -1))]
        @test cts3[4].polys == [(1 + x1 * x2, as_ee(-1, -1)), (1 + x2, as_ee(-1, -1))]
        @test cts3[1].mono == [em1, as_ee(-1, 2)]         # x1^{ε−1} x2^{2ε−1}

        for eq in (1//64, 3//128)
            e = AS_RR(QQ(eq))
            tot_exp = as_gamma(e)^2 * as_gamma(1 + 3e) / as_gamma(1 + 5e)
            fA_exp  = (as_B(e, 1 + 4e) - as_B(e, AS_RR(1))) / e
            fem_exp = tot_exp - 2 * as_B(e, AS_RR(1)) * as_B(e, 1 + 4e) +
                      as_B(e, AS_RR(1))^2
            @test as_close(as_face_value(cts3, [1],    e), fA_exp)
            @test as_close(as_face_value(cts3, [3],    e), fA_exp)
            @test as_close(as_face_value(cts3, [1, 3], e), 1 / e^2)
            @test as_close(as_face_value(cts3, Int[],  e), fem_exp)
            @test as_close(sum(as_ct_value(ct, e) for ct in cts3), tot_exp)
        end
        # face ∅ → ζ₂ (bounded): pole cancellation on the double-pole fixture
        e = AS_RR(QQ(1//1024))
        @test abs(as_face_value(cts3, Int[], e) -
                  const_pi(AS_RR)^2 / 6) < AS_RR(1//10)
    end

    # -----------------------------------------------------------------------
    @testset "ε-dependent jacobian exponent (G5 / RT-extra-jac)" begin
        # F1 ∅-face subface [1]: the (1+x1) jacobian must carry −(1+ε) —
        # ε-DEPENDENT EpsExp end-to-end into the C8-facing poly list
        c = cts1[2]
        jac = [pe for pe in c.polys if pe[1] == 1 + x1]
        @test length(jac) == 1
        @test jac[1][2] == as_ee(-1, -1)
        @test iseps(jac[1][2])
        # face [1] subface [2] (the paper-printed (1+x2)^{−1−ε}): same pin
        @test cts1[6].polys[1][2] == as_ee(-1, -1) && iseps(cts1[6].polys[1][2])
        # BEHAVIORAL: dropping the ε-part (exponent −1 instead of −1−ε)
        # must break the ∅-face value by a strictly-greater-than-slack margin
        e = AS_RR(QQ(1//64))
        mut = CounterTerm(c.face, c.subface, c.sign, c.vol_det, c.vol_trops,
                          c.J, c.face_vars, c.mono, c.u_sigma,
                          [(1 + x2, as_ee(-1, -3)), (1 + x1, as_ee(-1, 0))])
        @test abs(as_ct_value(mut, e) - as_ct_value(c, e)) > AS_RR(1//100)
    end

    @testset "synthetic-truth ×(1+δπ) control (G9; known pitfall: a rational δ aliases)" begin
        e = AS_RR(QQ(1//64))
        tot_exp = as_gamma(e)^2 * as_gamma(1 + e) / as_gamma(1 + 3e)
        tot     = sum(as_ct_value(ct, e) for ct in cts1)
        pert    = tot + (const_pi(AS_RR) * AS_RR(10)^-8) * as_ct_value(cts1[2], e)
        @test as_close(tot, tot_exp)                    # comparator passes clean
        @test !as_close(pert, tot_exp)                  # ... and FAILS perturbed
        @test abs(pert - tot_exp) > AS_RR(10)^-10
    end

    # -----------------------------------------------------------------------
    @testset "step 0 — locally finite early return (wl:10739-10746)" begin
        Ec = EulerIntegrand(Prefactor(), [as_ee(1, 1), as_ee(1, 1)],
                            [(1 + x1 + x2, as_ee(-3, -3))],
                            [:x1, :x2], Symbol[], R)
        for sd in (SigmaFace[], [SigmaFace(Int[], [1, 2, 3])])
            cts = subtraction_terms(Ec, td1, sd, Dict{Int,Vector{BigInt}}())
            @test length(cts) == 1
            c = cts[1]
            @test c.face == Int[] && c.subface == Int[] && c.sign == 1
            @test c.vol_det == 1 && isempty(c.vol_trops)
            @test c.J == [1, 2] && c.face_vars == [:x1, :x2]
            @test c.mono == [AS_EPS, AS_EPS]              # ν − 1: dlog → flat
            @test c.polys == Ec.polys
            # value = Γ(1+ε)²Γ(1+ε)/Γ(3+3ε)
            e = AS_RR(QQ(1//64))
            @test as_close(as_ct_value(c, e),
                           as_gamma(1 + e)^2 * as_gamma(1 + e) / as_gamma(3 + 3e))
        end
    end

    @testset "typed refusals (SPEC_B1 §6 — C7 defensive re-checks)" begin
        # power-divergent ray (a > 0): dlog ν1 = −1+ε ⇒ TropI(ρ1) = 1−ε
        Ep = EulerIntegrand(Prefactor(), [as_ee(-1, 1), AS_EPS],
                            [(1 + x1 + x2, as_ee(-1, -3))],
                            [:x1, :x2], Symbol[], R)
        @test_throws PowerDivergentRefusal subtraction_terms(Ep, td1, sd1, ws1)
        # ill-defined: TropI(ρ1) ≡ 0 (ν1 = 0) — wl:10948-10956 guard class
        Ei = EulerIntegrand(Prefactor(), [AS_ZERO, AS_EPS],
                            [(1 + x1 + x2, as_ee(-1, -3))],
                            [:x1, :x2], Symbol[], R)
        thrown = try subtraction_terms(Ei, td1, sd1, ws1); nothing
                 catch ex; ex end
        @test thrown isa B1Refusal && thrown.kind == :TropIllDefined
        @test haskey(thrown.detail, "ray")                # detail ALWAYS populated
        # w·ρ ≠ −1 (D2 pin): w1 = (2,0) gives product −2
        wsbad = Dict{Int,Vector{BigInt}}(1 => BigInt[2, 0], 2 => BigInt[0, 1])
        @test_throws WNotUnitNormalized subtraction_terms(E1, td1, sd1, wsbad)
        # missing w = the D1-corrected NotFound sentinel ⇒ GP refusal
        wsgp = Dict{Int,Vector{BigInt}}(2 => BigInt[0, 1])
        gp = try subtraction_terms(E1, td1, sd1, wsgp); nothing
             catch ex; ex end
        @test gp isa GeometricPropertyViolated
        @test gp.facet == 1 && gp.compatible == [2]
        # missing-typed value (DivergentFacet convenience path)
        dfs = [DivergentFacet(1, log_divergent, as_ee(0, -1), missing),
               DivergentFacet(2, log_divergent, as_ee(0, -1), BigInt[0, 1])]
        @test_throws GeometricPropertyViolated subtraction_terms(E1, td1, sd1, dfs)
        # refusal messages name the B2 continuation route (message law,
        # allow_continuation named)
        msg = sprint(showerror, PowerDivergentRefusal(1, as_ee(1, -1)))
        @test occursin("Nilsson-Passare", msg) && occursin("REFUSING", msg)
        @test occursin("allow_continuation", msg)
    end

    @testset "DivergentFacet convenience overload (happy path)" begin
        dfs = [DivergentFacet(1, log_divergent, as_ee(0, -1), BigInt[1, 0]),
               DivergentFacet(2, log_divergent, as_ee(0, -1), BigInt[0, 1])]
        @test subtraction_terms(E1, td1, sd1, dfs) == cts1
    end
end

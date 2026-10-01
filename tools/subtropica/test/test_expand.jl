# test_expand.jl — OWNER: unit-2. Unit tests for src/eps_expand.jl (C8).
# Self-contained (no cross-unit test deps). Run standalone:
#   ulimit -v 32505856; julia --project=tools/subtropica \
#     tools/subtropica/test/test_expand.jl
# or via test/runtests.jl (auto-included).
#
# Numeric gates in here:
#  - loggamma_series vs mpmath loggamma at >=50 dps, 3 rational offsets x 5
#    orders (RT-vii). mpmath derivatives via contour
#    integration (mp.diff method='quad', radius 1/8 — no finite-difference
#    noise), mp.dps=120.
#  - gamma_prefactor_series vs mpmath Gamma products (poles, half-integer
#    pairs, e^{gammaE*eps} normalization, eps_power shifts).
#  - eps_expand remainder RATIO gate: |S(eps)-D(eps)| scales like
#    eps^(maxorder+1) between two eps values (R2-style check), arb 600 bits.
# Symbolic atoms are evaluated with Nemo arb (independent of mpmath).

using Test
using SubTropica
using Nemo

const RBI = Rational{BigInt}

# ---------------------------------------------------------------------------
# helpers: numeric evaluation of unit-2's HlogExpr constants via Nemo arb
# (atoms emitted by eps_expand.jl only: EulerGamma, Pi, Log2, HZeta([k]),
#  HLogA(rational or "Pi"); coef strings are pure rationals)
# ---------------------------------------------------------------------------
function _parse_rat(s::AbstractString)
    if occursin("/", s)
        p, q = split(s, "/")
        return parse(BigInt, p) // parse(BigInt, q)
    end
    return parse(BigInt, s) // big(1)
end

function _atomval(RR::ArbField, a::SubTropica.HAtom)
    if a isa SubTropica.HConst
        a.name === :Pi         && return const_pi(RR)
        a.name === :EulerGamma && return const_euler(RR)
        a.name === :Log2       && return log(RR(2))
        error("test evaluator: unexpected HConst $(a.name)")
    elseif a isa SubTropica.HZeta
        length(a.idx) == 1 && a.idx[1] >= 2 || error("test evaluator: unexpected HZeta $(a.idx)")
        return zeta(RR(a.idx[1]))
    elseif a isa SubTropica.HLogA
        a.arg == "Pi" && return log(const_pi(RR))
        return log(RR(QQ(_parse_rat(a.arg))))
    end
    error("test evaluator: unexpected atom type $(typeof(a))")
end

function hval(RR::ArbField, e::SubTropica.HlogExpr)
    acc = RR(0)
    for t in e.terms
        v = RR(QQ(_parse_rat(t.coef)))
        for (a, p) in t.atoms
            v *= _atomval(RR, a)^p
        end
        acc += v
    end
    return acc
end

# mpmath reference values through a fresh python3 (ulimit prefix as everywhere in this package).
function run_python(code::String)
    io  = IOBuffer()
    cmd = pipeline(`bash -c "ulimit -v 32505856; exec python3 -"`; stdout = io)
    open(cmd, "w") do pin
        write(pin, code)
    end
    return String(take!(io))
end

const PYREF_CODE = """
from mpmath import mp, mpf, gamma, loggamma, exp, euler, factorial, diff
mp.dps = 120
r8 = mpf(1)/8
def tay(f, n):
    # contour-integration Taylor coefficients (no finite-difference noise);
    # all case singularities lie outside |t| = 1/8
    return [diff(f, 0, k, method='quad', radius=r8).real/factorial(k) for k in range(n+1)]
def pr(tag, cs):
    for i, c in enumerate(cs):
        print(tag, i, mp.nstr(c, 75))
# --- loggamma: 3 rational offsets x 5 orders (RT-vii gate) ---
pr('lg1', tay(lambda t: loggamma(1 + t), 4))
pr('lg2', tay(lambda t: loggamma(mpf(1)/2 - 2*t), 4))
pr('lg3', tay(lambda t: loggamma(mpf(7)/2 + 3*t), 4))
# --- gamma prefactor cases ---
pr('g1', tay(lambda t: mpf(2)/3*gamma(1 + t), 4))
pr('g2', tay(lambda t: t*gamma(t), 5))                    # Gamma(eps), pole -1
pr('g3', tay(lambda t: gamma(mpf(1)/2 + t)**2, 4))        # pi atom, even half-int pair
pr('g4', tay(lambda t: t*gamma(-2 + t), 5))               # Gamma(-2+eps), pole -1
pr('g5', tay(lambda t: exp(euler*t)*gamma(1 + t), 4))     # gammaE_eps normalization
pr('g6', tay(lambda t: gamma(mpf(5)/2 + t)/gamma(mpf(1)/2 - t), 4))  # half-int ratio, q_pi = 0
"""

function parse_pyref(out::String)
    ref = Dict{String,Dict{Int,BigFloat}}()
    for line in split(out, '\n')
        isempty(strip(line)) && continue
        tag, idx, val = split(strip(line))
        d = get!(ref, String(tag), Dict{Int,BigFloat}())
        d[parse(Int, idx)] = parse(BigFloat, val)
    end
    return ref
end

"floored matched digits between symbolic (arb) and mpmath reference"
function matched_digits(sym::BigFloat, ref::BigFloat)
    d = abs(sym - ref)
    iszero(d) && return 999
    scale = max(abs(ref), BigFloat(1))
    return Int(floor(-log10(d / scale)))
end

@testset "eps_expand (unit-2, C8)" begin

    RR = ArbField(600)
    pyout = run_python(PYREF_CODE)
    REF   = parse_pyref(pyout)
    @test haskey(REF, "lg1") && haskey(REF, "g6")   # subprocess sanity

    setprecision(BigFloat, 512)

    # -----------------------------------------------------------------
    # loggamma_series vs mpmath: 3 rational offsets x 5 orders, >=50 dps
    # -----------------------------------------------------------------
    @testset "loggamma_series vs mpmath (>=50 dps)" begin
        cases = [("lg1", SubTropica.EpsExp(1//1, 1//1)),
                 ("lg2", SubTropica.EpsExp(1//2, -2//1)),
                 ("lg3", SubTropica.EpsExp(7//2, 3//1))]
        for (tag, arg) in cases
            coeffs = loggamma_series(arg, 4)
            @test length(coeffs) == 5
            for k in 0:4
                sym = BigFloat(hval(RR, coeffs[k + 1]))
                @test matched_digits(sym, REF[tag][k]) >= 50
            end
        end
        # structural spot-checks: logGamma(1+eps) coeff of eps^1 is -EulerGamma,
        # of eps^3 is -zeta(3)/3; logGamma(1/2+...) constant is (1/2)Log[Pi]
        c = loggamma_series(SubTropica.EpsExp(1//1, 1//1), 3)
        @test c[1].terms == SubTropica.HTerm[]                     # logGamma(1) = 0
        @test length(c[2].terms) == 1 &&
              c[2].terms[1].coef == "-1" &&
              c[2].terms[1].atoms == [SubTropica.HConst(:EulerGamma) => 1]
        @test length(c[4].terms) == 1 &&
              c[4].terms[1].coef == "-1/3" &&
              c[4].terms[1].atoms == [SubTropica.HZeta([3]) => 1]
        ch = loggamma_series(SubTropica.EpsExp(1//2, 1//1), 0)
        @test length(ch[1].terms) == 1 &&
              ch[1].terms[1].coef == "1/2" &&
              ch[1].terms[1].atoms == [SubTropica.HLogA("Pi") => 1]
    end

    @testset "loggamma_series typed refusals" begin
        @test_throws EpsExpandRefusal loggamma_series(SubTropica.EpsExp(0//1, 1//1), 3)   # pole: caller's job
        @test_throws EpsExpandRefusal loggamma_series(SubTropica.EpsExp(-1//1, 1//1), 3)
        @test_throws EpsExpandRefusal loggamma_series(SubTropica.EpsExp(1//3, 1//1), 3)   # non-half-int base
        @test_throws EpsExpandRefusal loggamma_series(SubTropica.EpsExp(1//1, 1//1), -1)
    end

    # -----------------------------------------------------------------
    # gamma_prefactor_series vs mpmath (poles, pi atoms, normalization)
    # -----------------------------------------------------------------
    @testset "gamma_prefactor_series vs mpmath (>=50 dps)" begin
        # (tag, prefactor, expected minorder, #orders to check, ref index offset)
        # ref tag gN holds Taylor coeffs of the ANALYTIC product; for pole cases
        # the reference function is t*f(t): its coeff i == our coeff of eps^(i-1).
        gcases = [
            ("g1", SubTropica.Prefactor(2//3, 0, [(SubTropica.EpsExp(1//1, 1//1), 1)], 0//1), 0, 4),
            ("g2", SubTropica.Prefactor(1//1, 0, [(SubTropica.EpsExp(0//1, 1//1), 1)], 0//1), -1, 5),
            ("g3", SubTropica.Prefactor(1//1, 0, [(SubTropica.EpsExp(1//2, 1//1), 2)], 0//1), 0, 4),
            ("g4", SubTropica.Prefactor(1//1, 0, [(SubTropica.EpsExp(-2//1, 1//1), 1)], 0//1), -1, 5),
            ("g5", SubTropica.Prefactor(1//1, 0, [(SubTropica.EpsExp(1//1, 1//1), 1)], 1//1), 0, 4),
            ("g6", SubTropica.Prefactor(1//1, 0,
                   [(SubTropica.EpsExp(5//2, 1//1), 1), (SubTropica.EpsExp(1//2, -1//1), -1)], 0//1), 0, 4),
        ]
        for (tag, pf, expect_min, nord) in gcases
            mo, coeffs = gamma_prefactor_series(pf, nord)
            @test mo == expect_min
            @test length(coeffs) == nord + 1
            for i in 0:nord
                # our coeffs[i+1] is eps^(mo+i); ref index for pole cases is
                # taylor coeff (mo+i) - mo of t^{-mo}*f — i.e. plain i either way
                sym = BigFloat(hval(RR, coeffs[i + 1]))
                @test matched_digits(sym, REF[tag][i]) >= 50
            end
        end
        # structural: Gamma(1/2+eps)^2 leading coefficient is EXACTLY Pi
        _, c3 = gamma_prefactor_series(
            SubTropica.Prefactor(1//1, 0, [(SubTropica.EpsExp(1//2, 1//1), 2)], 0//1), 2)
        @test length(c3[1].terms) == 1 &&
              c3[1].terms[1].coef == "1" &&
              c3[1].terms[1].atoms == [SubTropica.HConst(:Pi) => 1]
        # eps_power shift: 5 * eps^2 * Gamma(1+eps) has minorder 2, leading 5
        mo5, c5 = gamma_prefactor_series(
            SubTropica.Prefactor(5//1, 2, [(SubTropica.EpsExp(1//1, 1//1), 1)], 0//1), 2)
        @test mo5 == 2
        @test length(c5[1].terms) == 1 && c5[1].terms[1].coef == "5" && isempty(c5[1].terms[1].atoms)
        # 1/Gamma(eps): ZERO of order 1 => minorder +1
        moz, _ = gamma_prefactor_series(
            SubTropica.Prefactor(1//1, 0, [(SubTropica.EpsExp(0//1, 1//1), -1)], 0//1), 2)
        @test moz == 1
    end

    @testset "gamma_prefactor_series typed refusals" begin
        # odd net half-integer Gamma count => pi^(1/2) constant
        @test_throws EpsExpandRefusal gamma_prefactor_series(
            SubTropica.Prefactor(1//1, 0, [(SubTropica.EpsExp(1//2, 1//1), 1)], 0//1), 2)
        # exact pole Gamma(0), b = 0
        @test_throws EpsExpandRefusal gamma_prefactor_series(
            SubTropica.Prefactor(1//1, 0, [(SubTropica.EpsExp(0//1, 0//1), 1)], 0//1), 2)
        # unsupported base Gamma(1/3+eps)
        @test_throws EpsExpandRefusal gamma_prefactor_series(
            SubTropica.Prefactor(1//1, 0, [(SubTropica.EpsExp(1//3, 1//1), 1)], 0//1), 2)
        # zero prefactor
        @test_throws EpsExpandRefusal gamma_prefactor_series(SubTropica.Prefactor(0//1), 2)
    end

    # -----------------------------------------------------------------
    # sign_split (RT-viii rule)
    # -----------------------------------------------------------------
    @testset "sign_split RT-viii" begin
        R2v, (x, y) = polynomial_ring(QQ, ["x", "y"])
        # tier 1: all-positive coefficients
        s = sign_split(x + y, [:x, :y], Symbol[])
        @test s.sign == 1 && s.content == 1 && s.evidence == [:coefficient_signs]
        @test s.factors == [(x + y, 1)]
        # tier 1, negative-definite: -x - y  => sign -1, positive factor stored
        s = sign_split(-x - y, [:x, :y], Symbol[])
        @test s.sign == -1 && s.factors == [(x + y, 1)]
        # content extraction: 6x + 6y
        s = sign_split(6x + 6y, [:x, :y], Symbol[])
        @test s.content == 6 && s.factors == [(x + y, 1)] && s.sign == 1
        # tier 2 probe: x^2 - x*y + y^2 positive-definite but mixed coefficient signs
        s = sign_split(x^2 - x*y + y^2, [:x, :y], Symbol[])
        @test s.sign == 1 && s.evidence == [:numeric_probe]
        @test haskey(s.provenance["probes"], string(x^2 - x*y + y^2))
        # tier 3 MANDATORY refusal: (x-y)^2 — factor (x-y) changes sign
        err = try
            sign_split(x^2 - 2x*y + y^2, [:x, :y], Symbol[])
            nothing
        catch e
            e
        end
        @test err isa SignIndefiniteError
        @test occursin("x", err.factor) && occursin("y", err.factor)   # names the factor
        # nprobe floor (R1: N >= 8)
        @test_throws EpsExpandRefusal sign_split(x + y, [:x, :y], Symbol[]; nprobe = 4)
        # domains: (1-u) indefinite on (0,inf), definite on (0,1)
        Ru, (u,) = polynomial_ring(QQ, ["u"])
        @test_throws SignIndefiniteError sign_split(1 - u, [:u], Symbol[])
        s = sign_split(1 - u, [:u], Symbol[]; domains = Dict(:u => (0//1, 1//1)))
        @test s.sign == 1 && s.factors == [(1 - u, 1)] && s.evidence == [:numeric_probe]
        # invalid domain refuses
        @test_throws EpsExpandRefusal sign_split(1 - u, [:u], Symbol[];
                                                 domains = Dict(:u => (-1//1, 1//1)))
    end

    # -----------------------------------------------------------------
    # Laurent pole guard (R2)
    # -----------------------------------------------------------------
    @testset "laurent_pole_guard" begin
        @test laurent_pole_guard(-2, 2, 4) === nothing
        @test_throws EpsExpandRefusal laurent_pole_guard(-2, 2, 3)   # pole eats top order
        @test laurent_pole_guard(1, 2, 2) === nothing                # zeros need no padding
    end

    # -----------------------------------------------------------------
    # MANDATORY fixture (RT-extra-jac): counterterm-style jacobian factor
    # (1-u)^(1-TropI) with eps-dependent exponent; TropI = -2 + 3*eps
    # => exponent 3 - 3*eps. Domain u in (0,1) (farm counterterm variable).
    # -----------------------------------------------------------------
    @testset "MANDATORY jacobian fixture (1-u)^(1-TropI), eps-dependent" begin
        Ru, (u,) = polynomial_ring(QQ, ["u"])
        quad = SubTropica.EulerIntegrand(
            SubTropica.Prefactor(1//1),
            [SubTropica.EpsExp(0//1, 0//1)],
            [(1 - u, SubTropica.EpsExp(3//1, -3//1))],       # (1-u)^(3-3eps)
            [:u], Symbol[], Ru)
        res = eps_expand(quad, 2; domains = Dict(:u => (0//1, 1//1)))
        @test res.minorder == 0 && res.maxorder == 2 && res.prefactor_minorder == 0
        @test length(res.integrands) == 3
        # eps^0: (1-u)^3, no logs
        @test length(res.integrands[1]) == 1
        li0 = res.integrands[1][1]
        @test li0.num == (1 - u)^3 && li0.den == one(Ru) && isempty(li0.logs)
        @test li0.vars == [:u]
        # eps^1: -3*(1-u)^3 * Log[1-u]
        @test length(res.integrands[2]) == 1
        li1 = res.integrands[2][1]
        @test li1.num == -3 * (1 - u)^3 && li1.logs == [(1 - u) => 1]
        # eps^2: (9/2)*(1-u)^3 * Log[1-u]^2
        @test length(res.integrands[3]) == 1
        li2 = res.integrands[3][1]
        @test li2.num == Ru(QQ(9//2)) * (1 - u)^3 && li2.logs == [(1 - u) => 2]
        # numeric flow-through: S(eps) matches (1-u)^(3-3eps) at u=1/3,
        # remainder RATIO gate ~ eps^(maxorder+1)
        pt   = QQ(1//3)
        base = RR(QQ(2//3))                                # (1-u) at u=1/3
        function S(epsv::ArbFieldElem)
            acc = RR(0)
            for (o, lst) in enumerate(res.integrands)
                for li in lst
                    v = RR(evaluate(li.num, [pt])) / RR(evaluate(li.den, [pt]))
                    for (g, k) in li.logs
                        v *= log(RR(evaluate(g, [pt])))^k
                    end
                    acc += v * epsv^(o - 1)
                end
            end
            return acc
        end
        D(epsv::ArbFieldElem) = exp((RR(3) - 3epsv) * log(base))
        e1, e2 = RR(QQ(1//10^6)), RR(QQ(1//10^8))
        r1, r2 = abs(S(e1) - D(e1)), abs(S(e2) - D(e2))
        @test BigFloat(r1) < 1e-15                       # remainder ~ eps^3 * O(1)
        # ratio gate: r2/r1 approx (e2/e1)^(maxo+1) = 1e-6 within a factor 10
        @test BigFloat(r2) < 10 * BigFloat(r1) * BigFloat(1e-6)
        @test BigFloat(r2) > BigFloat(r1) * BigFloat(1e-7) / 10
    end

    # -----------------------------------------------------------------
    # eps_expand: factorization + letter merge + prefactor pole padding
    # integrand: x^(2+eps/2) * y * (x^2*y + x*y^2)^(-1+2eps) * (x+s)^(1-eps)
    # prefactor: Gamma(eps) (pole => maxorder = order + 1)
    # -----------------------------------------------------------------
    @testset "eps_expand: factored polys, letter merge, pole padding" begin
        R3, (x, y, s) = polynomial_ring(QQ, ["x", "y", "s"])
        quad = SubTropica.EulerIntegrand(
            SubTropica.Prefactor(1//1, 0, [(SubTropica.EpsExp(0//1, 1//1), 1)], 0//1),  # Gamma(eps)
            [SubTropica.EpsExp(2//1, 1//2), SubTropica.EpsExp(1//1, 0//1)],
            [(x^2*y + x*y^2, SubTropica.EpsExp(-1//1, 2//1)),
             (x + s,         SubTropica.EpsExp(1//1, -1//1))],
            [:x, :y], [:s], R3)
        res = eps_expand(quad, 2)
        @test res.prefactor_minorder == -1
        @test res.maxorder == 3                       # order + pole padding (wl:11186-11187)
        @test length(res.integrands) == 4
        @test length(res.prefactor_coeffs) == 5       # covers eps^-1 .. eps^3
        # letters: x (1/2 + 2 = 5/2 after merge with factor x of P1), y (2), x+y (2), x+s (-1)
        # eps^0: single term, A = x^2*y * (x+s) / (x*y*(x+y))
        @test length(res.integrands[1]) == 1
        li0 = res.integrands[1][1]
        @test li0.num * (x*y*(x + y)) == li0.den * (x^2*y*(x + s))   # A as a fraction
        # eps^1: 4 letters => 4 terms; eps^2: C(2+3,3) = 10 terms
        @test length(res.integrands[2]) == 4
        @test length(res.integrands[3]) == 10
        # letter merge check: exactly one eps^1 term carries Log[x], with
        # coefficient 5/2 = 1/2 (from nu_x) + 2 (from the x factor of P1)
        logx = [li for li in res.integrands[2] if li.logs == [x => 1]]
        @test length(logx) == 1
        @test logx[1].den == li0.den
        @test logx[1].num == li0.num * QQ(5//2)
        # numeric remainder RATIO gate at (x,y,s) = (3/7, 2/5, 1/3)
        ptv  = [QQ(3//7), QQ(2//5), QQ(1//3)]
        function S(epsv::ArbFieldElem)
            acc = RR(0)
            for (o, lst) in enumerate(res.integrands)
                for li in lst
                    v = RR(evaluate(li.num, ptv)) / RR(evaluate(li.den, ptv))
                    for (g, k) in li.logs
                        v *= log(RR(evaluate(g, ptv)))^k
                    end
                    acc += v * epsv^(o - 1)
                end
            end
            return acc
        end
        function D(epsv::ArbFieldElem)
            xv, yv, sv = RR(QQ(3//7)), RR(QQ(2//5)), RR(QQ(1//3))
            return exp((RR(2) + epsv/2)*log(xv)) * yv *
                   exp((RR(-1) + 2epsv)*log(xv^2*yv + xv*yv^2)) *
                   exp((RR(1) - epsv)*log(xv + sv))
        end
        e1, e2 = RR(QQ(1//10^6)), RR(QQ(1//10^8))
        r1, r2 = abs(S(e1) - D(e1)), abs(S(e2) - D(e2))
        @test BigFloat(r1) < 1e-20                    # remainder ~ eps^4
        @test BigFloat(r2) < 10 * BigFloat(r1) * BigFloat(1e-8)   # (e2/e1)^4
    end

    # -----------------------------------------------------------------
    # eps_expand: padding, b=0 integer-power path, refusals
    # -----------------------------------------------------------------
    @testset "eps_expand: padding + b=0 path + refusals" begin
        R2v, (x, y) = polynomial_ring(QQ, ["x", "y"])
        # all-b=0 quad: higher orders are EXPLICIT empty lists (padFromSeries analogue)
        quad0 = SubTropica.EulerIntegrand(
            SubTropica.Prefactor(1//1), [SubTropica.EpsExp(1//1), SubTropica.EpsExp(0//1)],
            [(x + y, SubTropica.EpsExp(-2//1, 0//1))], [:x, :y], Symbol[], R2v)
        res0 = eps_expand(quad0, 2)
        @test length(res0.integrands) == 3
        @test length(res0.integrands[1]) == 1
        @test isempty(res0.integrands[2]) && isempty(res0.integrands[3])
        @test res0.integrands[1][1].num == x && res0.integrands[1][1].den == (x + y)^2
        # b=0 with an indefinite poly: NO refusal (integer power, no split needed)
        quadsq = SubTropica.EulerIntegrand(
            SubTropica.Prefactor(1//1), [SubTropica.EpsExp(0//1), SubTropica.EpsExp(0//1)],
            [(x^2 - 2x*y + y^2, SubTropica.EpsExp(2//1, 0//1))], [:x, :y], Symbol[], R2v)
        ressq = eps_expand(quadsq, 0)
        @test ressq.integrands[1][1].num == (x^2 - 2x*y + y^2)^2
        # b != 0 with the same indefinite poly: MUST refuse (RT-viii tier 3)
        quadbad = SubTropica.EulerIntegrand(
            SubTropica.Prefactor(1//1), [SubTropica.EpsExp(0//1), SubTropica.EpsExp(0//1)],
            [(x^2 - 2x*y + y^2, SubTropica.EpsExp(2//1, 1//1))], [:x, :y], Symbol[], R2v)
        @test_throws SignIndefiniteError eps_expand(quadbad, 1)
        # negative-on-domain poly with eps exponent: Euclidean refusal
        quadneg = SubTropica.EulerIntegrand(
            SubTropica.Prefactor(1//1), [SubTropica.EpsExp(0//1), SubTropica.EpsExp(0//1)],
            [(-x - y, SubTropica.EpsExp(1//1, 1//1))], [:x, :y], Symbol[], R2v)
        @test_throws EpsExpandRefusal eps_expand(quadneg, 1)
        # non-integer constant exponent part: typed refusal
        quadhalf = SubTropica.EulerIntegrand(
            SubTropica.Prefactor(1//1), [SubTropica.EpsExp(0//1), SubTropica.EpsExp(0//1)],
            [(x + y, SubTropica.EpsExp(1//2, 1//1))], [:x, :y], Symbol[], R2v)
        @test_throws EpsExpandRefusal eps_expand(quadhalf, 1)
        # nu length mismatch
        quadlen = SubTropica.EulerIntegrand(
            SubTropica.Prefactor(1//1), [SubTropica.EpsExp(0//1)],
            [(x + y, SubTropica.EpsExp(1//1, 0//1))], [:x, :y], Symbol[], R2v)
        @test_throws EpsExpandRefusal eps_expand(quadlen, 1)
    end

end

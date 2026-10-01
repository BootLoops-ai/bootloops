# test_graph.jl — OWNER: unit-G. Self-contained tests for src/graph.jl
# (Phase-C graph front-end). No cross-unit surfaces are called.
#
# Run standalone (DESIGN_B1 §3 — DEFAULT global Julia env, include-by-path;
# NOT --project=tools/subtropica):
#   ulimit -v 32505856; julia tools/subtropica/test/test_graph.jl
# Auto-joins the shared suite via test/runtests.jl at integration (no edits).
#
# What is tested (task charter):
#  * U,F EXACT match vs pySecDec LoopIntegralFromGraph on 5 graphs
#    (bubble, 2-mass sunrise, massless 1-loop box, the moller_box_eq41
#    fixture graph, the 3-loop 8-propagator lbl3se fixture graph) via a
#    small python subprocess helper — the live-tested glue pattern of
#    scripts/make_fixtures.py (pySecDec 1.6.6, probed live
#    the lbl3se graph reconstruction below was verified against
#    an independent momentum-space LoopIntegralFromPropagators
#    definition — the helper re-asserts
#    that equality on every run).
#  * graph_to_quadruple == the FROZEN moller_box_eq41 fixture (byte-level);
#    the synthetic 7-var fixture (synth_prod7_7var) is reconstructed
#    independently in-test and must match exactly.
#  * Gamma prefactor per paper eq 1.4 (subtropica.txt:236-249) incl. the
#    (−1)^{ν_tot} sign and the Γ(ν_e) denominators.
#  * Lee-Pomeransky variant flag (G = U + F, no gauge).
#  * Pinch (ν_e = 0) semantics; typed refusals (numerators, non-Euclidean,
#    scaleless, disconnected, tree-level, missing dot products).
#  * Fixture-loader round-trip: load -> emit (subtropica-euler-quad-v1) ->
#    reload over ALL six frozen fixture files; byte-parity of poly strings
#    for the make_fixtures.py-generated ones.

using Test, Nemo, JSON

include(joinpath(@__DIR__, "..", "src", "types.jl"))
include(joinpath(@__DIR__, "..", "src", "b1_types.jl"))
include(joinpath(@__DIR__, "..", "src", "graph.jl"))

const _TG_FIXDIR = normpath(joinpath(@__DIR__, "..", "fixtures"))
const _TG_SCRATCH = mktempdir(; prefix="subtropica_graph_test_")

# ---------------------------------------------------------------------------
# pySecDec subprocess helper (make_fixtures.py live-tested glue pattern).
# Written to scratch at test time — unit-G owns only src/graph.jl +
# test/test_graph.jl, so the helper is embedded here, not a repo file.
# Every invocation runs under the house memory cap (ulimit -v 32505856).
# ---------------------------------------------------------------------------
const _TG_PYHELPER = raw"""
#!/usr/bin/env python3
# pySecDec U/F helper for subtropica test_graph.jl (OWNER: unit-G; embedded in
# the test file, materialized to scratch at run time).
# argv[1]: JSON spec {internal_lines, external_lines, replacement_rules,
#   powerlist?, euclid_subs?, allvars, propagators?/loop_momenta?/
#   external_momenta? (optional LoopIntegralFromPropagators cross-check)}.
# stdout: JSON {U:[[coef,[exps...]],...], F:..., expU:[a,b], expF:[a,b],
#   gamma_arg:[a,b], powerlist, L} — coefs "p" or "p/q", exps over allvars.
import json, sys
import sympy as sp
from fractions import Fraction
from pySecDec.loop_integral import LoopIntegralFromGraph, LoopIntegralFromPropagators

spec = json.load(open(sys.argv[1]))
kwargs = {}
if spec.get("powerlist"):
    kwargs["powerlist"] = spec["powerlist"]
li = LoopIntegralFromGraph(
    internal_lines=[[m, list(vw)] for (m, vw) in spec["internal_lines"]],
    external_lines=[[p, v] for (p, v) in spec["external_lines"]],
    replacement_rules=[tuple(r) for r in spec["replacement_rules"]],
    **kwargs)

allvars = [sp.Symbol(v, positive=True) for v in spec["allvars"]]
loc = {str(v): v for v in allvars}
eps = sp.Symbol("eps")
loc["eps"] = eps

def euclid(expr):
    e = sp.expand(expr)
    for lhs, rhs in spec.get("euclid_subs", []):
        e = sp.expand(e.subs(sp.sympify(lhs, locals=dict(loc)),
                             sp.sympify(rhs, locals=dict(loc))))
    return e

U = euclid(sp.sympify(str(li.U), locals=dict(loc)))
F = euclid(sp.sympify(str(li.F), locals=dict(loc)))

# optional cross-check: the momentum-space build must agree with the graph
# (used for lbl3se — pins the graph reconstruction against its validated
# momentum-space propagator list)
if spec.get("propagators"):
    li2 = LoopIntegralFromPropagators(
        propagators=spec["propagators"],
        loop_momenta=spec["loop_momenta"],
        external_momenta=spec.get("external_momenta", []),
        replacement_rules=[tuple(r) for r in spec["replacement_rules"]])
    assert sp.expand(U - euclid(sp.sympify(str(li2.U), locals=dict(loc)))) == 0, \
        "graph-vs-propagator U mismatch"
    assert sp.expand(F - euclid(sp.sympify(str(li2.F), locals=dict(loc)))) == 0, \
        "graph-vs-propagator F mismatch"

def eps_pair(expr):
    p = sp.Poly(sp.expand(sp.sympify(str(expr), locals=dict(loc))), eps)
    assert p.degree() <= 1, f"exponent not eps-linear: {expr}"
    return [str(Fraction(str(p.coeff_monomial(1)))),
            str(Fraction(str(p.coeff_monomial(eps))))]

def terms(poly):
    p = sp.Poly(poly, *allvars)
    out = []
    for mono, coeff in sorted(p.terms(), key=lambda t: t[0], reverse=True):
        c = sp.Rational(coeff)
        out.append([f"{c.p}/{c.q}" if c.q != 1 else str(c.p), list(mono)])
    return out

# Gamma_factor = (+/-) gamma(arg): pySecDec bakes the (-1)^{nu_tot} of the
# Feynman parametrization INTO Gamma_factor (observed live: sunrise nu_tot=3
# gives "-gamma(2*eps - 1)") — split it off so the sign itself becomes an
# external oracle for the eq 1.4 prefactor sign c = (-1)^{nu_tot}.
gf = sp.sympify(str(li.Gamma_factor), locals=dict(loc, gamma=sp.Function("gamma")))
gsign, gf = gf.as_coeff_Mul()
assert gsign in (1, -1), str(li.Gamma_factor)
assert gf.func.__name__ == "gamma" and len(gf.args) == 1, str(li.Gamma_factor)

json.dump({"U": terms(U), "F": terms(F),
           "expU": eps_pair(li.exponent_U), "expF": eps_pair(li.exponent_F),
           "gamma_arg": eps_pair(gf.args[0]), "gamma_sign": int(gsign),
           "powerlist": [int(p) for p in li.powerlist], "L": int(li.L)},
          sys.stdout)
print()
"""

function _tg_pysecdec_available()
    return success(`bash -c "ulimit -v 32505856; python3 -c 'import pySecDec, sympy' 2>/dev/null"`)
end

let helper = joinpath(_TG_SCRATCH, "pysecdec_uf_helper.py")
    write(helper, _TG_PYHELPER)
end

function _tg_run_helper(spec::AbstractDict, tag::AbstractString)
    helper = joinpath(_TG_SCRATCH, "pysecdec_uf_helper.py")
    sfile = joinpath(_TG_SCRATCH, "spec_$(tag).json")
    ofile = joinpath(_TG_SCRATCH, "out_$(tag).json")
    write(sfile, JSON.json(spec))
    run(`bash -c "ulimit -v 32505856; python3 $helper $sfile > $ofile"`)
    return JSON.parsefile(ofile)
end

function _tg_poly_from_terms(R::QQMPolyRing, ts)
    B = MPolyBuildCtx(R)
    for t in ts
        push_term!(B, QQ(_graph_qparse(t[1])), Int[Int(e) for e in t[2]])
    end
    return finish(B)
end

_tg_eps(pair) = EpsExp(_graph_qparse(pair[1]), _graph_qparse(pair[2]))

_tg_kind(f) = try
    f()
    return :none
catch e
    e isa GraphRefusal && return e.kind
    rethrow()
end

_tg_pref_eq(a::Prefactor, b::Prefactor) =
    a.c == b.c && a.eps_power == b.eps_power && a.gammas == b.gammas &&
    a.gammaE_eps == b.gammaE_eps

# parent-independent poly comparison via the canonical fixture string
_tg_poly_eq(a::QQMPolyRingElem, b::QQMPolyRingElem) =
    _graph_poly_string(a) == _graph_poly_string(b)

function _tg_E_eq(a::EulerIntegrand, b::EulerIntegrand)
    a.vars == b.vars || return false
    a.kinvars == b.kinvars || return false
    a.nu == b.nu || return false
    length(a.polys) == length(b.polys) || return false
    for i in eachindex(a.polys)
        _tg_poly_eq(a.polys[i][1], b.polys[i][1]) || return false
        a.polys[i][2] == b.polys[i][2] || return false
    end
    return _tg_pref_eq(a.prefactor, b.prefactor)
end

# ---------------------------------------------------------------------------
# The 5 parity graphs: native spec + pySecDec spec.
# Dot-product tables use legs 1..nlegs−1 (leg nlegs eliminated by momentum
# conservation — src/graph.jl _graph_psq). Euclidean charts follow
# scripts/make_fixtures.py (s = −sb, t = −tb, p² = −psb; sb,tb,psb,msq > 0).
# ---------------------------------------------------------------------------

# (1) massive bubble: edges 2×(1,2), masses (msq, 0), p² = −psb
const _TG_BUB = (
    edges = [(1, 2), (1, 2)], nodes = [1, 2],
    masses = Any[:msq, 0],
    kin = Dict((1, 1) => [:psb => -1]),
    kinvars = Symbol[:msq, :psb],
    spec = Dict{String,Any}(
        "internal_lines" => Any[Any["m", [1, 2]], Any["0", [1, 2]]],
        "external_lines" => Any[Any["p1", 1], Any["p2", 2]],
        "replacement_rules" => Any[Any["p1*p1", "-psb"], Any["p2*p2", "-psb"],
                                   Any["p1*p2", "psb"], Any["m*m", "msq"]],
        "allvars" => Any["x0", "x1", "msq", "psb"]))

# (2) 2-mass sunrise: edges 3×(1,2), masses (m1sq, m2sq, 0), p² = −psb
const _TG_SUN = (
    edges = [(1, 2), (1, 2), (1, 2)], nodes = [1, 2],
    masses = Any[:m1sq, :m2sq, 0],
    kin = Dict((1, 1) => [:psb => -1]),
    kinvars = Symbol[:m1sq, :m2sq, :psb],
    spec = Dict{String,Any}(
        "internal_lines" => Any[Any["m1", [1, 2]], Any["m2", [1, 2]],
                                Any["0", [1, 2]]],
        "external_lines" => Any[Any["p1", 1], Any["p2", 2]],
        "replacement_rules" => Any[Any["p1*p1", "-psb"], Any["p2*p2", "-psb"],
                                   Any["p1*p2", "psb"],
                                   Any["m1*m1", "m1sq"], Any["m2*m2", "m2sq"]],
        "allvars" => Any["x0", "x1", "x2", "m1sq", "m2sq", "psb"]))

# (3) massless 1-loop box: cycle (1,2)(2,3)(3,4)(4,1); s=−sb, t=−tb
const _TG_BOX = (
    edges = [(1, 2), (2, 3), (3, 4), (4, 1)], nodes = [1, 2, 3, 4],
    masses = Any[0, 0, 0, 0],
    kin = Dict((1, 1) => 0, (2, 2) => 0, (3, 3) => 0,
               (1, 2) => [:sb => -1 // 2], (2, 3) => [:tb => -1 // 2],
               (1, 3) => [:sb => 1 // 2, :tb => 1 // 2]),
    kinvars = Symbol[:sb, :tb],
    spec = Dict{String,Any}(
        "internal_lines" => Any[Any["0", [1, 2]], Any["0", [2, 3]],
                                Any["0", [3, 4]], Any["0", [4, 1]]],
        "external_lines" => Any[Any["p1", 1], Any["p2", 2],
                                Any["p3", 3], Any["p4", 4]],
        "replacement_rules" => Any[
            Any["p1*p1", "0"], Any["p2*p2", "0"], Any["p3*p3", "0"],
            Any["p4*p4", "0"],
            Any["p1*p2", "s/2"], Any["p3*p4", "s/2"],
            Any["p2*p3", "t/2"], Any["p1*p4", "t/2"],
            Any["p1*p3", "-s/2-t/2"], Any["p2*p4", "-s/2-t/2"]],
        "euclid_subs" => Any[Any["s", "-sb"], Any["t", "-tb"]],
        "allvars" => Any["x0", "x1", "x2", "x3", "sb", "tb"]))

# (4) moller box — the moller_box_eq41 fixture graph (make_fixtures.py
# build_moller_box graph + Mandelstams, Euclidean chart s=−sb, t=−tb;
# u = 4msq − s − t ⇒ p1·p3 = msq + sb/2 + tb/2)
const _TG_MOL = (
    edges = [(1, 2), (2, 3), (3, 4), (4, 1)], nodes = [1, 2, 3, 4],
    masses = Any[:msq, 0, :msq, 0],
    kin = Dict((1, 1) => :msq, (2, 2) => :msq, (3, 3) => :msq,
               (1, 2) => [:sb => -1 // 2, :msq => -1],
               (2, 3) => [:tb => -1 // 2, :msq => -1],
               (1, 3) => [:msq => 1, :sb => 1 // 2, :tb => 1 // 2]),
    kinvars = Symbol[:sb, :tb, :msq],          # fixture kinvar ORDER (pinned)
    spec = Dict{String,Any}(
        "internal_lines" => Any[Any["m", [1, 2]], Any["0", [2, 3]],
                                Any["m", [3, 4]], Any["0", [4, 1]]],
        "external_lines" => Any[Any["p1", 1], Any["p2", 2],
                                Any["p3", 3], Any["p4", 4]],
        "replacement_rules" => Any[
            Any["p1*p1", "msq"], Any["p2*p2", "msq"], Any["p3*p3", "msq"],
            Any["p4*p4", "msq"],
            Any["p1*p2", "s/2-msq"], Any["p3*p4", "s/2-msq"],
            Any["p2*p3", "t/2-msq"], Any["p1*p4", "t/2-msq"],
            Any["p1*p3", "msq-s/2-t/2"], Any["p2*p4", "msq-s/2-t/2"],
            Any["m*m", "msq"]],
        "euclid_subs" => Any[Any["s", "-sb"], Any["t", "-tb"]],
        "allvars" => Any["x0", "x1", "x2", "x3", "sb", "tb", "msq"]))

# (5) lbl3se — 3-loop 8-propagator light-by-light self-energy. Graph
# reconstruction of the momentum-space definition (its validated
# PROPAGATORS list): vertices 1..6, ell-loop 1→2→3→4 with photons p1..p4, the
# k1/k2 subgraph on vertices 4,5,6:
#   e1 (1,2) m [ell+p1]     e2 (2,3) m [ell+p1+p2]  e3 (3,4) m [ell−p4]
#   e4 (4,6) m [k1]         e5 (5,1) m [k2]         e6 (5,6) m [ell−k1−k2]
#   e7 (4,5) 0 [ell−k1]     e8 (6,1) 0 [ell−k2]
# The helper re-verifies this reconstruction against the propagator list on
# every run (LoopIntegralFromPropagators cross-assert). Kinematics kept
# SYMBOLIC (s,t,m2 — not Euclidean-positive), so this entry exercises
# symanzik_UF only; the quadruple test uses the physical point below.
const _TG_LBL_EDGES = [(1, 2), (2, 3), (3, 4), (4, 6), (5, 1), (5, 6),
                       (4, 5), (6, 1)]
const _TG_LBL = (
    edges = _TG_LBL_EDGES, nodes = [1, 2, 3, 4],
    masses = Any[:m2, :m2, :m2, :m2, :m2, :m2, 0, 0],
    kin = Dict((1, 1) => 0, (2, 2) => 0, (3, 3) => 0,
               (1, 2) => [:s => 1 // 2], (1, 3) => [:t => 1 // 2],
               (2, 3) => [:s => -1 // 2, :t => -1 // 2]),
    kinvars = Symbol[:s, :t, :m2],
    spec = Dict{String,Any}(
        "internal_lines" => Any[Any["m", [1, 2]], Any["m", [2, 3]],
                                Any["m", [3, 4]], Any["m", [4, 6]],
                                Any["m", [5, 1]], Any["m", [5, 6]],
                                Any["0", [4, 5]], Any["0", [6, 1]]],
        "external_lines" => Any[Any["p1", 1], Any["p2", 2],
                                Any["p3", 3], Any["p4", 4]],
        "replacement_rules" => Any[
            Any["p1*p1", "0"], Any["p2*p2", "0"], Any["p3*p3", "0"],
            Any["p4*p4", "0"],
            Any["p1*p2", "s/2"], Any["p1*p3", "t/2"],
            Any["p2*p3", "(-s-t)/2"], Any["p1*p4", "(-s-t)/2"],
            Any["p2*p4", "t/2"], Any["p3*p4", "s/2"], Any["m*m", "m2"]],
        "allvars" => Any["x0", "x1", "x2", "x3", "x4", "x5", "x6", "x7",
                         "s", "t", "m2"],
        # momentum-space cross-check (extract_uf.py definition, verbatim)
        "propagators" => Any["(ell+p1)**2 - m2", "(ell+p1+p2)**2 - m2",
                             "(ell-p4)**2 - m2", "k1**2 - m2", "k2**2 - m2",
                             "(ell-k1-k2)**2 - m2", "(ell-k1)**2",
                             "(ell-k2)**2"],
        "loop_momenta" => Any["ell", "k1", "k2"],
        "external_momenta" => Any["p1", "p2", "p3", "p4"]))

const _TG_PARITY = [("bubble", _TG_BUB), ("sunrise2m", _TG_SUN),
                    ("box1l", _TG_BOX), ("moller", _TG_MOL),
                    ("lbl3se", _TG_LBL)]

# === testsets appended below ===

@testset "graph.jl (unit-G)" begin

    # -----------------------------------------------------------------------
    # Hand-pinned Symanzik polynomials (no external deps; textbook values)
    # -----------------------------------------------------------------------
    @testset "hand-pinned U,F" begin
        # massive bubble: U = x0+x1, F = psb·x0·x1 + msq·x0·(x0+x1)
        UF = symanzik_UF(_TG_BUB.edges, _TG_BUB.nodes, _TG_BUB.masses,
                         _TG_BUB.kin; kinvars=_TG_BUB.kinvars)
        x0, x1 = UF.xgens
        msq, psb = UF.kinmap[:msq], UF.kinmap[:psb]
        @test UF.L == 1 && UF.nv == 2
        @test UF.U == x0 + x1
        @test UF.F == psb * x0 * x1 + msq * x0 * (x0 + x1)

        # 2-mass sunrise: U = x0x1+x0x2+x1x2,
        # F = psb·x0x1x2 + (m1sq·x0 + m2sq·x1)·U
        UF = symanzik_UF(_TG_SUN.edges, _TG_SUN.nodes, _TG_SUN.masses,
                         _TG_SUN.kin; kinvars=_TG_SUN.kinvars)
        x0, x1, x2 = UF.xgens
        m1sq, m2sq, psb = UF.kinmap[:m1sq], UF.kinmap[:m2sq], UF.kinmap[:psb]
        @test UF.L == 2
        @test UF.U == x0 * x1 + x0 * x2 + x1 * x2
        @test UF.F == psb * x0 * x1 * x2 + (m1sq * x0 + m2sq * x1) * UF.U

        # massless box: U = Σx, F = sb·x1x3 + tb·x0x2
        UF = symanzik_UF(_TG_BOX.edges, _TG_BOX.nodes, _TG_BOX.masses,
                         _TG_BOX.kin; kinvars=_TG_BOX.kinvars)
        x0, x1, x2, x3 = UF.xgens
        sb, tb = UF.kinmap[:sb], UF.kinmap[:tb]
        @test UF.L == 1
        @test UF.U == x0 + x1 + x2 + x3
        @test UF.F == sb * x1 * x3 + tb * x0 * x2

        # 1-loop massive tadpole (self-loop edge): U = x0, F = msq·x0²
        UF = symanzik_UF([(1, 1)], Int[], Any[:msq], Dict();
                         kinvars=Symbol[:msq])
        @test UF.L == 1 && UF.U == UF.xgens[1]
        @test UF.F == UF.kinmap[:msq] * UF.xgens[1]^2
    end

    # -----------------------------------------------------------------------
    # Typed refusals (refuse LOUDLY, typed — B1 refusal style)
    # -----------------------------------------------------------------------
    @testset "typed refusals" begin
        gq(args...; kw...) = graph_to_quadruple(args...; kw...)
        # numerators: OPTIONAL this pass ⇒ typed refusal (charter;
        # STSymanzik CenterDot machinery wl:3546-3568/3606-3660 not ported)
        @test _tg_kind(() -> gq(_TG_BUB.edges, _TG_BUB.nodes, _TG_BUB.masses,
                _TG_BUB.kin; kinvars=_TG_BUB.kinvars,
                numerators=["l1*l1"])) == :NumeratorsNotSupported
        # negative exponent = numerator territory
        @test _tg_kind(() -> gq(_TG_BUB.edges, _TG_BUB.nodes, _TG_BUB.masses,
                _TG_BUB.kin; kinvars=_TG_BUB.kinvars,
                nu=[1, -1])) == :NumeratorExponent
        # disconnected graph (STSymanzik::zeroU analogue, wl:3101)
        @test _tg_kind(() -> symanzik_UF([(1, 2), (3, 4)], [1, 3],
                Any[:msq, :msq], Dict((1, 1) => [:psb => -1]);
                kinvars=Symbol[:msq, :psb])) == :Disconnected
        # tree-level graph
        @test _tg_kind(() -> symanzik_UF([(1, 2)], [1, 2], Any[:msq],
                Dict((1, 1) => [:psb => -1]);
                kinvars=Symbol[:msq, :psb])) == :NoLoop
        # scaleless: massless bubble at p² = 0 ⇒ F ≡ 0
        # (STSymanzik::zeroF wl:3527; vacuumPeriod shortcut NOT ported)
        @test _tg_kind(() -> symanzik_UF([(1, 2), (1, 2)], [1, 2], Any[0, 0],
                Dict((1, 1) => 0))) == :Scaleless
        # non-Euclidean chart: p² = +psb flips the F sign ⇒ typed refusal
        # at the quadruple layer (types.jl declared-Euclidean contract)
        @test _tg_kind(() -> gq([(1, 2), (1, 2)], [1, 2], Any[0, 0],
                Dict((1, 1) => [:psb => 1]);
                kinvars=Symbol[:psb])) == :NonEuclidean
        # missing dot product (needed pair not supplied)
        @test _tg_kind(() -> symanzik_UF(_TG_BOX.edges, _TG_BOX.nodes,
                _TG_BOX.masses,
                Dict((1, 1) => 0, (2, 2) => 0, (3, 3) => 0,
                     (1, 2) => [:sb => -1 // 2], (2, 3) => [:tb => -1 // 2]);
                kinvars=Symbol[:sb, :tb])) == :MissingDotProduct
        # dot product involving the conservation-eliminated last leg
        @test _tg_kind(() -> symanzik_UF(_TG_BUB.edges, _TG_BUB.nodes,
                _TG_BUB.masses, Dict((1, 2) => 0);
                kinvars=Symbol[:msq])) == :BadInput
        # Lee-Pomeransky is projectively complete: gauge must be nothing
        @test _tg_kind(() -> gq(_TG_MOL.edges, _TG_MOL.nodes, _TG_MOL.masses,
                _TG_MOL.kin; kinvars=_TG_MOL.kinvars,
                lee_pomeransky=true, gauge=1)) == :BadInput
        # gauge index must be a kept edge
        @test _tg_kind(() -> gq(_TG_MOL.edges, _TG_MOL.nodes, _TG_MOL.masses,
                _TG_MOL.kin; kinvars=_TG_MOL.kinvars, gauge=7)) == :BadInput
    end

    # -----------------------------------------------------------------------
    # pySecDec parity — U,F EXACT on the 5 charter graphs
    # -----------------------------------------------------------------------
    if _tg_pysecdec_available()
        @testset "pySecDec LoopIntegralFromGraph parity (5 graphs)" begin
            for (tag, G) in _TG_PARITY
                UF = symanzik_UF(G.edges, G.nodes, G.masses, G.kin;
                                 kinvars=G.kinvars)
                out = _tg_run_helper(G.spec, tag)
                @test out["L"] == UF.L
                # EXACT polynomial equality in OUR ring (helper terms are
                # coefficient/exponent lists — no string-parsing ambiguity)
                @test UF.U == _tg_poly_from_terms(UF.ring, out["U"])
                @test UF.F == _tg_poly_from_terms(UF.ring, out["F"])
                # exponent + Gamma conventions (DESIGN.md C2 / wl:2996-2997,
                # wl:3010): pySecDec exponent_U/exponent_F/Gamma_factor must
                # equal our eU/eF/Γ(−eF) with ν_e = 1, D = 4−2ε
                ne, L = UF.ne, UF.L
                @test _tg_eps(out["expU"]) == EpsExp(ne - (L + 1) * 2, L + 1)
                @test _tg_eps(out["expF"]) == EpsExp(L * 2 - ne, -L)
                @test _tg_eps(out["gamma_arg"]) == EpsExp(ne - L * 2, L)
                # pySecDec's Gamma_factor sign == our eq 1.4 c = (−1)^{ν_tot}
                # (external oracle for the prefactor sign; ν_e = 1 here)
                @test out["gamma_sign"] == (iseven(ne) ? 1 : -1)
                @test all(==(1), out["powerlist"])
            end
        end
    else
        println("=" ^ 72)
        println("LOUD SKIP: python3/pySecDec NOT importable — the 5-graph")
        println("pySecDec parity suite did NOT run. Do NOT gate integration")
        println("on this run (DESIGN_B1 §3: loud skip, never a silent pass).")
        println("=" ^ 72)
        @test_skip false
    end

    # -----------------------------------------------------------------------
    # graph_to_quadruple == FROZEN Phase-A moller fixture (byte + semantic)
    # -----------------------------------------------------------------------
    @testset "moller quadruple == frozen fixture" begin
        q = graph_to_quadruple(_TG_MOL.edges, _TG_MOL.nodes, _TG_MOL.masses,
                               _TG_MOL.kin; kinvars=_TG_MOL.kinvars)
        fxpath = joinpath(_TG_FIXDIR, "moller_box_eq41.json")
        fx = load_quadruple(fxpath)
        raw = JSON.parsefile(fxpath)

        @test q.vars == fx.E.vars == Symbol[:x0, :x1, :x2]   # Cheng-Wu x3=1
        @test q.kinvars == fx.E.kinvars == Symbol[:sb, :tb, :msq]
        @test q.nu == fx.E.nu                                 # flat ν_e−1 = 0
        @test length(q.polys) == 2
        for i in 1:2
            @test _tg_poly_eq(q.polys[i][1], fx.E.polys[i][1])
            @test q.polys[i][2] == fx.E.polys[i][2]
        end
        @test q.polys[1][2] == EpsExp(0, 2)      # U^{2ε}      (fixture pin)
        @test q.polys[2][2] == EpsExp(-2, -1)    # F^{−2−ε}    (fixture pin)
        @test _tg_pref_eq(q.prefactor, fx.E.prefactor)
        @test q.prefactor.gammas == [(EpsExp(2, 1), 1)]       # Γ(2+ε)
        @test q.prefactor.c == 1                 # (−1)^{ν_tot}, ν_tot = 4

        # byte-level: emitted schema fields == the frozen fixture file
        em = quadruple_json(q; id="moller_box_eq41", title=raw["title"])
        @test em["vars"] == raw["vars"] && em["kinvars"] == raw["kinvars"]
        @test em["nu"] == raw["nu"]
        @test em["prefactor"] == raw["prefactor"]
        for i in 1:2
            @test em["polys"][i]["poly"] == raw["polys"][i]["poly"]  # bytes
            @test em["polys"][i]["exp"] == raw["polys"][i]["exp"]
        end
    end

    # -----------------------------------------------------------------------
    # Synthetic 7-var fixture vs an INDEPENDENT in-test reconstruction
    # -----------------------------------------------------------------------
    @testset "synth 7-var fixture: independent reconstruction" begin
        fx = load_quadruple(joinpath(_TG_FIXDIR, "synth_prod7_7var.json"))
        @test fx.E.vars == Symbol[:x0, :x1, :x2, :x3, :x4, :x5, :x6]
        @test isempty(fx.E.kinvars)
        @test fx.E.nu == fill(EpsExp(0), 7)  # 7 × flat 0
        @test length(fx.E.polys) == 1
        @test fx.E.polys[1][2] == EpsExp(-2, -1)   # P^(−2−ε)
        # regenerate P = prod(1+x_i) here, independently of the Python
        # generator, and require EXACT equality with the frozen bytes'
        # parse — 2^7 = 128 monomials, every coefficient 1
        P = prod(1 + g for g in gens(parent(fx.E.polys[1][1])))
        @test fx.E.polys[1][1] == P
        @test length(fx.E.polys[1][1]) == 128
        # factorization is exactly why the closed form (1+ε)^(−7) in the
        # fixture notes holds: each factor integrates to 1/(1+ε)
        @test fx.E.prefactor.gammas == []
        @test fx.E.prefactor.c == 1
    end

    # -----------------------------------------------------------------------
    # Lee-Pomeransky variant flag
    # -----------------------------------------------------------------------
    @testset "Lee-Pomeransky variant" begin
        UF = symanzik_UF(_TG_MOL.edges, _TG_MOL.nodes, _TG_MOL.masses,
                         _TG_MOL.kin; kinvars=_TG_MOL.kinvars)
        q = graph_to_quadruple(_TG_MOL.edges, _TG_MOL.nodes, _TG_MOL.masses,
                               _TG_MOL.kin; kinvars=_TG_MOL.kinvars,
                               lee_pomeransky=true)
        @test q.vars == Symbol[:x0, :x1, :x2, :x3]        # NO gauge fixing
        @test length(q.polys) == 1
        @test _tg_poly_eq(q.polys[1][1], UF.U + UF.F)     # G = U + F
        @test q.polys[1][2] == EpsExp(-2, 1)              # G^{−D/2}, D=4−2ε
        # Γ(D/2) / Γ((L+1)D/2 − ν_tot):  Γ(2−ε) / Γ(−2ε)  (L=1, ν_tot=4)
        @test q.prefactor.gammas == [(EpsExp(2, -1), 1), (EpsExp(0, -2), -1)]
        @test q.prefactor.c == 1
        @test q.nu == fill(EpsExp(0), 4)
    end

    # -----------------------------------------------------------------------
    # Gamma prefactor per paper eq 1.4 (subtropica.txt:236-249)
    # -----------------------------------------------------------------------
    @testset "eq 1.4 prefactor bookkeeping" begin
        # sunrise ν=(2,1,1): ν_tot=4, L=2 ⇒ c=+1, Γ(−eF)=Γ(0+2ε), Γ(2)⁻¹
        q = graph_to_quadruple(_TG_SUN.edges, _TG_SUN.nodes, _TG_SUN.masses,
                               _TG_SUN.kin; kinvars=_TG_SUN.kinvars,
                               nu=[2, 1, 1])
        @test q.prefactor.c == 1
        @test q.prefactor.gammas == [(EpsExp(0, 2), 1), (EpsExp(2), -1)]
        @test q.polys[1][2] == EpsExp(-2, 3)   # eU = ν−(L+1)(D−2ε)/2
        @test q.polys[2][2] == EpsExp(0, -2)   # eF = −(ν−L(D−2ε)/2)
        @test q.nu == [EpsExp(1), EpsExp(0)]   # flat ν_e−1, gauge = edge 3
        # sunrise ν=(1,1,1): ν_tot=3 ⇒ c = (−1)^{ν_tot} = −1 (eq 1.4)
        q = graph_to_quadruple(_TG_SUN.edges, _TG_SUN.nodes, _TG_SUN.masses,
                               _TG_SUN.kin; kinvars=_TG_SUN.kinvars)
        @test q.prefactor.c == -1
        @test quadruple_json(q; id="t")["prefactor"]["c"] == "-1"
        # eq 1.3 e^{εLγ_E} normalization is opt-in
        @test q.prefactor.gammaE_eps == 0
        q = graph_to_quadruple(_TG_SUN.edges, _TG_SUN.nodes, _TG_SUN.masses,
                               _TG_SUN.kin; kinvars=_TG_SUN.kinvars,
                               gammaE_norm=true)
        @test q.prefactor.gammaE_eps == 2      # L = 2
    end

    # -----------------------------------------------------------------------
    # Pinch semantics (ν_e = 0 drops the edge, keeps 0-based naming —
    # STGetIntegrandData wl:2977-2979/2984)
    # -----------------------------------------------------------------------
    @testset "pinch: sunrise ν=(1,1,0) == bubble" begin
        qp = graph_to_quadruple(_TG_SUN.edges, _TG_SUN.nodes,
                                Any[:m1sq, :m2sq, 0], _TG_SUN.kin;
                                kinvars=_TG_SUN.kinvars, nu=[1, 1, 0])
        qb = graph_to_quadruple([(1, 2), (1, 2)], [1, 2], Any[:m1sq, :m2sq],
                                _TG_SUN.kin; kinvars=_TG_SUN.kinvars)
        @test _tg_E_eq(qp, qb)
        @test qp.vars == Symbol[:x0]           # x1 gauged (last kept edge)
    end

    # -----------------------------------------------------------------------
    # Fixture-loader round-trip over ALL six frozen fixture files
    # -----------------------------------------------------------------------
    @testset "fixture loader round-trip (6 files)" begin
        # poly-string BYTE parity holds where the file's strings came from
        # make_fixtures.py mma_poly (exponent-tuple-descending): moller U/F
        # and the smirnov polys. The synth/divergent files carry hand-written
        # literals ("1 + x", ascending) — semantic parity only there.
        byte_poly_files = Set(["moller_box_eq41.json", "smirnov_tst2_5var.json",
                               "synth_prod7_7var.json"])
        make_fixture_files = ["moller_box_eq41.json", "smirnov_tst2_5var.json",
                              "synth_log2_1var.json", "synth_zeta2_2var.json",
                              "synth_zeta3_beta_1var.json",
                              "synth_prod7_7var.json",
                              "divergent_xinv_1var.json"]
        for fname in make_fixture_files
            raw = JSON.parsefile(joinpath(_TG_FIXDIR, fname))
            fx = load_quadruple(joinpath(_TG_FIXDIR, fname))
            em = quadruple_json(fx.E; id=raw["id"], title=raw["title"],
                                divergent=raw["divergent"],
                                notes=raw["notes"],
                                provenance=raw["provenance"])
            # semantic round-trip: emit -> reload -> identical integrand
            re = load_quadruple(em)
            @test _tg_E_eq(re.E, fx.E)
            @test em["divergent"] == raw["divergent"]
            # byte parity on the schema fields with a pinned canonical form
            @test em["vars"] == raw["vars"] && em["kinvars"] == raw["kinvars"]
            @test em["nu"] == raw["nu"]
            @test em["prefactor"] == raw["prefactor"]
            @test [p["exp"] for p in em["polys"]] ==
                  [p["exp"] for p in raw["polys"]]
            if fname in byte_poly_files
                @test [p["poly"] for p in em["polys"]] ==
                      [p["poly"] for p in raw["polys"]]
            end
        end
        # emitted JSON is writable + reloadable from disk too (7-var case)
        raw = JSON.parsefile(joinpath(_TG_FIXDIR, "synth_prod7_7var.json"))
        fx = load_quadruple(joinpath(_TG_FIXDIR, "synth_prod7_7var.json"))
        em = quadruple_json(fx.E; id=raw["id"], title=raw["title"])
        p = write_quadruple_json(joinpath(_TG_SCRATCH, "rt_synth7.json"), em)
        @test _tg_E_eq(load_quadruple(p).E, fx.E)
    end
end

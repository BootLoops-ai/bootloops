#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision.
r"""selftest.py -- the Longhand battery: one entry point for both routes.

ROUTE A (the arbitrary-precision parametric evaluator), default tier: the
FAST code-check battery. Real checks with known answers on every code path,
single-process (no multiprocessing pools), a few seconds of wall on a laptop:

  * analytic loop-integration kernel _I2 on BOTH discriminant closed-form
    branches AND the near-double-root series branch (small |D| of either
    sign, plus D=0 exactly) vs direct adaptive quadrature;
  * dimension-reduced box through the fixed Gauss-Legendre product engine at
    dps=15 vs the pinned cross-validated reference (>= 12 d required);
  * box depth gate: box at dps=35 must agree with the pinned 40-d reference
    to >= 30 d; the independent cross-chart reduction box_crosschart at
    dps=20 must agree to >= 18 d; and an input-truncation control (a
    kinematic constant built at 15-digit ambient precision) must cap
    agreement in the 12-22 d band -- proving the gate detects the
    inexact-input failure mode;
  * nested tanh-sinh engine on a 2-D closed form (exact value 1/2);
  * check_spec positive check + two negative controls (inhomogeneous poly,
    wrong projective degree -- both must be rejected);
  * make_eval compiled integrand vs an independent exact-Fraction evaluation
    at a rational point (exact agreement required);
  * serial (nproc=1) CBC-lattice QMC on the 3-D box spec and on the shipped
    6-D closed-form fixture toy_qmc_spec.json (>= 3 d and <= 5 sigma vs the
    exact values);
  * reduce_linear stage 1 + reduce_linear_full (degenerate and generic
    stage-2 partial-fraction eliminations vs exact values) +
    reduce_one_quadratic vs the numeric kernel on both branches (named SKIP
    if sympy is absent);
  * assemble_integrand_expr general-(a,b) fallback: an a=-1,b=1 topology
    must emit the unreduced U^a/F^b Cheng-Wu integrand, checked exactly at a
    rational point (named SKIP if sympy is absent);
  * pySecDec U/F front end: the one-loop box propagator list must reproduce
    the pinned box spec's F polynomial exactly (named SKIP if pySecDec is
    absent -- it is only needed for this front end).

ROUTE B (longhand.disteval, the pySecDec disteval driver): the pytest battery
under tests/ (test_disteval.py), run in a subprocess with `-rs` so skipped
legs are named: the acceptance check over the worked example's pinned records
(examples/ndpent_top) must reproduce the record by name, the compare rule, the
chunk split, the per-chunk epsabs rule, the memory basis, every refusal by
name, the planted and foreign controls, the family / pins legs.  The two legs
that run disteval's stand-in form import pySecDec and SKIP BY NAME when
pySecDec (with pySecDecContrib) is not installed; every other leg replays the
shipped records and needs no pySecDec.  Requires pytest.

--full: ROUTE A's heavy tier instead of the fast one -- box benchmark at
dps=30 vs the pinned cross-validated 40-d reference (>= 25 d) + the 6-D QMC
farm on the fixture (N=8009, 16 shifts, 8 worker processes) checked against
the hand-provable exact value 13/1080 within its own error bar + the generic
UF-spec pathway's box control through an 8-process pool (the legs of
selftest_parametric.py).  Multi-process; wall depends on core count.  Route B
runs the same battery in both tiers.

--route parametric | disteval | all (default all): run one route's legs only.

Exits nonzero on failure.  Missing gmpy2 is a named hard failure (it is a
listed requirement of the engines), as is missing pytest for route B;
missing pySecDec/sympy are named skips.
"""
import argparse
import json
import os
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mpmath as mp

FAIL = 0


def _agree_digits(dev, ref):
    """Decimal digits of agreement; 99 flags exact-to-working-precision."""
    mp.mp.dps = 40
    rel = abs(mp.mpf(str(dev))) / abs(mp.mpf(str(ref)))
    return 99 if rel == 0 else -int(mp.log10(rel))


def _leg(tag, ok, detail):
    global FAIL
    print(f"[{tag}] {detail}")
    if ok:
        print("  PASS")
    else:
        print("  FAIL")
        FAIL += 1


def _skip(tag, why):
    print(f"[{tag}] SKIP ({why})")


# --------------------------------------------------------------------------- #
#  fast tier
# --------------------------------------------------------------------------- #
def fast_kernel():
    """_I2 vs direct adaptive quadrature: both discriminant closed-form
    branches plus the near-double-root series branch (|D| <= Q^2/8, both
    signs, and D = 0 exactly)."""
    from bench_box import _I2
    mp.mp.dps = 30
    worst = 99
    for (P, Q, R) in [(mp.mpf(2), mp.mpf(1), mp.mpf(3)),      # D = 23 > 0
                      (mp.mpf(1), mp.mpf(5), mp.mpf(2)),      # D = -17 < 0
                      (mp.mpf(1), mp.mpf(10), mp.mpf('25.2')),  # D = 0.8: series
                      (mp.mpf(1), mp.mpf(10), mp.mpf('24.9')),  # D = -0.4: series
                      (mp.mpf(2), mp.mpf(8), mp.mpf(8))]:       # D = 0 exactly
        closed = _I2(P, Q, R)
        direct = mp.quad(lambda x: 1 / (P * x * x + Q * x + R) ** 2, [0, mp.inf])
        worst = min(worst, _agree_digits(closed - direct, direct))
    _leg("kernel", worst >= 20,
         f"_I2 vs direct quad, closed + series branches: worst agree ~{worst} d "
         f"(need >=20)")


# 40-significant-digit box reference at s=-1, t=-1/3, m2=1, M5=2:
# cross-validated to 46+ d by dps=45-vs-60 runs of box AND the independent
# box_crosschart reduction (different Cheng-Wu chart, different eliminated
# parameter); a nested adaptive tanh-sinh run agrees to its own ~25 d floor.
BOX_REF_STR = '0.1245570572327092697104024441545542480922'


def fast_box_gauss():
    """Dimension-reduced box through the Gauss product engine at dps=15."""
    from fractions import Fraction
    from bench_box import box
    mp.mp.dps = 45
    REF = mp.mpf(BOX_REF_STR)
    v = box(-1, Fraction(-1, 3), 1, 2, dps=15)
    agree = _agree_digits(mp.mpf(str(v)) - REF, REF)
    _leg("box gauss", agree >= 12,
         f"M5=2 dps=15 -> agree vs pinned reference ~{agree} d (need >=12)")


def fast_box_depth():
    """Depth gate: >=30 d at dps=35, the independent cross-chart reduction,
    and an input-truncation control the gate must be able to see."""
    from fractions import Fraction
    from bench_box import box, box_crosschart
    mp.mp.dps = 45
    REF = mp.mpf(BOX_REF_STR)
    v35 = box(-1, Fraction(-1, 3), 1, 2, dps=35)
    mp.mp.dps = 45
    a35 = _agree_digits(mp.mpf(str(v35)) - REF, REF)
    vx = box_crosschart(-1, Fraction(-1, 3), 1, 2, dps=20)
    mp.mp.dps = 45
    ax = _agree_digits(mp.mpf(str(vx)) - REF, REF)
    # control: a kinematic input rounded at 15-digit ambient precision is a
    # different rational number; agreement must cap near the input precision
    # (proves this gate catches the inexact-input failure mode)
    mp.mp.dps = 15
    t_trunc = mp.mpf(-1) / 3
    vt = box(-1, t_trunc, 1, 2, dps=35)
    mp.mp.dps = 45
    at = _agree_digits(mp.mpf(str(vt)) - REF, REF)
    ok = (a35 >= 30) and (ax >= 18) and (12 <= at <= 22)
    _leg("box depth", ok,
         f"dps=35 ~{a35} d (need >=30); cross-chart dps=20 ~{ax} d (need >=18); "
         f"truncated-input control ~{at} d (need 12-22)")


def fast_tanhsinh():
    """Nested tanh-sinh engine on a 2-D closed form: int (1+x+y)^-3 = 1/2."""
    from hiprec import integrate_tanhsinh
    v = integrate_tanhsinh(lambda x, y: 1 / (1 + x + y) ** 3, 2, dps=12)
    agree = _agree_digits(v - mp.mpf(1) / 2, mp.mpf(1) / 2)
    _leg("tanh-sinh 2d", agree >= 8,
         f"int (1+x+y)^-3 over [0,inf)^2 -> agree vs 1/2 ~{agree} d (need >=8)")


def fast_spec_checks():
    """check_spec positive check + two negative controls."""
    import parametric
    E, degs = parametric.check_spec(parametric.box_spec())
    ok = (E == 4 and degs.get('F0') == 2)
    neg = 0
    bad_inhom = {'name': 'bad', 'n_den': 3,
                 'polys': {'F0': [[[1, 0, 0], [1, 1]], [[1, 1, 0], [1, 1]]]},
                 'numerator': {'gamma': 1, 'terms': [[[1, 1], [['F0', -2]]]]}}
    try:
        parametric.check_spec(bad_inhom)
    except AssertionError:
        neg += 1
    bad_deg = {'name': 'bad2', 'n_den': 3,
               'polys': {'F0': [[[1, 1, 0], [1, 1]], [[0, 1, 1], [1, 1]]]},
               'numerator': {'gamma': 1, 'terms': [[[1, 1], [['F0', -2]]]]}}
    try:
        parametric.check_spec(bad_deg)          # degree -4 != -3: must raise
    except AssertionError:
        neg += 1
    _leg("spec checks", ok and neg == 2,
         f"box spec E={E} deg(F0)={degs.get('F0')}; negative controls rejected {neg}/2")


def fast_make_eval():
    """Compiled sparse evaluator vs independent exact-Fraction evaluation."""
    from fractions import Fraction
    import gmpy2
    from gmpy2 import mpfr
    import parametric
    spec = parametric.box_spec()
    gmpy2.get_context().precision = 150
    f, dim = parametric.make_eval(spec)
    pt = [Fraction(1, 3), Fraction(2, 5), Fraction(7, 4)]
    got = f(*[mpfr(p.numerator) / p.denominator for p in pt])
    xs = pt + [Fraction(1)]                     # Cheng-Wu chart var = 1
    F = Fraction(0)
    for m, (n_, d_) in spec['polys']['F0']:
        term = Fraction(n_, d_)
        for xi, e in zip(xs, m):
            term *= xi ** e
        F += term
    want = F ** -2                              # gamma=1, single term F0^-2
    rel = abs(got - mpfr(want.numerator) / want.denominator) / abs(got)
    _leg("make_eval", dim == 3 and rel < mpfr('1e-30'),
         f"rational-point check vs Fraction reference: rel dev {float(rel):.1e} (need <1e-30)")


def fast_qmc_box3d():
    """Serial CBC-QMC (no pool) on the 3-D box spec vs the pinned reference."""
    import gmpy2
    from gmpy2 import mpfr
    import hiprec
    import parametric
    spec = parametric.box_spec()
    with tempfile.NamedTemporaryFile('w', suffix='.json', delete=False) as fh:
        json.dump({'box': spec}, fh)
        bpath = fh.name
    try:
        res = hiprec.integrate_qmc(parametric.spec_factory,
                                   {'json': bpath, 'family': 'box'},
                                   3, N=8009, n_shifts=4, korobov_p=3,
                                   dps=20, nproc=1, seed=2)
    finally:
        os.unlink(bpath)
    REF = mpfr(BOX_REF_STR)
    dev = abs(res['value'] - REF)
    pull = float(dev / res['error']) if res['error'] > 0 else float('inf')
    agree = -int(gmpy2.log10(dev)) if dev > 0 else 30
    _leg("qmc box 3d", agree >= 3 and pull <= 5,
         f"serial N=8009 4 shifts -> agree ~{agree} d pull {pull:.2f} sigma (need >=3 d, <=5 sigma)")


def fast_qmc_toy6d():
    """Serial CBC-QMC (no pool) on the shipped 6-D closed-form fixture."""
    import gmpy2
    from gmpy2 import mpfr
    import parametric
    SPEC = os.path.join(HERE, 'toy_qmc_spec.json')
    res = parametric.integrate_spec(SPEC, N=8009, n_shifts=4, korobov_p=3,
                                    dps=20, nproc=1, seed=3)
    en, ed = parametric.load_spec(SPEC)['exact']
    EXACT = mpfr(en) / mpfr(ed)
    dev = abs(res['value'] - EXACT)
    pull = float(dev / res['error']) if res['error'] > 0 else float('inf')
    agree = -int(gmpy2.log10(dev)) if dev > 0 else 30
    _leg("qmc toy 6d", agree >= 3 and pull <= 5,
         f"serial N=8009 4 shifts -> I vs exact {en}/{ed}: agree ~{agree} d "
         f"pull {pull:.2f} sigma (need >=3 d, <=5 sigma)")


def fast_reductions():
    """reduce_linear stage 1 + reduce_linear_full stage 2 + reduce_one_quadratic
    vs exact values / the numeric kernel."""
    try:
        import sympy as sp
    except ImportError:
        _skip("reductions", "sympy not installed -- symbolic reduction legs not exercised")
        return
    import feynman
    from bench_box import _I2
    x0, x1, x2 = sp.symbols('x0 x1 x2')
    # stage 1: F linear in x0 after Cheng-Wu x2=1 -> factors (x1+1)*(x1+1)
    F = x0 * (x1 + x2) + x1 * x2 + x2 ** 2
    den, rem = feynman.reduce_linear(F, [x0, x1, x2])
    ok_lin = (len(den) == 2 and rem == [x1]
              and sp.expand(den[0] * den[1] - (x1 + 1) ** 2) == 0)
    # stage 2 degenerate (proportional factors): full reduction to the exact
    # constant 1  (int_0^inf dx1/(x1+1)^2)
    e_deg, rem_deg = feynman.reduce_linear_full(F, [x0, x1, x2])
    ok_deg = (rem_deg == [] and sp.simplify(e_deg - 1) == 0)
    # stage 2 generic: Cheng-Wu Fc = (x1+2)x0 + (2x1+1); the iterated closed
    # forms give int dx1/((x1+2)(2x1+1)) = log(4)/3 exactly
    F2 = x0 * x1 + 2 * x0 * x2 + 2 * x1 * x2 + x2 ** 2
    e_gen, rem_gen = feynman.reduce_linear_full(F2, [x0, x1, x2])
    mp.mp.dps = 30
    ok_gen = (rem_gen == []
              and _agree_digits(mp.mpf(str(sp.N(e_gen, 30))) - mp.log(4) / 3,
                                mp.log(4) / 3) >= 25)
    # no stage 2 when a variable is linear in only one power-1 factor
    # (that integral diverges): expr stays 1/(A*B)
    F3 = x0 * x1 + x0 * x2 + x1 * x1 + x2 * x2
    e_one, rem_one = feynman.reduce_linear_full(F3, [x0, x1, x2])
    ok_one = (rem_one == [x1]
              and sp.simplify(e_one - 1 / ((x1 + 1) * (x1 ** 2 + 1))) == 0)
    # symbolic quadratic elimination must match the numeric kernel on both branches
    worst = 99
    for (P, Q, R) in [(2, 1, 3), (1, 5, 2)]:            # D>0 and D<0
        expr, rem2 = feynman.reduce_one_quadratic(
            sp.expand(P * x0 ** 2 + Q * x0 + R), [x0, x1], x0)
        v = mp.mpf(str(sp.N(expr, 25)))
        mp.mp.dps = 25
        ref = _I2(mp.mpf(P), mp.mpf(Q), mp.mpf(R))
        worst = min(worst, _agree_digits(v - ref, ref))
    _leg("reductions", ok_lin and ok_deg and ok_gen and ok_one and worst >= 10,
         f"stage1 ok={ok_lin}; full: degenerate ok={ok_deg}, generic ok={ok_gen}, "
         f"single-factor-stops ok={ok_one}; reduce_one_quadratic vs _I2 "
         f"worst ~{worst} d (need >=10)")


def fast_general_ab():
    """assemble_integrand_expr general-(a,b) fallback: unreduced U^a/F^b
    Cheng-Wu integrand, checked exactly at a rational point."""
    try:
        import sympy as sp
    except ImportError:
        _skip("general (a,b)", "sympy not installed -- fallback leg not exercised")
        return
    import feynman
    x0, x1, x2 = sp.symbols('x0 x1 x2')
    m2, s = sp.symbols('m2 s')
    U = x0 + x1 + x2
    F = m2 * (x0 + x1 + x2) ** 2 - s * x0 * x1
    # N=3 propagators, L=1 loop at d=4: a = -1, b = 1 -> no reduction, emit U^a/F^b
    expr, rem = feynman.assemble_integrand_expr(
        U, F, [x0, x1, x2], 1, 3, {m2: sp.Integer(1), s: sp.Integer(-1)})
    r = sp.Rational
    pt = {x0: r(1, 3), x1: r(2, 5)}
    Ucw = (x0 + x1 + 1).subs(pt)
    Fcw = ((x0 + x1 + 1) ** 2 + x0 * x1).subs(pt)
    want = Ucw ** (-1) * Fcw ** (-1)
    ok = (rem == [x0, x1] and sp.simplify(expr.subs(pt) - want) == 0)
    _leg("general (a,b)", ok,
         f"a=-1,b=1 fallback: vars {rem} (need [x0, x1]); rational-point value "
         f"exact match={sp.simplify(expr.subs(pt) - want) == 0}")


def fast_uf_frontend():
    """pySecDec U/F extraction must reproduce the pinned box spec exactly."""
    try:
        import pySecDec                              # noqa: F401
        import sympy as sp
    except ImportError:
        _skip("U/F front end",
              "pySecDec not installed -- propagator front end not exercised; "
              "pip install pySecDec to enable")
        return
    from fractions import Fraction
    import feynman
    import parametric
    props = ['k**2 - M5', '(k+p1)**2 - m2', '(k+p1+p2)**2 - m2', '(k-p4)**2 - m2']
    rr = [('p1*p1', 0), ('p2*p2', 0), ('p3*p3', 0), ('p4*p4', 0),
          ('p1*p2', 's/2'), ('p3*p4', 's/2'),
          ('p2*p3', '(-s-t)/2'), ('p1*p4', '(-s-t)/2'),
          ('p1*p3', 't/2'), ('p2*p4', 't/2')]
    U, F, xs, L, Nprop = feynman.build_UF(props, ['k'], ['p1', 'p2', 'p3', 'p4'], rr)
    s_, t_, m2_, M5_ = sp.symbols('s t m2 M5')
    Fnum = sp.expand(F.subs({s_: sp.Rational(-1), t_: sp.Rational(-1, 3),
                             m2_: 1, M5_: 2}))
    spec = parametric.box_spec(s=-1, t=Fraction(-1, 3), m2=1, M5=2)
    Fexp = sp.S(0)
    for m, (n_, d_) in spec['polys']['F0']:
        term = sp.Rational(n_, d_)
        for xi, e in zip(xs, m):
            term *= xi ** e
        Fexp += term
    ok = (L == 1 and Nprop == 4
          and sp.expand(U - sum(xs)) == 0 and sp.expand(Fnum - Fexp) == 0)
    _leg("U/F front end", ok,
         f"one-loop box: L={L} N={Nprop}, U == sum(x): "
         f"{sp.expand(U - sum(xs)) == 0}, F == pinned box spec F0: "
         f"{sp.expand(Fnum - Fexp) == 0}")


def run_fast():
    t0 = time.time()
    try:
        import gmpy2                                  # noqa: F401
    except ImportError:
        print("FAIL: requirement gmpy2 is not installed (pip install gmpy2) -- "
              "the engines cannot run")
        return 1
    fast_kernel()
    fast_box_gauss()
    fast_box_depth()
    fast_tanhsinh()
    fast_spec_checks()
    fast_make_eval()
    fast_qmc_box3d()
    fast_qmc_toy6d()
    fast_reductions()
    fast_general_ab()
    fast_uf_frontend()
    verdict = "FAIL" if FAIL else "PASS"
    print(f"parametric fast battery {verdict}  ({time.time() - t0:.1f}s wall; "
          f"single-process -- run --full for the multi-process farm)")
    return 1 if FAIL else 0


# --------------------------------------------------------------------------- #
#  --full tier: the original heavy battery, unchanged
# --------------------------------------------------------------------------- #
def run_full():
    """(1) box benchmark at dps=30 vs the pinned cross-validated 40-d
    reference -> expect ~29-30 d, far past the 8-9 d double-precision ceiling;
    (2) a 6-D QMC run on the shipped synthetic fixture toy_qmc_spec.json --
    every term has a Dirichlet closed form (exact value 13/1080; derivation
    in the fixture's note field), so the whole parametric+QMC pathway is
    checked against a hand-provable answer, within its own error bar."""
    global FAIL

    # (1) box, high precision (exact rational kinematics -- see bench_box's
    # precision contract)
    from fractions import Fraction
    from bench_box import box
    mp.mp.dps = 50
    REF = mp.mpf(BOX_REF_STR)   # 40 significant digits, cross-validated to 46+
    v = box(-1, Fraction(-1, 3), 1, 2, dps=30)
    mp.mp.dps = 50
    agree = -int(mp.log10(abs(mp.mpf(str(v)) - REF)))
    print(f"[box gauss] M5=2 dps=30 -> {mp.nstr(v, 32)}  agree vs ref ~{agree} d")
    if agree < 25:
        print("  FAIL: box below the 25 d agreement floor")
        FAIL += 1
    else:
        print("  PASS (>=25 d, far past the 8-9 d double-precision ceiling)")

    # (2) 6-D QMC farm on the synthetic closed-form fixture
    import gmpy2
    from gmpy2 import mpfr
    import parametric
    SPEC = os.path.join(HERE, 'toy_qmc_spec.json')
    res = parametric.integrate_spec(SPEC, N=8009, n_shifts=16, korobov_p=3,
                                    dps=25, nproc=8, seed=3)
    en, ed = parametric.load_spec(SPEC)['exact']
    EXACT = mpfr(en) / mpfr(ed)     # 13/1080, provable by hand (see fixture note)
    dev = abs(res['value'] - EXACT)
    pull = float(dev / res['error']) if res['error'] > 0 else float('inf')
    agree = -int(gmpy2.log10(dev)) if dev > 0 else 30
    print(f"[toy6d qmc] I={str(res['value'])[:14]} exact {en}/{ed}  "
          f"agree ~{agree} d  pull {pull:.2f} sigma")
    if agree < 3 or pull > 5:
        print("  FAIL: 6-D QMC off the closed form (>=3 d and <=5 sigma required)")
        FAIL += 1
    else:
        print("  PASS (matches the Dirichlet closed form within its error bar)")

    # (3) the generic UF-spec pathway's box control through an 8-process pool
    # (the legs of selftest_parametric.py, run as its own process so its pool
    # starts from a real __main__)
    sys.stdout.flush()
    rc = subprocess.run([sys.executable, os.path.join(HERE, 'selftest_parametric.py')]).returncode
    print(f"[parametric pool] selftest_parametric.py rc={rc}")
    if rc != 0:
        print("  FAIL: the generic-pathway box control (pool) did not pass")
        FAIL += 1
    else:
        print("  PASS")

    return 1 if FAIL else 0


# --------------------------------------------------------------------------- #
#  route B: the disteval driver's pytest battery (tests/test_disteval.py)
# --------------------------------------------------------------------------- #
def run_disteval():
    """Run the disteval route's battery in a subprocess: pytest over tests/
    with -rs (named skips) and no cache dir written into the tree.  The two
    stand-in stage legs skip by name without pySecDec; the record-replay legs
    always run.  Returns the pytest rc (0 = all passed or skipped by name)."""
    t0 = time.time()
    try:
        import pytest                                # noqa: F401
    except ImportError:
        print("[disteval] FAIL: requirement pytest is not installed "
              "(pip install pytest) -- the disteval route's battery cannot run")
        return 1
    tests = os.path.join(HERE, 'tests')
    cmd = [sys.executable, '-m', 'pytest', '-q', '-rs', '-p', 'no:cacheprovider', tests]
    print(f"[disteval] {' '.join(cmd[1:])}", flush=True)
    sys.stdout.flush()
    rc = subprocess.run(cmd, cwd=HERE).returncode
    verdict = "PASS" if rc == 0 else "FAIL"
    print(f"[disteval] pytest rc={rc}\n  {verdict}")
    print(f"disteval battery {verdict}  ({time.time() - t0:.1f}s wall)")
    return rc


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="python3 selftest.py",
        description="Longhand battery: route A (arbitrary-precision parametric "
                    "evaluator) + route B (pySecDec disteval driver)")
    ap.add_argument("--full", action="store_true",
                    help="route A: run the heavy tier (multi-process QMC farm + "
                         "the parametric pool control) instead of the fast "
                         "code-check tier")
    ap.add_argument("--route", choices=["parametric", "disteval", "all"], default="all",
                    help="run one route's legs only (default: all)")
    args = ap.parse_args(argv)
    t0 = time.time()
    rc_a = rc_b = 0
    if args.route in ("parametric", "all"):
        print("=== route A: arbitrary-precision parametric evaluator "
              f"({'full' if args.full else 'fast'} tier) ===", flush=True)
        rc_a = run_full() if args.full else run_fast()
    if args.route in ("disteval", "all"):
        print("=== route B: pySecDec disteval driver (tests/test_disteval.py) ===", flush=True)
        rc_b = run_disteval()
    verdict = "FAIL" if (rc_a or rc_b) else "PASS"
    print(f"longhand battery {verdict}  (route A rc={rc_a}, route B rc={rc_b}; "
          f"{time.time() - t0:.1f}s wall)")
    sys.exit(1 if (rc_a or rc_b) else 0)


if __name__ == '__main__':
    main()

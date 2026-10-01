#!/usr/bin/env python3
"""
test_conic_thirdkind.py — member battery for wayfinder/conic_thirdkind.py
(the conic third-kind reduction core).

Legs (all fast — the whole battery runs in seconds):
  * engine truth — hermite_conic_sqf on a small conic sig with the exact
    cancel==0 certificate PLUS an independent mpmath quadrature check of
    the reduced identity (no fabricated digits; reference computed
    in-test at two precisions);
  * deg-2 REFUSAL — the SCOPE leg: deg-3 q must raise ConicScopeError
    naming the conic deg-2 scope, from the engine AND from the builder
    (the Hermite tail's (j+1)*lc elimination is deg-2-only; lifting =
    second-kind basis extension, deliberately NOT done — by design);
  * truth case — rebuild the toy system from the family fixture via
    rows_from_connection + build_row_system (symbols/endpoints passed
    EXPLICITLY — the parameterization under test; the spectator pin
    exercises `subs`) and compare states/atoms/Q/numer/Amat ENTRY-EXACT
    against the banked toy fixture;
  * mutation control (comparator teeth) — a 1e-20-class exact rational
    perturbation of one banked Amat entry MUST be detected by the same
    comparator the truth case uses (exact sympy, never string/float);
  * mutation control (certificate teeth) — a corrupted Hermite triple
    must fail the exact cancel==0 certificate;
  * banked-row quadrature — INDEPENDENT numeric verification of the
    banked system: the L row and the J0 row of the banked Amat are
    checked against direct mpmath quadrature of their defining
    identities at a pinned rational t, at two working precisions
    (the reference values are recomputed in-test, never stored);
  * ExtSys numeric — the banked toy system wrapped as a duck-typed
    DESystem: A(x) Horner evaluation vs direct exact-sympy substitution
    at two working precisions; declared singular_points nonempty (the
    moving poles).

Fixtures: tests/fixtures_conic_thirdkind/ — a SMALL synthetic conic
third-kind family (toy_family.json: one fenced row with a moving pole
and a deg-2 irreducible third-kind factor, plus a spectator symbol
pinned at load) and its built system (toy_banked.json). The banked file
is regenerated in-run by the truth-case leg from the family fixture
through the shipped builder, and cross-checked numerically by the
quadrature leg — nothing in the fixtures is trusted without being
recomputed here.

Runnable directly (python3 tests/test_conic_thirdkind.py) or under pytest.
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import sympy as sp                                       # noqa: E402
import conic_thirdkind as ct                             # noqa: E402

FIX = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   'fixtures_conic_thirdkind')
SYM = {n: sp.Symbol(n, real=True) for n in ('x', 't', 'u')}
x, t, u = SYM['x'], SYM['t'], SYM['u']

ENDPOINTS = (-sp.Rational(1, 2), sp.Rational(1, 2))      # exact in binary
PIN = {u: sp.Integer(1)}                                 # spectator pin


# ---------------------------------------------------------------- comparator
def entries_equal(a, b):
    """Exact equality of two sympy expressions (srepr fast path, then exact
    cancel of the difference). NEVER string/float-compares values."""
    if sp.srepr(a) == sp.srepr(b):
        return True
    return sp.cancel(sp.together(a - b)) == 0


def compare_system(built, banked):
    """built = build_row_system dict; banked = toy_banked.json dict.
    Returns list of mismatch strings (empty = pass)."""
    bad = []
    if built['states'] != banked['states']:
        bad.append('states: %s != %s' % (built['states'], banked['states']))
    ba = [[sp.srepr(f), int(k)] for f, k in built['atoms']]
    if ba != banked['atoms']:
        bad.append('atoms differ')
    if not entries_equal(built['Q'], ct.load_expr(banked['Q'], SYM)):
        bad.append('Q differs')
    if built.get('numer') is not None and banked.get('numer'):
        if not entries_equal(built['numer'],
                             ct.load_expr(banked['numer'], SYM)):
            bad.append('numer differs')
    n = len(banked['states'])
    for i in range(n):
        for j in range(n):
            if not entries_equal(built['Amat'][i, j],
                                 ct.load_expr(banked['Amat'][i][j], SYM)):
                bad.append('Amat[%d,%d] differs' % (i, j))
    return bad


def _load_family():
    return json.load(open(os.path.join(FIX, 'toy_family.json')))


def _load_banked():
    return json.load(open(os.path.join(FIX, 'toy_banked.json')))


def _banked_sysd(banked):
    """Banked JSON -> the sysd dict shape ExtSys consumes."""
    n = len(banked['states'])
    Amat = sp.Matrix(n, n, lambda i, j: ct.load_expr(banked['Amat'][i][j],
                                                     SYM))
    return {'tag': banked['tag'], 'states': banked['states'],
            'Amat': Amat, 't': t}


def _rebuild_from_family():
    fam = _load_family()
    row, mrow, Q, numer = ct.rows_from_connection(
        fam, 'toy_row', 'toy_M', 'toy', 't', SYM, subs=PIN)
    return ct.build_row_system(row, mrow, Q, x, t, ENDPOINTS,
                               tag='toy_row', numer=numer, verbose=False)


# --------------------------------------------------------------- fast legs
def test_engine_truth_small_conic():
    """Certified Hermite reduction on a small conic + independent
    quadrature identity check at two precisions (no fabricated digits)."""
    import mpmath
    q = 1 - x**2 + sp.Rational(1, 3) * x
    sig = (x**3 + 2 * x - sp.Rational(1, 5)) / (x**2 + 1)
    H, c, T = ct.hermite_certified(sig, q, x)     # cert asserts inside
    # independent identity check: INT sig/Y == [H*Y] + c*INT 1/Y + INT T/Y
    A, B = -sp.Rational(1, 2), sp.Rational(1, 2)          # exact in binary
    cr = sp.Rational(c)                    # c is exact rational here (assert)
    worst = None
    for dps in (30, 50):
        with mpmath.workdps(dps):
            c_mp = mpmath.mpf(int(cr.p)) / mpmath.mpf(int(cr.q))
            f_l = sp.lambdify(x, sig / sp.sqrt(q), 'mpmath')
            lhs = mpmath.quad(f_l, [float(A), float(B)])
            HY = sp.lambdify(x, H * sp.sqrt(q), 'mpmath')
            one = sp.lambdify(x, 1 / sp.sqrt(q), 'mpmath')
            TY = (sp.lambdify(x, T / sp.sqrt(q), 'mpmath')
                  if T != 0 else None)
            rhs = HY(float(B)) - HY(float(A)) \
                + c_mp * mpmath.quad(one, [float(A), float(B)])
            if TY is not None:
                rhs += mpmath.quad(TY, [float(A), float(B)])
            err = abs(mpmath.mpc(lhs) - mpmath.mpc(rhs))
            worst = err if worst is None else max(worst, err)
            assert err < mpmath.mpf(10) ** (-(dps - 8)), \
                'identity residual %s at dps %d' % (err, dps)
    print('  engine truth: identity residual worst %s (dps 30/50)' % worst)


def test_refusal_deg3():
    """SCOPE leg: deg-3 q refuses loudly, naming the conic deg-2 scope."""
    tt = sp.Symbol('tt', real=True)
    q3 = x**3 - x + sp.Rational(1, 7)
    try:
        ct.hermite_conic_sqf(1 / (x - 3), q3, x)
        raise AssertionError('deg-3 q was ACCEPTED by hermite_conic_sqf')
    except ct.ConicScopeError as e:
        msg = str(e)
        assert 'deg-2' in msg and 'CONIC' in msg, \
            'refusal does not name the scope: %s' % msg
    # builder-level refusal too (the fence is at every entry)
    row = {'H': sp.Integer(0), 'c': sp.Integer(1), 'T': 1 / (x - 3)}
    mrow = {'H': sp.Integer(0), 'c': sp.Integer(0)}
    try:
        ct.build_row_system(row, mrow, q3 + tt * x, x, tt,
                            (-sp.Rational(2, 5), sp.Rational(3, 10)),
                            verbose=False)
        raise AssertionError('deg-3 Q was ACCEPTED by build_row_system')
    except ct.ConicScopeError:
        pass
    print('  refusal: deg-3 q refused by engine AND builder, scope named')


def test_truth_case_toy_rebuild():
    """Rebuild the toy system from the family fixture through the
    PARAMETERIZED member API (spectator pinned via `subs`) and compare
    entry-exact against the banked toy fixture."""
    built = _rebuild_from_family()
    banked = _load_banked()
    bad = compare_system(built, banked)
    assert not bad, 'truth case mismatches: %s' % bad[:10]
    ea = [sp.srepr(e) for e in built['endpoints']]
    assert ea == banked['endpoints'], 'endpoints differ'
    print('  truth case: toy_row rebuilt (%.1f s), n=%d states, '
          'Amat/atoms/Q/numer ENTRY-EXACT vs banked fixture'
          % (built['wall_s'], len(built['states'])))


def test_mutation_comparator_teeth():
    """A (1+1e-20)-class exact perturbation of ONE banked entry must be
    DETECTED by the truth-case comparator."""
    banked = _load_banked()
    mut = json.loads(json.dumps(banked))          # deep copy
    n = len(banked['states'])
    tgt = None
    for i in range(n):
        for j in range(n):
            if banked['Amat'][i][j] != 'Integer(0)':
                tgt = (i, j)
                break
        if tgt:
            break
    e = ct.load_expr(banked['Amat'][tgt[0]][tgt[1]], SYM)
    mut['Amat'][tgt[0]][tgt[1]] = sp.srepr(
        e * (1 + sp.Rational(1, 10**20)))
    sysd = _banked_sysd(mut)
    sysd['atoms'] = [(ct.load_expr(f, SYM), k) for f, k in mut['atoms']]
    sysd['Q'] = ct.load_expr(mut['Q'], SYM)
    sysd['numer'] = None
    bad = compare_system(sysd, banked)
    assert bad, 'mutated Amat[%d,%d] (x(1+1e-20)) NOT detected' % tgt
    # and the unmutated roundtrip is clean (the comparator is not just angry)
    sysd0 = _banked_sysd(banked)
    sysd0['atoms'] = [(ct.load_expr(f, SYM), k) for f, k in banked['atoms']]
    sysd0['Q'] = ct.load_expr(banked['Q'], SYM)
    sysd0['numer'] = None
    assert not compare_system(sysd0, banked), 'clean roundtrip flagged'
    print('  mutation control: comparator detected %s at Amat[%d,%d]; '
          'clean roundtrip green' % (bad[0].split()[0], *tgt))


def test_mutation_certificate_teeth():
    """A corrupted Hermite triple must fail the exact cancel==0 cert."""
    q = 1 - x**2
    sig = x / (x**2 + 2)
    H, c, T = ct.hermite_conic_sqf(sig, q, x)
    cert = sp.cancel(ct.odd_deriv(H, q, x) + c + T - sig)
    assert cert == 0, 'true triple failed its own certificate'
    Hbad = H + sp.Rational(1, 10**20)
    certbad = sp.cancel(ct.odd_deriv(Hbad, q, x) + c + T - sig)
    assert certbad != 0, '1e-20-corrupted H PASSED the exact certificate'
    print('  mutation control: exact certificate rejects corrupted H')


def test_banked_rows_quadrature():
    """INDEPENDENT numeric check of the banked system: the L row and the
    J0 row must reproduce direct quadrature of their defining identities
    at a pinned rational t, at two working precisions. The reference
    values are recomputed here from the family fixture's definitions —
    nothing stored is trusted."""
    import mpmath
    fam = _load_family()
    row, mrow, Q, _numer = ct.rows_from_connection(
        fam, 'toy_row', 'toy_M', 'toy', 't', SYM, subs=PIN)
    banked = _load_banked()
    states = banked['states']
    n = len(states)
    Amat = [[ct.load_expr(banked['Amat'][i][j], SYM) for j in range(n)]
            for i in range(n)]
    atoms = [(ct.load_expr(f, SYM), k) for f, k in banked['atoms']]
    t0 = sp.Rational(1, 10)
    A_, B_ = ENDPOINTS
    q0 = Q.subs(t, t0)                       # conic at the pinned point
    # sanity: q positive on the closed panel (endpoint + midpoint checks)
    for pt in (A_, B_, sp.Integer(0)):
        assert q0.subs(x, pt) > 0

    def _rat(e):
        r = sp.Rational(sp.cancel(e))
        return r

    worst = None
    for dps in (30, 50):
        with mpmath.workdps(dps):
            Y = sp.lambdify(x, sp.sqrt(q0), 'mpmath')

            def quad_over_Y(numer_expr):
                f_l = sp.lambdify(x, numer_expr / sp.sqrt(q0), 'mpmath')
                return mpmath.quad(f_l, [float(A_), float(B_)])

            # state values, each from its own definition
            val = {'L': None, 'F': None,     # never referenced by these rows
                   'M': quad_over_Y(sp.Integer(1)),
                   'BcA': Y(float(A_)), 'BcB': Y(float(B_))}
            for a_i, (f, k) in enumerate(atoms):
                val['J%d' % a_i] = quad_over_Y(x**k / f.subs(t, t0))

            def row_rhs(i):
                acc = mpmath.mpf(0)
                for j, s in enumerate(states):
                    ent = Amat[i][j]
                    if ent == 0:
                        continue
                    assert val[s] is not None, \
                        'row %s references un-evaluated state %s' % (i, s)
                    r = _rat(ent.subs(t, t0))
                    acc += (mpmath.mpf(int(r.p)) / mpmath.mpf(int(r.q))) \
                        * val[s]
                return acc

            # L row: dL/dt = [H_L*Y]_A^B + c_L*M + INT T_L/Y (family defs)
            HY = sp.lambdify(x, row['H'].subs(t, t0) * sp.sqrt(q0), 'mpmath')
            cL = _rat(row['c'].subs(t, t0))
            lhs_L = HY(float(B_)) - HY(float(A_)) \
                + (mpmath.mpf(int(cL.p)) / mpmath.mpf(int(cL.q))) * val['M'] \
                + quad_over_Y(row['T'].subs(t, t0))
            err_L = abs(lhs_L - row_rhs(states.index('L')))

            # J0 row: dJ0/dt = INT sig/Y with sig = -x^k f_t/f^2 - x^k Q_t/(2 f Q)
            f0, k0 = atoms[0]
            sig = (-x**k0 * sp.diff(f0, t) / f0**2
                   - x**k0 * sp.diff(Q, t) / (2 * f0 * Q)).subs(t, t0)
            lhs_J = quad_over_Y(sp.cancel(sig))
            err_J = abs(lhs_J - row_rhs(states.index('J0')))

            for name, err in (('L', err_L), ('J0', err_J)):
                worst = err if worst is None else max(worst, err)
                assert err < mpmath.mpf(10) ** (-(dps - 8)), \
                    '%s-row residual %s at dps %d' % (name, err, dps)
    print('  banked-row quadrature: L and J0 rows match recomputed '
          'quadrature, worst residual %s (dps 30/50)' % worst)


def test_extsys_numeric():
    """Banked system -> ExtSys: Horner A(x) vs direct sympy substitution at
    two precisions; singular points declared."""
    import mpmath
    banked = _load_banked()
    sysd = _banked_sysd(banked)
    es = ct.ExtSys(sysd)
    assert es.n == len(banked['states']) and es.var == 't'
    assert es.singular_points, 'no moving poles declared'
    pt = sp.Rational(7, 100)
    for dps in (30, 50):
        with mpmath.workdps(dps + 20):     # form the point AT working
            pt_mp = mpmath.mpf(7) / 100    # precision (ambient-15 footgun)
        rows = es.A(pt_mp, 0, dps)
        worst = 0.0
        for (i, j), _ in es._ent.items():
            ref = sp.N(sysd['Amat'][i, j].subs(t, pt), dps + 10)
            with mpmath.workdps(dps + 10):
                err = abs(rows[i][j] - mpmath.mpc(str(sp.re(ref)),
                                                  str(sp.im(ref))))
                rel = float(err / max(1, abs(mpmath.mpc(str(sp.re(ref)),
                                                        str(sp.im(ref))))))
            worst = max(worst, rel)
        assert worst < 10 ** (-(dps - 5)), \
            'ExtSys A() off by rel %g at dps %d' % (worst, dps)
    print('  ExtSys: A(7/100) matches exact substitution at dps 30/50; '
          '%d singular points declared' % len(es.singular_points))


if __name__ == '__main__':
    t_ = __import__('time').time()
    test_engine_truth_small_conic()
    test_refusal_deg3()
    test_truth_case_toy_rebuild()
    test_mutation_comparator_teeth()
    test_mutation_certificate_teeth()
    test_banked_rows_quadrature()
    test_extsys_numeric()
    print('test_conic_thirdkind: ALL LEGS PASS (%.1f s)'
          % (__import__('time').time() - t_))

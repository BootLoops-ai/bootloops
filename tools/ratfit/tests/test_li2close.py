# Tests for tools/li2close.py
# (a) known integrals (pi^2/6 partial-fraction pair, finite-interval atoms,
#     u-parametric family) with divergence-ledger behavior;
# (b) positive control: re-derive the reference (5,7,13) wall constants
#     K0(1/2), K1 from the archived atom lists (READ-ONLY) and match the shipped
#     kclosed.py evaluators to >= 45 digits.
import importlib.util
import os

import mpmath
import pytest
import sympy as sp
from sympy import Rational as R
from sympy import oo

import li2close as lc

KCLOSE = os.environ.get('RATFIT_KCLOSE_DIR', '')  # reference data dir (not shipped; env-pointed): unset/absent -> positive-control tests SKIP


def _agree_digits(a, b):
    a, b = mpmath.mpf(a), mpmath.mpf(b)
    if a == b:
        return mpmath.inf
    return float(-mpmath.log10(abs(a - b) / max(abs(a), mpmath.mpf('1e-30'))))


def _load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ------------------------------------------------------------ known integrals
def test_pi2_over_6_partial_fraction_pair():
    # int_0^oo ln(1+t)/(t(1+t)) dt = pi^2/6;  1/(t(1+t)) = 1/t - 1/(1+t):
    # each atom alone is ln^2(Lambda)-divergent; the pair cancels in the ledger.
    atoms = [lc.Atom(1, 0, 1, 1, 1, 0, oo, 'pf.t'),
             lc.Atom(-1, 1, 1, 1, 1, 0, oo, 'pf.1pt')]
    cl = lc.close_atoms(atoms)
    assert sp.simplify(cl.finite - sp.pi ** 2 / 6) == 0
    assert any(e['monom'] == 'LE^0 LL^2' for e in cl.ledger)
    assert all('FAIL' not in e['verdict'] for e in cl.ledger)
    f = lc.emitted_function(lc.emit_mpmath(cl.finite, 'val'), 'val')
    ref = mpmath.mp.pi ** 2 / 6
    with mpmath.workdps(60):
        assert _agree_digits(f(dps=50), mpmath.pi ** 2 / 6) >= 45


def test_single_divergent_atom_raises():
    with pytest.raises(lc.DivergenceError) as ei:
        lc.close_atoms([lc.Atom(1, 0, 1, 1, 1, 0, oo)])
    assert ei.value.ledger    # ledger names the surviving monomial


def test_finite_interval_atom_vs_quadrature():
    # int_0^1 ln(1+2t)/(1+t) dt (generic case; ReLi2 argument crosses 1)
    atoms = [lc.Atom(1, 1, 1, 1, 2, 0, 1, 'fin')]
    cl = lc.close_atoms(atoms)
    assert cl.ledger == []          # no endpoint divergences at all
    f = lc.emitted_function(lc.emit_mpmath(cl.finite, 'val'), 'val')
    with mpmath.workdps(60):
        ref = mpmath.quad(lambda t: mpmath.log(1 + 2 * t) / (1 + t), [0, 1])
        assert _agree_digits(f(dps=50), ref) >= 35


def test_param_family_ledger_and_eval():
    # K(u) = int_0^oo ln(1+u t)/(t(1+t)) dt, u in (0,1): u-parametric atoms,
    # symbolic LL^2 cancellation, sign-sampling on the domain for -u/(1-u).
    u = sp.symbols('u', positive=True)
    atoms = [lc.Atom(1, 0, 1, 1, u, 0, oo, 'p.t'),
             lc.Atom(-1, 1, 1, 1, u, 0, oo, 'p.1pt')]
    cl = lc.close_atoms(atoms, param=u, domain=(0, 1))
    assert all('FAIL' not in e['verdict'] for e in cl.ledger)
    src = lc.emit_mpmath(cl.finite, 'Kp', param=u)
    f = lc.emitted_function(src, 'Kp')
    for uq in ('1/3', '2/3'):
        with mpmath.workdps(60):
            uf = mpmath.mpf(sp.Rational(uq).p) / sp.Rational(uq).q
            ref = mpmath.quad(
                lambda t: mpmath.log(1 + uf * t) / (t * (1 + t)),
                [0, 1, 10, mpmath.inf])
            assert _agree_digits(f(uq, dps=50), ref) >= 25, uq


def test_param_without_domain_raises():
    u = sp.symbols('u', positive=True)
    with pytest.raises(ValueError):
        lc.close_atoms([lc.Atom(1, 1, 1, 1, u, 0, 1)], param=u)


def test_ilin_atoms_match_direct_integral():
    # _ilin structure: int_{lo}^{hi} dc/((A0+A1 c)(B0+B1 c)) as 4 log atoms
    t = sp.symbols('tau', positive=True)
    atoms = lc.ilin_atoms(1, t + 1, 2, 3, 1, 0, 1, t, (R(1, 2), 2), 'il')
    assert len(atoms) == 4
    # pointwise: sum of atom integrands == the closed _ilin at a few t values
    tv = R(3, 4)
    with mpmath.workdps(40):
        s = sum(sp.N(a.integrand(t).subs(t, tv), 30) for a in atoms)
        A0, A1, B0, B1 = float(tv + 1), 2.0, 3.0, 1.0
        direct = mpmath.quad(
            lambda c: 1 / ((A0 + A1 * c) * (B0 + B1 * c)), [0, 1])
        assert abs(float(s) - float(direct)) < 1e-12


# ------------------------------------------- positive control (reference K0/K1)
@pytest.mark.skipif(not os.path.isdir(KCLOSE),
                    reason='kclose reference data not present')
def test_positive_control_k0_k1_rederivation():
    """Re-derive K0(u), K1 from the archived atom lists (READ-ONLY import of
    assemble.py) and match the shipped kclosed.py to >= 45 digits."""
    asm = _load(os.path.join(KCLOSE, 'assemble.py'), 'kclose_assemble')
    kcl = _load(os.path.join(KCLOSE, 'kclosed.py'), 'kclose_kclosed')

    # K1 (pure constant, 28 atoms)
    atoms1, _m = asm.atoms_K1()
    cl1 = lc.close_atoms(atoms1)
    f1 = lc.emitted_function(lc.emit_mpmath(cl1.finite, 'K1_new'), 'K1_new')
    with mpmath.workdps(80):
        d = _agree_digits(f1(dps=70), kcl.K1_closed(70))
    assert d >= 45, f'K1 agreement only {d:.1f}d'

    # K0(u) (u-parametric, 18 atoms), evaluated at u = 1/2
    atoms0, _t = asm.atoms_K0()
    cl0 = lc.close_atoms(atoms0, param=asm.u, domain=(R(1, 5), 1))
    src = lc.emit_mpmath(cl0.finite, 'K0_new', param=asm.u)
    f0 = lc.emitted_function(src, 'K0_new')
    with mpmath.workdps(80):
        d = _agree_digits(f0('1/2', dps=70), kcl.K0_closed('1/2', 70))
    assert d >= 45, f'K0(1/2) agreement only {d:.1f}d'
    # ledger must show every divergent monomial cancelled
    for cl in (cl0, cl1):
        assert all('FAIL' not in e['verdict'] for e in cl.ledger)

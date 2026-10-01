"""Synthetic exact-data tests for tools/ratfit.py (production-close ports)."""
from fractions import Fraction

import flint
import pytest

import ratfit

# True function f = P/Q, P = 3 - 2x + x^2, Q = 1 + 5x^2 + 7x^3 (coprime,
# no rational roots, so Q != 0 at every rational sample node).
P_TRUE = flint.fmpq_poly([3, -2, 1])
Q_TRUE = flint.fmpq_poly([1, 0, 5, 7])
N_NODES = 24
XS = [Fraction(k, 97) for k in range(1, N_NODES + 1)]


def frac(v):
    return Fraction(int(v.p), int(v.q))


def feval(x):
    xq = ratfit.fq(x)
    return frac(P_TRUE(xq) / Q_TRUE(xq))


YS = [feval(x) for x in XS]
FRESH = [Fraction(5, 11), Fraction(-3, 13), Fraction(101, 103)]


def assert_matches_truth(P, Q):
    for x in FRESH:
        xq = ratfit.fq(x)
        assert P(xq) / Q(xq) == ratfit.fq(feval(x))


def test_thiele_blind_roundtrip():
    P, Q, ok, deg = ratfit.thiele_loo(XS, YS)
    assert ok and deg == 3
    assert_matches_truth(P, Q)


def test_thiele_screened_agrees_on_success():
    P1, Q1, ok1, d1 = ratfit.thiele_loo(XS, YS)
    P2, Q2, ok2, d2 = ratfit.thiele_loo_screened(XS, YS)
    assert ok1 and ok2 and d1 == d2
    assert_matches_truth(P2, Q2)


def test_thiele_screened_agrees_on_failure():
    ys = list(YS)
    ys[2] += 1  # node index 2 is a held-out test node for n=24, n_loo=6
    P1, Q1, ok1, d1 = ratfit.thiele_loo(XS, ys)
    P2, Q2, ok2, d2 = ratfit.thiele_loo_screened(XS, ys)
    assert not ok1 and not ok2
    # both must fall back to the SAME single full-node exact Thiele
    assert str(P1) == str(P2) and str(Q1) == str(Q2) and d1 == d2


def test_known_denominator_roundtrip():
    Pn, Qd, ok, deg = ratfit.ratrecon_with_denom(XS, YS, Q_TRUE)
    assert ok and deg == 3
    assert str(Pn) == str(P_TRUE) and str(Qd) == str(Q_TRUE)


def test_known_denominator_gcd_reduces_multiple():
    Qbig = Q_TRUE * flint.fmpq_poly([-2, 1])  # extra factor (x-2)
    Pn, Qd, ok, _ = ratfit.ratrecon_with_denom(XS, YS, Qbig)
    assert ok
    assert str(Pn) == str(P_TRUE) and str(Qd) == str(Q_TRUE)


def test_screen_accepts_true_Q_rejects_wrong_Q():
    assert ratfit.screen_modp(XS, YS, Q_TRUE) is True
    Qwrong = flint.fmpq_poly([2, 3, 1, 1])  # cubic, not a multiple of Q_TRUE
    assert ratfit.screen_modp(XS, YS, Qwrong) is False
    # screen verdicts must agree with the exact certify
    assert ratfit.ratrecon_with_denom(XS, YS, Qwrong)[2] is False


def _coeff_strs(poly):
    return [str(c) for c in poly.coeffs()]


Q0 = flint.fmpq_poly([1, 2])              # 1 + 2x
R = flint.fmpq_poly([1, 0, 1])            # 1 + x^2 (the stair ratio, deg 2)
QEX = {(5, 0): {0: _coeff_strs(Q0), 1: _coeff_strs(Q0 * R),
                "2": _coeff_strs(Q0 * R * R)},   # str key: json-loaded form
       (5, 1): {0: _coeff_strs(Q0)}}


def test_extract_stairs_recovers_ratio():
    Rlib = ratfit.extract_stairs(QEX)
    assert list(Rlib) == [5]
    assert len(Rlib[5]) == 1  # two consecutive pairs, one distinct ratio
    assert str(Rlib[5][0]) == str(R)


def test_row_lcm():
    L = ratfit.row_lcm(QEX)[5]
    top = Q0 * R * R
    assert L.degree() == top.degree() == 5
    for q in (Q0, Q0 * R, top):
        assert ratfit.exact_div(L, q) is not None


def test_degree_budget_sanity():
    bud = ratfit.degree_budget(QEX, n=40, n_loo=6)
    b5 = bud[5]
    assert b5["blind_ceiling"] == 17
    assert b5["certify_cap"] == 32
    assert b5["stair_step_degs"] == [2]
    assert b5["projected_degQ"][0] == 1
    assert b5["projected_degQ"][1] == 3
    assert b5["projected_degQ"][2] == 5
    assert b5["projected_degQ"][3] == 7  # projected: +max stair deg per step
    assert b5["projected_degQ"][6] == 13


def test_zero_function():
    zs = [Fraction(0)] * N_NODES
    P, Q, ok, _ = ratfit.thiele_loo(XS, zs)
    assert ok and P(ratfit.fq(Fraction(1, 2))) == 0
    P2, Q2, ok2, _ = ratfit.thiele_loo_screened(XS, zs)
    assert ok2 and P2(ratfit.fq(Fraction(1, 2))) == 0


# ===== parametric_rec_interp: multivariate iterated-Newton lift =====

# Truth: normalized recurrence coefficients that are exact polynomials in
# (a, b), hidden behind a param-dependent overall scale s(a,b) != 0 that the
# normalize_key coefficient carries (the fit must divide it out).
T1 = {(0, 0): Fraction(3), (1, 0): Fraction(2), (0, 1): Fraction(-1),
      (1, 1): Fraction(1)}                       # 3 + 2a - b + ab
T2 = {(2, 0): Fraction(1), (0, 2): Fraction(-5)}  # a^2 - 5b^2


def _scaled_rec(params):
    a, b = (Fraction(p) for p in params)
    s = 1 + a * a + b * b  # nonzero scale
    return {(0, 0): s,
            (1, 0): ratfit.mpoly_eval(T1, params) * s,
            (1, 1): ratfit.mpoly_eval(T2, params) * s}


GRID2 = [(Fraction(i), Fraction(j)) for i in range(3) for j in range(3)]


def test_parametric_rec_interp_roundtrip():
    fit = ratfit.parametric_rec_interp(_scaled_rec, GRID2, (0, 0))
    assert fit[(0, 0)] == {(0, 0): Fraction(1)}   # normalizer maps to 1
    assert fit[(1, 0)] == T1
    assert fit[(1, 1)] == T2
    fresh = (Fraction(7, 5), Fraction(-3, 11))    # off-grid evaluation
    assert ratfit.mpoly_eval(fit[(1, 0)], fresh) == ratfit.mpoly_eval(T1, fresh)


def test_parametric_rec_interp_scalar_grid():
    def rec(params):
        (a,) = (Fraction(p) for p in params)
        return {(0, 0): Fraction(1), (2, 1): a * a - 3}
    fit = ratfit.parametric_rec_interp(rec, [Fraction(k) for k in range(3)],
                                       (0, 0))
    assert fit[(2, 1)] == {(2,): Fraction(1), (0,): Fraction(-3)}


def test_parametric_rec_interp_crosscheck_catches_underresolved():
    # a^3 dependence on a 3-node axis: the grid fit exists and matches every
    # grid node, so ONLY the off-grid cross-check can catch it — must raise.
    def rec(params):
        (a,) = (Fraction(p) for p in params)
        return {(0, 0): Fraction(1), (1, 0): a ** 3}
    with pytest.raises(ValueError, match="cross-check FAILED"):
        ratfit.parametric_rec_interp(rec, [Fraction(k) for k in range(3)],
                                     (0, 0))


def test_parametric_rec_interp_zero_normalizer_raises():
    def rec(params):
        a, b = (Fraction(p) for p in params)
        return {(0, 0): a, (1, 0): a + b}  # normalizer vanishes at a=0
    with pytest.raises(ValueError, match="cannot normalize"):
        ratfit.parametric_rec_interp(rec, GRID2, (0, 0))


def test_parametric_rec_interp_incomplete_grid_raises():
    with pytest.raises(ValueError, match="tensor-product"):
        ratfit.parametric_rec_interp(_scaled_rec, GRID2[:-1], (0, 0))


def test_parametric_rec_interp_deg_bound_guard():
    def rec(params):
        (a,) = (Fraction(p) for p in params)
        return {(0, 0): Fraction(1), (1, 0): a}
    with pytest.raises(ValueError, match="deg_bound"):
        ratfit.parametric_rec_interp(rec, [Fraction(k) for k in range(14)],
                                     (0, 0), deg_bound=12)

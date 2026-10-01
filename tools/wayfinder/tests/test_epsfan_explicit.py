"""Tests for the additive epsfan extensions (explicit-node wiring):
explicit-node grids (eps_nodes1/eps_nodes2) + k_certify subwindow in
vandermonde_laurent_certified.  Legacy behavior covered by
test_epsfan_certified.py (unchanged)."""
import os
import sys

import mpmath as mp
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from epsfan import vandermonde_laurent_certified, EpsFanCertifyError  # noqa: E402


def f_vec(e):
    return [3 / e ** 2 - 1 / e + 2 + 5 * e, 7 / e ** 2 + mp.mpf('0.25') * e]


def _primes(lo, n, skip=()):
    out, x = [], lo

    def isp(m):
        if m < 2:
            return False
        i = 2
        while i * i <= m:
            if m % i == 0:
                return False
            i += 1
        return True

    while len(out) < n:
        if isp(x) and x not in skip:
            out.append(x)
        x += 1
    return out


Q1 = _primes(1013, 24)
Q2 = _primes(1091, 36, skip=set(Q1))
N1 = ['1/%d' % q for q in Q1]
N2 = ['1/%d' % q for q in Q2]


def test_explicit_nodes_known_fan():
    res = vandermonde_laurent_certified(f_vec, -2, 1, 30, R0='1/4',
                                        eps_nodes1=N1, eps_nodes2=N2,
                                        guard=8)
    with mp.workdps(30):
        c = res['coeffs']
        assert abs(c[-2][0] - 3) < mp.mpf('1e-25')
        assert abs(c[1][1] - mp.mpf('0.25')) < mp.mpf('1e-25')
    assert res['diag']['explicit_nodes'] is True
    assert res['diag']['n_grid1'] == 24 and res['diag']['n_grid2'] == 36


def test_shared_node_refused():
    with pytest.raises(ValueError, match='SHARE nodes'):
        vandermonde_laurent_certified(f_vec, -2, 1, 30, R0='1/4',
                                      eps_nodes1=N1,
                                      eps_nodes2=N1[:1] + N2[1:], guard=8)


def test_node_count_rule_refused():
    with pytest.raises(ValueError, match='node-count rule'):
        vandermonde_laurent_certified(f_vec, -2, 1, 30, R0='1/4',
                                      eps_nodes1=N1[:8], eps_nodes2=N2,
                                      guard=8)


def test_r0_clamp_refused():
    with pytest.raises(ValueError, match='R0/10'):
        vandermonde_laurent_certified(f_vec, -2, 1, 30, R0='1/500',
                                      eps_nodes1=N1, eps_nodes2=N2, guard=8)


def test_both_or_neither():
    with pytest.raises(ValueError, match='together'):
        vandermonde_laurent_certified(f_vec, -2, 1, 30, R0='1/4',
                                      eps_nodes1=N1, guard=8)


def test_grid2_corruption_raises():
    q1set = set(Q1)

    def router(e):
        q = int(round(1 / float(e)))
        v = f_vec(e)
        if q in q1set:
            return v
        return [v[0] * (1 + mp.mpf('1e-20')), v[1]]

    with pytest.raises(EpsFanCertifyError) as exc:
        vandermonde_laurent_certified(router, -2, 1, 30, R0='1/4',
                                      eps_nodes1=N1, eps_nodes2=N2, guard=8)
    assert exc.value.kind == 'vand-grid'


def test_k_certify_subwindow():
    q1set = set(Q1)

    def router(e):
        q = int(round(1 / float(e)))
        c3 = mp.mpf('1e-6') if q in q1set else mp.mpf('2e-6')
        return [3 / e ** 2 + 2 + c3 * e ** 3]

    # gate confined to (-2,1): k=3 disagreement reported but not raising
    res = vandermonde_laurent_certified(router, -2, 3, 12, R0='1/4',
                                        eps_nodes1=N1, eps_nodes2=N2,
                                        guard=8, k_certify=(-2, 1))
    assert res['diag']['k_certify'] == [-2, 1]
    assert res['diag']['grid_agree_digits'][3] < 10
    assert res['diag']['certify_min_agree_digits'] >= 12
    # full-window gate raises on the same data
    with pytest.raises(EpsFanCertifyError) as exc:
        vandermonde_laurent_certified(router, -2, 3, 12, R0='1/4',
                                      eps_nodes1=N1, eps_nodes2=N2, guard=8)
    assert exc.value.kind == 'vand-grid'


def test_k_certify_bounds_checked():
    with pytest.raises(ValueError, match='k_certify'):
        vandermonde_laurent_certified(f_vec, -2, 1, 30, R0='1/4',
                                      eps_nodes1=N1, eps_nodes2=N2,
                                      guard=8, k_certify=(-4, 1))

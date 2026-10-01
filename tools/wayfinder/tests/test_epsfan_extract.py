"""Extractor-mode tests (vandermonde_laurent_extract).

Charter: synthetic Laurent recovery with ON-DEMAND nodes at 3 dps values;
cond/wp machinery unchanged (delegation equality vs a direct certified
call on the same grids); f(dps) node-count rule enforced; on-demand
sampling actually happens through the caller's transport callback;
fail-closed refusals (no R0, under-supplying node_gen, bar-lowering
n1_request, corrupted sample -> certified-layer raise).
"""
import math
from fractions import Fraction

import mpmath as mp
import pytest

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from epsfan import (vandermonde_laurent_extract, vandermonde_laurent_certified,
                    extract_node_counts, EpsFanCertifyError)

KMIN, KMAX = -2, 3
P = KMAX - KMIN + 1


def _coeff_truth(dps):
    """Known irrational Laurent window (window-limited: probe measures pure
    extraction error)."""
    with mp.workdps(dps + 30):
        return {-2: +mp.pi, -1: -mp.sqrt(2), 0: mp.mpf(3) / 7,
                1: +mp.log(2), 2: -mp.e, 3: mp.zeta(3)}


def _make_f(dps, counter=None, corrupt_node=None):
    truth = _coeff_truth(dps + 60)

    def f(e):
        if counter is not None:
            counter.append(mp.mpf(e))
        v = mp.mpf(0)
        for k, c in truth.items():
            v += c * mp.mpf(e) ** k
        if corrupt_node is not None and abs(mp.mpf(e) - mp.mpf(corrupt_node)) \
                < mp.mpf('1e-30'):
            v += mp.mpf('1e-25') * max(abs(v), mp.mpf(1))
        return v
    return f


def _prime_node_gen(issued_log=None):
    """Caller node design: exact 1/q strings, q prime, geometric-ish band
    1013..~40000, disjoint from exclude (the production pattern)."""
    def is_prime(n):
        if n < 2:
            return False
        i = 2
        while i * i <= n:
            if n % i == 0:
                return False
            i += 1
        return True

    def gen(n, exclude):
        if issued_log is not None:
            issued_log.append((n, tuple(exclude)))
        ex = set(str(Fraction(x)) for x in exclude)
        out, q = [], 1013
        while len(out) < n:
            if is_prime(q):
                s = '1/%d' % q
                if str(Fraction(s)) not in ex:
                    out.append(s)
            q += 2 if q > 2 else 1
        return out
    return gen


@pytest.mark.parametrize('dps', [30, 120, 250])
def test_extract_synthetic_recovery_on_demand(dps):
    """Synthetic Laurent recovery with on-demand nodes at 3 dps values;
    node count follows the f(dps) rule; sampling goes through the
    callback."""
    seed = 6                     # small seed so the f(dps) term can bind
    n1_want, n2_want = extract_node_counts(KMIN, KMAX, dps, nodes_seed=seed)
    assert n1_want == max(seed, P + 4, P + math.ceil(dps / 25.0))
    counter, issued = [], []
    f = _make_f(dps, counter=counter)
    res = vandermonde_laurent_extract(
        f, _prime_node_gen(issued), KMIN, KMAX, dps,
        R0=Fraction(1, 4), nodes_seed=seed)
    d = res['diag']
    assert d['extractor'] is True
    assert d['n_demanded'] == [n1_want, n2_want]
    assert d['n_grid1'] == n1_want and d['n_grid2'] == n2_want
    # engine issued exactly two node requests: (n1, ()) then (n2, g1)
    assert issued[0][0] == n1_want and issued[0][1] == ()
    assert issued[1][0] == n2_want and len(issued[1][1]) == n1_want
    # on-demand sampling: every evaluation went through the callback;
    # each pass costs n1+n2 evals, wp escalation re-calls in multiples
    assert len(counter) >= n1_want + n2_want
    assert len(counter) % (n1_want + n2_want) == 0
    truth = _coeff_truth(dps + 60)
    with mp.workdps(dps + 10):
        for k in range(KMIN, KMAX + 1):
            err = abs(res['coeffs'][k] - truth[k]) / max(abs(truth[k]),
                                                         mp.mpf(1))
            assert err < mp.mpf(10) ** (-(dps - 2)), (k, dps, err)


def test_extract_node_count_grows_with_dps():
    """The f(dps) rule makes precision dps-dialable: node demand grows."""
    ns = [extract_node_counts(KMIN, KMAX, d, nodes_seed=6)[0]
          for d in (30, 120, 250)]
    assert ns == sorted(ns) and ns[0] < ns[-1]
    # production seed: 24-node-equivalent config at dps 40
    assert extract_node_counts(-4, 15, 40, nodes_seed=24) == (24, 36)


def test_extract_delegation_machinery_unchanged():
    """cond/wp machinery unchanged BY CONSTRUCTION: extractor result equals
    a direct vandermonde_laurent_certified call on the same grids,
    coefficient-for-coefficient, and carries the same cond/wp diag."""
    dps = 40
    gen = _prime_node_gen()
    n1, n2 = extract_node_counts(KMIN, KMAX, dps, nodes_seed=6)
    g1 = gen(n1, ())
    g2 = gen(n2, tuple(g1))
    f = _make_f(dps)
    r_x = vandermonde_laurent_extract(f, gen, KMIN, KMAX, dps,
                                      R0='1/4', nodes_seed=6)
    r_c = vandermonde_laurent_certified(f, KMIN, KMAX, dps, R0='1/4',
                                        nodes_seed=6, eps_nodes1=g1,
                                        eps_nodes2=g2)
    assert r_x['diag']['wp'] == r_c['diag']['wp']
    assert r_x['diag']['cond'] == r_c['diag']['cond']
    with mp.workdps(dps):
        for k in range(KMIN, KMAX + 1):
            assert r_x['coeffs'][k] == r_c['coeffs'][k]


def test_extract_refuses_without_R0():
    with pytest.raises(ValueError, match='REFUSING to assume'):
        vandermonde_laurent_extract(_make_f(30), _prime_node_gen(),
                                    KMIN, KMAX, 30, nodes_seed=6)


def test_extract_refuses_undersupplying_node_gen():
    def stingy(n, exclude):
        return _prime_node_gen()(max(n - 2, 1), exclude)
    with pytest.raises(ValueError, match='DEMANDS n1'):
        vandermonde_laurent_extract(_make_f(30), stingy, KMIN, KMAX, 30,
                                    R0='1/4', nodes_seed=6)


def test_extract_refuses_bar_lowering_n1_request():
    n1, _ = extract_node_counts(KMIN, KMAX, 30, nodes_seed=6)
    with pytest.raises(ValueError, match='BELOW the'):
        vandermonde_laurent_extract(_make_f(30), _prime_node_gen(),
                                    KMIN, KMAX, 30, R0='1/4', nodes_seed=6,
                                    n1_request=n1 - 1)


def test_extract_n1_request_dial_up():
    """n1_request above the rule is honored and n2 scales with it."""
    dps = 30
    n1_rule, _ = extract_node_counts(KMIN, KMAX, dps, nodes_seed=6)
    issued = []
    res = vandermonde_laurent_extract(
        _make_f(dps), _prime_node_gen(issued), KMIN, KMAX, dps, R0='1/4',
        nodes_seed=6, n1_request=n1_rule + 3)
    assert res['diag']['n_demanded'] == [n1_rule + 3,
                                         math.ceil(1.5 * (n1_rule + 3))]


def test_extract_sabotage_corrupted_sample_raises():
    """A 1e-25-relative corruption at ONE grid-2 node must trip the
    certified layer (power depth: 25 d perturbation vs two-grid tol at
    10^-(dps+guard//2); single-depth measurement)."""
    dps = 30
    gen = _prime_node_gen()
    n1, n2 = extract_node_counts(KMIN, KMAX, dps, nodes_seed=6)
    g2 = gen(n2, tuple(gen(n1, ())))
    f_bad = _make_f(dps, corrupt_node=g2[3])
    with pytest.raises(EpsFanCertifyError) as ei:
        vandermonde_laurent_extract(f_bad, gen, KMIN, KMAX, dps,
                                    R0='1/4', nodes_seed=6)
    assert ei.value.kind in ('vand-grid', 'vand-verify')

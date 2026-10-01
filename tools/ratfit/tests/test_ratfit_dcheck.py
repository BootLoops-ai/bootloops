"""The fit-point D check in ratfit: synthetic exact-data legs.

A fit must refuse a D that vanishes at a FIT sample, not only at the held-out
nodes -- otherwise that sample is absorbed.  Blind Thiele routes already verify
every fit node with Q != 0 (thiele_loo, thiele_loo_screened, thiele_gate.fit:
no hole -- the controls below pin that); the known-Q certify
`ratrecon_with_denom` certified the held-out nodes only, so a node on a root of
Qknown (y*Qknown = 0 there, whatever y) would be absorbed without the check.
The PLANTED legs fail by name ("PLANTED SAMPLE ABSORBED") if the check is
removed and pass with it (ok=False).
"""
from fractions import Fraction

import flint
import pytest

import ratfit
from ratfit import thiele_gate

P_TRUE = flint.fmpq_poly([3, -2, 1])
Q_TRUE = flint.fmpq_poly([1, 0, 5, 7])
N_NODES = 24
XS = [Fraction(k, 97) for k in range(1, N_NODES + 1)]
PLANT = Fraction(1, 7)


def frac(v):
    return Fraction(int(v.p), int(v.q))


def feval(x):
    xq = ratfit.fq(x)
    return frac(P_TRUE(xq) / Q_TRUE(xq))


YS = [feval(x) for x in XS]
TEST_IDX, FIT_IDX = ratfit._loo_split(N_NODES, 6)


def q_with_root_at(k):
    """Q_TRUE times (x - x_k): a known denominator vanishing at node k."""
    return Q_TRUE * flint.fmpq_poly([-ratfit.fq(XS[k]), 1])


def test_known_denominator_positive_control():
    Pn, Qd, ok, deg = ratfit.ratrecon_with_denom(XS, YS, Q_TRUE)
    assert ok and deg == 3
    assert str(Pn) == str(P_TRUE) and str(Qd) == str(Q_TRUE)


@pytest.mark.parametrize("which", ["fit-node", "held-out-node"])
def test_known_denominator_planted_sample_on_a_Q_root_refused(which):
    k = FIT_IDX[2] if which == "fit-node" else TEST_IDX[1]
    ys = list(YS)
    ys[k] = ys[k] + PLANT
    Pn, Qd, ok, deg = ratfit.ratrecon_with_denom(XS, ys, q_with_root_at(k))
    if ok:
        xq = ratfit.fq(XS[k])
        got = "Q=0" if Qd(xq) == 0 else str(frac(Pn(xq) / Qd(xq)) == ys[k])
        pytest.fail("PLANTED SAMPLE ABSORBED (%s %d): ratrecon_with_denom certified ok with Qknown vanishing at "
                    "that node; returned P/Q = %s / %s reproduces the planted sample: %s"
                    % (which, k, Pn, Qd, got))
    assert not ok


def test_known_denominator_root_off_the_grid_still_reduces():
    # a Qknown root at a non-node is legitimate (gcd-reduced), unchanged
    Qbig = Q_TRUE * flint.fmpq_poly([-2, 1])
    Pn, Qd, ok, _ = ratfit.ratrecon_with_denom(XS, YS, Qbig)
    assert ok and str(Pn) == str(P_TRUE) and str(Qd) == str(Q_TRUE)


def test_blind_routes_refuse_the_planted_sample_controls():
    k = FIT_IDX[2]
    ys = list(YS)
    ys[k] = ys[k] + PLANT
    assert ratfit.thiele_loo(XS, ys)[2] is False
    assert ratfit.thiele_loo_screened(XS, ys)[2] is False
    assert thiele_gate.fit(XS, ys, max_fit=18) is None
    zs = [Fraction(0)] * N_NODES
    zs[k] = PLANT
    assert ratfit.thiele_loo(XS, zs)[2] is False
    assert thiele_gate.fit(XS, zs, max_fit=18) is None

#!/usr/bin/env python3
"""
test_quad.py — unit self-checks for wayfinder/quad.py (E3 refine-until-gate
adaptive quadrature).

Discipline: NO fabricated reference digits — every reference is an mpmath
closed form computed in-test at higher working precision (dps+40). Controls
required by the E3 specification:
  * endpoint-singular:  int_0^1 log(x) dx      = -1
  * oscillatory:        int_0^1 cos(50 x) dx   = sin(50)/50
  * smooth:             int_0^1 exp(x) dx      = e - 1
each at dps 30 / 60 / 100, asserting (a) >= dps matched digits vs exact,
(b) the certificate holds: |value - exact| <= err_bound, (c) the gated raw
agreement beat tol = 10^-(dps+guard).
Fail-closed: impossible tol at a tiny depth cap -> QuadNonConvergence with
named diagnostics; divergent 1/x -> raise; nan integrand -> raise. Knobs
(depth0/max_depth) are speed seeds only — pinned by a value-invariance test.
Crank: dps 30 vs 70 (D+40) agree to >= 30 digits with genuinely different
internal depths. Global mp.dps must be untouched (import-dps footgun check).

Runnable directly (python3 tests/test_quad.py) or under pytest.
Runtime target: seconds (Gauss-Legendre smoke kept at dps 30: its node
computation is the expensive part and it is opt-in for consumers).
"""

import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from mpmath import mp, mpf, mpc  # noqa: E402
from quad import quad_refine, QuadCert, QuadNonConvergence  # noqa: E402

DPS_LIST = (30, 60, 100)


def _digits(a, b, wp):
    """Matched significant digits of a vs reference b (computed at wp) —
    same metric as test_transport.py / gate.py (rel to max(|b|,1))."""
    with mp.workdps(wp):
        a = mpc(a)
        b = mpc(b)
        if a == b:
            return wp
        rel = abs(a - b) / max(abs(b), mpf(1))
        if rel == 0:
            return wp
        return int(-mp.log10(rel))


def _check_control(f, exact_fn, dps, name):
    """Run one control at one dps; assert digits, certificate, gate."""
    cert = quad_refine(f, [0, 1], dps, full_output=True)
    assert isinstance(cert, QuadCert)
    wp = dps + 40
    with mp.workdps(wp):
        exact = exact_fn()
        err = abs(mpc(cert.value) - exact)
        # (a) >= dps correct digits
        d = _digits(cert.value, exact, wp)
        assert d >= dps, "%s dps=%d: only %d digits" % (name, dps, d)
        # (b) certified bound holds vs exact
        assert err <= cert.err_bound, \
            "%s dps=%d: |err|=%s > err_bound=%s" % (
                name, dps, mp.nstr(err, 5), mp.nstr(cert.err_bound, 5))
        # (c) the raw double-refinement agreement beat the gate tol
        assert cert.agreement <= cert.tol, \
            "%s dps=%d: agreement above tol" % (name, dps)
    # certificate line is printable and non-empty
    assert cert.bound_line(name)
    return cert


def test_smooth_exp():
    for dps in DPS_LIST:
        _check_control(lambda x: mp.exp(x), lambda: mp.e - 1, dps, "exp")


def test_endpoint_singular_log():
    # tanh-sinh must eat the x=0 log endpoint singularity
    for dps in DPS_LIST:
        _check_control(lambda x: mp.log(x), lambda: mpf(-1), dps, "log")


def test_oscillatory_cos50():
    for dps in DPS_LIST:
        _check_control(lambda x: mp.cos(50 * x),
                       lambda: mp.sin(50) / 50, dps, "cos50")


def test_fail_closed_cap():
    # impossible tol at a tiny depth cap -> named fail-closed raise
    try:
        quad_refine(lambda x: mp.cos(50 * x), [0, 1], 60,
                    depth0=2, max_depth=3)
    except QuadNonConvergence as e:
        assert e.kind == "cap"
        assert e.depth == 3 and e.cap == 3
        assert e.agreement is not None and e.agreement > e.tol
        msg = str(e)
        for token in ("fail-closed", "agreement", "tol", "cap", "dps=60"):
            assert token in msg, "diagnostic %r missing from message" % token
    else:
        raise AssertionError("cap cripple did not raise")


def test_fail_closed_divergent():
    # 1/x on [0,1] diverges: the agreement never settles -> cap raise
    # (never a silently truncated value)
    try:
        quad_refine(lambda x: 1 / x, [0, 1], 30)
    except QuadNonConvergence as e:
        assert e.kind in ("cap", "nonfinite")
    else:
        raise AssertionError("divergent integral did not raise")


def test_fail_closed_nonfinite():
    try:
        quad_refine(lambda x: mp.nan, [0, 1], 30)
    except QuadNonConvergence as e:
        assert e.kind == "nonfinite"
    else:
        raise AssertionError("nan integrand did not raise")


def test_knobs_are_seeds_only():
    # depth0 / max_depth change speed, never the certified value
    f = lambda x: mp.cos(50 * x)  # noqa: E731
    v_default = quad_refine(f, [0, 1], 40)
    v_seeded = quad_refine(f, [0, 1], 40, depth0=4, max_depth=20)
    d = _digits(v_default, v_seeded, 80)
    assert d >= 40, "seed knob moved the value (agree only %d digits)" % d


def test_crank_dps30_vs_70():
    # charter point 3: D and D+40 agree to D with genuinely different
    # internal refinement work
    f = lambda x: mp.cos(50 * x)  # noqa: E731
    c30 = quad_refine(f, [0, 1], 30, full_output=True)
    c70 = quad_refine(f, [0, 1], 70, full_output=True)
    assert (c30.depths != c70.depths) or (c30.evals != c70.evals), \
        "crank runs did identical internal work"
    with mp.workdps(110):
        exact = mp.sin(50) / 50
    d = _digits(c30.value, c70.value, 110)
    assert d >= 30, "crank agreement only %d digits" % d
    assert _digits(c70.value, exact, 110) >= 70


def test_waypoints_and_bound_sum():
    # explicit waypoint split: same value, per-segment agreements sum
    c = quad_refine(lambda x: mp.log(x), [0, mpf(1) / 3, 1], 40,
                    full_output=True)
    assert len(c.segments) == 2 and len(c.depths) == 2
    with mp.workdps(80):
        err = abs(mpc(c.value) - mpf(-1))
        assert err <= c.err_bound
    assert _digits(c.value, mpf(-1), 80) >= 40


def test_infinite_intervals():
    c = quad_refine(lambda x: mp.exp(-x), [0, mp.inf], 40, full_output=True)
    with mp.workdps(80):
        assert abs(mpc(c.value) - 1) <= c.err_bound
    assert _digits(c.value, mpf(1), 80) >= 40
    c = quad_refine(lambda x: mp.exp(-x * x), [mp.ninf, mp.inf], 40,
                    full_output=True)
    with mp.workdps(80):
        ref = mp.sqrt(mp.pi)
        assert abs(mpc(c.value) - ref) <= c.err_bound
        assert _digits(c.value, ref, 80) >= 40


def test_complex_segment():
    c = quad_refine(lambda z: mp.exp(z), [0, mpc(1, 1)], 40,
                    full_output=True)
    with mp.workdps(80):
        ref = mp.exp(mpc(1, 1)) - 1
        assert abs(mpc(c.value) - ref) <= c.err_bound
        assert _digits(c.value, ref, 80) >= 40


def test_gauss_legendre_smoke():
    # GL is opt-in (expensive nodes); dps 30 smoke on the smooth control
    _t = time.time()
    cert = quad_refine(lambda x: mp.exp(x), [0, 1], 30, method="gl",
                       full_output=True)
    assert cert.method == "gauss-legendre"
    with mp.workdps(70):
        ref = mp.e - 1
        assert abs(mpc(cert.value) - ref) <= cert.err_bound
    assert _digits(cert.value, ref, 70) >= 30
    assert time.time() - _t < 60


def test_bad_arguments():
    for bad in (dict(method="simpson"),
                dict(depth0=8, max_depth=4)):
        try:
            quad_refine(lambda x: x, [0, 1], 30, **bad)
        except ValueError:
            pass
        else:
            raise AssertionError("bad args %r accepted" % (bad,))
    try:
        quad_refine(lambda x: x, [1], 30)
    except ValueError:
        pass
    else:
        raise AssertionError("1-point interval accepted")
    try:
        quad_refine(lambda x: x, [0, 1], 0)
    except ValueError:
        pass
    else:
        raise AssertionError("dps=0 accepted")


def test_global_dps_untouched():
    before = mp.dps
    quad_refine(lambda x: mp.log(x), [0, 1], 45)
    try:
        quad_refine(lambda x: mp.cos(50 * x), [0, 1], 45, depth0=2,
                    max_depth=3)
    except QuadNonConvergence:
        pass
    assert mp.dps == before, "global mp.dps was mutated"


ALL_TESTS = [
    test_smooth_exp,
    test_endpoint_singular_log,
    test_oscillatory_cos50,
    test_fail_closed_cap,
    test_fail_closed_divergent,
    test_fail_closed_nonfinite,
    test_knobs_are_seeds_only,
    test_crank_dps30_vs_70,
    test_waypoints_and_bound_sum,
    test_infinite_intervals,
    test_complex_segment,
    test_gauss_legendre_smoke,
    test_bad_arguments,
    test_global_dps_untouched,
]


if __name__ == "__main__":
    n_fail = 0
    for t in ALL_TESTS:
        t0 = time.time()
        try:
            t()
            print("PASS  %-32s %.2fs" % (t.__name__, time.time() - t0))
        except AssertionError as e:
            n_fail += 1
            print("FAIL  %-32s %s" % (t.__name__, e))
    print("%d/%d passed" % (len(ALL_TESTS) - n_fail, len(ALL_TESTS)))
    sys.exit(1 if n_fail else 0)

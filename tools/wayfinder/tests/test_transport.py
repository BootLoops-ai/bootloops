#!/usr/bin/env python3
"""
test_transport.py — unit self-checks for wayfinder/transport.py.

Discipline: NO fabricated reference digits — every reference is an mpmath
closed form computed in-test at higher working precision. Two-precision
agreement (dps 80 vs 110) is asserted, and global mp.dps must be untouched
(import-dps footgun check).

Runnable directly (python3 tests/test_transport.py) or under pytest.
Runtime target: seconds (Verify phase owns the heavy controls).
"""

import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from mpmath import mp, mpf, mpc  # noqa: E402
from transport import transport_fixed_eps  # noqa: E402


def _agree_digits(a, b, wp):
    """Matched significant digits of a vs reference b (computed at wp)."""
    with mp.workdps(wp):
        a = mpc(a)
        b = mpc(b)
        if a == b:
            return wp
        rel = abs(a - b) / max(abs(b), mpf(1))
        if rel == 0:
            return wp
        return int(-mp.log10(rel))


class ScalarPow:
    """y' = (a/x) y with a = 3/2 + eps; solution y(x1) = y0 (x1/x0)^a."""
    n = 1
    var = "x"
    meta = {"singular_points": ["0"]}

    def A(self, x, eps, dps):
        with mp.workdps(dps):
            return [[(mpf(3) / 2 + eps) / mpc(x)]]


class ConstRot:
    """y' = [[0,1],[-1,0]] y; y(x) = [cos dx, -sin dx] from y0=[1,0]."""
    n = 2
    var = "x"
    meta = {}

    def A(self, x, eps, dps):
        return [[mpc(0), mpc(1)], [mpc(-1), mpc(0)]]


class SimplePole:
    """y' = y/(x-1); y(x) = C (x-1). Analytic continuation around the pole
    is single-valued, so the detour side does not matter."""
    n = 1
    var = "x"
    meta = {"singular_points": ["1"]}

    def A(self, x, eps, dps):
        with mp.workdps(dps):
            return [[1 / (mpc(x) - 1)]]


def test_scalar_power_law_two_precision():
    """Closed form (x1/x0)^a at dps 80 and 110 + cross-precision agreement."""
    dps_saved = mp.dps
    results = {}
    for dps in (80, 110):
        y = transport_fixed_eps(ScalarPow(), "0.1", "1", "3", ["1"], dps)
        results[dps] = y[0]
        with mp.workdps(dps + 20):
            a = mpf(3) / 2 + mpf(1) / 10
            ref = mp.power(mpf(3), a)
            d = _agree_digits(y[0], ref, dps + 20)
        assert d >= dps - 6, f"scalar dps={dps}: only {d} digits vs closed form"
    d2 = _agree_digits(results[80], results[110], 130)
    assert d2 >= 72, f"two-precision agreement only {d2} digits"
    assert mp.dps == dps_saved, "global mp.dps was mutated (import-dps footgun)"
    print(f"  scalar power law: dps80/dps110 vs closed form OK, "
          f"two-precision {d2}d")


def test_const_matrix_two_precision():
    """2x2 constant A vs [cos, -sin] closed form at dps 80 and 110."""
    dps_saved = mp.dps
    results = {}
    for dps in (80, 110):
        y = transport_fixed_eps(ConstRot(), "0", "0", "2", ["1", "0"], dps)
        results[dps] = y
        with mp.workdps(dps + 20):
            refs = [mp.cos(mpf(2)), -mp.sin(mpf(2))]
        for i in (0, 1):
            d = _agree_digits(y[i], refs[i], dps + 20)
            assert d >= dps - 6, f"rot dps={dps} comp {i}: only {d} digits"
    for i in (0, 1):
        d2 = _agree_digits(results[80][i], results[110][i], 130)
        assert d2 >= 72, f"rot two-precision comp {i}: only {d2} digits"
    assert mp.dps == dps_saved, "global mp.dps was mutated (import-dps footgun)"
    print("  2x2 constant-A rotation: closed form + two-precision OK")


def test_complex_detour_around_pole():
    """Pole at x=1 on the straight path 0.5 -> 2: a complex waypoint detour
    (MUM-crossing detour pattern) must give y0*(x1-1)/(x0-1); the straight
    path must raise (Zeno/step-collapse), not hang."""
    dps_saved = mp.dps
    dps = 80
    y = transport_fixed_eps(SimplePole(), "0", "0.5", "2", ["1"], dps,
                            path=["1+0.75j"])
    with mp.workdps(dps + 20):
        ref = (mpf(2) - 1) / (mpf("0.5") - 1)
    d = _agree_digits(y[0], ref, dps + 20)
    assert d >= dps - 6, f"detour: only {d} digits vs closed form"
    # detour BELOW the pole must give the same value (single-valued here)
    y2 = transport_fixed_eps(SimplePole(), "0", "0.5", "2", ["1"], dps,
                             path=["1-0.75j"])
    d12 = _agree_digits(y[0], y2[0], dps + 20)
    assert d12 >= dps - 8, f"detour side-independence: only {d12} digits"
    # straight path through the pole: must raise, quickly
    t0 = time.time()
    raised = False
    try:
        transport_fixed_eps(SimplePole(), "0", "0.5", "2", ["1"], 40)
    except (RuntimeError, AssertionError):
        raised = True
    assert raised, "straight path through pole did not raise"
    assert time.time() - t0 < 60, "Zeno guard too slow"
    assert mp.dps == dps_saved, "global mp.dps was mutated (import-dps footgun)"
    print(f"  complex detour: {d}d vs closed form, side-independent {d12}d, "
          f"straight path raised OK")


class ComplexPairPole:
    """y' = -2x/(x^2 + delta^2) y, delta = 1/8 (dyadic): the ONLY
    singularities are the complex pair +-i*delta just off the real path.
    Closed form from y(x0)=y0: y(x1) = y0 (x0^2+delta^2)/(x1^2+delta^2)."""
    n = 1
    var = "x"
    # 0.125 is dyadic-exact; declared as (re, im) pairs (transport._to_mpc)
    meta = {"singular_points": [(0, "0.125"), (0, "-0.125")]}

    def A(self, x, eps, dps):
        with mp.workdps(dps):
            xv = mpc(x)
            return [[-2 * xv / (xv * xv + mpf("0.125") ** 2)]]


def test_complex_pole_pair_off_axis():
    """G1 regression (transport_yline.py step-control parity).

    Two catalogued silent-wrong rules are pinned here at once:
      (a) REAL-pole-only clearance: both poles are complex (+-i/8), so a
          real-pole-only rule sees a singularity-free path 0 -> 1 and steps
          h_wrong = ratio*leg = 1/2 >> delta = radius of convergence. The
          degree-M Taylor polynomial of the TRUE solution about x=0 then
          misses the closed form by ~(h_wrong/delta)^M — demonstrated below
          against the closed form (provable failure, no engine needed).
      (b) single-term acceptance: about x=0 the solution is EVEN, so every
          odd Taylor term is ~0 (exactly zero up to circle roundoff) — the
          pre-port single-term early-break fired at order ~9 with a ~1e-3
          relative tail (dip-then-climb). The geometric two-term bound does
          not; the closed-form digit gate below kills any regression to
          either rule.
    """
    dps_saved = mp.dps
    results = {}
    d2 = None
    with mp.workdps(40):
        d2 = mpf("0.125") ** 2
    for dps in (60, 80):
        y, diag = transport_fixed_eps(ComplexPairPole(), "0", "0", "1",
                                      ["1"], dps, mtay=140, return_diag=True)
        results[dps] = y[0]
        with mp.workdps(dps + 20):
            ref = mpf("0.125") ** 2 / (1 + mpf("0.125") ** 2)
        d = _agree_digits(y[0], ref, dps + 20)
        assert d >= dps - 6, \
            f"complex-pair dps={dps}: only {d} digits vs closed form " \
            f"(real-pole-only clearance or single-term acceptance regression)"
        assert diag["steps"] >= 8, \
            f"only {diag['steps']} steps: complex poles did not throttle " \
            f"the march — real-pole-only clearance regression"
    d2p = _agree_digits(results[60], results[80], 110)
    assert d2p >= 54, f"two-precision agreement only {d2p} digits"
    # (a) made concrete: the wrong rule's first step from x=0 provably fails.
    # Taylor of y about 0 is sum (-1)^k (x/delta)^(2k); at h_wrong = 1/2 the
    # degree-140 partial sum is off by ~ (h/delta)^140 ~ 10^84.
    with mp.workdps(120):
        delta = mpf("0.125")
        h_wrong = mpf(1) / 2          # 0.5 x leg length: no real pole seen
        assert h_wrong > delta        # exceeds the true convergence radius
        psum = mpc(0)
        for k in range(0, 71):
            psum += (-1) ** k * (h_wrong / delta) ** (2 * k)
        truth = delta ** 2 / (h_wrong ** 2 + delta ** 2)
        rel = abs(psum - truth) / abs(truth)
        assert rel > mpf(10) ** 30, \
            f"wrong-rule partial sum unexpectedly close (rel={mp.nstr(rel, 3)})"
    assert mp.dps == dps_saved, "global mp.dps was mutated (import-dps footgun)"
    print(f"  complex pole pair +-i/8: closed form OK at 60/80 "
          f"(pair {d2p}d, steps={diag['steps']}); wrong-rule step provably "
          f"off by 1e{int(mp.log10(rel))}")


def test_singular_landing_refused():
    """Landing on a declared singular point must AssertionError up front
    (chassis rule: 'entry singular at y=0 — system not eligible')."""
    dps_saved = mp.dps
    raised = False
    try:
        transport_fixed_eps(ScalarPow(), "0.1", "1", "0", ["1"], 50)
    except AssertionError:
        raised = True
    assert raised, "singular landing point was not refused"
    assert mp.dps == dps_saved
    print("  singular landing refused OK")


def test_small_norm_y0_relative_guard():
    """F1 regression: the convergence guard must be RELATIVE to
    ||y0||. Before the fix the absolute early-break cut the Taylor series at
    order ~9 for small-norm states (measured on this exact probe at dps=65:
    y0=1e-30 -> 43.0 digits, 1e-60 -> 12.8, 1e-90 -> 5.6). Bar: every scale
    gives >= dps-2 RELATIVE digits vs the closed form y0*(x1/x0)^a."""
    dps_saved = mp.dps
    dps = 65
    for s in ("1", "1e-30", "1e-60", "1e-90"):
        y = transport_fixed_eps(ScalarPow(), "0.1", "1", "3", [s], dps)
        with mp.workdps(dps + 40):
            a = mpf(3) / 2 + mpf(1) / 10
            ref = mp.mpmathify(s) * mp.power(mpf(3), a)
            rel = abs(mpc(y[0]) - ref) / abs(ref)
            d = int(-mp.log10(rel)) if rel > 0 else dps + 40
        assert d >= dps - 2, \
            f"y0={s}: only {d} relative digits (F1 small-norm regression)"
    assert mp.dps == dps_saved, "global mp.dps was mutated"
    print("  small-norm y0 (1..1e-90): >= dps-2 relative digits at all scales")


class MonomialA:
    """y' = x^k y (entire A, no singular points): local Taylor support gap
    k+1 about x = 0 — the F1 sparse-lattice geometry (de_load monomial-cone
    shape). Closed form y(x) = y0 exp(x^(k+1)/(k+1))."""
    n = 1
    var = "x"
    meta = {}

    def __init__(self, k):
        self.k = k

    def A(self, x, eps, dps):
        with mp.workdps(dps):
            return [[mpc(x) ** self.k]]


def test_sparse_lattice_monomials():
    """Sparse-Taylor-lattice regression (
    FIXED same day): y' = x^k y launched from x = 0 puts the solution's
    Taylor support on the lattice (k+1)Z, so for k >= 2 the series has >= 2
    consecutive zero/roundoff terms; the pre-fix raw-last-two-terms early
    break saw _tail_bound(0, 0) = 0 and accepted a catastrophically
    truncated step while banking a false trunc_worst certificate (measured
    pre-fix on THIS probe, dps 60 leg 0->1: k=2 -> 6 correct digits, recorded
    trunc_worst = 9.9e-69 — false by ~63 orders; k=5 -> 5 digits).

    Bar: k in {2,3,5,10,12} at
    dps 60 AND 80, each >= dps-6 digits vs exp(1/(k+1)), certificate honest
    up to the output-rounding floor: rel_err <= 10*trunc_worst +
    10^-(dps-2). (The returned vector is ROUNDED to dps, so ~10^-dps
    relative error is a floor no truncation certificate below it can beat —
    the honesty clause catches any falsification beyond that floor; the
    pre-fix falsification was ~63 orders beyond it.) k = 1 is the
    dense-boundary control (period-2 solution lattice, worked pre-fix) and
    must stay at full precision."""
    dps_saved = mp.dps
    for k in (1, 2, 3, 5, 10, 12):
        for dps in (60, 80):
            y, diag = transport_fixed_eps(MonomialA(k), "0", "0", "1",
                                          ["1"], dps, return_diag=True)
            with mp.workdps(dps + 30):
                ref = mp.exp(mpf(1) / (k + 1))
                rel = abs(mpc(y[0]) - ref) / abs(ref)
                d = int(-mp.log10(rel)) if rel > 0 else dps + 30
                tw = mp.mpmathify(diag["trunc_worst"])
                honest = rel <= 10 * tw + mpf(10) ** (-(dps - 2))
            assert d >= dps - 6, \
                f"k={k} dps={dps}: only {d} digits vs exp(1/(k+1)) " \
                f"(F1 sparse-lattice regression)"
            assert honest, \
                f"k={k} dps={dps}: trunc_worst {diag['trunc_worst']} " \
                f"FALSIFIED (true rel err ~1e-{d} beyond 10x certificate + " \
                f"rounding floor)"
    assert mp.dps == dps_saved, "global mp.dps was mutated (import-dps footgun)"
    print("  sparse-lattice monomials y'=x^k y, k=1,2,3,5,10,12 @ dps 60/80: "
          ">= dps-6 digits, trunc_worst honest at all 12")


def test_diag_return_and_aseries_guard():
    """F7: return_diag=True returns (y, diag) with steps/trunc_worst/
    min_hfrac. F6: an A_series fast path returning fewer than M+1 orders
    must raise ValueError, never zero-pad."""
    dps_saved = mp.dps
    y, diag = transport_fixed_eps(ScalarPow(), "0.1", "1", "3", ["1"], 60,
                                  return_diag=True)
    assert isinstance(diag["steps"], int) and diag["steps"] > 0
    assert isinstance(diag["trunc_worst"], str)
    assert diag["trunc_worst_log10"] is None or \
        diag["trunc_worst_log10"] < -50, \
        f"trunc_worst_log10 {diag['trunc_worst_log10']} not below guard scale"
    assert 0 < diag["min_hfrac"] <= 0.5
    # F6: truncated fast path must raise, not silently zero-pad
    class ShortSeries:
        n = 1
        var = "x"
        meta = {}

        def A(self, x, eps, dps):
            return [[mpc(1)]]

        def A_series(self, x, eps, dps, M):
            return [[[mpc(1)]]]        # only order 0, whatever M was asked

    raised = False
    try:
        transport_fixed_eps(ShortSeries(), "0", "0", "1", ["1"], 50)
    except ValueError as e:
        raised = "zero-pad" in str(e)
    assert raised, "short A_series return was not refused"
    assert mp.dps == dps_saved, "global mp.dps was mutated"
    print(f"  diag return OK (steps={diag['steps']}, trunc_worst="
          f"{diag['trunc_worst']}, min_hfrac={diag['min_hfrac']}); "
          f"short A_series refused OK")


def main():
    t0 = time.time()
    print("test_transport.py self-checks:")
    test_scalar_power_law_two_precision()
    test_const_matrix_two_precision()
    test_complex_detour_around_pole()
    test_complex_pole_pair_off_axis()
    test_singular_landing_refused()
    test_small_norm_y0_relative_guard()
    test_sparse_lattice_monomials()
    test_diag_return_and_aseries_guard()
    print(f"ALL PASS ({time.time() - t0:.1f}s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
r"""
test_two_sector.py — synthetic unit test for two_sector.py (the fixed-eps
port of the reference two-sector Frobenius engine pieces; provenance in the
module header).

The reference is a GAUGE-CONSTRUCTED closed form (no pasted digits):
start from the decoupled base B = [[0, 0], [g0*s, -eps/s]] with solutions
    w1 = c1,   w2 = c2*(-s)^(-eps) + g0*c1*s^2/(2+eps),
and gauge y = T w with T = [[1, c], [0, s]] (det = s vanishing at 0 — that
is exactly how APPARENT poles arise). Then A = T B T^-1 + T' T^-1:

    A = [[ c*g0*s,            -c^2*g0 - c*eps/s^2 ],
         [ g0*s^2,            -c*g0*s + (1-eps)/s ]]

(A_01 has an apparent s^-2; the within-block coupling is one-directional:
p_01 = 2, A_10 analytic — the certified slot-schedule condition), and

    y1 = c1*(1 + c*g0*s^2/(2+eps)) + c*c2*(-s)^(-eps)
    y2 = g0*c1*s^3/(2+eps)         + c2*s*(-s)^(-eps)

i.e. EXACT terminating series
    a_0 = (c1, 0), a_2 = (c*g0*c1/(2+eps), 0), a_3 = (0, g0*c1/(2+eps))
    b_0 = (c*c2, 0), b_1 = (0, c2)                       (all others zero).

The B-branch level-1 slot is STRUCTURALLY RESONANT (the gauge pushes an
exponent to 1 - eps = (-eps) + 1), and b_{1,1} = c2 is determined by the
level-0 constraint row — the exact pattern solve_sector_autopin
handles on row 20. The closed form itself is verified against the DE
in-test (y' - A y == 0 at high precision) before being used as reference,
so a derivation slip here cannot silently pass.
"""

import os
import sys

from mpmath import mp, mpf, mpc

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from two_sector import two_sector_series, eval_two_sector  # noqa: E402

C_GAUGE = 2
G0 = 3


class GaugeApparent:
    """The gauge-constructed system above; only true singularity: s = 0."""
    n = 2
    var = "s"
    meta = {"singular_points": ["0"]}

    def A(self, s, eps, dps):
        with mp.workdps(dps):
            sv = mpc(s)
            e = mpc(eps)
            c = mpc(C_GAUGE)
            g = mpc(G0)
            return [[c * g * sv, -c * c * g - c * e / (sv * sv)],
                    [g * sv * sv, -c * g * sv + (1 - e) / sv]]


def _closed_form(s, eps, c1, c2, wp):
    """(y1, y2) and (y1', y2') of the closed form at s (inside workdps)."""
    with mp.workdps(wp):
        sv = mpc(s)
        e = mpc(eps)
        c = mpc(C_GAUGE)
        g = mpc(G0)
        pref = mp.exp(-e * mp.log(-sv))
        y1 = c1 * (1 + c * g * sv ** 2 / (2 + e)) + c * c2 * pref
        y2 = g * c1 * sv ** 3 / (2 + e) + c2 * sv * pref
        d1 = c1 * 2 * c * g * sv / (2 + e) + c * c2 * pref * (-e / sv)
        d2 = 3 * g * c1 * sv ** 2 / (2 + e) + c2 * pref * (1 - e)
        return (y1, y2), (d1, d2)


def _check_reference_satisfies_de():
    """Honesty gate on the reference itself: y' - A y == 0 to ~wp digits."""
    wp = 120
    desys = GaugeApparent()
    eps = mpf(1) / 3
    with mp.workdps(wp):
        for s in (mpf(-1) / 32, mpf(-3) / 64):
            (y1, y2), (d1, d2) = _closed_form(s, eps, mpc(1), mpc(1), wp)
            A = desys.A(s, eps, wp)
            r1 = d1 - (A[0][0] * y1 + A[0][1] * y2)
            r2 = d2 - (A[1][0] * y1 + A[1][1] * y2)
            scale = max(abs(d1), abs(d2), mpf(1))
            assert abs(r1) < scale * mpf(10) ** (-(wp - 15)), abs(r1)
            assert abs(r2) < scale * mpf(10) ** (-(wp - 15)), abs(r2)


def test_two_sector_autopin_closed_form():
    """Full engine vs the gauge closed form, two precisions."""
    dps_saved = mp.dps
    _check_reference_satisfies_de()
    eps = mpf(1) / 3
    NL = 8
    results = {}
    for dps in (50, 70):
        res = two_sector_series(GaugeApparent(), eps, "0",
                                seed_a=["1", "0"], seed_b=[str(C_GAUGE), "0"],
                                NL=NL, dps=dps, pmax=2, autopin=True)
        results[dps] = res
        # structure: apparent pole + delays detected
        assert res["poles"]["0,1"] == 2, res["poles"]
        assert res["blocks"] == [[0, 1]], res["blocks"]
        # autopin found the structurally resonant slot (component 1, level 1)
        assert list(res["pins_B"]) == [(1, 1)], res["pins_B"]
        assert res["pins_A"] == {}, res["pins_A"]
        with mp.workdps(dps + 10):
            tol = mpf(10) ** (-(dps - 8))
            # pin value t = c2 = 1
            assert abs(mpc(res["pins_B"][(1, 1)]) - 1) < tol, res["pins_B"]
            e = mpc(eps)
            want_A = [[mpc(0)] * (NL + 1) for _ in range(2)]
            want_B = [[mpc(0)] * (NL + 1) for _ in range(2)]
            want_A[0][0] = mpc(1)
            want_A[0][2] = C_GAUGE * G0 / (2 + e)
            want_A[1][3] = G0 / (2 + e)
            want_B[0][0] = mpc(C_GAUGE)
            want_B[1][1] = mpc(1)
            for i in range(2):
                for m in range(NL + 1):
                    ga = mpc(res["A"][i][m])
                    gb = mpc(res["B"][i][m])
                    assert abs(ga - want_A[i][m]) < tol, \
                        ("A", i, m, ga, want_A[i][m])
                    assert abs(gb - want_B[i][m]) < tol, \
                        ("B", i, m, gb, want_B[i][m])
            # constraint residuals: reported and tiny
            for k, v in list(res["residuals_A"].items()) + \
                    list(res["residuals_B"].items()):
                assert abs(mpc(v)) < mpf(10) ** (-(dps - 2)), (k, v)
        # evaluation against the closed form off the expansion point
        s_eval = mpf(-1) / 64
        vals = eval_two_sector(res, s_eval, dps)
        with mp.workdps(dps + 15):
            (y1, y2), _ = _closed_form(s_eval, eps, mpc(1), mpc(1), dps + 15)
            for got, want in zip(vals, (y1, y2)):
                assert abs(mpc(got) - want) < abs(want) * \
                    mpf(10) ** (-(dps - 8)), (got, want)
    # two-precision agreement of every coefficient
    with mp.workdps(90):
        for key in ("A", "B"):
            for i in range(2):
                for m in range(NL + 1):
                    d = abs(mpc(results[50][key][i][m]) -
                            mpc(results[70][key][i][m]))
                    assert d < mpf(10) ** -42, (key, i, m, d)
    assert mp.dps == dps_saved, "global mp.dps mutated"
    print("  two-sector autopin: closed form matched at dps 50/70, "
          "pin t=c2, residuals reported tiny, eval off-point OK")


def test_two_sector_resonance_raises_without_pin():
    """autopin=False on the resonant branch must raise loudly (reference
    'RESONANT slot needs pin' contract), naming the coordinate."""
    dps_saved = mp.dps
    eps = mpf(1) / 3
    try:
        two_sector_series(GaugeApparent(), eps, "0",
                          seed_a=["1", "0"], seed_b=["2", "0"],
                          NL=6, dps=40, pmax=2, autopin=False)
        raise RuntimeError("resonant slot did not raise")
    except AssertionError as e:
        assert "RESONANT" in str(e) and "(1, 1)" in str(e), str(e)
    assert mp.dps == dps_saved
    print("  two-sector: unpinned resonant slot raises with coordinate")


def test_two_sector_pole_order_tripwire():
    """pmax too small must trip the s^-k tripwire, not alias silently."""
    dps_saved = mp.dps
    try:
        two_sector_series(GaugeApparent(), mpf(1) / 3, "0",
                          seed_a=["1", "0"], seed_b=["2", "0"],
                          NL=6, dps=40, pmax=1, autopin=True)
        raise RuntimeError("pmax=1 did not trip on the s^-2 entry")
    except ValueError as e:
        assert "pole order" in str(e), str(e)
    assert mp.dps == dps_saved
    print("  two-sector: pole-order > pmax tripwire fires")


def main():
    import time
    t0 = time.time()
    print("test_two_sector.py self-checks:")
    test_two_sector_autopin_closed_form()
    test_two_sector_resonance_raises_without_pin()
    test_two_sector_pole_order_tripwire()
    print(f"ALL PASS ({time.time() - t0:.1f}s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())

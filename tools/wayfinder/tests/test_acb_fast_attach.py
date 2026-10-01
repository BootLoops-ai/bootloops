#!/usr/bin/env python3
"""test_acb_fast_attach.py — battery for the ACB FAST-ATTACH member
(acb_fast_attach.py: attach_fast_path + attach_acb_fast for exact-rational
scalar companion systems).

SCALED-DOWN EXACT CONTROL: the order-2 exact-Q operator

    L = (1-x) D^2 - D            (solutions {1, -log(1-x)})

marched as a scalar companion system through attach_fast_path /
attach_acb_fast + transport_fixed_eps, x=1/8 -> 3/4, eps=0, on BOTH
backends.  Gates (all script-emitted):
  G1 exact control  : each backend's landing agrees with the closed form
                      [-log(1/4), 4] to >= dps - 5 digits (dps=120).
  G2 backend parity : mpmath vs acb landings agree to >= dps - 5 digits.
  G3 identical chassis certificates: same accepted step count; the
      per-step worst geometric-tail certificates (trunc_worst) agree.
  G4 speed ratio    : wall ratio mpmath/acb RECORDED, not asserted (toy
      n=2 scale; the backend rule rests on a measured order-17 march
      at dps 430, mpmath 73 s -> acb 2.2 s, 33x — see GUIDE.md).
  G5 fences         : acb backend refuses eps != 0 and dps > table wp.

The mpmath G1 leg needs no python-flint; the acb legs (G1-acb, G2-G5)
SKIP without it.  The receipt JSON goes to a temp dir (pytest tmp_path),
never beside the test.  Runnable directly
(python3 tests/test_acb_fast_attach.py) or under pytest.
"""
import json
import os
import sys
import tempfile
import time
from fractions import Fraction as F

import mpmath as mp
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import acb_fast_attach as AFA                       # noqa: E402
from transport import transport_fixed_eps, HAVE_FLINT  # noqa: E402

DPS = 120
WP_TABLE = DPS + 40
X0, X1 = F(1, 8), F(3, 4)


class ToyCompanionSystem:
    """Minimal DESystem duck-type for the battery operator (contract per
    the acb_fast_attach module docstring: .n / .coeffs / .A / .meta)."""

    def __init__(self):
        # L = p2 D^2 + p1 D + p0, p2 = 1 - x, p1 = -1, p0 = 0
        self.coeffs = [[F(0)], [F(-1)], [F(1), F(-1)]]
        self.n = 2
        self.var = "x"
        self.meta = {"name": "battery_log_toy",
                     "singular_points": [1.0]}

    def A(self, x, eps, dps):
        if mp.mpc(eps) != 0:
            raise ValueError("pure-Q operator: eps must be 0")
        with mp.workdps(dps):
            z = mp.mpc(x)
            p2 = 1 - z
            if p2 == 0:
                raise ZeroDivisionError("lc vanishes")
            return [[mp.mpc(0), mp.mpc(1)],
                    [mp.mpc(0), 1 / p2]]   # y2' = -p1/p2 y2 = y2/(1-x)


def digits_agree(a, b):
    """Common digits of two mpc vectors (relative)."""
    d = max(abs(a[i] - b[i]) for i in range(len(a)))
    s = max(abs(v) for v in b)
    if d == 0:
        return 999.0
    return float(-mp.log10(d / s))


def _seed_and_truth():
    with mp.workdps(DPS + 60):
        y0 = [mp.log(mp.mpf(8) / 7), mp.mpf(8) / 7]
        truth = [mp.log(4), mp.mpf(4)]
    return y0, truth


def _run_backend(backend):
    sysm = ToyCompanionSystem()
    AFA.attach_fast_path(sysm)
    if backend == "acb":
        AFA.attach_acb_fast(sysm, WP_TABLE)
    y0, truth = _seed_and_truth()
    t0 = time.perf_counter()
    y, diag = transport_fixed_eps(
        sysm, 0, X0, X1, list(y0), DPS, return_diag=True,
        backend=(None if backend == "mpmath" else "acb"))
    wall = time.perf_counter() - t0
    with mp.workdps(DPS + 60):
        dtruth = digits_agree(y, truth)
    return {"y": [mp.nstr(v, DPS) for v in y], "diag": diag,
            "wall_s": round(wall, 4),
            "digits_vs_exact_control": round(dtruth, 2)}


def test_g1_mpmath_exact_control():
    """G1 (mpmath half): pure-mpmath fast path lands on the closed form.
    Needs no python-flint."""
    res = _run_backend("mpmath")
    assert res["digits_vs_exact_control"] >= DPS - 5, res


def _acb_gates(out_dir):
    """G1(acb)+G2-G5, receipt JSON script-emitted into out_dir."""
    receipt = {"PRODUCER": {
        "leg": "wayfinder acb fast-attach member battery",
        "script": "test_acb_fast_attach.py",
        "stamp_utc": time.strftime("%a %b %d %H:%M:%S UTC %Y", time.gmtime()),
        "control": ("SCALED-DOWN EXACT CONTROL: order-2 (1-x)D^2 - D, "
                    "solutions {1, -log(1-x)}")},
        "params": {"dps": DPS, "acb_table_wp": WP_TABLE,
                   "x0": str(X0), "x1": str(X1), "eps": 0}}

    results = {b: _run_backend(b) for b in ("mpmath", "acb")}

    with mp.workdps(DPS + 60):
        ya = [mp.mpmathify(v) for v in results["mpmath"]["y"]]
        yb = [mp.mpmathify(v) for v in results["acb"]["y"]]
        parity = digits_agree(ya, yb)

    da, db = results["mpmath"]["diag"], results["acb"]["diag"]
    g1 = (results["mpmath"]["digits_vs_exact_control"] >= DPS - 5 and
          results["acb"]["digits_vs_exact_control"] >= DPS - 5)
    g2 = parity >= DPS - 5
    g3_steps = (da["steps"] == db["steps"])
    tw_a, tw_b = da["trunc_worst_log10"], db["trunc_worst_log10"]
    g3_cert = (tw_a is None and tw_b is None) or \
        (tw_a is not None and tw_b is not None and abs(tw_a - tw_b) < 1.0)
    ratio = results["mpmath"]["wall_s"] / max(results["acb"]["wall_s"], 1e-9)

    # G5 fences on the acb backend
    fence_eps = fence_wp = False
    y0, _ = _seed_and_truth()
    sysm = ToyCompanionSystem()
    AFA.attach_fast_path(sysm)
    AFA.attach_acb_fast(sysm, WP_TABLE)
    try:
        transport_fixed_eps(sysm, mp.mpf("1e-30"), X0, X1, list(y0), DPS,
                            backend="acb")
    except ValueError:
        fence_eps = True
    try:
        transport_fixed_eps(sysm, 0, X0, X1, list(y0), WP_TABLE + 50,
                            backend="acb")
    except ValueError:
        fence_wp = True

    receipt["results"] = results
    receipt["gates"] = {
        "G1_exact_control_pass": bool(g1),
        "G2_backend_parity_digits": round(parity, 2),
        "G2_backend_parity_pass": bool(g2),
        "G3_steps_mpmath": da["steps"], "G3_steps_acb": db["steps"],
        "G3_identical_step_count": bool(g3_steps),
        "G3_trunc_worst_log10_mpmath": tw_a,
        "G3_trunc_worst_log10_acb": tw_b,
        "G3_certificate_agree": bool(g3_cert),
        "G4_wall_mpmath_s": results["mpmath"]["wall_s"],
        "G4_wall_acb_s": results["acb"]["wall_s"],
        "G4_speed_ratio_mpmath_over_acb": round(ratio, 2),
        "G4_note": "toy n=2 scale; recorded, not asserted",
        "G5_refuses_nonzero_eps": bool(fence_eps),
        "G5_refuses_dps_gt_table_wp": bool(fence_wp),
    }
    verdict = all([g1, g2, g3_steps, g3_cert, fence_eps, fence_wp])
    receipt["verdict"] = "PASS" if verdict else "FAIL"

    out = os.path.join(out_dir, "RECEIPT_acb_fast_attach.json")
    with open(out, "w") as f:
        json.dump(receipt, f, indent=1)
    assert g1, receipt["gates"]
    assert g2, receipt["gates"]
    assert g3_steps, receipt["gates"]
    assert g3_cert, receipt["gates"]
    assert fence_eps and fence_wp, receipt["gates"]
    return receipt


def test_acb_backend_gates(tmp_path):
    """G1(acb)+G2+G3+G4(recorded)+G5 — SKIP without python-flint."""
    if not HAVE_FLINT:
        pytest.skip("python-flint not installed: acb legs skipped")
    _acb_gates(str(tmp_path))


if __name__ == "__main__":
    test_g1_mpmath_exact_control()
    print("G1 (mpmath exact control): PASS")
    if HAVE_FLINT:
        r = _acb_gates(tempfile.mkdtemp(prefix="acb_fast_attach_"))
        print(json.dumps(r["gates"], indent=1))
        print("VERDICT:", r["verdict"])
        sys.exit(0 if r["verdict"] == "PASS" else 1)
    print("acb legs: SKIP (python-flint not installed)")
    sys.exit(0)

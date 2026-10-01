#!/usr/bin/env python3
"""Battery for common/verdict.py: gate-primitive units + sha-stamped
atomic receipt round-trip (no external inputs needed — pure chassis)."""
import json, os, sys, tempfile
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import verdict as V

m, r = V.parse_ball("[2.0371060933e-8 +/- 4.36e-148]")
assert float(m) == 2.0371060933e-8 and float(r) == 4.36e-148
assert V.matched_digits(m, m) == 999
V.ball_gate("unit like-for-like", "[1.0 +/- 1e-10]", "[1.0 +/- 1e-10]", 9)
V.route_gate("unit cross-route", "[1.0 +/- 1e-3]", "[1.0001 +/- 1e-12]", 3)
assert V.two_dps_gate("unit two-dps", "[3.14159265 +/- 1e-3]",
                      "[3.141592653589 +/- 1e-40]", 8) >= 8
ok = False
try:
    V.ball_gate("unit radius-bar must-fail", "[1.0 +/- 1.0]",
                "[1.0 +/- 1e-20]", 1)
except AssertionError:
    ok = True
assert ok, "ball_gate waived the radius-comparability bar"
print("[unit radius-bar] PASS  must-fail mutation control held")

with tempfile.TemporaryDirectory() as td:
    out = os.path.join(td, "r.json")
    rec = V.receipt({"x": 1}, out, code_paths=(os.path.join(HERE, "verdict.py"),))
    back = json.load(open(out))
    assert back["schema"] == "terrier-receipt-v1"
    assert back["code_sha256"] == V.code_sha(os.path.join(HERE, "verdict.py"))
    assert back["payload"] == {"x": 1} and not os.path.exists(out + ".tmp")
print("[receipt] PASS  atomic sha-stamped round-trip")
print("SELFTEST common/verdict.py: ALL PASS")

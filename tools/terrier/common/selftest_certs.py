#!/usr/bin/env python3
"""battery: certs.py — parse/matched-digits, gates, honesty, float64 guard."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import certs as X

n = [0]
def ok(name, cond):
    n[0] += 1
    print("[%s] %s" % ("PASS" if cond else "FAIL", name))
    assert cond, name

W = "[2.0371060933111834182083117862342977e-8 +/- 3e-150]"
WB = "[2.0371060933111834182083117999999999e-8 +/- 3e-150]"   # differs @ ~26 dg
m, r = X.parse_ball(W)
ok("parse ball", float(m) == 2.0371060933111834e-08 and r > 0)
ok("parse plain", X.parse_ball("1.5")[1] == 0)
ok("matched_digits identical", X.matched_digits(m, m) == 999)
d = X.matched_digits(m, X.parse_ball(WB)[0])
ok("matched_digits ~26", 24 <= d <= 28)
ok("is_ball_str law", X.is_ball_str(W) and not X.is_ball_str("2.037e-8"))
# two-dps (c1): pass at bar 20, fail at bar 40
ok("two-dps pass", X.two_dps_gate(WB, W, 20)["pass"])
ok("two-dps fail", not X.two_dps_gate(WB, W, 40)["pass"])
# ball gate: same-route regression — radius comparability enforced
ok("ball gate pass", X.ball_gate(W, W, 140)["pass"])
wide = "[2.0371060933111834182083117862342977e-8 +/- 3e-100]"
g = X.ball_gate(wide, W, 20)
ok("ball gate radius-ratio fail", not g["pass"] and g["radius_ratio"] > 16)
far = "[2.0372e-8 +/- 1e-20]"
ok("ball gate overlap fail", not X.ball_gate(far, W, 2)["pass"])
# route gate: no radius bar (C/R ratio ~2e9 must pass), digits still gate
ok("route gate wide-radius pass", X.route_gate(wide, W, 25)["pass"])
ok("route gate digits fail", not X.route_gate(WB, W, 40)["pass"])
# ball honesty: radius quoted + rel bar in decimal-string (never float64)
h = X.ball_honesty(W, "1e-140")
ok("honesty pass + quotes radius", h["pass"] and h["quoted_ball"] == W)
ok("honesty fail on wide ball", not X.ball_honesty(wide, "1e-140")["pass"])
try:
    X.ball_honesty("2.037e-8", "1e-3"); bad = False
except AssertionError:
    bad = True
ok("bare decimal refused as certified", bad)
# float64 trim guard: subnormal mid / underflowing radius must be flagged
ok("float64 safe (normal)", X.float64_safe(W))
ok("float64 unsafe: rad underflow", not X.float64_safe(
    "[2.03e-8 +/- 1e-330]"))
ok("float64 unsafe: subnormal mid", not X.float64_safe(
    "[1.5e-322 +/- 1e-340]"))
# published-value gate (H4): DKMM 2.037e-8 in-band; Broeckel 2.048e-8 out
ok("ulp 2.037e-8", abs(float(X.last_digit_ulp("2.037e-8")) - 1e-11) < 1e-26)
ok("published pass", X.published_value_gate(W, "2.037e-8")["pass"])
ok("published fail", not X.published_value_gate(W, "2.048e-8")["pass"])
ok("published int-mantissa", X.published_value_gate("[124 +/- 0]", "124")["pass"])
print("selftest_certs: %d/%d green" % (n[0], n[0]))

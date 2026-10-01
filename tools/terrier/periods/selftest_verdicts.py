#!/usr/bin/env python3
"""Battery (verdicts.py + common/verdict.py): primitive unit tests + the
full 13/13 DKMM scoreboard (run_pipe gates, regenerated in place)."""
import json, os, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__))
PIPE = os.path.join(HERE, "pipeline")
sys.path.insert(0, HERE)
from verdicts import parse_ball, matched_digits, route_gate, two_dps_gate

# unit tests of the gate primitives
m, r = parse_ball("[2.037e-8 +/- 4.4e-148]")
assert float(m) == 2.037e-8 and float(r) == 4.4e-148
m2, r2 = parse_ball("1.5")
assert float(m2) == 1.5 and r2 == 0
assert matched_digits(m, m) == 999
import mpmath as mp
assert matched_digits(mp.mpf("1.00000001"), mp.mpf("1.0")) == 8
route_gate("unit route_gate", "[1.0 +/- 1e-3]", "[1.0005 +/- 1e-9]", 3)
two_dps_gate("unit two_dps", "[1.23456 +/- 1e-2]", "[1.234561 +/- 1e-30]", 5)
try:    # radius-comparability must NOT be waived in ball_gate
    from verdicts import ball_gate
    ball_gate("unit ball_gate must-fail", "[1.0 +/- 1.0]", "[1.0 +/- 1e-20]", 1)
    raise SystemExit("ball_gate accepted a 1e20 radius ratio — FAIL")
except AssertionError:
    print("[unit ball_gate] PASS  radius bar enforced (mutation control)")

# full scoreboard, regenerated in place
rp = subprocess.run([sys.executable, "run_pipe.py", "gates"], cwd=PIPE,
                    capture_output=True, text=True, timeout=600)
sys.stdout.write(rp.stdout[-1400:])
assert rp.returncode == 0, f"gates FAILED:\n{rp.stderr[-1500:]}"
g = json.load(open(os.path.join(PIPE, "dkmm", "gates_pipe.json")))
assert g["G5_W0_oncurve"]["matched_digits"] == 999
assert g["G6_W0_vac"]["matched_digits"] == 999
assert "CONFIRMED" in g["verdict"]
print("SELFTEST verdicts.py: primitives + 13/13 scoreboard ALL PASS")

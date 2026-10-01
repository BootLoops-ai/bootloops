#!/usr/bin/env python3
"""Battery CEN-S (census/control_sampler): run the control_sampler battery
(exact-rational AD rejection sampling, dyadic coins, exact CP brackets,
Hoeffding p_U only-widens law). Gate: rc 0 + 'SUMMARY: 19 pass / 0 fail'."""
import os, re, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "control_sampler", "tests")
r = subprocess.run([sys.executable, "run_all.py"], cwd=D,
                   capture_output=True, text=True, timeout=1500)
out = r.stdout + r.stderr
print(out[-2000:])
assert r.returncode == 0, "control_sampler battery rc != 0"
m = re.search(r"SUMMARY: (\d+) pass / (\d+) fail", out)
assert m and m.group(2) == "0", "SUMMARY line missing or failures"
assert int(m.group(1)) >= 19, "battery shrank below the expected 19"
print(f"SELFTEST controlsampler: ALL PASS ({m.group(0)})")

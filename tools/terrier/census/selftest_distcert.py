#!/usr/bin/env python3
"""Battery CEN-C (census/dist_cert): run the dist_cert battery (exact
Fraction balls, certified sqrt/ln, z/WP charts, root/CM predicates,
Fincke-Pohst OFF certificates). Gate: rc 0 + '37/37 checks passed'."""
import os, re, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "dist_cert", "tests")
r = subprocess.run([sys.executable, "test_dist_cert.py"], cwd=D,
                   capture_output=True, text=True, timeout=1500)
out = r.stdout + r.stderr
print(out[-2000:])
assert r.returncode == 0, "dist_cert battery rc != 0"
m = re.search(r"(\d+)/(\d+) checks passed", out)
assert m and m.group(1) == m.group(2), "pass-count line missing/short"
assert int(m.group(1)) >= 37, "battery shrank below the expected 37"
print(f"SELFTEST distcert: ALL PASS ({m.group(0)})")

#!/usr/bin/env python3
"""Battery CEN-D (census/dedup): replay the two-stack dedup battery,
T1-T8 synthetic (fixed seeds; proof-level ground truth by
construction). Gates: rc 0 + 'BATTERY GREEN (8/8)'. Writes
census/dedup/test_receipts.json (volatile, regenerated). Distinct from
common/selftest_dedup.py (known-orbit dedup)."""
import os, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "dedup")
r = subprocess.run([sys.executable, "test_dedup.py"], cwd=D,
                   capture_output=True, text=True, timeout=1500)
out = r.stdout + r.stderr
print(out[-2000:])
assert r.returncode == 0, "dedup census battery rc != 0"
assert "BATTERY GREEN (8/8)" in out, "8/8 GREEN line missing"
assert "T7-red-reconciled-miss" in out and "T8-red-witness-fail" in out, \
    "drill rows missing"
print("SELFTEST dedup_census: ALL PASS (8/8 synthetic)")

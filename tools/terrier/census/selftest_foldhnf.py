#!/usr/bin/env python3
"""Battery CEN-F (census/fold_hnf): run the T1-T4 fold_hnf battery.
Gates: rc 0 + 'BATTERY GREEN: T1 T2 T3 T4' + T2 exact acceptance line present.
Needs the G4 anchor-set reference directory (not included in the package;
env TERRIER_G4_DIR)."""
import os, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "fold_hnf", "tests")
r = subprocess.run([sys.executable, "test_fold_hnf.py"], cwd=D,
                   capture_output=True, text=True, timeout=1500)
out = r.stdout + r.stderr
print(out[-2000:])
assert r.returncode == 0, "fold_hnf battery rc != 0"
assert "BATTERY GREEN: T1 T2 T3 T4" in out, "T1-T4 GREEN line missing"
assert "T2 PASS" in out and "T1 PASS" in out, "per-T PASS lines missing"
print("SELFTEST foldhnf: ALL PASS (T1-T4)")

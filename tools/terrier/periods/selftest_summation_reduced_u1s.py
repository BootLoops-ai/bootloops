#!/usr/bin/env python3
"""Battery (periods/summation): the reduced-size checks.
(1) b7_checks.py — B7 Cauchy floors: all 4 PASS lines EXACT (CS+jet
    majorant, Eulerian closed forms, tail corollary, log modulus).
(2) import-green — every summation module imports clean (boxops_u1 imports
    its generic engine from the box-operator bank, TERRIER_BOXOPS_BANK,
    which is not included in the package).
The full-size gates (B1 84-op termwise, B3 GF==lattice 210/210, ratio-law
200) and their receipts are not included in the package."""
import os, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "summation")

# (1) B7 Cauchy floors: 4/4 exact PASS lines
r = subprocess.run([sys.executable, "b7_checks.py"], cwd=D,
                   capture_output=True, text=True, timeout=300)
print(r.stdout, end="")
assert r.returncode == 0, "b7_checks rc != 0:\n" + r.stderr[-800:]
LINES = ("CS+jet majorant 20 trials: PASS",
         "Eulerian closed forms d=0..5: PASS",
         "tail corollary bound: PASS",
         "log modulus inequality: PASS")
for ln in LINES:
    assert ln in r.stdout, f"b7 floor line missing/failed: {ln!r}"
assert "FAIL" not in r.stdout, "b7_checks printed a FAIL"

# (2) import-green: all summation modules load in-tree
MODS = "towers6,boxops_u1,b3_entropy,b3_gates,b3_oracle,b3_pass,b7_numbers"
r2 = subprocess.run([sys.executable, "-c",
                     "import " + MODS.replace(",", ", ")],
                    cwd=D, capture_output=True, text=True, timeout=300)
assert r2.returncode == 0, "import-green FAIL:\n" + r2.stderr[-1200:]
print(f"SELFTEST summation_reduced_u1s: ALL PASS "
      f"(b7 Cauchy floors 4/4 + {len(MODS.split(','))} modules import-green)")

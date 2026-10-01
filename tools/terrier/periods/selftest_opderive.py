#!/usr/bin/env python3
"""Battery (opderive.py): run_pipe g3 — 2-prime x 2-truncation pin,
CRT+Wang exact, full-series annihilation proof, theta-form + recurrence ==
the stored pipeline/dkmm/operator_LS.json exactly (engine: tools/annihilator/;
the input period series is reference data not included in the package --
set TERRIER_KKLT_BANK)."""
import os, subprocess, sys
PIPE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "pipeline")
r = subprocess.run([sys.executable, "run_pipe.py", "g3"], cwd=PIPE,
                   capture_output=True, text=True, timeout=1500)
sys.stdout.write(r.stdout[-1200:])
assert r.returncode == 0, f"g3 battery FAILED:\n{r.stderr[-1500:]}"
assert "FAIL" not in r.stdout
print("SELFTEST opderive.py: g3 operator regression ALL PASS")

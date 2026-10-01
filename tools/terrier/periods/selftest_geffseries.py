#!/usr/bin/env python3
"""Battery geffseries.py: run_pipe g1 — DKMM Gamma-series/GV/w0/tower/
racetrack/mirror BYTE-EXACT vs the reference series (not included in the
package; set TERRIER_KKLT_BANK). g1 also rewrites series_curve/tower_curve/
gv_extracted/mirror_curve in pipeline/dkmm/."""
import os, subprocess, sys
PIPE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "pipeline")
r = subprocess.run([sys.executable, "run_pipe.py", "g1"], cwd=PIPE,
                   capture_output=True, text=True, timeout=1500)
sys.stdout.write(r.stdout[-1500:])
assert r.returncode == 0, f"g1 battery FAILED:\n{r.stderr[-1500:]}"
assert "FAIL" not in r.stdout
print("SELFTEST geffseries.py: g1 byte-exact regression ALL PASS")

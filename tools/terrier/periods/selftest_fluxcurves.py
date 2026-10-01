#!/usr/bin/env python3
"""Battery (fluxcurves.py): G2 exact flux/curve + stage1 frame gates —
G4 frame dictionary == reference certified-W0 frame exactly, f3 xi-mutation
MUST-FAIL control, extended towers exact.  Runs run_pipe {g2, stage1}
(stage1 consumes the g1/g3 outputs in pipeline/dkmm/; the reference data is
not included in the package -- set TERRIER_KKLT_BANK)."""
import os, subprocess, sys
PIPE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "pipeline")
for st in ("g2", "stage1"):
    r = subprocess.run([sys.executable, "run_pipe.py", st], cwd=PIPE,
                       capture_output=True, text=True, timeout=1500)
    sys.stdout.write(r.stdout[-1200:])
    assert r.returncode == 0, f"{st} FAILED:\n{r.stderr[-1500:]}"
    assert "FAIL" not in r.stdout
assert os.path.exists(os.path.join(PIPE, "dkmm", "frame_solution.json"))
print("SELFTEST fluxcurves.py: G2 + G4 frame/f3-mutation/towers ALL PASS")

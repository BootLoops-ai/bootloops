#!/usr/bin/env python3
"""Battery (transport.py): reproduce the reference DKMM certified balls
EXACTLY through this chassis.  Re-runs route R 150 + vac 150 R (consuming
the stage1/jets outputs in pipeline/dkmm/), then gates midpoint string-EXACT
(999) vs the reference certified-W0 receipts (not included in the package;
set TERRIER_KKLT_BANK)."""
import json, os, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__))
PIPE = os.path.join(HERE, "pipeline")
sys.path.insert(0, os.path.join(HERE, "..", "common"))
from verdict import parse_ball, matched_digits, verdict

for args in (["route", "R", "150"], ["vac", "150", "R"]):
    r = subprocess.run([sys.executable, "run_pipe.py"] + args, cwd=PIPE,
                       capture_output=True, text=True, timeout=1700)
    sys.stdout.write(r.stdout[-900:])
    assert r.returncode == 0, f"{args} FAILED:\n{r.stderr[-1500:]}"

KKLT = os.environ.get("TERRIER_KKLT_BANK")
if not KKLT:
    sys.exit("REFUSE: reference data not included in the package (certified "
             "W0 / W0-at-vacuum receipts of the two-modulus KKLT example) "
             "-- set TERRIER_KKLT_BANK to run")
mine = json.load(open(os.path.join(PIPE, "dkmm", "result_R_150.json")))
ref = json.load(open(os.path.join(KKLT, "cert_w0", "result_R_150.json")))
m1, r1 = parse_ball(mine["W0"]); m2, r2 = parse_ball(ref["W0"])
verdict("on-curve W0 midpoint EXACT", matched_digits(m1, m2) == 999,
        "on-curve ball midpoint string-identical to the reference receipt")
verdict("on-curve radius ratio", r2 > 0 and 0.9 < float(r1 / r2) < 1.12,
        f"ratio {float(r1/r2):.3f} (reference radius class ~4.6e-166)")
vm = json.load(open(os.path.join(PIPE, "dkmm", "result_vac_R_150.json")))
vp = json.load(open(os.path.join(KKLT, "cert_w0_vac", "result_vac_R_150.json")))
v1, s1 = parse_ball(vm["W0_vac"]); v2, s2 = parse_ball(vp["W0_vac"])
verdict("vacuum |W0|_vac midpoint EXACT", matched_digits(v1, v2) == 999,
        "Krawczyk-certified vacuum ball midpoint string-identical")
verdict("vacuum radius ratio", s2 > 0 and 0.9 < float(s1 / s2) < 1.12,
        f"ratio {float(s1/s2):.3f}")
print("SELFTEST transport.py: reference DKMM certified balls reproduced "
      "EXACTLY — ALL PASS")

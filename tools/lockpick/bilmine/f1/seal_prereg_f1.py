#!/usr/bin/env python3
# lockpick bilmine member (F1 ring stage) — seal ceremony emitter: pins the F1 prereg + every input of record (sha256) BEFORE any relation number exists.
"""seal_prereg_f1.py — BILMINE-F1 seal writer: pins the F1 leg's prereg +
every input of record (sha256) into work/SEAL_RECORD.json BEFORE any control
or mine run. Producer law: emits its own producer block; optional
producer-lint check.

Leg-dir fixture contract (every pinned file must exist before sealing; the
writer refuses loudly and lists what is missing):
  $BILMINE_F1_LEG/PREREG_BILMINE_F1.md, $BILMINE_F1_LEG/LAUNCH_LINE.txt
  $BILMINE_PARENT — the CLOSED parent leg: PREREG_BILMINE.md +
      work/{SEAL_RECORD,CONTROL_VERDICT,EXACT_INVARIANTS,MINE_P1,MINE_M1,
      MINE_W5,BILMINE_VERDICT}.json — plus the shipped parent library
      ../bilmine_lib.py (pinned by content)
  $BILMINE_T5/conn_{p1,m1}_{lo,hi}.json, linkb_{lo,hi}.json,
      CURE1_AGREEMENT_REEMIT.json
  $BILMINE_PFRAME, $BILMINE_C1_PI, $BILMINE_C1_N
"""
import hashlib, json, os, subprocess, sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from bilmine_lib import require_env, t5_dir  # noqa: E402
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bilmine_f1_lib import leg_dir_f1, parent_dir  # noqa: E402


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def main():
    LEG = leg_dir_f1()
    PARENT = parent_dir()
    T5 = t5_dir()
    parent_lib = os.path.join(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))), "bilmine_lib.py")
    inputs = {
        "conn_p1_lo": f"{T5}/conn_p1_lo.json",
        "conn_p1_hi": f"{T5}/conn_p1_hi.json",
        "conn_m1_lo": f"{T5}/conn_m1_lo.json",
        "conn_m1_hi": f"{T5}/conn_m1_hi.json",
        "linkb_lo": f"{T5}/linkb_lo.json",
        "linkb_hi": f"{T5}/linkb_hi.json",
        "cure1_caps": f"{T5}/CURE1_AGREEMENT_REEMIT.json",
        "exact_pframe": require_env("BILMINE_PFRAME",
                                    "the exact generator-frame JSON"),
        "c1_Pi": require_env("BILMINE_C1_PI",
                             "the C1 control's certified period 5-vector JSON"),
        "c1_N": require_env("BILMINE_C1_N",
                            "the C1 control's integer certificate JSON"),
        "parent_prereg": f"{PARENT}/PREREG_BILMINE.md",
        "parent_seal": f"{PARENT}/work/SEAL_RECORD.json",
        "parent_control_verdict": f"{PARENT}/work/CONTROL_VERDICT.json",
        "parent_exact_invariants": f"{PARENT}/work/EXACT_INVARIANTS.json",
        "parent_mine_p1": f"{PARENT}/work/MINE_P1.json",
        "parent_mine_m1": f"{PARENT}/work/MINE_M1.json",
        "parent_mine_w5": f"{PARENT}/work/MINE_W5.json",
        "parent_verdict": f"{PARENT}/work/BILMINE_VERDICT.json",
        "parent_lib": parent_lib,
        "prereg": f"{LEG}/PREREG_BILMINE_F1.md",
        "launch_line": f"{LEG}/LAUNCH_LINE.txt",
    }
    missing = [k for k, v in inputs.items() if not os.path.isfile(v)]
    if missing:
        print("seal_prereg_f1: REFUSING — missing input(s) of record:",
              ", ".join(f"{k} ({inputs[k]})" for k in missing))
        return 4
    stamp = subprocess.check_output(["date", "-u"]).decode().strip()
    rec = {
        "seal": ("BILMINE-F1 prereg seal — written before any control or mine "
                 "run; no relation number exists in this leg at this stamp"),
        "stamp_utc": stamp,
        "pins": {k: {"path": v, "sha256": sha(v)} for k, v in inputs.items()},
        "producer": {
            "leg": "BILMINE-F1",
            "script": "seal_prereg_f1.py",
            "script_sha256": sha(os.path.abspath(__file__)),
            "stamp_utc": stamp,
        },
    }
    out = f"{LEG}/work/SEAL_RECORD.json"
    json.dump(rec, open(out, "w"), indent=1)
    rc = 0
    lint = os.environ.get("BILMINE_PRODUCER_LINT")
    if lint:
        r = subprocess.run([sys.executable, lint, "--check", out],
                           capture_output=True, text=True)
        print(r.stdout.strip() or r.stderr.strip())
        rc = r.returncode
    print("SEALED", out, "prereg_sha256", rec["pins"]["prereg"]["sha256"])
    return rc


if __name__ == "__main__":
    sys.exit(main())

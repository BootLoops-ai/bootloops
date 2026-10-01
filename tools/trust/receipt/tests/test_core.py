#!/usr/bin/env python3
"""test_core.py — synthetic-system unit tests: verify/emit/detect semantics,
rc discipline (0/1/2), and the 6-class sabotage suite on a system small
enough to check by hand.

Sabotage class mapping to the archived sabotage battery (all caught at 2
primes there):
  S1 coefficient perturb      = archived S1
  S2 master-label swap        = archived S2 (the dictionary/change-of-basis bug)
  S3 global scale             = archived S3
  S4 lambda perturb           = archived S6
  S5 lambda support drop      (new here: truncated witness)
  S6 wrong-witness-for-target (new here: witness/target association bug)
The archived battery's S4 wrong-twist / S5 wrong-basis-element were
ENGINE-level sabotages; their table-level product is exactly the wrong-table
exhibit, covered in test_banked.py as the canonical MUST-FAIL.

Run: python3 test_core.py <out_dir>   (writes scratch witnesses under out_dir)
"""
import json
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
TOOL = os.path.dirname(HERE)
sys.path.insert(0, TOOL)

import core                                     # noqa: E402
from core import (RC_FAIL, RC_MALFORMED, RC_OK,  # noqa: E402
                  load_witness, verify_row, write_witness)
from detector import detect                     # noqa: E402
from emitter import FpSystem, emit_for_target   # noqa: E402

P = 2147483647
# scratch goes under argv[1], EXCEPT when that resolves into the package
# source tree (under `pytest <pkg>` argv[1] is the package dir itself) —
# then, or with no argument, a fresh mkdtemp dir is used instead.
PKG = os.path.realpath(TOOL)
_out = sys.argv[1] if len(sys.argv) > 1 else None
if _out and os.path.realpath(_out) != PKG \
        and not os.path.realpath(_out).startswith(PKG + os.sep):
    OUT = _out
    os.makedirs(OUT, exist_ok=True)
else:
    OUT = tempfile.mkdtemp(prefix="receipt_test_core_")
FAILED = []


def check(name, cond, extra=""):
    tag = "PASS" if cond else "FAIL"
    print(f"[{tag}] {name} {extra}")
    if not cond:
        FAILED.append(name)


# ---------------------------------------------------------------- fixture
# columns: 101,102 = masters; 201,202,203 = reducible; system (mod P):
#   R0:  203 - 5*101            = 0
#   R1:  202 - 3*203 - 7*102    = 0
#   R2:  201 - 2*202 - 101      = 0
#   R3:  201 + 202 + 203 - X    (X chosen to make it dependent-free filler:
#        use an independent combo: 201 - 11*102)
ROWS = [
    {203: 1, 101: (-5) % P},
    {202: 1, 203: (-3) % P, 102: (-7) % P},
    {201: 1, 202: (-2) % P, 101: (-1) % P},
]
# ground truth closures:
#   203 = 5*101
#   202 = 3*203 + 7*102 = 15*101 + 7*102
#   201 = 2*202 + 101   = 31*101 + 14*102
TRUTH = {203: {101: 5}, 202: {101: 15, 102: 7}, 201: {101: 31, 102: 14}}
PIVOTS = [201, 202, 203]
MASTERS = {101, 102}
ALLCOLS = sorted({c for r in ROWS for c in r})


def emit_witness(target, path):
    import numpy as np
    rng = np.random.default_rng(7)
    S = FpSystem([list(r.items()) for r in ROWS], PIVOTS, ALLCOLS, P)
    rec, lam, _ = emit_for_target(S, target, MASTERS, rng,
                                  claimed=TRUTH[target])
    assert rec["status"] == "CERTIFIED", rec
    nz = np.flatnonzero(lam)
    write_witness(path, p=P, target_col=target, c=rec["c"],
                  lam_idx=[int(i) for i in nz],
                  lam_val=[int(lam[i]) for i in nz], n_rows=len(ROWS),
                  system={"n_rows": len(ROWS),
                          "fingerprint": core.system_fingerprint(ROWS, P)})
    return path


# ---------------------------------------------------- emit + verify positive
wp = {t: emit_witness(t, os.path.join(OUT, f"w_{t}.json"))
      for t in (201, 202, 203)}
for t in (201, 202, 203):
    wit = load_witness(wp[t])
    ok, det = verify_row(ROWS, wit,
                         check_fingerprint=core.system_fingerprint(ROWS, P))
    check(f"positive_control_verify_{t}", ok, det.get("reason", ""))
    check(f"emit_c_matches_truth_{t}",
          wit["c"] == {k: v % P for k, v in TRUTH[t].items()})

# emit CLAIM_MISMATCH: wrong claimed row must be refused a witness
import numpy as np                                # noqa: E402
rng = np.random.default_rng(8)
S = FpSystem([list(r.items()) for r in ROWS], PIVOTS, ALLCOLS, P)
bad_claim = dict(TRUTH[202]); bad_claim[101] = 16   # noqa: E702
rec, lam, _ = emit_for_target(S, 202, MASTERS, rng, claimed=bad_claim)
check("emit_claim_mismatch_detected", rec["status"] == "CLAIM_MISMATCH", rec["status"])

# ------------------------------------------------------------ sabotage suite
wit0 = load_witness(wp[202])


def mutate(**kw):
    w = json.loads(json.dumps(json.load(open(wp[202]))))
    for k, v in kw.items():
        w[k] = v
    q = os.path.join(OUT, "w_mut.json")
    json.dump(w, open(q, "w"))
    return load_witness(q)


# S1 coefficient perturb
w = mutate(c={"101": 16, "102": 7})
ok, _ = verify_row(ROWS, w)
check("S1_coefficient_perturb_caught", not ok)
# S2 master-label swap
w = mutate(c={"101": 7, "102": 15})
ok, _ = verify_row(ROWS, w)
check("S2_label_swap_caught", not ok)
# S3 global scale x2
w = mutate(c={"101": 30, "102": 14})
ok, _ = verify_row(ROWS, w)
check("S3_global_scale_caught", not ok)
# S4 lambda perturb
raw = json.load(open(wp[202]))
raw["lam"]["val"][0] = (raw["lam"]["val"][0] + 1) % P
json.dump(raw, open(os.path.join(OUT, "w_s4.json"), "w"))
ok, _ = verify_row(ROWS, load_witness(os.path.join(OUT, "w_s4.json")))
check("S4_lambda_perturb_caught", not ok)
# S5 lambda support drop (truncated witness)
raw = json.load(open(wp[202]))
raw["lam"]["idx"] = raw["lam"]["idx"][:-1]
raw["lam"]["val"] = raw["lam"]["val"][:-1]
json.dump(raw, open(os.path.join(OUT, "w_s5.json"), "w"))
ok, _ = verify_row(ROWS, load_witness(os.path.join(OUT, "w_s5.json")))
check("S5_lambda_truncation_caught", not ok)
# S6 wrong witness for target (association bug): witness of 203 claimed for 202
raw = json.load(open(wp[203]))
raw["target"]["col"] = 202
json.dump(raw, open(os.path.join(OUT, "w_s6.json"), "w"))
ok, _ = verify_row(ROWS, load_witness(os.path.join(OUT, "w_s6.json")))
check("S6_wrong_witness_for_target_caught", not ok)

# table-mismatch mode: valid witness + wrong TABLE row
ok, det = verify_row(ROWS, wit0, table_row={101: 16, 102: 7})
check("table_mismatch_caught", not ok and "TABLE" in det.get("reason", ""))
# fingerprint mismatch
ok, det = verify_row(ROWS, wit0, check_fingerprint="deadbeef")
check("fingerprint_mismatch_caught", not ok)
# lambda length mismatch = FAIL (not malformed): claims wrong system size
raw = json.load(open(wp[202]))
raw["lam"]["n_rows"] = 4
json.dump(raw, open(os.path.join(OUT, "w_len.json"), "w"))
ok, det = verify_row(ROWS, load_witness(os.path.join(OUT, "w_len.json")))
check("n_rows_mismatch_fails", not ok)

# --------------------------------------------------------------- rc via CLI
sysj = os.path.join(OUT, "system.jsonl")
with open(sysj, "w") as fh:
    for r in ROWS:
        fh.write(json.dumps({str(k): v for k, v in r.items()}) + "\n")

env = dict(os.environ)


def cli(*a):
    return subprocess.run([sys.executable,
                           os.path.join(TOOL, "receipt.py"), *a],
                          capture_output=True, text=True, env=env)

r = cli("verify", "--system-jsonl", sysj, "--witness", wp[201], wp[202],
        wp[203], "-q")
check("rc0_all_pass", r.returncode == RC_OK, f"rc={r.returncode}")
r = cli("verify", "--system-jsonl", sysj, "--witness", wp[201],
        os.path.join(OUT, "w_s4.json"), "-q")
check("rc1_any_fail", r.returncode == RC_FAIL, f"rc={r.returncode}")
badj = os.path.join(OUT, "bad.json")
open(badj, "w").write("{not json")
r = cli("verify", "--system-jsonl", sysj, "--witness", badj, "-q")
check("rc2_malformed", r.returncode == RC_MALFORMED, f"rc={r.returncode}")
r = cli("verify", "--system-jsonl", sysj, "--witness",
        os.path.join(OUT, "nonexistent_*.json"), "-q")
check("rc2_missing_witness", r.returncode == RC_MALFORMED,
      f"rc={r.returncode}")
missing = os.path.join(OUT, "w_missing_fields.json")
json.dump({"receipt_version": "1.0", "p": P}, open(missing, "w"))
r = cli("verify", "--system-jsonl", sysj, "--witness", missing, "-q")
check("rc2_missing_fields", r.returncode == RC_MALFORMED,
      f"rc={r.returncode}")

# ------------------------------------------------------------------ detect
alarm, rep = detect(ROWS, P, MASTERS, targets=[201, 202, 203])
check("detect_clean_no_alarm", not alarm, json.dumps(rep)[:120])
# false master declared (202 is reducible): A1 or A2 must fire
alarm, rep = detect(ROWS, P, MASTERS | {202}, targets=[201, 203])
check("detect_false_master_alarm", alarm,
      f"A1={rep['A1_master_relations']['fired']} "
      f"A2={rep['A2_pseudo_masters']['fired']}")
# true master omitted from declaration: survivors alarm
alarm, rep = detect(ROWS, P, {101}, targets=[201])
check("detect_omitted_master_alarm",
      alarm and rep["A2_pseudo_masters"]["fired"])
# uncovered target
alarm, rep = detect(ROWS, P, MASTERS, targets=[999])
check("detect_uncovered_target_alarm",
      alarm and rep["A3_uncovered_targets"]["fired"])
# foreign support in claimed table
alarm, rep = detect(ROWS, P, MASTERS, table={201: {101: 31, 555: 1}})
check("detect_foreign_support_alarm",
      alarm and rep["A4_foreign_support"]["fired"])

print(f"\ntest_core: {'ALL PASS' if not FAILED else 'FAILURES: ' + str(FAILED)}")
sys.exit(0 if not FAILED else 1)

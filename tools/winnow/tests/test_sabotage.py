#!/usr/bin/env python3
"""Sabotage suite + the canonical wrong-table MUST-FAIL, through the
VENDORED receipt component (ibplapper.receipt).

Leg A — 6-class synthetic sabotage (6/6 at 2 primes):
  S1 coefficient perturb, S2 master-label swap, S3 global scale,
  S4 lambda perturb, S5 lambda truncation, S6 wrong-witness-for-target
  — every one must FAIL verify; the untampered witnesses must PASS.

Leg B — the recorded confident-wrong-table exhibit (the design reason receipts
are native; archived — this leg skips without WINNOW_BANKED_ROOT):
3 rows x 2 primes of D-ladder-stable,
rank-clean, 2-prime-CRT-consistent WRONG values. Witnesses emitted by
ibplapper's retrofit path; the wrong table rows MUST FAIL verify_row
(table_row mode) and the independently evaluated TRUE rows must PASS as
positive control. 6/6 + 6/6 required.
"""
import _common  # noqa: F401  MUST be first: sys.path bootstrap for the package

import json
import os

import ibplapper as lap
from ibplapper.adapters import kira
from ibplapper.adapters.kira import CoeffEvaluator
from ibplapper.receipt import core
from ibplapper.receipt.core import load_witness, verify_row, write_witness
from ibplapper.witness import retrofit
from _common import (COPAIR, D4, D0, E0, NOMAT, PRIMES, REGISTRY, check,
                     finish, out_dir, skip_without)
from _oracle import anchor_table, classify, eval_vec, parse_oracle_table

OUT = out_dir("/tmp/ibplapper_test_sabotage")

# ------------------------------------------------- Leg A: synthetic 6-class
P = 2147483647
ROWS = [{203: 1, 101: (-5) % P},
        {202: 1, 203: (-3) % P, 102: (-7) % P},
        {201: 1, 202: (-2) % P, 101: (-1) % P}]
TRUTH = {203: {101: 5}, 202: {101: 15, 102: 7}, 201: {101: 31, 102: 14}}

sys_ = lap.System(ROWS, P, {201: 3, 202: 2, 203: 1}, forbid={101, 102})
res = lap.eliminate(sys_, lap.Schedule(policy="B2FT"), witnesses="retrofit",
                    witness_targets=[201, 202, 203],
                    witness_dir=os.path.join(OUT, "wits"))
for t in (201, 202, 203):
    rec = res.witnesses[t]
    check(f"legA_positive_control::{t}",
          rec["status"] == "CERTIFIED" and rec["verify_pass"]
          and rec["c"] == TRUTH[t])

wp = {t: res.witnesses[t]["witness_path"] for t in (201, 202, 203)}


def mutated(base_path, fn, name):
    raw = json.load(open(base_path))
    fn(raw)
    q = os.path.join(OUT, f"w_{name}.json")
    json.dump(raw, open(q, "w"))
    ok, _ = verify_row(ROWS, load_witness(q))
    check(f"legA_{name}_caught", not ok)


mutated(wp[202], lambda w: w.update(c={"101": 16, "102": 7}), "S1_perturb")
mutated(wp[202], lambda w: w.update(c={"101": 7, "102": 15}), "S2_labelswap")
mutated(wp[202], lambda w: w.update(c={"101": 30, "102": 14}), "S3_scale")
mutated(wp[202], lambda w: w["lam"]["val"].__setitem__(
    0, (w["lam"]["val"][0] + 1) % P), "S4_lam_perturb")


def _truncate(w):
    w["lam"]["idx"] = w["lam"]["idx"][:-1]
    w["lam"]["val"] = w["lam"]["val"][:-1]


mutated(wp[202], _truncate, "S5_lam_truncate")
mutated(wp[203], lambda w: w["target"].__setitem__("col", 202),
        "S6_wrong_target")

# table-row mode + fingerprint mode
wit0 = load_witness(wp[202])
ok, det = verify_row(ROWS, wit0, table_row={101: 16, 102: 7})
check("legA_table_mismatch_caught", not ok and "TABLE" in det.get("reason", ""))
ok, _ = verify_row(ROWS, wit0, check_fingerprint="deadbeef")
check("legA_fingerprint_mismatch_caught", not ok)

# ----------------------------------- Leg B: recorded wrong-table exhibit
if COPAIR is None:
    print("SKIP leg B: the archived wrong-table exhibit needs WINNOW_BANKED_ROOT")
    finish("test_sabotage")
GEN = os.path.join(NOMAT, "gen_lbl3m2L_k2disp")
STG = os.path.join(D4, "lbl3m2L_k2disp")
FAM = "lbl3m2L_k2disp"
SEC119 = ["[1, 1, 1, 0, 1, 1, 2, 0, 0]", "[1, 1, 2, 0, 1, 1, 2, 0, 0]",
          "[1, 2, 1, 0, 1, 1, 2, 0, 0]"]
pair = json.load(open(os.path.join(COPAIR, "PAIR_RESULT.json")))
oracle_blocks = json.load(open(os.path.join(COPAIR, "ORACLE_BLOCKS.json")))
orc = {str(list(t)): terms for t, terms in parse_oracle_table(
    os.path.join(STG, "corpus", REGISTRY, "kira_target.m"))}

n_mustfail = n_failed = n_ctl = n_ctl_pass = 0
for p in PRIMES:
    system, sysd = kira.load(GEN, FAM, p, D0, E0)
    res = lap.eliminate(system, lap.Schedule(policy="B2F"))
    info = classify(sysd["sector_of"], STG, GEN, set(res.subs), FAM)
    masters_map = info["masters_map"]
    midx2w = {idx: w for w, idx in masters_map.items()}
    ev = CoeffEvaluator(p, D0, E0)
    # anchor the three exhibit targets
    orc3 = [(tuple(json.loads(t)), orc[t]) for t in SEC119]
    targets, skip = anchor_table(orc3, info, res.subs, sysd["sector_of"],
                                 ev, p)
    check(f"legB_exhibit_targets_anchored::{p}", len(targets) == 3,
          f"{len(targets)} skip={skip}")
    wits = retrofit(system, res.subs, [w for _t, _tm, w in targets],
                    masters=set(masters_map))
    sec_blocks = {s["sector"]: s for s in pair["sectors"]
                  if s["p"] == p and s["sector"] == 119}
    blocks = sec_blocks[119]["blocks"]
    for (tgt, terms, w_t) in targets:
        tstr = str(list(tgt))
        rec = wits[w_t]
        check(f"legB_witness_certified::{tstr}|{p}",
              rec["status"] == "CERTIFIED" and rec["verify_pass"])
        wit = {"p": p, "target_col": w_t, "c": rec["c"],
               "lam_idx": sorted(rec["lam"]),
               "lam_val": [rec["lam"][i] for i in sorted(rec["lam"])],
               "n_rows": len(system.rows), "meta": {}, "path": "<mem>"}
        # TRUE row: independent oracle evaluation (c = -vec convention)
        vec = eval_vec(terms, ev, midx2w, p)
        true_c = {w: (-v) % p for w, v in vec.items()}
        # convention gate vs the recorded in-sector oracle blocks
        blk_true = {midx2w[tuple(json.loads(m))]: (-v) % p
                    for m, v in oracle_blocks[f"{tstr}|{p}"].items()}
        check(f"legB_block_convention_gate::{tstr}|{p}",
              all(true_c.get(w) == v for w, v in blk_true.items()))
        # WRONG row: overwrite in-sector entries with the pairing's values
        wrong_c = dict(true_c)
        for m, v in blocks[tstr].items():
            wrong_c[midx2w[tuple(json.loads(m))]] = (-v) % p
        assert wrong_c != true_c, f"pairing row equals oracle row? {tstr}|{p}"
        ok, det = verify_row(system.rows, wit, table_row=wrong_c)
        n_mustfail += 1
        n_failed += (not ok)
        check(f"legB_MUSTFAIL_wrong_row::{tstr}|{p}", not ok,
              det.get("reason", "")[:60])
        ok2, det2 = verify_row(system.rows, wit, table_row=true_c)
        n_ctl += 1
        n_ctl_pass += bool(ok2)
        check(f"legB_control_true_row_passes::{tstr}|{p}", ok2,
              det2.get("reason", ""))

check("legB_all_6_wrong_rows_failed", n_mustfail == 6 and n_failed == 6,
      f"{n_failed}/{n_mustfail}")
check("legB_all_6_controls_passed", n_ctl == 6 and n_ctl_pass == 6,
      f"{n_ctl_pass}/{n_ctl}")

finish("test_sabotage")

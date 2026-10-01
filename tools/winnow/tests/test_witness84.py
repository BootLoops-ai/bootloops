#!/usr/bin/env python3
"""84/84 witness round-trip — the archived retrofit battery re-run through
the ibplapper API (archived record: RETROFIT_lbl3m2L_k2disp.json, 42
registry rows x 2 primes, oracle_match + independent-checker pass).
Replays archived artifacts — skips without WINNOW_BANKED_ROOT.

Per prime: load k2disp via the kira adapter -> eliminate (B2F) -> anchor the
42-row registry oracle table to closed rows (harness-side anchoring,
_oracle.py) -> check ORACLE MATCH (our closed row == independently evaluated
kira table row) -> retrofit-emit a lambda witness per row -> the vendored
RECEIPT core must verify every one against the loaded rows.
Bar: 84/84 oracle-match AND 84/84 CERTIFIED+verified. Witness files are
written to scratch and re-verified from disk (full round-trip).
"""
import _common  # noqa: F401  MUST be first: sys.path bootstrap for the package

import glob
import os
import time

import ibplapper as lap
from ibplapper.adapters import kira
from ibplapper.adapters.kira import CoeffEvaluator
from ibplapper.receipt import core
from ibplapper.witness import retrofit
from _common import (D4, D0, E0, NOMAT, PRIMES, REGISTRY, check, finish,
                     out_dir, skip_without)

skip_without("test_witness84", NOMAT, D4)
from _oracle import anchor_table, classify, eval_vec, parse_oracle_table

OUT = out_dir("/tmp/ibplapper_test_witness84")
GEN = os.path.join(NOMAT, "gen_lbl3m2L_k2disp")
STG = os.path.join(D4, "lbl3m2L_k2disp")
FAM = "lbl3m2L_k2disp"

orc = parse_oracle_table(os.path.join(STG, "corpus", REGISTRY,
                                      "kira_target.m"))
check("registry_table_42_rows", len(orc) == 42, f"n={len(orc)}")

targets_by_w = None
n_match = n_cert = n_verify = n_total = 0
t_solve_total = 0.0
for p in PRIMES:
    system, sysd = kira.load(GEN, FAM, p, D0, E0)
    res = lap.eliminate(system, lap.Schedule(policy="B2F"))
    check(f"no_leftover::{p}", not res.leftover)
    info = classify(sysd["sector_of"], STG, GEN, set(res.subs), FAM)
    check(f"partition_square_complete::{p}",
          len(info["pivots"]) == len(system.rows)
          and info["n_unclassified"] == 0,
          f"pivots={len(info['pivots'])} rows={len(system.rows)} "
          f"uncls={info['n_unclassified']}")
    ev = CoeffEvaluator(p, D0, E0)
    if targets_by_w is None:        # anchor once, at the anchor slice
        targets, skip = anchor_table(orc, info, res.subs,
                                     sysd["sector_of"], ev, p)
        check("anchored_42_of_42", len(targets) == 42,
              f"anchored={len(targets)} skip={skip}")
        targets_by_w = targets
    midx2w = {idx: w for w, idx in info["masters_map"].items()}
    # oracle match: our closed row == independently evaluated table row
    for tgt, terms, w_t in targets_by_w:
        vec = eval_vec(terms, ev, midx2w, p)
        n_total += 1
        if res.subs.get(w_t) == vec:
            n_match += 1
        else:
            check(f"ORACLE_MISMATCH::{tgt}|{p}", False)
    # retrofit witnesses (shell columns excluded from the coefficient basis)
    wdir = os.path.join(OUT, f"w_{p}")
    t0 = time.time()
    wits = retrofit(system, res.subs, [w for _t, _tm, w in targets_by_w],
                    masters=set(info["masters_map"]), write_dir=wdir,
                    family=FAM, point={"d": D0, "eta": E0})
    t_solve_total += time.time() - t0
    for w_t, rec in wits.items():
        n_cert += rec["status"] == "CERTIFIED"
    # full round-trip: re-load every witness FILE, re-verify via the
    # vendored solver-agnostic core against the adapter-loaded rows
    t0 = time.time()
    wfiles = sorted(glob.glob(os.path.join(wdir, "*.json")))
    for q in wfiles:
        ok, det = core.verify_row(system.rows, core.load_witness(q))
        n_verify += bool(ok)
    t_pure = time.time() - t0
    print(f"[measure] p={p}: verify {len(wfiles)} rows in {t_pure:.3f}s "
          f"({1000*t_pure/max(1,len(wfiles)):.2f} ms/row; recorded prior "
          f"0.96 ms/row)")

check("oracle_match_84_of_84", n_match == 84 and n_total == 84,
      f"{n_match}/{n_total}")
check("witness_certified_84_of_84", n_cert == 84, f"{n_cert}/84")
check("witness_verified_84_of_84", n_verify == 84, f"{n_verify}/84")
print(f"[measure] retrofit solve total {t_solve_total:.1f}s for 84 rows "
      f"({t_solve_total/84:.2f} s/row; recorded prior 2.2-2.3 s/row/prime "
      f"fresh, min-poly reuse 5.5x)")

finish("test_witness84")

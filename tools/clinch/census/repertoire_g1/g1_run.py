"""Gate-1 driver: bench / batch / assemble.

bench    — planted instrument controls (PLANT_BISTABLE truth=2, PLANT_MAC
           truth=1) + named grid cells inline, single core; phase timings =
           the T2 (census rate on the anchored family), T3 (pruning shares),
           T4 (certified-layer per-cell cost vs the 1 core-h bar), C0
           (measured support/equilibrium density) receipts. HALT rule: bench
           cell over T4 bar => no batch, re-budget.
batch    — all grid cells, multiprocessing width from --width, per-cell wall
           cap from --cellcap (basis: measured bench cells x slip).
assemble — read cells/*.json -> GATE1_VERDICT.json (decision rule verbatim
           from the study's plan) + K-control equality check.
"""
import argparse
import json
import multiprocessing as mp
import os
import sys
import time

import numpy as np

import g1_model as M
import g1_census as CEN

# Receipts land OUTSIDE the shipped tree: $G1_OUT, else ./g1_out under the
# caller's cwd — never inside the package dir.
OUT = os.environ.get("G1_OUT") or os.path.join(os.getcwd(), "g1_out")


def _write_cell(rec, disp, outdir):
    """Atomic per-cell receipt writes (checkpoint-restart discipline:
    a <cell>.json existing == cell done; tmp+os.replace makes that decidable
    after a seat/process death — no truncated receipt can be mistaken for a
    completed cell)."""
    cid = rec["cell_id"]
    os.makedirs(outdir, exist_ok=True)
    jtmp = os.path.join(outdir, f"{cid}.json.tmp")
    with open(jtmp, "w") as fh:
        json.dump(rec, fh, indent=1)
    ntmp = os.path.join(outdir, f"{cid}_disp.tmp.npz")
    np.savez_compressed(ntmp, disp=disp)
    os.replace(ntmp, os.path.join(outdir, f"{cid}_disp.npz"))
    os.replace(jtmp, os.path.join(outdir, f"{cid}.json"))


def _hard_class(cell):
    """The measured hard class (h=2 x lo-loss, m=8
    block, 16/140 cells) — carries the adopted 3.6 core-h bracket cap."""
    return cell["m"] == 8 and cell["h"] == 2 and cell["loss"] == "lo"


def run_one(args):
    cell, cap, outdir = args
    t0 = time.time()
    try:
        rec, disp = CEN.run_cell(cell, wall_cap_s=cap)
    except Exception as e:  # loud, honest failure -> UNCERTIFIED record
        rec, disp = ({"cell_id": cell["cell_id"], "status": "UNCERTIFIED",
                      "reason": f"EXCEPTION {type(e).__name__}: {e}"},
                     np.zeros(0, np.uint8))
    rec["wall_total_s"] = round(time.time() - t0, 2)
    _write_cell(rec, disp, outdir)
    return (rec["cell_id"], rec.get("status"), rec.get("n_stable"),
            rec.get("completeness_closed"), rec.get("wall_s"))


def cmd_bench(ns):
    outdir = os.path.join(OUT, "bench")
    cells = {c["cell_id"]: c for c in M.build_plants() + M.build_grid()}
    ids = ns.ids or ["PLANT_BISTABLE", "PLANT_MAC",
                     "m08_h1_w1.0_sp0_a0.0_K1.0_Lslice_Dhi_A",
                     "m08_h1_w0.1_sp1_a0.5_K1.0_Lslice_Dhi_A"]
    summary = []
    for cid in ids:
        cell = cells[cid]
        print(f"[bench] {cid} start", flush=True)

        def log(msg):
            print(f"[bench] {cid} {msg}", flush=True)
        t0 = time.time()
        rec, disp = CEN.run_cell(cell, wall_cap_s=ns.cellcap, log=log)
        rec["wall_total_s"] = round(time.time() - t0, 2)
        _write_cell(rec, disp, outdir)
        truth = cell.get("truth_stable")
        row = {"cell_id": cid, "status": rec.get("status"),
               "n_stable": rec.get("n_stable"),
               "truth_stable": truth,
               "control_pass": (None if truth is None
                                else rec.get("n_stable") == truth
                                and rec.get("completeness_closed")),
               "n_supports": rec.get("n_supports_enumerated"),
               "dispositions": rec.get("dispositions"),
               "n_equilibria": rec.get("n_equilibria"),
               "completeness_closed": rec.get("completeness_closed"),
               "wall_s": rec.get("wall_s"),
               "boxes_total": rec.get("boxes_total"),
               "timings": rec.get("timings")}
        print(json.dumps(row), flush=True)
        summary.append(row)
    hard = [r for r in summary if r["cell_id"].startswith("m08")]
    t4_bar_s = 3600.0
    bench = {"seed": M.SEED_BASE, "prec_bits": CEN.PREC,
             "rows": summary,
             "controls_pass": all(r["control_pass"] for r in summary
                                  if r["control_pass"] is not None),
             "T4_bar_core_h_per_cell": 1.0,
             "T4_worst_bench_cell_s": max((r["wall_s"] or 0) for r in hard)
             if hard else None,
             "T4_within_bar": all((r["wall_s"] or 0) <= t4_bar_s
                                  for r in hard) if hard else None}
    with open(os.path.join(outdir, "BENCH_RECEIPT.json"), "w") as fh:
        json.dump(bench, fh, indent=1)
    print(json.dumps({"BENCH": bench["controls_pass"],
                      "T4_within_bar": bench["T4_within_bar"]}), flush=True)


def cmd_batch(ns):
    outdir = os.path.join(OUT, "cells")
    cells = M.build_grid()
    done = set()
    if os.path.isdir(outdir):
        done = {f[:-5] for f in os.listdir(outdir) if f.endswith(".json")}
    todo = [c for c in cells if c["cell_id"] not in done]
    # v3: hard class (h2 x lo-loss) first — the 16 long poles pack the pool
    # head instead of straggling the tail; scheduling only, verdict-free.
    todo.sort(key=lambda c: (0 if _hard_class(c) else 1, c["cell_id"]))
    n_hard = sum(1 for c in todo if _hard_class(c))
    print(f"[batch] {len(cells)} cells, {len(done)} done, {len(todo)} to run "
          f"({n_hard} hard-class), width={ns.width}, cellcap={ns.cellcap}s, "
          f"cellcap_hard={ns.cellcap_hard}s", flush=True)
    args = [(c, ns.cellcap_hard if _hard_class(c) else ns.cellcap, outdir)
            for c in todo]
    t0 = time.time()
    with mp.get_context("spawn").Pool(ns.width, maxtasksperchild=4) as pool:
        for i, res in enumerate(pool.imap_unordered(run_one, args)):
            print(f"[batch] {i+1}/{len(args)} {res} "
                  f"elapsed={time.time()-t0:.0f}s", flush=True)
    print(f"[batch] DONE wall={time.time()-t0:.0f}s", flush=True)


def cmd_assemble(ns):
    outdir = os.path.join(OUT, "cells")
    recs = []
    for f in sorted(os.listdir(outdir)):
        if f.endswith(".json"):
            recs.append(json.load(open(os.path.join(outdir, f))))
    man = M.cells_manifest(M.build_grid())
    n_cells = len(man)
    ok = [r for r in recs if r.get("status") == "OK"]
    closed = [r for r in ok if r.get("completeness_closed")]
    m8 = [r for r in closed if r.get("axes", {}).get("m") == 8]
    tot_eq = sum(r["n_equilibria"] for r in ok)
    tot_unk = sum(r["n_unknown_sign"] for r in ok)
    stable_counts = {r["cell_id"]: r["n_stable"] for r in ok}
    multi = {cid: c for cid, c in stable_counts.items() if c > 1}
    coex_multi = {}
    for r in ok:
        ns_coex = sum(1 for e in r["equilibria"]
                      if e["stability"] == "STABLE" and len(e["S"]) > 0)
        if ns_coex > 1:
            coex_multi[r["cell_id"]] = ns_coex
    # K-control law (L2): counts must match K=1 partners exactly
    kc = {}
    for r in ok:
        if r["cell_id"].endswith("_KC"):
            ax = r["axes"]
            partner = (f"m08_h{ax['h']}_w{ax['w']}_sp{ax['sp']}_a{ax['a']}"
                       f"_K1.0_L{ax['loss']}_D{ax['div']}_{ax['seed_tag']}")
            kc[r["cell_id"]] = {"partner": partner,
                                "n_stable": r["n_stable"],
                                "partner_n_stable":
                                stable_counts.get(partner),
                                "match": stable_counts.get(partner)
                                == r["n_stable"]}
    uncert_cells = [{"cell_id": r["cell_id"],
                     "status": r.get("status"),
                     "reason": r.get("reason"),
                     "uncert": r.get("dispositions", {}).get("UNCERT"),
                     "excl_open": r.get("excl", {}).get("n_open"),
                     "capped": r.get("capped")}
                    for r in recs
                    if r.get("status") != "OK"
                    or not r.get("completeness_closed")]
    unk_frac = (tot_unk / tot_eq) if tot_eq else 0.0
    gate1 = {
        "bar_certificate_cells": len(closed),
        "bar_m8_certified_cells": len(m8),
        "bar_m8_required": 64,
        "bar_unknown_frac": unk_frac,
        "bar_unknown_required": 0.01,
        "GO": len(m8) >= 64 and unk_frac <= 0.01,
    }
    every_closed_one = all(r["n_stable"] == 1 for r in closed)
    verdict = {
        "date_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "n_cells_designed": n_cells, "n_cells_run": len(recs),
        "n_cells_ok": len(ok), "n_cells_certified_complete": len(closed),
        "uncertified_cells_verbatim": uncert_cells,
        "total_equilibria": tot_eq, "total_unknown_sign": tot_unk,
        "gate1_bars": gate1,
        "stable_count_by_cell": stable_counts,
        "cells_with_multiple_stable_states": multi,
        "cells_with_multiple_stable_COEXISTENCE_states": coex_multi,
        "one_everywhere_over_certified_cells": every_closed_one,
        "k_control_law": kc,
        "decision_rule": ("PLAN Sec.3/pitch Gate-1 verbatim: completeness "
                         "certificate + >=99% certified stability signs at "
                         "n>=20, m=8, >=64 cells, within 8 core-days; count "
                         "ONE everywhere => exclusion branch; >1 anywhere "
                         "=> alternative-repertoires branch"),
    }
    with open(os.path.join(OUT, "GATE1_VERDICT.json"), "w") as fh:
        json.dump(verdict, fh, indent=1)
    print(json.dumps({k: verdict[k] for k in
                      ("n_cells_certified_complete", "gate1_bars",
                       "cells_with_multiple_stable_states",
                       "one_everywhere_over_certified_cells")}, indent=1),
          flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("bench")
    b.add_argument("--ids", nargs="*", default=None)
    b.add_argument("--cellcap", type=float, default=3600.0)
    t = sub.add_parser("batch")
    t.add_argument("--width", type=int, default=24)
    t.add_argument("--cellcap", type=float, default=3600.0)
    t.add_argument("--cellcap-hard", dest="cellcap_hard", type=float,
                   default=12960.0)  # 3.6 core-h, adopted hard-class bracket
    sub.add_parser("assemble")
    ns = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    {"bench": cmd_bench, "batch": cmd_batch,
     "assemble": cmd_assemble}[ns.cmd](ns)

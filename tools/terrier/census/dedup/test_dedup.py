#!/usr/bin/env python3
"""Validation battery — dedup Engines A+B.

SYNTHETIC inputs only (no census data): fixed seeds in TEST_SEEDS.json;
ground truth is proof-level (budget invariant separates families;
membership is witnessed by construction). Battery:
  T1 generators exact-symplectic (all three stacks validate)
  T2 budget invariance + inverse-pair probes (fixed seed)
  T3 Engine A partition == ground truth (no truncation)
  T4 Engine B partition == ground truth (no truncation)
  T5 every merge witness verifies (referee stack witness.py)
  T6 reconcile GREEN-AGREE + plant own-class + RAW/DEDUP + hist columns
  T7 disagreement drill: crippled Engine A -> RED-RECONCILED-MISS, halt,
     reconciled partition == truth (the under-merge path)
  T8 witness-tamper drill -> RED-WITNESS-FAIL (integrity halt path)
Receipts -> test_receipts.json (volatile, regenerated). Exit 1 on any
failure.
"""
import json
import os
import random
import sys

import engine_a
import engine_b
import gens_synth
import reconcile
import witness

HERE = os.path.dirname(os.path.abspath(__file__))
CHECKS = []

def check(name, ok, detail=""):
    CHECKS.append({"check": name, "pass": bool(ok), "detail": str(detail)})
    print("%s %s %s" % ("PASS" if ok else "FAIL", name, detail))
    return ok

def main():
    cfgj = json.load(open(os.path.join(HERE, "TEST_SEEDS.json")))
    seeds, cfg = cfgj["seeds"], cfgj["config"]
    gens = gens_synth.make_gens(seeds["gens_synth"], cfg["n"], cfg["n_gens"])
    S = witness.sig(cfg["n"])
    check("T1-symplectic-gens",
          all(witness.sympl_ok(M, S) for M in gens) and len(gens) == 2,
          "gens=%s" % gens)
    # T2: budget invariance + inverse pairs, fixed probe seed
    rng = random.Random(seeds["invariance_probe"])
    tab = witness.make_table(gens, cfg["n"])
    toks = gens_synth.alphabet(len(gens))
    inv_ok = True
    for _ in range(300):
        st = (tuple(rng.randrange(-9, 10) for _ in range(cfg["n"])),
              tuple(rng.randrange(-9, 10) for _ in range(cfg["n"])))
        b0 = witness.budget(st, S)
        tok = rng.choice(toks)
        nx = witness.step(st, tok, tab)
        back = witness.step(nx, reconcile.inv_word((tok,))[0], tab)
        if witness.budget(nx, S) != b0 or back != st:
            inv_ok = False
    check("T2-budget-invariance-and-inverses", inv_ok, "300 probes")
    # Battery inputs (fixed seeds; plants inserted PRE-dedup)
    states, truth, tags, labels = gens_synth.build_battery(
        seeds["families"], seeds["shuffle"], gens, cfg)
    print("battery: N_raw=%d classes_truth=%d plants=%d"
          % (len(states), len(truth), tags.count("plant")))
    resA = engine_a.dedup(states, gens)
    resB = engine_b.dedup_b(states, gens)
    # Snapshot the battery results the receipt block below is built from
    # (the receipt reads ONLY these bindings, never the drill results).
    states_syn, resA_syn, resB_syn = states, resA, resB
    pa = sorted(sorted(c) for c in resA["partition"])
    pb = sorted(sorted(c) for c in resB["partition"])
    check("T3-engineA-ground-truth",
          pa == truth and not resA["truncated"],
          "N_dedup_A=%d nodes=%d" % (resA["N_dedup"], resA["nodes"]))
    check("T4-engineB-ground-truth",
          pb == truth and not resB["truncated"],
          "N_dedup_B=%d nodes=%d" % (resB["N_dedup"], resB["nodes"]))
    wit_ok = True
    for res in (resA, resB):
        for i, st in enumerate(states):
            if not witness.verify(st, res["words"][i], res["reps"][i],
                                  gens, tab):
                wit_ok = False
    check("T5-all-witnesses-verify", wit_ok,
          "%d witnesses x 2 engines" % len(states))
    rep = reconcile.reconcile(states, gens, resA, resB, tags, labels)
    check("T6-reconcile-green",
          rep["status"] == "GREEN-AGREE" and not rep["halt"]
          and rep["plant_own_class"] and rep["N_raw"] == len(states)
          and rep["N_dedup_A"] == rep["N_dedup_B"] == len(truth)
          and "orbit_size_hist_by_label" in rep,
          "status=%s hist=%s" % (rep["status"],
                                 rep.get("orbit_size_hist_by_label")))
    # T7: disagreement path — crippled Engine A under-merges
    resA_cr = engine_a.dedup(states, gens, max_nodes=0)
    rep7 = reconcile.reconcile(states, gens, resA_cr, resB, tags, labels)
    part7 = rep7["partition"]
    check("T7-red-reconciled-miss",
          rep7["status"] == "RED-RECONCILED-MISS" and rep7["halt"]
          and part7 == truth and len(rep7["engine_miss_pairs"]["A"]) > 0
          and not rep7["witness_failures"],
          "missA=%d" % len(rep7["engine_miss_pairs"]["A"]))
    # T8: witness-tamper drill -> integrity halt (RED-WITNESS-FAIL)
    resBt = dict(resB)
    resBt["words"] = list(resB["words"])
    pa_cr = sorted(sorted(c) for c in resA_cr["partition"])
    j = reconcile._extra_pairs(pb, pa_cr, len(states))[0][1]
    resBt["words"][j] = resBt["words"][j] + ("T",)
    rep8 = reconcile.reconcile(states, gens, resA_cr, resBt, tags, labels)
    check("T8-red-witness-fail",
          rep8["status"] == "RED-WITNESS-FAIL" and rep8["halt"]
          and len(rep8["witness_failures"]) >= 1,
          "failures=%d" % len(rep8["witness_failures"]))
    receipt = {
        "module": "census/dedup",
        "inputs": "SYNTHETIC fixed-seed battery only (no census data)",
        "seeds": seeds, "config": cfg, "gens": gens,
        "N_raw": len(states_syn), "N_classes_truth": len(truth),
        "n_plants": tags.count("plant"),
        "engine_A": {"N_dedup": resA_syn["N_dedup"],
                     "nodes": resA_syn["nodes"],
                     "truncated": resA_syn["truncated"]},
        "engine_B": {"N_dedup": resB_syn["N_dedup"],
                     "nodes": resB_syn["nodes"],
                     "truncated": resB_syn["truncated"]},
        "agreement": {"status": rep["status"], "halt": rep["halt"],
                      "plant_own_class": rep.get("plant_own_class"),
                      "plant_rows": rep.get("plant_rows"),
                      "orbit_size_hist_by_label":
                          rep.get("orbit_size_hist_by_label")},
        "drill_T7": {"status": rep7["status"], "halt": rep7["halt"],
                     "engine_miss_A": len(rep7["engine_miss_pairs"]["A"]),
                     "reconciled_equals_truth": part7 == truth},
        "drill_T8": {"status": rep8["status"], "halt": rep8["halt"],
                     "witness_failures": len(rep8["witness_failures"])},
        "checks": CHECKS,
    }
    out = os.path.join(HERE, "test_receipts.json")
    with open(out, "w") as fh:
        json.dump(receipt, fh, indent=1, sort_keys=True)
    print("receipt -> %s" % out)
    ok = all(c["pass"] for c in CHECKS)
    print("BATTERY %s (%d/%d)" % ("GREEN" if ok else "RED",
                                  sum(c["pass"] for c in CHECKS),
                                  len(CHECKS)))
    return 0 if ok else 1

if __name__ == "__main__":
    sys.exit(main())

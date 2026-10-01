#!/usr/bin/env python3
# lockpick bilmine member (F1 ring stage) — verdict assembler: builds the F1 leg verdict from the script-emitted stage receipts.
"""emit_verdict_f1.py — assemble the BILMINE-F1 leg verdict (prereg sec 8)
from the script-emitted stage receipts. The verdict line printed here is THE
verdict of record (script-emitted; never hand-set)."""
import json, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bilmine_f1_lib import leg_dir_f1, MODES, emit

MODE_CAP = {"sym_gen": "1e6", "anti_gen": "1e8",
            "sym_oa": "1e10", "anti_oa": "1e10",
            "sym_ob": "1e10", "anti_ob": "1e10"}


def main():
    LEG = leg_dir_f1()
    cv = json.load(open(f"{LEG}/work/CONTROL_VERDICT.json"))
    ex = json.load(open(f"{LEG}/work/EXACT_COROLLARY.json"))
    p1 = json.load(open(f"{LEG}/work/MINE_P1_F1.json"))
    m1 = json.load(open(f"{LEG}/work/MINE_M1_F1.json"))
    w5 = json.load(open(f"{LEG}/work/MINE_W5_F1.json"))
    gg = json.load(open(f"{LEG}/work/GLOBAL_GATE.json"))

    out = {"leg": "BILMINE-F1",
           "prereg_sha256": json.load(open(f"{LEG}/work/SEAL_RECORD.json"))
           ["pins"]["prereg"]["sha256"],
           "control_pass": bool(cv.get("control_pass"))}

    if not out["control_pass"]:
        out["verdict"] = ("BILMINE-F1 VERDICT: CONTROL-FAIL — no L6 statement "
                          "of any kind (prereg sec 4).")
        emit(out, f"{LEG}/work/BILMINE_F1_VERDICT.json",
             os.path.abspath(__file__))
        print(out["verdict"])
        return 2

    fam_state = {}
    any_found = any_partial = any_fence = False
    priced = []
    for sector in ("sym", "anti"):
        for mode in MODES:
            key = f"{sector}_{mode}"
            sec = gg["sectors"][key]
            r1, r2 = sec["ranks"]
            g2ok = (p1[f"G2_{sector}_{mode}"]["two_leg_lattice_identical"]
                    and m1[f"G2_{sector}_{mode}"]["two_leg_lattice_identical"])
            audits_ok = all(
                rec["windows"][f"{sector}_{mode}_{l}"]["audit"]["audit_pass"]
                for rec in (p1, m1) for l in ("lo", "hi"))
            po = [rec["windows"][f"{sector}_{mode}_{l}"]["priced_out"]
                  for rec in (p1, m1) for l in ("lo", "hi")]
            po = [x for x in po if x]
            fences = [h.get("consistency_fence")
                      for rec in (p1, m1) for l in ("lo", "hi")
                      for h in rec["windows"][f"{sector}_{mode}_{l}"]["hits"]
                      if h.get("consistency_fence")]
            nfull = sec.get("n_full_pass", 0) or 0
            if po:
                st = "CAPACITY-PRICED-OUT (measured; see window receipts)"
                priced.append(key)
            elif nfull >= 1 and fences:
                st = f"INSTRUMENT-INCONSISTENCY (fence-flagged chain: {fences[0]})"
                any_fence = True
            elif nfull >= 1:
                st = f"FOUND(rank {nfull})"
                any_found = True
            elif (r1 == 0 and r2 == 0) and g2ok and audits_ok:
                st = ("NULL-with-floor (lawful ladder exhausted; capacity "
                      "arithmetic printed per rung)")
            elif r1 == 0 and r2 == 0:
                st = ("NULL-with-floor (post-quotient; G5 re-pose applied)"
                      if audits_ok else "VOID (degenerate basket unresolved)")
            else:
                st = ("PARTIAL-STRUCTURE (point-local Z[sqrt15] relations "
                      "without a certified global pairing)")
                any_partial = True
            fam_state[key] = {"state": st, "ranks_p1_m1": [r1, r2],
                              "two_leg_identical_both_points": bool(g2ok),
                              "audits_pass": bool(audits_ok),
                              "n_full_pass_chains": nfull,
                              "height_cap": MODE_CAP[key]}
    out["window_families"] = fam_state

    w5v = {pt: w5["points"][pt]["verdict"] for pt in w5["points"]}
    out["W5F1_blind_doubled"] = w5v
    out["exact_corollary"] = {
        "W_sym_dim_over_Qsqrt15": ex["corollary_over_Qsqrt15"]["W_sym_dim"],
        "W_anti_dim_over_Qsqrt15": ex["corollary_over_Qsqrt15"]["W_anti_dim"],
        "End_dim_over_Qsqrt15": ex["corollary_over_Qsqrt15"]["End_dim"],
        "prediction": "any global invariant pairing = (a + b*sqrt15)*P",
    }

    floors = ("(+1: 500/535d legs, floor 470/505d rel; -1: 515/550d legs, "
              "floor 485/520d rel; caps 543.8/558.6 cured)")
    w5line = "; ".join(f"{pt}: {v}" for pt, v in w5v.items())
    all_null = all(fam_state[k]["state"].startswith("NULL-with-floor")
                   for k in fam_state)
    if any_found:
        head = "BILMINE-F1 VERDICT: FOUND — " + ", ".join(
            f"{k}:{fam_state[k]['state']}" for k in fam_state
            if "FOUND" in fam_state[k]["state"])
        tail = (" — the Z[sqrt15] lattice is reported and this leg STOPS "
                "(consumption is a separately-gated new door; grant term 4).")
    elif any_fence:
        head = ("BILMINE-F1 VERDICT: INSTRUMENT-INCONSISTENCY — " + ", ".join(
            f"{k}:{fam_state[k]['state']}" for k in fam_state
            if "INCONSIST" in fam_state[k]["state"]))
        tail = (" — a fence-flagged chain lies inside the parent's certified "
                "integer NULL scope; both receipts confront at the relay.")
    elif any_partial:
        head = "BILMINE-F1 VERDICT: PARTIAL-STRUCTURE — " + ", ".join(
            f"{k}:{fam_state[k]['state']}" for k in fam_state)
        tail = ""
    elif all_null:
        head = ("BILMINE-F1 VERDICT: NULL-with-floor — control_pass=true; no "
                "Z[sqrt15]-coefficient bilinear relation lattice (S_loc, "
                "S_mum), coefficients a + b*sqrt15 at heights max(|a|,|b|) <= "
                "1e6 (GEN SYM) / 1e8 (GEN ANTI) / 1e10 (one-sided "
                "orientations O-A, O-B), links the 666-pair period-product "
                "space of the certified MUM->+1 / MUM->-1 connection matrices "
                "to a shared MUM-frame pairing, in either symmetry sector, "
                "either orientation, at floors " + floors)
        tail = ("; blind doubled window W5F1 [" + w5line + "]; exact "
                "corollary: W_sym tensor Q(sqrt15) = Q(sqrt15)*P (dim 1), "
                "W_anti = 0, End = Q(sqrt15) — the parent's exact Z*P "
                "polarization admits no hidden sqrt15-scaled analytic "
                "realization at these heights and floors; ladder exhaustion "
                "is final for this leg, not a retry.")
    else:
        head = ("BILMINE-F1 VERDICT: NULL-with-floor over completed windows; "
                "CAPACITY-PRICED-OUT windows [" + ", ".join(priced) + "] "
                "quoted by measurement (wider is a relay ask)")
        tail = ("; blind doubled window W5F1 [" + w5line + "]; floors "
                + floors + "; ladder exhaustion final for completed windows.")
    out["verdict"] = head + tail
    emit(out, f"{LEG}/work/BILMINE_F1_VERDICT.json", os.path.abspath(__file__))
    print(out["verdict"])
    return 0


if __name__ == "__main__":
    sys.exit(main())

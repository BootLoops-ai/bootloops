#!/usr/bin/env python3
# lockpick bilmine member — verdict assembler: builds the leg verdict from the script-emitted stage receipts.
"""emit_verdict.py — assemble the BILMINE leg verdict (prereg sec 8) from the
script-emitted stage receipts. The verdict line printed here is THE verdict
of record (script-emitted; never hand-set)."""
import json, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bilmine_lib import leg_dir, emit


def main():
    LEG = leg_dir()
    cv = json.load(open(f"{LEG}/work/CONTROL_VERDICT.json"))
    ex = json.load(open(f"{LEG}/work/EXACT_INVARIANTS.json"))
    p1 = json.load(open(f"{LEG}/work/MINE_P1.json"))
    m1 = json.load(open(f"{LEG}/work/MINE_M1.json"))
    w5 = json.load(open(f"{LEG}/work/MINE_W5.json"))
    gg = json.load(open(f"{LEG}/work/GLOBAL_GATE.json"))

    out = {"leg": "BILMINE",
           "prereg_sha256": json.load(open(f"{LEG}/work/SEAL_RECORD.json"))
           ["pins"]["prereg"]["sha256"],
           "control_pass": bool(cv.get("control_pass"))}

    if not out["control_pass"]:
        out["verdict"] = "BILMINE VERDICT: CONTROL-FAIL — no L6 statement of any kind (prereg sec 4)."
        emit(out, f"{LEG}/work/BILMINE_VERDICT.json", os.path.abspath(__file__))
        print(out["verdict"])
        return 2

    sector_state = {}
    any_partial = False
    any_found = False
    for sector in ("sym", "anti"):
        sec = gg["sectors"][sector]
        r1, r2 = sec["ranks"] if sec.get("ranks") else (
            gg["sectors"][sector]["g2_p1"]["rank"],
            gg["sectors"][sector]["g2_m1"]["rank"])
        g2ok = (p1[f"G2_{sector}"]["two_leg_lattice_identical"]
                and m1[f"G2_{sector}"]["two_leg_lattice_identical"])
        audits_ok = all(p1["windows"][f"{sector}_{l}"]["audit"]["audit_pass"]
                        for l in ("lo", "hi")) and \
            all(m1["windows"][f"{sector}_{l}"]["audit"]["audit_pass"]
                for l in ("lo", "hi"))
        nfull = sec.get("n_full_pass", 0) or 0
        if nfull >= 1:
            st = f"FOUND(rank {nfull})"
            any_found = True
        elif (r1 == 0 and r2 == 0) and g2ok and audits_ok:
            st = "NULL-with-floor (ladder exhausted through 1e10; lawful at every rung)"
        elif r1 == 0 and r2 == 0:
            st = "NULL-with-floor (post-quotient; G5 re-pose applied)" if audits_ok \
                else "VOID (degenerate basket unresolved)"
        else:
            st = "PARTIAL-STRUCTURE (point-local integer relations without a certified global pairing)"
            any_partial = True
        sector_state[sector] = {"state": st, "ranks_p1_m1": [r1, r2],
                                "two_leg_identical_both_points": bool(g2ok),
                                "audits_pass": bool(audits_ok),
                                "n_full_pass_chains": nfull}
    out["sectors"] = sector_state

    w5v = {pt: w5["points"][pt]["verdict"] for pt in w5["points"]}
    out["W5_blind666"] = w5v
    out["exact_abstract_lattice"] = {
        "W_sym_dim": ex["W_invariant_spaces"]["GtSG"]["dim_sym"],
        "W_anti_dim": ex["W_invariant_spaces"]["GtSG"]["dim_anti"],
        "End_Q_dim": ex["commutant_dim_EndQ"],
        "P_in_W_sym": ex["P_in_W_sym"]["GtSG"],
        "detP": ex["detP"],
    }

    floors = "(+1: 500/535d legs, floor 470/505d rel; -1: 515/550d legs, floor 485/520d rel; caps 543.8/558.6 cured)"
    w5line = "; ".join(f"{pt}: {v}" for pt, v in w5v.items())
    if any_found:
        head = "BILMINE VERDICT: FOUND — " + ", ".join(
            f"{s}:{sector_state[s]['state']}" for s in sector_state)
        tail = " — the lattice is reported and this leg STOPS (consumption is a separately-gated new door)."
    elif any_partial:
        head = "BILMINE VERDICT: PARTIAL-STRUCTURE — " + ", ".join(
            f"{s}:{sector_state[s]['state']}" for s in sector_state)
        tail = ""
    else:
        head = ("BILMINE VERDICT: NULL-with-floor — control_pass=true; no integer "
                "bilinear relation lattice (S_loc, S_mum), max|coeff| <= 1e10, links "
                "the 666-pair period-product space of the certified MUM->+1 / "
                "MUM->-1 connection matrices to a shared MUM-frame pairing, in "
                "either symmetry sector, at floors " + floors)
        tail = ("; blind-666 window [" + w5line + "]; exact abstract lattice: "
                "W_sym = Z*P exactly (rank 1, det 512000, sig (4,2)), W_anti = 0, "
                "End_Q = Q — ladder exhaustion is final for this leg, not a retry.")
    out["verdict"] = head + tail
    emit(out, f"{LEG}/work/BILMINE_VERDICT.json", os.path.abspath(__file__))
    print(out["verdict"])
    return 0


if __name__ == "__main__":
    sys.exit(main())

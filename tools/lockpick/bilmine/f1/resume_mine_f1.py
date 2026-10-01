#!/usr/bin/env python3
# lockpick bilmine member (F1 ring stage) — fresh-process resume helper: re-runs only interrupted mine stages with unchanged machinery.
"""resume_mine_f1.py — fresh-process resume of the F1 mine after an
interrupted run in which the p1 receipt MINE_P1_F1.json was already emitted
by the unchanged run_mine_f1 and STANDS. This script re-runs ONLY the
interrupted stages with the same unchanged machinery: mine_point("m1") ->
MINE_M1_F1.json, then W5F1 for both points -> MINE_W5_F1.json.
Instrument-lineage note rides each receipt (the parent remine_points.py
pattern)."""
import os, sys, time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bilmine_f1_lib import leg_dir_f1, emit
import run_mine_f1


def main():
    LEG = leg_dir_f1()
    stamp = run_mine_f1.check_control()
    t0 = time.time()
    rec = run_mine_f1.mine_point("m1")
    rec["control_pass_stamp"] = stamp
    rec["wall_s"] = round(time.time() - t0, 1)
    rec["resumed_after"] = ("interrupted run (p1 receipt already emitted by "
                            "the unchanged script and stands); fresh-process "
                            "re-run of m1 only, machinery unchanged")
    emit(rec, f"{LEG}/work/MINE_M1_F1.json", os.path.abspath(__file__))
    for sector in ("sym", "anti"):
        for mode in run_mine_f1.MODES:
            g2 = rec[f"G2_{sector}_{mode}"]
            print("m1", sector, mode, "rank", g2["rank"],
                  "two_leg_identical", g2["two_leg_lattice_identical"])
    w5 = {"window": "W5F1 blind doubled basket (prereg sec 5)", "points": {},
          "resumed_after": rec["resumed_after"]}
    for pt in ("p1", "m1"):
        w5["points"][pt] = run_mine_f1.w5f1_window(pt)
        print("W5F1", pt, w5["points"][pt]["verdict"])
    emit(w5, f"{LEG}/work/MINE_W5_F1.json", os.path.abspath(__file__))


if __name__ == "__main__":
    main()

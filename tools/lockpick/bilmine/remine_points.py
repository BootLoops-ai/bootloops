#!/usr/bin/env python3
# lockpick bilmine member — fresh-process re-emit helper for the mine receipts (so a patched mine_window verifiably executes).
"""remine_points.py — re-emit MINE_P1/MINE_M1 from a FRESH process so a
patched run_control.mine_window (e.g. the G5 quotient re-pose) verifiably
executes: a long-lived process can hold a stale import of the pre-patch
mine_window (launch/edit race), so the re-emit runs the patched machinery in
a new process for one clean lineage.
"""
import os, sys, time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bilmine_lib import leg_dir, emit
import run_mine


def main():
    LEG = leg_dir()
    stamp = run_mine.check_control()
    for pt, out in (("p1", "MINE_P1.json"), ("m1", "MINE_M1.json")):
        t0 = time.time()
        rec = run_mine.mine_point(pt)
        rec["control_pass_stamp"] = stamp
        rec["wall_s"] = round(time.time() - t0, 1)
        rec["supersedes"] = ("fresh-process re-emit: the patched mine_window "
                             "verifiably executed in a new process")
        emit(rec, f"{LEG}/work/{out}", os.path.abspath(__file__))
        for sector in ("sym", "anti"):
            g2 = rec[f"G2_{sector}"]
            print(pt, sector, "rank", g2["rank"],
                  "two_leg_identical", g2["two_leg_lattice_identical"])
        for k, w in rec["windows"].items():
            q = w["audit"].get("quotient")
            print(" ", k, "audit_pass", w["audit"]["audit_pass"],
                  "quotient_dropped", (q["dropped_rows"] if q else None))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# lockpick bilmine member (F1 ring stage) — sealed-window production miner driver (point x sector x mode window families + blind doubled window W5F1).
"""run_mine_f1.py — the two-point Q(sqrt15) bilinear production mine (prereg
sec 5; runs ONLY after CONTROL-PASS). 12 sector window families (point x
sector x mode{GEN,O-A,O-B}, lo/hi legs) + W5F1 doubled blind window
(pilot-priced). Emits work/MINE_P1_F1.json, work/MINE_M1_F1.json,
work/MINE_W5_F1.json under $BILMINE_F1_LEG."""
import json, os, sys, time

import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bilmine_f1_lib import (leg_dir_f1, t5_dir, MODES, emit, parse_c_matrix,
                            hnf_basis,
                            capacity_line_f1, mine_homogeneous, residuals,
                            passes_floor, mine_window_f1, w5f1_products,
                            build_rows_f1, entry_list, sym_pairs, anti_pairs)


def slice_check(C, sector, use_im):
    """Structural assertion (not a sealed bar): the O-A / O-B windows are
    exact column slices of the GEN window — [a_loc|mum_b] and [b_loc|mum_a]
    respectively. Built at reduced dps; mpf equality is exact at fixed dps."""
    fit, _h = entry_list(sector, 0)
    n_loc = len(sym_pairs()) if sector == "sym" else len(anti_pairs())
    n_fit = len(fit)
    g, _t = build_rows_f1(C, sector, fit, use_im, 120, "gen")
    oa, _t = build_rows_f1(C, sector, fit, use_im, 120, "oa")
    ob, _t = build_rows_f1(C, sector, fit, use_im, 120, "ob")
    for r in range(len(g)):
        ga = g[r][:n_loc] + g[r][2 * n_loc + n_fit:]
        gb = g[r][n_loc:2 * n_loc] + g[r][2 * n_loc:2 * n_loc + n_fit]
        assert oa[r] == ga, f"OA slice mismatch {sector} row {r}"
        assert ob[r] == gb, f"OB slice mismatch {sector} row {r}"
    return True

RUNGS = [100, 10 ** 4, 10 ** 6, 10 ** 8, 10 ** 10]
POINTS = {
    "p1": {"cap": 543, "d": {"lo": 500, "hi": 535}, "use_im": False,
           "d_capacity": 500},
    "m1": {"cap": 558, "d": {"lo": 515, "hi": 550}, "use_im": True,
           "d_capacity": 515},
}


def check_control():
    cv = json.load(open(f"{leg_dir_f1()}/work/CONTROL_VERDICT.json"))
    if not cv.get("control_pass"):
        print("CONTROL-FAIL on record — refusing to mine (prereg sec 4).")
        sys.exit(2)
    return cv["producer"]["stamp_utc"]


def mine_point(pt):
    cfg = POINTS[pt]
    rec = {"point": pt, "cap_digits": cfg["cap"], "windows": {}}
    Cs = {}
    for leg in ("lo", "hi"):
        path = f"{t5_dir()}/conn_{pt}_{leg}.json"
        C, meta = parse_c_matrix(path, cfg["cap"])
        Cs[leg] = C
        rec[f"input_{leg}"] = meta
    rec["slice_check"] = {}
    for sector in ("sym", "anti"):
        rec["slice_check"][sector] = bool(
            slice_check(Cs["lo"], sector, cfg["use_im"]))
        for mode in MODES:
            for leg in ("lo", "hi"):
                d = cfg["d"][leg]
                w = mine_window_f1(Cs[leg], sector, mode, d,
                                   cfg["d_capacity"], RUNGS,
                                   use_im=cfg["use_im"],
                                   tag=f"{pt}-{sector}-{mode}-{leg}")
                rec["windows"][f"{sector}_{mode}_{leg}"] = w
                print(" ", w["tag"], "audit", w["audit"]["audit_pass"],
                      "hits", len(w["hits"]), "wall",
                      w["wall_s"], "priced_out", bool(w["priced_out"]),
                      flush=True)
    # G2 lattice-level two-leg identity per (sector, mode)
    for sector in ("sym", "anti"):
        for mode in MODES:
            lo = [list(h["vec"]) for h in
                  rec["windows"][f"{sector}_{mode}_lo"]["hits"]]
            hi = [list(h["vec"]) for h in
                  rec["windows"][f"{sector}_{mode}_hi"]["hits"]]
            rec[f"G2_{sector}_{mode}"] = {
                "hnf_lo": hnf_basis(lo), "hnf_hi": hnf_basis(hi),
                "two_leg_lattice_identical": hnf_basis(lo) == hnf_basis(hi),
                "rank": len(hnf_basis(lo)),
            }
    return rec


def w5f1_window(pt):
    cfg = POINTS[pt]
    d = cfg["d"]["lo"]
    C, meta = parse_c_matrix(f"{t5_dir()}/conn_{pt}_lo.json", cfg["cap"])
    rows, keep, dz, dd, pairs = w5f1_products(C, cfg["use_im"], d)
    n = len(keep)
    rec = {"point": pt, "leg": "lo", "d_leg": d,
           "n_products_doubled_total": 2 * len(pairs),
           "n_after_quotient": n, "n_zero_dropped": len(dz),
           "n_dup_dropped": len(dd),
           "capacity": [capacity_line_f1(n, H, cfg["d_capacity"])
                        for H in (2, 3, 10)]}
    lawful = [c for c in rec["capacity"] if c["lawful"]]
    if not lawful:
        rec["verdict"] = "no lawful rung at this dimension (capacity refused)"
        rec["hits"] = []
        return rec
    top = int(lawful[-1]["H"])
    sub = [r[:120] for r in rows]
    t0 = time.time()
    _c, backend, _w = mine_homogeneous(sub, 120, d, top, bkz_beta=0)
    pilot = time.time() - t0
    proj = pilot * (n / 120.0) ** 4
    rec["pilot"] = {"dim": 120, "wall_s": round(pilot, 2), "backend": backend,
                    "projected_full_s": round(proj, 1), "wall_cap_s": 1800,
                    "parent_integer_pilot_context":
                        "302.09 s (p1) / 558.25 s (m1) at dim 120 (reference run)"}
    if proj > 1800:
        rec["verdict"] = ("CAPACITY-PRICED-OUT: full dim-%d reduction projects "
                          "%.0f s > 1800 s desk wall (pilot %.2f s at dim 120); "
                          "wider is a relay ask" % (n, proj, pilot))
        rec["hits"] = []
        return rec
    t0 = time.time()
    cands, backend, _w = mine_homogeneous(rows, n, d, top, bkz_beta=0)
    hits = []
    for c in cands:
        res = residuals(rows, c, d + 60)
        ok, worst = passes_floor(res, d, margin=30)
        if ok:
            hits.append({"vec_on_kept": c, "kept_index": keep,
                         "height": max(abs(x) for x in c),
                         "worst_resid_log10": worst})
    rec["full"] = {"wall_s": round(time.time() - t0, 1), "backend": backend,
                   "n_candidates": len(cands)}
    rec["hits"] = hits
    rec["verdict"] = ("HITS:%d" % len(hits)) if hits else \
        ("NULL-with-floor at height %d over the quotiented %d-column doubled "
         "basket" % (top, n))
    return rec


def main():
    stamp = check_control()
    for pt, out in (("p1", "MINE_P1_F1.json"), ("m1", "MINE_M1_F1.json")):
        t0 = time.time()
        rec = mine_point(pt)
        rec["control_pass_stamp"] = stamp
        rec["wall_s"] = round(time.time() - t0, 1)
        emit(rec, f"{leg_dir_f1()}/work/{out}", os.path.abspath(__file__))
        for sector in ("sym", "anti"):
            for mode in MODES:
                g2 = rec[f"G2_{sector}_{mode}"]
                print(pt, sector, mode, "rank", g2["rank"],
                      "two_leg_identical", g2["two_leg_lattice_identical"])
    w5 = {"window": "W5F1 blind doubled basket (prereg sec 5)", "points": {}}
    for pt in ("p1", "m1"):
        w5["points"][pt] = w5f1_window(pt)
        print("W5F1", pt, w5["points"][pt]["verdict"])
    emit(w5, f"{leg_dir_f1()}/work/MINE_W5_F1.json", os.path.abspath(__file__))


if __name__ == "__main__":
    main()

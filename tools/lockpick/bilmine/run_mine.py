#!/usr/bin/env python3
# lockpick bilmine member — sealed-window production miner driver (per-point SYM/ANTI windows x two precision legs, height ladder, blind product window W5).
"""run_mine.py — the two-point bilinear production mine (prereg sec 5; runs
ONLY after CONTROL-PASS). Windows W1-W4 (per-point SYM/ANTI, lo/hi legs) + W5
blind-666 honesty window (pilot-priced). Emits work/MINE_P1.json,
work/MINE_M1.json, work/MINE_W5.json under $BILMINE_LEG."""
import json, os, sys, time

import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bilmine_lib import (leg_dir, t5_dir, emit, parse_c_matrix, hnf_basis,
                         capacity_line,
                         mine_homogeneous, residuals, passes_floor, audit_rows)
from run_control import mine_window

RUNGS = [100, 10 ** 4, 10 ** 6, 10 ** 8, 10 ** 10]

POINTS = {
    "p1": {"cap": 543, "d": {"lo": 500, "hi": 535}, "use_im": False},
    "m1": {"cap": 558, "d": {"lo": 515, "hi": 550}, "use_im": True},
}


def check_control():
    cv = json.load(open(f"{leg_dir()}/work/CONTROL_VERDICT.json"))
    if not cv.get("control_pass"):
        print("CONTROL-FAIL on record — refusing to mine (prereg sec 4).")
        sys.exit(2)
    return cv["producer"]["stamp_utc"]


def mine_point(pt):
    cfg = POINTS[pt]
    rec = {"point": pt, "cap_digits": cfg["cap"], "windows": {}}
    for leg in ("lo", "hi"):
        path = f"{t5_dir()}/conn_{pt}_{leg}.json"
        C, meta = parse_c_matrix(path, cfg["cap"])
        rec[f"input_{leg}"] = meta
        d = cfg["d"][leg]
        for sector in ("sym", "anti"):
            w = mine_window(C, sector, d, RUNGS, use_im=cfg["use_im"],
                            tag=f"{pt}-{sector}-{leg}")
            rec["windows"][f"{sector}_{leg}"] = w
    # G2 lattice-level two-leg identity per sector
    for sector in ("sym", "anti"):
        lo = [list(h["vec"]) for h in rec["windows"][f"{sector}_lo"]["hits"]]
        hi = [list(h["vec"]) for h in rec["windows"][f"{sector}_hi"]["hits"]]
        rec[f"G2_{sector}"] = {
            "hnf_lo": hnf_basis(lo), "hnf_hi": hnf_basis(hi),
            "two_leg_lattice_identical": hnf_basis(lo) == hnf_basis(hi),
            "rank": len(hnf_basis(lo)),
        }
    return rec


def w5_products(C, use_im, d):
    """666-pair product vector + degeneracy quotient (zero/dup columns)."""
    ents = [(i, j) for i in range(6) for j in range(6)]
    pairs = [(a, b) for a in range(len(ents)) for b in range(a, len(ents))]
    with mp.workdps(d + 60):
        vals = []
        for (a, b) in pairs:
            (i, j), (k, l) = ents[a], ents[b]
            vals.append(C[i][j] * C[k][l])
        scale = max(max(abs(v.real), abs(v.imag)) for v in vals)
        thr = scale * mp.mpf(10) ** (-(d - 30))
        keep, dropped_zero = [], []
        for t, v in enumerate(vals):
            if max(abs(v.real), abs(v.imag)) <= thr:
                dropped_zero.append(t)
            else:
                keep.append(t)
        # duplicate detection among kept columns (relative at 30d)
        seen, dropped_dup = {}, []
        keep2 = []
        for t in keep:
            key = mp.nstr(vals[t].real / scale, 25) + "|" + mp.nstr(vals[t].imag / scale, 25)
            if key in seen:
                dropped_dup.append((seen[key], t))
            else:
                seen[key] = t
                keep2.append(t)
        rows = [[vals[t].real for t in keep2]]
        if use_im:
            rows.append([vals[t].imag for t in keep2])
    return rows, keep2, dropped_zero, dropped_dup, pairs


def w5_window(pt):
    cfg = POINTS[pt]
    d = cfg["d"]["lo"]
    C, meta = parse_c_matrix(f"{t5_dir()}/conn_{pt}_lo.json", cfg["cap"])
    rows, keep, dz, dd, pairs = w5_products(C, cfg["use_im"], d)
    n = len(keep)
    rec = {"point": pt, "leg": "lo", "d_leg": d, "n_products_total": len(pairs),
           "n_after_quotient": n, "n_zero_dropped": len(dz),
           "n_dup_dropped": len(dd),
           "capacity": [capacity_line(n, H, d) for H in (2, 3, 10)]}
    # pilot pricing (prereg W5): dim-120 pilot, quartic projection, 30-min wall
    sub = [r[:120] for r in rows]
    t0 = time.time()
    _c, backend, wall = mine_homogeneous(sub, 120, d, 3, bkz_beta=0)
    pilot = time.time() - t0
    proj = pilot * (n / 120.0) ** 4
    rec["pilot"] = {"dim": 120, "wall_s": round(pilot, 2), "backend": backend,
                    "projected_full_s": round(proj, 1), "wall_cap_s": 1800}
    if proj > 1800:
        rec["verdict"] = ("CAPACITY-PRICED-OUT: full dim-%d reduction projects "
                          "%.0f s > 1800 s desk wall (pilot %.2f s at dim 120); "
                          "wider is a relay ask" % (n, proj, pilot))
        rec["hits"] = []
        return rec
    t0 = time.time()
    cands, backend, wall = mine_homogeneous(rows, n, d, 3, bkz_beta=0)
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
        "NULL-with-floor at height 3 over the quotiented %d-product basket" % n
    return rec


def main():
    stamp = check_control()
    for pt, out in (("p1", "MINE_P1.json"), ("m1", "MINE_M1.json")):
        t0 = time.time()
        rec = mine_point(pt)
        rec["control_pass_stamp"] = stamp
        rec["wall_s"] = round(time.time() - t0, 1)
        emit(rec, f"{leg_dir()}/work/{out}", os.path.abspath(__file__))
        for sector in ("sym", "anti"):
            g2 = rec[f"G2_{sector}"]
            print(pt, sector, "rank", g2["rank"],
                  "two_leg_identical", g2["two_leg_lattice_identical"])
    w5 = {"window": "W5 blind-666 (prereg sec 5)", "points": {}}
    for pt in ("p1", "m1"):
        w5["points"][pt] = w5_window(pt)
        print("W5", pt, w5["points"][pt]["verdict"])
    emit(w5, f"{leg_dir()}/work/MINE_W5.json", os.path.abspath(__file__))


if __name__ == "__main__":
    main()

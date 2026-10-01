#!/usr/bin/env python3
# lockpick bilmine member — the control battery (C1 planted bilinear mined BLIND; C2 negative; C3 full-shape 6x6 plant) + mine_window, the shared window miner.
"""run_control.py — BILMINE controls (ALL before any production number;
prereg sec 4).

C1 planted bilinear control: a certified period 5-vector (Pi, N supplied via
   BILMINE_C1_PI / BILMINE_C1_N) — Pi products mined BLIND (miner never reads
   N); recovery compared to the primitive integer quadric of N^-1 only after
   the mine emits. Two legs (40/50 dps) x two engines (lattice route +
   mp.pslq route).
C2 negative control: sha512-derived synthetic 5-vector, matched magnitudes —
   must NULL at every rung on both engines.
C3 full-shape 6x6 plant: M_plant = exp(K).R with K^T S_a + S_a K = 0 exactly;
   blind SYM-window production mine must CONTAIN the planted (S_a, S_b) and
   every extracted vector must verify at the floors.

control_pass := C1 and C2 and C3 — script-emitted to work/CONTROL_VERDICT.json.
"""
import hashlib, json, os, sys
from fractions import Fraction
from math import gcd

import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bilmine_lib import (leg_dir, require_env, emit, sha_file,
                         mine_homogeneous, residuals,
                         passes_floor, capacity_line, audit_rows, canonicalize,
                         build_rows, entry_list, mats_to_vec, hnf_basis)

P5 = [(a, b) for a in range(5) for b in range(a, 5)]  # 15 sym pairs on 5


def parse_pi(path):
    import re as rex
    d = json.load(open(path))
    out = []
    for s in d["Pi"]:
        m = rex.match(r"\(?\s*([+-]?[0-9.eE]+)\s*([+-]\s*[0-9.eE]+)j", s)
        re_s, im_s = m.group(1), m.group(2).replace(" ", "")
        with mp.workdps(len(s) + 50):
            out.append(mp.mpc(mp.mpf(re_s), mp.mpf(im_s)))
    return out


def product_rows(vec5):
    """15 sym products of a 5-vector -> 2 real constraint rows (Re, Im)."""
    prods = [vec5[a] * vec5[b] for (a, b) in P5]
    return [[z.real for z in prods], [z.imag for z in prods]], prods


def mine_quadric(vec5, d_leg, rungs, tag):
    """Blind quadric mine on a 5-vector: lattice engine + pslq engine."""
    with mp.workdps(d_leg + 60):
        rows, prods = product_rows(vec5)
    audit = audit_rows(rows, d_leg)
    ladder, found = [], []
    for H in rungs:
        cl = capacity_line(15, H, d_leg)
        ladder.append(cl)
        if not cl["lawful"]:
            cl["refused"] = True
            continue
        cands, backend, wall = mine_homogeneous(rows, 15, d_leg, H)
        hits = []
        for c in cands:
            res = residuals(rows, c, d_leg + 60)
            ok, worst = passes_floor(res, d_leg, margin=10)
            if ok:
                hits.append({"vec": c, "height": max(abs(x) for x in c),
                             "worst_resid_log10": worst})
        cl.update({"backend": backend, "wall_s": round(wall, 3),
                   "n_candidates": len(cands), "n_hits": len(hits)})
        for h in hits:
            if h["vec"] not in [f["vec"] for f in found]:
                found.append(h)
    # pslq engine on the Re row, verified on Im (second engine)
    pslq_hit = None
    with mp.workdps(d_leg):
        r = mp.pslq([x for x in rows[0]], maxcoeff=int(rungs[-1]), maxsteps=100000)
    if r:
        cc = list(canonicalize(r))
        res = residuals(rows, cc, d_leg + 60)
        ok, worst = passes_floor(res, d_leg, margin=10)
        pslq_hit = {"vec": cc, "passes_both_rows": ok, "worst_resid_log10": worst}
    return {"tag": tag, "d_leg": d_leg, "audit": audit, "ladder": ladder,
            "lattice_hits": found, "pslq_hit": pslq_hit}


def known_quadric_from_N():
    """Primitive integer c-vector of Pi^T N^-1 Pi in the P5 product ordering.
    READ ONLY AT COMPARISON TIME (after the blind mine)."""
    N_PATH = require_env("BILMINE_C1_N",
                         "the C1 control's integer certificate JSON")
    N = json.load(open(N_PATH))["N"]
    n = len(N)
    import sympy
    Ninv = sympy.Matrix(N).inv()
    cs = []
    for (a, b) in P5:
        v = Fraction(int(Ninv[a, b].p), int(Ninv[a, b].q))
        cs.append(v if a == b else 2 * v)
    den = 1
    for f in cs:
        den = den // gcd(den, f.denominator) * f.denominator
    iv = [int(f * den) for f in cs]
    return list(canonicalize(iv))


def synthetic_pi(mags):
    """sha512-derived deterministic 5-vector with matched magnitudes/precision."""
    out = []
    for k in range(5):
        digs = ""
        seed = f"BILMINE-NEG-{k+1}"
        h = seed.encode()
        while len(digs) < 140:
            h = hashlib.sha512(h).digest()
            digs += "".join(str(b % 10) for b in h)
        re_s = digs[:62]
        im_s = digs[62:124]
        with mp.workdps(120):
            re = mp.mpf("0." + re_s) * mags[k]
            im = mp.mpf("0." + im_s) * mags[k]
            out.append(mp.mpc(re, im))
    return out


def c3_plant():
    """Committed deterministic 6x6 plant (prereg C3)."""
    S_a = [[2, 1, 0, 0, 0, 1],
           [1, -2, 1, 0, 0, 0],
           [0, 1, 3, 1, 0, 0],
           [0, 0, 1, -1, 1, 0],
           [0, 0, 0, 1, 2, 1],
           [1, 0, 0, 0, 1, -3]]
    A0 = [[0, 1, -1, 0, 2, 0],
          [-1, 0, 1, -2, 0, 1],
          [1, -1, 0, 1, 0, -1],
          [0, 2, -1, 0, 1, 0],
          [-2, 0, 0, -1, 0, 2],
          [0, -1, 1, 0, -2, 0]]
    R = [[1, 0, 1, 0, 0, -1],
         [0, 1, 0, 1, 0, 0],
         [-1, 0, 1, 0, 1, 0],
         [0, 1, 0, 1, 1, 0],
         [0, 0, -1, 0, 1, 1],
         [1, 0, 0, 1, 0, 1]]
    import sympy
    Sa, A0m, Rm = sympy.Matrix(S_a), sympy.Matrix(A0), sympy.Matrix(R)
    assert Sa.det() != 0 and Rm.det() != 0
    assert (A0m.T + A0m).is_zero_matrix
    K = Sa.inv() * A0m
    assert (K.T * Sa + Sa * K).is_zero_matrix
    S_b = (Rm.T * Sa * Rm).tolist()
    d = 600
    with mp.workdps(d + 40):
        Km = mp.matrix(6, 6)
        for i in range(6):
            for j in range(6):
                fr = Fraction(int(K[i, j].p), int(K[i, j].q))
                Km[i, j] = mp.mpf(fr.numerator) / mp.mpf(fr.denominator)
        E = mp.expm(Km)
        Mm = E * mp.matrix([[mp.mpf(x) for x in row] for row in R])
        Mp = [[mp.mpc(Mm[i, j], 0) for j in range(6)] for i in range(6)]
    return S_a, [[int(x) for x in row] for row in S_b], Mp


def mine_window(C, sector, d_leg, rungs, use_im, tag, floor_margin=30):
    """One window: single LLL at the top lawful rung; ladder applied at
    extraction; holdout = near-integer gate on determined S_mum entries."""
    from bilmine_lib import holdout_eval
    fit, hold = entry_list(sector, holdout_index=0)
    rows, tags = build_rows(C, sector, fit, use_im, d_leg + 60)
    n_loc = len(sym_pairs_len(sector))
    n = n_loc + len(fit)
    audit = audit_rows(rows, d_leg)
    quotient = None
    if audit["duplicate_row_pairs"] and not audit["zero_rows"] \
            and not audit["zero_cols"]:
        # G5 re-pose: quotient the declared duplicates (drop the LATER row of
        # each duplicate pair; prereg sec 6 G5 — no new freedom introduced)
        drop = sorted({b for (_a, b) in audit["duplicate_row_pairs"]})
        rows = [r for t, r in enumerate(rows) if t not in drop]
        tags = [g for t, g in enumerate(tags) if t not in drop]
        re_audit = audit_rows(rows, d_leg)
        quotient = {"quotiented": True, "dropped_rows": drop,
                    "dropped_tags": [],
                    "re_audit": re_audit}
        audit = dict(audit)
        audit["audit_pass"] = re_audit["audit_pass"]
        audit["quotient"] = quotient
    ladder = [capacity_line(n, H, d_leg) for H in rungs]
    lawful = [cl for cl in ladder if cl["lawful"]]
    for cl in ladder:
        if not cl["lawful"]:
            cl["refused"] = True
    found = []
    backend = wall = None
    if lawful and audit["audit_pass"]:
        top = int(lawful[-1]["H"])
        cands, backend, wall = mine_homogeneous(rows, n, d_leg, top, bkz_beta=0)
        thr = mp.mpf(10) ** (-(d_leg - floor_margin))
        for c in cands:
            res_f = residuals(rows, c, d_leg + 60)
            ok_f, worst_f = passes_floor(res_f, d_leg, margin=floor_margin)
            if not ok_f:
                continue
            hev = holdout_eval(C, sector, hold, c[:n_loc], d_leg + 60)
            ok_h = all(h["dist_to_int_rel"] <= thr and h["im_rel"] <= thr
                       for h in hev)
            worst_h = max([max(h["dist_to_int_rel"], h["im_rel"]) for h in hev],
                          default=mp.mpf(0))
            if ok_h:
                found.append({
                    "vec": c, "height": max(abs(x) for x in c),
                    "worst_fit_resid_log10": worst_f,
                    "worst_holdout_resid_log10":
                        float(mp.log10(worst_h)) if worst_h > 0 else -9999.0,
                    "holdout_determined_S_mum":
                        {str(h["entry"]): h["nearest_int"] for h in hev}})
        lawful[-1].update({"n_candidates": len(cands)})
    for cl in ladder:
        if cl["lawful"]:
            cl["n_hits_at_or_below"] = sum(1 for f in found
                                           if f["height"] <= cl["H"])
    return {"tag": tag, "sector": sector, "d_leg": d_leg, "audit": audit,
            "n_unknowns": n, "fit_entries": fit, "holdout_entries": hold,
            "backend": backend, "wall_s": (round(wall, 2) if wall else None),
            "ladder": ladder, "hits": found}


def sym_pairs_len(sector):
    from bilmine_lib import sym_pairs, anti_pairs
    return sym_pairs() if sector == "sym" else anti_pairs()


def main():
    out = {"control": "BILMINE C1/C2/C3 (prereg sec 4)"}

    # ---- C1: planted bilinear control (BLIND mine, then compare) ----
    PI_PATH = require_env("BILMINE_C1_PI",
                          "the C1 control's certified period 5-vector JSON")
    Pi = parse_pi(PI_PATH)
    c1 = {"input": {"path": PI_PATH, "sha256": sha_file(PI_PATH)},
          "legs": {}}
    for leg, d in (("lo", 40), ("hi", 50)):
        c1["legs"][leg] = mine_quadric(Pi, d, [10, 100], f"C1-{leg}")
    # canonical agreement across legs + engines
    def lat_vecs(legrec):
        return sorted(tuple(h["vec"]) for h in legrec["lattice_hits"])
    agree_legs = lat_vecs(c1["legs"]["lo"]) == lat_vecs(c1["legs"]["hi"])
    lo_hits = lat_vecs(c1["legs"]["lo"])
    pslq_lo = c1["legs"]["lo"]["pslq_hit"]
    pslq_hi = c1["legs"]["hi"]["pslq_hit"]
    engines_agree = (pslq_lo and pslq_hi and pslq_lo["passes_both_rows"]
                     and pslq_hi["passes_both_rows"]
                     and tuple(pslq_lo["vec"]) in lo_hits
                     and tuple(pslq_hi["vec"]) in lo_hits)
    rank1 = len(lo_hits) == 1
    # comparison to the known lattice — read N only NOW
    target = known_quadric_from_N()
    c1["known_primitive_quadric_from_Ninv"] = target
    recovered = rank1 and agree_legs and bool(engines_agree) and \
        list(lo_hits[0]) == target
    c1.update({"two_leg_agree": agree_legs, "two_engine_agree": bool(engines_agree),
               "lattice_rank": len(lo_hits),
               "recovered_known_lattice_exactly": bool(recovered),
               "pass": bool(recovered)})
    out["C1"] = c1

    # ---- C2: negative control ----
    mags = [max(abs(z), mp.mpf(1) / 4) for z in Pi]
    Neg = synthetic_pi(mags)
    c2 = {"legs": {}}
    for leg, d in (("lo", 40), ("hi", 50)):
        c2["legs"][leg] = mine_quadric(Neg, d, [10, 100], f"C2-{leg}")
    n_hits = sum(len(c2["legs"][l]["lattice_hits"]) for l in c2["legs"])
    pslq_null = all(not c2["legs"][l]["pslq_hit"] or
                    not c2["legs"][l]["pslq_hit"]["passes_both_rows"]
                    for l in c2["legs"])
    c2.update({"n_lattice_hits": n_hits, "pslq_null": bool(pslq_null),
               "pass": n_hits == 0 and bool(pslq_null)})
    out["C2"] = c2

    # ---- C3: full-shape 6x6 plant through the production window ----
    S_a, S_b, Mp = c3_plant()
    from bilmine_lib import sym_pairs
    fit_e, hold_e = entry_list("sym", 0)
    plant_vec = list(canonicalize(
        [S_a[k][l] for (k, l) in sym_pairs()] + [S_b[a][b] for (a, b) in fit_e]))
    c3 = {"planted_pair_canonical": plant_vec, "legs": {}}
    for leg, d in (("lo", 500), ("hi", 535)):
        c3["legs"][leg] = mine_window(Mp, "sym", d, [100, 10**4, 10**6, 10**8, 10**10],
                                      use_im=False, tag=f"C3-{leg}")
    def vecs(rec):
        return [list(h["vec"]) for h in rec["hits"]]
    contains = all(plant_vec in vecs(c3["legs"][l]) for l in c3["legs"])
    # G2 is a LATTICE-level identity (prereg sec 6): canonical HNF of the
    # Z-span must agree — individual LLL representatives may differ.
    hnf_lo = hnf_basis(vecs(c3["legs"]["lo"]))
    hnf_hi = hnf_basis(vecs(c3["legs"]["hi"]))
    legs_agree = hnf_lo == hnf_hi
    c3.update({"contains_plant": bool(contains), "two_leg_agree": bool(legs_agree),
               "lattice_rank_lo": len(hnf_lo), "lattice_rank_hi": len(hnf_hi),
               "hnf_lattice": hnf_lo,
               "n_hits_lo": len(c3["legs"]["lo"]["hits"]),
               "pass": bool(contains and legs_agree)})
    out["C3"] = c3

    out["control_pass"] = bool(c1["pass"] and c2["pass"] and c3["pass"])
    out["verdict"] = ("CONTROL-PASS" if out["control_pass"] else "CONTROL-FAIL")
    emit(out, f"{leg_dir()}/work/CONTROL_VERDICT.json", os.path.abspath(__file__))
    print("CONTROL:", out["verdict"],
          "| C1", c1["pass"], "| C2", c2["pass"], "| C3", c3["pass"])
    return 0 if out["control_pass"] else 2


if __name__ == "__main__":
    sys.exit(main())

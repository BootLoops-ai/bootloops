#!/usr/bin/env python3
# lockpick bilmine member (F1 ring stage) — the ring control battery (C1' planted sqrt15 mined BLIND; C2' negative; C3' sqrt15-scaled full-shape plant; detector smoke).
"""run_control_f1.py — BILMINE-F1 controls (ALL before any production
number; prereg sec 4). Control-upgrade law: the planted control is itself a
sqrt15-COEFFICIENT relation mined BLIND (an integer plant cannot license a
ring instrument).

C1' planted sqrt15 control: Pi' = diag(s,1,s,1,s)*Pi (s = sqrt15, committed
    e = (1,0,1,0,1)); doubled 30-unknown basket mined blind (N read only at
    comparison time); two legs (52/56 dps) x two engines (LLL + mp.pslq).
C2' negative: sha512 "BILMINE-F1-NEG-k" synthetic 5-vector, doubled basket —
    must NULL everywhere.
C3' full-shape plant: C' = exp(K)*R*diag(1,1,s,1,1,1) — sqrt15-scaled pairing
    on the MUM side; blind production GEN SYM mine must CONTAIN the plant.
DETECTOR smoke test (sealed prereg sec 4): the G4a' Q(sqrt15) detector.

control_pass := C1' and C2' and C3' and DET — script-emitted to
work/CONTROL_VERDICT.json.
"""
import hashlib, json, os, sys
from fractions import Fraction
from math import gcd

import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bilmine_f1_lib import (leg_dir_f1, emit, sha_file, mine_homogeneous,
                            residuals, require_env,
                            passes_floor, audit_rows, canonicalize,
                            capacity_line_f1, mine_window_f1, hnf_basis,
                            entry_list, sym_pairs, detect_qsqrt15, s15)
sys.path.insert(0, os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))  # the parent bilmine member dir
from run_control import parse_pi, known_quadric_from_N, c3_plant  # noqa: E402

E_COMMITTED = (1, 0, 1, 0, 1)   # sealed sqrt15 exponent vector (prereg sec 4)
P5 = [(a, b) for a in range(5) for b in range(a, 5)]


def scaled_pi(Pi):
    with mp.workdps(90):
        s = s15()
        return [Pi[k] * (s ** E_COMMITTED[k]) for k in range(5)]


def doubled_product_rows(vec5, d):
    """15 sym products doubled by sqrt15 -> 2 real rows in 30 unknowns
    (block order [a-parts | b-parts], prereg sec 4)."""
    with mp.workdps(d + 60):
        s = s15()
        prods = [vec5[a] * vec5[b] for (a, b) in P5]
        cols = prods + [s * v for v in prods]
        return [[z.real for z in cols], [z.imag for z in cols]]


def mine_quadric_f1(vec5, d_leg, rungs, tag):
    rows = doubled_product_rows(vec5, d_leg)
    audit = audit_rows(rows, d_leg)
    ladder, found = [], []
    for H in rungs:
        cl = capacity_line_f1(30, H, d_leg)
        ladder.append(cl)
        if not cl["lawful"]:
            cl["refused"] = True
            continue
        cands, backend, wall = mine_homogeneous(rows, 30, d_leg, H, bkz_beta=0)
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
    # second engine: pslq on the Re row at the TOP LAWFUL rung, verified on both
    top_lawful = max((cl["H"] for cl in ladder if cl["lawful"]), default=0)
    pslq_hit = None
    if top_lawful:
        with mp.workdps(d_leg):
            r = mp.pslq([x for x in rows[0]], maxcoeff=int(top_lawful),
                        maxsteps=100000)
        if r:
            cc = list(canonicalize(r))
            res = residuals(rows, cc, d_leg + 60)
            ok, worst = passes_floor(res, d_leg, margin=10)
            pslq_hit = {"vec": cc, "passes_both_rows": ok,
                        "worst_resid_log10": worst}
    return {"tag": tag, "d_leg": d_leg, "audit": audit, "ladder": ladder,
            "lattice_hits": found, "pslq_hit": pslq_hit}


def expected_from_N():
    """The N-derived expected canonical Z[sqrt15] vector for Pi' (read N only
    at comparison time). Committed transform: coefficient at slot (p,q) of the
    known integer quadric c becomes c * sqrt15^(2 - e_p - e_q) after clearing
    by 15: esum=0 -> a += 15c; esum=1 -> b += c; esum=2 -> a += c."""
    iq = known_quadric_from_N()
    a = [0] * 15
    b = [0] * 15
    for t, (p, q) in enumerate(P5):
        c = iq[t]
        if c == 0:
            continue
        esum = E_COMMITTED[p] + E_COMMITTED[q]
        if esum == 0:
            a[t] += 15 * c
        elif esum == 1:
            b[t] += c
        else:
            a[t] += c
    return list(canonicalize(a + b))


def committed_expected():
    """The prereg sec-4 committed paragraph, as a vector: relation
    8*sqrt15*v(0,1) + 4*sqrt15*v(2,3) - v(4,4) = 0."""
    a = [0] * 15
    b = [0] * 15
    b[P5.index((0, 1))] = 8
    b[P5.index((2, 3))] = 4
    a[P5.index((4, 4))] = -1
    return list(canonicalize(a + b))


def synthetic_pi_f1(mags):
    out = []
    for k in range(5):
        digs = ""
        seed = f"BILMINE-F1-NEG-{k+1}"
        h = seed.encode()
        while len(digs) < 140:
            h = hashlib.sha512(h).digest()
            digs += "".join(str(x % 10) for x in h)
        re_s = digs[:62]
        im_s = digs[62:124]
        with mp.workdps(120):
            out.append(mp.mpc(mp.mpf("0." + re_s) * mags[k],
                              mp.mpf("0." + im_s) * mags[k]))
    return out


def c3prime_matrix():
    """C' = M * diag(1,1,sqrt15,1,1,1) with the parent's committed plant M."""
    S_a, S_b, Mp = c3_plant()
    with mp.workdps(640):
        s = s15()
        Cp = [[Mp[i][j] * (s if j == 2 else 1) for j in range(6)]
              for i in range(6)]
    return S_a, S_b, Cp


def c3prime_plant_vec(S_a, S_b):
    """Expected plant in GEN SYM coordinates: (S_a, D_r S_b D_r), D_r scaling
    MUM column 2 by sqrt15. Blocks [a_loc | b_loc | X_fit | Y_fit]."""
    fit, _hold = entry_list("sym", 0)
    a_loc = [S_a[k][l] for (k, l) in sym_pairs()]
    b_loc = [0] * len(sym_pairs())
    X, Y = [], []
    for (i, j) in fit:
        n2 = (1 if i == 2 else 0) + (1 if j == 2 else 0)
        v = S_b[i][j]
        if n2 == 0:
            X.append(v); Y.append(0)
        elif n2 == 1:
            X.append(0); Y.append(v)
        else:
            X.append(15 * v); Y.append(0)
    return list(canonicalize(a_loc + b_loc + X + Y))


def detector_smoke():
    tests = []
    with mp.workdps(400):
        big = mp.mpf(10) ** 60
        # (i) rational fractional part
        r1, e1 = detect_qsqrt15(big + mp.mpf(7) / 2000)
        tests.append({"case": "10^60 + 7/2000",
                      "detected": (str(r1[0]), str(r1[1])) if r1 else None,
                      "err_log10": e1,
                      "pass": bool(r1 and r1[0] == Fraction(10**60) +
                                   Fraction(7, 2000) and r1[1] == 0)})
        # (ii) genuine Q(sqrt15) fractional part
        s = s15()
        r2, e2 = detect_qsqrt15(big + (3 + 2 * s) / 7)
        tests.append({"case": "10^60 + (3+2*sqrt15)/7",
                      "detected": (str(r2[0]), str(r2[1])) if r2 else None,
                      "err_log10": e2,
                      "pass": bool(r2 and r2[0] == Fraction(10**60) +
                                   Fraction(3, 7) and r2[1] == Fraction(2, 7))})
        # (iii) pi must refuse
        r3, e3 = detect_qsqrt15(big + mp.pi)
        tests.append({"case": "10^60 + pi", "detected": bool(r3),
                      "err_log10": e3, "pass": r3 is None})
        # (iv) sha-seeded random 40-digit fractional must refuse
        h = hashlib.sha512(b"BILMINE-F1-DET-NEG").digest()
        digs = "".join(str(x % 10) for x in h)[:40]
        r4, e4 = detect_qsqrt15(big + mp.mpf("0." + digs))
        tests.append({"case": "10^60 + sha-frac(40d)", "detected": bool(r4),
                      "err_log10": e4, "pass": r4 is None})
        # (v) near-zero short-circuit
        r5, e5 = detect_qsqrt15(mp.mpf(10) ** -320, near_zero_log10=-300)
        tests.append({"case": "1e-320 (near-zero, thr -300)",
                      "detected": (str(r5[0]), str(r5[1])) if r5 else None,
                      "err_log10": e5,
                      "pass": bool(r5 and r5[0] == 0 and r5[1] == 0)})
    return {"tests": tests, "pass": all(t["pass"] for t in tests)}


def main():
    out = {"control": "BILMINE-F1 C1'/C2'/C3' + detector smoke (prereg sec 4)",
           "committed_exponents": list(E_COMMITTED)}

    # ---- C1': planted sqrt15-coefficient control (BLIND, then compare) ----
    PI_PATH = require_env("BILMINE_C1_PI",
                          "the C1 control's certified period 5-vector JSON")
    Pi = parse_pi(PI_PATH)
    PiP = scaled_pi(Pi)
    c1 = {"input": {"path": PI_PATH, "sha256": sha_file(PI_PATH)},
          "scaling": "Pi' = diag(s,1,s,1,s) * Pi, s = sqrt15 (committed)",
          "legs": {}}
    for leg, d in (("lo", 52), ("hi", 56)):
        c1["legs"][leg] = mine_quadric_f1(PiP, d, [10, 100], f"C1p-{leg}")

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
    # comparison targets — N read only NOW, after the blind mine
    target_N = expected_from_N()
    target_committed = committed_expected()
    commit_match = target_N == target_committed
    c1["expected_from_N"] = target_N
    c1["expected_committed_prereg"] = target_committed
    c1["commit_matches_N_derivation"] = bool(commit_match)
    recovered = (rank1 and agree_legs and bool(engines_agree)
                 and list(lo_hits[0]) == target_N and commit_match)
    c1.update({"two_leg_agree": agree_legs,
               "two_engine_agree": bool(engines_agree),
               "lattice_rank": len(lo_hits),
               "recovered_known_sqrt15_lattice_exactly": bool(recovered),
               "pass": bool(recovered)})
    out["C1p"] = c1

    # ---- C2': negative control on the doubled basket ----
    mags = [max(abs(z), mp.mpf(1) / 4) for z in PiP]
    Neg = synthetic_pi_f1(mags)
    c2 = {"legs": {}}
    for leg, d in (("lo", 52), ("hi", 56)):
        c2["legs"][leg] = mine_quadric_f1(Neg, d, [10, 100], f"C2p-{leg}")
    n_hits = sum(len(c2["legs"][l]["lattice_hits"]) for l in c2["legs"])
    pslq_null = all(not c2["legs"][l]["pslq_hit"] or
                    not c2["legs"][l]["pslq_hit"]["passes_both_rows"]
                    for l in c2["legs"])
    c2.update({"n_lattice_hits": n_hits, "pslq_null": bool(pslq_null),
               "pass": n_hits == 0 and bool(pslq_null)})
    out["C2p"] = c2

    # ---- C3': production-shape plant with sqrt15-scaled pairing ----
    S_a, S_b, Cp = c3prime_matrix()
    plant_vec = c3prime_plant_vec(S_a, S_b)
    c3 = {"planted_pair_canonical": plant_vec,
          "construction": "C' = exp(K)*R*diag(1,1,sqrt15,1,1,1); planted "
                          "(S_loc, S_mum) = (S_a, D_r S_b D_r)",
          "legs": {}}
    for leg, d in (("lo", 500), ("hi", 535)):
        c3["legs"][leg] = mine_window_f1(
            Cp, "sym", "gen", d, 500, [100, 10**4, 10**6, 10**8, 10**10],
            use_im=False, tag=f"C3p-{leg}")

    def vecs(rec):
        return [list(h["vec"]) for h in rec["hits"]]
    # direct residual verify of the plant on each leg's rows (arbiter)
    plant_direct = {}
    for leg, d in (("lo", 500), ("hi", 535)):
        from bilmine_f1_lib import build_rows_f1
        fit, _h = entry_list("sym", 0)
        rows, _t = build_rows_f1(Cp, "sym", fit, False, d + 60, "gen")
        res = residuals(rows, plant_vec, d + 60)
        okp, worstp = passes_floor(res, d, margin=30)
        plant_direct[leg] = {"passes_floor": bool(okp),
                             "worst_resid_log10": worstp}
    hnf_lo = hnf_basis(vecs(c3["legs"]["lo"]))
    hnf_hi = hnf_basis(vecs(c3["legs"]["hi"]))
    legs_agree = hnf_lo == hnf_hi
    contains_verbatim = all(plant_vec in vecs(c3["legs"][l]) for l in c3["legs"])
    contains_lattice = all(
        hnf_basis(vecs(c3["legs"][l])) ==
        hnf_basis(vecs(c3["legs"][l]) + [plant_vec]) and len(vecs(c3["legs"][l])) > 0
        for l in c3["legs"])
    c3.update({"contains_plant_verbatim": bool(contains_verbatim),
               "contains_plant_lattice_level": bool(contains_lattice),
               "plant_direct_residual": plant_direct,
               "two_leg_agree": bool(legs_agree),
               "lattice_rank_lo": len(hnf_lo), "lattice_rank_hi": len(hnf_hi),
               "hnf_lattice": hnf_lo,
               "pass": bool(contains_lattice and legs_agree
                            and all(plant_direct[l]["passes_floor"]
                                    for l in plant_direct))})
    out["C3p"] = c3

    # ---- detector smoke (sealed as part of the control battery) ----
    out["DET"] = detector_smoke()

    out["control_pass"] = bool(out["C1p"]["pass"] and out["C2p"]["pass"]
                               and out["C3p"]["pass"] and out["DET"]["pass"])
    out["verdict"] = ("CONTROL-PASS" if out["control_pass"] else "CONTROL-FAIL")
    emit(out, f"{leg_dir_f1()}/work/CONTROL_VERDICT.json",
         os.path.abspath(__file__))
    print("CONTROL:", out["verdict"], "| C1p", out["C1p"]["pass"],
          "| C2p", out["C2p"]["pass"], "| C3p", out["C3p"]["pass"],
          "| DET", out["DET"]["pass"])
    return 0 if out["control_pass"] else 2


if __name__ == "__main__":
    sys.exit(main())

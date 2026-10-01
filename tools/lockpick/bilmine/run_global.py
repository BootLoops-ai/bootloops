#!/usr/bin/env python3
# lockpick bilmine member — gates G3 (cross-point shared S_mum) + G4 (global invariance through the certified b-link) on the mined lattices.
"""run_global.py — BILMINE gates G3 (cross-point shared S_mum) + G4 (global
monodromy invariance through the stored b-link) on the mined lattices
(prereg sec 6). Emits work/GLOBAL_GATE.json."""
import json, os, sys
from fractions import Fraction
from math import gcd

import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bilmine_lib import (leg_dir, t5_dir, require_env, emit, parse_c_matrix,
                         sym_pairs, anti_pairs,
                         holdout_eval, canonicalize, exact_invariance_check,
                         detect_rational, sha_file)
from run_control import entry_list

CAPS = {"p1": 543, "m1": 558}
DLEG = {"p1": 500, "m1": 515}


def full_smum(pt, sector, vec):
    """HNF basis vector -> full 6x6 integer S_mum (fit part + determined
    holdout ints recomputed for THIS vector) + the S_loc 6x6."""
    loc = sym_pairs() if sector == "sym" else anti_pairs()
    n_loc = len(loc)
    fit, hold = entry_list(sector, 0)
    C, _ = parse_c_matrix(f"{t5_dir()}/conn_{pt}_lo.json", CAPS[pt])
    d = DLEG[pt]
    hev = holdout_eval(C, sector, hold, vec[:n_loc], d + 60)
    thr = mp.mpf(10) ** (-(d - 30))
    ok = all(h["dist_to_int_rel"] <= thr and h["im_rel"] <= thr for h in hev)
    Sl = [[0] * 6 for _ in range(6)]
    Sm = [[0] * 6 for _ in range(6)]
    for t, (k, l) in enumerate(loc):
        Sl[k][l] = int(vec[t])
        Sl[l][k] = int(vec[t]) if sector == "sym" else -int(vec[t])
    for t, (a, b) in enumerate(fit):
        Sm[a][b] = int(vec[n_loc + t])
        Sm[b][a] = int(vec[n_loc + t]) if sector == "sym" else -int(vec[n_loc + t])
    for h in hev:
        a, b = h["entry"]
        Sm[a][b] = h["nearest_int"]
        Sm[b][a] = h["nearest_int"] if sector == "sym" else -h["nearest_int"]
    return Sl, Sm, ok


def parse_A(tier):
    T5 = t5_dir()
    d = json.load(open(f"{T5}/linkb_{tier}.json"))
    A = d["A_cols_cycles_rows_mum"]
    out = [[None] * 6 for _ in range(6)]
    with mp.workdps(460):
        for i in range(6):
            for j in range(6):
                re_s, im_s = A[i][j]
                out[i][j] = mp.mpc(mp.mpf(re_s), mp.mpf(im_s))
    return out, {"dps": d["dps"], "resid": d["A_solve_rel_resid"],
                 "sha256": sha_file(f"{T5}/linkb_{tier}.json")}


def g4_gate(Sm, gens, W_sym_basis, P):
    """Pull S_mum through the b-link; rational-detect; exact invariance; W membership."""
    rec = {}
    fracs = None
    for tier in ("lo", "hi"):
        A, meta = parse_A(tier)
        with mp.workdps(380):
            Sc = [[mp.fsum(A[i][a] * Sm[i][j] * A[j][b]
                           for i in range(6) for j in range(6))
                   for b in range(6)] for a in range(6)]
            im_max = max(abs(Sc[a][b].imag) /
                         max(max(abs(Sc[x][y].real) for x in range(6)
                                 for y in range(6)), mp.mpf(1))
                         for a in range(6) for b in range(6))
            fr = [[None] * 6 for _ in range(6)]
            errs = []
            allok = True
            scale = max(abs(Sc[x][y].real) for x in range(6) for y in range(6))
            for a in range(6):
                for b in range(6):
                    f, e = detect_rational(
                        Sc[a][b].real, den_cap=10 ** 6, tol_log10=-40,
                        near_zero_log10=float(mp.log10(scale)) - 300)
                    fr[a][b] = f
                    errs.append(e)
                    if f is None:
                        allok = False
        rec[tier] = {"A_meta": meta, "rational_all": bool(allok),
                     "max_rel_im": mp.nstr(im_max, 4),
                     "worst_rational_err_log10":
                         (max(e for e in errs if e > -9000) if any(e > -9000 for e in errs) else -9999)}
        if allok and fracs is None:
            fracs = fr
        elif allok and fracs is not None:
            rec["tiers_agree"] = (fr == fracs)
    if fracs is None:
        rec["S_cyc_exact"] = None
        return rec, None
    den = 1
    for row in fracs:
        for f in row:
            den = den // gcd(den, f.denominator) * f.denominator
    Sz = [[int(f * den) for f in row] for row in fracs]
    g = 0
    for row in Sz:
        for x in row:
            g = gcd(g, abs(x))
    if g > 1:
        Sz = [[x // g for x in row] for row in Sz]
    rec["S_cyc_exact_primitive"] = Sz
    rec["denominator_cleared"] = den
    inv = {o: exact_invariance_check(Sz, gens, o) for o in ("GtSG", "GSGt")}
    rec["exact_invariance_all10"] = inv
    # membership in W_sym (dim-1: proportionality to the primitive W generator)
    wgen = W_sym_basis[0] if W_sym_basis else None
    prop = None
    if wgen:
        loc = sym_pairs()
        wmat = [[0] * 6 for _ in range(6)]
        for t, (k, l) in enumerate(loc):
            wmat[k][l] = wgen[t]
            wmat[l][k] = wgen[t]
        num = den2 = None
        okp = True
        for a in range(6):
            for b in range(6):
                if wmat[a][b] == 0 and Sz[a][b] == 0:
                    continue
                if wmat[a][b] == 0 or Sz[a][b] == 0:
                    okp = False
                    break
                f = Fraction(Sz[a][b], wmat[a][b])
                if num is None:
                    num = f
                elif f != num:
                    okp = False
                    break
            if not okp:
                break
        prop = {"proportional_to_W_generator": bool(okp and num is not None),
                "ratio": str(num) if okp and num is not None else None}
        # and to the stored P
        okP = True
        numP = None
        for a in range(6):
            for b in range(6):
                if P[a][b] == 0 and Sz[a][b] == 0:
                    continue
                if P[a][b] == 0 or Sz[a][b] == 0:
                    okP = False
                    break
                f = Fraction(Sz[a][b], P[a][b])
                if numP is None:
                    numP = f
                elif f != numP:
                    okP = False
                    break
            if not okP:
                break
        prop["proportional_to_banked_P"] = bool(okP and numP is not None)
        prop["ratio_to_P"] = str(numP) if okP and numP is not None else None
    rec["W_membership"] = prop
    passed = (rec["lo"]["rational_all"] and rec["hi"]["rational_all"]
              and rec.get("tiers_agree", False)
              and (inv["GtSG"] or inv["GSGt"])
              and prop and prop["proportional_to_W_generator"])
    return rec, bool(passed)


def main():
    LEG = leg_dir()
    p1 = json.load(open(f"{LEG}/work/MINE_P1.json"))
    m1 = json.load(open(f"{LEG}/work/MINE_M1.json"))
    ex = json.load(open(f"{LEG}/work/EXACT_INVARIANTS.json"))
    pf = json.load(open(require_env(
        "BILMINE_PFRAME", "the exact generator-frame JSON")))
    gens = [pf["G_integer"][k] for k in pf["frame_specification"]["GENS_order"]]
    P = pf["invariant_form"]["P"]
    W_sym = ex["W_invariant_spaces"]["GtSG"]["sym_basis_pairsorder"]

    out = {"gate": "BILMINE G3+G4 (prereg sec 6)", "sectors": {}}
    for sector in ("sym", "anti"):
        sec = {"g2_p1": p1[f"G2_{sector}"], "g2_m1": m1[f"G2_{sector}"]}
        r1, r2 = sec["g2_p1"]["rank"], sec["g2_m1"]["rank"]
        sec["ranks"] = [r1, r2]
        if r1 == 0 or r2 == 0:
            sec["G3_shared_S_mum"] = None
            sec["chain"] = "no candidate chain (rank 0 at a point)"
            out["sectors"][sector] = sec
            continue
        chains = []
        for v1 in sec["g2_p1"]["hnf_lo"]:
            Sl1, Sm1_, ok1 = full_smum("p1", sector, [int(x) for x in v1])
            for v2 in sec["g2_m1"]["hnf_lo"]:
                Sl2, Sm2_, ok2 = full_smum("m1", sector, [int(x) for x in v2])
                c1 = canonicalize([Sm1_[a][b] for (a, b) in sym_pairs()]) \
                    if sector == "sym" else \
                    canonicalize([Sm1_[a][b] for (a, b) in anti_pairs()])
                c2 = canonicalize([Sm2_[a][b] for (a, b) in sym_pairs()]) \
                    if sector == "sym" else \
                    canonicalize([Sm2_[a][b] for (a, b) in anti_pairs()])
                match = (list(c1) == list(c2)) and any(c1)
                ch = {"S_loc_p1": Sl1, "S_mum_p1": Sm1_, "holdout_ok_p1": ok1,
                      "S_loc_m1": Sl2, "S_mum_m1": Sm2_, "holdout_ok_m1": ok2,
                      "G3_S_mum_match": bool(match)}
                if match and ok1 and ok2:
                    g4, g4pass = g4_gate(Sm1_, gens, W_sym, P)
                    ch["G4"] = g4
                    ch["G4_pass"] = g4pass
                chains.append(ch)
        sec["chains"] = chains
        sec["n_full_pass"] = sum(1 for c in chains
                                 if c.get("G3_S_mum_match") and c.get("G4_pass"))
        out["sectors"][sector] = sec
    emit(out, f"{LEG}/work/GLOBAL_GATE.json", os.path.abspath(__file__))
    for sector, sec in out["sectors"].items():
        print(sector, "ranks", sec.get("ranks"),
              "full-pass chains", sec.get("n_full_pass"))


if __name__ == "__main__":
    main()

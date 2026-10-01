#!/usr/bin/env python3
# lockpick bilmine member (F1 ring stage) — gates G3 (cross-point shared S_mum over Z[sqrt15]) + G4 (global invariance with the sealed ring detector).
"""run_global_f1.py — BILMINE-F1 gates G3 (cross-point shared S_mum over
Z[sqrt15]) + G4 (global monodromy invariance through the stored b-link with
the sealed Q(sqrt15) detector) on the mined lattices (prereg sec 6).
Emits work/GLOBAL_GATE.json."""
import json, os, sys
from fractions import Fraction
from math import gcd

import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bilmine_f1_lib import (leg_dir_f1, t5_dir, require_env, MODES, emit,
                            parse_c_matrix, sha_file,
                            sloc_parts, smum_fit_parts, holdout_eval_f1,
                            full_smum_xy, smum_canonical, detect_qsqrt15,
                            exact_invariance_check, s15, entry_list)

CAPS = {"p1": 543, "m1": 558}
DLEG = {"p1": 500, "m1": 515}


def top_lawful(win):
    return max((cl["H"] for cl in win["ladder"] if cl["lawful"]), default=0)


def candidate_smum(pt, sector, mode, vec, win, C_cache):
    """HNF basis vector -> (ok, (X, Y) full 6x6 S_mum pair, S_loc parts)."""
    if pt not in C_cache:
        C_cache[pt], _ = parse_c_matrix(f"{t5_dir()}/conn_{pt}_lo.json",
                                        CAPS[pt])
    C = C_cache[pt]
    d = DLEG[pt]
    A_loc, B_loc = sloc_parts(vec, sector, mode)
    Xf, Yf = smum_fit_parts(vec, sector, mode)
    _fit, hold = entry_list(sector, 0)
    ok, worst, hev = holdout_eval_f1(C, sector, hold, A_loc, B_loc, d + 60,
                                     mode, top_lawful(win), d)
    hold_ab = {str(h["entry"]): [h["alpha"], h["beta"]] for h in hev}
    if not ok:
        return False, None, None, (A_loc, B_loc)
    X, Y = full_smum_xy(sector, Xf, Yf, hold_ab)
    return True, X, Y, (A_loc, B_loc)


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


def prop_to_P(M, P):
    """Exact proportionality of a rational matrix to integer P (zero counts)."""
    if all(x == 0 for row in M for x in row):
        return True, Fraction(0)
    ratio = None
    for a in range(6):
        for b in range(6):
            if P[a][b] == 0 and M[a][b] == 0:
                continue
            if P[a][b] == 0 or M[a][b] == 0:
                return False, None
            f = Fraction(M[a][b], P[a][b]) if isinstance(M[a][b], int) \
                else M[a][b] / P[a][b]
            if ratio is None:
                ratio = f
            elif f != ratio:
                return False, None
    return True, ratio


def g4_gate(X, Y, gens, P, sector):
    """Pull S_mum = X + sqrt15*Y through the b-link; detect over Q(sqrt15)
    (sealed absolute detector); exact componentwise invariance; membership."""
    rec = {}
    detected = None
    for tier in ("lo", "hi"):
        A, meta = parse_A(tier)
        with mp.workdps(380):
            s = s15()
            Sm = [[X[i][j] + s * Y[i][j] for j in range(6)] for i in range(6)]
            Sc = [[mp.fsum(A[i][a] * Sm[i][j] * A[j][b]
                           for i in range(6) for j in range(6))
                   for b in range(6)] for a in range(6)]
            scale = max(max(abs(Sc[a][b].real) for a in range(6)
                            for b in range(6)), mp.mpf(1))
            im_max = max(abs(Sc[a][b].imag) / scale
                         for a in range(6) for b in range(6))
            det = [[None] * 6 for _ in range(6)]
            errs = []
            allok = True
            nz = float(mp.log10(scale)) - 300
            for a in range(6):
                for b in range(6):
                    rs, e = detect_qsqrt15(Sc[a][b].real, den_cap=10 ** 6,
                                           tol_log10=-40, near_zero_log10=nz)
                    det[a][b] = rs
                    errs.append(e)
                    if rs is None:
                        allok = False
        rec[tier] = {"A_meta": meta, "detected_all": bool(allok),
                     "max_rel_im": mp.nstr(im_max, 4),
                     "worst_err_log10":
                         (max(e for e in errs if e > -9000)
                          if any(e > -9000 for e in errs) else -9999)}
        if allok and detected is None:
            detected = det
        elif allok and detected is not None:
            rec["tiers_agree"] = (det == detected)
    if detected is None:
        rec["S_cyc_exact"] = None
        return rec, None
    Xc = [[detected[a][b][0] for b in range(6)] for a in range(6)]
    Yc = [[detected[a][b][1] for b in range(6)] for a in range(6)]

    def clear(M):
        den = 1
        for row in M:
            for f in row:
                den = den // gcd(den, f.denominator) * f.denominator
        Z = [[int(f * den) for f in row] for row in M]
        g = 0
        for row in Z:
            for x in row:
                g = gcd(g, abs(x))
        if g > 1:
            Z = [[x // g for x in row] for row in Z]
        return Z, den

    Xz, denX = clear(Xc)
    Yz, denY = clear(Yc)
    rec["X_cyc_primitive"] = Xz
    rec["Y_cyc_primitive"] = Yz
    inv = {}
    for o in ("GtSG", "GSGt"):
        okX = exact_invariance_check(Xz, gens, o) if any(any(r) for r in Xz) \
            else True
        okY = exact_invariance_check(Yz, gens, o) if any(any(r) for r in Yz) \
            else True
        inv[o] = bool(okX and okY)
    rec["exact_invariance_all10_componentwise"] = inv
    okXp, ratX = prop_to_P(Xc, P)
    okYp, ratY = prop_to_P(Yc, P)
    if sector == "anti":
        member = all(x == 0 for row in Xc for x in row) and \
            all(y == 0 for row in Yc for y in row)
        rec["W_membership"] = {"anti_sector_requires_zero": bool(member)}
    else:
        member = okXp and okYp
        rec["W_membership"] = {
            "X_prop_P": bool(okXp), "ratio_X": str(ratX) if okXp else None,
            "Y_prop_P": bool(okYp), "ratio_Y": str(ratY) if okYp else None,
            "S_cyc_eq_u_plus_w_sqrt15_times_P": bool(member)}
    passed = (rec["lo"]["detected_all"] and rec["hi"]["detected_all"]
              and rec.get("tiers_agree", False)
              and (inv["GtSG"] or inv["GSGt"]) and bool(member))
    return rec, bool(passed)


def main():
    LEG = leg_dir_f1()
    p1 = json.load(open(f"{LEG}/work/MINE_P1_F1.json"))
    m1 = json.load(open(f"{LEG}/work/MINE_M1_F1.json"))
    pf = json.load(open(require_env(
        "BILMINE_PFRAME", "the exact generator-frame JSON")))
    gens = [pf["G_integer"][k] for k in pf["frame_specification"]["GENS_order"]]
    P = pf["invariant_form"]["P"]
    C_cache = {}

    out = {"gate": "BILMINE-F1 G3+G4 (prereg sec 6)", "sectors": {}}
    for sector in ("sym", "anti"):
        for mode in MODES:
            key = f"{sector}_{mode}"
            sec = {"g2_p1": p1[f"G2_{sector}_{mode}"],
                   "g2_m1": m1[f"G2_{sector}_{mode}"]}
            r1, r2 = sec["g2_p1"]["rank"], sec["g2_m1"]["rank"]
            sec["ranks"] = [r1, r2]
            if r1 == 0 or r2 == 0:
                sec["G3_shared_S_mum"] = None
                sec["chain"] = "no candidate chain (rank 0 at a point)"
                sec["n_full_pass"] = 0
                out["sectors"][key] = sec
                continue
            chains = []
            for v1 in sec["g2_p1"]["hnf_lo"]:
                ok1, X1, Y1, sl1 = candidate_smum(
                    "p1", sector, mode, [int(x) for x in v1],
                    p1["windows"][f"{sector}_{mode}_lo"], C_cache)
                for v2 in sec["g2_m1"]["hnf_lo"]:
                    ok2, X2, Y2, sl2 = candidate_smum(
                        "m1", sector, mode, [int(x) for x in v2],
                        m1["windows"][f"{sector}_{mode}_lo"], C_cache)
                    ch = {"holdout_ok_p1": bool(ok1), "holdout_ok_m1": bool(ok2)}
                    if ok1 and ok2:
                        c1 = smum_canonical(sector, X1, Y1)
                        c2 = smum_canonical(sector, X2, Y2)
                        match = (c1 is not None and list(c1) == list(c2))
                        ch.update({"S_mum_p1_XY": [X1, Y1],
                                   "S_mum_m1_XY": [X2, Y2],
                                   "S_loc_p1_AB": sl1, "S_loc_m1_AB": sl2,
                                   "G3_S_mum_match": bool(match)})
                        if match:
                            g4, g4pass = g4_gate(X1, Y1, gens, P, sector)
                            ch["G4"] = g4
                            ch["G4_pass"] = g4pass
                    chains.append(ch)
            sec["chains"] = chains
            sec["n_full_pass"] = sum(1 for c in chains
                                     if c.get("G3_S_mum_match")
                                     and c.get("G4_pass"))
            out["sectors"][key] = sec
    emit(out, f"{LEG}/work/GLOBAL_GATE.json", os.path.abspath(__file__))
    for key, sec in out["sectors"].items():
        print(key, "ranks", sec.get("ranks"),
              "full-pass chains", sec.get("n_full_pass"))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# lockpick bilmine member (F1 ring stage) — ring-window layer: Z[sqrt15] doubled-basket windows GEN/O-A/O-B over the parent miner library.
"""bilmine_f1_lib.py — BILMINE-F1: Q(sqrt15)-coefficient bilinear-lattice
windows on two-point certified period data.

REUSES the parent BILMINE instrument (../bilmine_lib.py): parse/trim laws,
the homogeneous LLL miner, audits, floors, HNF lattice
comparison, canonicalization, exact invariance. This module adds ONLY the
ring-window layer sealed by the F1 prereg ceremony:
  - doubled / one-sided row builders (modes GEN, OA, OB; sec 2);
  - the near-Z[sqrt15] holdout gate (sec 6 G1h; held-out coords never free);
  - the Q(sqrt15) fractional-part-ABSOLUTE detector (sec 6 G4a');
  - W5F1 doubled blind basket (sec 5);
  - F1 producer emit.
Backend law: LLL only (bkz_beta=0, fplll-cli "-m proved" chain); per-
reduction wall 3600 s via MPLLL_FPLLL_TIMEOUT (sealed sec 3).
"""
import hashlib, json, os, subprocess, sys
from fractions import Fraction

import mpmath as mp

os.environ.setdefault("MPLLL_FPLLL_TIMEOUT", "3600")  # sealed per-reduction wall: this leg's explicit opt-in (the tool default is unbounded)

sys.path.insert(0, os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))  # the parent bilmine member dir
from bilmine_lib import (  # noqa: E402  (parent instrument, reused verbatim)
    parse_c_matrix, sym_pairs, anti_pairs, entry_list, mine_homogeneous,
    residuals, passes_floor, audit_rows, hnf_basis, canonicalize,
    exact_invariance_check, sha_file, require_env, t5_dir)

MODES = ("gen", "oa", "ob")


def leg_dir_f1():
    return require_env("BILMINE_F1_LEG", "the F1 leg working directory — "
                       "receipts land in $BILMINE_F1_LEG/work")


def parent_dir():
    return require_env("BILMINE_PARENT", "the closed parent BILMINE leg "
                       "directory (read-only inputs of record)")


def emit(obj, out_path, script_path):
    stamp = subprocess.check_output(["date", "-u"]).decode().strip()
    obj["producer"] = {
        "leg": "BILMINE-F1",
        "script": os.path.basename(script_path),
        "script_sha256": sha_file(script_path),
        "stamp_utc": stamp,
    }
    json.dump(obj, open(out_path, "w"), indent=1, default=str)
    lint = os.environ.get("BILMINE_PRODUCER_LINT")
    if not lint:
        return None
    r = subprocess.run([sys.executable, lint, "--check", out_path],
                       capture_output=True, text=True)
    line = (r.stdout or r.stderr).strip()
    print(line)
    if r.returncode != 0:
        raise RuntimeError(f"producer_lint failed on {out_path}: {line}")
    return line


def s15():
    """sqrt15 at the CURRENT working precision (sealed: never cached)."""
    return mp.sqrt(15)


# ---------- window construction (prereg sec 2) ----------

def n_unknowns(sector, mode):
    n_loc = len(sym_pairs()) if sector == "sym" else len(anti_pairs())
    fit, _ = entry_list(sector, 0)
    n_fit = len(fit)
    return (2 * (n_loc + n_fit)) if mode == "gen" else (n_loc + n_fit)


def build_rows_f1(C, sector, entries, use_im, dps, mode):
    """Constraint rows for C^T S_loc C = S_mum with Z[sqrt15] coefficients.

    Unknown block order (sealed):
      gen: [S_loc a | S_loc b | S_mum a (fit) | S_mum b (fit)]
      oa:  [S_loc a | S_mum b (fit)]      (S_mum = sqrt15 * integer)
      ob:  [S_loc b | S_mum a (fit)]      (S_loc = sqrt15 * integer)
    Held-out S_mum entries are never unknowns (G1h gates them post-hoc).
    """
    loc = sym_pairs() if sector == "sym" else anti_pairs()
    rows, tags = [], []
    with mp.workdps(dps):
        s = s15()
        zero = mp.mpc(0, 0)
        for (i, j) in entries:
            locv = []
            for (k, l) in loc:
                if sector == "sym":
                    v = C[k][i] * C[l][j] + (C[l][i] * C[k][j] if k != l else zero)
                else:
                    v = C[k][i] * C[l][j] - C[l][i] * C[k][j]
                locv.append(v)
            mum_a = [mp.mpc(-1, 0) if (a, b) == (i, j) else zero
                     for (a, b) in entries]
            mum_b = [mp.mpc(-s, 0) if (a, b) == (i, j) else zero
                     for (a, b) in entries]
            if mode == "gen":
                coeffs = locv + [s * v for v in locv] + mum_a + mum_b
            elif mode == "oa":
                coeffs = locv + mum_b
            elif mode == "ob":
                coeffs = [s * v for v in locv] + mum_a
            else:
                raise ValueError(mode)
            rows.append([c.real for c in coeffs]); tags.append((i, j, "re"))
            if use_im:
                rows.append([c.imag for c in coeffs]); tags.append((i, j, "im"))
    return rows, tags


def sloc_parts(c, sector, mode):
    """Candidate vector -> (A_loc, B_loc) integer coefficient lists (loc order),
    S_loc = A + sqrt15*B."""
    n_loc = len(sym_pairs()) if sector == "sym" else len(anti_pairs())
    if mode == "gen":
        return list(c[:n_loc]), list(c[n_loc:2 * n_loc])
    if mode == "oa":
        return list(c[:n_loc]), [0] * n_loc
    if mode == "ob":
        return [0] * n_loc, list(c[:n_loc])
    raise ValueError(mode)


def smum_fit_parts(c, sector, mode):
    """Candidate vector -> (X_fit, Y_fit) integer lists over fit entries,
    S_mum(fit) = X + sqrt15*Y."""
    n_loc = len(sym_pairs()) if sector == "sym" else len(anti_pairs())
    fit, _ = entry_list(sector, 0)
    n_fit = len(fit)
    if mode == "gen":
        return (list(c[2 * n_loc:2 * n_loc + n_fit]),
                list(c[2 * n_loc + n_fit:2 * n_loc + 2 * n_fit]))
    if mode == "oa":
        return [0] * n_fit, list(c[n_loc:])
    if mode == "ob":
        return list(c[n_loc:]), [0] * n_fit
    raise ValueError(mode)


# ---------- holdout gate G1h (prereg sec 6; near-RING, never free) ----------

def holdout_eval_f1(C, sector, hold_entries, A_loc, B_loc, dps, mode,
                    H_top, d_leg, floor_margin=30):
    """For each held-out entry: v = (C^T (A + sqrt15 B) C)[i,j] at dps; gate
    v against the mode's allowed ring shape at the sealed floor.
    Returns (all_ok, worst_rel, [{entry, alpha, beta, ok, ...}])."""
    loc = sym_pairs() if sector == "sym" else anti_pairs()
    out, allok = [], True
    worst = mp.mpf(0)
    with mp.workdps(dps):
        s = s15()
        thr_rel = mp.mpf(10) ** (-(d_leg - floor_margin))
        for (i, j) in hold_entries:
            acc = mp.mpc(0, 0)
            scale = mp.mpf(0)
            for t, (k, l) in enumerate(loc):
                if sector == "sym":
                    co = C[k][i] * C[l][j] + (C[l][i] * C[k][j] if k != l else 0)
                else:
                    co = C[k][i] * C[l][j] - C[l][i] * C[k][j]
                acc += co * (A_loc[t] + s * B_loc[t])
                scale = max(scale, abs(co), abs(s * co))
            scale = scale or mp.mpf(1)
            thr = scale * thr_rel
            rec = {"entry": (i, j)}
            vr, vi = acc.real, acc.imag
            ok = abs(vi) <= thr
            alpha = beta = None
            if mode == "ob":
                alpha = int(mp.nint(vr))
                ok = ok and abs(vr - alpha) <= thr and abs(alpha) <= H_top
                beta = 0
            elif mode == "oa":
                beta = int(mp.nint(vr / s))
                ok = ok and abs(vr - beta * s) <= thr and abs(beta) <= H_top
                alpha = 0
            else:  # gen: 3-unknown mini-reduction on [v, 1, sqrt15]
                if abs(vr) > H_top * (1 + s) + 1:
                    ok = False
                else:
                    cands, _bk, _w = mine_homogeneous(
                        [[vr, mp.mpf(1), s]], 3, d_leg, 30 * H_top, bkz_beta=0)
                    best = None
                    for cv in cands:
                        if abs(cv[0]) != 1:
                            continue
                        al = -cv[1] * cv[0]
                        be = -cv[2] * cv[0]
                        if abs(al) > H_top or abs(be) > H_top:
                            continue
                        err = abs(vr - al - be * s)
                        if err <= thr and (best is None or err < best[2]):
                            best = (al, be, err)
                    if best is None:
                        ok = False
                    else:
                        alpha, beta = best[0], best[1]
            rec.update({"alpha": alpha, "beta": beta, "ok": bool(ok),
                        "im_rel": float(mp.log10(abs(vi) / scale))
                        if vi != 0 else -9999.0})
            if ok and alpha is not None:
                err = abs(vr - alpha - beta * s) / scale
                worst = max(worst, err, abs(vi) / scale)
            allok = allok and ok
            out.append(rec)
    return allok, (float(mp.log10(worst)) if worst > 0 else -9999.0), out


# ---------- capacity (sealed: gated on the point's d_lo for BOTH legs) ----------

def capacity_line_f1(n, H, d_capacity):
    need = float((n + 1) * mp.log10(mp.mpf(H)) + 20)
    return {"H": H, "need_digits": round(need, 1), "d_capacity": d_capacity,
            "lawful": need <= d_capacity}


# ---------- the F1 window miner (prereg secs 5-6) ----------

def mine_window_f1(C, sector, mode, d_leg, d_capacity, rungs, use_im, tag,
                   floor_margin=30):
    """One (point, sector, mode) window, one leg: G5 audit (+quotient re-pose),
    capacity-lawful single LLL at top lawful rung, extraction with floors +
    the near-ring holdout gate. Returns the window record."""
    fit, hold = entry_list(sector, holdout_index=0)
    n = n_unknowns(sector, mode)
    rows, tags = build_rows_f1(C, sector, fit, use_im, d_leg + 60, mode)
    audit = audit_rows(rows, d_leg)
    quotient = None
    if audit["duplicate_row_pairs"] and not audit["zero_rows"] \
            and not audit["zero_cols"]:
        drop = sorted({b for (_a, b) in audit["duplicate_row_pairs"]})
        dropped_tags = [list(tags[t]) for t in drop]
        rows = [r for t, r in enumerate(rows) if t not in drop]
        tags = [g for t, g in enumerate(tags) if t not in drop]
        re_audit = audit_rows(rows, d_leg)
        quotient = {"quotiented": True, "dropped_rows": drop,
                    "dropped_tags": dropped_tags,
                    "re_audit": re_audit}
        audit = dict(audit)
        audit["audit_pass"] = re_audit["audit_pass"]
        audit["quotient"] = quotient
    ladder = [capacity_line_f1(n, H, d_capacity) for H in rungs]
    lawful = [cl for cl in ladder if cl["lawful"]]
    for cl in ladder:
        if not cl["lawful"]:
            cl["refused"] = True
    found = []
    backend = wall = None
    priced_out = None
    if lawful and audit["audit_pass"]:
        top = int(lawful[-1]["H"])
        try:
            cands, backend, wall = mine_homogeneous(rows, n, d_leg, top,
                                                    bkz_beta=0)
        except subprocess.TimeoutExpired:
            priced_out = ("CAPACITY-PRICED-OUT: reduction exceeded the sealed "
                          "3600 s per-reduction wall (measured timeout); "
                          "wider is a relay ask")
            cands = []
        for c in cands:
            res_f = residuals(rows, c, d_leg + 60)
            ok_f, worst_f = passes_floor(res_f, d_leg, margin=floor_margin)
            if not ok_f:
                continue
            A_loc, B_loc = sloc_parts(c, sector, mode)
            ok_h, worst_h, hev = holdout_eval_f1(
                C, sector, hold, A_loc, B_loc, d_leg + 60, mode, top, d_leg,
                floor_margin)
            if not ok_h:
                continue
            fence = None
            if mode == "gen":
                Xf, Yf = smum_fit_parts(c, sector, mode)
                if not any(B_loc) and not any(Yf):
                    fence = "pure-integer content (inside parent NULL scope)"
                elif not any(A_loc) and not any(Xf):
                    fence = ("pure sqrt15-multiple content (inside parent NULL "
                             "scope after dividing by sqrt15)")
            found.append({
                "vec": c, "height": max(abs(x) for x in c),
                "worst_fit_resid_log10": worst_f,
                "worst_holdout_resid_log10": worst_h,
                "holdout_alpha_beta": {str(h["entry"]): [h["alpha"], h["beta"]]
                                       for h in hev},
                "consistency_fence": fence})
        if lawful and not priced_out:
            lawful[-1]["n_candidates"] = len(cands)
    for cl in ladder:
        if cl["lawful"]:
            cl["n_hits_at_or_below"] = sum(1 for f in found
                                           if f["height"] <= cl["H"])
    return {"tag": tag, "sector": sector, "mode": mode, "d_leg": d_leg,
            "d_capacity": d_capacity, "audit": audit, "n_unknowns": n,
            "fit_entries": fit, "holdout_entries": hold, "backend": backend,
            "wall_s": (round(wall, 2) if wall else None),
            "priced_out": priced_out, "ladder": ladder, "hits": found}


# ---------- W5F1 doubled blind basket (prereg sec 5) ----------

def w5f1_products(C, use_im, d):
    """Doubled 666-pair basket [v, sqrt15*v] interleaved + degeneracy quotient."""
    ents = [(i, j) for i in range(6) for j in range(6)]
    pairs = [(a, b) for a in range(len(ents)) for b in range(a, len(ents))]
    with mp.workdps(d + 60):
        s = s15()
        vals = []
        for (a, b) in pairs:
            (i, j), (k, l) = ents[a], ents[b]
            v = C[i][j] * C[k][l]
            vals.append(v)
            vals.append(s * v)
        scale = max(max(abs(v.real), abs(v.imag)) for v in vals)
        thr = scale * mp.mpf(10) ** (-(d - 30))
        keep, dropped_zero = [], []
        for t, v in enumerate(vals):
            if max(abs(v.real), abs(v.imag)) <= thr:
                dropped_zero.append(t)
            else:
                keep.append(t)
        seen, dropped_dup, keep2 = {}, [], []
        for t in keep:
            key = (mp.nstr(vals[t].real / scale, 25) + "|" +
                   mp.nstr(vals[t].imag / scale, 25))
            if key in seen:
                dropped_dup.append((seen[key], t))
            else:
                seen[key] = t
                keep2.append(t)
        rows = [[vals[t].real for t in keep2]]
        if use_im:
            rows.append([vals[t].imag for t in keep2])
    return rows, keep2, dropped_zero, dropped_dup, pairs


# ---------- Q(sqrt15) detector, fractional-part ABSOLUTE (G4a') ----------

def detect_qsqrt15(x, den_cap=10 ** 6, tol_log10=-40, near_zero_log10=None):
    """x (mpf) -> ((r, s) Fractions with r + s*sqrt15 ~ x, or None, err_log10).

    Sealed criterion (prereg sec 6 G4a'): near-zero short-circuit; split
    x = n + f (n = floor(x), f in [0,1)); mp.pslq([f, 1, sqrt15],
    maxcoeff 1e7); accept iff c0 != 0, |c0| <= den_cap, |c2| <= den_cap and
    |f - fhat| <= 10^tol ABSOLUTE (fhat = -(c1 + c2*sqrt15)/c0). A relative
    criterion is banned (the parent's caught G4 defect)."""
    if near_zero_log10 is not None and (
            x == 0 or abs(x) <= mp.mpf(10) ** near_zero_log10):
        return (Fraction(0), Fraction(0)), -9999.0
    with mp.workdps(mp.mp.dps + 20):
        n = int(mp.floor(x))
        f = x - n
    with mp.workdps(220):
        s = s15()
        try:
            rel = mp.pslq([mp.mpf(f), mp.mpf(1), s], maxcoeff=10 ** 7,
                          maxsteps=10 ** 6)
        except Exception:
            rel = None
    if not rel or rel[0] == 0 or abs(rel[0]) > den_cap or abs(rel[2]) > den_cap:
        return None, 0.0
    c0, c1, c2 = int(rel[0]), int(rel[1]), int(rel[2])
    with mp.workdps(320):
        s = s15()
        fhat = -(mp.mpf(c1) + mp.mpf(c2) * s) / c0
        err = abs(mp.mpf(f) - fhat)
        ok = err <= mp.mpf(10) ** tol_log10
    if not ok:
        return None, (float(mp.log10(err)) if err > 0 else -9999.0)
    r = Fraction(n) + Fraction(-c1, c0)
    t = Fraction(-c2, c0)
    return (r, t), (float(mp.log10(err)) if err > 0 else -9999.0)


# ---------- misc ----------

def pairs_of(sector):
    return sym_pairs() if sector == "sym" else anti_pairs()


def full_smum_xy(sector, Xf, Yf, hold_ab):
    """(X_fit, Y_fit) lists + holdout {entry: (alpha,beta)} -> 6x6 integer
    matrix pair (X, Y) with S_mum = X + sqrt15*Y."""
    fit, hold = entry_list(sector, 0)
    sgn = 1 if sector == "sym" else -1
    X = [[0] * 6 for _ in range(6)]
    Y = [[0] * 6 for _ in range(6)]
    for t, (a, b) in enumerate(fit):
        X[a][b] = int(Xf[t]); X[b][a] = sgn * int(Xf[t])
        Y[a][b] = int(Yf[t]); Y[b][a] = sgn * int(Yf[t])
    for (a, b) in hold:
        al, be = hold_ab[str((a, b))]
        X[a][b] = int(al); X[b][a] = sgn * int(al)
        Y[a][b] = int(be); Y[b][a] = sgn * int(be)
    return X, Y


def smum_canonical(sector, X, Y):
    """Canonical doubled vector of the full S_mum pair over sector pairs."""
    pr = pairs_of(sector)
    return canonicalize([X[a][b] for (a, b) in pr] + [Y[a][b] for (a, b) in pr])

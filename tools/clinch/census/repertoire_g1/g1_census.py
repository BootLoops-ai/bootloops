"""Gate-1 — per-cell certified census (enumeration + B&B driver).

Support universe per cell: all S with |S| <= min(m, n) in itertools
.combinations order, k ascending (the canonical order of the disposition
array); |S| > m covered by the EXCL layer (see g1_cert docstring).
"""
import itertools
import time

import numpy as np
from flint import arb, arb_mat, ctx

import g1_model as M
from g1_cert import (CellCert, ball, amat, _lo, _hi, PREC, BOX_CAP, MINREL,
                     D_V, D_INV, D_BB_EMPTY, D_ROOTS, D_UNCERT)


def conflict_pairs(cc, cell, live, phimin):
    """CERT-INV: ordered pairs (i,l), both live: certified M_il > 0.

    v3 ENGINE AMENDMENT (INV-R fold): the v2 box-form margin is kept
    verbatim and OR-ed with the h-free ratio-form bound INV-R:
      ratio form fires iff b_lj*d_i - d_l*b_ij > 0 certified in arb for
      EVERY niche j with b_ij > 0 (equivalently M_il = r_il*d_i - d_l > 0
      with r_il = certified min_j b_lj/b_ij; division-free, exact-dyadic
      inputs, no float pre-rounding).
    SOUNDNESS (same certificate class, GATE1_RECEIPT Sec. 4 INV semantics):
    at any equilibrium x* whose support contains i, x*_i G_i(x*) = 0 forces
    the resident identity sum_j b_ij phi_j(y*) = d_i; with b >= 0 and
    phi >= 0 the certified rowwise bound b_lj >= r_il b_ij (all j) gives
    G_l(x*) = sum_j b_lj phi_j(y*) - d_l >= r_il*d_i - d_l > 0
    (block-triangular transversal law, cc.stability, already of record)
    => no STABLE state on any support containing i and missing l. Uses NO
    property of phi beyond nonnegativity: valid at every h and K.
    Box path note: q_j = b_lj - b_ij exact dyadic => sign
    decidable in float; the float subtraction pre-rounds before arb
    (an epsilon-thin wart that the independent re-derivation gate covers;
    the ratio path has no such wart)."""
    b, d = cell["b"], cell["d"]
    req = {i: [] for i in live}
    margins = []
    for i in live:
        di = arb(float(d[i]))
        for l in live:
            if l == i:
                continue
            # --- v2 box form, verbatim ---
            s = arb(0)
            for j in range(cc.m):
                q = b[l][j] - b[i][j]
                if q > 0:
                    s += arb(float(q)) * phimin[j]
                elif q < 0:
                    s += arb(float(q))
            s -= arb(float(d[l])) - di
            fired, form, mlo = bool(s > 0), "box", _lo(s)
            # --- ratio form (INV-R, h-free) ---
            if not fired:
                dl = arb(float(d[l]))
                ok, mmin, seen = True, None, False
                for j in range(cc.m):
                    bij = float(b[i][j])
                    if bij > 0.0:
                        seen = True
                        t = arb(float(b[l][j])) * di - dl * arb(bij)
                        if not (t > 0):
                            ok = False
                            break
                        tl = _lo(t)
                        mmin = tl if mmin is None else min(mmin, tl)
                if ok and seen:
                    fired, form, mlo = True, "ratio", mmin
            if fired:
                req[i].append(l)
                margins.append({"i": int(i), "l": int(l),
                                "margin_lo": mlo, "form": form})
    return req, margins


def farkas_infeasible(cc, cell, S):
    """Certified infeasibility of B_S p = d_S, p in [0,1]^m, via a left-null
    witness y of B_S (float svd; verified in arb): for feasible p,
    y.d = (y.B) p, so |y.d| <= sum_j ub|(y.B)_j|; a certified strict
    violation => infeasible. Covers exactly-singular augmented dets
    (e.g. sparse cells with an all-zero niche column)."""
    m = cc.m
    Bf = cell["b"][list(S)]
    df = cell["d"][list(S)]
    U, sv, Vt = np.linalg.svd(Bf)
    ncols = U.shape[1]
    for idx in range(ncols - 1, -1, -1):
        if idx < len(sv) and sv[idx] > 1e-9 * sv[0]:
            break
        y = U[:, idx]
        yd = arb(0)
        for i, si in enumerate(S):
            yd += arb(float(y[i])) * cc.dv[si]
        tot = arb(0)
        for j in range(m):
            s = arb(0)
            for i, si in enumerate(S):
                s += arb(float(y[i])) * cc.B[si, j]
            tot += abs(s.mid()) + s.rad()
        if abs(yd.mid()) - yd.rad() > tot:
            return True
    return False


def excl_layer(cc, cell, live):
    """CERT-EXCL: |S|=m+1 over live clones, certified det[B_S|d_S] != 0;
    Farkas fallback on undischarged dets."""
    m = cc.m
    if len(live) <= m:
        return {"needed": False, "n_checked": 0, "n_open": 0, "open": []}
    n_checked, n_open, n_farkas, opens = 0, 0, 0, []
    min_absdet = None
    for S in itertools.combinations(live, m + 1):
        rows = [[cc.B[i, j] for j in range(m)] + [cc.dv[i]] for i in S]
        det = arb_mat(rows).det()
        n_checked += 1
        if not (det != 0):
            if farkas_infeasible(cc, cell, S):
                n_farkas += 1
            else:
                n_open += 1
                opens.append(list(map(int, S)))
        else:
            a = abs(_lo(abs(det)))
            min_absdet = a if min_absdet is None else min(min_absdet, a)
    return {"needed": True, "n_checked": n_checked, "n_open": n_open,
            "n_farkas_closed": n_farkas, "open": opens[:20],
            "min_absdet_lo": min_absdet}


def census_support(cc, cell, S, xmax, deadline):
    """Certified-complete root census for one support. Returns
    (disposition, roots, stats). roots: list of dicts."""
    Sl = list(S)
    k = len(Sl)
    Bs, Cs, ds = cc.sub(Sl)
    xm = [xmax[i] for i in Sl]
    key = f"fr|{cell['cell_id']}|{','.join(map(str, Sl))}"
    verified = []   # (boxlo, boxhi, encl_lo, encl_hi)

    # k = m: exact linear decision (unique p* -> unique y* -> unique x*)
    if k == cc.m:
        lst, lout = cc.lin_route(Sl, Bs, Cs, ds, xm)
        if lst == "empty":
            return D_BB_EMPTY, [], {"boxes": 0, "cert": "LIN"}
        if lst == "root":
            elo = [_lo(v) for v in lout]
            ehi = [_hi(v) for v in lout]
            return D_ROOTS, [{"S": [int(i) for i in Sl], "encl_lo": elo,
                              "encl_hi": ehi, "positive": True,
                              "cert": "LIN"}], {"boxes": 0, "cert": "LIN"}
        # 'fall' -> full B&B below

    def try_verify(lo, hi):
        st, out = cc.krawczyk(cell, Sl, Bs, Cs, ds, lo, hi)
        if st != "contained":
            return False
        blo, bhi = list(lo), list(hi)
        elo, ehi = out
        for _ in range(3):
            st2, out2 = cc.krawczyk(cell, Sl, Bs, Cs, ds, elo, ehi)
            if st2 == "contained":
                elo, ehi = out2
            else:
                break
        for (_, _, plo, phi_) in verified:
            if all(ehi[i] >= plo[i] and elo[i] <= phi_[i] for i in range(k)):
                return True  # same root (enclosures intersect): keep first
        verified.append((blo, bhi, elo, ehi))
        return True

    for r in M.float_roots(cell, Sl, xm, key):
        for frac in (0.08, 0.02, 0.004):
            lo = [max(0.0, r[i] - frac * xm[i]) for i in range(k)]
            hi = [min(xm[i], r[i] + frac * xm[i]) for i in range(k)]
            if try_verify(lo, hi):
                break

    stack = [([0.0] * k, list(xm))]
    nboxes = 0
    while stack:
        if time.monotonic() > deadline:
            return D_UNCERT, [], {"reason": "CELL_WALL_CAP", "boxes": nboxes}
        nboxes += 1
        if nboxes > BOX_CAP:
            return D_UNCERT, [], {"reason": "BOX_CAP", "boxes": nboxes}
        lo, hi = stack.pop()
        absorbed = False
        for (blo, bhi, _, _) in verified:
            if all(lo[i] >= blo[i] - 1e-300 and hi[i] <= bhi[i] + 1e-300
                   for i in range(k)):
                absorbed = True
                break
        if absorbed:
            continue
        if cc.corner_excludes(Bs, Cs, ds, lo, hi):
            continue
        st, out = cc.krawczyk(cell, Sl, Bs, Cs, ds, lo, hi)
        if st == "empty":
            continue
        if st == "contained":
            elo, ehi = out
            dup = False
            for (_, _, plo, phi_) in verified:
                if all(ehi[i] >= plo[i] and elo[i] <= phi_[i]
                       for i in range(k)):
                    dup = True
                    break
            if not dup:
                verified.append((lo, hi, elo, ehi))
            continue
        if out is not None:
            lo, hi = out  # Krawczyk-tightened
        widths = [(hi[i] - lo[i]) / max(xm[i], 1e-300) for i in range(k)]
        wi = int(np.argmax(widths))
        if widths[wi] < MINREL:
            return D_UNCERT, [], {"reason": "MINREL_FEATURE", "boxes": nboxes,
                                  "at": lo}
        mid = 0.5 * (lo[wi] + hi[wi])
        hi_a = list(hi)
        hi_a[wi] = mid
        lo_b = list(lo)
        lo_b[wi] = mid
        stack.append((list(lo), hi_a))
        stack.append((lo_b, list(hi)))
    roots = []
    for (blo, bhi, elo, ehi) in verified:
        pos = all(elo[i] > 0.0 for i in range(k))
        roots.append({"S": [int(i) for i in Sl], "encl_lo": elo,
                      "encl_hi": ehi, "positive": pos})
    if not roots:
        return D_BB_EMPTY, [], {"boxes": nboxes}
    return D_ROOTS, roots, {"boxes": nboxes}


def enumerate_prefilter(cell, live, dead, req):
    """Yield (k, combos_array, disp_array) with V/INV pre-dispositions."""
    n, m = cell["n"], cell["m"]
    deadv = np.zeros(n, bool)
    deadv[dead] = True
    reqmask = np.zeros(n, np.int64)
    for i, ls in req.items():
        mm = 0
        for l in ls:
            mm |= (1 << l)
        reqmask[i] = mm
    bit = np.array([1 << i for i in range(n)], np.int64)
    out = []
    for k in range(1, min(m, n) + 1):
        combos = np.array(list(itertools.combinations(range(n), k)),
                          np.int64).reshape(-1, k)
        disp = np.zeros(len(combos), np.uint8)
        isdead = deadv[combos].any(axis=1)
        disp[isdead] = D_V
        smask = np.bitwise_or.reduce(bit[combos], axis=1)
        needs = np.bitwise_or.reduce(reqmask[combos], axis=1)
        viol = (needs & ~smask) != 0
        disp[(~isdead) & viol] = D_INV
        out.append((k, combos, disp))
    return out


def run_cell(cell, wall_cap_s=3600.0, log=None):
    t0 = time.monotonic()
    deadline = t0 + wall_cap_s
    ctx.prec = PREC
    cc = CellCert(cell)
    tm = {}
    g, dead, viable, vunk = cc.viability()
    live = sorted(viable + vunk)
    tm["viability"] = time.monotonic() - t0
    # certified equilibrium bounds + global phi floor
    xmax = {}
    bound_fail = []
    for i in live:
        fx = M.f_xmax(cell, i)
        xb = cc.xmax(i, fx)
        if xb is None:
            bound_fail.append(i)
        else:
            # 1.5x slack keeps roots off the search-box boundary (k=1 roots
            # ARE the bound's defining root); any inflation of a certified
            # upper bound stays a certified upper bound (h_i monotone).
            xmax[i] = 1.5 * xb + 1e-12
    tm["xmax"] = time.monotonic() - t0 - tm["viability"]
    if bound_fail:
        return ({"cell_id": cell["cell_id"], "status": "UNCERTIFIED",
                 "reason": f"no equilibrium bound for clones {bound_fail}"},
                np.zeros(0, np.uint8))
    phimin = []
    for j in range(cc.m):
        ym = arb(0)
        for i in live:
            ym += cc.C[i, j] * arb(xmax[i])
        phimin.append(cc.phi(ball(0, _hi(ym))))  # ball covers [phi(ymax), 1]
    # NOTE: phi decreasing => lower bound of phi over feasible y is phi(ymax);
    # using ball(0,ymax) through phi gives a rigorous enclosure of the range.
    req, inv_margins = conflict_pairs(cc, cell, live, phimin)
    tm["conflicts"] = time.monotonic() - t0 - sum(tm.values())
    excl = excl_layer(cc, cell, live)
    tm["excl"] = time.monotonic() - t0 - sum(tm.values())

    levels = enumerate_prefilter(cell, live, dead, req)
    equilibria = []
    counts = {"V": 0, "INV": 0, "BB_EMPTY": 0, "ROOTS": 0, "UNCERT": 0}
    uncert_reasons = []
    boxes_total = 0
    disp_arrays = []
    capped = False
    for (k, combos, disp) in levels:
        hard = np.nonzero(disp == 0)[0]
        for idx in hard:
            if capped or time.monotonic() > deadline:
                capped = True
                disp[idx] = D_UNCERT
                continue
            S = tuple(int(v) for v in combos[idx])
            dcode, roots, stats = census_support(cc, cell, S, xmax, deadline)
            boxes_total += stats.get("boxes", 0)
            disp[idx] = dcode
            if dcode == D_UNCERT:
                uncert_reasons.append({"S": list(S), **{k2: v for k2, v
                                       in stats.items() if k2 != "boxes"}})
            elif dcode == D_ROOTS:
                others = [l for l in live if l not in S]
                for r in roots:
                    if not r["positive"]:
                        uncert_reasons.append({"S": list(S),
                                               "reason": "BOUNDARY_ROOT"})
                        disp[idx] = D_UNCERT
                        continue
                    verdict, inv_signs, internal = cc.stability(
                        cell, S, r["encl_lo"], r["encl_hi"], others)
                    equilibria.append({
                        "S": r["S"], "encl_lo": r["encl_lo"],
                        "encl_hi": r["encl_hi"], "stability": verdict,
                        "internal": internal,
                        "n_inv_pos": sum(1 for v in inv_signs.values()
                                         if v > 0),
                        "n_inv_unk": sum(1 for v in inv_signs.values()
                                         if v == 0)})
        for code, name in ((D_V, "V"), (D_INV, "INV"),
                           (D_BB_EMPTY, "BB_EMPTY"), (D_ROOTS, "ROOTS"),
                           (D_UNCERT, "UNCERT")):
            counts[name] += int((disp == code).sum())
        disp_arrays.append((k, disp))
        if log:
            log(f"  k={k}: {len(combos)} supports, hard={len(hard)}, "
                f"t={time.monotonic()-t0:.1f}s")
    # empty support (x=0): J = diag(g) exactly
    if not vunk:
        empty_st = "STABLE" if all(g[i] < 0 for i in range(cc.n)) else \
                   "UNSTABLE"
    else:
        empty_st = "UNSTABLE" if any(g[i] > 0 for i in range(cc.n)) else \
                   "UNKNOWN"
    equilibria.append({"S": [], "encl_lo": [], "encl_hi": [],
                       "stability": empty_st, "internal": empty_st,
                       "n_inv_pos": 0, "n_inv_unk": 0,
                       "tag": "extinction_state"})
    n_stable = sum(1 for e in equilibria if e["stability"] == "STABLE")
    n_unstable = sum(1 for e in equilibria if e["stability"] == "UNSTABLE")
    n_unk = sum(1 for e in equilibria if e["stability"] == "UNKNOWN")
    stable_supports = [e["S"] for e in equilibria if e["stability"] == "STABLE"]
    complete = (counts["UNCERT"] == 0 and excl.get("n_open", 0) == 0
                and not capped)
    wall = time.monotonic() - t0
    disp_concat = np.concatenate([d for (_, d) in disp_arrays]) \
        if disp_arrays else np.zeros(0, np.uint8)
    return {"cell_id": cell["cell_id"], "status": "OK",
            "axes": {a: cell[a] for a in
                     ("m", "h", "w", "sp", "a", "K", "loss", "div",
                      "seed_tag")},
            "n_dead": len(dead), "n_viable": len(viable), "n_vunk": len(vunk),
            "dispositions": counts, "excl": excl,
            "n_inv_pairs": len(inv_margins),
            "n_supports_enumerated": int(sum(len(c) for (_, c, _) in levels)),
            "n_equilibria": len(equilibria), "n_stable": n_stable,
            "n_unstable": n_unstable, "n_unknown_sign": n_unk,
            "stable_supports": stable_supports,
            "extinction_stable": empty_st == "STABLE",
            "completeness_closed": bool(complete), "capped": bool(capped),
            "uncert_reasons": uncert_reasons[:50],
            "equilibria": equilibria, "inv_margins": inv_margins,
            "boxes_total": int(boxes_total), "wall_s": round(wall, 2),
            "timings": {k2: round(v, 2) for k2, v in tm.items()},
            "prec_bits": PREC}, disp_concat

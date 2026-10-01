#!/usr/bin/env python3
"""mplll.py — multi-point LLL analytic-regression fitter (v2 core).
PRIMARY qrbkz: E=[I;F]→LQ→col-equilibrate L→[I|10^d L]→BKZ-β (kills κ-plateau).
BASELINE rawlll: B=[W·I|10^d·F].  rawbkz = raw lattice + BKZ-β (use when
K+1 > N_pts so qrbkz's LQ is inapplicable AND dim>~200 where plain LLL loses
planted relations (K=236 control)).
GRADED/CVP/BATCH: mplll_graded / mplll_batch.
Gate = held-out ≥30d, two-precision (d,d−20), capacity, +/− controls.
CLI: mplll.py --basis F.json --target I.json --dps D [--method qrbkz|rawlll|rawbkz|cvp]"""
import math, os, re, sys
import mpmath as mp
from .mplll_lattice import reduce_rows, row_norm2, gcd_list   # g3 package-relative
from .pslq_gate import canonicalize                           # (was flat sys.path.insert)
GATE_D, MARGIN = 30, 15
_CRX = re.compile(r'\(?\s*([+\-]?[\d.eE+\-]+)\s*([+\-]\s*[\d.eE+\-]+)j\s*\)?$')

def _mpc(x):
    if isinstance(x, (mp.mpc, mp.mpf, complex, float, int)): return mp.mpc(x)
    s = str(x); m = _CRX.match(s)
    return mp.mpc(m[1], m[2].replace(' ', '')) if m else mp.mpc(mp.mpf(s), 0)


def _cols(I, F, fit, d):
    """Return (col_fn, N_eff): real or [re,im]-split columns per fit index."""
    thr = mp.mpf(10) ** (-d)
    cx = any(abs(_mpc(v[j]).imag) > thr * max(abs(_mpc(v[j]).real), mp.mpf(1))
             for v in [I] + F for j in fit)
    if cx:
        return (lambda v, j: (lambda z: (z.real, z.imag))(_mpc(v[j]))), 2 * len(fit)
    return (lambda v, j: (_mpc(v[j]).real,)), len(fit)


def _extract(R, K, hmax, c0slot=0):
    """→ (relation, basis_rels); basis_rels = short c0==0 rows with negligible
    right-block ⇒ integer relation among F only (redundant basis)."""
    cand, brels = [], []
    Rmax = max((abs(x) for r in R for x in r[K + 1:]), default=1) or 1
    for r in R:
        blk = r[:K + 1]
        g = abs(gcd_list([x for x in blk if x]) or 1)   # rawlll: identity block
        blk = [x // g for x in blk]                     # carries W ⇒ de-scale
        h = max(abs(x) for x in blk); rmax = max((abs(x) for x in r[K + 1:]), default=0)
        if blk[c0slot] != 0 and h <= hmax:
            cand.append((row_norm2(r), blk))
        elif blk[c0slot] == 0 and 0 < h <= hmax and rmax <= hmax and rmax * hmax < Rmax:
            brels.append(canonicalize(blk[1:K + 1]))
    if not cand: return None, brels
    c = list(min(cand)[1]); g = gcd_list(c) or 1; c = [x // g for x in c]
    return (tuple(-x for x in c) if c[0] < 0 else tuple(c)), brels


def _rawlll_lattice(I, F, fit, d, W):
    K = len(F)
    with mp.workdps(d + 35):    # NEVER read ambient dps: API calls at dps=15
        S = mp.mpf(10) ** d     # (an ambient dps=15 silently produced junk relations)
        col, _ = _cols(I, F, fit, d)
        return ([[(W if j == i else 0) for j in range(K + 1)]
                 + [int(mp.nint(x * S)) for j in fit for x in col(v, j)]
                 for i, v in enumerate([I] + F)], None)


def _qrbkz_lattice(I, F, fit, d, guard=30):
    """LQ on rows of E, column-equilibrate L, embed [I|10^d·L].  Near-zero
    L-columns (< 10^{-d}·max) are DROPPED (null-space signal, not constraints).
    Returns (rows, span_residual_d)."""
    K = len(F)
    with mp.workdps(d + guard):     # _cols parses strings: keep off ambient dps
        col, Ne = _cols(I, F, fit, d)
        E = mp.matrix(K + 1, Ne)
        for i, v in enumerate([I] + F):
            cs = [x for j in fit for x in col(v, j)]
            for jj, x in enumerate(cs): E[i, jj] = x
        Qp, Rp = mp.qr(E.T)                       # E^T = Q'R'  ⇒  E = R'^T Q'^T
        L = Rp.T                                  # (K+1)×(K+1) lower-triangular
        nrm = [mp.sqrt(mp.fsum(L[i, j] ** 2 for i in range(K + 1)))
               for j in range(K + 1)]
        span_res = float(-mp.log10(min(nrm) / max(nrm))) if min(nrm) > 0 else float(d + guard)
        thr = max(nrm) * mp.mpf(10) ** (-d)
        keep = [j for j in range(K + 1) if nrm[j] > thr]
        S = mp.mpf(10) ** d
        return ([[(1 if j == i else 0) for j in range(K + 1)]
                 + [int(mp.nint(L[i, j] / nrm[j] * S)) for j in keep]
                 for i in range(K + 1)], span_res)


def heldout_digits(rel, I, F, idx, d):
    with mp.workdps(d + 20):
        c0, ds = mp.mpf(rel[0]), []
        for j in idx:
            pred = -mp.fsum(rel[i + 1] * _mpc(F[i][j]) for i in range(len(F))) / c0
            e = abs(pred - _mpc(I[j])) / max(abs(_mpc(I[j])), mp.mpf(1))
            ds.append(float(-mp.log10(e)) if e > 0 else float(d))
    return min(ds), ds


def capacity(rel, K, N_fit, d):
    h = max(abs(x) for x in rel); need = (K + 1) * math.log10(max(h, 2)) + MARGIN
    return {"height": h, "need_d": round(need, 1), "have_d": N_fit * d, "ok": need < N_fit * d}


def _fit(method, I, F, names, fit, ho, d, W, hmax, two_prec, backend, beta):
    K = len(F)
    build = (_qrbkz_lattice if method == "qrbkz" else
             lambda I, F, fit, d: _rawlll_lattice(I, F, fit, d, W))
    use_beta = beta if method in ("qrbkz", "rawbkz") else 0
    out = {"method": method, "K": K, "N_fit": len(fit), "N_ho": len(ho), "dps": d}
    rows, span_res = build(I, F, fit, d)
    R, out["backend"], out["lll_wall_s"] = reduce_rows(
        rows, backend, bkz_beta=use_beta)
    if span_res is not None: out["span_residual_d"] = round(span_res, 2)
    rel, brels = _extract(R, K, hmax)
    out["relation"] = rel
    if brels:
        out["basis_relations"] = [list(b) for b in sorted(set(brels))]
        sys.stderr.write(f"[mplll] WARN: {len(brels)} basis-internal relation(s) "
                         f"detected (redundant basis functions)\n")
    if rel is None:
        out["status"] = "NULL"; return out
    out["capacity"] = capacity(rel, K, len(fit), d)
    out["heldout_min_d"], out["heldout_per_pt"] = heldout_digits(rel, I, F, ho, d)
    out["insample_min_d"], _ = heldout_digits(rel, I, F, fit, d)
    if two_prec:
        d2 = max(10, d - 20)
        rows2, _ = build(I, F, fit, d2)
        R2, _, out["lll_wall_s_leg2"] = reduce_rows(
            rows2, backend, bkz_beta=use_beta)
        rel2, _ = _extract(R2, K, hmax)
        out["two_prec_stable"] = rel2 is not None and canonicalize(rel) == canonicalize(rel2)
    out["c0"], out["coeffs"] = rel[0], {names[i]: -rel[i + 1] for i in range(K)}
    out["status"] = ("HIT" if out["heldout_min_d"] >= GATE_D
                     and out.get("two_prec_stable", True) else "PARTIAL")
    return out


def mplll_fit(I, F, names, fit, ho, d, W=10 ** 8, hmax=10 ** 30,
              two_prec=True, backend=None, method="qrbkz", beta=None):
    beta = beta or min(len(F) + 1, 40)
    return _fit(method, I, F, names, fit, ho, d, W, hmax, two_prec, backend, beta)


def _rank_exact(rows):
    """Rank over Q of integer row vectors: fraction-free Gaussian elimination on
    Python ints (rows gcd-reduced as they go), never a floating rank."""
    M = [[int(x) for x in r] for r in rows if any(r)]
    if not M: return 0
    rank, ncols = 0, len(M[0])
    for c in range(ncols):
        piv = next((i for i in range(rank, len(M)) if M[i][c]), None)
        if piv is None: continue
        M[rank], M[piv] = M[piv], M[rank]
        p = M[rank][c]
        for i in range(rank + 1, len(M)):
            if M[i][c]:
                q = M[i][c]
                M[i] = [p * a - q * b for a, b in zip(M[i], M[rank])]
                g = gcd_list(M[i]) or 1
                if g > 1: M[i] = [x // g for x in M[i]]
        rank += 1
        if rank == len(M): break
    return rank


def house_verdict(got, want, relations):
    """The synthetic house control's verdict WITH THE SPAN TEST.  got / want are
    canonical relation vectors (c0, c_1..c_K) for the same target; relations =
    the member set's basis-internal relations (integer vectors over the K
    members vanishing at the fit points: the fit's own `basis_relations`).  On
    a linearly DEPENDENT member set the coefficient vector is unique only
    modulo those relations, so the lattice may return a SHORTER representative
    of the SAME function — got != want is then not a failure.  The difference
    want[0]*got - got[0]*want (c0 slot 0; = got - want when the c0 agree) is
    tested for membership in the rational span of the relations by the exact
    rank test rank(R) == rank(R + [difference]).  reason: 'exact' (got == want),
    'equivalent_representative' (in the span: ok), 'not_in_span' (ok False),
    'null' (no relation found: ok False)."""
    if got is None:
        return {"ok": False, "reason": "null"}
    got, want = tuple(int(x) for x in got), tuple(int(x) for x in want)
    if got == want:
        return {"ok": True, "reason": "exact"}
    diff = [want[0] * g - got[0] * w for g, w in zip(got, want)][1:]
    R = [[int(x) for x in r] for r in (relations or [])]
    r0, r1 = _rank_exact(R), _rank_exact(R + [diff])
    in_span = r1 == r0
    return {"ok": in_span,
            "reason": "equivalent_representative" if in_span else "not_in_span",
            "difference": diff, "relations_n": len(R), "relations_rank": r0,
            "rank_with_difference": r1, "height_got": max(abs(x) for x in got),
            "height_want": max(abs(x) for x in want)}


def controls(I, F, names, fit, ho, d, W, backend=None, method="qrbkz"):
    with mp.workdps(d + 30):
        Isyn = [mp.nstr((3 * _mpc(F[0][j]) - 7 * _mpc(F[1][j])) / 5, d + 20)
                for j in range(len(F[0]))]
        Irand = [mp.nstr(mp.rand() * 10 - 5, d + 20) for _ in I]
    want = canonicalize([5, -3, 7] + [0] * (len(F) - 2))
    r = mplll_fit(Isyn, F, names, fit, ho, d, W, two_prec=False, backend=backend, method=method)
    rn = mplll_fit(Irand, F, names, fit, ho, d, W, two_prec=False, backend=backend, method=method)
    got = canonicalize(r["relation"]) if r["relation"] else None
    rels = sorted(set(tuple(int(x) for x in b)                # the member set's basis-internal
                      for b in r.get("basis_relations", []) + rn.get("basis_relations", [])))
    syn = {"control": "synthetic (3f0-7f1)/5", "want": want, "got": got,
           "relations_n": len(rels),
           "relations_source": "basis_relations of the synthetic and negative fits "
                               "(the lattice's c0 == 0 short rows; same backend, dps, points)"}
    syn.update(house_verdict(got, want, rels))           # ok + reason (+ the span-test fields)
    det = [syn,
           {"control": "negative (random target)", "got_status": rn["status"],
            "cap": rn.get("capacity"),
            "ok": rn["relation"] is None or rn.get("heldout_min_d", 0) < 3}]
    return {"ok": all(x["ok"] for x in det), "detail": det}


if __name__ == "__main__": from .mplll_cli import main; main()  # noqa (via shim / python -m lockpick.mplll)

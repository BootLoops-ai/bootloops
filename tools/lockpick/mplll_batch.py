#!/usr/bin/env python3
"""mplll_batch.py — v2 helpers: multi-target LQ cache, CVP/Babai, d_min sweep."""
import os, sys
import mpmath as mp
from .mplll_lattice import reduce_rows, reduce_rows_batch, row_norm2  # g3 package-relative
from . import mplll
from .pslq_gate import canonicalize


def _basis_LQ(F, fit, d, guard=30):
    """LQ of the K×N basis block once. Returns (L_F K×K, Q_F K×N, col_norms)."""
    K = len(F)
    with mp.workdps(d + guard):     # _cols parses strings: keep off ambient dps
        col, Ne = mplll._cols(F[0], F[1:], fit, d)   # detect complex on F only
        E = mp.matrix(K, Ne)
        for i, v in enumerate(F):
            cs = [x for j in fit for x in col(v, j)]
            for jj, x in enumerate(cs): E[i, jj] = x
        Qp, Rp = mp.qr(E.T)              # F = L_F · Q_F,  L_F=Rp.T[:,:K], Q_F=Qp.T[:K,:]
        L_F = mp.matrix([[Rp[j, i] for j in range(K)] for i in range(K)])
        Q_F = mp.matrix([[Qp[jj, k] for jj in range(Ne)] for k in range(K)])
        cn = [mp.sqrt(mp.fsum(L_F[i, j] ** 2 for i in range(K))) for j in range(K)]
    return L_F, Q_F, cn, col, Ne


def _target_row(I, Q_F, col, fit, Ne, d, guard=30):
    """Project target onto Q_F rowspace: returns (p[K], iperp_norm)."""
    K = Q_F.rows
    with mp.workdps(d + guard):
        Iv = mp.matrix(1, Ne)
        cs = [x for j in fit for x in col(I, j)]
        for jj, x in enumerate(cs): Iv[0, jj] = x
        p = [mp.fsum(Iv[0, jj] * Q_F[k, jj] for jj in range(Ne)) for k in range(K)]
        r = Iv - mp.matrix([p]) * Q_F
        iperp = mp.sqrt(mp.fsum(r[0, jj] ** 2 for jj in range(Ne)))
    return p, iperp


def _embed(L_F, cn, p, iperp, d, K):
    """Assemble (K+1)×(K+1) L=[[iperp,p];[0,L_F]], equilibrate, embed [I|10^d L]."""
    with mp.workdps(d + 30):
        nrm = [max(iperp, mp.mpf(10) ** (-d - 20))] + \
              [mp.sqrt(p[j] ** 2 + cn[j] ** 2) for j in range(K)]
        span_res = float(-mp.log10(min(nrm) / max(nrm)))
        thr = max(nrm) * mp.mpf(10) ** (-d + 5)      # 5-digit margin (rank-1 iperp)
        keep = [j for j in range(K + 1) if nrm[j] > thr]
        S = mp.mpf(10) ** d
        def Lij(i, j):
            if j == 0: return iperp if i == 0 else mp.mpf(0)
            return p[j - 1] if i == 0 else L_F[i - 1, j - 1]
        rows = [[(1 if j == i else 0) for j in range(K + 1)]
                + [int(mp.nint(Lij(i, j) / nrm[j] * S)) for j in keep]
                for i in range(K + 1)]
    return rows, span_res


def qrbkz_batch(I_dict, F, names, fit, ho, d, hmax=10 ** 30, backend=None, beta=None):
    """Multi-target qrbkz: LQ(F) once, rank-1 update per target, batch reduce."""
    K = len(F); beta = beta or min(K + 1, 40)
    L_F, Q_F, cn, col, Ne = _basis_LQ(F, fit, d)
    latt, spans, keys = [], [], list(I_dict)
    for k in keys:
        p, ip = _target_row(I_dict[k], Q_F, col, fit, Ne, d)
        rows, sr = _embed(L_F, cn, p, ip, d, K)
        latt.append(rows); spans.append(sr)
    Rs, bknd, wall = reduce_rows_batch(latt, backend, bkz_beta=beta)
    out = {}
    for k, R, sr in zip(keys, Rs, spans):
        rel, brels = mplll._extract(R, K, hmax)
        o = {"method": "qrbkz-batch", "K": K, "dps": d, "backend": bknd,
             "span_residual_d": round(sr, 2), "relation": rel,
             "basis_relations": [list(b) for b in sorted(set(brels))] if brels else None}
        if rel:
            o["heldout_min_d"], _ = mplll.heldout_digits(rel, I_dict[k], F, ho, d)
            o["insample_min_d"], _ = mplll.heldout_digits(rel, I_dict[k], F, fit, d)
            o["capacity"] = mplll.capacity(rel, K, len(fit), d)
            o["c0"], o["coeffs"] = rel[0], {names[i]: -rel[i + 1] for i in range(K)}
            o["status"] = "HIT" if o["heldout_min_d"] >= mplll.GATE_D else "PARTIAL"
        else:
            o["status"] = "NULL"
        out[k] = o
    out["_meta"] = {"lll_wall_s": round(wall, 3), "n_targets": len(keys), "backend": bknd}
    return out


def cvp_fit(I, F, names, fit, ho, d, denom, c_LS=None, hmax=10 ** 30, backend=None, beta=None):
    """CVP/Babai: LLL-reduce basis-only lattice L_F, Babai nearest-plane on target
    t = D·c_LS ⊕ 10^d·D·p/cn → nearest integer coeff vector.  c_LS from
    L_F^T c = p (LS via LQ) if not supplied."""
    K = len(F); beta = beta or min(K, 40)
    L_F, Q_F, cn, col, Ne = _basis_LQ(F, fit, d)
    p, _ = _target_row(I, Q_F, col, fit, Ne, d)
    with mp.workdps(d + 30):
        if c_LS is None:
            c_LS = mp.lu_solve(L_F.T, mp.matrix(p)); c_LS = [c_LS[i] for i in range(K)]
        else:
            c_LS = [mp.mpf(str(x)) for x in c_LS]
        S = mp.mpf(10) ** d
        rowsF = [[(1 if j == i else 0) for j in range(K)]
                 + [int(mp.nint(L_F[i, j] / cn[j] * S)) for j in range(K)] for i in range(K)]
        R, bknd, wall = reduce_rows(rowsF, backend, bkz_beta=beta)
        t = [denom * c_LS[i] for i in range(K)] + [denom * p[j] / cn[j] * S for j in range(K)]
        M = 2 * K; G = [[mp.mpf(x) for x in r] for r in R]
        for i in range(K):                                   # Gram–Schmidt of R
            for j in range(i):
                m = mp.fsum(G[i][k] * G[j][k] for k in range(M)) / \
                    mp.fsum(G[j][k] ** 2 for k in range(M))
                for k in range(M): G[i][k] -= m * G[j][k]
        v = list(t)
        for i in range(K - 1, -1, -1):                       # Babai nearest-plane
            m = int(mp.nint(mp.fsum(v[k] * G[i][k] for k in range(M)) /
                            mp.fsum(G[i][k] ** 2 for k in range(M))))
            for k in range(M): v[k] -= m * R[i][k]
        c = [int(mp.nint(t[i] - v[i])) for i in range(K)]    # nearest pt, left block
    rel = tuple([denom] + [-x for x in c])
    g = mplll.gcd_list(rel) or 1; rel = tuple(x // g for x in rel)
    if rel[0] < 0: rel = tuple(-x for x in rel)
    hd, _ = mplll.heldout_digits(rel, I, F, ho, d)
    isd, _ = mplll.heldout_digits(rel, I, F, fit, d)
    return {"method": "cvp", "K": K, "dps": d, "denom": denom, "backend": bknd,
            "lll_wall_s": round(wall, 3), "relation": rel,
            "c0": rel[0], "coeffs": {names[i]: -rel[i + 1] for i in range(K)},
            "heldout_min_d": hd, "insample_min_d": isd,
            "capacity": mplll.capacity(rel, K, len(fit), d),
            "status": "HIT" if hd >= mplll.GATE_D else "PARTIAL"}


def dmin_sweep(I, F, names, fit, ho, dps, W=10 ** 4, method="qrbkz", backend=None):
    """Bisect on d∈[4,dps]: smallest d giving two-prec-stable rel == d=dps rel."""
    ref = mplll.mplll_fit(I, F, names, fit, ho, dps, W, method=method,
                          two_prec=True, backend=backend)
    if not ref["relation"]:
        return {"d_min": None, "reference": None, "d_tested": []}
    want = canonicalize(ref["relation"]); tested = [dps]
    def ok(d):
        tested.append(d)
        r = mplll.mplll_fit(I, F, names, fit, ho, d, W, method=method,
                            two_prec=True, backend=backend)
        return (r["relation"] and canonicalize(r["relation"]) == want
                and r.get("two_prec_stable", False))
    lo, hi = 4, dps
    if ok(lo): hi = lo
    while hi - lo > 1:
        mid = (lo + hi) // 2; lo, hi = (lo, mid) if ok(mid) else (mid, hi)
    return {"d_min": hi, "reference": list(want), "d_tested": sorted(set(tested))}

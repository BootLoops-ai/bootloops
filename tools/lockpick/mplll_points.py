#!/usr/bin/env python3
"""mplll_points.py — D-optimal sample-point selector for multipoint lattice regression.
Answers arXiv:2507.17815 §5.3 (κ gradient-descent on target side) with the standard
D-optimal design on the BASIS side: pick x_j maximizing |det(E E^T)|, E_{ij}=B_i(x_j).
Greedy: k<K = Leja/Fekete (column-pivoted MGS pivot); k≥K = leverage argmax
b^T G^{-1} b with Sherman–Morrison downdate.  Basis is cheap, target (AMFlow) is
expensive → select points BEFORE any target eval."""
import argparse, json, sys
import mpmath as mp


def _greedy_dopt(B, n_select, tol):
    """B: list of M columns, each len-K list[mpf].  Returns (sel_idx, gains) where
    gains[k] = log10 multiplicative gain in sqrt(det E E^T) at step k."""
    K, M = len(B[0]), len(B)
    R = [list(c) for c in B]                          # residuals (phase-1)
    nrm2 = [mp.fsum(v * v for v in c) for c in R]
    sel, gains, Q, S = [], [], [], set()
    # ---- phase 1: k < K, column-pivoted MGS ----
    for k in range(min(n_select, K, M)):
        j = max((i for i in range(M) if i not in S), key=lambda i: nrm2[i])
        p = mp.sqrt(max(nrm2[j], mp.mpf(0)))
        if gains and p < tol * mp.mpf(10) ** gains[0]:
            break
        sel.append(j); S.add(j); gains.append(mp.log10(p)); q = [v / p for v in R[j]]
        Q.append(q)
        for i in range(M):
            if i in S: continue
            d = mp.fsum(q[r] * R[i][r] for r in range(K))
            for r in range(K): R[i][r] -= d * q[r]
            nrm2[i] -= d * d
            if nrm2[i] < 0: nrm2[i] = mp.fsum(v * v for v in R[i])
    # ---- phase 2: k ≥ K, leverage argmax with Sherman–Morrison ----
    if len(sel) < n_select and len(sel) < M:
        G = mp.zeros(K)
        for j in sel:
            for a in range(K):
                for b in range(K): G[a, b] += B[j][a] * B[j][b]
        Ginv = G ** -1
        while len(sel) < min(n_select, M):
            best, jb, vb = mp.mpf(-1), None, None
            for i in range(M):
                if i in S: continue
                v = Ginv * mp.matrix(B[i])
                lev = mp.fsum(B[i][r] * v[r] for r in range(K))
                if lev > best: best, jb, vb = lev, i, v
            sel.append(jb); S.add(jb)
            gains.append(mp.log10(1 + best) / 2)
            den = 1 + best
            for a in range(K):
                for b in range(K): Ginv[a, b] -= vb[a] * vb[b] / den
    return sel, [float(g) for g in gains]


def select_from_matrix(B_cols, n_select, dps=60):
    """B_cols: M columns × K basis values (mpf-able). Returns (sel_idx, log10_det, gains)."""
    with mp.workdps(dps):
        B = [[mp.mpf(v) for v in c] for c in B_cols]
        sel, gains = _greedy_dopt(B, n_select, mp.mpf(10) ** (-dps // 2))
    return sel, sum(gains), gains


def select_points(basis_fn, candidate_xs, K, n_select, dps=60):
    """basis_fn: x -> list[K mpf].  Evaluates basis at all candidates then greedy-picks
    n_select D-optimal points.  Returns (selected_xs, log10_det_L, per_step_gain).
    Early stop in phase-1 if pivot < 10^{-dps/2}·pivot_0 (rank-deficient basis)."""
    xs = list(candidate_xs)
    with mp.workdps(dps):
        B = [[mp.mpf(v) for v in basis_fn(x)] for x in xs]
    if any(len(c) != K for c in B):
        raise ValueError(f"basis_fn must return length-{K} vectors")
    sel, logdet, gains = select_from_matrix(B, n_select, dps)
    return [xs[j] for j in sel], logdet, gains


def cond_E(B_cols, dps=60):
    """κ(E)=σ_max/σ_min for K×N E with columns B_cols (via eig of E E^T)."""
    with mp.workdps(dps):
        K, N = len(B_cols[0]), len(B_cols)
        E = mp.matrix(K, N)
        for j, c in enumerate(B_cols):
            for i in range(K): E[i, j] = mp.mpf(c[i])
        ev = mp.mp.eigsy(E * E.T)[0]
        s = sorted(mp.sqrt(max(e, mp.mpf(0))) for e in ev)
        return float(s[-1] / s[0]) if s[0] > 0 else float('inf')


def _cli():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--from-json", required=True,
                    help='{"candidates":[x..], "basis_values":{name:[B_i(x_j)..]}}')
    ap.add_argument("--n-select", type=int, required=True)
    ap.add_argument("--dps", type=int, default=60)
    ap.add_argument("--out", default="-")
    a = ap.parse_args()
    d = json.load(open(a.from_json))
    xs = d["candidates"]; names = list(d["basis_values"].keys())
    K, M = len(names), len(xs)
    B_cols = [[d["basis_values"][n][j] for n in names] for j in range(M)]
    sel, logdet, gains = select_from_matrix(B_cols, a.n_select, a.dps)
    Bsel = [B_cols[j] for j in sel]
    out = {"K": K, "M": M, "n_selected": len(sel), "selected_idx": sel,
           "selected_xs": [xs[j] for j in sel], "log10_det_L": round(logdet, 4),
           "per_step_gain": [round(g, 4) for g in gains],
           "kappa_E_selected": cond_E(Bsel, a.dps) if len(sel) >= K else None}
    js = json.dumps(out, indent=1, default=str)
    (sys.stdout if a.out == "-" else open(a.out, "w")).write(js + "\n")


if __name__ == "__main__":
    _cli()

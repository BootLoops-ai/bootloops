#!/usr/bin/env python3
# lp_kill_rank_agnostic.py — Venkov 2-design LP kill per (K class, Coxeter h),
# rank-agnostic complements: N = len(GK) is read from the input Gram via
# set_rank() (IJ symmetric entries, the LP row count, and the Farkas dim
# nIJ+1 are all sized from it). Exact Fractions det (no float det).
# The s_K>=2 Leech guard (no root-free complement class) is enforced HERE, in
# analyze_class; an independent verifier must enforce it too. Kills are
# exact rational Farkas certificates; no float step load-bearing.
# Generalization contract, input format and modes:
# lattice/LP_KILL_RANK_GENERALIZATION.md. Exercised at complement ranks
# c=3, c=4, c=6, c=7 and the original c=5 Niemeier case.
# short_vectors/cholesky_frac: a Python implementation of Fincke-Pohst enumeration
# (U. Fincke and M. Pohst, Math. Comp. 44 (1985) 463-471; H. Cohen, A Course in
# Computational Algebraic Number Theory, Algorithm 2.7.5). The exact-rational
# Cholesky step, the search window isqrt(floor(T_i/q_ii)) + 2 about -U_i, and the
# output contract (one (vector, norm) pair per +-v with 0 < norm <= bound) follow
# Hecke.jl, src/QuadForm/Enumeration.jl, by the Hecke contributors (Claus Fieker,
# Tommy Hofmann and others; https://github.com/thofma/Hecke.jl), BSD 2-Clause
# License; the recursive control flow, the sign deduplication and everything else
# in this file are our own.
# NIEMEIER_H: the 18 distinct Coxeter numbers of the 23 Niemeier lattices with
# roots (H.-V. Niemeier, J. Number Theory 5 (1973) 142-178; Conway and Sloane,
# Sphere Packings, Lattices and Groups, Chapter 16, Table 16.1).
import sys, json, math
from fractions import Fraction as Fr
N = IJ = NROWS = None                     # set per genus by set_rank(len(GK))
NIEMEIER_H = [2,3,4,5,6,7,8,9,10,12,13,14,16,18,22,25,30,46]  # 23 rooted Niemeiers
def set_rank(n):
    global N, IJ, NROWS
    N = n
    IJ = [(i, j) for i in range(N) for j in range(i, N)]      # n(n+1)/2 sym entries
    NROWS = len(IJ) + 1                                        # Farkas dim nIJ+1
def cholesky_frac(A):
    n = len(A); Q = [[Fr(A[i][j]) for j in range(n)] for i in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            Q[j][i] = Q[i][j]; Q[i][j] = Q[i][j] / Q[i][i]
        for k in range(i+1, n):
            for l in range(k, n):
                Q[k][l] = Q[k][l] - Q[k][i] * Q[i][l]
    return Q
def short_vectors(A, bound):
    n = len(A); Q = cholesky_frac(A); out = []
    def rec(i, x, T):
        if T < 0: return
        ub = math.isqrt(int(T / Q[i][i])) + 2
        c = -sum(Q[i][j] * x[j] for j in range(i+1, n))
        xi = math.floor(c - ub)
        while xi <= c + ub:
            d = Q[i][i] * (Fr(xi) - c)**2
            if d <= T:
                if i == 0:
                    v = [xi] + x[1:]
                    nv = sum(v[a]*A[a][b]*v[b] for a in range(n) for b in range(n))
                    if 0 < nv <= bound: out.append((v, nv))
                else:
                    x2 = list(x); x2[i] = xi; rec(i-1, x2, T - d)
            xi += 1
    rec(n-1, [0]*n, Fr(bound))
    seen = set(); res = []
    for v, nv in out:
        t = tuple(v)
        if tuple(-a for a in v) in seen: continue
        seen.add(t); res.append((list(v), int(nv)))
    return res
def det_frac(A):
    n = len(A); M = [[Fr(x) for x in row] for row in A]; d = Fr(1)
    for c in range(n):
        p = next((r for r in range(c, n) if M[r][c] != 0), None)
        if p is None: return Fr(0)
        if p != c: M[c], M[p] = M[p], M[c]; d = -d
        d *= M[c][c]
        for r in range(c+1, n):
            f = M[r][c] / M[c][c]
            M[r] = [a - f*b for a, b in zip(M[r], M[c])]
    return d
def matinv_frac(A):
    n = len(A); M = [[Fr(A[i][j]) for j in range(n)] + [Fr(int(i==j)) for j in range(n)] for i in range(n)]
    for c in range(n):
        p = next(r for r in range(c, n) if M[r][c] != 0)
        M[c], M[p] = M[p], M[c]; piv = M[c][c]
        M[c] = [e / piv for e in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c]; M[r] = [a - f*b for a, b in zip(M[r], M[c])]
    return [row[n:] for row in M]
def analyze_class(GK):
    set_rank(len(GK))                     # documented delta: rank from input
    d = det_frac(GK); assert d.denominator == 1 and d > 0
    detG = int(d)
    Ginv = matinv_frac(GK)
    adj = [[Ginv[i][j] * detG for j in range(N)] for i in range(N)]
    assert all(x.denominator == 1 for row in adj for x in row)
    adj = [[int(x) for x in row] for row in adj]
    roots = [v for v, nv in short_vectors(GK, 2) if nv == 2]
    sK = 2 * len(roots)
    assert sK >= 2, "LEECH GUARD FAIL: root-free complement class (s_K=0)"
    S = [[0]*N for _ in range(N)]
    for x in roots:
        Gx = [sum(GK[i][j]*x[j] for j in range(N)) for i in range(N)]
        for i in range(N):
            for j in range(N):
                S[i][j] += 2 * Gx[i] * Gx[j]      # both signs +-x
    shells = []
    for u, a in short_vectors(adj, 2*detG - 1):   # 0 < a = D*|v|^2 < 2D
        shells.append((u, a, 2 * ((2*detG) // a)))  # cap 2*floor(2/t)
    return detG, sK, S, shells
def lp_test(GK, detG, sK, S, shells, h):
    import numpy as np
    from scipy.optimize import linprog
    Nv = len(shells)
    b = [2*h*GK[i][j] - S[i][j] for (i, j) in IJ] + [24*h - sK]
    A = np.zeros((NROWS, Nv))
    for k, (u, a, cap) in enumerate(shells):
        for r, (i, j) in enumerate(IJ):
            A[r, k] = u[i]*u[j]
        A[NROWS-1, k] = 1
    caps = [c for (_, _, c) in shells]
    res = linprog(np.zeros(Nv), A_eq=A, b_eq=np.array(b, dtype=float),
                  bounds=list(zip([0]*Nv, caps)), method='highs')
    if res.status == 0:
        return "FEASIBLE", None, res.x
    cvec = np.concatenate([np.array(b, dtype=float), np.array(caps, dtype=float)])
    A_ub = np.hstack([-A.T, -np.eye(Nv)]); b_ub = np.zeros(Nv)
    bnds = [(-1, 1)]*NROWS + [(0, None)]*Nv
    r2 = linprog(cvec, A_ub=A_ub, b_ub=b_ub, bounds=bnds, method='highs')
    if r2.status != 0 or r2.fun > -1e-9:
        return "INFEASIBLE-NOCERT", None, None
    yf = r2.x[:NROWS]
    for den in (10**3, 10**5, 10**7, 10**9, 10**12):
        y = [Fr(v).limit_denominator(den) for v in yf]
        val = sum(yi*bi for yi, bi in zip(y, b))
        for k, (u, a, cap) in enumerate(shells):
            col = sum(y[r]*u[i]*u[j] for r, (i, j) in enumerate(IJ)) + y[NROWS-1]
            if col < 0: val += cap * (-col)
        if val < 0:
            return "KILLED", [str(v) for v in y], None
    return "INFEASIBLE-NOCERT", None, None
def exact_feasible(shells, b, x):
    # try to lift float-feasible point to exact rational witness (Ax=b, 0<=x<=cap)
    sup = [k for k in range(len(shells)) if x[k] > 1e-7]
    cols = {k: [Fr(shells[k][0][i]*shells[k][0][j]) for (i, j) in IJ] + [Fr(1)] for k in sup}
    for den in (1, 2, 3, 4, 6, 12, 24):
        M = [[cols[k][r] for k in sup] + [Fr(b[r])] for r in range(NROWS)]
        piv = []  # (row, colidx-in-sup)
        rr = 0
        for ci in range(len(sup)):
            p = next((r for r in range(rr, NROWS) if M[r][ci] != 0), None)
            if p is None: continue
            M[rr], M[p] = M[p], M[rr]
            M[rr] = [e / M[rr][ci] for e in M[rr]]
            for r in range(NROWS):
                if r != rr and M[r][ci] != 0:
                    f = M[r][ci]; M[r] = [a - f*c for a, c in zip(M[r], M[rr])]
            piv.append((rr, ci)); rr += 1
        if any(all(M[r][c] == 0 for c in range(len(sup))) and M[r][-1] != 0
               for r in range(NROWS)):
            return None
        val = {}
        for ci in range(len(sup)):
            if ci not in [c for _, c in piv]:
                val[ci] = Fr(round(x[sup[ci]]*den), den)
        for r, ci in reversed(piv):
            val[ci] = M[r][-1] - sum(M[r][c]*val[c] for c in range(len(sup))
                                     if c != ci and c in val)
        m = [Fr(0)]*len(shells)
        for ci, k in enumerate(sup): m[k] = val[ci]
        ok = all(0 <= m[k] <= shells[k][2] for k in range(len(shells)))
        if ok:
            for r in range(NROWS):
                if r < NROWS - 1:
                    i, j = IJ[r]
                    lhs = sum(m[k]*shells[k][0][i]*shells[k][0][j] for k in sup)
                else:
                    lhs = sum(m[k] for k in sup)
                if lhs != Fr(b[r]): ok = False; break
        if ok:
            return {str(k): str(m[k]) for k in sup}
    return None
def parse_hecke(path, key):
    classes = []
    for line in open(path):
        p = line.strip().split("|")
        if p[0] == "HCLS" and p[1] == key:
            rows = p[6].strip()[1:-1].split(";")
            GK = [[int(x) for x in r.split(",")] for r in rows]
            classes.append({"idx": int(p[2]), "GK": GK, "aut": int(p[4])})
        if p[0] == "HKEY" and p[1] == key:
            info = {"det": p[3], "mass": p[5], "n": int(p[7]), "sum": p[9],
                    "match": p[11]}
    return classes, info
def run_genus(key, hecke_out, dst):
    classes, info = parse_hecke(hecke_out, key)
    assert info["match"] == "true", "MASS GATE FAIL: enumeration incomplete"
    rec = {"key": key, "hecke": info, "nclasses": len(classes), "H": NIEMEIER_H,
           "classes": [], "closed": False, "closing_class": None}
    for cl in classes:
        detG, sK, S, shells = analyze_class(cl["GK"])
        cres = {"idx": cl["idx"], "GK": cl["GK"], "aut": cl["aut"], "det": detG,
                "sK": sK, "nshells": len(shells), "perh": {}, "certs": {},
                "feas_witness": {}}
        allkill = True
        for h in NIEMEIER_H:
            status, cert, xf = lp_test(cl["GK"], detG, sK, S, shells, h)
            if status == "FEASIBLE":
                b = [2*h*cl["GK"][i][j] - S[i][j] for (i, j) in IJ] + [24*h - sK]
                w = exact_feasible(shells, b, xf)
                status = "FEASIBLE-EXACT" if w else "FEASIBLE-FLOAT"
                if w: cres["feas_witness"][str(h)] = w
            cres["perh"][str(h)] = status
            if status == "KILLED" and cert: cres["certs"][str(h)] = cert
            if status != "KILLED": allkill = False
        cres["all_h_killed"] = allkill
        rec["classes"].append(cres)
        if allkill and not rec["closed"]:
            rec["closed"] = True; rec["closing_class"] = cl["idx"]
    json.dump(rec, open(dst, "w"), indent=1)
    st = "CLOSED-ALLROOTED" if rec["closed"] else "NOT-CLOSED"
    print(f"{key}|{st}|classes={rec['nclasses']}|" +
          ";".join(f"K{c['idx']}(sK={c['sK']}):" + ("ALLKILL" if c["all_h_killed"]
          else ",".join(h for h, s in c["perh"].items() if s != "KILLED") + "-open")
          for c in rec["classes"]))
def control_a16():
    import numpy as np
    from scipy.optimize import linprog
    GK = [[2*int(i==j) for j in range(6)] for i in range(6)]  # A1^6 control fixed at c=6
    detG, sK, S, shells = analyze_class(GK)
    assert (detG, sK) == (64, 12), (detG, sK)
    st, cert, _ = lp_test(GK, detG, sK, S, shells, 2)
    b = [2*2*GK[i][j] - S[i][j] for (i, j) in IJ] + [24*2 - sK]
    assert all(x == 0 for x in b[:-1]) and b[-1] == 36
    Nv = len(shells)
    A = np.zeros((NROWS, Nv+1))
    for k, (u, a, cap) in enumerate(shells):
        for r, (i, j) in enumerate(IJ): A[r, k] = u[i]*u[j]
        A[NROWS-1, k] = 1
    A[NROWS-1, Nv] = 1  # zero-projection column restored
    res = linprog(np.zeros(Nv+1), A_eq=A, b_eq=np.array(b, dtype=float),
                  bounds=[(0, None)]*(Nv+1), method='highs')  # caps dropped
    print(f"CONTROL-A1^6|det={detG}|sK={sK}|strict h=2: {st}|relaxed: "
          f"{'FEASIBLE' if res.status == 0 else 'INFEASIBLE'}|n0={res.x[Nv]:.6f} (expect 36)")
if __name__ == "__main__":
    if sys.argv[1] == "control": control_a16()
    else: run_genus(sys.argv[1], sys.argv[2], sys.argv[3])

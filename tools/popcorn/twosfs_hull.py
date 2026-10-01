#!/usr/bin/env python3
"""twosfs_hull — exact Euclidean projection of a normalized 2-SFS onto the
convex hull of finitely many others (support-restricted mixture test).

Register: EXACT. fractions.Fraction end to end; standard library only.

Setting. A symmetric (n-1) x (n-1) moment matrix is stored as a dict over the
pairs (i, j), 1 <= i <= j <= n-1 (pairs(n)), and compared in the
multiplicity-weighted Frobenius inner product
    <A, B> = sum_{i<=j} mult_ij A_ij B_ij,   mult = 1 on the diagonal, 2 off it
(dot), i.e. the ordinary Frobenius product of the full symmetric matrices.

project_hull(n, x, comps) computes the EXACT Euclidean projection of x onto
conv(comps) by active-set enumeration: for every subset S of the components
it solves the equality-constrained least-squares KKT system on the Gram
matrix (solve_exact, Fraction Gaussian elimination), keeps solutions with all
weights >= 0 that satisfy the global variational inequality
<x - p, P_k - p> <= 0 for EVERY component k (which characterizes the
projection onto a convex set), and returns the one of least squared distance:
    {'S': active subset, 'lam': exact convex weights, 'p': projection,
     'w': x - p, 'd2': exact squared distance}.
d2 == 0 iff x is in the hull; when d2 > 0, w = x - p is a separating
functional (<w, P_k - p> <= 0 < <w, x - p> = d2 for all k) and
d2 / max|w_ij| is a lower bound on the l1 distance from x to the hull.
Cost is 2^m subsets for m components: meant for m up to about a dozen.

Use with twosfs_engine: x = the normalized exact 2-SFS q of a Lambda-
coalescent (normalize_matrix), comps = the normalized single-time Kingman
moments component_moments(n, t_k) at a finite support {t_k}. Because an atom
history's 1-SFS scales with its dwell D and its 2-SFS with D^2, the 2-SFS
mixtures achievable by pooling atom histories on a FIXED support, with the
per-atom weights free, sweep exactly conv{normalized m2(t_k)}; d2 > 0 then
says no pooling of atom histories on that support reproduces q. This is the
weakest of the three registers (it says nothing about other supports or about
genuinely time-varying histories: the constant-size Kingman coalescent itself
sits outside the hull of any finite atom support); the support-free and
full-class certificates are in twosfs_certificates.

Also: solve_exact(A, b) (square system, returns None if singular; note that
popcorn.certificates.region.solve_exact is a different, rectangular routine),
q_of_moments(n, v) == twosfs_engine.normalize_matrix(v, n), mult(i, j).
"""
from fractions import Fraction as F
from itertools import combinations
from twosfs_engine import normalize_matrix

def pairs(n):
    return [(i, j) for i in range(1, n) for j in range(i, n)]

def mult(i, j):
    return 1 if i == j else 2

def dot(n, a, b):
    return sum(mult(i, j) * a[(i, j)] * b[(i, j)] for (i, j) in pairs(n))

def solve_exact(Amat, bvec):
    """Gaussian elimination, Fractions. Returns x or None if singular."""
    m = len(Amat)
    M = [row[:] + [bvec[r]] for r, row in enumerate(Amat)]
    for c in range(m):
        piv = next((r for r in range(c, m) if M[r][c] != 0), None)
        if piv is None:
            return None
        M[c], M[piv] = M[piv], M[c]
        pv = M[c][c]
        M[c] = [x / pv for x in M[c]]
        for r in range(m):
            if r != c and M[r][c] != 0:
                f = M[r][c]
                M[r] = [x - f * y for x, y in zip(M[r], M[c])]
    return [M[r][m] for r in range(m)]

def project_hull(n, x, comps):
    """Exact projection of x onto conv(comps) (dicts over pairs).

    Enumerate subsets S; solve min ||x - sum lam_k P_k||^2, sum lam = 1 via
    KKT on the Gram system; keep solutions with lam >= 0 that satisfy the
    global variational inequality <x - p, P_k - p> <= 0 for ALL k."""
    m = len(comps)
    best = None
    for r in range(1, m + 1):
        for S in combinations(range(m), r):
            # KKT: for k in S: sum_l lam_l <P_k, P_l> + nu = <P_k, x>; sum lam = 1
            k_ = len(S)
            A = [[dot(n, comps[S[a]], comps[S[b]]) for b in range(k_)] + [F(1)]
                 for a in range(k_)]
            A.append([F(1)] * k_ + [F(0)])
            b = [dot(n, comps[S[a]], x) for a in range(k_)] + [F(1)]
            sol = solve_exact(A, b)
            if sol is None:
                continue
            lam = sol[:k_]
            if any(l < 0 for l in lam):
                continue
            p = {ij: sum(lam[a] * comps[S[a]][ij] for a in range(k_))
                 for ij in pairs(n)}
            w = {ij: x[ij] - p[ij] for ij in pairs(n)}
            if all(dot(n, w, {ij: comps[k][ij] - p[ij] for ij in pairs(n)}) <= 0
                   for k in range(m)):
                d2 = dot(n, w, w)
                if best is None or d2 < best['d2']:
                    best = {'S': S, 'lam': lam, 'p': p, 'w': w, 'd2': d2}
    return best

def q_of_moments(n, v):
    q, tot = normalize_matrix(v, n)
    return q, tot


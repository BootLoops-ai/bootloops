#!/usr/bin/env python3
# frame_rank5.py (DIAG-SLICE) — derive, from the 29x29 integer frame
# matrices of the 6-parameter HV4 family, that the S6-invariant sector of
# the 29-frame is EXACTLY rank 5 and MUM with unipotency index 5 under
# T = T1..T6 — the frame-side pin that the diagonal invariant period system
# has order 5.  Output: FRAME_RANK5.json (shipped, sha-pinned; mkcard.py
# copies it into the card).
# Input: the directory TERRIER_HV4_FRAME_DIR (not included in the package)
# holding the 29x29 integer matrices MSWAP (S6 transposition), M6 (S6
# 6-cycle), T1..T6 (large-complex-structure monodromies) and SIGMA
# (intersection form) of the arXiv:2404.12422 sec 5 frame, one JSON array per
# matrix in a file <name>.json.
# Note: M6 satisfies (5.59) as displayed = transpose of the printed (5.58)
# block form; M6 is a slot PERMUTATION, so the +1-eigenspace is identical
# for M6 and its transpose — the convention delta is immaterial here
# (asserted below).
import json, os
from fractions import Fraction as Fr
# frame matrix dir from env (required; fail-closed when unset).
FD = os.environ.get("TERRIER_HV4_FRAME_DIR")
if not FD:
    raise SystemExit("REFUSED: TERRIER_HV4_FRAME_DIR is not set (frame matrix dir)")
L = lambda f: json.load(open(os.path.join(FD, f)))
RS, RC = L("MSWAP.json"), L("M6.json")
TS = [L(f"T{k}.json") for k in range(1, 7)]
SIG = L("SIGMA.json")
n = 29
def mm(A, B):
    return [[sum(A[i][k]*B[k][j] for k in range(len(B)))
             for j in range(len(B[0]))] for i in range(len(A))]
def rank_frac(M):
    M = [[Fr(x) for x in r] for r in M]; r = 0; rows = len(M); cols = len(M[0])
    for c in range(cols):
        p = next((i for i in range(r, rows) if M[i][c] != 0), None)
        if p is None: continue
        M[r], M[p] = M[p], M[r]
        M[r] = [x / M[r][c] for x in M[r]]
        for i in range(rows):
            if i != r and M[i][c] != 0:
                M[i] = [a - M[i][c]*b for a, b in zip(M[i], M[r])]
        r += 1
    return r
# permutation check + transpose-equivalence of the +1 eigenspace
assert all(sorted(r) == [0]*28+[1] for r in RC) and \
       all(sorted(c) == [0]*28+[1] for c in zip(*RC)), "M6 not a permutation"
# invariant sector: kernel of stacked (RS - I; RC - I); dim = 29 - rank
I29 = [[1 if i == j else 0 for j in range(n)] for i in range(n)]
stack = [[RS[i][j] - I29[i][j] for j in range(n)] for i in range(n)] + \
        [[RC[i][j] - I29[i][j] for j in range(n)] for i in range(n)]
inv_dim = n - rank_frac(stack)
# orbit-sum basis: the 5 S6 orbits of the 29 frame slots
# (slot 0; slots 1-6; slots 7-21; slots 22-27; slot 28)
orbits = [[0], list(range(1, 7)), list(range(7, 22)), list(range(22, 28)), [28]]
B = [[1 if j in orb else 0 for orb in orbits] for j in range(n)]  # 29x5
ok_basis = all(all(sum(R[i][j]*B[j][b] for j in range(n)) == B[i][b]
               for i in range(n)) for R in (RS, RC) for b in range(5))
# T = T1..T6 restricted to the invariant sector
T = TS[0]
for k in range(1, 6):
    T = mm(T, TS[k])
TB = mm(T, B)                       # 29x5
reps = [orb[0] for orb in orbits]   # representative slot per orbit
N5 = [[TB[reps[i]][b] for b in range(5)] for i in range(5)]
ok_closed = all(TB[i][b] == sum(B[i][j]*N5[j][b] for j in range(5))
                for i in range(n) for b in range(5))
U = [[N5[i][j] - (1 if i == j else 0) for j in range(5)] for i in range(5)]
U2 = mm(U, U); U4 = mm(U2, U2); U5 = mm(U4, U)
mum_ok = any(any(r) for r in U4) and not any(any(r) for r in U5)
# Sigma restricted to the invariant basis (5x5 Gram) + signature
SB = mm([[SIG[i][j] for j in range(n)] for i in range(n)], B)
G5 = [[sum(B[i][a]*SB[i][b] for i in range(n)) for b in range(5)] for a in range(5)]
sym_ok = all(G5[i][j] == G5[j][i] for i in range(5) for j in range(5))
# char poly by Faddeev-LeVerrier (exact), then Descartes (all roots real)
def charpoly(A):
    m = len(A); Mk = [r[:] for r in I5] if False else None
    c = [Fr(1)]; Mk = [[Fr(0)]*m for _ in range(m)]
    for i in range(m): Mk[i][i] = Fr(1)
    Ak = [[Fr(x) for x in r] for r in A]
    Mcur = [r[:] for r in Mk]
    for k in range(1, m + 1):
        AM = mm(Ak, Mcur)
        ck = -Fr(sum(AM[i][i] for i in range(m)), k)
        c.append(ck)
        Mcur = [[AM[i][j] + (ck if i == j else 0) for j in range(m)] for i in range(m)]
    return c  # lambda^m + c1 lambda^(m-1) + ... + cm
cp = charpoly(G5)
def signchanges(cs):
    s = [x for x in cs if x != 0]
    return sum(1 for a, b in zip(s, s[1:]) if (a > 0) != (b > 0))
pos = signchanges(cp)                                  # roots all real (sym)
neg = signchanges([c if i % 2 == 0 else -c for i, c in enumerate(cp)])
detG5 = cp[-1] * (1 if len(cp) % 2 == 1 else -1)       # (-1)^5 * cm ... fix:
# det(A) = (-1)^m * cp_m  with cp monic char poly of A: p(0) = cm = (-1)^m det
detG5 = (-1)**5 * cp[-1]
rec = {"invariant_sector_dim": inv_dim, "orbit_basis_ok": ok_basis,
       "Tprod_preserves_sector": ok_closed,
       "N5_restricted": N5, "mum_index5": mum_ok,
       "U4_nonzero": any(any(r) for r in U4),
       "gram5": [[str(x) for x in r] for r in G5], "gram5_symmetric": sym_ok,
       "gram5_det": str(detG5), "gram5_signature": [pos, neg],
       "charpoly_gram5": [str(x) for x in cp]}
json.dump(rec, open("FRAME_RANK5.json", "w"), indent=1)
print("FRAME-RANK5", {k: rec[k] for k in ("invariant_sector_dim", "orbit_basis_ok",
      "Tprod_preserves_sector", "mum_index5", "gram5_det", "gram5_signature")})
assert inv_dim == 5 and ok_basis and ok_closed and mum_ok and sym_ok

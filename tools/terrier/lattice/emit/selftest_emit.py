#!/usr/bin/env python3
"""selftest_emit.py — battery for the emit/ glue-extraction predicate chain on
the shipped fixtures (exact arithmetic, no external data needed).

The Julia engines (harness.jl / s11.jl) regenerate the full emission lists;
they need an Oscar/Hecke project env and are exercised separately (README.md
run law). This battery checks the MATHEMATICAL CONTENT of the glue extraction
on shipped fixtures, in exact integer/Fraction arithmetic:

E1 known-answer table: all 44 rows of bank/sound_minvecs.json verify through
   a Python port of the tau-extraction identities (invertible-glue case):
   with A = T1, At = T2, B = N*T2, M = B*At^-1 (= N), Bt = At^-1 * B^T * A,
   S0 = B*Bt — dual-integrality (M, Bt integral), Gram symmetry/evenness/
   positive-definiteness, trace identity nmin_lat == tr(S0) == 2*nflux_min,
   even trace, and the det-chain identity |det S0| == det(M)^2 |det A det At|.
E2 closed-form synthetic case: T1 = T2 = A2 Gram [[2,1],[1,2]] glued by
   N = [[0,1],[-1,0]] gives S0 = 3*I exactly (trace 6, det-chain 9 == 9),
   derived by hand in this file — the port must reproduce it.
E3 MUTATION (must-fail): a tampered glue entry and a tampered nmin_lat value
   must each be caught by the E1 checks.

Exit 0 iff all gates pass.  Run: ulimit -v 32505856; nice -n 5 python3
selftest_emit.py
"""
import json
import os
import sys
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
GATES = []


def gate(name, ok, detail=""):
    GATES.append((name, bool(ok)))
    print("  %-30s %s  %s" % (name, "PASS" if ok else "FAIL", detail))


def det2(M):
    return M[0][0] * M[1][1] - M[0][1] * M[1][0]


def mm(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(2)) for j in range(2)]
            for i in range(2)]


def tp(M):
    return [[M[j][i] for j in range(2)] for i in range(2)]


def inv2(M):
    d = Fr(det2(M))
    return [[M[1][1] / d, -M[0][1] / d], [-M[1][0] / d, M[0][0] / d]]


def check_row(row):
    """Exact predicate chain on one glued pair; returns (ok, reason)."""
    N, T1, T2 = row["N"], row["T1"], row["T2"]
    if det2(N) == 0:
        return False, "singular glue (this port covers det N != 0 only)"
    A, At = T1, T2                      # invertible glue: saturation trivial
    for G in (A, At):
        if G != tp(G):
            return False, "Gram not symmetric"
        if any(G[i][i] % 2 for i in range(2)):
            return False, "Gram not even"
        if not (G[0][0] > 0 and det2(G) > 0):
            return False, "Gram not positive definite"
    B = mm(N, T2)
    M = mm([[Fr(x) for x in r] for r in B], inv2(At))
    Bt = mm(mm(inv2(At), tp(B)), [[Fr(x) for x in r] for r in A])
    if not all(x.denominator == 1 for r in M for x in r):
        return False, "M not integral (dual-integrality fails)"
    if not all(x.denominator == 1 for r in Bt for x in r):
        return False, "Bt not integral (dual-integrality fails)"
    S0 = mm(B, [[int(x) for x in r] for r in Bt])
    t = S0[0][0] + S0[1][1]
    if t % 2:
        return False, "odd trace on even lattice"
    if t != row["nmin_lat"]:
        return False, "trace != recorded nmin_lat (%d != %d)" % (
            t, row["nmin_lat"])
    if t != 2 * row["nflux_min"]:
        return False, "trace != 2 * nflux_min"
    Mi = [[int(x) for x in r] for r in M]
    chi0 = det2(Mi) ** 2 * abs(det2(A) * det2(At))
    if abs(det2(S0)) != chi0:
        return False, "det-chain identity broken"
    return True, ""


# E1 — the shipped 44-row known-answer table verifies exactly
rows = json.load(open(os.path.join(HERE, "bank", "sound_minvecs.json")))["rows"]
bad = [(r["d1"], r["k1"], r["d2"], r["k2"], why)
       for r in rows for ok, why in [check_row(r)] if not ok]
gate("E1.table_44_rows_exact", len(rows) == 44 and not bad,
     bad[:2] or "44/44 rows: dual-integrality + trace + det-chain exact")

# E2 — closed-form synthetic case, derived by hand:
#   A2 x A2 with N = [[0,1],[-1,0]]:  B = N*T2 = [[1,2],[-2,-1]],
#   Bt = N^T*T1 = [[-1,-2],[2,1]],  S0 = B*Bt = [[3,0],[0,3]].
toy = {"T1": [[2, 1], [1, 2]], "T2": [[2, 1], [1, 2]],
       "N": [[0, 1], [-1, 0]], "nmin_lat": 6, "nflux_min": 3}
ok, why = check_row(toy)
B = mm(toy["N"], toy["T2"])
Bt = mm(mm(inv2(toy["T2"]), tp(B)), [[Fr(x) for x in r] for r in toy["T1"]])
S0 = mm(B, [[int(x) for x in r] for r in Bt])
gate("E2.closed_form_a2xa2", ok and S0 == [[3, 0], [0, 3]],
     why or "S0 == 3*I, trace 6, det-chain 9 == 9")

# E3 — mutation controls MUST be caught
mut1 = json.loads(json.dumps(rows[0]))
mut1["N"][0][1] += 1                       # tampered glue entry
mut2 = json.loads(json.dumps(rows[0]))
mut2["nmin_lat"] += 2                      # tampered recorded value
c1, _ = check_row(mut1)
c2, _ = check_row(mut2)
gate("E3.mutations_caught", (not c1) and (not c2),
     "tampered glue + tampered nmin_lat both rejected")

npass = sum(1 for _n, okk in GATES if okk)
print("selftest_emit: %s (%d/%d)" % ("PASS" if npass == len(GATES) else "FAIL",
                                     npass, len(GATES)))
sys.exit(0 if npass == len(GATES) else 1)

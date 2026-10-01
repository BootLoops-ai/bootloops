#!/usr/bin/env python3
"""b3_oracle.py — B3 full-M direct-lattice oracle (independent algorithm):
brute enumeration of ALL C(M+6,6) support points |k|<=M with incremental arb
term products (prec 200), at the worst Stage-A corner (19..29)/1024 (real).
Computes jet-0 tower T_0 and the rho_1 jet T_{rho1} = sum c(k)phi^k 2(H_n-H_k1)
and gates them against the separable-GF pass balls (difference must contain 0).
This is the B3 pricing pass (3.2e7 pts at M=50); prints measured pts/s
vs the 14k pts/s p8d reference rate. Timed pilot: M=20 probe + extrapolation printed first;
checkpoint line per outer k1 shell (2h law). Usage: b3_oracle.py [M]"""
import sys, time
from fractions import Fraction as Fr
from math import comb, factorial
from flint import arb, ctx
import towers6 as tw
from towers6 import CQ

ctx.prec = 200
C0 = (11, 13, 15, 17, 19, 21)
PHI = [Fr(c + 8, 1024) for c in C0]


def tables(M):
    A = []
    for I in range(6):
        col, pw = [], Fr(1)
        for b in range(M + 1):
            q = pw / factorial(b) ** 2
            col.append(arb(q.numerator) / q.denominator)
            pw *= PHI[I]
        A.append(col)
    w2 = [arb(factorial(n) ** 2) for n in range(M + 1)]
    Ht = tw.harmonic_tables(M + 1)
    D = [[arb((2 * (Ht[1][n] - Ht[1][k])).numerator) /
          (2 * (Ht[1][n] - Ht[1][k])).denominator if n != k else arb(0)
          for k in range(M + 1)] for n in range(M + 1)]
    return A, w2, D


def lattice_pass(M, A, w2, D, t0):
    S0, S1, npts = arb(0), arb(0), 0
    for k1 in range(M + 1):
        p1 = A[0][k1]
        for k2 in range(M + 1 - k1):
            p2 = p1 * A[1][k2]
            n2 = k1 + k2
            for k3 in range(M + 1 - n2):
                p3 = p2 * A[2][k3]
                n3 = n2 + k3
                for k4 in range(M + 1 - n3):
                    p4 = p3 * A[3][k4]
                    n4 = n3 + k4
                    for k5 in range(M + 1 - n4):
                        p5 = p4 * A[4][k5]
                        n5 = n4 + k5
                        A6, DK = A[5], D
                        for k6 in range(M + 1 - n5):
                            t = p5 * A6[k6] * w2[n5 + k6]
                            S0 += t
                            S1 += t * DK[n5 + k6][k1]
                            npts += 1
        print(f"[ckpt] k1={k1} pts={npts} elapsed={time.time()-t0:.0f}s",
              flush=True)
    return S0, S1, npts


def main():
    M = int(sys.argv[1]) if len(sys.argv) > 1 else 50
    t0 = time.time()
    # timed-pilot probe at M=20, extrapolate by support-point count
    A, w2, D = tables(20)
    tp = time.time()
    _s0, _s1, np20 = lattice_pass(20, A, w2, D, tp)
    dt20 = time.time() - tp
    npM = comb(M + 6, 6)
    est = dt20 * npM / np20
    print(f"[pilot-probe] M=20: {np20} pts in {dt20:.1f}s "
          f"({np20/dt20:.0f} pts/s) -> M={M}: {npM} pts ~ {est/60:.1f} min",
          flush=True)
    assert est < 7200, "timed-pilot: extrapolates past 2h envelope, FREEZE"
    A, w2, D = tables(M)
    t1 = time.time()
    S0, S1, npts = lattice_pass(M, A, w2, D, t1)
    dt = time.time() - t1
    print(f"[rate] M={M}: {npts} pts in {dt:.1f}s = {npts/dt:.0f} pts/s "
          f"(anchor 14k pts/s p8d exact oracle)", flush=True)
    # gate vs separable GF at same M (independent algorithm, same ring)
    Ht = tw.harmonic_tables(M + 1)
    phis = [CQ(p) for p in PHI]
    import b3_pass as bp
    T = bp.eval_point(phis, True, M, Ht, "rho")
    for mono, S in ((0, S0), (tw.PB[0], S1)):
        diff = T[mono] - S
        assert diff.contains(arb(0)), f"oracle-vs-GF FAIL mono {mono}: {diff}"
        print(f"[gate] mono {mono}: |GF-oracle| ball contains 0 "
              f"(rad {float(diff.rad()):.2e}) PASS", flush=True)
    print(f"[done] wall {time.time()-t0:.1f}s  S0={S0}", flush=True)


if __name__ == "__main__":
    main()

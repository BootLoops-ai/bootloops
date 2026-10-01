#!/usr/bin/env python3
"""17-component analytic-solution jets on the curve (PA-2 oracle).
Y_j(s) = (theta^{B17[j]} w0)|curve = sum_cone c(n) prod_a n_a^{alpha_a} s^m.
Same cone walk as the reference series generator (verified EXACT against
geff_series.py); weights added per B17 monomial.  Emits jets_<M>.json (str ints).
usage: jets_series.py M"""
import json, sys, time
from math import factorial

B17 = [(0,0,0,0,0),(1,0,0,0,0),(0,1,0,0,0),(0,0,1,0,0),(0,0,0,1,0),
       (0,0,0,0,1),(1,1,0,0,0),(1,0,1,0,0),(0,1,1,0,0),(0,0,0,2,0),
       (0,0,0,1,1),(1,1,1,0,0),(2,0,0,0,0),(3,0,0,0,0),(1,0,0,1,0),
       (2,1,0,0,0),(2,0,1,0,0)]
M = int(sys.argv[1])
t0 = time.time()
Y = [[0] * (M + 1) for _ in range(17)]
F = [factorial(k) for k in range(4 * (M // 10) + 5)]
npts = 0
for n1 in range(M // 9 + 1):
    r1 = M - 9 * n1
    for n2 in range(n1, r1 // 10 + 1):
        r2 = r1 - 10 * n2
        for n3 in range(n1, r2 // 10 + 1):
            r3 = r2 - 10 * n3
            for n4 in range(max(n2, n3), r3 // 10 + 1):
                r4 = r3 - 10 * n4
                lo5 = max(n4, n2 + n3 - n1)
                hi5 = min(n2 + n3, n1 + n4, r4 // 10)
                for n5 in range(lo5, hi5 + 1):
                    npts += 1
                    m = 9 * n1 + 10 * (n2 + n3 + n4 + n5)
                    c = F[4 * n4] // (
                        F[2 * n4] * F[n1 + n4 - n5] * F[n5 - n4] *
                        F[n2 - n1] * F[n4 - n2] * F[n3 - n1] *
                        F[n4 - n3] * F[n1 - n2 - n3 + n5] * F[n2 + n3 - n5])
                    nn = (n1, n2, n3, n4, n5)
                    Y[0][m] += c
                    if n1:
                        Y[1][m] += c * n1
                        cn1 = c * n1
                        Y[12][m] += cn1 * n1
                        Y[13][m] += cn1 * n1 * n1
                        if n2:
                            Y[6][m] += cn1 * n2
                            Y[15][m] += cn1 * n1 * n2
                        if n3:
                            Y[7][m] += cn1 * n3
                            Y[16][m] += cn1 * n1 * n3
                        if n4:
                            Y[14][m] += cn1 * n4
                        if n2 and n3:
                            Y[11][m] += cn1 * n2 * n3
                    if n2:
                        Y[2][m] += c * n2
                        if n3:
                            Y[8][m] += c * n2 * n3
                    if n3:
                        Y[3][m] += c * n3
                    if n4:
                        Y[4][m] += c * n4
                        Y[9][m] += c * n4 * n4
                        if n5:
                            Y[10][m] += c * n4 * n5
                    if n5:
                        Y[5][m] += c * n5
print(f"M={M}: {npts} pts ({time.time()-t0:.0f}s)", flush=True)
import os
HERE = os.path.dirname(os.path.abspath(__file__))
json.dump({"M": M, "B17": [list(b) for b in B17],
           "Y": [[str(x) for x in row] for row in Y]},
          open(os.path.join(HERE, f"jets_{M}.json"), "w"))
print("saved jets")

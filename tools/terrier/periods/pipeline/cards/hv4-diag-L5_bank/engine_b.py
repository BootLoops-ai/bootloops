#!/usr/bin/env python3
# Part of the hv4-diag-L5 card battery: writes SERIES_B.json (cross-check verdicts vs SERIES_A.json), which is sha-pinned in regression_manifest.json and read by battery_l5.py; edits that change the output break the pin.
# engine_b.py (DIAG-SLICE) — Engine B: INDEPENDENT route. A3 by direct
# composition enumeration (math.comb; no shared code with Engine A), then
# the 3+3 multinomial split c6(n) = sum_j C(n,j)^2 A3(j) A3(n-j) and the
# 3+1 split A4(n) = sum_j C(n,j)^2 A3(j). Engine C anchor: direct 6-letter
# composition enumeration to n=40. Full-range cross-check vs SERIES_A.json.
import json, time
from math import comb
N = 400
t0 = time.time()
A3 = []
for n in range(N + 1):
    s = 0
    for a in range(n + 1):
        ca = comb(n, a)
        for b in range(n - a + 1):
            m = ca * comb(n - a, b)
            s += m * m
    A3.append(s)
c6 = [sum(comb(n,j)**2 * A3[j] * A3[n-j] for j in range(n+1)) for n in range(N+1)]
A4 = [sum(comb(n,j)**2 * A3[j] for j in range(n+1)) for n in range(N+1)]
NC = 40; cC = []
for n in range(NC + 1):
    s = 0
    for a in range(n + 1):
     for b in range(n - a + 1):
      for c in range(n - a - b + 1):
       for d in range(n - a - b - c + 1):
        for e in range(n - a - b - c - d + 1):
            m = (comb(n,a)*comb(n-a,b)*comb(n-a-b,c)
                 *comb(n-a-b-c,d)*comb(n-a-b-c-d,e))
            s += m * m
    cC.append(s)
SA = json.load(open("SERIES_A.json"))
rec = {"N": N, "agree_c6_full": [str(x) for x in c6] == SA["c6"],
       "agree_domb_full": [str(x) for x in A4] == SA["c4_domb"],
       "agree_c3_full": [str(x) for x in A3] == SA["c3"],
       "agree_direct6_n40": cC == c6[:NC+1],
       "elapsed_s": round(time.time() - t0, 3),
       "engine": "B: direct A3 enum + 3+3/3+1 splits; C: direct 6-comp n<=40"}
json.dump(rec, open("SERIES_B.json", "w"))
print("ENGINE-B", rec)
assert all(v for k, v in rec.items() if k.startswith("agree"))

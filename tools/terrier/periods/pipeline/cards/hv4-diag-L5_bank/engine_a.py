#!/usr/bin/env python3
# Part of the hv4-diag-L5 card battery: writes SERIES_A.json, which is sha-pinned in regression_manifest.json and read by engine_b.py, battery_l5.py and mkcard.py; edits that change the output break the pin.
# engine_a.py (DIAG-SLICE) — Engine A: k-letter abelian-square counts by
# letter adjunction c^{(k+1)}(n) = sum_j C(n,j)^2 c^{(k)}(j)  (multinomial
# factorization; the add-one-variable torus-period construction). Own
# Pascal triangle (Engine B uses math.comb).
# c^{(6)} = diagonal HV4 fundamental period coeffs (eq 5.18 restricted);
# c^{(4)} = Domb / K3-diagonal series (tex eq K3pi0diag). Exact ints.
import json, time
N = 400
t0 = time.time()
C = [[1]]
for n in range(1, N + 1):
    C.append([1] + [C[-1][j-1] + C[-1][j] for j in range(1, n)] + [1])
def adjoin(a):
    return [sum(C[n][j]*C[n][j]*a[j] for j in range(n+1)) for n in range(N+1)]
ks = {1: [1]*(N+1)}
for k in range(2, 7):
    ks[k] = adjoin(ks[k-1])
assert ks[4][:5] == [1, 4, 28, 256, 2716]      # Domb anchor (printed in tex)
assert ks[6][:3] == [1, 6, 66]                 # 6-letter hand anchor
out = {"N": N, "c6": [str(x) for x in ks[6]],
       "c4_domb": [str(x) for x in ks[4]], "c3": [str(x) for x in ks[3]],
       "elapsed_s": round(time.time()-t0, 3),
       "engine": "A: letter-adjunction convolution, own Pascal triangle"}
json.dump(out, open("SERIES_A.json", "w"))
print("ENGINE-A OK", out["elapsed_s"], "s; digits(c6[400]) =", len(out["c6"][400]))

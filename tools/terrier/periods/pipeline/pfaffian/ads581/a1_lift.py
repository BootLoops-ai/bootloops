#!/usr/bin/env python3
"""A-1 CRT/Wang lift: interp_<P>.npz (>=2 primes) -> A_exact.json.
Per-coeff CRT then Wang rational reconstruction; stability gate = lift
from all-but-last prime must agree on a 5% subsample (design A-1 (i)).
usage: a1_lift.py P1 P2 [P3 ...]"""
import sys, os, json, time, random
import numpy as np
from fractions import Fraction as Fr
from math import isqrt, gcd

HERE = os.path.dirname(os.path.abspath(__file__))
primes = [int(x) for x in sys.argv[1:]]
assert len(primes) >= 2
banks = [np.load(os.path.join(HERE, f"interp_{p}.npz")) for p in primes]
for b in banks[1:]:
    assert np.array_equal(b["degs"], banks[0]["degs"]), "deg table mismatch"
    assert int(b["dden"]) == int(banks[0]["dden"])
degs = banks[0]["degs"]; dden = int(banks[0]["dden"]); RK = degs.shape[0]

def crt_pair(a1, m1, a2, m2):
    d = pow(m1, -1, m2)
    return (a1 + ((a2 - a1) * d % m2) * m1) % (m1 * m2), m1 * m2

def wang(r, M):
    B = isqrt(M // 2)
    r0, r1 = M, r % M
    t0, t1 = 0, 1
    while r1 > B:
        q = r0 // r1
        r0, r1, t0, t1 = r1, r0 - q * r1, t1, t0 - q * t1
    if abs(t1) > B or t1 == 0 or gcd(r1, abs(t1)) != 1 or gcd(t1, M) != 1:
        return None
    return Fr(r1, t1) if t1 > 0 else Fr(-r1, -t1)

def lift_arr(vals_by_prime, plist):
    """vals_by_prime: list of int lists (same length). -> (fracs, nfail)."""
    n = len(vals_by_prime[0])
    out, nfail = [], 0
    for k in range(n):
        a, m = int(vals_by_prime[0][k]) % plist[0], plist[0]
        for t in range(1, len(plist)):
            a, m = crt_pair(a, m, int(vals_by_prime[t][k]) % plist[t],
                            plist[t])
        f = wang(a, m)
        if f is None:
            nfail += 1
            out.append(None)
        else:
            out.append(f)
    return out, nfail

t0 = time.time()
DENf, dfail = lift_arr([b["DEN"] for b in banks], primes)
entries, nfail_tot, ncoef = {}, dfail, len(DENf)
hmax = max((abs(f.numerator).bit_length() + f.denominator.bit_length()
            for f in DENf if f is not None), default=0)
for i in range(RK):
    for j in range(RK):
        d = int(degs[i, j])
        if d < 0:
            continue
        fr, nf = lift_arr([b["Anum"][i, j, :d + 1] for b in banks], primes)
        nfail_tot += nf; ncoef += d + 1
        if nf == 0:
            hmax = max(hmax, max(abs(f.numerator).bit_length()
                                 + f.denominator.bit_length() for f in fr))
        entries[f"{i},{j}"] = fr
print(f"[lift k={len(primes)}] {ncoef} coeffs, {nfail_tot} Wang fails, "
      f"max height {hmax} bits ({time.time()-t0:.0f}s)", flush=True)
if nfail_tot:
    print("LIFT INCOMPLETE — need more primes"); sys.exit(2)

# stability gate: all-but-last prime lift on 5% subsample must agree
random.seed(20260711)
sub = random.sample([(i, j, k) for i in range(RK) for j in range(RK)
                     if degs[i, j] >= 0 for k in range(int(degs[i, j]) + 1)],
                    max(1, ncoef // 20))
plist2 = primes[:-1]
mism = 0
if len(plist2) >= 2:
    for (i, j, k) in sub:
        a, m = None, None
        for t, p in enumerate(plist2):
            v = int(banks[t]["Anum"][i, j, k]) % p
            a, m = (v, p) if a is None else crt_pair(a, m, v, p)
        f = wang(a, m)
        if f != entries[f"{i},{j}"][k]:
            mism += 1
    print(f"[stability] k-1 subsample: {mism}/{len(sub)} mismatch "
          f"-> {'PASS' if mism == 0 else 'FAIL'}", flush=True)
else:
    print("[stability] k-1 == 1 prime: SKIPPED (need k>=3)"); mism = -1

out = {"primes": primes, "dden": dden, "rank": RK,
       "den": [str(f) for f in DENf],
       "deg_table": degs.tolist(),
       "entries": {k: [str(f) for f in v] for k, v in entries.items()},
       "max_height_bits": hmax, "stability_pass": (mism == 0)}
json.dump(out, open(os.path.join(HERE, "A_exact.json"), "w"))
json.dump({"dden": dden, "deg_table": degs.tolist(),
           "max_num_deg": int(degs.max())},
          open(os.path.join(HERE, "deg_table.json"), "w"))
print(f"[saved] A_exact.json ({os.path.getsize(os.path.join(HERE, 'A_exact.json'))//1024} KB)")

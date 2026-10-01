#!/usr/bin/env python3
"""PA-1: exact A(s) reproduces ALL 289 entries at >=256 nodes of a FRESH
prime never used in the lift.  usage: gate_pa1.py P_fresh"""
import sys, os, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
# reference-data root from the environment (no machine-local default;
# fail-closed). Layout: $TERRIER_REFERENCE_BANKS/cards/ads-5-81-3213_bank/
# with the flux frame, the s^3000 series and dmodule/ (pmodp, boxops,
# connection samples) -- not included in the package, see pfaffian/README.md
_CB = os.environ.get("TERRIER_REFERENCE_BANKS", "")
if not _CB:
    raise SystemExit("REFUSE: TERRIER_REFERENCE_BANKS is not set — point it "
                     "at the reference-data root holding "
                     "cards/ads-5-81-3213_bank/ (not included in the "
                     "package; see periods/pipeline/pfaffian/README.md)")
DM = os.path.join(_CB, "cards", "ads-5-81-3213_bank", "dmodule")
P = int(sys.argv[1])
AX = json.load(open(os.path.join(HERE, "A_exact.json")))
assert P not in AX["primes"], "fresh prime was used in the lift!"
from fractions import Fraction as Fr

def red(fs):  # exact rational coeff list -> mod-p int64 array
    out = np.zeros(len(fs), dtype=np.int64)
    for k, s in enumerate(fs):
        f = Fr(s)
        out[k] = f.numerator % P * pow(f.denominator % P, P - 2, P) % P
    return out

D = np.load(os.path.join(DM, f"conn_samples_{P}.npz"))
sig, A = D["sig"].astype(np.int64), D["A"].astype(np.int64)
N = len(sig); RK = A.shape[0]
denp = red(AX["den"])

def pev(pl, xs):
    acc = np.zeros(len(xs), dtype=np.int64)
    for cf in pl[::-1]:
        acc = (acc * xs + int(cf)) % P
    return acc

dv = pev(denp, sig)
assert np.all(dv % P != 0), "den vanishes at a fresh node"
dinv = np.array([pow(int(v), P - 2, P) for v in dv], dtype=np.int64)
ndiff = 0
for i in range(RK):
    for j in range(RK):
        key = f"{i},{j}"
        if key in AX["entries"]:
            v = pev(red(AX["entries"][key]), sig) * dinv % P
        else:
            v = np.zeros(N, dtype=np.int64)
        ndiff += int(np.count_nonzero((v - A[i, j]) % P))
tot = RK * RK * N
res = {"gate": "PA-1", "fresh_prime": P, "nodes": N,
       "entry_values_checked": tot, "diffs": ndiff,
       "pass": ndiff == 0}
json.dump(res, open(os.path.join(HERE, "gate_PA1.json"), "w"), indent=1)
print(f"[PA-1] fresh p={P}: {ndiff}/{tot} diffs over {N} nodes -> "
      f"{'PASS' if ndiff == 0 else 'FAIL'}")
sys.exit(0 if ndiff == 0 else 1)

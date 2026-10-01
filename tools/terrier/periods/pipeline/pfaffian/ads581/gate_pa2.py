#!/usr/bin/env python3
"""PA-2: exact A(s) annihilates the held-out analytic jets over Q:
DEN * theta Y_j - sum_i Num[i][j] Y_i == 0 through s^M EXACTLY (flint).
Jets Y from jets_series.py (held out: never used in reconstruction).
Component 0 == the reference s^3000 fundamental-period series, re-asserted
here (that series is read from the reference-data root)."""
import os, json, gzip, time
from fractions import Fraction as Fr
import flint

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
BANK = os.path.join(_CB, "cards", "ads-5-81-3213_bank")
AX = json.load(open(os.path.join(HERE, "A_exact.json")))
J = json.load(open(os.path.join(HERE, "jets_3000.json")))
M = J["M"]; RK = AX["rank"]
W = json.load(gzip.open(os.path.join(BANK, "w0_series_3000.json.gz")))
assert all(J["Y"][0][m] == W["a"][m] for m in range(M + 1)), \
    "jets comp-0 != reference w0 series"

def qpoly(strs):
    return flint.fmpq_poly([flint.fmpq(Fr(s).numerator, Fr(s).denominator)
                            for s in strs])

t0 = time.time()
DEN = qpoly(AX["den"])
Y = [flint.fmpq_poly([int(v) for v in J["Y"][i]]) for i in range(RK)]
thY = [flint.fmpq_poly([m * int(J["Y"][i][m]) for m in range(M + 1)])
       for i in range(RK)]
NUM = {}
for k, v in AX["entries"].items():
    i, j = map(int, k.split(","))
    NUM[(i, j)] = qpoly(v)
worst_ok = M
allzero = True
for j in range(RK):
    R = DEN * thY[j]
    for i in range(RK):
        if (i, j) in NUM:
            R -= NUM[(i, j)] * Y[i]
    cs = R.coeffs()
    bad = [m for m in range(min(len(cs), M + 1)) if cs[m] != 0]
    if bad:
        allzero = False
        worst_ok = min(worst_ok, bad[0] - 1)
        print(f"[PA-2] row {j}: first nonzero at m={bad[0]}")
res = {"gate": "PA-2", "M": M, "rows": RK,
       "checked_through": M, "pass": allzero,
       "note": "DEN*thetaY - Num^T-contraction == 0 through s^M exactly "
               "(Y complete to M, so products exact to M)"}
json.dump(res, open(os.path.join(HERE, "gate_PA2.json"), "w"), indent=1)
print(f"[PA-2] all {RK} rows zero through s^{M}: "
      f"{'PASS' if allzero else 'FAIL'} ({time.time()-t0:.0f}s)")

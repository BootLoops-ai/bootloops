#!/usr/bin/env python3
"""A-1 per-prime step: conn_samples_<P>.npz -> monic-den normalized
(Anum, DEN) mod p via ratrecon2 per entry (prod2 interp block, reused).
Saves interp_<P>.npz: Anum (17,17,<=575) padded, DEN, per-entry degs.
usage: a1_perprime.py P [NVAL]"""
import sys, os, time
import numpy as np
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
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, DM)
import pmodp

P = int(sys.argv[1])
NVAL = int(sys.argv[2]) if len(sys.argv) > 2 else 64
T = pmodp.make(P)
D = np.load(os.path.join(DM, f"conn_samples_{P}.npz"))
sig, A = D["sig"].astype(np.int64), D["A"].astype(np.int64)
N = len(sig); RK = A.shape[0]
x, xv = sig[:N - NVAL], sig[N - NVAL:]
n = len(x)
t0 = time.time()
G = T["prodtree"](x)
num = [[None] * RK for _ in range(RK)]
den = [[None] * RK for _ in range(RK)]
for i in range(RK):
    for j in range(RK):
        y = A[i, j]
        if np.all(y[:n] % P == 0):
            num[i][j] = np.array([0], dtype=np.int64)
            den[i][j] = np.array([1], dtype=np.int64)
            continue
        r, s_ = T["ratrecon2"](y[:n], x, y[N - NVAL:], xv, G)
        assert r is not None, f"A[{i}][{j}] interp fails at p={P}"
        c = pow(int(s_[-1]), P - 2, P)
        num[i][j] = r * c % P
        den[i][j] = s_ * c % P
print(f"[{P}] all {RK*RK} entries interp ({time.time()-t0:.0f}s)", flush=True)
lcm = np.array([1], dtype=np.int64)
for i in range(RK):
    for j in range(RK):
        g = T["pgcd"](lcm.copy(), den[i][j].copy())
        q, _ = T["pdivmod"](den[i][j], g)
        lcm = T["pmulc"](lcm, q)
DEN = lcm  # monic by construction
dden = len(DEN) - 1
Anum = np.zeros((RK, RK, 704), dtype=np.int64)
degs = np.zeros((RK, RK), dtype=np.int64) - 1
for i in range(RK):
    for j in range(RK):
        q, rem = T["pdivmod"](DEN, den[i][j])
        assert len(rem) == 1 and rem[0] == 0
        nm = T["pmulc"](num[i][j], q) if num[i][j].any() else \
            np.array([0], dtype=np.int64)
        assert len(nm) <= 704
        if nm.any():
            Anum[i, j, :len(nm)] = nm
            degs[i, j] = len(nm) - 1
np.savez(os.path.join(HERE, f"interp_{P}.npz"), Anum=Anum, DEN=DEN,
         degs=degs, p=P, dden=dden)
print(f"[{P}] den deg {dden}, max num deg {int(degs.max())} "
      f"({time.time()-t0:.0f}s) -> interp_{P}.npz", flush=True)

# rank_ladder.py — certified-transport cost vs holonomic rank, algorithm='naive'.
# Run under Sage's Python (with ore_algebra installed), one rank per process:
#   RL_RANK=24 RL_STAGES=ord20,ord40,sing20,sing40 sage -python rank_ladder.py
#
# Test operator (synthetic hypergeometric family — no external data):
#   L_R = 2^R * theta^R - x * (2*theta+1)^R,   theta = x*Dx
# = the pFq operator of f_R(x) = RF_{R-1}(1/2,...,1/2; 1,...,1; x):
#   rank R, regular singular points exactly {0, 1, oo},
#   leading coefficient 2^R x^R (1-x),
#   series c_n = ((n-1/2)^R / n^R) c_{n-1}, c_0 = 1  (known solution).
# Built directly from Stirling numbers (theta^m = sum_k S(m,k) x^k Dx^k):
#   Dx^k coefficient p_k = x^k * (2^R S(R,k) - W_k x),
#   W_k = sum_{m>=k} binom(R,m) 2^m S(m,k).
#
# Paths (fixed representatives across the whole ladder):
#   ordinary:  [1/3, 1/2]                       (min dist to sing = 1/2)
#   near-sing: [1/2, 1 + I/32, 3/2]             (passes x=1 at dist 1/32)
# Precisions: eps = 1e-20 and 1e-40.
import json
import os
import sys
import time

from sage.all import (QQ, ZZ, PolynomialRing, RealBallField, QQbar, I,
                      binomial)
from sage.version import version as sage_version
from ore_algebra import OreAlgebra

R = int(os.environ.get("RL_RANK", "24"))
STAGES = os.environ.get("RL_STAGES", "ord20,ord40,sing20,sing40").split(",")
OUT = os.environ.get("RL_OUT", ".")

Rx = PolynomialRing(QQ, "x")
x = Rx.gen()
A = OreAlgebra(Rx, "Dx")
Dx = A.gen()
RBF = RealBallField(430)

rep = {"rank": R, "sage": str(sage_version), "algorithm": "naive",
       "stages": {}}


def stamp():
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def maxrss_mb():
    try:
        import resource
    except ImportError:          # non-POSIX: report nothing rather than lie
        return 0.0
    ru = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    # ru_maxrss is BYTES on macOS, KB on Linux
    return ru / (1024.0 * 1024.0) if sys.platform == "darwin" else ru / 1024.0


def ball_stats(M):
    """max |mid| and max radius over the entries of a ball matrix."""
    mx, mr = 0.0, 0.0
    for i in range(M.nrows()):
        for j in range(M.ncols()):
            z = M[i, j]
            try:
                m = z.abs().mid()
                r = z.rad()
            except AttributeError:
                m = abs(z.mid())
                r = max(z.real().rad(), z.imag().rad())
            mx = max(mx, float(m))
            mr = max(mr, float(r))
    return mx, mr


# ---------------- build L_R from the Stirling triangle ------------------------
t0 = time.time()
S = [ZZ(1)]                      # row m of S(m,k), k=0..m; start m=0
W = [ZZ(0)] * (R + 1)
for m in range(1, R + 1):
    Snew = [ZZ(0)] * (m + 1)
    Snew[0] = ZZ(0)
    for k in range(1, m):
        Snew[k] = k * S[k] + S[k - 1]
    Snew[m] = ZZ(1)
    if m < R:
        c = binomial(R, m) * ZZ(2) ** m
        for k in range(0, m + 1):
            W[k] += c * Snew[k]
    S = Snew
W[0] += ZZ(1)                    # m=0 term: binom(R,0) 2^0 S(0,0)
cR = ZZ(2) ** R                  # m=R term of W, and the 2^R prefactor
for k in range(0, R + 1):
    W[k] += cR * S[k]
coeffs = []
for k in range(0, R + 1):
    Sk = S[k] if k <= R else ZZ(0)
    coeffs.append(x ** k * (cR * Sk - W[k] * x))
L = A(coeffs)
t_build = time.time() - t0
hmax = max(max(ZZ(c.numerator()).nbits() for c in p.coefficients())
           for p in coeffs if p != 0)
rep["build"] = {"seconds": round(t_build, 3),
                "order": int(L.order()),
                "degree": int(max(p.degree() for p in coeffs)),
                "coeff_height_bits": int(hmax),
                "maxrss_mb_after": round(maxrss_mb(), 1)}
assert L.order() == R
lead = L.leading_coefficient()
assert lead == cR * x ** R * (1 - x), "leading coefficient sanity failed"
print(f"[{stamp()}] rank {R}: built (order {R}, deg {rep['build']['degree']}, "
      f"height {hmax} bits, {t_build:.2f}s)")

PATHS = {"ord": [QQ(1) / 3, QQ(1) / 2],
         "sing": [QQ(1) / 2, 1 + QQbar(I) / 32, QQ(3) / 2]}
EPS = {"20": 1e-20, "40": 1e-40, "60": 1e-60, "80": 1e-80}

for st in STAGES:
    st = st.strip()
    if st == "validate":
        continue
    pkey = "ord" if st.startswith("ord") else "sing"
    ekey = st[-2:]
    t0 = time.time()
    M = L.numerical_transition_matrix(PATHS[pkey], EPS[ekey],
                                      algorithm='naive')
    dt = time.time() - t0
    mx, mr = ball_stats(M)
    rep["stages"][st] = {"path": pkey, "eps": EPS[ekey],
                         "wall_s": round(dt, 2),
                         "max_entry_abs": mx, "max_entry_rad": mr,
                         "maxrss_mb_after": round(maxrss_mb(), 1)}
    print(f"[{stamp()}] rank {R} {st}: {dt:.2f}s  max|entry|={mx:.3e}  "
          f"max rad={mr:.3e}  maxrss={rep['stages'][st]['maxrss_mb_after']}MB")
    json.dump(rep, open(f"{OUT}/ladder_R{R}.json", "w"), indent=1,
              default=str)

# ---------------- validation stage (small ranks only) -------------------------
if "validate" in [s.strip() for s in STAGES]:
    # (i) IC convention determination on Dx^3 (solutions 1, x, x^2)
    x0, x1 = QQ(1) / 3, QQ(1) / 2
    Mc = (Dx ** 3).numerical_transition_matrix([x0, x1], 1e-30,
                                               algorithm='naive')
    vder = [RBF(x0 ** 2), RBF(2 * x0), RBF(2)]
    vtay = [RBF(x0 ** 2), RBF(2 * x0), RBF(1)]
    tgt = RBF(x1 ** 2)
    od = sum(Mc[0, j] * vder[j] for j in range(3))
    ot = sum(Mc[0, j] * vtay[j] for j in range(3))
    conv = ("derivatives" if abs((od.real() - tgt).mid()) < 1e-20 else
            "taylor" if abs((ot.real() - tgt).mid()) < 1e-20 else None)
    assert conv is not None
    rep["validate"] = {"ic_convention": conv}

    # (ii) round trip: Id must lie INSIDE the ball matrix Mb*Mf (containment
    # is the certified statement; ball WIDTH degrades with rank — measured)
    Mf = L.numerical_transition_matrix([x0, x1], 1e-40, algorithm='naive')
    Mb = L.numerical_transition_matrix([x1, x0], 1e-40, algorithm='naive')
    P = Mb * Mf
    dev, ok = 0.0, True
    for i in range(R):
        for j in range(R):
            d = P[i, j] - (1 if i == j else 0)
            dev = max(dev, float(d.abs().upper()))
            ok = ok and d.contains_zero()
    rep["validate"]["roundtrip_contains_identity"] = bool(ok)
    rep["validate"]["roundtrip_ball_width"] = dev
    assert ok, "identity NOT contained in roundtrip ball"

    # (iii) known-solution gate: series of f_R at 0 -> ICs at 1/3 ->
    #       transport -> f_R(1/2) vs direct series
    N = 260
    c = [QQ(1)]
    for n in range(1, N):
        c.append(c[-1] * (QQ(2 * n - 1) / QQ(2 * n)) ** R)

    def f_ders(xq, nder):
        xb = RBF(xq)
        outs = []
        for d in range(nder + 1):
            v = RBF(0)
            for n in range(d, N):
                term = RBF(c[n]) * ZZ(n).factorial() / ZZ(n - d).factorial()
                v += term * xb ** (n - d)
            # geometric tail bound: |c_n| <= 1, ratio <= x  (c_n decreasing)
            tail = RBF(c[N - 1]) * ZZ(N - 1).factorial() \
                / ZZ(N - 1 - d).factorial() * xb ** (N - 1 - d) \
                / (1 - xb) * 2
            outs.append(v.add_error(tail.abs().upper()))
        return outs
    ics = f_ders(x0, R - 1)
    if conv == "taylor":
        v0 = [ics[i] / ZZ(i).factorial() for i in range(R)]
    else:
        v0 = list(ics)
    ft = sum(Mf[0, j] * v0[j] for j in range(R))
    ft = ft.real() if hasattr(ft, "real") else ft
    fd = f_ders(x1, 0)[0]
    d = ft - fd
    rep["validate"]["known_solution_contains_zero"] = bool(d.contains_zero())
    rep["validate"]["known_solution_ball_width"] = float(d.abs().upper())
    rep["validate"]["known_solution_mid_dev"] = float(abs(d.mid()))
    assert d.contains_zero(), "transported value inconsistent with series"
    print(f"[{stamp()}] rank {R} validate: conv={conv}  roundtrip contains "
          f"Id (width {dev:.2e})  known-sol contains 0 "
          f"(width {rep['validate']['known_solution_ball_width']:.2e}, "
          f"mid dev {rep['validate']['known_solution_mid_dev']:.2e})  PASS")

rep["maxrss_mb_final"] = round(maxrss_mb(), 1)
rep["stamp"] = stamp()
json.dump(rep, open(f"{OUT}/ladder_R{R}.json", "w"), indent=1, default=str)
print(f"[{stamp()}] rank {R}: receipt -> {OUT}/ladder_R{R}.json")

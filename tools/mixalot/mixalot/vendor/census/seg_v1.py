#!/usr/bin/env python3
"""seg_v1: contiguous-segment author-mixture evidence for ordered count data.

Model and notation: see the docstrings below. Channel/count matrices: the
counts_*_{A,B}.json files (two independent extractions, A and B).

Modes:
  exact  — Python Fraction end-to-end (rational alpha), moderate sizes.
  float  — log-domain lgamma + logsumexp, full-vocabulary scale.

Public surface:
  prefix_sums(X)                       -> P (n+1 x k int lists)
  Z_exact(X, g, a) / Z_float(X, g, a)  -> (log)evidence, uniform seg prior
  dp_exact(X, g, a) / dp_float(X, g, a)-> raw sum_T prod W (no 1/C factor)
  enum_exact(X, g, a)                  -> same by composition enumeration
  boundary_posterior_exact/float       -> P(boundary after t | g, data)
  free2_exact(X, a) / free2_float(X,a) -> free-assignment g=2 evidence
  selfcheck()                          -> writes seg_v1_selfcheck.txt

Run selfcheck:  python3 seg_v1.py --selfcheck
"""
from fractions import Fraction
from math import lgamma, log, exp, comb, inf
from itertools import combinations, product
import json, sys, random

# ---------------------------------------------------------------- primitives

def rising_exact(a, m):
    """(a)_m = a (a+1) ... (a+m-1) as an exact Fraction."""
    out = Fraction(1)
    for r in range(m):
        out *= a + r
    return out

def W_exact(X, a):
    """Collapsed segment weight prod_c (a)_{X_c} / (k a)_N, exact."""
    k = len(X)
    N = sum(X)
    num = Fraction(1)
    for x in X:
        if x:
            num *= rising_exact(a, x)
    return num / rising_exact(k * a, N)

class LogTables:
    """T[x] = lgamma(a+x)-lgamma(a); U[N] = lgamma(ka)-lgamma(ka+N)."""
    def __init__(self, a, k, nmax):
        af = float(a)
        kaf = float(k * a)
        la, lka = lgamma(af), lgamma(kaf)
        self.T = [lgamma(af + x) - la for x in range(nmax + 1)]
        self.U = [lka - lgamma(kaf + N) for N in range(nmax + 1)]

def logW_float(X, tab):
    N = 0
    s = 0.0
    for x in X:
        N += x
        s += tab.T[x]
    return s + tab.U[N]

def logsumexp(vals):
    m = max(vals)
    if m == -inf:
        return -inf
    return m + log(sum(exp(v - m) for v in vals))

def prefix_sums(X):
    """X: n x k int matrix -> P: (n+1) x k with P[j] = sum of rows < j."""
    n = len(X)
    k = len(X[0])
    P = [[0] * k for _ in range(n + 1)]
    for i in range(n):
        Pi, Pi1, Xi = P[i], P[i + 1], X[i]
        for c in range(k):
            Pi1[c] = Pi[c] + Xi[c]
    return P

def seg_counts(P, i, j):
    """Pooled counts of units (i, j]."""
    return [b - a for a, b in zip(P[i], P[j])]

# ------------------------------------------------- contiguous evidence exact

def enum_exact(X, g, a):
    """sum over segmentations of prod_s W(segment s): direct composition
    enumeration via boundary positions (O(C(n-1,g-1)) terms)."""
    n = len(X)
    P = prefix_sums(X)
    tot = Fraction(0)
    for bnds in combinations(range(1, n), g - 1):
        cuts = (0,) + bnds + (n,)
        term = Fraction(1)
        for s in range(g):
            term *= W_exact(seg_counts(P, cuts[s], cuts[s + 1]), a)
        tot += term
    return tot

def dp_exact(X, g, a):
    """Same quantity by DP: F(j,m) = sum_i F(i,m-1) W((i,j])."""
    n = len(X)
    P = prefix_sums(X)
    W = {}
    for i in range(n):
        for j in range(i + 1, n + 1):
            W[(i, j)] = W_exact(seg_counts(P, i, j), a)
    F = [[Fraction(0)] * (g + 1) for _ in range(n + 1)]
    F[0][0] = Fraction(1)
    for m in range(1, g + 1):
        for j in range(m, n + 1):
            acc = Fraction(0)
            for i in range(m - 1, j):
                if F[i][m - 1]:
                    acc += F[i][m - 1] * W[(i, j)]
            F[j][m] = acc
    return F[n][g]

def Z_exact(X, g, a):
    """Evidence with the uniform prior over the C(n-1,g-1) segmentations."""
    return dp_exact(X, g, a) / comb(len(X) - 1, g - 1)

# ------------------------------------------------- contiguous evidence float

def _logW_table(X, a):
    n = len(X)
    P = prefix_sums(X)
    tab = LogTables(a, len(X[0]), sum(P[n]))
    LW = {}
    for i in range(n):
        for j in range(i + 1, n + 1):
            LW[(i, j)] = logW_float(seg_counts(P, i, j), tab)
    return LW

def dp_float(X, g, a, LW=None):
    """log of sum_T prod W via forward DP in log domain."""
    n = len(X)
    if LW is None:
        LW = _logW_table(X, a)
    F = [[-inf] * (g + 1) for _ in range(n + 1)]
    F[0][0] = 0.0
    for m in range(1, g + 1):
        for j in range(m, n + 1):
            vals = [F[i][m - 1] + LW[(i, j)] for i in range(m - 1, j)]
            F[j][m] = logsumexp(vals) if vals else -inf
    return F[n][g]

def Z_float(X, g, a):
    return dp_float(X, g, a) - log(comb(len(X) - 1, g - 1))

def _backward_float(X, g, a, LW):
    n = len(X)
    B = [[-inf] * (g + 1) for _ in range(n + 1)]
    B[n][0] = 0.0
    for m in range(1, g + 1):
        for i in range(n - m, -1, -1):
            vals = [LW[(i, j)] + B[j][m - 1] for j in range(i + 1, n - m + 2)]
            B[i][m] = logsumexp(vals) if vals else -inf
    return B

def boundary_posterior_float(X, g, a):
    """P(boundary after unit t | g, data), t = 1..n-1 (log-float mode)."""
    n = len(X)
    LW = _logW_table(X, a)
    Ff = [[-inf] * (g + 1) for _ in range(n + 1)]
    Ff[0][0] = 0.0
    for m in range(1, g + 1):
        for j in range(m, n + 1):
            vals = [Ff[i][m - 1] + LW[(i, j)] for i in range(m - 1, j)]
            Ff[j][m] = logsumexp(vals) if vals else -inf
    B = _backward_float(X, g, a, LW)
    Zlog = Ff[n][g]
    out = {}
    for t in range(1, n):
        vals = [Ff[t][m] + B[t][g - m] for m in range(1, g)]
        out[t] = exp(logsumexp(vals) - Zlog) if vals else 0.0
    return out

def boundary_posterior_exact(X, g, a):
    """Exact-Fraction boundary posterior (moderate sizes)."""
    n = len(X)
    P = prefix_sums(X)
    W = {(i, j): W_exact(seg_counts(P, i, j), a)
         for i in range(n) for j in range(i + 1, n + 1)}
    F = [[Fraction(0)] * (g + 1) for _ in range(n + 1)]
    F[0][0] = Fraction(1)
    for m in range(1, g + 1):
        for j in range(m, n + 1):
            F[j][m] = sum((F[i][m - 1] * W[(i, j)]
                           for i in range(m - 1, j)), Fraction(0))
    B = [[Fraction(0)] * (g + 1) for _ in range(n + 1)]
    B[n][0] = Fraction(1)
    for m in range(1, g + 1):
        for i in range(n - m, -1, -1):
            B[i][m] = sum((W[(i, j)] * B[j][m - 1]
                           for j in range(i + 1, n - m + 2)), Fraction(0))
    Z = F[n][g]
    return {t: sum((F[t][m] * B[t][g - m] for m in range(1, g)),
                   Fraction(0)) / Z for t in range(1, n)}

# ------------------------------------------------------ free-assignment g=2

def free2_direct_exact(X, a):
    """2^{-n} sum over all label vectors of W(group0) W(group1); naive
    per-assignment pooling. Reference implementation for the selfcheck."""
    n = len(X)
    k = len(X[0])
    tot = Fraction(0)
    for lab in product((0, 1), repeat=n):
        c0 = [0] * k
        c1 = [0] * k
        for i, l in enumerate(lab):
            tgt = c1 if l else c0
            for c in range(k):
                tgt[c] += X[i][c]
        tot += W_exact(c0, a) * W_exact(c1, a)
    return tot / (Fraction(2) ** n)

def free2_gray_exact(X, a):
    """Same sum in Gray-code order: one unit moves per step, O(2^n k)."""
    n = len(X)
    k = len(X[0])
    total = [sum(X[i][c] for i in range(n)) for c in range(k)]
    c1 = [0] * k
    in1 = [False] * n
    tot = W_exact(total, a) * W_exact(c1, a)     # mask 0
    for step in range(1, 2 ** n):
        u = (step & -step).bit_length() - 1      # Gray code flips bit u
        sgn = -1 if in1[u] else 1
        in1[u] = not in1[u]
        for c in range(k):
            c1[c] += sgn * X[u][c]
        c0 = [t - x for t, x in zip(total, c1)]
        tot += W_exact(c0, a) * W_exact(c1, a)
    return tot / (Fraction(2) ** n)

def free2_float(X, a, chunk_bits=16):
    """log Z_free by chunked numpy over 2^n assignments (n <= 25).
    Meet-in-the-middle pooled counts + lgamma lookup tables."""
    import numpy as np
    X = np.asarray(X, dtype=np.int64)
    n, k = X.shape
    h1 = min(n, chunk_bits)                      # low bits enumerated per chunk
    lo, hi = X[:h1], X[h1:]
    nlo, nhi = 2 ** h1, 2 ** (n - h1)
    cLo = np.zeros((nlo, k), dtype=np.int64)
    for u in range(h1):                          # doubling construction
        blk = 1 << u
        cLo[blk:2 * blk] = cLo[:blk] + X[u]
    Ntot = int(X.sum())
    tab = LogTables(a, k, Ntot)
    T = np.array(tab.T)
    U = np.array(tab.U)
    total = X.sum(axis=0)
    parts = []
    chi = np.zeros(k, dtype=np.int64)
    for j in range(nhi):
        if j:
            u = (j & -j).bit_length() - 1        # Gray order on high bits
            g = j ^ (j >> 1)
            sgn = 1 if (g >> u) & 1 else -1
            chi += sgn * hi[u]
        c1 = cLo + chi
        c0 = total - c1
        lw = (T[c1].sum(axis=1) + U[c1.sum(axis=1)]
              + T[c0].sum(axis=1) + U[c0.sum(axis=1)])
        m = lw.max()
        parts.append(m + log(np.exp(lw - m).sum()))
    return logsumexp(parts) - n * log(2.0)

# ------------------------------------------------------------- data loading

def load_matrix(json_path, channel="closed", vocab=None):
    """Build (n x k) count matrix from a counts_*.json. Channel keys are the
    union over chapters (sorted) unless an explicit vocab list is given.
    Returns (X, vocab, chapters)."""
    obj = json.load(open(json_path))
    chs = obj["chapters"]
    if vocab is None:
        keys = set()
        for c in chs:
            keys.update(c[channel])
        vocab = sorted(keys)
    idx = {v: i for i, v in enumerate(vocab)}
    X = [[0] * len(vocab) for _ in chs]
    for r, c in enumerate(chs):
        for kk, v in c[channel].items():
            if kk in idx:
                X[r][idx[kk]] += v
    return X, vocab, [c["chapter"] for c in chs]

# ---------------------------------------------------------------- selfcheck

def _rand_table(n, k, seed, hi=6):
    rng = random.Random(seed)
    return [[rng.randint(0, hi) for _ in range(k)] for _ in range(n)]

def _brute_contiguous(X, g, a):
    """Sum prod W over ALL g^n label assignments that realize a contiguous
    g-segmentation (non-decreasing labels, onto {0..g-1}). Returns
    (sum, number_of_kept_assignments)."""
    n = len(X)
    k = len(X[0])
    tot = Fraction(0)
    kept = 0
    for lab in product(range(g), repeat=n):
        if any(lab[i + 1] < lab[i] for i in range(n - 1)):
            continue
        if len(set(lab)) != g:
            continue
        kept += 1
        pooled = [[0] * k for _ in range(g)]
        for i, l in enumerate(lab):
            for c in range(k):
                pooled[l][c] += X[i][c]
        term = Fraction(1)
        for s in range(g):
            term *= W_exact(pooled[s], a)
        tot += term
    return tot, kept

def selfcheck(out_path):
    L = []
    ok_all = True

    def rec(name, ok, detail):
        nonlocal ok_all
        ok_all &= ok
        L.append(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}")

    # S1 collapse vs brute g^n-restricted enumeration, exact Fractions
    for (n, g, k), seed in zip([(6, 2, 3), (7, 3, 4), (8, 2, 5)], (11, 12, 13)):
        X = _rand_table(n, k, seed)
        for a in (Fraction(1, k), Fraction(1, 2), Fraction(1)):
            brute, kept = _brute_contiguous(X, g, a)
            en = enum_exact(X, g, a)
            dp = dp_exact(X, g, a)
            ok = (brute == en == dp) and kept == comb(n - 1, g - 1)
            rec(f"S1 n={n} g={g} k={k} alpha={a}", ok,
                f"kept={kept} C={comb(n-1,g-1)} brute==enum==dp: "
                f"{brute == en == dp}")
    # S2 free-assignment n=10: direct vs Gray exact; float paths in log
    n, k = 10, 4
    X = _rand_table(n, k, 21)
    a = Fraction(1, 2)
    zd = free2_direct_exact(X, a)
    zg = free2_gray_exact(X, a)
    rec("S2a free2 direct==gray exact (n=10,k=4)", zd == zg, f"equal={zd == zg}")
    lz = log(zd.numerator) - log(zd.denominator)
    lf = free2_float(X, a)
    rec("S2b free2 float vs exact log", abs(lf - lz) < 1e-9,
        f"delta={abs(lf - lz):.3e}")
    # S3 g=1 closed form
    X = _rand_table(9, 5, 31)
    a = Fraction(1, 5)
    ok = dp_exact(X, 1, a) == W_exact(prefix_sums(X)[len(X)], a)
    rec("S3 g=1 closed form (n=9,k=5)", ok, f"equal={ok}")
    # S4 exact vs float DP log evidence
    X = _rand_table(12, 6, 41)
    for g in (2, 3, 4):
        a = Fraction(1, 2)
        ze = Z_exact(X, g, a)
        lze = log(ze.numerator) - log(ze.denominator)
        lzf = Z_float(X, g, a)
        rec(f"S4 g={g} exact-vs-float log Z (n=12,k=6)",
            abs(lze - lzf) < 1e-9, f"delta={abs(lze - lzf):.3e}")
    # S5 boundary posterior sums to g-1 (exact), float matches exact
    X = _rand_table(8, 3, 51)
    g, a = 3, Fraction(1, 3)
    be = boundary_posterior_exact(X, g, a)
    ok = sum(be.values()) == g - 1
    rec("S5a sum_t P(boundary t)==g-1 exact (n=8,g=3)", ok, f"equal={ok}")
    bf = boundary_posterior_float(X, g, a)
    md = max(abs(float(be[t]) - bf[t]) for t in be)
    rec("S5b boundary exact vs float", md < 1e-9, f"maxdelta={md:.3e}")

    L.append(f"OVERALL: {'ALL PASS' if ok_all else 'FAILURES PRESENT'}")
    with open(out_path, "w") as f:
        f.write("\n".join(L) + "\n")
    print("\n".join(L))
    return ok_all

if __name__ == "__main__":
    if "--selfcheck" in sys.argv:
        out = sys.argv[sys.argv.index("--selfcheck") + 1] \
            if len(sys.argv) > sys.argv.index("--selfcheck") + 1 \
            else "seg_v1_selfcheck.txt"
        sys.exit(0 if selfcheck(out) else 1)
    print(__doc__)

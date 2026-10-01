"""
B3: PF/contiguity recurrence order — the GO/NO-GO gate.

For a 1-parameter slice  Z(n) = ∫_{[0,1]^5} p_A(x)^n · R(x) dx  with
R(x) = ∏_{k≠A} p_k(x)^{u_k} (small fixed exponents), find the minimal
linear recurrence  Σ_{j=0}^{r} c_j(n) Z(n+j) = 0  with c_j ∈ ℚ[n].

Method (mod-p nullspace):
  1. Represent 256·p_k as a 2^5 int kernel (multilinear).
  2. Iteratively convolve to build Q_n = (256·p_A)^n · (256^|u'|·R) as a
     dense (N+1)^5 int64 array mod a large prime p.
  3. ∫ Q_n dx  =  Σ_i Q_n[i] · ∏_e 1/(i_e+1)   (mod p).
  4. Z(n) ≡ 256^{-N} ∫Q_n  (mod p),  N = n + |u'|.
  5. Nullspace-search over (order r, coeff-degree d) for the recurrence.
  6. Repeat at a second prime; report stable (r,d).

Decision rule (the recurrence-order screen of this bench):
  r ≤ 4  → GO     |  5 ≤ r ≤ 8 → CONDITIONAL  |  r > 8 → NO-GO
"""
import os, sys, time, json, itertools
import numpy as np
from fractions import Fraction
import sympy as sp

# module resolution: the kit directory itself ships every module imported
# below; no environment override is consulted.
_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)
from jc4_patterns import PATTERNS, EDGES

# ---- multilinear kernels: 256*p_k as 2x2x2x2x2 int arrays ----------------
x1, x2, x3, x4, x5 = EDGES
KERNEL = {}
for lab, _, _, mult, p in PATTERNS:
    poly = sp.Poly(256 * p, x1, x2, x3, x4, x5, domain='ZZ')
    K = np.zeros((2, 2, 2, 2, 2), dtype=np.int64)
    for mono, c in poly.terms():
        K[mono] = int(c)
    KERNEL[lab] = K

PRIMES = [2_147_483_629, 2_147_483_587]  # < 2^31


def modinv(a, p):
    return pow(int(a) % p, p - 2, p)


def conv5d_mod(Q, K, p):
    """5D polynomial multiplication of dense Q (shape s) by 2^5 kernel K, mod p.
    Returns array of shape s+1 in each dim."""
    s = tuple(d + 1 for d in Q.shape)
    out = np.zeros(s, dtype=np.int64)
    for eps in itertools.product((0, 1), repeat=5):
        c = int(K[eps]) % p
        if c == 0:
            continue
        sl_out = tuple(slice(e, e + Q.shape[i]) for i, e in enumerate(eps))
        out[sl_out] = (out[sl_out] + c * Q) % p
    return out


def integrate_mod(Q, p):
    """∫_{[0,1]^5} Q(x) dx  mod p, where Q[i1..i5] is coeff of ∏ x_e^{i_e}."""
    s = Q.shape
    # weight[i1..i5] = ∏ 1/(i_e+1) mod p ; build as outer product
    inv = [np.array([modinv(j + 1, p) for j in range(s[d])], dtype=np.int64)
           for d in range(5)]
    W = inv[0].reshape(-1,1,1,1,1)
    acc = (Q * W) % p
    for d in range(1, 5):
        shape = [1]*5; shape[d] = s[d]
        acc = (acc * inv[d].reshape(shape)) % p
    return int(acc.sum() % p)


def Z_sequence_mod(varying, fixed, nmax, p, alpha=1, log=None):
    """Compute Z(n) mod p for n=0..nmax where
       Z(n) = ∫ p_{varying}^n · ∏ p_k^{fixed[k]} · ∏ x_e^{alpha-1} dx.
    Only alpha=1 supported here (uniform prior); general alpha → shift weight.
    Returns list of ints mod p, plus timing."""
    assert alpha == 1
    Kv = KERNEL[varying]
    # Build R = ∏ (256 p_k)^{fixed[k]}
    R = np.ones((1, 1, 1, 1, 1), dtype=np.int64)
    nfix = 0
    for lab, e in fixed.items():
        for _ in range(e):
            R = conv5d_mod(R, KERNEL[lab], p)
            nfix += 1
    inv256 = modinv(256, p)
    Zs = []
    Q = R.copy()
    t0 = time.time()
    for n in range(nmax + 1):
        N = n + nfix
        I = integrate_mod(Q, p)
        Zn = (I * pow(inv256, N, p)) % p
        Zs.append(Zn)
        if log and (n % 5 == 0 or n == nmax):
            dt = time.time() - t0
            mem = Q.nbytes / 1e6
            print(f"[{log}] n={n:3d} N={N:3d} shape={Q.shape} mem={mem:.0f}MB "
                  f"t={dt:.1f}s", flush=True)
        if n < nmax:
            Q = conv5d_mod(Q, Kv, p)
    return Zs


def find_recurrence_mod(Zs, p, rmax=12, dmax=12, overdet=3):
    """Find minimal (r,d) with a poly-coeff recurrence
       Σ_{j=0}^r Σ_{i=0}^d a_{ij} n^i Z(n+j) = 0  for all valid n,
    via mod-p nullspace. Search by increasing ncols=(r+1)(d+1).
    Returns (r, d, nullvec) or None."""
    M = len(Zs)
    cands = sorted(((r, d) for r in range(1, rmax + 1) for d in range(dmax + 1)),
                   key=lambda rd: ((rd[0] + 1) * (rd[1] + 1), rd[0]))
    for r, d in cands:
        ncols = (r + 1) * (d + 1)
        nrows = M - r
        if nrows < ncols + overdet:
            continue
        A = np.empty((nrows, ncols), dtype=object)
        for n in range(nrows):
            for j in range(r + 1):
                zj = Zs[n + j]
                npow = 1
                for i in range(d + 1):
                    A[n, j * (d + 1) + i] = (npow * zj) % p
                    npow = (npow * n) % p
        null = nullspace_mod(A, p)
        if null is not None and verify_recurrence(Zs, r, d, null, p):
            return r, d, null
    return None


def nullspace_mod(A, p):
    """Return a single nullspace vector of A over GF(p), or None if full rank."""
    A = np.array(A, dtype=object)
    m, n = A.shape
    A = A % p
    piv_cols = []
    row = 0
    for col in range(n):
        # find pivot
        sel = None
        for r in range(row, m):
            if A[r, col] % p != 0:
                sel = r; break
        if sel is None:
            continue
        A[[row, sel]] = A[[sel, row]]
        inv = modinv(A[row, col], p)
        A[row] = (A[row] * inv) % p
        for r in range(m):
            if r != row and A[r, col] % p != 0:
                A[r] = (A[r] - A[r, col] * A[row]) % p
        piv_cols.append(col)
        row += 1
        if row == m:
            break
    free = [c for c in range(n) if c not in piv_cols]
    if not free:
        return None
    # build a nullvector with the LAST free var = 1
    f = free[-1]
    v = [0] * n
    v[f] = 1
    for i, pc in enumerate(piv_cols):
        v[pc] = (-int(A[i, f])) % p
    return v


def verify_recurrence(Zs, r, d, v, p):
    """Check Σ_{j,i} v[j*(d+1)+i] n^i Z[n+j] == 0 mod p for all n."""
    M = len(Zs)
    for n in range(M - r):
        s = 0
        for j in range(r + 1):
            for i in range(d + 1):
                s = (s + v[j * (d + 1) + i] * pow(n, i, p) * Zs[n + j]) % p
        if s % p != 0:
            return False
    return True


# ---- exact rational Z(n) via sympy (small-n cross-check) -----------------
def Z_exact_small(varying, fixed, n):
    """Brute sympy.integrate for tiny n (cross-check)."""
    expr = sp.Integer(1)
    for lab, e in fixed.items():
        idx = [i for i, P in enumerate(PATTERNS) if P[0] == lab][0]
        expr *= PATTERNS[idx][4] ** e
    idx = [i for i, P in enumerate(PATTERNS) if P[0] == varying][0]
    expr *= PATTERNS[idx][4] ** n
    for xe in EDGES:
        expr = sp.integrate(sp.expand(expr), (xe, 0, 1))
    return sp.nsimplify(expr)


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--varying", default="xxxx")
    ap.add_argument("--fixed", default="", help="e.g. xxyy:2,xyzw:1")
    ap.add_argument("--nmax", type=int, default=25)
    ap.add_argument("--rmax", type=int, default=12)
    ap.add_argument("--dmax", type=int, default=12)
    ap.add_argument("--check", type=int, default=3, help="cross-check Z(0..k) vs sympy")
    args = ap.parse_args()

    fixed = {}
    if args.fixed:
        for tok in args.fixed.split(","):
            lab, e = tok.split(":")
            fixed[lab] = int(e)

    print(f"=== B3 recurrence probe: varying={args.varying}, fixed={fixed}, "
          f"nmax={args.nmax} ===", flush=True)

    # cross-check small n
    if args.check:
        p = PRIMES[0]
        Zs = Z_sequence_mod(args.varying, fixed, args.check, p)
        for n in range(args.check + 1):
            zex = Z_exact_small(args.varying, fixed, n)
            zmod = Fraction(zex).numerator * modinv(Fraction(zex).denominator, p) % p
            ok = "OK" if zmod == Zs[n] else "FAIL"
            print(f"  cross-check n={n}: exact={zex}  mod-p match: {ok}")
            assert zmod == Zs[n], f"mismatch at n={n}"

    results = {}
    for ip, p in enumerate(PRIMES):
        t0 = time.time()
        Zs = Z_sequence_mod(args.varying, fixed, args.nmax, p,
                            log=f"p{ip}" if ip == 0 else None)
        t1 = time.time()
        rec = find_recurrence_mod(Zs, p, rmax=args.rmax, dmax=args.dmax)
        t2 = time.time()
        if rec is None:
            print(f"[p={p}] NO recurrence found with r<={args.rmax}, d<={args.dmax}")
            results[p] = None
        else:
            r, d, v = rec
            print(f"[p={p}] recurrence FOUND: order r={r}, coeff-deg d={d} "
                  f"(Z-seq {t1-t0:.1f}s, null {t2-t1:.1f}s)")
            results[p] = (r, d)
    # 2-prime stability
    rs = set(v for v in results.values() if v)
    if len(rs) == 1:
        r, d = rs.pop()
        verdict = ("GO" if r <= 4 else "CONDITIONAL" if r <= 8 else "NO-GO")
        print(f"\n>>> STABLE across 2 primes: order={r}, deg={d} → {verdict}")
    else:
        print(f"\n>>> UNSTABLE across primes: {results}")

#!/usr/bin/env python3
"""zseries.py — BAYES-MIXTURES Phase 0b: ray series z(n) = Z(n*U0) engines.

Engines:
  * z_ray_1var(k, U0, nmax)          : exact Fractions, one k-state variable
                                       (A=I collapse: Z = prod U_v! (k-1)!^2 *
                                       sum_m W_U(m) m!(N-m)!/[(N+1)!(m+k-1)!(N-m+k-1)!],
                                       W_U(m) = #{0<=x<=U, |x|=m}), O(N^2)/term.
  * z_ray_1var_modp(k, U0, nmax, p)  : same mod p (int64 numpy poly-DP).
  * z_ray_table_modp(U0, nmax, p)    : general 2-block table, dense phi-DP mod p
                                       + per-m assembly mod p.  p < 2^25 required
                                       (int64 accumulate-without-overflow margin).
Cross-gates: mod-p engines vs exact engines / lsx_direct.Z_phi at small n.
"""

import os
import sys
from fractions import Fraction
from math import comb, factorial

import numpy as np

# sibling imports resolve to the vendored copies beside this file
# (MIXALOT_PILOT_DIR overrides the default sibling path)
sys.path.insert(0, os.environ.get("MIXALOT_PILOT_DIR",
                                  os.path.dirname(os.path.abspath(__file__))))
from lsx_direct import Model, Z_phi  # noqa: E402


# ----------------------------------------------------------------------------
# Tier 1: one k-state variable (A = I collapse)
# ----------------------------------------------------------------------------

def z_1var_exact(U):
    k = len(U)
    N = sum(U)
    # W_U(m) via bounded-parts poly DP (exact ints)
    W = [1]
    for u in U:
        newW = [0] * (len(W) + u)
        for m, w in enumerate(W):
            if w:
                for x in range(u + 1):
                    newW[m + x] += w
        W = newW
    pref = 1
    for u in U:
        pref *= factorial(u)
    pref *= factorial(k - 1) ** 2
    Z = Fraction(0)
    for m in range(N + 1):
        if W[m]:
            Z += Fraction(W[m] * factorial(m) * factorial(N - m),
                          factorial(m + k - 1) * factorial(N - m + k - 1))
    return Z * Fraction(pref, factorial(N + 1))


def z_ray_1var(k, U0, nmax, progress=False):
    assert len(U0) == k
    return [z_1var_exact([n * u for u in U0]) for n in range(nmax + 1)]


def z_ray_1var_modp(k, U0, nmax, p):
    """Same series mod p.  Requires p > N_max + k."""
    Nmax = nmax * sum(U0)
    assert p > Nmax + k + 1, "prime too small for factorial arguments"
    # factorial and inverse factorial tables mod p
    F = np.ones(Nmax + k + 1, dtype=np.int64)
    for i in range(1, len(F)):
        F[i] = (F[i - 1] * i) % p
    Finv = np.ones_like(F)
    Finv[-1] = pow(int(F[-1]), p - 2, p)
    for i in range(len(F) - 1, 0, -1):
        Finv[i - 1] = (Finv[i] * i) % p
    out = []
    for n in range(nmax + 1):
        U = [n * u for u in U0]
        N = sum(U)
        W = np.zeros(N + 1, dtype=np.int64)
        W[0] = 1
        top = 0
        for u in U:
            newW = np.zeros(N + 1, dtype=np.int64)
            for x in range(u + 1):
                newW[x:top + x + 1] = (newW[x:top + x + 1] + W[:top + 1]) % p
            W = newW
            top += u
        pref = 1
        for u in U:
            pref = (pref * int(F[u])) % p
        pref = (pref * int(F[k - 1]) ** 2) % p
        pref = (pref * int(Finv[N + 1])) % p
        m = np.arange(N + 1)
        terms = (W * F[m] % p) * F[N - m] % p
        terms = (terms * Finv[m + k - 1]) % p
        terms = (terms * Finv[N - m + k - 1]) % p
        out.append(int(np.sum(terms % p) % p * pref % p))
    return out


# ----------------------------------------------------------------------------
# Tier 2: general 2-block table mod p (dense DP, p < 2^25)
# ----------------------------------------------------------------------------

def Z_table_modp(U, p):
    """Z(U) mod p for a 2-way table (uniform prior).  p < 2^25."""
    assert p < (1 << 25), "need p < 2^25 for int64 accumulation margins"
    model = Model.from_table(U)
    N = model.N
    t1 = model.block_len[0] - 1
    t2 = model.block_len[1] - 1
    assert p > N + max(t1, t2) + 2, "prime too small"
    caps1 = list(model.AU[0][1:])
    caps2 = list(model.AU[1][1:])
    dims = [N + 1] + [c + 1 for c in caps1] + [c + 1 for c in caps2]
    phi = np.zeros(dims, dtype=np.int64)
    phi[(0,) * len(dims)] = 1
    for a_v, Uv in model.cells:
        d = [1] + list(a_v[0][1:]) + list(a_v[1][1:])
        new = np.zeros_like(phi)
        for x in range(Uv + 1):
            C = comb(Uv, x) % p
            src = tuple(slice(0, dims[kk] - x * d[kk]) if d[kk] else slice(None)
                        for kk in range(len(dims)))
            dst = tuple(slice(x * d[kk], dims[kk]) if d[kk] else slice(None)
                        for kk in range(len(dims)))
            new[dst] = (new[dst] + C * phi[src]) % p
        phi = new
    # factorial tables
    top = N + max(t1, t2) + 2
    F = np.ones(top, dtype=np.int64)
    for i in range(1, top):
        F[i] = (F[i - 1] * i) % p
    Finv = np.ones_like(F)
    Finv[-1] = pow(int(F[-1]), p - 2, p)
    for i in range(top - 1, 0, -1):
        Finv[i - 1] = (Finv[i] * i) % p
    # per-m assembly: block weight grids W1[m,kept1], W2[m,kept2]
    n1 = int(np.prod([c + 1 for c in caps1])) if caps1 else 1
    n2 = int(np.prod([c + 1 for c in caps2])) if caps2 else 1
    grids1 = np.indices([c + 1 for c in caps1]).reshape(t1, n1)
    grids2 = np.indices([c + 1 for c in caps2]).reshape(t2, n2)
    sum1 = grids1.sum(axis=0)
    sum2 = grids2.sum(axis=0)
    Z = 0
    for m in range(N + 1):
        phim = phi[m].reshape(n1, n2)
        if not phim.any():
            continue
        # W1: F[b0]*prod_j F[b_j]*F[cap_j-b_j]*F[c0], b0=m-sum1, c0=(N-m)-(sumcaps1-sum1)
        b0 = m - sum1
        c0 = (N - m) - (sum(caps1) - sum1)
        ok1 = (b0 >= 0) & (c0 >= 0)
        W1 = np.where(ok1, F[np.clip(b0, 0, None)] * F[np.clip(c0, 0, None)] % p, 0)
        for j in range(t1):
            W1 = W1 * F[grids1[j]] % p * F[caps1[j] - grids1[j]] % p
        b0 = m - sum2
        c0 = (N - m) - (sum(caps2) - sum2)
        ok2 = (b0 >= 0) & (c0 >= 0)
        W2 = np.where(ok2, F[np.clip(b0, 0, None)] * F[np.clip(c0, 0, None)] % p, 0)
        for j in range(t2):
            W2 = W2 * F[grids2[j]] % p * F[caps2[j] - grids2[j]] % p
        # S_m = W1 . phim . W2 mod p   (p<2^25: products <2^50, safe row sums)
        v = phim @ (W2 % p)            # sums of n2 terms, each < p^2 < 2^50
        v %= p
        S = int((W1 % p) @ v % p)
        pref = int(F[m]) * int(F[N - m]) % p * int(Finv[N + 1]) % p
        pref = pref * int(F[t1]) % p * int(F[t1]) % p
        pref = pref * int(F[t2]) % p * int(F[t2]) % p
        pref = pref * int(Finv[m + t1]) % p * int(Finv[N - m + t1]) % p
        pref = pref * int(Finv[m + t2]) % p * int(Finv[N - m + t2]) % p
        Z = (Z + S * pref) % p
    return Z


def z_ray_table_modp(U0, nmax, p, progress=None):
    out = []
    for n in range(nmax + 1):
        U = [[n * x for x in row] for row in U0]
        if n == 0:
            out.append(1 % p)   # empty product: Z(0)=1
            continue
        out.append(Z_table_modp(U, p))
        if progress and n % progress == 0:
            print(f"  n={n}", file=sys.stderr, flush=True)
    return out


# ----------------------------------------------------------------------------
# selftest
# ----------------------------------------------------------------------------

def selftest():
    ok = True

    def chk(name, cond):
        nonlocal ok
        ok &= cond
        print(f"  [{'PASS' if cond else 'FAIL'}] {name}")

    # Tier-1 exact vs lsx_direct engines
    from lsx_direct import model_1var, Z_xsum
    for U in [(1, 1), (2, 1), (3, 2)]:
        chk(f"1var exact collapse U={U} == x-sum",
            z_1var_exact(list(U)) == Z_xsum(model_1var(U)))
    m3 = Model([1], [(((1, 0, 0),), 2), (((0, 1, 0),), 1), (((0, 0, 1),), 1)])
    chk("1var k=3 collapse == x-sum", z_1var_exact([2, 1, 1]) == Z_xsum(m3))
    # Tier-1 modp vs exact
    p = (1 << 24) - 3  # 16777213, prime
    zs = z_ray_1var(2, (1, 1), 8)
    zp = z_ray_1var_modp(2, (1, 1), 8, p)
    match = all((z.numerator * pow(z.denominator, p - 2, p)) % p == zp[n]
                for n, z in enumerate(zs))
    chk("1var k=2 ray mod-p == exact (n<=8)", match)
    # Tier-2 table modp vs exact Z_phi
    for U in [[[1, 1], [1, 0]], [[2, 1], [1, 1]], [[2, 2], [2, 2]]]:
        ze = Z_phi(Model.from_table(U))
        zm = Z_table_modp(U, p)
        chk(f"table {U} mod-p == exact",
            (ze.numerator * pow(ze.denominator, p - 2, p)) % p == zm)
    # known value gate: Tier-0 recurrence residues (independently verified)
    zs = z_ray_1var(2, (1, 1), 10)
    rec_ok = all(
        (n + 1) * (n + 2) ** 2 * zs[n]
        - 2 * (n + 2) * (4 * n * n + 18 * n + 21) * zs[n + 1]
        + 4 * (n + 3) * (2 * n + 5) ** 2 * zs[n + 2] == 0
        for n in range(9))
    chk("Tier-0 recurrence annihilates exact k=2 ray", rec_ok)
    print("SELFTEST:", "PASS" if ok else "FAIL")
    return ok


if __name__ == "__main__":
    sys.exit(0 if selftest() else 1)


# ----------------------------------------------------------------------------
# Veronese coin model (k=1, s_1=4, t_1=1) mod p: states (m, ell)
# ----------------------------------------------------------------------------

def Z_coin_modp(U, p):
    """Z(U) mod p for the reduced coin model (columns (4-v, v), UNSCALED).
    b = (4m-ell, ell); 2D DP over (m, ell). p < 2^25."""
    assert p < (1 << 25) and len(U) == 5
    N = sum(U)
    L = sum(v * u for v, u in enumerate(U))
    assert p > 4 * N + 2
    import numpy as np
    phi = np.zeros((N + 1, L + 1), dtype=np.int64)
    phi[0, 0] = 1
    from math import comb
    for v, u in enumerate(U):
        new = np.zeros_like(phi)
        for x in range(u + 1):
            C = comb(u, x) % p
            if v == 0:
                new[x:, :] = (new[x:, :] + C * phi[:N + 1 - x, :]) % p
            else:
                new[x:, x * v:] = (new[x:, x * v:]
                                   + C * phi[:N + 1 - x, :L + 1 - x * v]) % p
        phi = new
    top = 4 * N + 2
    F = np.ones(top, dtype=np.int64)
    for i in range(1, top):
        F[i] = (F[i - 1] * i) % p
    Finv = np.ones_like(F)
    Finv[-1] = pow(int(F[-1]), p - 2, p)
    for i in range(top - 1, 0, -1):
        Finv[i - 1] = (Finv[i] * i) % p
    Z = 0
    ells = np.arange(L + 1)
    for m in range(N + 1):
        row = phi[m]
        if not row.any():
            continue
        b1 = ells
        b0 = 4 * m - ells
        c1 = L - ells
        c0 = 4 * (N - m) - c1
        ok = (b0 >= 0) & (c0 >= 0)
        w = np.where(ok, F[np.clip(b0, 0, None)] * F[b1] % p, 0)
        w = np.where(ok, w * F[np.clip(c0, 0, None)] % p * F[c1] % p, 0)
        S = int(np.sum(row * w % p) % p)
        pref = int(F[m]) * int(F[N - m]) % p * int(Finv[N + 1]) % p
        pref = pref * int(Finv[4 * m + 1]) % p * int(Finv[4 * (N - m) + 1]) % p
        Z = (Z + S * pref) % p
    return Z

#!/usr/bin/env python3
r"""
annihilator.py — ANNIHILATOR: find the annihilating (Picard–Fuchs / holonomic)
operator of a series.  The ansatz sweep frontend and the GF(p) numpy kernel
live here; the series_pf and nullspace_fast paths are import shims.  Related
names, for grep: series_pf, nullspace_fast, pf_from_series,
find_recurrence_fast, nullspace_mod_fast.

DATA-route sibling of tools/pf_rank.jl: pf_rank.jl needs the integrand polynomial
(Griffiths–Dwork on the Baikov/Feynman rep); this needs ONLY the exact moment
sequence a[0..N].  Use this when you have a fast series generator (LGF moments,
maxcut Taylor, diagonal coefficients) but no clean integrand.

THIN FRONTEND to tools/nullspace_fast (numpy int64 mod-p RREF, p < 2^31), adding:
  * ODE-order-first search (s outer, r inner) — the piece this frontend adds (cf. Chen–Jaroschek–Kauers–Singer on order-degree curves).  The
    default nullspace_fast sort is by #unknowns, which returns inflated ODE
    orders carrying apparent-singularity indicial factors.
    s-first guarantees the FIRST hit is the true minimal θ-order.
  * mode='auto': mod-p s-first scan pins (r,s), then ONE exact ℚ-nullspace at
    that single (r,s) reconstructs coefficients — avoids the Fraction-elimination
    blowup of a full ℚ scan (the modp-first-then-ℚ pattern).
  * rec_to_theta: recurrence → θ-form ODE.

Note: numpy engine requires p < 2^31 (int64 overflow); default p = 2^31−1.
"""
from __future__ import annotations
import sys, os
from fractions import Fraction as F
from math import gcd

P31 = (1 << 31) - 1  # Mersenne prime, int64-safe for numpy engine


# ===== GF(p) kernel (nullspace_fast; tools/annihilator/nullspace_fast/ forwards here) =====
import numpy as np


def nullspace_mod_fast(A, p):
    """A: (m,n) int64 array, entries in [0,p). Return one nullvector or None."""
    A = np.array(A, dtype=np.int64) % p
    m, n = A.shape
    piv_cols = []
    row = 0
    for col in range(n):
        # find pivot
        sel = -1
        for r in range(row, m):
            if A[r, col] != 0:
                sel = r; break
        if sel < 0:
            continue
        if sel != row:
            A[[row, sel]] = A[[sel, row]]
        inv = pow(int(A[row, col]), p - 2, p)
        A[row] = (A[row] * inv) % p
        # eliminate (vectorized)
        mask = np.ones(m, dtype=bool); mask[row] = False
        factors = A[mask, col].copy()
        nz = factors != 0
        idx = np.where(mask)[0][nz]
        if len(idx):
            A[idx] = (A[idx] - factors[nz, None] * A[row]) % p
        piv_cols.append(col)
        row += 1
        if row == m:
            break
    free = [c for c in range(n) if c not in piv_cols]
    if not free:
        return None
    f = free[-1]
    v = np.zeros(n, dtype=np.int64)
    v[f] = 1
    for i, pc in enumerate(piv_cols):
        v[pc] = (-int(A[i, f])) % p
    return v.tolist()


def find_recurrence_fast(Zs, p, rmax=40, dmax=40, overdet=5, verbose=False,
                         prefer='unknowns'):
    """prefer='unknowns' (default, legacy): sort by (r+1)(d+1) — cheapest system
    first, but can inflate ODE order via apparent-singularity indicial factors.
    prefer='ode_order': d-first (d outer, r inner) — first hit is the true
    MINIMAL ODE order (d ↔ θ-order under rec→θ). Use this for PF operators."""
    M = len(Zs)
    Zs = np.array(Zs, dtype=np.int64)
    if prefer == 'ode_order':
        cands = sorted(((r, d) for r in range(1, rmax + 1) for d in range(dmax + 1)),
                       key=lambda rd: (rd[1], rd[0]))
    else:
        cands = sorted(((r, d) for r in range(1, rmax + 1) for d in range(dmax + 1)),
                       key=lambda rd: ((rd[0] + 1) * (rd[1] + 1), rd[0]))
    for r, d in cands:
        ncols = (r + 1) * (d + 1)
        nrows = M - r
        if nrows < ncols + overdet:
            continue
        # build matrix int64
        ns = np.arange(nrows, dtype=np.int64)
        npows = np.ones((nrows, d + 1), dtype=np.int64)
        for i in range(1, d + 1):
            npows[:, i] = (npows[:, i - 1] * ns) % p
        A = np.empty((nrows, ncols), dtype=np.int64)
        for j in range(r + 1):
            zj = Zs[j:j + nrows].reshape(-1, 1)
            A[:, j * (d + 1):(j + 1) * (d + 1)] = (npows * zj) % p
        v = nullspace_mod_fast(A, p)
        if v is not None:
            # verify
            ok = True
            for nn in range(nrows):
                s = 0
                for j in range(r + 1):
                    cj = sum(v[j*(d+1)+i] * pow(nn, i, p) for i in range(d + 1)) % p
                    s = (s + cj * int(Zs[nn + j])) % p
                if s != 0:
                    ok = False; break
            if ok:
                if verbose:
                    print(f"    found r={r} d={d} ncols={ncols}")
                return r, d, v
    return None


# ===== ansatz sweep (series_pf; tools/annihilator/series_pf.py forwards here) =====


def _to_modp(x, p):
    if isinstance(x, F):
        return x.numerator % p * pow(x.denominator % p, p - 2, p) % p
    return int(x) % p


def _nullspace_Q(rows, ncols):
    """Exact ℚ nullspace via Gaussian elimination."""
    M = [list(r) for r in rows]; nrows = len(M); piv = []; rlead = 0
    for c in range(ncols):
        pr = next((r for r in range(rlead, nrows) if M[r][c] != 0), None)
        if pr is None: continue
        M[rlead], M[pr] = M[pr], M[rlead]
        pv = M[rlead][c]; M[rlead] = [x / pv for x in M[rlead]]
        for r in range(nrows):
            if r != rlead and M[r][c] != 0:
                f = M[r][c]
                M[r] = [M[r][i] - f * M[rlead][i] for i in range(ncols)]
        piv.append(c); rlead += 1
        if rlead == nrows: break
    free = [c for c in range(ncols) if c not in piv]; basis = []
    for fc in free:
        v = [F(0)] * ncols; v[fc] = F(1)
        for i, pc in enumerate(piv): v[pc] = -M[i][fc]
        L = 1
        for x in v: L = L * x.denominator // gcd(L, x.denominator)
        v = [x * L for x in v]
        g = 0
        for x in v:
            if x != 0: g = gcd(g, abs(int(x)))
        if g > 1: v = [x / g for x in v]
        basis.append(v)
    return basis


def pf_from_series(a, rmax=30, smax=8, nverify=15, mode='auto', p=P31):
    """
    Minimal-order Picard-Fuchs recurrence Σ_j c_j(n) a_{n+j}=0 from a series a[0..N].
    mode='modp': fast rank-only scan over F_p, returns (r,s,None,0) — order certified, no coeffs.
    mode='exact': ℚ-nullspace at the modp-pinned (r,s); returns (r,s,coeffs,n_heldout).
    mode='auto': modp scan to pin (r,s), then ℚ-reconstruct at that single (r,s).
    a can be list[Fraction] (exact) or list[int] (interpreted mod p in modp mode).
    Returns (r, s, coeffs_dict_or_None, n_heldout_verified) or None if none in range.
    """
    assert p < (1 << 31), "numpy int64 engine requires p < 2^31"
    N = len(a)
    ap = [_to_modp(x, p) for x in a]
    hit = find_recurrence_fast(ap, p, rmax=rmax, dmax=smax, overdet=5,
                               prefer='ode_order')
    if hit is None:
        return None
    r, s, _ = hit
    if mode == 'modp':
        return r, s, None, 0
    # exact ℚ reconstruct at the single pinned (r,s)
    aQ = [x if isinstance(x, F) else F(x) for x in a]
    nunk = (r + 1) * (s + 1)
    neq = N - r - nverify
    if neq < nunk + 3:
        raise ValueError(f"need ≥{nunk + 3 + r + nverify} terms for (r,s)=({r},{s})")
    rows = []
    for n in range(neq):
        row = []
        for j in range(r + 1):
            an = aQ[n + j]; nk = F(1)
            for k in range(s + 1):
                row.append(an * nk); nk *= n
        rows.append(row)
    ker = _nullspace_Q(rows, nunk)
    for vec in ker:
        ok = True
        for n in range(neq, N - r):
            val = F(0); idx = 0
            for j in range(r + 1):
                an = aQ[n + j]
                for k in range(s + 1):
                    val += vec[idx] * an * F(n) ** k; idx += 1
            if val != 0: ok = False; break
        if ok:
            coeffs = {(j, k): vec[j * (s + 1) + k]
                      for j in range(r + 1) for k in range(s + 1)
                      if vec[j * (s + 1) + k] != 0}
            return r, s, coeffs, N - r - neq
    return r, s, None, 0  # modp hit but ℚ-verify failed (unlucky prime / insufficient N)


def rec_to_theta(r, s, coeffs):
    """Σ_j c_j(n) a_{n+j}=0  →  ODE  Σ_j x^{r-j} c_j(θ-j) f = 0,  θ=x∂_x, f=Σ a_n x^n.
    Returns {(xpow, thetapow): Fraction}.  (sympy only here, not in hot path.)"""
    from sympy import symbols, Poly, expand, Rational
    n = symbols('n'); ode = {}
    for j in range(r + 1):
        cj = sum(Rational(coeffs.get((j, k), F(0))) * n ** k for k in range(s + 1))
        pcs = Poly(expand(cj.subs(n, n - j)), n).all_coeffs()
        for k, ck in enumerate(reversed(pcs)):
            ode[(r - j, k)] = ode.get((r - j, k), F(0)) + F(ck.p, ck.q)
    return ode


def pretty_rec(r, s, coeffs):
    terms = []
    for j in range(r + 1):
        poly = [f"{c}" if k == 0 else f"{c}*n" if k == 1 else f"{c}*n^{k}"
                for k in range(s, -1, -1) if (c := coeffs.get((j, k), F(0))) != 0]
        if poly: terms.append(f"({'+'.join(poly)})*a[n+{j}]")
    return " + ".join(terms) + " = 0"


def moments_sc(d, N):
    """Return-probability moments μ_d(2n) for the d-dim simple-cubic lattice, n=0..N.
    Exact Fractions.  Inlined from lgf_oracle_v2.py (d-fold conv of 1/(n!)²)."""
    b = [F(1)]; fac = 1
    for n in range(1, N + 1):
        fac *= n; b.append(F(1, fac * fac))
    conv = [F(0)] * (N + 1); conv[0] = F(1)
    for _ in range(d):
        nxt = [F(0)] * (N + 1)
        for i in range(N + 1):
            if conv[i] == 0: continue
            for j in range(N + 1 - i): nxt[i + j] += conv[i] * b[j]
        conv = nxt
    mu = []; fac2 = 1
    for n in range(N + 1):
        if n > 0: fac2 *= (2 * n - 1) * (2 * n)
        mu.append(F(fac2) * conv[n] / (F(4) ** n * F(d) ** (2 * n)))
    return mu


if __name__ == '__main__':
    import argparse, time
    ap = argparse.ArgumentParser(
        description="ANNIHILATOR: PF/holonomic operator of a series (see GUIDE.md).")
    ap.add_argument('--selftest', action='store_true',
                    help="run the built-in battery (the default action); "
                         "exit 0 only if every rung passes")
    ap.parse_args()
    print("[series_pf self-test] d-dim simple-cubic LGF, mode='auto' (ode_order-first)",
          flush=True)
    # Known minimal (r,s) for the d-dim simple-cubic LGF moment recurrences.
    expected = {2: (1, 2), 3: (2, 3), 4: (2, 4), 5: (3, 5)}
    all_ok = True
    for d in (2, 3, 4, 5):
        t0 = time.time()
        a = moments_sc(d, 100)
        res = pf_from_series(a, rmax=14, smax=8, nverify=15, mode='auto')
        if res is None:
            print(f"  d={d}: NO recurrence in range  FAIL", flush=True)
            all_ok = False
            continue
        r, s, coeffs, nver = res
        ok = (r, s) == expected[d] and coeffs is not None
        if not ok:
            all_ok = False
        tag = "PASS" if ok else (
            f"FAIL (want (r,s)={expected[d]})" if (r, s) != expected[d]
            else "FAIL (no exact coefficients)")
        print(f"  d={d}: (r,s)=({r},{s})  held-out={nver}  [{time.time()-t0:.2f}s]  {tag}",
              flush=True)
    print("ALL PASS" if all_ok else "SELFTEST FAILED", flush=True)
    sys.exit(0 if all_ok else 1)

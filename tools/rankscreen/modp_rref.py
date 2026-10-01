#!/usr/bin/env python3
"""Mod-p RREF / rank + inconsistency detection for exact-Q closure row systems.

ALGORITHM CONTRACT: mirrors the certified exact eliminators of the pm_core /
pm_ext engine lineage (window_core.gauss_rank == elim.gauss_rank_mpq, verified
cross-implementations): rows processed in list order; incoming row reduced
against all existing pivot rows; pivot column = MIN nonzero column of the
reduced row; a row reducing to zero coefficients with nonzero rhs is
INCONSISTENT.  Under this contract the reduced form of every incoming row is
the unique representative of row + span(pivots) vanishing on all pivot
columns, so rank, pivot set and inconsistent-row set are what the exact
algorithm computes, imaged through Q -> F_p (valid whenever no pivot leading
coefficient hits 0 mod p; cross-prime agreement is the certificate — see
screen.py and the manual for the false-agreement bounds).

Two backends, IDENTICAL verdict semantics (gate-verified identical output):
  rref_modp_sparse : dict-of-int rows, pure Python ints.  Primes < 2^62.
  rref_modp_dense  : numpy int64, fully-reduced (RREF) pivot block, one
                     matvec per incoming row.  Primes < 2^25 REQUIRED
                     (overflow bound: ncols_max * p^2 < 2^63 needs
                     p < sqrt(2^63/ncols); enforced for ncols <= 8192).
Backend choice by measured speed on the real 3064-col hard-twin system:
dense wins (see manual); sparse kept as the wide-prime and
low-memory fallback and as the semantics reference.

A Fraction coefficient num/den maps to num * den^-1 mod p; if p | den the
prime is UNUSABLE for the system: BadPrime is raised and the caller draws a
replacement prime (screen.py does this automatically and records it).

Returns ModpResult(rank, pivots, incon_idx, prime, backend, secs):
  pivots    : tuple of (pivot_col, row_index) in pivot-creation order
  incon_idx : tuple of row indices (into the input row list) that reduced to
              0 = nonzero mod p
"""
from collections import namedtuple
import time

ModpResult = namedtuple('ModpResult',
                        'rank pivots incon_idx prime backend secs')


class BadPrime(Exception):
    """p divides a denominator (or leading-structure) — replace the prime."""


# ---------------------------------------------------------------- primes ----
# Fixed pools (deterministic receipts).  DENSE pool: 25-bit primes — largest
# primes below 2^25; ncols*p^2 < 2^63 holds for ncols <= 8192.
# SPARSE pool: 30-bit primes (largest below 2^30).
PRIMES_DENSE = (33554393, 33554383, 33554371, 33554347, 33554341,
                33554317, 33554291, 33554273, 33554267, 33554249,
                33554171, 33554167, 33554159, 33554137, 33554123,
                33554093)
PRIMES_SPARSE = (1073741789, 1073741783, 1073741741, 1073741723,
                 1073741719, 1073741717, 1073741689, 1073741671,
                 1073741663, 1073741651, 1073741621, 1073741567,
                 1073741527, 1073741477, 1073741467, 1073741441)
DENSE_NCOLS_MAX = 8192  # overflow guard for the 25-bit dense pool


def _frac_modp(num, den, p):
    """num/den mod p (Fraction or int inputs already split)."""
    if den % p == 0:
        raise BadPrime(f"p={p} divides denominator {den}")
    n = num % p
    if n == 0:
        return 0
    if den == 1:
        return n
    return (n * pow(den % p, p - 2, p)) % p


def rows_to_modp(rows, p):
    """[( {col: Fraction}, rhs )] -> [( {col:int}, int )] mod p, zeros dropped.
    Raises BadPrime if p divides any denominator."""
    out = []
    for rdict, rhs in rows:
        d = {}
        for c, q in rdict.items():
            v = _frac_modp(q.numerator, q.denominator, p)
            if v:
                d[c] = v
        out.append((d, _frac_modp(rhs.numerator, rhs.denominator, p)))
    return out


# ---------------------------------------------------------------- sparse ----
def rref_modp_sparse(rows, p, ncols=None):
    """rows: [({col: Fraction}, Fraction rhs)] (or pre-mapped ints via
    rows_int=True path in screen worker).  Mirrors gauss_rank_mpq literally."""
    t0 = time.time()
    mrows = rows_to_modp(rows, p)
    piv = {}          # col -> index into red
    pivorder = []     # (col, original row idx) in creation order
    red = []; redr = []
    incon = []
    for idx, (rdict, rv) in enumerate(mrows):
        rdict = dict(rdict)
        for col, ri in piv.items():
            v = rdict.get(col)
            if v:
                # pivot rows are stored monic at their pivot column
                row_p = red[ri]
                for c2, v2 in row_p.items():
                    nv = (rdict.get(c2, 0) - v * v2) % p
                    if nv:
                        rdict[c2] = nv
                    elif c2 in rdict:
                        del rdict[c2]
                rv = (rv - v * redr[ri]) % p
        if rdict:
            c0 = min(rdict)
            inv = pow(rdict[c0], p - 2, p)
            row_m = {c: (v * inv) % p for c, v in rdict.items()}
            rv_m = (rv * inv) % p
            piv[c0] = len(red)
            red.append(row_m); redr.append(rv_m)
            pivorder.append((c0, idx))
        elif rv:
            incon.append(idx)
    return ModpResult(len(pivorder), tuple(pivorder), tuple(incon), p,
                      'sparse', round(time.time() - t0, 3))


# ----------------------------------------------------------------- dense ----
def rref_modp_dense(rows, p, ncols):
    """numpy int64 backend; fully-reduced pivot block so one fused
    reduction per incoming row.  Verdict-identical to sparse (the reduced
    row is the unique span-representative in both variants)."""
    import numpy as np
    if ncols > DENSE_NCOLS_MAX:
        raise ValueError(f"dense backend capped at ncols<={DENSE_NCOLS_MAX}")
    if p >= 1 << 25:
        raise ValueError("dense backend requires p < 2^25 (overflow bound)")
    t0 = time.time()
    mrows = rows_to_modp(rows, p)
    W = ncols + 1                      # last column = rhs
    cap = min(len(mrows), ncols) + 1
    P = np.zeros((cap, W), dtype=np.int64)   # monic, mutually-RREF pivot rows
    pivcols = np.empty(cap, dtype=np.int64)
    r = 0
    pivorder = []
    incon = []
    buf = np.zeros(W, dtype=np.int64)
    for idx, (rdict, rv) in enumerate(mrows):
        buf[:] = 0
        if rdict:
            cols = np.fromiter(rdict.keys(), dtype=np.int64, count=len(rdict))
            vals = np.fromiter(rdict.values(), dtype=np.int64,
                               count=len(rdict))
            buf[cols] = vals
        buf[ncols] = rv
        if r:
            f = buf[pivcols[:r]] % p          # pivot block is monic
            nz = np.nonzero(f)[0]
            if nz.size:
                buf -= f[nz] @ P[nz]          # |sum| <= r*p^2 < 2^63
                buf %= p
        nzc = np.nonzero(buf[:ncols])[0]
        if nzc.size:
            c0 = int(nzc[0])
            inv = pow(int(buf[c0]), p - 2, p)
            row_m = (buf * inv) % p
            # maintain full RREF: clear column c0 from existing pivot rows
            if r:
                colv = P[:r, c0].copy()
                nz2 = np.nonzero(colv)[0]
                if nz2.size:
                    P[nz2] -= np.outer(colv[nz2], row_m)
                    P[nz2] %= p
            P[r] = row_m
            pivcols[r] = c0
            pivorder.append((c0, idx))
            r += 1
        elif buf[ncols] % p:
            incon.append(idx)
    return ModpResult(r, tuple(pivorder), tuple(incon), p,
                      'dense', round(time.time() - t0, 3))


def rref_modp(rows, p, ncols, backend='auto'):
    if backend == 'sparse':
        return rref_modp_sparse(rows, p, ncols)
    if backend == 'dense':
        return rref_modp_dense(rows, p, ncols)
    if backend == 'auto':
        # measured on the reference 3064-col hard-twin system: dense wins
        # whenever it is applicable (see manual); sparse otherwise.
        if ncols is not None and ncols <= DENSE_NCOLS_MAX and p < (1 << 25):
            return rref_modp_dense(rows, p, ncols)
        return rref_modp_sparse(rows, p, ncols)
    raise ValueError(backend)

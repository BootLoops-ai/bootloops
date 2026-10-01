#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision.
"""vopclose.py — VoP function-level closure engine for graded 1-D path-DEs.

The reusable, sector-agnostic core of a FULL-A closure run; sector-specific
drivers are outside the scope of this package (everything sector-specific goes
through VopConfig).

Problem it solves: FULL-A function-level closure of a UT/graded path-DE
    y'_{i,k}(t) = sum_{j,m} A^{(m)}_{ij}(t) y_{j,k-m}(t)
along a 1-D path t in [0,1], when value-level PSLQ fails (mixed-weight layers)
and a global canonical eps-form is blocked (e.g. half-integer-orbit radical
blocks). Layer by layer, block by block (topological order), it builds
variation-of-parameters closed forms in the exact function model
    F(t) = sum_terms  C_tag * sqrt(rad) * R(t) * G_word(t)
(named constants; radical monomials continued along the detour contour; exact
flint fmpq rational prefactors; iterated-integral words over dlog / rational /
algebraic kernels), closed under +, Rat*, sqrt-mul, exact d/dt, exact int_0^t.

Machinery: exact rational homogeneous-solution solver with eigenvalue-bounded
denominators; block strategies FULL (full fundamental basis, incl. radical
sqrt-letter columns — linear through quartic letters, e.g. a production
sqrt(3721t^4+...+7040) Gram determinant block with indicials (1/2,1/2,0)), CONJ
(conjugation by the radical), TRI (triangularization by certified rational
columns), all certified EXACTLY (col' = A0 col as a polynomial identity);
zero-remainder symbolic DE certification of every closed layer; det-basis
apparent-pole detection (zeros of det(fundamental basis) on (0,1) are poles of
C' = M^{-1}S that cancel only in Y = M C — the contour MUST detour them; this
was an honest gate FAIL in production before the fix); spectral Chebyshev-cumint
evaluation on the detour contour (arb/acb); held-out endpoint gate harness.

Rational fitting of the exact connection entries from path samples is
delegated to the sibling package tools/ratfit (known-Q Newton + mod-p screen +
blind Thiele) — see exact_entry_from_samples(). A node-table numeric mirror (a
tier-2 layered route) and exactA repair passes (eps-recurrence, d-space 2-var
recon) are outside the scope of this package.

Validation: an equal-mass nonplanar five-point top sector, all 15 TOP
coefficients (3 masters x eps^-4..eps^0) held out 59.57-62.46d vs a deep
independent oracle. Self-test: python3 tools/vopclose/vopclose.py --test  (<60 s).
"""
import os
import sys
import time
from dataclasses import dataclass, field
from fractions import Fraction
from itertools import combinations

import flint
import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

FQ = flint.fmpq
FP = flint.fmpq_poly

def fq(x):
    if isinstance(x, FQ): return x
    f = Fraction(x)
    return FQ(f.numerator, f.denominator)

def fp(coeffs):  # low-to-high
    return FP([fq(c) for c in coeffs])

ZERO_P = FP([0]); ONE_P = FP([1])

# ───────────────────────── exact rational functions ─────────────────────────
class Rat:
    __slots__ = ('n', 'd')

    def __init__(self, n, d=None):
        if d is None: d = ONE_P
        if not isinstance(n, FP): n = fp(n if isinstance(n, (list, tuple)) else [n])
        if not isinstance(d, FP): d = fp(d if isinstance(d, (list, tuple)) else [d])
        if d == ZERO_P: raise ZeroDivisionError
        if n == ZERO_P: self.n, self.d = ZERO_P, ONE_P; return
        g = n.gcd(d)
        if g.degree() > 0: n, d = n // g, d // g
        lc = d.coeffs()[-1]
        if lc != 1: n, d = n * FP([1 / lc]), d * FP([1 / lc])
        self.n, self.d = n, d

    def is_zero(self): return self.n == ZERO_P
    def __add__(a, b): return Rat(a.n * b.d + b.n * a.d, a.d * b.d)
    def __sub__(a, b): return Rat(a.n * b.d - b.n * a.d, a.d * b.d)

    def __mul__(a, b):
        if isinstance(b, FQ): return Rat(a.n * FP([b]), a.d)
        return Rat(a.n * b.n, a.d * b.d)

    def __neg__(a): return Rat(-a.n, a.d)
    def inv(a): return Rat(a.d, a.n)
    def deriv(a): return Rat(a.n.derivative() * a.d - a.n * a.d.derivative(), a.d * a.d)
    def __eq__(a, b): return a.n == b.n and a.d == b.d
    def key(a): return (tuple(str(c) for c in a.n.coeffs()), tuple(str(c) for c in a.d.coeffs()))
    def __hash__(a): return hash(a.key())
    def ev_fq(a, x): x = fq(x); return a.n(x) / a.d(x)

    def ev_mpc(a, z):
        num = mp.mpc(0)
        for c in reversed(a.n.coeffs()): num = num * z + mp.mpf(int(c.p)) / int(c.q)
        den = mp.mpc(0)
        for c in reversed(a.d.coeffs()): den = den * z + mp.mpf(int(c.p)) / int(c.q)
        return num / den

    def __repr__(a): return f"({a.n})/({a.d})"

R_ZERO = Rat(ZERO_P); R_ONE = Rat(ONE_P)

# ─────────────────── registries (letters / kernels / constants) ───────────────────
LETTERS = []; _LK = {}
RKERNS = []; _RK = {}
AKERNS = []; _AK = {}          # (Rat, rad)
CONSTS = {'1': mp.mpc(1)}
_HCACHE = {}

def reset_registries():
    """Fresh registries (call between independent closures in one process)."""
    global CONSTS
    LETTERS.clear(); _LK.clear(); RKERNS.clear(); _RK.clear()
    AKERNS.clear(); _AK.clear(); _HCACHE.clear()
    CONSTS = {'1': mp.mpc(1)}

def reg_letter(f):
    lc = f.coeffs()[-1]
    if lc != 1: f = f * FP([1 / lc])
    k = tuple(str(c) for c in f.coeffs())
    if k in _LK: return _LK[k]
    LETTERS.append(f); _LK[k] = len(LETTERS) - 1
    return len(LETTERS) - 1

def reg_rkern(r):
    k = r.key()
    if k in _RK: return _RK[k]
    RKERNS.append(r); _RK[k] = len(RKERNS) - 1
    return len(RKERNS) - 1

def reg_akern(r, rad):
    k = (r.key(), rad)
    if k in _AK: return _AK[k]
    AKERNS.append((r, rad)); _AK[k] = len(AKERNS) - 1
    return len(AKERNS) - 1

def reg_const(tag, val):
    if tag in CONSTS:
        assert abs(CONSTS[tag] - val) < mp.mpf(10) ** (-mp.mp.dps + 12), f"const clash {tag}"
    CONSTS[tag] = val
    return tag

def rad_mul(rad1, rad2):
    """sqrt(rad1)*sqrt(rad2) = Rfac * sqrt(rad).  Rfac = prod_common l(t)/l(0)."""
    s1, s2 = set(rad1), set(rad2)
    rad = tuple(sorted(s1 ^ s2))
    R = R_ONE
    for li in s1 & s2:
        f = LETTERS[li]
        R = R * Rat(f, FP([f(fq(0))]))
    return R, rad

def rad_dlog_half(rad):
    """(1/2) sum dlog l_a as Rat."""
    R = R_ZERO
    for li in rad:
        f = LETTERS[li]
        R = R + Rat(f.derivative(), f * FP([2]))
    return R

def kern_eval_parts(kk):
    """kernel key -> (Rat, rad)."""
    if kk[0] == 'L': f = LETTERS[kk[1]]; return Rat(f.derivative(), f), ()
    if kk[0] == 'R': return RKERNS[kk[1]], ()
    return AKERNS[kk[1]]

# ─────────────────── Hermite split (exact rational integration) ───────────────────
def hermite(r):
    """r = H' + sum c_a dlog(l_a) + sum c_k K_k (registered simple-pole kernels)."""
    ck = r.key()
    if ck in _HCACHE: return _HCACHE[ck]
    H = R_ZERO; logs = {}; rems = []
    q, rem = divmod(r.n, r.d)
    if q != ZERO_P:
        qc = q.coeffs()
        H = H + Rat(FP([FQ(0)] + [qc[i] / (i + 1) for i in range(len(qc))]))
    if rem != ZERO_P:
        c, facs = r.d.factor()
        num = rem * FP([1 / c])
        pieces = []
        rest = list(facs)
        while len(rest) > 1:
            f, e = rest.pop()
            fe = f ** e
            D2 = ONE_P
            for (g, eg) in rest: D2 = D2 * g ** eg
            g0, a, b = fe.xgcd(D2)
            assert g0.degree() == 0
            a, b = a * FP([1 / g0.coeffs()[0]]), b * FP([1 / g0.coeffs()[0]])
            pieces.append(((num * b) % fe, f, e))
            num = (num * a) % D2
        f, e = rest[0]
        pieces.append((num % (f ** e), f, e))
        for (nf, f, e) in pieces:
            fp_ = f.derivative()
            g0, u, v = fp_.xgcd(f)
            u = u * FP([1 / g0.coeffs()[0]])
            while e > 1:
                s = (nf * u) % f
                w = (nf - s * fp_) // f
                H = H + Rat(-s * FP([FQ(1, e - 1)]), f ** (e - 1))
                nf = w + s.derivative() * FP([FQ(1, e - 1)])
                e -= 1
            if nf == ZERO_P: continue
            dn, df = nf.degree(), f.degree()
            cl = FQ(0)
            if dn == df - 1: cl = nf.coeffs()[-1] / (fq(df) * f.coeffs()[-1])
            if cl != 0:
                li = reg_letter(f)
                logs[li] = logs.get(li, FQ(0)) + cl
                nf = nf - fp_ * FP([cl])
            if nf != ZERO_P:
                rr = Rat(nf, f)
                lead = rr.n.coeffs()[-1]
                rk = reg_rkern(Rat(rr.n * FP([1 / lead]), rr.d))
                rems.append((lead, rk))
    out = (H, logs, rems)
    _HCACHE[ck] = out
    return out

# ─────────────── containers: {(ctag, rad, word): Rat} — closed function algebra ───────────────
def cf_zero(): return {}

def cf_const(ctag, r=None):
    return {(ctag, (), ()): (r if r is not None else R_ONE)}

def cf_add(A, B):
    out = dict(A)
    for k, r in B.items():
        if k in out:
            s = out[k] + r
            if s.is_zero(): del out[k]
            else: out[k] = s
        elif not r.is_zero(): out[k] = r
    return out

def cf_neg(A): return {k: -r for k, r in A.items()}

def cf_scale(A, c):
    c = fq(c)
    if c == 0: return {}
    return {k: r * c for k, r in A.items()}

def cf_mulrat(A, R):
    if R.is_zero(): return {}
    return {k: r * R for k, r in A.items()}

def cf_mulalg(A, R, rad):
    """A * R*sqrt(rad)."""
    if R.is_zero(): return {}
    out = cf_zero()
    for (ctag, rd, word), r in A.items():
        Rf, nrd = rad_mul(rd, rad)
        out = cf_add(out, {(ctag, nrd, word): r * R * Rf})
    return out

def cf_diff(A):
    out = cf_zero()
    for (ctag, rad, word), r in A.items():
        rp = r.deriv()
        if rad: rp = rp + r * rad_dlog_half(rad)
        if not rp.is_zero(): out = cf_add(out, {(ctag, rad, word): rp})
        if word:
            K, radk = kern_eval_parts(word[0])
            Rf, nrd = rad_mul(rad, radk)
            out = cf_add(out, {(ctag, nrd, word[1:]): r * K * Rf})
    return out

def cf_int(A):
    out = cf_zero()
    for (ctag, rad, word), r in A.items(): out = cf_add(out, _int_term(ctag, rad, word, r))
    return out

def _int_term(ctag, rad, word, r):
    if rad:
        lead = r.n.coeffs()[-1]
        ai = reg_akern(Rat(r.n * FP([1 / lead]), r.d), rad)
        return {(ctag, (), (('A', ai),) + word): Rat(FP([lead]))}
    H, logs, rems = hermite(r)
    out = cf_zero()
    if not H.is_zero():
        out = cf_add(out, {(ctag, (), word): H})
        if word == ():
            h0 = H.ev_fq(0)
            if h0 != 0: out = cf_add(out, {(ctag, (), ()): Rat(FP([-h0]))})
        else:
            K, radk = kern_eval_parts(word[0])
            out = cf_add(out, cf_neg(_int_term(ctag, radk, word[1:], H * K)))
    for li, cl in logs.items(): out = cf_add(out, {(ctag, (), (('L', li),) + word): Rat(FP([cl]))})
    for (cl, rk) in rems: out = cf_add(out, {(ctag, (), (('R', rk),) + word): Rat(FP([cl]))})
    return out

# ─────────────── exact rational solutions of Y' = (A - 1/2 dlog g) Y ───────────────
def solve_rat_solutions(Ain, extra_rad=None, dN_extra=30, eig_fn=None):
    """Exact rational solution vectors. A: nb x nb Rat matrix. eig_fn(f) ->
    numeric residue-matrix eigenvalues at a root of irreducible factor f (used
    only for the denominator BOUND; results are certified exactly downstream)."""
    nb = len(Ain)
    A = [row[:] for row in Ain]
    if extra_rad:
        g = ONE_P
        for li in extra_rad: g = g * LETTERS[li]
        half = Rat(g.derivative(), g * FP([2]))
        for a in range(nb): A[a][a] = A[a][a] - half
    facs = {}
    for a in range(nb):
        for b in range(nb):
            if A[a][b].is_zero(): continue
            c, fl = A[a][b].d.factor()
            for f, e in fl:
                k2 = tuple(str(x) for x in f.coeffs())
                if k2 not in facs: facs[k2] = [f, e]
                else: facs[k2][1] = max(facs[k2][1], e)
    D = ONE_P
    for k2, (f, emax) in facs.items():
        ev = eig_fn(f) if (eig_fn is not None and emax == 1) else []
        if ev != [] and emax == 1:
            lam = min(float(mp.re(x)) for x in ev)
            pf = max(0, int(mp.ceil(-lam - 0.25)))
        else: pf = emax + 1
        if extra_rad: pf = pf + 1
        if pf > 0: D = D * f ** pf
    dN = D.degree() + dN_extra
    Vv = []; Wv = []; Uv = []
    for a in range(nb):
        La = ONE_P
        for b in range(nb):
            if A[a][b].is_zero(): continue
            gg = La.gcd(A[a][b].d)
            La = La * (A[a][b].d // gg)
        Vv.append(La * D); Wv.append(-(La * D.derivative()))
        U = []
        for b in range(nb):
            if A[a][b].is_zero(): U.append(None); continue
            U.append(-(D * A[a][b].n * (La // A[a][b].d)))
        Uv.append(U)
    maxdeg = max(Vv[a].degree() for a in range(nb)) + dN + 2
    ncols = nb * (dN + 1)
    rowsM = {}

    def addpoly(a, P, col):
        for kk, c in enumerate(P.coeffs()):
            if c == 0: continue
            key = a * (maxdeg + 1) + kk
            rowsM.setdefault(key, {})
            rowsM[key][col] = rowsM[key].get(col, FQ(0)) + c

    tpow = [FP([0] * k2 + [1]) for k2 in range(dN + 2)]
    for b in range(nb):
        for m2 in range(dN + 1):
            col = b * (dN + 1) + m2
            for a in range(nb):
                acc = ZERO_P
                if Uv[a][b] is not None: acc = acc + Uv[a][b] * tpow[m2]
                if a == b:
                    if m2 > 0: acc = acc + Vv[a] * FP([fq(m2)]) * tpow[m2 - 1]
                    acc = acc + Wv[a] * tpow[m2]
                if acc != ZERO_P: addpoly(a, acc, col)
    rowkeys = sorted(rowsM)
    Mz = flint.fmpq_mat(len(rowkeys), ncols)
    for ri, rk in enumerate(rowkeys):
        for col, c in rowsM[rk].items(): Mz[ri, col] = c
    R2, rank = Mz.rref()
    piv = []; ci = 0
    for ri in range(rank):
        while ci < ncols and R2[ri, ci] == 0: ci += 1
        piv.append(ci)
    free = [c2 for c2 in range(ncols) if c2 not in piv]
    sols = []
    for fcol in free:
        vec = [FQ(0)] * ncols
        vec[fcol] = FQ(1)
        for ri in range(rank - 1, -1, -1):
            pc = piv[ri]
            s = FQ(0)
            for c2 in range(pc + 1, ncols):
                if R2[ri, c2] != 0 and vec[c2] != 0: s += R2[ri, c2] * vec[c2]
            vec[pc] = -s / R2[ri, pc]
        Nvecs = [Rat(fp(vec[b * (dN + 1):(b + 1) * (dN + 1)]), D) for b in range(nb)]
        if all(x.is_zero() for x in Nvecs): continue
        sols.append(Nvecs)
    return sols

# ─────────────────── block strategies (FULL / TRI / CONJ) ───────────────────
def rat_matinv(Rm):
    """Exact inverse + determinant of an n x n Rat matrix (Gauss-Jordan).
    Raises ZeroDivisionError on a singular matrix."""
    nb = len(Rm)
    A = [list(Rm[i]) + [R_ONE if i == j else R_ZERO for j in range(nb)] for i in range(nb)]
    det = R_ONE
    for c in range(nb):
        p = next((r for r in range(c, nb) if not A[r][c].is_zero()), None)
        if p is None: raise ZeroDivisionError("singular basis matrix")
        if p != c: A[c], A[p] = A[p], A[c]; det = -det
        piv = A[c][c]; det = det * piv
        pi = piv.inv()
        A[c] = [x * pi for x in A[c]]
        for r in range(nb):
            if r == c or A[r][c].is_zero(): continue
            f = A[r][c]
            A[r] = [A[r][j] - f * A[c][j] for j in range(2 * nb)]
    return [row[nb:] for row in A], det

def matmul_rat(X, Y):
    n1 = len(X); n2 = len(Y[0])
    out = [[R_ZERO for _ in range(n2)] for _ in range(n1)]
    for a in range(n1):
        for b in range(n2):
            acc = R_ZERO
            for c in range(len(Y)):
                if X[a][c].is_zero() or Y[c][b].is_zero(): continue
                acc = acc + X[a][c] * Y[c][b]
            out[a][b] = acc
    return out

def eig_generic(Amat):
    """Numeric residue-matrix eigenvalues of Amat at one root of factor f
    (denominator-bound helper for solve_rat_solutions; never load-bearing)."""
    def ef(f):
        rts = mp.polyroots([mp.mpf(int(c.p)) / int(c.q) for c in reversed(f.coeffs())],
                           maxsteps=300, extraprec=80)
        r0 = rts[0]
        nb = len(Amat)
        Mm = mp.matrix(nb, nb)
        K = 8; eps = mp.mpf('1e-18')
        for a in range(nb):
            for b in range(nb):
                if Amat[a][b].is_zero(): continue
                s = mp.mpc(0)
                for kk in range(K):
                    z = r0 + eps * mp.exp(2j * mp.pi * kk / K)
                    s += (z - r0) * Amat[a][b].ev_mpc(z)
                Mm[a, b] = s / K
        try:
            return mp.eig(Mm, left=False, right=False)
        except Exception:
            return []
    return ef

def certify_generic(Amat, vec, rad):
    """EXACT check: (vec*sqrt(rad))' = Amat (vec*sqrt(rad)) as Rat identity."""
    nb = len(Amat)
    half = rad_dlog_half(rad) if rad else R_ZERO
    for a in range(nb):
        resid = vec[a].deriv() + vec[a] * half
        for b in range(nb):
            if Amat[a][b].is_zero(): continue
            resid = resid - Amat[a][b] * vec[b]
        if not resid.is_zero(): return False
    return True

def scalar_strategy(b):
    if b.is_zero(): return ('FULL', [[R_ONE]], [()])
    H, logs, rems = hermite(b)
    if not H.is_zero() or rems:
        raise RuntimeError(f"scalar mu non-elementary (H={not H.is_zero()}, rems={len(rems)})")
    num = ONE_P; den = ONE_P; rad = []
    for li, cl in logs.items():
        c2 = Fraction(int(cl.p), int(cl.q))
        if c2.denominator == 1:
            if c2 >= 0: num = num * LETTERS[li] ** c2.numerator
            else: den = den * LETTERS[li] ** (-c2.numerator)
        elif c2.denominator == 2:
            fl = (c2.numerator - 1) // 2
            if fl >= 0: num = num * LETTERS[li] ** fl
            else: den = den * LETTERS[li] ** (-fl)
            rad.append(li)
        else: raise RuntimeError(f"scalar non-half-integer exponent {c2} letter {li}")
    col = Rat(num, den)
    rad = tuple(sorted(rad))
    assert certify_generic([[b]], [col], rad), "scalar cert fail"
    return ('FULL', [[col]], [rad])

def find_solutions(Amat):
    """Certified homogeneous columns: rational first, then radical (sqrt of the
    letters whose residue matrix carries half-integer eigenvalues)."""
    nb = len(Amat)
    facs = {}
    for a in range(nb):
        for b in range(nb):
            if Amat[a][b].is_zero(): continue
            c, fl = Amat[a][b].d.factor()
            for f, e in fl: facs[tuple(str(x) for x in f.coeffs())] = f
    ef = eig_generic(Amat)
    halfset = []
    for k2, f in facs.items():
        for x in ef(f):
            fr = float(mp.re(x)) * 2
            if abs(fr - round(fr)) < 1e-5 and round(fr) % 2 != 0: halfset.append(f); break
    cols = []; rads = []
    for v in solve_rat_solutions(Amat, None, eig_fn=ef):
        if len(cols) < nb and certify_generic(Amat, v, ()): cols.append(v); rads.append(())
    if len(cols) < nb and halfset:
        radli = tuple(sorted(reg_letter(f) for f in halfset))
        for v in solve_rat_solutions(Amat, radli, eig_fn=ef):
            if len(cols) < nb and certify_generic(Amat, v, radli):
                cols.append(v); rads.append(radli)
    return cols, rads

def strategy_for(Amat, tag=''):
    """Exact solve strategy for Y' = Amat Y + S.
    ('FULL', cols, rads): full homogeneous basis (col c = sqrt(rad_c)*Rat-vec).
    ('CONJ', radli, sub): conjugate by the radical, Y = sqrt(g/g0) Z.
    ('TRI', T, Tinv, B, f, sub): T = [rational solutions | unit cols],
      B = Tinv(AT - T'), first f columns of B exactly zero; sub = corner strategy."""
    nb = len(Amat)
    if all(Amat[a][b].is_zero() for a in range(nb) for b in range(nb)):
        return ('FULL', [[R_ONE if a == c else R_ZERO for a in range(nb)]
                         for c in range(nb)], [()] * nb)
    if nb == 1: return scalar_strategy(Amat[0][0])
    cols, rads = find_solutions(Amat)
    if len(cols) == nb: return ('FULL', cols, rads)
    ratcols = [cols[c] for c in range(len(cols)) if rads[c] == ()]
    if not ratcols:
        radlis = [rads[c] for c in range(len(cols)) if rads[c] != ()]
        if radlis:
            radli = radlis[0]
            g = ONE_P
            for li in radli: g = g * LETTERS[li]
            half = Rat(g.derivative(), g * FP([2]))
            Ag = [[(Amat[a][b] - half if a == b else Amat[a][b]) for b in range(nb)]
                  for a in range(nb)]
            return ('CONJ', radli, strategy_for(Ag, tag + "/conj"))
        raise RuntimeError(f"{tag}: no rational homogeneous column for triangularization "
                           f"({len(cols)}/{nb} found)")
    f = len(ratcols)
    T = Tinv = None
    for unit_idx in combinations(range(nb), nb - f):
        Tc = [[(ratcols[c][a] if c < f else (R_ONE if a == unit_idx[c - f] else R_ZERO))
               for c in range(nb)] for a in range(nb)]
        try:
            Ti, det = rat_matinv(Tc)
            if det.is_zero() or det.ev_fq(0) == 0: continue
            T, Tinv = Tc, Ti
            break
        except ZeroDivisionError:
            continue
    if T is None: raise RuntimeError(f"{tag}: no invertible T padding found")
    Tp = [[T[a][c].deriv() for c in range(nb)] for a in range(nb)]
    AT = matmul_rat(Amat, T)
    B = matmul_rat(Tinv, [[AT[a][c] - Tp[a][c] for c in range(nb)] for a in range(nb)])
    for c in range(f):
        for a in range(nb):
            assert B[a][c].is_zero(), f"{tag}: B col {c} not zero — bad rational solution"
    corner = [[B[f + a][f + b] for b in range(nb - f)] for a in range(nb - f)]
    return ('TRI', T, Tinv, B, f, strategy_for(corner, tag + f"/corner{nb-f}"))

def solve_lin(strat, Svec, w0cf):
    """Containers Y solving Y' = A Y + S with Y(0) given by containers w0cf."""
    if strat[0] == 'CONJ':
        _, radli, sub = strat
        invfac = R_ONE
        for li in radli:
            fpoly = LETTERS[li]
            invfac = invfac * Rat(FP([fpoly(fq(0))]), fpoly)
        Smod = [cf_mulalg(S, invfac, radli) if S else cf_zero() for S in Svec]
        Z = solve_lin(sub, Smod, w0cf)
        return [cf_mulalg(z, R_ONE, radli) if z else cf_zero() for z in Z]
    if strat[0] == 'FULL':
        _, cols, rads = strat
        nb = len(cols)
        Rm = [[cols[c][a] for c in range(nb)] for a in range(nb)]
        Rinv, det = rat_matinv(Rm)
        Ys = [cf_zero() for _ in range(nb)]
        for c in range(nb):
            Tc = cf_zero()
            for b in range(nb):
                if Rinv[c][b].is_zero() or not Svec[b]: continue
                Tc = cf_add(Tc, cf_mulrat(Svec[b], Rinv[c][b]))
            if rads[c] and Tc:
                invfac = R_ONE
                for li in rads[c]:
                    fpoly = LETTERS[li]
                    invfac = invfac * Rat(FP([fpoly(fq(0))]), fpoly)
                Tc = cf_mulalg(Tc, invfac, rads[c])
            Ic = cf_int(Tc) if Tc else cf_zero()
            C0 = cf_zero()
            for b in range(nb):
                if Rinv[c][b].is_zero() or not w0cf[b]: continue
                rv = Rinv[c][b].ev_fq(0)
                if rv != 0: C0 = cf_add(C0, cf_scale(w0cf[b], rv))
            inner = cf_add(C0, Ic)
            if not inner: continue
            for a in range(nb):
                if Rm[a][c].is_zero(): continue
                Ys[a] = cf_add(Ys[a], cf_mulalg(inner, Rm[a][c], rads[c]))
        return Ys
    _, T, Tinv, B, f, sub = strat
    nb = len(T)
    TS = []; W0 = []
    for c in range(nb):
        acc = cf_zero()
        for b in range(nb):
            if Tinv[c][b].is_zero() or not Svec[b]: continue
            acc = cf_add(acc, cf_mulrat(Svec[b], Tinv[c][b]))
        TS.append(acc)
        acc0 = cf_zero()
        for b in range(nb):
            if Tinv[c][b].is_zero() or not w0cf[b]: continue
            rv = Tinv[c][b].ev_fq(0)
            if rv != 0: acc0 = cf_add(acc0, cf_scale(w0cf[b], rv))
        W0.append(acc0)
    wsub = solve_lin(sub, TS[f:], W0[f:])
    w = [None] * nb
    for c in range(nb - f): w[f + c] = wsub[c]
    for a in range(f):
        integ = TS[a]
        for b in range(f, nb):
            if B[a][b].is_zero() or not w[b]: continue
            integ = cf_add(integ, cf_mulrat(w[b], B[a][b]))
        w[a] = cf_add(W0[a], cf_int(integ) if integ else cf_zero())
    Ys = [cf_zero() for _ in range(nb)]
    for a in range(nb):
        for c in range(nb):
            if T[a][c].is_zero() or not w[c]: continue
            Ys[a] = cf_add(Ys[a], cf_mulrat(w[c], T[a][c]))
    return Ys

def det_roots_of_strategies(strats, lo=1e-6, hi=None):
    """APPARENT poles of the VoP representation: real zeros on (lo,hi) of
    det(fundamental basis) / det(T) for every strategy in the tree. C'=M^{-1}S
    has poles there that cancel only in Y=MC — the contour MUST detour them
    (integrating through these points is an honest gate FAIL — they are
    excluded by construction)."""
    hi = 1 - lo if hi is None else hi
    polys = []

    def walk(st):
        if st[0] == 'FULL':
            cols = st[1]; nb = len(cols)
            try:
                _, det = rat_matinv([[cols[c][a] for c in range(nb)] for a in range(nb)])
            except ZeroDivisionError:
                return
            if det.n.degree() >= 1: polys.append(det.n)
        elif st[0] == 'CONJ': walk(st[2])
        else:
            try:
                _, det = rat_matinv(st[1])
                if det.n.degree() >= 1: polys.append(det.n)
            except ZeroDivisionError:
                pass
            walk(st[5])

    for st in strats: walk(st)
    roots = set()
    for p in polys:
        cs = [mp.mpf(int(c.p)) / int(c.q) for c in reversed(p.coeffs())]
        try:
            rr = mp.polyroots(cs, maxsteps=300, extraprec=80)
        except Exception:
            continue
        for r in rr:
            if abs(mp.im(r)) < 1e-9 and lo < mp.re(r) < hi: roots.add(round(float(mp.re(r)), 8))
    return sorted(roots)

# ─────────────────── numeric spectral contour engine ───────────────────
def _mpf2arb(x):
    s, m, e, _ = mp.mpf(x)._mpf_
    v = flint.arb(-int(m) if s else int(m))
    return v * (flint.arb(2) ** int(e)) if e else v

def _mpc2acb(z):
    z = mp.mpc(z); return flint.acb(_mpf2arb(z.real), _mpf2arb(z.imag))

def _acb2mpc(z):
    def f(a):
        m, e = a.man_exp(); return mp.mpf(int(m)) * mp.mpf(2) ** int(e)
    return mp.mpc(f(z.real.mid()), f(z.imag.mid()))

class Engine:
    """Chebyshev spectral cumulative-integration engine on the detour contour
    0 -> 1 in the upper half plane (real poles AND det-basis apparent poles
    passed in real_poles are detoured at height delta)."""

    def __init__(self, real_poles, NC=120, delta='0.08', pad='0.04'):
        self.NC = NC
        flint.ctx.prec = int(3.33 * mp.mp.dps) + 20
        delta = mp.mpf(delta); pad = mp.mpf(pad)
        rp = sorted(set(round(float(x), 10) for x in real_poles if 1e-8 < float(x) < 1 - 1e-8))
        clusters = []
        for x in rp:
            if clusters and x - clusters[-1][1] < 0.08: clusters[-1] = (clusters[-1][0], x)
            else: clusters.append((x, x))
        verts = [mp.mpc(0)]
        for (a, b) in clusters:
            a_ = max(mp.mpf(a) - pad, mp.mpf(float(verts[-1].real)) + mp.mpf('0.005'))
            b_ = min(mp.mpf(b) + pad, mp.mpf(1) - mp.mpf('0.005'))
            if a_ >= b_: continue
            verts += [mp.mpc(a_), mp.mpc(a_, delta), mp.mpc(b_, delta), mp.mpc(b_)]
        verts.append(mp.mpc(1))
        self.segs = [(verts[k], verts[k + 1]) for k in range(len(verts) - 1)]
        self.uch = [(1 - mp.cos(mp.pi * k / NC)) / 2 for k in range(NC + 1)]
        self._COS = [[mp.cos(mp.pi * m * k / NC) for k in range(NC + 1)] for m in range(NC + 2)]
        self.zg = [[za + u * (zb - za) for u in self.uch] for (za, zb) in self.segs]
        self.zga = [[_mpc2acb(z) for z in row] for row in self.zg]
        self.dz = [zb - za for (za, zb) in self.segs]
        self.dza = [_mpc2acb(d) for d in self.dz]
        self._wcache = {}; self._radcache = {}; self._kcache = {}
        N = NC
        Cm = flint.acb_mat(N + 1, N + 1)
        for m2 in range(N + 1):
            for k in range(N + 1):
                w = mp.mpf(1) / 2 if k in (0, N) else mp.mpf(1)
                Cm[m2, k] = _mpc2acb(2 * w * self._COS[m2][k] / N)
        Wm = flint.acb_mat(N + 1, N + 1)
        for k in range(N + 1):
            Wm[k, 0] = _mpc2acb(self.uch[k] / 2)
            Wm[k, 1] = _mpc2acb((1 - self._COS[2][k]) / 8)
            for m2 in range(2, N + 1):
                v = ((1 - self._COS[m2 + 1][k]) / (m2 + 1) - (1 - self._COS[m2 - 1][k]) / (m2 - 1)) / 4
                if m2 == N: v = v / 2
                Wm[k, m2] = _mpc2acb(v)
        self._CUM = Wm * Cm

    def _cumcheb_acb(self, vals_acb):
        N = self.NC
        V = flint.acb_mat(N + 1, 1)
        for k in range(N + 1): V[k, 0] = vals_acb[k]
        O = self._CUM * V
        return [O[k, 0] for k in range(N + 1)]

    def kern_vals(self, kk):
        if kk in self._kcache: return self._kcache[kk]
        K, radk = kern_eval_parts(kk)
        rv = self.rad_vals(radk) if radk else None
        ncf = [_mpc2acb(mp.mpf(int(c.p)) / int(c.q)) for c in reversed(K.n.coeffs())]
        dcf = [_mpc2acb(mp.mpf(int(c.p)) / int(c.q)) for c in reversed(K.d.coeffs())]
        out = []
        for s in range(len(self.segs)):
            row = []
            for k in range(self.NC + 1):
                z = self.zga[s][k]
                num = flint.acb(0)
                for c in ncf: num = num * z + c
                den = flint.acb(0)
                for c in dcf: den = den * z + c
                v = num / den
                if rv is not None: v = v * _mpc2acb(rv[s][k])
                row.append(v)
            out.append(row)
        self._kcache[kk] = out
        return out

    def rad_vals(self, rad):
        """Continued sqrt(prod l(z)/l(0)) at all contour nodes."""
        if rad in self._radcache: return self._radcache[rad]
        w = []
        for s in range(len(self.segs)):
            row = []
            for k in range(self.NC + 1):
                z = self.zg[s][k]
                v = mp.mpc(1)
                for li in rad:
                    f = LETTERS[li]
                    num = mp.mpc(0)
                    for c in reversed(f.coeffs()): num = num * z + mp.mpf(int(c.p)) / int(c.q)
                    f0 = f(fq(0))
                    v *= num / (mp.mpf(int(f0.p)) / int(f0.q))
                row.append(v)
            w.append(row)
        out = []; prev = mp.mpc(1); prevw = mp.mpc(1)
        for s in range(len(self.segs)):
            row = []
            for k in range(self.NC + 1):
                cur = prev * mp.sqrt(w[s][k] / prevw)
                row.append(cur); prev = cur; prevw = w[s][k]
            out.append(row)
        self._radcache[rad] = out
        return out

    def word_vals_acb(self, word):
        if word in self._wcache: return self._wcache[word]
        if word == ():
            one = flint.acb(1)
            sv = [[one] * (self.NC + 1) for _ in self.segs]
            self._wcache[word] = sv
            return sv
        inner = self.word_vals_acb(word[1:])
        kv = self.kern_vals(word[0])
        sv = []; cum = flint.acb(0)
        for s in range(len(self.segs)):
            dz = self.dza[s]
            g = [kv[s][k] * inner[s][k] * dz for k in range(self.NC + 1)]
            segcum = self._cumcheb_acb(g)
            sv.append([cum + segcum[k] for k in range(self.NC + 1)])
            cum += segcum[-1]
        self._wcache[word] = sv
        return sv

    def word_vals(self, word):
        key = ('mpc', word)
        if key in self._wcache: return self._wcache[key]
        out = [[_acb2mpc(x) for x in row] for row in self.word_vals_acb(word)]
        self._wcache[key] = out
        return out

    def locate(self, z):
        for s, (za, zb) in enumerate(self.segs):
            d = self.dz[s]
            u = ((z - za) / d).real
            offax = abs(za + u * d - z)
            if -1e-12 <= u <= 1 + 1e-12 and offax < 1e-10: return s, mp.mpf(u)
        raise ValueError(f"point {z} not on contour")

    def interp(self, vals, u):
        num = mp.mpc(0); den = mp.mpc(0)
        for k in range(self.NC + 1):
            d = u - self.uch[k]
            if abs(d) < mp.mpf(10) ** -60: return vals[k]
            w = (-1) ** k * (mp.mpf(1) / 2 if k in (0, self.NC) else mp.mpf(1)) / d
            num += w * vals[k]; den += w
        return num / den

    def eval_cf(self, A, z):
        z = mp.mpc(z)
        s, u = self.locate(z)
        tot = mp.mpc(0)
        for (ctag, rad, word), r in A.items():
            gv = self.interp(self.word_vals(word)[s], u) if word else mp.mpc(1)
            rv = self.interp(self.rad_vals(rad)[s], u) if rad else mp.mpc(1)
            tot += CONSTS[ctag] * r.ev_mpc(z) * gv * rv
        return tot

    def eval_cf_end(self, A): return self.eval_cf(A, mp.mpc(1))

# ─────────────────── boundary constant closing (PSLQ) ───────────────────
_LOGB_CACHE = {}

def _log_basis(primes):
    key = (mp.mp.dps, tuple(primes))
    if key not in _LOGB_CACHE: _LOGB_CACHE[key] = [mp.log(p) for p in primes]
    return _LOGB_CACHE[key]

def try_rational(x, maxc=10 ** 15, tol=None):
    tol = tol or mp.mpf(10) ** -52
    if abs(mp.im(x)) > tol: return None
    if mp.re(x) == 0: return Fraction(0)   # mp.pslq rejects zero vectors (exactly-real/zero inputs)
    r = mp.pslq([mp.re(x), mp.mpf(1)], maxcoeff=maxc, maxsteps=10 ** 6)
    if not r or r[0] == 0: return None
    v = Fraction(-r[1], r[0])
    if abs(mp.re(x) - mp.mpf(v.numerator) / v.denominator) < tol: return v
    return None

def try_logring(x, primes=(2, 3, 5, 7, 11, 13, 17), tol=None):
    """x = q0 + sum q_p log(p) + q_pi*i*pi with small rational q's, or None."""
    tol = tol or mp.mpf(10) ** -48
    lb = _log_basis(primes)
    vec = [mp.re(x)] + lb + [mp.mpf(1)]
    r = mp.pslq(vec, maxcoeff=10 ** 4, maxsteps=10 ** 6)
    out = None
    if r and r[0] != 0:
        den = -r[0]
        co = {f"log{primes[a]}": Fraction(r[a + 1], den) for a in range(len(primes))}
        co['1'] = Fraction(r[len(primes) + 1], den)
        chk = sum((mp.mpf(f.numerator) / f.denominator) * (lb[primes.index(int(t2[3:]))] if t2 != '1' else 1)
                  for t2, f in co.items())
        if abs(mp.re(x) - chk) < tol: out = co
    if out is None: return None
    qpi = try_rational(mp.im(x) / mp.pi, maxc=10 ** 6)
    if qpi is None and abs(mp.im(x)) > tol: return None
    if qpi: out['ipi'] = qpi
    return {t2: f for t2, f in out.items() if f != 0}

# ─────────────────── the closure driver ───────────────────
@dataclass
class VopConfig:
    """Sector-agnostic configuration (everything sector-specific lives here).
    RatA:          {(i,j,m): Rat} exact graded connection A^{(m)}_{ij}(t)
    blocks:        {block_id: [row indices]} diagonal blocks of the system
    boundary_vals: {(i,k): mpc} values at t=0 (layer k = graded eps-layer)
    kmin,kmax:     layer range to close (e.g. -4..0)
    L:             optional per-row grading offsets (eps-power = k - L[i]),
                   bookkeeping only
    closed_sources:optional (j,k) -> container for rows NOT in blocks (e.g.
                   sunrise-type subsectors with a priori closed forms)
    top_rows:      rows reported by certify()/gate() by default
    scope_rows:    if set, the LAST layer only solves the upward closure of
                   these rows (a production top-row optimization)
    poles:         known real path poles in (0,1) (contour detours); det-basis
                   apparent poles are added automatically by engine()
    NC:            Chebyshev nodes per contour segment
    boundary_mode: 'pslq' (rational -> log-ring -> tagged const) or 'tag'
    log_primes:    primes for the log-ring boundary basket"""
    RatA: dict
    blocks: dict
    boundary_vals: dict
    kmin: int
    kmax: int
    L: dict = None
    closed_sources: object = None
    top_rows: tuple = ()
    scope_rows: tuple = ()
    poles: tuple = ()
    NC: int = 140
    boundary_mode: str = 'pslq'
    log_primes: tuple = (2, 3, 5, 7, 11, 13, 17)

class VopClosure:
    """Layered VoP closure of a graded 1-D path-DE. Usage:
        vc = VopClosure(cfg); NCF = vc.close()
        cert = vc.certify()          # exact zero-remainder DE certificates
        eng  = vc.engine()           # spectral engine on the detoured contour
        gate = vc.gate(refs)         # held-out endpoint digits vs refs"""

    def __init__(self, cfg: VopConfig):
        self.cfg = cfg
        self.COLS = {}; self.SUPP = {}
        for (i, j, m) in cfg.RatA:
            self.COLS.setdefault(i, set()).add(j)
            self.SUPP.setdefault((i, j), []).append(m)
        self.row_block = {i: s for s, rows in cfg.blocks.items() for i in rows}
        self.block_rows = set(self.row_block)
        self.BAS_CACHE = {}
        self.NCF = {}
        self.boundary_record = {}
        self._eng = None

    def topo_order(self):
        edges = {s: set() for s in self.cfg.blocks}
        for (i, j) in self.SUPP:
            si = self.row_block.get(i); sj = self.row_block.get(j)
            if si is None or sj is None or si == sj: continue
            edges[si].add(sj)
        order = []; seen = set()

        def visit(s):
            if s in seen: return
            seen.add(s)
            for p in sorted(edges[s]): visit(p)
            order.append(s)

        for s in sorted(self.cfg.blocks): visit(s)
        return order

    def block_strategy(self, s):
        if s in self.BAS_CACHE: return self.BAS_CACHE[s]
        rows = self.cfg.blocks[s]
        nb = len(rows)
        Amat = [[self.cfg.RatA.get((rows[a], rows[b], 0), R_ZERO) for b in range(nb)]
                for a in range(nb)]
        st = strategy_for(Amat, tag=f"blk{s}{rows}")
        self.BAS_CACHE[s] = st
        return st

    def upward_closure(self, rows):
        cl = set(); stack = list(rows)
        while stack:
            i = stack.pop()
            if i in cl: continue
            cl.add(i)
            for j in self.COLS.get(i, ()):
                if j in self.block_rows and j not in cl: stack.append(j)
        return cl

    def getn(self, j, k):
        if j in self.block_rows: return self.NCF.get((j, k), cf_zero())
        if self.cfg.closed_sources is not None: return self.cfg.closed_sources(j, k) or cf_zero()
        return cf_zero()

    def _boundary_cf(self, i, k):
        x = self.cfg.boundary_vals.get((i, k))
        if x is None or abs(x) < mp.mpf(10) ** -65: return cf_zero()
        if self.cfg.boundary_mode == 'pslq':
            q = try_rational(x)
            if q is not None:
                self.boundary_record[(i, k)] = f"exact {q}"
                return {('1', (), ()): Rat(fp([q]))}
            lr = try_logring(x, self.cfg.log_primes)
            if lr is not None:
                C = cf_zero()
                for t2, f in lr.items():
                    if t2 == 'ipi': reg_const('ipi', mp.pi * 1j)
                    elif t2 != '1': reg_const(t2, mp.mpc(mp.log(int(t2[3:]))))
                    C = cf_add(C, {(t2, (), ()): Rat(fp([f]))})
                self.boundary_record[(i, k)] = "logring " + str(lr)
                return C
        tag = reg_const(f"B{i}k{k}", mp.mpc(x))
        self.boundary_record[(i, k)] = f"tagged {mp.nstr(mp.mpc(x), 35)}"
        return cf_const(tag)

    def close(self, log=None):
        order = self.topo_order()
        scope = self.upward_closure(self.cfg.scope_rows) if self.cfg.scope_rows else None
        for k in range(self.cfg.kmin, self.cfg.kmax + 1):
            t0 = time.time()
            for s in order:
                rows = self.cfg.blocks[s]
                if scope is not None and k == self.cfg.kmax and not any(i in scope for i in rows):
                    continue
                Svec = []
                for i in rows:
                    S = cf_zero()
                    for j in self.COLS.get(i, ()):
                        for m in self.SUPP.get((i, j), ()):
                            if m == 0 and j in rows: continue
                            src = self.getn(j, k - m)
                            if not src: continue
                            S = cf_add(S, cf_mulrat(src, self.cfg.RatA[(i, j, m)]))
                    Svec.append(S)
                v0cf = [self._boundary_cf(i, k) for i in rows]
                if all(not c for c in v0cf) and all(not S for S in Svec): continue
                Ys = solve_lin(self.block_strategy(s), Svec, v0cf)
                for a, i in enumerate(rows):
                    if Ys[a]: self.NCF[(i, k)] = Ys[a]
            if log:
                nnz = sum(1 for (i2, k2) in self.NCF if k2 == k)
                log(f"[vopclose] layer {k}: {nnz} nonzero rows ({time.time()-t0:.1f}s)")
        return self.NCF

    def certify(self, rows=None, layers=None):
        """EXACT zero-remainder symbolic DE check per closed (row, layer):
        cf_diff(Y_{i,k}) - sum_{j,m} A^{(m)}_{ij} y_{j,k-m} == 0 as containers."""
        rows = list(rows if rows is not None else (self.cfg.top_rows or sorted(self.block_rows)))
        layers = list(layers if layers is not None else range(self.cfg.kmin, self.cfg.kmax + 1))
        out = {}
        for i in rows:
            for k in layers:
                Y = self.NCF.get((i, k))
                if Y is None: continue
                D = cf_diff(Y)
                for j in self.COLS.get(i, ()):
                    for m in self.SUPP.get((i, j), ()):
                        src = self.getn(j, k - m)
                        if not src: continue
                        D = cf_add(D, cf_neg(cf_mulrat(src, self.cfg.RatA[(i, j, m)])))
                bad = [kk for kk, r in D.items() if not r.is_zero()]
                out[(i, k)] = (len(bad) == 0, len(bad))
        return out

    def detour_roots(self):
        return det_roots_of_strategies([self.BAS_CACHE[s] for s in sorted(self.BAS_CACHE)])

    def engine(self, NC=None):
        """Spectral engine on the contour detouring BOTH the configured path
        poles and the det-basis apparent poles of the strategies (built-in)."""
        if self._eng is None:
            pts = sorted(set(float(p) for p in self.cfg.poles) | set(self.detour_roots()))
            self._eng = Engine(pts, NC=NC or self.cfg.NC)
        return self._eng

    def gate(self, refs, z=None):
        """Held-out gate: refs {(i,k): mpc reference at path end (default t=1)}.
        Returns {(i,k): relative agreed digits}; read the held-out oracle ONLY here."""
        eng = self.engine()
        z = mp.mpc(1) if z is None else mp.mpc(z)
        out = {}
        for (i, k), ref in refs.items():
            A = self.getn(i, k)
            if not A: continue
            v = eng.eval_cf(A, z)
            dd = abs(v - ref)
            out[(i, k)] = 75.0 if dd == 0 else float(-mp.log10(dd / abs(ref)))
        return out

# ─────────────── exact connection entries from path samples (via ratfit) ───────────────
def exact_entry_from_samples(xs, ys, Qknown=None, n_loo=6):
    """One graded connection entry a_m(t) as an exact Rat from exact node
    samples, delegating to tools/ratfit: known-Q route (GF(p) screen +
    exact Newton certify) when a validated denominator is supplied, else
    blind screened Thiele-LOO. Returns (Rat, ok, deg). This is the
    sector-agnostic core of the production exact-A build route."""
    try:
        import ratfit
    except ImportError:
        # vendored layout: ratfit is the sibling package tools/ratfit/
        import os as _os, sys as _sys
        _sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
        import ratfit
    if Qknown is not None:
        if ratfit.screen_modp(xs, ys, Qknown, n_loo=n_loo):
            P, Q, ok, deg = ratfit.ratrecon_with_denom(xs, ys, Qknown, n_loo=n_loo)
            if ok: return Rat(P, Q), True, deg
    P, Q, ok, deg = ratfit.thiele_loo_screened(xs, ys, n_loo=n_loo)
    return Rat(P, Q), ok, deg

# ─────────────────────────── self-test ───────────────────────────
def _selftest():
    mp.mp.dps = 60
    t0 = time.time()
    results = []

    def check(name, cond):
        results.append((name, bool(cond)))
        print(f"  [{'PASS' if cond else 'FAIL'}] {name}", flush=True)

    def digits(v, ref):
        dd = abs(v - ref)
        return 75.0 if dd == 0 else float(-mp.log10(dd / abs(ref)))

    reset_registries()
    # 1: rational int/diff round trip + endpoint value
    A = {('1', (), ()): Rat(fp([1]), fp([2, 1]))}
    I = cf_int(A)
    dr = cf_add(cf_diff(I), cf_neg(A))
    check("T1 exact d/dt(int) - id zero-remainder", all(v.is_zero() for v in dr.values()))
    eng = Engine([], NC=96)
    check("T2 spectral endpoint log(3/2) >=45d",
          digits(eng.eval_cf_end(I), mp.log(mp.mpf(3) / 2)) >= 45)
    # 3: nested word vs quadrature
    I2 = cf_int(cf_mulrat(I, Rat(fp([1]), fp([1, 1]))))
    ref2 = mp.quad(lambda s: mp.log((s + 2) / 2) / (s + 1), [0, 1])
    check("T3 nested iterated integral >=40d", digits(eng.eval_cf_end(I2), ref2) >= 40)
    # 4: multi-factor Hermite (repeated + irreducible quadratic)
    r3 = Rat(fp([1, 0, 3]), fp([2, 1]) ** 2 * fp([5, 1, 1]))
    I3 = cf_int({('1', (), ()): r3})
    dr3 = cf_add(cf_diff(I3), cf_neg({('1', (), ()): r3}))
    ref3 = mp.quad(lambda s: (3 * s ** 2 + 1) / ((s + 2) ** 2 * (s ** 2 + s + 5)), [0, 1])
    check("T4 Hermite multi-factor DE + value >=40d",
          all(v.is_zero() for v in dr3.values()) and digits(eng.eval_cf_end(I3), ref3) >= 40)
    # 5: radical algebra + continuation (linear and quartic letters)
    li = reg_letter(fp([3, 1]))
    Y = {('1', (li,), ()): R_ONE}
    dr4 = cf_add(cf_diff(Y), cf_neg(cf_mulalg(Y, Rat(fp([1]), fp([6, 2])), ())))
    li4 = reg_letter(fp([1, 0, 0, 0, 1]))
    Y4 = {('1', (li4,), ()): R_ONE}
    check("T5 radical d/dt exact + sqrt values (linear 50d, quartic 50d)",
          all(v.is_zero() for v in dr4.values())
          and digits(eng.eval_cf(Y, mp.mpc(1)), mp.sqrt(mp.mpf(4) / 3)) >= 50
          and digits(eng.eval_cf(Y4, mp.mpc(1)), mp.sqrt(2)) >= 50)
    # 6: algebraic-kernel integration
    I5 = cf_int(cf_mulrat(Y, Rat(fp([1]), fp([1, 1]))))
    ref5 = mp.quad(lambda s: mp.sqrt((s + 3) / 3) / (s + 1), [0, 1])
    check("T6 algebraic kernel int >=40d", digits(eng.eval_cf_end(I5), ref5) >= 40)
    # 7: negative-letter radical branch on the detour contour
    li2 = reg_letter(fp([Fraction(2, 5), -1]))
    Y2 = {('1', (li2,), ()): R_ONE}
    eng2 = Engine([0.4], NC=96)
    ref6 = mp.sqrt(mp.mpf(3) / 2) * mp.exp(-1j * mp.pi / 2)
    check("T7 branch continuation through detour >=45d",
          digits(eng2.eval_cf(Y2, mp.mpc(1)), ref6) >= 45)
    # 8: det-basis apparent-pole detection (root at t=1/2)
    st = ('FULL', [[R_ONE, Rat(fp([0, 2]))], [Rat(fp([0, 2])), R_ONE]], [(), ()])
    check("T8 det-basis apparent-pole detection", det_roots_of_strategies([st]) == [0.5])

    # 9-14: graded 3x3 toy closure (rational block + 2x2 radical block with
    # half-integer indicials at t=-3; exercises FULL/CONJ/TRI + layer recursion)
    reset_registries()
    Rr = lambda n, d=None: Rat(fp(n), fp(d) if d is not None else ONE_P)
    RatA = {(0, 0, 0): Rr([1], [2, 1]),          # 1/(t+2)
            (0, 0, 1): Rr([1], [2, 1]),
            (1, 1, 0): Rr([1], [6, 2]),          # 1/(2(t+3)) — indicial 1/2
            (1, 2, 0): Rr([1], [1, 1]),          # 1/(t+1)
            (2, 2, 0): Rr([1], [6, 2]),
            (2, 0, 1): R_ONE}
    cfg = VopConfig(RatA=RatA, blocks={0: [0], 1: [1, 2]},
                    boundary_vals={(0, 0): mp.mpc(1), (1, 0): mp.mpc(1), (2, 0): mp.mpc(2)},
                    kmin=0, kmax=1, top_rows=(0, 1, 2), NC=96)
    vc = VopClosure(cfg)
    NCF = vc.close()
    check("T9 toy closure: 6 nonzero (row,layer) containers", len(NCF) == 6)
    st1 = vc.block_strategy(1)
    check("T10 radical block strategy = CONJ(sqrt letter) -> TRI",
          st1[0] == 'CONJ' and st1[2][0] == 'TRI')
    cert = vc.certify()
    check("T11 exact DE certification 6/6 zero-remainder",
          len(cert) == 6 and all(ok for (ok, nb) in cert.values()))
    w1 = mp.sqrt(mp.mpf(4) / 3)
    refs = {(0, 0): mp.mpc(3) / 2,
            (0, 1): mp.mpc(3) / 2 * mp.log(mp.mpf(3) / 2),
            (1, 0): w1 * (1 + 2 * mp.log(2)),
            (2, 0): 2 * w1,
            (2, 1): mp.mpc(4) / 3}
    refs[(1, 1)] = w1 * mp.quad(
        lambda s: mp.sqrt(3 / (s + 3)) * ((s + 3) ** 2 / 3 - (s + 3)) / (s + 1), [0, 1])
    g = vc.gate(refs)
    worst = min(g.values())
    check("T12 held-out endpoint gate 6/6 >=40d (worst %.1fd)" % worst,
          len(g) == 6 and worst >= 40)
    check("T13 boundary PSLQ closed all 3 boundaries rational",
          all(str(v).startswith("exact") for v in vc.boundary_record.values())
          and len(vc.boundary_record) == 3)
    # 14: exact_entry_from_samples via tools/ratfit (blind + known-Q routes)
    tgt = Rr([1, 3], [2, 0, 1])                  # (1+3t)/(2+t^2)
    xs = [Fraction(k, 37) for k in range(1, 25)]
    ys = [Fraction(int(tgt.n(fq(x)).p), int(tgt.n(fq(x)).q))
          / Fraction(int(tgt.d(fq(x)).p), int(tgt.d(fq(x)).q)) for x in xs]
    Rb, okb, _ = exact_entry_from_samples(xs, ys)
    Rk, okk, _ = exact_entry_from_samples(xs, ys, Qknown=fp([2, 0, 1]))
    check("T14 ratfit delegation (blind + known-Q) exact", okb and okk and Rb == tgt and Rk == tgt)

    npass = sum(1 for _, ok in results if ok)
    print(f"vopclose selftest: {npass}/{len(results)} PASS ({time.time()-t0:.1f}s)", flush=True)
    return npass, len(results)

if __name__ == "__main__":
    if "--test" in sys.argv:
        n, N = _selftest()
        sys.exit(0 if n == N else 1)
    print(__doc__)

#!/usr/bin/env python3
# galois member — motivic Galois derivation engine (D3/D5/D7 coaction cuts).
"""core.py — motivic Galois derivations D3/D5/D7 on this package's
single-valued HPL spaces (the coaction-cut engine).

Mathematical content
--------------------
Work in the Goncharov iterated-integral convention
    I(a0; a1..an; a_{n+1}) = int_{a0<t1<..<tn<a_{n+1}} prod dt_i/(t_i - a_i)
with letters in {0,1}, shuffle regularization (zeta_sh(1) -> 0, log z -> 0 at
the 0-basepoint), and the motivic coproduct
    Delta I = sum_{subseq S} I(a0; A_S; a_{n+1}) (x) prod_p I(arc_p).
The weight-(2k+1) Galois derivation is D_m = (id (x) pi_m) Delta, m=3,5,7,
where pi_m projects the de Rham right factor onto the coefficient of the
single letter f_m (equivalently: the coefficient of the pure generator
monomial z_m in the double-shuffle-reduced zpoly — products and even zetas
project to 0).  D_m is a derivation (pi_m is an infinitesimal character).

Conventions (locked, package-wide):
  - H-words (Remiddi-Vermaseren, mzv_formal/svmpl): H_u(z), leftmost letter
    outermost;  H_u(z) = (-1)^{#1s(u)} I(0; rev(u); z).
  - zeta_sh(w) = shuffle-regularized H_w(1)  (t0deep mzv_tables_w8.json,
    zpoly over generators (p=i*pi, z3, z5, z7, z53), IKZ dual-gated).
  - Single-valued: calL_w(z,zb) = sum_{(u,s,c) in T[w]} c * H_u(z) H_s(zb)
    (t0deep svt_tables_w8.json, p-free), the I[z,w,0] = (-1)^{#1s} calL_w.
  - Ansatz zeta tokens: f3=z3, f5=z5, f7=z7, f33=z3^2/2, f3*f3=z3^2,
    f53=zeta(5,3) (=z53, word 00001001), z3*z5; f35=z3*z5-z53 (stuffle).

Arc values used by the coproduct (x,y in {0,1}, B a word):
  I(x; B; x) = 0 for |B|>=1, 1 for B empty.
  I(0; B; 1) = (-1)^{#1s(B)} zeta_sh(rev(B)).
  I(1; B; 0) = (-1)^{|B|} I(0; rev(B); 1) = (-1)^{|B|+#1s(B)} zeta_sh(B).
Final arc of a FUNCTION word I(0; A; z): only its z-constant part survives
pi_m (generic z kills every z-dependent de Rham factor):
  constpart I(x; B; z) = I(x; B; 0)  (0 if x=0; the reversal value if x=1).

Everything exact (Fraction / zpoly); no numerics except the G4 oracle gate.
"""
import json
import os
import sys
from fractions import Fraction
from functools import lru_cache
from itertools import combinations

T0DEEP = os.environ.get("GALOIS_T0DEEP")
if not T0DEEP:
    raise RuntimeError(
        "GALOIS_T0DEEP is not set. galois.core needs the t0deep MZV/sv table "
        "bank (mzv_formal.py + mzv_tables_w8.json + svt_tables_w8.json) — a "
        "reference table bank that does not ship with this repo. "
        "Point GALOIS_T0DEEP at a directory holding those three files.")
sys.path.insert(0, T0DEEP)

from mzv_formal import (all_words, zp_add, zp_mul, zp_scale, zp_const,
                        zp_weight_check)

W = 8

# ---------------- table load (SHA-verified by the driver) -------------------

def _parse_zpoly(d):
    return {tuple(int(x) for x in k.split(",")): Fraction(v)
            for k, v in d.items()}


def load_tables():
    """returns (ZSH, T): ZSH word->zpoly (all 510 words); T word->[(u,s,c)]"""
    mz = json.load(open(f"{T0DEEP}/mzv_tables_w8.json"))
    ZSH = {tuple(int(a) for a in w): _parse_zpoly(d)
           for w, d in mz["zeta_sh"].items()}
    sv = json.load(open(f"{T0DEEP}/svt_tables_w8.json"))
    T = {}
    for w, terms in sv["T"].items():
        key = tuple(int(a) for a in w)
        T[key] = [(tuple(int(a) for a in u), tuple(int(a) for a in s),
                   _parse_zpoly(c)) for (u, s, c) in terms]
    return ZSH, T


ZSH, T = load_tables()

# ---------------- zpoly helpers --------------------------------------------

GEN_KEY = {3: (0, 1, 0, 0, 0), 5: (0, 0, 1, 0, 0), 7: (0, 0, 0, 1, 0)}
Z53_KEY = (0, 0, 0, 0, 1)


def zp_weight(k):
    ep, e3, e5, e7, e53 = k
    return ep + 3 * e3 + 5 * e5 + 7 * e7 + 8 * e53


def pi_m(zp, m):
    """infinitesimal character: coefficient of the pure generator z_m."""
    return zp.get(GEN_KEY[m], Fraction(0))


# ---------------- I-convention word values ----------------------------------

def nones(B):
    return sum(1 for a in B if a == 1)


def I01(B):
    """I(0; B; 1) as zpoly (shuffle-regularized)."""
    if not B:
        return zp_const(1)
    z = ZSH[tuple(reversed(B))]
    return zp_scale(z, Fraction((-1) ** nones(B)))


def arc_const(x, B, y):
    """I(x; B; y) as zpoly, x,y in {0,1} (constant arc)."""
    if x == y:
        return zp_const(1) if not B else {}
    if not B:
        return zp_const(1)
    if x == 0:                       # (0 -> 1)
        return I01(B)
    # (1 -> 0): path reversal
    return zp_scale(I01(tuple(reversed(B))), Fraction((-1) ** len(B)))


def final_arc_constpart(x, B):
    """z-constant part of I(x; B; z) (shuffle-reg at the 0 basepoint):
    = I(x; B; 0)."""
    if not B:
        return zp_const(1)
    if x == 0:
        return {}
    return zp_scale(I01(tuple(reversed(B))), Fraction((-1) ** len(B)))


# ---------------- Goncharov D_m on I-words ----------------------------------

def Dm_iword(A, m, endpoint):
    """D_m I(0; A; endpoint), endpoint in {'one','z'}.
    Returns dict {kept_subword: Fraction} (endpoint='z': kept I(0;.;z) words)
    or, for endpoint='one', dict {kept_subword: Fraction} with kept words to
    be reduced via I01 by the caller."""
    n = len(A)
    if n < m:
        return {}
    out = {}
    for S in combinations(range(n), n - m):
        # product of arc constants
        prod = zp_const(1)
        prev = None          # index of previous kept letter
        left_val = 0         # a0
        ok = True
        for idx in S:
            gap = A[(prev + 1 if prev is not None else 0):idx]
            xval = left_val if prev is None else A[prev]
            arc = arc_const(xval, gap, A[idx])
            if not arc:
                ok = False
                break
            prod = zp_mul(prod, arc)
            prev = idx
        if not ok:
            continue
        tail = A[(S[-1] + 1 if S else 0):]
        xval = A[S[-1]] if S else 0
        if endpoint == 'one':
            arc = arc_const(xval, tail, 1)
        else:
            arc = final_arc_constpart(xval, tail)
        if not arc:
            continue
        prod = zp_mul(prod, arc)
        c = pi_m(prod, m)
        if c == 0:
            continue
        kept = tuple(A[i] for i in S)
        out[kept] = out.get(kept, Fraction(0)) + c
    return {k: v for k, v in out.items() if v != 0}


def Dm_word_value(w, m):
    """D_m zeta_sh(w) (H-convention word w) as zpoly, via the Goncharov
    coproduct on I(0; rev(w); 1)."""
    sign = Fraction((-1) ** nones(w))
    A = tuple(reversed(w))
    acc = {}
    for kept, c in Dm_iword(A, m, 'one').items():
        acc = zp_add(acc, I01(kept), sign * c)
    return acc


def Dm_hword(u, m):
    """D_m H_u(z) = sum_{u'} q * H_{u'}(z): dict {u': Fraction}."""
    if len(u) < m:
        return {}
    sign = Fraction((-1) ** nones(u))
    A = tuple(reversed(u))
    out = {}
    for kept, c in Dm_iword(A, m, 'z').items():
        up = tuple(reversed(kept))
        q = sign * c * Fraction((-1) ** nones(up))
        out[up] = out.get(up, Fraction(0)) + q
    return {k: v for k, v in out.items() if v != 0}


# ---------------- derivation on zpolys (generator images) -------------------

def fit_z53_images():
    """derive D_m(z53) directly: Goncharov on the zeta(5,3) word 00001001."""
    images = {}
    for m in (3, 5, 7):
        images[m] = Dm_word_value((0, 0, 0, 0, 1, 0, 0, 1), m)
    return images


Z53_IMG = fit_z53_images()


def Dm_zpoly(zp, m):
    """derivation on Q[p,z3,z5,z7,z53]: D_m p = 0, D_m z_j = delta_{mj},
    D_m z53 = Z53_IMG[m]; Leibniz on monomials."""
    out = {}
    gk = GEN_KEY[m]
    gi = gk.index(1)
    for k, q in zp.items():
        ep, e3, e5, e7, e53 = k
        exps = [ep, e3, e5, e7, e53]
        # generator z_m factor
        if exps[gi] > 0:
            k2 = list(exps)
            mult = k2[gi]
            k2[gi] -= 1
            key = tuple(k2)
            out[key] = out.get(key, Fraction(0)) + q * mult
        # z53 factor
        if e53 > 0 and Z53_IMG[m]:
            k2 = list(exps)
            mult = k2[4]
            k2[4] -= 1
            base = {tuple(k2): q * mult}
            for kk, vv in zp_mul(base, Z53_IMG[m]).items():
                out[kk] = out.get(kk, Fraction(0)) + vv
    return {k: v for k, v in out.items() if v != 0}


# ---------------- single-valued space: coordinates + decomposition ----------
# weight-w sv space basis: pairs (word, mono) with |word| + wt(mono) = w,
# mono a p-free zpoly monomial key in the ODD-zeta subring (e53 included).

def sv_monos(wt):
    """p-free zeta monomial keys of weight wt (e3,e5,e7,e53 only)."""
    out = []
    for e53 in range(0, wt // 8 + 1):
        for e7 in range(0, (wt - 8 * e53) // 7 + 1):
            for e5 in range(0, (wt - 8 * e53 - 7 * e7) // 5 + 1):
                r = wt - 8 * e53 - 7 * e7 - 5 * e5
                if r % 3 == 0:
                    out.append((0, r // 3, e5, e7, e53))
    return out


def sv_basis(wt):
    """ordered list of (word, mono) for the weight-wt sv space."""
    basis = []
    for lw in range(wt, -1, -1):
        for mono in sv_monos(wt - lw):
            for w in all_words(lw):
                if len(w) == lw:
                    basis.append((w, mono))
    return basis


class SVDecomposer:
    """exact decomposition of bilinear objects {(u,s): zpoly} into the
    weight-wt sv basis  sum_j x_j * mono_j * calL_{w_j}."""

    def __init__(self, wt):
        self.wt = wt
        self.basis = sv_basis(wt)
        # coordinate index: (u, s, monokey) -> row
        coords = {}
        cols = []
        for (w, mono) in self.basis:
            col = {}
            for (u, s, c) in T[w]:
                cm = zp_mul(c, {mono: Fraction(1)})
                for k, q in cm.items():
                    key = (u, s, k)
                    col[key] = col.get(key, Fraction(0)) + q
            col = {k: v for k, v in col.items() if v != 0}
            for k in col:
                coords.setdefault(k, len(coords))
            cols.append(col)
        self.coords = coords
        # RREF of the basis matrix (columns = basis elements)
        # rows: coordinate index -> {basis col: coeff}
        rows = {}
        for j, col in enumerate(cols):
            for k, v in col.items():
                rows.setdefault(self.coords[k], {})[j] = v
        # Gaussian elimination over coordinates to build a solver:
        # we store the matrix column-wise and do elimination on demand via
        # exact LU on the (ncoord x nbasis) system each solve -- but nbasis
        # is small (<=37), so precompute a pivot decomposition once.
        self.cols = cols
        self.n = len(cols)
        # build pivot structure: forward elimination on columns
        piv = []                     # list of (coordkey_index, colvec, j)
        reduced_cols = []
        for j, col in enumerate(cols):
            r = dict(col)
            xrep = {j: Fraction(1)}  # representation in original basis cols
            for (pk, pcol, prep) in piv:
                if pk in r and r[pk] != 0:
                    f = r[pk] / pcol[pk]
                    for k2, v2 in pcol.items():
                        nv = r.get(k2, Fraction(0)) - f * v2
                        if nv:
                            r[k2] = nv
                        elif k2 in r:
                            del r[k2]
                    for k2, v2 in prep.items():
                        nv = xrep.get(k2, Fraction(0)) - f * v2
                        if nv:
                            xrep[k2] = nv
                        elif k2 in xrep:
                            del xrep[k2]
            r = {self.coords[k] if not isinstance(k, int) else k: v
                 for k, v in r.items() if v != 0} if False else \
                {k: v for k, v in r.items() if v != 0}
            if not r:
                raise AssertionError(
                    f"sv basis at weight {self.wt} linearly dependent (col {j})")
            pk = min(r, key=lambda k: self.coords[k] if not isinstance(k, int)
                     else k) if False else min(r, key=lambda k: str(k))
            piv.append((pk, r, xrep))
            reduced_cols.append(r)
        self.piv = piv

    def target_coords(self, bil):
        """bilinear {(u,s): zpoly} -> {(u,s,monokey): Fraction}"""
        t = {}
        for (u, s), c in bil.items():
            for k, q in c.items():
                key = (u, s, k)
                t[key] = t.get(key, Fraction(0)) + q
        return {k: v for k, v in t.items() if v != 0}

    def decompose(self, bil):
        """solve sum x_j basis_j == bil exactly; returns vector list
        [(basis_index, Fraction)]; raises if not in span."""
        t = self.target_coords(bil)
        # any coordinate not seen in the basis span => not in span
        for k in t:
            if k not in self.coords:
                raise AssertionError(
                    f"sv-closure FAIL at weight {self.wt}: coordinate {k} "
                    f"outside basis span")
        r = dict(t)
        sol = {}
        for (pk, pcol, prep) in self.piv:
            if pk in r and r[pk] != 0:
                f = r[pk] / pcol[pk]
                for k2, v2 in pcol.items():
                    nv = r.get(k2, Fraction(0)) - f * v2
                    if nv:
                        r[k2] = nv
                    elif k2 in r:
                        del r[k2]
                for k2, v2 in prep.items():
                    sol[k2] = sol.get(k2, Fraction(0)) + f * v2
        if r:
            raise AssertionError(
                f"sv-closure FAIL at weight {self.wt}: nonzero residual "
                f"({len(r)} coords, sample {list(r.items())[:3]})")
        return [(j, q) for j, q in sorted(sol.items()) if q != 0]


_DECOMP = {}


def decomposer(wt):
    if wt not in _DECOMP:
        _DECOMP[wt] = SVDecomposer(wt)
    return _DECOMP[wt]


# ---------------- D_m on calL words and on sv vectors -----------------------

@lru_cache(maxsize=None)
def Dm_calL(w, m):
    """D_m calL_w as an exact vector over the weight-(|w|-m) sv basis.
    Returns tuple of (basis_index, Fraction)."""
    wt2 = len(w) - m
    if wt2 < 0:
        return ()
    bil = {}
    for (u, s, c) in T[w]:
        # constants
        dc = Dm_zpoly(c, m)
        if dc:
            bil[(u, s)] = zp_add(bil.get((u, s), {}), dc)
        # holomorphic function factor
        for u2, q in Dm_hword(u, m).items():
            bil[(u2, s)] = zp_add(bil.get((u2, s), {}), zp_scale(c, q))
        # antiholomorphic function factor
        for s2, q in Dm_hword(s, m).items():
            bil[(u, s2)] = zp_add(bil.get((u, s2), {}), zp_scale(c, q))
    bil = {k: v for k, v in bil.items() if v}
    if not bil:
        return ()
    dec = decomposer(wt2)
    return tuple(dec.decompose(bil))


def Dm_vector(vec, wt, m):
    """D_m on an sv vector {(word,mono): Fraction} of pure weight wt.
    Returns {(word,mono): Fraction} at weight wt-m."""
    wt2 = wt - m
    out = {}
    if wt2 < 0:
        return out
    dec2 = decomposer(wt2) if wt2 >= 0 else None
    for (w, mono), q in vec.items():
        # term q * mono * calL_w
        # (a) derivation on the constant monomial
        dmono = Dm_zpoly({mono: Fraction(1)}, m)
        for k2, v2 in dmono.items():
            assert zp_weight(k2) == zp_weight(mono) - m
            key = (w, k2)
            out[key] = out.get(key, Fraction(0)) + q * v2
        # (b) derivation on the function
        if len(w) >= m:
            for (j, c) in Dm_calL(w, m):
                (w2, mono2) = decomposer(len(w) - m).basis[j]
                mm = zp_mul({mono2: Fraction(1)}, {mono: Fraction(1)})
                assert len(mm) == 1
                k3 = next(iter(mm))
                key = (w2, k3)
                out[key] = out.get(key, Fraction(0)) + q * c * mm[k3]
    return {k: v for k, v in out.items() if v != 0}


# ---------------- conjugation (z <-> zb) ------------------------------------

@lru_cache(maxsize=None)
def conj_calL(w):
    """conj calL_w as vector over the weight-|w| sv basis (transpose of the
    (u,s) bilinear expansion)."""
    bil = {}
    for (u, s, c) in T[w]:
        bil[(s, u)] = zp_add(bil.get((s, u), {}), c)
    dec = decomposer(len(w))
    return tuple(dec.decompose(bil))


def conj_vector(vec):
    """conjugation on an sv vector {(word,mono): Fraction} (pure weight)."""
    out = {}
    for (w, mono), q in vec.items():
        for (j, c) in conj_calL(w):
            (w2, mono2) = decomposer(len(w)).basis[j]
            mm = zp_mul({mono2: Fraction(1)}, {mono: Fraction(1)})
            k3 = next(iter(mm))
            key = (w2, k3)
            out[key] = out.get(key, Fraction(0)) + q * c
    return {k: v for k, v in out.items() if v != 0}


# ---------------- ansatz-element conversion ---------------------------------

ZTOK = {
    "1": {(0, 0, 0, 0, 0): Fraction(1)},
    "f3": {(0, 1, 0, 0, 0): Fraction(1)},
    "f5": {(0, 0, 1, 0, 0): Fraction(1)},
    "f7": {(0, 0, 0, 1, 0): Fraction(1)},
    "f33": {(0, 2, 0, 0, 0): Fraction(1, 2)},
    "f3*f3": {(0, 2, 0, 0, 0): Fraction(1)},
    "f53": {(0, 0, 0, 0, 1): Fraction(1)},
    "f35": {(0, 1, 1, 0, 0): Fraction(1), (0, 0, 0, 0, 1): Fraction(-1)},
    "z3*z5": {(0, 1, 1, 0, 0): Fraction(1)},
    "f3*f5": {(0, 1, 1, 0, 0): Fraction(1)},
}


def element_vector(el):
    """parsed ansatz element (list of {coeff, word, zeta}) -> sv vector
    {(word,mono): Fraction} with the the I[z,w,0] -> calL sign."""
    vec = {}
    for t in el:
        q = Fraction(t["coeff"])
        w = tuple(t["word"])
        sign = Fraction((-1) ** nones(w))
        zt = ZTOK[t["zeta"]]
        for k, f in zt.items():
            key = (w, k)
            vec[key] = vec.get(key, Fraction(0)) + q * sign * f
    return {k: v for k, v in vec.items() if v != 0}

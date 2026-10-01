#!/usr/bin/env python3
"""
heldout_cv.py — leave-one-out cross-validation + PSLQ rational certifier for
oracle per-master-per-point dumps.

Honest value-fit gate: a master m at ε-order w is CERTIFIED only if, for
*every* held-out point P_k, the rational fit on the remaining N-1 points
predicts the value at P_k to ≥ `gate` digits.

Honesty fields reported per (master, weight):
    n_points_used, n_duplicates_excluded, basis_dimension, condition_number,
    heldout_digits[], min, mean, certified, pslq_rationals (or null).

Dependencies: stdlib + mpmath only; fpylll optional (fast LLL).

Conditioning layer:
    A recurring wall on boundary fits is that the *basis-evaluation
    matrix* can reach cond ~1e36, so LSQ+PSLQ floors on conditioning, not data.
    This module now ACTS on the condition number:

        lyndon_basis()    canonical well-conditioned shuffle-algebra basis
        auto_condition()  LLL-reduce columns of the basis matrix → unimodular U
        pslq_with_lll()   LLL-preconditioned integer-relation finder
                          (robust for >20 basis elements where raw PSLQ stalls)

    heldout_certify() applies an optional UT rotation (raw masters → canonical
    leading-singularity basis), then auto-conditions, then fits, then
    rationalises with pslq_with_lll.  New honesty fields:
        cond_before, cond_after, lll_applied, pslq_method.

CLI:
    python heldout_cv.py --oracle DIR [DIR ...] \
        --basis basis_values.json --weight 2 \
        [--gate 30] [--integrity integrity_report.json] \
        [--report cv_report.json] [--dps 80]

basis_values.json schema (one of):
  {"points": {"<s>|<t>": {"B0": "<mpf-string>", "B1": "...", ...}, ...},
   "names": ["B0","B1",...]}
or per-weight:
  {"weights": {"2": {"points": {...}, "names": [...]}, ...}}

A basis evaluator may also be supplied programmatically:
    heldout_certify(values, basis_fn=lambda s,t: [mpf,...], names=[...], ...)
"""
from __future__ import annotations
import argparse
import glob
import json
import math
import os
import sys
from collections import defaultdict
from dataclasses import dataclass, field, asdict

import mpmath as mp

try:                                  # optional fast LLL
    from fpylll import IntegerMatrix, LLL as _fpylll_LLL
    _HAVE_FPYLLL = True
except (ImportError, ValueError):     # pure-Python fallback below
    # ValueError: a binary-incompatible fpylll (numpy ABI drift) raises
    # ValueError, not ImportError, at import — an uncaught raise here
    # would kill the whole CV/PSLQ certifier + its consumers; nothing
    # recorded cites the dead path
    _HAVE_FPYLLL = False

try:                                  # allow both `python heldout_cv.py` and
    from .oracle_integrity import scan as _scan_integrity, IDENTITY_KEYS
except ImportError:                   # `from gatekeeper import ...`
    from oracle_integrity import scan as _scan_integrity, IDENTITY_KEYS


# ====================================================================== I/O
def _tag_of(d):
    for k in ("tag", "k", "master", "name"):
        if k in d:
            return str(d[k])
    return "?"


def _ptkey(d):
    return f"{d.get('s')}|{d.get('t')}"


def _to_mpf(x):
    if isinstance(x, (int, float)):
        return mp.mpf(x)
    s = str(x)
    if "/" in s:
        a, b = s.split("/")
        return mp.mpf(a) / mp.mpf(b)
    return mp.mpf(s)


def load_oracle_values(dirs, exclude_files=None):
    """Walk DIRS for oracle JSONs, return
        values[tag][ptkey][order] = mpc
        kinematics[ptkey] = (s_mpf, t_mpf)
    Files in `exclude_files` (a set of abs paths) are skipped."""
    if isinstance(dirs, (str, os.PathLike)):
        dirs = [dirs]
    exclude_files = set(exclude_files or ())
    values = defaultdict(lambda: defaultdict(dict))
    kinematics = {}
    n_excluded = 0
    for d in dirs:
        pattern = os.path.join(d, "**", "*.json")
        for f in sorted(glob.glob(pattern, recursive=True)):
            if os.path.abspath(f) in exclude_files or f in exclude_files:
                n_excluded += 1
                continue
            try:
                rec = json.load(open(f))
            except Exception:
                continue
            tag = _tag_of(rec)
            pk = _ptkey(rec)
            try:
                kinematics[pk] = (_to_mpf(rec.get("s", 0)),
                                  _to_mpf(rec.get("t", 0)))
            except Exception:
                pass
            coeffs = rec.get("coeffs") or rec.get("coefficients") or {}
            for o, v in coeffs.items():
                try:
                    re = mp.mpf(str(v.get("re", "0")))
                    im = mp.mpf(str(v.get("im", "0")))
                except Exception:
                    continue
                values[tag][pk][int(o)] = mp.mpc(re, im)
    return dict(values), kinematics, n_excluded


# =============================================================== basis input
def load_basis_json(path, weight=None):
    """Return (basis_fn, names) from a basis_values.json. basis_fn maps a
    point-key string -> list[mpf]."""
    spec = json.load(open(path))
    if "weights" in spec:
        if weight is None:
            raise ValueError("basis file is per-weight; pass --weight")
        spec = spec["weights"][str(weight)]
    names = spec.get("names")
    pts = spec["points"]
    if names is None:
        # infer from first point
        names = sorted(next(iter(pts.values())).keys())
    table = {pk: [mp.mpf(str(pts[pk][n])) for n in names] for pk in pts}

    def basis_fn(ptkey, s=None, t=None):
        return table[ptkey]

    return basis_fn, names


# ============================================================ Lyndon words
def _duval_lyndon_words(alphabet, n):
    """Duval's algorithm: yield all Lyndon words of length EXACTLY n over
    `alphabet` (a sorted sequence of letters), in lexicographic order."""
    k = len(alphabet)
    idx = {a: i for i, a in enumerate(alphabet)}
    w = [alphabet[0]]
    while w:
        if len(w) == n:
            yield tuple(w)
        # Duval step: repeat w to length n, then increment last letter
        m = len(w)
        w = [w[i % m] for i in range(n)]
        while w and idx[w[-1]] == k - 1:
            w.pop()
        if w:
            w[-1] = alphabet[idx[w[-1]] + 1]


def _shuffle(u, v):
    """Shuffle product of two words (tuples) → dict {word: integer coeff}."""
    if not u:
        return {v: 1}
    if not v:
        return {u: 1}
    out = defaultdict(int)
    for w, c in _shuffle(u[1:], v).items():
        out[(u[0],) + w] += c
    for w, c in _shuffle(u, v[1:]).items():
        out[(v[0],) + w] += c
    return out


def _shuffle_many(words):
    """Iterated shuffle of a list of words → dict {word: int coeff}."""
    acc = {(): 1}
    for w in words:
        nxt = defaultdict(int)
        for a, ca in acc.items():
            for b, cb in _shuffle(a, w).items():
                nxt[b] += ca * cb
        acc = nxt
    return acc


def lyndon_basis(alphabet, weight):
    """Canonical Lyndon-word basis for the weight-`weight` graded piece of the
    shuffle algebra over `alphabet`.

    By Radford's theorem the shuffle algebra is the FREE polynomial algebra on
    Lyndon words, so the set
        { L_{i1} ⧢ L_{i2} ⧢ ... ⧢ L_{ik}  :  L_{i1} >= ... >= L_{ik} Lyndon,
                                              |L_{i1}|+...+|L_{ik}| = weight }
    is a basis of the same dimension (|alphabet|^weight) as the naive
    "all words of length `weight`" basis, but is the well-conditioned one for
    iterated-integral value fits (it kills the shuffle redundancies that make
    the naive basis nearly-degenerate numerically).

    Returns
    -------
    lyn_names : list[str]
        Human-readable names "L1⧢L2⧢..." for each Lyndon basis element.
    all_words : list[tuple]
        All words of length `weight` (the naive basis), lexicographically.
    T : list[list[Fraction]]
        Transition matrix, shape (n_lyn, n_words), with RATIONAL entries:
            lyn_basis[i]  =  Σ_j  T[i][j] · all_words[j]
        (so basis_values_lyn = T · basis_values_words).
    lyndon_words_by_len : dict[int, list[tuple]]
        The raw Lyndon words at each length 1..weight (Duval output).
    """
    from fractions import Fraction
    from itertools import product, combinations_with_replacement
    alphabet = sorted(alphabet)
    # Lyndon words of each length 1..weight
    lyn = {L: list(_duval_lyndon_words(alphabet, L)) for L in range(1, weight + 1)}
    # all multisets of Lyndon words with total length == weight
    pool = [(L, w) for L in range(1, weight + 1) for w in lyn[L]]
    factored = []
    for k in range(1, weight + 1):
        for combo in combinations_with_replacement(range(len(pool)), k):
            if sum(pool[i][0] for i in combo) == weight:
                # store in non-increasing order (Radford convention)
                factored.append(tuple(sorted((pool[i][1] for i in combo),
                                             reverse=True)))
    factored = sorted(set(factored), reverse=True)
    # naive basis
    all_words = ["".join(w) for w in product(alphabet, repeat=weight)]
    word_index = {w: j for j, w in enumerate(all_words)}
    # transition matrix (rational; entries are in fact integers here)
    T = [[Fraction(0)] * len(all_words) for _ in factored]
    lyn_names = []
    for i, fac in enumerate(factored):
        lyn_names.append("⧢".join("".join(w) for w in fac))
        for word, c in _shuffle_many(list(fac)).items():
            T[i][word_index["".join(word)]] += Fraction(c)
    return lyn_names, all_words, T, lyn


# ====================================================================== LLL
def _int_lll(B, delta=0.99):
    """Pure-Python integer LLL on the ROWS of B (list of list[int]).
    Lovász condition parameter δ ∈ (1/4, 1).  Returns the reduced basis
    in place (B is mutated) and also returns B for convenience.

    ~50 lines; O(n^2 d) Gram–Schmidt per swap.  Fine for n ≲ 100; for
    bigger lattices install fpylll."""
    n = len(B)
    if n == 0:
        return B
    m = len(B[0])

    # working precision must cover the integer magnitudes (entries can be
    # 10^prec-scaled values, i.e. hundreds of digits)
    max_bits = 1
    for row in B:
        for x in row:
            bl = int(x).bit_length() if x else 1
            if bl > max_bits:
                max_bits = bl
    _wdps = max(mp.mp.dps, int(max_bits * 0.3011) + 30)

    # Gram–Schmidt (rational, via mp.mpf at high prec for stability)
    def gram():
        Bs = [[mp.mpf(x) for x in row] for row in B]
        mu = [[mp.mpf(0)] * n for _ in range(n)]
        Bstar = [list(Bs[0])]
        beta = [sum(x * x for x in Bstar[0])]
        for i in range(1, n):
            v = list(Bs[i])
            for j in range(i):
                bj = beta[j]
                mu[i][j] = (sum(Bs[i][k] * Bstar[j][k] for k in range(m)) / bj
                            if bj != 0 else mp.mpf(0))
                for k in range(m):
                    v[k] -= mu[i][j] * Bstar[j][k]
            Bstar.append(v)
            beta.append(sum(x * x for x in v))
        return mu, beta

    with mp.workdps(_wdps):
        mu, beta = gram()
        k = 1
        while k < n:
            # size-reduce
            for j in range(k - 1, -1, -1):
                q = int(mp.nint(mu[k][j]))
                if q:
                    for t in range(m):
                        B[k][t] -= q * B[j][t]
                    mu[k][j] -= q
                    for t in range(j):
                        mu[k][t] -= q * mu[j][t]
            # Lovász condition
            if beta[k] >= (mp.mpf(delta) - mu[k][k - 1] ** 2) * beta[k - 1]:
                k += 1
            else:
                B[k], B[k - 1] = B[k - 1], B[k]
                mu, beta = gram()           # cheap & correct; recompute
                k = max(1, k - 1)
    return B


def _lll_rows(M):
    """Dispatch: fpylll if available, else pure-Python.  M is list[list[int]].
    Returns reduced rows as list[list[int]]."""
    if _HAVE_FPYLLL:
        A = IntegerMatrix(len(M), len(M[0]))
        for i, row in enumerate(M):
            for j, x in enumerate(row):
                A[i, j] = int(x)
        _fpylll_LLL.reduction(A)
        return [[int(A[i, j]) for j in range(A.ncols)] for i in range(A.nrows)]
    return _int_lll([list(r) for r in M])


# ============================================================= conditioning
def _cond_svd(A):
    """SVD-based condition number of an mp.matrix (real)."""
    try:
        _, S, _ = mp.svd_r(A, full_matrices=False)
        sv = [abs(S[i]) for i in range(min(A.rows, A.cols))]
    except Exception:
        # fall back: sqrt of eig(AtA) extremes
        AtA = A.T * A
        ev = mp.mp.eig(AtA, left=False, right=False)
        sv = [mp.sqrt(abs(e)) for e in ev]
    sv = [s for s in sv if s > 0]
    return (max(sv) / min(sv)) if sv else mp.inf


def auto_condition(basis_matrix, target_cond=1e6, prec=None, _log=True,
                   drop_tol_digits=None):
    """If the basis-evaluation matrix is ill-conditioned, LLL-reduce its
    COLUMN lattice to find an integer unimodular transform U with
    cond(basis_matrix · U) ≪ cond(basis_matrix).

    The construction is the standard one: embed column j of A (scaled to a
    `prec`-digit integer) together with the j-th unit vector as row j of an
    n×(m+n) lattice, LLL-reduce the rows, and read U off the trailing n×n
    block.  Reduced columns whose norm is ≤ 10^{-drop_tol_digits} · max_norm
    are DROPPED — these are the (near-)exact integer dependencies among the
    original basis functions (shuffle relations, in the polylog application),
    so the returned U is n×r with r = effective rank.

    Parameters
    ----------
    basis_matrix : mp.matrix, shape (N_points, n_basis)
    target_cond  : float — skip LLL if already ≤ this.
    prec         : int   — digits to scale by before rounding to integers
                           (default: 2/3 of mp.mp.dps, capped).
    drop_tol_digits : near-zero reduced columns are dropped if their norm is
                      this many digits below the largest (default: prec//2).

    Returns
    -------
    U            : mp.matrix (n_basis × r), integer entries, with
                   det(full U)=±1 before any column drop; or None if no
                   reduction was applied.
    cond_before  : mp.mpf
    cond_after   : mp.mpf  (of AU on the r retained columns)
    AU           : mp.matrix (N_points × r) — the reduced basis matrix,
                   read DIRECTLY from the LLL leading block (so accurate
                   to ~prec digits; do NOT recompute as basis_matrix*U at
                   the caller's dps, which can cancel catastrophically when
                   |U| is large).
    """
    A = basis_matrix
    mrows, n = A.rows, A.cols
    cond0 = _cond_svd(A)
    if cond0 <= target_cond:
        if _log:
            print(f"[cv] basis cond {mp.nstr(cond0, 3)} ≤ target "
                  f"{target_cond:.0e}; no LLL.")
        return None, cond0, cond0, A

    if prec is None:
        prec = max(20, min((2 * mp.mp.dps) // 3, 300))
    if drop_tol_digits is None:
        drop_tol_digits = max(8, prec // 2)
    Nscale = mp.mpf(10) ** prec
    # Weight on the trailing identity block: large enough that LLL is
    # penalised for using transform coefficients larger than ~10^max_u_digits.
    # Without this, near-exact column dependencies create rounding-noise
    # lattice vectors that LLL multiplies by ~Nscale, polluting U.
    max_u_digits = 8
    Amax = max(abs(A[i, j]) for i in range(mrows) for j in range(n))
    W = max(1, int(Nscale * max(Amax, 1) / mp.mpf(10) ** max_u_digits))
    # n × (mrows + n) integer lattice:  row j = [ round(N·col_j(A)) | W·e_j ]
    M = []
    for j in range(n):
        row = [int(mp.nint(Nscale * A[i, j])) for i in range(mrows)]
        row += [W if jj == j else 0 for jj in range(n)]
        M.append(row)
    R = _lll_rows(M)
    # full unimodular transform (columns of U = trailing-n of each reduced row,
    # divided by the weight W — exact since every reduced row is an integer
    # combination of the originals, whose trailing blocks are all W·e_j)
    # and reduced basis AU read off the LEADING block (exact-int / Nscale)
    Ufull = mp.matrix(n, n)
    AU = mp.matrix(mrows, n)
    for k in range(n):
        for j in range(n):
            Ufull[j, k] = R[k][mrows + j] // W
        for i in range(mrows):
            AU[i, k] = mp.mpf(R[k][i]) / Nscale
    # drop near-zero reduced columns (= integer dependencies in the basis)
    col_norms = [mp.sqrt(sum(AU[i, k] ** 2 for i in range(mrows)))
                 for k in range(n)]
    nm_max = max(col_norms) if col_norms else mp.mpf(1)
    keep = [k for k in range(n)
            if col_norms[k] > nm_max * mp.mpf(10) ** (-drop_tol_digits)]
    r = len(keep)
    if r < n:
        U = mp.matrix(n, r)
        AUr = mp.matrix(mrows, r)
        for kk, k in enumerate(keep):
            for j in range(n):
                U[j, kk] = Ufull[j, k]
            for i in range(mrows):
                AUr[i, kk] = AU[i, k]
    else:
        U, AUr = Ufull, AU
    cond1 = _cond_svd(AUr)
    if _log:
        drop_msg = f", dropped {n - r} dependent col(s)" if r < n else ""
        print(f"[cv] basis cond {mp.nstr(cond0, 3)} → LLL → "
              f"{mp.nstr(cond1, 3)}{drop_msg}; transform stored.")
    return U, cond0, cond1, AUr


# ================================================== LLL-preconditioned PSLQ
def pslq_with_lll(values, prec=None, max_coeff_bits=64, _verify_tol_digits=None):
    """Find a small integer relation  Σ a_i · values[i] ≈ 0  by LLL on the
    standard relation lattice, falling back to mpmath.pslq.

    Much more robust than raw PSLQ once len(values) ≳ 20, because LLL is a
    *basis* reduction (simultaneous), whereas PSLQ is a sequential nearest-
    plane walk that stalls on ill-conditioned inputs.

    Parameters
    ----------
    values : sequence of mp.mpf / mp.mpc (real parts used).
    prec   : working digits for the integer scaling (default mp.mp.dps - 8).
    max_coeff_bits : reject relations whose max |a_i| exceeds 2**this.

    Returns
    -------
    (relation, method)  where relation is a list[int] or None, and
    method ∈ {"lll", "pslq", None}.
    """
    vals = [mp.mpf(v.real if isinstance(v, mp.mpc) else v) for v in values]
    n = len(vals)
    if prec is None:
        prec = max(10, mp.mp.dps - 8)
    if _verify_tol_digits is None:
        _verify_tol_digits = max(5, prec // 2)
    Nscale = mp.mpf(10) ** prec
    cap = 2 ** max_coeff_bits
    # lattice rows:  [ e_i  |  round(N · v_i) ]
    M = []
    for i in range(n):
        row = [0] * n
        row[i] = 1
        row.append(int(mp.nint(Nscale * vals[i])))
        M.append(row)
    R = _lll_rows(M)
    # candidate = shortest row by ℓ∞ on the coeff block with small residual
    best = None
    for row in R:
        coeffs, resid = row[:n], row[n]
        if all(c == 0 for c in coeffs):
            continue
        h = max(abs(c) for c in coeffs)
        if h > cap:
            continue
        # true residual at full precision
        r = abs(sum(c * v for c, v in zip(coeffs, vals)))
        scale = max(abs(v) for v in vals) * h
        ok = (r / scale if scale > 0 else r) < mp.mpf(10) ** (-_verify_tol_digits)
        if ok and (best is None or h < best[1]):
            best = (list(coeffs), h)
    if best is not None:
        return best[0], "lll"
    # fallback
    try:
        rel = mp.pslq(vals, tol=mp.mpf(10) ** (-_verify_tol_digits),
                      maxcoeff=cap)
    except Exception:
        rel = None
    if rel is not None:
        return list(rel), "pslq"
    return None, None


def _rational_from_relation(values, rel, target_index=0):
    """Given Σ a_i v_i = 0 with v[target_index] the unknown and the rest a
    basis, return the rational string for v[target] in that basis when the
    basis is just [1] (i.e. n=2, p/q recognition), else return the integer
    coeff vector with the target normalised out."""
    a0 = rel[target_index]
    if a0 == 0:
        return None
    # simple p/q case: values = [x, 1]  →  a0·x + a1 = 0  →  x = -a1/a0
    if len(values) == 2 and target_index == 0:
        p, q = -rel[1], rel[0]
        g = math.gcd(abs(p), abs(q))
        if g:
            p, q = p // g, q // g
        if q < 0:
            p, q = -p, -q
        return f"{p}/{q}" if q != 1 else f"{p}"
    return rel


# ============================================================ linear algebra
def _lstsq(A, b):
    """Least-squares solve A c = b via normal equations at current mp.dps.
    Returns (c, condition_number_estimate)."""
    m, n = A.rows, A.cols
    AtA = mp.matrix(n, n)
    Atb = mp.matrix(n, 1)
    for i in range(n):
        for j in range(n):
            s = mp.mpf(0)
            for k in range(m):
                s += A[k, i] * A[k, j]
            AtA[i, j] = s
        s = mp.mpf(0)
        for k in range(m):
            s += A[k, i] * b[k]
        Atb[i] = s
    c = mp.lu_solve(AtA, Atb)
    # crude condition estimate: ratio of diag extremes of AtA
    diags = [abs(AtA[i, i]) for i in range(n)]
    cond = (max(diags) / min(diags)) if min(diags) > 0 else mp.inf
    return c, cond


def _pslq_rational(x, tol_digits, max_coeff_bits=48):
    """Recognise x as p/q via LLL-preconditioned relation search on [x, 1],
    falling back to mpmath.pslq.  Returns ('p/q' or 'p', method) or
    (None, None)."""
    if x == 0:
        return "0", "exact"
    rel, method = pslq_with_lll([mp.mpf(x), mp.mpf(1)],
                                prec=max(tol_digits, 20),
                                max_coeff_bits=max_coeff_bits,
                                _verify_tol_digits=tol_digits)
    if rel is None:
        return None, None
    s = _rational_from_relation([x, 1], rel, target_index=0)
    if s is None:
        return None, None
    # verify
    if "/" in s:
        p, q = s.split("/")
        val = mp.mpf(p) / mp.mpf(q)
    else:
        val = mp.mpf(s)
    if abs(val - x) > mp.mpf(10) ** (-(tol_digits - 2)):
        return None, None
    return s, method


# =================================================================== result
@dataclass
class CVResult:
    master: str
    weight: int
    n_points_used: int
    n_duplicates_excluded: int
    basis_dimension: int
    basis_names: list
    condition_number: str
    cond_before: str = "n/a"
    cond_after: str = "n/a"
    lll_applied: bool = False
    pslq_method: str | None = None
    ut_rotation_applied: bool = False
    heldout_digits: list = field(default_factory=list)
    min: float = 0.0
    mean: float = 0.0
    certified: bool = False
    pslq_rationals: list | None = None
    note: str | None = None

    def to_json(self):
        return asdict(self)


# ============================================================ UT rotation
def _apply_ut_rotation(values, ut_rotation):
    """Rotate raw-master values to the canonical UT (leading-singularity)
    basis BEFORE any fitting.  `ut_rotation` is one of

        {ut_tag: {raw_tag: 'p/q' | number}, ...}
        (matrix, row_names, col_names)   — matrix is list[list] or mp.matrix

    Returns a NEW values dict {ut_tag: {ptkey: {order: mpc}}} containing only
    the rotated masters; raw tags not appearing as inputs are passed through
    untouched."""
    if ut_rotation is None:
        return values
    if isinstance(ut_rotation, dict):
        spec = ut_rotation
    else:
        Mx, row_names, col_names = ut_rotation
        if not isinstance(Mx, mp.matrix):
            Mx = mp.matrix(Mx)
        spec = {row_names[i]: {col_names[j]: Mx[i, j]
                               for j in range(Mx.cols) if Mx[i, j] != 0}
                for i in range(Mx.rows)}
    consumed = set()
    for combo in spec.values():
        consumed.update(combo.keys())
    new = {t: v for t, v in values.items() if t not in consumed}
    for ut_tag, combo in spec.items():
        coeffs = {rt: _to_mpf(c) if not isinstance(c, (mp.mpf, mp.mpc))
                  else c for rt, c in combo.items()}
        # all points where EVERY input raw master is available
        common_pks = None
        for rt in coeffs:
            pks = set(values.get(rt, {}).keys())
            common_pks = pks if common_pks is None else (common_pks & pks)
        if not common_pks:
            continue
        new[ut_tag] = {}
        for pk in common_pks:
            orders = None
            for rt in coeffs:
                o = set(values[rt][pk].keys())
                orders = o if orders is None else (orders & o)
            new[ut_tag][pk] = {}
            for o in orders:
                s = mp.mpc(0)
                for rt, c in coeffs.items():
                    s += c * values[rt][pk][o]
                new[ut_tag][pk][o] = s
    return new


# ==================================================================== core
def heldout_certify(values, kinematics, *, weight, basis_fn, basis_names,
                    gate=30, dps=80, n_duplicates_excluded=0,
                    masters=None, use_real=True,
                    ut_rotation=None, cond_target=1e6, auto_cond=True):
    """Leave-one-out CV over all masters in `values` at ε-order `weight`.

    values      : {tag: {ptkey: {order: mpc}}}
    kinematics  : {ptkey: (s_mpf, t_mpf)}  (only needed if basis_fn wants s,t)
    basis_fn    : callable(ptkey, s, t) -> list[mpf]  OR  callable(ptkey)
    basis_names : list[str]
    gate        : min held-out digits to CERTIFY
    use_real    : fit the real part (Euclidean default); set False to fit |.|
    ut_rotation : optional rational rotation raw→UT basis, applied FIRST
                  (see _apply_ut_rotation for accepted formats).
    cond_target : if SVD cond of basis matrix exceeds this, auto_condition()
                  LLL-reduces the column lattice and the fit is done in the
                  transformed basis (coefficients mapped back via U).
    auto_cond   : set False to disable conditioning (reproduce old behaviour).

    Returns {tag: CVResult}.
    """
    mp.mp.dps = dps
    if ut_rotation is not None:
        values = _apply_ut_rotation(values, ut_rotation)
    out = {}
    targets = masters if masters is not None else sorted(values.keys())
    nB = len(basis_names)

    for tag in targets:
        pts = values.get(tag, {})
        # collect the rows that have this weight AND a basis evaluation
        rows, rhs, pks = [], [], []
        for pk, orders in pts.items():
            if weight not in orders:
                continue
            s, t = kinematics.get(pk, (None, None))
            try:
                try:
                    bv = basis_fn(pk, s, t)
                except TypeError:
                    bv = basis_fn(pk)
            except KeyError:
                continue
            rows.append([mp.mpf(x) for x in bv])
            v = orders[weight]
            rhs.append(v.real if use_real else abs(v))
            pks.append(pk)

        N = len(rows)
        res = CVResult(master=tag, weight=weight,
                       n_points_used=N,
                       n_duplicates_excluded=n_duplicates_excluded,
                       basis_dimension=nB, basis_names=list(basis_names),
                       condition_number="n/a")
        if N <= nB:
            res.note = (f"underdetermined: N={N} <= basis_dim={nB} "
                        "(need at least basis_dim+1 for held-out CV)")
            out[tag] = res
            continue

        # ---- conditioning ------------------------------------------------
        A = mp.matrix(rows)
        b = mp.matrix([[x] for x in rhs])
        cond0 = _cond_svd(A)
        res.cond_before = mp.nstr(cond0, 4)
        res.ut_rotation_applied = ut_rotation is not None
        U = None
        nB_eff = nB
        if auto_cond and cond0 > cond_target:
            U, _, cond1, A = auto_condition(A, target_cond=cond_target,
                                            prec=max(30, dps - 10))
            if U is not None:
                nB_eff = U.cols
                rows = [[A[i, j] for j in range(nB_eff)] for i in range(N)]
                res.lll_applied = True
                res.cond_after = mp.nstr(cond1, 4)
                if nB_eff < nB:
                    res.note = (f"LLL exposed {nB - nB_eff} integer "
                                f"dependence(s); fitting on rank={nB_eff}")
        if not res.lll_applied:
            res.cond_after = res.cond_before
        # legacy field (kept for backward compat with existing reports)
        res.condition_number = res.cond_after

        # ---- full fit (for PSLQ) ----------------------------------------
        try:
            c_full, _ = _lstsq(A, b)
        except ZeroDivisionError:
            c_full = None
            res.note = ((res.note + "; " if res.note else "")
                        + "normal equations numerically singular "
                          f"(cond={res.cond_after})")

        # ---- leave-one-out ----------------------------------------------
        digits = []
        for k in range(N):
            ridx = [i for i in range(N) if i != k]
            Ak = mp.matrix([rows[i] for i in ridx])
            bk = mp.matrix([[rhs[i]] for i in ridx])
            try:
                ck, _ = _lstsq(Ak, bk)
            except Exception:            # singular sub-system
                digits.append(0.0)
                continue
            pred = sum(ck[j] * rows[k][j] for j in range(nB_eff))
            err = abs(pred - rhs[k])
            ref = max(abs(rhs[k]), mp.mpf(1))
            d = float(-mp.log10(err / ref)) if err > 0 else float(dps)
            digits.append(d)

        res.heldout_digits = [round(d, 2) for d in digits]
        res.min = round(min(digits), 2)
        res.mean = round(sum(digits) / len(digits), 2)
        res.certified = res.min >= gate

        # ---- rationalise (LLL-preconditioned) ---------------------------
        if res.certified and c_full is not None:
            from fractions import Fraction
            tol = max(10, int(res.min * 0.6))
            methods = set()
            if U is not None:
                # Rationalise in the REDUCED basis (well-conditioned, so
                # c_full is accurate), then map back via U in EXACT integer/
                # Fraction arithmetic — avoids the catastrophic cancellation
                # that U·c_full suffers at float precision when |U| is large.
                Umax = max(abs(U[i, j]) for i in range(nB)
                           for j in range(nB_eff))
                bits = int(mp.log(max(Umax, 2), 2)) + 24
                cU_rat = []
                for j in range(nB_eff):
                    s, m = _pslq_rational(c_full[j], tol,
                                          max_coeff_bits=max(48, 2 * bits))
                    if s is None:
                        cU_rat = None
                        break
                    methods.add(m)
                    cU_rat.append(Fraction(s))
                if cU_rat is not None:
                    Uint = [[int(mp.nint(U[i, j])) for j in range(nB_eff)]
                            for i in range(nB)]
                    c_orig = [sum(Uint[i][j] * cU_rat[j]
                                  for j in range(nB_eff)) for i in range(nB)]
                    rats = [(f"{c.numerator}/{c.denominator}"
                             if c.denominator != 1 else str(c.numerator))
                            for c in c_orig]
                else:
                    rats = [None] * nB
            else:
                pairs = [_pslq_rational(c_full[j], tol) for j in range(nB)]
                rats = [p[0] for p in pairs]
                methods = {p[1] for p in pairs if p[1]}
            res.pslq_rationals = rats if all(r is not None for r in rats) \
                else None
            res.pslq_method = ("lll" if methods <= {"lll", "exact"} and "lll" in methods
                               else "/".join(sorted(m for m in methods if m))
                               or None)
        out[tag] = res
    return out


# ====================================================================== CLI
def _main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--oracle", nargs="+", required=True,
                   help="oracle data dir(s)")
    p.add_argument("--basis", required=True,
                   help="basis_values.json (see module docstring)")
    p.add_argument("--weight", type=int, required=True,
                   help="ε-order to certify")
    p.add_argument("--gate", type=int, default=30)
    p.add_argument("--dps", type=int, default=80)
    p.add_argument("--integrity", default=None,
                   help="integrity_report.json from oracle_integrity.py; "
                        "if absent, a scan is run on the fly")
    p.add_argument("--report", default="cv_report.json")
    p.add_argument("--masters", nargs="*", default=None,
                   help="restrict to these master tags")
    a = p.parse_args(argv)

    mp.mp.dps = a.dps

    # 1. integrity → exclusion set
    if a.integrity and os.path.exists(a.integrity):
        irep = json.load(open(a.integrity))
        excl = set()
        for g in irep.get("duplicates", []):
            excl.update(g["files"][1:])
        for s in irep.get("suspicious_zero", []):
            excl.add(s["file"])
        for u in irep.get("unreadable", []):
            excl.add(u["file"])
    else:
        rep = _scan_integrity(a.oracle)
        excl = rep.quarantine_set()
    n_excl = len(excl)

    # 2. load oracle values
    values, kin, _ = load_oracle_values(a.oracle, exclude_files=excl)

    # 3. basis
    basis_fn, names = load_basis_json(a.basis, weight=a.weight)

    # 4. CV
    results = heldout_certify(values, kin, weight=a.weight,
                              basis_fn=basis_fn, basis_names=names,
                              gate=a.gate, dps=a.dps,
                              n_duplicates_excluded=n_excl,
                              masters=a.masters)

    # 5. report
    out = {tag: r.to_json() for tag, r in results.items()}
    summary = {
        "weight": a.weight, "gate": a.gate, "dps": a.dps,
        "n_masters": len(out),
        "n_certified": sum(1 for r in results.values() if r.certified),
        "n_duplicates_excluded": n_excl,
        "basis_dimension": len(names),
    }
    with open(a.report, "w") as fh:
        json.dump({"summary": summary, "results": out}, fh, indent=1)

    print(f"[heldout_cv] weight={a.weight}  gate={a.gate}d  dps={a.dps}")
    print(f"  excluded {n_excl} duplicate/zero file(s) before fitting")
    print(f"  basis dim = {len(names)}  ({', '.join(names)})")
    for tag in sorted(results):
        r = results[tag]
        flag = "CERT " if r.certified else "     "
        rat = ("  rats=" + ",".join(r.pslq_rationals)
               if r.pslq_rationals else "")
        if r.note:
            print(f"  {flag}{tag:10s}  {r.note}")
        else:
            cond_str = (f"cond={r.cond_before}→{r.cond_after} [LLL]"
                        if r.lll_applied else f"cond={r.cond_after}")
            print(f"  {flag}{tag:10s}  N={r.n_points_used:3d}  "
                  f"min={r.min:6.1f}d  mean={r.mean:6.1f}d  "
                  f"{cond_str}{rat}")
    print(f"  CERTIFIED: {summary['n_certified']}/{summary['n_masters']}")
    print(f"  wrote {a.report}")
    return 0


if __name__ == "__main__":
    sys.exit(_main())

#!/usr/bin/env python3
r"""
ansatzer.py — sequential-discontinuity / coproduct COLLAPSE predictor.

PURPOSE.  BEFORE oracle-farming a graph, predict the residual symbol-ansatz
dimension per weight per master after the standard big structural cuts.  The
residual dimension IS the number of unknowns the value-fit must determine, hence
(× safety factor) the number of high-precision oracle points the farm needs.
A 7-letter alphabet can collapse to a couple of unknowns per master at w≤3
("few points"), while a 2-mass 5pt non-planar family (19-letter,
K=[9,28,62,120]) leaves ≈79 points to farm.

THE MATH.
At weight w the raw tensor symbol space over an n-letter alphabet A has dim n^w.
Six linear constraint families cut it down, IN ORDER:
  1. INTEGRABILITY (dA=A∧A): for every adjacent slot k, every environment, the
     wedge Σ c · dlog L_{a} ∧ dlog L_{b} vanishes as a 2-form → an exact-ℚ (or
     mod-p) nullspace.  This is the engine: an exact-ℚ or mod-p nullspace over the wedge
     constraint rows.
  2. FIRST-ENTRY: leading letter ∈ {Mandelstams/masses/physical thresholds}.
  3. LAST-ENTRY: trailing letter ∈ {restricted set from the connection / disc}.
  4. STEINMANN: double discontinuities in overlapping channels vanish, so two
     ADJACENT entries may not carry disjoint channel-tag sets.  Letters sharing
     a tag are same-channel (iterated discs in one channel are allowed);
     untagged letters are unconstrained.  Active whenever the alphabet supplies
     "channel" tags (like the coproduct cut activating on a supplied connection).
  5. PARITY/REALITY: if a per-letter ℤ₂ parity is supplied (radicals, σ-symmetry),
     project to the master's target sector.
  6. COPRODUCT / Δ^{(w-1,1)} (sequential disc): if the canonical DE is
        dF = ε (Σ_a A_a dlog ℓ_a) · F
     then  Disc_a F^{(w)} ∝ A_a · F^{(w-1)} ; for master m the LAST letter a is
     allowed only if row m of A_a is non-zero, AND the (w-1)-prefix lies in the
     already-cut weight-(w-1) space.  Recursive.

The predictor reports the surviving dimension after EACH cut and the final
RESIDUAL → pts_needed → verdict {COLLAPSES | FULL}.

ENGINES:
  integrability engine — recursive integrable-symbol builder
  exact-Q engine       — exact rational collapse counts
  mod-p backend        — sampled modular collapse counts
  alphabet layer       — alphabet + first/last-entry classification

CLI:
  python ansatzer.py --alphabet alphabet.json [--connection conn.json] \
                               --wmax 4 --out collapse_report.json
"""
from __future__ import annotations
import argparse, itertools, json, os, random, signal, sys, time
from fractions import Fraction
from dataclasses import dataclass, field

import numpy as np
import sympy as sp


# =========================================================================
#  Timeout / honesty machinery (per-cut hard cap)
# =========================================================================
class CutTimeout(Exception):
    pass

def _alarm(signum, frame):  # pragma: no cover
    raise CutTimeout()


# =========================================================================
#  Alphabet model
#
#  Input JSON schema (the "collapse" schema).  Any landau_alphabet.py output
#  (alphabet/first_entry/last_entry/graph.kinematics.invariants) is auto-adapted.
#
#  {
#    "name": "myfamily",
#    "vars": ["s","t","m2"],
#    "letters": {                           # name -> sympy-parsable poly/expr
#       "s": "s", "t": "t", "u": "-s-t", "sm4": "s-4*m2", ...
#    },
#    "first_entry": ["s","t","sm4","tm4"],  # names; default = ALL
#    "last_entry":  ["s","t","sm4","tm4","rs","rt"],  # names; default = ALL
#    "parity": { "rs":[1], "rt":[1] },      # optional ℤ₂^k vector per letter (default 0)
#    "parity_target": [0],                  # optional; default all-even
#    "channel": { "s":["S"], "sm4":["S"], ... },  # optional Steinmann channel-tags
#    "substitutions": {"m2":"1"}            # optional; applied to letters before dlog
#  }
# =========================================================================
@dataclass
class Alphabet:
    name: str
    vars: list[sp.Symbol]
    letter_names: list[str]
    letter_exprs: list[sp.Expr]
    first_entry: set[int]
    last_entry: set[int]
    parity: list[tuple[int, ...]]
    parity_target: tuple[int, ...]
    channel: list[set[str]]
    raw: dict = field(default_factory=dict)

    @property
    def n(self):
        return len(self.letter_names)


def load_alphabet(path_or_dict) -> Alphabet:
    if isinstance(path_or_dict, dict):
        d = dict(path_or_dict)
    else:
        with open(path_or_dict) as fh:
            d = json.load(fh)
    # ---- adapt landau_alphabet.py output -------------------------------
    if "alphabet" in d and "letters" not in d:
        inv = (d.get("graph", {}).get("kinematics", {}).get("invariants")
               or d.get("invariants") or [])
        d = {
            "name": d.get("name", "graph"),
            "vars": inv if inv else sorted({str(s) for L in d["alphabet"]
                                            for s in sp.sympify(L).free_symbols}),
            "letters": {L: L for L in d["alphabet"]},
            "first_entry": d.get("first_entry") or list(d["alphabet"]),
            "last_entry": d.get("last_entry") or list(d["alphabet"]),
            "_source": "landau_alphabet.py output (auto-adapted)",
        }
    # ---- parse ----------------------------------------------------------
    name = d.get("name", "graph")
    var_syms = [sp.Symbol(v, real=True) for v in d["vars"]]
    locs = {str(v): v for v in var_syms}
    subs = {sp.Symbol(k): sp.sympify(v, locals=locs)
            for k, v in d.get("substitutions", {}).items()}
    Lnames = list(d["letters"].keys())
    Lexprs = [sp.expand(sp.sympify(d["letters"][n], locals=locs).subs(subs))
              for n in Lnames]
    idx = {n: i for i, n in enumerate(Lnames)}
    fe = set(idx[n] for n in d.get("first_entry", Lnames) if n in idx)
    le_raw = d.get("last_entry") or []
    le = set(idx[n] for n in le_raw if n in idx) or set(range(len(Lnames)))
    # parity
    par_raw = d.get("parity", {})
    krad = max((len(v) for v in par_raw.values()), default=0)
    parity = [tuple(par_raw.get(n, [0]*krad)) if krad else ()
              for n in Lnames]
    ptarget = tuple(d.get("parity_target", [0]*krad)) if krad else ()
    # channel tags (for the Steinmann cut; optional).  A tag entry may be a
    # list of tags or a single bare tag string.
    chan_raw = d.get("channel", {})
    channel = []
    for n in Lnames:
        tags = chan_raw.get(n, [])
        if isinstance(tags, str):
            tags = [tags]
        channel.append(set(tags))
    return Alphabet(name, var_syms, Lnames, Lexprs, fe, le, parity, ptarget,
                    channel, raw=d)


# =========================================================================
#  Connection model (optional).
#
#  conn.json schema:
#    {
#      "masters": ["m1","m2",...],                        # length M
#      "A": { "<letter_name>": [[i,j,"q"],...], ... }     # sparse M×M, rationals
#    }
#  or  {"A": {"<letter>": [[row list of M rationals]...]}}  (dense)
#
#  The COPRODUCT cut uses ONLY the SUPPORT structure of A_a (which masters feed
#  which masters with last-letter a) — it does not need the full numeric DE.
# =========================================================================
@dataclass
class Connection:
    masters: list[str]
    A: dict[str, np.ndarray]       # letter -> M×M boolean support matrix

    @property
    def M(self):
        return len(self.masters)


def load_connection(path_or_dict, letter_names) -> Connection | None:
    if path_or_dict is None:
        return None
    if isinstance(path_or_dict, dict):
        d = path_or_dict
    else:
        with open(path_or_dict) as fh:
            d = json.load(fh)
    masters = d.get("masters") or [f"m{i}" for i in range(len(next(iter(d["A"].values()))))]
    M = len(masters)
    A: dict[str, np.ndarray] = {}
    for ln, mat in d["A"].items():
        if ln not in letter_names:
            continue
        S = np.zeros((M, M), dtype=bool)
        if mat and isinstance(mat[0], (list, tuple)) and len(mat[0]) == 3 \
           and not isinstance(mat[0][2], (list, tuple)):
            for (i, j, v) in mat:
                if str(v) not in ("0", "0/1"):
                    S[int(i), int(j)] = True
        else:  # dense
            arr = np.array([[Fraction(str(x)) != 0 for x in row] for row in mat])
            S[:arr.shape[0], :arr.shape[1]] = arr
        A[ln] = S
    return Connection(masters, A)


# =========================================================================
#  EXACT-ℚ backend (Fraction RREF)
# =========================================================================
def _rref_nullspace_Q(rows: list[dict[int, Fraction]], N: int):
    """rows = list of sparse {col: Fraction} → (rank, basis:list[list[Fraction]])."""
    if not rows:
        return 0, [[Fraction(int(j == c)) for c in range(N)] for j in range(N)]
    A = [[Fraction(0)]*N for _ in rows]
    for r, row in enumerate(rows):
        for c, v in row.items():
            A[r][c] = v
    nrows = len(A); pivots = []; r = 0
    for col in range(N):
        piv = None
        for rr in range(r, nrows):
            if A[rr][col] != 0:
                piv = rr; break
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        pv = A[r][col]; A[r] = [x/pv for x in A[r]]
        for rr in range(nrows):
            if rr != r and A[rr][col] != 0:
                f = A[rr][col]; A[rr] = [a - f*b for a, b in zip(A[rr], A[r])]
        pivots.append(col); r += 1
        if r == nrows:
            break
    pivset = set(pivots); free = [c for c in range(N) if c not in pivset]
    basis = []
    for fc in free:
        v = [Fraction(0)]*N; v[fc] = Fraction(1)
        for ri, pc in enumerate(pivots):
            v[pc] = -A[ri][fc]
        basis.append(v)
    return len(pivots), basis


def _intersect_Q(basis, rows, N):
    """Intersect span(basis) with the nullspace of `rows`. Returns new basis."""
    if not basis or not rows:
        return basis
    k = len(basis); proj = []
    for row in rows:
        pr = {}
        for j in range(k):
            acc = Fraction(0)
            bj = basis[j]
            for c, val in row.items():
                acc += val * bj[c]
            if acc != 0:
                pr[j] = acc
        if pr:
            proj.append(pr)
    if not proj:
        return basis
    _, sub = _rref_nullspace_Q(proj, k)
    out = []
    for y in sub:
        v = [Fraction(0)]*N
        for j in range(k):
            if y[j] != 0:
                bj = basis[j]
                for c in range(N):
                    v[c] += y[j]*bj[c]
        out.append(v)
    return out


# =========================================================================
#  MOD-P backend (numpy int64)
#  Used when N or #rows is large (auto-switch); reports rank_method='modp'.
# =========================================================================
PRIME = 2147483647  # 2^31-1

def _row_echelon_modp(A: np.ndarray) -> tuple[np.ndarray, list[int]]:
    """In-place row-echelon over F_p. Returns (pivot rows only, pivot cols)."""
    nrows, N = A.shape
    pr = 0; pivcols = []
    for col in range(N):
        piv = -1
        for r in range(pr, nrows):
            if A[r, col] % PRIME != 0:
                piv = r; break
        if piv < 0:
            continue
        A[[pr, piv]] = A[[piv, pr]]
        inv = pow(int(A[pr, col]), PRIME - 2, PRIME)
        A[pr] = (A[pr] * inv) % PRIME
        colvals = A[:, col].copy(); colvals[pr] = 0
        mask = colvals % PRIME != 0
        if mask.any():
            A[mask] = (A[mask] - colvals[mask].reshape(-1, 1) * A[pr]) % PRIME
        pivcols.append(col); pr += 1
        if pr >= nrows:
            break
    return A[:pr].copy(), pivcols


def _rank_modp(rows: list[dict[int, int]], N: int, batch: int = 4096) -> int:
    """Incremental mod-p rank: process rows in batches, carry forward the
    row-echelon pivot rows.  Memory O((batch+rank)×N) instead of O(nrows×N);
    enables nrows >> N (the w≥3 large-alphabet regime)."""
    if not rows:
        return 0
    carry = np.zeros((0, N), dtype=np.int64)
    for start in range(0, len(rows), batch):
        chunk = rows[start:start+batch]
        B = np.zeros((carry.shape[0] + len(chunk), N), dtype=np.int64)
        B[:carry.shape[0]] = carry
        for r, d in enumerate(chunk):
            for c, v in d.items():
                B[carry.shape[0] + r, c] = int(v) % PRIME
        carry, _ = _row_echelon_modp(B)
        if carry.shape[0] >= N:
            return N
    return carry.shape[0]


# =========================================================================
#  dlog forms + wedge table
# =========================================================================
def _dlog_forms(alpha: Alphabet):
    """forms[i] = tuple of len(vars) sympy rationals: ∂_v log ℓ_i."""
    forms = []
    for e in alpha.letter_exprs:
        forms.append(tuple(sp.cancel(sp.together(sp.diff(e, v)/e)) for v in alpha.vars))
    return forms


def _wedge_exact(alpha: Alphabet, forms):
    """Exact wedge on a SINGLE global cleared denominator D (the symbols.jl
    `D = ∏F²` trick): W[(a,b)] = list[{mon→Fraction}] over the (i<j) 2-form
    components, where every entry is the polynomial numerator of
    D · (dlog ℓ_a ∧ dlog ℓ_b)_{ij}.  Because D is COMMON to all pairs, the
    row-builder can compare monomial coefficients directly (no per-group lcm)."""
    V = alpha.vars
    pairs = [(i, j) for i in range(len(V)) for j in range(i+1, len(V))]
    n = alpha.n
    # global denominator: lcm over all (a,b,component) of the wedge denominator.
    # symbols.jl uses ∏F² which is an upper bound; we take the actual lcm.
    raw: dict[tuple[int, int], list[tuple[sp.Expr, sp.Expr]]] = {}
    D = sp.Integer(1)
    for a in range(n):
        for b in range(a+1, n):
            comps = []
            for (i, j) in pairs:
                J = sp.cancel(sp.together(forms[a][i]*forms[b][j]
                                          - forms[a][j]*forms[b][i]))
                num, den = sp.fraction(J)
                comps.append((sp.expand(num), sp.expand(den)))
                if den != 1:
                    D = sp.lcm(D, den)
            raw[(a, b)] = comps
    D = sp.expand(D)
    W: dict[tuple[int, int], list[dict]] = {}
    for ab, comps in raw.items():
        out = []
        for (num, den) in comps:
            if num == 0:
                out.append({}); continue
            scaled = sp.expand(num * sp.cancel(D/den))
            pol = sp.Poly(scaled, *V)
            d = {}
            for mon, c in pol.terms():
                r = sp.Rational(c)
                d[mon] = Fraction(int(r.p), int(r.q))
            out.append(d)
        W[ab] = out
    return W, len(pairs)


def _wedge_modp(alpha: Alphabet, forms, npts: int = 6, seed: int = 11):
    """Multi-point mod-p wedge evaluation; returns W_pts, ncomp where
    W_pts[k][(a,b)] = list[int] of length ncomp (one per 2-form component)."""
    V = alpha.vars
    pairs = [(i, j) for i in range(len(V)) for j in range(i+1, len(V))]
    ncomp = len(pairs)
    n = alpha.n
    rng = random.Random(seed)
    W_pts = []
    for _ in range(npts):
        pt = {v: rng.randint(2, PRIME - 2) for v in V}
        D = []
        for fr in forms:
            vec = []
            for c in fr:
                num, den = sp.fraction(sp.together(c))
                nn = int(num.subs(pt)) % PRIME
                dd = int(den.subs(pt)) % PRIME
                vec.append(nn * pow(dd, PRIME - 2, PRIME) % PRIME if dd else 0)
            D.append(vec)
        W = {}
        for a in range(n):
            for b in range(n):
                if a == b:
                    continue
                W[(a, b)] = [(D[a][i]*D[b][j] - D[a][j]*D[b][i]) % PRIME
                             for (i, j) in pairs]
        W_pts.append(W)
    return W_pts, ncomp


# =========================================================================
#  1.  symbol_space_dims — raw n^w
# =========================================================================
def symbol_space_dims(alpha: Alphabet, weight_max: int) -> dict[int, int]:
    return {w: alpha.n ** w for w in range(weight_max + 1)}


# =========================================================================
#  2.  integrability_cut(alpha, w, *, method) → (basis_or_None, dim_after, info)
#      Exact-ℚ path returns the explicit basis; mod-p path returns dim only
#      (basis=None) and downstream cuts re-impose all rows together.
# =========================================================================
def integrability_cut(alpha: Alphabet, w: int, *, method: str = "auto",
                      timeout_s: int = 0, _cache: dict | None = None):
    n = alpha.n
    words = list(itertools.product(range(n), repeat=w))
    widx = {wd: i for i, wd in enumerate(words)}
    N = len(words)
    if w <= 1:
        basis = [[Fraction(int(j == c)) for c in range(N)] for j in range(N)]
        return words, basis, N, {"rank_method": "trivial"}
    if method == "auto":
        method = "exact_Q" if N <= 700 else "modp"
    forms = (_cache or {}).get("forms")
    if forms is None:
        forms = _dlog_forms(alpha)
        if _cache is not None:
            _cache["forms"] = forms

    if timeout_s:
        signal.signal(signal.SIGALRM, _alarm); signal.alarm(timeout_s)
    try:
        if method == "exact_Q":
            W, ncomp = (_cache or {}).get("Wexact") or _wedge_exact(alpha, forms)
            if _cache is not None:
                _cache["Wexact"] = (W, ncomp)
            rows = _integ_rows_exact(words, widx, W, ncomp, w, alpha.vars)
            _, basis = _rref_nullspace_Q(rows, N)
            return words, basis, len(basis), {"rank_method": "exact_Q",
                                              "n_rows": len(rows)}
        else:  # mod-p
            Wp, ncomp = (_cache or {}).get("Wmodp") or _wedge_modp(alpha, forms)
            if _cache is not None:
                _cache["Wmodp"] = (Wp, ncomp)
            rows = _integ_rows_modp(words, widx, Wp, ncomp, w)
            rk = _rank_modp(rows, N)
            return words, None, N - rk, {"rank_method": "modp",
                                         "n_rows": len(rows), "npts": len(Wp)}
    except CutTimeout:
        return words, None, N, {"rank_method": method, "timeout": True,
                                "bound": "≥ (raw; integrability TIMED OUT)"}
    finally:
        if timeout_s:
            signal.alarm(0)


def _integ_rows_exact(words, widx, W, ncomp, w, V):
    """For each (slot k, environment, 2-form component, monomial), one ℚ-row.
    Wedge numerators in W are already on a common global denominator, so
    comparison is direct monomial bookkeeping (mirrors symbols.jl)."""
    rows = []
    for k in range(w - 1):
        groups: dict[tuple, list[int]] = {}
        for ci, wd in enumerate(words):
            groups.setdefault(wd[:k] + wd[k+2:], []).append(ci)
        for env, cols in groups.items():
            buckets: dict[tuple, dict[int, Fraction]] = {}
            for ci in cols:
                a, b = words[ci][k], words[ci][k+1]
                if a == b:
                    continue
                lo, hi, sgn = (a, b, 1) if a < b else (b, a, -1)
                comps = W[(lo, hi)]
                for cc in range(ncomp):
                    for mon, cf in comps[cc].items():
                        key = (cc, mon)
                        d = buckets.setdefault(key, {})
                        d[ci] = d.get(ci, Fraction(0)) + sgn * cf
            for d in buckets.values():
                d = {c: v for c, v in d.items() if v != 0}
                if d:
                    rows.append(d)
    return rows


def _integ_rows_modp(words, widx, W_pts, ncomp, w):
    rows = []
    for Wp in W_pts:
        for k in range(w - 1):
            groups: dict[tuple, list[int]] = {}
            for ci, wd in enumerate(words):
                groups.setdefault(wd[:k] + wd[k+2:], []).append(ci)
            for env, cols in groups.items():
                comp = [dict() for _ in range(ncomp)]
                for ci in cols:
                    a, b = words[ci][k], words[ci][k+1]
                    if a == b:
                        continue
                    wv = Wp[(a, b)]
                    for cc in range(ncomp):
                        if wv[cc]:
                            comp[cc][ci] = (comp[cc].get(ci, 0) + wv[cc]) % PRIME
                for cr in comp:
                    if cr:
                        rows.append(cr)
    return rows


# =========================================================================
#  3.  entry_cut(basis, words, allowed_first, allowed_last) — exact-ℚ projection
# =========================================================================
def entry_cut(basis, words, allowed_first: set[int], allowed_last: set[int]):
    N = len(words)
    fe_rows = [{i: Fraction(1)} for i, wd in enumerate(words)
               if wd[0] not in allowed_first]
    le_rows = [{i: Fraction(1)} for i, wd in enumerate(words)
               if wd[-1] not in allowed_last]
    b1 = _intersect_Q(basis, fe_rows, N)
    b2 = _intersect_Q(b1, le_rows, N)
    return b1, b2


# =========================================================================
#  3b. steinmann_cut — adjacent-entry channel veto from the "channel" tags.
#
#  Steinmann relations kill double discontinuities in overlapping channels:
#  a word is forbidden iff some ADJACENT pair (slot k, k+1) carries two
#  non-empty channel-tag sets with EMPTY intersection.  Letters sharing a tag
#  are discontinuities in the same channel (iterated discs allowed); letters
#  with no tag are unconstrained.  Vacuous when the alphabet has no tags.
# =========================================================================
def steinmann_ok(wd, channel) -> bool:
    """True iff word `wd` (tuple of letter indices) has no adjacent pair of
    letters whose channel-tag sets are both non-empty and disjoint."""
    for a, b in zip(wd, wd[1:]):
        ca, cb = channel[a], channel[b]
        if ca and cb and not (ca & cb):
            return False
    return True


def steinmann_cut(basis, words, alpha: Alphabet):
    if not any(alpha.channel):
        return basis, "vacuous (no channel tags)"
    rows = [{i: Fraction(1)} for i, wd in enumerate(words)
            if not steinmann_ok(wd, alpha.channel)]
    b = _intersect_Q(basis, rows, len(words))
    ntag = sum(1 for c in alpha.channel if c)
    return b, (f"adjacent overlapping-channel veto; "
               f"{ntag}/{alpha.n} letters tagged")


# =========================================================================
#  4.  coproduct_cut — Δ^{(w-1,1)} via the connection.
#
#  For each master m, letter a is allowed in the LAST slot iff row m of A_a is
#  non-zero (Disc_a F_m^{(w)} = (A_a · F^{(w-1)})_m ≠ 0).  Recursively, the
#  (w-1)-prefix of every word in the support of the symbol must lie in the
#  already-cut weight-(w-1) FAMILY space (this is integrability re-applied to
#  the prefix; in the exact-ℚ path the integrable basis already enforces it).
#
#  If `conn` is None → return basis unchanged with note 'SKIPPED'.
# =========================================================================
def coproduct_cut(basis, words, alpha: Alphabet, conn: Connection | None,
                  master: str | None):
    if conn is None or master is None:
        return basis, "SKIPPED (no connection supplied)"
    if master not in conn.masters:
        return basis, f"SKIPPED (master '{master}' not in connection)"
    m = conn.masters.index(master)
    allowed_last = set()
    for ln, S in conn.A.items():
        if ln in alpha.letter_names and S[m, :].any():
            allowed_last.add(alpha.letter_names.index(ln))
    # any letter NOT keyed in conn.A but in last_entry stays allowed (honesty:
    # missing A_a ⇒ no constraint from the connection on that letter)
    for i in range(alpha.n):
        if alpha.letter_names[i] not in conn.A:
            allowed_last.add(i)
    rows = [{i: Fraction(1)} for i, wd in enumerate(words)
            if wd[-1] not in allowed_last]
    b = _intersect_Q(basis, rows, len(words))
    note = (f"per-master last-letter from rows of A_a; "
            f"{len(allowed_last)}/{alpha.n} letters allowed for {master}")
    return b, note


# =========================================================================
#  Parity cut (optional; ℤ₂^k)
# =========================================================================
def parity_cut(basis, words, alpha: Alphabet):
    if not alpha.parity_target and not any(any(p) for p in alpha.parity):
        return basis, "vacuous"
    k = len(alpha.parity_target) if alpha.parity_target else \
        (len(alpha.parity[0]) if alpha.parity and alpha.parity[0] else 0)
    if k == 0:
        return basis, "vacuous"
    tgt = alpha.parity_target or tuple([0]*k)
    rows = []
    for i, wd in enumerate(words):
        p = [0]*k
        for a in wd:
            for j in range(k):
                p[j] ^= alpha.parity[a][j]
        if tuple(p) != tgt:
            rows.append({i: Fraction(1)})
    return _intersect_Q(basis, rows, len(words)), f"target={list(tgt)}"


# =========================================================================
#  mod-p path — for large alphabets.  Monomial cuts (first/last/parity/coproduct
#  last-letter) are applied as a WORD PRE-FILTER before integrability, so the
#  matrix is built only on the surviving words (19⁴=130k columns is
#  infeasible otherwise).  Returns the per-cut counts.
# =========================================================================
def _frac_mod_p(q: Fraction) -> int:
    return int(q.numerator) * pow(int(q.denominator), PRIME - 2, PRIME) % PRIME


def _modp_layered(alpha: Alphabet, w: int, _cache: dict, *,
                  conn: Connection | None = None, master: str | None = None,
                  steinmann: bool = True, timeout_s: int = 0):
    n = alpha.n
    forms = _cache.get("forms") or _dlog_forms(alpha)
    _cache["forms"] = forms
    # Wedge backend: EXACT (global-denom polynomial monomials) when n_vars ≤ 3
    # — degree stays bounded so the row set is complete.  For many-variable
    # alphabets (5pt: 7 vars) the global-denom blows up, so use MULTI-POINT
    # mod-p sampling with a 6→12-point convergence check (validated
    # point/prime-stable).
    wedge_mode = _cache.get("wedge_mode")
    if wedge_mode is None:
        # COLLAPSE_WEDGE_MODE=sample|exact overrides the auto rule: _wedge_exact's
        # global-denom sympy expansion is unbounded for large alphabets (32-letter
        # T4 hung >10min, OUTSIDE the per-weight alarm) — sampled wedge is the
        # validated sampled-wedge regime and converges via npts doubling.
        wedge_mode = os.environ.get("COLLAPSE_WEDGE_MODE") or \
            ("exact" if len(alpha.vars) <= 3 else "sample")
        _cache["wedge_mode"] = wedge_mode
    if wedge_mode == "exact":
        W, ncomp = _cache.get("Wexact") or _wedge_exact(alpha, forms)
        _cache["Wexact"] = (W, ncomp)
    else:
        # adaptive sampling: npts ∈ {6,12,24,48,...} until rank stabilizes.
        # Cache the largest set; smaller sets are prefixes (same seed).
        npts_max = int(os.environ.get("COLLAPSE_NPTS_MAX", "12"))
        Wmax, ncomp = _cache.get("Wmodp_max") or _wedge_modp(alpha, forms,
                                                             npts=npts_max, seed=11)
        _cache["Wmodp_max"] = (Wmax, ncomp)

    # per-master coproduct last-letter set
    cop_allowed = None
    if conn is not None and master is not None and master in conn.masters:
        m = conn.masters.index(master)
        cop_allowed = set()
        for ln, S in conn.A.items():
            if ln in alpha.letter_names and S[m, :].any():
                cop_allowed.add(alpha.letter_names.index(ln))
        for i in range(n):
            if alpha.letter_names[i] not in conn.A:
                cop_allowed.add(i)

    # parity predicate
    have_par = bool(alpha.parity_target) or any(any(p) for p in alpha.parity)
    if have_par:
        k = len(alpha.parity_target) or len(alpha.parity[0])
        tgt = alpha.parity_target or tuple([0]*k)
        def par_ok(wd):
            p = [0]*k
            for a in wd:
                for j in range(k):
                    p[j] ^= alpha.parity[a][j]
            return tuple(p) == tgt
    else:
        par_ok = lambda wd: True

    # Steinmann adjacency predicate (vacuous without channel tags)
    if steinmann and any(alpha.channel):
        chan = alpha.channel
        st_ok = lambda wd: steinmann_ok(wd, chan)
    else:
        st_ok = lambda wd: True

    # word-layer counts.  When n is large (≳50) the naïve n^w enumeration is the
    # bottleneck (e.g. 437^4 ≈ 4×10^10).  Iterate FE × n^{w-2} × LE directly so
    # the loop body runs |survivors| times, not n^w.  Intermediate counts (n_fe,
    # n_le) are computed combinatorially since they're just monomial filters.
    raw = n ** w
    fe = sorted(alpha.first_entry); le = sorted(alpha.last_entry)
    n_fe = len(fe) * n ** (w - 1)
    n_le = len(fe) * (n ** (w - 2) if w >= 2 else 1) * len(le) if w >= 2 \
        else len(set(fe) & set(le))
    n_st = n_par = n_cop = 0
    survivors = []
    le_cop = [a for a in le if (cop_allowed is None or a in cop_allowed)]
    if w == 1:
        for a in sorted(set(fe) & set(le)):
            wd = (a,)
            n_st += 1          # no adjacent pair at w=1
            if not par_ok(wd):
                continue
            n_par += 1
            if cop_allowed is not None and a not in cop_allowed:
                continue
            n_cop += 1
            survivors.append(wd)
    else:
        mid_iter = itertools.product(range(n), repeat=w - 2) if w >= 3 else [()]
        for mid in mid_iter:
            for a in fe:
                for b in le:
                    wd = (a,) + mid + (b,)
                    if not st_ok(wd):
                        continue
                    n_st += 1
                    if not par_ok(wd):
                        continue
                    n_par += 1
                    if cop_allowed is not None and b not in cop_allowed:
                        continue
                    n_cop += 1
                    survivors.append(wd)
    N = len(survivors)
    info = {"rank_method": f"modp/{wedge_mode}-wedge", "n_survivor_words": N}
    if w <= 1 or N == 0:
        return dict(raw=raw, first=n_fe, last=n_le, steinmann=n_st,
                    parity=n_par, coproduct=n_cop, integrability=N,
                    RESIDUAL=N), info

    if timeout_s:
        signal.signal(signal.SIGALRM, _alarm); signal.alarm(timeout_s)
    try:
        widx = {wd: i for i, wd in enumerate(survivors)}
        if wedge_mode == "exact":
            rows_Q = _integ_rows_exact(survivors, widx, W, ncomp, w, alpha.vars)
            rows = [{c: _frac_mod_p(v) for c, v in r.items()} for r in rows_Q]
            rk = _rank_modp(rows, N)
            residual = N - rk
        else:
            # adaptive: double npts until rank stable or npts_max hit
            npts_seq = [6]
            while npts_seq[-1] < len(Wmax):
                npts_seq.append(min(2*npts_seq[-1], len(Wmax)))
            rk_prev = -1; rk = 0
            for np_ in npts_seq:
                rk = _rank_modp(_integ_rows_modp(survivors, widx,
                                                 Wmax[:np_], ncomp, w), N)
                info.setdefault("rank_by_npts", {})[np_] = rk
                if rk == rk_prev:
                    info["sample_converged"] = True
                    info["npts_converged"] = np_
                    break
                rk_prev = rk
            else:
                info["sample_converged"] = (len(npts_seq) == 1)
                if not info["sample_converged"]:
                    info["bound"] = (f"≤{N-rk} at {npts_seq[-1]}pt; sampling NOT "
                                     f"converged (set COLLAPSE_NPTS_MAX higher)")
            residual = N - rk
    except CutTimeout:
        info["timeout"] = True
        info["bound"] = f"≥? (integrability TIMED OUT on {N} words)"
        residual = N
    finally:
        if timeout_s:
            signal.alarm(0)
    # NOTE: in mod-p mode the "integrability" column reports the residual AFTER
    # all monomial pre-filters + integrability (i.e. the same number as RESIDUAL);
    # the per-cut intermediate dims are word-counts, not subspace dims.
    return dict(raw=raw, first=n_fe, last=n_le, steinmann=n_st, parity=n_par,
                coproduct=n_cop, integrability=residual, RESIDUAL=residual), info


# =========================================================================
#  5.  predict_collapse — the driver
# =========================================================================
SAFETY = 1.3            # pts_needed = ceil(SAFETY * residual)
COLLAPSE_PER_MASTER = 12   # ≤ this at every w ⇒ "COLLAPSES"
FULL_THRESHOLD = 60        # ≥ this at any w ⇒ "FULL"


def predict_collapse(alphabet_json, *, connection=None, weight_max: int = 4,
                     masters: list[str] | None = None, method: str = "auto",
                     steinmann: bool = True, timeout_s: int = 120,
                     verbose: bool = False) -> dict:
    alpha = alphabet_json if isinstance(alphabet_json, Alphabet) else load_alphabet(alphabet_json)
    conn = connection if (connection is None or isinstance(connection, Connection)) \
        else load_connection(connection, alpha.letter_names)
    masters = masters or (conn.masters if conn else [None])
    raw = symbol_space_dims(alpha, weight_max)
    cache: dict = {}
    table = {"name": alpha.name, "n_letters": alpha.n,
             "letters": alpha.letter_names,
             "first_entry": sorted(alpha.letter_names[i] for i in alpha.first_entry),
             "last_entry": sorted(alpha.letter_names[i] for i in alpha.last_entry),
             "weight_max": weight_max,
             "per_weight": {}, "per_master": {}, "honesty": {}}

    # rank backend chosen PER WEIGHT under method="auto"
    table["honesty"]["rank_method"] = method
    table["honesty"]["coproduct"] = ("applied" if conn is not None
                                     else "SKIPPED (no connection supplied)")
    have_chan = any(alpha.channel)
    if not steinmann:
        table["honesty"]["steinmann"] = "SKIPPED (disabled)"
    elif have_chan:
        ntag = sum(1 for c in alpha.channel if c)
        table["honesty"]["steinmann"] = (f"applied (adjacent overlapping-channel "
                                         f"veto; {ntag}/{alpha.n} letters tagged)")
    else:
        table["honesty"]["steinmann"] = "vacuous (no channel tags)"
    apply_st = steinmann and have_chan
    EXACT_Q_WORD_CAP = 500      # exact-ℚ RREF wall-clock limit; mod-p above

    for w in range(weight_max + 1):
        if verbose:
            sys.stderr.write(f"[w={w}] raw={raw[w]} ...\n")
        rec = {"raw": raw[w]}
        if w == 0:
            rec.update(integrability=1, first=1, last=1, steinmann=1, parity=1,
                       coproduct=1, RESIDUAL=1, pts_needed=2,
                       rank_method="trivial")
            table["per_weight"][w] = rec
            continue
        t0 = time.time()
        # per-weight backend choice
        meth_w = method
        if method == "auto":
            meth_w = "exact_Q" if alpha.n ** w <= EXACT_Q_WORD_CAP else "modp"
        if meth_w == "exact_Q":
            words, basis, dim_int, info = integrability_cut(
                alpha, w, method="exact_Q", timeout_s=timeout_s, _cache=cache)
            rec["integrability"] = dim_int
            rec["rank_method"] = info["rank_method"]
            if info.get("timeout"):
                rec["bound"] = f"≥{dim_int}"
                rec.update(first=dim_int, last=dim_int, steinmann=dim_int,
                           parity=dim_int, coproduct=dim_int, RESIDUAL=dim_int,
                           pts_needed=int(np.ceil(SAFETY*dim_int)))
                table["per_weight"][w] = rec
                table["honesty"].setdefault("timeouts", []).append(w)
                continue
            b_fe, b_le = entry_cut(basis, words, alpha.first_entry, alpha.last_entry)
            rec["first"] = len(b_fe)
            rec["last"] = len(b_le)
            if apply_st:
                b_st, stnote = steinmann_cut(b_le, words, alpha)
            else:
                b_st = b_le
            rec["steinmann"] = len(b_st)
            b_par, parnote = parity_cut(b_st, words, alpha)
            rec["parity"] = len(b_par)
            # coproduct: report the FAMILY (master=None) and per-master below
            b_cop, copnote = coproduct_cut(b_par, words, alpha, conn, masters[0])
            rec["coproduct"] = len(b_cop)
            rec["RESIDUAL"] = len(b_cop)
            rec["secs"] = round(time.time() - t0, 2)
            # per-master
            for m in masters:
                bm, _ = coproduct_cut(b_par, words, alpha, conn, m)
                table["per_master"].setdefault(str(m), {})[w] = len(bm)
        else:  # mod-p: word pre-filter + integrability on survivors
            counts, info2 = _modp_layered(alpha, w, cache, conn=conn,
                                          master=masters[0], steinmann=steinmann,
                                          timeout_s=timeout_s)
            rec.update(counts)
            rec["rank_method"] = info2["rank_method"]
            if "sample_converged" in info2:
                rec["sample_converged"] = info2["sample_converged"]
            rec["secs"] = round(time.time() - t0, 2)
            if info2.get("timeout"):
                rec["bound"] = info2["bound"]
                table["honesty"].setdefault("timeouts", []).append(w)
            if conn is not None:
                for m in masters:
                    cm, _ = _modp_layered(alpha, w, cache, conn=conn, master=m,
                                          steinmann=steinmann,
                                          timeout_s=timeout_s)
                    table["per_master"].setdefault(str(m), {})[w] = cm["RESIDUAL"]
            else:
                table["per_master"].setdefault("None", {})[w] = rec["RESIDUAL"]
        rec["pts_needed"] = int(np.ceil(SAFETY * rec["RESIDUAL"]))
        table["per_weight"][w] = rec

    # ---- verdict --------------------------------------------------------
    res = [table["per_weight"][w]["RESIDUAL"] for w in range(1, weight_max+1)]
    rmax = max(res)
    table["RESIDUAL_max"] = rmax
    table["pts_needed_max"] = int(np.ceil(SAFETY * rmax))
    if rmax <= COLLAPSE_PER_MASTER:
        table["verdict"] = f"COLLAPSES (≤{rmax}/weight; few-point fit)"
    elif rmax >= FULL_THRESHOLD:
        table["verdict"] = f"FULL (~{table['pts_needed_max']} points needed at w≤{weight_max})"
    else:
        table["verdict"] = f"INTERMEDIATE ({rmax}/weight; ~{table['pts_needed_max']} pts)"
    return table


# =========================================================================
#  6.  compare(A, B) — side-by-side
# =========================================================================
def compare(report_A: dict, report_B: dict) -> str:
    wA = report_A["weight_max"]; wB = report_B["weight_max"]
    wmax = min(wA, wB)
    lines = []
    h = (f"{'':<14}| {report_A['name']:^28} || {report_B['name']:^28}")
    lines.append(h); lines.append("-"*len(h))
    lines.append(f"{'n_letters':<14}| {report_A['n_letters']:^28} || "
                 f"{report_B['n_letters']:^28}")
    lines.append(f"{'rank_method':<14}| {report_A['honesty']['rank_method']:^28} || "
                 f"{report_B['honesty']['rank_method']:^28}")
    lines.append(f"{'coproduct':<14}| {report_A['honesty']['coproduct']:^28} || "
                 f"{report_B['honesty']['coproduct']:^28}")
    lines.append("")
    cols = ["raw", "integrability", "first", "last", "steinmann", "coproduct",
            "RESIDUAL", "pts_needed"]
    head = "w | " + " | ".join(f"{c:>6}" for c in cols)
    lines.append(f"{head}   ||   {head}")
    for w in range(1, wmax+1):
        rA = report_A["per_weight"][w]; rB = report_B["per_weight"][w]
        rowA = f"{w} | " + " | ".join(f"{rA.get(c, '-'):>6}" for c in cols)
        rowB = f"{w} | " + " | ".join(f"{rB.get(c, '-'):>6}" for c in cols)
        lines.append(f"{rowA}   ||   {rowB}")
    lines.append("")
    lines.append(f"VERDICT: {report_A['name']:>12} → {report_A['verdict']}")
    lines.append(f"         {report_B['name']:>12} → {report_B['verdict']}")
    return "\n".join(lines)


# =========================================================================
#  pretty-print one report
# =========================================================================
def format_table(rep: dict) -> str:
    cols = ["raw", "integrability", "first", "last", "steinmann", "parity",
            "coproduct", "RESIDUAL", "pts_needed"]
    out = [f"=== {rep['name']}  ({rep['n_letters']} letters; "
           f"rank={rep['honesty']['rank_method']}; "
           f"steinmann={rep['honesty'].get('steinmann', '-')}; "
           f"coproduct={rep['honesty']['coproduct']}) ==="]
    out.append("weight | " + " | ".join(f"{c:>13}" for c in cols))
    out.append("-"*(9 + 16*len(cols)))
    for w in sorted(rep["per_weight"]):
        r = rep["per_weight"][w]
        out.append(f"  {w:>4} | " + " | ".join(f"{r.get(c,'-'):>13}" for c in cols))
    out.append("")
    out.append(f"VERDICT: {rep['verdict']}   "
               f"(max RESIDUAL = {rep['RESIDUAL_max']}, "
               f"pts_needed ≈ {rep['pts_needed_max']})")
    if "timeouts" in rep["honesty"]:
        out.append(f"HONESTY: integrability TIMED OUT at w={rep['honesty']['timeouts']}; "
                   f"residuals at those weights are ≥-bounds.")
    return "\n".join(out)


# =========================================================================
#  --per-orbit path — recursive integrable-symbol builder (per-orbit K_w).
#
#  Where predict_collapse reports the global residual per weight, this path
#  builds I_w orbit-by-orbit and reports the converged K_w per orbit with a
#  mod-p sample-wedge convergence certificate:
#     I_1 = FE;  I_w = ker(integ_{slot w-1}) on I_{w-1} ⊗ MID, last slot → LE at
#     the target weight.  Works at large n because dim(I_w) stays small even when
#     n^w is huge.  Input goes through load_alphabet, so the collapse schema,
#  the landau-alphabet-output auto-adaptation, the substitutions, and the
#  first/last-entry defaults all match the global path.  Channel tags are NOT
#  applied on this path (global table only).
#  Validated against predict_collapse on hexabox sub8 (K_w0..3 = [1,2,3,4]).
#  PRIME is the shared 2^31-1 modulus.
# =========================================================================
def dlog_forms_modp(letters, varsyms, npts, seed=11):
    """D[k][i] = (∂_v log ℓ_i)(pt_k) mod p, for i in letters, v in vars."""
    rng = random.Random(seed)
    locs = {str(v): v for v in varsyms}
    exprs = [sp.sympify(e, locals=locs) for e in letters.values()]
    forms = [[sp.together(sp.diff(e, v) / e) for v in varsyms] for e in exprs]
    pts = [{v: rng.randint(2, PRIME - 2) for v in varsyms} for _ in range(npts)]
    D = []
    for pt in pts:
        row = []
        for fr in forms:
            vec = []
            for c in fr:
                num, den = sp.fraction(c)
                nn = int(num.subs(pt)) % PRIME
                dd = int(den.subs(pt)) % PRIME
                vec.append(nn * pow(dd, PRIME - 2, PRIME) % PRIME if dd else 0)
            row.append(vec)
        D.append(row)
    return D  # D[pt][letter] = [dlog components]


def wedge_at_pt(Dpt, a, b, ncomp_pairs):
    da, db = Dpt[a], Dpt[b]
    return [(da[i] * db[j] - da[j] * db[i]) % PRIME for (i, j) in ncomp_pairs]


def _row_echelon(A):
    A = A.copy() % PRIME
    m, n = A.shape
    pr = 0; pivcols = []
    for col in range(n):
        piv = -1
        for r in range(pr, m):
            if A[r, col] % PRIME:
                piv = r; break
        if piv < 0:
            continue
        A[[pr, piv]] = A[[piv, pr]]
        inv = pow(int(A[pr, col]), PRIME - 2, PRIME)
        A[pr] = (A[pr] * inv) % PRIME
        cv = A[:, col].copy(); cv[pr] = 0
        mask = cv % PRIME != 0
        if mask.any():
            A[mask] = (A[mask] - cv[mask].reshape(-1, 1) * A[pr]) % PRIME
        pivcols.append(col); pr += 1
        if pr >= m:
            break
    return A[:pr], pivcols


def nullspace_modp(A, batch=2048):
    """A: (m, n) int64 mod p. Returns basis (k, n) of ker(A).
    When m >> n, first incrementally row-reduce in batches so the final
    elimination works on ≤n rows (O(n³) instead of O(m·n²))."""
    m_in, n = A.shape
    if m_in > n + batch:
        carry = np.zeros((0, n), dtype=np.int64)
        for s in range(0, m_in, batch):
            blk = np.vstack([carry, A[s:s+batch] % PRIME])
            carry, _ = _row_echelon(blk)
            if carry.shape[0] >= n:
                break
        A = carry
    else:
        A = A.copy() % PRIME
    m, n = A.shape
    pr = 0; pivcols = []
    for col in range(n):
        piv = -1
        for r in range(pr, m):
            if A[r, col] % PRIME:
                piv = r; break
        if piv < 0:
            continue
        A[[pr, piv]] = A[[piv, pr]]
        inv = pow(int(A[pr, col]), PRIME - 2, PRIME)
        A[pr] = (A[pr] * inv) % PRIME
        cv = A[:, col].copy(); cv[pr] = 0
        mask = cv % PRIME != 0
        if mask.any():
            A[mask] = (A[mask] - cv[mask].reshape(-1, 1) * A[pr]) % PRIME
        pivcols.append(col); pr += 1
    pivset = set(pivcols)
    free = [c for c in range(n) if c not in pivset]
    basis = np.zeros((len(free), n), dtype=np.int64)
    for k, fc in enumerate(free):
        basis[k, fc] = 1
        for ri, pc in enumerate(pivcols):
            basis[k, pc] = (-A[ri, fc]) % PRIME
    return basis


def recursive_Kw(alpha_json, wmax=4, npts_max=96):
    """`alpha_json`: path, dict, or Alphabet — anything load_alphabet accepts
    (collapse schema or landau-alphabet output; first/last-entry defaults and
    substitutions handled there, identically to the global path)."""
    alpha = alpha_json if isinstance(alpha_json, Alphabet) \
        else load_alphabet(alpha_json)
    varsyms = alpha.vars
    nvars = len(varsyms)
    pairs = [(i, j) for i in range(nvars) for j in range(i+1, nvars)]
    names = alpha.letter_names
    n = alpha.n
    FE = sorted(alpha.first_entry)
    LE = sorted(alpha.last_entry)
    LEset = set(LE)

    t0 = time.time()
    D = dlog_forms_modp(dict(zip(names, alpha.letter_exprs)), varsyms, npts_max)
    setup_s = time.time() - t0

    # I_1 = FE letters: basis over word-space {a : a∈FE}, last_letters = FE
    # Represent I_w as: (basis: K×len(words), words: list[tuple], last_of_word)
    words = [(a,) for a in FE]
    B = np.eye(len(words), dtype=np.int64)
    Kw = {0: 1, 1: len(set(FE) & set(LE))}  # K_w1 = |FE∩LE| (same as predictor)
    Kw_full = {1: len(FE)}  # before LE projection
    conv = {}

    for w in range(2, wmax + 1):
        # Tensor B (K × |words|) with letters → new words = words × tensor_set.
        # At w<wmax tensor with ALL n (carries forward); at w=wmax tensor with LE
        # only (we only need K_w^{LE}, and Ncand = K_prev×|LE| << K_prev×n).
        K_prev = B.shape[0]
        tensor_set = LE if w == wmax else list(range(n))
        nT = len(tensor_set)
        new_words = [wd + (m,) for wd in words for m in tensor_set]
        Ncand = K_prev * nT
        # Integrability on last slot only: for each pt, each 2-form component,
        # row over candidates: Σ_{old_word, ℓ} c · B[k, old_word] · wedge(last(old_word), ℓ)
        # i.e. constraint indexed by (old_word, pt, comp) acting on (k, ℓ):
        #   Σ_k Σ_ℓ c_{k,ℓ} · B[k, old_word] · W_pt(last(old_word), ℓ)[comp] = 0
        # This factorizes: for fixed old_word ow, last=a=ow[-1]:
        #   Σ_k B[k,ow] · ( Σ_ℓ c_{k,ℓ} · W_pt(a,ℓ)[comp] ) = 0  ∀ ow, pt, comp
        # Build constraint matrix M: rows = (ow, pt, comp), cols = (k, ℓ) flattened.
        last_of = np.array([wd[-1] for wd in words])
        # group old words by env = ow[:-1] (integrability acts within each env)
        env_groups = {}
        for ow_i, ow in enumerate(words):
            env_groups.setdefault(ow[:-1], []).append(ow_i)
        distinct_a = sorted(set(int(x) for x in last_of))
        ncomp = len(pairs)
        # precompute W[p][a] = (ncomp, nT) array of wedge components
        Wcache = []
        for p in range(npts_max):
            Wp = {}
            for a in distinct_a:
                arr = np.zeros((ncomp, nT), dtype=np.int64)
                for li, l in enumerate(tensor_set):
                    wv = wedge_at_pt(D[p], a, l, pairs)
                    for c in range(ncomp):
                        arr[c, li] = wv[c]
                Wp[a] = arr
            Wcache.append(Wp)

        def build_rows(p_lo, p_hi):
            rows = []
            for p in range(p_lo, p_hi):
                Wp = Wcache[p]
                for env, ow_idx_list in env_groups.items():
                    acc = np.zeros((ncomp, K_prev, nT), dtype=np.int64)
                    for ow_i in ow_idx_list:
                        a = int(last_of[ow_i])
                        acc = (acc + np.einsum('k,cn->ckn', B[:, ow_i],
                                               Wp[a])) % PRIME
                    for c in range(ncomp):
                        if acc[c].any():
                            rows.append(acc[c].reshape(-1))
            return rows

        # adaptive: build rows incrementally, append, recompute nullspace
        npts_seq = [6]
        while npts_seq[-1] < npts_max:
            npts_seq.append(min(2 * npts_seq[-1], npts_max))
        rk_prev = -1; null = None; all_rows = []
        p_done = 0
        for np_ in npts_seq:
            all_rows.extend(build_rows(p_done, np_))
            p_done = np_
            M = np.array(all_rows, dtype=np.int64) if all_rows \
                else np.zeros((0, Ncand), dtype=np.int64)
            null = nullspace_modp(M)
            rk = Ncand - null.shape[0]
            if rk == rk_prev:
                conv[w] = (True, np_)
                break
            rk_prev = rk
        else:
            conv[w] = (False, npts_seq[-1])

        # null is (K_new, K_prev*nT) — coefficients c_{k_new}[k_old, ℓ]
        K_new = null.shape[0]
        null3 = null.reshape(K_new, K_prev, nT)
        new_B = np.einsum('abc,bd->adc', null3, B).reshape(K_new, len(words)*nT) % PRIME
        if w == wmax:
            # tensored with LE only ⇒ K_new IS K_w^{LE}
            Kw[w] = K_new
            Kw_full[w] = f"(skipped at wmax; Ncand={Ncand})"
        else:
            le_mask = np.array([wd[-1] in LEset for wd in new_words])
            nonLE_cols = np.where(~le_mask)[0]
            if len(nonLE_cols):
                ech, piv = _row_echelon(new_B[:, nonLE_cols].copy())
                Kw[w] = K_new - len(piv)
            else:
                Kw[w] = K_new
            Kw_full[w] = K_new
            # Compress for next iteration: keep words with non-zero column
            nz_cols = np.where(new_B.any(axis=0))[0]
            words = [new_words[j] for j in nz_cols]
            B = new_B[:, nz_cols]

    return {"Kw_LE": Kw, "Kw_full_noLE": Kw_full, "converged": conv,
            "setup_s": round(setup_s, 2), "n_letters": n}


def per_orbit_report(alphabet_json, *, wmax=4, npts=96):
    """`--per-orbit` driver: per-orbit residual K_w via recursive_Kw, with a
    wallclock stamp.  `alphabet_json` is a path/dict/Alphabet accepted by
    load_alphabet (collapse schema or landau-alphabet output)."""
    t0 = time.time()
    r = recursive_Kw(alphabet_json, wmax=wmax, npts_max=npts)
    r["wallclock_s"] = round(time.time() - t0, 2)
    return r


# =========================================================================
#  CLI
# =========================================================================
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--alphabet", required=True,
                    help="alphabet JSON (collapse-schema or landau_alphabet.py output)")
    ap.add_argument("--connection", default=None,
                    help="optional connection JSON {masters, A:{letter:[[i,j,v],..]}}")
    ap.add_argument("--wmax", type=int, default=4)
    ap.add_argument("--method", choices=["auto", "exact_Q", "modp"], default="auto")
    ap.add_argument("--no-steinmann", action="store_true",
                    help="disable the Steinmann adjacent-channel cut even when "
                         "the alphabet supplies channel tags")
    ap.add_argument("--timeout", type=int, default=120,
                    help="per-weight integrability timeout (s)")
    ap.add_argument("--out", default=None)
    ap.add_argument("--compare", default=None,
                    help="second alphabet JSON for side-by-side")
    ap.add_argument("--per-orbit", action="store_true",
                    help="per-orbit residual K_w via the recursive integrable-symbol "
                         "builder; reports "
                         "Kw_LE / Kw_full_noLE / per-orbit mod-p convergence instead "
                         "of the global collapse table")
    ap.add_argument("--npts", type=int, default=96,
                    help="[--per-orbit] max mod-p sample points (adaptive 6→npts)")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args(argv)

    if args.per_orbit:
        rep = per_orbit_report(args.alphabet, wmax=args.wmax, npts=args.npts)
        print(json.dumps(rep, indent=2))
        if args.out:
            with open(args.out, "w") as fh:
                json.dump(rep, fh, indent=2)
            print(f"\nwrote → {args.out}")
        return 0

    rep = predict_collapse(args.alphabet, connection=args.connection,
                           weight_max=args.wmax, method=args.method,
                           steinmann=not args.no_steinmann,
                           timeout_s=args.timeout, verbose=args.verbose)
    print(format_table(rep))
    if args.compare:
        rep2 = predict_collapse(args.compare, connection=None,
                                weight_max=args.wmax, method="auto",
                                steinmann=not args.no_steinmann,
                                timeout_s=args.timeout, verbose=args.verbose)
        print(); print(format_table(rep2)); print()
        print(compare(rep, rep2))
    if args.out:
        with open(args.out, "w") as fh:
            json.dump(rep, fh, indent=2, default=str)
        print(f"\nwrote → {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

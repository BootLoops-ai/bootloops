"""trust.fraction_oracle — the pure-stdlib exact leg (design decision R3).

WHY IT EXISTS (SCOPE independence caveat): FLINT is a common SUBSTRATE of
kira128/FireFly and lp_syz's fmpq RREF (not shared house code, but shared
bytes under both). Check 2 (receipt core) is already FLINT-free; this module
makes the lp_syz EXACT route cross-gated by arithmetic that shares NOTHING
with FLINT: python stdlib `fractions.Fraction` only. No FLINT, no gmpy2, no
numpy, no sympy — `import fractions` and `hashlib` is the entire dependency
surface.

DESIGN POINT (R3, verbatim): "spot-check scale: a bounded row/point subsample
is the design point, it must be FOREIGN, not fast."
  - bounded: the oracle refuses blocks above max_rows x max_cols instead of
    silently grinding (CensoredError-style loud refusal);
  - FOREIGN: subsample selection is sha256-seeded from the system fingerprint
    + a caller tag — never "the first k rows", never the rows the main solve
    found convenient; verification points (d0,eta0) must be points the flint
    leg did NOT use for its own fit (caller declares; recorded in output).

Precedent: the winnow m2 battery cross-gates adapter values against an
independent Fraction-arithmetic oracle — this module is that pattern promoted
to a package leg.
"""
import hashlib
from fractions import Fraction


class OracleBoundExceeded(RuntimeError):
    """Block larger than the declared spot-check bound — refuse loudly,
    never grind (timed-pilot law: this leg is a spot check, not a solver)."""


def _seed_indices(n, k, seed_material):
    """Deterministic FOREIGN subsample: k of n indices, sha256-chain-seeded."""
    if k >= n:
        return list(range(n))
    out, seen = [], set()
    h = hashlib.sha256(seed_material.encode()).digest()
    while len(out) < k:
        for i in range(0, len(h) - 3, 4):
            v = int.from_bytes(h[i:i + 4], "big") % n
            if v not in seen:
                seen.add(v)
                out.append(v)
                if len(out) == k:
                    break
        h = hashlib.sha256(h).digest()
    return sorted(out)


def rref(rows, max_rows=600, max_cols=600):
    """Pure-Fraction row reduction. rows: list of dict col->Fraction.
    Returns (pivots: dict pivot_col -> reduced row dict, rank).
    Refuses blocks above the bound (OracleBoundExceeded)."""
    cols = sorted({c for r in rows for c in r})
    if len(rows) > max_rows or len(cols) > max_cols:
        raise OracleBoundExceeded(
            f"block {len(rows)}x{len(cols)} exceeds spot-check bound "
            f"{max_rows}x{max_cols} — this leg is a bounded oracle, not a solver")
    work = [dict(r) for r in rows]
    pivots = {}
    for col in cols:
        best = None
        for i, r in enumerate(work):
            if r.get(col):
                best = i
                break
        if best is None:
            continue
        prow = work.pop(best)
        pv = prow[col]
        prow = {c: v / pv for c, v in prow.items() if v}
        for c, prow_v in list(prow.items()):
            if prow_v == 0:
                del prow[c]
        # reduce existing pivot rows and remaining rows
        for pc, pr in pivots.items():
            f = pr.get(col)
            if f:
                for c, v in prow.items():
                    nv = pr.get(c, Fraction(0)) - f * v
                    if nv:
                        pr[c] = nv
                    elif c in pr:
                        del pr[c]
        for r in work:
            f = r.get(col)
            if f:
                for c, v in prow.items():
                    nv = r.get(c, Fraction(0)) - f * v
                    if nv:
                        r[c] = nv
                    elif c in r:
                        del r[c]
        pivots[col] = prow
    return pivots, len(pivots)


def reduce_column(rows, target_col, master_cols, max_rows=600, max_cols=600):
    """Express e_target over master_cols using the given exact relation rows
    (pure Fractions). Returns dict master_col -> Fraction, or None if the
    target column is not pivoted by the block. Bound-refusing."""
    pivots, _rank = rref(rows, max_rows=max_rows, max_cols=max_cols)
    if target_col not in pivots:
        return None
    masters = set(master_cols)

    memo = {}

    def expand(col, depth=0):
        if depth > 10000:
            raise RecursionError("expansion depth")
        if col in memo:
            return memo[col]
        if col in masters or col not in pivots:
            memo[col] = {col: Fraction(1)}
            return memo[col]
        out = {}
        row = pivots[col]  # col + sum tail = 0  =>  col = -tail
        for c, v in row.items():
            if c == col:
                continue
            for cc, vv in expand(c, depth + 1).items():
                nv = out.get(cc, Fraction(0)) - v * vv
                if nv:
                    out[cc] = nv
                elif cc in out:
                    del out[cc]
        memo[col] = out
        return out

    return expand(target_col)


def verify_witness_combination(rows_exact, lam_idx, lam_val_exact, target_col,
                               c_exact):
    """FLINT-free exact check of a lambda identity over Q:
        sum_i lam[i]*R_i == e_t - sum_m c[m]*e_m   (exact Fractions)
    rows_exact: list of dict col->Fraction (the caller's OWN parse);
    lam_val_exact: Fractions parallel to lam_idx. Returns (ok, residual_nnz)."""
    r = {}
    for i, li in zip(lam_idx, lam_val_exact):
        if not li:
            continue
        for c, v in rows_exact[i].items():
            nv = r.get(c, Fraction(0)) + li * v
            if nv:
                r[c] = nv
            elif c in r:
                del r[c]
    expect = {target_col: Fraction(1)}
    for m, cm in c_exact.items():
        if cm:
            expect[m] = expect.get(m, Fraction(0)) - cm
    expect = {c: v for c, v in expect.items() if v}
    return r == expect, len(set(r) ^ set(expect))


def spot_check_rows(rows, k, tag):
    """FOREIGN bounded row subsample: k rows chosen sha256-seeded from the
    row-set content + caller tag (never 'first k'). Returns (indices, rows)."""
    material = tag + "|" + hashlib.sha256(
        repr(sorted((i, sorted(r.items())) for i, r in enumerate(rows))).encode()
    ).hexdigest()
    idx = _seed_indices(len(rows), k, material)
    return idx, [rows[i] for i in idx]

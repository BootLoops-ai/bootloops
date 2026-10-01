#!/usr/bin/env python3
"""Deterministic synthetic fixtures for the rankscreen battery.

Every fixture is built in code, with a ground truth that follows from the
construction itself — no external data files. The construction is a planted
staircase:

  * choose pivot columns c_1 < ... < c_r;
  * generator row i puts its MINIMUM-column entry at c_i, has zero
    coefficient on every other pivot column, and carries seeded
    exact-Fraction entries on non-pivot columns greater than c_i;
  * because a generator has no support on any earlier pivot column,
    reduction against the earlier pivot rows leaves it unchanged, so an
    eliminator with the pinned semantics (row order, min-column pivot)
    pivots each generator at its own c_i: rank = r and the pivot sequence
    is exactly [(c_1, 0), ..., (c_r, r-1)];
  * a DEPENDENT row is an exact Fraction combination of generator rows with
    the right-hand side combined identically: it reduces to 0 = 0;
  * an INCONSISTENT row is the same combination with the right-hand side
    shifted by a nonzero constant: it reduces to 0 = shift — a planted
    0 = 1 row, visible at every prime that maps the system faithfully.

Builders return (rows, labels, ncols, nunk, truth). rows is the tool's
input form [({col: Fraction}, Fraction rhs)]; truth is a dict holding the
planted rank / pivot_cols / pivot_rows / incon_idx / closure_candidate
(plus builder-specific extras documented on each function).
"""
import random
from fractions import Fraction


def _frac(rng, num_mag=10**6, den_mag=10**4):
    return Fraction(rng.randint(1, num_mag) * rng.choice((1, -1)),
                    rng.randint(1, den_mag))


def _staircase(rng, ncols, pivot_cols):
    """Generator rows for the planted staircase (see module docstring)."""
    pivset = set(pivot_cols)
    gens = []
    for c in pivot_cols:
        extras = [j for j in range(c + 1, ncols) if j not in pivset]
        take = rng.sample(extras, max(1, len(extras) // 2)) if extras else []
        row = {c: _frac(rng)}
        for j in take:
            row[j] = _frac(rng)
        gens.append((row, _frac(rng)))
    return gens


def _combo(rng, gens, picks):
    """Exact Fraction combination of generator rows (coeffs and rhs alike)."""
    row, rhs = {}, Fraction(0)
    for gi in picks:
        a = _frac(rng, num_mag=100, den_mag=10)
        grow, grhs = gens[gi]
        for c, v in grow.items():
            row[c] = row.get(c, Fraction(0)) + a * v
            if row[c] == 0:
                del row[c]
        rhs += a * grhs
    return row, rhs


def planted_control(seed=20260826, ncols=64, nunk=48, n_dependent=8):
    """Closure-candidate-shaped known truth: pivots on ALL nunk pure
    columns, trailing dependent rows, no inconsistency. truth extras: none."""
    rng = random.Random(seed)
    pivot_cols = list(range(nunk))
    gens = _staircase(rng, ncols, pivot_cols)
    rows = list(gens)
    for _ in range(n_dependent):
        picks = rng.sample(range(len(gens)), rng.randint(2, 5))
        rows.append(_combo(rng, gens, picks))
    labels = [f"row{i:03d}" for i in range(len(rows))]
    truth = {"rank": len(pivot_cols),
             "pivot_cols": pivot_cols,
             "pivot_rows": list(range(len(pivot_cols))),
             "incon_idx": [],
             "closure_candidate": True}
    return rows, labels, ncols, nunk, truth


def planted_inconsistent(seed=20260827, ncols=64, nunk=48, r=40,
                         n_dependent=6):
    """Rank-deficit staircase plus planted 0 = 1 rows: one literal empty row
    with rhs 1 and two masked combinations with shifted rhs. truth extras:
    incon_idx lists all three planted rows (in row order)."""
    rng = random.Random(seed)
    pivot_cols = sorted(rng.sample(range(nunk), r))
    gens = _staircase(rng, ncols, pivot_cols)
    rows = list(gens)
    for _ in range(n_dependent):
        picks = rng.sample(range(len(gens)), rng.randint(2, 5))
        rows.append(_combo(rng, gens, picks))
    incon_idx = []
    for shift in (Fraction(1), Fraction(-7, 3)):
        picks = rng.sample(range(len(gens)), rng.randint(2, 5))
        row, rhs = _combo(rng, gens, picks)
        incon_idx.append(len(rows))
        rows.append((row, rhs + shift))
    incon_idx.append(len(rows))
    rows.append(({}, Fraction(1)))          # the literal 0 = 1 row
    labels = [f"row{i:03d}" for i in range(len(rows))]
    truth = {"rank": r,
             "pivot_cols": pivot_cols,
             "pivot_rows": list(range(r)),
             "incon_idx": incon_idx,
             "closure_candidate": False}
    return rows, labels, ncols, nunk, truth


def planted_prime_disagreement(sensitive_prime, seed=20260828, ncols=40,
                               nunk=32, r=24):
    """Staircase plus one appended row {c_free: P} = P with P =
    sensitive_prime and c_free a non-pivot pure column. Over Q (and mod any
    prime not dividing P) the row is a genuine pivot: rank r+1. Mod P the
    row maps to 0 = 0 and vanishes: rank r. Verdicts therefore DIFFER
    between a prime panel containing P and one avoiding it — the correct
    screen outcome is escalation. truth extras: rank_generic, rank_at_p."""
    rng = random.Random(seed)
    pivot_cols = sorted(rng.sample(range(nunk), r))
    gens = _staircase(rng, ncols, pivot_cols)
    rows = list(gens)
    c_free = min(j for j in range(nunk) if j not in set(pivot_cols))
    P = Fraction(sensitive_prime)
    rows.append(({c_free: P}, P))
    labels = [f"row{i:03d}" for i in range(len(rows))]
    truth = {"rank_generic": r + 1, "rank_at_p": r,
             "sensitive_prime": sensitive_prime, "c_free": c_free}
    return rows, labels, ncols, nunk, truth


def planted_bad_denominator(bad_den, seed=20260829, ncols=40, nunk=32, r=24):
    """Staircase with one coefficient denominator equal to bad_den (a pool
    prime, or a product of panel primes). Any prime dividing bad_den cannot
    map the system to F_p at all — the classic false-witness trap. The
    correct behavior is detection (BadPrime) and receipted replacement,
    never a silent wrong verdict. Planted truth (at any usable prime) is
    the staircase truth."""
    rng = random.Random(seed)
    pivot_cols = sorted(rng.sample(range(nunk), r))
    gens = _staircase(rng, ncols, pivot_cols)
    grow, grhs = gens[r // 2]
    poison_col = max(c for c in grow if c != pivot_cols[r // 2])
    grow[poison_col] = Fraction(3, bad_den)
    rows = list(gens)
    for _ in range(4):
        picks = rng.sample(range(len(gens)), rng.randint(2, 4))
        rows.append(_combo(rng, gens, picks))
    labels = [f"row{i:03d}" for i in range(len(rows))]
    truth = {"rank": r,
             "pivot_cols": pivot_cols,
             "pivot_rows": list(range(r)),
             "incon_idx": [],
             "closure_candidate": False,
             "bad_den": bad_den, "poison_col": poison_col}
    return rows, labels, ncols, nunk, truth

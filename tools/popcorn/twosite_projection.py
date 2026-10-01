#!/usr/bin/env python3
"""twosite_projection — exact hypergeometric projection of a two-site
diploid genotype table to a fixed subsample size.

WHAT IT COMPUTES
  Two biallelic sites typed in the same N diploid individuals are summarized
  by their 3x3 genotype table  cells[a][b] = #{individuals with alt-allele
  dosage a at site 1 and b at site 2},  a, b in {0, 1, 2}  (nine_cell builds
  it from the two dosage vectors). For a subsample size of m individuals
  (n = 2m chromosomes),

      project_grid(cells, N, m)[i][j]
        = P(alt count at site 1 = i and alt count at site 2 = j
            in a uniformly random subset of m of the N individuals)
        = E over ALL C(N, m) subsets of the indicator {counts = (i, j)},

  i, j = 0..n, returned as an (n+1) x (n+1) grid of fractions.Fraction that
  sums to exactly 1. The subset draws the nine cells jointly, i.e. the cell
  occupancies (k_ab) of the subsample are multivariate hypergeometric,

      P(k) = prod_ab C(cells[a][b], k_ab) / C(N, m),   sum_ab k_ab = m,

  and (i, j) = (sum_ab a k_ab, sum_ab b k_ab). The grid is accumulated by an
  integer dynamic program over the nine cells (state = individuals chosen so
  far, running i, running j; weights are products of binomial coefficients),
  with ONE exact division by C(N, m) at the end. This is the two-site
  analogue of the hypergeometric projection of a single-site frequency
  spectrum to a smaller sample (Marth et al. 2004; as used by dadi,
  Gutenkunst et al. 2009), taken over individuals rather than chromosomes so
  that both sites are projected onto the SAME subsample; a genotype (dosage)
  route, no phase needed. Summing the grid over pairs of sites gives the
  expected two-site (joint) frequency spectrum of the panel at sample size n
  (Hudson 2001 for the two-locus sampling distribution it estimates).

  brute_grid(d1, d2, N, m)  the same grid by explicit enumeration of all
      C(N, m) subsets: the reference route, O(C(N, m)) per pair.
  seg_matrix(grid)          the unfolded segregating block i, j = 1..n-1
      of a derived-allele-oriented grid ((n-1) x (n-1)).
  flip_grid(grid, fr, fc)   reorient: reverse rows (site 1) and/or columns
      (site 2), i -> n-i, for sites whose ALT allele is ancestral.
  fold_matrix(grid)         per-site minor-allele fold of the segregating
      cells, i -> min(i, n-i), giving an (n/2) x (n/2) matrix; folding
      after projection equals projecting folded counts.
  add_into(acc, mat), add_scaled(acc, mat, k), zeros(d), cells_key(cells)
      exact accumulation helpers (acc += mat; acc += k*mat with k an int;
      Fraction arithmetic is exact, so k successive add_into == add_scaled).
  GridCache                 memo of project_grid / fold_matrix / oriented
      seg_matrix keyed on (cells_key, N, m): equal tables give identical
      rationals, so a panel scan computes each distinct table once.
  brute_check(dosages, N, n_list)  project_grid vs brute_grid on every pair
      of a small panel, cell by cell with exact ==, plus sum-to-one.

REGISTER: EXACT. Every grid/matrix entry is a fractions.Fraction; equality
tests against the enumeration route are exact `==`, no tolerance anywhere.
GridCache.report() hit rates are float bookkeeping, never results. Standard
library only.

COST: project_grid is polynomial (at most (m+1) x (n+1)^2 DP states per
cell, nine cells); N = 113, m = 4 takes about a millisecond. brute_grid is
C(N, m) subset visits per pair (C(113, 4) = 6.4 million: seconds for ONE
pair) and exists only as the checker; brute_check refuses panels above
30 individuals or 200 sites.

REFERENCES
  Marth, Czabarka, Murvai & Sherry (2004), Genetics 166:351-372
    (hypergeometric projection of the frequency spectrum to a smaller
    sample).
  Gutenkunst, Hernandez, Williamson & Bustamante (2009), PLoS Genet.
    5:e1000695 (dadi; projection as standard practice).
  Hudson (2001), "Two-locus sampling distributions and their application",
    Genetics 159:1805-1817.

USAGE
  from twosite_projection import nine_cell, project_grid, brute_grid
  cells = nine_cell(d1, d2)              # d1, d2: alt dosages (0/1/2) per individual
  G = project_grid(cells, len(d1), 3)    # n = 6 chromosomes; sum == 1 exactly
  assert G == brute_grid(d1, d2, len(d1), 3)
  U = seg_matrix(flip_grid(G, alt1_is_ancestral, alt2_is_ancestral))
  F = fold_matrix(G)
  python3 twosite_projection.py [--geno tiny_geno.tsv] [--n 4,6,8]   # brute check
"""
import argparse
import itertools
import math
import os
import sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
REFERENCE_DIR = os.path.join(HERE, "reference", "twosite")

__all__ = ["nine_cell", "project_grid", "brute_grid", "flip_grid",
           "seg_matrix", "fold_matrix", "add_into", "add_scaled", "zeros",
           "cells_key", "GridCache", "brute_check", "main"]


# ------------------------------------------------ exact projection operator
def nine_cell(d1, d2):
    """3x3 genotype table of a pair of sites from two alt-dosage vectors
    (entries 0/1/2, one per individual, same individuals in the same order):
    c[a][b] = number of individuals with dosage a at site 1 and b at site 2."""
    c = [[0] * 3 for _ in range(3)]
    for a, b in zip(d1, d2):
        c[a][b] += 1
    return c


def project_grid(cells, n_ind, m):
    """EXACT projection: expectation over ALL C(N,m) subsamples of m
    individuals (n=2m chromosomes). DP over the 9 genotype cells of the
    multivariate hypergeometric; integer weights, one exact Fraction
    division at the end. Returns (n+1)x(n+1) grid of Fractions:
    grid[i][j] = P(alt count at site1 = i, at site2 = j) in the subsample."""
    n = 2 * m
    dp = {(0, 0, 0): 1}  # (t chosen, i, j) -> integer weight
    for a in range(3):
        for b in range(3):
            c = cells[a][b]
            if c == 0:
                continue
            ndp = {}
            for (t, i, j), w in dp.items():
                for k in range(0, min(c, m - t) + 1):
                    key = (t + k, i + a * k, j + b * k)
                    if key[0] > m:
                        break
                    ndp[key] = ndp.get(key, 0) + w * math.comb(c, k)
            dp = ndp
    denom = math.comb(n_ind, m)
    grid = [[Fraction(0)] * (n + 1) for _ in range(n + 1)]
    for (t, i, j), w in dp.items():
        if t == m:
            grid[i][j] += Fraction(w, denom)
    return grid


def brute_grid(d1, d2, n_ind, m, combos=None):
    """Explicit enumeration of ALL C(N,m) subsamples (the reference route
    for project_grid; O(C(N,m)) per pair). d1, d2 are the dosage VECTORS,
    not the table; `combos` may pass a precomputed list of index tuples."""
    n = 2 * m
    grid = [[Fraction(0)] * (n + 1) for _ in range(n + 1)]
    denom = math.comb(n_ind, m)
    unit = Fraction(1, denom)
    if combos is None:
        combos = itertools.combinations(range(n_ind), m)
    for combo in combos:
        i = sum(d1[x] for x in combo)
        j = sum(d2[x] for x in combo)
        grid[i][j] += unit
    return grid


def flip_grid(grid, flip_rows, flip_cols):
    """Reorient a grid to derived-allele counts: reverse the row index
    (site 1) and/or the column index (site 2), i -> n - i. Returns a new
    outer list (rows are shared when nothing is flipped along them)."""
    n = len(grid) - 1
    g = grid
    if flip_rows:
        g = g[::-1]
    if flip_cols:
        g = [row[::-1] for row in g]
    return g if (flip_rows or flip_cols) else [row[:] for row in grid]


def seg_matrix(grid):
    """Unfolded segregating cells i,j = 1..n-1 (derived-oriented grid)."""
    n = len(grid) - 1
    return [[grid[i][j] for j in range(1, n)] for i in range(1, n)]


def fold_matrix(grid):
    """Per-site minor-allele fold of segregating cells -> nf x nf,
    nf = n//2; folding after projection == projecting folded counts."""
    n = len(grid) - 1
    nf = n // 2
    out = [[Fraction(0)] * nf for _ in range(nf)]
    for i in range(1, n):
        for j in range(1, n):
            out[min(i, n - i) - 1][min(j, n - j) - 1] += grid[i][j]
    return out


def add_into(acc, mat):
    """acc += mat, entrywise. Mutates its FIRST argument only."""
    for r, row in enumerate(mat):
        for c, v in enumerate(row):
            acc[r][c] += v


def add_scaled(acc, mat, mult):
    """acc += mult * mat with mult an exact int. Fraction * int and
    Fraction + Fraction are exact, so this equals mult successive
    add_into(acc, mat) calls identically; zero entries are skipped (adding
    0 is the identity). Mutates its FIRST argument only."""
    for r, row in enumerate(mat):
        arow = acc[r]
        for c, v in enumerate(row):
            if v:
                arow[c] += v * mult


def zeros(d):
    """d x d matrix of Fraction(0)."""
    return [[Fraction(0)] * d for _ in range(d)]


def cells_key(cells):
    """Canonical immutable key form of a 3x3 contingency table."""
    return (tuple(cells[0]), tuple(cells[1]), tuple(cells[2]))


# ---------------------------------------------------------------- memo cache
class GridCache:
    """Memo of the projection keyed on (cells_key(cells), N, m).

    project_grid is a pure function of (cells, N, m), so equal keys give THE
    SAME rationals; caching changes which work is done, never what value
    comes out. Each slot holds [grid, fold_matrix(grid) or None,
    {(flip_rows, flip_cols): seg_matrix(flip_grid(grid, ...))} or None],
    the derived matrices filled lazily. Slots are never evicted and never
    mutated after first fill; every returned object is SHARED and must be
    treated as read-only (pass it only as the SECOND argument of add_into /
    add_scaled). max_entries > 0 stops admitting new slots past the cap
    (values are still computed and returned); 0 = unbounded."""

    def __init__(self, max_entries=0):
        self.max_entries = int(max_entries)
        self._slots = {}
        self.stats = {"grid_lookups": 0, "grid_misses": 0,
                      "grid_not_admitted": 0, "fold_misses": 0,
                      "umat_lookups": 0, "umat_misses": 0}

    def entry(self, ckey, n_ind, m):
        """Slot for (table key, N, m); computes the exact grid on a miss.
        ckey is a cells_key tuple (it indexes like the table itself)."""
        key = (ckey, n_ind, m)
        self.stats["grid_lookups"] += 1
        ent = self._slots.get(key)
        if ent is None:
            self.stats["grid_misses"] += 1
            ent = [project_grid(ckey, n_ind, m), None, None]
            if self.max_entries == 0 or len(self._slots) < self.max_entries:
                self._slots[key] = ent
            else:
                self.stats["grid_not_admitted"] += 1
        return ent

    def grid(self, cells, n_ind, m):
        """Memoized project_grid(cells, n_ind, m) (shared, read-only)."""
        return self.entry(cells_key(cells), n_ind, m)[0]

    def fold(self, ent):
        """fold_matrix(ent grid), memoized in the same slot."""
        fmat = ent[1]
        if fmat is None:
            self.stats["fold_misses"] += 1
            fmat = ent[1] = fold_matrix(ent[0])
        return fmat

    def umat(self, ent, flip_rows, flip_cols):
        """seg_matrix(flip_grid(ent grid, fr, fc)), memoized in the slot."""
        self.stats["umat_lookups"] += 1
        um = ent[2]
        if um is None:
            um = ent[2] = {}
        fk = (flip_rows, flip_cols)
        umat = um.get(fk)
        if umat is None:
            self.stats["umat_misses"] += 1
            umat = um[fk] = seg_matrix(flip_grid(ent[0], flip_rows,
                                                 flip_cols))
        return umat

    def clear(self):
        self._slots.clear()
        for k in self.stats:
            self.stats[k] = 0

    def __len__(self):
        return len(self._slots)

    def report(self):
        """Counters + hit rates (floats; bookkeeping only)."""
        s = dict(self.stats)
        s["distinct_entries"] = len(self._slots)
        lk, ms = s["grid_lookups"], s["grid_misses"]
        s["grid_hit_rate"] = (lk - ms) / lk if lk else None
        ul, um = s["umat_lookups"], s["umat_misses"]
        s["umat_hit_rate"] = (ul - um) / ul if ul else None
        s["max_entries"] = self.max_entries
        return s


# --------------------------------------------------------------- checker
def brute_check(dosages, n_ind, n_list=(4, 6, 8), max_ind=30, max_sites=200):
    """project_grid vs brute_grid on every pair of a small panel.

    dosages  list of per-site alt-dosage tuples (each of length n_ind).
    n_list   projected chromosome counts n = 2m to check.
    For every site pair and every n, the DP grid and the explicit
    enumeration over ALL C(N, m) subsets are compared cell by cell with
    exact Fraction ==, and each DP grid is checked to sum to exactly 1.
    Returns {"n_individuals", "n_sites", "per_n": {"n4": {m, n_subsamples,
    n_pairs, cells_checked, mismatches, sum_to_one_failures, exact_equal},
    ...}, "EXACT_EQUAL_ALL": bool}. Refuses panels beyond max_ind
    individuals or max_sites sites (the enumeration is O(C(N, m)))."""
    S = len(dosages)
    if n_ind > max_ind or S > max_sites:
        raise ValueError(f"brute_check limits: <={max_ind} individuals, "
                         f"<={max_sites} sites (got N={n_ind}, S={S})")
    for d in dosages:
        if len(d) != n_ind:
            raise ValueError("inconsistent individual count across sites")
    n_list = sorted(int(x) for x in n_list)
    for n in n_list:
        if n % 2 or n < 2 or n // 2 > n_ind:
            raise ValueError(f"n must be even, >=2, with m=n/2 <= N={n_ind}; "
                             f"got {n}")
    all_equal = True
    per_n = {}
    n_pairs = S * (S - 1) // 2
    for n in n_list:
        m = n // 2
        combos = list(itertools.combinations(range(n_ind), m))
        cells_checked = 0
        mismatches = 0
        sum_one_fail = 0
        for i in range(S):
            for j in range(i + 1, S):
                d1, d2 = dosages[i], dosages[j]
                g_dp = project_grid(nine_cell(d1, d2), n_ind, m)
                g_br = brute_grid(d1, d2, n_ind, m, combos)
                if sum(sum(r) for r in g_dp) != 1:
                    sum_one_fail += 1
                for r in range(n + 1):
                    for c in range(n + 1):
                        cells_checked += 1
                        if g_dp[r][c] != g_br[r][c]:  # EXACT Fraction ==
                            mismatches += 1
        eq = (mismatches == 0 and sum_one_fail == 0)
        all_equal &= eq
        per_n[f"n{n}"] = {"m": m, "n_subsamples": len(combos),
                          "n_pairs": n_pairs, "cells_checked": cells_checked,
                          "mismatches": mismatches,
                          "sum_to_one_failures": sum_one_fail,
                          "exact_equal": eq}
    return {"n_individuals": n_ind, "n_sites": S, "per_n": per_n,
            "EXACT_EQUAL_ALL": bool(all_equal)}


# -------------------------------------------------------------------- CLI
def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Exact two-site hypergeometric projection: check the "
                    "DP operator against explicit enumeration of all "
                    "subsamples on a small genotype panel.")
    ap.add_argument("--geno", default=os.path.join(REFERENCE_DIR,
                                                   "tiny_geno.tsv"),
                    help="genotype TSV 'pos ref alt gts' or VCF(.gz) "
                         "(default: the shipped 12-individual, 24-site "
                         "fixture)")
    ap.add_argument("--format", choices=("auto", "vcf", "tsv"),
                    default="auto")
    ap.add_argument("--chrom", default=None, help="VCF chromosome filter")
    ap.add_argument("--n", default="4,6,8",
                    help="projected chromosome counts n=2m, comma-separated")
    a = ap.parse_args(argv)
    from twosite_spectrum import load_geno   # reader lives with the estimator
    sites, n_ind, _counts, _fmt = load_geno(a.geno, a.format, a.chrom)
    r = brute_check([s[3] for s in sites], n_ind,
                    [int(x) for x in a.n.split(",")])
    for k, v in r["per_n"].items():
        print(f"{k}: m={v['m']} subsamples={v['n_subsamples']} "
              f"pairs={v['n_pairs']} cells={v['cells_checked']} "
              f"mismatches={v['mismatches']} "
              f"sum_to_one_failures={v['sum_to_one_failures']}")
    print(f"N={r['n_individuals']} sites={r['n_sites']} "
          f"EXACT_EQUAL_ALL = {r['EXACT_EQUAL_ALL']}")
    return 0 if r["EXACT_EQUAL_ALL"] else 1


if __name__ == "__main__":
    sys.exit(main())

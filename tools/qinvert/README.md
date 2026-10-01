# qinvert — exact consistency certificates for published summary statistics over partially released entity-level tables

## The problem it solves

An agency publishes statistics computed on more rows than it releases
(e.g. entities terminated or withheld before publication, such as insurance
plan contracts).  qinvert certifies whether the published quantiles, outlier
fences or means are arithmetically consistent with the released rows plus k
unreleased entries, and returns the witness families; a family is a
consistency certificate, never an identification of the true unreleased rows
(families are non-unique by construction).  Given

1. the published statistic values, taken EXACT (they are exact rationals of
   the hidden input under the publisher's stated convention), and
2. the released partial data,

qinvert solves for the **minimal delta families** — additions with exact
value windows, and removals drawn from the released rows — that make the
recomputed statistics match the published values exactly.  Deltas are
reported as witness FAMILIES with parity constraints per quantile
definition; families can be **stacked** across many datasets that share the
suppressed entities (one added entity carries one value per dataset it
participates in), under a **hard no-regression gate**: every dataset whose
published targets already recompute exactly must remain exact.

qinvert generalizes the profile-solver and cross-measure delta-set
machinery validated in a CMS Star Ratings audit.
Everything is exact `fractions.Fraction` arithmetic: **no float participates
in any accept/reject decision anywhere in this package** (the publisher
conventions handled so far — the documented QNTLDEF quantile definitions
applied to decimal-string scores — have exact-rational semantics, so there
are zero float-semantic sites to document; if a future target's semantics
are floating-point, each such site must be added and documented here).

## Package

```
qinvert/
  targets.py   statistic definitions (value + order-stat sensitivity)
  solver.py    profile-enumeration engine (single dataset, many targets)
  stacker.py   cross-dataset stacking with shared entities + no-regression
  tests/       executable receipts (planted truth, negative controls,
               parity units, stacking e2e, origin-tree regression)
```

## Supported statistics (`targets.py`)

- `Quantile(p, definition)` — all **5 SAS QNTLDEF** definitions
  (`"sas1".."sas5"`; Base SAS 9.4 Procedures Guide: Statistical Procedures,
  "The UNIVARIATE Procedure", section *Calculating Percentiles*; PROC MEANS
  shares them, default QNTLDEF=5) and all **9 R types**
  (`"r1".."r9"`; Hyndman & Fan 1996, *Sample Quantiles in Statistical
  Packages*, The American Statistician 50(4) 361-365).  Documented
  equivalences SAS1=R4, SAS2=R3, SAS3=R1, SAS4=R6, SAS5=R2 are asserted by
  the test battery, along with R's published `quantile(1:10, .25, type=i)`
  values and a read-only cross-check against the audit's validated
  `qntldef5`.
- `TukeyFence(side, mult, definition, cap)` — outer/inner fences with the
  CMS-style cap convention (lo-fence maxed with cap, hi-fence minned).
- `OrderStat(j)` / `MinStat()` / `MaxStat()`, `RangeStat()`.
- `Mean()` and `LinearCombo([(coeff, stat), ...])`.
- Published targets are exact values or `Interval(lo, hi)` (capped /
  inequality targets).  A composite target is simply a list
  `[(statistic, published), ...]` over ONE dataset.

Every statistic exposes `value(xs_sorted)` (exact) and `support(n)` — the
1-based order-stat indices and rational weights it depends on at size n.
`value` is derived from `support`, so the solver's sensitivity structure and
the evaluator cannot disagree; the tests pin both to external references.

## The algebra

**Fence → quartile inversion.** With multiplier m (m=3 outer):
`lo = Q1 − m·IQR = (1+m)·Q1 − m·Q3`, `hi = Q3 + m·IQR = (1+m)·Q3 − m·Q1`.
When both fences are published and uncapped:

```
Q1 = ((1+m)·lo + m·hi) / (1+2m)        Q3 = (m·lo + (1+m)·hi) / (1+2m)
```

(m=3: `Q1 = (4lo+3hi)/7`, `Q3 = (3lo+4hi)/7` — the audit's validated
form).  A side whose published value EQUALS its cap
contributes only an inequality (`raw_lo ≤ cap`), leaving a one-parameter
family: the solver grids the free quartile over an index-shift window of
candidate values and solves the other from the surviving equation
(modes `both` / `lo_only` / `hi_only` / `none`, `invert_fence_pair`).

**Index-shift lemma.** A delta of k = k_add + k_rem entries changes the rank
of any fixed value by at most k.  Hence the order statistic at position j of
the new array is either an added value or an original value whose original
index lies in `[j − k_add, j + k_rem]`.  Proof: each addition below a value
raises its rank by exactly 1, each removal below lowers it by 1; additions
and removals elsewhere leave it unchanged.  All candidate enumeration is
windowed by this lemma (pad k+4), which is what keeps k ≈ 30 tractable.

**Parity / averaging constraints.** For p = a/b in lowest terms, `n·p` is an
integer iff `b | n`.  For QNTLDEF=5 (and R type 2) this is the averaging
branch: at p = 1/4, 3/4 the quartile is `(x_(np) + x_(np+1))/2` **iff
n ≡ 0 (mod 4)** — the "n mod 4" parity argument, generalized per
definition through `support(n')`: a required value v at a two-index support
`w_j·x'_(j) + w_(j+1)·x'_(j+1) = v` is realized either by **both** order
stats equal to v (added copies of v when v is scarce in the base) or by a
**straddle pair** (s, s') with `w_j·s + w_(j+1)·s' = v`, `s < v < s'`, and
no merged value strictly between — with exact counting constraints pinning
how many adds sit below s.  Interpolating definitions (R4–R9, SAS1, SAS4)
are the same machinery with unequal weights.

**Equal pinned values.**  Blocks pinned to the SAME value v (adjacent order
stats, non-adjacent equal published quantiles, Q1 = Q3 on tie-heavy data —
the CMS complaint-measure shape) are merged into one **span constraint**
`x'_(i) = v` over every covered index: order-statistic monotonicity forces
all indices between two equal-pinned blocks to v, and an averaging block in
such a run can only realize its both-equal branch.  Silently pruning these
scenarios is the failure mode pinned by the planted-truth receipts
(`tests/test_all.py::TestSameValueBlocks`).

**Dense targets.** An exact Mean is a sum constraint
`sum(adds) − sum(removals) = mean·n' − sum(base)`; the solver satisfies it
by sliding filler adds inside their windows on the register lattice (exact
integer feasibility), after the order-stat placement.

**Verification gate.** Every witness is materialized and EVERY target
re-evaluated from scratch (`targets.evaluate_targets`); anything inexact is
dropped.  The gate makes the search sound regardless of enumeration
shortcuts; completeness is bounded by the caps below.

## Quickstart

```python
from fractions import Fraction as F
from qinvert import (Quantile, TukeyFence, Mean, Interval,
                     solve, Dataset, stack)

# one dataset, published Q1/Q3 (SAS default definition) + capped fences
targets = [
    (TukeyFence("lo", 3, "sas5", cap=0), F(0)),        # published "0" = cap
    (TukeyFence("hi", 3, "sas5"),        F("1.41")),   # exact published
]
res = solve(values, targets, register=F(1, 100),  # data are 2dp
            value_lo=0, k_cap=8, rem_cap=2)
res["k"]          # minimal delta size, or None if unreachable within caps
res["families"]   # [{k_add, k_rem, removed, adds:[{value, window, role}]}]
res["caps_hit"]   # loud: every truncation leaves a note; a witness-free
                  # search always carries "k_cap-exhausted" (emptiness is
                  # certified only up to the echoed k_cap / rem_cap)

# cross-dataset stacking with shared suppressed entities
d1 = Dataset("D1", rows1, targets1, register=F(1, 100))   # rows: (eid, value)
d2 = Dataset("D2", rows2, targets2, register=F(1, 100))
out = stack([d1, d2], classes=("MA-PD", "MA-only"),
            participation=lambda cls, ds: not (cls == "MA-only" and ds == "D2"),
            couplings=[("D1", "D2")],   # entity in D2 carries equal value in D1
            k_cap=8, rem_cap=1)
out["entities"], out["removed_entities"], out["battery"], out["regressions"]
```

`window` on an add spec is `(lo, hi, lo_strict, hi_strict)`: any register
point inside keeps every order-stat target exact (tests re-verify an
alternate in-window point — the window claim is itself an executable
receipt).  `stack` raises `StackerError` rather than ever returning a delta
that regresses an already-exact dataset.

## Tests (executable receipts)

`python -m pytest tests/test_all.py -q` from this directory (or
`python tests/test_all.py`).  The battery:

- evaluator receipts against R's documented type-1..9 outputs, SAS hand
  computations, the SAS↔R equivalences, and the audit's `qntldef5`
  (read-only), plus a negative control that the 14 definitions do disagree
  somewhere (no aliasing);
- planted-truth recovery at ALL 14 definitions (plant an add, publish the
  perturbed quartiles, require the plant among the minimal families), and
  planted-removal recovery;
- negative controls: unreachable / inverted targets return empty;
- parity units (n mod 4 generalized per definition) and solver-level
  averaging: existing straddle, added straddle partner (pinned crit), and
  the two-copy both-equal branch;
- equal-pinned-value (span merge) planted truths: adjacent order stats,
  equal published quantiles on tie-heavy data, Q1 == Q3 with a k=1 witness,
  non-adjacent equal blocks, wavg merged with a point block — plus negative
  controls (non-monotone pins and unreachable spans stay empty);
- out-of-range OrderStat semantics: solvable at larger k without crashing,
  loud when k_cap is too small, not-exact under evaluation, safe in stack();
- loudness receipts: `k_cap-exhausted` on every witness-free search (absent
  on success and on exact-already), `register-none-grid-data-values-only`
  in register-None one-sided fence modes, and
  `touching-realizations-skipped` on the straddle-touching structural skip
  (with the planted k=1 truth the note is honest about missing);
- capped-fence (hi_only) recovery, Min/Range/Mean/Mean+median, a
  three-block Q1+median+Q3 assembly, interval-target probe path, the
  in-window swap receipt, a forced minimal k=30 family, and a spawn-pool
  equality check;
- stacking end-to-end: 3 datasets sharing 2 planted entities (one coupled
  super/sub pair) plus an untouched exact dataset; unreachable dataset
  reported unsolved; and a removal blocked by the no-regression gate;
- REGRESSION vs the audit (read-only): star-year 2024, D02/MA-PD,
  published fences (0, 1.41), n=534 — the recorded unique +1-addition family
  with witness 21/50 = 0.42 (receipt: the recorded
  `membership_families.json`, sha pinned),
  reproduced by qinvert from the audit's parsed CSVs and double-gated
  through the audit's own `tukey.tukey_fences`.

## Honest limits

- **Not for person-level data** or any table suppressed to protect the
  confidentiality of individuals (small-cell / statistical-disclosure-control
  suppression); the tool operates on organization-level releases and cannot
  identify any true suppressed record.
- **Witness families are non-unique.**  A family is a constructive
  certificate that the published values are consistent with SOME delta of
  that size and shape; it is never an identification of the true suppressed
  rows.  Distinct families (and whole windows of add values) are typically
  compatible with the same publications.
- **Completeness is bounded by the enumeration caps**, all loud (recorded in
  `caps_hit` when hit): search depth — every witness-free search records
  `k_cap-exhausted` (emptiness is only ever certified UP TO the echoed
  `k_cap` / `rem_cap`, never absolutely); `max_witnesses` (families returned
  per minimal k), `grid_cap` = 800 register multiples per candidate window,
  `scenario_cap` = 20000 pinned value scenarios per (k_add, k_rem),
  `assembly_budget` = 300000 iterations per add-solve, `rem_cands_cap` = 48
  removal candidate values, `rem_bases_cap` = 20000 removal multisets,
  `probe_budget` / `probe_k_max` for the interval-only fallback.  An empty
  result whose `caps_hit` contains MORE than `k_cap-exhausted` is NOT a
  proof of unreachability at the searched depths; an empty result with
  `caps_hit == ["k_cap-exhausted"]` is exhaustive for every k ≤ k_cap **up
  to the structural limits below**.
- **Structural limits**: one-sided (capped) fence modes grid the free
  quantile over index-shift windows on the data/register lattice — solutions
  whose free quantile falls outside every window sub-lattice are missed;
  interpolated-value sub-lattices are gridded only up to weight denominator
  8 (`grid-sublattice-*` cap note); with `register=None` the grid holds data
  values and their interpolations only (loud:
  `register-none-grid-data-values-only`).  Overlapping order-stat supports
  (two targets pinning the same index) are skipped with a cap note
  (`overlapping-supports-skipped`).  Blocks pinned to EQUAL values (Q1 = Q3,
  coincident published quantiles) are handled exactly via the span merge
  (see "Equal pinned values" above) — but a straddle end that merely
  TOUCHES a neighboring different pinned value is skipped loudly
  (`touching-realizations-skipped`).  An `OrderStat` whose index exceeds the
  array size at some k is an infeasible SIZE, skipped silently (exact, not a
  truncation), and evaluates as not-exact rather than raising.  Removal
  candidates are complete up to region-equivalence for pure order-stat
  targets (removing any value in a region shifts pinned indices
  identically); for dense targets (Mean) every distinct value matters and
  the candidate list is truncated loudly at `rem_cands_cap`.  `LinearCombo`
  targets drive the search only in the two-order-stat case (Range,
  midrange); otherwise they are verification-only.  Windows attached to
  filler adds guarantee ORDER-STAT targets only; when a Mean target is
  active the mean-adjustment pins filler values instead.
- **The stacker's assembly is greedy** (slot reuse, first-fit classes,
  first matching coupling sub-multiset): it certifies feasibility and
  no-regression of the emitted delta set, not global minimality of the
  entity count, and its unsolved list means "not found within caps", not
  impossibility.
- Published values must be EXACT rationals of the hidden input under the
  publisher's convention.  If the agency rounds for display, feed the
  display interval as an `Interval` target instead of an exact value.

## Provenance

| Convention | Source |
|---|---|
| SAS QNTLDEF 1-5 | Base SAS 9.4 Procedures Guide: Statistical Procedures, "The UNIVARIATE Procedure", *Calculating Percentiles* (QNTLDEF=/PCTLDEF=) |
| R types 1-9 | Hyndman & Fan (1996), Am. Stat. 50(4) 361-365; R `?quantile` |
| QNTLDEF=5 end-to-end | CMS Star Ratings Technical Notes 2024 final pp.149-150, validated vs published K-5/K-6 fence tables |
| Fence inversion, capped modes | planted-truth batteries (`tests/test_all.py::TestFenceAlgebra`) |
| Profile solver, realization counting, 800-grid cap | planted-truth batteries (`tests/test_all.py`, all 14 definitions) |
| Equal-pinned-value span merge (Q1 = Q3, coincident quantiles) | order-statistic monotonicity; planted-truth receipts `tests/test_all.py::TestSameValueBlocks` |
| Out-of-range OrderStat = infeasible size / not-exact (never a crash) | our convention; receipts `tests/test_all.py::TestOrderStatOutOfRange` |
| spawn-only multiprocessing | fork+BLAS deadlocks under load; spawn is mandatory (receipt `test_solve_many_spawn_matches_serial`) |
| Regression anchor | the recorded `membership_families.json` receipt, 2024 D02/MA-PD entry |

Full per-convention sourcing lives in `PROVENANCE.md`; the package `REGISTER`
dict in `__init__.py` states each module's validation register.

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).

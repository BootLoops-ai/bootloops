# qinvert — GUIDE

Tool page: https://bootloops.ai/tools/index.html (registry index)

KIND: package (stdlib-only; all-exact Fraction arithmetic)

## PURPOSE

Exact consistency certificates for published summary statistics over partially
released entity-level tables: given published
statistics (quantiles / Tukey fences / order stats / means, exact rationals under the
publisher's convention) plus released partial data, solves for minimal delta families
(adds with exact value windows, removals from released rows) that make recomputation
match exactly; stacks families across datasets sharing suppressed entities under a
hard no-regression gate. All-exact Fraction arithmetic — no float in any
accept/reject decision.

## USE-WHEN

- Testing whether published quantiles/fences/means recompute from the released rows,
  and certifying the minimal number of unreleased entries needed (validated on CMS
  Part C/D Star Ratings cut points, where the unreleased entities are terminated plan
  contracts, i.e. organizations).
- Cross-dataset entity stacking where one added entity must carry one value per
  participating dataset.

## NOT-FOR

Not for person-level data or any table suppressed to protect the confidentiality of
individuals (small-cell / statistical-disclosure-control suppression); the tool
operates on organization-level releases and cannot identify any true suppressed
record. Identification of the TRUE suppressed rows — families are constructive consistency
certificates only, non-unique. Emptiness is certified only up to echoed
k_cap/rem_cap (never absolute); empty result with caps_hit beyond "k_cap-exhausted"
proves nothing. Stacker assembly is greedy: certifies feasibility + no-regression,
not minimal entity count. Float-semantic pipelines undocumented (only exact-rational
replays covered so far). Rounded published values must enter as Interval targets,
not exact.

## INVOKE

`from qinvert import Quantile, TukeyFence, OrderStat, MinStat, MaxStat, RangeStat,
Mean, LinearCombo, Interval, solve, Dataset, stack`; `solve(values, targets,
register=F(1,100), value_lo=0, k_cap=8, rem_cap=2)`; `stack([Dataset(...), ...],
classes=..., participation=..., couplings=..., k_cap=, rem_cap=)`.

Quantile definitions: 5 SAS QNTLDEFs ("sas1".."sas5") + 9 R types ("r1".."r9").
Caps (all loud in caps_hit): k_cap/rem_cap, max_witnesses, grid_cap=800,
scenario_cap=20000, assembly_budget=300000, rem_cands_cap=48, rem_bases_cap=20000,
probe_budget/probe_k_max. Multiprocessing: spawn-only (fork+BLAS deadlocks).

## GATES

Built-in verification gate: every witness materialized and every target re-evaluated
from scratch; inexact dropped. Suite 52 tests + 14 subtests green from clean cwd;
audit regression anchor (unique +1 family, witness 21/50=0.42) double-gated
through the source audit's independent fence code (not shipped — those tests
SKIP on public copies).

## FOOTGUNS

- Structural misses: one-sided capped-fence modes grid on window sub-lattices
  (out-of-sublattice free quantiles missed; interpolated sub-lattices only to weight
  denom 8); register=None grids data values + interpolations only; overlapping
  order-stat supports skipped; straddle end touching a different pinned value
  skipped — all loud via cap notes.
- Equal-pinned blocks (Q1=Q3, coincident quantiles) handled via span merge — a
  strict-chaining solver silently prunes these (TestSameValueBlocks
  receipts pin the guard).
- Out-of-range OrderStat = infeasible size, evaluates not-exact, never raises (incl.
  inside stack).
- LinearCombo drives search only in the two-order-stat case (Range/midrange);
  otherwise verification-only. Filler-add windows guarantee order-stat targets only;
  active Mean pins filler values instead.

## What runs publicly (battery: PARTIAL — verified from this copy)

- `python3 tests/test_all.py` — 52 tests, all green from a scratch cwd. The
  audit cross-check tests SKIP when the source audit's code is absent
  (marked in the test file); everything else is self-contained.
- `PROVENANCE.md` / docstring receipt citations point at artifacts not shipped in this repo.

CREDIT: quantile definitions follow the Base SAS 9.4 Procedures Guide (QNTLDEF/PCTLDEF)
and R. J. Hyndman & Y. Fan (1996, Am. Stat. 50:361); the end-to-end anchor is the CMS
Star Ratings Technical Notes. See README 'Sources'.

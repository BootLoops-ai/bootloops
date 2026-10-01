# rankscreen — multi-prime parallel rank screen for exact-Q Gaussian eliminations

Compute rank + inconsistency detection mod k primes IN PARALLEL (one process
per prime), agree across all k = the verdict; pay the exact-Q Gaussian elimination ONLY for closure
candidates (rank-full + 0 inconsistent), seeded with the mod-p pivot
structure.  Verdict semantics pinned to a specific certified eliminator convention: row order, min-column
pivot, insertion-order reduction.

Files: `modp_rref.py` (dense-int64/sparse mod-p RREF backends),
`screen.py` (k-prime parallel screen + CLI + receipts),
`shim.py` (engine drop-in: serialize/load rows_all,
`gauss_rank_screened`, `exact_seeded`),
`battery.py` + `fixtures.py` (the in-tree synthetic battery: planted
known-truth control, planted 0=1 rows, engineered per-prime disagreement
that must escalate, p-divides-denominator detection/refusal — plus a
`--mutation-controls` mode that sabotages each tested mechanism on a
scratch copy and requires that leg to fail).

Battery (self-contained, no external data): `python3 selftest.py` — runs
`battery.py --mutation-controls` from a scratch directory; per-leg
PASS/FAIL, report JSON, one OVERALL line, exit 0 only if every leg passes
and every sabotage is caught.

Authority docs: manual `GUIDE.md` (this directory)
(invocation, gates, false-agreement bounds, footguns).  Gate receipts are
summarized in GUIDE.md; full receipts are not shipped.

Never retrofit a LIVE chain — adopt at the next natural run.

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).

# ratfit — GUIDE

Tool page: https://bootloops.ai/tools/ratfit.html

KIND: package (python-flint + fractions; no sympy in hot path). Layout: `import
ratfit` (package), submodules `ratfit.thiele_gate`, `ratfit.li2close`,
`ratfit.degree_alias` (pure standard library; its law receipt
`DEGREE_ALIAS_DISCRIMINATOR.json` ships beside it as package data);
`coincidence_loci`/`design_grid` (the pathgrid surface). Flat shims ship at the tools root
(`tools/thiele_gate.py` — forwards the CLI too, `tools/li2close.py`,
`tools/pathgrid.py`) so flat imports also work. Regression tests ship in
`tests/` here; one command runs them all (see the battery section below).

## PURPOSE

exact rational-function reconstruction and gating suite for sampled path-DE
transport — the fitting backend of vopclose and of any driver that samples a
connection along a path — plus the family-level
byte-exact node gate (thiele_gate), the weight-2 rational×log closed-form closer
(li2close), the degeneracy-safe grid designer (design_grid), and the failure-axis
triage for reconstructions that stop converging (degree_alias: grid alias vs
CRT height).

## USE-WHEN

- Reconstructing DE-matrix entries / ε-layer coefficients exactly from sampled
  points (blind Thiele, LOO-certified, known-Q Newton, mod-p screened).
- BEFORE farming: `degree_budget` (samples-needed analytics) and `design_grid`
  (degeneracy-safe 1-D grids).
- A per-node derivation/verification stage dominates farm cost → thiele_gate
  config-forcing pattern (~9×/node measured): fit the family CF-exact on verified
  nodes, byte-exact gate every config-forced node.
- A 1-d rational×log integral (wall-term kernel, corner constant, endpoint ε/Λ-
  regularized) needs a CLOSED FORM with the divergence web checked → li2close.
- Adding prime fields to a modular (CRT) rational reconstruction never lifts the
  failing entries, pass after pass, and you must choose between buying another
  prime and extending the interpolation grid → degree_alias (section below).

## NOT-FOR

numeric fitting or approximate regression (everything here is exact or fail-loud);
multivariate RATIONAL reconstruction beyond the bivariate patterns delegated by
numkin (`parametric_rec_interp` lifts polynomial param-dependence only);
li2close does NOT re-check interior W-root GROUP cancellation (see FOOTGUNS).

## INVOKE

core API — `thiele` / `thiele_loo` (blind + leave-one-out), `thiele_loo_screened`
(GF(p) prescreen, ~2× on blind-fail tails), `ratrecon_with_denom` (known-Q Newton),
`screen_modp` (screen-then-certify: no false rejects, certify stays exact),
`extract_stairs` (exact ratios R=Q_{m+1}/Q_m), `row_lcm`, `degree_budget`,
`coincidence_loci`/`design_grid`, `exact_entry_from_samples` (the vopclose
entry point), `parametric_rec_interp`/`mpoly_eval` (multivariate iterated-Newton
lift of parametric recurrence coefficients over a tensor-product param grid,
cross-checked at 2 off-grid param-tuples). `ratfit.thiele_gate`: `fit/validate/gate` API + `gate_family`
(parallel over components) + JSON CLI (also via the flat shim). `ratfit.li2close`:
atom-list → ReLi₂ + log-bilinear closed form + standalone pure-mpmath evaluator
emitter (no CLI).

INPUTS: exact rational sample points (Fraction); screen primes 2⁶¹−1 + a second
61-bit prime; known denominators when blind fitting saturates. li2close atoms:
`c·ln|P0+P1·t|/(W0+W1·t)` on (lo,hi), lo=0 ε-regularized, hi=∞ Λ-regularized;
coefficients rational or rational-in-one-parameter sympy exprs.

## GATES

16 core tests; degree_alias 11 (proof sides, negatives, receipt cross-check,
both CLI forms); thiele_gate 16/16 (positive control reproduces a reference family gate
EXACTLY — 88/88 entries, 10032/10032 new-node checks); li2close 16/16 incl. π²/6-
class known integrals (positive control re-derived reference wall constants
byte-identical at dps 70). Known-Q window law: `ratrecon_with_denom` reaches degree
n−8, NOT the blind ceiling (n−6)/2 — this decides farm sizes. Row-LCM family proven
31/31.

## FOOTGUNS

- Symptom split (the diagnostic law): LOO-fails pinned AT the degree ceiling =
  saturation (go known-Q or farm more samples); fails BELOW the ceiling, especially
  on previously-passing m=0 = corrupt node (→ design_grid / preflight /
  verify-consistency), NOT a fitting problem.
- DENOM-v2 pattern for any sampled-DE ε-layer reconstruction: blind-fit m0/m1 →
  build fresh per-entry denominator stairs Q_m ≈ Q_{m0}·R^Δ from the CURRENT run's
  qex (stale caches lack the needed denominators) + row-LCM family + ADAPTIVE degree
  cap (a hard cap silently zeroes deep-m fixes) + mod-p screen.
- Degenerate sample slices return a deterministic wrong-but-valid reduction,
  byte-identical on rerun, silent — ALWAYS design farm grids through `design_grid`;
  never hand-pick rational grid points (they can land on real poles of the
  connection: pre-compute the path-pole locus, use irrational-offset grids).
- thiele_gate saturation detection: the CF depth is an honest degree certificate
  verified on the CF's own support — but a zero-tail degeneracy can silently break
  interpolation there; CFDead → colliding-point eviction retry is built in, trust
  its verdicts not a manual refit.
- li2close raises `DivergenceError` if the lnε/lnΛ web does not cancel — an
  incomplete atom list is a WRONG ANSWER, not a regularization choice. Interior
  W-root GROUP cancellation is NOT re-checked: verify Σc·ln|P(t0)|=0 per shared-W
  group yourself (spurious-pole minefield class).
- `parametric_rec_interp` demands a COMPLETE duplicate-free tensor-product
  param grid and always cross-checks at 2 deterministic OFF-GRID param-tuples
  (`rec_at_point` is called there too — it must be evaluable off the grid). A
  cross-check ValueError means the grid under-resolves the parameter
  dependence (or a slice is degenerate): add axis nodes, never catch-and-trust
  the grid-only fit.
- The fit-point D check: a fit must refuse a D that vanishes at a
  FIT sample, not only at the held-out nodes -- otherwise that sample is absorbed
  (its row holds as 0/0; a gcd reduction then strips the root, so the returned
  D shows nothing).  `thiele_loo` / `thiele_loo_screened` / `thiele_gate.fit`
  verify every fit node with Q != 0 (no hole).  `ratrecon_with_denom` certified
  the held-out nodes only, so a node on a root of Qknown was absorbed (y*Q = 0
  there whatever y): it now returns ok=False when Qknown vanishes at ANY node.
  `screen_modp` may still pass that form -- the exact certify decides, as
  documented.  Battery: tests/test_ratfit_dcheck.py (the planted legs fail by
  name if the check is removed).

## Member: `ratfit.degree_alias` — grid alias vs CRT height

PURPOSE: when a rational reconstruction over prime fields keeps failing on the
same entries no matter how many primes are added, decide which axis is actually
short. A fit-degree probe on the FAILING entries separates the two causes:
numerator degrees pinned at the interpolation-grid ceiling (nd−1 / ne−1) with a
trivial denominator (0,0), identical across every probed prime, mean the grid is
too small (GRID ALIAS → extend the grid along one axis, never buy primes);
degrees stable below the ceiling with a nontrivial denominator mean the CRT
modulus is too small (HEIGHT → add prime passes, never touch the grid);
cross-prime disagreement is INCONCLUSIVE (probe more, never pick a fix).
Classifiers that only look at WHERE a failure sits (interior vs grid edge) are
blind to aliasing, because an alias is consistent across all primes and reads as
an interior failure.

INVOKE: library — `from ratfit.degree_alias import classify, classify_cell`;
`classify(grid=(nd, ne), cells=[per-cell list of per-prime ((num_d, num_e),
(den_d, den_e)) fit degrees])` → `{'verdict': ALIAS|HEIGHT|INCONCLUSIVE,
'axis': 'd'|'eta'|None, 'counts': {...}}`; `discriminator()` returns the parsed
law receipt (`DEGREE_ALIAS_DISCRIMINATOR.json`, package data beside the module).
CLI — `python3 -m ratfit.degree_alias --selftest` (tools/ on PYTHONPATH) or
`python3 tools/ratfit/degree_alias.py --selftest`. Inputs: the per-prime fit
degrees of ~5 failing entries and the grid point counts per axis; nothing else.

GATES: verdict ALIAS requires EVERY fitted cell alias-pinned; one below-ceiling
cell makes it HEIGHT, one disagreeing cell makes it INCONCLUSIVE (fail closed).
The selftest replays the receipt's two proof sides (side A: grid (50,100), fits
(49,99)/(0,0) identical at 7 primes → ALIAS, corroborated by 0/272 entries lifted
after two more primes; side B: grid (29,101), fits well below ceiling with
nontrivial denominators at 5 primes → HEIGHT, corroborated by all 1031 entries
lifted by one more prime), two negatives, and a cross-check that the shipped
receipt still carries those grids. Standard positive control before a live
verdict: run it on functions of the same batch that already closed, where the
number of primes each needed is known; a single mismatch voids the verdict.

LIMITS / FOOTGUNS: it decides the axis, not the extension mechanics; the probe is
on failing entries only (passing entries say nothing about the residual set); a
pinned numerator WITH a nontrivial denominator is not the alias signature (mixed
→ INCONCLUSIVE, do not extend the grid on it); the analogous but different law
for exact-sample fits lives in the core FOOTGUNS above (LOO-fails at the ceiling
= saturation, below it = corrupt node).

## What runs without the reference data (battery: PARTIAL)

- `PYTHONPATH=<repo>/tools python3 -m pytest tools/ratfit/tests -q` → 53 passed, 2
  skipped. The 2 skips are positive controls against reference data (not shipped; marked
  in the test files; they SKIP when that data is absent). Requires python-flint.


## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).

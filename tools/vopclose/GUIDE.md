# vopclose — GUIDE

Tool page: https://bootloops.ai/tools/vopclose.html

KIND: script (single file `vopclose.py`, used as CLI selftest + library). Its fitting
backend is the sibling vendored package `tools/ratfit/` (found automatically from this
layout). Companion passes (numeric node-table mirror, exactA repair) are outside the
scope of this package.

## PURPOSE

FULL-A function-level closure engine for a UT/graded 1-D path-DE, for when value-PSLQ
fails (mixed-weight layers) and a canonical ε-form is blocked or unavailable
(half-integer-orbit radical blocks).
Exact function model Σ C_tag·√rad·R(t)·G_word(t), closed under +, Rat·, √-mul, exact
d/dt and ∫₀ᵗ (Hermite split); layered VoP recursion with exactly-certified block
strategies.

## USE-WHEN

- Need function-level (not value-level) closure of a graded path-DE with mixed-weight
  layers.
- Radical (√-letter) blocks in the homogeneous flow: FULL/CONJ/TRI strategies incl.
  √quartic conjugation (linear→quartic letters; indicials (1/2,1/2,0) class).
- A log-solution finder on an augmented [[A,B],[0,0]] system is needed
  (`solve_rat_solutions` — also usable standalone).

## NOT-FOR

multi-variable closure (1-D path scope — one closure per kinematic line); LOO-saturated
entries (those need hybrid certified-numeric kernels — same epistemic status as a recorded
numeric result: documented, not exact).

## INVOKE

`python3 tools/vopclose/vopclose.py --test` (selftest 14/14, <2 s); as library: build a
`VopConfig` (ALL sector specifics live there: RatA entries, blocks, boundary values,
layer range, path poles, NC) and drive the layered closure;
`solve_rat_solutions(Ain, extra_rad=None, dN_extra=30, eig_fn=None)` for rational/log
solution hunting; rational fitting of connection entries delegated to
`ratfit.exact_entry_from_samples`.

INPUTS: exact rational connection entries (or samples → ratfit); block structure +
boundary values via VopConfig; path-pole list for the detour logic.
OUTPUTS: per-layer closures with zero-remainder symbolic DE certificate; Chebyshev
spectral-cumint evaluations (arb/acb); held-out endpoint gate results.

## GATES

zero-remainder symbolic DE certificate per closed layer is mandatory; held-out endpoint
gate harness built in. Validations: a production top sector, 15/15 TOP
coefficients (3 masters × ε⁻⁴..ε⁰) at 59.6–62.5 held-out digits; production
radical blocks (5 per-block layered closures).

## FOOTGUNS

- Det-basis apparent poles: zeros of det(fundamental basis) on the path are poles of
  C′=M⁻¹S that cancel only in Y=MC — integrating through them is a gate FAIL; the
  built-in detection + contour detour must stay on.
- Whole-triangle tag-mode closures grow ~12×/layer in container terms — split per
  target block (upward-closure subsystems) BEFORE committing deep layer ranges.
- Log-ring boundary PSLQ with >~20-prime baskets at ~110d inputs false-positives (dense
  height-1e4 fits) — sane-height rational closures only.
- `try_rational` zero-vector guard: exactly-real/zero inputs return Fraction(0)
  (mp.pslq refuses zero vectors) — in place; don't bypass it.

## Battery

- `python3 vopclose.py --test` → 14/14 PASS in under 2 s (any cwd), including exact DE certification 6/6 zero-remainder, held-out
  endpoint gate 6/6 ≥40d, and the ratfit-delegation leg against the sibling vendored
  `tools/ratfit/` package. Needs mpmath (+ python-flint for the ratfit known-Q route).

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).

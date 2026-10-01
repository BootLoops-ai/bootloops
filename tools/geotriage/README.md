# geotriage

One-command per-graph maximal-cut **geometry + value-fittability** auto-triage.
A single entry point that emits the closure route ∈ {A,B,C,D,E} with mandatory honesty flags.

## Usage

```bash
python geotriage.py graph_specs/icc.json --out icc_classify.json
```

```python
from geotriage import classify
rep = classify("graph_specs/sunrise_eq.json")
rep['route']            # 'B'
rep['elliptic']         # {j, kodaira, level_N, on_X1N, congruence_group, ...}
rep['honesty']          # {genus_method, cm_evidence, pf_order_certified, caveats}
```

## Input spec

Same JSON schema as `tools/landau-alphabet/graph_specs/` plus optional:
- `"maxcut": {"poly":..., "vars":[...], "modulus":..., "source":...}` —
  feed the known Baikov maximal-cut polynomial directly (proxy mode for
  graphs whose loop-by-loop localization is computed elsewhere — see the
  crossedbox3l and icc fixtures for worked examples).
- `"fiber_point": {...}` — kinematic specialization for the CM-j / a_p test.

## Decision tree → route

| variety | sub-case | route | ring |
|---|---|---|---|
| genus 0 | — | **A** | MZV |
| genus 1 | j ∈ CM-13 | **A** | Chowla–Selberg Γ |
| genus 1 | on X₁(N), N≤10,12 | **B** | Γ₁(N), L(χ_D,k) |
| genus 1 | non-congruence | **C** | — (DE-transport) |
| K3 | rank-2 (CM, Livné) | **B** | ℚ(√D), L(f3,s) |
| K3 | rank≥3 | **C/D** | Hilbert/Siegel frontier |
| CY3+ | — | **D** | CY transport |
| genus≥2 | — | **E** | Siegel (hard limit) |

## Honesty flags (always present)

- `genus_method`: `"exact"` (square-free degree / hyperelliptic R–H) |
  `"degree-bound"` (hypersurface degree only) | `"numeric-probe"`
- `cm_evidence`: list of `(p, a_p, χ_D(p))` tuples or `"j=... ∈ CM-13"`
- `pf_order_certified`: `False` unless the PF op was actually built (it isn't here)
- `caveats`: every assumption is listed

## Tests (6 known cases + unit/control legs)

```bash
python test_geotriage.py
```

| graph | type | genus | route | certified / proxy |
|---|---|---|---|---|
| box1l | polylog | 0 | A | proxy (LS=1/st) |
| sunrise_eq | elliptic | 1 | B | **derived** from Symanzik F (Γ₁(6), config (1,2,3,6)) |
| icc | elliptic | 1 | C | proxy (Källén quartic) |
| banana3l_eq | K3 | — | B | **derived** F + **certified** Livné a_p⇔χ_{−15} |
| crossedbox3l | polylog | 0 | A | proxy (residual quadratics) |
| banana3l_1112 | K3 | — | B | **derived** F + **probe** background-centred a_p⇔χ_{−3}, cross-checked against the Hecke trace 2(a²−3b²) |

Unit/control legs: wired a_p counter supersingular patterns (D=−4, D=−3, wrong-D
control), synthetic-truth K3 centring bank (planted background + level-15 CM a_p,
plus mutation / degenerate-fiber / no-lock controls). Member legs: exactj
known-answer + independent-oracle COPAIR cases (skipped by name without
`python-flint`); fiberstack a_p vs naive count and the a_{p²} identity with a
mutation control.

## Members

- **exactj** (`exactj.py`): exact minimal polynomial over ℚ of j(λ), λ = x₀²,
  for an algebraic fiber x₀ given by any irreducible integer minpoly — two
  independent flint-exact routes (literal resultant + multiplication-matrix
  minpoly) COPAIRed per call, failing loudly on disagreement. `j_minpoly_lc != 1`
  is an exact char-0 **non-CM certificate** (a CM j is an algebraic integer);
  never read lc = 1 as a CM claim. Requires `python-flint` (the battery skips
  these legs by name when it is absent).
- **fiberstack** (`fiberstack.py`): exact a_p engines for Legendre curves —
  `chi_table`/`legendre_ap_fast` (character tabled once per prime),
  `ap2_charsum` (a_{p²} over F_p² by norm-character sum; the identity
  a_{p²} = a_p² − 2p is a built-in COPAIR when λ ∈ F_p; O(p²) — short menus
  only), and `rational_heights_menu`. For fibers of field degree > 1: split
  primes read a_p through the roots of the minpoly mod p (λ = r²); inert primes
  read supersingularity via `ap2_charsum` (p | a_{p²}, p ≥ 5).

## What is heuristic vs certified

**Certified** (exact symbolic / theorem):
- genus from square-free degree of a 1-var Baikov polynomial (hyperelliptic R–H)
- j-invariant via binary-quartic I,J (validated on y²=x⁴−1 → 1728, y²=x³−1 → 0)
- Kodaira fiber config from j-pole orders; Beauville/Sebbar level match
- equal-mass-banana K3 CM via η-quotient a_p ⇔ χ_{−15}(p)=−1 (Livné)
- elliptic CM at a rational fiber: CM-13 table (complete for j ∈ ℚ) corroborated
  by the wired a_p point counter — evidence carries (p, a_p, χ_D(p)) triples and
  the twist-invariant supersingular pattern a_p = 0 ⇔ χ_D(p) = −1

**Probe** (`certified: false` in the report — a verdict, not a theorem):
- generic-K3 CM via the toric point count with exact background subtraction
  (`k3_ap_probe`): the algebraic background c₂p²+c₁p+c₀ is locked in exact
  arithmetic on the a_p = 0 primes (Deligne-bounded residuals, never a
  least-squares fit), then Livné-scanned over candidate χ_D.  A unique match at
  zero-density ≈ ½ ⇒ `is_cm: true` as a PROBE verdict with a mandatory caveat;
  no exact lock ⇒ honest `is_cm: null` (non-CM zeros are too sparse to centre,
  or the background is character-twisted — the probe refuses rather than guesses)

**Heuristic / degree-bound** (flagged in `caveats`):
- PF order = #residual-vars + 1 (banana-ladder rule; op NOT built)
- non-congruence verdict from "fiber config ∉ Beauville list ∧ deg(j) ∉ index table"
  (a Hauptmodul search on higher N is NOT run)
- `baikov_maxcut` via Symanzik-F is birational to the true Baikov maximal cut
  only for banana/vertex topologies; other graphs should supply `spec["maxcut"]`
  (and the K3 probe's toric model is a banana-shaped PROXY for non-banana K3s)

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).

# terrier

Unified string-vacua / F-theory toolkit: certified Calabi–Yau period
point-values, flux-lattice enumeration with Nikulin discriminant-form
certificates, and a certified flux-census pipeline, all built on a shared
verification chassis (receipts, controls, ball certificates, verdicts).

Tool page: https://bootloops.ai/tools/terrier.html

The design principle throughout is *certified, fail-closed computation*:
results are interval/ball values with explicit radii, every stage writes
sha-stamped receipts, mutation and planted-error controls must fail when they
should, and anything the suite cannot certify is refused loudly rather than
approximated silently.

## Package layout

```
terrier/
├── periods/        certified CY period wing
│   ├── pipeline/     pipeline modules and fixture data (cards/, conipfv/, dkmm/,
│   │                 ffp/, pfaffian/, sunit/) — sha-pinned; do not edit in place
│   └── summation/    direct-summation route for multi-parameter point-values
├── lattice/        flux-lattice wing (Julia + Python)
│   └── emit/         stratum emission engines + generic sweep chassis
├── census/         certified flux-census pipeline
│   ├── fold_hnf/  dedup/  dist_cert/  control_sampler/  isd_cutoff/  tile_engine/
├── common/         shared chassis library (both wings import it)
│   └── data/         packaged data (known-point orbit table)
├── selftest.py     suite-level regression gate (see "Verification")
├── build_manifest.py  rebuilds regression_manifest.json after a full green run
├── regression_manifest.json  sha256 pins for artifacts + registered test batteries
├── MODULE_MAP.md   per-module index: what each module does + its test battery
└── GUIDE.md        condensed use-when / not-for / invoke reference
```

Per-directory READMEs (`common/README.md`, `census/README.md`,
`lattice/emit/README.md`, `periods/summation/README.md`) document each area
in depth.

The tree ships a compact set of verification fixtures — toy and known-answer
data the batteries check real results against (synthetic control cards,
closed-form tables, small published samples). These are what the sha-pinned
selftest verifies — do not edit them; the manifest pins them and any edit
breaks verification. Full-scale reference data sets are **not** included in
the package: batteries that need them print a named SKIP naming the
environment variable that supplies the data (see "Environment"). The complete
replay of the flux-vacua computations of the paper *Exact methods for string
flux vacua*, with the full reference banks, is published at
https://www.bootloops.ai/files/kklt/terrier/.

## The three wings

### periods/ — certified CY period point-values

Certified period/transport chassis: card-based target specification,
Γ-series construction, exact flux contraction, stepped MUM + conifold
transport with a *proven* tail majorant, Krawczyk existence/uniqueness
certification of vacua, and two-dps verdicts with mutation controls.

| module | what it does |
|---|---|
| `cards.py` | Card contract for period targets, fail-closed PENDING-CARD semantics. The shipped fixture cards (`pipeline/cards/`) validate green, including `hv4-diag-L5.json`: an order-5 diagonal Picard–Fuchs operator derived uniquely from series data and matching the printed literature operator. `ads-5-81-3213.json` is derived from the KKLT vacuum 5-81-3213 of M. Demirtas, M. Kim, L. McAllister, J. Moritz and A. Rios-Tascon, JHEP 12 (2021) 136 [arXiv:2107.09064] (flux vectors, D3-brane count, W_0 normalization and polytope data from the paper and its arXiv ancillary set; journal version CC BY 4.0, adapted; no ancillary file is reproduced) — see THIRD_PARTY.md section B. |
| `geffseries.py` | Effective Γ-series facade (series + GV extraction, non-simplicial cones handled). Battery, with the two-modulus KKLT reference data supplied via `TERRIER_KKLT_BANK`: 901 w0 coefficients exact, 48 GV invariants exact. |
| `fluxcurves.py` | Exact flux contraction / frame layer; the ξ-sign convention is a code toggle (`xi_flipped`), not a comment. |
| `transport.py` | MUM-transport chassis with the Krawczyk vacuum layer (`transport.vac_layer()`). Battery, with `TERRIER_KKLT_BANK` set: 13/13 gates on the two-modulus KKLT example of Demirtas–Kim–McAllister–Moritz, result balls byte-match the manifest pins. |
| `envelope_certified.py` | Certified tail majorant for stepped transport, with a proved lemma (`PROVEN_MAJORANT.md`). Replaces an empirical ×4 envelope that was falsified twice — see the lemma's §5 and `envelope_comparison.md`. |
| `pfaffian.py` | Matrix (Pfaffian-system) transport: from system samples plus a GKZ ideal and curve to a certified period vector at the target. Front door is `route_choice`, which encodes a scalar-elimination refusal law (hand-eliminating to a scalar operator can blow up: one measured example produced order 13 with z-degree > 4052). CRT/Wang lift with fresh-prime and stability gates, exact annihilation certificate, closed-form tail balls. With the reference receipts supplied via `TERRIER_PFAFFIAN_BANK`, the battery checks the ads-5-81 result ball byte-exact at dps 60 and 90, the fresh-prime gate at 0 diffs over 147,968 values, and a mutation control that must be caught. See `PFAFFIAN.md`. |
| `conifold.py` | Conifold-chart transport (CP0–CP10 gate chart + local towers). Leg-1 scalar transport is fail-closed by design. Selftest 6/6 on the shipped fixtures: the DKMM routing control, the exact tower regression and caught branch-flip/branch-swap mutations. Honesty notes: `s_star`/C2-tail/clearance quantities are labeled estimates; verdicts pinned to float64 reference values are CONDITIONAL with a ~1e-13 relative floor. |
| `gvderive.py` | Derives and certifies Gopakumar–Vafa working sets from toric data (rather than trusting a repository table). No battery is registered for it in this release: the worked example it was validated on is not included in the package. GV-band labeling law: completeness is certified only through the margin level; beyond it the ×5 band is *labeled*, never folded into a certified radius. |
| `opderive.py` | Operator-derivation wrapper (the engine is `tools/annihilator/`): CRT/Wang staging plus an annihilation-proof gate. Battery, with `TERRIER_KKLT_BANK` set: θ-form matches the reference operator exactly; 838 series windows annihilated with 316 held out. |
| `ffp.py` | Finite-field fingerprint engines (all-z Frobenius scan for attractor hunting). Refuses p ≥ 2³¹ (`BigPrimeRefusal`) unless explicitly overridden — the int64 arithmetic overflows there. Validated envelope: HV p ≤ 127, BCM p ≤ 337. Selftest 12/12 on the shipped fixtures, including brute-force validation at tiny q and regeneration of the shipped reference pass-sets. |
| `sunit.py` | Refined S-integrality candidate generator (quadratic-field enumeration, distances-to-singular-values law). Truncation statements are emitted as machine-readable output: a bounded net, **not** a finiteness proof. Selftest 9/9 on a shipped synthetic control card whose acceptance values are derived by hand in the card itself, including planted-point recovery under a reduced net and out-of-box/singular must-fail controls. |
| `summation/` | Direct-summation route: evaluate Π and its jets as the multivariate Γ-series at the target point — no operator derivation, no path marching. The route of choice when scalar operators explode (a 6-parameter example has no practical scalar operator) and the target lies inside the proven convergence polydisc. Applicability is gated by an exact entropy radius law (s < 1, per-order ratio s²; if s² ≥ 0.9 the route is refused). Verified on 84 operators × 15,625 points of termwise identities; the battery replays these gates when the box-operator reference data is supplied via `TERRIER_BOXOPS_BANK` (a named SKIP otherwise). |
| `verdicts.py` | Verdict facade; delegates to `common/verdict.py`. |

Route-choice law: point-value targets go through the `pfaffian.py` matrix
route; derive a scalar operator only if a caller needs L_s itself, and price
it first via `route_choice` and a mod-p probe.

### lattice/ — flux-lattice enumeration and certificates

| module | what it does |
|---|---|
| `lattice.jl` + `validate.jl` | Exact SL(2,ℤ) reduction and class enumeration for definite binary quadratic forms. Fills two Hecke gaps: (a) definite-BQF `reduction`/`representatives` are `NotImplemented` for d < 0 in Hecke — this is an independent exact implementation, cross-checked against Hecke ideal class groups; (b) the proper-vs-improper trap — SL2 (proper) form classes and GL2 lattice classes are kept distinct, and imprimitive classes are real (e.g. T(23), T(47), T(71): 2, 3, 4 lattice classes vs form class numbers h = 3, 5, 7). Validation battery: 18/18 (needs `TERRIER_JULIA_PROJECT`). |
| `nikulin.jl` | Nikulin discriminant-form certificates. Turnkey `primitive_embeddings` at sig(3,0) hangs (> 25 min/class on a definite rank-19 complement genus); this module routes through the discriminant form at sig(0,3) (Milgram sign flip, 3 mod 8) with glue rigidity via `image_in_Oq` — 14.9 s/class. The certificate covers glue-level uniqueness only; residual multiplicity equals complement-class multiplicity. |
| `genus_enum.jl` | Aut-free Kneser genus enumeration with PARI `qfauto` mass checks (avoids a Hecke OOM). |
| `lp_kill_rank_agnostic.py` | Rank-agnostic Venkov LP-kill engine with exact `Fraction` Farkas certificates (exercised at ranks 3, 4, 6, 7). See `LP_KILL_RANK_GENERALIZATION.md`. |
| `emit/` | Stratum emission engines: `harness.jl` (emission core with a pluggable `PREDICATE_HOOKS` interface — swap in a different exact predicate chain), `s11.jl`, `predicates.jl`, and `sweepchassis.py`, a generic (entry stream, predicate chain) sweep driver with blind controls, planted-error and halt/resume semantics. Shipped battery: 3/3 — exact-arithmetic verification of the glue predicate chain on a 44-row known-answer table plus a hand-derived closed-form case and mutation must-fails. See `emit/README.md`. |
| `lattice_theta.py` | Jacobi theta nullwerte and lattice theta series (Z^n, D_n, D_n^+, E8, Leech), GKP flatness factors, and vectorized closest-vector routines. Standalone NumPy module; it does not inherit this wing's frozen K3 conventions (theta argument is the nome q, not tau). |

### census/ — certified flux census

Pipeline: fold → dedup (two independent stacks) → adjudicate → distance
certificates → controls → verdicts, with a tile route for *continuous*
exclusion over moduli boxes.

| stage | what it does |
|---|---|
| `fold_hnf/` | HNF + det-U complete SL(2,ℤ)-orbit key: fold the raw flux window before any counting, with witnessed merges and tamper must-fail tests. |
| `dedup/` | Two independent dedup engines plus a third-stack witness/reconcile layer. Law: any disagreement between the two engines is a RED result — file it with witnesses; the reconciled count is audit-only and never quoted as the answer. Battery 8/8 on shipped synthetic seeds. |
| `dist_cert/` | Exact `Fraction` balls, certified sqrt/ln, z/WP chart handling, exact root/CM predicates, and Fincke–Pohst exhaustive "OFF" certificates. With the planted-point reference data supplied via `TERRIER_PLANTS_DIR`, the battery checks 37/37. |
| `control_sampler/` | Exact-rational rejection sampling (dyadic k/2⁶⁴ coins, zero floats), exact CP brackets, Hoeffding p_U labels that can only widen. With the seed file supplied via `TERRIER_SEEDS_JSON`, the battery checks 19/19. |
| `isd_cutoff/` | τ-eliminated ISD finiteness cutoff on n₂ = f·f + h·h, with the rule that the formula is fixed *before* any count is computed under it. Derivation (`ISD_CUTOFF_DERIVATION.md`) and the two referee probes (`isd_ref.py`, `sign_probe.py`) included. |
| `tile_engine/` | Continuous-exclusion route over moduli tiles; fail-closed (a missing certificate means OPEN, never EXCLUDED). Requires the summation card (`periods/summation/`) and external reference data (`TERRIER_TILE_BANK`). |

A synthetic end-to-end smoke battery (`census/selftest_census_smoke.py`)
runs fold → two-stack dedup → dist_cert → receipt on a shipped 2×2 window.
`common/dedup.py` complements the census wing: known-point pullback-orbit
dedup with a newform backstop, run *before* any escalation compute — a
newform match is a rediscovery, never a new point. Selftest 5/5, including
a real pullback rediscovery case as a must-fire fixture.

### common/ — shared chassis

Small by design: only patterns both wings use.

| module | what it owns | selftest |
|---|---|---|
| `receipts.py` | Receipt schema, atomic tmp+fsync+rename writes, sha-streaming, code-sha stamps, receipt-as-checkpoint resume | 12/12 |
| `controls.py` | Hash-placed blind controls, mutation-must-fail harness, planted-error drills (plants are loudly synthetic) | 14/14 |
| `frames.py` | Convention registry template: sign/normalization conventions as pinned, tamper-evident rows with fail-closed lookup | 9/9 |
| `certs.py` | Two-dps agreement, ball-honesty radius reporting, published-value gates (last-digit ulp band), float64 trim guard | 22/22 |
| `verdict.py` | Asserting scoreboard/verdict primitives (the wings' `verdicts` facades delegate here) | green |
| `dedup.py` | Known-point orbit dedup + newform backstop (see census section) | 5/5 |

## Running it

The Python wings are modules to import (each `periods/*.py` facade re-exports its engine's front door; the `selftest_*.py` scripts beside them show complete calls) plus a few command-line entry points; the lattice wing is Julia scripts:

```sh
# periods wing (Python 3)
python3 -c 'import sys; sys.path.insert(0, "periods"); import cards, transport, pfaffian'
                                    # facades: cards, geffseries, fluxcurves, transport,
                                    # verdicts, opderive, pfaffian, conifold,
                                    # gvderive, ffp, sunit
# common chassis
python3 -c 'import sys; sys.path.insert(0, "common"); import receipts, certs'

# lattice wing (Julia 1.10 with an Oscar/Hecke project environment)
julia --project=$TERRIER_JULIA_PROJECT lattice/lattice.jl
julia --project=$TERRIER_JULIA_PROJECT lattice/nikulin.jl
```

Inputs: toric/GKZ data and target cards (fail-closed PENDING-CARD schema),
flux vectors, lattice Gram data, census entry streams. Outputs: certified
balls (two-dps convention), sha-stamped atomic receipts, exact annihilation
certificates, exact Farkas certificates, census counts with witness files,
and machine-readable truncation statements.

### Environment

The package reads one family of environment variables, `TERRIER_<MEANING>`.
None has a machine-local default: an unset variable that a leg needs produces
a named SKIP (in `selftest.py`) or a loud refusal (in the module), never a
silent pass.

| variable | read by | meaning |
|---|---|---|
| `TERRIER_JULIA_PROJECT` | `selftest.py`, all Julia legs | path to a Julia 1.10 project with Oscar/Hecke installed. Two Hecke operations the suite deliberately routes around: definite-BQF reduction (`NotImplemented` for d < 0) and sig(3,0) `primitive_embeddings` (hangs). |
| `TERRIER_TOOLS_ROOT` | `periods/pipeline/conipfv/coniop/routeA_coni.py` | the repository `tools/` directory (default: resolved relative to the package), used to import `tools/annihilator/`. |
| `TERRIER_KKLT_BANK` | `periods/pipeline/run_pipe.py`, `periods/selftest_{transport,envelope}.py` and the geffseries/fluxcurves/opderive/verdicts batteries | reference data for the two-modulus KKLT example (restricted operators, certified W0 and W0-at-vacuum receipts). |
| `TERRIER_PFAFFIAN_BANK` | `periods/selftest_pfaffian.py` | reference receipts for the Pfaffian transport battery (ads-5-81). |
| `TERRIER_REFERENCE_BANKS` | `periods/pipeline/pfaffian/ads581/*.py` | root of the reference banks the Pfaffian route-A scripts resolve their D-module inputs from (layout: `periods/pipeline/pfaffian/README.md`; same root layout as `TERRIER_PFAFFIAN_BANK`). |
| `TERRIER_SKIP_JETS` | `periods/pipeline/conipfv/coniop/gseries.py` | if set, skip the jet computation in the conifold Γ-series. |
| `TERRIER_HV4_FRAME_DIR` | `periods/pipeline/cards/hv4-diag-L5_bank/frame_rank5.py` | frame matrices from which the rank-5 S6-invariant frame is derived. |
| `TERRIER_SUNIT_ATLAS`, `TERRIER_SUNIT_OUT`, `TERRIER_KNOWN_ORBIT` | `periods/sunit.py`, `periods/pipeline/sunit/*.py`, `common/dedup.py` | operator atlas directory, candidate output directory, known-orbit table (defaults: the packaged `atlas/`, `candidates2/`, `common/data/KNOWN_ORBIT.jsonl`). |
| `TERRIER_SHARD`, `TERRIER_NSHARD` | `periods/pipeline/sunit/gen_refined.py` | shard index / shard count for a split census run (defaults 0 / 1). |
| `TERRIER_BOXOPS_BANK` | `periods/summation/boxops_u1.py` | box-operator bank for the full-size summation gates. |
| `TERRIER_SWEEP_BANK` | `lattice/emit/s11.jl` | root holding the reference list of the rank-2 (U,U) cell (`<root>/glue/R2S11_UU_PILOT.jsonl`) that the s11 V1 gate reproduces byte-exactly. |
| `TERRIER_SWEEP_SHARD_RECEIPT` | `lattice/emit/selftest_sweepchassis.py` | a reference shard receipt for the byte-identical replay legs. |
| `TERRIER_TRANSVERSAL3` | `lattice/nikulin.jl` | path to the transversal-lattice JSON input (Gram matrices listed by determinant; not included) for the `validate03` and `cert:` modes. |
| `TERRIER_GENUS_OUT` | `lattice/genus_enum.jl` | output directory for genus-enumeration receipts (default: the module directory). |
| `TERRIER_GENUS_CONTROL_RECEIPT` | `lattice/selftest_genus_control_replay.jl` | reference receipt of the det-92 control genus. |
| `TERRIER_R18_DIR` | `lattice/selftest_lpkill_a1n_mustfail.py` | reference receipts for the LP-kill A1^6 control and the d4 genus replay. |
| `TERRIER_G4_DIR` | `census/fold_hnf/tests/test_fold_hnf.py` | reference fold data (G4 anchor set) for the fold_hnf acceptance legs. |
| `TERRIER_PLANTS_DIR` | `census/dist_cert/tests/test_dist_cert.py` | planted-point reference data for the dist_cert battery. |
| `TERRIER_SEEDS_JSON` | `census/control_sampler/tests/*.py` | seed file for the control_sampler battery. |
| `TERRIER_TILE_BANK` | `census/tile_engine/tile_engine.py` | reference bank root (u1s card + μ₂ pass module) for the tile engine. |

## Verification

`selftest.py` is the suite-level regression gate:

```sh
python3 selftest.py            # full run
python3 selftest.py --only cards,receipts,certs   # subset of batteries
```

It runs three sections, in order, and exits nonzero on any failure:

1. **M. manifest** — sha256 of every pinned file in
   `regression_manifest.json` (59 shipped fixtures). Every missing or
   changed file gets its own named FAIL row.
2. **B. batteries** — the 30 registered test batteries, each replayed in its
   own directory under resource caps (Python: `ulimit -v 32505856`, nice 5;
   Julia: single-threaded). Batteries run **only if the verified manifest
   sections passed** — a corrupted or missing pin refuses all replays rather
   than silently skipping anything. A battery whose inputs are a reference
   data set not included in the package (or a Julia/Oscar environment) is a
   named SKIP row saying exactly which env var to set; with the var set, it
   runs for real.
3. **X. coverage** — every `selftest_*.{py,jl}` on disk must be registered in
   the manifest, and the wing-regression map must point at registered
   batteries (no silently added or dropped tests).

**What this means on a checkout of this distribution:** `python3 selftest.py`
runs green (exit 0) as shipped, in well under a minute: the 59 shipped pins
verify, the 14 batteries whose inputs ship with the tree replay and PASS, and
every other leg prints a named SKIP row — the reference-data batteries name
the env var that would supply the missing data and the Julia legs name
`TERRIER_JULIA_PROJECT`. Nothing is skipped silently, and supplying wrong
data still fails loudly. A subset can be selected with `--only` (sections M
and X always run in full), e.g.:

```sh
python3 selftest.py --only cards,conifold,ffp,sunit,emit_L2
```

The sabotage behavior is by design: a single corrupted pin produces a named
`CHANGED artifact` failure, refuses all battery replays, and exits rc=1.

`build_manifest.py` rebuilds the pins; only do this after a full green run.

## Scope and caveats

- **Not for Feynman fixed-ε DE transport** — that is a different problem;
  use the DE-transport tool in this repository (`tools/wayfinder/`).
- `dist_cert`, `control_sampler`, and `tile_engine` carry criterion constants
  tuned to their original application; do not reuse them elsewhere without an
  independent validation receipt.
- `sunit.py` truncation output is a bounded net, not a finiteness proof (the
  output says so in machine-readable form).
- `ffp.py` refuses p ≥ 2³¹ for a reason (int64 overflow); do not bypass
  `allow_big_prime` without cause. Validated envelope: HV p ≤ 127,
  BCM p ≤ 337.
- `nikulin.jl` certifies glue-level uniqueness only.
- GV completeness is certified only through the margin level; the band beyond
  it is labeled, never a certified radius.
- `conifold.py` leg-1 scalar transport is fail-closed by design.
- Census dedup disagreements are RED results; the reconciled count is
  audit-only.
- The modules and fixture data under `periods/pipeline/` are sha-pinned; do
  not edit them in place — change them deliberately, re-run the corresponding
  battery, and rebuild the manifest (see `MODULE_MAP.md`).

## Related tools in this repository

- `tools/annihilator/` — the series→PF operator engine that
  `periods/opderive.py` wraps; also imported directly by
  `periods/pipeline/restrict_op.py` (`run_pipe.py g3`) and
  `periods/pipeline/conipfv/coniop/routeA_coni.py`, so those legs need the
  repository checkout, not a stand-alone copy of `tools/terrier`.
- `tools/pf_rank.jl` — mod-p holonomic rank probes.
- `tools/wayfinder/` — Feynman fixed-ε DE transport; not
  certified-PF, never force-fit one route into the other.
- `tools/pslq_gate.py` — integer-relation closure gates.

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).

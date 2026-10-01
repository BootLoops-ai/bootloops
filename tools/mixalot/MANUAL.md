# MIXALOT — exact Bayesian evidence for mixture models, in a box

**What.** A set of validated mixture-evidence engines, vendored byte-exact and
sha-pinned behind one front door. A user with a mixture question — *how many components? where are the
seams? is component X present?* — at large k, large g, or large N computes
with the original instruments instead of rebuilding them from the
paper's description.

**Where.** `tools/mixalot/` (vendored snapshot; the original source
trees remain the record of origin).
Receipts: `battery/battery.py` writes its run record to `battery/BATTERY.txt`
on every run; `mixalot/_pins.py` (sha pins).

## Quickstart
```python
import sys; sys.path.insert(0, "<your-checkout>/tools/mixalot")
import mixalot
mixalot.verify()                    # fail-closed integrity (run it first)
mixalot.Z([4,1,3,2], 2)             # exact evidence, Fraction
mixalot.gstar([4,1,3,2], gmax=3)    # exact component-count posterior
from mixalot.engines import zseries, seg_v1, frozen_comp_v1, w4_blind_gf
```
CLI: `python3 -m mixalot.cli counts.json [--gmax 4]`.

**Worked example — the regimes table (exact vs MC):**
`python3 examples/large_kgn_table.py` writes `examples/TABLE_LARGE_KGN.md`
live: large-N exact-vs-nested/ss with measured walls (exact logZ at N=1e5 in
~12 s; sampler biases printed per cell), large-k exact-only rows to k=1000 at
N=1e5 (~140 s — the regime no sampler reaches), and large-g + the exact
DP-limit -1/48 law. No sampler claim is made beyond the live cells.
Seeded MC cells are deterministic on a given stack (environment fence:
bit-reproducibility is stack-dependent; the battery's hard gate is
calibration |z|<4 against the exact truth, not bit equality).

## Capability map (engine -> what it computes -> battery receipt)
| need | engine (vendored) | battery leg |
|---|---|---|
| exact evidence, k-state counts, small/medium | `bigg` (front door `mixalot.Z/gstar`) | M3 |
| 1-var closed forms / large-N series | `closed_form_1var`, `zseries` | M3 |
| 2-way tables at LSX (Lin–Sturmfels–Xu 2009) scale (mod-p + CRT) | `zseries.Z_table_modp` (+ `lsx55_exact` driver) | M4c/M4d |
| LSX flagship direct expansion | `mixalot.core.lsx_example('swiss'\|'coin10'\|'coin242')` (runnable front door) | M4a/M4b |
| large k (thousands of states) | `bigk2.Z_fast2` (`bigk` dep) | M8 |
| large g / generic-g | `biggen`, `bigg` | M7b grid |
| exact Dirichlet-process limit | `w4_blind_gf.Zdpm_gf` (+ `w4_dpm_limit` gates) | M7 |
| float route to 1e6 components (k=2) | `swap_route` (numpy; certified vs exact at N=100) | M17 (selftest; fences still ride) |
| contiguous-process census (seams, g*) | `seg_v1` (exact + float DP), `seg_blind` twin | M5 |
| planted-seam power probes | `mic_power` (plant + detect harness) | M9 |
| blend / frozen-component evidence + presence BF | `frozen_comp_v1` (+ `_blind` twin) | M6 |
| holonomic closure helpers | `formula_emitter{,_dirichlet}`, `annihilator` | M18 (dirichlet emitter, incl. b0>A/b1>B polynomial branch) |
| Etienne (2005) neutral-biodiversity likelihood, exact + certified-ball | `etienne_evaluate` (vendored eco/ member) | M14 |
| exact box-moment Z(u) via certified contiguity walk (the fourth evidence route) | `mixalot.boxwalk` (subpackage `mixalot/boxwalk/`; CLI `python3 -m mixalot.boxwalk`) | M16 (its 14-leg selftest) |
| MC comparison suite (hm/bridge/ss/chib/nested; JAX HMC + Gibbs) | `estimators.run(family,k,counts,estimator,seed)` | M11 |
| Etienne-formula product-tree DP, exact + certified-ball (quasi-linear, overflow-impossible) | `ball_engine` (vendored eco/pilot/ member) | M19 (gate_engines at load) |
| certified ML by interval Newton (log coords, ball-ND Hessian) | `phase2_certify` (eco/pilot/) | — |
| equal-I multisample factorization (M = I^A·M~) | `multisample`, `multisample_ball`, `phase2_multisample_eqI` | — |
| Kronecker-packed 3-level hierarchical certified evaluator | `hier_kron` (+ `hierarchical3`, `hier_ball`, `hier_eval`) | — |
| independent Hoppe/CRP urn oracle (naive convolution) | `etienne_oracle` | M19 (gate_engines G1) |

### Member — Etienne pilot engine cluster (`vendor/eco/pilot/`, 16 files)
The engine cluster behind the single-file `etienne_evaluate` member: balanced
product-tree DP over fmpz_poly/arb_poly (quasi-linear, overflow-impossible);
certified ML by interval Newton on the exact gradient in log coords (running
digamma/trigamma; uniqueness + ball-ND Hessian = certified maximum); equal-I
multisample factorization (M = I^A·M~, parameter-free, built once); and a
Kronecker-packed 3-level certified evaluator.
- **Loadable engines:** ball_engine, etienne_oracle, phase2_certify,
  phase2_multisample_eqI, multisample, multisample_ball, hier_kron,
  hierarchical3, hier_ball, hier_eval, rf_engine, gate_engines
  (`gate_engines` EXECUTES its G1/G2 gate at load — exact == independent urn
  oracle + ball containment/width — ~0.2 s, print-only; battery leg M19).
- **Pinned but NOT loadable** (top-level scripts; run from a scratch cwd
  beside their input data): gate_rf.py (its `../data/volkov2005/` read is at
  module TOP LEVEL — it runs on import, not only under `__main__`),
  gate_worked_example.py, verify_hierarchical3.py, verify_recurrence.py.
- **cwd-write note:** `phase2_multisample_eqI.py` run as a script writes
  `PHASE2_MS_EQI_PANAMA.json` into the CURRENT directory and reads its
  multi-sample input by relative path (same law as lsx55_exact — script runs
  belong in a scratch dir; loading it as an engine is side-effect free).
- **Cross-ref clinch:** two interval-Newton certified-ML implementations
  ship — clinch (block-arrow Krawczyk + interval Schur, hierarchical/block
  axis) and this member's phase2_certify (scalar/low-dim log-coord interval
  Newton for the Etienne likelihood). Route hierarchical/block-scale
  certification to clinch.

### Member — boxwalk (`mixalot/boxwalk/`, the exact box-moment route)
Exact rational Z(u) = ∫_{[0,1]^n} ∏_k P_k(x)^{u_k} dx for integer
polynomials P_k and integer exponents u, by a mod-p contiguity walk on the
exponent lattice + CRT + rational reconstruction: Singular finds
boundary-safe relation modules, python-flint certifies them over Z, a
recorded pivot program refills the moving moment window (or the window is
raised by slab steps), per-prime replay closes. The fourth evidence route
beside the allocation sum, the lattice DP and the closed forms, with a
different scaling profile (window cells ~ side^n; n ≤ 6 in practice).
- **API:** `from mixalot import boxwalk`; `load_spec`, `plan`, `walk`,
  `produce` (gated: G1 crossover vs dense oracle, G2 mutation detected, G3
  two disjoint prime sets agree; writes `MANIFEST.json`),
  `emit_fiber_recurrence` (exact scalar contiguity recurrence along one
  exponent ray, fitted with the pinned `annihilator` engine, 2-prime-stable
  order, held-out-prime verified) + `fiber.gcrd/lclm/gcrd_reduce/...`
  shift-operator calculus, `verify_manifest` (independent replay: spec
  hash, residues, local CRT+ratrec, fresh-prime dense replay with a cost
  budget, gates, program hash).
- **CLI:** `python3 -m mixalot.boxwalk {selftest|plan|produce|fiber|verify}`
  from `tools/mixalot/`; `mixalot/boxwalk/cli.py` and `selftest.py` also
  run as plain scripts. Spec JSON: `variables`, `polynomials`
  {label: {"e1 .. en": coeff}}, `target_u`, `order`, `options`; example
  `mixalot/boxwalk/examples/jc_quartet.json`.
- **Env:** `BOXWALK_SELFTEST_OUT` (selftest record path; default beside the
  subpackage), `BOXWALK_ANNIHILATOR` (alternate annihilator.py),
  `PHYLO_GKZ` (kit dir for `examples/make_jc_quartet.py`).
- **Deps:** numpy, sympy, python-flint, Singular on PATH (selftest leg E
  skips by name without it).
- **Battery:** M16 runs the member's 14-leg selftest from a clean cwd
  (record routed to scratch); reference record `mixalot/boxwalk/SELFTEST.json`;
  full caveats in `mixalot/boxwalk/README.md`.

### MEMBER SECTION — authorship members: nullcal, seam, dcm

**Null calibration (`mixalot/nullcal.py`).** `calibrate(claim_bf, known_bfs)`
locates a mixed-vs-pure log10 Bayes factor in the distribution of the same
statistic over known-pure objects (leave-one-out wherever profiles are
trained on them) and returns counts, quantiles and a verdict
(INSIDE-NULL / TAIL / ABOVE-NULL). `blend_vs_pure_bf(U, gmix, Z=None)` and
`calibrate_objects(...)` compute the factors exactly for exchangeable count
vectors — inside the package via `engines.load("bigg")`, standalone via an
explicit `Z`. A mixture's extra freedom also absorbs object-to-object rate
variation the iid idealisation leaves out: a factor inside the known-pure
tail is not evidence of mixture. Reference fixture: the Federalist 65-paper
lists under four word lists (`fixtures/fixtures_nullcal.json`).

**Contiguous-block evidence (`mixalot/seam.py`).** `z_seam(seq, pA, pB)` —
exact evidence that profile A wrote a prefix and B the suffix, seam uniform
over the N+1 positions, with the seam posterior summary; `z_seam_brute` is
the enumeration reference. Named refusal `[seam.profile-normalization]` on a
non-normalised profile.

**Burstiness arm (`mixalot/dcm.py`).** Dirichlet-compound-multinomial
profiles `alpha = kappa * theta` with add-half `theta` over V listed + 1
residual categories and ONE shared `kappa` (fitted on known-pure objects,
then treated as an exact rational, `kn/1000`). `dcm_pure`, `dcm_mixture`
(token mixture, f ~ Beta(1,1): exact evidence, exact posterior-mean of f,
the three Bayes factors, float quantiles), `dcm_seam` (change-point, both
orientations pooled — the reference convention), `fit_kappa` (exact
comparison on a rational grid). Exact throughout; scaled-integer lattice of
the same class as the frozen-profile blend.

Coverage rows (battery):
| member | leg | instances | check |
|---|---|---|---|
| nullcal | S1 | 3 lists x 65 known BFs + No. 55 claims | counts + verdicts exact |
| seam | S2/S3 | pinned + 25 seeded; all 16 length-4 seqs | Fraction equality; sum == 1 |
| dcm | S5 | tiny instance, all 8 assignments | mixture + change-point == direct enumeration |
| dcm | acceptance leg A | reference instance MW30 (the 30 Mosteller–Wallace function words), V=31, N=2043 | 10 certified values to 1e-30, wall 167 s |
| all | S0/planted | fixture-sha gate; digit flip; bad profile | refusals BY NAME, rc 3 / rc 1 / exception |

## The applications, as recipes
- **Bible-class census** (the Isaiah instrument): per-unit count
  vectors -> `seg_v1.dp_float` over g with the `1/C(n-1,g-1)` correction ->
  posterior over g + boundary posteriors. Battery M5 reproduces the frozen
  Isaiah primary cell from public-text-derived counts: N=11637, g*=3,
  log10BF 95.998, exact-vs-float 3.6e-12. BOUNDARY CONVENTION: the recorded
  "ch39 boundary 0.989" is the g=2-CONDITIONED boundary posterior (condition
  on g=2, then locate the cut); the unconditioned g*=3 cut at ch39 is 1.000 —
  state the conditioning or the numbers will not match.
- **Federalist / Shakespeare-class presence questions**: frozen component
  profiles + `frozen_comp_v1.presence_bayes_factor` (the original
  Federalist / Shakespeare analyses used exactly this engine
  family; M6 reproduces the gate cases).
- **Patient-table-class problems** (the LSX 5.5 table, 3x3): 2-way tables via
  `Z_table_modp` + CRT; M4d reproduces the flagship's held-out-prime residue
  against the recorded 143/262-digit fraction (Lin, Sturmfels & Xu reported 16 days of Maple
  for this integral in 2009; this route: minutes per prime).
- **Planted-truth discipline** (plants before claims):
  `mic_power.sample_units` + `mixalot.core.p_g1` (the GUARDED front door —
  the raw engine's g=1..6 grid needs n>=6 units), and
  `mixalot.plant/recover_blind` for count-vector g-inference (its docstring
  states the MEASURED small-N power: 0/5, 0/5, 4/5 at the defaults — that is
  the weak-identifiability fence in action, not a defect).

**Worked examples (runnable as shipped, print-only, any cwd):**
`python3 examples/worked_examples.py` — large-N exact,
large-k with safe output handling (print lnZ — a
raw `repr()` of these rationals exceeds Python's 4300-digit int-str limit),
large-g + the exact DP −1/48 law via `mixalot.core.Zdpm` (front door coerces
alpha to Fraction — the raw engine returns float on float alpha), the
runnable LSX route, and one live exact-vs-MC cell.

## Input rules (the public front door)
- Counts must be INTEGERS: non-integer counts raise `NonIntegerCountError`
  — MIXALOT never truncates. Priors must be strictly positive.
- `measure=` accepts exactly `"dirichlet"` (default) or `"lebesgue"`.
- Front-door `Z/gstar` carry a cost guard at N=2,000,000 (k=1 answers are
  returned instantly); beyond it, use the scale engines per the capability
  map.
- CLI `counts.json` schema is printed by `--help` (list form, or object with
  required key `"U"`).
- `frozen_comp_v1`'s docstring cites a derivation note
  that is not distributed with this package; the battery's M6 gate
  rationals are the shipped evidence.

## Fences (ride every use)
- **Convention**: Dir(1) mixture evidence, the gated-identity convention
  (`lebesgue` switch on the front door). LSX rows carry their own per-example
  convention pins (bare integral vs full ML) inside the engine.
- **Weak identifiability of iid category counts**: for single-draw count
  vectors a g-mixture is marginally another categorical; g-inference is
  prior-geometry dominated at small N (measured: log BF(2v1) ~ +0.8 on
  1-component uniform-ish data, k=4 N=60). The exact posterior QUANTIFIES
  this honestly — where sampler noise manufactures certainty, the exact
  number refuses to (measured; recorded control shape-direction
  rows). Structure (units/segments as in seg_v1, frozen profiles
  as in blend) is what buys power; plant at your own shape before trusting a
  verdict.
- **Exact-frontier advisories** (`mixalot.advise`): quoted from the engine
  docstrings' measured practical frontiers — where the exact routes reach,
  and the instruction to validate any sampler against them (the reason the
  exact routes exist).
- **swap_route** is k=2, balanced-U, certified at N=100 against exact — the
  1e6-component float route is a demonstrated ROUTE, not a general-purpose
  evaluator (caveat recorded in the DP pin).
- **lsx55/Z_table_modp prime discipline**: p < 2^25 (int64 accumulate);
  ratrec needs ASYMMETRIC bounds (num << den for evidence fractions).
- **Engine bytes are never edited**; re-vendors restamp `_pins.py`
  in the same write. `mixalot.verify()` is fail-closed on the vendor set,
  advisory-loud on upstream-source drift — and the LOADER re-hashes source
  bytes per load and executes `exec(compile(source))`, never Python's
  bytecode cache (`__pycache__` is purged on every load; the pyc-shadow
  attack is a battery mutation control, M2b).
- **lsx55 cwd-write law**: `lsx55_exact.py` run as a script writes
  `LSX55_EXACT_RESULT.json` into the CURRENT directory — always run it from
  a scratch dir (receipt-clobber hazard).
  The battery's cheap row uses the single-prime residue check instead.
- **Checkpoint-writing hunt engines are EXCLUDED from the vendor set**
  (run_rayfarm2/ray_lift default-write into their source tree; gdred_n and
  rescan_16 are cwd-sensitive readers) — out of the tool's front-door
  scope by design; vendoring them would first require wrappers that
  parameterize their checkpoint/read dirs.

## Costs (measured)
coin10 sub-second; Isaiah cell ~20 s; blend/DP legs instant; bigk2 k=100
N=1e4 ~0.5 s; lsx55 one prime ~6.5 min + 3.6 GB (full 66-prime CRT ~7
core-hours; Lin, Sturmfels & Xu reported 16 days of Maple in 2009); bigg g=3 memory frontier N~3000
(~1 GB) .. N~1e4 (~25-30 GB, do not attempt under a 31 GB limit).

## References

- S. Lin, B. Sturmfels & Z. Xu, JMLR 10 (2009) 1611-1631, arXiv:0805.3602 ("LSX").
- Evans, Gilula & Guttman, Biometrika 76 (1989) 557-563 (the 3x3 table of LSX Ex. 5.5).
- Ferguson, Ann. Stat. 1 (1973) 209; Antoniak, Ann. Stat. 2 (1974) 1152.
- R. S. Etienne, Ecol. Lett. 8 (2005) 253-260 and 10 (2007) 608-618; Etienne & Haegeman, Theor. Ecol. 4 (2011) 87-109 ("EH2011").
- R. K. S. Hankin, J. Stat. Softw. 22(12) (2007) (the untb package; "SA5" = Etienne's sampling-formula example as shipped there).
- Volkov, Banavar, He, Hubbell & Maritan, Nature 438 (2005) 658-661.
- Hoppe, J. Math. Biol. 20 (1984) 91-94.
- J. Peccoud & B. Ycart, Theor. Popul. Biol. 48 (1995) 222-234; Larsson et al., Nature 565 (2019) 251-254 (txburst).
- Madsen, Kauchak & Elkan, ICML 2005, 545-552.
- F. Mosteller & D. Wallace, JASA 58 (1963) 275 ("MW30" = their 30 function words).
- Newton & Raftery, JRSS-B 56 (1994) 3; Meng & Wong, Stat. Sinica 6 (1996) 831-860; Xie, Lewis, Fan, Kuo & Chen, Syst. Biol. 60 (2011) 150; Chib, JASA 90 (1995) 1313; Neal (1999); Skilling, Bayesian Anal. 1 (2006) 833; Speagle, MNRAS 493 (2020) 3132 (dynesty).
- R. E. Moore, SIAM J. Numer. Anal. 14 (1977) 611; R. Krawczyk, Computing 4 (1969) 187; S. M. Rump, Acta Numer. 19 (2010) 287.
- F. Johansson, IEEE Trans. Comput. 66 (2017) 1281 (Arb).

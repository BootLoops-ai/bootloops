# tools/ — the BootLoops toolkit packages

One directory per package. Every package carries a `GUIDE.md` — purpose, when to
reach for it, its acceptance gates and known hazards, and a link to its page on
[bootloops.ai](https://www.bootloops.ai). The thematic roster (what belongs to what) is `../toolkit/ours/README.md`; the
working recipes are `../toolkit/ours/RECIPES.md`; external engines the packages drive are catalogued in
`../toolkit/external/TOOLS.md`. Our own engines (Eichler.jl, SOFIA.jl, Leviathan) ship in full
under `../upgrades/`; the patched forks of Kira, Blade and AMFlow.cpp that several packages drive
live in their own repositories, listed in [`../upgrades/ENGINES.md`](../upgrades/ENGINES.md)
together with the environment variables through which the packages find the built binaries.

**Scope.** These are the general-purpose packages. Code specific to an individual
problem (per-problem evaluators, results-table scripts, single-result engines) is not
in this repository — it ships with the problem pages on
[bootloops.ai](https://www.bootloops.ai), beside the results it produced.

**Time standard.** A package's default battery completes in ~20 s on a laptop; heavier tiers go behind flags.

**Two terms used throughout.** A *banked* value, fixture or receipt is one that passed its acceptance
battery and was then stored read-only as reference data; batteries compare fresh output against banked
objects. The *value of record* (form, figure, point of record) is the single banked object a package
treats as authoritative when several candidates exist; anything else is a candidate or a control.

**Flat files at this level** are compatibility shims and members that provide the
packages' flat import names (`mplll*.py` and `pslq_gate.py` →
lockpick; `thiele_gate.py`, `li2close.py`, `pathgrid.py` → ratfit;
`quad_probe.py`, `ell_subst.py`, `oracle_registry.py` → gatekeeper;
`numkin_harvest.py`, `numkin_sweep.sh`,
`shift_opt.py` → numkin; `parallel_kira_gen*.py` → kira-stack; `maxcut_*.py`,
`deriv_m2.py`, `sparse_cascade.py`, `pf_rank.jl` → maxcut/dipstick).
`fixtures/` is the shared fixture home.

## Verification classes

Every package ships behind an acceptance battery. What runs from a cold
clone differs by package:

- **selftest** — the battery or selftest runs green from the clone as shipped.
- **partial** — the public legs run green from the clone; legs that read internal
  reference data, or need an external engine you must build first, are marked in the
  GUIDE and skip or fail closed with a named error.
- **smoke** — no battery of record; smoke checks (worked examples, refusal paths,
  positive/negative controls) verified.
- **data-gated** — the package operates on data you supply (or a reference-data
  root that does not ship with the clone); without it, it refuses loudly by design.

## Package index

| package | what it is | battery |
|---|---|---|
| `abacus/` | certified point counts, Frobenius charpolys, and L-factors for abelian fourfolds from a period lattice | selftest |
| `amflow-kit/` | the house kit around the AMFlow.cpp fork: three gates that validate `amflow_cli solve_integrals` output (lint, smoke, A/B compare), the offline IBP-cache key predictor hashing the canonical jobs.yaml form the fork hashes (`amflow_kit.keypred`), and the shared launch memory fence with post-spawn `/proc/<pid>/environ` readback and abort by pid (`amflow_kit.memfence`) | partial |
| `annihilator/` | Picard–Fuchs / holonomic operator from an exact series alone (minimal P-finite recurrence + θ-form ODE) | selftest |
| `ansatzer/` | Ansatzer — predicts residual symbol-ansatz dimension per weight after structural cuts, before oracle-farming | selftest |
| `baller/` | arbitrary-precision ball arithmetic for science — the certified instruments under one front door; its hygiene wing owns the dps_lint mpmath import-time-dps AST linter (`python3 tools/baller/baller/dps_lint.py FILE_OR_DIR ...` / `python3 -m baller.hygiene lint ...`, battery leg L21) and the purity-scan literal scanner | partial |
| `blade/` | Wolfram-free Blade block-triangular IBP pipeline (BT search, per-point probes, symbolic tables) | partial |
| `clinch/` | certified local interval-Newton convergence certificates for fitted hierarchy optima | partial |
| `coalescer/` | connection coefficients of a Fuchsian operator onto its fractional-power (finite-monodromy) Frobenius branches by spectral projection of the local monodromy; generic door --op FILE.json (your operator + exact seed; template examples/gauss_2f1.json), banana threshold = the worked example (--masses) | selftest |
| `cosmoflow/` | FRW/dS wavefunction and correlator integrands and letter alphabets for an arbitrary site graph (chains, rings, stars, multi-loop polygons) via `alphabet_graph(nv, edges)` / CLI kind `graph`, plus the Cayley-Menger Baikov polynomial of any 1-loop n-gon; the 1-loop triangle numerics (oracle, maxcut periods, two-site chain) ship as the worked example under `examples/triangle/` | selftest |
| `counterweight/` | Counterweight — no-Mathematica ε-factorization of Kira-built connections (with `canonical_form/` member) | selftest |
| `dipstick/` | pre-compute triage: PF order / GKZ rank / count / cycle-type probes before committing a farm, plus the expansion-by-regions completeness certifier (`regions` member: Newton-polytope facet census, two independent checks; julia + Oscar, skips by name without them); carries `run_msolve_capped.sh`, the mandatory address-space-cap wrapper for every msolve exec | selftest |
| `dogtag/` | Dogtag — integral-family identity checks before compute is spent on a family file: isomorphism to a catalog family up to loop-momentum relabeling, planar or crossed in the physical leg order, cut signature, family-versus-drawn-graph compare (catches family/topology mislabels; formerly `topology-audit`, script `topology_audit.py`); member loomcheck/: Yangian/loom/fishnet applicability screen for position-space conformal integrals | selftest |
| `eichler/` | members of the Eichler.jl family (engine in `../upgrades/Eichler.jl`): `mirror6`, an exact mirror-map / instanton-number fingerprint engine for MUM points of D-finite operators with a published-number control battery; `theta9`, certified Siegel theta-constant observables; and `genus2`, genus-2 curve identification and recognition hygiene (Igusa–Clebsch invariants, Mestre reconstruction, recognition-free candidate arbitration, the PSLQ recognition law as structure) | partial |
| `ellipticus/` | Ellipticus: certified evaluator for elliptic multiple polylogarithms and iterated integrals on an elliptic curve from an exact curve/family spec (letters, word, point, precision in; certified values out); members ellred (standard-elliptic reduction) and gmtel (Gauss-Manin q-telescoper); formerly empl_eval | selftest |
| `emitall/` | claims-integrity harness: re-emit every quoted headline number from its receipt, fail on drift; document-side members `battery.py` (paper claims rows) and `paper_seams.py` (comment-swallow seam gate over a paper's tex tree) | selftest |
| `eras/` | ERAS — dependency-tamed certified enclosure forms in a shared parameter, with the consumer-pluggable adversarial verification battery | selftest |
| `famhar/` | multi-sheet family evaluator harness over a 2-var dlog connection (one JSON config per family; affine + polynomial letters) | selftest |
| `ffcapital/` | exact salvage of FireFly `ff_save` reconstruction state from dead or capped runs, fail-closed | data-gated |
| `formglue/` | glue around FORM 5: traces, color, term-stream rewriting, MZV evaluation, GRACE integration | partial |
| `frobenius-boundary/` | family-agnostic Frobenius-branch boundary machinery at var=∞ for non-analytic anchors | data-gated |
| `galois/` | motivic Galois coaction-cut engine + the exact MZV relation ring it reduces against | data-gated |
| `gatekeeper/` | Gatekeeper — data-integrity + held-out-CV honesty layer that gates every closure | selftest |
| `geotriage/` | GeoTriage — per-graph maximal-cut geometry and value-fittability triage with mandatory honesty gates | selftest |
| `gpl-eval/` | GPLEval.jl — GPL/HPL evaluator matching the fit's basis conventions at arbitrary precision | partial (GiNaC cross-check legs need the external GiNaC library; see GUIDE) |
| `holonomic/` | certified analytic continuation of holonomic ODEs with Arb ball enclosures (house recipe over the ore_algebra engine; the independent oracle the house transport evaluators check against) | partial |
| `kira-stack/` | IBP reduction to masters + raw DE, with the operational recipes that keep it alive at scale | smoke |
| `landau-alphabet/` | Landau Alphabet — face-by-face Landau / symbol-alphabet derivation from a graph spec, with a completeness-honesty contract; members: the capped msolve/Singular Gröbner escalation backend for heavy faces, and `LandauAlphabet.jl/`, the Julia iterated-integral bootstrap engine (weight-graded ansätze, exact constraint assembly over Q, LLL/PSLQ acceptance gates) behind GPLEval.jl | selftest |
| `lockpick/` | the integer-relation / lattice value-fit family: mplll + pslq_gate + the curated ring-basis member | partial |
| `longhand/` | Longhand — independent numerical ground truth for Feynman integrals, done the slow honest way: (A) the arbitrary-precision Feynman-parametric evaluator (fixed Gauss-Legendre product to 40+ digits at effective dimension ≤3, nested tanh-sinh, parallel CBC-lattice QMC, closed-form elimination of F-linear parameters, the UF-spec JSON front end) and (B) `longhand.disteval`, the pySecDec disteval driver (one compiled package at one Euclidean point in three checkpointed stages, a three-way lattice compare and an acceptance check with planted-fail and foreign controls; the generic-mass five-point double-pentagon top sector as the worked example under `examples/ndpent_top`) | selftest |
| `maxcut/` | maximal-cut kit: cut-geometry classification, canonical UT rotation, maxcut DE assembly | partial |
| `membound/` | soft-region (ω→0) boundary constants as frequency-space Bessel-kernel integrals for worldline/post-Minkowskian-type expansions: ε-layer moments of user-supplied two- and three-frequency ω-cores from a JSON spec (`python3 -m membound --spec core.json`), assembled into ω^(−kε)-layer boundary vectors, with retarded/advanced i0⁺ per frequency factor and a phase-only 'fey' diagnostic flag; worked example and gate: the PM memory family (c_M = 1) | partial |
| `mixalot/` | exact Bayesian evidence for mixture models, in a box (sha-pinned assembly package); members include the exact box-moment contiguity walker `mixalot/boxwalk/` (mod-p walk + CRT, `python3 -m mixalot.boxwalk`) | selftest |
| `nestor/` | certified-quadrature front door: node-farmed nested tanh-sinh with per-level precision ladder, plus singularity-subtracted dispersion quadrature (nestor.dispersion) with a worked one-loop-box kernel example | selftest |
| `numkin/` | freeze-one-scale rescue for stuck multi-var FireFly reductions (with basisland member) | smoke |
| `pmflow/` | self-consistent auxiliary-mass-flow CLI for cut-eikonal (PM) integral families | smoke |
| `popcorn/` | certified population-genetics likelihoods: exact/ball-certified selection-SFS vectors + gradients, certified DFE-mixing kernel, dominance; plus exact Lambda-coalescent expected SFS, float-validated two-locus moments and transient selected SFS under piecewise-constant N(t), an EPO ancestral-state join utility; exact linked two-site frequency spectra E[L_i L_j] under Lambda-coalescents with exact in/out certificates against every variable-size Kingman history, exact minor-allele folding of two-site spectra with class-membership deciders, and an exact two-site genotype projection + two-site spectrum estimator (VCF/TSV, cM or bp bins, block jackknife); certified (Arb-ball, model-conditional) enclosures of the transient selected SFS under piecewise-constant N(t) at small n, and msprime simulation harnesses — planted SFS/two-site truths with seeded synthetic genotype + polarization tables, and a two-window/two-locus Monte Carlo with jackknife errors | selftest |
| `posq/` | certified two-sided Bayesian evidence / Bayes-factor quadrature (BALLER member) | partial |
| `qinvert/` | exact consistency certificates for published summary statistics over partially released entity-level tables | selftest |
| `rankscreen/` | multi-prime parallel rank + inconsistency screen for exact-Q eliminations | selftest |
| `ratfit/` | exact rational-function reconstruction and gating suite for sampled path-DE data (Thiele/Newton/mod-p core, thiele_gate, li2close, design_grid), now with degree_alias: failure-axis triage for reconstruction failures, grid alias vs CRT height | partial |
| `seedling/` | pre-reduction front door for kira: certified predictive staging (pin/run/preflight/certify/census) + family preflight (`audit` symmetry group, `identity` transcription check) | partial |
| `subtropica/` | SubTropica — no-Wolfram pipeline over the HyperFLINT hyperlogarithm engine (one GPL-3.0 file, see its NOTICE) | partial (build the engine first) |
| `surd/` | Surd — exact hyperlogarithm integration of three-fold Cheng-Wu integrals on a one-parameter kinematic slice with radical (square- and cube-root) last-variable letters carried exactly: a-priori-alphabet fibration in Q[x,t,j]/(j²+1), function-level third step in Q(i)(t)[ρ]/(L), weight ≤ 3 hyperlogs with exact coefficients, structure tools and certified Arb reference evaluators (LO QCD/N=4 collinear four-point energy correlator as the worked case) | partial |
| `terrier/` | unified string-vacua/F-theory suite: certified CY period values, lattice wing, census wing | partial |
| `tropical-sampler/` | Borinsky–Sattelberger–Sturmfels–Telen (BSST) tropical importance sampling (arXiv:2204.06414 Algorithm 1, exact-fraction reimplementation) + `cegm_gj/` certified-value leg (tropical-cone Gauss-Jacobi, 30-60d two-run gates) | selftest |
| `trust/` | the reduction-verification triad — three checks with disjoint core lineages (vendored lpsyz + the strata streamed-IBP certifier + the receipt per-(point, prime) lambda-multiplier certificate tool for IBP/Laporta reduction tables, all inside); one command `sh tools/trust/selftest.sh` runs the import/pin gates, the receipt member battery and the strata member gates | partial |
| `vopclose/` | full-A function-level closure engine for a UT/graded 1-D path-DE | selftest |
| `wayfinder/` | Wayfinder — the eps-graded DE-transport pipeline: one evaluator contract, certified transport, Frobenius landings, the dlog-connection sampler wing, the engine-agnostic eps->0 limit extractor (epslimit: Laurent/Richardson extrapolation of an external engine's eps-grid) and the FLINT exact rational-function exporter (flintexport) | partial |
| `winnow/` | upgraded-Laporta F_p eliminator library with native lambda-witness receipts — Seedling's exact elimination engine (`import ibplapper`) | partial |

49 packages (`wayfinder` includes the sampler wing; `counterweight`
includes `canonical_form/`; `landau-alphabet` includes the Gröbner escalation
backend and `LandauAlphabet.jl/`, the Julia engine `gpl-eval` builds on).
The one operations package, Turnstile (admission control for long jobs on a
shared machine), lives under [`../ops/`](../ops/README.md) rather than here; it
is outside this index and outside `run_selftests.py`, and its self-test is
`bash ops/turnstile/selftest.sh` from the repository root.

## Not in this repository

- **Problem-specific code** — per-problem evaluators, results-table scripts, and
  single-result engines ship with their problem pages on
  [bootloops.ai](https://www.bootloops.ai), not here. The per-row closed-form
  evaluators are an instance: they ship with
  their problem pages — the sunrise, kite, and banana bundles on the site.
  Three packages the site's tool pages describe — Recount, Actuary and
  Trireme — are in the same position: as built they are tied to the projects
  they served, so they are not part of this release; their general kernels
  return in a later one.
- **External engines** — Kira, FireFly, FORM, AMFlow, Blade, FLINT, msolve,
  Singular, OSCAR, and the rest are catalogued with versions, licenses, and
  obtain-upstream notes in `../toolkit/external/TOOLS.md`. Our own engines
  ship in full under `../upgrades/`; the patched forks of Kira, Blade and
  AMFlow.cpp live in their own repositories, listed in
  [`../upgrades/ENGINES.md`](../upgrades/ENGINES.md).
- **JaCK & Jill** — the exact and certified Bayesian-evidence package for
  phylogenetic models (Python package `phyloexact`) lives in its own
  repository, `jackandjill` (published beside this one),
  with its own self-certification suite. The Mixalot `boxwalk` member's example that builds a
  Jukes–Cantor quartet vendors the two small kit modules it needs under
  `mixalot/mixalot/boxwalk/examples/jc_kit/`.

Some packages keep sha-pinned records (pin ledgers, battery logs) in-tree; these are the packages' verification chain and are not meant
to be edited — the selftests check the pins first and refuse to run on a modified
tree.

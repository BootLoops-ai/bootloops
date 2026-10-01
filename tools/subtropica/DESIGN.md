# subtropica — DESIGN

Julia rewrite (no Wolfram dependency) of SubTropica's front end, driving
the HyperFLINT C++ engine. This file records the component layout, the
design corrections and the acceptance gates. Source comments cite it by
tag: **[RT-i] .. [RT-x]** and **[RT-extra-*]** are design corrections —
places where a naive reading of the upstream `.wl` source is wrong and
the port deliberately does something else, with the wrong reading spelled
out; **R1 .. R11** are entries of the hazard register (section 4);
**C2 .. C15** are the components (section 1); **Gate A / Gate B** items
are the acceptance checks (section 3). Wire and handoff contracts are in
`CONTRACTS.md`; the divergent-route notes are in `DESIGN_B1.md`.

Source pins: `.wl` at upstream commit `adac2f72` (vendored copies in
`reference/`, see `reference/PROVENANCE.md`); engine binary v1.2.8 at
upstream commit `ead8c6e`; the engine tree itself (the `hyperflint.sh`
LD_LIBRARY_PATH wrapper and `data/mzv_reductions.json`) is built by the
user and treated as read-only. Note the skew between the two commits
(PROVENANCE.md). Engine paths are configured through `SUBTROPICA_HF_BIN`
and `SUBTROPICA_MZV_DATA` (defaults in `scripts/env.sh` and `HFConfig` in
`src/types.jl`). Every engine or helper process is started under
`ulimit -v 32505856` and is ended, when it has to be, by its own PID; run
output goes to a per-run directory outside the package
(`subtropica_runs/<tag>` by default).

Capability tiers, as named throughout the code:
- **A** — convergent driver: proven linearly-reducible Euler integrals
  from a raw quadruple, Euclidean region, `check_divergences` hard-ON.
- **B1** — log-divergent tropical subtraction (geometric property holds).
- **B2** — Nilsson–Passare continuation for power-divergent inputs and
  geometric-property violations, engaged only with `allow_continuation=true`.
- **C** — graph front end (shipped: `src/graph.jl`) and Method of Regions
  (not included in this release).

====================================================================
## 1. COMPONENTS
====================================================================
Edges: C2→C3→{C4,C10}; C4→C5→{C6,C7}; C6→C7→C8→C9; C10→C9; C9→C11→C12;
C13 gates everything; C14 hangs off C3+C4. The graph front end (C1) is not
part of tier A [RT-x]; it ships as the tier-C module `src/graph.jl`.

- C2 Euler integrand — `EulerIntegrand` in `src/types.jl`, constructor
  `euler_integrand` in `src/SubTropica.jl`. Tier-A input is the RAW EULER
  QUADRUPLE only: (prefactor, polys with a+b·ε exponents, ν, vars)
  [RT-x]. Fixtures are produced by the pySecDec `LoopIntegralFromGraph`
  glue in `scripts/make_fixtures.py`. Graph-layer conventions
  (STGetIntegrandData 2972 / STSymanzik 3104): eU=ν−(L+1)(D−2ε)/2,
  eF=−(ν−L(D−2ε)/2), prefactor Γ(−eF)/∏Γ(ν_i), the (−F)^eF sign tracked
  explicitly.
- C3 `src/tropical.jl` — STtoCoeffMonPols (multivariate factorization
  over Q via Nemo `factor` on `QQMPolyRingElem`), STforgetCoeffs, Trop I
  max-plus evaluation (eq 3.12/3.20), STEvalRay, initial-form restriction
  (eq 3.30), STPuiseux rescale series.
- C4 `src/polytope.jl` — [RT-iv] written against Oscar's FIRST-CLASS
  API: `newton_polytope`, `minkowski_sum`, `facets`, `vertices`,
  `normal_fan`, `affine_hull`; NOT raw Polymake property strings (the
  FACETS/VERTICES/MAXIMAL_CONES/AFFINE_HULL munging with the
  `rays := −FACETS[:,2:]` sign convention is the same Polymake underneath
  but re-adds the most error-prone string and sign handling). One raw
  Polymake escape hatch stays behind a flag for queries Oscar does not
  surface. Cache keyed by (coefficient-free poly product, vars). The
  1-variable dummy-projection case is included. Oscar is in-process and
  not thread-safe: the layer runs serially.
- C5 divergent-ray classification (in `src/tropical.jl`, predicates in
  `src/b1_types.jl`) — exact `fmpq` sign comparison of the factored
  Trop I at regulators→0; Σ_div via vertex-set argmax intersection
  (eq 3.24, STGetFaces 10255); w-vector search STProduceUs 10587 (exact
  rational nullspace, score {max|w|, count, sum}, w·ρ_f=−1 normalization,
  u = 1 − 1/(1+x^w) — the CODE convention, not the paper's v; checked
  once against the worked eq 3.37 example).
- C6 `src/continue.jl` — Nilsson–Passare continuation (eqs 3.22–3.23;
  STContinueRays 10545, STfindFirstNPContinuation 10852 subset search,
  STExpandIntegral 10906/11029 orchestration; FullSimplify → Nemo).
  [RT-ix] This is the tier-B2 core, engaged only via
  `subtropica_integrate(...; allow_continuation=true)`
  (test_b2_continue.jl, test_tworeg_driver.jl). On the B1 route
  (`allow_continuation=false`, the default) an input that needs
  continuation gets a LOUD typed refusal (`PowerDivergentRefusal` /
  `GeometricPropertyViolated`) — never a silent skip, never a partial
  answer.
- C7 `src/subtract.jl` — Möbius/inclusion–exclusion counter-term
  assembly (eqs 3.29–3.36; STSubtractionFormula 10723).
  [RT-i] **SIGN CORRECTION — MANDATORY**: the NET counter-term sign is
  **(−1)^|subface|**. The upstream code builds `Times[-a[subface],...]`
  with `a[s]=(−1)^(|s|+1)`, so net = −(−1)^(|s|+1) = (−1)^|s| (empty
  subface → +1 = identity term, single facet → −1 = subtracted). The
  "(−1)^(|subface|+1)" reading is WRONG — a transliteration of it ADDS
  single-facet counter-terms instead of subtracting them. See
  CONTRACTS.md WARNING-1.
  [RT-ii] **I/J CORRECTION — MANDATORY**: I = FIRST index subset with
  det(rays[fc,I]) ≠ 0; **J = Complement(1:n, I)**. The Kronecker-delta
  completion rows in Vol, the /∏x_J division AND the gauge-fix (vars ∉ J
  set to 1) all use **J**, the complement. Reading J = I yields finite,
  numerically plausible, wrong-volume wrong-gauge results, silent until
  the face-sum comparator. See CONTRACTS.md WARNING-2. Jacobian
  ∏(1−u_j)^(1−TropI(ρ_j)), Vol=|det|/∏(−TropI), x^(u_σ) monomial. Output
  grouped by face like Σ_div.
  [RT-extra-jac] The jacobian exponent uses TropI WITH regulators, i.e.
  it is ε-DEPENDENT — (1−u_j) factors enter C8's poly list with a+bε
  exponents; a unit fixture where a counter-term jacobian carries an ε
  exponent is mandatory in the subtraction tests.
- C8 `src/eps_expand.jl` — closed-form Laurent expansion of
  pref·∏x^m·∏P_j^(a_j+b_j ε) (STfastEpsSeries 17380): coefficients =
  rational functions × polynomials in log(P_j) (explicit exp(εb·logP)
  multinomial structure, no general Series); minOrder = |vars(last face)|
  − |vars(first face)| (SubTropica.wl 11188), maxOrder = requested +
  prefactor pole offset; Laurent pole-padding bookkeeping (known pitfall:
  a truncated Laurent series must be guarded — a ratio-style check that
  requested order + padding covers the deepest pole, R2). [RT-vii] **the
  Γ-prefactor symbolic ε-series is C8's job**: no library call gives the
  Γ(a+bε) series symbolically in γ_E/ζ(k) in Julia, so
  `loggamma_series(a+bε, order)` is hand-coded from the standard
  ψ-expansion, emitting EulerGamma and zeta tokens as `HlogExpr` atoms
  (closed recursion).
- C9 `src/hf_bridge.jl` (+ `src/hlogexpr.jl`, `src/serialize.jl`) — the
  engine boundary (CONTRACTS.md §a): request builder (flat JSON, ONE
  object, stdin), response parser (stdout LAST line only), in-band error
  taxonomy, env contract, coefficient-string → `HlogExpr` token parser
  (mzv_a_b_c → zeta tokens with 'm' = negative argument, Log2, delta[v],
  Wm_i/Wp_i/sqrt_disc_i via the returned algebraic_letters table).
  [RT-v] **Wm/Wp handling routes through the engine's OWN ops**:
  `simplify_with_vieta`, `back_substitute`, `combine_wm_wp_ratios` — the
  port does NOT reimplement the Vieta/Sqrt[disc] rewrite of upstream .wl
  12930–13005 (a semantics-drift surface). Only the token→HlogExpr
  PARSING is local. Residual ZeroInfPeriod keys are evaluated via the
  `zero_inf_period` op; [RT-extra-period] every residual key is ALSO
  cross-evaluated via ginac (GPLEval.jl/ginac_gpl) to ≥60 dps as a unit
  gate — this is our own code with no upstream to compare against.
  [RT-vi] **VERSION GATE = STAMP + BEHAVIORAL CANARY**: source builds of
  the engine stamp `hf_version=1.2.0` unless HF_VERSION is set at cmake
  time, so a stamp-only "≥1.2.8" gate fails or misleads on a source-built
  C-ABI library. Gate = (stamp: `hyperflint.sh --version` → `HF_VERSION:`
  line and/or response-envelope `hf_version`) AND (canary: the
  unary-minus precedence fix that motivated the pin — `-x^2` must mean
  −(x²)). Canary, verified on the v1.2.8 binary:
  `{"op":"eval","a":"-x^2","vars":["x"],"values":["3"]}` → `"result":"-9"`
  (pre-fix builds silently drop the sign → `9`). Both checks run at
  bridge init; either failing is fatal.
- C10 `src/lr_dispatch.jl` — thin wrappers on `find_lr_orders`
  (MULTI-GROUP `groups:[[..]]` per-addend semantics, MANDATORY for
  counter-term sums — a NOLR on a counter-term sum means multi-group was
  forgotten; `groups` is asserted for sums), `factor_table`; per-face
  best_order persistence + REPLAY on resume (R3). PIN CORRECTION: there is
  NO standalone `verify_order` op in the CLI dispatch — order
  verification is the OPTIONAL `"verify_order":[v1,...]` FIELD of the
  `find_lr_orders` request (handlers.cpp:741). stDeriveGaugeFromHomogeneousLR
  port. `find_lr_orders_scan` (Cheng-Wu gauge scan) is an engine op
  (CONTRACTS.md) without a wrapper in this release. `src/lr_refine.jl`
  adds a compatibility-graph refinement of the letter bound.
- C11 counter-term farm — not included in this release: the face loop in
  `src/b1_driver.jl` runs serially in-process with one engine call per
  face × order (see DESIGN_B1.md §4b). The design for a process pool
  (setsid-launched CLI workers, `ulimit -v 32505856` on every worker,
  FindRoots escalation strict → algebraic_letters:true → carry_discharge,
  narrow-ctx retry with the flag off once and no recursion,
  `{divergent:true}` fatal) is recorded here for the scale-out seam
  `_b1_face_result`. Code paths that would need the farm refuse with a
  message naming C11.
- C12 `src/assemble.jl` — sums counter-term results per face per ε
  order, letter-table union, NP prefactors, normalization (e^(LγEε)), final
  ε-series to the output order; `LaurentSeries{HlogExpr}` result with
  provenance metadata.
- C13 `src/verify.jl` + `scripts/make_fixtures.py` — verification harness.
  [RT-iii] **ORACLE SET CORRECTION — MANDATORY**: the ≥30-digit gate
  oracles are **{an AMFlow-class numerical evaluator, a high-precision
  sector-decomposition evaluator (tools/longhand, route A; oracle key
  `:hiprec_sectordecomp`), mpmath
  tanh-sinh quadrature (low-dimensional convergent cases)}**. pySecDec is
  SANITY-ONLY (QMC/MC, ~1e-3..1e-10 — it cannot reach 30 digits on
  anything nontrivial) and NEVER a gate oracle; FORM is not a
  Feynman-integral evaluator and is off the oracle list. **Fixture-value
  comparison is FORBIDDEN in gates**: comparing against a value stored in
  an upstream fixture (the engine's test/Smirnov cases, the 22-case Long
  suite's embedded references — upstream verified them only against
  pySecDec at relErr<1e-3) validates plumbing, not correctness; such
  comparisons are labeled REGRESSION, never VALIDATION. Also here: the
  rational-point picker avoiding letter degeneracies (STVerify 6258
  logic), the Laurent comparator, and the synthetic-truth control
  injector (mandatory: known-answer positive control + coefficient-
  perturbed negative control; perturbations NON-RATIONAL, ×(1+δπ),
  because rational perturbations alias into exact relations).
- C14 Method of Regions (STGetRegionVectors + STMoRExpand; would reuse
  C3/C4) — not included in this release.
- C15 `src/SubTropica.jl` — module shell, top-level API, option
  validation, tier routing.

Upstream material not ported: GUI/web-server/library layer; Maple/ginsh/
FIESTA/AMFlow wrappers; HyperIntica.wl + mzv.wl (already in the C++
engine); notebook/paclet/LibraryLink plumbing; the PV one-loop library
fetch; AutoRationalize M1/M2/M3 (off by default upstream; R5); the .wl LR
search STFasterFubini2 (superseded by the engine's `find_lr_orders`); the
benchmarking harness (its cases are extracted as fixtures); the Normaliz
triangulation path (legacy/diagnostic).

====================================================================
## 2. INTERFACE CONTRACTS
====================================================================
Pinned in `CONTRACTS.md`, with file:line citations pointing at the
vendored `reference/` copies [RT-extra-vendor]. Summary: (a) the engine
eval-json wire contract; (b) the locally-finite integrand handoff; (c)
the two corrected pins as WARNINGS; (d) the HlogExpr serialization
grammar. Transport is CLI-first (`hyperflint.sh eval-json`); the OPT-IN
in-process library (`deps/libhf_cabi_static_f2.so`, statically isolated,
built by `deps/build_cabi.sh` from a patched engine tree; see
`deps/BUILD.md`) covers only the C-ABI ops and MUST pass the dlopen gate
(strict 1.2.8 stamp + envelope + canary + pair probe) because source
builds mis-stamp 1.2.0.

====================================================================
## 3. ACCEPTANCE GATES
====================================================================
GATE A (tier A, exact):
(i) reproduce ≥2 published convergent benchmarks — the eq 4.1
Møller-box-class integral and the eq 4.12 Smirnov period — with the
symbolic output evaluated via GPLEval/ginac_gpl + mpmath and matched
against an INDEPENDENT oracle (never used in any fit) to ≥30 genuine
decimal digits at a non-degenerate rational kinematic point, the floored
integer digit count recorded. **Comparison against a fixture's stored
value is FORBIDDEN in this gate** [RT-iii].
(ii) synthetic-truth control: one integrand with a known closed form
passes; one coefficient-perturbed (×(1+δπ), non-rational) negative control
FAILS the comparator.
(iii) MZV provenance gate (R7): 25 random rules of mzv_reductions.json
verified LHS=RHS to ≥60 dps via ginac (an independent implementation).
(iv) version gate wired (stamp + canary [RT-vi]); one deliberately
divergent standalone input returns `{"divergent":true}` and is surfaced
as fatal.
(v) residual-key unit gate [RT-extra-period]: every residual ZeroInfPeriod
key in gate outputs cross-evaluated, `zero_inf_period` op vs ginac, to
≥60 dps.

GATE B (tiers B1/B2, exact):
(i) reproduce the published divergent benchmarks eq 4.3 (eikonal, 1/ε⁴)
and eq 4.11 (small-x, 1/ε³): ALL pole orders present with correct
multiplicity, every Laurent coefficient ε^k (deepest pole through ε⁰
minimum) matched ≥30 digits against an independent numerical evaluator at
a rational Euclidean point (floored digits recorded per order).
(ii) one further divergent target of the class where sector
decomposition stalls, cross-checked ≥30 digits through ε⁰.
(iii) Möbius-layer unit gate: pole-part CANCELLATION — the sum of
counter-terms per face is finite where the theorem says so (numeric probe
at the boundary). This is the gate that catches the sign error [RT-i] and
the I/J inversion [RT-ii].
(iv) synthetic divergent control (hand-computable 1/ε² Laurent) +
perturbed (×(1+δπ)) negative control.
(v) leave-one-out [RT-iii]: at least one benchmark verified ≥30 digits by
an oracle NEVER consulted while debugging that benchmark (sector
decomposition or direct mpmath quadrature).
(vi) jacobian-ε unit fixture [RT-extra-jac] present and passing.

Regression suite (tier C): extracted upstream benchmark cases pass
end-to-end against their embedded references — **labeled REGRESSION, not
validation** [RT-iii].

====================================================================
## 4. HAZARD REGISTER
====================================================================
R1 (a·b)^(cε) power-splitting / positivity. [RT-viii] **CONCRETE RULE —
MANDATORY** (a sign ledger for numeric leading coefficients ALONE is
insufficient: F Euclidean-positive overall does not make each Q-factor
coefficient sign-definite; (x−y)²-type factors exist; Nemo has no
Refine): after Nemo `factor`, for EACH factor run (1) the per-factor
coefficient-sign test via `coefficients`/`exponent_vectors` — all
coefficients of one sign under positive-kinematics substitution ⇒ sign-
definite, split allowed; (2) otherwise probe at N random positive-rational
points (N≥8, fresh randomness per run) — consistent sign ⇒ split allowed
with the probed sign recorded in provenance; (3) ANY sign change ⇒
REFUSE-WITH-ERROR (typed, names the factor). Implemented in eps_expand.jl,
unit-tested on an (x−y)² factor. The ledger covers POLYNOMIAL factors,
not just numeric leading coefficients. The (−F)^eF branch convention is
unit-tested with a numeric probe. Physical-region input is refused (R4).
R2 Laurent bookkeeping (the truncated-Laurent pole guard): carry
minOrder from the face codimension analytically; check that requested
order + padding covers the deepest pole; property tests vs sympy series.
R3 LR-order variance: `find_lr_orders` is a time-budgeted branch-and-bound
(HF_LR_TIME_BUDGET_S=180 default) — different budgets pick different
orders, changing intermediate alphabets (NOT the answer). best_order is
persisted per face (`best_order.json`) and REPLAYED via the
`verify_order` FIELD of find_lr_orders; the ≥30-digit numeric gate
re-verifies after any order change; NOLR on a counter-term sum ⇒
assert-fail (multi-group forgotten).
R4 iπ / physical region: not supported in this release
(`region=:physical` refuses). If added, delta-residue signs would be fixed
by matching a high-precision numeric oracle at ≥2 kinematic points,
recorded as `:numerically_matched_sign` + evidence and never silently
promoted to analytic; [RT-extra-r4] those ≥2 points MUST be disjoint from
at least one HELD-OUT verification point — otherwise the "evidence" is
the fit. Gates A/B are Euclidean-only and never touch this path.
R5 AutoRationalize: not ported (off by default upstream, genuinely open);
the engine's algebraic_letters covers quadratics; FKV/Källén out of scope.
R6 Divergence-check ambiguity: hard-ON standalone, OFF only inside face
sums; Gate A(iv) is the positive divergence-detection test; a multi-pole
CONVERGENT integrand detects false positives.
R7 MZV table provenance: Gate A(iii); table path pinned in scripts/env.sh
+ config; any table update re-triggers the provenance test.
R8 Ownership: the upstream engine tree is read-only; both commits + the
binary sha256 are recorded (PROVENANCE.md); source builds only under
`deps/`.
R9 Polymake only via Oscar/Polymake.jl in-process (no standalone binary;
the jll is broken standalone); Oscar-native fallback.
R10 Upstream churn: pins + skew recorded (PROVENANCE.md); stamp+canary
runtime gate [RT-vi]; ANY upstream update re-runs the full fixture suite
[RT-extra-vendor].
R11 Wrapper corner cases: the FindRoots cascade is an explicit small state
machine with unit tests.

====================================================================
## 5. SCOPE
====================================================================
Tier A alone is a symbolic evaluator for finite linearly-reducible
Euler/Feynman integrals in the Euclidean region: an independent symbolic
oracle beside GPLEval for PSLQ-free closure of finite masters,
energy-correlator-class Euler integrands, finite UT-basis members where a
numerical evaluator gives only numbers, MZV/period identities on demand.
B1 [RT-ix] adds log-divergent Laurent-in-ε results symbolically — the
class where sector decomposition stalls. B2 completes the eikonal/small-x
divergent classes (eq 4.3/4.11 patterns).

====================================================================
## 6. TOP-LEVEL API (src/SubTropica.jl)
====================================================================
```julia
subtropica_integrate(quad::EulerIntegrand; order=:auto, dps=60,
                  gauge=:auto, region=:euclidean, check_divergences=:auto,
                  findroots=:cascade, workers=1, hf=HFConfig(),
                  allow_continuation=false, continuation_order=16)
    -> LaurentSeries{HlogExpr}
euler_integrand(prefactor::Prefactor, polys, nu, vars; kinvars=Symbol[])
    -> EulerIntegrand           # raw quadruple constructor (tier-A input)
verify_laurent(L::LaurentSeries, oracle::Symbol; point=:auto,
               min_digits=30) -> GateReport
    # oracle ∈ (:amflow, :hiprec_sectordecomp, :mpmath)   [RT-iii]
    # :pysecdec accepted ONLY with sanity=true, never in gates
hf_call(op::String, fields::Dict; cfg=HFConfig()) -> Dict
hf_gate(cfg=HFConfig()) -> Bool   # stamp + unary-minus canary [RT-vi]
loggamma_series(arg::EpsExp, order::Int) -> Vector{HlogExpr}  [RT-vii]
```
`region=:physical` refuses (R4). `check_divergences=:auto` = hard-ON for
standalone calls, OFF inside face sums only. `workers` is accepted and
recorded; the face loop is serial in this release (C11).

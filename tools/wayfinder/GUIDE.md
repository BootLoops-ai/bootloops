# Wayfinder — guide

Wayfinder — `tools/wayfinder` — https://bootloops.ai/tools/wayfinder.html

KIND: package (de_load, transport, frobenius, frob_scalar, boundary_branches,
two_sector, epsfan, ratfun, boundary, gate, manifest, quad, tailcut, conic_thirdkind,
acb_fast_attach, sampler + targets/ — the dlog-connection sampler, epslimit — the
engine-agnostic ε→0 limit extractor, flintexport — the FLINT exact rational-function
exporter, tests/).

PURPOSE: Eps-graded DE transport: load a saved DE system into ONE evaluator
contract, transport boundary/anchor data along it (fixed-ε Taylor march / η-march),
land on Frobenius bases at regular-singular points, extract ε-Laurent by Cauchy fan or
real-grid Vandermonde, seed AMFlow-convention vacuum boundaries, and gate everything
with two-precision digit gates + per-literal FEED tables + sha256 manifests.
Throughout the package's code and docs, a **saved** (or reference) system means a
sha-pinned reference DE system on disk.

USE-WHEN:
- Shared exact-rational transport on a saved DE file: fixed-ε endpoint transport,
  η-march, Frobenius landing (matrix or scalar resonant towers), ε-Laurent
  extraction, FEED gate.
- Build canonical DE + transport + gate vs oracle for a new target family
  (invariant-derivative ops + integrability check).
- ε-graded m²-connection transport from the m²→∞ massless limit boundary.
- Fixed-ε boundary transported along a kinematic line (`transport_fixed_eps`).
- Fixed-ε nodes → ε-Laurent + oracle gate (`epsfan`).
- Third-kind-fenced CONIC connection row (dL/dt = [H·Y]_A^B + c·M + ∫_A^B T/Y dx over
  Y²=q(x,t) QUADRATIC in x, T≠0, t-poles MOVING over the family) needs an exact
  rational t-connection for the march → `conic_thirdkind` member (below).
- An external engine (AMFlow, pySecDec, …) already produced values on an ε-grid and
  you need the ε→0 Laurent coefficients, or a series must be summed AT its radius
  of convergence with known tail exponents → `epslimit` member (below; not
  `epsfan`, which owns the evaluator and chooses its own nodes).
- sympy `cancel`/`together` is the wall in exact DE-connection assembly or bulk
  rational-function linear algebra over Q(x1..xk) → `flintexport` member (below).

MEMBER — conic_thirdkind (engine core): exact extended t-system builder for
third-kind fenced conic rows — adjoins third-kind atoms J_{f,k}=∫_A^B x^k/(f·Y) dx per
(irreducible x-dependent factor f of den T, k<deg f) + M + endpoint carriers
Bc_X=√q(X,t) + accumulator F'=L (n~9-14); rows closed EXACTLY (partial fractions
residual==0 cert; Hermite cancel==0 certs); output = exact rational A(t) →
`transport_fixed_eps` with declared singular points = denominator roots = the MOVING
POLES; anchors/gates = `quad.quad_refine`. FULLY PARAMETERIZED:
fiber/connection symbols, panel endpoints, and the connection-input schema
(`rows_from_connection` adapter + `subs` spectator pins — pin FIRST) are all
arguments. Import: `from wayfinder import conic_thirdkind` (LAZY module — package
import stays sympy-free). SCOPE (named, fail-closed — do NOT silently generalize):
(1) CONIC ONLY — deg(q,x)==2; deg>2 (elliptic/hyperelliptic) is a basis EXTENSION,
not a parameterization; `ConicScopeError` raises at engine AND builder entry;
(2) T-den x-dependent factors squarefree (mult≥2 = the unbuilt Hermite loop
extension; loud assert). Battery: `tests/test_conic_thirdkind.py` — all legs fast
(seconds): engine identity vs independent 2-precision quadrature, deg-3 refusal,
comparator+certificate mutation controls, toy truth case (a small synthetic fenced
row in `tests/fixtures_conic_thirdkind/` rebuilt through the parameterized API,
compared entry-exact against its saved reference build, and cross-checked row-by-row
against recomputed quadrature), ExtSys numeric leg.

MEMBER — epslimit (`epslimit.py`; formerly the standalone eps-extrapolator
package): engine-agnostic ε→0 limit extraction. PURPOSE: given any engine's
(ε_i, value_i) grid — values as python numbers, decimal strings, Arb-ball strings
`'[mid +/- rad]'` or `{'re':..,'im':..}` dicts; AMFlow result JSONs via
`load_amflow_grid`, which drops planted failure markers — return the ε-Laurent
coefficients c_{k,j} of Σ c_{k,j} ε^k log(ε)^j with the finite part c_{0,0} the usual
target. The pole order k_min is auto-detected from the log-log slope; a log(ε)
resonance is accepted only on a decisive leave-one-out improvement; the working
precision is set from the grid's dynamic range and node clustering; four methods
(rescaled generalized Vandermonde, Neville/Romberg power elimination, Bulirsch–Stoer
rational extrapolation, ε-FFT on a circle of nodes) are run as the geometry allows and
CROSS-VALIDATED — the reported achieved-digit figure is method disagreement, not a
self-reported residual. Second verb: `richardson_boundary(partial_sums, Ms,
tail_exponents=None)` sums a series AT its radius of convergence when the boundary
singularity type is known (√-type → half-integer M^{−α} tails, the default), by
Richardson-eliminating the listed tail powers from the partial sums (37.7 digits
measured at a lattice-Green-function threshold; this is the verb
tools/frobenius-boundary points at for disk-edge sums). INVOKE: `from wayfinder.epslimit
import extrapolate, richardson_boundary, load_amflow_grid` (or `from wayfinder import
epslimit`, a LAZY module); `r = extrapolate(samples, max_log=2)` → `r.finite_part`,
`r.c(k, j)`, `r.k_min`, `r.log_order`, `r.achieved_digits`, `r.diagnostics`; CLI
`python3 tools/wayfinder/epslimit.py out1.json out2.json --part re` (or, from tools/,
`python3 -m wayfinder.epslimit …`). RELATION TO epsfan: epsfan extracts the Laurent fan
from a DE system Wayfinder can evaluate at ε-nodes of its own design, with certified
alias bounds; epslimit sits downstream of an engine you do not control and takes the
grid as given — the two share no code and answer different questions. LIMITS: achievable
digits are information-theoretic in the grid — 3 ε-samples cap at ~6–8 digits whatever
the algorithm (3 points resolve only ε^0..ε^2, so the finite part carries an O(c_3·Πε_i)
truncation error); for ≥30-digit finite parts take ≥6–8 samples on a geometric ladder
(ratio 2–4) at engine working precision ≥2–3× the target digits. Against real engine
grids the remaining error tracks the engine's internal precision, not the extrapolator.
It does not duplicate engine-side ε=0 resonance fixes (e.g. in the AMFlow.cpp fork,
the sibling repository amflow-cpp). Precision handling differs from the rest of the package:
`extrapolate()` sets the global mpmath precision to its working value for the duration
of the solve and restores it in a `finally` block. Requirements: mpmath (gmpy2 backend
recommended, optional). Battery legs (`pytest tools/wayfinder`, seconds):
`tests/test_epslimit_synthetic.py` — 6 legs / 12 graded checks against known exact
expansions (finite part >30d, pole c_{−1}=ζ(3) and c_0=π² >25d, ε·log ε resonance with
the log order auto-detected, ε-FFT on a circle >30d, narrow 3-point exact quadratic
>60d, and the truncation-floor demonstration: 3 pts ~8d → 8 pts >60d);
`tests/test_epslimit_fixtures.py` — the file pipeline on committed synthetic fixtures
(`tests/fixtures_epslimit/`: a 3-file exact-quadratic grid, c_0=7/2 to >40d; an 8-point
pole ladder, c_{−1}=ζ(3) >25d and c_0=π² >20d, with a planted engine-failure sample the
loader must drop).

MEMBER — flintexport (`flintexport.py`; formerly the standalone flint-export
script): FLINT `fmpq_mpoly` fraction walker for bulk EXACT rational-function algebra over
Q(x1..xk) — DE-connection assembly, reduction-row contraction, big substitute-and-compare
— and the exporter behind the monomial-json connection format `load_monomial_json`
reads. PURPOSE: walk each sympy expression tree ONCE with values held as (num, den)
`fmpq_mpoly` pairs, gcd-reducing on add and cross-cancelling on multiply, with leaf
substitutions (d → 4−2ε, kinematic pins) applied through an env; sympy supplies only the
input tree, all arithmetic is FLINT (measured ~1000× over sympy `cancel`/`together`:
285.6 s for a full 22×22 three-variable connection where sympy took 128 s for one 7-entry
row). INVOKE: `from wayfinder.flintexport import mk_ctx, gen, const, poly, make_env,
walk, canonical_terms, Frac` (or `from wayfinder import flintexport`, LAZY);
`ctx = mk_ctx(('eps','s2','s12'))`; `env = make_env(ctx, {sym: gen/poly/const …})`;
`fr = walk(sympy_expr, env, cache={}, ctx)`; `num_terms, den_terms =
canonical_terms(fr)` (exact integer-coprime export form, denominator lead positive;
`None` for an exact zero). Coefficients int/Fraction/fmpq; any number of variables;
INTEGER exponents only — a fractional or symbolic power raises `ValueError`, never
silently floors. LIMITS: not a general CAS; keep every operation inside FLINT and cancel
once at export — dropping back to sympy arithmetic mid-walk reintroduces the wall.
Requirements: python-flint (pip) and sympy; the rest of Wayfinder treats python-flint as
an optional accelerator, so this member is exposed lazily and its battery leg skips by
name without it. Julia sibling: tools/counterweight `src/MpolyFeed.jl`
(`feed_matrix_from_chunks`). Battery leg: `tests/test_flintexport.py` — the module's five
exact-identity checks (exact zero after expansion, gcd cancel, a three-variable
(d,s2,s12) pin identity, canonical export shape + positive-lead law, non-integer-exponent
refusal) plus a fresh-interpreter check that package import stays sympy-free and does
not load either lazy member; standalone: `python3 tools/wayfinder/flintexport.py`
→ `wayfinder.flintexport selftest: PASS`, rc=0.

NOT-FOR: deferred (all loud raises) — kira loader cuts / symbolic invariants /
multi-family; `vacuum_ending_seed` partial-η + non-builtin ending masters (ending
IBP); complex-ε Frobenius matching (Im-linear flag raises); two_sector ε-Laurent
window + ε-offset resonance pivoting (fixed-ε by design). Not an IBP or
ε-factorization tool.

INVOKE: `from wayfinder import load_monomial_json, load_amatrix_json,
load_kira_targets, transport_fixed_eps, frobenius_basis, land, cauchy_laurent,
vandermonde_laurent, anchor_roundtrip, tadpole, singlemass_vacuum,
vacuum_ending_seed, gate_strings, feed_gate_table, write_manifest, check_manifest,
RF, two_sector_series, eval_two_sector, quad_refine` (+ lazy modules `frob_scalar`,
`boundary_branches`). Contract: `DESystem` with `.n/.var/.A(x,eps,dps)/.meta`;
optional `.singular_points` (loaders do NOT populate — set it when you know the
poles) and `.A_series` via `DESystem.enable_fast_path(eps, wp)` (exact-rational
ratfun.RF layer; kills the ~N A-evals/step generic cost — a 79×79 system is
minutes/step without it). BIG-n FAST PATH: `enable_fast_path(eps, wp, backend="acb")`
+ `transport_fixed_eps(..., backend="acb")` — flint acb ball step kernel
(acbfast.py), chassis semantics/certificates identical, state mid-trimmed per
accepted step (NOT ball-certified — same honesty class as mpmath, gates unchanged);
measured n=378/dps38: 10.7–11.8 s/step vs ≥60–120 mpmath; small-n parity 47–49d.
Default backend untouched. `entry_exact` is a debug hook, NOT contract API. Certified
Laurent: `epsfan.vandermonde_laurent_certified(..., eps_nodes1=, eps_nodes2=,
k_certify=(ka,kb))`.

INPUTS: three accepted formats — monomial-json (exact ε-monomials), amatrix-json
(rational-in-d y-poly fits), AMFlow `kira_target.m` (η-DE from reduction rules;
masters-file order = AMFlow boundary order; tadpole-scaling convention check built
in). Explicit dps everywhere; global mp.dps never touched. Pure mpmath+stdlib;
python-flint OPTIONAL accelerator (HAVE_FLINT); sympy only lazily inside
frob_scalar/boundary_branches/conic_thirdkind.

OUTPUTS: transported vectors + `trunc_worst` truncation certificates; Frobenius
landing matrices; ε-Laurent tables with `min_agree_digits`; FEED verdict tables at
charter bars (zero-vendored refs gate absolutely; 1-digit refs refused); sha256
manifests (bytes, never mtime).

MEMBER — sampler (`sampler.py` + `targets/`): samples exact ∂J/∂x_i data at
exact-rational kinematic points (via amflow-cpp `solve_integrals` finite differences
or the `diffeq` backend's exact IBP-derivative mode) and fits a dlog connection A_a
per alphabet letter, with held-out digit verification and honesty fields
(`n_points_used`, `cond_number`, `fit_residual_digits`, `min_heldout_digits`,
`failed_to_rationalize`). Needs `heldout_cv` (gatekeeper) and `canonical_form`
(Counterweight member) importable — sibling packages under tools/ resolve
automatically; clear ImportError otherwise. The amflow legs need an amflow-cpp
`amflow_cli` build (from the fork, the sibling repository amflow-cpp); set `AMFLOW_CLI`.
A target spec may declare its own `ibp_cache_dir`; `AMFLOW_IBP_CACHE` sets the
default cache dir. Battery: `python3 tests/test_sampler.py` (or the registered
`pytest tools/wayfinder`, which collects the `test_t*` functions) — T1
(synthetic round-trip), T4 (exact A_a from 2 pts) and T5 (negative control: one
perturbed ∂J value must fail to rationalize or drop the held-out digits) run as
shipped; T2/T3 (box1l end-to-end, ~10-15 min) are gated by
`WAYFINDER_RUN_AMFLOW=1` and skip otherwise.

MEMBER — acb fast-attach (`acb_fast_attach.py`): attacher pair for
scalar-operator companion systems with EXACT rational coefficients — the
scalar-companion route `de_load` does not build (NOT a replacement for
`DESystem.enable_fast_path` on matrix systems; that path already exists).
Contract: a companion duck-type exposing `.n` / `.coeffs` (n+1 Fraction lists,
`coeffs[j][k]` = coeff of x^k in p_j, operator Σ_j p_j(x)·D^j) /
`.A(x, eps, dps)` / `.meta['singular_points']`. `attach_fast_path(sysm)` binds
a duck-typed mpmath `A_series` (numeric Taylor shift of the exact p_j +
truncated series division, per `transport._local_a_series`) — also what makes
the mpmath control route comparable step-for-step with acb; part of the member,
keep it. `attach_acb_fast(sysm, wp)` builds the `acbfast.AcbFastTable` at
working precision wp DIRECTLY from the companion entries -p_j/p_n
(shared-denominator dedup inside), stores `sysm._acb_fast`, wires
`sysm._acb_check_eps` so nonzero eps and dps > wp are REFUSED; requires
python-flint (ImportError otherwise). Then march:
`transport_fixed_eps(sysm, 0, x0, x1, y0, dps, backend="acb")`. Fixed at
eps=0 (pure-Q operators). Import: `from wayfinder import attach_fast_path,
attach_acb_fast` (or flat with the package dir on sys.path). BACKEND LAW
(measured): backend="acb" for every large-order high-dps march — same short
march, same steps, same truncation certificate; an order-17 march at dps 430
measured mpmath 73 s -> acb 2.2 s (33x). Battery:
`tests/test_acb_fast_attach.py` (exact log-toy control; the acb legs skip
without python-flint).
- F1 — PRECISION LAW (measured; binds callers of frob_scalar/transport on
  EITHER backend): dps_work ≥ ~2.2·Nmax + target digits, where 2.2 is the
  measured worst-case per-term digit loss when apparent singularities crowd
  the expansion point — the coefficient recurrence carries parasitic modes
  ~(1/d)^N for d the distance to the nearest leading-coefficient root.
  RECOMPUTE the local rate as log10(1/d_nearest-lc-root) at EACH expansion
  point; never reuse one point's rate.
- F2 — GAUGE-RADIUS LAW: a saved gauge representative can have GENUINE
  solution singularities inside its apparent-singularity cloud — check the
  SOLUTION's convergence radius at the anchor BEFORE trusting a Frobenius
  series there. The nearest-equation-singularity distance is only a LOWER
  bound (Cauchy); the solution's actual radius is a MEASURED dominant-shell
  selection among the discrete spectrum of characteristic roots (reciprocals
  of the nonzero leading-symbol roots — exact identity). Say "measured
  radius"/"certificate", never "proved". Per-gauge forward-stability figures
  (digits/term) are gauge-specific measurements, not laws — re-measure per
  gauge before pricing by them.
- F2 noise fence (dps-scaling): a fixed-dps forward recurrence shows a
  spurious slope-break toward the nearest-shell rate at break-N ~
  dps/(per-term loss); a break that MOVES OUT with working precision is
  rounding noise — never read a radius at one fixed dps.

ENV: OPENBLAS_NUM_THREADS≤2 MANDATORY for jet-transport-class jobs (else thrashes).
SERIALIZER-WIDTH POLICY: serializers write working-precision-derived
widths — never fixed nstr(...,N); carry honest_cap + serialized_digits fields;
dps-keyed caches. `WAYFINDER_REFERENCE_FIXTURES` — points the optional reference-fixture
test legs (large production DE files, not distributed with the package) at a local
dir; without it those legs SKIP.

GATES: suite 57/57 green; tadpole control ~152d vs direct Gamma
oracle; 2F1 fixed-ε control 101d generic / 148d fast path vs held-out hyp2f1; a
26-key two-loop five-point production chain at full stored ~70d; a resonant-tower
production constant at 260d two-precision. Always run `anchor_roundtrip` (extraction-error isolator, boundary
cancels) BEFORE transporting. ε-graded-transport gate ceiling is the ε₀-alias ∝ε₀⁷:
1e-7→22d, 1e-9→36d, 1e-10→~43d — tighten eps0, it is NOT a Vandermonde-dps issue.
PUBLIC-RUNNABLE test legs (verified in this repo): test_de_load
(synthetic groups), test_manifest, test_gate, test_quad (14/14), test_tailcut
(12/12), test_boundary (B1-B5 + refusals), test_conic_thirdkind (all legs),
test_acb_fast_attach (exact log-toy control; acb legs skip without
python-flint), test_epslimit_synthetic + test_epslimit_fixtures (8 legs, ~4 s),
test_flintexport (skips by name without python-flint/sympy) — all green. REFERENCE-FIXTURE legs (SKIP without `WAYFINDER_REFERENCE_FIXTURES`):
large production-file loads (a monomial-json cone connection and its closure, a 79x79
amatrix-json fit, a kira_target.m chain, a triple-tadpole B4 cross-check) and the heavier controls
(control_hypergeom 101d, full-suite timing).

FOOTGUNS:
- Frobenius shear path at scale (frobenius.py): `_solve_sylvester` n>6
  = Bartels–Stewart on ONE cached mp.schur (a dense LU would be O(n⁶) bignum per Taylor order; n≤6 takes a direct path).
- Defective clusters: `cluster_tol=` kwarg on `frobenius_basis`/`land` (default 1e-(wp/2)). Jordan pseudospectrum rings sit AT the default tol and shatter
  defective clusters into singletons whose projectors then rightly fail the quality
  gates. Use cluster_tol ~1e-40: above any ring at wp≥120, below physical O(δ)
  exponent gaps. Per-cluster scalar towers CANNOT represent integer-ladder classes
  with defective rungs (cross-rung log mixing) — use the shear engine and fix its
  gates; never rebuild driver-side landing.
- `frobenius` is REAL-ε0 only (Im-linear flag raises on complex ε).
- epsfan `cauchy_laurent`: trust `min_agree_digits`, not dps; FULL circle unless the
  row is Schwarz-real — `schwarz_real` mirror is OPT-IN, complex rows silently
  corrupted otherwise.
- `vandermonde_laurent` node DESIGN RULE (both mandatory, measured): geometric
  q-ratio ≈10 (narrow band ⇒ Vandermonde cap ~25d) AND ε_max ≈1e-3 (tail alias ⇒
  ~29.5d cap at 5e-3); `n_verify=0` refused. Certified variant: two-grid RAISE gate
  on the quoted subwindow — higher fitted orders are tail-absorbers, full-window
  gates unpassable; claim cap = measured node capacity.
- transport step control: real-pole-only step control is SILENTLY WRONG with tight
  arb radii — complex-pole clearance + geometric tail-bound halving.
- Sparse-Taylor-lattice acceptance: zero-run-aware (`_tail_bound_nz`) —
  `trunc_worst` is trustworthy on monomial-cone/sparse-support legs. Advisory: the
  3/4 ratio clamp soundness argument needs ratio≤0.5 (all in-repo callers comply);
  `two_sector` near-resonant slots just below the SVD trip carry silent
  amplification — bank per-branch worst smin/smax; two-precision gates cover.
- de_load fast path: unreduced RF denominators may declare APPARENT singular points
  (conservative, safe).
- On-shell derivative ops for DE assembly: use Oₖ=pₖ·∂pₖ and assert O[pᵢ²]=0 for
  EVERY leg incl. conserved.

RELATED: `boundary_branches` delegates to tools/frobenius-boundary (var=∞), which in
turn points at `epslimit.richardson_boundary` for disk-edge series sums;
frob_scalar pulls the annihilator tool; ratfit backs exact entry fits; the
sampler (sampler.py) needs gatekeeper + counterweight/canonical_form; the eMPL
evaluator tools/ellipticus is the same contract family (builds its own systems instead of
loading saved DEs).

CREDIT: differential equations for master integrals after Kotikov [Kotikov], Remiddi
[Rem] and Gehrmann–Remiddi [GR], in Henn's ε-form [Henn] where available; transport by
generalized series expansions around singular points in the line of Lee, Smirnov and
Smirnov [LSS], Moriello [Moriello], Hidding's DiffExp [DiffExp],
Armadillo–Bonciani–Devoto–Rana–Vicini's SeaSyde [SeaSyde] and Liu–Ma's AMFlow [AMFlow],
whose vacuum-boundary conventions (Liu, Ma, Wang [AMF1]; AMFlow.cpp contributors
[AMFcpp]; see the AMFlow.cpp fork, the sibling repository amflow-cpp) are ported
verbatim; local theory after Wasow and Moser/Turrittin (classical), Sylvester solves by
Bartels–Stewart [BS]; Ising-class test operators from Assis et al. [Ising11]. The
epslimit member is Richardson extrapolation (L. F. Richardson 1911; Richardson and Gaunt
1927 [Richardson]) on ε-grid inputs from AMFlow [AMFlow] and pySecDec [pySecDec], with
arithmetic by mpmath (F. Johansson and contributors); the flintexport member's
arithmetic is all FLINT's fmpq_mpoly (multivariate polynomials by D. Schultz; FLINT by
W. Hart, F. Johansson, A. Ahlbäck and the FLINT developers; flintlib.org) through
python-flint (F. Johansson, O. Benjamin and contributors) — that module is a thin walker
over their work. Bracketed keys resolve in REFERENCES.md at the repository root.
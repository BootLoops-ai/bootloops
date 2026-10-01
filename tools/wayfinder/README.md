# Wayfinder — anchor transport and connection fitting for saved DE systems

Anchor-transport package. Loads saved differential-
equation systems (three on-disk formats) into ONE evaluator contract, then
transports boundary/anchor data along the DE: fixed-eps Taylor endpoint
transport, Frobenius landing at regular-singular points, Cauchy eps-fan
Laurent extraction, AMFlow-convention vacuum boundary seeds, two-precision
digit gates, and sha256 input manifests. Two further members sit beside the
DE route: `epslimit`, the engine-agnostic eps->0 limit extractor for an
external engine's (eps, value) grid (plus Richardson summation of a series
at its radius of convergence), and `flintexport`, the FLINT exact
rational-function walker that writes the monomial-json connections
`de_load` reads. Both are lazy modules (`from wayfinder import epslimit,
flintexport`); see the two member sections near the end.

```
from wayfinder import (load_monomial_json, load_amatrix_json,
                         load_kira_targets,                    # + DESystem.enable_fast_path
                         transport_fixed_eps,
                         frobenius_basis, land,                # + frob_scalar (module)
                         cauchy_laurent, vandermonde_laurent, anchor_roundtrip,
                         tadpole, singlemass_vacuum, vacuum_ending_seed,
                         gate_strings, feed_gate_table,
                         write_manifest, check_manifest,
                         RF, two_sector_series, eval_two_sector,
                         quad_refine)                          # + QuadCert, QuadNonConvergence
from wayfinder import frob_scalar, boundary_branches         # lazy modules
from wayfinder import epslimit, flintexport                  # lazy members (flintexport needs python-flint)
```

Module map:

| module | role | headline caveat |
|---|---|---|
| `de_load` | 3 saved formats -> `DESystem`; `enable_fast_path(eps, wp)` = exact-rational `.A_series` + `.singular_points` | kira grammar edges raise; apparent sings from unreduced RF denominators |
| `transport` | fixed-eps Taylor march (endpoint / eta-march), NO Laurent window | zero-run-aware step acceptance (see the G1 bullet) |
| `frobenius` | matrix Frobenius basis + LS landing; integer-gap Moser/Turrittin shearing | REAL eps0 only (Im-linear flag) |
| `frob_scalar` | SCALAR resonant integer towers, value-at-singularity (validated 260d two-precision) | lazy sympy + annihilator.py |
| `boundary_branches` | delegator to `tools/frobenius-boundary` (var = infinity) | eig(D1) caveat verbatim |
| `two_sector` | two-branch series at apparent higher poles | fixed-eps; k=1 autopin; NL truncation |
| `epsfan` | Cauchy fan R=1e-3 + real-grid Vandermonde + anchor round-trip | full circle unless `schwarz_real` (opt-in) |
| `ratfun` | fixed-eps RF layer; backs the fast path | num/den unreduced |
| `boundary` | AMFlow-convention tadpole/vacuum seeds (verbatim ports) | ending IBP not implemented |
| `gate` | two-precision digit gate + per-literal FEED table | feed it MEASURED numbers |
| `manifest` | sha256 input manifests | bytes+sha256, never mtime |
| `quad` | E3 refine-until-gate adaptive quadrature (`quad_refine`) | agreement is a two-depth estimate, not a proof — guard digits absorb the gap; GL nodes expensive, opt-in |
| `epslimit` | eps->0 Laurent extraction from an EXTERNAL engine's (eps, value) grid (auto pole order, log(eps) resonance by leave-one-out, cross-method digit estimate) + `richardson_boundary` disk-edge series summation | digits are capped by the grid, not the algorithm (3 samples ~6-8 d); sets/restores global mp.dps around the solve; not `epsfan` |
| `flintexport` | FLINT `fmpq_mpoly` fraction walker: bulk exact rational-function algebra for connection assembly; exporter of the monomial-json format | requires python-flint + sympy (lazy member); integer exponents only (ValueError otherwise) |

Pure `mpmath` + stdlib. `python-flint` is an OPTIONAL accelerator behind
try/import (`HAVE_FLINT`); nothing requires it. No sympy anywhere in the
matrix/transport paths; the two exceptions are LAZY and outside package
import: `frob_scalar` pulls sympy + `tools/annihilator/annihilator.py` inside its entry
points (exact indicial-root factoring), and `boundary_branches` delegates to `tools/frobenius-boundary`
(sympy) behind an import guard. `flintexport` is the one member that
REQUIRES python-flint (and sympy for its input trees); it is a lazy module,
so package import stays flint/sympy-free and its battery leg skips by name
without them. `epslimit` is pure mpmath, also lazy.

---

## The shared contract

```
class DESystem:
    n: int                       # system size
    var: str                     # transport variable ('s', 'y', 'eta', ...)
    A(x, eps, dps) -> n x n list[list[mpc]]
    meta: dict                   # {'source', 'sha256', 'format', ...}
```

Every consumer (`transport`, `frobenius`) takes any object with this shape —
the loaders produce it, and inline systems (see
`tests/control_hypergeom.py`) satisfy it in a dozen lines. Two optional
duck-typed extensions are honored:

* `.singular_points` (or `meta['singular_points']`): list of singular
  x-points. Transport works without it (two-radius cross-check discipline)
  but is faster and alias-provably-safe with it. The loaders do not populate
  it at load time; the opt-in `DESystem.enable_fast_path(eps, wp)` attaches
  it (denominator roots, deduped — may include APPARENT roots, see Known
  gaps), or set it on the instance yourself when you know the poles.
* `.A_series(x, eps, dps, M) -> [A_0..A_M]` Taylor coefficients of A about
  x: fast path replacing the generic per-step Cauchy circle (~N A-evals per
  step). Strongly recommended for big loaded systems (a 79x79 production
  matrix measures 0.77 s/A-eval at dps 50 — generic-path transport is
  minutes/step without this hook).

`DESystem.entry_exact(i, j, x, eps)` is a documented debug/test hook
(exact Gaussian-rational entry value) — NOT contract API; do not build
production code on it.

## Modules

### de_load — saved formats -> DESystem

| loader | format | pinned against |
|---|---|---|
| `load_monomial_json(path, dps_check=50, var=None, closure_json=None)` | exact rational monomials `'i,j' -> {num, den}` in (eps, var); already in eps | exporter: `flintexport.py` (this package, `wayfinder.flintexport`) |
| `load_amatrix_json(path, n=None, var="y")` | rational-function-of-d y-polynomial fits, `entries['i,j'].{N,Q}`; file in d, evaluated at d=4-2eps | a 79x79 production `A_of_y.json` fit file |
| `load_kira_targets(path, masters_path=None, family_yaml=None, var="eta")` | AMFlow `kira_target.m` (eta,d) reduction rules -> assembled eta-DE | AMFlow-produced `kira_target.m` files (a two-loop five-point family); convention check: tadpole row = (d-2)/(2 eta) = exact eta-scaling of the massive tadpole |

Entries evaluate EXACTLY (`Fraction` / flint `fmpq` / Gaussian rationals for
complex points) and round to dps only at the end inside `mp.workdps`.
Unstored entries are exact zeros. Basis order for kira systems = the
`masters.final` file order = AMFlow's boundary-vector order (boundary
builders MUST keep that convention). `load_amatrix_json` refuses files with
a non-empty `failures` list (incomplete matrix); n is inferred as 1+max
index (A_of_y.json does not store n — pass `n=` if your boundary says
otherwise). Unhandled kira constructs (cut propagators, non-trailing eta,
symbolic invariants, multi-family files, missing dotted targets, ...) raise
`NotImplementedError` naming the unsupported construct — partial-but-correct by
design, never guessy.

### transport — `transport_fixed_eps(desys, eps0, x0, x1, y0, dps, path=None, ...)`

Scalar-mpc Taylor march of y' = A(x, eps0) y at FIXED eps0 (complex eps0
allowed):
adaptive step = ratio x dist-to-singularity, halve on non-convergence,
exact final landing, landing-regularity assert, 400-step/leg cap. Local
A-Taylor coefficients come from `.A_series` if present, else a Cauchy
circle (radix-2 FFT) with either a declared-singularity radius bound or a
two-radius (r, r/2) cross-check. Near-singular straight paths raise (Zeno
abort / step-fraction assert) with the detour fix named — pass complex
waypoints via `path=[...]` (e.g. `[x0/4, 0.03j, x1]`).

Contract-wording deviation (deliberate, documented at `_MIN_HFRAC`): the
"assert if h/|x1-x0| < 1/64" rule is enforced on the step FRACTION of the
distance-to-singularity proxy (the chassis invariant), not on absolute |h| —
the literal reading would forbid legitimate Frobenius-landing approach legs
where |h| scales with the distance to the singular point.

Guard rules:

* **Relative convergence guard.** Both the early-break and
  the step-acceptance test are relative to the state scale
  `||y||_inf` (`tail < guard * max(||y||, tiny)`). An absolute guard is
  accidentally relative only for a `W0=[I|0]` state with norm >= 1;
  for small-norm states (eps-fan seeds ~ eps^k, eta^lambda
  Tradition seeds, Frobenius feeds) it silently loses
  ~log10(1/||y0||) digits (measured at dps=65: y0=1e-60 -> 12.8 digits;
  under the relative guard all of y0 = 1, 1e-30, 1e-60, 1e-90 give >= 66 relative
  digits — regression `test_small_norm_y0_relative_guard`).
* **Zeno threshold scales with REMAINING leg distance (F5)** — a
  leg-length-relative threshold would spuriously abort legitimate approaches on
  long legs (e.g. `land()` launched far from x_sing with a second declared
  singularity capping the match radius).
* **Short `A_series` returns raise (F6)** — fewer than M+1 orders is a
  ValueError, never a silent zero-pad.
* **Honesty diagnostics returned (F7)** — `return_diag=True` returns
  `(y, {'steps', 'trunc_worst', 'trunc_worst_log10', 'min_hfrac'})`; callers should record these too.
* **Geometric tail-bound acceptance (G1), zero-run-aware.**
  Step acceptance uses the geometric tail rule
  generalized to sparse Taylor lattices: tail bounded
  by `(t_a+t_b)/(1 - min(3/4, (t_b/t_a)^(1/(b-a))))` over the last two
  ABOVE-FLOOR term magnitudes at indices `a < b` (floor = `10^-(wp-10)` x
  running max term — Cauchy-circle roundoff sits ~`10^-wp` relative, so
  genuine threshold-scale terms sit 15+ orders above it), never a single
  small term and never raw consecutive terms that may sit on the
  zero/roundoff lattice. A complex pole pair symmetric about the path
  zeroes alternate Taylor terms (dip-then-climb), so a single-term
  early-break silently loses ~55 digits on that geometry; the pole
  clearance is COMPLEX-safe (min |z-s| over
  declared sings — Euclidean clearance over ALL denominator roots).
  `trunc_worst` banks the geometric bound. Regressions:
  `test_complex_pole_pair_off_axis` (closed form + a provable wrong-rule
  failure demonstration), `test_sparse_lattice_monomials` (below).
  **Zero-run-aware step acceptance (the sparse-lattice rule).**
  A sliding `m > 8` early break on the raw
  last two terms has `_tail_bound(0, 0) = 0`, so a local Taylor lattice
  with gaps >= 2 orders (monomial-type A stepped from a symmetry/anchor
  point — exactly the de_load monomial-cone shape) zeroes two consecutive
  terms and ACCEPTS a catastrophically truncated step while banking a
  false `trunc_worst` (measured on the raw rule: `y' = x^2 y` from x=0 at
  dps 60 -> 6 correct digits, recorded bound 9.9e-69, ~63 orders false).
  The rule here
  (`_tail_bound_nz` + freshness gate): bound from the last two ABOVE-FLOOR
  terms with the gap-normalized ratio; early break only while that pair is
  FRESH (trailing zero-run since `b` strictly shorter than the largest
  observed inter-term gap — dense series reduce exactly to the raw
  behavior); an all-zero/roundoff trailing window NEVER accepts early and
  falls through to the fixed top order `mtay`.
  At the top order the zero-run-aware bound is used whenever
  the pair is current through the top (closing the raw-top-two
  lattice-misalignment hole too); the `t_prev==0 -> r=1/2` fallback
  applies ONLY there (A==0 / terminated polynomial series stay correct).
  Regression evidence (`test_sparse_lattice_monomials`):
  `y'=x^k y` vs `exp(x^(k+1)/(k+1))`, k in
  {1,2,3,5,10,12} x dps {60,80} -> 61-82 digits (bars 54/74), certificate
  honest at all 12 (rel err <= 10 x `trunc_worst` + the `10^-(dps-2)`
  output-rounding floor); the raw-rule positive control reproduces 6d/5d at
  k=2/5. `trunc_worst` is trustworthy on sparse-support legs again.
  Advisory: the 3/4-clamp soundness argument assumes
  `ratio <= 0.5`; nothing asserts `ratio <= 0.75` — keep callers at <= 0.5
  (all in-repo callers are).

Fast path: `desys.enable_fast_path(eps, wp=250)` (any de_load DESystem)
attaches `.A_series` (exact-rational fixed-eps `ratfun.RF` Taylor — no
Cauchy-circle aliasing, no two-radius cost) and `.singular_points`
(denominator roots) (see `de_load.py`). Plain
`.A()` systems behave exactly as before; `A_series` refuses wrong-eps
(beyond 10^-dps, exact-rational compare) and dps > wp calls. Caveat: RF
num/den are unreduced, so `singular_points` may include APPARENT roots —
conservative for stepping; prune manually if a landing point is flagged.

### frobenius — `frobenius_basis(desys, eps0, x_sing, dps, kmax, direction=+1)`, `land(...)`

Full matrix-Frobenius fundamental basis Y(u) = P(u) exp(M log u) at a
regular-singular point (Sylvester recursion for P_q; Jordan/equal-eigenvalue
log layers included), plus least-squares landing of a transported vector
onto that basis at three matching points (overdetermined; residual
REPORTED, never silently trusted). Double pole at
x_sing -> ValueError (not regular-singular).

Nonzero-INTEGER indicial gaps are
removed by a Moser/Turrittin SHEARING chain. The top cluster
of the gapped congruence class is split off via its resolvent-contour
spectral projector (trapezoid circle at sep/2; idempotency/commutation/
rank/trace quality gates), and the gauge Y = T diag(u I_G, I) Z lowers that
eigengroup by exactly 1; repeat until gap-free, then the Sylvester
recursion applies. Basis dict gains `sheared`, `n_shears`, `shear_factors`,
`indicial_reduced` (`indicial` stays = eigenvalues of the ORIGINAL residue;
`M`/`P` are the reduced system's). Gap-free systems take the IDENTICAL old
path. Verified against four hand-derived 2x2 closed forms (gap 1 diagonal /
resonant-log / conjugated, gap 2 double-shear) plus `land()` through the
shearing path in `tests/test_frobenius.py`. Any shear quality-gate failure
raises RuntimeError naming the alternatives (frob_scalar / driver-side
landing / frobenius_boundary).

The LS solve records `solver` ('qr' vs
'normal_eq') and a svd condition-number estimate `cond_est` of the matching
matrix in the returned dict. If `mp.qr_solve` fails (narrow except — 
programming errors propagate), `land()` FAILS LOUDLY by default: the
normal-equations fallback squares the condition number of a matrix whose
columns scale as u^lambda_i, and a small residual does not bound the kappa
error in near-degenerate directions. Opt in with `allow_normal_eq=True`.

### frob_scalar — `transport_value(rec_json, wfr, dps, s0=1/2, s_check=2/5)` (module import: `from wayfinder import frob_scalar`)

SCALAR Frobenius transport-to-singular-point: value of a holonomic power
series AT the regular-singular point x=1 of its ODE, from an exact
P-recurrence for its Taylor coefficients (annihilator-fit format).
Validated at
260-digit two-precision on a resonant-tower production constant plus
two oracle controls. Handles the FULL resonant integer towers
the matrix module shears away (per-root eps-jet basis, e.g. indicials
{0..6} u {3/2,3/2}); this is the scalar-route answer to integer gaps.
Two design warnings are documented in the module:
(1) sequence-recurrence fits do NOT annihilate the generating function —
the exact boundary polynomial q must be built and `D^{deg q+1} o L`
composed (done automatically in `transport_value`); (2) the Frobenius
eps-trick needs the PER-ROOT derivative window `[P_rho, P_rho+mult)` else
degenerate non-log resonances silently drop basis elements. Gates
reported: independent-check-point residual, divergent-branch coefficients,
resonance-drop residuals, series tail. Lazily imports sympy +
`tools/annihilator/annihilator.py` inside its entry points (package import stays
sympy-free). Fast hand-derived ports of both control classes:
`tests/test_frob_scalar.py` (dps <= 120; NOT the full 260d production run).

### boundary_branches — thin delegator to `tools/frobenius-boundary` (var = INFINITY)

Import-guarded, lazy: delegates `build_poly_DE / spectrum / classify /
exclude_strata / branch_series / phi_matrix / seed_vectors / eps_to_d`
straight to `frobenius_boundary.core` — NOTHING vendored (that package has
its own 43-check selftest: `cd tools/frobenius-boundary && python3 -m frobenius_boundary.selftest`). Routing: `frobenius` = landing at a FINITE
regular-singular point; `frob_scalar` = scalar operators with resonant
towers; `frobenius_boundary` = branch CLASSIFICATION at infinity for
cut-block DEs. The wrapper docstring carries the FLAGGED
caveat verbatim: "exponents = eig(D1)" is exact only when the N-frame is
Fuchsian — with a genuine nilpotent D_0 block the TRUE exponent multiset is
the eig of the shear-reduced residue (classify/exclusion still run on
eig(D1), regression-pinned), plus the `kill_integer_eps_indep` vacuum-
theorem precondition. dps hygiene: the delegate raises global mp.dps at
import; the wrapper imports it inside `mp.workdps(>=50)` so ambient
`mp.dps` is NEVER mutated (tripwired in `tests/test_boundary_branches.py`).

### epsfan — `cauchy_laurent(f, kmin, kmax, dps, R='1e-3', nodes=32, check='R2', schwarz_real=False)`, `vandermonde_laurent(...)`, `anchor_roundtrip(...)`

Laurent coefficients of f(eps) around eps=0 on a Cauchy circle, with a
SECOND extraction (2R, or half nodes) and per-k measured agreement digits.
`diag['min_agree_digits']` is the honest digit count — trust it, not dps.
Working precision auto-covers the R^-k pole leverage.

`schwarz_real=True` (OPT-IN) enables the conjugate-mirror
shortcut: only the upper half circle is evaluated (17 of 32 nodes) and the
rest is filled by conjugation. **Valid ONLY for Schwarz-real rows** (ALL
Laurent coefficients real). On a complex-Laurent row the mirror silently
destroys the imaginary content, the corruption at the affected order is
R-independent so the two-pass diag does NOT reliably catch it there, and
the measured half-circle rescue (17-node LS) was truncation-limited to ~13
digits — catastrophically worse, not degraded (measured). The default
`schwarz_real=False` = full circle, always safe. `diag` reports
`schwarz_real` and the measured `n_evals`.

`vandermonde_laurent(f, kmin, kmax, dps, eps_max='1e-3', span='10',
n_verify=7, eps_nodes=None)` is the SECOND validated
extraction mode:
REAL nodes eps=1/q on a geometric grid with total span ~10 (one decade,
per-step ratio ~1.105) and eps_max ~1e-3, square Vandermonde solve on the
P=window fit nodes, INTERLEAVED held-out verify spares.
`diag['self_digits']` is the MEASURED held-out extraction error — the
honest count. `n_verify=0` is refused (unmeasured = unfabricatable).
Working precision covers the measured solve loss (~3.6*P digits,
P-dominated, window-offset independent — calibrated).

`anchor_roundtrip(f, coeffs, dps, eps_anchor='6.1e-4')` resums extracted
coefficients at a fresh anchor eps and compares against a direct
evaluation (callable f, or a precomputed transported value). Both sides
carry the same upstream boundary content — the boundary CANCELS, so the
residual isolates extraction error (plus window truncation unless f is
window-limited). Unit-tested on a known Laurent: clean coefficients
round-trip at ~dps; a 1e-20-relative coefficient corruption is measured at
~20 digits.

### ratfun — fixed-eps rational-function layer (`RF`)

The polynomial/RF layer (gates passed at full stored ~70 digits through
it).
mpc-coefficient num/den polynomials at FIXED eps: `taylor` (shifted-Horner
+ series division), `eval`, `series_inf`, `laurent0`. Backs
`DESystem.enable_fast_path`; import directly for driver-style use. All
calls inside `mp.workdps` (caller-owned).

### two_sector — `two_sector_series(desys, eps0, x_sing, seed_a, seed_b, NL, dps, pmax=4, pins_a/pins_b=None, autopin=False)`, `eval_two_sector(result, s, dps=None)`

Two-branch Frobenius series `y ~ sum a_n s^n + (-s)^(-eps) sum b_n s^n` at
a point with APPARENT higher poles (entries to `s^-pmax`), fixed eps.
Implements: `scc_blocks` +
`block_delays` (Bellman `phi_i = max(0, max_j(phi_j + p_ij -
1))`, positive pole cycles refused), the slot schedule `tau = n + phi_i`
with small dense slot solves, level<=0 constraint rows REPORTED as
residuals (never fudged), and the k=1 autopin (resonant coordinate
determined from the constraint rows via particular+homogeneous runs,
`t = -r_p/r_h`). Deliberately omitted: the eps-Laurent WINDOW
arithmetic and eps-offset resonance pivoting — package is fixed-eps by
design; resonance is SVD-detected per slot and must be pinned (loud raise
otherwise). Entry Laurent tables come from one Cauchy circle on the
pointwise `.A` contract with an on-circle pole-order significance test and
a `pole order > pmax` tripwire. References beyond NL truncate to zero: the
top ~pmax levels of downstream components carry truncation-order error
— keep NL above what you consume. Seeds
are the caller's (region-derived). Unit test: a
gauge-constructed closed form with an apparent `s^-2`, a structurally
resonant B-slot and a self-verifying DE check
(`tests/test_two_sector.py`).

### quad — `quad_refine(f, interval, dps, guard=10, method="tanh-sinh", depth0=None, max_depth=None, wp_extra=20, scale=1, full_output=False)`

THE standardized adaptive-quadrature refine-until-gate helper. Integrates
f over `[a, (waypoints...,) b]` (finite complex endpoints; mp.inf/mp.ninf
supported) on mpmath's nested tanh-sinh levels (Gauss-Legendre opt-in) and
accepts ONLY when the double-refinement agreement `|I_{d+1} - I_d|` beats
`10^-(dps+guard) * scale`; escalation is exact continuation of the nested
rule; non-convergence at the cap raises `QuadNonConvergence` NAMING the
measured agreement, tol, depth, cap, wp, method and segment — a value the
loop did not certify is never returned. `depth0`/`max_depth`/`wp_extra` are
STARTING seeds / resource caps only (speed, never the answer; the seed
scales with dps via mpmath `guess_degree`). `full_output=True` returns a
`QuadCert` (value, gated raw agreement, `err_bound` covering the returned
dps-rounded value = agreement + measured dps rounding + wp roundoff floor,
depths, evals, per-segment breakdown, printable `bound_line()`). Endpoint
singularities (log-class) and oscillatory integrands are unit-tested at dps
30/60/100 against in-test closed forms; fail-closed cap/divergent/nan
controls, seed-invariance, dps 30-vs-70 crank and the import-dps check live
in `tests/test_quad.py`. Calibration (measured): default runs
accept at the seeded depth (8/9/9 at dps 30/60/100), worst raw agreement
10^-(dps+34) vs tol 10^-(dps+10) — ~24 orders of headroom; this toolkit's
default gate runs should never escalate. Wire it by IMPORT from wave
scripts; this module touches no row script itself.

### boundary — AMFlow-convention boundary values

`tadpole(m2, nu, d, dps)`, `vacuum_known(loops, props)`,
`singlemass_vacuum(loops, props, eps, dps)` (builtin table (1,1) (2,3)
(3,4) (3,5) (4,5), ported VERBATIM from amflow-cpp
(`src/qft/vacuum.cpp` in the sibling repository amflow-cpp), mirroring AMFlow.m lines 819-824),
`single_mass_prefactor(n, m0, m_eps_coef, eps, dps)` (amfsystem.cpp ~line
2118), `vacuum_ending_seed(desys, scheme, params, dps)` for
Tradition/SingleMass endings (leading eta-power constants only; the DE owns
all subleading orders). Normalization: measure `d^d k/(i pi^(d/2))` per
loop, propagator `1/(k^2 - m^2 + i0)^nu`, d = 4-2eps, NO exp(eps*gammaE).
Exact Gamma poles raise with the nonzero-eps-grid fix named. Ending
masters outside the builtin vacuum table need the ending IBP reduction ->
`NotImplementedError` (the message names what is not supported).

### gate — `gate_strings(computed, vendored, dps_pair, zero_scale=None, zero_digits=None)`, `write_gate_report(report, out_json)`, `feed_gate_table(entries, digit_bar=30)`

Two-precision digit gate: PASS iff matched digits vs the FULL vendored
string >= stored-1 AND the two computed precisions agree to >= stored
digits (a value that moves with working precision is noise, not a result).
Pure function; `write_gate_report` is the module's only writer.

`feed_gate_table` formats the per-literal FEED gate table
— one row per literal: name, stored digits,
recomputed-vs-stored matched digits, two-precision agreement, verdict.
Pure formatting + charter-bar thresholds (PASS iff stored >= 30, matched >=
stored-1, pair-agree >= stored; unmeasured literals FAIL; empty table
FAILs). Feed it MEASURED numbers — `gate_strings` per_key entries compose
directly (see tests/test_gate.py F6). Writes nothing; returns the table
string + per-literal report.

Degenerate-reference rules (the relative-digit metric is meaningless
on degenerate references and must never pass vacuously):

* **Zero-vendored strings gate ABSOLUTELY.** A vendored '0'/'0.0' (the
  common vanishing-imaginary-part reference) passes only if BOTH pair
  members satisfy `|computed| <= scale_ref * 10^-zero_digits`; scale_ref
  defaults to the max |value| over nonzero sibling keys (override with
  `zero_scale=`, scalar or per-key dict), zero_digits to the max sibling
  digit claim. A lone un-scaled zero HARD-FAILS. (Without this rule,
  computed ('5.0','5.0') vs vendored '0' would PASS.)
* **Nonzero 1-digit vendored strings HARD-FAIL** — `matched >= n_stored-1
  = 0` is vacuous, the key gates nothing; store >= 2 digits.
* **`dps_pair` must be strictly increasing** (suggest dhi >= dlo + 10);
  equal/reversed precisions raise ValueError — a (v, v) caller would
  silently degrade the gate to one precision. Byte-identical pair strings
  are flagged per key (`pair_identical_strings`).

### manifest — `sha256_file(path)`, `write_manifest(paths, out_json)`, `check_manifest(manifest_json)`

sha256 manifests for input artifacts, keyed relative to the manifest file
so they travel with the run directory. `check_manifest`
compares bytes+sha256 (never mtime). Verify manifests before reusing artifacts after any interruption.

---

## Design rules

1. **No Laurent window.** Windowed Laurent endpoint transport
   is UNSOUND pointwise: the connection's
   (d-4)-poles make the Laurent depth unbounded and the window aliases.
   The calibrated replacement is
   fixed-eps scalar transport + a Cauchy eps-fan. That is the ONLY eps
   handling in this package.
2. **R = 1e-3 Cauchy default (32 nodes).** Measured calibration: R = 1e-17
   loses ~50 digits to 1/eps pole leverage; R ~ 1e-3 costs ~3 digits per
   Laurent order and the second-pass (r, delta) diagnostics MEASURE the
   alias/noise floor (never trust a single extraction). The real-grid Vandermonde mode uses the
   matching measured design rule: geometric decade
   grid (total span ~10), eps_max ~1e-3, held-out verify spares.
2b. **Full circle unless Schwarz-real.** The conjugate-mirror shortcut
   (17-of-32 nodes) is valid ONLY for rows with all-real Laurent
   coefficients; complex-Laurent rows need the full circle (measured:
   17-node LS ~13d — catastrophic).
   `schwarz_real` is opt-in and defaults to False.
3. **Real-eps Frobenius matching only.** Complex-eps Frobenius matching carries a known Im-linear
   residual defect, root cause not closed. `frobenius_basis`/`land` raise on complex eps0
   until that flag is lifted. Complex eps is fine in `transport_fixed_eps`
   itself (the eps-fan circle feeds transport, not the Frobenius matcher).
4. **Explicit dps everywhere; global `mp.dps` is never touched.** All
   string parsing happens inside `mp.workdps` (parsing at ambient precision silently truncates). Every module obeys this; the controls trip-wire it.
5. **Partial-but-correct, never guessy.** Anything not pinned against real
   producing/consuming code raises `NotImplementedError` naming the
   unsupported construct (kira grammar edge cases, partial-eta endings, ending IBP
   reduction); integer indicial gaps are handled by the pinned shearing
   implementation (hand-derived closed-form tests).
6. **Measured honesty handles.** Fan agreement digits, Frobenius landing
   residuals, gate matched digits, Taylor truncation norms are REPORTED —
   nothing is silently trusted, no digits are fabricated.

## Convention pins

* Vacuum table, SingleMass prefactor and ending-value rules mirror AMFlow.m
  (lines 819-824 / 1039-1061; the amflow-cpp fork lives in its own repository,
  the sibling repository amflow-cpp).
* T_def convention: vendored 60d Laurent coefficients re-derived in
  tests/test_boundary.py.
* Real-grid Vandermonde node design: geometric decade q-grid 1013..10037,
  eps_max ~1e-3, interleaved verify spares, digits = min(oracle,
  self-consistency) — `vandermonde_laurent`.
* Conjugate-mirror-is-Schwarz-only + the measured 17-node-LS ~13d failure
  mode sit behind the `schwarz_real` opt-in default.
* MUM-crossing Zeno detour rule — a straight path through/near a singular
  point collapses the adaptive step; use an explicit complex detour.

## Tests and controls

```
python3 -m pytest tests/ -q            # 71 unit tests (~810 s under load; the 14-test quad suite alone ~8 s)
python3 tests/control_hypergeom.py     # end-to-end positive control, ~6 min
python3 tests/control_tadpole.py       # boundary vs direct Gamma, seconds
python3 tests/test_quad.py             # quadrature self-test (also under pytest), ~8 s
python3 tests/test_epslimit_synthetic.py; python3 tests/test_epslimit_fixtures.py   # epslimit legs (also under pytest), ~4 s
python3 tests/test_flintexport.py      # flintexport exact identities (also under pytest; skips without python-flint), <1 s
```

(The suite covers: 32 integration tests, the
small-norm-y0 regression, the diag/A_series-guard test, the pytest
wrapper for the gate degenerate-reference controls, vandermonde realgrid,
anchor round-trip,
schwarz_real mirror (incl. the complex-row corruption must-fail), the
FEED-table controls inside the gate suite, a test_epsfan pytest
enforcement wrapper — check() failures in test_epsfan must be
visible to pytest, not only the __main__ exit code — plus the
Frobenius-route (shearing/frob_scalar/boundary_branches), fastpath and
two_sector suites listed below, and the sparse-lattice regression
`test_sparse_lattice_monomials`.)

Frobenius-route suites: 
`tests/test_frob_scalar.py` — fast hand-derived ports of the two scalar-transport
control CLASSES (3F2 half-integer vs closed form pi/Gamma(3/4)^4;
log-resonant class with the D^2 o L boundary composition asserted, vs a
two-precision tanh-sinh quad oracle; recurrences verified EXACTLY on every
term before transport; dps <= 120, ~55 s); `tests/test_frobenius.py` gained
the 4 shearing closed-form checks + land()-through-shearing (~35 s extra);
`tests/test_boundary_branches.py` — delegator contract (lazy import, global
mp.dps never mutated, synthetic 4-eps classify/exclude round-trip; skips if
frobenius_boundary is absent).

`control_hypergeom.py` runs the FULL pipeline on 2F1(eps,-2eps;1-3eps;x)
(inline DESystem -> fixed-eps transport 1/2 -> -0.6i -> -1 around the x=0
singularity -> gate vs held-out `mpmath.hyp2f1` at 90 stored digits at dps
100 AND 140 -> eps-fan -> Laurent coefficients vs an independent
`mpmath.taylor` expansion, with the exact anchor c2 = pi^2/6).
Control results (measured, integration run):

* unit tests: 32/32 pass;
* tadpole control: PASS — 151.5/151.9/152.4 matched digits vs the direct
  Gamma oracle at the 150d bar on all 3 (m2, nu, d) triples, two-precision
  stable;
* hypergeometric fixed-eps transport (integration probes, measured):
  generic path dps=100 -> 101.4 matched digits vs held-out
  mpmath.hyp2f1(0.1,-0.2;0.7;-1); A_series fast path dps=105 -> 106.6d,
  dps=145 -> 147.5d. The full packaged control (gate + eps-fan + taylor
  comparison) is `tests/control_hypergeom.py`; run it and quote only its
  own measured verdict — do NOT quote numbers that are not in its log.

## Known limitations

* `DESystem.enable_fast_path(eps, wp)` exports `singular_points` and the
  exact-rational `A_series`
  (ratfun RF layer; opt-in per eps, all
  three loader formats). Caveats: unreduced RF denominators can declare
  APPARENT singular points; the numpy root fallback for high-degree /
  repeated-factor denominators is float64 (fine
  for clearance steering, the transport tail bound still gates steps).
* `two_sector` autopin implements k=1 pins only;
  >1 resonant coordinates per branch raise. Cross-level truncation: series
  references beyond NL are dropped (truncation-order effect on the top
  ~pmax levels).
* `load_kira_targets`: cut families, symbolic kinematic invariants,
  multi-family files raise NotImplementedError (lift when such a file shows
  up).
* `vacuum_ending_seed`: partial-eta insertion modes and non-builtin ending
  masters (ending IBP reduction), SingleMass FactorizeFamily derivation of
  gamma_params/sub_values — all NotImplementedError with the missing piece named.
* `frobenius`: nonzero-integer indicial gaps are handled by the
  Moser/Turrittin shearing chain (see the module section above); the
  RuntimeError quality gates name the frob_scalar / driver-side /
  frobenius_boundary alternatives if a shear split ever fails.
* Complex-eps Frobenius matching blocked on the upstream Im-linear flag.

## The sampler — the dlog-connection-once route

Fits the constant ℚ dlog-connection matrices `{A_a}` of a canonical
differential equation

    d J  =  ε · ( Σ_a  A_a · d log a ) · J

from cheap **numeric** IBP samples, never running symbolic multivariate
FireFly reconstruction.

### Why

Symbolic IBP over ≳9 kinematic variables reconstructs every
entry of `A(x,ε)` as a multivariate rational function — ~30 min/node with
FireFly.  But in the canonical basis the unknowns are just `n_letters × M²`
small rationals.  So: alphabet (Landau Alphabet) + K cheap numeric points
(cached amflow-cpp evaluations) + one linear solve + PSLQ → done.

### Usage

```bash
python sampler.py targets/box1l_conn.json \
    --n-fit 3 --n-heldout 2 --goal-digits 40 --eps-power 3 \
    --parallel 6 --out out/box1l/connection.json
```

or programmatically

```python
from wayfinder.sampler import connection_once
conn = connection_once("targets/box1l_conn.json")
```

The output `connection.json` is the schema `tools/ansatzer/`
consumes via `--connection` (sparse `A: {letter: [[i,j,"p/q"],...]}`).

### Target schema

Same `family` / `masters` block the pipeline already uses, plus

| key        | meaning |
|------------|---------|
| `kinvars`  | list of kinematic invariants to differentiate in |
| `alphabet` | list of sympy strings (letters); falls back to Landau Alphabet on `graph` |
| `rotation` | `{tag: {"coeff": "<sympy in kinvars,eps>"}}` — the diagonal UT rotation.  May be polynomial in `eps` (e.g. `eps*(1-2*eps)` for massless bubbles).  Defaults to `canonical_form.ut_rotation`. |
| `points`   | optional explicit kinematic points (else auto-generated) |

### Backend status

| backend | status |
|---|---|
| `finite_diff` | **working** — central difference of `mode:"solve_integrals"` Laurent values at exact-rational shifted points; ~⅔·`goal_digits` effective derivative digits. |
| `amflow_diffeq` | **working**.  `libp_deriv` carries the upstream-MMA momentum-derivative chain rule (`LIBPDerivivative`, `Kira/interface.m`) for replacement-defined invariants (`src/ibp/libp_deriv.cpp`, `build_kinematic_chain`), wired through `mode:"diffeq"` in `src/api/run_json.cpp` and `_sample_de_amflow_diffeq` here.  One Kira IBP per point, exact ∂J/∂x_i as rational functions of `d`.  Validation (off-shell 1L triangle, s1=-1, s2=-18/25, s12=-13/25): nonzero d/ds1,d/ds2,d/ds12 in 6.7 s; bubble diagonals = (d-4)/(2s_i) exact; Euler contraction Σ s_a A_a = (d-6)/2 (tri) with vanishing off-diagonals; eps^-1 pole cancellation in the triangle rows exact at 120-digit precision; Richardson central-difference cross-check vs `mode:solve_integrals` at s2,s12 ± {1/1000, 1/2000} agrees to 11-12 digits on every eps order (plain h=1/1000 central diff gives the h²-limited 5-6 d). |

### Honesty fields

`connection.json["honesty"]` always reports `n_points_used`, `cond_number`,
`fit_residual_digits`, `min_heldout_digits`, `failed_to_rationalize`
(`[letter, row, col, float]`), `ut_check` (is the rotated leading order
constant?), and `elliptic_masters_flagged`.

### box1l reference

3 fit + 2 held-out points, 25 amflow calls, ~130 s wall (parallel-6), cached
re-run 0.1 s.  Recovers the Henn 1-loop massless box connection (basis
`ε(1-2ε)·bub_s`, `ε(1-2ε)·bub_t`, `ε²st·box`):

    A[s]   = [[-1,0,0],[0,0,0],[2,0,-1]]
    A[t]   = [[0,0,0],[0,-1,0],[0,2,-1]]
    A[s+t] = [[0,0,0],[0,0,0],[-2,-2,1]]

fit residual 31 d, held-out 32 d, 0 failed.  `Σ_a A_a = -𝟙` (Euler scaling
check).

### Tests

```bash
python tests/test_sampler.py                    # T1 (synth) + T4 (sparse synth) + T5 (tampered-sample negative control), ~3 s
WAYFINDER_RUN_AMFLOW=1 python tests/test_sampler.py   # + T2/T3 (amflow), ~10-20 min
python3 -m pytest tools/wayfinder/tests/test_sampler.py -q   # from the repo root: same legs; T2/T3 skip with their env reason
```

The functions are `test_t1_synthetic` … `test_t5_tampered_sample`, so the
registered wayfinder battery (`python3 -m pytest tools/wayfinder -x -q`)
collects all five.  T5 is the negative control: one derivative value of the
T1 table is perturbed by a relative 1e-6 (19 orders above the 25-digit
rationalization tolerance) and the test asserts that the fit reports
`failed_to_rationalize` entries or the held-out digits fall below 30, while
the same table unperturbed passes in the same process.

## Member — epslimit: eps->0 limits from an external engine's grid

`epslimit.py` (formerly the standalone `eps_extrapolator`). Generic,
engine-agnostic eps->0 Laurent-coefficient extractor: given a list of
`(eps_i, value_i)` pairs — each `value_i` a high-precision numeric sample of a
master integral (or any function) at a chosen `eps_i`, produced by AMFlow,
pySecDec or anything else — it returns the coefficients of

    f(eps) = sum_{k>=k_min} sum_{j>=0} c_{k,j} eps^k log(eps)^j

with the finite part `c_{0,0}` the usual target. It sits DOWNSTREAM of the
engine and does not duplicate engine-side eps=0 resonance fixes.

**Relation to `epsfan`.** `epsfan` extracts the Laurent fan from a DE system
this package can evaluate wherever it likes, so it DESIGNS its eps-nodes
(Cauchy circle R=1e-3, geometric real grids) and certifies alias bounds.
`epslimit` is for the other situation: the grid already exists, chosen by an
engine you do not control, and you want the most digits it honestly carries.
The two share no code.

Why naive solves cap at ~6 digits, and what this handles:

1. *Conditioning.* An eps-Vandermonde with all `|eps| ~ 1e-3` is
   catastrophically ill-conditioned. Power columns are rescaled to
   `(eps/scale)^k`, and the mpmath working precision is auto-set from the
   grid's modulus dynamic range AND a node-clustering condition estimate.
2. *Unknown / negative leading power* (genuine poles): `k_min` is
   auto-detected from the log-log slope of `|value|` vs `|eps|`.
3. *Resonances / Jordan blocks* (a `log(eps)` appears): detected by
   leave-one-out — each candidate log order gets its full power budget and a
   log term is accepted only on a DECISIVE improvement (a true resonance jumps
   from ~7 to ~220 digits; engine round-off wobble does not); handled by the
   generalized `eps^k log(eps)^j` Vandermonde.

| method | when | notes |
|---|---|---|
| `vandermonde_fit` | always | rescaled generalized (power x log) solve; exact (square) or normal-equations (over-determined) |
| `eps_fft` | eps on a circle `\|eps\|=r` | Cauchy/DFT recovery of the full Laurent series; immune to Vandermonde conditioning |
| `neville_zero` / `romberg_ladder` | real/geometric ladder, finite point | Richardson power elimination; condition-robust, matrix-free |
| `bulirsch_stoer_zero` | nearby complex eps-singularities | rational (Pade-type) extrapolation to eps=0 |

`extrapolate()` auto-routes by grid geometry, leading power and resonance
order, runs the applicable methods, and CROSS-VALIDATES them: the reported
achieved-digit figure is method disagreement, not a self-reported residual.

```python
from wayfinder.epslimit import extrapolate, load_amflow_grid, richardson_boundary
samples = load_amflow_grid(['out_eps1.json', 'out_eps2.json', ...])  # AMFlow result JSONs; failure markers dropped
r = extrapolate(samples, max_log=2)
print(r.finite_part)        # c_0
print(r.c(-1))              # 1/eps coefficient, if present
print(r.c(1, 1))            # eps*log(eps) coefficient, if a resonance was found
print(r.k_min, r.log_order, r.achieved_digits, r.diagnostics)

# second verb: sum a series AT its radius of convergence, known tail exponents
S_inf = richardson_boundary(partial_sums, Ms, tail_exponents=None)   # default tails 1/2, 3/2, 5/2, ...
```

CLI: `python3 tools/wayfinder/epslimit.py out1.json out2.json --part re`
(`--kmin`, `--maxlog`, `--dps`), or from `tools/`:
`python3 -m wayfinder.epslimit ...`. Inputs accept python numbers, strings,
Arb-ball strings `'[mid +/- rad]'`, or `{'re':..,'im':..}` dicts. Pure
`mpmath`; the gmpy2 backend is recommended but optional. Unlike every other
module here, `extrapolate()` sets the GLOBAL `mp.dps` to its working value
for the duration of the solve and restores it in a `finally` block.

**Why sparse grids cap low.** A finite-part fit from only 3 eps samples caps
at ~6-8 digits no matter the algorithm: 3 points resolve only
`{eps^0, eps^1, eps^2}`, so the finite part carries an `O(c_3 * prod eps_i)`
truncation error. The cap is information-theoretic, not an algorithm defect
— `test_epslimit_truncation_floor` demonstrates it (3 pts ~8 d, 5 and 8 pts
full precision on an exactly-known function). Against real engine grids the
remaining error is set by the engine's internal working precision.
Operational rule for >=30-digit finite parts: >=6-8 eps samples on a
geometric ladder (ratio 2-4) at engine working precision >= 2-3x the target
digits.

Battery legs (collected by `pytest tools/wayfinder`; ~4 s together):
`tests/test_epslimit_synthetic.py` (6 legs, 12 graded checks against known
exact expansions: clean finite part, pole with `c_{-1}=zeta(3)`, `c_0=pi^2`,
`eps log eps` resonance with the log order auto-detected, eps-FFT on a
circle, narrow 3-point exact quadratic, truncation floor) and
`tests/test_epslimit_fixtures.py` (the file pipeline `load_amflow_grid ->
extrapolate` on the committed fixtures in `tests/fixtures_epslimit/`: a
3-file exact-quadratic grid, `c_0 = 7/2`; an 8-point pole ladder,
`c_{-1} = zeta(3)`, `c_0 = pi^2`, including a planted engine-failure sample
the loader must drop).

## Member — flintexport: bulk exact rational-function algebra on FLINT

`flintexport.py` (formerly the standalone `flint_export`). A fraction walker
over Q(x1..xk) on FLINT `fmpq_mpoly`: DE-connection assembly, reduction-row
contraction, big substitute-and-compare — the jobs where sympy
`cancel`/`together` is the wall (measured: sympy 128 s for ONE 7-entry row;
this walker 285.6 s for the full 22x22 three-variable connection, ~1000x).
It is the exporter behind the monomial-json connection format that
`load_monomial_json` reads.

Method: walk each sympy expression tree ONCE with values held as `(num, den)`
`fmpq_mpoly` pairs (`Frac`), gcd-reducing on add and cross-cancelling on
multiply; leaf substitutions (`d -> 4-2*eps`, kinematic pins) are applied
through the env; sympy supplies only the input tree, ALL arithmetic is FLINT.

```python
from wayfinder.flintexport import mk_ctx, gen, const, poly, make_env, walk, canonical_terms
ctx = mk_ctx(('eps', 's2', 's12'))
env = make_env(ctx, {sym_s2: gen(ctx, 1), sym_s12: gen(ctx, 2),
                     sym_d: poly(ctx, {(0,0,0): 4, (1,0,0): -2})})   # d -> 4 - 2 eps
fr = walk(sympy_expr, env, cache={}, ctx)          # Frac (num, den), lazily reduced
num_terms, den_terms = canonical_terms(fr)         # exact integer-coprime export form (den lead positive); None for exact zero
```

Coefficients as int/Fraction/fmpq; any number of variables; INTEGER exponents
only — a fractional or symbolic power raises `ValueError` (never silently
floors). Not a general CAS. Footgun: keep every operation inside FLINT and
cancel once at export; falling back to sympy arithmetic mid-walk reintroduces
the wall. Requires python-flint (pip) and sympy — hence a LAZY member
(`from wayfinder import flintexport`). Julia sibling:
`tools/counterweight/src/MpolyFeed.jl` (`feed_matrix_from_chunks`).

Battery leg: `tests/test_flintexport.py` — five exact identities (exact zero
after expansion, gcd cancel, a three-variable `(d,s2,s12)` pin identity,
canonical export shape + positive-lead law, non-integer-exponent refusal)
plus a fresh-interpreter check that `import wayfinder` loads neither lazy
member nor sympy; skips by name without python-flint/sympy. Standalone:
`python3 tools/wayfinder/flintexport.py` -> `wayfinder.flintexport selftest: PASS`.

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).

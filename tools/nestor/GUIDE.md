# nestor — guide

Tool page: https://bootloops.ai/tools/nestor.html

KIND: package

PURPOSE: NESTOR (NESted Tanh-sinh ORacle) — one certified-quadrature front door over
the standing layer: node-farmed nested tanh-sinh + wp-scaled per-level precision ladder
+ load-time singularity verification + fail-closed refuse-don't-guess.

USE-WHEN:
- Certified high-precision quadrature ORACLE behind a gate.
- Nested / farmed outer folds; a serial mp.quad leaves the machine idle — farm an explicit
  tanh-sinh grid over a Pool instead.

NOT-FOR: general singular LEFT endpoints at a!=0 need the wants_dists (u,d_a,d_b)
integrand form for full relative accuracy (the a=0 source case is covered). Not a
symbolic integrator; not for integrands whose singularity spec you cannot declare.

INVOKE: `nestor.integrate(f, a, b, dps=, spec=, md=, nproc=)` -> mpf at >=dps agreed
digits OR a NAMED refusal (exit 0/2/3/4). `nestor.oracle(...)` -> full JSON receipt
(receipt_path= writes atomically). `python3 -m nestor.selftest` (13/13: nine
oracle controls + four dispersion legs). `nestor.farm_integrate` /
`nestor.farm_campaign` (checkpointed, resume=skip). `nestor.dispersion` — the
singularity-subtracted dispersion quadrature member (section below).
`nestor.probe` IS quad_probe (delegates by identity to
tools/gatekeeper/probes/quad_probe.py — one source of truth). b=mp.inf allowed
(serial cert-ladder tail).

INPUTS: integrand f + [a,b]; spec = endpoint singularity declaration (type
algebraic|log|none per endpoint) verified at LOAD; dps, md (tanh-sinh depth), nproc.

OUTPUTS: mpf or named refusal (EvalRefusal / SpecRefusal / CertFail / AbortRefusal)
carrying a needed-dps clause; JSON receipt: per-level walls, agreement digits,
singcheck record, refusal record.

ENV: nproc for the farm; no absolute paths in the package.

GATES: validation by parity with the original implementation (2447-node farm grid,
75.9 digits; nested 42.59 digits); selftest 13/13 in-tree (6 positive incl. nested
two-fold + log-endpoint P6 + farm smoke, 3 mutation, 4 dispersion legs D1–D4 —
see the dispersion section), ~3 s.

FOOTGUNS:
- SIGTERM mid-run -> AbortRefusal recorded in the receipt; pool workers reset SIGTERM
  to SIG_DFL so terminate() does not deadlock (do NOT remove _winit's reset — removing it hangs the pool).
- singcheck refuses undeclared log/algebraic growth at an endpoint (slope < -0.05) —
  declare the endpoint type or it fails closed.
- consistency != control: a fixed-parameter dps rerun cannot certify a truncation
  (Nch/grid) bias.
- inverse-sqrt endpoints cap an mpmath tanh-sinh leg at ~wdps/2 (endpoint-node
  cancellation): ANY half-precision plateau is this footgun — substitute the sqrt away
  (t=sin^2(theta) for dt/sqrt(t(1-t)R): analytic integrand, roundoff-limited after,
  <1 s fix).

RELATED: tools/gatekeeper (quad_probe is delegated by identity; held-out
certification downstream); tools/lockpick (PSLQ closure downstream).

CREDIT: tanh-sinh (double-exponential) quadrature is due to H. Takahasi & M. Mori (1974,
Publ. RIMS Kyoto 9:721); its use as the workhorse of very-high-precision integration
follows D. H. Bailey and collaborators (Bailey, Jeyabalan & Li 2005, Exp. Math. 14:317);
the per-level rule nestor farms is mpmath's implementation (F. Johansson and the mpmath
developers). The node-farm, wp-ladder and nesting-pattern consolidation are this
program's.

## dispersion member (`nestor.dispersion`)

KIND: module folder `dispersion/` (`disp_sub.py`, mpmath only; module README beside
it), worked kernel `examples/box1_dilog.py`, battery legs
`tests/test_dispersion_moments.py`, `tests/test_dispersion_surrogate.py`.

PURPOSE: exponentially convergent, singularity-subtracted quadrature for a dispersion
integral I=(1/π)∫ρ(w′)K(w′)dw′ in which the spectral density ρ=Im Σ carries known
threshold non-analyticities (a turn-on (w′−w*)^α[log]^m at the lowest threshold,
finite jumps at higher cuts) and the kernel K is smooth and analytic on each panel —
the shape of every two-particle-reducible self-energy insertion. On a short sub-panel
next to the threshold the known singular model S is subtracted so ρ−S is analytic and
tanh-sinh converges exponentially; ∫S·K is added back in CLOSED FORM from exact
power-log moments times the kernel's local Taylor coefficients.

USE-WHEN:
- A dispersion integral with a threshold singularity converges only polynomially or
  floors at a few digits under Gauss–Legendre (reference case: GL floored at ~6
  digits with 180 nodes; tanh-sinh with panel breaks at the thresholds alone reached
  15/29/99 digits at level 3/4/5; subtraction takes the rest at the lowest node
  budget).
- ρ(w′) is expensive across many quadrature nodes: pair with the DE-transport ρ
  pattern of tools/wayfinder — one high-precision boundary value, then Taylor steps
  of the sub-integral's differential equation along the node lattice instead of a
  fresh quadrature per node (reference: 1669 nodes at 99.9 digits, cross-validated
  at 3 points).

NOT-FOR: insertions that are not two-particle-reducible; kernels with singularities
inside the dispersion range (K must be analytic on each sub-panel); singular models
beyond power×log² (moment limits below). This member is plain quadrature, not the
certified refuse-don't-guess oracle above: gate its output with a two-level or
two-precision agreement, or feed the panel integrands to `nestor.integrate`.

INVOKE: `from nestor.dispersion import disp_subtracted, moment_powerlog,
kernel_taylor, addback_endpoint, tanhsinh_panel` (tools/ on the path; also
`nestor.disp_subtracted`). Main entry `disp_subtracted(rho, K, panels, dps=40,
kernel_nmax=None, maxdegree=None)`; each panel is `{'a','b','thresh': None |
('left'|'right', w*, [(c,α,m),…]), 'sub_width', 'radius', 'map':'tail'}`.
`moment_powerlog(p, m, L)` = ∫₀^L u^p (log u)^m du in closed form;
`kernel_taylor(K, wstar, nmax, radius=None, dps=40, method='cheb'|'quad')` (the
default Chebyshev collocation uses real evaluations only); `addback_endpoint(coeffs,
Kc, L, side='left'|'right')`.

INPUTS: callables ρ and K; the panel list with thresholds and singular-model
coefficients per threshold side; mpmath dps. OUTPUTS: the dispersion integral at the
target dps (mpf).

GATES: two-precision / two-level agreement per panel; for a DE-transported ρ,
multi-point cross-validation against direct quadrature at a few nodes before trusting
the walk. Battery (selftest legs D1–D4, same assertions as the pytest files): closed
forms ≥30 digits vs adaptive quadrature and add-back ≥25 digits (measured 33); level-5
subtracted assembly on a two-threshold analytic surrogate ≥30 digits vs a tanh-sinh
reference (measured 46); an add-back-sign mutation control that must collapse the
agreement to ≤5 digits, proving the closed-form term is load-bearing.
`python3 tools/nestor/tests/test_dispersion_surrogate.py` also prints the full
Gauss–Legendre / tanh-sinh / subtracted convergence ladder (level 3/4/5 → ~16/29/46
digits).

FOOTGUNS:
- The sub-panel width MUST be smaller than the kernel's analyticity radius at the
  threshold — the add-back is a kernel Taylor series; violating this silently
  corrupts the closed-form term.
- `moment_powerlog` supports log powers m∈{0,1,2} and p>−1 only; outside that,
  restructure the model.
- A finite JUMP needs NO subtraction model: place a panel break exactly at it and
  pass `thresh=None` on both sides — each side is then analytic.

### examples/box1_dilog.py (worked dispersion kernel)

PURPOSE: the ε⁰ part of a one-loop box with three equal internal masses and one large
corner mass (pySecDec normalization, factor 1), reduced from 2D to 1D by two analytic
integrations (four-root partial fractions) plus ONE fixed Gauss–Legendre rule in v.
The singularities of the v-integrand sit at v=−1 and |v|=1 independently of the corner
mass M5, so convergence is uniform in M5; the small quadratic roots are taken by
Vieta (r_k/(M5·large root), no subtraction), so floating-point accuracy is also
uniform in M5 (measured flat 66 digits at dps=50 across M5∈[1e8,1e26]). ≥dps digits
at ~6–8 ms per call, ~40–48× the 2D reference.

USE-WHEN: the box kernel K(w′)=Box(s,t;m²,M5=w′) is needed to ≥45–50 digits per call
at millisecond cost inside dispersion panels (`make_K` memoized factory).

NOT-FOR: anything outside the deep-Euclidean region s<0, t<0, m²>0, M5>0, u=−s−t<4m²
(ValueError at or past the pseudo-threshold). This particular box only, not a general
four-point function; the full dilogarithm closed form is a documented path, not
assembled.

INVOKE: `from box1_dilog import box_closed_dilog, make_K` with `tools/nestor/examples`
on the path — `box_closed_dilog(s, t, m2, M5, dps=40)` (alias `box_closed`);
`make_K(s, t, m2, dps)`. `python3 tools/nestor/examples/box1_dilog.py` reproduces the
validation table against the reference 2D oracle, which is not included in this
repository: set `BOX1_DILOG_ORACLE_DIR` to a directory holding `box1_closed.py`,
otherwise the run prints SKIP and exits 0 (vendor validation: table reproduced at
≥50.95 digits on the 5-point M5 gate, uniformity probe 51.16–51.5 digits).

FOOTGUNS: without the Vieta stabilization ~log10(M5/(2m²)) digits are lost to
small-root cancellation, invisible under the +18 guard digits until M5~1e20 at
dps=50; the dv-first variant is singular at x0=−m²/M5, which collapses onto the
endpoint as M5→∞ (6 digits at M5=500) — dx0-first form only; restores mp.mp.dps on
exit including error paths.

RELATED: tools/wayfinder (the DE-transport ρ walk this member pairs with).
CREDIT: singularity subtraction after Kahaner, Monegato and Davis–Rabinowitz,
specialized to the dispersion-plus-analytic-kernel setting; the box1_dilog example
specializes the scalar one-loop box D0 of Denner, Nierste and Scharf [DNS] in the
conventions of Denner's review [Denner]; the quadrature is Takahasi–Mori tanh-sinh [TM]
via mpmath (keys in REFERENCES.md at the repository root).

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).

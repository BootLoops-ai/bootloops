# posq — certified evidence by exact positive quadrature

posq computes a Bayesian evidence integral as a certified positive value-sum:
a degree-matched exact Gauss sum on the raw likelihood product, plus log-space
patch enclosures for what survives — never expanding coefficients, never
subtracting — so two-sided machine-width certificates come from positivity
alone. The deliverable is the evidence value, and evidence differences / Bayes
factors between models on the same data, with certified two-sided intervals.

Tool page: https://bootloops.ai/tools/posq.html. posq also ships as a member of
the [`baller`](../baller/) package (`quad.posq` imports these bytes). A
paper-register description of the method is in
[`DESCRIPTION.md`](DESCRIPTION.md) (literature positioning, with references,
is on the tool page); the one-page operational summary is
[`GUIDE.md`](GUIDE.md).

## How it works

For the production object (a binary continuous-time Markov chain on a rooted
4-leaf tree with a strict clock), substitute `u_k = exp(-beta*s_k)`. The
Exp(lambda) clock-increment priors then marginalize exactly into

    Z_s(p) = (c^3/2) * INT_{[0,1]^3} u^(c-1) F(u; p) du,   c = 2*lambda*p*(1-p),

where `F` is the raw positive product of pattern probabilities — a polynomial
of known multidegree. A certified Gauss–Jacobi tensor rule built at that degree
(interval-Newton node certification at 2048 bits, with receipt gates: weight
balls strictly positive, weight sum reproduces the exact moment, monomial
exactness verified) makes the clock quadrature exact: zero s-error, no cells,
no remainder engine. Positive weights times a positive integrand give a
sign-positive streamed ball sum, hence certified lower and upper endpoints at
machine width — measured relative width `2^-179` at 192-bit working precision
on the production object. The one remaining dimension (`p`) is integrated by a
non-adaptive Cauchy-ellipse register described under "The p-layer" below.

## What posq is for

Likelihoods of CTMCs on small state graphs with clock priors and sign-positive
products:

- a small state graph (here: binary CTMC on a rooted 4-leaf tree, strict
  clock),
- prior marginalization that lands on polynomial-weight integrals (Exp/Gamma
  clock increments give `u^(c-1)`-class weights after the exp substitution;
  uniform or polynomial-density priors on rate/frequency parameters),
- an integrand that is a raw product of pattern/site probabilities — positive
  on the domain, polynomial of known multidegree in the transformed variables,
  never coefficient-expanded.

Model computed for the results quoted here: binary CTMC, strict clock,
Exp(lambda=10) increments, `pi1 ~ U(0,1)`, uncorrected likelihood.

## Measured results

All numbers below are measured on the production 318-count four-taxon data
set (not shipped in this repository — see "Not included" below). What ships
is the engine itself plus a LIVE verification battery
(`verify_posq_adaptation.py`) that re-runs the kernel and rule builder from
this checkout.

- Exact clock sweep: `Z_s(p0)` over 48,717,760 nodes, `ln Z_s` enclosure width
  about `2e-18` nats (relative width `2^-179.4`), 0.201 CPU-h per sweep
  (14.9 us/node in Python; the C kernel runs the same chain). An
  exact-rational N=8 surrogate is contained at relative difference
  `<= 4.8e-56`, and a planted count corruption is detected.
- First complete certified two-sided evidence row:
  `lnZ in [-833.5351, -833.4725]` (width 0.0627 nats, certified positive lower
  endpoint), 6 panels plus exact endpoint bands, 27.4 CPU-h per topology.
- Full 15-topology table: all 15 topologies certified two-sided, widths
  0.0095–1.157 nats, 93.25 CPU-h total. Cross-checks all pass: a barycentric
  edge-interpolation referee (worst residual `3.0e-4` nats), 3 planted count
  corruptions detected, monotone flank checks. Certified decision sentences:
  maximum-evidence topology `(((A,D),C),B)`, certified `lnBF >= 0.4675` over
  the runner-up and `>= 5.28 / 5.21` over the losing split classes, agreeing
  with an independently computed reference model's ranking.

Decision sentences follow a conservative display rule (implemented in
`derive_posq_sentences.py`): a certified gap is the winner's lower endpoint
minus the rival's upper endpoint; lower endpoints and Bayes-factor floors
round down, upper endpoints round up.

## Usage limits

These are measured limits, not style preferences.

1. **Four-leaf only.** The exact clock collapse is quartet-native: proven and
   priced for rooted 4-leaf shapes only (3 clock dimensions; the 12
   caterpillar topologies share one fiber stage covering 88.47% of the work;
   the 3 balanced shapes go through the min/diff decomposition below). Six
   leaves are not supplied by this tool, and corrected/ascertained likelihoods
   are not covered. A log-space order-k patch route (k=4 measured at `1.7e-3`
   nats/patch, interval inflation factor 11.85) is a measured candidate for
   objects outside this wall, but it is not part of the shipped engine.

2. **The p-layer is the only engineered error dimension**, and its route menu
   is priced by measurement (one sweep-equivalent = 724.6 CPU-s):
   - *Two-exponent sandwich with corner-patch lemma*: sound (corner-patch
     defect certified at `e^-1248`, `2.5e-14` nats on the lower bound) but far
     too expensive — 0.172 nats/panel at radius `1e-4`, 262 core-h per
     topology in Python. Use only if nothing else applies.
   - *Taylor jets through the rule*: center jets are sound and cheap
     (certified `dlnZ_s/dp` at relative width `2^-173.8` for 2.55
     sweep-equivalents) and remain a usable primitive. The jet remainder over
     a fat panel fails at production degrees — interval Newton rides the
     explosive second solution of the three-term recurrence at endpoint
     clusters (even n=17 fails at c-radius `1e-12`) — and the sound
     complete-monotonicity variant costs 1103–1748 sweep-equivalents for
     1.0–0.1 nats. Not viable as the production layer.
   - *Cauchy-ellipse register (shipped)*: certified Gauss–Legendre p-nodes per
     panel, one Bernstein-basis absolute-majorant sweep per panel (the
     pointwise absolute value `|b_k(u)|` is not a polynomial, but the
     same-degree Bernstein-basis absolute majorant is exactly integrable —
     that substitution is part of the register), and exact endpoint bands
     (`P_x <= 8p` / `8(1-p)`; no band sweeps needed). Measured cost: 186–190
     sweep-equivalents, i.e. 37–38 CPU-h per topology, nearly independent of
     whether the width target is 0.1, 0.5, or 1.0 nats (the majorant drives
     the cost) — if you can afford it at all, it is a 0.1-nat instrument.
     Production came in at 27.4 CPU-h per topology at width 0.0627 nats. The
     absolute majorant pays a measured penalty over signed sensitivity
     (effective rate ~990 per unit-p vs ~230 signed): panel counts are set by
     the majorant — do not re-derive them from signed sensitivity.

3. **Balanced-topology exactness gates.** Balanced shapes are computed via a
   min/diff decomposition (`m ~ Exp(2*lambda)` giving a `w1^(2c-1)` rule;
   `d, s3 ~ Exp(lambda)`; two ordered cases summed). Its exactness gates must
   pass before any sweep: per-character degrees (4,3,1) with `KMAX_BAL = 4`;
   exact w3-parity; exact sum-to-one; the dual (index-arithmetic vs
   string-permutation) exponent-vector constructions asserted equal; measure
   factorization; the p-Taylor-shift identity checked exactly on `Fraction`;
   per-panel pole caps — the nearest pole is `c = -1/2` (at
   `p = (1 - sqrt(11/10))/2`), so each panel's ellipse parameter must sit
   strictly below the pole's Bernstein image on both flanks; and per-panel
   dominance trials. The production driver re-asserted all of these on every
   run and refused to sweep on failure; any new driver must do the same.

4. **The assembly must stay cancellation-free.** Everything is a
   positive-weighted sum of positive balls by construction, and everything
   that could break that is screened, not assumed: the only mixed-sign step is
   the per-node value `P_y = A_y + B_y*v`; every `P` ball and every assembled
   `Z` ball must certify strictly positive; per-sweep relative width worse
   than `2^-100` is a hard failure (the driver raises, it does not warn); and
   boundary-clustered nodes are screened explicitly — corner nodes sit within
   about `n^-2` of the boundary, where the mixed-sign evaluation is tightest,
   and `F` is certified positive at relative width `2^-166` to `2^-179` at
   192-bit at all 8 rule corners. Never expand coefficients: the two failure
   mechanisms this tool exists to avoid are expansion cost and expansion
   cancellation (a 560-bit-precision class); value evaluation at 192-bit is
   the rule.

## Files

Engine (do not edit in place — any change re-runs the live battery first):

- `posq_kernel.c` / `posq_kernel` — C (FLINT/Arb) multi-topology certified
  sweep kernel: caterpillar/balanced/multi-sweep modes, shared squaring
  tables, exponent vectors over 16 shared pattern values, exact ball
  serialization in and out. No binary is shipped: `kernel_io` compiles the
  kernel on first use, or build it yourself with
  `gcc -O2 -o posq_kernel posq_kernel.c -lflint -lmpfr -lgmp -lm`.
- `kernel_io.py` — exact ball serialization between Python and the kernel
  (radius rounding is upward, hence sound).
- `bench_rung1_gauss.py` — the certified Gauss–Jacobi rule builder
  (interval-Newton at 2048 bits with receipt gates), the Python reference
  sweep, and the exact-rational N=8 surrogate.
- `topos.py`, `exact_polys.py`, `exact_polys_bal.py`, `bal_fiber.py` — exact
  structure: 15-topology enumeration, pattern polynomials, the balanced
  min/diff decomposition.
- `derive_posq_sentences.py` — the certified-sentence derivation and display
  rule.

Packaging and verification surface:

- `posq.py` — the package surface: `SurrogateExact` (the N=8 exact-rational
  cross-check object as an exact closed form with exact-derivative jets) and
  `PosqKernelAdaptation` (the verification-battery adaptation; see below).
  `python3 posq.py --selftest` runs the exact-surrogate pin plus one live
  kernel-containment leg in seconds.
- `verify_posq_adaptation.py` — the battery driver (see Verification).

Closure interface (a separate, read-only evaluation surface for plug-in
log-likelihood closures): `closure_interface.py`,
`selftest_closure_interface.py`, `closure_fixtures/`, and
`CLOSURE_INTERFACE_MANIFEST.sha256` — documented in
[`CLOSURE_INTERFACE_README.md`](CLOSURE_INTERFACE_README.md). Note its numbers
are float-class estimates, not certified enclosures.

Not included in this repository: the production driver (`run_quartet.py`),
the production result JSONs (they live with the problem pages on
bootloops.ai), and the production pattern-count data
(`QUARTET_COLLAPSE.json`). Some scripts here (`derive_posq_sentences.py`,
`topos.py`'s exponent vectors, parts of `bench_rung1_gauss.py`) read
production input data through environment variables and refuse with a named
error when those point nowhere; the engine, the surrogate objects, and the
verification battery run from this checkout alone.

## Requirements

- Python 3 with `python-flint` (Arb ball arithmetic) for all engine paths.
- A C compiler with the FLINT, MPFR, and GMP development headers to build the
  C kernel (built automatically on first use; no binary is included).
- `numpy` and `scipy` for the closure-interface selftest only.

## Verification

The packaged engine was run as a consumer adaptation through an adversarial
verification battery (the `eras` verifier, `verify_eras_adversarial.py`):
parameter vector `(p, lambda)`; the C kernel as the point engine, with the
`(17,13,5)` rules rebuilt and receipt-gated at each point's exact rational `c`
(the same rebuild-per-p discipline as production); a monotone c-bracket times
ball-q closed-box enclosure form; and exact closed-form derivative jets as an
independent calculus path with no quadrature. The battery object is the N=8
exact-rational surrogate — the same rule builder at the same spec and the same
kernel in the same mode — because production-scale points cost 0.2 CPU-h each
and the battery needs about 700 of them; production-scale validation rests on
the containment and planted-corruption controls of the measured results
above.

Result: pass at 3 seeds (`987654321202607`, `20260718`, `424242424242`), each
with: domain gates passed (the exact closed form reproduces the pinned exact
value by exact `Fraction` identity; the kernel ball contains the exact
rational at relative width `2^-248`; a planted count corruption is detected at
`dln = 0.195`; 9/9 exact-value containments); 224/224 closed-box containments
per shell across 3 shells (corners at the radius, faces, a near-face band, two
fresh batches); 3/3 informative coverage cells (width/spread 6–33x, well under
the 1000x vacuousness threshold); all 6 planted corruptions detected
(narrowing by `1-1e-6` and by 0.5; absolute shifts at `width*1e-3` and
scale-adaptive; dropped-dimension for both dimensions); a precision probe
caught (a 53-bit rule rebuild fails interval-Newton certification and returns
no-information rather than a wrong number); finite-difference vs
exact-calculus derivatives agreeing to a worst error of `3.1e-23` against a
`1.1e-12` gate; zero fatals. The battery is live in this directory — re-run
`python3 verify_posq_adaptation.py` and the run you make is the record.

Honest limits of that claim:

- A single-seed pass is necessary, not sufficient. Three seeds are recorded;
  vary `--seed` further before relying on new adaptations.
- The battery exercises the engine at surrogate scale (`n = (17,13,5)`,
  `N = 8`); the production-scale claims rest on the production controls, not
  on this battery.
- The c-bracket enclosure form in `posq.py` is battery scaffolding (valid only
  on p-boxes inside `(0, 1/2)`; it returns no-information outside), not a
  production p-layer — the production p-layer is the ellipse register above.
- Containment batteries have residual limits: thin-set corruptions are
  detected with probability scaling in points times seeds; a corruption
  smaller than the enclosure slack is undetectable by containment in
  principle; the battery defends against honest bugs, not against edits to
  the verifier itself.

`verify_posq_adaptation.py` (usage:
`python3 verify_posq_adaptation.py [--seed N] [--points N]`; exit 0 pass /
1 fail / 2 indeterminate) imports the `eras` verifier from the `tools/eras/`
directory beside this one and runs end-to-end from this checkout.

## Adapting posq to a new model

Small-graph CTMC evidence under clock-style priors (for example,
coalescent-class priors) is the same object class, but nothing transfers by
analogy — only by receipts:

- Re-derive the degree bookkeeping and the prior-to-weight mapping for the new
  model; rebuild the Gauss rule at the object's own multidegree and re-assert
  the dominance and parity gates.
- Run a measured trial first (real counts through a production-style sweep;
  record widths and CPU costs).
- Rerun the verification battery against the adaptation, following the
  `verify_posq_adaptation.py` pattern; check `result.ok`; use several seeds.
- The four-leaf limit is not lifted by adaptation.

CREDIT: the kernel is written directly against FLINT/Arb (W. Hart, F. Johansson, A.
Ahlbäck and the FLINT developers; Arb: Johansson 2017, IEEE Trans. Comput. 66:1281);
certified Gauss–Legendre/Gauss–Jacobi rules are obtained by interval Newton in ball
arithmetic, in the manner of Johansson & Mezzarobba (2018, SIAM J. Sci. Comput. 40:C726)
for the Gauss–Legendre case, seeded from SciPy's double-precision Gauss–Jacobi nodes
(scipy.special.roots_jacobi); the likelihood integrated is Felsenstein's (1981, J. Mol.
Evol. 17:368) CTMC tree likelihood; the production counts derive from DravLex v1.0
(Kolipakam, Jordan, Dunn, Greenhill, Bouckaert, Gray & Verkerk 2018, R. Soc. Open Sci.
5:171504; CC BY 4.0; not shipped here).

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).

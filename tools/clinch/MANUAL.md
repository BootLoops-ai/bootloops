# clinch — manual

**CLINCH** — Certified Local Interval-Newton Convergence for Hierarchies.
Given a fitted MAP optimum theta* and model oracles evaluable in ball
arithmetic, CLINCH either issues a CERTIFICATE — inside the box around the
candidate a true stationary point of the penalized objective EXISTS, is
UNIQUE, and the Hessian is positive definite throughout (a certified strict
local minimum) — or REFUSES with the failed contraction NAMED. Verdict set: CERTIFIED-INTERIOR / CERTIFIED-CORNER /
REFUSED(named) / OUT-OF-SCOPE. No third outcome; a fat ball is a refusal.

SCOPE: plain local interval certification
exists in sibling packages. CLINCH's axis is the
HIERARCHICAL/BLOCK one — block-arrow Krawczyk + interval Schur at
latent-group scale (5592-9 coordinate certificates measured in ~90 s) —
and every outward sentence scopes it so.

LAYOUT: package clinch/ + battery/ + adapters/ (etienne, certpass) +
census/ + this manual.
KIND: package.

MODEL SPACES AND REFERENCE DATA. The shipped model oracles are named after
the layout labels of the hierarchical demographic model they were built for:
`v3.1` is the uncentered dial layout (`oracle_v31`, `adapter_v31`), `v3.6` the
centered layout and `v3.7` the centered layout with the advance-level pin
(both in `oracle_v36`). These labels name model spaces, not versions of this
package. The case studies below were run on four reference fits of that
model, labeled `fit13`, `fit17` (v3.1 layout) and `fit16`, `fit18a`
(v3.6/v3.7 layout). The reference model code, the fit artifacts and their
certificate receipts are not shipped; every loader and battery leg that
needs them is gated on the environment variables `CLINCH_CODE_DIR` (the
reference model code) and `CLINCH_REFERENCE_DIR` (the reference fits and
their stored results) and refuses or skips by name when they are unset. Everything else (engine, jets, kkt_pd, mask, cert,
the synthetic fixtures and battery legs L1, L2, L5, L7, L8) is
self-contained.

## Quickstart (self-contained)
```python
import sys
sys.path.insert(0, "tools/baller"); sys.path.insert(0, "tools/clinch")
sys.path.insert(0, "tools/clinch/battery")
import clinch                        # verify() runs fail-closed at first use
from clinch import engine
from synthetic import SynthOracle, exact_optimum   # planted hierarchical problem
zs, g = exact_optimum()              # exact rational optimum (Fraction solve)
theta = [float(v) for v in zs[0] + zs[1] + g]
res = engine.certify(SynthOracle(), theta, prec=128)
print(res["verdict"], res["radius_min"])            # CERTIFIED-INTERIOR
```

## Quickstart (a reference-model fit; needs the reference data, see above)
```python
from clinch import adapter_v31 as A, engine, cert
from clinch.oracle_v31 import ModelV31Oracle
arrs, consts, theta, prov = A.load("fit13")     # reference loader + npz cache
o = ModelV31Oracle(arrs, consts)                # check_space-guarded
gate = A.gate_gradient_adjudicated(o, theta, "fit13")   # BEFORE any cert
res  = engine.certify(o, theta)                 # direct certificate
thp, prcpt = engine.newton_polish(o, theta)     # labelled second product
resp = engine.certify(o, thp)                   # polished-center cert
cert.write_cert(res, out_dir, "fit13_null", gates=...)  # CERT.json
```
For your own model: write an oracle with the `SynthOracle` interface
(block-arrow partition, `F`/`H` over arb balls via `clinch.jets`) and call
`engine.certify(oracle, theta)`.

## Architecture (two pieces at the package boundary)
- The generic certification kernel is **baller.certify.block_krawczyk**
  (blockwise approximate inverse + interval Schur complement on the
  border; strict-interior contraction; per-coordinate radius contract;
  global block-arrow congruence PD leg), consumed BY IDENTITY — never
  copied. Its battery is BALLER L16.
- CLINCH owns: the model ORACLE (gradient + block-arrow Hessian in arb via
  second-order jets), the identity gates, the radius policy, the corner
  (KKT) and out-of-scope semantics, and the products (CERT.json +
  human statement).

## Modules
- **jets** — N-var second-order jets over arb balls (exp/log/log1p/expm1/
  expit/lgamma/pow_const; order-1 mode for gradient-only passes). Carries
  the certified MONOTONE ENDPOINT HULLS dig_ball/trig_ball/lgam_ball —
  arb's own ball digamma NaNs at rad ~ 0.85 and is garbage-wide before
  that (measured with python-flint 0.8); psi/psi' are monotone on (0, inf)
  so endpoint hulls are rigorous AND tight.
- **oracle_v31** — the oracle of the v3.1 (uncentered) model space
  (51 globals + 3x(S+NG+NF) latents). check_space REFUSES any other slot set cold (OUT-OF-SCOPE at
  the engine level). Partition: (process, family-component) blocks with
  VERIFIED genus->family consistency (inconsistent genera merge; the
  block-diagonality of the latent Hessian is then exact — v3.1 has no
  centering). LOG-SPACE likelihood forms throughout (log-sigmoid, expm1,
  shifted-exponent fecundity A, branch-stable log(phi+mu)) — naive
  expit/product balls lose positivity over wide boxes (measured at
  r_unit ~ 0.06 sigma). CLIP GUARDS prove every reference clip strictly
  inactive over the box or raise a NAMED OracleRefusal.
- **adapter_v31** — reference data via the runners' own load_data (refbridge by
  identity; REGISTRY keys asserted) + sha-pinned npz cache; the ORACLE
  IDENTITY GATES: gate_gradient (registered float64-register metric,
  receipted VERBATIM) and gate_gradient_adjudicated (the top-K
  disagreement coordinates re-derived in mpmath at 50 dps via the
  reference's own hand-derived chain — derivation-independent;
  bar 1e-12 at matched precision), plus gate_fd_hessian (assembled arb
  Hessian vs central FD of the arb gradient, exact power-of-two step).
- **oracle_v36** — the CENTERED-space oracle (v3.6 unpinned / v3.7
  pinned) via the AUXILIARY-BORDER KKT EXTENSION: hard mean-zero
  centering makes the latent Hessian family-block + rank-2, so ubar
  becomes an auxiliary BORDER variable per centered channel with its
  constraint multiplier (border = 53 globals + 4 aux; the v3.7
  advance-level pin enters as ONE more KKT border row, lam, touching
  only border coords + m_a — the pin statistic T in extended variables).
  Zeros of the extended system correspond exactly to (constrained)
  stationary points of the true composed objective; the theta-theta
  block is family-block-arrow again and baller's Krawczyk applies
  unchanged. The extended Jacobian is a KKT matrix (indefinite by
  construction), so the KRAWCZYK leg certifies the KKT point; the PD
  statement for centered spaces is the REDUCED Hessian's — the rank-4
  congruence extension, carried by **kkt_pd** (below); the oracle
  declares its multiplier rows (kkt_multiplier_idx) and the engine
  routes the PD leg there. Survival carries the
  v3.6 colimited surface (per-site f_sat enclosed from the log_m0b
  ball with exact dfs/d2fs chains). Identity-gated at the two centered
  reference fits (fit16, fit18a: 3.0e-10 / 4.2e-10, reference-noise class;
  aux rows resolve to 1e-16; FD gate incl. aux columns <= 9.4e-11).
- **kkt_pd** — the reduced-Hessian PD certificate for CENTERED
  (extended-KKT) spaces: the rank-4 congruence extension. The extended
  Jacobian K = [[H_L, A^T], [A, 0]] is indefinite by construction, so
  the plain PD target does not apply; pd_reduced_congruence builds the
  SAME global block-arrow float factor as the PD leg but with a SIGNED
  border eigenfactor, encloses M = X^T K X in arb against the signed
  target diag(+1…, ±1 border pattern), and every row sum of |M - T|
  provably < 1 certifies the INERTIA of every matrix in the enclosure
  (Weyl + Sylvester; M nonsingular forces X nonsingular, so no rigorous
  inverse anywhere). Exactly one certified negative direction per
  multiplier row then gives, by the classical KKT inertia identity
  (Gould 1985; nonsingularity of K certifies the constraint rows' full
  rank for free), positive definiteness of the reduced Hessian
  Z^T H_L Z on the constraint tangent space — with the Krawczyk KKT
  point, second-order sufficiency: a certified strict local minimum of
  the equality-constrained problem. Refusals: kkt-signature (wrong
  negative count — the indefinite-reduced-Hessian verdict, counts
  named), midpoint_cholesky:block_<i>, block/border congruence rows.
  Consumed by engine when the oracle declares kkt_multiplier_idx;
  battery leg L8 proves it on the toy centered fixture at rank 2.
- **mask** — MaskedOracle: fixes the ACTIVE SET at declared bounds and
  presents the reduced block-arrow system (the corner path).
- **engine** — the certificate driver; all policy constants CODE-PINNED
  (no knobs): PURE curvature-scaled box r_k = r_unit/sqrt(mid H_kk)
  (uniform in the Newton metric — the only geometry that survives, see
  FOOTGUNS), r_unit ladder self-seeded at STEP_COVER x the scaled Newton
  residual (climb-only), FAT-BALL LAW (radii capped at 1 curvature-sigma:
  a certificate around a poorly-converged candidate would be true but
  WEAK — refused instead, with the residual named), wall budget 1800 s
  (checkpoint-and-report, never silent grinding). newton_polish +
  polished-center certificate = the SEPARATE, LABELLED product
  (certify at the polished center, receipt the distance to the
  published point). certify_guarded -> OUT-OF-SCOPE.
- **cert** — CERT.json (ball-center sha, per-coordinate radius range,
  contraction margins per NAMED block (process:family), PD margins,
  oracle identity-gate receipts, arb precision, rung history, wall
  times) + the one-paragraph human statement. LOCAL statements only —
  the engine makes no global-optimality claim anywhere, including prose.

## GATES (battery/battery.py — run from a scratch cwd; refuses the tool
## tree and the protected roots)
L1 planted optimum (exact Fraction truth) certified + contained; corner
   candidate CERTIFIED-CORNER with strict KKT margins; undeclared-bound
   control -> fat-ball refusal; out-of-scope typed · L2 planted saddle:
   Krawczyk certifies the stationary point, PD REFUSES naming block +
   direction · L3 flat-valley control GREEN (results below) · L4 tamper:
   math sign-flip refused at fixed radius; bit-flipped scratch copy -> verify()
   refuses typed; poisoned pyc inert; sign flip in the real reference-fit
   oracle inputs moves nll 6.6e4 vs the receipted gate value · L5 dps 30/50/80
   agree; prec=8 refuses · L6 flagship receipts verified (below) ·
   L7 engine exhaustion path: an
   oracle that raises on EVERY radius rung gets the TYPED refusal
   REFUSED(all-rungs-failed) with the per-rung failure classes listed —
   never a best=None crash · L8 centered-space reduced-Hessian PD: the
   toy CENTERED fixture (rank-2 auxiliary-border KKT extension) gets
   the FULL certificate — Krawczyk KKT point + PD_REDUCED_CERTIFIED
   with the certified inertia, exact rational KKT point contained;
   flipped-border control -> typed kkt-signature refusal (counts
   named); flipped-latent control -> midpoint-Cholesky refusal naming
   the block.
OVERALL: PASS, all 8 legs (with the reference data present; without it
L3, L6 and L4's real-oracle sub-leg skip by design, see `selftest.py`).
Receipt producers: battery/flagship.py (L6 certificates), run_gates.py /
run_fdgate.py (oracle identity gates), l3_gates.py / l3_run.py (L3).

## CASE STUDY 1 — two uncentered reference fits (fit13, fit17)
A worked application of the tool to a real ~5,600-coordinate production
fit. The receipts (clinch_fit13_null/ and clinch_fit17_null/, CERT.json
each, sha-pinned) and the underlying fit data are reference data (not
shipped, see above); battery leg L6 re-verifies the certificates where they
are present.
- Identity gates: fit13 registered-metric 4.97e-10 / fit17 2.53e-10
  (float64 register, receipted verbatim; the reference's OWN rounding —
  PROVEN by mpmath adjudication: at every top-5 disagreement coordinate
  |arb - mp50| <= 1.4e-13 while |ref_float64 - mp50| equals the
  registered metric). FD-Hessian gates ok (<= 1.7e-10 / 7.1e-9 rel).
- fit13_null DIRECT: REFUSED — the receipted theta's scaled Newton
  residual is 0.0286 curvature-sigma (raw gmax 9.6e-5 LOOKS converged;
  in the Newton metric it is 4 orders too loose): the covering radius
  forces Hessian-variation margins ~ 4.7e10 and no radius wins (margin
  ~ 2.9e4 even at 1e-6 sigma, c-part-dominated; superlinear width above
  1e-3 sigma). THE FINDING THE TOOL EXISTS FOR: "the optimizer
  converged" at the fit's own bar is not certifiable as stated.
- fit13_null POLISHED-CENTER: 2 structured-Newton steps take the
  residual 0.0286 -> 5.0e-11 sigma; theta+ at distance 4.64e-5 (inf-norm)
  from the receipted theta. CERTIFIED-INTERIOR: 5592 coordinates, 582
  blocks + 51-border, contraction max margin 0.735 < 1, PD margins >=
  0.98, radii [1.8e-13, 1e-9] (curvature-scaled), prec 192 bits, ~90 s.
- fit17_null: same shape — DIRECT REFUSED (border, best margin 3.3e5);
  POLISHED-CENTER (1 Newton step, distance 1.83e-6) CERTIFIED-INTERIOR:
  5955 coordinates, margin 0.834, PD >= 0.978.
- Bounds: fit13/fit17 are UNCONSTRAINED fits (no bounds kwarg anywhere in
  the runners/polish) — no dial sits at a bound; the corner register is
  exercised by battery L1 and stands ready.

## CASE STUDY 2 — the flat-valley control (centered fits fit16, fit18a)
The flat-valley control, run on the two centered reference fits (not
shipped):
- The v3.6 stall optimum (the fit's stage-1 checkpoint, gmax 0.026999):
  REFUSED. The Newton correction there is
  0.97 CURVATURE-SIGMA and concentrates on the ADVANCE-LEVEL machinery —
  a_adv intercepts + the sp_adv pooling scale log_sig[3] (the exact
  level-redistribution pair) + advance slopes: the
  refusal points at the advance-level direction. At its coverage radius
  the enclosure NaNs (1000+ rail crossings, census receipted) — the
  fail-closed register.
- The SAME model under the v3.7 pin (fit18a, extended-KKT system,
  deep-polished center at distance 2.6e-5 from the receipted init;
  residual 7.7e-14 sigma): Krawczyk CERTIFIED, margin 0.626 at
  r_unit 3e-10 sigma. BORDER CONTRACTION DEFECT AT THE COMMON RADIUS:
  fit16-stall 3.23e9 vs fit18a 0.626 — ratio 5.16e9 (~9.7 orders of
  magnitude). The pin removes the valley.
- Flat-valley geometry, measured: the two receipted optima sit 0.77
  apart along the advance-level direction (a_adv shift +0.769 uniform;
  ubar_adv -0.001 -> +0.778) at ONLY 0.651 nll cost — the near-flat
  direction, quantified; the final fit16 theta's gauge gradient is
  8.0e-11 (prior-settled).
FINDINGS (receipted, nothing weakened): (1) the fit16 stall optimum is
the stage-1 checkpoint; the fit's final theta is a later tail-polished
point (gmax 7.6e-6, gauge-settled; its own direct refusal is the ordinary
under-polish class, margin 1.2e5 at coverage 5.9e-5). (2) At fit18a's
optimum the reference's q/stay5 clip rail is genuinely ACTIVE in 4
advance cells with (n-y)=0 / ~1e-10 value effect — which forced the
CLIP-CENSUS contract (below). (3) The extended KKT Jacobian is
indefinite by construction, so L3's certificate register is the Krawczyk
leg (existence+uniqueness of the KKT point); the PD statement for
centered spaces (reduced-Hessian rank-4 congruence) is kkt_pd's own leg,
proven self-contained by battery leg L8.

## CLIP-CENSUS CONTRACT
The certified target is the UNCLIPPED penalized objective — the
objective whose gradient the reference's analytic chain computes
(clips-as-identity — the reference port's documented semantics). Reference clip
bounds the box cannot be proven clear of are COUNTED per cell class and
receipted (oracle.last_clip_census -> CERT), never silently crossed and
never a hard refusal: at every v3.1 certificate the census is ZERO
(equivalent to the old guards passing); fit18a's census (av_q_stay5: 4)
is the honest record of the active rail. Enclosure-domain failures still
refuse fail-closed (NaN -> nonfinite_oracle).

## FOOTGUNS (measured on the fit13 reference fit)
- RADIUS GEOMETRY: a raw-uniform box suffers Bauer-Skeel width
  amplification ~ cond(H) (~6 decades of curvature; margin 7.7e6); mixing
  per-coordinate Newton-step floors into the radii creates radius RATIOS
  spanning decades that amplify the cross-block terms by 1/r_k (margin
  8.9e3 with the own-row width term at 1.0 exactly). ONLY the pure
  curvature scaling survives. Both wrong geometries are documented in
  engine.py's header.
- WIDE-BALL ENCLOSURES: never expit/products of wide positive balls
  (mid+/-rad products include 0; log then NaNs) — the oracle is log-space
  end to end; never x**n on a possibly-zero ball (baller law); arb ball
  digamma NaNs at rad ~ 0.85 and is garbage-wide before that — use the
  jets' monotone endpoint hulls.
- PD LEG (in BALLER, consumed here): per-block-Schur routes FAIL on
  collinear monospecific-genus latents (cond ~ 1e13): the Neumann bound
  rank-collapses (eps ~ 1e-3 x ||W||_1^2 ~ 1e2) and arb's certified
  interval solve fattens by ~1e2 — the GLOBAL block-arrow congruence
  (one float block-triangular factor, exactly nonsingular; no rigorous
  inverse anywhere) certifies the same enclosure with eps ~ 1e-2.
- The identity gate's registered 1e-12 bar is UNATTAINABLE against a
  float64 reference in the max(1,|g|) register (the reference's own
  rounding is ~5e-10 on FIA-scale coordinates) — the identity claim rides
  the matched-precision mpmath adjudication register; BOTH numbers are
  receipted verbatim, neither is weakened.
- Certificates are AT THE CANDIDATE AS GIVEN. The polished-center
  certificate is a second, labelled product; its statement includes the
  measured distance. Never blur the two.

LIMITATION: certified model comparison (dAIC/LRT; max-vs-max with both
optima enclosed) is out of scope; the verdict register is per optimum.

RELATED: baller (the certification kernel; tools/baller/MANUAL.md, battery leg L16);
eras (mesh sibling).
CREDIT: the certification kernel is the Krawczyk operator (R. Krawczyk 1969, Computing
4:187) with Moore's interior test (R. E. Moore 1977, SIAM J. Numer. Anal. 14:611) in the
verification semantics of S. M. Rump (2010, Acta Numerica 19:287), assembled blockwise
in baller.certify.block_krawczyk over Arb (F. Johansson 2017) via python-flint.
adapters/certpass certifies optima of the ETAS model (Y. Ogata 1988, JASA 83:9) as
fitted by the open `etas` package of L. Mizrahi, S. Nandan & S. Wiemer (2021, Seismol.
Res. Lett. 92:2333; github.com/lmizrahi/etas, MIT), whose EM follows Veen & Schoenberg
(2008, JASA 103:614) — we thank its authors for a codebase clean enough to certify
against; adapters/etienne certifies optima of Etienne's (2005, 2007) neutral sampling
formula. The hierarchical reference model and its fit code are an in-house study,
consumed read-only.

## MEMBERS (absorbed external certified-optimum instruments)
Two externally built certified-optimum instruments are absorbed into
this tree, code-only; each member dir carries PROVENANCE.md +
PROVENANCE.sha256.
- adapters/certpass/ — certified point-process optima (incomplete-data
  ETAS reference fits, coded-Poisson multipliers, conditional-spatial
  multinomial), the second external adapter after adapters/etienne/
  (promote-on-reuse). The originating study's fit artifacts and helper
  modules are not shipped: the run_*/assemble_* scripts resolve them
  through CERTPASS_REFERENCE_DIR and refuse loudly when it is unset; the
  oracles and gates are self-contained.
- census/repertoire_g1/ — Gate-1 complete-census certificates:
  V/INV/EXCL/BB disposition classes, the disposition array itself is the
  completeness certificate, UNKNOWN recorded never forced, 128-bit pinned.
  rederive_lib (the non-arb fixed-point verifier) is declared under
  baller, not here. Receipts land in $G1_OUT (default ./g1_out under the
  caller's cwd).

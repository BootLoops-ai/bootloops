# tropical-sampler — GUIDE

Tool page: https://bootloops.ai/tools/index.html (no dedicated page; listed in the index)

KIND: script (pure python/numpy/scipy — no Polymake, no Wolfram). Single file
`tropical_sampler.py`, imported as a library.

MEMBER: `cegm_gj/` — the CEGM tropical-cone Gauss–Jacobi evaluator, the
certified-value leg on the same fan machinery (section below).

REQUIREMENTS: `numpy` + `scipy` (`scipy.spatial.ConvexHull` builds the hulls);
scipy is not in the common core — `pip install scipy`.

## PURPOSE

Borinsky–Sattelberger–Sturmfels–Telen (BSST) tropical importance sampling (Algorithm 1 of
arXiv:2204.06414, reimplemented in
pure python/numpy/scipy) for Euler-type integrals
∫_{ℝ₊^E} ∏_k Q_k(t)^{u_k} ∏_e (1+t_e)^{−c} dt with subtraction-free Q_k. Exact-integer
Minkowski normal fan (Bareiss determinants, Fraction sector integrals,
fan-refinement-invariant I^tr) + TWO additions over plain BSST: per-cone Neyman
stratification and pilot-fitted per-cone exponential tilting (proposal exact, weights
bounded, unbiased).

## USE-WHEN

- Numeric oracle needed for a toric/GKZ/Bayesian-evidence Euler integral with
  subtraction-free presentation (e.g. phylogenetic evidence integrals; generic in E
  — candidate oracle for toric Bayesian/GKZ integrals broadly).
- Data-sized exponents make plain tropical MC useless → the tilted variant (see
  FOOTGUNS).

## NOT-FOR

Integrands that are not subtraction-free (positivity-gated Q_k presentation is the input
contract). A few digits in minutes, not a high-precision engine (~3–4 digits in ~5
min/topology at 7D / 1670–2540 cones measured) — use as ranking/oracle leg, PSLQ-grade
digits come from elsewhere. For certified 30–60-digit values on the cone
decomposition, use the `cegm_gj` member (below).

## INVOKE

`import tropical_sampler` — build the fan from supp Q_k, then sample: cone ~ c_σ/I^tr,
y = R·λ with λ_i ~ Exp(β_i), Ẑ = prefac · I^tr · mean(w). Key internals: `idet`
(Bareiss), `primitive_normal`, `hull_vertices`; per-cone c_σ = |det R| / ∏β_i exact
Fraction with β_i = −⟨g, r_i⟩ > 0.

INPUTS: supports/exponents of the Q_k (integer point sets), exponents u_k, edge-decay
exponent c; works in log coords y = log t with tropical PL exponent Φ(y).
OUTPUTS: unbiased estimate of the integral + MC error; exact-rational I^tr and per-cone
weights (fan-refinement invariant — an internal consistency check).

## GATES

validated vs exact rational values (mutation control 273σ separation on the
planted-wrong vs −0.51σ on truth) and vs an independent exact leg on all 15 five-taxon
topologies, max pull 2.23σ, ranking exact on every distinct level — reproduce that
pattern (positive + mutation control) when adapting to a new integrand family.

## FOOTGUNS

- Plain Alg-1 is coefficient-blind and STALLS (~2.5 d projected) when data-sized
  exponents make the weights span e^50+ — the pilot-fitted exponential tilting is what
  restores 1/√n convergence. Do not run the untilted sampler on data-sized u_k.
- qhull flat-simplex crash in non-simple high-D cones — handled in-file (rank guard + QJ
  fallback); keep the guard if you fork.
- Sector integrals must stay exact Fractions until the final float — the unbiasedness
  argument leans on the exact fan.

## What runs publicly (battery: `python3 tropical_sampler.py --selftest`)

- Import + fan-construction smoke, built into the module as `--selftest`
  (deterministic seeded rng, seconds-scale, exit code = number of failed
  legs): E=2 three-group fan builds (4 cones, exact I^tr = 25/16), `tiling_probe` exact
  (500 directions each covered exactly once), and a weight-0 `refine_extra` summand
  (7 cones) reproduces the SAME exact Fraction I^tr — the documented
  fan-refinement-invariance consistency check.
- The full validation pattern (15 five-taxon topologies vs an exact leg +
  mutation control) was run on a 15-topology phylogenetic data set (not included);
  reproduce the positive + mutation-control pattern on your own integrand family per
  GATES.

## Member: `cegm_gj/` — CEGM tropical-cone Gauss–Jacobi evaluator (certified-value leg)

KIND: per-family evaluator package (~55K py + 9K C): `tropical_cones.py` (fan
build + 20k-direction coverage), `cone_engine.py` (production cone evaluator),
`x36_engine.py` (GJ machinery + the plain tensor engine, kept because it
documents the measured failure), `setup_x36.py` / `run_pilot.py`
(parametrization verification, kinematics, oracle, QMC cross-check),
`cone_identity_check.py`, `control_oracle_check.py`, `fan_x36.json` (sha256
04807e002de2..., the certified fan fixture), `t3_prod.py` (production driver:
sanity/cert/run), `cchunk.c` (C-MPFR inner loop), `gates.py` (two-run
floored-digit gate), `bank.py` (value-bank assembler; data-gated — it reads
the run JSONs a production run emits, set `CEGM_GJ_RUNDIR`).

REQUIREMENTS: `gmpy2`, `mpmath`, `sympy` (numpy/scipy shared with the
sampler). The optional C engine is not shipped compiled — build it in place
first (`gcc -O2 -shared -fPIC -o cchunk.so cchunk.c -lmpfr -lgmp`, needs
libmpfr/libgmp dev headers); `--engine c` fails to load without it, the
default engine is `plain`.

PURPOSE: certified high-precision values of AHL Grassmannian string integrals
(X(3,6) = G+(3,6)/T production instance): log-space decomposed into the
maximal cones of the common-refinement Newton fan (48 maximal cones,
reproducing AHL's P(3,6) f-vector as an independent structural check; 52
unimodular simplicial subcones), per-cone tensor Gauss–Jacobi (geometric
convergence), support-hoisting, C-MPFR inner loop.

VS THE SAMPLER (why both legs): same fan machinery, different quadrature.
`tropical_sampler.py` = unbiased MC oracle, ~3–4 digits in minutes, generic
in E. `cegm_gj` = certified-value leg, 30–60+ digits with two-run/
two-precision gates, but per-family setup (fan + kinematics + certificates
are X(3,6)-specific; the construction applies verbatim to X(3,7) at cost
n^6 — measure a small probe before committing).

INVOKE (from any run dir — artifacts land in `$CEGM_GJ_OUT` or cwd, never
the tools tree):

- `python3 <path>/t3_prod.py sanity [--procs P] [--engine plain|hoist|c]` —
  gated selftest, rc=0 iff rel-vs-pilot-n10 < 1e-38.
- `python3 <path>/t3_prod.py cert --point {1,2}` — kinematic-point
  certificate (exact conservation, AHL Claim-1 LP + 8 perturbations,
  kappa>0 all subcones).
- `python3 <path>/t3_prod.py run --point {1,2} --n N --prec BITS [--procs P]
  [--engine c]` — production value + per-run JSON with provenance.
- `python3 <path>/gates.py A B [--bar D] [--relmax R]` — floored-integer
  digit agreement between two runs; `bank.py` (set `CEGM_GJ_RUNDIR`)
  reassembles the value bank, fails loudly on gate regression.

MEASURED: plain 4-dim tensor GJ REFUTED at measured n^-1.37 (corner
singularities); cone-route control digits = 1.158n − 0.08; generic corrected
fit digits = 0.7825n + 4.64 (residual ±0.12d; the pilot fit 0.875n + 2.28 is
low-n-biased — overpredicts ~4d at n=68). Control gate 34.6d vs the
split-kinematics Gamma oracle; per-subcone integrand identity worst rel dev
6.9e-40; 20,000-direction fan coverage each direction in exactly one subcone.
Engines: hoist = exact reordering, bit-identical/1.2e-77 vs plain, x3.5–3.8
single-core; c bit-identical to hoist, additional x1.16–1.36; single-core CPU
(n=24/256b) plain 73 us/pt, hoist 19 us/pt, c 15.9 us/pt. Reference values:
I_{3,6}(s*) 60 digits certified (n=72@256 vs n=78@272, rel 1.107e-61; walls
1688 s and 5340 s at 64 processes); point-2 32 digits (n=36@192 vs n=44@224,
rel 3.26e-33). Display policy: papers print <=10 digits.

GATES: reproduce this pattern on any new use — two-run/two-precision
floored-digit agreement (`gates.py`), certificate stage before runs, and the
sanity stage (which reproduces the pilot value at the 1.49e-45 pilot-string
floor, engines c AND plain).

FOOTGUNS:

- mpfr_str value strings carry a TRAILING EXPONENT ("...e2") — truncated
  console displays dropped it and briefly misreported 60.4 as 0.604; parse
  full strings, as the gates/JSONs do.
- Rate fits from low n overpredict high-n digits (~4d at n=68) — budget gap
  legs.
- mp.quad silent 10–15d footgun in semi-infinite oracle forms — the
  production gate never relies on mp.quad.
- 100-digit and X(3,7) commitments: rethink the route and measure a small
  probe first (levers: cone-orbit symmetry, sum-factorization, further C
  work, a DE-transport leg).

CREDIT: the sampler reimplements Algorithm 1 of Michael Borinsky, Anna-Laura
Sattelberger, Bernd Sturmfels and Simon Telen [BSST], itself built on Borinsky's
tropical Monte Carlo [Bor20] (see also feyntrop by Borinsky, Munch and Tellander
[feyntrop]); stratification and tilting are our additions. Stringy canonical forms /
Grassmannian string integrals after Arkani-Hamed, He and Lam [AHL]; CEGM amplitudes
after Cachazo, Early, Guevara and Mizera [CEGM]; positive parametrization and the
split-kinematics oracle after Bruno Giménez Umbert and Bernd Sturmfels [GUS]. Cone
decomposition, Gauss–Jacobi hoisting and the C-MPFR kernel are in-house; hulls by
SciPy/Qhull. Bracketed keys resolve in REFERENCES.md at the repository root.

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).

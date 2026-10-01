# ERAS — dependency-tamed certified enclosures

Tool page: https://www.bootloops.ai/tools/eras.html

KIND: package (`eras.py` the enclosure substrate, `verify_eras_adversarial.py`
the consumer-pluggable adversarial battery, `QUARTET_COLLAPSE.json` the pinned
318-count reference data, shared with `tools/posq`). ERAS = **E**xact
**R**emainder-**A**ware **S**eries.

PURPOSE: certified enclosure forms in a shared parameter `p` for arb/acb
integrands where `p` enters every factor — the regime where naive ball
evaluation explodes by interval dependency. Measured on the reference
integrand shipped here: naive-ball enclosure relative width 8.84e22 at p-ball
radius 1e-4, on an object whose true sensitivity is O(1e2–1e3). The centered
forms restore usable widths: `F(p0) + remainder-from-derivative-on-ball`
(mean-value at first order; order-K Taylor–Lagrange when first order is
rigorous but too wide — measured, K=1 misses practical width targets, and
order-K with panel half-width <= ~3e-3 is the working configuration). The
derivative comes from forward-mode dual numbers over the same assembly (~2x
the value-only cost, measured), the order-K coefficients from truncated series
arithmetic with ball coefficients, so every form is rigorous, not heuristic.

USE-WHEN:

- A shared parameter (a rate inside hundreds of per-observation terms is the
  canonical case) makes naive interval evaluation vacuously wide and you need
  a certified enclosure with usable width instead.
- You need forward-mode `d/dp` duals or order-K ball-coefficient series over
  flint/arb without hand-rolling the product-rule bookkeeping (`D`, `S`,
  `mv_form`, `taylor_encl`, `taylor_ladder` are the reusable pieces).
- You are building ANY adaptation of these enclosure forms — including a
  multivariate one — and need the acceptance battery that assumes the
  adaptation is buggy or adaptive until proven otherwise.

NOT-FOR: general ball-arithmetic instruments (that is `baller`, which consumes
this package rather than absorbing it); out-of-the-box multivariate enclosure
forms (the substrate is univariate; the battery accepts and judges
consumer-built vector adaptations, with a worked vector-mode example in
`MultivariateSelfTest`).

INVOKE: `python3 eras.py selftest` (anchors both arms at gate 1e-30,
forward-mode vs central FD at gate 1e-20, planted-error controls);
`python3 eras.py blowup` / `taylork` (the naive-vs-centered tables with
containment checks) / `timeeval`. Battery:
`python3 verify_eras_adversarial.py all` from scratch space (~50 s; results
JSON go under `$TMPDIR`, override `ERAS_RESULTS`). Consumers: import the
verifier, implement `enclose`/`point_eval` (optionally `deriv`, `degraded`),
call `run_battery(...)` and CHECK the returned `BatteryResult.ok` — exiting 0
without checking it is a consumer bug. The battery must be rerun against any
adaptation before that adaptation is used; a single-seed pass is necessary,
not sufficient — rerun with several `--seed` values.

Dependency: `python-flint` (arb/acb). `QUARTET_COLLAPSE.json` is looked up
beside `eras.py` (`ERAS_BASE` overrides the directory); a missing data file
raises at import — there is no fallback. Data attribution: the 318 site-pattern
counts in `QUARTET_COLLAPSE.json` are derived from the DravLex v1.0 Dravidian
lexical database (Kolipakam, Jordan, Dunn, Greenhill, Bouckaert, Gray & Verkerk
2018), distributed under CC BY 4.0; see
THIRD_PARTY.md at the repository root.

GATES: the shipped reference battery is the acceptance gate — VERIFY PASS with
82 checks, 19 informative coverage cells, 0 undeclared no-information cells,
59 declared-expected blowup-regression cells, 0 fatals, in under a minute. The
reference numbers are regression-gated every run: naive relwidth at r=1e-4
must reproduce the 8.84e22 class, the working configuration (K=16, r=1e-4)
the 4.7e-3 class, and the 53-bit precision probe the ~2e9 degradation
(x4.3e11). Verified consumer in this repository: posq
(`tools/posq/verify_posq_adaptation.py` runs its kernel engine through this
battery; baller's battery leg L6 drives that same run).

FOOTGUNS:

- A non-finite enclosure is NO-INFORMATION, never a pass; a vacuous enclosure
  (width > 1e3x the sampled midpoint spread) is likewise INDET, never a pass.
- Build input constants at working precision: a 53-bit default-precision
  constant degrades a 4.7e-3 enclosure to ~2e9 (measured; the precision probe
  gates it). `p0` must be exactly representable or its representation error
  absorbed outward.
- `point_eval` tightness is quantitative: max point-ball radius <=
  enclosure_width/1e4 at defaults, or the cell is POINT-PRECISION INDET and
  plants cannot be verified from it.
- Corruptions confined to thin sets are detected with probability scaling in
  points x seeds, not deterministically (measured rates in the verifier
  docstring) — hence the several-seeds law.
- The battery defends against honest-but-buggy and adaptive adaptations, not
  against a consumer who edits the verifier or fakes its report.

See the verifier's docstring for the full consumer protocol, the battery laws
(enclosure commitment ordering, closed-box sampling, vacuousness and
point-precision gates, plant channels, FD and precision probes), and the exit
classes (0 PASS / 1 FAIL / 2 INDET).

CREDIT: the enclosure forms are the classical remedies for interval dependency — Moore's
mean-value form (R. E. Moore, Interval Analysis, Prentice-Hall 1966), the centered forms
of Krawczyk & Neumaier (1985, SIAM J. Numer. Anal. 22:604), the Taylor models of Makino
& Berz (2003, Int. J. Pure Appl. Math. 4:379) and the Taylor forms of Neumaier (2003,
Reliable Computing 9:43) — carried in Arb ball arithmetic (Johansson 2017) via
python-flint; ERAS's contribution is the forward-mode assembly and the adversarial
battery, not the forms. Data: DravLex v1.0 (Kolipakam, Jordan, Dunn, Greenhill,
Bouckaert, Gray & Verkerk 2018, R. Soc. Open Sci. 5:171504; lexibank/dravlex,
doi:10.5281/zenodo.5121580), CC BY 4.0.

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).

# longhand — guide

Longhand: independent numerical ground truth for Feynman integrals, done the slow honest
way. Two routes in one package — (A) an arbitrary-precision evaluator of the
Feynman-parametric representation, ours end to end, and (B) a driver that runs pySecDec's
`disteval` on a compiled sector-decomposition package as a checkpointed, self-checking
chain — so that a value obtained by a fast method (differential equations, auxiliary mass
flow, a fitted ansatz) can be confronted with a number that shares none of its machinery.

NAME: `longhand` (directory `tools/longhand/`; Python import name `longhand` with the
repository's `tools/` on `sys.path`; route B is the member `longhand.disteval`).

Tool page: https://bootloops.ai/tools/longhand.html (tools index: https://www.bootloops.ai/tools/,
numerical-oracle group; the external engine route B drives and route A cross-checks:
https://bootloops.ai/tools/pysecdec.html).

KIND: package (python). Route A, flat modules in this directory: `hiprec.py` (engines),
`feynman.py` (U/F builder + closed-form linear elimination), `parametric.py` (generic
UF-spec pathway), `bench_box.py`, `bench_se3l.py` (benchmarks), `cbc_cache.json`
(vendored generating-vector cache, read-only), `toy_qmc_spec.json` (synthetic 6-D
battery fixture with a hand-provable closed form; also a worked example of the UF-spec
schema), `selftest_parametric.py`. Route B, the member directory `disteval/`:
`pysecdec_point.py` (the command line), `stage_chunks.py`, `assemble_compare.py`,
`gate.py` (the acceptance check), `memory_basis.py`, `PACKAGE_PINS.json`. Worked example
of route B: `examples/ndpent_top/` (the pinned records of one complete run; see
`examples/README.md`). Battery: `selftest.py` (both routes; route B's legs live in
`tests/test_disteval.py`). Manuals: `PARAMETRIC.md` (route A, with the benchmark
numbers), `DISTEVAL.md` (route B, with every flag, refusal and measured wall).

REQUIREMENTS: the common core (`mpmath`, `sympy`, `numpy`) plus `gmpy2` — `pip install
gmpy2` — imported at load by route A's engines and by the battery (a missing gmpy2 is a
named hard failure). `pytest` for route B's battery. GNU time at `/usr/bin/time` (the
`-v` form, Linux) for route B's stages, which record each chunk's wall and peak memory
through it. pySecDec is OPTIONAL and user-installed (`pip install pySecDec`; version
1.6.6 is the one the records were made with): route A needs it only for the
propagator-list front end (`feynman.build_UF`; the UF-spec pathway and both engines run
without it), route B needs it for the disteval stages and for the stand-in form `--dry`
(the compare, the acceptance check, the pins and the whole record replay run without it).
pySecDec is licensed GPL-3.0 (not bundled here; importing it at run time places the
running combination under its terms, not this MIT source).

PURPOSE: a SECOND numerical value for a Feynman integral, independent of the method that
produced the first one, with an error statement that is measured rather than assumed.
- Route A (arbitrary precision, ours): finite (eps^0) scalar integrals — and, through the
  UF-spec pathway, any finite projective parametric integrand built from Symanzik-type
  polynomials with a polynomial numerator — in the deep-Euclidean region (s,t<0, m²>0:
  a smooth positive integrand, no contour deformation). It evaluates
  I = Γ(N−L·d/2) ∫ δ(1−Σx) U^{N−(L+1)d/2} / F^{N−L·d/2} directly in arbitrary precision
  and breaks the ~8-9 digit double-precision ceiling of pySecDec disteval / FIESTA double
  mode: 40+ digits when the effective dimension after exact reduction is ≤3 (fixed
  Gauss-Legendre product rule, spectral), ~17 digits at dimension ≤2 by nested tanh-sinh,
  and a few digits with a statistical error bar at irreducible dimension 4-12 (CBC rank-1
  lattice QMC with Korobov periodization and random shifts, multi-process).
- Route B (double precision, pySecDec's numerics, our staging): any integral pySecDec has
  compiled into a `loop_package` + `make disteval`, all regulator orders the package
  carries, at one kinematic point. Three disteval stages (two lattice sizes with standard
  rank-1 lattices, then the median-lattice second integrator at the target precision),
  each split into chunks of sectors that are resumable checkpoints, then a per-order
  pairwise compare and a conservative acceptance check whose digit count is a minimum
  over the B-C pair and the B, C error bars — never a lattice setting — with a
  PLANTED-FAIL control (the compare must see a 1e-3 shift) and a FOREIGN control (coarse
  independent evaluations of the same library must agree with stage A). Established on
  Euclidean points with a package compiled without contour deformation (the worked
  example); nothing in the driver is specific to that, but no deformed package has been
  run through it.

USE-WHEN:
- You need an independent oracle value for a FINITE Euclidean integral where the
  auxiliary-mass-flow route is hard or its eps→0 extrapolation is precision-limited, and
  you want more digits than double precision gives → route A. Effective dimension after
  analytic reduction ≤3 → `integrate_gauss_product` (40-digit box references cross-
  validated to 46+ digits by two independent reductions at dps=45/60); ≤2 →
  `integrate_tanhsinh` (~17 d); irreducible 4-7 (up to 12 measured) → `integrate_qmc`
  with its median-of-shifts error estimate (~3 d on the 6-D benchmark).
- You have ANY finite projective parametric integrand in the UF-spec schema (polynomials
  U0/U1/F0/F1/... + numerator terms — numerators via the Mellin/derivative trick
  upstream) → `parametric.integrate_spec` feeds the QMC engine directly (validated on
  12-dimensional irreducible 4-loop integrands as an independent 2-digit cross-check).
- The integral has poles, or many sectors, or you want every eps order at once, and
  double precision is enough → route B on a compiled pySecDec package: it turns one long
  disteval call into resumable chunks, gives you a second integrator on the same kernels,
  and states digits by object with controls.
- You are about to quote a disteval number and want to know how many of its digits two
  lattice constructions actually agree on → route B's compare + check over your three
  stage results.

NOT-FOR: route A — eps poles, non-Euclidean kinematics, thresholds; inexact kinematic
inputs (see FOOTGUNS); a ≥30-digit value at irreducible dimension ≥6 (QMC error decays
polynomially there: measured ~N^-0.5..-1 at 12-D, rel ~1e-2 class in ~10 min-class wall
— a genuine independent 2-digit cross-check, NOT a 10-digit pin; a ≥30 d value there must
come from an analytic dimension-collapse route, e.g. a nested self-energy + low-dim outer
integral, or a dispersive reduction). Route B — digits beyond disteval's double-precision
ceiling (~8-9); Minkowski points unless your package was compiled with contour
deformation (untested here); building the package (that is pySecDec's `loop_package` +
`make disteval`, done by you beforehand; the compiled package is a build of its host and
never ships with this repository).

INVOKE (route A, Python; `tools/` on `sys.path`):
`from longhand import box, box_crosschart` (the 40-digit one-loop box oracle and its
independent validation route; pass EXACT kinematics: int/Fraction/decimal string —
`box(-1, Fraction(-1, 3), 1, 2, dps=35)`); `from longhand import
integrate_gauss_product, integrate_tanhsinh, integrate_qmc` (the engines;
`integrate_qmc(factory, factory_arg, dim, N=64007, n_shifts=16, korobov_p=3, dps=40,
nproc=None, seed=1)` returns `{'value','error','digits','shifts','N','n_shifts','gvec'}`);
`from longhand import integrate_spec, check_spec, make_eval, load_spec, box_spec`
(the UF-spec pathway: `integrate_spec(json_path, family=None, N=8009, n_shifts=8,
korobov_p=3, dps=25, nproc=None)`); `from longhand.feynman import build_UF,
reduce_linear, reduce_linear_full, reduce_one_quadratic, build_integrand_expr,
assemble_integrand_expr, make_factory` (sympy; `build_UF` needs pySecDec). The same
modules answer to their bare names with `tools/longhand` itself on `sys.path` (`from
hiprec import integrate_qmc`, `from bench_box import box`), which is how the benches and
the battery import them. Benches: `python3 bench_box.py [dps]`, `python3 bench_se3l.py
--N 64007 --shifts 24 --p 3 --dps 40 --nproc 24 --out result.json` (run inside
`tools/longhand`). Any pool must start under `if __name__ == '__main__'` on platforms
whose multiprocessing start method is spawn (macOS): the package's own entry points do.

INVOKE (route B, command line; `cd tools/longhand/disteval` or give the path):
`python3 pysecdec_point.py --disteval-dir PKG/disteval --stage A --workers W --work WORK`
then `--stage B`, `--stage C` over the same `--work` (B and C read their per-chunk
`epsabs` from `WORK/stage_A/RESULT_A.json`); `--stage compare --results
A=WORK/stage_A/RESULT_A.json B=... C=... --out WORK/COMPARE.json`; `--stage check --work
WORK --compare WORK/COMPARE.json --leaf-stamp LEAF.stamp [--smoke NAME=smoke.json ...]
--out WORK/GATE.json` (`--stage gate` is the older spelling of `check`, kept working);
`--stage check --example` replays the worked example's records and must reproduce them
(`--fixtures` is the older spelling); `--planted-fail-control`; `--basis W` (the memory
figure to set on the run's cgroup); `--check-package --disteval-dir D [--pins FILE]`;
`--emit-pins OUT --disteval-dir D --family NAME [--generate-receipt FILE]`. Every line
takes `--family NAME` (default `ndpent_top`, the example's family; `--name` is the alias)
and `--pins FILE` (default `disteval/PACKAGE_PINS.json`, the example's). Common flags:
`--point 'name=value ...'` (every real parameter of the package exactly once; default the
example's point), `--workers W`, `--points`, `--shifts`, `--epsrel`, `--epsabs` (floor),
`--epsabs-from`, `--lattice-candidates`, `--chunk-sectors`, `--resume`, `--max-chunks N`,
`--dry [--dry-result FILE]`, `--unpinned-package`. From Python: `from longhand.disteval
import gate, stage_chunks, assemble_compare, memory_basis` (`check` is an alias of
`gate`; `gate.run_gate(results, rcs, smokes, stamp, out_path, ...)`). Full flag, refusal
and exit-code tables: `DISTEVAL.md`.

INPUTS:
- Route A, propagator list: pySecDec `LoopIntegralFromPropagators` conventions and
  normalization (values compare digit-for-digit with disteval): `build_UF(propagators,
  loop_momenta, external_momenta, replacement_rules)` → `(U, F, xs, L, Nprop)`.
- Route A, UF-spec (a dict, or a JSON file holding one spec or `{family: spec}`):
  `{"name": str, "n_den": E, "polys": {tag: [[monomial_exponents (E ints), [num, den]],
  ...]}, "numerator": {"gamma": int, "terms": [[[cn, cd], [[tag, power], ...]], ...]}}`,
  meaning I = gamma · Σ_t (cn/cd) · T[Π_k P_tag^power] with T[.] the projective simplex
  integral. `check_spec` enforces: every polynomial homogeneous with monomials of length
  E, every numerator term of total projective degree −E (else AssertionError naming the
  offender). Evaluation is in the Cheng-Wu chart x_last=1 over [0,inf)^(E−1); rational
  coefficients are exact. `toy_qmc_spec.json` and `parametric.box_spec()` are worked
  examples. Optional `"exact": [num, den]` is read by the battery only.
- Route A, kinematics: a deep-Euclidean point as EXACT numbers; `dps`, lattice size `N`
  (prime; vectors for N=8009 at dim 3 and 6, N=64007 at dim 12 and N=101 at dim 5 are
  vendored, any other (N, dim) is constructed once — minutes at N=64007 — and cached in
  the user cache), `n_shifts`, `korobov_p`, `nproc`.
- Route B: a compiled package's `disteval/` directory (`NAME.json`, `NAME_integral.json`,
  `NAME_integral.so`, `builtin.so`, `coefficients/`) matching `--pins` (or run
  `--unpinned-package`, labelled); the point; per stage the lattice flags above; for the
  check, the three `RESULT_<S>.json`, optionally FOREIGN smokes (`disteval --format=json`
  stdout of a coarse run of the same `.so`) and the cgroup's end-of-run stamp
  (`key=value` lines: memory.peak, memory.max, oom_kill, wall_s, ...).

OUTPUTS: route A — an mpmath value (`box`, the Gauss and tanh-sinh engines) or a dict
with a gmpy2 value, the standard error of the shift mean, a digit estimate and the
per-shift estimates (QMC); `bench_se3l.py --out FILE.json` writes a JSON record of a
run. Route B — under `--work`: `stage_<S>/chunks/chunk_NNN/{result.json, disteval.log,
time.txt, rc, DONE}`, `stage_<S>/{CHUNKS_<S>.json, SETTINGS_<S>.json, RESULT_<S>.json,
STAGE_DONE, stage.rc}`, `progress.jsonl` (one line per chunk: wall, maxrss, disteval's
statistics), `COMPARE.json`, `GATE.json` (the acceptance check: `gate` PASS/FAIL, the
per-order digit table, the controls, the not-established list, a PRODUCER block); every
JSON written carries a `date -u` stamp. Exit codes 0 PASS / 1 FAIL by name / 2 usage or
refusal by name / 3 pin mismatch / 4 a pinned file missing.

ENV: `LONGHAND_CBC_CACHE` (older name `HIPREC_CBC_CACHE`, still honored) — the user
cache file for newly constructed CBC generating vectors; default
`$XDG_CACHE_HOME/longhand/cbc_cache.json` (else `~/.cache/...`); the cache written by
this module's earlier packaging (`.../hiprec-sectordecomp/cbc_cache.json`) is still read.
`LONGHAND_DISTEVAL_PACKAGE_DIR` (older name `PYSECDEC_POINT_PACKAGE_DIR`) — if set to a
compiled example package's `disteval/` dir, the battery also verifies it against the pins
on this host (otherwise that leg skips by name: the package is not a fixture).
`PYTHONUSERBASE` is carried into route B's subprocesses so a scratch `HOME` never hides a
user-site pySecDec. gmpy2 precision and `mp.dps` are set per call.

GATES (acceptance criteria the battery enforces): route A — dimension-reduce FIRST
(F-linear parameters integrated out exactly in the Cheng-Wu chart, iterated by
`reduce_linear_full`: 1/(A B), then the partial-fraction log form — valid because the
coefficient polynomials are >0 in the Euclidean region); two-dps monotone convergence
before quoting digits; the box at dps=35 ≥30 d vs the 40-d cross-validated reference,
the independent cross-chart reduction ≥18 d at dps=20, and a truncated-input control
that MUST land in the 12-22 d band (proving the inexact-input failure mode is
detectable); QMC values within 5σ of exact references and ≥3 d. Route B — the check over
the example's records reproduces the record by name (21 checks: digit class 7, per-order
9/7/7/7/7, zero B-C lines over 3σ, A-consistency MARGINAL at max σ 3.0519681777552745,
PLANTED-FAIL raised at σ 9.99707344790171 vs A, FOREIGN max σ 1.717724232013687 and
2.9485284820009676); a planted digit re-pinned FAILS the check by name, un-re-pinned is a
pin mismatch; a shifted smoke fails FOREIGN; no smoke is never a pass; every refusal
asserts its named line AND its exit code, never a traceback.

BATTERY (SELF-CONTAINED; run from the package directory): `python3 selftest.py` → route
A's FAST code-check tier (single-process; eleven legs with checked answers on every code
path: the analytic loop-integration kernel on both discriminant branches and the
near-double-root series branch vs direct quadrature; the dimension-reduced box through
the Gauss-product engine at dps=15 ≥12 d; the box DEPTH legs above; tanh-sinh on a 2-D
closed form; `check_spec` positive check + two negative controls; `make_eval` vs an independent
exact-Fraction evaluation; serial CBC-QMC on the 3-D box spec AND the 6-D fixture (≥3 d,
≤5σ vs the exact values); `reduce_linear` / `reduce_linear_full` / `reduce_one_quadratic`
vs exact values and the numeric kernel; the general-(a,b) fallback checked exactly at a
rational point; the pySecDec U/F front end reproducing the pinned box-spec F exactly —
named SKIP when pySecDec is absent), then route B's pytest battery (`tests/`, 32 legs:
the record replay, the compare rule on a synthetic pair, the chunk split on the example's
package JSON, the per-chunk epsabs rule, the memory basis, the refusals, the controls,
the family / pins / alias legs; the two stand-in stage legs import pySecDec and SKIP BY
NAME without it; the on-host package leg skips unless `LONGHAND_DISTEVAL_PACKAGE_DIR` is
set). Measured on a Linux server: 63 s with pySecDec 1.6.6 installed (route A 5 s; route
B 31 passed, 1 skipped, 57 s — the two stand-in legs dominate), 13 s without pySecDec
(route B 29 passed, 3 skipped by name). rc=0. `python3 selftest.py --full` → route A's
heavy tier instead (multi-process: box Gauss dps=30 ≥25 d vs the pinned reference,
measured ~33 d; the 6-D QMC farm on `toy_qmc_spec.json`, N=8009, 16 shifts, 8 worker
processes, within its error bar of the exact 13/1080 — every fixture term has a
Dirichlet closed form, derivation in the fixture's note field; and the generic
pathway's box control through an 8-process pool, ≥4 d ≤5σ, = `python3
selftest_parametric.py`) plus route B as before. `--route parametric|disteval` runs one
route's legs. All entries are spawn-safe (pools start only under `if __name__ ==
'__main__'`). Nothing is written into the package tree (pytest runs with `-p
no:cacheprovider`; temporaries go to the system temp dir).

FOOTGUNS:
- (A) INEXACT KINEMATIC INPUTS silently cap every engine: a constant built at low
  ambient precision (e.g. mpf(-1)/3 at the mpmath default dps) is a different rational
  number at the 1e-16 level, and the spectral rules then converge — to full requested
  depth, two-dps checks included — on the integral of THAT number. Two implementations
  sharing one truncated input "cross-validate" each other down to the same wrong tail.
  Pass int/Fraction/decimal strings; `bench_box` converts all inputs at working
  precision, and the battery's truncated-input control proves the failure is detectable.
- (A) A Gauss-Legendre node cache keyed without dps produces a stale-node false plateau
  — the engine keys its cache on (nodes, dps); keep it that way in derived code.
- (A) The Korobov periodization + CBC lattice is honest about its error bar — do not
  quote QMC central values past the bar. Korobov p=3 (the default) is measured 3.7x WORSE
  than p=1 at 12-D on heavy-tailed integrands (rel 1.5e-1 vs 4.2e-2 at identical cost):
  at high dimension sweep p at pilot size before scaling.
- (A) HIGH-DIM BIAS (measured at 12-D): at small N the shift-spread bar can UNDERCOVER —
  the lattice misses integrable face-spike mass (F0→0 corners) in a way correlated
  across shifts (an N=8009 p=1 pilot read 185.0±7.7 vs true ~223.9, a 5σ-low bias the
  bar did not see). Scale N and check drift against an independent route before trusting
  the bar at dim >> 6.
- (A) `bench_se3l.py` caches its one-time symbolic reduction in `se3l_reduced.json`
  beside the script (needs pySecDec + sympy to build); delete it when changing the
  propagator list.
- (B) A lattice setting is never a digit count: digits come from the compare of two
  constructions and their bars, per order; stage A alone establishes nothing.
- (B) `--unpinned-package` makes the example's smokes meaningless as a FOREIGN control
  (different library): give `--smoke NAME=smoke.json` of YOUR library or the check fails
  on `foreign_ok` by design. A smoke without `sums.NAME` is a failing control, not a
  skipped one.
- (B) `--dry` still imports pySecDec (it is a library-presence rehearsal); without
  pySecDec the stand-in stage cannot run — the record replay (`--stage check --example`)
  can.
- (B) Chunk directories hard-link the `.so` files (symlink fallback): keep `--work` on
  the same filesystem as the package for the cheap form. There is no clock anywhere: a
  stage runs until its chunks land or it is stopped by pid, and `--resume` continues; do
  not wrap it in a timeout and call the result a measurement.
- (B) The tool sets no memory or CPU limit itself: `--basis W` prints the figure
  (page_up(1.3 x (driver + W x per_worker)) from measured maxrss; the nominal per-worker
  form alone would OOM at spawn) and creating the cgroup with it is the caller's job, as
  is reading it back: hand the end-of-run readback to the check with `--leaf-stamp FILE`
  (`key=value` lines; a non-zero `oom_kill` fails the check). The example's
  `FENCE_AT_SPAWN_record.json` is the readback the record run's wrapper made at spawn.
- (B) GNU time must be `/usr/bin/time` with `-v` (Linux); the BSD `time` on macOS has no
  `-v`, so the stages are Linux-only as shipped (the record replay and the compare are
  not).

RELATED: the AMFlow.cpp fork in the sibling repository amflow-cpp and
`tools/amflow-kit` (the primary oracle route A cross-checks; amflow-kit's `memfence` is
the cgroup helper a caller can use around route B); `tools/nestor` (singularity-subtracted
dispersion quadrature — the analytic dimension-collapse route when QMC cannot reach the
digits); `tools/baller` (its hygiene member `dps_lint` for precision hygiene on mpmath
code); `tools/subtropica` (whose verification harness lists this evaluator among its
oracle slots).

CREDIT: route A's U/F extraction and comparison targets come from pySecDec by Borowka,
Heinrich, Jahn, Jones, Kerner, Schlenk and Zirke [pySecDec] and its disteval evaluator by
Heinrich, Jones, Kerner, Magerya, Olsson and Schlenk [disteval] (sector decomposition:
Binoth and Heinrich [BH]); rank-1 lattice QMC for Feynman integrals after Li, Wang, Yan
and Zhao [LWYZ] and Borowka et al. [pySecDecQMC], generating vectors by the
component-by-component construction of Nuyens and Cools [NC] with Korobov periodization;
FIESTA (Smirnov, Shapurov, Vysotsky [FIESTA5]) as second cross-check. Route B orchestrates
pySecDec and disteval [pySecDec, disteval]; its lattice rules are those of Borowka et al.
[pySecDecQMC] and the median-lattice construction of Goda and L'Ecuyer [GL]; all of route
B's numerical work is pySecDec's, the staging, compare and check are ours. Bracketed keys
resolve in REFERENCES.md at the repository root.

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).

# frobenius-boundary — guide

Tool page: see the tools index, https://bootloops.ai/tools/ ("Frobenius" row in the
periods/arithmetic-geometry group).

KIND: package (`frobenius_boundary/` python package: core.py, poly_de.py, de_core.py,
infinity.py, classify_mod.py, branches_mod.py, kira_parse.py, fast_lcd.py, io_util.py,
cli.py, selftest.py, config_example.json) + two member dirs: `asdminer/` and
`deform_transport/` (see MEMBER sections below)

PURPOSE: Family-agnostic Frobenius-branch boundary machinery at var=∞: for a
single-variable transport anchored at infinity of a mass-like variable where the
solution is NOT analytic there (cut blocks whose vacuum boundary vanishes), pins
which Frobenius branches can carry the physical solution BEFORE
any quadrature/oracle work (reference family: 30→27 allowed local exponent classes = 22
survivors + 5 half-integer FLAGGED + 3 killed; a naive "30→12" eig(D1) reading is
WRONG as a counting claim).

USE-WHEN:
- Cut-block boundary at var=∞ where the naive vacuum boundary VANISHES on the maximal
  cut; which branches carry the physical solution.
- Any transport pipeline needing the exponent spectrum, ε-strata survivor cut, and branch
  series Φ(Mstart) at an infinity anchor.

NOT-FOR: finite-point Frobenius landing (that is wayfinder.frobenius / frob_scalar);
analytic-at-∞ boundaries (plain vacuum seeds suffice).

COMPANION: `wayfinder.epslimit` (tools/wayfinder/epslimit.py; formerly the standalone
eps-extrapolator) — `from wayfinder.epslimit import richardson_boundary`;
`richardson_boundary(S,Ms,tail_exponents)`: boundary-of-disk series summation at the
radius of convergence (√-type singularity → half-integer M^{−α} tails eliminated);
37.7d measured at a lattice-Green-function threshold point. Use for summing a
Frobenius/boundary series AT its convergence-disk edge.

MEMBER: `asdminer/` — Atkin–Swinnerton-Dyer congruence miner: the series-side
arithmetic complement.  Fed an exact integer series, it solves every three-term
row a(np) − γ_p·a(n) + p^(k−1)·a(n/p) ≡ 0 (mod p^s) per prime, CRT-merges with
witness extraction, and emits a script-verdicted table (CONSISTENT /
INCONSISTENT with named witnesses / UNDERDETERMINED / P-DENOMINATOR-WALL —
walls are verdicts, never silent garbage), with Weil-box integer candidates,
Eisenstein pins and a mechanical γ classification.  Validated against three
families of proven congruence laws (Honda 16/16 vs exact point counts; Hecke
eigenvalues 5/5; Beukers congruences 12/12) with fault and null controls
caught.  GAUGE LESSON (load-bearing; it bit the builders first): a holomorphic
period at a family-degenerate point mines to the TRIVIAL unit-root class — a
CONSISTENT verdict there looks like "modularity found" but is the unit root
everything has; `classify_gamma()` separates the trivial and genuinely modular
classes mechanically, and the member's first production run was a decisive
witnessed null (0/33 mines consistent).  Files: `asdlib.py` (core, stdlib-only),
`run_control.py` (positive-control battery, gates G1–G6; needs `gp`/PARI on
PATH; ~6 s; the member's fast battery leg), `run_l6_gauge.py` /
`run_l6_mine.py` / `run_l6_topup.py` (operator-driven miners; need
FROB_L6_OPERATOR, see ENV), `fixtures/` (reference gp comparator script+log
pairs).  Outputs land in `asdminer/work/` (override: FROB_ASDMINER_WORK).

MEMBER: `deform_transport/` — the matrix wing: where the trace/series routes
read arithmetic data a prime at a time, this member transports the FULL matrix
of Frobenius along a deformation and lands the complete local factor — the
whole integer characteristic polynomial, certified by functional-equation and
weight bounds and stable under two independent truncations.  The step up is
measured, not asserted: on the rank-6 Sym^5 control every trace-grade route
fails outright while the matrix transport lands all six eigenvalues exactly at
both test primes, and a fault planted to be invisible at rank 2
(`battery/fault_harness.py`, the historical formal-dual index bug) is caught
by the rank-6 instrument on both routes.  Its negative control travels with
it: a measured refusal boundary — two independent recipes both ending in a
licensed, receipted refusal — reproduced (scaled) in the battery so every
consumer watches the instrument refuse correctly before trusting its
acceptances.  Layout: `code/` (deform_kit.py matrix engine; lfact_kit.py
single-series wing; vendored engines r2_kit.py / f2_lib.py / emit_receipt.py;
build_sym5.py operator builder; run_cal_a/b/d/d2 calibration chain;
preflight_invariance.py; run_license/license2; emit_wall),
`battery/` (run_battery_deform.py stages S0–S5 + fixtures/deform_controls.json
checked-answer comparators).  Deps: python-flint + gmpy2 (guarded — named
refusal at first use, not at import); numpy/scipy only inside the sym5 series
control.  Battery: `python3 battery/run_battery_deform.py --stages S0,S1`
(fast tier, seconds-to-minutes class) — full `--stages S0,S1,S2,S3,S4` is
hour-class (S3/S4 measured ~30–40 min each) and S4 needs
FROB_L6_OPERATOR.  The asdminer member above is its congruence-side partner:
the cheap first strike before the matrix machinery is paid for.

INVOKE: config-driven `cli.py` (see frobenius_boundary/config_example.json for the
schema); API pipeline: `build_poly_DE` (exact D(x)·dM/dx=G(x)·M from a per-point Kira
table; scalar LCD at rational d + DFT fit of G, fit relerr ~1e-77 reference) →
`spectrum` (exact Laurent of x·A at ∞, D1=R_Δ−diag(L·d/2+ai), exponents = eig(D1)) →
`classify` (pair spectra across ≥2 ε, rationalize λ=a+b·ε, Jordan via SVD nullity,
trace rule tr D1 rational-linear in ε = build-integrity gate) → `exclude_strata`
(survivor argument; killed branches hardened by no-rational-within-1e-25 (q≤1e6) scan
at 2 ε; complex pairs by |Im|>tol) → `branch_series` / `phi_matrix(Mstart)` (D_p
recursion seeded by eigvectors → N×J branch matrix for the anchor solve
R_rows·Φ·c=O_rows). Selftest: `python3 -m frobenius_boundary.selftest` (run from this
directory) = the reference-family regression, 73 checks; `FROB_CONFIG=<path>` env override
for the config. NOTE: the reduction table the regression targets belong to does not
ship; fill in a copy of `config_example.json` (`kira_targets_m`/`masters` → your own
per-point Kira table + preferred_masters) and point FROB_CONFIG at it — without a
table the selftest refuses loudly (exit 3). Wall ~15-25 min.

INPUTS: per-point Kira reduction table (kira_parse.py); ε list (multi-ε classification
wants ≥2, the 4-ε pattern is the validated one); convention M_k(var) ~
var^(L·d/2 + ai_k + λ_j), ai = (#ISP powers) − (Σ dots), ε=(4−d)/2, L=n_loops;
N-frame N=diag(var^−(L·d/2+ai))·M ~ var^λ.

OUTPUTS: exact poly-DE (D,G); exponent spectrum + rational classification; Jordan data
+ trace-rule check; survivor strata list; branch series + Φ anchor matrix; diagnostics
fwd_reduced / orbit_multiset / refine_rfinal / LOG_BRANCH_REQUIRED / reduce_failed.

ENV: FROBENIUS_POLYDE_CACHE=<dir> disk-caches the fitted poly-DE;
FROBENIUS_SAMPLE_NCPU=<n> fork-parallel D*A sampling. Both default OFF; numerics
unchanged when on.  Member env (all refuse loudly when a needed one is unset):
FROB_L6_OPERATOR=<json> — the order-6 reference operator ({"coeffs": [[...],...]},
7 rows; NOT distributed) consumed by the asdminer L6 drivers, deform op_L6() and
battery stage S4; FROB_ASDMINER_WORK=<dir> — asdminer output dir (default
asdminer/work/); FROB_PRODUCER_LINT=<script> — optional producer-lint hook
(receipt emission skips lint with a printed note when unset);
DEFORM_R2_GRID=<json> + DEFORM_F1_RECEIPTS=<dir> — comparator inputs for the
deform run_cal_d/run_cal_d2 gates (not distributed).

GATES: selftest (73 checks) mandatory after any edit; trace rule =
build-integrity gate; PROVEN vs ASSUMED ledger in README — PROVEN: the DE and D1, the
exponent spectrum + rational/non-rational classification at the given ε set, (reference family)
vacuum-sector emptiness killing integer ε-independent branches. ASSUMED (flagged in
output): survivor (a,b) rationalizations are multi-ε fits at ≤1e-26, not symbolic;
ε-incommensurate kills rest on the single-Laurent-ε argument; strata kill ASSUMES
regions completeness + quasi-unipotency.

FOOTGUNS:
- `kill_integer_eps_indep` ONLY with a per-family vacuum IBP check — never blanket.
  NOTE: config_example.json therefore ships it FALSE, but the selftest's regression
  targets (3 KILLED rows) require TRUE for the reference family — selftest.py forces
  the flag itself, so a filled-in copy of config_example.json reproduces the targets;
  set it deliberately in your own production configs.
- Integer resonances on the backward path: compatible/min-norm SVD solve
  (incompatibility surfaces as refine_rfinal/LOG_BRANCH_REQUIRED; span of Φ
  unchanged). Genuine nilpotent D_0 forward blocks (fit-noise D_{−k} dropped) are
  removed by exact shear reduction → backward recursion in the Fuchsian reduced frame
  + map-back; the old forward-coupled dense lstsq is fdepth>1 fallback ONLY (flagged
  reduce_failed).
- Do not read strata multiplicities off eig(D1) when D_0≠0.

BATTERY: the selftest is fully scripted BUT its input (the reference family's Kira table) is
not included in this repo, and the run walls at ~15-25 min; on a clean checkout it
refuses loudly (exit 3) — that refusal is the expected core-battery outcome.
`selftest_out/` is regenerated by the selftest. Public-runnable today: the API
pipeline on your own per-point Kira table (kira_parse.py reads standard kira2math
output). Member fast legs (runnable on a clean checkout): asdminer —
`python3 asdminer/run_control.py` (gates G1-G6, ~6 s, needs gp/PARI on PATH);
deform_transport — `python3 deform_transport/battery/run_battery_deform.py
--stages S0,S1` (Sym^5 operator rebuild + rank-2 theorem control; needs
python-flint + gmpy2). The full deform battery (--stages S0,S1,S2,S3,S4) is
hour-class and S4 additionally needs FROB_L6_OPERATOR.

RELATED: consumed by wayfinder.boundary_branches (thin delegator, restates the
eig(D1) caveat); upstream = per-point Kira tables; downstream = anchor solves in
transport pipelines.

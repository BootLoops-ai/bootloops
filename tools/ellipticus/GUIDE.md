# Ellipticus — certified evaluator for elliptic multiple polylogarithms and iterated integrals on an elliptic curve

Ellipticus evaluates elliptic multiple polylogarithms (eMPLs) and iterated
integrals on an elliptic curve to a requested precision and returns each value
with a certificate. Given an exact curve or family spec (a quartic or cubic
y^2 = F(x;z) with coefficients rational in z, exact rationals only), a set of
letters, a word, an evaluation point and a number of digits, it returns the
values together with the evidence that they are right to that precision:
proven per-step Taylor tails on the transport engine, proven q-tail bounds on
the Kronecker-Eisenstein engine, and two-precision plus method agreement
throughout. Python import name `ellipticus` (formerly `empl_eval`; the former
EMPL_EVAL_* environment names are still read).

Tool page: https://bootloops.ai/tools/ellipticus.html (tools index:
https://bootloops.ai/tools/). Related: tools/eichler for modular-form periods
and iterated Eisenstein integrals (the Eichler.jl engine); tools/gpl-eval for
ordinary MPLs.

KIND: package (curve.py, extend.py, march.py, qengine.py, api.py, selftest.py,
battery.py; README.md + SUMMARY.md in-tree; member subpackages `ellred/` and
`gmtel/` — see MEMBERS below). LAYOUT: supporting engines ship VENDORED here — `vendor/kronecker.py` + `vendor/itint.py`
(Kronecker-Eisenstein kernels + AGW tau-words), `vendor/rowscripts_empl.py` (BMSW
sunrise rows machinery) and
`vendor/row{09,10,11}_data.json` (held-out engine gate points); `fixtures/` holds
the synthetic demonstration families + their pinned exact systems the selftest
compares against, `fixtures/battery/` the independent-quadrature reference
values for battery legs L1-L3 + the pinned modulus/chart strings for leg L5
(`ellred_controls.json`), and `fixtures/gmtel/L5_theta.txt` the reference L5
operator the gmtel member gates against.

PURPOSE: General certified eMPL / elliptic-iterated-integral evaluator. Given an
exact curve/family spec (quartic or cubic y^2=F(x;z), coefficients rational in z,
exact rationals only), letters, a word, an evaluation point and dps: returns values
with a certificate (proven per-step tails / proven q-tail bounds / two-precision +
method-agreement, per engine).

USE-WHEN:
- eMPL / elliptic-polylog evaluation is slow or ad hoc — window-moment kernels
  I_k(z)=int x^k dx/sqrt(F(x;z)) and outer dlog-letter words on them at 100d+
  certified.
- Weight-1/period layer: real-oval Legendre-K/AGM periods, cycle moments
  oint x^k dx/y, or COMPLEX-pair contour cycle values on an exact quartic.
- Kronecker-Eisenstein g^(n)(z,tau)/omega_k values with a proven tail bound; BMSW
  unequal-mass sunrise eps^0 (E^(0), J^(0)) at any Euclidean t.
- Building an extended inhomogeneous-PF first-order system (moment block +
  algebraic boundary states with rational log-derivative) for ANY squarefree
  cubic/quartic family — no IBP needed.

NOT-FOR: Gamma_1(6) iterated-Eisenstein/modular-form words with ball-certified
enclosures — route to the Eichler package (upgrades/Eichler.jl + tools/gpl-eval's
empl bridge; certified, majorant tails; NOT duplicated here). Equal-mass kite Ebar
words — not implemented here. v1 scope refusals
(loud, fail-closed): dlog(z) letter allowed OUTERMOST only (single-log
tangential-base reg; deeper log letters need the L-algebra — engine q has it for
tau-words); analytic seed bases must be perfect-square or square-x-linear fibres at
z=0 chart; nested letter seeds beyond depth-1 not wired; no z=1-class
gauge-singular endpoint bridge (transport to 1-delta only).

INVOKE:
- import: put this directory's parent (tools/) on sys.path or PYTHONPATH, then
  `from ellipticus import MomentFamily, QuarticCurve, evaluate, qengine`.
- march (engine i): `fam = MomentFamily(coeffs_x_desc_strings_in_z,
  window=(Fraction,Fraction), kmax=2|1, param='z')`; `fam.polyform()` (exact
  system); `fam.eval_words([[a_w,...,a_2,('mom',k)],...], endpoints, dps,
  waypoints=None)` → {endpoint: {kernels, boundary, words}} at wdps=dps+80.
  Letters: exact rationals or mp values from exact data; a=0 outermost = reg dlog.
- periods (engine agm): `cur = QuarticCurve([c4..c0] exact)`;
  `cur.period_oval(dps)`; `cur.cycle_moments(ks, dps)`;
  `cur.cycle_pair_contour((i,j), ks, dps)` for complex-pair cycles.
- q engine: `qengine.g_n_certified(n, z, tau, dps)` → (value, PROVEN bound);
  `qengine.sunrise_E0(t, masses_sq, dps)`, `qengine.sunrise_J(...)`,
  `qengine.sunrise_equal_mass_control(t, dps)`; `qengine.tau_word_frame` +
  `tau_words`.
- self-gate wrapper: `evaluate(verb_callable, dps, extra=20)` → dict(value,
  shared_digits, lo, hi).
- selftest: `python3 selftest.py` from this directory (or `PYTHONPATH=<tools
  dir> python3 -m ellipticus.selftest`); battery: `python3 battery.py` here, or
  `PYTHONPATH=<tools dir> python3 -m ellipticus.battery` (add `--pilot` for a
  low-dps timing pass).

INPUTS: EXACT rationals everywhere (Fraction/int/str; floats REFUSED — purity
bar). Family coefficients as sympy-parseable strings in the parameter;
rational-in-z dependence allowed. Window endpoints exact rationals. dps explicit
per call; ambient mp.dps never left modified.

OUTPUTS: mp values at working precision (wdps=dps+80 march; dps+30 curve verbs);
battery receipt written to $ELLIPTICUS_BATTERY_OUT, default
battery_results.json in the package directory (per-leg rows/digits/pass/secs);
march per-step tail asserts raise AssertionError on failure (fail-closed).

ENV: ELLIPTICUS_CACHE (optional seed-series disk cache dir);
ELLIPTICUS_BATTERY_OUT (battery receipt path). ELLIPTICUS_KRON_DIR /
ELLIPTICUS_ROWS_DIR (override the vendor/ engine copies), ELLIPTICUS_W3
(override fixtures/), ELLIPTICUS_ARTIFACTS (root of the reference values for
battery legs L1-L3 + L5; when unset they are read from the shipped
fixtures/battery/, and legs whose files are absent print a named SKIP.
Layout under the root: march_words.json + cycle_controls.json +
contour_controls.json + ellred_controls.json — each file's `comment`/
`_comment` documents its provenance). Former names: each variable is also
read as EMPL_EVAL_CACHE, EMPL_BATTERY_OUT, EMPL_EVAL_KRON_DIR,
EMPL_EVAL_ROWS_DIR, EMPL_EVAL_W3, EMPL_EVAL_ARTIFACTS respectively; when both
are set the ELLIPTICUS_* value wins (`ellipticus.env`). Member env: FAMGEN_L5 /
FAMGEN_TOOLS (gmtel; see MEMBERS).

GATES: selftest 4/4 (S2 = EXACT string match of the generic-built quartic 5x5
system vs the pinned fixtures/system_polyform.json, a snapshot the same call
regenerates; S4 tail-bound must COVER measured truncation error). Battery =
reference-value reproduction end-to-end against INDEPENDENT oracles (direct
tanh-sinh quadrature of the defining integrals / documented letter
conventions, plus a closed-form 2*K(1/4) period — never the code path under
test). Rerun battery after ANY engine edit.
Measured on this copy: selftest S1-S4 ALL PASS (~15 s; S1 80.0d, S3 80.0d, S4
routes 49.5d + equal-mass 50d); battery ALL FIVE legs PASS in ~21 s — L1
march kernels ~70d / words ~47d vs quadrature oracles + EXACT polyform pins
(both families), L2 period 69.9d vs closed form + moments 47d vs quadrature,
L3 contour 47d vs the segment-integral oracle (orientation sign recorded),
L4 sunrise rows 09/10/11 ~100d vs the bundled held-out engine gate points +
equal-mass classical control, L5 ellred quartic F/E/Pi atoms ~60d + cubic
J-moments >=39d vs direct quadrature + 1r2c Carlson-route agreement 61d +
J1 derivative certificate 8e-62 (~1 s).

MEMBERS:
- `ellred/` — standard-elliptic (Legendre/Carlson) reduction of one-fold
  elliptic integrals int rat(x) dx / sqrt(cubic|quartic): quartic Legendre
  chart + closed F/E/Pi atoms + J-moment recurrence (`ellred_engine`), cubic
  3-real-root and 1-real-2-complex Legendre moments (`ellred_cubic`), x-space
  reducer to F/E/Pi + algebraic (`ellred_reduce`), analytic 1r2c reducer with
  exact-derivative certification (`ellred_1r2c`), Carlson R_F/R_D/R_J route
  cross-check (`ellred_carlson`), consolidated re-exports + demonstration
  P4_at family (`ellred_kernels`), vendored chart layer (`ell_normal`).
  Invoke: `from ellipticus import ellred` (or `from ellipticus.ellred import
  ellred_engine, ...`). Atom-derivation gate script: `PYTHONPATH=<tools dir>
  python3 -m ellipticus.ellred.derive_atoms`. Battery: leg L5 (above) gates
  chart moduli vs the pinned strings in fixtures/battery/ellred_controls.json
  and every atom class vs direct quadrature; measured ~1 s at dps 45.
- `gmtel/` — certified Gauss-Manin q-telescoper: exact Picard-Fuchs operators
  of K3 pencil periods from a Weierstrass model (exact 2-parameter fiber GM +
  creative telescoping over the base q with an explicit certificate
  row-vector). Route A = geometric telescoper (`famgen_pf_routeA`), route B =
  data route via tools/annihilator pf_from_series (`famgen_pf_routeB`),
  shared exact operator algebra (`famgen_common`). Script-style modules (they
  sys.path themselves); deps: python-flint + numpy (routeA imports flint at
  module top). Env: FAMGEN_L5 = reference L5 operator file (default the
  packaged fixtures/gmtel/L5_theta.txt), FAMGEN_TOOLS = dir holding the
  annihilator engine module (default the sibling tools/annihilator/).
  Member check (GUIDE-documented, not a separate BATTERIES.json row — the
  full control is minutes-scale): `python3
  tools/ellipticus/gmtel/gmtel_selftest.py --pilot` (degeneration smoke,
  measured 53 s) or full (positive control (1,4,16): two routes ==
  reference L5 exactly up to the exhibited x^(-1/2) twist, exact certificate
  at 8/8 rational s points; measured 268 s on this copy).

FOOTGUNS:
- Seed-series generation dominates cold march cost (exact fmpq partial fractions,
  N~wdps/1.7 terms); set ELLIPTICUS_CACHE — warm runs skip it entirely.
- cycle_pair_contour sign convention: principal sqrt at the theta=0 basepoint,
  counterclockwise; a stored value from a different contour/basepoint can differ
  by overall cycle orientation — the battery records the sign explicitly; never
  silently flip.
- Some families are real only inside an exact frontier (endpoint-root
  crossing); transport toward it makes steps shrink (dist-to-sing control) — stop
  short of the frontier, don't push through.
- Letters at a=0 anywhere but outermost raise NotImplementedError by design; the
  tau-word engine (itint L-algebra) is the door for inner log letters.
- The packaged kronecker qbar route near a horizontal lattice line auto-falls back
  to the theta route; g_n_certified then certifies by route agreement (not a
  series bound) — the returned "bound" is the measured route difference there.

CREDIT: Kronecker–Eisenstein kernels and τ-iterated integrals implement the conventions
of Bogner, Müller-Stach and Weinzierl [BMSW] and Adams and Weinzierl [AW]; elliptic
multiple polylogarithms in the sense of Brown and Levin [BL] and Broedel, Duhr, Dulat,
Penante and Tancredi [BDDT, BDDPT]; independent numerics cross-checked against
Walden–Weinzierl's GiNaC implementation [WW] where applicable;
moment_system/transport/seed charts are validated against independent quadrature oracles
in the battery; sunrise closed forms = the rows machinery (vendor/rowscripts_empl.py).
Certification layers are ours. Bracketed keys resolve in REFERENCES.md at the repository
root.

RELATED: tools/eichler + upgrades/Eichler.jl (modular-form periods and the
Gamma_1(6)/iterated-Eisenstein certified engine — route, don't rebuild);
tools/gpl-eval (ordinary MPLs; its empl_bridge.jl fronts Eichler.jl's
iterated.jl); tools/wayfinder (stored-DE marcher — same contract family;
Ellipticus builds its own systems instead of loading stored DEs); famhar
(multi-sheet family evaluator; elliptic slots declared there route here);
nestor (quadrature oracle layer used as independent gates).

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).

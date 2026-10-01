# famhar — guide

Tool page: https://bootloops.ai/tools/famhar.html

KIND: package (python engine `famhar.py` + the ladder calibration family: reference plugin
`ladder_reference.py`, config emitter `build_ladder_config.py`, `families/ladder_L4.json`,
`PROVENANCE.json`, `tests/`; + README)

PURPOSE: Multi-sheet FAMILY evaluator harness over a 2-var dlog connection: a
family = ONE JSON config stamp {connection A, anchor, alphabet/branch loci,
masters<->targets map}; harness supplies certified eval(path_spec, dps) -> value +
du,dv (A.J + prefactor/Jacobian rule) + sheet record + monodromy loop driver. Thin
layer over UNCHANGED wayfinder.

USE-WHEN:
- Need value + ∂u,∂v + sheet/monodromy for a family of masters as a config stamp
  instead of a bespoke evaluator.
- Second-sheet / through-cut evaluations where naive principal-value logs are wrong by
  a disc (G5 class).

LETTERS: affine (const + per-var coeffs) or POLYNOMIAL ({"terms": [[[exps],"p/q"],...]}
canonical-terms schema, total degree >= 2) — polynomial letters get per-segment EXACT
Gaussian-rational composition + mp.polyroots (dlog w(t) = sum_roots dt/(t-r)),
t-root clearance certification against the same clearance_min bar, per-factor
continuous winding; root-finding failure REFUSES the segment; polynomial marches
REFUSE floating waypoints (exact rationals only — evaluate with rational-string
waypoints already qualifies).

NOT-FOR: connection entries eps-independent exact rationals
(eps-dependent stored systems go through wayfinder loaders on their own variable);
kinematic map is 2x2 only; ELLIPTIC families fit the same config slots (Eichler period transport / coalescer pinned-convention projectors) but no elliptic family ships in this release;
higher jets through discs (out of scope here).

INVOKE (public): `import famhar; fam = famhar.Family.from_json("families/<fam>.json")`;
anchor values via config `anchor.values_plugin` (module.func);
`res = fam.evaluate(path_spec, y0, dps)`; loops via `fam.march(closed_loop, J, dps)`.
famhar finds `wayfinder` by inserting its parent tools/ dir on sys.path —
both packages ship in this repo's tools/.

INPUTS: JSON config stamp (schema in README: chart_vars, eps0, masters, letters
(affine or polynomial), dlog connection, anchor point+plugin, kinematics
map/jacobian, targets coeffs+prefactor, clearance_min); ALL numbers exact rational
strings parsed inside workdps; path_spec = piecewise-straight waypoints in the
family CHART.

OUTPUTS: `res["targets"][name]` -> value / d_chart / d_kin (this sheet);
`res["sheet"]` = per-letter continuous winding + continued log(letter) + per-segment
wayfinder diags (trunc_worst recorded); full path retained as homotopy datum.

GATES: two-pass dps law — every quoted digit marched at dps 50 AND 80,
quoted = min(vs-reference, pair agreement). Gate battery pattern per family (ladder:
G0 A.J-vs-FD before any transport trusted, G1 values vs 100d reference values, G2
derivatives, G3 identity loop, G4 branch-locus loop vs INDEPENDENT exact monodromy —
exp(2*pi*i*Res A) is NOT the reference, path-ordering — G5 second sheet, R1/R2
regressions); ladder PASS_ALL 7/7 at 50d class. Run a timed pilot before any dps>=100 march
(wayfinder ~dps^1.7 cost law). The ladder gate battery itself (the record
`evaluator_harness/gates_ladder.py` sha256 752b454378b66386…, GATES_LADDER.json,
gates_ratmarch.py, the 100d reference values) is not included in this repo; the harness
code, the config schema and the calibration family (plugin + config emitter + emitted
config, pinned in PROVENANCE.json) ship whole.

MEMBERS: `famhar.py` (engine) · `ladder_reference.py` (closed-form reference plugin of the
ladder family; `phi_value`, `phi_derivs_chart`, `basis_values`, `anchor_values`,
`default_ell_k`, `ell_sheet`) · `build_ladder_config.py` (emits `families/ladder_L4.json`
from the plugin's basis, exact Fractions, no reference data) · `families/ladder_L4.json` ·
`PROVENANCE.json` · `tests/test_ladder_sheet.py` + `tests/fixtures/ladder_sheet_receipt.json`.

FOOTGUNS:
- Segment with letter-clearance < clearance_min is REFUSED with detour named — never
  march silently (seg-gating law, regression R1); path endpoint ON a locus raises the
  named REFUSAL (zero-endpoint guard in `famhar._delta_arg`).
- Naive principal-value logs across a branch crossing are WRONG by exactly the disc —
  the tracker's arg-unwrap (bisection to |d arg|<pi/2) is the fix.
- Ambient-dps input-conversion trap: config numbers must be exact rational strings;
  harness parses inside workdps for exactly this reason.
- Polynomial-letter marches need EXACT rational waypoints (the letter t-polynomial
  is composed exactly before root-finding) — float waypoints raise the named
  REFUSAL, never a silent approximate compose.
- THE SHEET OF ell in `ladder_reference` (sheet-selected default): cont=None means the
  anchored PRINCIPAL branch ell = principal log(a*b), realized as log(a)+log(b)+2*pi*i*k
  with k = `default_ell_k(a, b, dps)` (k=0 on the Euclidean slice b=conj(a) and at the
  anchor; k=-1 on the real lambda>0 sheet, where a=w and b=wbar are BOTH negative real
  and the naive principal sum carries an extra 2*pi*i — the old default returned
  -383.41…+113.78i-class garbage for a real Phi^(3) = 117.915…). An explicit
  cont={'ell_k': k} stays RELATIVE TO log(a)+log(b) (so ell_k=-1 == the default there,
  ell_k=0 == the old wrong branch); `ell_sheet` reports the k in force.
  `phi_derivs_chart` uses the same default.

BATTERY: `python3 tools/famhar/famhar.py --selftest` — self-contained synthetic
battery (no reference data; seconds). Verdict contract: every leg registers itself
as it opens and the verdict line's executed/skipped counts are read from that
registry, never typed (`legs 6 executed / 0 skipped` on the shipped tree);
`PASS_ALL` rc 0 only at 0 skipped and 0 failed; `FAIL` rc 1 on any failed check
(named); `PASS_EXECUTED_WITH_SKIPS` rc 2 when no check failed but a leg was skipped
— each skip printed as `SKIP: <reason>` under its `[Bn]` banner and the skipped legs
named on the verdict line; no expected-fail legs. The legs: B1 straight
march values + A.J chart/kin derivatives vs closed forms; B2 identity loop (zero
windings, J returns); B3 branch-locus loop (winding = 2*pi, monodromy vs the exact
closed-form continuation — function-level reference, never exp(2*pi*i*Res A)); B4
polynomial letter (exact composition + root-finding: value/derivative/monodromy);
B5 refusal laws (near-locus, through-locus, endpoint-on-locus, float waypoints);
B6 ladder sheet fixture (`tests/fixtures/ladder_sheet_receipt.json`, Phi^(3) at
(u,v)=(1/4,1/10) on the real lambda>0 sheet at dps 50: the record's ell_k=-1 call and
the default agree with the record to its dps-40 cap (bar 38 digits, measured 40.3) and
with each other to >= 48; the wrong-sheet control ell_k=0 reproduces the old branch and
fails by name if the default returned it; Euclidean-slice and anchor k=0 controls) —
SKIPPED BY NAME (`SKIP: absent ladder_reference.py`, verdict
`PASS_EXECUTED_WITH_SKIPS`, rc 2) if `ladder_reference.py` or the fixture is absent.
`python3 -m pytest tools/famhar/tests -q` — the full sheet battery (both record points,
tampered-fixture control, derivatives vs finite differences, two-precision
self-consistency, README/GUIDE membership with a planted control, PROVENANCE pins,
builder drift guard, config load + one-point A.J-vs-FD, and the registered self-test
run twice: on the tree, where the verdict line's counts must equal the number of
`[Bn]` banners printed with the `PASS_ALL` token and rc 0, and on a plugin-removed
copy, where B6 is skipped by name with the distinct token and rc 2). The full gate LADDER consumes
the 100d reference values (not shipped); use it as the pattern for gating your own
family JSON. Measured: ~3.0 s/step at n=65 dps40 (sparse-aware, 156 nonzero entries);
ladder battery 2182 s.

RELATED: wayfinder (transport core, used as-is: transport_fixed_eps, gate, manifest,
stored-system loaders); declared elliptic extensions: Eichler period_matrix, coalescer.

CREDIT: calibration family = the Ussyukina–Davydychev ladder integrals Phi^(1..4) [UD]
in the conventions of Drummond, Duhr, Eden, Heslop, Pennington and Smirnov [DDEHPS];
transport by wayfinder. Bracketed keys resolve in REFERENCES.md at the repository root.
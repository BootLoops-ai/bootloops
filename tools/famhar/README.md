# famhar — multi-sheet family evaluator harness

Certified numerical evaluation for a *family* of master integrals connected by
a two-variable dlog differential equation. A family is one JSON config: the
connection matrix, an anchor point with a values plugin, the letter alphabet
(whose zeros are the branch loci), and a masters→targets map. Given that
config, famhar provides

- `evaluate(path_spec, y0, dps)` — values of the targets at the path endpoint,
  on the sheet selected by the path, with partial derivatives in both the
  chart variables and the kinematic invariants;
- a sheet record — per-letter continuous winding and continued logarithms,
  the abelianized homotopy data of the path;
- `march(closed_loop, y0, dps)` — monodromy: the action of a closed loop on
  the transported basis vector.

famhar is a thin layer over [`wayfinder`](../wayfinder/README.md), which
ships beside it in `tools/` and does the actual Taylor transport.
Single-variable systems stored on disk still route through wayfinder's own
loaders; famhar's job is families given as an in-config connection.

Tool page: <https://bootloops.ai/tools/famhar.html>. Condensed reference card:
[`GUIDE.md`](GUIDE.md).

## When to use it

- You need value + ∂u,∂v + sheet/monodromy for a family of masters, and would
  rather write one JSON config than a bespoke evaluator.
- Second-sheet or through-cut evaluations, where naive principal-value
  logarithms are wrong by exactly the discontinuity.

## Requirements

Python 3 with `mpmath`. famhar locates `wayfinder` by putting its parent
`tools/` directory on `sys.path`, so no install step is needed beyond having
both packages in place as shipped.

## Usage

```python
import famhar, importlib

fam = famhar.Family.from_json("families/my_family.json")

# anchor values come from the plugin named in the config
mod, fn = fam.cfg["anchor"]["values_plugin"].rsplit(".", 1)
anchor = fam.parse_point(fam.cfg["anchor"]["point"], dps + 15)
y0 = getattr(importlib.import_module(mod), fn)(anchor, dps + 15, fam.masters)

res = fam.evaluate(path_spec, y0, dps)   # path_spec = [anchor_pt, ..., end_pt]
res["targets"]["Phi4"]["value"]   # value at the end point, on this sheet
res["targets"]["Phi4"]["d_chart"] # d/da, d/db  (A.J + prefactor rule)
res["targets"]["Phi4"]["d_kin"]   # d/du, d/dv  (config's 2x2 Jacobian)
res["sheet"]                      # per-letter winding + continued logs
                                  # + per-segment transport diagnostics
```

Monodromy loops go through `fam.march(closed_loop, y0, dps)` with a closed
`path_spec`.

## How evaluation works

- **Paths and sheets.** The path is given in the family chart as
  piecewise-straight waypoints; the sheet is determined by the path's homotopy
  class. The harness returns the abelianized record — per-letter continuous
  winding plus the continued log of each letter — which is exact bookkeeping
  for dlog alphabets whose branch loci are letter zeros. The full path stays
  in the result as the complete homotopy datum.
- **Certified marches.** Each segment is pulled back to an inline DE system in
  t on [0,1], with exact partial-fraction `A_series` and declared
  `singular_points`. Affine letters contribute their pole in closed form;
  polynomial letters contribute one unit pole per root of their exactly
  composed t-polynomial (see the letters section below). wayfinder's
  complex-pole clearance and geometric tail-bound acceptance do the
  stepping; the worst truncation certificate per segment is kept in the
  diagnostics.
- **Fail-closed segments.** A segment whose letter clearance falls below the
  config's `clearance_min` is refused with the detour named — nothing is ever
  marched silently across an uncertified region. A path endpoint sitting
  exactly on a branch locus raises a named refusal rather than dividing by
  zero.
- **Branch crossings.** The sheet tracker unwraps principal arguments to the
  accumulated angle (bisecting to |Δarg| < π/2 per sub-step), never
  re-principalizing through a crossing. Naive principal-value evaluation
  across a cut is wrong by exactly the discontinuity; the tracker is the fix.
- **Monodromy references.** When gating a loop against an exact reference,
  derive the reference from the function representation itself. Because of
  path-ordering, exp(2πi·Res A) is *not* the monodromy of the transported
  system in general (see the note in `famhar.py`).
- **eps.** Marches are fixed-eps (no Laurent window; put an eps-fan on top if
  a family needs the Laurent expansion). Frobenius landings, if used, inherit
  wayfinder's real-eps0-only rule.

## Config schema (a family = one JSON file)

```jsonc
{
 "name": "my_family",
 "chart_vars": ["a","b"],          // the connection's variables
 "eps0": "0",                      // fixed-eps march value (rational string)
 "n_masters": 65,
 "masters": [ {"type":"lA","r":0,"s":1}, ... ],   // basis bookkeeping (free-form dicts)
 "letters": [ {"name":"a","const":"0","a":"1","b":"0"},          // affine letter
              {"name":"w","terms":[[[0,0],"1"],[[1,1],"-1"]]},   // polynomial letter (1 - a*b)
              ... ],
 "connection": { "a": {"i,j": "p/q", ...}, ... }, // dJ = [sum M_letter dlog letter] J
 "anchor": { "point": {"a":"1/8","b":"1/9"},
             "values_plugin": "module.func" },     // -> J(anchor) at dps
 "kinematics": { "vars": ["u","v"], "map": {...}, "jacobian": {...} },
 "targets": { "Phi4": { "coeffs": {"masterindex":"p/q", ...},
                        "prefactor": "((a-1)*(b-1))/(a-b)",
                        "dprefactor": {"a": "...", "b": "..."} } },
 "clearance_min": "1/50"
}
```

All numbers are exact rational strings, parsed inside the working precision —
never converted at ambient precision, which would silently truncate them. The
same rule is why floats are refused wherever exactness matters.

### Letters: affine and polynomial

Letters may be affine in the chart variables (the `const`/per-variable
schema) or polynomial (the `terms` schema: a list of
`[[exponent tuple], "p/q"]` terms over the chart variables, the same
canonical-terms format as the rational connection block below). Affine
letters cover dlog/polylog alphabets in a rational chart (e.g. an (a,b)
chart in which an algebraic square root has been rationalized away); a
`terms` letter of total degree ≤ 1 is folded into the affine fast path.

Polynomial letters (total degree ≥ 2) switch the segment machinery to
per-segment root-finding: the letter is composed exactly along each
straight chart segment (Gaussian-rational arithmetic — degree drops from
special segment directions are detected exactly, never by threshold), its
t-polynomial is rooted with `mp.polyroots`, and since
w(t) = lc·∏(t−r_k) with lc constant in t,
dlog w(t) = Σ_k dt/(t−r_k): each root contributes a unit pole term to the
same inline DE system, clearance is certified as the nearest root's
distance to the real segment [0,1] in t-units (the rational-block
convention, same `clearance_min` bar), and the continuous winding is the
exact per-factor sum of the affine tracker applied to each linear factor.
Fail-closed: root-finding failure refuses the segment (clearance cannot be
certified without the roots), and polynomial-letter marches refuse
floating waypoints — pass exact rationals (strings, ints, `Fraction`s), as
`evaluate` already does when given rational-string waypoints.

### Scope (further limits)
- `connection` entries are eps-independent exact rationals. eps-dependent
  systems stored on disk go through wayfinder's loaders on their own
  variable.
- The kinematic derivative map is 2×2 (two chart variables, two invariants).
- Elliptic families use the same config slots — the anchor/values plugin
  carries the period machinery, so the harness contract does not change. This
  is the point of the config design; no elliptic family ships in this
  release.
- First derivatives come exactly from A·J; higher jets through discontinuities
  are out of scope.

## Rational connection blocks

Beyond the dlog form, famhar can march an exact *rational* connection declared
in a config's `connection_rational` slot: one sha256-pinned file per operator,
entries `"i,j" -> [[num_terms],[den_terms]]` with each term
`[[exponents], coeff]` over the file's `vars`.

- `Family.load_connection_rational(spec=None, verify_sha=True)` →
  `RationalConnection`. File paths resolve relative to the config's directory
  (or its parent if the config lives in a `families/` directory); every file's
  sha256 is verified against its pin, fail-closed.
- `Family.march_rational(path_spec, y0, dps, eps0=None, subs=None, mtay=None,
  ratconn=None)` → `(J_end, record)`. Straight chart segments are pulled back
  exactly: waypoints must be exact rationals (`parse_point_exact`; floats
  refuse), entries are composed into univariate rational functions of t with
  exact Gaussian-rational coefficients, and each segment serves wayfinder's
  DE-system contract (`.A`, exact-shift `.A_series`, `.singular_points` = the
  denominator roots). `subs` pins every non-chart file variable exactly, e.g.
  `{"d": "28/5"}` for a fixed-eps march with d = 4−2eps; the caller owns the
  eps→variable convention. The sheet record tracks continuous winding and
  continued log per distinct denominator, with per-segment clearances,
  diagnostics, walls, the closure census, and file provenance.
- `Family.A_rational_at(point, dps, subs=None, ratconn=None)` — pointwise
  sparse evaluation `{op: {(i,j): mpc}}`. This is legal even on a block that
  is not closed; only marching requires closure.

Refusals (fail-closed, never worked around):

- **`SYSTEM_NOT_CLOSED`** — `march_rational` refuses unless every entry column
  is covered by the block's own DE rows: marching an open block would
  fabricate a solution. Because the sparse format cannot distinguish a zero
  row from a missing row, coverage credit for zero rows requires an explicit
  `zero_rows` declaration in the files; a declared-zero row that carries
  entries refuses at load.
- Any denominator root in t within `clearance_min` of the segment refuses,
  naming the detour; if root-finding fails, the segment refuses too, since
  clearance cannot be certified without the roots.
- Non-exact waypoints, unpinned non-chart variables, sha mismatches,
  empty or zero denominators, and negative exponents all refuse.

Measured on the calibration family below: a rational-block march at dps 30
took 9.77 s/step versus 7.32 s/step for the equivalent dlog march (1.33× the
cost, same 22 steps), with 47-digit agreement between the two engines.

## Validation: the four-loop ladder family

The harness was calibrated on the Usyukina–Davydychev ladder integrals
Phi^(1..4) (conventions of arXiv:1303.6909), a family whose closed forms make
every check free. Chart: a = z/(z−1), b = z̄/(z̄−1) treated as independent,
with u = ab/((a−1)(b−1)), v = 1/((a−1)(b−1)) and prefactor
−1/(z−z̄) = (a−1)(b−1)/(a−b). Branch loci in the chart are the letter zeros
{a, 1−a, b, 1−b}; a = b is a prefactor pole only, not a DE singularity.
Masters: {ℓ^r Li_s(a), ℓ^r Li_s(b), ℓ^r} with ℓ = log(ab), r+s ≤ 8, n = 65.
The connection is derived analytically (dlog form, integer entries) and
validated against finite differences of the closed-form basis before any
transport is trusted.

Quoted digits use a two-precision rule: every number is marched at dps 50
*and* 80, and the quoted count is the minimum of agreement-with-reference and
pair agreement. The full battery (acceptance bar 30 digits; 25 for G0):

| gate | quoted digits | wall |
|---|---|---|
| G0 — connection A·J vs finite differences (260 checks, 2 fresh points, both vars) | 38 worst | 0.2 s |
| G1 — values at two points, Phi^(1..4), vs an independent 100-digit closed-form reference | 50–51 (all 8) | 766 s |
| G2 — dPhi4/du, dv vs independent termwise closed-form derivatives | 50 (all 4) | 0.2 s |
| G3 — identity loop (windings 0) | 50/80 componentwise | 17 s |
| G4 — branch-locus loop vs exact independent monodromy; winding exactly 2π | 50/80; disc 49/80 | 776 s |
| G5 — second-sheet point through the a ∈ (1,∞) cut vs a continuation-safe parametric reference, componentwise on the continued masters and on Phi4 | 50/81; naive-principal asserted wrong by exactly the disc, 47/78; crossing-sign pin 40 | 623 s |
| refusal regressions — through-locus and near-locus paths refuse; tracker vs independent dense winding | 46 vs dense sampler | 0.0 s |

All gates pass. The 100-digit reference values were independently
cross-checked against the one-fold parametric integral representation, and
the G5 crossing sign was pinned by an independent elementary continuation.

The ladder gate battery itself (`gates_ladder.py` and its receipts,
`GATES_LADDER.json`) consumes the 100-digit reference values and writes
receipts beside them; it is a record of the calibration run and is not
included in this repository (`evaluator_harness/gates_ladder.py`, sha256
`752b454378b66386…`, 434 lines). What ships here is the engine
(`famhar.py`), the config schema, the calibration family itself — the
closed-form reference plugin `ladder_reference.py`, its config emitter
`build_ladder_config.py` and the emitted config `families/ladder_L4.json`
(`PROVENANCE.json` pins each vendored file to its record) — and a
self-contained selftest battery:

```
python3 tools/famhar/famhar.py --selftest
```

runs in seconds on synthetic families whose closed forms make every check
free — a straight march gated on values and A·J chart/kinematic derivatives
against closed forms, an identity loop (zero windings, J returns), a
branch-locus loop gated on winding = 2π and on the exact continuation of
the closed-form basis (a function-level reference, never exp(2πi·Res A)),
a polynomial-letter leg (value, derivative, locus-loop monodromy through
the per-segment root-finding path), the refusal laws (near-locus,
through-locus, endpoint-on-locus, floating waypoints on a polynomial
family), and the ladder sheet fixture below (skipped by name — `SKIP:
absent <file>` under its `[B6]` banner — if the plugin or its fixture is
absent). Each leg registers itself as it opens, and the verdict line's
counts are read from that registry (`legs 6 executed / 0 skipped` on the
shipped tree), never typed. The verdict is `PASS_ALL` (exit 0) only when
every leg executed and every check passed; `FAIL` (exit 1) when any check
failed, the failing checks printed by name; `PASS_EXECUTED_WITH_SKIPS`
(exit 2) when no check failed but a leg was skipped — a skipped leg is
never an all-pass, and a consumer reading the token or the exit code alone
sees it. The battery has no expected-fail legs. To gate
your own family you supply its JSON config; the schema and refusal
behavior above are the contract, and the ladder battery is the pattern to
copy for gating a new family against whatever independent reference it
has. The pytest battery `tests/test_ladder_sheet.py` carries the full
sheet fixture, the planted controls, the builder drift guard and a
one-point A·J-vs-finite-difference check of the shipped config.

### The sheet of ell in the reference plugin

The plugin's master ell is `log(a·b)` on the anchored principal branch
(real at the anchor 0 < a0, b0 < 1). `phi_value(L, a, b, dps, cont=None)`
realizes it as `log(a) + log(b) + 2πi·k` with k chosen so the sum lands on
the principal branch of the product (`default_ell_k`): k = 0 on the
Euclidean slice b = conj(a) and at the anchor, k = −1 on the real λ > 0
sheet (z, z̄ real, u + v < 1), where the chart point a = w, b = w̄ has both
coordinates negative real and the naive principal sum overshoots the real
master by one turn. An explicit `cont={'ell_k': k}` keeps its meaning
relative to `log(a) + log(b)`; `ell_sheet(a, b, dps, cont)` reports the k
in force and whether it was selected or given. `phi_derivs_chart` uses the
same default. The fixture `tests/fixtures/ladder_sheet_receipt.json` pins
Φ^(3) at (u, v) = (1/4, 1/10) on that sheet from an identity receipt
(record sha256 `3c50aeda13a5c9a2…`, computed at dps 40): the default and
the explicit `ell_k = −1` call agree with it to the record's precision
(measured 40.3 digits at dps 50, the record's own cap), agree with each
other and with a dps-80 run to ≥ 48 digits, and `ell_k = 0` reproduces the
old wrong branch (real part −383.41…, |Im| 113.78) as a control that fails
by name if the default ever returned it.

## Performance

- Ladder family, n = 65, dps 40: ~3.0 s/step (the transport recursion is
  sparse-aware; this connection has 156 nonzero entries). A single-segment
  23-step march cost 69.7 s.
- Cost grows roughly as dps^1.7 (wayfinder's measured scaling). Time a
  short pilot march before committing to anything at dps ≥ 100.

## Files

- `famhar.py` — the engine. Generic: it contains no family content.
- `ladder_reference.py` — the closed-form reference plugin of the ladder
  calibration family (continuation-safe; the sheet-selected default for ell).
- `build_ladder_config.py` — emits `families/ladder_L4.json` (exact
  Fraction entries) from the plugin's basis; no reference data needed.
- `families/ladder_L4.json` — the emitted calibration family config.
- `PROVENANCE.json` — record sha256, vendored sha256 and the listed changes
  of every vendored file, plus the fixture pin.
- `tests/test_ladder_sheet.py`, `tests/fixtures/ladder_sheet_receipt.json`
  — the ladder sheet battery and its record fixture.
- `GUIDE.md` — condensed reference card (purpose, invocation, footguns).

Not shipped (records of the calibration run, named for the reader):
`evaluator_harness/gates_ladder.py` (sha256 `752b454378b66386…`) with
`GATES_LADDER.json` and the 100-digit reference values `phi_values.json`.

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).

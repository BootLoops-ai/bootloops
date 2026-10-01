# CONI-PFV EXTENSION — END-TO-END VALIDATION

Scope: validation gates v1-v4 of the conifold-PFV extension (DESIGN.md), run on
the conifold-PFV example cards of arXiv:2406.13751. None of those cards is
included in the package; the example quoted below, the 'Lórien' candidate of
that paper (card name ds-lorien), and the non-simplicial 7-ray AdS example
ads-5-113-4627-main record what the code reported on cards that do not ship.
Compute: single core, ulimit -v 31G, nice 5, minutes-class.

## Verdict table

| gate | statement | verdict |
|---|---|---|
| v1 | dS certified ball at dps 60, verdict-grade width | **FAIL — STRUCTURAL (designed fail-closure, not a code defect)** |
| v2 | value vs published/pinned values | **PARTIAL — estimate grade only** (no certified margin possible without v1) |
| v3 | log-branch mutation caught | **PASS** (CP9 cross-gate, both example cards) — one documented escape + one classification limitation, below |
| v4 | non-simplicial cone fix g1 PASS + DKMM regression byte-exact | **PASS** (4627-main battery ALL PASS; dkmm g1 outputs sha256-identical) |

Consequence: the extension validates as designed at g2'/g3' (estimate grade),
but no dS card reaches a certified W0 through the transport legs (see v1);
that requires the coni_op derivation route (coniop/).

## v1 — certified transport: structural diagnosis

The target chain is curve -> series -> certified transport (dps 60) -> W0 on
one dS vacuum. It is not exercisable on a card that lacks the transport
inputs, by design:

* No shipped card carries `coni_op` / `coni_towers` / `coni_routes` / `s_star` /
  `coni_jet_ideal` / `coni_C2`.
* `coni_transport.transport_ready` = False on every conifold-PFV card tested;
  `run_leg1` raises NotImplementedError until a coni_op-bearing card exists —
  DESIGN.md §7 scope.
* Therefore none of the DKMM transport-class gates (two-dps, dual-route,
  ball honesty, Krawczyk, FD-Jacobian) has a coni-PFV analog that can run on
  these cards. Their analogs are specified in DESIGN.md §3 (legs 1-3 reuse the
  parametric pipe_transport/pipe_vac majorants unchanged) and are unvalidated
  until a card carries coni_op. No transport is faked.

What the coni_op derivation route (coniop/) must deliver per card before v1 can
close: the bulk theta-form operator L_s (z_cf->0 elimination of the card's
UNGATED GKZ boxes, or annihilation of the class-graded restricted series — the
conifold-ray resummation (D5) makes the graded series a research step, not a
mechanical one), exact MUM log-towers, envelope-scan routes, exact s_star, the
2-parameter jet ideal (LEG 2), and a certified C2 bound (LEG 3).

## v2 — estimates vs published/pinned values (ESTIMATE grade, labeled)

End-to-end g2'+g3' runs (byte-identical on repetition):

| card | quantity | here (PFV step-1 est.) | published/pinned | margin |
|---|---|---|---|---|
| ds-lorien | Im tau | 18.1822 | 17.5173 (tau_pin) | +3.8% (gate ok, bar 20%) |
| ds-lorien | lam | 93/19 exact | published Table-3 K' = 93/19 | EXACT |

These are window-limited PFV step-1 estimates gating CONSISTENCY only (the
2406 step-1/step-2 gap dominates); the certified statement remains open with
v1.

## v3 — mutation control (log-branch convention)

* **Mutation B (the flip proper):** log(-2 pi i z_cf) -> log(+2 pi i z_cf),
  i.e. the i*pi branch shift absorbed into the F-term magnitude equation
  (joint_vev zmag line). CAUGHT on both example cards by the CP9 cross-gate
  "F-term vev vs Q_throat vev disagree beyond Li-correction level"
  (gap = pi = 3.14 e-folds > 1.5 bar; the Q_throat form is an independent
  pinned formula, so the gate is a true cross-check, not a tautology).
  Nonzero exit; run_conipfv reports the card CONI-FAIL.
* **ESCAPE (documented):** mutation A, the phase-only variant (z* = -i zmag
  entering the tau F-term with flipped sign, block coefficients untouched):
  NOT caught — t moves 0.12%, W0_est moves ~1e-5 rel;
  all estimate-grade bars (20% / 0.75 / 1.5) are orders of magnitude wider.
  This is inherent to estimate-grade gates: an O(zmag) ~ 1e-5 phase defect is
  invisible. The certified transport legs (v1) are the layer that must
  catch phase errors; this is an acceptance test for that leg.
* **Classification limitation of run_conipfv.run_card:** under mutation B
  ds-lorien is reported CONI-BLOCKED (exit contribution 0), not CONI-FAIL,
  because run_card classifies ANY AssertionError as "blocked" when the card
  carries a `blocker` field (ds-lorien does). A genuine regression on a
  blocker-bearing card is therefore not visible in the sweep exit code; read
  the gate_fail text of such cards.
* Restoration check: with the pristine coni_frobenius.py a clean re-run is
  byte-identical to the reference gates JSONs.

## v4 — non-simplicial cone fix (routed design)

* ads-5-113-4627-main battery (a 7-ray non-simplicial h=5 card, fresh code
  path): [route] eff_route non-simplicial, 7 extreme rays / 8 facets exact, all
  pinned + 88 gv_table classes inside the cone; [failclose] gv_extract raises
  the structural message; [frame/gv-table] sigma=(-1,1,1,1,-1) eps=+1
  verified, 88 integer classes, 7 leader pins exact; [peel] MN=44 (exercises
  the nu_a>MN rule, nu_0=184): pure orders 14/29/44 peel == pins AND naive
  m-graded sums == ALL 5 pins; composite orders 28/43 peel-contaminated
  (m=28: 7172 vs -4) — peel is a diagnostic, not a gate, on non-simplicial
  cones; [w0-prefix] support(w0, s<=100) == semigroup<14,72>.
* DKMM regression through the patched geff_series: run_pipe.py g1 ALL PASS
  (48 GV, lit extras, 901 w0 coeffs, tower to s^240, N_3=3/N_4=540, mirror
  map; wall ~5 min) and the four g1 outputs (GV table, series, tower, mirror
  map) sha256 BYTE-IDENTICAL before/after the fix.
* Reproducibility control: repeated run_conipfv single-card runs (dkmm +
  4627-main NOT-CONI controls; ds-lorien and two other dS examples
  CONI-PENDING) byte-identical.
* **Scoreboard note:** a single-card invocation of run_conipfv.py rewrites the
  scoreboard JSON in out/ with only the cards given on the command line; the
  full board is the one written by `--all`.

## Bottom line

g2'/g3' + routing + controls + the cone fix validate cleanly and
reproducibly; the mutation control passes with one documented estimate-grade
escape and the classification limitation above. The certified end-to-end leg
(v1/v2) is structurally out of reach for every dS card until a card carries a
derived coni_op (coniop/ route).

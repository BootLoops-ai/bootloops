# Closure interface — READ-ONLY evaluation tier [bench]

A boundary surface for plug-in log-likelihood closures: a consumer writes a
closure to the harness spec and gets register-labeled evidence evaluations
back. Consumers make no direct calls into ungated internals — **this
interface is the boundary**.

READ-ONLY EVALUATION TIER; INTERNALS REMAIN GATED. The module imports the
pinned R3 upstream build privately, sha-pins all 14 upstream files at load
(drift fails loud, G1), returns only plain JSON-able dicts, and re-exports
nothing from behind the boundary.

## Files

- `closure_interface.py` — the surface (the ONLY sanctioned entry point)
- `closure_fixtures/fixtures.json` + `closure_fixtures/FIXTURES_README.md`
  — bench fixtures: exact closed-form truths, recorded measured tolerances,
  plants
- `selftest_closure_interface.py` — consumer-runnable, no input data needed
  (numpy+scipy, ~1 CPU-min); receipt below
- `CLOSURE_INTERFACE_MANIFEST.sha256` — shas of everything staged + the
  pinned upstream

Sanctioned build behind the boundary: the pinned R3 build (the
`certified_evidence` package; env `POSQ_R3_ROOT` points at it). That build
is not distributed with this repository, so the interface is reference
material here: the contract, the fixtures and the selftest ship, but no
number is produced without the pinned upstream. The pins were taken against
a build already exercised end-to-end (`PlugInLoglik` wired with no code
changes).

## Closure objects — object-by-object coverage

Coverage, also available mechanically as `closure_interface.coverage()`:

| object | status | where |
|---|---|---|
| loglik-closure harness spec | **STAGED** | `HARNESS_SPEC` + `check_closure()` — the contract a plug-in closure is written to; covariates/offsets/heterogeneity live inside the closure |
| MC-register evaluation | **STAGED** | `evaluate_evidence()` -> logZ at register `R3-certified-calibrated`, certificate attached; evidence VALUES only |
| register law / labeling | **STAGED** | `MC_REGISTER`, `REGISTER_LAW`; every result self-labels and carries the calibration-transfer caveat |
| hazard-field loglik closure | **MISSING-BY-DESIGN** | the closure is the consumer's own contribution; fixture F3 demonstrates a hazard-shaped closure wiring through unchanged |
| posterior over {h_b} (rate re-inference posteriors) | **MISSING** | no sanctioned posterior-sampling surface exists behind the boundary — the validated packages ship evidence values, not draws; the pilot MH is ungated scaffolding and stays behind the boundary. `infer_hazard_field()` / `posterior_draws()` raise `ClosureObjectMissing` with reason + unlock |
| posterior-draw MC runner | **MISSING** | no sanctioned posterior-draw runner exists behind the boundary |

Consequence stated plainly: texture statistics that need POSTERIOR DRAWS
(e.g. burstiness as the CV of a rate field {h_b}, or credible intervals on
likelihood-ratio statistics) cannot be served by this tier. What IS live: a
consumer can write a log-likelihood closure to `HARNESS_SPEC`, validate it
with `check_closure`, and get register-labeled evidence evaluations (model
comparison, closure debugging, semi-synthetic calibration of the closure
itself) at bench cost. If a calibrated, gated draw surface is ever built,
this interface grows the call and consumers re-wire nothing.

## Usage (bench)

```python
import sys; sys.path.insert(0, '<your-checkout>/tools/posq')
import closure_interface as ci

print(ci.HARNESS_SPEC)              # the closure contract
ci.check_closure(my_loglik, dim=9)  # spec compliance, read-only
res = ci.evaluate_evidence(my_loglik, dim=9, budget=100000, seed=1,
                           name='hazard-A')
res['logZ'], res['register'], res['calibration_transfer']
ci.coverage()                        # the table above, mechanically
```

Binding caveat on every number: register is `R3-certified-calibrated` =
float-class; the shipped calibration was measured on JC anchors and does
NOT transfer to plug-in closures — re-grade against exact/enclosure anchors
of the consumer's own object before leaning on error bars. R3 output never
feeds a certified claim.

## Selftest

`python3 selftest_closure_interface.py` -> **PASS 14/14 gates**: G1 pin
(14 upstream shas), G2 spec x3 fixtures + plants P1 (shape) and P2 (NaN)
fired, G3 evidence vs exact truths (worst 1.6e-5 / 6.8e-6 / 1.6e-4 nats vs
tols 5e-4 / 5e-4 / 5e-3; 4 seeds x budget 1e5 each), G4 register labeling
on every result, G5 missing-object law (both raisers fire with
reason+unlock). Exit 0. (Needs the pinned upstream build at `POSQ_R3_ROOT`;
without it, G1 fails loud and no number is produced.)

# closure_interface bench fixtures

Three small, self-contained loglik closures with analytically exact truth
values, plus three plants that must fire. They exist so any consumer can
verify the closure
interface end-to-end WITHOUT any phylogenetic context: no alignment, no
tree, no stored dataset — just `python3 ../selftest_closure_interface.py`.

Truth values are EXACT closed forms (that is the exact register; they were
computed at 40-dps mpmath precision, quoted to 30 digits). The estimator
outputs gated against them are FLOAT-CLASS R3 numbers; the gate tolerances
are recorded in `fixtures.json` and were derived from a measured 8-seed sweep
(2026-07-20 UTC, numpy 2.4.1 / scipy 1.17.1), with >= 20x margin over the
worst measured deviation at HALF the gate budget. Measured numbers ride in
`fixtures.json` under `measured`.

## The fixtures

- **F1-separable-power** (d=3): `loglik(u) = sum_i a_i ln u_i`, a=(1,2,3).
  Z = prod 1/(a_i+1) = **1/24 exact rational**; lnZ = -ln 24.
- **F2-constant** (d=2): `loglik(u) = -7/2`. Z = e^{-7/2} exactly; the
  degenerate control (constant IS weights).
- **F3-hazard7-toy** (d=7, hazard-field-SHAPED): per bin b in 1..7, event count
  d_b with exposure E_b and hazard h_b = u_b (uniform prior):
  `loglik(h) = sum_b [ d_b ln h_b - h_b E_b ]`, d = (3,2,9,2,3,8,2),
  E_b = 40. Exact truth per bin:
  `Integral_0^1 h^d e^{-E h} dh = E^{-(d+1)} * gamma_lower(d+1, E)`
  (substitute t = E h), so lnZ = sum_b [ ln gamma_lower(d_b+1, E_b)
  - (d_b+1) ln E_b ]. Re-derive with mpmath:
  `sum(mp.log(mp.gammainc(d+1, 0, E)/mp.mpf(E)**(d+1)) for d,E in ...)`
  at `mp.mp.dps = 40`.

  F3 is a hazard-field-SHAPED demonstration ONLY: 7 time bins with
  elevated counts in bins 3 and 6 (a "crisis" texture). It is
  NOT a full hazard-field model — no per-slot covariates, no per-bin
  exposure offset, no gamma slot heterogeneity. Those belong inside
  the consumer's own closure;
  F3 shows that a hazard-field-style closure wires through HARNESS_SPEC
  unchanged.

## The plants (must fire)

- **P1** wrong output shape `(n,2)` -> `check_closure` raises.
- **P2** NaN in output -> `check_closure` raises.
- **P3** `infer_hazard_field()` / `posterior_draws()` -> raise
  `ClosureObjectMissing` with non-empty reason AND unlock (missing
  objects are refused mechanically, never improvised).

## Run

```bash
python3 <your-checkout>/tools/posq/selftest_closure_interface.py
```

Exit 0 = all gates pass. Requires numpy + scipy (the pinned upstream uses
both); no input data, no network, ~1 CPU-min. The selftest also re-verifies
the upstream sha pins (gate G1) — a consumer on a drifted build fails loud
before any number is produced.

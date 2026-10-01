# periods/pipeline/pfaffian — scripts behind `periods/pfaffian.py`

`DESIGN.md` is the method note (why a first-order matrix system, gate
tables, cost model). `ads581/` holds the route-A scripts for the card
ads-5-81-3213. The route-B scripts (the conifold-curve example) are not
included in this release. The library front door and the
replay battery are one level up: `periods/pfaffian.py`, `periods/PFAFFIAN.md`,
`periods/selftest_pfaffian.py`.

The scripts are the stage-by-stage drivers that produced the reference
results the battery replays. Each one is run from its own directory, reads
the outputs of the previous stage from that directory, and writes its own
output next to itself. Python 3 with numpy and python-flint.

## Reference data (not included in the package)

Two kinds of input do not ship:

* The ads-5-81-3213 D-module bank: the generic box-operator module
  (`boxops`: reduce_sides, poly_from_factors, theta_mul), the big-int-safe
  mod-p polynomial module (`pmodp`: prodtree, ratrecon2, pgcd, pdivmod,
  pmulc), the 17x17 connection samples per prime, the exact fundamental-period
  series to s^3000 and the flux frame (F, H) of the card. Set
  `TERRIER_REFERENCE_BANKS` to a root laid out as the scripts open it:

      $TERRIER_REFERENCE_BANKS/
        cards/
          ads-5-81-3213_bank/
            gate_g2.json              flux frame: {"F": [...], "H": [...]} (12 entries each)
            w0_series_3000.json.gz    {"a": [a_0, ..., a_3000]} exact fundamental-period coefficients
            dmodule/
              boxops.py  pmodp.py     the two modules above
              conn_samples_<p>.npz    keys sig (nodes), A (17,17,nodes) mod p

  Scripts that need it refuse with a one-line message when the variable is
  unset: `ads581/a1_perprime.py`, `gate_pa1.py`, `gate_pa2.py`, `verdict.py`.
* `ads581/verdict.py` additionally reads a GV-invariant table for the card
  placed next to the script (`gv_full.json` = {"d1,d2,d3,d4,d5": n_d}); the
  card's pinned GV values are asserted against it.

Everything else (cards, `family.py`, `geff_series.py`, `pipe_lib.py`,
`conipfv/coniop/gseries.py`) is in the package.

## Run-time files

Outputs have fixed names and are written next to the script that makes them;
the next stage reads them from there. Per-prime and per-order files carry the
parameter in the name (`interp_<P>.npz`, `jets_<M>.json`, `tower_<M>.json`).
The fixed names are:
`A_exact.json` (exact connection: den, entries, deg table, primes used, max height in bits, stability flag) and `deg_table.json` from the lift; `gate_PA1.json`, `gate_PA2.json` from the gates; `f3_frame.json` (frame-fit coefficients) and `verdict_raw.json` (certified W0 balls at dps 60 and 90) from the route-A verdict.
Below, "connection" means that exact connection file and "bank" means the
reference-data root above.

## Route A — `ads581/` (ads-5-81-3213, rank 17)

| step | command | reads | writes |
|---|---|---|---|
| per-prime interpolation | `a1_perprime.py P [NVAL]` | bank samples at P, bank `pmodp` | `interp_<P>.npz` (Anum, DEN, degs mod P) |
| exact lift | `a1_lift.py P1 P2 [P3 ...]` | `interp_<P>.npz` | connection + degree table |
| gate PA-1 | `gate_pa1.py P_fresh` | connection, bank samples at P_fresh | PA-1 result |
| held-out jets | `jets_series.py M` | — | `jets_<M>.json` (17 analytic-solution jets to s^M; PA-2 uses M=3000) |
| gate PA-2 | `gate_pa2.py` | connection, `jets_<M>.json` at M=3000, bank s^3000 series | PA-2 result |
| verdict | `verdict.py` | card, bank flux frame + s^3000 series, GV table; caches `tower_<M>.json` | frame-fit coefficients, certified W0 balls |

## Route B (conifold-curve cards)

The route-B stage scripts (box-operator supply check, closure certificates,
connection sampler, per-prime interpolation, exact lift, gates PB-2 / PB-3,
exact series oracle, entropy certificate with certified landing, and the
ball marcher) are not included in this release; `DESIGN.md` section 3
describes the route and the library functions in `periods/pfaffian.py`
(`entropy_landing`) implement its landing step.

`periods/selftest_pfaffian.py` replays the route-A fresh-prime gate, a
mutation control and the certified route-A verdict against stored copies
of these outputs (`TERRIER_PFAFFIAN_BANK`; see `periods/PFAFFIAN.md`).

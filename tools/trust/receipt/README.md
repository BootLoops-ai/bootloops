# RECEIPT — every reduction table ships with a receipt

The receipt member of the trust package (`tools/trust/receipt`): a
self-contained directory that also runs standalone from any location.

**Positioning:** a solver-agnostic standalone verifier/emitter CLI — one witness contract for any producer.

Per-(kinematic point, prime) lambda-multiplier certificates for IBP/Laporta
reduction tables: **verify** any claimed table row against the generating
system in ~1 ms/row, **emit** certificates for any existing table (retrofit
Wiedemann driver), **detect** false-master/structural defects.

## The point: solvers become drop-in replaceable

The certificate layer makes ANY future reduction algorithm a drop-in
replacement. Correctness is verified **per row, independent of the
solver**: a witness is an exact linear identity

```
sum_i lam[i] * R_i  =  e_t - sum_m c[m] * e_m     (mod p)
```

over the staged system's rows, checkable by an independent parser in a
millisecond. Trust attaches to the receipt, not to the engine — so an
experimental engine can enter production without a trust wall: its rows
either carry receipts or they don't. The verify path is strictly
solver-agnostic by design invariant (WITNESS_FORMAT.md): `core.py` imports
no solver code and makes zero assumptions about the table's
producer. Witness format = the contract.

## The exhibit: a confident wrong table (why this layer is non-optional)

The wrong-table exhibit — an
experimental cohomology-pairing engine produced reduction rows for two
sectors of a two-loop family that were **internally beautiful**: D-ladder
stable, rank-clean, and 2-prime CRT-consistent, reconstructing to structured
rationals at the node (`(d-4)/(eta+5)` and `-1/(2(eta+5))`-class values).
Every internal consistency check a reduction pipeline normally runs —
including the "compute at two primes and compare" folklore gate — PASSED on
a wrong answer. Only the certificate layer caught it: the archived oracle rows carry
lambda witnesses; the pairing's rows are provably NOT in the system's row
span with those coefficients (0/16 certified, two structural root causes,
measured). Without receipts this ships as a confident wrong table
(the same failure class appeared in a second independent pipeline: 25 wrong rows with a complete master basis, rejected only by two independent oracles).

That exhibit is this repo's canonical MUST-FAIL regression
(`tests/test_banked.py` leg 3): the 6 wrong (row, prime) pairs must fail
verification forever; the archived oracle rows must pass.

## Quickstart

```bash
T=<path to tools/trust/receipt>
RL=<your run root>
ART=$RL/gen_myfam                                      # kira staging artifacts
STG=$RL/staging/myfam                                  # staged config

# 1. EMIT lambda-witnesses for an existing kira table (retrofit; both primes)
python3 $T/receipt.py emit --adapter strata \
  --system-dir $ART --staging $STG --family myfam \
  --table $STG/corpus/<id>/kira_target.m \
  --slices 2147483647:1234577:87654321,2147483629:1234577:87654321 \
  --out /path/to/witnesses

# 2. VERIFY witnesses (independent parse; ~1 ms/row after system load)
python3 $T/receipt.py verify --adapter strata \
  --system-dir $ART --staging $STG --family myfam \
  --witness '/path/to/witnesses/w_*.json'
# rc: 0 all pass | 1 any row fails | 2 malformed input

# 3. VERIFY a claimed table against its receipts (wrong-table catch mode)
python3 $T/receipt.py verify --adapter strata ... \
  --witness '...' --table claimed_rows.json     # rc=1 if the table lies

# 4. DETECT structural defects (false masters, master relations,
#    uncovered targets)
python3 $T/receipt.py detect --adapter strata \
  --system-dir $ART --staging $STG --family myfam \
  --slices 2147483647:1234577:87654321
```

Generic (non-kira) systems: `--system-jsonl` (one JSON row-object per line)
replaces the adapter everywhere; see WITNESS_FORMAT.md for the contract.

## What's measured (this README makes no new claims)

| number | value | source |
|---|---|---|
| verify cost | 0.69 ms/row CANONICAL (median of 3 on a pinned quiet core; 0.57–0.96 on a loaded machine; this repo's regression re-measures and prints it) | archived timing receipt |
| emit cost | ~2.2–2.3 s/(row,prime) at the 7,135-eq reference family; min-poly reuse floor 0.86–1.0n matvecs | archived probe receipts |
| retrofit record | 84/84 rows (42 registry rows x 2 primes) oracle-match + independent-checker pass | archived retrofit receipt (replayed by tests/test_banked.py) |
| witness size | KB-class: density mean 3.3%, <= 32 KB/row/prime flat | archived probe receipts |
| sabotage battery | 6 classes caught at both primes (S1 perturb, S2 label swap, S3 scale, S4/S6 lambda tamper + engine-level wrong-twist/wrong-basis) | archived sabotage receipts |
| wrong-table exhibit | 0/16 pairing rows certified; CRT-consistent wrong values caught only by receipts | archived exhibit verdict (replayed by tests/test_banked.py) |

## Scope honesty

- Certificates are **per-(kinematic point, prime)**. No cross-point
  (Schwartz-Zippel) claim; symbolic tables out of scope; 2-loop
  measurements do not derisk 3-loop claims (pre-registered scope).
- A witness certifies row-span membership in the GIVEN system. System
  completeness/correctness is separate: `detect` (structural audit) +
  staging provenance gates. Measured per-row solver economics run
  >= 4,500x the whole-table eliminator's — this tool ships the
  certificate layer, not a solver replacement.
- Emission uses a column partition; the strata adapter stages it via one
  eliminator pass at the anchor slice (seconds-class at the 7,135-eq
  reference scale).

## Layout

```
receipt.py            CLI (verify / emit / detect); TOOL_NAME constant
core.py               solver-agnostic verify core + witness IO (the invariant)
emitter.py            Wiedemann retrofit driver (probe-harness port)
detector.py           false-master/structural audit + self-contained eliminator
adapters/strata.py    ALL kira/STRATA specifics (loader bridge, labels, legacy)
WITNESS_FORMAT.md     the versioned contract (v1.0)
tests/                test_core.py (synthetic + 6-sabotage + rc discipline)
                      test_banked.py (30 archived pairs; 84/84 retrofit;
                                      wrong-table MUST-FAIL exhibit; throughput)
                      mutation_test.py (mutates the CHECKER on scratch
                                        copies; suite must catch every mutant)
                      run_tests.sh
```

## Integration (any live run)

A solve+oracle pipeline can consume receipts directly:

- **emit-side**: after the oracle-certification pass, call
  `adapters.strata` + `emitter` per certified row (or the `emit` CLI on the
  run's staged dir + table) to store witnesses next to the run ledger;
- **verify-side**: any later consumer replays `receipt.py verify` against
  the staged system — no kira relaunch, no solver trust.

This tool does not write into live run directories.

Dependencies: numpy + the strata toolchain for the ADAPTER only — found by
default at the sibling strata member of the trust package (`tools/trust/strata`,
i.e. `../strata`);
`RECEIPT_STRATA_PATH` points it at an external strata tree instead.
The core has no dependency beyond the standard library (hashlib/json) — by
design, so a referee can reimplement the checker from WITNESS_FORMAT.md
alone.

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).

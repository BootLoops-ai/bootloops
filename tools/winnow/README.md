# Winnow (`ibplapper`) — upgraded-Laporta F_p eliminator with native receipts

Winnow is the exact elimination engine behind Seedling: a standalone
eliminator LIBRARY callable by Kira or any IBP pipeline. One sentence:
**sparse identity rows over F_p plus an ordering/stratum policy in; closed
reduced rows, per-row lambda-witnesses, and a signed audit ledger out —
everything upstream of the rows and downstream of the per-(point,prime)
tables stays caller-side.** The Python package/import name is always
`ibplapper`, never bare `lapper` (rust-lapper/nim-lapper own that name in
genomics); `import winnow` is a verbatim re-export alias.

Certificates are **per-(kinematic point, prime)**. No cross-point claim, no
symbolic-table claim. 2-loop measured; anything 3-loop is labeled
extrapolation.

## Quickstart — library

```python
import ibplapper as lap

sys = lap.System(
    rows,             # iterable of sparse rows: dict[col:int -> val mod p]
    p,                # ONE prime per System (per-(point,prime) scope)
    order,            # dict col -> rank; HIGHER rank eliminated FIRST
    forbid={...},     # master columns: never pivoted on
    stratum_of=None,  # optional dict col -> stratum (higher = earlier);
                      # default single stratum
)
sched = lap.Schedule(
    policy="B2FT",              # B2 | B3 | B2F | B2FT (answer-identical;
                                # policy is a COST knob only)
    bank=lap.Bank(path),        # optional interface table T (append-only,
                                # sha256-verified read-back)
    caps={"stratum_wall_s": 60},  # loud CensoredError, never silent
)
res = lap.eliminate(sys, sched, witnesses="retrofit",
                    witness_targets=[cols...])

res.subs        # pivot -> CLOSED row (master/forbidden cols only)
res.leftover    # relations among masters — NEVER silently absorbed
res.witnesses   # per target: lambda vector, c, verify status, file path
res.ledger      # audit ledger; res.save_ledger(path) appends JSONL

from ibplapper import receipt
ok, detail = receipt.verify(my_own_rows_parse, "w_2147483647_....json")
```

Column ids are opaque ints. `order`/`stratum_of`/`forbid`/`meta` are the
entire physics surface.

## Quickstart — Kira sidecar (Route A)

Run Kira with `run_initiate: true` only (Generate+Select, no back-sub),
then:

```python
from ibplapper.adapters import kira
system, sysd = kira.load(art_dir, family, p, d0, eta0,
                         expect_counters={"eqs": 7135, "terms": 42286})
res = lap.eliminate(system, lap.Schedule(policy="B2FT"))
```

## Quickstart — user-defined system (Route B)

Route B reads the DOCUMENTED, stable kira `user_defined_system` input
format (the equation files a `reduce_user_defined_system` job consumes:
one `integral*coefficient` term per line, blank-line equation splits,
integer-weight or `name[indices]` notation, plain or gzipped, file or
directory of `.kira`/`.kira.gz`). No Kira installation is needed — these
are files a user writes or any IBP generator emits:

```python
from ibplapper.adapters import user_system
system, sysd = user_system.load(
    "eqs/", p, values={"d": d0, "s": s0, "m2": m20},   # one numeric slice
    forbid=[("T", (1, 0, 0))],                         # masters, by label
    expect_counters={"eqs": 7135, "terms": 42286})
res = lap.eliminate(system, lap.Schedule(policy="B2FT"))
```

Because the format is documented there is nothing empirical to pin:
Route B has NO version guard and NO weight decode — the
`KiraVersionError` / `ABSENT-UNGUARDED` risk class of Route A does not
exist here. In integral notation the adapter assigns deterministic
Laporta ranks (ascending `(t, r, s, sector, family, indices)`; higher
rank eliminated first) and strata default to `t`; in weight notation the
file's own integers are the ranks. Order/strata/forbid remain the
caller's physics surface either way.

Validated end-to-end scales: 18,553 eqs (vacuum 3-loop family); 7,135
(the reference dispersion family); 16,116; 179,447 (full-family
selection); 454k eqs / 838.7 s / 20 GB per slice at the largest measured
sector.

**Honest caveats, out of our control:**
- `SYSTEM_*.gz` is an UNDOCUMENTED internal kira dump. The weight decode is
  EMPIRICAL and version-pinned (kira 3.1 / kira.db VERSION '2.1' /
  integral_ordering 5 / WEIGHTBITS (5,4,17,13)). The adapter REFUSES a
  kira.db with different constants (`KiraVersionError`) and records
  `ABSENT-UNGUARDED` when no kira.db ships with the tree.
- Kira never exports its weight<->integral dictionary. Ordinal-master
  labels are asserted empirics — "best-effort, receipt-checked", never
  guaranteed; the downstream cross-check is MANDATORY.
- No path feeds results back into Kira (no import API); consumers take
  Winnow tables/receipts directly.
- These caveats are Route A's alone. Route B (`adapters/user_system.py`,
  the documented stable `user_defined_system` format — see its quickstart
  above) has no version-pinned empirics at all. Route C (patching Kira to
  call Winnow) remains NOT REALISTIC: C++ with no eliminator plugin
  interface.

## What stays caller-side, and why (measured reasons — this table is binding)

| stays out | why (measured) |
|---|---|
| Seeding / identity generation | 6.1–24.4% of needed pivot columns can have no covering identity in any staged box — the box must exist and is domain/tool-specific. Library consumes rows. |
| Symmetries / sector maps / dictionaries | Where silent wrongness lives: mixed-convention symmetry bugs, wrong dictionary routes (a measured 25/545 deep-dot defect class, QUARANTINED). Library works on opaque cols, hands `meta` back untouched. |
| Multi-prime orchestration + CRT + rational reconstruction | Certificates are per-(point, prime) by pre-registered scope. Stacking primes/points is orchestration (FireFly-class); putting it inside would smuggle in a symbolic-correctness claim the witnesses do not make. |
| Kinematics / coefficient evaluation | Rows arrive already in F_p. The slice evaluator ships in the Kira ADAPTER, not the core. |
| Choice of order/weights | Hook. The kira adapter supplies kira weights; any caller can supply its own Laporta rank. |

**EXCLUDED (binding):** the `weightdict` dict route (measured-WRONG
on a closed deep-dot class: 25/545 rows rejected by two independent oracles
— quarantined until the mechanism is pinned); `incell`/`syzbound` (syzygy
closure certification — separate tool); `gridmode` (grid-splice assembler —
separate tool).

**Witness modes:** `witnesses="retrofit"` (post-hoc Wiedemann lambda solve
— any policy, any backend) and `witnesses="native"` (in-elimination lambda
tracking via the engine transcript hooks, `certify.py` — B2F/B2FT on the
cpu backend only, the paths with hooks; B2/B3 and the dense backend refuse
loudly). Native tracking multiplies fill, so it measures itself on every
run — registered rows, lambda ops, total/peak lambda nnz land in
`ledger["witnesses"]["transcript_fill"]`, and the battery
(`tests/test_native_witness.py`) re-measures and prints the numbers on its
toy fixtures each time it runs. Either way the certificate is verified
through the vendored receipt core against the caller's rows — trust
attaches to the receipt, never to the recording.

## Why receipts are native, not optional

The recorded confident-wrong-table exhibit: reduction rows that were D-ladder
stable, rank-clean, and 2-prime CRT-consistent — every folklore internal
consistency check PASSED on a wrong answer. Only the lambda certificate
caught it (0/16 rows in the row span). Corroborated by a wrong dictionary
route (25 wrong rows with a complete master basis) and a large sector with
+20 false masters, all checks green. Moral, stated plainly: internal
consistency is not correctness; every reduction table ships with a receipt.
A witness verifies against the CALLER's OWN parse of the rows —
`ibplapper.receipt` is a BYTE-IDENTICAL vendor of this repository's
`tools/trust/receipt` (the receipt member of the trust package; core imported, never forked; `tests/test_vendor_integrity`
enforces sha256 equality).

## Measured record

| item | measured | bar |
|---|---|---|
| Interface-table (T) bank replays, reference family, 4 slices | BYTE-EQUAL via API, B2F AND B2FT | reference banks |
| T replay from a FRESH kira Generate | BYTE-EQUAL (also witnesses Generate determinism) | reference banks |
| retrofit battery | 84/84 oracle-match + certified + verified | 84/84 |
| verify throughput | 0.69 ms/row canonical (0.85–0.86 ms/row loaded host) | ~1 ms/row |
| retrofit solve | ~2.2–3.1 s/row/prime at the 7,135-eq reference family | measured envelope |
| sabotage S1–S6 | 6/6 caught | 6/6 |
| recorded wrong-table exhibit | 6/6 MUST-FAIL failed, 6/6 true-row controls passed | canonical regression |
| mutation test | 4/4 engine mutants CAUGHT, control clean | all caught |

This release ships the pure-python engine with the measured envelope above; a
compiled core port is out of scope.

## Layout

```
pyproject.toml            package `ibplapper` (deps: numpy; test: python-flint;
                          optional dense backend: torch)
ibplapper/
  api.py                  System / Schedule / Bank / Result / eliminate / caps
  engine.py               stratified F_p eliminator (B2/B3/B2F/B2FT,
                          stratified_solve, rref_canonical, b2b3_compare)
  bank.py                 interface table T (append-only, sha256-verified)
  coverage.py             coverage-counted fresh-identity gate
  ledger.py               audit ledger
  witness.py              retrofit lambda emission (glue over vendored emitter)
  certify.py              native lambda tracking: Transcript (engine hook
                          recorder) + extract_lambdas + native emission,
                          per-run fill measurement into the ledger
  receipt/                VENDORED tools/trust/receipt: core/emitter/detector +
                          WITNESS_FORMAT.md + VENDOR_MANIFEST.json
  backends/               optional dense (torch) backend, byte-equal contract
  adapters/kira.py        Route A — SYSTEM_*.gz sidecar: loader + weights +
                          version guard + provenance gate (all kira
                          empirics HERE)
  adapters/weights.py     the version-pinned empirical weight decode
  adapters/user_system.py Route B — documented user_defined_system reader
                          (no version guard: nothing empirical to pin)
winnow/                   `import winnow` alias (re-exports ibplapper)
tests/
  run_tests.sh            full battery (any FAIL = nonzero exit)
  test_vendor_integrity   byte-identity of the receipt vendor (drift is loud)
  test_gates              G1-G4 vs an independent flint RREF oracle + API gates
  test_user_system        Route B reader gates (self-authored fixtures,
                          exact-row + oracle + refusal + provenance legs)
  test_native_witness     native-witness gates (certify + verify + retrofit
                          agreement + fill measurement + MUST-FAIL legs)
  test_m2_coeffs          m2-coefficient evaluator gates vs a Fraction oracle
  test_sabotage           6-class sabotage MUST-FAIL suite (+ archived leg)
  test_dense_backend      dense-backend contract gates (torch; + archived leg)
  test_kira_adapter       version guard MUST-REFUSE + provenance gates (archived)
  test_tbanks             BYTE-EQUAL T-bank replays (archived)
  test_witness84          84/84 retrofit round-trip (archived)
```

Run: `tests/run_tests.sh [scratch_dir]`. The self-contained suites run
as-is; the legs marked "archived" replay reference artifacts not shipped
with this tree and SKIP cleanly unless `WINNOW_BANKED_ROOT` /
`WINNOW_T1_ROOT` point at a copy.

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).

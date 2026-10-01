# winnow — GUIDE

Tool page: https://bootloops.ai/tools/index.html (registry index)

KIND: package (Python library `ibplapper`; `import winnow` is a verbatim
re-export alias). Deps: numpy; test: python-flint; optional dense backend:
torch.

## PURPOSE

Standalone upgraded-Laporta F_p eliminator LIBRARY — Seedling's exact
elimination engine. Sparse identity rows over F_p + an ordering/stratum
policy in; closed reduced rows, per-row lambda-witnesses, and a signed
audit ledger out. Trust attaches to the receipt, not the engine: every
reduced row can carry a lambda-witness verifiable against the caller's own
parse of the input rows.

## USE-WHEN

- Reducing an IBP system where every table must ship a receipt (internal
  consistency is not correctness: the recorded wrong-table exhibit passed
  D-ladder, rank, and 2-prime CRT checks on a wrong table; only the lambda
  certificate caught it, 0/16 rows in span).
- Kira sidecar Route A: run Kira with `run_initiate: true` (Generate+Select,
  no back-sub), then eliminate the `SYSTEM_*.gz` dump yourself.
- Route B: eliminate equation files in kira's documented
  `user_defined_system` input format (`integral*coefficient` term lines,
  blank-line splits; weight or `name[indices]` notation) — no Kira
  install, no version pinning.
- When leftover master-relations must surface, never be silently absorbed.

## NOT-FOR

Seeding/identity generation, symmetries/sector maps/dictionaries,
multi-prime CRT/rational reconstruction, kinematics/coefficient evaluation,
order/weight choice — all caller-side (binding table with measured reasons
in README.md). Certificates are per-(kinematic point, prime) ONLY — no
cross-point or symbolic-table claim; 2-loop measured, 3-loop =
extrapolation. EXCLUDED: the `weightdict` dict route (measured-WRONG,
quarantined), `incell`/`syzbound`, `gridmode`. `witnesses="native"`
(in-elimination lambda tracking, certify.py) runs on B2F/B2FT cpu only —
B2/B3 and the dense backend refuse loudly; retrofit is the
policy/backend-independent mode.

## INVOKE

```python
import ibplapper as lap
sys = lap.System(rows, p, order, forbid={...}, stratum_of=None)
res = lap.eliminate(sys,
                    lap.Schedule(policy="B2FT", bank=lap.Bank(path),
                                 caps={"stratum_wall_s": 60}),
                    witnesses="retrofit", witness_targets=[...])
```

Kira sidecar (Route A): `from ibplapper.adapters import kira; system,
sysd = kira.load(art_dir, family, p, d0, eta0, expect_counters=...)`.
User-defined system (Route B): `from ibplapper.adapters import
user_system; system, sysd = user_system.load(src, p, values={"d": d0,
...}, forbid=[...], expect_counters=...)` — src = file or directory of
`.kira`/`.kira.gz`; every coefficient symbol must be declared in
`values` (one numeric slice; undeclared symbols refuse loudly).
Witnesses: `witnesses="retrofit"` (any policy/backend) or `"native"`
(B2F/B2FT cpu; per-run fill measurement lands in
`ledger["witnesses"]["transcript_fill"]`).
Verify downstream: `ibplapper.receipt.verify(my_own_rows_parse,
"w_....json")` — against YOUR OWN parse of the rows, mandatory.
Tests: `tests/run_tests.sh [scratch_dir]`.

INPUTS: sparse rows `dict[col->val mod p]`, ONE prime per System; `order` =
dict col->rank (HIGHER rank eliminated FIRST); `forbid` = master cols;
optional `stratum_of` (higher = earlier). Column ids are opaque ints;
order/stratum_of/forbid/meta are the entire physics surface. Policies
B2|B3|B2F|B2FT are answer-identical (cost knob only).

OUTPUTS: `res.subs` (pivot -> closed row over master/forbidden cols),
`res.leftover` (master relations, never absorbed), `res.witnesses` (lambda
vector, c, verify status, file path per target), `res.ledger` (+
`save_ledger` JSONL append). Bank = append-only interface table T,
sha256-verified read-back.

ENV: caps produce a loud `CensoredError` with a partial ledger at stratum
boundaries, never silence. The Kira adapter is version-pinned (kira 3.1 /
kira.db VERSION '2.1' / integral_ordering 5 / WEIGHTBITS (5,4,17,13)):
it refuses a mismatch (`KiraVersionError`) and records `ABSENT-UNGUARDED`
when kira.db is missing.

## Independent check (standard library only)

The point of a lambda-witness is that anyone can check it with a parser that
imports none of the solver. The normative semantics are in
`ibplapper/receipt/WITNESS_FORMAT.md`; this is the whole check in plain
Python, on the three-row toy system that `tests/test_sabotage.py` (leg A)
certifies and then sabotages:

```python
import json

def verify(rows, w, table_row=None):
    """rows: ordered list of {col: val mod p}; w: a v1 lambda-witness dict."""
    p, t = int(w["p"]), int(w["target"]["col"])
    assert int(w["lam"]["n_rows"]) == len(rows), "lam n_rows mismatch"
    claimed = table_row if table_row is not None else w["c"]
    r = {}
    for i, lv in zip(w["lam"]["idx"], w["lam"]["val"]):
        for col, val in rows[int(i)].items():
            r[int(col)] = (r.get(int(col), 0) + int(lv) * int(val)) % p
    expect = {int(m): (-int(v)) % p for m, v in claimed.items()}
    expect[t] = 1
    strip = lambda d: {k: v for k, v in d.items() if v % p}
    return strip(r) == strip(expect)          # exact; no tolerances

p = 2147483647
rows = [{203: 1, 101: -5 % p}, {202: 1, 203: -3 % p, 102: -7 % p},
        {201: 1, 202: -2 % p, 101: -1 % p}]            # I_203 = 5 I_101, ...
w = {"receipt_version": "1.0", "kind": "lambda-witness", "p": p,
     "target": {"col": 202}, "c": {"101": 15, "102": 7},
     "lam": {"n_rows": 3, "idx": [0, 1], "val": [3, 1]}}
assert verify(rows, w)                                        # true witness: PASS
assert not verify(rows, dict(w, c={"101": 30, "102": 14}))    # tampered c: FAIL
assert not verify(rows, w, table_row={"101": 16, "102": 7})   # wrong claimed row: FAIL
```

One wrong coefficient, one wrong lambda entry or one stray residual column
fails the comparison. A witness file written by `eliminate(...,
witness_dir=...)` is the same object (`json.load` it); the rows must come from
your own parse of the system, in the system's row order.

## GATES

Battery record: 84/84 retrofit round-trip (2 primes); sabotage S1-S6 6/6
caught; wrong-table MUST-FAIL 6/6 wrong pairs FAIL + 6/6 controls PASS; engine
mutation test 4/4 caught; interface-table replays BYTE-EQUAL to the
reference banks via the API. Measured envelope: 454k eqs / 838.7 s / 20 GB
per slice at the largest measured sector; verify 0.69 ms/row canonical;
retrofit solve ~2.2-3.1 s/row/prime. Rerun tests/run_tests.sh after any
change.

## FOOTGUNS

- Do not trust folklore consistency checks (D-ladder, rank, CRT): they all
  passed on measured wrong tables — the receipt is the only correctness
  claim.
- `SYSTEM_*.gz` is an UNDOCUMENTED kira internal; the weight decode is
  empirical and version-pinned — never bypass the adapter guard.
- Ordinal-master labels are asserted empirics (best-effort,
  receipt-checked): witness verification against the caller's OWN parse of
  the rows is MANDATORY downstream.
- Price big witness batches first (retrofit solve is seconds/row/prime).

## What runs publicly (battery: PARTIAL — verified from this copy)

- Self-contained suites, ALL PASS from a scratch cwd: `test_vendor_integrity`
  (the vendored `ibplapper/receipt/` files byte-identical to this
  repository's `tools/trust/receipt` — sha256-enforced), `test_gates` (G1-G4 vs an
  independent python-flint RREF oracle + API contract gates),
  `test_user_system` (Route B reader: exact hand-computed rows, rank/label
  gates, refusal + provenance legs, end-to-end witness verify against an
  independent re-parse — fixtures authored by the test itself),
  `test_native_witness` (native certificates verified through the receipt
  core, retrofit agreement, tamper MUST-FAILs, per-run fill measurement),
  `test_m2_coeffs` (evaluator vs a Fraction-arithmetic oracle; its archived
  corpus leg skips), `test_sabotage` leg A (6-class MUST-FAIL suite), and
  the synthetic legs of `test_dense_backend` (skips wholly without torch).
- Archived-replay legs (`test_kira_adapter`, `test_tbanks`,
  `test_witness84`, `test_sabotage` leg B, the dense golden gate) replay
  reference artifacts not in this repo and SKIP cleanly unless
  `WINNOW_BANKED_ROOT` / `WINNOW_T1_ROOT` point at them.

## RELATED

Seedling (tools/seedling — the certified predictive-staging front door;
winnow is its elimination engine); receipt (tools/trust/receipt — the receipt member of the trust package, the
solver-agnostic standalone certificate CLI; winnow vendors its core
byte-identically); the Kira fork in the sibling repository kira (clone it beside this one: ../kira).

## CREDIT

The elimination is Laporta's algorithm [Lap] on Chetyrkin–Tkachov IBP identities [CT],
run over prime fields in the finite-field tradition of von Manteuffel–Schabinger [vMS],
Peraro [Per16] and Kira 2 [Kira2], with CRT and Wang/Monagan rational reconstruction
[Wang81, Monagan]. Input formats and Generate/Select artifacts are Kira's (Maierhöfer,
Usovitsch, Uwer; Klappert, Lange; Wu [Kira1, Kira2, Kira3]) — we thank the Kira
developers for documenting the user-defined-system format that makes the sidecar route
possible. The SYSTEM_*.gz term-weight packing that `ibplapper/adapters/weights.py` decodes
is Kira's internal convention (Kira 3.1, WEIGHTBITS as recorded in kira.db) [Kira3]; the
decode was established empirically from run artifacts and is version-guarded
(`adapters/kira.py`); no Kira source is included. The λ-witness receipts are ours.
Bracketed keys resolve in REFERENCES.md at the repository root.
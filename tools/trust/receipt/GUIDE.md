# receipt — GUIDE (member of the trust package, `tools/trust/receipt`)

Tool page: https://bootloops.ai/tools/trust.html (the trust package page)

KIND: package member (`receipt.py` CLI; `shim.py` programmatic API; `core.py` verify core;
`emitter.py`; `detector.py`; `adapters/strata.py`; `WITNESS_FORMAT.md` v1.0).
Solver-agnostic standalone CLI; core is stdlib-only — reimplementable from
WITNESS_FORMAT.md alone.

## PURPOSE

Per-(kinematic point, prime) lambda-multiplier certificates for IBP/Laporta
reduction tables: verify any claimed row against the generating system (~1 ms/row),
emit witnesses retroactively (Wiedemann retrofit), detect false-master/structural
defects. Trust attaches to the receipt, not the engine — any solver becomes drop-in
replaceable.

## USE-WHEN

- Any reduction table you must TRUST (the default reduction-table audit: per-(point,prime) row certificates).
- Retrofitting certificates onto an existing table; auditing a system for false
  masters / master relations / uncovered targets / foreign support (A1-A4).
- Catching confident-wrong tables: the archived wrong-table exhibit — internally
  consistent, 2-prime CRT-consistent wrong rows passed every folklore gate; only
  receipts caught it (0/16 certified). Canonical MUST-FAIL regression in
  tests/test_banked.py.

## NOT-FOR

Per-(point,prime) ONLY — no cross-point (Schwartz-Zippel) claim; symbolic tables
out of scope; 2-loop measurements do not derisk 3-loop. Witness certifies row-span
membership in the GIVEN system only — system completeness = detect + provenance
gates, not the witness. Not a solver (per-row solver economics measured >=4,500x
per-row vs whole-table eliminator). Production-scale emit = pre-registered row samples
only (>=24 min/(target,prime) lower bound).

## INVOKE

`python3 receipt.py {verify|emit|detect} --adapter strata --system-dir <kira
artifacts> --staging <staged config> --family <fam> [--table ...] [--witness
'w_*.json'] [--slices p:d:eta,...] [--out ...]`; **generic systems (the public
front): `--system-jsonl` replaces the adapter** (JSONL of sparse rows dict col->val
mod p). Programmatic: `from shim import ReceiptSession`.

Witness JSON v1.0 (p, point, target col, coeff dict c, sparse lam) — full spec in
`WITNESS_FORMAT.md`. Column ids opaque; solver meaning in optional labels only.
rc: 0 clean / 1 any-row-fail-or-alarm / 2 malformed.

ENV: `RECEIPT_STRATA_PATH` — adapter only; by default the adapter uses the
sibling strata member of the trust package (`tools/trust/strata`, `../strata`), and the env var
points it at an external strata tree instead; everything else
(core/verify/emit/detect on JSONL systems) is self-contained.

## GATES

Battery record: verify 0.69 ms/row canonical (0.57-0.96 loaded); retrofit 84/84 (42
rows x 2 primes); emit ~2.2-2.3 s/(row,prime) at the 7,135-eq reference family;
wrong-table must-fail 6/6 wrong pairs FAIL + 6/6 controls PASS; sabotage suite 6 classes caught;
checker mutation-tested 5/5 mutants caught; detect on the real 7,135-pivot system 0
alarms 0.37 s. Rerun tests/run_tests.sh after any change.

## FOOTGUNS

- Passing verify does NOT certify the system itself — always pair with detect +
  staging provenance when the system's completeness is in question.
- Emit uses a column partition staged via one eliminator pass at the anchor
  slice; do not treat emit as production-scale for full
  production tables.
- Tool never writes into live run directories (coordination is read-only).

## What runs without the reference data (battery: PARTIAL)

- `tests/test_core.py` (synthetic + sabotage + rc discipline): ALL PASS from a
  scratch cwd. `tests/mutation_test.py`: ALL MUTANTS CAUGHT. These two are the
  self-contained legs (`run_tests.sh` runs them; its default OUT is a tmp dir).
- `tests/test_banked.py` (30 archived pairs + 84/84 retrofit + wrong-table MUST-FAIL) is a
  reference-data leg: it replays archived artifacts not in this repo and SKIPS
  cleanly unless `RECEIPT_BANKED_ROOT` points at them.
- `adapters/strata.py` pairs with the strata toolchain, shipped as the sibling
  strata member of the trust package (`tools/trust/strata`) and found there by default;
  systems without kira artifacts use `--system-jsonl`.

CREDIT: the tables certified are Laporta-algorithm reductions (S. Laporta 2000, Int. J.
Mod. Phys. A 15:5087) over finite fields (von Manteuffel & Schabinger 2015,
arXiv:1406.4513; Peraro 2016, JHEP 12:030) as produced by Kira (Maierhöfer, Usovitsch &
Uwer 2018; Klappert, Lange, Maierhöfer & Usovitsch 2021; Lange, Usovitsch & Wu,
arXiv:2505.20197); the retrofit emitter uses Wiedemann's sparse solver (D. Wiedemann
1986, IEEE Trans. Inf. Theory 32:54); the witness-per-row design is in the spirit of
certifying algorithms (McConnell, Mehlhorn, Näher & Schweitzer 2011, Comput. Sci. Rev.
5:119).

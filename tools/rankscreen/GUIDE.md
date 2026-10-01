# rankscreen — GUIDE

Tool page: https://bootloops.ai/tools/index.html (registry index)

KIND: package (`modp_rref.py`, `screen.py`, `shim.py`, plus the in-tree
synthetic battery `battery.py` + `fixtures.py` with its `selftest.py` entry;
stdlib + multiprocessing)

## PURPOSE

Multi-prime parallel rank + inconsistency screen for exact-Q Gaussian eliminations of
large row systems. k mod-p RREFs run in parallel (one process per prime); identical
(rank, pivot sequence, inconsistent-row set) across all k primes = CERTIFIED-SCREEN. The exact-Q solve runs ONLY on
closure candidates (rank-full pure block + 0 inconsistent), seeded with the agreed
pivot structure.

## SCOPE

A rankscreen verdict is a screen, not a proof of closure: CERTIFIED-SCREEN
certifies only that k independent primes reproduced identical rank, pivot
sequence, and inconsistent-row set under the pinned eliminator convention
(rows in the order given, minimum-column pivot, insertion-order reduction)
— enough to certify refutations and to decide and seed the exact solve —
while any cross-prime disagreement escalates to the exact solve
unconditionally (there is no vote), and a closure is only ever finalized
by the exact solve itself.

## USE-WHEN

- A closure-engine elimination (3000+ cols, exact-Q, hours single-core) needs its
  rank/inconsistency verdict in minutes.
- You want refutation verdicts (rank deficit, inconsistency sets) certified without
  paying the exact solve at all.
- You have idle cores and a serial exact elimination in front of you.

## NOT-FOR

finalizing CLOSURES — a closure candidate ALWAYS goes to the exact-Q solve (the
screen only decides WHETHER to pay for it and seeds it). Not a general
sparse-linear-algebra library; verdict semantics are pinned to a specific certified
eliminator (row order, min-column pivot, insertion-order reduction). Dense backend
capped at ncols<=8192 and p<2^25 (overflow bounds); sparse backend is the
wide-prime/low-memory fallback.

## INVOKE

- CLI: `python3 tools/rankscreen/screen.py ROWS.json.gz [--k 2] [--backend
  auto|dense|sparse] [--nunk N] [--out RECEIPT.json] [--primes p1,p2,...] [--tag t]`
- Drop-in (adopt at the next natural run, NEVER retrofit a live chain): `from shim
  import gauss_rank_screened` — same call shape and 6-tuple return as an exact
  `gauss_rank_mpq`-style eliminator; screened refutation branch returns
  red=redr=None and piv={pivot_col: original_row_index}.
- Pieces: `shim.serialize_rows / load_rows` (rankscreen-rows-v1 gzip JSON);
  `shim.screen_rows_labeled(rows_all, ncols, ...)`; `shim.exact_seeded(rowlist,
  labels, pivot_rows)` (pivot-rows-first exact with prediction verification +
  automatic plain-order fallback); `screen.screen_rows` (API core);
  `modp_rref.rref_modp{,_dense,_sparse}`.

INPUTS: rows as [({col:int -> Fraction}, Fraction rhs)] (+ parallel labels), or the
serialized rankscreen-rows-v1 form; fixed deterministic prime pools; k>=2 enforced.
OUTPUTS: verdict dict / receipt JSON (rankscreen-verdict-v1) with per-prime receipts.

## GATES (summary — full receipts not shipped)

Positive control (2348x304): exact 304/304 rank / 0 inc reproduced; screen == exact
rank + pivot cols + empty incon set on BOTH backends, closure_candidate=True; exact
77 s vs dense k=4 screen 0.75 s. Recorded refutation controls: identical pivot cols +
incon sets vs exact on multiple real systems (241x like-for-like measured on 864-col
eliminations); planted negatives (perturbed rhs, 0=1 row, p|denominator
auto-replacement, engineered per-prime disagreement → ESCALATE-TO-EXACT) all caught.
Speed on 3064-col systems (5706/5891 rows): screen k=2 wall 35-38 s vs measured
exact lower bound >= 3104 s (>=82-90x).

FALSE-AGREEMENT RISK (honest statement): a wrong verdict at one prime requires p | D
for a specific nonzero minor-class invariant D; a false CERTIFIED-SCREEN needs ALL k
primes bad WITH identical wrong outputs; a false CLOSURE is impossible to finalize
(exact always runs on closure candidates) — the exposed surface is a false
REFUTATION only. Mitigation: raise --k and/or rerun with a disjoint --primes panel;
rank-full at any single prime is already a rank-full proof over Q (rank can only
DROP mod p).

## FOOTGUNS

- NEVER retrofit into a running chain — adopt at the next natural run.
- The screened refutation branch returns red=redr=None: any downstream consumer of
  exact reduced rows must be inside the closure-candidate/exact path.
- Inconsistent-row LABEL SETS are order-dependent; compare verdicts only on
  identically-ordered row lists — the tool always preserves engine row order.
- exact_seeded measured ~parity with plain exact on the 304-col control (its value
  is the prediction cross-check + fallback, not speed); don't sell it as a speedup.
- Dense backend: p<2^25 and ncols<=8192 are hard overflow caps — auto pool selection
  handles this; only --primes overrides can violate it (guarded, raises).
- multiprocessing fork start method (rows shared copy-on-write); nice the CLI
  yourself on shared boxes.

## Battery (in-tree, synthetic, self-contained)

One command: `python3 selftest.py` — runs `battery.py --mutation-controls`
from a scratch directory. Fixtures are built deterministically in code
(`fixtures.py`) with ground truths that follow from their construction, so the
battery needs no external data and runs anywhere the tool runs (~3 s). Legs:

- **planted_control** — known-truth staircase system: CERTIFIED-SCREEN with
  exactly the planted rank, pivot sequence, empty inconsistent set and the
  closure flag, identically on BOTH backends and in agreement with the exact
  eliminator on the same rows.
- **zero_eq_one** — planted inconsistent system (a literal 0 = 1 row plus
  masked shifted combinations): every planted row flagged, on both backends
  and in agreement with exact; never a closure candidate.
- **prime_disagree** — a row engineered so different primes see different
  ranks: the verdict must be ESCALATE-TO-EXACT — including at k = 3 where two
  of the three primes agree (escalation on ANY dissent; there is no majority
  vote).
- **bad_denominator** — a coefficient denominator divisible by a pool prime:
  the prime must be detected, replaced and receipted in `primes_replaced`
  (never a silent wrong verdict), and a fully-unusable prime panel must be
  refused.

`--mutation-controls` copies the tool to scratch, applies ONE targeted
sabotage per leg (wrong pivot choice, disabled inconsistency detection,
forced-agree verdict, disabled bad-prime detection), and requires exactly the
sabotaged leg to FAIL — the mechanisms are proven live, not assumed.

CREDIT: multi-prime modular rank screening is standard exact linear algebra (as in the
finite-field IBP programmes of von Manteuffel–Schabinger [vMS] and Kira 2 [Kira2]); the
certification contract is ours. Bracketed keys resolve in REFERENCES.md at the
repository root.

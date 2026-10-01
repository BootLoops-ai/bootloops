# blade — python library for the Wolfram-free Blade C pipeline

Status: **tool-complete**. Plain python3 + sympy (lazy, assembly
only) + flint (scheme RREF); no heavy deps. C binaries:
`{redg1,fitrel,ssolve,recmod,dynamicrr,dumppoints,fflowcli}` built from the
`no-wolfram-port` blade fork, which lives in the sibling repository
`blade` (clone it beside this one: `../blade`); build
it and point `BLADE_BIN_DIR` at the resulting bin dir (default:
`../blade/bin` relative to this repository's root). No deployment wiring
in this package.

## Modules

- **`formats.py`** — read/write dataclasses for every contract file: `red_*`
  databases, `kin_*` kinematics, `sch_*` schemes, `config/<k>/{in_,tmp_,out_}*`,
  fit tables, ssolve system txt, evallist, `points/eval/degrees .fflow` (flat LE
  uint64), `rec_*_{coeff,mono,part}`, `rrres`, `state`/`flag`, `manifest.json`.
  Layouts cross-checked against the C readers (file:line citations in the
  docstrings; where docs and C disagreed, C won). Writers are byte-exact per
  producer. `classify(path)` + `round_trip(path, scratch)`;
  `big_uint_primes()` regenerates the fflow prime table;
  `assemble_rational_functions()` rebuilds sympy expressions.

- **`pipeline.py`** — gated stage runner (redg1 → fitrel → dumppoints →
  ssolve → recmod → dynamicrr). Measured gate discipline: **rc is recorded
  but never gates** (redg1 always exits 0; fitrel writes flag=1 over an empty
  fit nullspace; ssolve exits 0 with a wrong-size eval; recmod prints
  "failed" with rc 0). The gates that fire are artifact gates: flag content,
  fit-table dimensions, eval byte sizes, recmod coeff/mono/part consistency +
  cross-prime identity, dynamicrr state==2 + rrres count. Failures raise
  `BladeGateError` naming the artifact.

- **`ratrec.py`** — standalone: reconstruct multivariate rational functions
  from a black-box modular evaluator through recmod + dynamicrr. Needs
  fflow's own structured sample points (dumppoints) and EXACT total num/den
  degrees (over-estimates make recmod fail — measured); the default path
  learns degrees via held-out-validated line scans.

- **`extension.py`** — integral-set extension (OperatorExtendG staircase,
  `divide_sub_family_basic`); `scheme.py` — block-scheme construction
  (G1/G2 via flint nmod RREF with the reversed-column pivot convention,
  Tarjan SCC magic loops, SymMap dedup); `ansatz.py` — kin_table +
  cumulative `config/<k>/in_*` generation (ten ordering contracts O1–O10,
  byte-gated 38/38 vs the reference trees); `search.py` — the search loop
  (escalation with measured STOP lines, fitnumber bisection, Kira-expression
  GF(p) database grower, raw-system fflow grower — `FflowSystem` +
  `make_fflow_grower`: a caller-supplied sparse system in the fflow JSON
  grammar, reduced per point by `fflowcli learn`/`evalmany`, no Kira
  anywhere — growing-CRT rational export with typed
  `ReconstructionOverflow` refusal, held-out relation oracle); `adaptive.py`
  — AdaptiveSearch probes (power limits → weights, single-variable probes).

- **`search.py` nana mode** (reduced-analytic search, upstream
  BLAdaptiveSearch semantics): search relations analytic in `nana < n`
  parameters with the rest numerically pinned; every pinned pool carries
  `nana_pins.json` and exports carry a `reduced_analytic` marker —
  relations are valid ONLY at the pins, fail-closed at the producer.
  Validated: gate 1 — `nana=full` is
  byte-path-identical to the classic search (213/213 artifacts identical,
  relations == the reference export); gate 2 — db at nana=1 closes at L10
  vs L12 full, `reduced_analytic` marker present, monomials free of the
  pinned parameter. Gate 3 (dbox-w0 at nana=2, pin msq): **w0 CLOSES at
  L18 in 2.4 s** (vs 1829.6 s at L18 with zero relations through L20 for
  the full 3-var search), pin-generic at 4 distinct pins; GF(p) fit-table
  annihilation vs the reference Kira reduction 256/256 across 8 primes x 2
  pins; rational EXPORT AT PINS emitted at msq=3/2 (4 relations, 27 primes,
  802-bit heights, annihilation 32/32 + off-pin control). CAVEATS
  (measured): slice coefficients are intrinsically tall — at pin 41/23 the
  export refused typed at 12/24/44-prime ceilings (the per-pin GF(p) fit
  tables are the production artifact; export optional); the trimmed template
  support is pin-specific for small pins, so the production recipe is
  search ONCE at a GENERIC (large-height) pin then `fit_at_pin` per
  consumer pin (validated: 13-18 s per pin instance).

- **`semibl.py`** — SolveSemiBL port: closed BT form → FULL SYMBOLIC
  reduction table (kira_table-comparable; open blocks are explicit PENDING
  rows, never silent partial numbers). I/S/R orderings with an executable
  permutation spec; `ptsord` digits only at the C boundary, typed refusal
  >9 params (the `evaluator=` lane covers those, gate-proven
  byte-equivalent); grow-on-demand dynamicrr handshake.

- **`probe.py`** — BTOracle: per-point numeric reduction oracle over a
  stored, sha256-pinned BT form. Exact backend at ANY prime (works on
  partial families: closed blocks solved, declared extras as independents);
  ssolve backend on stored fit tables at fflow primes, self-gated against
  the exact backend per batch. Corrupted family directories are refused by the
  provenance gate.

- **`blade_family.py`** — turnkey CLI (`python3 -m blade.blade_family
  prepare|probe|status`): family spec JSON → full chain → self-describing
  stored family dir (btform.json). Database lanes: kira2math table,
  discovered Kira config dir, stored `red_*` pools, or a raw fflow JSON
  sparse system (`fflow_system` + `fflow_columns` spec keys — no Kira
  anywhere; the held-out gate then replays the system's own sparse solve).
  Partial closures are stored honestly with
  measured stop reasons; typed failures leave nothing half-written.

- **`selftest.py`** — see next section.

## Selftest and the --full regression battery

    cd tools; ulimit -v 32505856
    python3 -m blade.selftest            # quick: 11 rows, ~7 s
    python3 -m blade.selftest --full     # battery: 29 rows

Quick mode: format round-trips on the reference tree (114 files
byte-identical), pipeline rerun on a copy (rrres/rec_* byte-identical to the
reference), ratrec plants + two mutation controls — plus the SELF-CONTAINED
raw-system section (needs no reference artifacts): typed lane-resolution
refusals, then, where the `redg1`/`fitrel`/`fflowcli` binaries are built,
`blade_family prepare` end-to-end on a planted-truth toy sparse system, probe
values from the stored family vs the planted coefficients at a fresh prime,
and a mutated-truth control (without the binaries the end-to-end legs SKIP by
name).  With no reference bank configured the selftest runs just that
section, so a public clone gets a real battery.

`--full` replays every decisive gate from a reference artifact bank that is
not distributed with this repository (set `BLADE_PORT_ROOT` to it; missing
artifacts are NAMED and the replay aborts) — no new searches, everything
deterministic:

| section | replayed gate |
|---|---|
| M | sha256 of all **242** reference artifacts in the bank-side `regression_manifest.json`; missing/changed artifacts each get a loud FAIL row and the replays abort (never a silent skip) |
| A/B/C | the quick sections above |
| D | gate-mutation trio (truncated eval / empty fit with flag left 1 / dynamicrr state tamper) via wrapped binaries — the NAMED BladeGateError must fire |
| E | three-route database md5 agreement (synthetic == Kira-feed == fflow-native) vs the reference `md5_all.txt` |
| F | search: db self-generated search replay (scheme+ansatz+escalation, seed 20260710): sch_* byte-identity vs the reference tree, exported relations content == reference `db_relations.json`, annihilation vs Kira rows at 2 fresh points |
| G | probe: BTOracle probe at 12 points (6 db dense 192/192 + 6 dbox-w1 compositional 384/384) vs the reference Kira tables; typed `ReconstructionOverflow` refusal on the reference w1 tree at 3 primes, nothing emitted |
| H | semibl: identity run on a fresh copy of the reference tree (rec_*/rrres byte-identical, 10 files); reference symbolic dbox table spot-eval 256/256 vs the Kira table; w0 rows still explicit PENDING |

Measured: **29/29 reference-replay rows PASS, 12.1–12.6 s wall** (three runs,
loaded machine); the self-contained raw-system section adds 5 rows,
measured 3.0–4.1 s.
Harness controls, both PASS:
(a) one flipped byte in a scratch mirror of the bank → the battery names the
exact artifact and refuses to replay, exit 1;
(b) a 1-char bug planted in a scratch copy of `scheme.py` → caught by three
named search rows (wrong g1, sch_* diff, escalation stall), exit 1.
The manifest pins reference artifacts only — the C binaries and the blade
modules are deliberately unpinned (they are the objects under test).

## Measured envelope (all walls from one loaded many-core machine; re-measure on your own hardware before quoting)

- **db family: fully closed.** Warm per-point reduction 0.30 ms (python
  exact) / 0.05–0.10 ms (C ssolve, batched) vs 19.7–22.0 s per-point Kira —
  ~7×10⁴ (independently re-measured). Scope: per-point
  re-evaluation over a stored BT form vs a from-scratch Kira run.
- **dbox: symbolic parity 4/8 targets** (w1+w2 closed, exported, annihilate
  fresh Kira exactly; G2 three-way table gate 1536/1536 + 128/128).
  **w0 OPEN**: measured ladder walls L6→L20 =
  0.9 / 7.7 / 29.7 / 67.5 / 356 / 608 / 1830 / 4035 s (8 threads, loadavg
  285–375), zero relations through L20, so all 4 w0 relations live at
  L ≥ 22; the L22 round projected 9054 s > the 9000 s STOP line and was not
  launched. Levers measured, not guessed: divide_sub_family refuted (+420
  out-of-database integrals), adaptive weights (1,1,1) = no-op, 32-thread
  redg1 no help (redg1 search is ~2-core bound); single-variable probes
  collapse the projected close from L28 to ~L22–24. The nana lane above is
  the production route to w0.
- **Deployment probes (both GO):** a production family's sector-452
  gate block — BTOracle == stored exact rows 1872/1872 and == a FRESH
  independent Kira table 288/288, controls detected, ~9 min total; a second family's
  maximal-cut rungs 3p/5p — 252/252 + 3×4032/4032 vs stored and independent
  companions, controls detected, ~6 min total. Cuts enter as zero-sector
  declarations; no semantics change needed.

## Limits (typed refusals or documented, not silent)

- Multi-group (multi-family) operation untested end-to-end; both real
  families run single-group.
- ncons>0 partial-analytic slice solving on the ssolve lane is a typed
  refusal (per-slice database refits not ported); nana pools + per-pin
  instances are the supported route.
- `ptsord` breaks at ≥10 parameters in the shipped C (single-digit parser);
  the `evaluator=` lane covers those systems.
- dynamicrr needs ≥2 primes with consecutive fflow ids 0..k−1; the reserved
  oracle primes BIG_UINT[3]/[4] are never fitted (gate-enforced).
- `scan_degrees` is O(maxdeg³) per line — supply exact `degrees=` beyond
  maxdeg ~50; over-estimated TOTAL degrees make recmod fail (measured).
- Byte-identical writers are guaranteed for the producers seen in the banks;
  Wolfram-written trees may differ in whitespace (readers are
  fscanf-equivalent).

## Footgun ledger

The standing footguns: rc lies everywhere (gate on artifacts); flags lie over empty fit
tables (dimensions gate); recmod rejects over-estimated total degrees; a
corrupt degree-scan point is ABSORBED (the scan is not a corruption
detector — the recmod/dynamicrr gates are); stale dynamicrr `state` desyncs
the handshake (the runner deletes it first).

## Upstream notice

Portions of this package are a Python translation of Blade
(https://gitee.com/multiloop-pku/blade), Copyright (c) 2022 Xin Guan, Xiao
Liu, Yan-Qing Ma and Wen-Hao Wu, released under the MIT License; their
copyright and permission notice is reproduced in [`NOTICE`](NOTICE) in this
directory and in `LICENSE.md` of the fork (the sibling repository `blade`),
and applies to those portions. The FiniteFlow-derived C interface in the fork
is Copyright (c) 2019 Tiziano Peraro, MIT License (`THIRD_PARTY.md` in that
repository).

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).

# trust — GUIDE

Tool page: https://bootloops.ai/tools/trust.html

KIND: package (assembly: sha-pinned vendored lp-syz set + the strata member +
the receipt member + identity fronts + battery).
Battery: 95 checks ALL PASS, 19/19 mutations caught, 3 named skips (reference
receipts not included), ~30 s on a fast workstation.
`vendor/lpsyz` IS lp-syz's home, `strata/` IS the streamed-IBP certifier's
home and `receipt/` IS the certificate-format tool's home — all three ship
INSIDE this package. One command tests the package: `sh selftest.sh`.

## PURPOSE

The reduction-verification triad under one roof — three checks, CORE lineages disjoint
(battery T7C import-graph + sha census; "no shared code" holds at the core level ONLY —
see the substrate caveats), on ONE reduction table: (1) kira crosscheck [strata; Kira
C++/Fermat lineage; kira-side coeff eval = strata CoeffEvaluator], (2) per-row
lambda-witness audit vs the caller's re-parse [receipt CORE stdlib-only; AS-EXECUTED row
re-parse = strata loader/CoeffEvaluator — the named STRATA-LOADER SUBSTRATE caveat,
mitigated by the battery's independent stdlib re-parse spot-check], (3) lp_syz syzygy
oracle [Singular + flint fmpq at rational (d0,eta0); FLINT-substrate caveat →
fraction_oracle leg] — the class 1+2 cannot see: basis rank-deficiency/non-uniqueness
(sec39/sec45/LP(431)).

## USE-WHEN

- A streamed/winnow/strata-lineage eliminator table must move downstream and needs the
  full triad, not one check.
- lp_syz certificates need to travel as WITNESS v1.0 (trust.witness_bridge — the schema
  is used as is, provenance in system.source, verified by the receipt member's core and by a
  second, independently written receipt core that is not part of this repo).
- A FLINT-free exact spot-check of the lp_syz route is required (trust.fraction_oracle —
  pure-stdlib Fraction, bounded FOREIGN subsample/point by design; the head-to-head FLINT
  cross-gate is the Fraction-vs-flint RREF on the same exact rows; the foreign-point leg
  is a Fraction RE-ELIMINATION of the tower's exact rows, tower build + kira eval shared
  with the flint route).

## NOT-FOR

production-scale lp_syz (8-var top-node syz wall: >1200 s, measured;
certifier scope only); cross-point/symbolic-table claims; presenting
strata-vs-winnow agreement as independent (ONE lineage — engine/bank/weights verbatim
lifts; caveat text in trust/strata.py); presenting check 2 AS EXECUTED as "no shared
code" with the eliminator (STRATA-LOADER SUBSTRATE caveat, trust/receipt.py docstring —
the independence claim stops at the receipt CORE).

## INVOKE

`sys.path.insert(0, "<repo>/tools/trust"); import trust; trust.verify()` (import itself
refuses PYTHONPATH/PYTHONSTARTUP/PYTHONHOME worlds typed EnvPoisonError unless the
interpreter ran -E/-I — MEMBERSHIP form, set-but-empty refuses too;
AND refuses typed StdlibShadowError when a security-critical stdlib module resolved
from outside the STDLIB roots (sysconfig stdlib/platstdlib, BASE-installation vars) or
through ANY site-packages/dist-packages segment — the env-FREE sys.path[0]/cwd shadow
class AND the site-.pth reorder class; honest-import perimeter:
__file__ spoofable by an already-hostile interpreter. Child mode -E -P -s -B closes
env + script-dir + user-site-.pth vectors; NAMED residuals: root-owned system-site .pth
(root write required) + the host's own site processing (same-user user-site write =
outside the pin perimeter)); child mode `trust.run_vendored('lp_syz_431', [...])`
(fail-closed TRUST_OUT_ROOT, VALUE-validated; engines reach sympy/flint via an inert
tail append — user site never site-processed);
`trust.oracle_k2disp.reduce_and_compare(table, d0, eta0)` (2-loop fixture adapter;
table text parsed under a LOCKED rational-function grammar, never sympify —
out-of-grammar coefficients refuse typed; FULL parse coverage —
unconsumed records/lines refuse typed NAMING line numbers, never silently dropped;
raw-byte ALPHABET asserted BEFORE parsing — ASCII-printable + newline
only, lookalike bytes refuse naming offset/line; resource bombs refuse at
grammar/build/eval time with AGGREGATE whole-coefficient budgets — literal >10^4
digits, |exponent| >10^4 incl. MERGED symbolic exponents, per-intermediate AND
cumulative bits >2*10^6, op-units >10^5, div0 — zoo never built, eval is a budgeted
Fraction walk, never a hang);
`trust.witness_bridge.emit_lpsyz_witness(...)` (path must resolve under TRUST_OUT_ROOT;
provenance keys producer/schema_note RESERVED); battery
`python3 battery/battery.py <scratch>` (clean env + clean cwd).

OUTPUTS: exact-QQ reductions vs kira tables; v1.0 witnesses.
ENV: Singular + python-flint (vendored legs); receipt cores stdlib-only;
TRUST_LPSYZ_KTAB overrides the kira table path (the built-in default points at a table not shipped in this repo — set it explicitly).

## THE STRATA MEMBER

`strata/` beside the package is the streamed-IBP certifier that check 1
fronts — a self-contained toolchain over kira artifacts:

- `corpus_bench.py` — the end-to-end certification driver: kira Generate of
  the solve box, stratified Route-B F_p solve with per-stratum banking,
  weight<->integral dictionary matched at the anchor slice and held out on
  the others (doubling as the full-table kira crosscheck), fresh identities
  from the independent vacuum IBP generator with a coverage-counted residual
  gate + planted-fault trials, and the corpus oracle row-match.
- `fp_eliminate.py` — kira-free stratified F_p eliminator (B2 Laporta order,
  B2F fast bucketed forward-sub — contract-identical to B2 — B3 Markowitz
  fill-aware; every policy finishes in canonical RREF, so outputs are
  identical by construction).
- `certify.py` — lambda-certificate check (sum_i lam_i R_i == {w:1}+rest
  mod p by sparse accumulation over the ORIGINAL rows, independent of the
  eliminator).
- `loader.py`/`weights.py` — F_p-exact SYSTEM_*.gz loader + the empirical
  64-bit weight decode; `weightdict.py` — the dictionary kira runs;
  `ibp_gen.py` — the independent vacuum IBP identity generator;
  `coverage.py` — the coverage gate; `interface_table.py` — stored T tables.
- `tests/test_certify.py` — strata's own test leg (pytest, self-contained).

`trust.strata.load()` imports it in place by identity; TRUST_STRATA_ROOT
points the front at an external strata tree instead. All stored strata
values are slice-numeric (NOT-VALID-OFF-GRID). Check 1 end-to-end needs
kira (Fermat backend) on the machine; numpy for `certify.py`.

## THE RECEIPT MEMBER

`receipt/` beside the package is the certificate-format tool that check 2
fronts (`trust.receipt.core()` imports its `core.py` in place by identity;
TRUST_RECEIPT_ROOT points the front at an external receipt tree instead).
It is a self-contained directory and also runs standalone from any location.

- PURPOSE: per-(kinematic point, prime) lambda-multiplier certificates for
  IBP/Laporta reduction tables — verify any claimed row against the
  generating system (~1 ms/row), emit witnesses for an existing table
  (Wiedemann retrofit), detect false-master/structural defects. Trust
  attaches to the certificate, not to the engine, so any solver becomes a
  drop-in replacement. The verify core (`core.py`) is standard-library only
  and reimplementable from `receipt/WITNESS_FORMAT.md` (v1.0) alone; Winnow
  vendors byte-identical copies of `core.py`/`emitter.py`/`detector.py`/
  `WITNESS_FORMAT.md` (its `tests/test_vendor_integrity` enforces equality
  against this member).
- INVOKE: `python3 receipt/receipt.py {verify|emit|detect} --adapter strata
  --system-dir <kira artifacts> --staging <staged config> --family <fam>
  [--table ...] [--witness 'w_*.json'] [--slices p:d:eta,...] [--out ...]`;
  generic systems use `--system-jsonl` (one sparse row dict col->val mod p
  per line) in place of the adapter. Programmatic: add `receipt/` to
  `sys.path`, then `from shim import ReceiptSession`. rc: 0 clean / 1 any
  row fails or an alarm fires / 2 malformed input. `adapters/strata.py` finds
  the strata member at `../strata` by default (RECEIPT_STRATA_PATH overrides).
- LIMITS: certificates are per-(point, prime) only — no cross-point or
  symbolic-table claim; a witness certifies row-span membership in the GIVEN
  system, so system completeness is the `detect` audit plus staging
  provenance, never the witness; not a solver; production-scale emit is for
  pre-registered row samples only.
- TESTS: `sh receipt/tests/run_tests.sh [scratch]` — `test_core` (synthetic +
  six sabotage classes + rc discipline) and `mutation_test` (5/5 checker
  mutants caught) are self-contained; `test_banked` (30 archived pairs, 84/84
  retrofit, the wrong-table MUST-FAIL exhibit) SKIPS by name unless
  RECEIPT_BANKED_ROOT points at the archived reference artifacts. Full member
  guide: `receipt/GUIDE.md`; positioning and measured numbers:
  `receipt/README.md`.

## GATES

pins verify fail-closed (4 files, source shas advisory; LINKED_SHAS live-front
baselines compared ADVISORY-LOUD — linked_drift reported, never raises);
battery T0-T8 + T7C census, mutation-controlled, OVERALL PASS from a scratch directory
(the run's numbers are written to BATTERY_REPORT.json); battery AND `import trust` refuse
PYTHONPATH/PYTHONSTARTUP/PYTHONHOME worlds (membership form) AND shadowed-stdlib worlds
(StdlibShadowError — script-dir/cwd shadow, site-.pth reorder);
TRUST_OUT_ROOT is allowlist-validated and the witness-emit path honors it too;
linked tools' own batteries certify their own bytes.

## FOOTGUNS

inherited lp-syz set verbatim (descent sign, Aut(G) quotient mandatory, sec39-class
quotient-before-compare, the wrong-table two-primes fallacy); tower_relation_rows_modp refuses
primes hitting coefficient denominators; witness col ids are tower J-space colids —
meaning lives in labels/provenance only. The kira TABLE is semi-trusted
input — a sympify call would be a code-execution surface, so parsing is
locked-grammar (rational functions
in d,eta only; anything else refuses without evaluation); never point TRUST_OUT_ROOT at
a pinned reference-receipt store (legal for the allowlist, clobbers the battery's
truth — own-tree guard protects only the importing copy); the battery T4 re-parse
shares the trivialsector sidecar with the loader (common-mode term-drop; backstopped by
the fingerprint-vs-pinned-rollup leg — code-path independent, not input-independent).
The stdlib-shadow class needs NO env var — sys.path[0] (script dir; cwd under
-c/REPL/stdin) sits ahead of the stdlib; import refuses via post-import stdlib
identity, but that is the HONEST-import perimeter (__file__ spoofable by an
already-hostile interpreter) — don't park driver scripts in dirs writable by untrusted
agents; a nonconforming-but-honest table refuses on parse coverage rather than
certifying its parsed subset (fix the table, don't relax the parser); eval_coeff
refuses typed when a denominator vanishes at the evaluation point (pick another
rational point). .pth files execute arbitrary 'import' lines
at EVERY interpreter startup, before any user code — the gate refuses stdlib modules
from any site dir (never a stdlib home) and children run -s (user-site .pth dead;
sympy/flint via inert tail append), but a .pth in the ROOT-owned system site-packages
remains the named residual (root write required), and the HOST interpreter's own site
processing precedes the import gates (same-user user-site write = outside the pin
perimeter); resource budgets are AGGREGATE per coefficient (node/op/cumulative-bit caps
+ merged-exponent walk + Mul-head accounting — chains rebuilt from in-cap pieces refuse
<1 s; eval is a budgeted Fraction walk, never sympy .subs); the raw table must be pure
ASCII-printable + newline — lookalike bytes (U+2192, fullwidth digits, lookalike family
names) refuse BEFORE parsing naming offset/line instead of silently vanishing from the
certified set.

## What runs here (battery: PARTIAL — verified from this copy)

The one command for the whole package (any cwd; writes only under the
scratch dir):

```sh
sh selftest.sh [scratch]     # leg 1: import gates + 4/4 vendor pins + member homes
                             # leg 2: receipt member battery (receipt/tests/run_tests.sh)
                             # leg 3: strata member certify gates (pytest; skips by name without it)
python3 -E -c "import trust; trust.verify()"   # leg 1 alone, from this directory
```

- leg 1 exits 0 on import gates + pins verify clean (vendor 4/4 OK; both
  member homes present; linked-baseline drift, if any, is reported ADVISORY,
  never fatal); leg 2 prints `RECEIPT TEST BATTERY: ALL PASS` (the archived
  `test_banked` leg SKIPS by name); leg 3 passes 5/5 — verified from this copy.
- `battery/battery.py <scratch>` (clean env + clean cwd) → OVERALL: ALL PASS — 95
  checks, 19/19 mutations caught, 3 named skips, ~30 s — verified from this copy.
  The named skips are legs whose reference receipts are not included (the sec63
  value-level leg and the descent-sign catcher; the tadpole value-level control
  skips when TRUST_LPSYZ_KTAB is unset). The
  lp_syz_431 planted controls (one planted violation found; clean-port zero
  violations) run and grade from this tree as-is.
  Elsewhere, legs that consume pinned reference receipts and the pinned kira
  table report their sources missing (reference receipts not included); the
  source-sha comparison is ADVISORY by design and the fail-closed refusal behavior
  (banned-cwd/OUT_ROOT allowlist, env-poison, stdlib-shadow) is fully exercisable
  anywhere.
- The kira-crosscheck leg (check 1) fronts the `strata/` member beside the
  package (THE STRATA MEMBER above) and runs end-to-end wherever kira is
  installed — `strata/corpus_bench.py` is its driver. Checks 2 (receipt core,
  stdlib-only) and 3 (lp_syz, vendored in `vendor/lpsyz` — needs Singular +
  python-flint; FLINT-free spot-check via `trust.fraction_oracle`) run from this tree as-is.
- The banned-cwd / output-root refusal texts name directory classes — they are the
  boundary of what the tool REFUSES to write over.
- Licensing note: `vendor/lpsyz` holds in-house engines (MIT, sha-pinned), not
  third-party code; the directory name only marks that the package executes them
  as pinned, vendored bytes.


CREDIT: check 1 runs Kira (P. Maierhöfer, J. Usovitsch & P. Uwer, CPC 230 (2018) 99; J.
Klappert, F. Lange, P. Maierhöfer & J. Usovitsch, CPC 266 (2021) 108024; F. Lange, J.
Usovitsch & Z. Wu, Kira 3, CPC 322 (2026) 109999, arXiv:2505.20197) with R. H. Lewis's
Fermat (home.bway.net/lewis) and FireFly (J. Klappert & F. Lange, CPC 247 (2020) 106951;
J. Klappert, S. Y. Klein & F. Lange, CPC 264 (2021) 107968); check 3 computes syzygy
modules with Singular (W. Decker, G.-M. Greuel, G. Pfister, H. Schönemann;
singular.uni-kl.de) and exact linear algebra with FLINT (W. Hart, F. Johansson, A.
Ahlbäck and the FLINT developers) via python-flint; the elimination order is Laporta's
(2000). We are grateful to the Kira, FireFly, Fermat, Singular and FLINT developers: the
triad is only meaningful because their engines are independently excellent.

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).

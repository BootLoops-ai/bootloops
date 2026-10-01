# TRUST — streamed-IBP certification triad (assembly package)

Layout: package `trust` + vendor/lpsyz (lp_syz's home) + strata (the
streamed-IBP certifier, a member of this package) + receipt (the
certificate-format tool, a member of this package) + battery + selftest.sh.
KIND: package (assembly: sha-pinned vendored lp-syz set + the strata member +
the receipt member + identity fronts to linked tools + battery).
Battery: 95 checks ALL PASS, 19/19 mutations caught, 3 named skips
(reference receipts not included), ~30 s on a fast workstation.

## PURPOSE

Three checks on ONE reduction table (the streamed winnow/strata-lineage
eliminator is the system under test). CORE lineages disjoint — censused by
the battery (T7C); "no shared code" holds at the core level ONLY, see the
two substrate caveats below for where the checks as executed share bytes:

| check | tool | lineage | catches |
|---|---|---|---|
| 1. kira crosscheck | strata full-table + dict held-out | Kira 3.1 C++ + Fermat (third-party); kira-side coeff eval = strata CoeffEvaluator | engine-vs-engine disagreement |
| 2. lambda-witness audit | receipt core vs the caller's re-parse | CORE: stdlib-only pure-python ints; AS-EXECUTED row re-parse: strata loader/CoeffEvaluator (see substrate caveat below) | confident-wrong rows (the wrong-table exhibit class) |
| 3. lp_syz syzygy oracle | Singular `syz` + flint fmpq RREF at rational (d0,eta0) | no code shared with kira/eliminator/receipt (FLINT substrate caveat below) | basis rank-deficiency / non-uniqueness (sec39, sec45, LP(431)) — structurally invisible to 1+2 |

ARBITRATION: strata = end-to-end
certification; receipt = per-(point,prime) row audit; lp_syz = independent
oracle when kira itself is in question. (strata ships as the member tree
`strata/` beside the package — trust/strata.py is its identity front; check 1
runs end-to-end wherever kira is installed, checks 2 and 3 from this tree
alone.)

ONE-LINEAGE CAVEAT: strata and winnow are ONE lineage (engine/bank/weights
verbatim lifts) — strata-vs-winnow agreement is NEVER independent evidence.
Full text: trust/strata.py docstring.

FLINT SUBSTRATE CAVEAT: FLINT underlies both kira128/FireFly and lp_syz's
RREF. Mitigations: check 2's core is FLINT-free; trust.fraction_oracle
is the pure-stdlib Fraction leg — bounded FOREIGN spot-checks of the exact
route (measured: flint-vs-Fraction exact equality head-to-head on the same
exact rows [the FLINT cross-gate proper] + a Fraction RE-ELIMINATION of the
tower's exact rows matching kira at a foreign rational point — tower build +
kira-side eval SHARED with the flint route in that leg).

STRATA-LOADER SUBSTRATE CAVEAT (check 2 AS EXECUTED — named parallel to
the FLINT caveat): the receipt CORE is stdlib-only, but the
verify-path rows in the pilots and battery come from
receipt/adapters/strata.py::load_rows (the receipt member) = strata's OWN loader.load_system
+ loader.CoeffEvaluator — the same strata-lineage parse/eval that stages the
eliminator's rows — and that adapter imports fp_eliminate at module level.
Check 1's kira-side coefficients are ALSO CoeffEvaluator-evaluated, so a
systematic loader/CoeffEvaluator bug is common-mode across checks 1 AND 2;
only check 3 escapes it (real example: the fixture family's float-division
exactness bug, recorded in CoeffEvaluator's own docstring). Never present check 2 as
executed as "no shared code" with the eliminator. Wired mitigations: the
battery's T4 INDEPENDENT stdlib re-parse (sha-seeded FOREIGN SYSTEM-row
subsample, ast-based modular arithmetic, no loader import) and the T7C
import-graph/sha census that names this edge.

## LAYOUT

    trust/            _core (pins verify / exec-loader / child runner /
                      identity alias), _pins (PINS+LINKED),
                      strata.py + receipt.py (identity fronts),
                      oracle_k2disp.py (2-loop lp_syz family adapter),
                      witness_bridge.py (lp_syz -> WITNESS v1.0),
                      fraction_oracle.py (stdlib leg)
    strata/           the streamed-IBP certifier (check 1's engine), a
                      member of this package: corpus_bench end-to-end
                      driver, fp_eliminate stratified F_p eliminator,
                      certify lambda-certificate check, loader/weights/
                      weightdict/ibp_gen/coverage/interface_table, own
                      tests/ — the identity front trust/strata.py resolves
                      here by default (TRUST_STRATA_ROOT overrides)
    receipt/          the certificate-format tool (check 2's core), a
                      member of this package: receipt.py CLI (verify /
                      emit / detect), core.py stdlib-only verify core +
                      witness IO, emitter.py Wiedemann retrofit driver,
                      detector.py structural audit, shim.py programmatic
                      API, adapters/strata.py (finds ../strata by default;
                      RECEIPT_STRATA_PATH overrides), WITNESS_FORMAT.md
                      v1.0, own tests/ — the identity front
                      trust/receipt.py resolves here by default
                      (TRUST_RECEIPT_ROOT overrides); winnow vendors
                      byte-identical copies of core/emitter/detector/
                      WITNESS_FORMAT (test-enforced on its side)
    selftest.sh       one command for the whole package: leg 1 import
                      gates + vendor pins + member homes, leg 2 the
                      receipt member battery, leg 3 the strata member
                      certify gates
    vendor/lpsyz/     the vendored certifier set (in-house engines, MIT,
                      sha-pinned — not third-party code):
                      lp_syz.py (KTAB via env
                      TRUST_LPSYZ_KTAB; int_max_str_digits raised), lp_syz_431.py,
                      lp_syz_prod.py, bessel_oracle2.py
    battery/          adversarial battery battery.py:
                      T0 package pins pre-import / T1
                      vendor pins / T2 tamper copies / T3 restore-
                      verify + planted lp-syz traps / T4 joint fixture bars
                      + live bounded legs + INDEPENDENT stdlib re-parse /
                      T5 witness round-trip (CPU-floor stub detector) / T6
                      Fraction-oracle isolated child / T7 identity fronts +
                      caveat asserts / T7C lineage census (import-graph +
                      sha, as-executed) / T8 refusal battery (allowlist
                      posture, relocated-copy probes, typed fixture
                      refusals; cwd/script-dir stdlib-shadow
                      refusal, parse-coverage refusals naming line numbers,
                      resource-bomb caps — each with positive controls AND
                      scratch-copy mutation controls); mutation-controlled
                      throughout; refuses PYTHONPATH/PYTHONSTARTUP/
                      PYTHONHOME worlds at the preamble (MEMBERSHIP form —
                      set-but-empty refuses too); U1-U7
                      declared-unwired with measured numbers

## INVOKE

    import sys; sys.path.insert(0, "<path-to>/tools/trust")
    import trust
    trust.verify()                                   # fail-closed pins gate
    trust.run_vendored("lp_syz_431", [...])          # child mode; TRUST_OUT_ROOT required (VALUE-validated)
    from trust import oracle_k2disp as ok            # check-3 on the fixture family
    r = ok.reduce_and_compare(table, d0, eta0)       # exact-QQ vs kira table
    from trust import witness_bridge as wb           # lp_syz -> v1.0 witness
    wb.emit_lpsyz_witness(path, rows_modp, tcol, c, p, family=..., point=...,
                          provenance={...})          # provenance in system.source
                                                     # path must resolve under TRUST_OUT_ROOT;
                                                     # producer/schema_note provenance keys RESERVED
    python3 battery/battery.py <scratch>             # battery (~30 s quiet, minutes loaded; clean env + clean cwd — refuses PYTHONPATH worlds AND shadowed-stdlib worlds)
    sh selftest.sh [scratch]                         # whole-package public legs (verify + receipt member tests + strata member tests)
    python3 receipt/receipt.py verify|emit|detect ... # the receipt member CLI (receipt/GUIDE.md)

`import trust` itself REFUSES typed (EnvPoisonError) in a PYTHONPATH/
PYTHONSTARTUP/PYTHONHOME world unless the interpreter ran -E/-I (the
refusal applies at package import — in-process verify()/load_vendored/
driver code paths are enforced, not just the battery; membership test,
set-but-EMPTY vars refuse too), and REFUSES
typed (StdlibShadowError, naming the offending path) when any
security-critical stdlib module it rides resolved from OUTSIDE the STDLIB
roots — the env-FREE sys.path[0] shadow class: script dir in script mode,
cwd under -c/REPL/stdin — or through ANY site-packages/
dist-packages dir (a 2-line .pth reorders sys.path at
interpreter startup, before any user code; site dirs are never a stdlib
home). The identity check authenticates __file__ against sysconfig
stdlib/platstdlib under BASE-installation vars (builtin/frozen modules
exempt); it is the HONEST-import perimeter — __file__ is spoofable by an
already-hostile interpreter. Child mode (run_vendored: -E -P -s -B, engines
reach sympy/flint via an inert tail append) closes the env, script-dir AND
user-site-.pth vectors; the NAMED residual is a .pth in the ROOT-owned
system site-packages (arbitrary code at every interpreter startup — root
write access required), plus the host process's own site processing, which
runs before the import gates (same-user user-site write = outside the pin
perimeter).

## WITNESS SCHEMA

WITNESS v1.0 (receipt/WITNESS_FORMAT.md, the receipt member), used as is. Syzygy provenance
rides the extension slot: `system.source` (JSON-encoded producer/tower/point) +
optional family/point/labels. lp_syz witnesses verify through the receipt member's
core AND a second, independently written receipt core (not part of this repo)
with ZERO schema changes (14/14 through both cores, 3/3 planted faults caught).

## VALIDATION

- Restore-verify of the vendored engines: 4/4 pinned controls reproduced
  (tadpole control; sec63 tower field-equal incl. per-target diffs; ctrl452
  negative; ctrl45 positive, violation support matches the pinned record).
- Fixture family: k2disp.
- Check 3: the lp_syz 2-loop adapter exact-matches the pinned kira reference
  table — P0 {6,7}; P1 {2,5,6,7} 7/7; FULL {1..7} ALL 42 targets 42/42 EXACT
  at TWO rational points (sigma=-1, 0 violations; the 2-loop envelope is NOT
  the 3-loop 8-var wall).
- Three checks on one table: emit table-mode 84/84 CERTIFIED, 0
  CLAIM_MISMATCH at 2 primes; dual-core verify vs own re-parse 84/84 + 84/84;
  tamper MUST-FAIL caught.
- Witness round-trip: 14 lp_syz v1.0 witnesses (7 targets x 2 primes), 14/14
  PASS through BOTH cores; planted faults 3/3 caught; Fraction leg:
  flint-vs-Fraction exact equality head-to-head on the same exact rows
  (437x182 bounded block — the FLINT cross-gate proper) + Fraction
  RE-ELIMINATION of the tower's exact rows EQUAL to kira at a foreign 3rd
  rational point (tower build + kira eval shared with the flint route in that
  leg).
- 8-var production-adapter limit (one bounded probe, nothing built on it):
  sec503 top-node syz did not finish within 1800 s at RSS ~0.11 G — a time
  wall, not memory.

## NOT-FOR

- Production-scale lp_syz reduction (8-var top-node syz wall;
  certifier scope only; >30-min single points refused).
- Cross-point / symbolic-table claims (witnesses are per-(point,prime)).
- Presenting strata-winnow agreement as independent (one lineage).
- Presenting check 2 AS EXECUTED as "no shared code" with the eliminator
  (STRATA-LOADER SUBSTRATE CAVEAT above — the independence claim stops at
  the receipt CORE).
- 2-loop fixture measurements do not derisk 3-loop scale (receipt scope pin).

## FOOTGUNS

Inherited verbatim from the vendored certifier (lp-syz manual): descent sign
J_S[R] = -sum J_child (tadpole controls are BLIND to it); Aut(G) quotient
MANDATORY; sec39-class rank-deficient bases make reductions NON-UNIQUE —
quotient before coefficient-level comparison; two-primes-agree does not
certify a table (the wrong-table exhibit). Specific to this package: tower_relation_rows_modp refuses primes that
hit coefficient denominators; witness col ids are TOWER colids (J-space) —
solver meaning lives in labels/provenance only, per the schema.
Env perimeter: NEVER run the battery (or any pin-trusting code) in a
PYTHONPATH/PYTHONSTARTUP/PYTHONHOME world — a pins-aware stdlib shadow
can green a tampered vendor; the battery refuses such worlds and children
launch -E -P -B env-scrubbed. TRUST_OUT_ROOT is ALLOWLIST-validated
(scratch-class roots only) — the configured deny list refuses everything
else; the own-tree guard derives from the package's actual
location (travels with relocated copies). reduce_and_compare refuses junk
tables and absent targets TYPED. Env-less lp_syz KTAB default is sha-pinned.
Import + input hardening: `import trust` refuses the SAME poisoned worlds
typed (EnvPoisonError; -E/-I exempt) — the in-process package code paths are
enforced, not just the battery. The kira TABLE
is SEMI-TRUSTED INPUT — sympify on it would be a code-execution surface: parse_kira_table/
eval_coeff never sympify table text — a locked rational-function grammar
(int literals, d, eta, + - * / and integer-literal powers) builds sympy
objects directly; ANY out-of-grammar coefficient refuses ValueError without
evaluation. LINKED-tool drift is compared (advisory
LINKED_SHAS baselines in _pins.py; verify() reports linked_drift loudly,
never raises — re-baseline only with a recorded reason). The witness
EMIT path honors the output-root rule (path must resolve under the validated
TRUST_OUT_ROOT; env-less/out-of-root emission refuses OutputRootError before
anything is written) and the bridge's system.source producer/schema_note
stamps are RESERVED (colliding provenance keys refuse — forged-lineage guard).
ALLOWLIST NUANCE (by design but sharp): an allowlisted
scratch prefix can CONTAIN trust-critical bytes — the pinned reference-receipt
store and any sibling package
copy under it are LEGAL TRUST_OUT_ROOT values; the own-tree guard protects
only the IMPORTING copy's tree. Never point TRUST_OUT_ROOT at a pinned
receipt store: a misdirected run can clobber the battery's truth — detection
is loud after the fact (pin mismatch bricks T1/T4) but the receipt bytes may
be unrecoverable. T4 SIDECAR COMMON-MODE (backstopped): the
battery's "independent" stdlib re-parse reads the SAME
sectormappings/<fam>/trivialsector sidecar the strata loader uses for term
drops — a wrong/tampered sidecar moves both parsers identically and the
equality comparator cannot see it. Backstop one leg earlier: changed rows
change the system fingerprint, and T4 verifies the rollup-PINNED pilot-B
witnesses against that fingerprint. Read T4's independence as CODE-PATH
independent, not fully input-independent.
Stdlib-shadow class: it needs NO env var — sys.path[0]
(script dir; cwd under -c/REPL/stdin) sits ahead of the stdlib, so a
hashlib.py planted where you run from can green a tampered vendor with
a clean environment. `import trust` refuses typed StdlibShadowError via
post-import stdlib IDENTITY vs the interpreter-owned roots — an
authentication of the HONEST-import case (__file__ is spoofable by an
already-hostile interpreter; that residual is NAMED, not claimed away —
child mode from a clean cwd has no such caveat). Scratch dirs are
exactly where drivers get run from and may be writable by others:
do not put driver scripts next to untrusted files. Table parse asserts
FULL COVERAGE: any '->' record or entry-body line the locked grammar cannot
consume refuses typed NAMING the line numbers — out-of-grammar bytes are
never silently dropped, so a nonconforming-but-honest table refuses rather
than certifying its parsed subset (fix the table, don't relax the parser).
In-grammar resource bombs refuse typed at grammar/build time (literal
>10^4 digits, |exponent| >10^4, numeric-power result >2*10^6 bits,
identically-zero denominators, 0^negative — sympy zoo is never built);
eval_coeff refuses typed when a denominator vanishes AT the evaluation
point (pick another rational point). Env gates use MEMBERSHIP: a
set-but-empty PYTHONPATH refuses too (empty-var path-inertness is
undocumented interpreter behavior and is not relied on).
Site/.pth + resource hardening: the stdlib-
identity gate authenticates against STDLIB roots ONLY — site.py executes
.pth 'import' lines at interpreter startup BEFORE any user code, so a
2-line .pth in site-packages can reorder sys.path and let a planted
hashlib win import from a whitelisted root; site-packages/dist-packages
are never a stdlib home, whoever owns them, and a venv's own lib dir is
not one either (the venv sysconfig scheme would prefix-cover the venv's
site-packages — the gate uses BASE-installation vars). Children launch
-E -P -s -B: the user-writable user site is never site-processed (its
.pth files never execute); sympy/flint ride the launch bootstrap's inert
TAIL append. NAMED residuals: a .pth in the ROOT-owned system
site-packages still executes at every interpreter startup (root write
required); the HOST process's site processing runs before the import
gates, so same-user user-site write remains outside the pin perimeter
(child mode narrows the exposure; it does not sandbox the host). The
per-Pow resource caps are AGGREGATE: whole-coefficient ast-node/
op-unit/cumulative-bit budgets, Mul-numeric-head accounting (sympy merges
coefficients there and eagerly distributes integer powers over them), and
a merged-exponent walk — Mult/Div-chain and symbolic-merge rebuilds of the
refused bombs refuse typed in <1 s instead of grinding 6-60+ s at
multi-GiB; eval_coeff is a budgeted stdlib-Fraction tree-walk (never
sympy .subs on the table path) — typed refusal, never a hang, and ~2x
faster on the real table. The table byte ALPHABET is asserted on the raw
bytes BEFORE parsing: ASCII-printable + newline only, any other byte
refuses typed naming offset/line — lookalike records (U+2192 arrows,
fullwidth digits, lookalike family names) would match neither coverage
counter and vanish silently; a non-ASCII-but-honest table refuses loudly
(fix the table, don't relax the parser; CRLF still parses via explicit
normalization).

## GATES

Battery OVERALL PASS from a scratch directory: 95 checks, 19/19 mutations
caught, 3 named skips (reference receipts not included), ~30 s on a fast
workstation (the run's numbers are written to BATTERY_REPORT.json, not this file).

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

# BALLER — arbitrary-precision ball arithmetic for science

Tool page: https://www.bootloops.ai/tools/baller.html

KIND: package (`baller/` + `vendor/` + `battery/`). Deep manual: `MANUAL.md` in this
directory.

PURPOSE: Arbitrary precision ball arithmetic for science — the certified ball instruments behind ONE front door: the user-facing run/solve/render layer,
certified quadrature (posq), certified period/ODE transport (kklt engines),
Krawczyk/interval-Newton certification, dual-path halt-not-average tripwires +
mutation-tested-gate harness, and precision-hygiene lints/discipline. Ball-arithmetic
ONLY. BALLER never aliases/shadows Arb's namespace, but it OWNS the user-facing
ball-arithmetic layer that builds with it.

USE-WHEN:

- You want a user calculation (expression or iteration) executed in ball arithmetic,
  a value certified to N digits with automatic precision escalation, or fail-closed
  digit printing — WITHOUT hand-rolling the radius-watch/escalate loop (the front
  door: `baller.run` / `baller.solve` / `baller.render`).
- You need a certified enclosure instrument (quadrature evidence, ODE endpoint
  transport, unique-zero proof) without rebuilding bespoke machinery.
- You are wiring a two-route computation and want the halt-not-average tripwire or a
  mutation-tested gate instead of hand-rolling either.
- You want the documented ball-footgun catalog ENFORCED (ambient-dps, fixed-wdps,
  unkeyed caches, unrestored flint ctx, radius blowup, mpf re-round bridge).

NOT-FOR: Taylor-form enclosures; wrapping/renaming Arb;
Fraction/Decimal/str inputs to march (raw TypeErrors — use float/arb/acb/fmpq).
The vendored certlane member is a research re-implementation used to study numerical
certification of a published hazard computation (Selva et al. 2021, INGV matPTF). It
is not a tsunami-warning system, is not endorsed by INGV, and must not be used
operationally.

INVOKE: `sys.path.insert(0, ".../tools/baller"); import baller; baller.verify()`
then `from baller import quad, transport, certify, contract, hygiene`. FRONT DOOR
(top level): `baller.run(f, x0, n, dps)` / `baller.solve(f, target_digits, x0=, n=)`
/ `baller.render(x, digits=, strict=)` / `baller.certified_digits(x)`. Battery:
`python3 battery/battery.py` from scratch space — it REFUSES in-tree and
deny-listed working directories (fail-closed refusal, ships as-is).

Dependency: `python-flint` (Arb). Vendored engines and members (`vendor/kklt`,
`vendor/geo`, `vendor/certlane`, `vendor/ling_onesided`, `vendor/lgf_monodromy`)
are sha256-pinned — 22 files; `baller.verify()` is fail-closed on any mismatch.

## Aliased members (identity imports — code beside this package)

BALLER aliases two registered tools IN PLACE by identity (never copies): it
imports them from the tools root beside this package. In this repository that root
is `tools/` (this package's parent); `BALLER_TOOLS_ROOT` overrides it. The members
and their registry rows:

- **posq** — `tools/posq/` · certified evidence by exact positive quadrature:
  compute a Bayesian evidence integral as a certified positive VALUE-sum (degree-
  matched exact Gauss sum on the raw likelihood product + log-space patch
  enclosures), never expanding coefficients, never subtracting — two-sided
  machine-width certificates from positivity alone. Aliased as `quad.posq` /
  `quad.kernel_io`.
- **wayfinder** (display name Wayfinder) — `tools/wayfinder/` · the eps-graded
  DE-transport layer: packaged DE systems into one evaluator contract, fixed-ε /
  η-march transport, Frobenius landings, ε-Laurent extraction, acceptance gates + manifests. Aliased under `contract` (gate/manifest).
## Hygiene member: dps_lint (mpmath precision linter)

`baller/dps_lint.py` — this package is its only home (exposed as
`baller.hygiene.dps_lint`, the same module object as `baller.dps_lint`).

PURPOSE: AST-based linter for mpmath precision footguns — module-level mp
constants evaluated at import time at the default dps=15 (the class that capped
a production evaluator at 17.6 digits, twice), plus cross-module dps-reset
rules over the whole linted file set, plus complex-step `mp.diff` on
real-branched callables. Stdlib-only; the linted files are parsed, never
executed, so mpmath need not be installed to lint.

USE-WHEN: before any ≥30-digit gate on an mpmath evaluator; when results change
with import order or a module import silently resets `mp.mp.dps` mid-run; when
held-out digits plateau as dps grows.

NOT-FOR: runtime guards inside functions (the fixed-wdps class below is a grep
discipline plus `hygiene.lint_fixed_wdps`, not this linter). Does not lint Julia.

INVOKE (three equivalent doors):
`python3 tools/baller/baller/dps_lint.py FILE_OR_DIR [FILE_OR_DIR ...]` ·
`python3 -m baller.hygiene lint FILE_OR_DIR [...]` (with `tools/baller` on
`sys.path`/`PYTHONPATH`) · `from baller import hygiene;
hygiene.dps_lint.main([...])` / `hygiene.dps_lint.lint_file(path)` /
`hygiene.dps_lint.lint_cross(files)`. Exit 1 if any findings, exit 2 on a
nonexistent path or a non-.py file argument, 0 clean. Cross-module rules run
over the whole given file set, so lint the import closure together, not files
one at a time. Directories recurse.

BATTERY: `python3 tools/baller/baller/dps_lint.py --selftest` (or
`python3 -m baller.hygiene selftest`) plants one finding per rule in toy
fixtures under a temp dir (honors `TMPDIR`), asserts per-rule counts, the
exit-code contract (1 planted / 0 clean / 2 bad path / 2 non-.py) and a
fixed-fixture control; exit 0 all-pass. The package battery runs it as leg L21.

RULES:
- `[module-level]` (and `default-arg in f()` / `decorator on ...` / `class-body
  ...` contexts) — `mp.mpf(...)`/`mp.mpc(...)`/numeric calls/lazy constants
  (pi, euler, catalan, ...) evaluated at import time before the first
  module-level dps assignment (e.g. `T = mp.mpf(-1)/3` at module level).
- `[imported-sets-dps]` — module-level `mp.mp.dps =` in a module imported by
  another linted file: the import RESETS the importer's dps; re-set dps after
  every such import.
- `[import-frozen-const]` — module-level mpf/mpc in an imported module: frozen
  at import-time dps even though the module sets its own dps; a main that
  raises dps later still inherits the frozen value.
- `[pre-main-const]` — module-level literal whose only protecting dps
  assignment sits under `if __name__ == '__main__'`: as an import the guard
  never runs (dps=15).
- `[mpdiff-real-branched]` — `mp.diff(f, ...)` defaults to a complex step:
  garbage (~1e96) on real-branched integrands (asserts/min/max/sorted/ordered
  `if` inside f), silently. Use an explicit-h central difference. Heuristic;
  false positives accepted.

LIMITS (adjacent class the linter does NOT catch — sweep manually or use
`hygiene.lint_fixed_wdps`): fixed working-precision guards/thresholds
(`workdps(<literal>)`, `mpf('1e-K')`) that must scale with `mp.mp.dps`
(cancellation guards need ~2·dps). Symptom: held-out digits plateau as dps
grows, cross-dps agreement constant (e.g. 40v48 = 48v60), or ZeroDivision only
at high dps. Two-precision self-agreement cannot catch it (shared plateau);
diagnose the integrand pointwise at two dps before re-running nodes. Deep
tanh-sinh nodes also round x to endpoint roots exactly — express vanishing
sqrt factors analytically in the substitution variable. Reach limits shared
with the other lints: named-constant dataflow, `setattr`, alias token games.

## Hygiene member: purity_scan (finite-precision-literal scanner)

`baller/purity_scan.py`, the standalone finite-precision-literal scanner for a
numerics purity bar ("no finite-precision numbers on the value path"). It
enumerates and pre-classifies every literal — MPF-OF-FLOAT (`mp.mpf(<float>)`,
the classic disease signature; exact dyadics are false positives), DEC-STRING
(decimal strings ≥8 sig digits), FLOAT-EXPR (float literal in a value
expression), PURE-DPSREL / PURE-CONV (auto-cleared convergence context), plus
FLOAT-EXPR-SHORT and PARSE-ERROR for C sources — so a human adjudicates only
the residual (RAT-STRING and PURE-INT-RAT are recognized-exact classes, not
reported). Report-only, rc=0; python via AST walk + regex, `.c`/`.h` via
regex. CLI (positional, no flags): `python3 baller/purity_scan.py
<list-of-abs-paths.txt> <out.json>`. The scan cannot see build-time constant
freezes behind clean sources — pair it with re-evaluation of stored constants
at unusual dps vs their stored values.

## Vendored members (sha-pinned, loaded as plain modules)

Beside the kklt/geo engines, three member families ship under `vendor/`
(all 22 vendor files pinned in `baller/_pins.py`; battery leg L20 keeps
them live): **certlane/** — certified arb-ball enclosure primitives for the
matPTF alert computation, with a 26-case containment pytest suite and a
four-surface adversarial harness; **ling_onesided/** — the
positive-integrand one-sided evidence kit (centered forms, priced-lever
gates, certified tail sups, peeling-circuit evidence engine, one-sided
fixed-grid route) with its `QUARTET_COLLAPSE.json` fixture vendored beside
it; **lgf_monodromy/** — `valmono.py` validated ball-Taylor Fuchsian
monodromy transport with the pinned `BC_cache.pkl` coefficient cache
(KNOWN GAP: loop geometry not vendored — pass explicit waypoints to
`transport_loop`). Member outputs default to the CURRENT working directory
(env-overridable), never `vendor/` — unpinned files there trip
`baller.verify()`. Full contracts: MANUAL.md, "Vendored members".

GATES: acceptance battery 21 legs OVERALL PASS, ~2 min (L17
front-door Muller gates + L18 front-door mutation controls 4/4 + L19
purity-scan planted 3-class smoke + L20 vendored fold members: certlane
pytest 26/26, ling/lgf import + null-loop identity probes + L21 dps_lint
member selftest 11/11 and `-m baller.hygiene lint` exit codes) incl.
tamper/pyc/shadow mutation controls, march mutation-tested 3/3, and the
L12-L14 integrity/minted-check/tripwire guard legs.

The battery is 21/21 OVERALL PASS against this repository's tools tree as
shipped (posq, wayfinder and eras beside this package). Leg L6 shells out to posq's OWN adaptation
battery, which imports the eras battery beside it
(`tools/eras/verify_eras_adversarial.py`). posq's C kernel is a prebuilt
binary: on a machine where it does not execute (different FLINT
shared-library version or CPU architecture), posq rebuilds it from source at
first use, and with no C toolchain L6 reports a loud named SKIP instead —
the battery stays green either way (see `tools/posq/GUIDE.md`,
Requirements). The tools root in `baller/_core.py`
defaults to this package's grandparent with `BALLER_TOOLS_ROOT` override,
and `battery/battery.py` pins the real
root into the environment so its scratch-copy probe subprocesses inherit it.

FOOTGUNS:

- Run baller.verify() first; vendored loads are fail-closed + thread-locked;
  bytecode cache never an input.
- tripwire REJECTS balls/complex typed — compare mids/parts, watch radii with
  RadiusWatch; NaN/inf never pass a bar.
- ftrim only after ACCEPTED steps; kklt engines mutate flint ctx without restore —
  wrap in hygiene.ctx_guard.
- bridge trap: mp.mpf() on an existing mpf re-rounds (both direct and tuple form) —
  hygiene.bridge_mpf/make_mpf is the raw path.
- PL/PT configure law: always transport.configure_all, never PL.configure alone.
- KNOWN RESIDUAL: importing quad leaves tools/posq on sys.path (posq's own
  receipt-pinned bytes) — generic names there (kernel_io) can shadow later user
  imports.
- Lints are advisory instruments with documented reach limits (named-constant
  dataflow, setattr, alias token games are out of reach).

See MANUAL.md for the full module contracts (quad.integral_certified LAWS,
certify.block_krawczyk oracle contract, front-door semantics and the documented
divisor-ball deviation).


CREDIT: BALLER is a front door over Arb, the ball-arithmetic library of Fredrik
Johansson (2017, IEEE Trans. Comput. 66:1281), now part of FLINT (W. Hart, F. Johansson,
A. Ahlbäck and the FLINT developers; flintlib.org), reached through python-flint (F.
Johansson, O. Benjamin and contributors); every certified digit in this package is
theirs first, and we are grateful to them. Certification: Krawczyk (1969, Computing
4:187), Moore (1977, SIAM J. Numer. Anal. 14:611), Rump (2010, Acta Numerica 19:287).
The front-door reference figure is J.-M. Muller's recurrence (Muller et al., Handbook of
Floating-Point Arithmetic, Birkhäuser 2010). certlane re-derives, in ball arithmetic,
layers of INGV's matPTF (Selva, Lorito, Volpe, Romano, Tonini, Perfetti, Bernardi et al.
2021, Nat. Commun. 12:5677; code github.com/INGV/matPTF) — no matPTF code is copied and
INGV does not endorse this re-implementation; ling_onesided's counts derive from DravLex
v1.0 (Kolipakam, Jordan, Dunn, Greenhill, Bouckaert, Gray & Verkerk 2018, R. Soc. Open
Sci. 5:171504; CC BY 4.0). The vendored kklt/geo engines are in-house and ship
sha-pinned; the minted checks encode the documented footgun catalog.

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).

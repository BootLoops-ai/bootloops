# formglue — guide

Tool page: https://bootloops.ai/tools/formglue.html — the external FORM engine's page:
https://bootloops.ai/tools/form.html

KIND: external-wrapper + package (form_oracle.py, form_route.py, form_io.py,
form_hyper.py, tornheim.py, tornheim_w5.py; tests/ with vendored fixtures;
examples/gg2qq/ — the gg -> q qbar FORM worked example)

PURPOSE: glue around FORM 5 — Dirac traces, color algebra, out-of-RAM term-stream
rewriting, arbitrary-precision MZV/euler-sum evaluation (#StartFloat/Evaluate),
integrated GRACE/grcc diagram generator — plus tornheim, the independent
Mordell–Tornheim–Witten / depth-2 MZV cross-check oracle.

EXTERNAL ENGINES (called, never absorbed — obtain upstream):
- FORM 5 with float support: build from github.com/form-dev/form with
  `--enable-float`; point env `FORM5_BIN` at the binary (default: `form` on PATH).
  Distro FORM 4.3 has NO float support and cannot back the oracle (it remains the
  trace cross-check engine).
- HyperFORM (Kardos–Moch–Schnetz, arXiv:2607.01163): clone
  github.com/adamkardos/HyperFORM and point env `HYPERFORM_SRC` at its src/ dir —
  form_hyper refuses loudly when unset.
- GRACE/grcc (T. Kaneko): external diagram generator for the form_route front-end.

USE-WHEN:
- Dirac traces / color algebra / giant-expression rewriting that never holds the
  expression in RAM → form/tform.
- Independent zeta/MZV/euler-sum oracle (PSLQ basket cross-check, ζ-ring gate) →
  `form_oracle.zeta(n,dps)` / `mzv(idx,dps)` / `euler(idx,dps)` → mpmath-ready strings
  ≥dps digits (runs FORM at dps+10; weight cap 22 = highest probed).
- FORM 5 diagrams_ output has NO momentum routing (every internal edge its own label)
  → `form_route.py`: spanning-tree glue, chords keep own labels as loop momenta, tree
  edges solved by exact leaf-strip; emits routing table / FORM `#procedure route<i>()`
  id-block / JSON; 5 rc-coded self-checks EVERY run (per-vertex conservation, Euler,
  externals-once, union-find Betti recount, chord-unit); `--mutate-edge` control hook.
- Any nested-harmonic / Tornheim-type constant (massless propagators, banana ε-coeffs,
  lattice sums) → `tornheim`: `mzv2(s,t)` (Richardson-accelerated nsum, hard-wired —
  no method keyword; the 'r+s+e' route was 100× slower and is not exposed),
  `eulerH(s)`, `T2(a,b,c)`, `MT3(a,b,c,d)` (head + asymptotic tail), signed
  lattice `W3,W4`. All caches dps-keyed.
- 5-leg signed lattice sum W₅(p;q) or 4-fold Mordell–Tornheim MT4(a,b,c,d;e) →
  `tornheim_w5` (member — imports tornheim by identity): `W5(p,q)`
  (30=2⁵−2 sign sectors via lone-leg MT4 + pair conv23), `MT4(...)` (recursion to
  min-index 0, exact-S₃ boundary + EM asymptotic tails), `_clear_caches()`. Walls
  3.6 s (dps60) / 7.4 s (dps120) per W5 value cold.

NOT-FOR: HPL oracles via FORM — `lin_`/`hpl_`/`mpl_` are reserved names with NO
numeric evaluation in 5.0.1; `mzvhalf_` unexercised. FORM is never absorbed into
python — binaries only.

INVOKE: API as above. Read ALL FORM output through `formglue.form_io.read_form_output`
— never raw. Tests: `python3 tests/test_form_io.py` and `tests/test_form_route.py`
(pure python vs vendored fixtures, no FORM needed); `tests/test_form_oracle.py`
(engine legs need a float-capable FORM 5 via `FORM5_BIN`; on a box where FORM
is absent OR present but float-incapable — e.g. distro FORM 4.x, no
`#StartFloat` — those legs SKIP with a named reason stating the version found
and the capability missing, via a once-per-session probe
(`form_oracle.probe_float_support()`); the engine-free legs still run).
Worked example: `examples/gg2qq/` — tree-level |M|^2 for g g -> q qbar in FORM
(`form gg2qq.frm`; `python3 check_gg2qq.py` compares it symbolically with the
textbook result and requires two planted mutations to fail), the FORM 5
`diagrams_` one-loop enumeration (`enum_1loop.frm`, 3 / 109 / 30 / 7), one box
numerator (`box_numerator.frm`), and the routed topology data + renderer
(`gg2qq_topologies.json`, `render_gg2qq_diagrams.py`, matplotlib only). See its
README.md; not collected by the battery.

INPUT CONVENTIONS (empirically pinned):
- mzv_ is outer-index-first (mzv_(2,1)=ζ(3)); euler_ negative index = (−1)^m on its
  own summation variable.
- form_route: node_ args = INFLOWING momenta; edge_(id, f(m), v1, v2) carries m v1→v2;
  labels appear twice with opposite signs. Masses are NOT in the grcc model language —
  distinct named particles + user-side mass table. Quartic couplings must be declared
  squared (grcc coupling power = nlegs−2).

OUTPUTS: mpmath-ready constant strings; routing tables / FORM id-blocks / JSON.

ENV: `FORM5_BIN` (read at call time; no ambient-dps use); `HYPERFORM_SRC`.

GATES: traces byte-identical to FORM 4.3 (the untouched cross-check engine); MZV legs
60-100d vs mpmath/tornheim/CVZ; strict Run/Parse exceptions; selftest with
sabotaged-copy rc=1 proof. Tornheim analytic controls: W₃((1,1,1);1)=π⁴/60,
W₄((1,1,1,1);1)=30ζ₅−12ζ₂ζ₃. HONESTY: the assembled-XA/XF sub-check of the trace
battery was SKIPPED (per-piece check is the claim). tornheim_w5 controls (reference
values not distributed): published D₅=C_{1,1,1,1,1} Laurent (DGV 1502.06698 eq powerd5, taken
verbatim) reproduced 78–79d on all 7 coefficients at dps80 through the full chain; MT4
vs certified brute 15–24d (truncation-limited); W5 vs independent box-lattice brute
9–18d; cross-dps 59d (=dps60 floor); script-mode smoke = reference-value regression,
rc-coded.

FOOTGUNS:
- #write with ~102-char symbol names CRASHES 5.0.1 (rc=134 after a partial line) —
  keep long-symbol fixtures ≤~80 chars.
- #write wraps MID-TOKEN (indented at token boundaries, backslash mid-word) — hence
  read_form_output, always.
- Don't lean on reserved lin_/hpl_/mpl_ names for numerics (silent no-evaluation).

BATTERY (PARTIAL; run from a scratch cwd):
public-runnable WITHOUT FORM: tests/test_form_io.py → 5/5 PASS;
tests/test_form_route.py → 6/6 PASS (both against the vendored tests/fixtures/);
tests/test_form_oracle.py engine-free legs (validation, parse strictness, error
paths) → 4/7 PASS with the 3 engine legs a named SKIP (probe reports the FORM
version found and the missing float capability — an incapable FORM 4.x on PATH
is a SKIP, not a fail).
Runnable WITH a float-capable FORM 5 (`FORM5_BIN`): tests/test_form_oracle.py → 7/7
PASS (verified with a FORM 5.0.1 --enable-float build). Reference-only legs: the
FORM 4.3 trace-identity battery, the HyperFORM check-suite/examples runs (need the
upstream engine), tornheim_w5's reference-value regression (reference values not distributed; the
analytic W3/W4 controls above are public).

RELATED: consumers — lockpick PSLQ baskets (independent oracle leg); MZV cross-checks
vs the galois ring (tools/galois). Cross-check engine: distro FORM 4.3.

CREDIT: FORM is by Jos Vermaseren [FORM0], FORM 4 with Jan Kuipers, Takahiro Ueda and
Jens Vollinga [FORM4], 4.2/4.3 with Ben Ruijl [FORM42], TFORM with M. Tentyukov [TFORM],
and FORM 5 — whose float/Evaluate support this package relies on — by Josh Davies,
Toshiaki Kaneko, Coenraad Marinissen, Takahiro Ueda and Jos Vermaseren [FORM5],
maintained by the form-dev community (github.com/form-dev/form, GPL-3.0); we are
grateful to the maintainers. grc/GRACE diagram generation is Toshiaki Kaneko's [grc].
The router test fixtures (tests/fixtures/) are `diagrams_` output generated with FORM 5.0.1
from the reference manual's example QCD model (gg->qqbar demo); no FORM source is included.
HyperFORM is by Adam Kardos, Sven-Olaf Moch and Oliver Schnetz [HyperFORM] (GPL-3.0;
obtain upstream, not bundled; the form_hyper job skeleton follows HyperFORM's shipped
examples), a FORM port of the core of Erik Panzer's HyperInt [HyperInt]. Because its
driver template follows HyperFORM's shipped example drivers, form_hyper.py is distributed
under GPL-3.0-only (see License below, and NOTICE and LICENSE-GPL-3.0 in this directory).
Mordell–Tornheim–Witten sums are classical. Bracketed keys resolve in REFERENCES.md at
the repository root.

## License

MIT License, except `form_hyper.py`, which is GPL-3.0-only (its driver template follows the example drivers shipped with HyperFORM by Adam Kardos, Sven-Olaf Moch and Oliver Schnetz, GPL-3.0, github.com/adamkardos/HyperFORM; see `NOTICE` and `LICENSE-GPL-3.0` in this directory). Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).

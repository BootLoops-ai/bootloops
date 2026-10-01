# emitall — declarative emit-all harness for project headline claims

## The problem it solves

A project files a results report (or paper section, or blog table) quoting dozens to
hundreds of headline numbers, each derived from a stored receipt. Numbers rot
in three ways: the quote was transcribed instead of emitted, the receipt was
regenerated after the sentence was filed, or the quoted count silently
depends on a universe/register choice the sentence does not state. The
mechanical fix: re-emit EVERY quoted value from the receipts and compare at
the quote's own printed precision, listing every discrepancy as a finding.
Catch classes this surfaces include:

- **an un-universed count**: a report quotes "19/21 cases one-signed"; the
  canonical universe (one census per year, the final year at a single
  register) emits **16/18** — 19/21 arises only when one year's census is
  counted at BOTH of its registers. The number was not wrong; it was
  un-universed, and only a mechanical recount of everything noticed.
- **a printed-precision drift**: the stored skew −0.6845 was quoted as −0.69
  (prints −0.68 at the report's own 2-dp register).

emitall makes that loop standing machinery — the
catch-your-own-filing-in-10-minutes loop every project should be running.
Instead of a per-project script, a project writes one declarative
**claims spec** and runs one engine.

## What it does

```
python3 tools/emitall/run.py SPEC.json [--out report.json] [--quiet]
```

The spec maps each quoted headline to (receipt, field expression, comparison
mode). The engine loads every receipt (JSON or CSV), emits every value,
compares each against its quoted string, prints the table, writes a JSON
findings report, and **exits nonzero on any finding** — mismatch,
outward-band violation, evaluation error, missing receipt, or stale receipt
(mtime ordering). Findings are listed, never smoothed. Report timestamps are
machine-read (`date -u` subprocess), never typed.

## Spec format

JSON (YAML accepted when PyYAML is installed):

```json
{
  "campaign": "...", "register": "...", "quote_source": "which doc is quoted",
  "base": "/abs/path",
  "receipts": {
    "cc":  {"path": "runs/census.json"},
    "dol": {"path": "runs/dollars.json", "newer_than": ["runs/flips.json"]},
    "tbl": {"path": "work/table.csv", "format": "csv"}
  },
  "claims": [
    {"section": "S2", "label": "census rows differing", "quoted": "213",
     "expr": "f('cc','/total_rows_differ')"},
    {"section": "S2", "label": "per-group differing", "quoted": "55/68/90",
     "format": "{a}/{b}/{c}",
     "exprs": {"a": "f('cc','/groups/g3/rows_differ')", "...": "..."}},
    {"section": "BANDS", "label": "gross band $M", "mode": "outward-band",
     "quoted_lo": 66, "quoted_hi": 100, "granularity": 1, "scale": "1e-6",
     "lo_expr": "f('dol','/lo_usd')", "hi_expr": "f('dol','/hi_usd')"}
  ]
}
```

A claim with no `"quoted"` is emitted for the record (status EMITTED, never a
finding). `lets` (ordered name→expression) factor repeated subexpressions;
`note` / `note_format`+`note_exprs` attach riders that travel with the row.

## Field expressions

Python-syntax strings evaluated by a whitelisted AST walker (never `eval`;
no attributes, no subscripts, no lambdas — receipts are reachable only
through the accessors). Full grammar in `expr.py`'s docstring. Highlights:

| need | expression |
|---|---|
| nested key / list index | `f('cc','/peer_groups/peer3/n')`, `f('r','/rows/0/id')` |
| collect + aggregate | `sum(vals('sr','/peer_groups/{peer3\|peer4\|peer5}/ties'))`, `max(...)`, `len(...)` |
| arithmetic combos | `100 * f('cc','/differs') / f('cc','/rated')`, `f('a','/x') + f('b','/y')` |
| count-where | `cw('gaps','','[["off_optimum","==",true]]')` — ops `== != > >= < <= in not_in contains not_contains`, optional `"int"`/`"float"` field coercion, conditions ANDed in order with short-circuit (guard `not_in ["NA",""]` before an `int` coercion), value may be `{"$ref": "alias:/pointer"}` for cross-receipt membership |
| filtered aggregate | `min(where(L, '>', 0))` with `lets: {"L": "vals(...)"}` |
| string surgery | `int(split(f('cc','/s'),'/',0))`, `join(xs,'/')`, `jsondump(x)`, `joinrows('dc','/mixed','{year}-{peer} +{up}/-{down}','; ')` |
| row-wise zip of two list fields | `zipjoin(vals('r','/yrs/*'), vals('r','/regs/*'), '{a}={b}', '; ')` — LOUD on length mismatch |
| substring predicate | `len(where(L, 'contains', 'restated'))`, cw condition `["tag","contains","restated"]` — LOUD on a non-container element |
| exact interval widths | `str(frac('232.4') - frac('154.8'))` → `388/5` — `frac()` parses `"p/q"`, decimal, and scientific strings to exact `Fraction`; flows through arithmetic; `printed` mode formats Fractions at the quote's register, `exact` mode compares `str()` = `p/q` |

Anything not expressible here should be computed into a receipt field by the
project's own emitting script and then quoted — derivation belongs in
receipts, verification in specs.

## Comparison modes

- **`printed`** (default): compare at the quote's own printed precision
  (numbers formatted to the quote's decimal count, `str(int(round(v)))` for
  integer quotes; strings compared after stripping `,` and trailing `%`).
  Optional `"decimals"` overrides the inferred precision. Python booleans
  count as numeric — wrap them in `str(...)`.
- **`exact`**: `str(emitted) == str(quoted)`, no normalization.
- **`outward-band`**: enforces the rounding law mechanically. The printed
  interval must be the OUTWARD rounding of the computed interval at the
  stated granularity: `quoted_lo == floor(lo/g)·g` and
  `quoted_hi == ceil(hi/g)·g`, exact `Fraction` arithmetic (no float
  participates in the accept/reject). Violations are graded **inward-edge**
  (a printed edge cuts into the computed interval) vs **outward-loose**
  (encloses but wider than floor/ceil). Failure shape: a nearest-rounding
  recount prints **3 of 6** one-decimal band edges INWARD ([154.9, 232.3]
  against computed [154.8676, 232.3014]). The rule: display values are
  computed outward by the emitting script, never transcribed from a
  nearest-rounding recount. The case is replayed as an executable receipt in
  the test suite and fires both inward edges.

Receipt-level checks: a missing receipt file and a stale receipt (older than
a declared `newer_than` input, or than `max_age_days`) are findings in their
own right; claims on a stale receipt still evaluate (the staleness itself is
the finding), claims on a missing one are marked MISSING-RECEIPT. Evaluation
errors (schema drift, bad pointer) are LOUD findings (`eval-error`), and the
rest of the spec still runs.

## Example deployment (`specs/`)

- **`toy_stars.json`** over **`toy_receipts/`** — a complete working
  deployment of the spec format on an INVENTED star-rating program (every
  number synthetic): printed-register claims, a percent claim, multi-field
  `format`+`exprs`, a `vals()` pattern recount, a `lets`+`where` count, a
  CSV `cw()` with the NA-guard-before-coercion pattern, two outward-band
  checks (integer and 0.1 granularity), an `exact`-mode claim, a boolean
  wrapped in `str(...)`, and an emitted-only record row — **11 MATCH + 1
  emitted-only, 0 findings**. Executable receipts:
  `tests/test_all.py::TestToySpecReference` runs it clean AND proves it can
  fail — a planted transcription drift fires a `mismatch`, and
  nearest-rounded band edges fire two `inward-edge` outward-band findings.
  Copy this spec as the starting point for your own deployment.

## Template lint: `lint.py`

```
python3 tools/emitall/lint.py <script-or-dir> [...] [--json report.json] [--quiet]
```

Enforces the emission-integrity rule: digit literals in emission prose are
findings, sampled fields must be sorted/seeded — it turns this
error class into a machine-caught one for every project.
Failure shape: an emitted headline carries the
hardcoded prose **"(max 5 cents)"** while the stored receipt's maximum is
**$0.09** — every other number in that string was interpolated from the
receipt, so the one typed digit survives eyeball review (defect class
"typed-not-emitted"; such defects can UNDERSTATE a correct receipt). A
second class:
emission order depending on PYTHONHASHSEED (set/dict-hash iteration feeding
emitted lists). Both classes are reproduced near-verbatim in the test
corpus and must keep firing.

What fires (exit 1 on any finding; 0 clean):

- **digit-literal** — any digit run in the literal text of a string feeding
  an emission sink (variables/keys/calls matching `emit|headline|summar`,
  case-insensitive; assignment, append/extend/add, emit-function args,
  keyword args, dict keys, returns from emit-named functions). Interpolation
  placeholders are exempt: f-string `{...}` fields incl. format specs,
  `.format()` `{...}` fields (`{{ }}` escapes stay literal), %-placeholders
  incl. width/precision digits. Escape hatch:
  `# emitall-lint: allow-literal <reason>` on a line of the flagged string —
  the reason is MANDATORY (a bare pragma is itself a finding,
  `pragma-missing-reason`, and suppresses nothing); allowed literals are
  reported with their reasons, never silent.
- **unordered-iteration** — set iteration (`set()`/literals/comprehensions,
  `&|^-` combinations, `.union()`-family) or dict-view iteration
  (`.keys()/.values()/.items()`) feeding an emission sink, including inside
  f-string interpolations and statement-level `for` loops whose body writes
  a sink — unless wrapped in `sorted()`. Order-insensitive reducers
  (`sum/min/max/len/any/all/set/frozenset/sorted`) are exempt consumers.
- **unseeded-sampling** — `random.sample/choice/choices/shuffle` (module or
  `*.random`) with no preceding `.seed()` in the file, or the same methods
  on an RNG constructed without a seed (`random.Random()`,
  `np.random.default_rng()`); seeded instances and post-seed calls pass.

Honest limits: sink detection is by NAME pattern — an emission list called
`lines` is invisible; quote-transcription args to an emit-named function
are flagged even though
they are targets, not prose (fix: pragma with reason, or migrate the quotes
into a claims spec — spec JSON is data, not linted Python); templates held
in variables are not resolved; a bare dict NAME iterated cannot be
classified; an unresolvable `x.sample()` receiver is not flagged; seeding is
detected per file, not per RNG stream.

## Paper claims battery: `battery.py`

```
python3 tools/emitall/battery.py SPEC.json [--quiet]
python3 tools/emitall/battery.py --selftest
```

The document-side verb: a FROZEN/NEVER/NEG/ABSTRACT assertion battery over a
paper's tex sources and rendered text, driven by one declarative rows spec
(JSON; shape in `battery.py`'s docstring). Row classes:

- **FROZEN** — a pinned headline number/phrase must appear verbatim
  (whitespace-normalized) in its OWNING section file or the rendered pdf
  text, optionally at an exact count.
- **NEVER** — forbidden phrasing appears NOWHERE in rendered prose (comments
  stripped); optional context whitelist for scoped uses.
- **NEG** — pins an ABSENCE: the string must land exactly 0 times.
- **ABSTRACT** — sentinel phrases present in (or dropped from) the abstract
  block (tex `\begin{abstract}`, or a rendered-text span between two
  markers).

Prose is compared comment-stripped (a commented-out sentence can neither
satisfy nor trip a row), and rendered text is wrap-proof joined (form feeds
and page-number lines dropped, hyphen wraps rejoined — a pinned phrase that
straddles a page turn still counts as contiguous). Rendered text comes from
`pdftotext` (the `pdf` spec key; requires the `pdftotext` binary) or from a
pre-rendered `txt` file (no binary needed). Fail-closed freshness: with
`pdf_fresh_vs`, a render older than the tex it claims to render REFUSES
(rc 2) rather than passing against a stale dump. Exit codes: 0 all rows
pass / 1 any FAIL / 2 refusal or malformed spec.

`--selftest` runs a gated fixture battery — 8 gates: every row class
catches its planted failure and passes its clean twin, comment stripping
holds in both directions, the wrap-proof join and the stale-render refusal
each fire. It builds its fixtures in a temp dir and renders from `txt`, so
it needs no pdftotext; the test suite runs it under both `python3` and
`python3 -O`.

## Comment-swallow seam gate: `paper_seams.py`

```
python3 tools/emitall/paper_seams.py <master.tex> [--snapshot <dir>] [--marker <label> ...] [--allow <phrase> ...] [--quiet]
```

The battery's sibling on the document side. Where `battery.py` pins what the
prose says, `paper_seams.py` checks that an edit did not break a sentence
across a `%` comment. A **comment swallow** is a comment block dropped INSIDE
a sentence: the head of the sentence is stranded on the comment's last line
and the rendered text reads "...end. lowercase fragment...". The cure is a
move of the comment to its own line boundary, never a rewording. Four rows,
over the master's `\input`/`\include` roster (one level of nesting):

- **S1** — a live line that ends a sentence, then a comment block, then a
  live line that begins lowercase.
- **S1b** — a comment block whose last line ends on a capitalized bare word
  (no terminal punctuation) before a lowercase live line.
- **S2** — rendered `word. lowercase` fragments (`pdftotext -layout` of the
  PDF beside the master; abbreviations and the bibliography excluded) that
  the de-commented tex also carries adjacently; reading-order artifacts
  (dropped section glyphs, ligatures, table layout) are INFO, not findings.
- **S3** — with `--snapshot <pre-fold dir>`: a fold comment block (one that
  carries a `--marker <label>`; default label `FOLD`) that ends on a bare word,
  or carries text that was live before the fold, ahead of a lowercase live line.

Tabular rows (`&`-separated) and macro lines legitimately begin lowercase and
are skipped; `--allow <phrase>` names a legitimate lowercase sentence head (a
tool name such as `txburst, the backend`) and suppresses strikes on live lines
that begin with it. Every finding names its file and line. Exit 0 clean / 1
at least one tex-confirmed swallow / 2 usage. Without a PDF beside the master
(or without `pdftotext`) S2 is skipped and the INFO line says so; the tex rows
still run. Fixtures under `specs/seams/` — a planted S1 swallow, its cured
twin, and a head-rule case — are run by `tests/test_paper_seams.py` as
subprocesses from a foreign cwd, as a script and as
`python3 -m emitall.paper_seams`.

## Tests

```
python3 tests/test_all.py        # from tools/emitall/ ; pytest -q also works
```

66 tests (57 in `tests/test_all.py` + 9 in `tests/test_paper_seams.py`):
pointer/pattern planted truths,
expression-language battery with safety rejections (attribute access,
subscripts, lambdas, keywords, unknown names all refused), count-where with
coercion guards and $ref membership, printed-mode byte-compatibility pins
(including the −0.6845 → "−0.68" F1 class), outward-band planted truths +
negative controls (inward lo, inward hi, loose band, on-grid boundary,
0.1-granularity Fraction exactness, malformed inputs), the 3-of-6 planted
replay, engine end-to-end planted specs (match/mismatch/exact/missing/
stale/max-age/eval-error/note templates), CLI exit codes (0 clean / 1
findings / 2 malformed spec), the shipped example deployment run clean and
with planted corruptions, and the lint batteries:
interpolated-literal planted truths (typed runs caught, interpolated
receipt fields exempt, fully-interpolated spec clean),
unseeded-sampler planted truths (set-intersection and dict-view
emission feeds caught; sorted/seeded variants clean; seed-after-draw still
fires), sink/template planted truths (emit-args, dict keys, keyword args,
returns, augassign/extend, %-and-format placeholder exemptions, for-loop
sinks, reducer exemptions, non-emission strings untouched), pragma honored/
reason-mandatory, lint CLI exit codes + dir scan + parse-error loudness,
and zipjoin/contains/frac planted truths with loud negative controls, and
the paper-claims battery selftest (`battery.py --selftest` run via
subprocess, 8/8 gates green, under both `python3` and `python3 -O`), and
the seam-gate fixture battery (planted S1 swallow exits 1 naming file:line,
cured twin exits 0, head rule struck bare and suppressed by `--allow`, an
unrelated `--allow` phrase suppresses nothing, `--quiet` keeps the verdict
and exit code, script and `-m` invocations, usage rc 2 — all from a foreign
cwd). Suite green under `python3` and `python3 -O`; no bare asserts in
package code.

## Honest limits

- **The spec cannot certify completeness.** It checks every claim it
  contains; whether the claims cover every number the report quotes is the
  spec author's responsibility. Port the document's numbers exhaustively —
  a real deployment runs to 100+ claims for a single report.
- **Quotes are transcribed once.** The quoted strings in the spec are typed
  from the document being audited; a typo in the spec's quote produces a
  spurious finding (loud), never a spurious match — but a quote omitted is
  silence.
- **`printed` mode verifies at the printed register only.** It inherits
  Python float formatting (round-half-even at the last digit), exactly like
  the emitting scripts it audits; it does not certify the underlying value,
  only that the receipt reproduces the quote at the quote's precision.
- **Staleness is mtime ordering, not content hashing.** A regenerated-then-
  backdated receipt would evade it; the check targets the honest failure
  mode (sentence filed before its receipt), not adversaries.
- **`eval-error` conflates spec bugs with receipt schema drift** — the
  engine cannot tell them apart; both demand attention and both are loud.
- The expression language is deliberately small (see the table above); it is
  a verification language, not a computation language.

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).

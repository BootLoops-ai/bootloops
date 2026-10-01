# emitall — guide

Tool page: https://bootloops.ai/tools/emitall.html

KIND: package (engine.py, expr.py, run.py, lint.py, battery.py, paper_seams.py;
specs/ example claim spec + toy receipts + seam-gate tex fixtures, tests/). NOTE:
the sha-pin writer/verifier pairs Emitall standardizes are project-side artifacts
and do NOT ship here.

PURPOSE: Claims-integrity emit-all harness — a project writes ONE declarative
claims spec (JSON/YAML) mapping each quoted headline to (receipt, field expression,
comparison mode); the engine re-emits every value from stored receipts, compares at
the quote's own register, and exits nonzero on any finding.

USE-WHEN:
- Filing/checking a report section or table quoting numbers derived from stored
  receipts — catch transcription drift, regenerated-receipt staleness, un-universed
  counts (e.g. a 19/21 quote vs the canonical 16/18), printed-precision drift
  (−0.6845 quoted −0.69 vs printed −0.68).
- Enforcing the outward-rounding law on printed band edges (motivating case: 3-of-6
  inward dollar band edges).
- Linting emission scripts for typed-not-emitted digits, hash-order emission,
  unseeded sampling (lint.py).
- Gating a paper's tex tree for comment swallows before the built document goes
  out — a `%` comment block dropped inside a sentence, stranding the sentence
  head on the comment's last line (paper_seams.py).

NOT-FOR: Spec completeness is the author's burden — it checks only the claims it
contains, not whether every quoted number is covered. printed mode verifies the
printed register only. mtime ordering is NOT content hashing. Lint sinks are
name-pattern (`emit|headline|summar`) — an emission list called `lines` is
invisible; variable-held templates unresolved; seeding detected per file, not per
RNG stream. Derivation belongs in receipts, not spec expressions.

INVOKE: `python3 run.py SPEC.json [--out report.json] [--quiet]`; lint:
`python3 lint.py <script-or-dir> [--json out] [--quiet]`; paper claims battery:
`python3 battery.py SPEC.json [--quiet]` (selftest: `python3 battery.py
--selftest`); comment-swallow seam gate: `python3 tools/emitall/paper_seams.py
<master.tex> [--snapshot <pre-fold dir>] [--marker <label> ...] [--allow <phrase> ...]
[--quiet]` (also
`python3 -m emitall.paper_seams …` with tools/ on the import path); tests:
`python3 -m pytest tools/emitall -q` from the repo root (or `python3
tests/test_all.py` + `python3 tests/test_paper_seams.py` from tools/emitall/).

INPUTS: Claims spec JSON (YAML if PyYAML present): campaign/register/quote_source/
base, receipts {alias: {path, format(json|csv), newer_than, max_age_days}}, claims
[{section,label,quoted,expr | exprs+format | outward-band quoted_lo/hi+granularity+
scale+lo_expr/hi_expr, mode, lets, note}]. Expressions: whitelisted-AST (never eval;
no attributes/subscripts/lambdas) — f() JSON-pointer, vals() with {a|b|c}/*
patterns, cw() count-where with coercion guards + $ref cross-receipt membership,
where/split/join/jsondump/joinrows, zipjoin/contains/frac (exact Fraction,
never through a double).

OUTPUTS: Findings table to stdout + JSON report (machine-read UTC date stamp, never
typed). Exit 0 clean / 1 any finding (mismatch, outward-band violation, eval-error,
missing receipt, stale receipt) / 2 malformed spec. Lint: exit 0/1/2; findings
digit-literal, unordered-iteration, unseeded-sampling, pragma-missing-reason,
parse-error. Seam gate: one `FAIL <row> <file>:<line> …` line per finding, an
INFO line (S2 candidate/artifact counts or "S2 skipped", roster size), a verdict
line `paper_seams: CLEAN | N swallow(s) -- <master>`; exit 0 clean / 1 any
tex-confirmed swallow / 2 usage.

ENV: no knobs beyond the in-source pragma `# emitall-lint: allow-literal <reason>`
(reason MANDATORY, bare pragma is itself a finding). The `pdftotext` binary is
needed only for pdf-mode battery specs and for the seam gate's S2 row (skipped,
with an INFO line, when no PDF sits beside the master or the binary is absent);
`battery.py --selftest` renders from txt fixtures and needs no binary.

GATES: Comparison modes: printed (byte-compatible fmt_like semantics; optional
"decimals"), exact (str==str), outward-band (quoted_lo == floor(lo/g)·g, quoted_hi
== ceil(hi/g)·g, exact Fraction accept/reject; inward-edge vs outward-loose
grading). Paper claims battery rows: FROZEN (verbatim pin in the owning file
or rendered text, optional exact count), NEVER (forbidden phrasing, optional
context whitelist), NEG (exact-0 absence), ABSTRACT (sentinel present/dropped);
comment-stripped prose, wrap-proof rendered join, fail-closed stale-render
refusal (rc 2). Seam-gate rows (paper_seams.py): S1 sentence-ending live line ->
`%` block -> lowercase live line; S1b block ending on a capitalized bare word
before a lowercase live line; S2 rendered "word. lowercase" fragments (pdftotext
-layout) confirmed adjacent in the de-commented tex (reading-order artifacts =
INFO); S3 (--snapshot) fold block (labeled by --marker, default FOLD) ending on a bare word or carrying
pre-fold live text before a lowercase live line; tabular rows and macro lines
skipped; --allow <phrase> = a legitimate lowercase sentence head. Battery:
66 tests OK under python3 AND python3 -O, no skips — including the
shipped example deployment (specs/toy_stars.json over specs/toy_receipts/, an
invented program with synthetic numbers: 11 MATCH + 1 emitted-only clean;
planted transcription drift and nearest-rounded band edges fire as findings),
the battery selftest (8/8 gates: every row class catches its planted failure
and passes its clean twin), and the seam-gate fixtures (specs/seams/: planted
S1 swallow exits 1 naming file:line, cured twin exits 0, head-rule case struck
bare and suppressed by --allow; run as subprocesses from a foreign cwd).

PUBLIC-RUNNABLE: everything — the example spec's receipts ship beside it.
Project-specific specs are not included; write your own claims spec against
your own receipts (specs/toy_stars.json is the format example).

FOOTGUNS:
- A claim with no "quoted" is EMITTED-only, never a finding — don't mistake it for
  a passing check.
- Python booleans count as numeric in printed mode — wrap in str(...).
- Stale receipts still evaluate (staleness is the finding); missing receipts mark
  claims MISSING-RECEIPT.
- Lint flags quote-transcription args to emit-named functions even though they're
  comparison targets — fix: pragma with reason, or migrate quotes into a claims
  spec (spec JSON is data, not linted).
- cw(): guard `not_in ["NA",""]` BEFORE an int/float coercion (conditions ANDed in
  order, short-circuit).
- paper_seams.py walks ONE level of `\input` nesting below the master; a section
  file that `\input`s a file that `\input`s another is not reached. S2 needs the
  built PDF beside the master (same stem) — with no PDF the run can only be as
  clean as the tex rows; the INFO line says S2 was skipped.

RELATED: intended audience = any project filing quoted results.
RELATED PRIOR ART: reproducible-manuscript pipelines that regenerate quoted numbers from
source (e.g. showyourwork; Luger, Bedell, Foreman-Mackey, Crossfield, Zhao & Hogg 2021,
arXiv:2110.06271) address the same rot from the build side; emitall works from the
receipt side. PDF text via poppler's pdftotext.

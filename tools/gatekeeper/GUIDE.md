# Gatekeeper — guide

Gatekeeper — `tools/gatekeeper` — https://bootloops.ai/tools/gatekeeper.html

KIND: package (stdlib + mpmath only)

PURPOSE: data-integrity + held-out-CV honesty layer for amflow/oracle
per-master-per-point JSON dumps, canonical AMFlow result parsing, and pre-farm
quadrature probes. Gates every closure.

USE-WHEN:
- Suspect duplicated/oracle-contaminated input → hash-dedup scan (the trap:
  byte-identical numerics under two master labels ⇒ the fit
  "certifies" its own training data).
- Certify a value-fit → `heldout_certify` (leave-one-out CV + PSLQ rational fit).
- Parse AMFlow `solve_integrals` JSON → `amf_laurent`/`amf_arb`. Always parse via this,
  never inline.
- Near-singular fit matrix / PSLQ won't close → `auto_condition` + `pslq_with_lll`
  (LLL pre-reduction).
- About to farm a multi-level quadrature / farm digits cap despite level bumps →
  `quad_probe` (3-min probe BEFORE any farm).
- 1D chain `prod P_i^(-1/2)` between two roots — which rule? → `ell_subst`.
- "Do towers exist at point X for masters Y?" / duplicate-oracle hazard →
  `gatekeeper.registry` (env-configured roots; flat `oracle_registry.py` = shim).
- Held-out digits look dps/boundary-limited — is it really boundary precision? →
  `ctrl_per_order.py`.

NOT-FOR: producing oracle values (that's AMFlow/evaluators); ell_subst outside its three
supported root configurations returns honest-unsupported by design (interior simple
root ⇒ the integral is not real as posed; endpoint multiplicity >1 ⇒ divergent; ≥2
pinching roots on one side ⇒ no rule built).

INVOKE:
- `python -m gatekeeper.oracle_integrity DIR... --report integrity_report.json
  [--quarantine] [--strict-zero-balls]` — exit 0 CLEAN / 2 QUARANTINE_NEEDED. Drops
  identity fields, SHA256 of normalized payload; flags duplicates, empty/all-zero coeffs,
  unreadable. Values are read through `amf_ball` (plain decimals, Arb balls
  `[mid +/- rad]` / `mid +/- rad` / `[+/- rad]`, `[a, b]` intervals); the zero test is on
  the MIDPOINT (mpmath, dps from the string length). A ball with rad > |mid| is
  `consistent_with_zero` — its own report field, REPORT-ONLY by default (verdict and exit
  code unchanged); `--strict-zero-balls` promotes those files to quarantine.
- Two dump forms are read: the tagged per-file form (top-level `tag`/`indices` beside
  `coeffs` = `{"<order>": {"re", "im"}}`, one record per file) and the AMFlow
  `solve_integrals` out.json form (`result[i].coefficients` = `[{order, value: {re, im}}]`
  with `result[i].integral.indices`): every result entry is one record tagged
  `<file-tag>#<indices>` (file-tag = the file's `tag`, else `integral.family`, else the
  basename stem), coefficients re-keyed `{str(order): value}` for the classifier and the
  duplicate hash; a file with neither key, or an empty `result` list, reads as
  `empty coeffs`; `n_files` counts files, `n_records` records; `--quarantine` stays
  file-granular (a flagged entry moves its whole file).
- `python -m gatekeeper.heldout_cv --oracle DIRS --basis basis_values.json
  --weight W --gate 30 --integrity integrity_report.json --report cv_report.json`.
- API: `scan_integrity`, `heldout_certify`, `load_oracle_values`, `lyndon_basis`,
  `auto_condition`, `pslq_with_lll`, `amf_laurent`, `amf_arb`, `amf_ball` (the one
  ball-string splitter: `"[X +/- E]"` → `("X","E")` as decimal strings).
- Probes: `python3 tools/gatekeeper/probes/quad_probe.py` (per-chain convergence:
  vary ONE chain's level, others frozen so their error cancels exactly in consecutive
  diffs; low dps 25-30, 1-2 points, parallel subprocess configs; env-var levels
  LEVEL_OUT/LEVEL_K pattern). `python3 tools/gatekeeper/probes/ell_subst.py`
  (2 endpoint roots → Gauss-Chebyshev, affine τ; 2 endpoint + 2 outside-pinching one
  per side → Jacobi τ=sn(u,m), Möbius-symmetrized outer pair ±σ, m=1/σ² from the
  cross-ratio, trapezoid-on-period rate exp(−πNK′/K); 2 endpoint + ONE
  outside-pinching root on a single side → τ=sn²(u,m) (the G-R 3.131 map,
  m=(b−a)/(o−a) mirrored per side, all three sqrt-roots absorbed exactly,
  midpoint-on-(0,K) = 2N-point trapezoid on the period, rate exp(−2πNK′/K);
  no `log_frame()` for this kind — it raises); smooth factors folded into
  weights cancellation-free; `log_frame()` = cheb_log_table hook for ε-tower log
  moments). Acceptance: `python3 tools/gatekeeper/probes/ell_subst_test.py` (T4
  cross-check leg needs a module not shipped in this repo; it SKIPs loudly).
- `python3 tools/gatekeeper/ctrl_per_order.py --selftest` (fixtures are not shipped in this repo —
  it SKIPs loudly without them; CTRL_TOWERS/CTRL_ORACLE point at your
  own pair). `python -m gatekeeper.registry --point -2,-3 --masters 13..40` /
  `--dupes`.

INPUTS: per-(master, point) JSON dumps; AMFlow solve_integrals result entries;
basis_values.json for CV.

OUTPUTS: integrity_report.json (duplicates/suspicious_zero/unreadable/consistent_with_zero
[{file,tag,order,part,mid,rad}] + strict_zero_balls, verdict CLEAN|QUARANTINE_NEEDED;
`--quarantine` moves offenders to `<dir>_QUARANTINE/`);
cv_report.json with per-master per-order certification; amf_laurent → ε-Laurent dict;
oracle/REGISTRY.json.

ENV (registry only): `ORACLE_REG_ROOTS` (colon-separated root globs),
`ORACLE_REG_REGISTRY` (registry json path), `ORACLE_REG_COMBOS` (combos dir, default =
dirname(registry)) — points the engine at your own oracle tree.

GATES (worst-held-out-point certification discipline): a master is CERTIFIED only if
EVERY held-out point matches the N−1-point fit to ≥ gate digits — the verdict is the
WORST held-out point, never the average. Duplicates from the integrity report are
excluded BEFORE fitting (`n_duplicates_excluded` reported). Run the integrity scan
before the CV; a fit certified on unscanned data is void.

FOOTGUNS:
- `amf_laurent`: `coefficients[].order` is the ABSOLUTE ε-power; `leading_order` is
  redundant (= min order). The `{lo+c['order']: ...}` pattern DOUBLE-COUNTS and is
  silent when lo==0 (most 2L) — this class can read 1d as 139d. Parser warns if
  leading_order≠min(order).
- AMFlow writes files even when it returned nothing (empty/all-zero coeffs) —
  integrity scan flags them; never consume unscanned dumps.
- Never `float()` an Arb ball string: `"[0.25… +/- 1e-115]"` fails `float()`, and a
  strip-the-sign fallback sees `[` and reads a valid nonzero value as ZERO (the false
  QUARANTINE_NEEDED class). Parse through `amf_ball`/`amf_arb`; the integrity scan does.
  The printed radius is goal-bounded (tracks working_pre), so `consistent_with_zero`
  is a precision statement, not a mis-write — hence report-only by default.
- quad_probe: <5 digits/level ⇒ METHOD-MISMATCHED flag (switch to classical-weight
  Gauss). tanh-sinh is superexponential — probe where consecutive diffs ≲1e-8 before
  trusting a ts flag.
- ell_subst: never silent tanh-sinh on unsupported configurations — honest-unsupported
  is the correct output.
- ctrl_per_order exists because a scalar min-over-orders "floor" misdiagnoses:
  per-(master, ε-order) table + digit-loss slope vs order compared to |log₁₀ε₀| gives
  an ALIAS-LIMITED verdict + usable-order cutoff (≥30d); a jet truncation-alias ladder
  (~8.5 d/order) masquerades as a boundary-precision problem. Skips |oracle|<1e-30,
  caps 99d.

TESTS: `python3 gatekeeper/test_gatekeeper.py` (T1 integrity duplicate trap,
LLL conditioning recovery, pslq_with_lll, T6-T9 ball forms: verbatim AMFlow ball fixture
`fixtures/gatekeeper/amflow_ball_nonzero.json` not flagged, six true-zero forms flagged,
consistent_with_zero status + strict promotion, CLI end-to-end exit codes — ALL TESTS
PASS, self-contained; this is the registered battery together with ctrl_per_order);
`python3 gatekeeper/probes/ell_subst_test.py` (acceptance ≥40d gates,
self-contained but T4).

RELATED: tools/lockpick (PSLQ engines + discipline); tools/nestor (nestor.probe IS
quad_probe by identity). Run quad_probe before ANY multi-level farm.

CREDIT: PSLQ after Ferguson, Bailey and Arno [PSLQ] via mpmath (F. Johansson and
contributors), LLL pre-reduction after Lenstra–Lenstra–Lovász [LLL]; the AMFlow result
format parsed here is Liu and Ma's [AMFlow] as emitted by AMFlow.cpp [AMFcpp]. Bracketed
keys resolve in REFERENCES.md at the repository root.

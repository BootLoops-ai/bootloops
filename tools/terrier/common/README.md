# common/ — shared chassis (both wings import; never re-grow these)

Small by design: only patterns BOTH wings demonstrably use.
All batteries are `selftest_*.py` next to each module (the root selftest.py
runs them); none needs external data.

| module | what it owns | battery |
|---|---|---|
| receipts.py | receipt schema (schema/status/pins.code_sha256/reason), atomic tmp+fsync+rename, sha-stream, code-sha stamp, keyed jsonl merge, fsync'd append-log, receipt-as-checkpoint resume | selftest_receipts.py 12/12 |
| controls.py | CTRL-POS/ID-V1 hash placement, blind injection order, VOID_CONTROL_FAIL rule, mutation-must-fail harness + battery, planted-error drill (loudly-synthetic law) | selftest_controls.py 14/14 |
| frames.py | convention-registry template: rows w/ status vocab (PINNED-V/-N/-V+N, OURS, PIN-REQ), evidence classes [SRC]/[PIL]/[REF], executed-toggle tables, tamper-evident freeze/load, require_pinned fail-closed gate, markdown dump | selftest_frames.py 9/9 |
| certs.py | non-asserting dict gates (two-dps, ball, route), ball-honesty radius-reporting law, published-value gate (last-digit-ulp band), float64 trim guard (subnormal pitfall), is_ball_str serialization law | selftest_certs.py 22/22 |
| verdict.py | asserting print-scoreboard gate primitives (the wings' verdicts facades delegate here), parse_ball / matched_digits, sha-stamped atomic receipt | selftest_verdict.py |
| dedup.py | known-point pullback-orbit dedup + newform backstop over the packaged data/KNOWN_ORBIT.jsonl; explicit-receipt law (no ledger write without a path) | selftest_dedup.py 5/5 |

## Delegation map (DELEGATION LAW: call, never duplicate)
- certs.py imports parse_ball + matched_digits FROM verdict.py (single ball
  primitive). Split of surfaces: verdict.py asserts and prints (scoreboard
  semantics); certs.py returns {"pass": ...} dicts for receipts and
  expected-fail batteries.
- frames.py reuses receipts.atomic_write_json/utc_now.
- code-sha: receipts.code_sha(files, root) is the chassis stamp;
  verdict.code_sha(*paths) delegates to it (sorted concatenation; do not add
  a third).

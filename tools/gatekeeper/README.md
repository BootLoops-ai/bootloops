# gatekeeper

Two cheap honesty wins for amflow/oracle per-master-per-point JSON dumps,
packaged as one reusable tool. **stdlib + mpmath only.**

## Why this exists

An oracle farm writes one JSON per `(master, point)`. Collect two masters with
the wrong propagator indices and their files come out **byte-identical** to two
others — the fit then "certifies" the duplicates, recognising *its own training
data under a second label*, and the certified-master count silently inflates.

A 30-second hash scan catches it. This package is that scan, plus
the held-out-CV gate that makes a value-fit honest.

## 1. `oracle_integrity.py` — dedup auditor (FAST, hash-only)

```bash
python -m gatekeeper.oracle_integrity mdata/ mdata_nt/ \
       --report integrity_report.json [--quarantine]
```

What it does
* Walk all `*.json`, parse, **drop the identity fields** (`tag`/`k`/`indices`),
  re-serialize sorted, SHA256.
* Any hash group with **>1 distinct (tag, indices)** ⇒ DUPLICATE — same
  numerics under two master labels.
* Also flags **empty / all-zero `coeffs`** (the `eps_order` footgun: amflow
  returned nothing, file still written) and unreadable files.
* Two dump forms are read: the tagged per-file form (top-level `tag`/`indices`
  beside `coeffs` = `{"<order>": {"re", "im"}}`, one record per file) and the
  AMFlow `solve_integrals` out.json form (`result[i].coefficients` =
  `[{order, value: {re, im}}]` with `result[i].integral.indices`): every
  result entry is one record tagged `<file-tag>#<indices>` (file-tag = the
  file's `tag`, else `integral.family`, else the basename stem), its
  coefficients re-keyed `{str(order): value}` for the classifier and for the
  duplicate hash. A file with neither key, or an empty `result` list, reads
  as `empty coeffs`. `n_files` counts files, `n_records` records;
  `--quarantine` stays file-granular (a flagged entry moves its whole file).
* Coefficient strings are read through `amf_result.amf_ball`, so every form
  AMFlow's Arb printer emits parses: plain decimals, `[mid +/- rad]` balls
  (also `mid +/- rad` and the zero-midpoint `[+/- rad]`), and `[a, b]`
  intervals. The zero test is applied to the **midpoint** with mpmath at a
  precision set from the string length — never `float()`.
* A ball whose **radius exceeds |mid|** is **consistent with zero**: listed
  in its own field `consistent_with_zero`, never folded into
  `suspicious_zero`. Report-only by default (verdict and exit code
  unchanged); `--strict-zero-balls` promotes those files to quarantine.
* `--quarantine` moves offending files to `<dir>_QUARANTINE/` (keeps the first
  per hash).

Output `integrity_report.json`
```json
{"dirs":[...], "n_files":N, "n_records":N, "ok_count":N,
 "duplicates":[{"hash":"...","files":[...],"tags":[...],"indices":[...],"n":k}],
 "suspicious_zero":[{"file":"...","tag":"...","reason":"..."}],
 "unreadable":[...],
 "consistent_with_zero":[{"file":"...","tag":"...","order":"-1","part":"re",
                          "mid":"1.0e-40","rad":"1.0e-30"}],
 "strict_zero_balls": false,
 "verdict":"CLEAN" | "QUARANTINE_NEEDED"}
```

Exit code: `0` if CLEAN, `2` if QUARANTINE_NEEDED (`consistent_with_zero`
entries change the exit code only under `--strict-zero-balls`).

## 2. `heldout_cv.py` — leave-one-out CV + PSLQ certifier

```bash
python -m gatekeeper.heldout_cv \
       --oracle mdata/ mdata_nt/ \
       --basis basis_values.json --weight 2 --gate 30 \
       --integrity integrity_report.json --report cv_report.json
```

For each master `m` at ε-order `w`:
1. collect values at all N points (duplicates from the integrity report are
   **excluded before fitting** — `n_duplicates_excluded` is reported);
2. for each held-out point `P_k`: least-squares fit on the other N−1, predict
   at `P_k`, record `-log10|pred − actual|`;
3. **CERTIFIED** ⇔ `min_k heldout_digits ≥ gate` (default 30);
4. only if certified: `mpmath.pslq` each fitted coefficient to a rational.

Honesty fields per `(master, weight)`: `n_points_used`,
`n_duplicates_excluded`, `basis_dimension`, `condition_number`,
`heldout_digits[]`, `min`, `mean`, `certified`, `pslq_rationals|null`.

`basis_values.json` schema:
```json
{"names":["B0","B1",...],
 "points":{"-3|-5":{"B0":"<mpf>","B1":"<mpf>",...}, ...}}
```
or per-weight under a top-level `"weights":{"2":{...}}`.

## Python API

```python
from gatekeeper import scan_integrity, heldout_certify, load_oracle_values

rep = scan_integrity(["mdata/", "mdata_nt/"])
excl = rep.quarantine_set()                       # files to drop from a fit

vals, kin, _ = load_oracle_values(["mdata/"], exclude_files=excl)
res = heldout_certify(vals, kin, weight=2,
                      basis_fn=lambda pk,s,t: [...], basis_names=[...],
                      gate=30, dps=80,
                      n_duplicates_excluded=len(excl))
```

## Tests

```bash
python tools/gatekeeper/test_gatekeeper.py
# T1 integrity-dedup: PASS
# T2 heldout-CV:      PASS  (min≈90d, rats=['1','3'])
# ...
# T6 ball-form nonzero NOT flagged: PASS   (verbatim AMFlow balls, fixtures/gatekeeper/)
# T7 true zeros MUST flag: PASS
# T8 consistent_with_zero: PASS
# T9 CLI end-to-end: PASS
```

## 3. `probes/` — Gatekeeper's probe kit

Pre-farm quadrature honesty. The flat shims `tools/ell_subst.py` and
`tools/quad_probe.py` forward here (`python3 tools/quad_probe.py ...` runs the
CLI).

* `probes/ell_subst.py` — weight-matched 1D substitution selector for
  sqrt-product chains: Gauss-Chebyshev / Jacobi-sn via cross-ratio /
  one-sided Jacobi-sn² (G-R 3.131 map for a single outside-pinching root,
  either side), honest `unsupported` instead of the silent tanh-sinh
  ~3 d/level trap. Needs **sympy** on top of the package's stdlib+mpmath
  baseline.
* `probes/quad_probe.py` — per-chain convergence prober: vary ONE chain's
  level (others cancel exactly in consecutive diffs); <5 digits/level =
  METHOD-MISMATCHED. Probe BEFORE farming.

```bash
python -m gatekeeper.probes.quad_probe --workdir <dir> --module <mod> \
    --func <fn> --base DPS=26,LEVEL_K=4 --chains LEVEL_K --span 2
python tools/gatekeeper/probes/ell_subst_test.py   # acceptance battery (T1–T8)
```

```python
from gatekeeper.probes import select, cubic_under_sqrt_K, probe, EnvEvaluator
```

## Dependencies

Python standard library and `mpmath`; `sympy` for `probes/ell_subst.py`.
`heldout_cv.py` uses `fpylll` (GPL-2.0-or-later) for LLL when it is importable
and falls back to a pure-Python reduction otherwise; `fpylll` is optional and is
not bundled.

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).

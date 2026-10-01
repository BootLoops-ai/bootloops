# fibrate — sound fibration engine (subtropica gap-tool)

Fibrates ZIP_reg words with ze-dependent letters into
`G({letters}; ze) × level-2 constants`, replacing HyperFLINT's
**unsound** `fibration_basis` op (which silently returned wrong values on
631/635 reg-divergent split jobs and corrupted on-contour words in a
large production run). Upstream's own fibration route
(`reference/SubTropica.wl:11411`, `:15565 STToFibrationBasis`) shells to
HyperInt's `fibrationBasis` with no per-word validation — same trust gap.

Validated in a production pilot: 370/370 words fibrated + validated;
the assembled formula passed the 38–42-digit held-out oracle gate.

## Method

Goncharov total differential on the regularized iterated integral
`I(0; b1..bn; ∞)` (`b` = reversed ZIP word): differentiation peels
letters via dlog one-forms (paper eqs 3.40/3.41; parameter-dependent
letters per the total-differential remark, paper:1423–1426), followed by
recursive **exact** sympy integration in `ze` (antiderivative
construction = paper §3.3.2 step 2), an integration constant fixed
numerically at `ze = 1/2` and identified by PSLQ over the level-2
weight ≤ 5 constant basis, and **mandatory** per-word numeric validation
against the trusted step-1 ZIP semantics.

## Hardening (vs the pilot prototype)

1. **Adaptive ZIP series lengths** (`zip_num.py`): per-word, keyed to the
   minimum mapped Hoelder-letter magnitude. A fixed length would be a
   dps-independent accuracy ceiling;
   the constructor accepts **no** fixed-length override, and
   `--selfcheck` verifies the invariant against the exact identity
   `ZIP_reg([-9/26]) = log(26/9)` (≥ 120 of 130 digits required).
2. **PSLQ control battery** on the integration-constant step
   (`Engine.identify_const`): dual-precision stability (`dps` and
   `dps−40`), height cap `maxcoeff = 10^10`, residual gate
   `10^-(dps−30)`, and a **positive control** (a known 4-term rational
   combination must be recovered exactly) that runs at every `Engine`
   construction — the engine refuses to start if it fails. A negative
   control (`e` must NOT be identified) runs in `--selfcheck`.
3. **Provenance stamping**: every word whose representation
   (transitively, through memoized subword fibrations) contains a
   PSLQ-identified nonzero constant carries `"pslq_constants"` in its
   `provenance` list, and the top-level `provenance` carries it if any
   word does. Downstream `LaurentSeries` builders MUST consume this and
   not claim `:analytic` for such coefficients.
   Zero constants are established by the numeric-zero gate
   `|C| < 10^-(dps−25)` plus the validation below (never by fiat).
4. **Mandatory 2-point validation**: every fibrated word (including all
   internal subwords) is evaluated at ≥ 2 rational points (default
   `2/5`, `3/8`; all distinct from the fit point `1/2`) and compared to
   the trusted step-1 ZIP semantics (`zip_num.ZipEval` — verified correct
   in all 869 step-1 calls of the pilot) at relative
   tolerance `10^-(dps−40)`. There is deliberately **no flag to disable
   this**; supplying fewer than 2 distinct points is a typed refusal.

Failure anywhere (PSLQ, validation, on-contour letters, algebraic poles)
is a loud `FibrateError` → CLI exit code 2, **no partial output**.

## CLI

```bash
ulimit -v 32505856   # standing memory cap, every invocation
python3 scripts/fibrate/fibrate.py WORDLIST.json [--out OUT.json] \
    [--dps 130] [--validation-points 2/5 3/8]
python3 scripts/fibrate/fibrate.py --selfcheck
```

- `WORDLIST.json`: list of words; each word is a list of letter strings
  in `ze`, e.g. `[["-1", "-ze", "1/(ze - 1)", "0"], ...]` (the pilot's
  `unres_words.json` format). Letters must evaluate ≤ 0 at the fit and
  validation points (on-contour letters are refused — that is exactly
  the class where the upstream op is unsound).
- `--dps` ≥ 60 (default 130). Constants are built at the requested
  precision (never import-time ambient — mplll-ambient-dps footgun
  class).
- Exit codes: `0` success, `2` any refusal/failure.

### Output document (JSON)

```
{ "tool": "subtropica/fibrate", "version": "1.0.0",
  "dps": 130, "fit_point": "1/2",
  "validation": {"points": ["2/5","3/8"], "rel_tol": "1e-90",
                  "mandatory": true, "oracle": "zip_num.ZipEval ..."},
  "pslq_controls": {"positive_control": "pass",
                     "dual_precision": [130, 90],
                     "height_cap": 10000000000,
                     "residual_gate": "1e-100", "zero_gate": "1e-105"},
  "provenance": ["goncharov_differential","validated_2pt",
                  "pslq_constants"?],
  "words": [ {"word": [...], "weight": n,
              "terms": [{"gword": ["0","1",...],   # target G letters
                          "consts": ["z3",...],     # constant multiset
                          "coef": "<sympy ratfn in ze>"}, ...],
              "used_pslq": bool, "provenance": [...],
              "validated_digits": {"2/5": 129, "3/8": 129}}, ... ] }
```

A word's value is
`ZIP_reg(word)(ze) = Σ_terms coef(ze) · G(gword; ze) · Π_{c∈consts} c`,
with `G` in the convention of eq (3.40) (trailing-zero shuffle strip,
`G(0^k;x) = log^k(x)/k!`) and the constant names from
`fibrate.build_consts` (`z2, z3, z4, z5, l2, li4, ..., z2z3 = ζ2·ζ3`).
Terms are canonically ordered (`fibrate.serialize_rep`).

## Python API

```python
import sys; sys.path.insert(0, ".../tools/subtropica/scripts")
from fibrate import Engine, run_wordlist, selfcheck, parse_letter

doc = run_wordlist([["-1", "-ze"]], dps=130)          # full document
eng = Engine(dps=130)                                  # positive control runs here
rep = eng.fib(tuple(parse_letter(x) for x in ["-1", "-ze"]))
# rep: {(gword tuple of sympy numbers, const-name tuple): sympy ratfn}
eng.meta[...]   # {'used_pslq': bool, 'digits': {...}} per word
```

`Engine` construction sets the ambient mpmath precision (single-purpose
process model). Both `import fibrate` (this dir on `sys.path`) and
`from fibrate import ...` (parent `scripts/` on `sys.path`, package
style) work.

## Regression corpus / tests

`fixtures/regression_words.json` + `fixtures/regression_expected.json`:
10 words (weights 1–5, with and without PSLQ constants) whose expected
canonical terms are **copied from the pilot's validated memo**
`engine_memo.pkl` (regenerate with `fixtures/make_regression.py`, which
needs `SUBTROPICA_PILOT_DIR` pointing at the pilot dir — not shipped).

Julia test (standalone, default global env, memory cap inside):

```bash
ulimit -v 32505856; julia <checkout>/test/test_fibrate.jl
```

covers: selfcheck battery, exact fixture equality on the 10-word corpus,
≥ 85 validated digits per word per point, provenance stamps, a
comparator negative control (perturbed terms must mismatch —
print-loud-rc0 class), and typed refusals (on-contour letter, single
validation point, malformed JSON). Last run: 94/94 pass, ~2.5 min wall.

## Scope / limits

- Letters: rational functions of `ze` that are ≤ 0 at all evaluation
  points; poles arising in integration must be RATIONAL (linear
  denominator factors) — algebraic poles are a typed refusal.
- Constants: level-2 (alternating) basis through weight 5. Words of
  weight > 5 would need basis extension (`build_consts`).
- The engine memoizes per-process; batch words per Engine/process. It is
  serial by design — parallelism belongs at the farm layer (one process
  per worker, as with the C11 farm).

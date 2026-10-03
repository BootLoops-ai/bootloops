# cohortgate - guide

KIND: package (stdlib only; no third-party dependencies)

PURPOSE: a cohort and **group-label** integrity gate. A group label must be
derivable from a *declared rule*; it may never be inferred from a file name, a
sample-name prefix, or a hand-typed column. The row count is additionally
re-derived by a second, independent route and the two routes must agree.

USE-WHEN:
- A cohort table with a group / phenotype column is about to be compared across
  groups -> run this gate first (the sourcing and label steps of a
  data-integrity gate).
- You suspect the grouping was produced from file names or name prefixes. A
  prefix reflects a naming habit, not a group definition, and the two can be
  completely decoupled: a batch of files can all carry the same prefix while
  every subject belongs to one group.
- The `group` column was hand-edited, or several scripts each pasted their own
  labels -> `label_mismatch` (column vs rule) names every affected row.
- Counts do not add up and you cannot find why -> `count_invariant` (n_rows vs
  the sum of label counts, two independent routes) together with `unclassified`
  and `ambiguous` points at the row.
- A pipeline needs a negative control that can actually fire -> this package's
  battery is the pattern.

NOT-FOR:
- No statistical modelling, no adjustment, no guessing of missing-value
  semantics. Missing = matches no rule -> reported as `unclassified`, never
  imputed with a default.
- It does not judge whether the grouping rule is biologically sensible - the
  rule is declared in the spec and the gate only checks that table and rule are
  self-consistent.
- It does not do literature or metadata provenance (that is the sourcing step).

INVOKE:
- `python3 cohortgate.py SPEC.json DATA.csv [--quiet]`
  Exit codes: 0 CLEAN / 2 FINDINGS / 3 SPEC_OR_IO (missing spec fields, an empty
  rule set, unreadable files - it refuses to guess).
- Battery: `python3 selftest.py` (8 legs, ~0.1 s, no external dependencies)

SPEC shape:
```json
{ "id_column": "subject_id", "declared_n": 6,
  "group_rule": { "check_declared_label": true,
                  "rules": { "balanced":    [["lactate", "<=", 2]],
                             "glycolytic":  [["lactate", ">", 2], ["bili", "<=", 1.2]],
                             "hepatic":     [["lactate", ">", 2], ["bili", ">", 1.2]] } } }
```
A condition is `[column, operator, value]` with the operator drawn from
`< <= > >= == !=`; all conditions of a rule must hold for it to match; a
non-numeric cell or a missing column never matches.

ACCEPTANCE GATES (this package's own standard):
- A clean fixture must read CLEAN, with `n_rows == sum(label counts)` (the
  second, independent route).
- Three planted defects (label mismatch / duplicate id / a row outside every
  rule) must each be caught, with the finding naming file:row and id.
- A wrong `declared_n` must be caught (fail-closed).
- An empty rule set must be **refused** (rc 3), never silently passed.
- With `check_declared_label=false` the label check must switch off - proving
  the check is spec-driven and not hard-coded.
- Closing law: the clean fixture passes and all three planted fixtures trip.
  **A gate that has never failed anything certifies nothing.**

VERIFICATION CLASS: selftest (the battery runs green from a cold clone: no
external engine, no gated data)

LIMITS (stated, not hidden):
- CSV with numeric conditions only; categorical conditions, non-interval logic
  and multi-table joins are not implemented.
- Overlapping rules are reported as `ambiguous`; the gate does not resolve the
  overlap - that is a modelling decision.
- The shipped fixtures are synthetic and exist only for the battery to
  self-certify. Before use on a real cohort, register the data through the
  project's usual intake.

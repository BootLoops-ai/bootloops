#!/usr/bin/env python3
"""cohortgate - a cohort and group-label integrity gate (stdlib only).

LAW: a group label must be derivable from a declared rule. It may never be
inferred from a file name, a sample-name prefix, or a hand-typed column: a
prefix reflects a naming habit, not a group definition, and the two can be
completely decoupled.

Six findings are possible; any one of them is a finding:
  1. label_mismatch      - the hand-typed `group` disagrees with the rule-derived label
  2. unclassified        - the row matches no group rule
  3. ambiguous           - the row matches more than one group rule
  4. count_invariant     - n_rows != sum(label counts) (the second, independent route)
  5. declared_n_mismatch - n_rows disagrees with the cohort size declared in the spec (fail-closed)
  6. duplicate_id        - the same id appears more than once

Usage:
    python3 cohortgate.py SPEC.json DATA.csv [--quiet]
Exit codes:
    0 = CLEAN      - no finding
    2 = FINDINGS   - findings present (each printed as kind / file:row / id / note)
    3 = SPEC_OR_IO - spec or data unreadable, or spec missing fields (refuses to guess)
"""
import csv
import json
import sys

OPS = {
    "<=": lambda a, b: a <= b,
    ">=": lambda a, b: a >= b,
    "<": lambda a, b: a < b,
    ">": lambda a, b: a > b,
    "==": lambda a, b: a == b,
    "!=": lambda a, b: a != b,
}


class SpecError(Exception):
    pass


def load_spec(path):
    try:
        with open(path, encoding="utf-8") as fh:
            spec = json.load(fh)
    except (OSError, ValueError) as e:
        raise SpecError(f"cannot read spec: {e}") from None
    if "group_rule" not in spec or "rules" not in spec["group_rule"]:
        raise SpecError("spec has no group_rule.rules - refusing to guess the grouping")
    if not spec["group_rule"]["rules"]:
        raise SpecError("spec group_rule.rules is empty - an empty rule set would pass every row")
    return spec


def load_rows(path):
    try:
        with open(path, encoding="utf-8", newline="") as fh:
            return list(csv.DictReader(fh))
    except OSError as e:
        raise SpecError(f"cannot read data: {e}") from None


def _num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def _match(row, conds):
    """conds: [[col, op, value], ...]; all conditions must hold.

    A non-numeric cell or a missing column never matches.
    """
    for c in conds:
        if len(c) != 3:
            raise SpecError(f"a condition must be [col, op, value], got {c!r}")
        col, op, val = c
        if op not in OPS:
            raise SpecError(f"unknown comparison operator {op!r}")
        x = _num(row.get(col))
        if x is None:
            return False
        if not OPS[op](x, float(val)):
            return False
    return True


def derive_label(row, rules):
    """Return (label, None), or (None, 'unclassified'|'ambiguous')."""
    hits = [lab for lab, conds in rules.items() if _match(row, conds)]
    if len(hits) == 1:
        return hits[0], None
    return None, ("unclassified" if not hits else "ambiguous")


def audit(rows, spec):
    """Return (findings, counts, n_rows)."""
    findings = []
    idc = spec.get("id_column", "subject_id")
    rules = spec["group_rule"]["rules"]
    check_lbl = spec["group_rule"].get("check_declared_label", True)

    counts = {lab: 0 for lab in rules}
    seen = {}
    for i, row in enumerate(rows, start=2):          # line numbers include the header row
        sid = row.get(idc)
        if sid in seen:
            findings.append({"kind": "duplicate_id", "row": i, "id": sid,
                             "note": f"id repeated (first seen on line {seen[sid]})"})
        else:
            seen[sid] = i

        lab, why = derive_label(row, rules)
        if why:
            findings.append({"kind": why, "row": i, "id": sid,
                             "note": "row matches no group rule" if why == "unclassified"
                                     else "row matches more than one group rule"})
            continue
        counts[lab] += 1
        if check_lbl:
            declared = row.get("group")
            if declared != lab:
                findings.append({"kind": "label_mismatch", "row": i, "id": sid,
                                 "note": f"column group={declared!r} but the rule derives {lab!r}"})

    n_rows = len(rows)
    n_sum = sum(counts.values())
    if n_rows != n_sum:
        findings.append({"kind": "count_invariant", "row": None, "id": None,
                         "note": f"n_rows={n_rows} != sum(label counts)={n_sum} "
                                 f"(the two independent routes disagree)"})

    declared = spec.get("declared_n")
    if declared is not None and int(declared) != n_rows:
        findings.append({"kind": "declared_n_mismatch", "row": None, "id": None,
                         "note": f"declared_n={declared} != n_rows={n_rows} (fail-closed)"})

    return findings, counts, n_rows


def main(argv):
    quiet = "--quiet" in argv
    args = [a for a in argv if a != "--quiet"]
    if len(args) != 2:
        print("usage: python3 cohortgate.py SPEC.json DATA.csv [--quiet]")
        return 3
    spec_path, data_path = args
    try:
        spec = load_spec(spec_path)
        rows = load_rows(data_path)
        findings, counts, n_rows = audit(rows, spec)
    except SpecError as e:
        print(f"SPEC-ERROR: {e}")
        return 3

    if not quiet:
        print(f"data: {data_path}  n_rows={n_rows}")
        print("counts: " + ", ".join(f"{k}={v}" for k, v in sorted(counts.items())))
        print(f"invariant: n_rows={n_rows} == sum(labels)={sum(counts.values())}")
        for f in findings:
            loc = f"{data_path}:{f['row']}" if f["row"] else data_path
            print(f"  FINDING {f['kind']:22s} {loc}  id={f['id']}  {f['note']}")
    if findings:
        if not quiet:
            print(f"verdict: FINDINGS ({len(findings)})")
        return 2
    if not quiet:
        print("verdict: CLEAN")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

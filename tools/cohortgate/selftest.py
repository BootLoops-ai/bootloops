#!/usr/bin/env python3
"""cohortgate acceptance battery.

Eight legs: a clean fixture passes, three planted defects are each caught, two
negative controls (wrong declared_n, an empty rule set), one control proving the
label check is spec-driven rather than hard-coded, and one closing law - the gate
can fail. A gate that has never failed anything certifies nothing.

Exit code: 0 = every leg as expected; 1 = a leg did not behave as expected.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cohortgate as cg  # noqa: E402

FIX = os.path.join(HERE, "fixtures")
LEGS = []


def run(spec, data):
    """Run the gate once; return (rc, findings, counts, n_rows).

    rc == 3 means the spec was refused, so the spec is not re-read.
    """
    rc = cg.main([os.path.join(FIX, spec), os.path.join(FIX, data), "--quiet"])
    try:
        spec_obj = cg.load_spec(os.path.join(FIX, spec))
        rows = cg.load_rows(os.path.join(FIX, data))
    except cg.SpecError:
        return rc, [], {}, 0
    findings, counts, n_rows = cg.audit(rows, spec_obj)
    return rc, findings, counts, n_rows


def leg(name, ok, detail=""):
    LEGS.append((name, ok, detail))
    print(f"  {'PASS' if ok else 'FAIL'}  {name}" + (f"  - {detail}" if detail else ""))


def kinds(findings):
    return {f["kind"] for f in findings}


# --- L1 clean fixture: passes, and the two independent routes agree --------
rc, f, c, n = run("spec.json", "clean.csv")
leg("L1 clean fixture passes, invariant holds",
    rc == 0 and not f and c == {"balanced": 2, "glycolytic": 2, "hepatic": 2} and n == 6,
    f"rc={rc} counts={c} n_rows={n} sum(labels)={sum(c.values())}")

# --- L2 planted defect: hand-typed label disagrees with the rule -----------
rc, f, _, _ = run("spec.json", "planted_label.csv")
lm = [x for x in f if x["kind"] == "label_mismatch"]
leg("L2 planted label mismatch is caught",
    rc == 2 and len(lm) == 1 and lm[0]["row"] == 5 and lm[0]["id"] == "S4",
    f"rc={rc} label_mismatch={len(lm)}" + (f" @row{lm[0]['row']} id={lm[0]['id']}" if lm else ""))

# --- L3 planted defect: duplicate id --------------------------------------
rc, f, _, _ = run("spec.json", "planted_dup.csv")
dup = [x for x in f if x["kind"] == "duplicate_id"]
leg("L3 planted duplicate id is caught",
    rc == 2 and len(dup) == 1 and dup[0]["id"] == "S5",
    f"rc={rc} duplicate_id={len(dup)}")

# --- L4 planted defect: a row outside every rule (the invariant fires too) --
rc, f, _, _ = run("spec.json", "planted_unclassified.csv")
k = kinds(f)
leg("L4 unclassified row caught, count invariant fires independently",
    rc == 2 and "unclassified" in k and "count_invariant" in k,
    f"rc={rc} kinds={sorted(k)}")

# --- L5 negative control: wrong declared_n (fail-closed) -------------------
rc, f, _, _ = run("spec_n7.json", "clean.csv")
leg("L5 wrong declared_n is caught (fail-closed)",
    rc == 2 and "declared_n_mismatch" in kinds(f),
    f"rc={rc} kinds={sorted(kinds(f))}")

# --- L6 negative control: an empty rule set must be refused ----------------
rc, _, _, _ = run("spec_empty.json", "clean.csv")
leg("L6 empty rule set is refused, not silently passed",
    rc == 3, f"rc={rc} (3 = SPEC_OR_IO: refuses to guess)")

# --- L7 the label check is spec-driven, not hard-coded ---------------------
rc, f, _, _ = run("spec_nolabelcheck.json", "planted_label.csv")
leg("L7 label check is spec-driven (turning it off makes the planted case clean)",
    rc == 0 and not f,
    f"rc={rc} (with check_declared_label off, the same planted data reads CLEAN)")

# --- L8 the closing law: the gate can fail --------------------------------
r = {}
for name in ("clean.csv", "planted_label.csv", "planted_dup.csv", "planted_unclassified.csv"):
    r[name], _, _, _ = run("spec.json", name)
can_fail = all(v == 2 for k_, v in r.items() if k_ != "clean.csv") and r["clean.csv"] == 0
leg("L8 the gate can fail: every planted fixture trips it, the clean twin passes",
    can_fail, " ".join(f"{k_}={v}" for k_, v in r.items()))

# --- summary --------------------------------------------------------------
fails = [name for name, ok, _ in LEGS if not ok]
print()
print(f"cohortgate battery: {len(LEGS) - len(fails)}/{len(LEGS)} legs pass"
      + (f" - FAIL: {fails}" if fails else " - OVERALL PASS"))
sys.exit(1 if fails else 0)

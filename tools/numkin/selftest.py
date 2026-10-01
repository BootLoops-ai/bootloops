#!/usr/bin/env python3
"""numkin smoke battery — what runs without a kira+FireFly stack:

the basisland member (the basis-solve backend) end-to-end on a
synthetic-truth system, plus its documented fail-closed refusal:

  1. rank: 5 sample rows over 3 basis symbols -> rank 3, multi-prime agree.
  2. solve: overdetermined 5x3 system with planted x = (2, -3/7, 5)
     -> exact vector recovered, exact_verified 5/5 rows.
  3. fail-closed: a pivot column with no support -> RANK DEFICIT naming the
     column, exit 2 (never a silent pad).

The sweep/harvest/shift_opt/etarerun members need staged kira SYSTEM files
and a kira+FireFly+Fermat stack; they are not run here.  Exit 0 only if all
three legs behave as documented.
"""
import json
import os
import subprocess
import sys
import tempfile
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
BL = os.path.join(HERE, "basisland.py")
fails = []


def grade(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name} {detail}")
    if not ok:
        fails.append(name)


tmp = tempfile.mkdtemp(prefix="numkin_selftest_")
# synthetic truth: f = 2*b1 - 3/7*b2 + 5*b3, sampled at 5 generic points
pts = [(1, 2, 3), (5, 1, 1), (2, 7, 1), (3, 3, 4), (1, 1, 9)]
x_true = (Fr(2), Fr(-3, 7), Fr(5))
rows = []
for b1, b2, b3 in pts:
    w = {"b1": str(Fr(b1, 3)), "b2": str(Fr(b2, 2)), "b3": str(Fr(b3))}
    rhs = x_true[0] * Fr(b1, 3) + x_true[1] * Fr(b2, 2) + x_true[2] * b3
    rows.append((w, str(rhs)))

# 1. rank verb
rank_in = os.path.join(tmp, "rank_in.json")
json.dump({"rows": [w for w, _ in rows]}, open(rank_in, "w"))
rank_out = os.path.join(tmp, "rank_out.json")
r = subprocess.run([sys.executable, BL, "rank", rank_in, rank_out],
                   capture_output=True, text=True)
d = json.load(open(rank_out)) if r.returncode == 0 else {}
grade("basisland rank (multi-prime agreement)",
      r.returncode == 0 and d.get("rank") == 3
      and d.get("pivot_symbols") == ["b1", "b2", "b3"], str(d))

# 2. solve verb: exact landing + mandatory exact verification
solve_in = os.path.join(tmp, "solve_in.json")
json.dump({"cols": ["b1", "b2", "b3"],
           "rows": [{"w": w, "rhs": [rhs]} for w, rhs in rows]},
          open(solve_in, "w"))
solve_out = os.path.join(tmp, "solve_out.json")
r = subprocess.run([sys.executable, BL, "solve", solve_in, solve_out, "8"],
                   capture_output=True, text=True)
d = json.load(open(solve_out)) if r.returncode == 0 else {}
grade("basisland solve lands the planted vector exactly",
      r.returncode == 0 and d.get("x") == [["2", "-3/7", "5"]]
      and d.get("exact_verified") == [5],
      str({k: d.get(k) for k in ("x", "exact_verified")}))
if r.returncode != 0:
    print((r.stdout + r.stderr)[-800:])

# 3. fail-closed: no support for pivot column b3 -> exit 2, column named
bad_in = os.path.join(tmp, "bad_in.json")
json.dump({"cols": ["b1", "b2", "b3"],
           "rows": [{"w": {"b1": "1"}, "rhs": ["1"]},
                    {"w": {"b2": "1"}, "rhs": ["1"]},
                    {"w": {"b1": "1", "b2": "2"}, "rhs": ["3"]}]},
          open(bad_in, "w"))
r = subprocess.run([sys.executable, BL, "solve", bad_in,
                    os.path.join(tmp, "bad_out.json"), "2"],
                   capture_output=True, text=True)
grade("solve refuses rank deficit (exit 2, names the column)",
      r.returncode == 2 and "RANK DEFICIT" in r.stdout and "b3" in r.stdout)

if fails:
    print(f"OVERALL: FAIL {fails}")
    sys.exit(1)
print("OVERALL: PASS (smoke; sweep/harvest/shift_opt legs need a kira+FireFly stack)")

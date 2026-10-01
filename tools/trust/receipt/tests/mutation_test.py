#!/usr/bin/env python3
"""mutation_test.py — mutation-test the CHECKER itself (scratch copies ONLY,
never in place).

A broken checker that passes everything is the worst outcome (the whole
point of the tool is that verification is louder than the solver). So: copy
core.py to a scratch dir, apply one checker-breaking mutation at a time, and
require that the sabotage suite CATCHES the mutant — i.e. at least one
sabotaged witness that the genuine checker rejects is ACCEPTED by the mutant
(proving the suite is sensitive to that breakage), or the positive control
breaks (mutant fails loudly either way -> counted as caught-by-suite).

Mutations (checker-killing classes):
  M1 verdict short-circuit: `ok = (r == expect)` -> `ok = True`
  M2 residual ignored:      `r == expect` -> `True or (r == expect)`
  M3 c ignored: expected residual drops the -c_m entries (only e_t checked)
  M4 modulus dropped in accumulation (accumulate over Z, no % p)
  M5 table comparison short-circuit: `claimed != wit["c"]` -> `False`

Run: python3 mutation_test.py [scratch_parent]
"""
import importlib.util
import json
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
TOOL = os.path.dirname(HERE)
sys.path.insert(0, TOOL)

from core import load_witness, verify_row as verify_row_genuine  # noqa: E402

# Scratch is ALWAYS a fresh mkdtemp dir — an argument only names its PARENT,
# and is ignored if it resolves into the package source tree (under
# `pytest <pkg>` argv[1] is the package dir itself; cleanup must never be
# able to delete the package).
PKG = os.path.realpath(TOOL)


def _inside_pkg(path):
    rp = os.path.realpath(path)
    return rp == PKG or rp.startswith(PKG + os.sep)


_parent = sys.argv[1] if len(sys.argv) > 1 else None
if _parent and not _inside_pkg(_parent):
    os.makedirs(_parent, exist_ok=True)
else:
    _parent = None
SCRATCH = tempfile.mkdtemp(prefix="receipt_mutation_", dir=_parent)
FAILED = []

MUTATIONS = {
    "M1_verdict_shortcircuit": (
        "    ok = (r == expect)",
        "    ok = True  # MUTANT"),
    "M2_residual_ignored": (
        "    ok = (r == expect)",
        "    ok = True or (r == expect)  # MUTANT"),
    "M3_c_ignored": (
        "    for m, cm in wit[\"c\"].items():\n"
        "        cm %= p\n"
        "        if cm:\n"
        "            expect[m] = (-cm) % p",
        "    pass  # MUTANT: claimed coefficients never enter the check"),
    "M4_modulus_dropped": (
        "            r[ccol] = (r.get(ccol, 0) + li * v) % p",
        "            r[ccol] = r.get(ccol, 0) + li * v  # MUTANT"),
    "M5_table_shortcircuit": (
        "        if claimed != wit[\"c\"]:",
        "        if False:  # MUTANT"),
}

# ---------------------------------------------------------------- fixture
P = 2147483647
ROWS = [
    {203: 1, 101: (-5) % P},
    {202: 1, 203: (-3) % P, 102: (-7) % P},
    {201: 1, 202: (-2) % P, 101: (-1) % P},
]
import numpy as np                               # noqa: E402
from emitter import FpSystem, emit_for_target    # noqa: E402
from core import write_witness                   # noqa: E402

rng = np.random.default_rng(7)
S = FpSystem([list(r.items()) for r in ROWS], [201, 202, 203],
             sorted({c for r in ROWS for c in r}), P)
rec, lam, _ = emit_for_target(S, 202, {101, 102}, rng)
assert rec["status"] == "CERTIFIED"
nz = np.flatnonzero(lam)
wp = os.path.join(SCRATCH, "w_good.json")
write_witness(wp, p=P, target_col=202, c=rec["c"],
              lam_idx=[int(i) for i in nz],
              lam_val=[int(lam[i]) for i in nz], n_rows=len(ROWS))
GOOD = load_witness(wp)


def sabotage_vectors():
    """(name, witness, table_row) triples the GENUINE checker must reject."""
    out = []
    for name, cmut in [("S1_perturb", {101: 16, 102: 7}),
                       ("S2_swap", {101: 7, 102: 15}),
                       ("S3_scale", {101: 30, 102: 14})]:
        w = json.load(open(wp))
        w["c"] = {str(k): v for k, v in cmut.items()}
        q = os.path.join(SCRATCH, f"w_{name}.json")
        json.dump(w, open(q, "w"))
        out.append((name, load_witness(q), None))
    w = json.load(open(wp))
    w["lam"]["val"][0] = (w["lam"]["val"][0] + 1) % P
    q = os.path.join(SCRATCH, "w_S4.json")
    json.dump(w, open(q, "w"))
    out.append(("S4_lambda_perturb", load_witness(q), None))
    w = json.load(open(wp))
    w["lam"]["idx"], w["lam"]["val"] = w["lam"]["idx"][:-1], w["lam"]["val"][:-1]
    q = os.path.join(SCRATCH, "w_S5.json")
    json.dump(w, open(q, "w"))
    out.append(("S5_lambda_truncation", load_witness(q), None))
    out.append(("S6_wrong_table", GOOD, {101: 16, 102: 7}))
    return out


VECTORS = sabotage_vectors()

# genuine checker: all sabotages rejected, control passes
ok, _ = verify_row_genuine(ROWS, GOOD)
assert ok, "genuine checker rejects the positive control — fixture broken"
for name, wit, trow in VECTORS:
    ok, _ = verify_row_genuine(ROWS, wit, table_row=trow)
    assert not ok, f"genuine checker fails to reject {name} — suite broken"
print("[PASS] genuine checker: control passes, 6/6 sabotages rejected")


def load_mutant(tag, old, new):
    src = open(os.path.join(TOOL, "core.py")).read()
    if old not in src:
        raise AssertionError(f"{tag}: mutation anchor not found in core.py")
    mdir = os.path.join(SCRATCH, tag)
    os.makedirs(mdir, exist_ok=True)
    with open(os.path.join(mdir, "core_mutant.py"), "w") as fh:
        fh.write(src.replace(old, new, 1))
    spec = importlib.util.spec_from_file_location(
        f"core_mutant_{tag}", os.path.join(mdir, "core_mutant.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


for tag, (old, new) in MUTATIONS.items():
    mod = load_mutant(tag, old, new)
    control_ok, _ = mod.verify_row(ROWS, GOOD)
    accepted = []
    for name, wit, trow in VECTORS:
        ok, _ = mod.verify_row(ROWS, wit, table_row=trow)
        if ok:
            accepted.append(name)
    # the suite catches the mutant iff it accepts a sabotage OR kills control
    caught = bool(accepted) or not control_ok
    tagr = "PASS" if caught else "FAIL"
    print(f"[{tagr}] {tag}: mutant "
          f"{'accepts ' + ','.join(accepted) if accepted else 'breaks control'}"
          f" -> suite would catch it")
    if not caught:
        FAILED.append(tag)

# cleanup: rmtree ONLY the mkdtemp scratch dir, never anything in the package tree
if _inside_pkg(SCRATCH):
    print(f"[guard] refusing to delete scratch inside the package tree: {SCRATCH}")
else:
    shutil.rmtree(SCRATCH, ignore_errors=False)
print(f"\nmutation_test: "
      f"{'ALL MUTANTS CAUGHT' if not FAILED else 'MISSED: ' + str(FAILED)}")
sys.exit(0 if not FAILED else 1)

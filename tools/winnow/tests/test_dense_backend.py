#!/usr/bin/env python3
"""Dense (torch) backend gates — backends/dense.py + backends/stratified.py.

Legs:
  a. synthetic stratified fixture (3 strata, masters forbidden, crafted
     master-relation leftover + pass-down rows): dense vs
     engine.stratified_solve — subs dict equality (keys, values, INSERTION
     ORDER), leftover multiset equality, bank-table byte equality via a
     mock table; plus a mixed dense/fallback budget run (same answer on
     every path mix) and a single-block eliminate_dense vs eliminate_fast
     contract check.
  b. GOLDEN GATE (archived; skips without WINNOW_BANKED_ROOT): the
     k2disp reference T bank replayed through backend='dense'
     (device cpu) for BOTH policies — produced .bin BYTE-EQUAL to the
     archived golden T; same for lp1disp (B2FT, fresh-stagepack Generate,
     load measured seconds-class); final subs/subs_raw identical to the
     stock backend on a slice.
  c. mutation control: one arithmetic op in the dense kernel corrupted
     (+1 on one value) -> the dense-vs-engine gate MUST fail loudly.
  d. typed refusals: p >= 2^31 (ValueError), transcript kwarg (TypeError —
     not a parameter by design), witnesses='native'+dense and non-fastpath
     policy+dense (NotImplementedError), unknown backend (ValueError);
     caps honored on the dense path (CensoredError via the same _CapBank
     proxy); witnesses='retrofit' works with backend='dense' (post-hoc);
     ledger carries "backend" only on dense runs (default path unchanged).
"""
import _common  # noqa: F401  MUST be first: sys.path bootstrap for the package

import os
import random
import sys

try:
    import torch
except ImportError:
    print("SKIP: the dense backend needs torch"
          " (optional extra: pip install 'ibplapper[dense]')")
    sys.exit(0)

import ibplapper as lap
from ibplapper import engine
from ibplapper.adapters import kira
from ibplapper.backends import dense as dense_mod
from ibplapper.backends.dense import eliminate_dense
from ibplapper.backends.stratified import stratified_solve_dense
from ibplapper.bank import InterfaceTable
from _common import NOMAT, SLICES4, check, finish, out_dir

torch.set_num_threads(8)        # tests stay <= 8 cores
OUT = out_dir("/tmp/ibplapper_test_dense")
P = 1048573


def canon_left(leftover):
    return sorted(tuple(sorted(r.items())) for r in leftover)


class MockTable:
    """Bank-protocol mock capturing the exact bytes the real Bank would
    append (InterfaceTable._encode is the lifted serializer)."""

    def __init__(self):
        self.blobs = []

    def append_rows(self, cell_key, node_key, rows, family=None):
        self.blobs.append((cell_key, node_key, InterfaceTable._encode(rows)))


# ---------------------------------------------------- a. synthetic fixture
rng = random.Random(42)
N_COLS = 90
ORDER = {c: c for c in range(N_COLS)}
STRAT = {c: c // 30 for c in range(N_COLS)}
MASTERS = {1000, 1001, 1002, 1003}
ROWS = []
for t in (2, 1, 0):
    for _ in range(24):
        r = {rng.randrange(30 * t, 30 * t + 30): rng.randrange(1, P)}
        for _ in range(4):
            r[rng.randrange(0, 30 * t + 30)] = rng.randrange(1, P)
        if rng.random() < 0.5:
            r[rng.choice(sorted(MASTERS))] = rng.randrange(1, P)
        ROWS.append(r)
# crafted: master-only relation (leftover) ...
ROWS.append({50: 7, 1000: 11})
ROWS.append({50: 14, 1000: 22, 1001: 5})
# ... and a stratum-2 row reducing to lower-stratum content (pass-down)
ROWS.append({70: 3, 5: 9})
ROWS.append({70: 6, 5: 4, 1002: 8})

ref = {}
for pol in ("B2F", "B2FT"):
    ta, tb = MockTable(), MockTable()
    subs_a, left_a, st_a = engine.stratified_solve(
        [dict(r) for r in ROWS], ORDER, P, STRAT, policy=pol,
        forbid=MASTERS, table=ta, cell_key="tiny", node_key="nX")
    subs_b, left_b, st_b = stratified_solve_dense(
        [dict(r) for r in ROWS], ORDER, P, STRAT, policy=pol,
        forbid=MASTERS, table=tb, cell_key="tiny", node_key="nX")
    ref[pol] = (subs_a, left_a, ta.blobs)
    check(f"synth_subs_equal::{pol}", subs_a == subs_b)
    check(f"synth_subs_insertion_order::{pol}",
          list(subs_a) == list(subs_b))
    check(f"synth_leftover_multiset::{pol}",
          canon_left(left_a) == canon_left(left_b))
    check(f"synth_bank_byte_equal::{pol}", ta.blobs == tb.blobs)
    check(f"synth_all_strata_dense_path::{pol}",
          all(s["path"] == "dense" for s in st_b.values()))
# fixture honesty: the crafted rows actually exercise both mechanisms
check("synth_fixture_has_leftover", len(ref["B2F"][1]) >= 1,
      f"n_leftover={len(ref['B2F'][1])}")
check("synth_fixture_has_passdown",
      any(s["rows_passed_down"] >= 1 for s in st_a.values()))

# mixed dense/fallback budget: same answer on every path mix
sizes = sorted(s["dense_bytes"] for s in st_b.values())
thr = sizes[len(sizes) // 2]
subs_c, left_c, st_c = stratified_solve_dense(
    [dict(r) for r in ROWS], ORDER, P, STRAT, policy="B2FT",
    forbid=MASTERS, dense_budget_bytes=thr)
check("synth_mixed_paths_present",
      {s["path"] for s in st_c.values()} == {"dense", "fallback-cpu"},
      {s["path"] for s in st_c.values()})
check("synth_mixed_answer_identical",
      subs_c == ref["B2FT"][0] and list(subs_c) == list(ref["B2FT"][0])
      and canon_left(left_c) == canon_left(ref["B2FT"][1]))

# single-block contract: eliminate_dense vs eliminate_fast directly
blk_rows = [dict(r) for r in ROWS]
s_f, l_f, t_f = engine.eliminate_fast([dict(r) for r in blk_rows], ORDER, P,
                                      forbid=MASTERS)
s_d, l_d, t_d = eliminate_dense([dict(r) for r in blk_rows], ORDER, P,
                                forbid=MASTERS)
check("block_subs_equal_and_ordered",
      s_f == s_d and list(s_f) == list(s_d))
check("block_leftover_multiset", canon_left(l_f) == canon_left(l_d))
check("block_pivots_match", t_f["pivots"] == t_d["pivots"],
      f"{t_f['pivots']} vs {t_d['pivots']}")
check("block_stats_backend_tag", t_d["backend"] == "dense-torch/cpu",
      t_d["backend"])

# ------------- b. GOLDEN GATE (shared impl: _dense_golden.golden_gates;
# a cuda harness can run the SAME gates with --device cuda)
from _dense_golden import golden_gates                        # noqa: E402

if NOMAT is None:
    print("SKIP leg b: the archived golden banks need WINNOW_BANKED_ROOT")
else:
    GEN_K2 = os.path.join(NOMAT, "gen_lbl3m2L_k2disp")
    for name, ok in golden_gates(OUT, device="cpu"):
        check(name, ok)

    # final subs identical to the stock path (k2disp slice 0, both
    # policies)
    system, _ = kira.load(GEN_K2, "lbl3m2L_k2disp", *SLICES4[0])
    for pol in ("B2F", "B2FT"):
        res_s = lap.eliminate(system, lap.Schedule(policy=pol))
        res_d = lap.eliminate(system, lap.Schedule(policy=pol),
                              backend="dense")
        check(f"k2disp_final_subs_identical::{pol}",
              res_s.subs == res_d.subs
              and list(res_s.subs) == list(res_d.subs))
        check(f"k2disp_subs_raw_identical::{pol}",
              res_s.subs_raw == res_d.subs_raw
              and list(res_s.subs_raw) == list(res_d.subs_raw))
        check(f"k2disp_leftover_identical::{pol}",
              canon_left(res_s.leftover) == canon_left(res_d.leftover))
    check("ledger_backend_key_dense_only",
          res_d.ledger.get("backend") == "dense-torch/cpu"
          and "backend" not in res_s.ledger)

# ----------------------------------------------------- c. mutation control
_honest_axpy = dense_mod._axpy_block
_armed = {"fire": True, "fired": False}


def _corrupt_axpy(A, sel, jpos, pivrow, p):
    n = _honest_axpy(A, sel, jpos, pivrow, p)
    if _armed["fire"]:
        _armed["fire"] = False
        _armed["fired"] = True
        A[sel[0], -1] = (A[sel[0], -1] + 1) % p     # +1 on ONE value
    return n


dense_mod._axpy_block = _corrupt_axpy
try:
    tm = MockTable()
    subs_m, left_m, _ = stratified_solve_dense(
        [dict(r) for r in ROWS], ORDER, P, STRAT, policy="B2FT",
        forbid=MASTERS, table=tm, cell_key="tiny", node_key="nX")
finally:
    dense_mod._axpy_block = _honest_axpy
check("mutation_actually_fired", _armed["fired"])
subs_a, left_a, blobs_a = ref["B2FT"]
detected = not (subs_m == subs_a and canon_left(left_m) == canon_left(left_a)
                and tm.blobs == blobs_a)
check("mutation_control_gate_fails_loudly", detected)

# ------------------------------------------------------ d. typed refusals
try:
    eliminate_dense([{3: 1, 1: 2}], {3: 3}, 2**31 + 11, forbid={1})
    check("refusal_p_ge_2p31", False, "no raise")
except ValueError as e:
    check("refusal_p_ge_2p31", "int64" in str(e), str(e)[:70])

try:
    stratified_solve_dense([{3: 1}], {3: 3}, P, {3: 0}, transcript=object())
    check("refusal_transcript_not_a_parameter", False, "accepted transcript")
except TypeError:
    check("refusal_transcript_not_a_parameter", True)

sys_small = lap.System([{9: 1, 1: 1}], P, {9: 9}, forbid={1})
try:
    lap.eliminate(sys_small, lap.Schedule(policy="B2F"),
                  witnesses="native", backend="dense")
    check("refusal_native_witness_dense", False)
except NotImplementedError as e:
    check("refusal_native_witness_dense", "retrofit" in str(e))

try:
    lap.eliminate(sys_small, lap.Schedule(policy="B2"), backend="dense")
    check("refusal_policy_B2_dense", False)
except NotImplementedError:
    check("refusal_policy_B2_dense", True)

try:
    lap.eliminate(sys_small, backend="tpu")
    check("refusal_unknown_backend", False)
except ValueError:
    check("refusal_unknown_backend", True)

# caps honored on the dense path (same _CapBank proxy, loud + partial ledger)
try:
    lap.eliminate(sys_small,
                  lap.Schedule(policy="B2F", caps={"stratum_fill": 0}),
                  backend="dense")
    check("caps_honored_on_dense", False, "no exception")
except lap.CensoredError as e:
    check("caps_honored_on_dense",
          e.ledger is not None
          and e.ledger["censored"]["cap"] == "stratum_fill")

# retrofit witnesses work with backend='dense' (post-hoc by construction).
# NOTE: retrofit's claim check needs a span-UNIQUE representation, so this
# leg uses a FULL-RANK square fixture — on the rank-deficient fixture above
# CLAIM_MISMATCH fires on the stock cpu backend too (measured; master
# relations make the master-coefficient representation ambiguous), i.e.
# that is a fixture property, not a backend property.
rng2 = random.Random(7)
ORDER2 = {c: c for c in range(30)}
STRAT2 = {c: c // 10 for c in range(30)}
MASTERS2 = {500, 501}
ROWS2 = []
for c in range(29, -1, -1):
    r = {c: rng2.randrange(1, P)}
    if c > 0:
        r[rng2.randrange(0, c)] = rng2.randrange(1, P)
    r[500 if c % 2 else 501] = rng2.randrange(1, P)
    ROWS2.append(r)
sys2 = lap.System([dict(r) for r in ROWS2], P, ORDER2, forbid=MASTERS2,
                  stratum_of=STRAT2)
targets = [29, 15]
res_w = lap.eliminate(sys2, lap.Schedule(policy="B2FT"), backend="dense",
                      witnesses="retrofit", witness_targets=targets)
check("retrofit_on_dense_certified_verified",
      all(res_w.witnesses[t]["status"] == "CERTIFIED"
          and res_w.witnesses[t]["verify_pass"] for t in targets),
      {t: res_w.witnesses[t]["status"] for t in targets})
res_w_cpu = lap.eliminate(sys2, lap.Schedule(policy="B2FT"),
                          witnesses="retrofit", witness_targets=targets)
check("retrofit_dense_matches_cpu_records",
      all(res_w.witnesses[t]["c"] == res_w_cpu.witnesses[t]["c"]
          and res_w.witnesses[t]["lam"] == res_w_cpu.witnesses[t]["lam"]
          for t in targets))

finish("test_dense_backend")

#!/usr/bin/env python3
"""test_banked.py — regression against the ARCHIVED probe artifacts.

Internal-reference leg: needs the archived artifact store (not shipped).
Set RECEIPT_BANKED_ROOT to its root; unset, this test SKIPS cleanly.

Legs (all measured priors from the archived probe runs; no new claims):
  1. Archived witness pairs: the 30 probe-era (c_*.json, lam_*.npy) pairs
     (three archived families, 10 pairs each) must verify PASS through the
     legacy reader — the independent-checker record they carry is 30/30
     (CHECK_*.json).
  2. Retrofit regression: emit over the archived full 42-row registry table
     at both archived primes must reproduce 84/84 CERTIFIED (the archived
     retrofit record: 42+42 oracle_match + checker_pass), and all 84 emitted
     v1 witnesses must verify PASS. Verify throughput is measured here
     (archived prior: 0.96 ms/row mean).
  3. CANONICAL MUST-FAIL — the archived confident-wrong-table exhibit:
     3 rows x 2 primes of D-ladder-stable, rank-clean, 2-prime-CRT-consistent
     WRONG values. Every one must FAIL against the lambda certificates; the
     archived oracle rows (ORACLE_BLOCKS.json) must PASS as positive control.
     rc discipline is exercised through the CLI on both tables.

Run: python3 test_banked.py <out_dir>
"""
import json
import os
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
TOOL = os.path.dirname(HERE)
sys.path.insert(0, TOOL)

RL = os.environ.get("RECEIPT_BANKED_ROOT")  # archived reference artifacts (not shipped)
if not RL:
    print("SKIP: RECEIPT_BANKED_ROOT unset — the archived-reference leg needs the"
          " reference artifact store (not shipped)")
    sys.exit(0)

import core                                      # noqa: E402
from core import load_witness, verify_row        # noqa: E402
from adapters import strata                      # noqa: E402

NOMAT = os.path.join(RL, "newmove_nomat_p0")
COPAIR = os.path.join(RL, "copair_day1")
D3 = os.path.join(RL, "strata_day3_tinyrung", "staging")
D4 = os.path.join(RL, "strata_day4_smallrung", "staging")
REGISTRY = "e2d8acded38e6885"
PRIMES = [2147483647, 2147483629]
D0, E0 = 1234577, 87654321

# scratch goes under argv[1], EXCEPT when that resolves into the package
# source tree (under `pytest <pkg>` argv[1] is the package dir itself) —
# then, or with no argument, a fresh mkdtemp dir is used instead.
PKG = os.path.realpath(TOOL)
_out = sys.argv[1] if len(sys.argv) > 1 else None
if _out and os.path.realpath(_out) != PKG \
        and not os.path.realpath(_out).startswith(PKG + os.sep):
    OUT = _out
    os.makedirs(OUT, exist_ok=True)
else:
    OUT = tempfile.mkdtemp(prefix="receipt_test_banked_")
FAILED = []
CK = os.path.join(OUT, "checkpoints.jsonl")


def ck(stage, obj):
    with open(CK, "a") as fh:
        fh.write(json.dumps({"t": time.strftime("%H:%M:%S"), "stage": stage,
                             **obj}) + "\n")


def check(name, cond, extra=""):
    tag = "PASS" if cond else "FAIL"
    print(f"[{tag}] {name} {extra}", flush=True)
    ck("check", {"name": name, "pass": bool(cond), "extra": str(extra)[:200]})
    if not cond:
        FAILED.append(name)


FAMS = {
    "vac3": (os.path.join(NOMAT, "gen_vac3"), os.path.join(D3, "vac3")),
    "lbl3m2L_k2disp": (os.path.join(NOMAT, "gen_lbl3m2L_k2disp"),
                       os.path.join(D4, "lbl3m2L_k2disp")),
    "vac3_b0_r2_f2_sm0": (os.path.join(NOMAT, "gen_vac3_b0_r2_f2_sm0"),
                          os.path.join(D3, "vac3_b0_r2_f2_sm0")),
}

# ---------------------------------------------- Leg 1: archived pairs, 30/30
t_leg1 = time.time()
n_pairs = 0
for fam, (art, stg) in FAMS.items():
    for p in PRIMES:
        import glob as _g
        pairs = sorted(_g.glob(os.path.join(NOMAT, f"c_{fam}_{p}_*.json")))
        if not pairs:
            continue
        sysd = strata.load_rows(art, fam, p, D0, E0)
        _s, info, _c = strata.stage_partition(art, stg, fam, p, D0, E0,
                                              sysd["sector_of"])
        for q in pairs:
            lw = strata.read_legacy_pair(q, info["masters_map"],
                                         len(sysd["rows"]))
            ok, det = verify_row(sysd["rows"], lw)
            check(f"banked_pair::{os.path.basename(q)}", ok,
                  det.get("reason", ""))
            n_pairs += 1
check("leg1_pair_count_30", n_pairs == 30, f"n={n_pairs}")
ck("leg1_done", {"wall_s": round(time.time() - t_leg1, 1), "n": n_pairs})

# ------------------------------- Leg 2: 84/84 emit regression + throughput
t_leg2 = time.time()
art, stg = FAMS["lbl3m2L_k2disp"]
emit_out = os.path.join(OUT, "emit_k2disp")
slices = ",".join(f"{p}:{D0}:{E0}" for p in PRIMES)
reuse = (os.environ.get("RECEIPT_TEST_REUSE_EMIT") == "1"
         and os.path.exists(os.path.join(emit_out, "EMIT_REPORT.json")))
if reuse:
    print("[note] RECEIPT_TEST_REUSE_EMIT=1: reusing existing emit output "
          "(dev shortcut; the reference record must come from a fresh run)")


class _R:
    returncode = 0
    stderr = ""


r = _R() if reuse else subprocess.run(
    [sys.executable, os.path.join(TOOL, "receipt.py"), "emit",
     "--adapter", "strata", "--system-dir", art, "--staging", stg,
     "--family", "lbl3m2L_k2disp",
     "--table", os.path.join(stg, "corpus", REGISTRY, "kira_target.m"),
     "--slices", slices, "--out", emit_out],
    capture_output=True, text=True)
check("emit_rc0", r.returncode == 0,
      f"rc={r.returncode} {r.stderr[-300:] if r.returncode else ''}")
rep = json.load(open(os.path.join(emit_out, "EMIT_REPORT.json")))
n_cert = sum(sl["n_certified"] for sl in rep["slices"])
n_rows = sum(sl["n_rows"] for sl in rep["slices"])
check("emit_84_of_84", n_cert == 84 and n_rows == 84,
      f"certified={n_cert}/{n_rows} anchored={rep['n_anchored']}")
ck("leg2_emit", {"wall_s": round(time.time() - t_leg2, 1),
                 "n_certified": n_cert,
                 "solve_totals": [sl["solve_wall_total_s"]
                                  for sl in rep["slices"]]})

# verify all 84 emitted witnesses + measure throughput
n84 = n84_pass = 0
t_pure = 0.0
for p in PRIMES:
    sysd = strata.load_rows(art, "lbl3m2L_k2disp", p, D0, E0)
    import glob as _g
    wfiles = sorted(_g.glob(os.path.join(emit_out,
                                         f"w_lbl3m2L_k2disp_{p}_*.json")))
    wits = [load_witness(q) for q in wfiles]
    t0 = time.time()
    for lw in wits:
        ok, det = verify_row(sysd["rows"], lw)
        n84 += 1
        n84_pass += bool(ok)
    t_pure += time.time() - t0
check("verify_84_of_84", n84 == 84 and n84_pass == 84,
      f"{n84_pass}/{n84}")
ms_per_row = 1000 * t_pure / max(n84, 1)
check("verify_throughput_ms_class", ms_per_row < 20,
      f"{ms_per_row:.2f} ms/row (archived prior 0.96 ms/row mean)")
ck("leg2_verify", {"ms_per_row": round(ms_per_row, 3), "n": n84})

# --------------------------- Leg 3: COPAIR wrong-table exhibit (MUST FAIL)
# The pairing engine produced sector-BLOCK coefficients (in-sector masters
# only; convention = eval_vec-negated, c = -block). The exhibit table row =
# the independently-evaluated TRUE oracle row with its in-sector entries
# overwritten by the pairing's values — the exact table a pairing pipeline
# would have shipped, differing ONLY in the wrong values.
t_leg3 = time.time()
pair = json.load(open(os.path.join(COPAIR, "PAIR_RESULT.json")))
oracle_blocks = json.load(open(os.path.join(COPAIR, "ORACLE_BLOCKS.json")))
SEC119 = ["[1, 1, 1, 0, 1, 1, 2, 0, 0]", "[1, 1, 2, 0, 1, 1, 2, 0, 0]",
          "[1, 2, 1, 0, 1, 1, 2, 0, 0]"]

# label -> emitted witness path, per prime
wit_by_label = {}
import glob as _g
for p in PRIMES:
    for q in _g.glob(os.path.join(emit_out, f"w_lbl3m2L_k2disp_{p}_*.json")):
        meta = json.load(open(q))
        wit_by_label[(str(meta["target"].get("label")), p)] = q

# independent oracle rows (registry table) + masters map (structural)
from adapters.strata import CoeffEvaluator, parse_oracle_table  # noqa: E402
sysd1, info, _closed = strata.stage_partition(
    art, stg, "lbl3m2L_k2disp", PRIMES[0], D0, E0,
    strata.load_rows(art, "lbl3m2L_k2disp", PRIMES[0], D0, E0)["sector_of"])
masters_map = info["masters_map"]
midx2w = {idx: w for w, idx in masters_map.items()}
orc = {str(list(t)): terms for t, terms in parse_oracle_table(
    os.path.join(stg, "corpus", REGISTRY, "kira_target.m"))}

n_mustfail = n_failed_correctly = n_control = n_control_pass = 0
wrong_tables, true_tables = {p: {} for p in PRIMES}, {p: {} for p in PRIMES}
for sec in pair["sectors"]:
    p = sec["p"]
    if sec["sector"] != 119 or p not in PRIMES:
        continue
    sysd = strata.load_rows(art, "lbl3m2L_k2disp", p, D0, E0)
    ev = CoeffEvaluator(p, D0, E0)
    for tstr in SEC119:
        wq = wit_by_label.get((tstr, p))
        check(f"exhibit_witness_exists::{tstr}|{p}", wq is not None)
        if wq is None:
            continue
        lw = load_witness(wq)
        # TRUE row: independent oracle evaluation (full row, all masters)
        true_c = strata.oracle_c(orc[tstr], ev, midx2w, p)
        assert not isinstance(true_c, str), (tstr, true_c)
        # convention gate: in-sector part of the true row must equal the
        # archived ORACLE_BLOCKS under c = -block
        blk_true = {midx2w[tuple(json.loads(m))]: (-v) % p
                    for m, v in oracle_blocks[f"{tstr}|{p}"].items()}
        check(f"block_convention_gate::{tstr}|{p}",
              all(true_c.get(w) == v for w, v in blk_true.items()))
        # WRONG row: overwrite in-sector entries with the pairing's values
        wrong_c = dict(true_c)
        for m, v in sec["blocks"][tstr].items():
            wrong_c[midx2w[tuple(json.loads(m))]] = (-v) % p
        assert wrong_c != true_c, f"pairing row equals oracle row? {tstr}|{p}"
        ok, det = verify_row(sysd["rows"], lw, table_row=wrong_c)
        n_mustfail += 1
        n_failed_correctly += (not ok)
        check(f"MUSTFAIL_copair_wrong_row::{tstr}|{p}", not ok,
              det.get("reason", "")[:60])
        # positive control: the independently-evaluated true row must PASS
        ok2, det2 = verify_row(sysd["rows"], lw, table_row=true_c)
        n_control += 1
        n_control_pass += bool(ok2)
        check(f"control_oracle_row_passes::{tstr}|{p}", ok2,
              det2.get("reason", ""))
        wrong_tables[p][lw["target_col"]] = {str(k): v
                                             for k, v in wrong_c.items()}
        true_tables[p][lw["target_col"]] = {str(k): v
                                            for k, v in true_c.items()}

check("exhibit_all_6_wrong_rows_failed",
      n_mustfail == 6 and n_failed_correctly == 6,
      f"{n_failed_correctly}/{n_mustfail}")
check("exhibit_all_6_controls_passed",
      n_control == 6 and n_control_pass == 6, f"{n_control_pass}/{n_control}")

# rc discipline through the CLI on the exhibit (one prime suffices for rc)
p = PRIMES[0]
wt = os.path.join(OUT, "copair_wrong_table.json")
tt = os.path.join(OUT, "copair_true_table.json")
json.dump(wrong_tables[p], open(wt, "w"))
json.dump(true_tables[p], open(tt, "w"))
wits_p = [wit_by_label[(t, p)] for t in SEC119]
base = [sys.executable, os.path.join(TOOL, "receipt.py"), "verify",
        "--adapter", "strata", "--system-dir", art, "--staging", stg,
        "--family", "lbl3m2L_k2disp", "--witness", *wits_p, "-q"]
r = subprocess.run(base + ["--table", wt], capture_output=True, text=True)
check("cli_rc1_on_wrong_table", r.returncode == 1, f"rc={r.returncode}")
r = subprocess.run(base + ["--table", tt], capture_output=True, text=True)
check("cli_rc0_on_true_table", r.returncode == 0,
      f"rc={r.returncode} {r.stderr[-200:] if r.returncode else ''}")
ck("leg3_done", {"wall_s": round(time.time() - t_leg3, 1)})

print(f"\ntest_banked: "
      f"{'ALL PASS' if not FAILED else 'FAILURES: ' + str(FAILED)}")
print(f"verify throughput: {ms_per_row:.2f} ms/row over {n84} rows")
sys.exit(0 if not FAILED else 1)

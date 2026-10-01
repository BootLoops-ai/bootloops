#!/usr/bin/env python3
"""Unit gates G1-G4 against the ibplapper API (gridmode is EXCLUDED from
v1) + API contract gates (caps loudness, leftover reporting,
native-witness smoke + policy refusal, Bank protocol immutability) +
coverage/planted-fault gates.

Independent oracle: dense flint nmod_mat rref (python-flint).
"""
import _common  # noqa: F401  MUST be first: sys.path bootstrap for the package

import os
import random

import flint

import ibplapper as lap
from ibplapper.coverage import coverage_gate, planted_fault_trials
from ibplapper.engine import b2b3_compare, eliminate, rref_canonical
from _common import check, finish, out_dir

OUT = out_dir("/tmp/ibplapper_test_gates")
PRIMES = [1048573, 2097143]


def random_system(rng, n_rows, n_cols, density, p):
    rows = []
    for _ in range(n_rows):
        r = {}
        for c in range(n_cols):
            if rng.random() < density:
                r[c] = rng.randrange(1, p)
        if r:
            rows.append(r)
    return rows


def flint_rref_rows(rows, order, p):
    cols = sorted(order, key=lambda c: -order[c])
    pos = {c: i for i, c in enumerate(cols)}
    nr, nc = len(rows), len(cols)
    flat = [0] * (nr * nc)
    for i, r in enumerate(rows):
        for c, v in r.items():
            flat[i * nc + pos[c]] = v % p
    M = flint.nmod_mat(nr, nc, flat, p)
    rr = M.rref()
    R = rr[0] if isinstance(rr, tuple) else rr
    out = []
    for i in range(nr):
        row = [(j, int(str(R[i, j]))) for j in range(nc) if R[i, j] != 0]
        if not row:
            continue
        lead_j, lead_v = row[0]
        assert lead_v == 1
        terms = tuple(sorted((cols[j], v) for j, v in row[1:]))
        out.append((cols[lead_j], terms))
    return sorted(out)


def recon_rows(subs, leftover=()):
    out = []
    for pc, pr in subs.items():
        row = {pc: 1}
        row.update(pr)          # int keys: never the dict(**) keyword form
        out.append(row)
    out.extend(dict(r) for r in leftover)
    return out


# ---------------------------------------------------- G1 exactness vs flint
n_ok = 0
for p in PRIMES:
    for seed in (11, 12, 13):
        rng = random.Random(seed)
        n_cols = 60
        rows = random_system(rng, 45, n_cols, 0.10, p)   # underdetermined
        order = {c: c for c in range(n_cols)}
        oracle = flint_rref_rows(rows, order, p)
        assert sum(1 for _, t in oracle if t) >= 10, "degenerate G1 system"
        mine = rref_canonical([dict(r) for r in rows], order, p)
        if not check(f"G1_rref_vs_flint::p{p}_s{seed}", mine == oracle):
            continue
        for pol, a in (("B2", 0), ("B3", 0), ("B3", 1), ("B3", 4)):
            shadow = {c: (c % 7) for c in range(n_cols)}
            subs, left, _ = eliminate([dict(r) for r in rows], order, p,
                                      policy=pol, alpha=a, shadow=shadow)
            canon = rref_canonical(recon_rows(subs, left), order, p)
            n_ok += check(f"G1_policy_canonical::{pol}_a{a}_p{p}_s{seed}",
                          canon == oracle)
check("G1_all_policy_runs", n_ok == 24, f"{n_ok}/24")

# ------------------------------ G2/G3 stratified == monolithic == flint + T
tbase = os.path.join(OUT, "T_g2")
for f in (tbase + ".bin", tbase + ".idx.json"):
    if os.path.exists(f):
        os.remove(f)
bank = lap.Bank(tbase, meta={"family": "synth", "gate": "G2"})
for p in PRIMES:
    rng = random.Random(101 + p)
    n_cols = 90
    stratum_of = {c: c // 30 for c in range(n_cols)}
    order = {c: c for c in range(n_cols)}
    rows = []
    for t in (2, 1, 0):
        for _ in range(18):
            r = {rng.randrange(30 * t, 30 * t + 30): rng.randrange(1, p)}
            for _ in range(5):
                r[rng.randrange(0, 30 * t + 30)] = rng.randrange(1, p)
            rows.append(r)
    oracle = flint_rref_rows(rows, order, p)
    mono = rref_canonical([dict(r) for r in rows], order, p)
    check(f"G2_monolithic_vs_flint::p{p}", mono == oracle)
    for pol in ("B2F", "B2FT"):
        sys_ = lap.System(rows, p, order, stratum_of=stratum_of)
        res = lap.eliminate(sys_, lap.Schedule(
            policy=pol, bank=bank if pol == "B2F" else None,
            bank_cell="cellA", bank_node=f"p{p}"))
        strat = rref_canonical(recon_rows(res.subs, res.leftover), order, p)
        check(f"G2_{pol}_vs_flint::p{p}", strat == oracle)
        check(f"G2_{pol}_ledger_strata::p{p}",
              set(res.ledger["per_stratum"]) == {"0", "1", "2"})
# G3: reopen + hash-verified round-trip + append-only immutability
bank2 = lap.Bank(tbase)
rows_back = bank2.get("cellA", f"p{PRIMES[0]}/t2", family="synth")
check("G3_bank_roundtrip_nonempty", bool(rows_back))
try:
    bank2.append_rows("cellA", f"p{PRIMES[0]}/t2", rows_back, family="synth")
    check("G3_bank_append_only", False, "re-append was ACCEPTED")
except AssertionError:
    check("G3_bank_append_only", True)

# ------------------------------------------------------------- G4 b2b3 hook
rng = random.Random(7)
p = PRIMES[0]
n_cols = 120
order = {c: c for c in range(n_cols)}
shadow = {c: (c * 13) % 9 for c in range(n_cols)}
rows = []
for lead in range(119, 20, -2):
    r = {lead: rng.randrange(1, p)}
    for dcol in (1, 2, 3, 7, 11):
        if lead - dcol >= 0 and rng.random() < 0.8:
            r[lead - dcol] = rng.randrange(1, p)
    rows.append(r)
masters = set(range(10))
jl = os.path.join(OUT, "b2b3_table.jsonl")
if os.path.exists(jl):
    os.remove(jl)
table = b2b3_compare(rows, order, p, shadow, masters, "banded_synth", jl)
check("G4_b2b3_hook", len(table) == 4 and os.path.exists(jl),
      f"{[r['policy'] for r in table]}")

# ------------------------------------------- API contract: leftover surfaced
p = PRIMES[0]
sys_ = lap.System([{5: 1, 1: 2}, {5: 2, 1: 4}, {1: 3, 2: 5}],
                  p, {5: 5}, forbid={1, 2})
res = lap.eliminate(sys_, lap.Schedule(policy="B2F"))
check("leftover_master_relation_reported",
      len(res.leftover) == 1 and set(res.leftover[0]) <= {1, 2},
      f"leftover={res.leftover}")
check("leftover_in_ledger", res.ledger["leftover"]["n"] == 1)
# master-only INPUT row also surfaces (pre-leftover path)
res = lap.eliminate(lap.System([{5: 1, 1: 2}, {1: 7}], p, {5: 5},
                               forbid={1}), lap.Schedule(policy="B2FT"))
check("masteronly_input_row_reported", len(res.leftover) == 1)

# order/forbid overlap refused
try:
    lap.System([{1: 1}], p, {1: 1}, forbid={1})
    check("order_forbid_overlap_refused", False)
except ValueError:
    check("order_forbid_overlap_refused", True)

# caps: CensoredError is loud and carries a partial ledger
try:
    lap.eliminate(lap.System([{9: 1, 1: 1}], p, {9: 9}, forbid={1}),
                  lap.Schedule(policy="B2F", caps={"stratum_fill": 0}))
    check("caps_censored_error_loud", False, "no exception")
except lap.CensoredError as e:
    check("caps_censored_error_loud",
          e.ledger is not None and e.ledger["censored"]["cap"]
          == "stratum_fill")

# witnesses='native' (in-elimination transcript) works on the fast path
# and refuses B2/B3 loudly (full battery: test_native_witness.py)
res = lap.eliminate(lap.System([{9: 1, 1: 1}], p, {9: 9}, forbid={1}),
                    witnesses="native", witness_targets=[9])
check("native_witness_smoke",
      res.witnesses[9]["status"] == "CERTIFIED"
      and res.witnesses[9]["verify_pass"]
      and res.ledger["witnesses"]["transcript_fill"]["n_registered"] >= 1)
try:
    lap.eliminate(lap.System([{9: 1, 1: 1}], p, {9: 9}, forbid={1}),
                  lap.Schedule(policy="B2"), witnesses="native")
    check("native_witness_b2_refused", False)
except NotImplementedError:
    check("native_witness_b2_refused", True)

# -------------------------------------- coverage gate + planted-fault trial
p = PRIMES[0]
rows = [{203: 1, 101: p - 5}, {202: 1, 203: p - 3, 102: p - 7},
        {201: 1, 202: p - 2, 101: p - 1}]
order = {201: 3, 202: 2, 203: 1}
res = lap.eliminate(lap.System(rows, p, order, forbid={101, 102}),
                    lap.Schedule(policy="B2FT"))
banked = res.subs_raw
masters = {101, 102}
fresh = [dict(r) for r in rows]     # original identities re-serve as fresh
rep, touch = coverage_gate(fresh, res.subs, banked, masters, p, c_min=1)
check("coverage_gate_pass", rep["pass"], rep)
pf = planted_fault_trials(fresh, banked, order, masters, p, n_trials=3)
check("planted_fault_all_detected", pf["all_detected"],
      f"rate={pf['detection_rate']}")

finish("test_gates")

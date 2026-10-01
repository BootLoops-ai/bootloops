#!/usr/bin/env python3
"""Kira adapter gates: version guard (GUARDED / MUST-REFUSE / ABSENT),
provenance gate, row-index integrity gate, and the T1 provenance counters
vs the recorded reference prior."""
import _common  # noqa: F401  MUST be first: sys.path bootstrap for the package

import os
import shutil
import sqlite3

import ibplapper as lap
from ibplapper.adapters import kira
from _common import SLICES4, T1GEN, check, finish, out_dir, skip_without

skip_without("test_kira_adapter", T1GEN)

OUT = out_dir("/tmp/ibplapper_test_kira_adapter")
p, d0, e0 = SLICES4[0]
FAM = "vac3_b0_r2_f2_sm0"

# guard: GUARDED on a pinned-kira artifact tree
g = kira.version_guard(T1GEN)
check("guard_status_guarded", g["status"] == "GUARDED", g)
check("guard_pins_echoed", tuple(g["WEIGHTBITS"]) == (5, 4, 17, 13)
      and g["INTEGRALORDERING"] == 5)

# loader counters vs the recorded reference prior (provenance invariant)
sysd = kira.load_system(T1GEN, FAM, p, d0, e0,
                        expect_counters={"eqs": 6001, "terms": 28089,
                                         "files": 10})
check("t1_counters_match_bench_prior", True,
      {k: sysd['counters'][k] for k in ('eqs', 'terms', 'files')})
check("t1_guard_in_counters",
      sysd["counters"]["version_guard"]["status"] == "GUARDED")

# provenance gate MUST refuse a wrong declaration
try:
    kira.load_system(T1GEN, FAM, p, d0, e0, expect_counters={"eqs": 9999})
    check("provenance_gate_refuses", False)
except kira.ProvenanceError:
    check("provenance_gate_refuses", True)

# full pipeline: System -> eliminate -> engine stats == BENCH prior
system, _ = kira.load(T1GEN, FAM, p, d0, e0)
res = lap.eliminate(system, lap.Schedule(policy="B2F"))
st = list(res.stats_per.values())
check("t1_solve_prior_counts",
      len(res.subs) == 6001 and not res.leftover
      and sum(s["pivots"] for s in st) == 6001,
      f"pivots={len(res.subs)}")

# version guard MUST-REFUSE: scratch artifact tree with alien WEIGHTBITS
alien = os.path.join(OUT, "alien_gen")
os.makedirs(os.path.join(alien, "results"), exist_ok=True)
shutil.rmtree(os.path.join(alien, "results"), ignore_errors=True)
os.makedirs(os.path.join(alien, "results"))
con = sqlite3.connect(os.path.join(alien, "results", "kira.db"))
con.execute("CREATE TABLE VERSION (number TEXT)")
con.execute("INSERT INTO VERSION VALUES ('9.9')")
con.execute("CREATE TABLE WEIGHTBITS (A INT, B INT, C INT, D INT)")
con.execute("INSERT INTO WEIGHTBITS VALUES (6, 4, 17, 13)")   # alien A
con.execute("CREATE TABLE INTEGRALORDERING (ID INT)")
con.execute("INSERT INTO INTEGRALORDERING VALUES (5)")
con.commit()
con.close()
try:
    kira.version_guard(alien)
    check("guard_refuses_alien_weightbits", False)
except kira.KiraVersionError:
    check("guard_refuses_alien_weightbits", True)

# absent kira.db -> loud ABSENT-UNGUARDED record, not a crash
empty = os.path.join(OUT, "empty_gen")
os.makedirs(empty, exist_ok=True)
g = kira.version_guard(empty)
check("guard_absent_recorded", g["status"] == "ABSENT-UNGUARDED", g)

# row-index integrity gate (checker lineage): to_system refuses a load
# where empty rows were dropped (lambda indices would misalign)
doctored = dict(sysd)
doctored["counters"] = dict(sysd["counters"], cancelled_to_empty=1)
try:
    kira.to_system(doctored, p)
    check("row_index_integrity_gate", False)
except AssertionError:
    check("row_index_integrity_gate", True)

finish("test_kira_adapter")

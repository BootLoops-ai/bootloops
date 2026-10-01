#!/usr/bin/env python3
"""Route B battery (adapters/user_system.py — the documented
user_defined_system input format).

Fixtures are authored by the test itself (toy scale, seconds): a
weight-notation directory holding one plain .kira file, one gzipped
.kira.gz file, and one decoy file that the documented extension filter
must ignore; and an integral-notation file with symbolic (d, s, m2)
coefficients.

Gates:
  U1  weight notation parses to EXACTLY the hand-computed F_p rows (file
      order sorted by name, blank-line equation splits, first-`*` term
      split, gz and plain text both read, decoy file ignored);
  U2  integral notation: deterministic Laporta ranks (t, r, s, sector,
      family, indices ascending), labels/col_of round-trip, stratum = t;
  U3  provenance gate: matching expect_counters passes, mismatching
      raises ProvenanceError;
  U4  refusals are loud: mixed notation, undeclared symbol, missing `*`,
      ill-formed integral, unknown forbid label, missing/empty source;
  U5  end-to-end Route B quickstart: load -> eliminate B2FT -> native +
      retrofit witnesses verify against an INDEPENDENT second load of the
      same files (the caller's-own-parse contract), and the closed rows
      match a dense in-test oracle solve;
  U6  master relations surface as leftover, never absorbed (a dependent
      equation planted in the fixture).
"""
import _common  # noqa: F401  MUST be first: sys.path bootstrap for the package

import gzip
import os

import ibplapper as lap
from ibplapper import receipt
from ibplapper.adapters import kira as kira_adapter
from ibplapper.adapters import user_system
from _common import check, finish, out_dir

OUT = out_dir("/tmp/ibplapper_test_user_system")
P = 1048573
D0, S0, M20 = 12347, 87651, 4242


def wtxt(path, text):
    with open(path, "w") as fh:
        fh.write(text)


# --------------------------------------------------------------- fixtures
WDIR = os.path.join(OUT, "sys_w")
os.makedirs(WDIR, exist_ok=True)
# a.kira: two equations (weight notation, one term per line, blank split)
wtxt(os.path.join(WDIR, "a.kira"),
     "5*(2)\n"
     "3*(-1)\n"
     "\n"
     "4*(d)\n"
     "3*(1)\n")
# b.kira.gz: one equation, gzipped; trailing equation without final blank
with gzip.open(os.path.join(WDIR, "b.kira.gz"), "wt") as fh:
    fh.write("5*(1)\n4*(1)\n2*(s)\n")
# decoy: must be ignored by the documented extension filter
wtxt(os.path.join(WDIR, "notes.txt"), "9*(1)\n")

IFILE = os.path.join(OUT, "sys_i.kira")
wtxt(IFILE,
     "T[1,1,0]*(d-2)\n"
     "T[1,0,0]*(s)\n"
     "T[0,1,0]*(-1)\n"
     "\n"
     "T[1,1,-1]*(1)\n"
     "T[1,1,0]*(2*m2)\n"
     "T[1,0,0]*(1)\n"
     "\n"
     # dependent equation: 2x the first -> must surface as rank deficiency
     "T[1,1,0]*(2*d-4)\n"
     "T[1,0,0]*(2*s)\n"
     "T[0,1,0]*(-2)\n")

# --------------------------------------------- U1: weight notation, exact
sysd_w = user_system.load_system(WDIR, P, {"d": D0, "s": S0})
expect_rows = [
    {5: 2 % P, 3: P - 1},           # a.kira eq 1
    {4: D0 % P, 3: 1},              # a.kira eq 2
    {5: 1, 4: 1, 2: S0 % P},        # b.kira.gz eq 1 (files sorted by name)
]
check("U1_weight_rows_exact", sysd_w["rows"] == expect_rows,
      f"rows={sysd_w['rows']}")
check("U1_counters", sysd_w["counters"]["eqs"] == 3
      and sysd_w["counters"]["terms"] == 7
      and sysd_w["counters"]["files"] == 2
      and sysd_w["counters"]["notation"] == "weight",
      sysd_w["counters"])
check("U1_order_is_weight",
      sysd_w["order"] == {c: c for c in (2, 3, 4, 5)})

# ------------------------------------------- U2: integral notation, ranks
sysd_i = user_system.load_system(IFILE, P, {"d": D0, "s": S0, "m2": M20})
# ascending (t, r, s, sector, fam, idx):
#   rank 0: T[1,0,0]  t1 r1 s0 sector 1
#   rank 1: T[0,1,0]  t1 r1 s0 sector 2
#   rank 2: T[1,1,0]  t2 r2 s0 sector 3
#   rank 3: T[1,1,-1] t2 r2 s1 sector 3
check("U2_ranks", sysd_i["col_of"] == {("T", (1, 0, 0)): 0,
                                       ("T", (0, 1, 0)): 1,
                                       ("T", (1, 1, 0)): 2,
                                       ("T", (1, 1, -1)): 3},
      sysd_i["col_of"])
check("U2_labels_roundtrip",
      all(sysd_i["col_of"][sysd_i["labels"][c]] == c
          for c in sysd_i["labels"]))
check("U2_strata_are_t",
      sysd_i["stratum_of"] == {0: 1, 1: 1, 2: 2, 3: 2})
check("U2_rows_exact", sysd_i["rows"] == [
    {2: (D0 - 2) % P, 0: S0 % P, 1: P - 1},
    {3: 1, 2: (2 * M20) % P, 0: 1},
    {2: (2 * D0 - 4) % P, 0: (2 * S0) % P, 1: P - 2}],
    sysd_i["rows"])

# --------------------------------------------------- U3: provenance gate
ok3 = True
try:
    user_system.load_system(IFILE, P, {"d": D0, "s": S0, "m2": M20},
                            expect_counters={"eqs": 3, "terms": 9})
except kira_adapter.ProvenanceError:
    ok3 = False
check("U3_expect_counters_pass", ok3)
try:
    user_system.load_system(IFILE, P, {"d": D0, "s": S0, "m2": M20},
                            expect_counters={"eqs": 4})
    check("U3_expect_counters_mismatch_raises", False, "no exception")
except kira_adapter.ProvenanceError:
    check("U3_expect_counters_mismatch_raises", True)

# -------------------------------------------------------- U4: loud refusals
def refuses(name, text, values=None, exc=user_system.UserSystemFormatError):
    path = os.path.join(OUT, f"bad_{name}.kira")
    wtxt(path, text)
    try:
        user_system.load_system(path, P, values or {"d": D0, "s": S0})
        return check(f"U4_{name}_refused", False, "no exception")
    except exc:
        return check(f"U4_{name}_refused", True)


refuses("mixed_notation", "5*(2)\nT[1,0]*(1)\n")
refuses("undeclared_symbol", "5*(2*x7)\n3*(1)\n")
refuses("missing_star", "5*(2)\n3\n")
refuses("illformed_integral", "T[1,q]*(1)\n")
refuses("index_count_drift", "T[1,0]*(1)\nT[1,0,0]*(1)\n")
refuses("empty_input", "\n\n")
try:
    user_system.to_system(sysd_i, P, forbid=[("T", (9, 9, 9))])
    check("U4_unknown_forbid_label_refused", False, "no exception")
except KeyError:
    check("U4_unknown_forbid_label_refused", True)
try:
    user_system.load_system(os.path.join(OUT, "no_such_source"), P,
                            {"d": D0})
    check("U4_missing_source_refused", False, "no exception")
except FileNotFoundError:
    check("U4_missing_source_refused", True)

# ------------------------------------- U5/U6: end-to-end Route B quickstart
system, sysd = user_system.load(
    IFILE, P, {"d": D0, "s": S0, "m2": M20},
    forbid=[("T", (1, 0, 0)), ("T", (0, 1, 0))],
    expect_counters={"eqs": 3, "terms": 9})
targets = [2, 3]
res = lap.eliminate(system, lap.Schedule(policy="B2FT"),
                    witnesses="native", witness_targets=targets,
                    witness_dir=os.path.join(OUT, "wit"))
check("U5_native_certified_verified",
      all(res.witnesses[t]["status"] == "CERTIFIED"
          and res.witnesses[t]["verify_pass"] for t in targets))
res_r = lap.eliminate(system, lap.Schedule(policy="B2FT"),
                      witnesses="retrofit", witness_targets=targets)
check("U5_retrofit_agrees",
      all(res_r.witnesses[t]["status"] == "CERTIFIED"
          and res_r.witnesses[t]["c"] == res.witnesses[t]["c"]
          for t in targets))

# independent second load = the caller's own parse (mandatory downstream)
system2, _ = user_system.load(IFILE, P, {"d": D0, "s": S0, "m2": M20},
                              forbid=[0, 1])
ok5 = all(receipt.verify([dict(r) for r in system2.rows],
                         res.witnesses[t]["witness_path"])[0]
          for t in targets)
check("U5_verify_against_independent_reparse", ok5)

# in-test oracle: hand-solve the 2x2 pivot block (mod P) in the subs
# convention  pivot + rest = 0.
#   eq0: a*c2 + b*c0 - c1 = 0        with a = d-2, b = s
#     -> c2 + (b/a)*c0 - (1/a)*c1 = 0
#   eq1: c3 + 2*m2*c2 + c0 = 0, substitute c2:
#     -> c3 + (1 - 2*m2*b/a)*c0 + (2*m2/a)*c1 = 0
a, b = (D0 - 2) % P, S0 % P
inv_a = pow(a, P - 2, P)
sub2 = {0: b * inv_a % P, 1: (P - inv_a) % P}
sub3 = {0: (1 - (2 * M20 * b % P) * inv_a) % P,
        1: (2 * M20 % P) * inv_a % P}
check("U5_closed_rows_match_oracle",
      res.subs[2] == {k: v for k, v in sub2.items() if v}
      and res.subs[3] == {k: v for k, v in sub3.items() if v},
      f"subs={res.subs}")
# the planted dependent equation (2x the first) must cancel to nothing:
# exactly two pivots, no spurious leftover
check("U6_dependent_row_cancels_cleanly",
      len(res.subs) == 2 and res.leftover == [],
      f"leftover={res.leftover}")

finish("test_user_system")

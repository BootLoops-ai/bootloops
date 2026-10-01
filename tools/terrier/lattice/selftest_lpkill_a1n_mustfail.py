#!/usr/bin/env python3
"""Battery LAT-LP (lattice/lp_kill_rank_agnostic.py): two replays against
reference receipts.
(1) A1^6 must-fail control: strict h=2 KILLED / relaxed FEASIBLE n0=36 —
    the output line must BYTE-match the reference control line.
(2) d4 genus kill replay from a reference Hecke genus listing of the rank-18
    case: the regenerated record must be JSON-equal to the reference d4
    record (1 class, all 18 h KILLED, exact Farkas certificates).
The three reference files (A1^6 control line, rank-18 Hecke genus listing,
d4 LP record; names in REF_FILES below) are not included in the package:
point TERRIER_R18_DIR at a directory holding them (required, no default).
Contract: lattice/LP_KILL_RANK_GENERALIZATION.md (also documents the Hecke
listing format the d4 mode reads)."""
import json, os, subprocess, sys, tempfile
HERE = os.path.dirname(os.path.abspath(__file__))
F_CTRL, F_HECKE, F_D4 = REF_FILES = ("control_a16_out.txt", "hecke_r18_out.txt", "lp_r18_d4.json")
R18 = os.environ.get("TERRIER_R18_DIR")
if not R18:
    sys.exit("SKIP[lpkill]: TERRIER_R18_DIR is unset — the LP-kill reference "
             "receipts (A1^6 control line, rank-18 Hecke genus listing, d4 "
             "record) are not included in the package; point TERRIER_R18_DIR "
             "at a directory holding them to run")
for f in REF_FILES:
    assert os.path.exists(os.path.join(R18, f)), f"reference file MISSING: {R18}/{f}"

# (1) A1^6 control — byte-match vs the reference line
r = subprocess.run([sys.executable, "lp_kill_rank_agnostic.py", "control"],
                   cwd=HERE, capture_output=True, text=True, timeout=600)
want_line = open(os.path.join(R18, F_CTRL)).read()
print(r.stdout, end="")
assert r.returncode == 0, "control rc != 0:\n" + r.stderr[-800:]
assert r.stdout == want_line, (f"A1^6 control line != reference\n got: {r.stdout!r}\n"
                               f"want: {want_line!r}")

# (2) d4 genus kill replay — JSON-equal vs the reference d4 record
dst = os.path.join(tempfile.mkdtemp(prefix="lpkill_d4_"), "lp_d4.json")
r2 = subprocess.run([sys.executable, "lp_kill_rank_agnostic.py", "d4",
                     os.path.join(R18, F_HECKE), dst],
                    cwd=HERE, capture_output=True, text=True, timeout=900)
print(r2.stdout, end="")
assert r2.returncode == 0, "d4 replay rc != 0:\n" + r2.stderr[-800:]
got = json.load(open(dst))
want = json.load(open(os.path.join(R18, F_D4)))
assert got == want, "d4 replay JSON != reference d4 record"
assert "d4|CLOSED-ALLROOTED|classes=1|" in r2.stdout, "d4 summary line missing"
print("SELFTEST lpkill_a1n_mustfail: ALL PASS "
      "(A1^6 byte-match + d4 JSON-equal vs reference)")

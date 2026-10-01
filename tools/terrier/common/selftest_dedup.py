#!/usr/bin/env python3
"""selftest_dedup.py — battery for common/dedup.py (known-point
pullback-orbit dedup; packaged data common/data/KNOWN_ORBIT.jsonl).
D1 packaged orbit intact: 9 rows, schema, T1-T4 identity rows
D2 orbit control battery:
   5 MUST-FIRE (the AESZ118 pullback candidate + T1-T4 themselves, exact
   verdicts) + 3 MUST-NOT-FIRE (census_h30 non-orbit points)
D3 newform backstop controls: 4 known labels fire with exact verdicts,
   unknown label (54.4.a.c) must NOT fire
D4 receipt law: explicit-path emit appends the ledger row; NO write
   without a receipt path; CLI rc 10 = rediscovery / 0 = no match
Resource law: ulimit -v 32505856, nice >= 5. Writes only volatile
out_selftest files in this directory."""
import json, os, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import dedup

FAIL = []


def gate(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name} {detail}", flush=True)
    if not ok:
        FAIL.append(name)


CONTROLS = [
    # (op, spec, expect_fire, expect_verdict, note)
    ("AESZ118", "z=-1/32", True, "REDISCOVERY-OF-T3",
     "the AESZ118 candidate: card-documented pullback of T3"),
    ("AESZ34", "z=-1/7", True, "REDISCOVERY-OF-T1", "T1 itself"),
    ("AESZ34", "minpoly=1,-66,1", True, "REDISCOVERY-OF-T2",
     "T2 pair itself (quadratic minpoly spec)"),
    ("AESZ22", "z=-1", True, "REDISCOVERY-OF-T3", "T3 itself"),
    ("AESZ17", "z=-1", True, "REDISCOVERY-OF-T4", "T4 itself"),
    ("AESZ34", "minpoly=3,-1,-1", False, None, "census non-orbit d13"),
    ("AESZ22", "minpoly=11,11,4", False, None, "census non-orbit d-55"),
    ("AESZ118", "minpoly=18,-21,-23", False, None,
     "census non-orbit d2097"),
]

NEWFORMS = [("14.4.a.a", "REDISCOVERY-OF-T1"),
            ("33.4.a.b", "REDISCOVERY-OF-T3"),
            ("34.4.b.a", "REDISCOVERY-OF-T2"),
            ("14.4.a.b", "REDISCOVERY-OF-T4"),
            ("54.4.a.c", None)]


def run_D1():
    rows = dedup._load_orbit()
    idn = {(r["operator_key"], r["match_key"]) for r in rows
           if r["relation"] == "identity"}
    want = {("AESZ34", "rat:-1/7"), ("AESZ34", "alg:1,-66,1"),
            ("AESZ22", "rat:-1"), ("AESZ17", "rat:-1")}
    schema = all({"operator_key", "z_image", "source_point", "relation",
                  "provenance", "match_key"} <= set(r) for r in rows)
    gate("D1 packaged KNOWN_ORBIT intact", len(rows) == 9 and schema
         and want <= idn, f"{len(rows)} rows, identity {len(idn)}")


def run_D2():
    ok, det = True, []
    for op, spec, want_fire, want_v, _note in CONTROLS:
        v, hits = dedup.orbit_match(op, spec)
        good = ((v is not None) == want_fire) \
            and (not want_fire or v == want_v) \
            and (bool(hits) == want_fire)
        ok &= good
        if not good:
            det.append(f"{op} {spec} -> {v}")
    gate("D2 orbit controls 5 must-fire + 3 must-not-fire", ok,
         "; ".join(det) or f"{len(CONTROLS)} controls")


def run_D3():
    ok = all(dedup.newform_backstop(lab) == want
             for lab, want in NEWFORMS)
    gate("D3 newform backstop controls (incl. unknown must-not-fire)",
         ok, f"{len(NEWFORMS)} labels")


def run_D4():
    rp = os.path.join(HERE, "out_selftest_dedup_receipt.jsonl")
    if os.path.exists(rp):
        os.remove(rp)
    v = dedup.check("AESZ118", "z=-1/32")          # no receipt path
    ok_nowrite = v == "REDISCOVERY-OF-T3" and not os.path.exists(rp)
    v2 = dedup.check("AESZ118", "z=-1/32", receipt=rp)
    row = json.loads(open(rp).readline())
    gate("D4a explicit-receipt law", ok_nowrite
         and v2 == "REDISCOVERY-OF-T3"
         and row["verdict"] == "REDISCOVERY-OF-T3"
         and row["stage"] == "dedup-pre-escalation",
         "no write without path; ledger row on explicit path")
    r1 = subprocess.run([sys.executable, "dedup.py", "check", "AESZ34",
                         "z=-1/7"], cwd=HERE, capture_output=True,
                        text=True)
    r2 = subprocess.run([sys.executable, "dedup.py", "check", "AESZ34",
                         "z=1/2"], cwd=HERE, capture_output=True,
                        text=True)
    r3 = subprocess.run([sys.executable, "dedup.py", "newform",
                         "14.4.a.a"], cwd=HERE, capture_output=True,
                        text=True)
    gate("D4b CLI rc semantics (10 rediscovery / 0 no-match)",
         r1.returncode == 10 and "REDISCOVERY-OF-T1" in r1.stdout
         and r2.returncode == 0 and "NO-MATCH" in r2.stdout
         and r3.returncode == 10)
    os.remove(rp)


if __name__ == "__main__":
    T0 = time.time()
    run_D1()
    run_D2()
    run_D3()
    run_D4()
    n = 5
    print(f"[selftest_dedup] {n - len(FAIL)}/{n} passed "
          f"({time.time()-T0:.0f}s wall)"
          + ("" if not FAIL else f"  FAILED: {FAIL}"))
    sys.exit(1 if FAIL else 0)

#!/usr/bin/env python3
"""amflow-kit battery entry — the public legs (no built amflow_cli needed):

  1. amflow_output_lint.py on the shipped vacuum fixture   -> must PASS
  2. compare_amflow_json.py ref-vs-ref on the same fixture -> must PASS
     (the goal from the fixture's own job config, in_vac2Bprobe.json, via
     --goal-from-config: the output form carries no goal key and the
     comparator refuses to run without a derivable goal)
  3. amflow_smoke.sh with no usable binary                 -> must fail closed
     with a named FATAL (the full smoke needs a built amflow_cli)
  4. amflow_kit.keypred --selftest                         -> 16/16, 0 fail
     (fixture, mutation and salt battery of the IBP-cache key predictor;
     its scratch copies go under $TMPDIR)
  5. `import amflow_kit` exposes keypred and memfence, and memfence's policy
     arithmetic on a planted no-leaf fixture reads 48 -> 120 GiB (2.5x)

Exit 0 only if all five legs behave as documented. The memfence process legs
and the comparator's rule battery are pytest: `python3 -m pytest tests -q -rs
-p no:cacheprovider` from this directory (see GUIDE.md, Battery).
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
FIX = os.environ.get("AMFLOW_SMOKE_FIX",
                     os.path.join(HERE, "..", "fixtures", "amflow_smoke"))
OUT = os.path.join(FIX, "out_vac2Bprobe.json")
IN = os.path.join(FIX, "in_vac2Bprobe.json")     # the job config of that output (carries the goal)
fails = []


def leg(name, argv, want_fail=False, want_text=None):
    r = subprocess.run(argv, capture_output=True, text=True, cwd=HERE)
    out = r.stdout + r.stderr
    ok = (r.returncode != 0) if want_fail else (r.returncode == 0)
    if want_text and want_text not in out:
        ok = False
    print(f"[{'PASS' if ok else 'FAIL'}] {name} (rc={r.returncode})")
    if not ok:
        print(out[-800:])
        fails.append(name)


leg("lint of the shipped vacuum fixture",
    [sys.executable, os.path.join(HERE, "amflow_output_lint.py"), OUT,
     "--class", "vacuum"])
leg("compare ref-vs-ref",
    [sys.executable, os.path.join(HERE, "compare_amflow_json.py"), OUT, OUT,
     "--min-digits", "20", "--goal-from-config", IN], want_text="OVERALL: PASS")
leg("smoke fails closed without a vendor binary (named FATAL)",
    ["bash", os.path.join(HERE, "amflow_smoke.sh"),
     os.path.join(HERE, "no_such_amflow_cli")],
    want_fail=True, want_text="FATAL")
leg("keypred selftest (fixture + mutation + salt battery)",
    [sys.executable, "-m", "amflow_kit.keypred", "--selftest"],
    want_text=" 0 fail")

# leg 5: the package imports and memfence's fuse arithmetic on a planted
# fixture (no /proc or cgroup reads: the leaf record is handed in)
name = "import amflow_kit; memfence fuse arithmetic 48 -> 120 GiB on a planted no-leaf fixture"
try:
    sys.path.insert(0, HERE)
    import amflow_kit
    mf = amflow_kit.memfence
    assert amflow_kit.keypred.numd_salt("-3/7") == "AMFLOW_DIFFEQ_NUMD=-3/7"
    leaf = {"cgroup_path": "/planted/no_leaf.scope", "memory_max_raw": "max",
            "memory_max_bytes": None, "ancestor_min_bytes": None}
    pol = mf.fuse_policy("48", leaf=leaf)
    assert pol["mode"] == "fuse", pol["mode"]
    assert pol["fuse_bytes"] == 120 * mf.GIB, pol["fuse_bytes"]
    assert pol["ulimit_v_kb"] == 125829120, pol["ulimit_v_kb"]
    assert pol["export_cap_gb"] == "120", pol["export_cap_gb"]
    assert "address-space fuse" in pol["words"]
    leaf["memory_max_bytes"], leaf["memory_max_raw"] = 68719476736, "68719476736"
    pol = mf.fuse_policy("48", leaf=leaf)
    assert pol["mode"] == "leaf" and pol["ulimit_v_kb"] is None
    assert pol["export_cap_gb"] == mf.NO_FUSE_EXPORT == "0"
    print(f"[PASS] {name}")
except Exception as e:  # noqa: BLE001 — any failure is a named FAIL
    print(f"[FAIL] {name}: {type(e).__name__}: {e}")
    fails.append(name)

if fails:
    print(f"OVERALL: FAIL {fails}")
    sys.exit(1)
print("OVERALL: PASS (public legs; the full smoke needs a built amflow_cli; "
      "the memfence process legs and the comparator rules run under pytest tests/)")

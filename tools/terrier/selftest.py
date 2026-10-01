#!/usr/bin/env python3
"""terrier suite selftest — sha-pinned manifest FIRST, then
replay every registered module's battery in this tree from its pinned inputs.
Missing or changed pinned artifacts are NAMED loudly and the battery replays
refuse to run (never a silent skip).

    ulimit -v 32505856; nice -n 5 python3 selftest.py [--only ID[,ID..]]

Sections:
  M. manifest  — sha256 of every artifact in regression_manifest.json
     ("artifacts" = the shipped in-tree fixtures; every pin verified FIRST).
  B. batteries — every battery row in the manifest runs in its own dir under
     the resource caps below (python: ulimit -v 32505856, nice 5; julia: +1.10
     -t 1, nice 5). Missing battery file = NAMED FAIL. Batteries whose inputs
     are reference data not included in the package (or a Julia env) are a
     loud NAMED SKIP naming the env var to set; set it and they run for
     real. --only filters
     by battery id but section M always runs in full.
  X. coverage  — every selftest_*.{py,jl} on disk must be registered in the
     manifest (no silently-added batteries) and the wing-regression map
     (DKMM balls, lattice 18/18, UU cell, shard replay, planted-error rc,
     mutation-control trio) must point at registered batteries.
Exit nonzero on any FAIL. Overview: README.md; module index: MODULE_MAP.md.
"""
import argparse, hashlib, json, os, subprocess, sys, time
ROOT = os.path.dirname(os.path.abspath(__file__))
MANIFEST = os.path.join(ROOT, "regression_manifest.json")
# julia project root comes from the environment (no machine-local default)
JULIA_PROJECT = os.environ.get("TERRIER_JULIA_PROJECT", "")
WINGS = ["periods", "common", "lattice", "census"]
RESULTS = []

# Batteries whose inputs are reference data sets not included in the
# package. Unset env var -> loud NAMED SKIP (each row names exactly what to
# set); point the var at the data and the battery runs for real.
# batteries that import an engine shipped as a SIBLING package of the
# repository (tools/<name>/), looked up the way the importing module does:
# the sibling directory first, then PYTHONPATH.
NEEDS_MODULE = {
    "opderive": ("annihilator", "tools/annihilator"),
}


def _module_available(mod, _label):
    import importlib.util
    sib = os.path.abspath(os.path.join(ROOT, "..", mod))
    sys.path.insert(0, sib)
    try:
        return importlib.util.find_spec(mod) is not None
    finally:
        sys.path.remove(sib)


NEEDS_DATA = {
    "geffseries": "TERRIER_KKLT_BANK",
    "fluxcurves": "TERRIER_KKLT_BANK",
    "opderive": "TERRIER_KKLT_BANK",
    "transport": "TERRIER_KKLT_BANK",
    "verdicts": "TERRIER_KKLT_BANK",
    "envelope": "TERRIER_KKLT_BANK",
    "pfaffian": "TERRIER_PFAFFIAN_BANK",
    "sweepchassis_L2b": "TERRIER_SWEEP_SHARD_RECEIPT",
    "foldhnf": "TERRIER_G4_DIR",
    "distcert": "TERRIER_PLANTS_DIR",
    "controlsampler": "TERRIER_SEEDS_JSON",
    "summation": "TERRIER_BOXOPS_BANK",
    "lpkill": "TERRIER_R18_DIR",
    "genusenum": "TERRIER_GENUS_CONTROL_RECEIPT",
}


def record(sec, name, ok, detail=""):
    """ok True = PASS, False = FAIL, None = named SKIP (input not included)."""
    RESULTS.append((sec, name, ok, detail))
    mark = ("SKIP (named)" if ok is None
            else "PASS" if ok else "FAIL <<<<<<<<<<")
    print(f"  [{mark}] {sec} :: {name}" + (f" -- {detail}" if detail else ""))


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def section_manifest(man):
    """Verify the shipped fixture pins. Every missing/changed file gets its
    own NAMED FAIL row. Returns True iff every pin is intact."""
    missing, changed = [], []
    for rel, want in sorted(man["artifacts"].items()):
        full = os.path.join(ROOT, rel)
        if not os.path.exists(full):
            missing.append(rel)
        elif sha256(full) != want:
            changed.append(rel)
    ok = not missing and not changed
    record("M.manifest", f"artifacts: {len(man['artifacts'])} sha256-pinned",
           ok, "" if ok else f"{len(missing)} missing, {len(changed)} changed")
    for rel in missing:
        record("M.manifest", f"MISSING artifact: {rel}", False)
    for rel in changed:
        record("M.manifest", f"CHANGED artifact: {rel}", False, "pin mismatch")
    return ok


def run_battery(b):
    """One battery, in its own directory, under the resource caps."""
    path = os.path.join(ROOT, b["path"])
    sec, name = f"B.{b['wing']}", f"{b['id']} [{b['path']}]"
    if not os.path.exists(path):
        record(sec, name, False, "BATTERY FILE MISSING")
        return
    d, f = os.path.dirname(path), os.path.basename(path)
    if b["runner"] == "julia" and not JULIA_PROJECT:
        record(sec, name, None, "needs a Julia Oscar/Hecke project env "
               "-- set TERRIER_JULIA_PROJECT to run")
        return
    need = NEEDS_DATA.get(b["id"])
    if need and not os.environ.get(need):
        record(sec, name, None,
               f"needs reference data not included in the package -- "
               f"set {need} to run")
        return
    ext = NEEDS_MODULE.get(b["id"])
    if ext and not _module_available(*ext):
        record(sec, name, None,
               f"needs the external engine {ext[1]} (a sibling package in "
               f"the BootLoops repository) -- run inside the repository "
               f"checkout or put it on PYTHONPATH")
        return
    if b["runner"] == "julia":
        cmd = ["nice", "-n", "5", "julia", "+1.10", "-t", "1",
               f"--project={JULIA_PROJECT}", f]
    else:
        cmd = ["bash", "-c",
               f"ulimit -v 32505856; exec nice -n 5 {sys.executable} {f}"]
    t0 = time.time()
    try:
        r = subprocess.run(cmd, cwd=d, capture_output=True, text=True,
                           timeout=b.get("timeout_s", 1800))
        ok, tail = r.returncode == 0, (r.stdout + r.stderr)[-1200:]
    except subprocess.TimeoutExpired:
        ok, tail = False, f"TIMEOUT > {b.get('timeout_s', 1800)}s"
    wall = time.time() - t0
    record(sec, name, ok, f"wall {wall:.1f}s" if ok
           else f"wall {wall:.1f}s\n{tail}")


def section_coverage(man):
    """No silent additions/omissions: on-disk selftest files == registered."""
    on_disk = set()
    for wing in WINGS:
        for base, _dirs, files in os.walk(os.path.join(ROOT, wing)):
            for f in files:
                if f.startswith("selftest_") and f.endswith((".py", ".jl")):
                    on_disk.add(os.path.relpath(os.path.join(base, f), ROOT))
    reg = {b["path"] for b in man["batteries"]}
    record("X.coverage", "every on-disk selftest registered",
           on_disk <= reg, ", ".join(sorted(on_disk - reg)))
    record("X.coverage", "every registered battery on disk",
           reg <= on_disk, ", ".join(sorted(reg - on_disk)))
    ids = {b["id"] for b in man["batteries"]}
    bad = {k: v for k, v in man["wing_regressions"].items()
           if not set(v["batteries"]) <= ids}
    record("X.coverage", "wing-regression map -> registered batteries",
           not bad, ", ".join(bad))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default=None,
                    help="comma-separated battery ids (M + X still run full)")
    args = ap.parse_args()
    t0 = time.time()
    if not os.path.exists(MANIFEST):
        record("M.manifest", "regression_manifest.json present", False,
               f"{MANIFEST} missing -- rebuild via build_manifest.py")
        man = None
    else:
        man = json.load(open(MANIFEST))
    pins_ok = section_manifest(man) if man else False
    if man and pins_ok:
        only = set(args.only.split(",")) if args.only else None
        for b in man["batteries"]:
            if only is None or b["id"] in only:
                run_battery(b)
        section_coverage(man)
    elif man:
        record("M.manifest", "battery replays REFUSED: pins not intact "
               "(restore the files or rebuild the manifest first)", False,
               "no silent skips -- every battery above is unexercised")
    wall = time.time() - t0
    nfail = sum(1 for _s, _n, ok, _d in RESULTS if ok is False)
    nskip = sum(1 for _s, _n, ok, _d in RESULTS if ok is None)
    print("\n" + "=" * 74)
    print(f"{'section':<12}{'test':<50}{'result'}")
    print("-" * 74)
    for s, n, ok, _d in RESULTS:
        print(f"{s:<12}{n[:49]:<50}" + ("SKIP (named)" if ok is None
                                        else "PASS" if ok else "FAIL  <<<<<<"))
    print("-" * 74)
    print(f"{len(RESULTS) - nfail - nskip}/{len(RESULTS)} passed"
          + (f", {nskip} named skip(s)" if nskip else "")
          + f" in {wall:.1f}s wall"
          + ("" if nfail == 0 else f"  --  {nfail} FAILURE(S)"))
    sys.exit(1 if nfail else 0)


if __name__ == "__main__":
    main()

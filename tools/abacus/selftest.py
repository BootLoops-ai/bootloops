#!/usr/bin/env python3
"""abacus selftest — the shipped battery in one command.

Runs the documented battery from a scratch working directory:
  1. run_s0.py -> expects "S0-PASS 7/7"
  2. run_s3.py -> expects "BATTERY-PASS" (the full harness, 22/22 checks)

run_s3's theta/period bridge legs (the S2 replay and the B1 receipts) need
the Eichler julia project: when ABACUS_EICHLER_PROJECT is unset, this
wrapper points it at the copy shipped in this repository
(upgrades/Eichler.jl), per GUIDE.md, and readies it (instantiate + load)
before entering the battery — a named, fail-closed refusal when it cannot.
Stage receipts land in the scratch directory. The harness writes its result
receipts beside itself, as designed — so the battery runs a byte-identical
scratch COPY of the package tree and the receipts land in that copy; the
shipped, sha-pinned tree is never written.
Exit 0 only if both stages pass.
"""
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
# ABACUS_JULIA overrides the julia executable (an absolute path or a name
# looked up on PATH); default: `julia` on PATH — mirrors abcount_s1.JULIA.
JULIA = os.environ.get("ABACUS_JULIA", "julia")


def ensure_eichler_ready(env):
    """The S3 battery's theta/period bridges run julia under the Eichler
    project (ABACUS_EICHLER_PROJECT).  A clone is complete, but the julia
    package depot may not be: on a machine that has never installed the
    project's pinned packages, every bridge call dies with "Package Arblib
    ... not installed" — the battery's first bridge call (the S2 replay)
    then banks an opaque STOP.  Settle that here, once, loudly: instantiate
    the project (a no-op on a ready depot; downloads and precompiles on
    first use) and load Arblib so the one-time precompilation cost lands
    here rather than inside the battery's timed calls.  Any failure is a
    named, fail-closed refusal — the battery is never entered."""
    proj = env.get("ABACUS_EICHLER_PROJECT", "")
    if not proj:
        print("ABACUS-REFUSAL: ABACUS_EICHLER_PROJECT is unset and the "
              "in-repo copy (upgrades/Eichler.jl) was not found — the S3 "
              "battery's julia bridges cannot run. Point the variable at an "
              "Eichler.jl project root (see GUIDE.md Dependencies).")
        return 1
    if not os.path.isdir(proj):
        print("ABACUS-REFUSAL: ABACUS_EICHLER_PROJECT points at %s, which "
              "is not a directory — the S3 battery's julia bridges cannot "
              "run (see GUIDE.md Dependencies)." % proj)
        return 1
    if shutil.which(JULIA) is None:
        print("ABACUS-REFUSAL: julia executable %r does not resolve — set "
              "ABACUS_JULIA to a julia executable (an absolute path or a "
              "name on PATH) or put julia on PATH (abcount_s1.JULIA "
              "resolves the same way); required for the S3 battery's "
              "theta/period bridges." % JULIA)
        return 1
    print("[eichler] readying the julia project (first run may download and "
          "precompile packages — this can take minutes) ...", flush=True)
    try:
        r = subprocess.run(
            [JULIA, "--project=%s" % proj, "-e",
             "import Pkg; Pkg.instantiate(); using Arblib"],
            env=env, capture_output=True, text=True, timeout=1800)
    except (subprocess.TimeoutExpired, OSError) as e:
        print("ABACUS-REFUSAL: the Eichler julia project could not be "
              "readied (%s). The S3 battery's bridges would fail "
              "identically; not entering the battery." % e)
        return 1
    if r.returncode != 0:
        print("ABACUS-REFUSAL: the Eichler julia project failed to "
              "instantiate/load (julia rc=%d). The S3 battery's bridges "
              "would fail identically; not entering the battery. "
              "julia said:" % r.returncode)
        for line in (r.stderr or r.stdout).strip().splitlines()[-12:]:
            print("  " + line)
        return 1
    return 0


def main():
    try:
        import cypari2                              # noqa: F401
    except ImportError:
        print("ABACUS-REFUSAL: python module `cypari2` is not installed — "
              "the S0 sanity run and the battery's PARI cross-paths import "
              "it (pip install cypari2; see GUIDE.md Dependencies).")
        return 1
    env = dict(os.environ)
    if not env.get("ABACUS_EICHLER_PROJECT"):
        cand = os.path.join(REPO, "upgrades", "Eichler.jl")
        if os.path.isdir(cand):
            env["ABACUS_EICHLER_PROJECT"] = cand
    rc = ensure_eichler_ready(env)
    if rc:
        return rc
    scratch = tempfile.mkdtemp(prefix="abacus_selftest_")
    # A registered battery never writes into the tracked tree.  The harness
    # writes its receipts beside itself by design (S3_RUN_OUTPUT.json,
    # S2_RUN_OUTPUT.json replay publish, S2_RUN_OUTPUT_s2stage.json,
    # build/S3_RESULT.json, battery/BATTERY_RESULT.json, theta bridge txt) —
    # so this wrapper runs a byte-identical COPY of the package tree under the
    # scratch directory (the same isolation idiom run_s3 itself uses for its
    # S2 replay): every receipt lands in the copy, the shipped sha-pinned tree
    # is never written.
    # The sha pins verify the copied bytes, which copytree preserves.
    pkg = os.path.join(scratch, "pkg")
    shutil.copytree(HERE, pkg, ignore=shutil.ignore_patterns("__pycache__"))
    for script, needle in (("run_s0.py", "S0-PASS 7/7"),
                           ("run_s3.py", "BATTERY-PASS")):
        r = subprocess.run([sys.executable, os.path.join(pkg, script)],
                           cwd=scratch, env=env, capture_output=True,
                           text=True, timeout=1800)
        out = r.stdout + r.stderr
        for line in out.strip().splitlines()[-2:]:
            print(f"[{script}] {line}")
        if r.returncode != 0 or needle not in r.stdout:
            print(f"FAIL: {script} rc={r.returncode} "
                  f"(expected {needle!r} in its output)")
            return 1
    print("abacus selftest PASS (S0 7/7 + S3 full battery)")
    return 0


if __name__ == "__main__":
    sys.exit(main())

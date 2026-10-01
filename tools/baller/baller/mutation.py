"""baller.mutation — the mutation-tested-gate harness: one reusable component
for the mutation controls that tools otherwise hand-write per test.

Law: mutants are built on SCRATCH
COPIES ONLY — the harness refuses to mutate a file in place. A gate is
mutation-tested iff it PASSES on the pristine bytes and FAILS on EVERY
mutant; a mutant that slips through means the gate is vacuous there.
"""
import os
import shutil
import subprocess
import sys
import tempfile


class MutationHarnessError(RuntimeError):
    pass


def run_mutation_gate(source_path, mutants, gate_cmd, scratch_dir=None,
                      timeout=300):
    """source_path: file the gate exercises. mutants: list of (old, new)
    unique-substring substitutions. gate_cmd: argv list; every occurrence of
    the literal token '{FILE}' is replaced by the path under test; rc==0 =
    gate passes. The gate must PASS on a pristine scratch copy and FAIL on
    every mutant scratch copy. Returns the report dict; raises
    MutationHarnessError on a vacuous gate (or a broken pristine)."""
    if not any("{FILE}" in a for a in gate_cmd):
        raise MutationHarnessError(
            "gate_cmd contains no '{FILE}' placeholder — every run would "
            "silently re-test the ORIGINAL file, certifying nothing")
    src = open(source_path, "rb").read().decode()
    own = scratch_dir is None
    scratch_dir = scratch_dir or tempfile.mkdtemp(prefix="baller_mut_")
    os.makedirs(scratch_dir, exist_ok=True)
    if os.path.realpath(scratch_dir) == os.path.dirname(os.path.realpath(source_path)):
        raise MutationHarnessError(
            "scratch_dir equals the source's own directory — the in-place "
            "mutation footgun; refusing (scratch copies ONLY)")

    def run_on(text, tag):
        p = os.path.join(scratch_dir, f"{tag}_{os.path.basename(source_path)}")
        if os.path.lexists(p):
            # the scratch name must not pre-exist: a pre-planted SYMLINK
            # there would pass the write THROUGH to its target (clobbering
            # the source the harness promises never to touch) — remove the
            # entry first
            os.remove(p)
        open(p, "w").write(text)
        cmd = [a.replace("{FILE}", p) for a in gate_cmd]
        try:
            r = subprocess.run(cmd, capture_output=True, text=True,
                               timeout=timeout)
        except subprocess.TimeoutExpired as e:
            raise MutationHarnessError(
                f"gate TIMED OUT ({timeout}s) on {tag} — a hanging gate "
                f"certifies nothing; fix the gate or raise timeout") from e
        return r.returncode, (r.stdout + r.stderr)[-400:]

    report = {"source": source_path, "mutants": [], "pristine_rc": None}
    try:
        rc, out = run_on(src, "pristine")
        report["pristine_rc"] = rc
        if rc != 0:
            raise MutationHarnessError(
                f"gate FAILS on pristine bytes (rc={rc}) — fix the gate or the "
                f"source before mutation-testing: {out}")
        caught = 0
        for i, (old, new) in enumerate(mutants):
            if old == new:
                raise MutationHarnessError(
                    f"mutant {i}: old == new (no-op mutation certifies nothing)")
            if src.count(old) != 1:
                raise MutationHarnessError(
                    f"mutant {i}: anchor {old!r} occurs {src.count(old)}x in "
                    f"source (must be unique — ambiguous mutation)")
            mutated = src.replace(old, new, 1)
            if source_path.endswith(".py"):
                # a syntax-breaking mutant makes ANY gate "fail" for
                # the wrong reason — the harness then certifies vacuous gates.
                # Mutants must be SEMANTIC: the mutated file must still compile.
                try:
                    compile(mutated, source_path, "exec")
                except SyntaxError as e:
                    raise MutationHarnessError(
                        f"mutant {i} breaks syntax ({e.msg} at line {e.lineno})"
                        f" — not a semantic mutation; a compile-only gate "
                        f"would 'catch' it while asserting nothing") from e
                # a compile-clean but LOAD-FATAL mutant (NameError at
                # module top level) lets an import-only gate 'catch' it while
                # asserting nothing — the mutated module must still IMPORT.
                probe = os.path.join(scratch_dir,
                                     f"_loadcheck_{os.path.basename(source_path)}")
                open(probe, "w").write(mutated)
                lc = subprocess.run(
                    [sys.executable, "-B", "-c",
                     "import importlib.util, sys\n"
                     "spec = importlib.util.spec_from_file_location('m', sys.argv[1])\n"
                     "m = importlib.util.module_from_spec(spec)\n"
                     "spec.loader.exec_module(m)\n", probe],
                    capture_output=True, text=True, timeout=min(timeout, 60))
                if lc.returncode != 0:
                    raise MutationHarnessError(
                        f"mutant {i} is LOAD-FATAL, not semantic "
                        f"({lc.stderr.strip().splitlines()[-1][:100] if lc.stderr.strip() else 'load failed'})"
                        f" — an import-only gate would 'catch' it while "
                        f"asserting nothing; supply value-changing mutants")
            rc, out = run_on(mutated, f"mut{i}")
            ok = rc != 0
            caught += ok
            report["mutants"].append(
                {"i": i, "old": old[:60], "new": new[:60], "caught": ok, "rc": rc})
        report["caught"] = f"{caught}/{len(mutants)}"
        if caught != len(mutants):
            missed = [m["i"] for m in report["mutants"] if not m["caught"]]
            raise MutationHarnessError(
                f"VACUOUS GATE: mutants {missed} not caught ({caught}/"
                f"{len(mutants)}); report={report}")
        return report
    finally:
        if own:
            shutil.rmtree(scratch_dir, ignore_errors=True)

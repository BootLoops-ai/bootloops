"""Planted controls for the cgroup.procs write discipline in
runner.build_cgroup_script / runner.run_phase / cli --cgroup parsing.

No sudo anywhere. Every test but the last uses a temporary directory as its
"leaf": the child is a plain `sleep`/`bash`, the shell path guard runs on a
script that exits before its first sudo, and the SIGKILL / move-to-parent
operations are recorded by mocks through runner.run_phase's launch_ops
injection. The last test (test_positive_control_live_leaf_size_zero_with_N_pids)
is the one LIVE control: it makes a real cgroup v2 leaf under the calling
process's OWN cgroup, moves only children it spawned itself into it, and
tears it down in `finally`; where the kernel refuses the mkdir or the write it
SKIPS by name with the errno (never a silent pass — under plain pytest the
skip is visible with -rs and the exit status stays 0; `python3
test_cgroup_law.py` exits 2 on it).
"""
import contextlib
import errno
import io
import json
import os
import shlex
import subprocess
import sys
import tempfile
import time

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from seedling import cli, runner  # noqa: E402

BAD_SPEC = ",12G,400000 100000"           # empty NAME -> used to mean the ROOT
CFG = {"name": "testcap", "mem_max": "12G", "cpu_max": "400000 100000"}


def _toy(td):
    cfg = os.path.join(td, "config")
    os.makedirs(cfg)
    with open(os.path.join(cfg, "integralfamilies.yaml"), "w") as fh:
        fh.write("integralfamilies:\n - name: \"toy\"\n"
                 "   loop_momenta: [l]\n   top_level_sectors: [3]\n"
                 "   propagators:\n    - [\"l^2\", 0]\n"
                 "    - [\"(l+p1)^2\", 0]\n")
    tpath = os.path.join(td, "targets")
    with open(tpath, "w") as fh:
        fh.write("toy[1,1]\n")
    nts = os.path.join(td, "nts")
    with open(nts, "w") as fh:
        fh.write("3 2\n1 1\n")
    return cfg, tpath, nts


def _raises(fn, *args):
    try:
        fn(*args)
    except ValueError:
        return True
    return False


def _cli_refused(argv):
    """Run the CLI with runner.run_phase booby-trapped; return (rc, stderr)."""
    orig = runner.run_phase

    def _trap(*a, **k):
        raise AssertionError("run_phase reached: a script would have been built")

    runner.run_phase = _trap
    err = io.StringIO()
    try:
        with contextlib.redirect_stderr(err), contextlib.redirect_stdout(io.StringIO()):
            rc = cli.main(argv)
    finally:
        runner.run_phase = orig
    return rc, err.getvalue()


def test_validate_cgroup_name():
    for good in ("testcap", "seed.1_a-b", "A", "x..y", "-"):
        assert runner.validate_cgroup_name(good) == good
    for bad in ("", ".", "..", "...", "../x", "a/b", "/a", "a b", "a;b",
                "a,b", "$CG", "a\n", None):
        assert _raises(runner.validate_cgroup_name, bad), repr(bad)
    # build_cgroup_script refuses the same names before emitting anything
    assert _raises(runner.build_cgroup_script, "", "1G", "max", "true", "/")
    assert _raises(runner.build_cgroup_script, "..", "1G", "max", "true", "/")


def test_cgroup_spec_fields():
    cfg, why = cli._cgroup_spec("leaf1,12G,400000 100000")
    assert why is None and cfg == {"name": "leaf1", "mem_max": "12G",
                                   "cpu_max": "400000 100000"}
    assert cli._cgroup_spec("leaf1,max,max")[1] is None
    for spec in (BAD_SPEC, "leaf1,,max", "leaf1,12G,", "leaf1,12G",
                 "leaf1,12G;rm,max", "leaf1,12G,max; rm -rf x", "../x,12G,max",
                 "..,12G,max", ",,"):
        cfg, why = cli._cgroup_spec(spec)
        assert cfg is None and why.startswith("--cgroup "), spec


def test_control_a_empty_name_rc2_nothing_built():
    """--cgroup ',12G,400000 100000' -> rc 2, one-line reason, NOTHING built:
    no run dir, no jobs.yaml, no ledger, run_phase never reached (pin AND run)."""
    with tempfile.TemporaryDirectory() as td:
        cfg, tpath, nts = _toy(td)
        out = os.path.join(td, "run")
        rc, err = _cli_refused(["pin", "--config", cfg, "--targets", tpath,
                                "--out", out, "--execute", "--cgroup", BAD_SPEC,
                                "--kira-cmd", "true"])
        assert rc == 2 and not os.path.exists(out)
        assert err.count("\n") == 1 and "refusing" in err and "empty" in err
        rc, err = _cli_refused(["run", "--config", cfg, "--targets", tpath,
                                "--out", out, "--nts", nts, "--margins", "const",
                                "--execute", "--cgroup", BAD_SPEC])
        assert rc == 2 and not os.path.exists(out)
        assert err.count("\n") == 1 and "refusing" in err and "empty" in err
        # a bad spec is refused at the rehearsal too (no --execute)
        rc, err = _cli_refused(["pin", "--config", cfg, "--targets", tpath,
                                "--out", out, "--cgroup", BAD_SPEC])
        assert rc == 2 and not os.path.exists(out)


def test_control_b_dotdot_name_rc2():
    """name '../x' -> rc 2 (would have written under /sys/fs/x)."""
    with tempfile.TemporaryDirectory() as td:
        cfg, tpath, nts = _toy(td)
        out = os.path.join(td, "run")
        for spec in ("../x,12G,400000 100000", "..,12G,max", "a/b,12G,max"):
            rc, err = _cli_refused(["pin", "--config", cfg, "--targets", tpath,
                                    "--out", out, "--execute", "--cgroup", spec,
                                    "--kira-cmd", "true"])
            assert rc == 2 and not os.path.exists(out), spec
            assert err.count("\n") == 1 and "refusing" in err
            rc, err = _cli_refused(["run", "--config", cfg, "--targets", tpath,
                                    "--out", out, "--nts", nts, "--margins", "const",
                                    "--execute", "--cgroup", spec])
            assert rc == 2 and not os.path.exists(out), spec


def test_control_c_dry_run_carries_guard_tokens():
    """The dry-run script text carries the guard tokens, the path guard
    precedes the first write, the self-check follows the cgroup.procs write,
    the LEAF line precedes exec — and the script parses (bash -n only)."""
    with tempfile.TemporaryDirectory() as d:
        rec, script = runner.run_phase("stage", d, "kira --parallel=4 jobs.yaml",
                                       CFG, None, dry_run=True)
    assert rec["exit_status"] is None and rec["dry_run"]
    for tok in ("set -eu", "${CG:?}", 'case "$CG"', "grep -qx", "cgroup.procs",
                "/proc/self/cgroup", "LEAF $CG pid $$", "sudo -n",
                f"exit {runner.RC_CG_PATH_REFUSED}", f"exit {runner.RC_CG_NOT_ENTERED}",
                "cmd=kira --parallel=4 jobs.yaml"):
        assert tok in script, tok
    assert "prlimit" not in script and "pgrep" not in script \
        and "pstree" not in script and "pgid" not in script
    assert script.index('case "$CG"') < script.index("mkdir")
    assert script.index('tee "$CG/cgroup.procs"') < script.index("grep -qx")
    assert script.index("LEAF $CG") < script.index("\nexec ")
    assert subprocess.run(["bash", "-n", "-c", script]).returncode == 0


def test_shell_path_guard_refuses_root_and_dot_components():
    """The emitted guard, run on its own with $CG forced: root / bare root /
    empty / dot components / foreign prefix -> exit 97 (empty -> the ${CG:?}
    abort), a proper child -> passes. No sudo is reached."""
    def run(cg):
        p = subprocess.run(["bash", "-c", f"set -eu\nCG={shlex.quote(cg)}\n"
                            f"{runner.CG_PATH_GUARD}\necho PASSED"],
                           capture_output=True, text=True)
        return p.returncode, p.stdout
    for cg in ("/sys/fs/cgroup/", "/sys/fs/cgroup", "/sys/fs/cgroup/../x",
               "/sys/fs/cgroup/a/./b", "/sys/fs/cgroup/..", "/var/x", "/", "x"):
        rc, out = run(cg)
        assert rc == runner.RC_CG_PATH_REFUSED and "PASSED" not in out, cg
    rc, out = run("")
    assert rc != 0 and "PASSED" not in out
    rc, out = run("/sys/fs/cgroup/ok.leaf-1")
    assert rc == 0 and out.strip() == "PASSED"


def test_descends_from_synthetic_chain():
    tree = {50: 40, 40: 30, 30: 1, 1: 0, 77: 1, 88: 99}
    pp = tree.get
    assert runner.descends_from(50, 30, pp) and runner.descends_from(30, 30, pp)
    assert not runner.descends_from(77, 30, pp)      # reparented to PID 1
    assert not runner.descends_from(1, 30, pp)       # PID 1 itself
    assert not runner.descends_from(88, 30, pp)      # chain vanishes (None)
    cyc = {5: 6, 6: 5}
    assert not runner.descends_from(5, 9, cyc.get)   # bounded walk


def _leaf_fixture(td, name):
    leaf = os.path.join(td, name)
    os.makedirs(leaf)
    return os.path.join(leaf, "cgroup.procs")


def test_control_d_foreign_pid_in_leaf_aborts():
    """Post-launch ppid check: a fixture cgroup.procs listing the child AND a
    foreign pid (PID 1 — the specimen class) -> the abort path fires: SIGKILL
    (mocked) on the child only, every leaf pid moved (mocked) to the parent,
    ledger `abort` record, exit_status 99."""
    with tempfile.TemporaryDirectory() as td:
        name = "seedling-test-leaf"
        procs = _leaf_fixture(td, name)
        killed, moved, echoed = [], [], []

        def launch(script):
            p = subprocess.Popen(["sleep", "30"])
            with open(procs, "w") as fh:
                fh.write(f"{p.pid}\n1\n")
            return p

        def kill(proc):            # the mock records; it ends the plain `sleep`
            killed.append(proc.pid)
            proc.kill()

        def move(pid, dest):       # the mock records; no cgroup is touched
            moved.append((pid, dest))
            return True

        ledger = os.path.join(td, "ledger.jsonl")
        rec, _ = runner.run_phase(
            "stage", td, "sleep 30", dict(CFG, name=name), ledger,
            launch_ops={"cg_root": td, "launch": launch, "kill": kill,
                        "move": move, "in_leaf": lambda pid, n: True,
                        "poll_s": 0.05, "echo": echoed.append})
        assert rec["exit_status"] == runner.RC_CG_FOREIGN_PID != 0
        child = killed[0]
        assert killed == [child] and child != 1
        assert rec["abort"]["child_pid"] == child
        assert [r["pid"] for r in rec["abort"]["foreign"]] == [1]
        assert rec["abort"]["foreign"][0]["ppid"] == 0
        assert rec["abort"]["foreign"][0]["cmd"]          # cmdline printed, not None
        assert [m[0] for m in moved] == [child, 1]        # the leaf is emptied ...
        assert all(m[1] == td for m in moved)             # ... into the parent
        assert rec["leaf"]["guard"] == "foreign_pid"
        led = [json.loads(l) for l in open(ledger)]
        assert [r["phase"] for r in led] == ["abort", "stage"]
        assert led[0]["child_pid"] == child and led[0]["foreign"][0]["pid"] == 1
        assert led[0]["moved"][1]["pid"] == 1
        assert any(l.startswith("LEAF-CENSUS") and "foreign=1" in l for l in echoed)
        assert any(l.startswith("ABORT") and "SIGKILL child" in l for l in echoed)
        assert rec["wall_s"] < 10


def test_positive_control_own_descendants_do_not_abort():
    """A leaf holding only the child and its grandchild (ppid chain -> child)
    is lawful: no kill, no move, no abort record, the real exit status."""
    with tempfile.TemporaryDirectory() as td:
        name = "seedling-test-leaf"
        procs = _leaf_fixture(td, name)
        pidfile = os.path.join(td, "grandchild.pid")
        touched, echoed = [], []

        def launch(script):
            p = subprocess.Popen(["bash", "-c", f"sleep 1 & echo $! > "
                                  f"{shlex.quote(pidfile)}; wait"])
            for _ in range(500):
                if os.path.exists(pidfile) and open(pidfile).read().strip():
                    break
                time.sleep(0.01)
            gc = int(open(pidfile).read().strip())
            with open(procs, "w") as fh:
                fh.write(f"{p.pid}\n{gc}\n")
            return p

        ledger = os.path.join(td, "ledger.jsonl")
        rec, _ = runner.run_phase(
            "stage", td, "bash -c 'sleep 1 & wait'", dict(CFG, name=name), ledger,
            launch_ops={"cg_root": td, "launch": launch,
                        "kill": lambda p: touched.append("kill"),
                        "move": lambda pid, d: touched.append("move"),
                        "in_leaf": lambda pid, n: True,
                        "poll_s": 0.05, "echo": echoed.append})
        assert rec["exit_status"] == 0 and touched == [] and "abort" not in rec
        assert rec["leaf"]["n_foreign"] == 0 and rec["leaf"]["n_own"] == 2
        assert rec["leaf"]["launch_line"].startswith("LEAF-CENSUS") \
            and "foreign=0" in rec["leaf"]["launch_line"] \
            and "own=2" in rec["leaf"]["launch_line"]
        led = [json.loads(l) for l in open(ledger)]
        assert [r["phase"] for r in led] == ["stage"]


def test_guard_exit_codes_named_in_record():
    """A leaf script that exits 97/98 (path guard / not entered) is named in
    the record's leaf.guard; nothing is aborted (no leaf to empty)."""
    with tempfile.TemporaryDirectory() as td:
        name = "seedling-test-leaf"
        _leaf_fixture(td, name)
        for code, guard in ((runner.RC_CG_PATH_REFUSED, "path_refused"),
                            (runner.RC_CG_NOT_ENTERED, "not_entered"),
                            (0, None)):
            touched = []
            rec, _ = runner.run_phase(
                "stage", td, f"exit {code}", dict(CFG, name=name), None,
                launch_ops={"cg_root": td,
                            "launch": lambda s, c=code: subprocess.Popen(
                                ["bash", "-c", f"exit {c}"]),
                            "kill": lambda p: touched.append("kill"),
                            "move": lambda pid, d: touched.append("move"),
                            "poll_s": 0.05, "echo": lambda line: None})
            assert rec["exit_status"] == code and rec["leaf"]["guard"] == guard
            assert touched == [] and "abort" not in rec


N_LIVE = 3                       # children moved into the live leaf
SYSFS_CGROUP = "/sys/fs/cgroup"  # cgroup v2 mount (the same root the leaf script guards)


def _own_cgroup_dir():
    """/sys/fs/cgroup<path> of THIS process from the '0::' line of
    /proc/self/cgroup; None when there is no cgroup v2 line."""
    try:
        with open("/proc/self/cgroup") as fh:
            for line in fh:
                if line.startswith("0::"):
                    return SYSFS_CGROUP + (line[3:].strip() or "/")
    except OSError:
        return None
    return None


def _size_keyed_liveness(path):
    """The INERT check: `test -s cgroup.procs` / os.path.getsize(...) > 0.
    On a cgroup v2 interface file this is False whatever the file holds."""
    return os.path.getsize(path) > 0


def _read_pids(path):
    """The honest check: read the file, one pid per line."""
    with open(path) as fh:
        return [int(tok) for tok in fh.read().split()]


def _errno_words(ex):
    return f"errno {ex.errno} {errno.errorcode.get(ex.errno, '?')} ({ex.strerror})"


def _teardown_live_leaf(leaf, children, procs):
    """Kill the children BY PID (terminate, wait; kill if one lingers), then
    rmdir the leaf with a bounded retry; a leaf that will not go is reported
    with what its cgroup.procs still lists — never left silently."""
    for p in children:
        try:
            p.terminate()
        except ProcessLookupError:
            pass
    for p in children:
        try:
            p.wait(timeout=10)
        except subprocess.TimeoutExpired:
            p.kill()
            p.wait(timeout=10)
    if not os.path.isdir(leaf):
        return
    last = None
    for _ in range(50):
        try:
            os.rmdir(leaf)
            return
        except OSError as ex:
            last = ex
            time.sleep(0.1)
    try:
        lingering = _read_pids(procs)
    except OSError:
        lingering = "<unreadable>"
    raise AssertionError(f"live leaf {leaf} not removed: {_errno_words(last)}; "
                         f"cgroup.procs still lists {lingering}")


def test_positive_control_live_leaf_size_zero_with_N_pids():
    """POSITIVE CONTROL, LIVE: cgroup v2 interface files stat at size 0 while
    carrying content. A real leaf made under this process's own cgroup and
    holding N = 3 children we spawned ourselves reports
    os.path.getsize(cgroup.procs) == 0 AND len(read().split()) == N with the
    pid set equal to the children's — so the `-s cgroup.procs` liveness form
    (_size_keyed_liveness) answers False on a populated leaf while the
    read-based form answers N. Both are asserted: if a kernel ever reports a
    non-zero st_size for cgroup.procs this test FAILS by name, and the
    size-keyed check would then be the one to re-examine. The cgroup safety
    law is kept to the letter: the leaf is a child of the CALLING process's
    own cgroup (never a parent's, never a foreign path), only pids of children
    this test spawned are written (never its own pid, never a foreign one),
    and teardown runs in `finally` (terminate + wait by pid, then rmdir). A
    refused mkdir or write SKIPS by name with the errno."""
    own = _own_cgroup_dir()
    if own is None:
        pytest.skip("no cgroup v2 '0::' line in /proc/self/cgroup: no live leaf")
    leaf = os.path.join(own, f"seedling-cgsize-{os.getpid()}")
    procs = os.path.join(leaf, "cgroup.procs")
    children = []
    try:
        try:
            os.mkdir(leaf)
        except OSError as ex:
            pytest.skip(f"live leaf refused: mkdir {leaf}: {_errno_words(ex)}")
        for _ in range(N_LIVE):
            children.append(subprocess.Popen(["sleep", "30"], stdout=subprocess.DEVNULL,
                                             stderr=subprocess.DEVNULL))
        for p in children:                       # one write per pid, own children only
            try:
                with open(procs, "w") as fh:
                    fh.write(f"{p.pid}\n")
            except OSError as ex:
                pytest.skip(f"live leaf refused: write pid {p.pid} to {procs}: "
                            f"{_errno_words(ex)}")
        size = os.path.getsize(procs)
        pids = _read_pids(procs)
        print(f"LEAF-SIZE-CONTROL leaf={leaf} n={N_LIVE} "
              f"pids={','.join(str(p.pid) for p in children)} getsize={size} "
              f"read_n={len(pids)}")
        # THE FACT: size 0, content N
        assert size == 0, f"cgroup.procs st_size {size} != 0 (kernel semantics changed?)"
        assert len(pids) == N_LIVE, pids
        assert set(pids) == {p.pid for p in children}, (pids, [p.pid for p in children])
        # THE INERT CHECK vs THE HONEST ONE, side by side on the same live leaf
        assert _size_keyed_liveness(procs) is False   # `test -s` says "empty"
        assert len(_read_pids(procs)) == N_LIVE       # the read says N
    finally:
        _teardown_live_leaf(leaf, children, procs)


if __name__ == "__main__":
    # direct run: every test in file order; a pytest.skip is reported by name
    # and the exit status is 2 (ALL PASS only at 0 skipped)
    tests = [obj for name, obj in list(globals().items())
             if name.startswith("test_") and callable(obj)]
    skipped = []
    for fn in tests:
        try:
            fn()
        except pytest.skip.Exception as ex:
            skipped.append(fn.__name__)
            print(f"SKIP {fn.__name__}: {ex}")
    if skipped:
        print(f"test_cgroup_law: PASS_EXECUTED_WITH_SKIPS ({len(tests) - len(skipped)} "
              f"executed / {len(skipped)} skipped: {', '.join(skipped)})")
        sys.exit(2)
    print(f"test_cgroup_law: ALL PASS ({len(tests)})")

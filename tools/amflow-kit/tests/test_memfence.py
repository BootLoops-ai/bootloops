#!/usr/bin/env python3
"""memfence battery (pytest; every leg is a test_* function so it collects).

Planted controls, each one measured on a real process where a process is
involved; nothing is asserted from a docstring.

  leaf fixture     memory.max = 68719476736 on the launcher's own cgroup ->
                   mode 'leaf', no ulimit, export '0' (the binary's no-cap
                   value)
  max fixture      memory.max = 'max' -> mode 'fuse'; env cap 48 -> fuse
                   120 GiB (2.5x) = 125829120 kB named 'address-space fuse';
                   env cap 15 -> 37.5 GiB = 39321600 kB
  binary semantics unset / '' / '0' / '-3' / 'abc' parse as no cap, '48' and
                   the strtod prefix '48abc' as 48 GiB
  max rule         env cap 100, budget 20 GiB -> fuse = max(100, 50) = 100
  readback match   a real `sleep 30` child spawned with the intended env
                   reads back equal from /proc/<pid>/environ
  comment swallow  `# export AMFLOW_X=1` inside a bash -c script -> readback
                   mismatch -> abort_by_name kills THAT pid (resolved at
                   spawn, identity-checked by cmdline) and it is gone
  identity refusal abort_by_name with the wrong argv REFUSES and the child
                   survives (then the test reaps it)
  live read        leaf_memory_max() on the running machine returns a readable
                   record (no claim about which mode that machine is in)
"""
import os
import subprocess
import sys
import time

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))  # tools/amflow-kit, the kit dir
from amflow_kit import memfence  # noqa: E402

GIB = 1024 ** 3


def _fixture_reads(monkeypatch, memory_max_text, cg="/fixture/leaf.scope"):
    """Serve /proc/self/cgroup and the leaf's memory.max from memory."""
    table = {
        "/proc/self/cgroup": f"0::{cg}\n",
        f"/sys/fs/cgroup{cg}/memory.max": memory_max_text + "\n",
    }

    def fake_read(path):
        return table.get(path)  # ancestors: unreadable -> None

    monkeypatch.setattr(memfence, "_read_text", fake_read)
    return cg


def test_leaf_fixture_no_ulimit_and_no_fuse_export(monkeypatch):
    cg = _fixture_reads(monkeypatch, "68719476736")
    leaf = memfence.leaf_memory_max()
    assert leaf["cgroup_path"] == cg
    assert leaf["memory_max_bytes"] == 68719476736
    assert leaf["decision"] == "leaf"
    pol = memfence.fuse_policy("48", leaf=leaf)
    assert pol["mode"] == "leaf"
    assert pol["ulimit_v_kb"] is None
    assert pol["export_cap_gb"] == "0" == memfence.NO_FUSE_EXPORT
    assert memfence.ulimit_prefix(pol) == ""
    assert pol["exports"][memfence.CAP_ENV] == "0"
    assert "the leaf is the fence" in pol["words"]
    assert "68719476736" in pol["words"]
    # the exported value is one the binary treats as no cap
    assert memfence.parse_cap_gb(pol["export_cap_gb"]) is None


def test_max_fixture_fuse_is_2p5x_budget(monkeypatch):
    _fixture_reads(monkeypatch, "max")
    leaf = memfence.leaf_memory_max()
    assert leaf["memory_max_bytes"] is None and leaf["decision"] == "max"
    pol = memfence.fuse_policy("48", leaf=leaf)
    assert pol["mode"] == "fuse"
    assert pol["fuse_bytes"] == 120 * GIB
    assert pol["ulimit_v_kb"] == 120 * GIB // 1024 == 125829120
    assert pol["export_cap_gb"] == "120"
    assert memfence.ulimit_prefix(pol) == "ulimit -v 125829120; "
    assert ("address-space fuse (RLIMIT_AS) = 120 GiB = 2.5 x the memory "
            "budget 48 GiB") in pol["words"]
    assert "not the memory cap" in pol["words"]
    pol15 = memfence.fuse_policy("15", leaf=leaf)
    assert pol15["fuse_bytes"] == int(37.5 * GIB)
    assert pol15["ulimit_v_kb"] == 39321600
    assert pol15["export_cap_gb"] == "37.5"
    assert memfence.parse_cap_gb(pol15["export_cap_gb"]) == 37.5


def test_max_rule_env_cap_above_scaled_budget(monkeypatch):
    _fixture_reads(monkeypatch, "max")
    pol = memfence.fuse_policy("100", budget_bytes=20 * GIB,
                               leaf=memfence.leaf_memory_max())
    assert pol["mode"] == "fuse"
    assert pol["fuse_bytes"] == 100 * GIB
    assert pol["export_cap_gb"] == "100"
    assert "max(env cap 100 GiB, 2.5 x budget 20 GiB = 50 GiB)" in pol["arithmetic"]


def test_unfenced_when_nothing_to_fence(monkeypatch):
    _fixture_reads(monkeypatch, "max")
    pol = memfence.fuse_policy(None, leaf=memfence.leaf_memory_max())
    assert pol["mode"] == "unfenced"
    assert pol["ulimit_v_kb"] is None and pol["export_cap_gb"] is None
    env = {memfence.CAP_ENV: "48", "KEEP": "1"}
    memfence.apply_policy_to_env(env, pol)
    assert memfence.CAP_ENV not in env and env["KEEP"] == "1"
    assert env[memfence.MODE_ENV] == "unfenced"


def test_binary_no_cap_semantics():
    # kira_run.cpp: cap && *cap, strtod, gib > 0.0 — else no RLIMIT_AS
    for v in (None, "", "0", "0.0", "-3", "abc", " ", "nan"):
        assert memfence.parse_cap_gb(v) is None, v
    assert memfence.parse_cap_gb("48") == 48.0
    assert memfence.parse_cap_gb("37.5") == 37.5
    assert memfence.parse_cap_gb("48abc") == 48.0   # strtod prefix
    assert memfence.parse_cap_gb(" 15 ") == 15.0
    assert memfence.format_gib(120 * GIB) == "120"
    assert memfence.format_gib(int(37.5 * GIB)) == "37.5"


def _spawn(script, env, expect=("sleep", "30")):
    """bash -c child that exec's `sleep 30`; returns the Popen once the pid's
    cmdline IS the exec'd sleep (the same pid, after bash's exec), so the
    identity check in abort_by_name sees ['sleep','30']."""
    p = subprocess.Popen(["bash", "-c", script], env=env,
                         stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    deadline = time.time() + 5.0
    while time.time() < deadline and memfence._cmdline(p.pid) != list(expect):
        time.sleep(0.02)
    return p


def _reap(p):
    try:
        p.kill()
    except ProcessLookupError:
        pass
    try:
        p.wait(timeout=5)
    except Exception:
        pass


def test_readback_match_on_real_child():
    env = {"PATH": os.environ.get("PATH", "/usr/bin:/bin"),
           memfence.CAP_ENV: "120", memfence.MODE_ENV: "fuse",
           memfence.FUSE_GIB_ENV: "120", memfence.BUDGET_GIB_ENV: "48"}
    p = _spawn("exec sleep 30", env)
    try:
        intended = memfence.intended_from_env(env)
        rb = memfence.env_readback(p.pid, intended)
        assert rb["verified"] is True
        assert rb["match"] is True, rb
        assert rb["diff"] == []
        assert rb["observed"][memfence.CAP_ENV] == "120"
        assert rb["cmdline"] == ["sleep", "30"]
        assert "match" in memfence.readback_line(rb)
    finally:
        _reap(p)


def test_comment_swallowed_export_mismatch_aborts_by_name():
    env = {"PATH": os.environ.get("PATH", "/usr/bin:/bin")}
    # the export line is swallowed by a comment: the child never sees AMFLOW_X
    p = _spawn("# export AMFLOW_X=1\nexec sleep 30", env)
    try:
        rb = memfence.env_readback(p.pid, {"AMFLOW_X": "1"})
        assert rb["verified"] is True
        assert rb["match"] is False
        assert rb["diff"] == [{"key": "AMFLOW_X", "intended": "1", "observed": None}]
        rec = memfence.abort_by_name(p.pid, ["sleep", "30"],
                                     reason="env readback mismatch")
        assert rec["refused"] is None
        assert rec["killed"] is True and rec["gone"] is True, rec
        p.wait(timeout=5)
        with pytest.raises(ProcessLookupError):
            os.kill(p.pid, 0)
    finally:
        _reap(p)


def test_abort_by_name_refuses_wrong_identity():
    env = {"PATH": os.environ.get("PATH", "/usr/bin:/bin")}
    p = _spawn("exec sleep 30", env)
    try:
        rec = memfence.abort_by_name(p.pid, ["some-other-program", "--flag"])
        assert rec["refused"] and "identity mismatch" in rec["refused"]
        assert rec["killed"] is False
        assert p.poll() is None          # still alive: the refusal held
    finally:
        _reap(p)


def test_readback_unresolved_pid_is_not_a_mismatch():
    rb = memfence.env_readback(None, {memfence.CAP_ENV: "120"})
    assert rb["verified"] is False and rb["match"] is False and rb["diff"] == []
    assert rb["error"] == "no pid resolved"


def test_find_spawned_pid_by_unique_tail():
    env = {"PATH": os.environ.get("PATH", "/usr/bin:/bin")}
    tag = f"--memfence-test-{os.getpid()}-{int(time.time())}"
    code = "import time; time.sleep(30)"
    t0 = time.time()
    # a wrapper shell exec's the real program; the unique tag is the identity
    p = subprocess.Popen(["bash", "-c", f"exec {sys.executable} -c '{code}' {tag}"],
                         env=env, stdout=subprocess.DEVNULL,
                         stderr=subprocess.DEVNULL)
    try:
        pid = memfence.find_spawned_pid(["-c", code, tag], t0, window=5.0)
        assert pid == p.pid
    finally:
        _reap(p)


def test_find_spawned_pid_exclude_skips_the_wrapper_stage(tmp_path):
    """A wrapper script whose cmdline ends with the same tail as the program
    it exec's: with exclude=(wrapper,) only the exec'd program is returned,
    never the wrapper stage (whose environ is still the launcher's)."""
    wrap = tmp_path / "wrap.sh"
    wrap.write_text("#!/bin/bash\nsleep 0.3\nexec \"$@\"\n")   # a slow wrapper stage
    wrap.chmod(0o755)
    env = {"PATH": os.environ.get("PATH", "/usr/bin:/bin")}
    tag = f"--memfence-excl-{os.getpid()}-{int(time.time())}"
    code = "import time; time.sleep(30)"
    t0 = time.time()
    p = subprocess.Popen([str(wrap), sys.executable, "-c", code, tag], env=env,
                         stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    try:
        # during the wrapper stage the tail already matches the wrapper's cmdline
        early = memfence.find_spawned_pid(["-c", code, tag], t0, window=5.0, retries=1)
        assert early == p.pid and memfence._cmdline(p.pid)[0].endswith("bash")
        pid = memfence.find_spawned_pid(["-c", code, tag], t0, window=5.0,
                                        exclude=(str(wrap),))
        assert pid == p.pid                      # same pid, after its exec
        assert memfence._cmdline(pid)[0] == sys.executable
        assert str(wrap) not in memfence._cmdline(pid)
    finally:
        _reap(p)


def test_tree_readback_catches_a_wrapper_that_drops_the_cap(tmp_path):
    """The direct child (a wrapper script) carries the cap; it unsets the
    variable before exec'ing the real program. The direct readback matches;
    only the tree readback sees the exec'd descendant without it — and the
    tree abort kills wrapper + descendant by pid."""
    wrap = tmp_path / "wrap.sh"
    wrap.write_text("#!/bin/bash\n# export AMFLOW_KIRA_MEM_CAP_GB=120  (swallowed)\n"
                    "unset AMFLOW_KIRA_MEM_CAP_GB\nsleep 30 & wait $!\n")
    wrap.chmod(0o755)
    env = {"PATH": os.environ.get("PATH", "/usr/bin:/bin"),
           memfence.CAP_ENV: "120", memfence.MODE_ENV: "fuse",
           memfence.FUSE_GIB_ENV: "120", memfence.BUDGET_GIB_ENV: "48"}
    p = subprocess.Popen([str(wrap)], env=env, stdout=subprocess.DEVNULL,
                         stderr=subprocess.DEVNULL)
    try:
        intended = memfence.intended_from_env(env)
        direct = memfence.env_readback(p.pid, intended)
        assert direct["match"] is True            # the wrapper itself has the cap
        tree = memfence.env_readback_tree(p.pid, intended, settle=0.5)
        assert tree["match"] is False, tree
        assert any(d["key"] == memfence.CAP_ENV and d["observed"] is None
                   and d.get("pid") for d in tree["diff"])
        kids = memfence.descendants(p.pid)
        assert kids, "the wrapper's sleep child was not found by lineage"
        recs = memfence.abort_tree_by_pid(p.pid, [str(wrap)], reason="control")
        p.wait(timeout=5)
        time.sleep(0.2)
        assert all(r["gone"] for r in recs), recs
        for k in kids:
            assert not memfence._alive(k)
    finally:
        _reap(p)


def test_nonce_identity_separates_same_argv_launches_in_the_window():
    """Two `sleep 30` children from the same cwd inside one start-time window:
    argv + time cannot tell them apart; the launch nonce in each environ can.
    The abort with the wrong nonce REFUSES; with the right one it kills only
    its own child."""
    base = {"PATH": os.environ.get("PATH", "/usr/bin:/bin")}
    n1, n2 = memfence.new_launch_id(), memfence.new_launch_id()
    t0 = time.time()
    p1 = _spawn("exec sleep 30", dict(base, **{memfence.LAUNCH_ID_ENV: n1}))
    p2 = _spawn("exec sleep 30", dict(base, **{memfence.LAUNCH_ID_ENV: n2}))
    try:
        found2 = memfence.find_spawned_pid(["sleep", "30"], t0, window=5.0,
                                           require_env={memfence.LAUNCH_ID_ENV: n2})
        found1 = memfence.find_spawned_pid(["sleep", "30"], t0, window=5.0,
                                           require_env={memfence.LAUNCH_ID_ENV: n1})
        assert found1 == p1.pid and found2 == p2.pid and found1 != found2
        # wrong nonce at kill time -> refused, p1 survives
        rec = memfence.abort_by_name(p1.pid, ["sleep", "30"],
                                     require_env={memfence.LAUNCH_ID_ENV: n2})
        assert rec["refused"] and "environ does not carry" in rec["refused"]
        assert p1.poll() is None
        # right nonce -> p2 killed, p1 untouched
        rec = memfence.abort_by_name(p2.pid, ["sleep", "30"],
                                     require_env={memfence.LAUNCH_ID_ENV: n2})
        assert rec["killed"] and rec["gone"], rec
        p2.wait(timeout=5)
        assert p1.poll() is None
    finally:
        _reap(p1)
        _reap(p2)


def test_live_cgroup_read_is_a_record():
    leaf = memfence.leaf_memory_max()
    assert set(leaf) >= {"cgroup_path", "memory_max_bytes", "ancestors", "decision"}
    if leaf["cgroup_path"] is not None:
        assert leaf["decision"] in ("leaf", "max", "unreadable")
    pol = memfence.fuse_policy("48", leaf=leaf)
    assert pol["mode"] in ("leaf", "fuse")
    assert pol["words"]

#!/usr/bin/env python3
"""memfence positive control: cgroup v2 interface files stat at size 0 while
carrying content (pytest; every leg is a test_* function so it collects).

memfence reads `memory.max` and `/proc/<pid>/cgroup` through `_read_text`
(open + read) and never keys on a file's size. This module pins WHY that is
the only honest reader: a size-keyed liveness check of the form
`test -s cgroup.procs` / `os.path.getsize(...) > 0` answers "empty" on every
cgroup v2 interface file whatever it holds. Two LIVE legs, no monkeypatch:

  own scope        the calling process's own `memory.max` stats at size 0
                   while its text parses (`max` or an integer) through
                   memfence.leaf_memory_max — no leaf needed
  live leaf        a real leaf made under the calling process's OWN cgroup,
                   holding N = 3 `sleep 30` children this test spawned:
                   getsize(cgroup.procs) == 0 AND the read lists exactly the
                   N pids; memfence.own_cgroup_path(pid) of each child names
                   the leaf; the size-keyed check answers False on it

Both assert the size AND the content, so a kernel that ever reports a
non-zero st_size for these files FAILS the legs by name. Where the mkdir or
the cgroup.procs write is refused the leg SKIPS by name with the errno
(pytest -rs shows it; the exit status stays 0 under plain pytest). Safety:
the leaf is a child of the calling process's own cgroup, only pids of
children spawned here are written (one write per pid), teardown runs in
`finally` (terminate + wait by pid, rmdir with a bounded retry).
"""
import errno
import os
import subprocess
import sys
import time

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))  # tools/amflow-kit, the kit dir
from amflow_kit import memfence  # noqa: E402

N_LIVE = 3


def _size_keyed_says_present(path):
    """The inert liveness check as a function: `test -s FILE` / getsize > 0."""
    return os.path.getsize(path) > 0


def _errno_words(ex):
    return f"errno {ex.errno} {errno.errorcode.get(ex.errno, '?')} ({ex.strerror})"


def test_own_scope_memory_max_is_size_zero_with_parsable_text():
    """memfence.leaf_memory_max reads the launcher's own memory.max by
    content; that file stats at 0 bytes while the text parses. A launcher
    keyed on `-s memory.max` would call every fence unset."""
    leaf = memfence.leaf_memory_max()
    if leaf["cgroup_path"] is None:
        pytest.skip("no cgroup v2 '0::' line in /proc/self/cgroup")
    path = leaf["memory_max_path"]
    try:
        size = os.path.getsize(path)
    except OSError as ex:
        pytest.skip(f"own scope has no readable memory.max ({path}): {_errno_words(ex)}")
    raw = leaf["memory_max_raw"]
    if raw is None:
        pytest.skip(f"own scope memory.max unreadable by memfence._read_text: {path}")
    print(f"MEMORY-MAX-SIZE-CONTROL path={path} getsize={size} text={raw}")
    assert size == 0, f"memory.max st_size {size} != 0 (kernel semantics changed?)"
    assert _size_keyed_says_present(path) is False
    assert raw == "max" or int(raw) > 0            # the content memfence keys on
    assert leaf["decision"] in ("leaf", "max")     # never 'unreadable' on a live scope
    assert memfence._read_text(path).strip() == raw


def _teardown(leaf, children, procs):
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
    left = memfence._read_text(procs)
    raise AssertionError(f"live leaf {leaf} not removed: {_errno_words(last)}; "
                         f"cgroup.procs still holds {left!r}")


def test_live_leaf_cgroup_procs_size_zero_with_N_pids():
    """A live leaf under this process's own cgroup with N = 3 own children:
    st_size 0, N pids on read, the set equal to the children's, and
    memfence.own_cgroup_path(child) == the leaf. The size-keyed check answers
    False on the populated leaf; the read answers N. Refusals skip by name."""
    cg = memfence.own_cgroup_path()
    if cg is None:
        pytest.skip("no cgroup v2 '0::' line in /proc/self/cgroup: no live leaf")
    name = f"memfence-cgsize-{os.getpid()}"
    leaf = os.path.join(memfence.SYSFS_CGROUP, cg.lstrip("/"), name)
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
        for p in children:
            try:
                with open(procs, "w") as fh:
                    fh.write(f"{p.pid}\n")
            except OSError as ex:
                pytest.skip(f"live leaf refused: write pid {p.pid} to {procs}: "
                            f"{_errno_words(ex)}")
        want = cg.rstrip("/") + "/" + name
        for p in children:                         # memfence's own reader sees the move
            assert memfence.own_cgroup_path(p.pid) == want, (p.pid, want)
        size = os.path.getsize(procs)
        text = memfence._read_text(procs)
        pids = [int(tok) for tok in text.split()]
        print(f"LEAF-SIZE-CONTROL leaf={leaf} n={N_LIVE} "
              f"pids={','.join(str(p.pid) for p in children)} getsize={size} "
              f"read_n={len(pids)}")
        assert size == 0, f"cgroup.procs st_size {size} != 0 (kernel semantics changed?)"
        assert len(pids) == N_LIVE, pids
        assert set(pids) == {p.pid for p in children}, (pids, [p.pid for p in children])
        assert _size_keyed_says_present(procs) is False   # `test -s` says "empty"
        assert len(memfence._read_text(procs).split()) == N_LIVE
    finally:
        _teardown(leaf, children, procs)

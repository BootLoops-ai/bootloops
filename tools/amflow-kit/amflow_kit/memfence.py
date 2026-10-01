#!/usr/bin/env python3
"""amflow_kit.memfence — the shared memory-fence policy for the AMFlow/Kira launchers.
Member of amflow-kit, the house kit around the AMFlow.cpp fork (../MEMFENCE.md).

One helper, several call sites (pmflow/pmflow.py, numkin/shift_opt.py, and
any launcher that follows the same pattern). It decides, for one
spawn, how the child's memory is fenced and it verifies after the spawn that
the child really received the environment the launcher intended.

The defect this replaces. AMFLOW_KIRA_MEM_CAP_GB used to be turned into a
fixed RLIMIT_AS (``ulimit -v``) sized as the memory budget. Under jemalloc the
virtual set (VmSize) of a Kira reduction runs 2-5x its resident set, so an
address-space limit "sized as the budget" kills the job at roughly a third of
the memory it was meant to be allowed: std::bad_alloc at 23.8 GB RSS under a
48 GiB cap inside a 64 GiB cgroup leaf (memory.events oom 0). The amflow-cpp binary applies the same variable a second
time, as RLIMIT_AS on the Kira child it forks (kira_run.cpp, the
AMFLOW_KIRA_MEM_CAP_GB block), so exporting the 1x budget reinstates the fuse
one level down even when the launcher sets no ulimit.

Policy (fuse_policy):

  * inside a cgroup v2 leaf with a finite memory.max  -> mode 'leaf'
      the leaf is the fence: no ulimit -v, and AMFLOW_KIRA_MEM_CAP_GB is
      exported as '0', the value the binary treats as "no cap" (it applies
      RLIMIT_AS only when strtod(value) > 0; unset, empty, '0', a negative
      number and a non-numeric string all mean no cap).
  * no finite memory.max on the leaf                   -> mode 'fuse'
      RLIMIT_AS = max(env cap, FUSE_FACTOR x budget), where the budget is the
      caller's intended memory budget (today the env cap IS the budget, so the
      fuse is 2.5x the env cap). The launch line and the receipt name it an
      ADDRESS-SPACE fuse, never "the memory cap". The same fuse value is
      exported as AMFLOW_KIRA_MEM_CAP_GB so the binary applies the same fuse,
      not the 1x budget, to its Kira child.
  * no leaf fence and no budget at all                 -> mode 'unfenced'
      nothing exported, no ulimit (the pre-existing behaviour of a launcher
      that had no cap to apply).

The decision reads the LEAF's own memory.max (the cgroup named by
/proc/self/cgroup), not an ancestor's. Ancestors are recorded for the receipt
but do not decide: on a shared box the user slice above every login carries
a finite memory.max equal to the whole machine, and treating that as "inside a
leaf" would remove the fuse from every launch on the box, which is exactly the
runaway-symbolic-swell case the fuse exists for.

Readback (env_readback): after the spawn the launcher reads
/proc/<pid>/environ of the child it resolved and compares the cap fields with
what it intended. A mismatch means the child is running under a different
fence than the receipt claims (the classic cause: an ``export`` line swallowed
by a comment in a launch script), and the launcher ABORTS BY NAME — it kills
the pid it resolved at spawn, after checking that pid's cmdline is the argv it
launched, never by pattern (abort_by_name).

Stdlib only. Linux /proc and cgroup v2 assumed; every read is best-effort and
an unreadable file is reported, never guessed.
"""
import json
import os
import re
import signal
import time
import uuid

__all__ = [
    "CAP_ENV", "CAP_FIELDS", "FUSE_FACTOR", "NO_FUSE_EXPORT", "GIB",
    "parse_cap_gb", "format_gib", "own_cgroup_path", "leaf_memory_max",
    "fuse_policy", "apply_policy_to_env", "intended_from_env",
    "read_environ", "env_readback", "env_readback_tree", "descendants",
    "readback_line", "abort_by_name", "abort_tree_by_pid",
    "find_spawned_pid", "ulimit_prefix", "new_launch_id", "pid_env_matches",
    "LAUNCH_ID_ENV",
]

CAP_ENV = "AMFLOW_KIRA_MEM_CAP_GB"
# The informational twins the launcher exports beside the cap so the child's
# environment states the fence it runs under (inert to the binary).
MODE_ENV = "MEMFENCE_MODE"
FUSE_GIB_ENV = "MEMFENCE_FUSE_GIB"
BUDGET_GIB_ENV = "MEMFENCE_BUDGET_GIB"
CAP_FIELDS = (CAP_ENV, MODE_ENV, FUSE_GIB_ENV, BUDGET_GIB_ENV)
# Per-launch nonce (new_launch_id) a launcher may put in the child's env so
# the pid it resolves and aborts is provably its own (see new_launch_id).
LAUNCH_ID_ENV = "MEMFENCE_LAUNCH_ID"

FUSE_FACTOR = 2.5
# kira_run.cpp applies RLIMIT_AS only when strtod(value) > 0, so '0' is the
# explicit "no cap" export (unset would also do, but an explicit value is
# visible in the child's environ and in the readback).
NO_FUSE_EXPORT = "0"
GIB = 1024 ** 3

PROC = "/proc"
SYSFS_CGROUP = "/sys/fs/cgroup"

_STRTOD_PREFIX = re.compile(
    r"^\s*[+-]?(?:(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?|inf(?:inity)?|nan)",
    re.IGNORECASE)


def _read_text(path):
    """Read a small file; None when unreadable. Tests monkeypatch this."""
    try:
        with open(path, "r") as fh:
            return fh.read()
    except OSError:
        return None


def parse_cap_gb(value):
    """AMFLOW_KIRA_MEM_CAP_GB as the binary parses it: float GiB when
    strtod(value) > 0, else None (unset / '' / '0' / negative / non-numeric
    all mean "no cap" in kira_run.cpp)."""
    if value is None:
        return None
    s = str(value)
    if not s:
        return None
    m = _STRTOD_PREFIX.match(s)
    if not m:
        return None
    try:
        gib = float(m.group(0))
    except ValueError:
        return None
    if gib != gib or gib <= 0.0 or gib == float("inf"):
        return None
    return gib


def format_gib(nbytes):
    """Bytes -> GiB string the binary re-parses: '120', '37.5', '0.25'."""
    g = nbytes / GIB
    if abs(g - round(g)) < 1e-9:
        return str(int(round(g)))
    return ("%.6f" % g).rstrip("0").rstrip(".")


def own_cgroup_path(pid="self"):
    """cgroup v2 path of <pid> from /proc/<pid>/cgroup ('0::<path>'), or None."""
    txt = _read_text(f"{PROC}/{pid}/cgroup")
    if txt is None:
        return None
    for ln in txt.splitlines():
        if ln.startswith("0::"):
            return ln[3:].strip() or "/"
    return None


def _memory_max_of(cg_path, sysfs_root):
    f = os.path.join(sysfs_root, cg_path.lstrip("/"), "memory.max")
    raw = _read_text(f)
    if raw is None:
        return f, None, None
    raw = raw.strip()
    if raw == "max":
        return f, raw, None
    try:
        return f, raw, int(raw)
    except ValueError:
        return f, raw, None


def leaf_memory_max(pid="self", sysfs_root=SYSFS_CGROUP):
    """The fence the launching process sits in.

    Returns a dict:
      cgroup_path        '0::' path from /proc/<pid>/cgroup (None if unreadable)
      memory_max_path    the leaf's own memory.max file
      memory_max_raw     its text ('max', a number, or None if unreadable)
      memory_max_bytes   int when finite, None for 'max' / unreadable
      ancestors          [{path, raw, bytes}] for every parent up to the root,
                         informational only (see the module docstring for why
                         the ancestors do not decide)
      ancestor_min_bytes the smallest finite ancestor value, or None
      decision           'leaf' | 'max' | 'unreadable'
    """
    cg = own_cgroup_path(pid)
    out = {"cgroup_path": cg, "memory_max_path": None, "memory_max_raw": None,
           "memory_max_bytes": None, "ancestors": [], "ancestor_min_bytes": None,
           "decision": "unreadable"}
    if cg is None:
        return out
    f, raw, nbytes = _memory_max_of(cg, sysfs_root)
    out.update({"memory_max_path": f, "memory_max_raw": raw,
                "memory_max_bytes": nbytes})
    out["decision"] = ("leaf" if nbytes is not None
                       else ("max" if raw == "max" else "unreadable"))
    parent = cg
    while parent not in ("", "/"):
        parent = os.path.dirname(parent.rstrip("/")) or "/"
        if parent == "/":
            break
        pf, praw, pbytes = _memory_max_of(parent, sysfs_root)
        out["ancestors"].append({"path": parent, "raw": praw, "bytes": pbytes})
    finite = [a["bytes"] for a in out["ancestors"] if a["bytes"] is not None]
    out["ancestor_min_bytes"] = min(finite) if finite else None
    return out


def fuse_policy(env_cap_gb, budget_bytes=None, leaf=None, factor=FUSE_FACTOR):
    """Decide the fence for one spawn.

    env_cap_gb   the AMFLOW_KIRA_MEM_CAP_GB value the caller would export
                 (string or None); parsed with the binary's own semantics.
    budget_bytes the caller's intended memory budget; None means "the env cap
                 is the budget" (today's convention at every call site).
    leaf         a leaf_memory_max() dict; None reads the live cgroup.

    Returns {'mode': 'leaf'|'fuse'|'unfenced', 'ulimit_v_kb': None|int,
             'export_cap_gb': str|None, 'exports': {...}, 'words': '...',
             'arithmetic': '...', plus the numbers behind them}.
    """
    if leaf is None:
        leaf = leaf_memory_max()
    cap_gib = parse_cap_gb(env_cap_gb)
    cap_bytes = int(cap_gib * GIB) if cap_gib is not None else None
    if budget_bytes is None:
        budget_bytes = cap_bytes
    leaf_bytes = leaf.get("memory_max_bytes")
    cg = leaf.get("cgroup_path")
    pol = {"mode": None, "ulimit_v_kb": None, "export_cap_gb": None,
           "exports": {}, "words": "", "arithmetic": "",
           "env_cap_gb": None if env_cap_gb is None else str(env_cap_gb),
           "env_cap_bytes": cap_bytes, "budget_bytes": budget_bytes,
           "budget_gib": None if budget_bytes is None else budget_bytes / GIB,
           "fuse_bytes": None, "fuse_gib": None, "factor": factor,
           "leaf_memory_max_bytes": leaf_bytes, "cgroup_path": cg,
           "leaf_memory_max_raw": leaf.get("memory_max_raw"),
           "ancestor_min_bytes": leaf.get("ancestor_min_bytes")}
    if leaf_bytes is not None:
        pol["mode"] = "leaf"
        pol["export_cap_gb"] = NO_FUSE_EXPORT
        pol["exports"] = {CAP_ENV: NO_FUSE_EXPORT, MODE_ENV: "leaf",
                          FUSE_GIB_ENV: "0",
                          BUDGET_GIB_ENV: (format_gib(budget_bytes)
                                           if budget_bytes is not None else "0")}
        pol["arithmetic"] = (f"memory.max = {leaf_bytes} bytes = "
                             f"{format_gib(leaf_bytes)} GiB on {cg}; "
                             f"fuse = none")
        pol["words"] = (
            f"inside cgroup leaf {cg} with memory.max = {leaf_bytes} bytes "
            f"({format_gib(leaf_bytes)} GiB): the leaf is the fence; no "
            f"RLIMIT_AS (ulimit -v) is set and {CAP_ENV} is exported as "
            f"'{NO_FUSE_EXPORT}', which the amflow-cpp binary treats as no cap "
            f"(kira_run.cpp applies RLIMIT_AS only when strtod(value) > 0)"
            + (f"; the caller's memory budget {format_gib(budget_bytes)} GiB "
               f"is recorded, not enforced by a fuse"
               if budget_bytes is not None else ""))
        return pol
    if budget_bytes is None and cap_bytes is None:
        pol["mode"] = "unfenced"
        pol["exports"] = {}
        pol["arithmetic"] = "no env cap, no budget: fuse = none"
        pol["words"] = (
            f"no finite memory.max on cgroup {cg} "
            f"(memory.max = {leaf.get('memory_max_raw')}) and no memory "
            f"budget given ({CAP_ENV} unset): nothing exported, no RLIMIT_AS")
        return pol
    budget_bytes = int(budget_bytes if budget_bytes is not None else cap_bytes)
    scaled = int(round(factor * budget_bytes))
    fuse_bytes = max(cap_bytes or 0, scaled)
    fuse_kb = -(-fuse_bytes // 1024)  # ceil
    fuse_gib_s = format_gib(fuse_bytes)
    budget_gib_s = format_gib(budget_bytes)
    pol.update({"mode": "fuse", "ulimit_v_kb": fuse_kb,
                "export_cap_gb": fuse_gib_s, "fuse_bytes": fuse_bytes,
                "fuse_gib": fuse_bytes / GIB, "budget_bytes": budget_bytes,
                "budget_gib": budget_bytes / GIB})
    pol["exports"] = {CAP_ENV: fuse_gib_s, MODE_ENV: "fuse",
                      FUSE_GIB_ENV: fuse_gib_s, BUDGET_GIB_ENV: budget_gib_s}
    pol["arithmetic"] = (
        f"fuse = max(env cap {format_gib(cap_bytes) if cap_bytes else 'unset'} "
        f"GiB, {factor} x budget {budget_gib_s} GiB = {format_gib(scaled)} GiB)"
        f" = {fuse_gib_s} GiB = {fuse_bytes} bytes = {fuse_kb} kB (ulimit -v)")
    pol["words"] = (
        f"address-space fuse (RLIMIT_AS) = {fuse_gib_s} GiB = {factor} x the "
        f"memory budget {budget_gib_s} GiB; not the memory cap (under jemalloc "
        f"VmSize runs 2-5x RSS); ulimit -v {fuse_kb} in the launch line and "
        f"{CAP_ENV} exported as '{fuse_gib_s}' so the binary applies the same "
        f"fuse, not the 1x budget, to its Kira child; no finite memory.max on "
        f"cgroup {cg} (memory.max = {leaf.get('memory_max_raw')})")
    return pol


def ulimit_prefix(policy):
    """The bash prefix for the launch line: 'ulimit -v N; ' in fuse mode,
    '' otherwise."""
    kb = policy.get("ulimit_v_kb")
    return f"ulimit -v {kb}; " if kb else ""


def apply_policy_to_env(env, policy):
    """Write the policy's exports into an env dict (in place; returned).
    In 'unfenced' mode the cap field is removed so the child cannot inherit a
    stray 1x fuse from the launcher's own shell."""
    if policy["mode"] == "unfenced":
        env.pop(CAP_ENV, None)
        for k in (MODE_ENV, FUSE_GIB_ENV, BUDGET_GIB_ENV):
            env.pop(k, None)
        env[MODE_ENV] = "unfenced"
        return env
    env.update(policy["exports"])
    return env


def intended_from_env(env, extra_keys=()):
    """The cap fields (plus extra keys) as the launcher intends the child to
    see them; a None value means 'must be absent'."""
    keys = list(CAP_FIELDS) + [k for k in extra_keys if k not in CAP_FIELDS]
    return {k: env.get(k) for k in keys}


def read_environ(pid):
    """/proc/<pid>/environ as a dict. Raises OSError when unreadable; returns
    {} for a process that has exited (a zombie's environ is empty)."""
    with open(f"{PROC}/{pid}/environ", "rb") as fh:
        raw = fh.read()
    out = {}
    for item in raw.split(b"\0"):
        if not item:
            continue
        k, _, v = item.partition(b"=")
        out[k.decode("utf-8", "replace")] = v.decode("utf-8", "replace")
    return out


def _cmdline(pid):
    try:
        with open(f"{PROC}/{pid}/cmdline", "rb") as fh:
            return [c.decode("utf-8", "replace")
                    for c in fh.read().split(b"\0") if c]
    except OSError:
        return None


def new_launch_id():
    """A per-launch nonce the launcher puts into the child's environment as
    LAUNCH_ID_ENV. argv + cwd + a start-time window is NOT an identity: two
    launches of the same command from the same directory within the window
    resolve to the same pid, and an abort keyed on that would kill the other
    launch's child. The nonce, read back from /proc/<pid>/environ, is."""
    return uuid.uuid4().hex


def pid_env_matches(pid, require_env):
    """True when /proc/<pid>/environ carries every key=value in require_env
    (None value = key must be absent). False when unreadable or empty."""
    if not require_env:
        return True
    try:
        env = read_environ(pid)
    except OSError:
        return False
    if not env:
        return False
    return all(env.get(k) == v for k, v in require_env.items())


def env_readback(pid, intended, retries=25, interval=0.2):
    """Read the spawned child's real environment and compare the intended
    fields. Retries briefly until /proc/<pid>/environ is readable and
    non-empty (the child may still be between fork and exec).

    Returns {'pid', 'observed', 'intended', 'match', 'diff', 'verified',
             'error', 'read_utc', 'cmdline', 'child_cgroup_path'}.
    'match' is True only when every intended field was read back equal;
    'verified' is False when the environ could not be read at all (exited or
    foreign) — that is not a mismatch, and 'diff' is then empty."""
    rb = {"pid": pid, "intended": dict(intended), "observed": {}, "match": False,
          "diff": [], "verified": False, "error": None,
          "read_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
          "cmdline": None, "child_cgroup_path": None}
    if pid is None:
        rb["error"] = "no pid resolved"
        return rb
    env = None
    err = None
    for _ in range(max(1, int(retries))):
        try:
            env = read_environ(pid)
            if env:
                break
            err = "environ empty (process exited before it was read?)"
        except OSError as e:
            err = f"{type(e).__name__}: {e}"
        time.sleep(interval)
    if not env:
        rb["error"] = err or "environ empty"
        return rb
    rb["verified"] = True
    rb["cmdline"] = _cmdline(pid)
    rb["child_cgroup_path"] = own_cgroup_path(pid)
    for k, want in intended.items():
        got = env.get(k)
        rb["observed"][k] = got
        if got != want:
            rb["diff"].append({"key": k, "intended": want, "observed": got})
    rb["match"] = not rb["diff"]
    return rb


def _ppid(pid):
    st = _read_text(f"{PROC}/{pid}/stat")
    try:
        return int(st.rsplit(")", 1)[1].split()[1])
    except (AttributeError, ValueError, IndexError):
        return None


def descendants(pid):
    """Every live process whose ppid chain leads to <pid>, parents before
    children (breadth-first). Identity by lineage, read from /proc/<n>/stat —
    never a name or pattern match."""
    kids = {}
    for p in os.listdir(PROC):
        if p.isdigit():
            pp = _ppid(p)
            if pp is not None:
                kids.setdefault(pp, []).append(int(p))
    out, queue = [], [int(pid)]
    while queue:
        cur = queue.pop(0)
        for c in sorted(kids.get(cur, [])):
            out.append(c)
            queue.append(c)
    return out


def env_readback_tree(pid, intended, settle=0.75, retries=25, interval=0.2):
    """env_readback on the spawned child AND on every descendant present
    after `settle` seconds (a `timeout`/wrapper-script chain exec's the real
    binary a few ms after the spawn; a wrapper that drops or rewrites an
    export is visible only in the exec'd binary's environ). Returns the direct
    child's readback with 'tree': [per-descendant readbacks] and its 'diff'
    extended by every descendant's diff (annotated with the pid), so 'match'
    is True only when every process in the tree carries the intended env."""
    rb = env_readback(pid, intended, retries=retries, interval=interval)
    rb["tree"] = []
    if pid is None or not rb["verified"]:
        return rb
    time.sleep(settle)
    for d in descendants(pid):
        sub = env_readback(d, intended, retries=1, interval=0.0)
        rb["tree"].append(sub)
        for item in sub["diff"]:
            rb["diff"].append(dict(item, pid=d, cmdline=sub["cmdline"]))
    rb["match"] = not rb["diff"]
    return rb


def abort_tree_by_pid(pid, argv=None, reason="", require_env=None):
    """abort_by_name on the direct child (identity-checked against argv and
    require_env) and then on every descendant it still has, deepest first,
    each by the pid read from /proc lineage. Returns the list of abort
    records; if the direct child's identity check REFUSES, nothing is killed."""
    recs = []
    kids = descendants(pid) if pid is not None else []
    head = abort_by_name(pid, argv, reason=reason, require_env=require_env)
    recs.append(head)
    if head["refused"]:
        return recs
    for k in reversed(kids):
        recs.append(abort_by_name(k, None, reason=reason + " (descendant)"))
    return recs


def readback_line(rb, prefix="memfence readback"):
    """One JSON line for a log file."""
    return prefix + " " + json.dumps(
        {"pid": rb.get("pid"), "match": rb.get("match"),
         "verified": rb.get("verified"), "diff": rb.get("diff"),
         "observed": rb.get("observed"), "error": rb.get("error"),
         "read_utc": rb.get("read_utc")}, sort_keys=True)


def _alive(pid):
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    # a zombie answers signal 0; treat 'Z' state as gone
    st = _read_text(f"{PROC}/{pid}/stat")
    if st and ") Z" in st:
        return False
    return True


def abort_by_name(pid, argv=None, reason="", grace=2.0, poll=0.05,
                  require_env=None):
    """Kill the child the launcher itself resolved at spawn — by pid, never by
    pattern. When argv is given the pid's cmdline must equal it (or end with
    it, for a wrapper-exec'd child); when require_env is given the pid's
    environ must carry it (the launch nonce — see new_launch_id); a mismatch
    of either REFUSES the kill and says so. SIGTERM, wait up to `grace`, then
    SIGKILL. Also reaps the pid if it is our own child.
    Returns {'pid','killed','gone','refused','reason','cmdline'}."""
    rec = {"pid": pid, "killed": False, "gone": False, "refused": None,
           "reason": reason, "cmdline": _cmdline(pid) if pid else None,
           "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    if pid is None:
        rec["refused"] = "no pid"
        return rec
    if require_env and not pid_env_matches(pid, require_env):
        rec["refused"] = (f"identity mismatch: pid {pid} environ does not carry "
                          f"{require_env!r}; not killed")
        rec["gone"] = not _alive(pid)
        return rec
    if argv is not None:
        cl = rec["cmdline"]
        argv = list(argv)
        if cl is None:
            rec["refused"] = "cmdline unreadable (already gone?)"
            rec["gone"] = not _alive(pid)
            return rec
        if not (cl == argv or (len(cl) >= len(argv) and cl[-len(argv):] == argv)):
            rec["refused"] = (f"identity mismatch: pid {pid} cmdline {cl!r} is "
                              f"not the launched argv {argv!r}; not killed")
            return rec
    for sig in (signal.SIGTERM, signal.SIGKILL):
        try:
            os.kill(pid, sig)
            rec["killed"] = True
        except ProcessLookupError:
            rec["gone"] = True
            return rec
        deadline = time.time() + grace
        while time.time() < deadline:
            try:
                os.waitpid(pid, os.WNOHANG)  # reap if it is our child
            except ChildProcessError:
                pass
            if not _alive(pid):
                rec["gone"] = True
                return rec
            time.sleep(poll)
    rec["gone"] = not _alive(pid)
    return rec


def _boot_epoch():
    txt = _read_text(f"{PROC}/stat") or ""
    for ln in txt.splitlines():
        if ln.startswith("btime"):
            return int(ln.split()[1])
    return 0


def find_spawned_pid(argv_tail, t0, window=10.0, retries=25, interval=0.2,
                     cwd=None, require_env=None, exclude=()):
    """Resolve a child that a wrapper (setsid/nohup/a shell script) exec'd on
    the launcher's behalf: the process whose cmdline ENDS WITH argv_tail,
    whose start time is within `window` seconds of t0 (optionally with a
    matching cwd) and — when require_env is given — whose environ carries it
    (the launch nonce, new_launch_id). `exclude` lists tokens (the wrapper
    path, 'setsid') that must NOT appear in the cmdline, so the pre-exec
    stages of the chain — which still carry the launcher's env — are skipped
    and only the exec'd program (whose environ is what the wrapper actually
    passed on) is returned. Without the nonce, argv_tail must carry a
    per-launch unique element (an output path); a repeated launch of the same
    command within the window is otherwise ambiguous. Returns the pid or
    None."""
    tail = list(argv_tail)
    excl = set(exclude or ())
    hz = os.sysconf("SC_CLK_TCK")
    boot = _boot_epoch()
    cwd = os.path.realpath(cwd) if cwd else None
    for _ in range(max(1, int(retries))):
        for pid in os.listdir(PROC):
            if not pid.isdigit():
                continue
            cl = _cmdline(pid)
            if not cl or len(cl) < len(tail) or cl[-len(tail):] != tail:
                continue
            if excl and any(tok in excl for tok in cl):
                continue
            try:
                if cwd and os.path.realpath(f"{PROC}/{pid}/cwd") != cwd:
                    continue
                st = _read_text(f"{PROC}/{pid}/stat")
                start_ticks = int(st.rsplit(")", 1)[1].split()[19])
            except (OSError, ValueError, IndexError, AttributeError):
                continue
            if abs(boot + start_ticks / hz - t0) > window:
                continue
            if require_env and not pid_env_matches(pid, require_env):
                continue
            return int(pid)
        time.sleep(interval)
    return None


if __name__ == "__main__":
    import sys
    pol = fuse_policy(os.environ.get(CAP_ENV))
    json.dump({"leaf": leaf_memory_max(), "policy": pol}, sys.stdout, indent=1)
    print()

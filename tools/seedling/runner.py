"""Cgroup-wrapped Kira invocation + telemetry + the MANDATORY two-slice value gate.

Every stage/solve runs in a FRESH cgroup v2 with cpu.max and memory.max from config.
NEVER prlimit --as (VSZ-binding: silently wrong cap type — measured footgun).

Telemetry per run (JSON lines ledger, one record per phase):
  {phase, wall_s, peak_rss_bytes (cgroup memory.peak), oom_kills, exit_status,
   rundir, jobs_yaml_sha256, cmd}
The ledger is the supervised-training substrate for learned seed prediction.

TWO-SLICE DIFFERENTIAL GATE (mandatory after every pinned stage; measured
exhibit): a pinned system can drop super-sector-chain equations and still SOLVE
CLEANLY with wrong values (false-master termination) — leftover-target/Regenerate
signals do NOT catch it. Protocol (as executed live on the exhibit family):
  1. solve one slice on the pinned staging; compare row values against a stored
     independent reference slice at the SAME prime -> FAIL rows?
  2. if FAIL: solve a SECOND fresh slice at a DIFFERENT prime; compare again.
     - second slice FAILS identically  -> PIN AT FAULT (escalate: widen s first
       here, since chain loss is an s-phenomenon), quarantine all outputs;
     - second slice PASSES             -> first slice was a run artifact
       (stale-db/ff_save class); rerun slice 1 fresh.
Only a two-slice PASS unlocks farming.
"""
import hashlib
import json
import os
import re
import shlex
import subprocess
import sys
import time

GATE_PASS = "PASS"
GATE_FAIL_PIN = "FAIL_PIN"          # both slices fail: pinned system invalid
GATE_RERUN_SLICE = "RERUN_SLICE"    # first failed, second passed: slice artifact

# ---------------------------------------------------------------------------
# cgroup leaf launch — the cgroup.procs write discipline.
# Clause-by-clause statement in build_cgroup_script's docstring.
# ---------------------------------------------------------------------------
CG_ROOT = "/sys/fs/cgroup"
CG_NAME_RE = re.compile(r"[A-Za-z0-9_.-]+")
RC_CG_PATH_REFUSED = 97   # leaf script: $CG failed the path guard; nothing was written
RC_CG_NOT_ENTERED = 98    # leaf script: $$ not in the leaf after the write; moved back
RC_CG_FOREIGN_PID = 99    # launcher: a leaf pid does not descend from the child; aborted
_GUARD_NAMES = {RC_CG_PATH_REFUSED: "path_refused",
                RC_CG_NOT_ENTERED: "not_entered",
                RC_CG_FOREIGN_PID: "foreign_pid"}

# Shell guard emitted BEFORE the first write of the leaf script: an empty $CG
# aborts (`${CG:?}`), the path must be a proper child of the cgroup root, and
# no `.`/`..` component may resolve it elsewhere. Exit 97, nothing touched.
CG_PATH_GUARD = f"""\
: "${{CG:?}}"
case "$CG" in
  {CG_ROOT}/?*) ;;
  *) echo "seedling: cgroup path guard refused '$CG'" >&2; exit {RC_CG_PATH_REFUSED} ;;
esac
case "/$CG/" in
  */./*|*/../*) echo "seedling: cgroup path guard refused '$CG' (dot component)" >&2; exit {RC_CG_PATH_REFUSED} ;;
esac"""


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def validate_cgroup_name(name):
    """Leaf name the launcher accepts: non-empty, ^[A-Za-z0-9_.-]+$, and not
    made of dots alone ('.' or '..' would resolve outside the leaf).
    Returns the name; raises ValueError naming the defect."""
    if not isinstance(name, str) or not name:
        raise ValueError("cgroup name is empty")
    if not CG_NAME_RE.fullmatch(name):
        raise ValueError(f"cgroup name {name!r} is not ^[A-Za-z0-9_.-]+$")
    if not name.strip("."):
        raise ValueError(f"cgroup name {name!r} is a dot-only component")
    return name


def build_cgroup_script(cg_name, mem_max, cpu_max, cmd, cwd):
    """Return the launch script text (fresh cgroup v2, pid entered before exec).

    The script obeys the cgroup.procs write discipline, five clauses. For
    each clause: how a reference shell launcher does it / how this module
    does it:

    1. A write to cgroup.procs names ONLY a pid the writer itself forked,
       verified by the ppid chain — never a pgid, a pstree grep, a pattern
       match, or the $$ of a shell you did not create.
       Reference: `descends()` walks /proc/<pid>/stat field 4 up to the
       launcher shell; PID 1 in the leaf is fatal.
       Here: the script writes its own `$$` and nothing else; `run_phase`
       reads the leaf and requires every live pid's ppid chain to end at the
       child it forked (`descends_from`); any other pid aborts the launch.
    2. Enter the leaf FIRST and fork INSIDE it — children inherit; no
       post-hoc moves.
       Reference: `echo $ME | tee cgroup.procs`, then `grep /proc/$ME/cgroup`,
       then fork.
       Here: `$$` is written, then `grep -qx "$$" cgroup.procs` AND
       `/proc/self/cgroup` must name the leaf, then `exec` replaces the shell
       with the worker. Failure moves `$$` back to the root and exits 98.
    3. Print every moved pid with its cmdline in the launch line; after the
       move the leaf count must equal the spawn count; a mismatch ABORTS and
       empties the leaf.
       Reference: `LEAF_PROCS n (expected k)` with pid/ppid/cmdline per row;
       mismatch -> ABORT exit 7, every pid moved back to the root.
       Here: the script prints `LEAF <cgroup> pid <$$> cmd=<cmd>`; `run_phase`
       prints the leaf census (pid, ppid, cmdline, verdict) once the leaf is
       populated. The spawn count is one exec'd pid whose descendants are the
       worker's own helpers, so the count law takes the form "zero
       non-descendants": a foreign pid -> SIGKILL the child, move every leaf
       pid to the parent, ledger `phase=abort`, exit_status 99.
    4. sudo is the lever: a root-level leaf's cgroup.procs accepts any pid, so
       every root-leaf write is a machine-wide action.
       Reference: `sudo -n` only; an empty argument -> exit 3; leaf `/` or
       `..` -> exit 3.
       Here: the name is validated in Python before anything is built
       (`validate_cgroup_name`; rc 2 at the CLI), the script runs `set -eu`,
       `: "${CG:?}"` and the `/sys/fs/cgroup/?*` + no-dot-component path guard
       before its first write (exit 97), and `sudo -n` never prompts.
    5. timeout/setsid/nohup re-group children — the pgid is never the process
       you launched; clause 1 stands regardless.
       Reference: forks with `setsid nohup` and still verifies by ppid chain.
       Here: no pgid is ever read; the child is launched directly (Popen) and
       verified by ppid chain.

    Exit codes: 97 path guard refused, 98 not in the leaf after the write,
    99 foreign pid (module constants RC_CG_*).
    """
    validate_cgroup_name(cg_name)
    leaf_line = shlex.quote(f"cmd={cmd}")
    return f"""\
set -eu
CG={CG_ROOT}/{cg_name}
{CG_PATH_GUARD}
sudo -n mkdir -p "$CG"
echo {mem_max} | sudo -n tee "$CG/memory.max" >/dev/null
echo 0 | sudo -n tee "$CG/memory.swap.max" >/dev/null
echo "{cpu_max}" | sudo -n tee "$CG/cpu.max" >/dev/null
echo $$ | sudo -n tee "$CG/cgroup.procs" >/dev/null
if ! grep -qx "$$" "$CG/cgroup.procs" || ! grep -qx "0::/{cg_name}" /proc/self/cgroup; then
  echo "seedling: pid $$ not in leaf $CG after the write; moving it back to {CG_ROOT}" >&2
  echo $$ | sudo -n tee "{CG_ROOT}/cgroup.procs" >/dev/null || true
  exit {RC_CG_NOT_ENTERED}
fi
echo "LEAF $CG pid $$ "{leaf_line} >&2
cd {cwd}
exec {cmd}
"""


def _proc_ppid(pid):
    """Parent pid from /proc/<pid>/stat (field 4, after the comm) — the same
    source `ps -o ppid= -p <pid>` prints; None once the pid is gone."""
    try:
        with open(f"/proc/{pid}/stat") as fh:
            stat = fh.read()
    except OSError:
        return None
    try:
        return int(stat.rsplit(")", 1)[1].split()[1])
    except (IndexError, ValueError):
        return None


def _proc_cmdline(pid, width=120):
    try:
        with open(f"/proc/{pid}/cmdline", "rb") as fh:
            raw = fh.read()
    except OSError:
        return None
    return raw.replace(b"\0", b" ").decode("utf-8", "replace").strip()[:width]


def _proc_in_leaf(pid, cg_name):
    """True iff /proc/<pid>/cgroup names the leaf (unified line `0::/<name>`)."""
    try:
        with open(f"/proc/{pid}/cgroup") as fh:
            return any(line.strip() == f"0::/{cg_name}" for line in fh)
    except OSError:
        return False


def _read_leaf_pids(procs_path):
    try:
        with open(procs_path) as fh:
            return [int(tok) for tok in fh.read().split()]
    except (OSError, ValueError):
        return []


def _move_pid_to(pid, dest_dir):
    """Move one pid into dest_dir's cgroup.procs (one write per pid; the
    kernel accepts a single pid per write). True on success."""
    r = subprocess.run(["sudo", "-n", "tee", os.path.join(dest_dir, "cgroup.procs")],
                       input=f"{pid}\n".encode(), stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL)
    return r.returncode == 0


def descends_from(pid, ancestor, ppid_of=_proc_ppid, max_depth=64):
    """True iff the ppid chain of pid reaches ancestor before pid 0/1 or a
    vanished process (clause 1: the only proof a pid is one we forked)."""
    p = pid
    for _ in range(max_depth):
        if p == ancestor:
            return True
        if p is None or p in (0, 1):
            return False
        p = ppid_of(p)
    return False


# Injection points for run_phase's launch leg (tests replace them; nothing
# in the test battery touches a real cgroup or sudo).
DEFAULT_LAUNCH_OPS = {
    "cg_root": CG_ROOT,
    "launch": lambda script: subprocess.Popen(["bash", "-c", script]),
    "ppid_of": _proc_ppid,
    "cmdline_of": _proc_cmdline,
    "in_leaf": _proc_in_leaf,
    "kill": lambda proc: proc.kill(),          # SIGKILL, the child only
    "move": _move_pid_to,                      # (pid, dest_dir) -> bool
    "poll_s": 0.5,
    "echo": lambda line: print(line, file=sys.stderr, flush=True),
}


def leaf_census(cg_name, child_pid, ops):
    """Classify every pid listed in the leaf's cgroup.procs: `own` (its ppid
    chain ends at child_pid), `foreign` (it does not), or `gone` (it exited or
    left the leaf between the read and the lookup)."""
    procs = os.path.join(ops["cg_root"], cg_name, "cgroup.procs")
    rows = []
    for pid in _read_leaf_pids(procs):
        ppid = ops["ppid_of"](pid)
        if ppid is None or not ops["in_leaf"](pid, cg_name):
            rows.append({"pid": pid, "ppid": ppid, "cmd": None, "verdict": "gone"})
            continue
        own = descends_from(pid, child_pid, ops["ppid_of"])
        rows.append({"pid": pid, "ppid": ppid, "cmd": ops["cmdline_of"](pid),
                     "verdict": "own" if own else "foreign"})
    return rows


def _census_line(cg_name, child_pid, rows, ops):
    live = [r for r in rows if r["verdict"] != "gone"]
    n_own = sum(r["verdict"] == "own" for r in live)
    n_for = len(live) - n_own
    body = "; ".join(f"pid {r['pid']} ppid {r['ppid']} {r['verdict']} [{r['cmd']}]"
                     for r in live)
    return (f"LEAF-CENSUS {os.path.join(ops['cg_root'], cg_name)} child {child_pid} "
            f"n={len(live)} own={n_own} foreign={n_for}: {body}")


def _abort_launch(proc, cg_name, rows, launch_line, ops):
    """Clause 3 abort: SIGKILL the child we forked, move EVERY pid still in the
    leaf back to the parent (foreign pids are moved, never signalled), print
    each one with its cmdline, and report exit_status 99."""
    child = proc.pid
    leaf = os.path.join(ops["cg_root"], cg_name)
    foreign = [r for r in rows if r["verdict"] == "foreign"]
    desc = "; ".join(f"pid {r['pid']} ppid {r['ppid']} [{r['cmd']}]" for r in foreign)
    ops["echo"](f"ABORT leaf {leaf}: {len(foreign)} pid(s) do not descend from child "
                f"{child} — {desc}; SIGKILL child {child}, moving every leaf pid "
                f"back to {ops['cg_root']}")
    ops["kill"](proc)
    try:
        proc.wait(timeout=5)
    except subprocess.TimeoutExpired:
        pass
    moved = []
    known_cmd = {r["pid"]: r["cmd"] for r in rows if r["cmd"]}
    for pid in _read_leaf_pids(os.path.join(leaf, "cgroup.procs")):
        ok = ops["move"](pid, ops["cg_root"])
        cmd = ops["cmdline_of"](pid) or known_cmd.get(pid)   # census cmdline if it is gone
        moved.append({"pid": pid, "cmd": cmd, "moved": ok})
        ops["echo"](f"ABORT leaf {leaf}: moved pid {pid} [{moved[-1]['cmd']}] "
                    f"-> {ops['cg_root']}: {'ok' if ok else 'FAILED'}")
    left = _read_leaf_pids(os.path.join(leaf, "cgroup.procs"))
    ops["echo"](f"ABORT leaf {leaf}: {len(moved)} pid(s) moved, {len(left)} left in the leaf")
    abort = {"reason": "foreign pid in leaf (cgroup.procs write discipline, clauses 1/3)",
             "child_pid": child, "child_exit": proc.returncode,
             "foreign": foreign, "census": rows, "moved": moved,
             "left_in_leaf": left}
    return {"exit_status": RC_CG_FOREIGN_PID, "abort": abort,
            "leaf": {"launch_line": launch_line, "guard": _GUARD_NAMES[RC_CG_FOREIGN_PID],
                     "n_own": sum(r["verdict"] == "own" for r in rows),
                     "n_foreign": len(foreign)}}


def _watch_leaf(proc, cg_name, ops):
    """Poll the leaf while the child runs. The first populated census is the
    launch line (clause 3); a foreign pid at any poll aborts (clause 1). A
    census taken while the child was exiting is discarded: its leftovers are
    the child's own orphans, not an entry."""
    child = proc.pid
    launch_line = None
    n_own = 0
    while True:
        if proc.poll() is not None:
            break
        rows = leaf_census(cg_name, child, ops)
        if proc.poll() is not None:
            break
        live = [r for r in rows if r["verdict"] != "gone"]
        if live:
            n_own = max(n_own, sum(r["verdict"] == "own" for r in live))
            if launch_line is None:
                launch_line = _census_line(cg_name, child, rows, ops)
                ops["echo"](launch_line)
            if any(r["verdict"] == "foreign" for r in live):
                return _abort_launch(proc, cg_name, rows, launch_line, ops)
        time.sleep(ops["poll_s"])
    return {"exit_status": proc.returncode, "abort": None,
            "leaf": {"launch_line": launch_line,
                     "guard": _GUARD_NAMES.get(proc.returncode),
                     "n_own": n_own, "n_foreign": 0}}


def run_phase(phase, rundir, cmd, cgroup_cfg, ledger_path, dry_run=False,
              launch_ops=None):
    """Run one phase in a fresh cgroup; append a ledger record; return it.

    cgroup_cfg: dict(name=, mem_max='12G', cpu_max='400000 100000').
    dry_run=True: build everything, execute nothing (returns record with
    exit_status=None) — used by tests and by budget-blocked planning.
    launch_ops: overrides for DEFAULT_LAUNCH_OPS (test injection only).

    Launch leg (see build_cgroup_script for the law it implements): the child
    is started with Popen — never `timeout`/`setsid`/`nohup` — and while it
    runs the leaf's cgroup.procs is read; every live pid must descend from the
    child by ppid chain. A foreign pid -> SIGKILL the child, every leaf pid
    moved back to the parent, a `phase=abort` ledger record, exit_status 99.
    The record carries `leaf` (launch line, guard name, own/foreign counts)
    and, on abort, `abort` (foreign pids with cmdlines, moved pids).
    """
    script = build_cgroup_script(cgroup_cfg["name"], cgroup_cfg["mem_max"],
                                 cgroup_cfg["cpu_max"], cmd, rundir)
    jy = os.path.join(rundir, "jobs.yaml")
    rec = {
        "phase": phase, "rundir": rundir, "cmd": cmd,
        "jobs_yaml_sha256": _sha256(jy) if os.path.exists(jy) else None,
        "cgroup": dict(cgroup_cfg), "dry_run": dry_run,
        "wall_s": None, "peak_rss_bytes": None, "oom_kills": None,
        "exit_status": None, "ts": time.time(),
    }
    if not dry_run:
        ops = dict(DEFAULT_LAUNCH_OPS)
        ops.update(launch_ops or {})
        t0 = time.time()
        proc = ops["launch"](script)
        watch = _watch_leaf(proc, cgroup_cfg["name"], ops)
        rec["wall_s"] = time.time() - t0
        rec["exit_status"] = watch["exit_status"]
        rec["leaf"] = watch["leaf"]
        if watch["abort"]:
            rec["abort"] = watch["abort"]
            if ledger_path:
                with open(ledger_path, "a") as fh:
                    fh.write(json.dumps({"phase": "abort", "rundir": rundir,
                                         "cmd": cmd, "cgroup": dict(cgroup_cfg),
                                         "ts": time.time(), **watch["abort"]})
                             + "\n")
        cg = os.path.join(ops["cg_root"], cgroup_cfg["name"])
        for key, field in (("memory.peak", "peak_rss_bytes"),):
            try:
                rec[field] = int(open(os.path.join(cg, key)).read())
            except OSError:
                pass
        try:
            ev = dict(l.split() for l in open(os.path.join(cg, "memory.events")))
            rec["oom_kills"] = int(ev.get("oom_kill", 0))
        except OSError:
            pass
    if ledger_path:
        with open(ledger_path, "a") as fh:
            fh.write(json.dumps(rec) + "\n")
    return rec, script


def compare_slice(fresh_rows, banked_rows):
    """Row-value comparison at one prime. Inputs: dict row_key -> value (int/str).
    Returns (n_common, fail_keys) — fail = present in both, values differ."""
    common = set(fresh_rows) & set(banked_rows)
    fails = sorted(k for k in common if fresh_rows[k] != banked_rows[k])
    return len(common), fails


def two_slice_gate(slice1, banked1, slice2=None, banked2=None):
    """Mandatory post-stage gate. slice1/banked1: row dicts at prime p1;
    slice2/banked2: row dicts at a DIFFERENT prime (only needed if slice1 fails).

    Returns dict(verdict=GATE_*, fails1=[...], fails2=[...]|None).
    """
    n1, fails1 = compare_slice(slice1, banked1)
    if n1 == 0:
        raise ValueError("gate has no rows in common with the reference")
    if not fails1:
        return {"verdict": GATE_PASS, "fails1": [], "fails2": None}
    if slice2 is None or banked2 is None:
        return {"verdict": "NEED_SECOND_SLICE", "fails1": fails1, "fails2": None}
    _, fails2 = compare_slice(slice2, banked2)
    if fails2:
        return {"verdict": GATE_FAIL_PIN, "fails1": fails1, "fails2": fails2}
    return {"verdict": GATE_RERUN_SLICE, "fails1": fails1, "fails2": []}


# ---------------------------------------------------------------------------
# DE-census gate (measured protocol lesson): overlap gates MISS
# DE-feeding rows. The eta-DE of an n-master system consumes the first
# eta-derivative images {nu + e_i} of every master nu (i over eta-carrying
# propagators); for ISP masters (negative indices) these images keep s>=1 —
# exactly the pin's zero-margin direction, invisible to stage-overlap gates.
# Measured census on the reference 54-master system: 317 rows = 276
# overlap-gated + 33 master-self + 8 uncovered (certified from-scratch in
# 451 s / 4.7G with seeds ONE UNIT past target support in every direction).
# Protocol: BEFORE consuming a pinned table in a DE solve, compute the census,
# partition it, and certify the uncovered subset with a from-scratch mini-kira.
# ---------------------------------------------------------------------------

def de_row_census(masters, eta_props):
    """First eta-derivative row images of the master set.

    masters: iterable of index tuples; eta_props: 0-based positions of
    eta-carrying propagators. Returns set of index tuples.
    """
    census = set()
    for nu in masters:
        for i in eta_props:
            img = list(nu)
            img[i] += 1
            census.add(tuple(img))
    return census


def de_census_gate(masters, eta_props, gated_rows):
    """Partition the DE census into covered/uncovered.

    gated_rows: rows already value-gated (e.g. stage-overlap set).
    Returns dict(census_n=, gated_n=, master_self_n=, uncovered=[tuples]).
    The 'uncovered' subset MUST be certified from-scratch before the DE solve.
    """
    masters_set = {tuple(m) for m in masters}
    gated = {tuple(g) for g in gated_rows}
    census = de_row_census(masters_set, eta_props)
    master_self = census & masters_set
    covered = (census & gated) | master_self
    uncovered = sorted(census - covered)
    return {
        "census_n": len(census),
        "gated_n": len(census & gated),
        "master_self_n": len(master_self),
        "uncovered": uncovered,
    }


def minikira_spec(family, top_sectors, uncovered_rows, margins=(1, 1, 1)):
    """Build the from-scratch certification job for the uncovered subset:
    target file text + jobs.yaml text, seeds pinned to the uncovered rows'
    support plus one unit in EVERY direction (validated recipe; ~451 s class on the
    reference family). Run with numeric eta at a generic point, no FireFly.
    """
    from . import support as _support, pinner as _pinner
    targets = [_support.Target(family, tuple(r)) for r in uncovered_rows]
    if not targets:
        raise ValueError("no uncovered rows — census gate already satisfied")
    g = _support.global_support(targets)
    r, s, d = g["r"] + margins[0], g["s"] + margins[1], g["d"] + margins[2]
    target_text = "\n".join(
        f"{family}[{', '.join(str(i) for i in t.indices)}]" for t in targets) + "\n"
    jobs_text = _pinner.render_jobs_yaml(family, top_sectors, r, s, d)
    return {"target": target_text, "jobs_yaml": jobs_text, "cell": (r, s, d)}

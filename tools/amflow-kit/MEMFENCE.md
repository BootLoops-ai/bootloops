# amflow-kit member manual: memfence — memory fence policy + environment readback for the AMFlow/Kira launchers

`amflow_kit.memfence` (`tools/amflow-kit/amflow_kit/memfence.py`), a member of
amflow-kit, the house kit around the AMFlow.cpp fork (overview: `GUIDE.md`).
One helper, imported by every launcher that spawns `amflow_cli` (and so Kira,
Fermat, FireFly under it): `pmflow/pmflow.py` and `numkin/shift_opt.py` in this
repository, and any launcher of your own that follows the same pattern. Stdlib
only; Linux `/proc` and cgroup v2.

## The defect it replaces

`AMFLOW_KIRA_MEM_CAP_GB` was turned into a fixed address-space limit
(`ulimit -v`, RLIMIT_AS) sized as the memory budget. Under jemalloc the
virtual set of a Kira reduction runs 2–5x its resident set, so a limit "sized
as the budget" kills the job at about a third of the memory it was meant to
have: `std::bad_alloc` at 23.8 GB RSS under a 48 GiB cap inside a 64 GiB
cgroup leaf whose `memory.events` showed `oom 0`. A second copy of the same fuse lives in the binary:
`amflow-cpp/src/ibp/kira_run.cpp` reads the variable and calls
`setrlimit(RLIMIT_AS, N GiB)` on the Kira child it forks. That binary is not
edited here; the helper drives it through the value it exports.

What the binary does with the variable (read from the source, not assumed):

| exported value | binary's action on its Kira child |
|---|---|
| unset | no `setrlimit` ("Unset (default) => no cap") |
| `""` (empty) | no `setrlimit` (`cap && *cap` fails) |
| `0`, `-3`, `abc`, ` ` | no `setrlimit` (`strtod` gives 0 / no digits; only `gib > 0.0` sets it) |
| `48` | `RLIMIT_AS = 48 GiB` (soft = hard) |
| `37.5` | `RLIMIT_AS = 37.5 GiB` (`strtod`, so fractions are fine) |

So `'0'` is the explicit "no fuse" export, and any value the helper exports in
fuse mode is applied by the binary as the Kira child's RLIMIT_AS too — which
is why the fuse value, never the 1x budget, is what gets exported.

## Policy — `fuse_policy(env_cap_gb, budget_bytes=None, leaf=None)`

| where the launcher runs | mode | `ulimit -v` | export `AMFLOW_KIRA_MEM_CAP_GB` |
|---|---|---|---|
| cgroup leaf with finite `memory.max` | `leaf` | none — the leaf is the fence | `'0'` (binary: no cap) |
| `memory.max = max` (or unreadable), budget known | `fuse` | `max(env cap, 2.5 x budget)` in kB | the fuse in GiB (`'120'`, `'37.5'`) |
| no fence and no budget at all | `unfenced` | none | removed (nothing to fence) |

The budget is the caller's intended memory budget; today the env cap IS the
budget at every call site, so the fuse is 2.5 x the env cap: 48 → 120 GiB =
125,829,120 kB; 15 → 37.5 GiB = 39,321,600 kB; 60 (the pmflow
default) → 150 GiB. The launch line and every receipt name it an
**address-space fuse (RLIMIT_AS) = N GiB = 2.5 x the memory budget M GiB; not
the memory cap** — `policy['words']` carries that sentence, `policy['arithmetic']`
the numbers.

The decision reads the launcher's **own leaf's** `memory.max` (the path in
`/proc/self/cgroup`, `0::<path>`), not an ancestor's. Ancestors are recorded
in `leaf_memory_max()['ancestors']` for the receipt but do not decide: on a
shared machine the user slice above every login carries a finite `memory.max`
equal to the whole machine, and treating that as "inside a leaf" would remove
the fuse from every launch on the machine — the runaway-symbolic-swell case the
fuse exists for. Enter the leaf first, then launch (write the pid into cgroup.procs before the spawn);
the launcher's leaf is then the child's leaf.

Beside the cap the helper exports three inert twins so the child's environment
states its fence: `MEMFENCE_MODE`, `MEMFENCE_FUSE_GIB`, `MEMFENCE_BUDGET_GIB`.

## Readback — `env_readback(pid, intended)` / `env_readback_tree`

After the spawn the launcher reads `/proc/<pid>/environ` of the child it
resolved and compares the intended fields (`intended_from_env(env)` = the cap
fields, plus whatever the launcher injected). `match` is True only when every
field reads back equal; `verified` is False when the environ could not be read
at all (exited, foreign) — that is reported, not treated as a mismatch.
`env_readback_tree` also reads every descendant present after a short settle
(`timeout` → jemalloc wrapper → the exec'd `amflow_cli`), because a wrapper
that drops or rewrites an export is visible only in the exec'd binary's
environ.

A mismatch means the job is running under a different fence than the receipt
says (the classic cause: an `export` line swallowed by a comment in a launch
script — only `/proc` shows it). The launcher then **aborts by name**:
`abort_by_name(pid, argv)` kills the pid it resolved at spawn after checking
that pid's cmdline is the argv it launched (a wrong identity REFUSES the
kill), SIGTERM then SIGKILL; `abort_tree_by_pid` does the same for the
descendants, deepest first, each by the pid read from `/proc` lineage. Never a
pattern, never `pkill -f`.

`find_spawned_pid(argv_tail, t0)` resolves a child that a wrapper exec'd on
the launcher's behalf (the `setsid` fork-parent footgun): the process whose
cmdline ends with `argv_tail` and started within the window. The tail must
carry a per-launch unique element (an output path) so it is an identity.

## Call sites and their receipt surfaces

| launcher | policy applied in | readback recorded in | on mismatch |
|---|---|---|---|
| `pmflow/pmflow.py` | `_env` (every subcommand's env) | `<log>.memfence.json` beside the run's log (no other per-run receipt) | `_run_fenced` kills the tree, RuntimeError |
| `numkin/shift_opt.py` | `probe` | `<log>.memfence.json` beside the probe log | kills the resolved amflow pid, RuntimeError |

A launcher of your own applies the policy where it builds the child's
environment (the `ulimit -v` prefix plus the exports), records the readback in
its launch receipt (`fence_mode`, `ulimit_v_kb`, `address_space_fuse_words`,
`fence_arithmetic`, `leaf_memory_max_bytes`, `cgroup_path`, `launch_line`,
`env_readback`, `aborted`), and on a mismatch kills the pid it resolved and
exits non-zero.

Import line used by the call sites (works in the repo and in a mirror tree laid
out as `tools/<pkg>/`; the name bound is `memfence`, so call sites read
`memfence.fuse_policy(...)`):

```python
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "amflow-kit")); from amflow_kit import memfence
```

`python3 -m amflow_kit.memfence` (from `tools/amflow-kit`) prints the calling
shell's own leaf record and the policy it would get, as JSON.

## Battery

The member's legs are part of the kit's one-command battery (`GUIDE.md`,
Battery); alone, from `tools/amflow-kit`:

```
python3 -m pytest tests/test_memfence.py tests/test_cgroup_procs_size.py -q -rs -p no:cacheprovider
```

Planted controls, every one measured on a real process where a process is
involved: leaf fixture (`memory.max = 68719476736` → leaf, no ulimit, export
`'0'`), max fixture (48 → 120 GiB named "address-space fuse", 15 → 37.5 GiB),
the binary's no-cap parse table, the max rule, a real `sleep 30` reading back
equal, a `# export AMFLOW_X=1` comment-swallowed export reading back as a
mismatch and the abort killing exactly that pid, a wrong-identity abort
refusing, a wrapper that drops the cap caught only by the tree readback, and
the unique-tail resolver.

`tests/test_cgroup_procs_size.py` adds two LIVE legs (no monkeypatch) that
pin why `_read_text` reads by content and never by size: cgroup v2 interface
files stat at size 0 whatever they hold. `test_own_scope_memory_max_is_size_zero_with_parsable_text`
— the calling process's own `memory.max` has `os.path.getsize(...) == 0`
while `leaf_memory_max` parses its text (`max` or an integer);
`test_live_leaf_cgroup_procs_size_zero_with_N_pids` — a real leaf made under
the calling process's own cgroup, holding 3 `sleep` children the test
spawned (one `cgroup.procs` write per pid, never its own pid), reports
`getsize(cgroup.procs) == 0` while the read lists exactly those 3 pids and
`own_cgroup_path(child)` names the leaf; so a `test -s cgroup.procs`
liveness check is shown inert beside the read-based one, and a kernel that
ever reports a non-zero size fails both legs by name. No sudo; the leaf is
torn down in `finally` (terminate + wait by pid, rmdir). Where the kernel
refuses the mkdir or the write the leg SKIPS by name with the errno. Skip
contract: the battery line above (and the kit's registered BATTERIES.json line) is
plain pytest, so a skip leaves the exit status 0 and is named only under
`-rs` (`-q -rs`); run_selftests.py then reports the kit PASS with the
live leg skipped — there is no separate runner that returns 2.

## Footguns

- The ancestors' `memory.max` never decides (see above); if a job wants the
  leaf to be the fence, the launcher must run inside the leaf.
- `RLIMIT_AS` is per process: the launcher's `ulimit -v` binds `amflow_cli`
  and everything it forks; the binary's own `setrlimit` binds only the Kira
  child. Both now carry the same fuse value.
- A readback that could not be read (`verified: False`) is not a mismatch; a
  launcher that needs a hard guarantee must resolve a long-lived pid (the
  `find_real_pid` pattern) — `env_readback(None, ...)` says "no pid resolved".

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).

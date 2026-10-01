# amflow-kit — the house kit around the AMFlow.cpp fork (lint, smoke, compare, keypred, memfence)

NAME: `amflow-kit` (directory `tools/amflow-kit/`, Python import name `amflow_kit`)
is the kit that goes around the AMFlow.cpp fork: everything this repository runs
before, beside and after an `amflow_cli` call that is not the solver itself —
three gates on the solver's output, an offline predictor of the fork's IBP-cache
key, and the memory fence every launcher applies to the solver process. The fork
itself lives in the sibling repository amflow-cpp (clone it beside this one: ../amflow-cpp); this
kit builds nothing and its public battery needs no binary.

Tool page: https://www.bootloops.ai/tools/amflow.html (the AMFlow page; the kit
rides the AMFlow entry, with the Kira page https://bootloops.ai/tools/kira-stack.html
for the cache-key member).

KIND: kit = three command-line gates (flat scripts in this directory:
`amflow_output_lint.py`, `amflow_smoke.sh`, `compare_amflow_json.py`) + the
importable package `amflow_kit/` with two stdlib-only modules, each also a
command line: `amflow_kit.keypred` and `amflow_kit.memfence`. Member manuals:
`KEYPRED.md`, `MEMFENCE.md`; the comparator's rules: `README.md`.

PURPOSE: One home for the checks and guards around `amflow_cli solve_integrals`:
validate every output before it is consumed (lint), validate every rebuilt
binary before it is used (smoke), gate every solver-lineage change on matched
digits (compare), predict where the fork's IBP cache will look before a run is
staged (keypred), and fence the solver's memory the same way from every
launcher, with a readback that proves the child got the fence the receipt
claims (memfence).

## Output gates — lint, smoke, A/B compare

- `amflow_output_lint.py out_*.json [--class vacuum|kinematic] [--allow-shallow]
  [--min-nonzero-frac F] [--in job.json] [-q]` — per-file sanity: nonzero-coef
  fraction, Laurent depth >=2, vacuum lead<0, arb-ball parse; loud FAIL + nonzero
  exit. Lint EVERY output before consuming.
- `amflow_smoke.sh BINARY_PATH [JOB_TIMEOUT_SECONDS]` — post-rebuild smoke: 2
  known-answer jobs (2L vacuum probe + 1L triangle), lint + >=20-digit compare vs
  fixtures, per-binary IBP cache isolation. Run after EVERY rebuild (catches
  double-precision Frobenius-resonance underflows).
- `compare_amflow_json.py ref.json test.json [--min-digits 25] [--n-integrals N]
  [--goal G|GREF,GTEST | --goal-from-config REF_CFG[,TEST_CFG] | --radii-only]
  [--strict-window]` — order-by-order matched-digit A/B gate for solver-lineage
  changes. A pair of coefficients both below the zero threshold
  T = max(radii, 10^-(goal-2)) is a NUMERICAL ZERO (named, excluded, counted); an
  order present in one file only is "not compared" (the eps window a solve returns
  depends on its goal), a FAIL only under `--strict-window`; OVERALL carries the two
  counts. The goal has three sources, in order: `--goal`; `--goal-from-config`
  reading the amflow_cli input config's top-level `goal_digits` (or `goal`); with
  neither flag the file's own key, else the config found BY NAME beside the out
  (`out/<name>.json` <-> `configs/<name>.json`), every location tried named on the
  `goal:` line. No derivable goal = a refusal by name (exit 2, nothing compared);
  `--radii-only` is the explicit opt-in to the radii-alone threshold (the previous
  behaviour); a missing `--goal-from-config` path exits 4 by name. The
  `solve_integrals` output form carries no goal key.

Fixtures ship in `../fixtures/amflow_smoke/` (see that README for reconstruction
provenance). NOTE: fixture `in_vac2Bprobe.json` carries `eps_parallel` keys in both `options`
and `blackbox`, so hygiene-checking binaries print two config warnings during
smoke — expected, harmless, and doubles as a live warning regression.
`amflow_smoke.sh` calls `amflow_output_lint.py` from this same directory
(resolved via `BASH_SOURCE`), so the two stay together.

## keypred — offline IBP-cache key predictor (`amflow_kit.keypred`)

PURPOSE: Predicts, without running anything, the 16-hex key under which the
fork's IBP cache will store or look up the reduction of a Kira input directory.
It hashes the CANONICAL `jobs.yaml` form — the same form the fork's own
`ibp_cache_key` hashes (thread-count lines dropped, the write-time `, d: N` seed
cap stripped, the result-determining `d:` of Masters mode kept) — so a
prediction made from either byte convention of `jobs.yaml` lands on the key the
live run computes. A raw-byte prediction from one convention misses keys
computed from the other; that is the miss class this member removes. Numeric-d
Reduce legs salt their key with the exact string `AMFLOW_DIFFEQ_NUMD=P/Q`;
`--numd` builds that salt, `--salt=` folds any string verbatim.

INVOKE (from this directory): `python3 -m amflow_kit.keypred key DIR [--raw]
[--explain] [--salt=S | --numd[=P/Q]]` (the canonical key; `--raw` = the
byte-exact key, for forensics only); `python3 -m amflow_kit.keypred predict DIR
[...]` (a hand-built mirror dir: CRLF / missing-final-newline repaired onto the
fork's writer convention first, then the canonical key); `python3 -m
amflow_kit.keypred --selftest` (<1 s). As a library: `from amflow_kit import
keypred; keypred.ibp_cache_key(DIR, raw=False, pin=False, salt="")`,
`keypred.numd_salt("P/Q")`. Inputs: a Kira input dir
(`config/{integralfamilies,kinematics}.yaml`, `jobs.yaml`, `jobs_kira2math.yaml`,
`preferred`, `target`; missing files skipped exactly as the fork does). Output:
the key on stdout; `--explain` adds the pin repairs and the canonical bytes as
hashed on stderr.

LIMITS: removes the byte-convention miss class only. A staged directory whose
CONTENT differs from what the fork writes (other sectors or rank box, another
reduction mode, another target set) still keys differently, correctly; read the
key the live run logged in that case. Never build a cache entry from a `--raw`
key. Never salt a Masters-mode prediction (its result does not depend on d, so
numeric-d and symbolic runs share the entry on purpose). The canonicalizer
tracks the fork's `ibp_cache.cpp` byte for byte: when that file, `kira_yaml.cpp`,
`reduce.cpp` or `numd_subst.cpp` changes, re-run the selftest in the same cycle.
Convention, fixtures and footguns: `KEYPRED.md`.

## memfence — launch memory fence + environment readback (`amflow_kit.memfence`)

PURPOSE: One policy, shared by every launcher that spawns `amflow_cli` (and so
Kira, Fermat and FireFly under it): `pmflow/pmflow.py`, `numkin/shift_opt.py`,
and any launcher of your own. It replaces a defect: the memory budget
`AMFLOW_KIRA_MEM_CAP_GB` used to be applied as an address-space limit
(`ulimit -v`, RLIMIT_AS) of the same size, and under jemalloc a Kira
reduction's virtual size runs 2-5x its resident size, so the job died at about
a third of the memory it was meant to have. The policy: inside a cgroup v2 leaf
with a finite `memory.max` the leaf is the fence (no `ulimit -v`; the variable
exported as `'0'`, which the fork reads as "no cap"); with no finite leaf limit
an address-space FUSE of 2.5x the budget is set and exported, and every receipt
names it "address-space fuse = N GiB = 2.5 x the memory budget; not the memory
cap"; with no leaf and no budget nothing is fenced. After the spawn the
launcher reads `/proc/<pid>/environ` of the child (and of every descendant a
wrapper exec'd) and compares it with what it intended; on a mismatch it kills
exactly the pid it resolved, after an identity check, never by pattern.

INVOKE (library): `from amflow_kit import memfence` (from a sibling package:
`sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
"..", "amflow-kit")); from amflow_kit import memfence`). `pol =
memfence.fuse_policy(env.get(memfence.CAP_ENV))`;
`memfence.apply_policy_to_env(env, pol)`; prefix the launch line with
`memfence.ulimit_prefix(pol)`; after the spawn `rb =
memfence.env_readback_tree(pid, memfence.intended_from_env(env))`; on `not
rb["match"]`: `memfence.abort_tree_by_pid(pid, argv, reason=...)` and exit
non-zero. `python3 -m amflow_kit.memfence` prints the calling shell's own leaf
record and the policy it would get, as JSON. Stdlib only; Linux `/proc` and
cgroup v2.

LIMITS: the decision reads the launcher's OWN leaf, never an ancestor's (on a
shared machine the slice above every login carries a finite limit equal to the
whole machine; treating that as a leaf would remove the fuse everywhere), so a
job that wants the leaf as its fence must be launched from inside the leaf.
RLIMIT_AS is per process: the launcher's `ulimit -v` binds `amflow_cli` and its
children, the fork's own `setrlimit` binds only its Kira child; both now carry
the same fuse value. A readback that could not be read (`verified: False`, the
child already gone) is reported, not treated as a mismatch. Policy tables, the
fork's parse table for the variable, call sites and receipt fields:
`MEMFENCE.md`.

## Configuration

The smoke script's paths are env-overridable; defaults assume this repo layout and a
scratch root under `$TMPDIR`:

- `AMFLOW_SMOKE_FIX` — fixtures dir (default: `../fixtures/amflow_smoke` beside the kit)
- `AMFLOW_SMOKE_ROOT` — scratch root for ref/cache/run dirs (default `$TMPDIR/amflow_smoke`;
  `AMFLOW_SMOKE_REF_DIR` / `AMFLOW_SMOKE_CACHE_ROOT` / `AMFLOW_SMOKE_RUN_ROOT` split it)
- `AMFLOW_SMOKE_VENDOR_BIN` — reference amflow_cli binary, needed only the first time
  the triangle reference is generated
- `AMFLOW_SMOKE_WRAP` — optional launch wrapper (e.g. a jemalloc wrapper); unset = run the binary directly
- `FERMATPATH` — must point at a `fer64` binary for amflow_cli; inherited from your
  environment

The two library members read the fork's own variables and add none that a user
must set: `AMFLOW_DIFFEQ_NUMD` (keypred's bare `--numd` reads it, as the live run
does) and `AMFLOW_KIRA_MEM_CAP_GB` (memfence reads it as the memory budget and
re-exports the fence value; it also exports the informational
`MEMFENCE_MODE` / `MEMFENCE_FUSE_GIB` / `MEMFENCE_BUDGET_GIB` and, when a
launcher asks, the per-launch nonce `MEMFENCE_LAUNCH_ID`).

Build the binary from the amflow-cpp fork (the sibling repository amflow-cpp).

## Battery

One command from this directory tests the whole kit:

```
python3 selftest.py && python3 -m pytest tests -q -rs -p no:cacheprovider
```

`python3 selftest.py` (<2 s, stdlib only) — lint of the shipped vacuum fixture
PASS (5 ints, 43/44 nonzero, class=vacuum); compare ref-vs-ref OVERALL PASS
(the goal via `--goal-from-config` from the fixture's own job config
`in_vac2Bprobe.json`, key `goal`); the smoke script resolves the shipped
fixtures and fails closed with a named FATAL when no vendor binary is
configured; `amflow_kit.keypred --selftest` 16/16 (5 fixture positive controls,
4 mutation predicts, 6 salt legs, the numd format gate; fixtures in
`tests/fixtures/keypred/`, scratch copies under `$TMPDIR`); `import amflow_kit`
exposes both members and memfence's fuse arithmetic on a planted no-leaf record
reads 48 -> 120 GiB. The full smoke (two solver jobs) needs a built
`amflow_cli`.

`python3 -m pytest tests -q -rs -p no:cacheprovider` (~10 s; pytest, stdlib +
mpmath) — `test_compare_amflow_json.py`: the compare gate's two rules on the
pinned reference fixture pairs (goals 50 / 70 with its two re-cut configs; goals
30 / 60), the three goal sources and the refusal, and synthetic pairs (the
earlier-tool positive control skips by name unless
`COMPARE_AMFLOW_JSON_PREVIOUS` names its file); `test_memfence.py`: the leaf and
no-leaf fixtures, the fork's no-cap parse table, the max rule, a real `sleep`
child reading back equal, a comment-swallowed export reading back as a mismatch
and the abort killing exactly that pid, a wrong-identity abort refusing, a
wrapper that drops the cap caught only by the tree readback, the unique-tail
resolver; `test_cgroup_procs_size.py`: two live legs pinning why cgroup v2
interface files are read by content and never by size (they SKIP by name with
the errno where the kernel refuses the leaf mkdir or the `cgroup.procs` write;
the exit status stays 0 and `-rs` names the skip).

Source of truth for the gates' flags: the kit's own `README.md`.

CREDIT: these gates validate output of AMFlow.cpp (the AMFlow.cpp contributors,
maintainer GitHub @chang18; github.com/chang18/amflow-cpp, doi:10.5281/zenodo.20087172;
our patched fork: the sibling repository amflow-cpp), a C++17 reimplementation of AMFlow by
Xiao Liu and Yan-Qing Ma [AMFlow] implementing the auxiliary-mass-flow method of Liu, Ma
and Chen-Yu Wang [AMF1, AMF2]. All algorithmic credit for the solver is theirs; the
lint/smoke/compare kit is ours. The ibp-cache key convention that the keypred member
predicts is AMFlow.cpp's (MIT), and the job files it canonicalizes are Kira's (P.
Maierhöfer, J. Usovitsch & P. Uwer [Kira1]; J. Klappert, F. Lange, P. Maierhöfer & J.
Usovitsch [Kira2]; F. Lange, J. Usovitsch & Z. Wu [Kira3]); keypred only predicts
hashes, the design is theirs. Bracketed keys resolve in REFERENCES.md at the repository
root.

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).

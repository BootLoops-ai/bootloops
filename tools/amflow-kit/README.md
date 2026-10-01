# amflow-kit — the house kit around the AMFlow.cpp fork

One home for what this repository runs around `amflow_cli` (the AMFlow.cpp fork,
the sibling repository amflow-cpp): three gates that validate
`solve_integrals` output, the offline IBP-cache key predictor, and the launch
memory fence. Directory `tools/amflow-kit/`; the two library members import as
`amflow_kit.keypred` and `amflow_kit.memfence`.

| member | what it does | invoke |
|--------|--------------|--------|
| `amflow_output_lint.py` | per-file sanity: nonzero-coef fraction, Laurent depth ≥2, vacuum lead<0, arb-ball parse; loud FAIL + nonzero exit | `python3 amflow_output_lint.py out_*.json [--class vacuum\|kinematic] [--allow-shallow] [--min-nonzero-frac F] [--in job.json] [-q]` |
| `amflow_smoke.sh` | post-rebuild smoke: 2 known-answer jobs through a binary, lint + ≥20-digit compare vs fixtures, per-binary IBP cache | `amflow_smoke.sh BINARY_PATH [JOB_TIMEOUT_SECONDS]` |
| `compare_amflow_json.py` | A/B two outputs order-by-order, matched relative decimal digits; a NUMERICAL ZERO pair (both members below T = max(radii, 10^-(goal-2))) is named, excluded and counted; an order present in one file only is "not compared" (the eps window is goal-dependent), never a FAIL unless `--strict-window`; the goal from `--goal`, from the amflow_cli config (`--goal-from-config`, or found by name beside the out), else the run refuses by name | `python3 compare_amflow_json.py REF.json TEST.json [--min-digits 25] [--n-integrals N] [--goal G\|GREF,GTEST \| --goal-from-config REF_CFG[,TEST_CFG] \| --radii-only] [--strict-window]` |
| `amflow_kit/keypred.py` | offline IBP-cache key predictor: hashes the CANONICAL `jobs.yaml` form the fork's own `ibp_cache_key` hashes (thread lines dropped, the write-time `, d: N` seed cap stripped, the Masters-mode `d:` kept), so a prediction from either byte convention lands on the live key; `predict` repairs CRLF / missing final newline on a hand-built mirror first; `--numd` / `--salt=` fold the numeric-d salt; `--raw` = the byte-exact key, forensics only | `python3 -m amflow_kit.keypred key DIR [--raw] [--explain] [--salt=S \| --numd[=P/Q]]`; `... predict DIR [...]`; `... --selftest`; library `from amflow_kit import keypred` |
| `amflow_kit/memfence.py` | shared memory-fence policy for the AMFlow/Kira launchers: inside a finite-`memory.max` cgroup leaf no `ulimit -v` (the leaf is the fence, `AMFLOW_KIRA_MEM_CAP_GB` exported as `'0'`), otherwise an address-space fuse of 2.5x the budget named as such; post-spawn `/proc/<pid>/environ` readback of the child and its exec'd descendants, abort by pid (identity-checked) on a mismatch | library: `from amflow_kit import memfence`; `fuse_policy`, `apply_policy_to_env`, `ulimit_prefix`, `env_readback_tree`, `abort_tree_by_pid`; `python3 -m amflow_kit.memfence` prints the calling shell's leaf record and policy |

`amflow_smoke.sh` calls `amflow_output_lint.py` from this same directory (resolved
via `BASH_SOURCE`), so the two stay together. The smoke fixtures are at `tools/fixtures/amflow_smoke/`;
keypred's synthetic Kira input dirs at `tests/fixtures/keypred/`. Member manuals: `KEYPRED.md`
(convention, footguns), `MEMFENCE.md` (policy tables, readback, call sites); overview and the
one-command battery: `GUIDE.md`.

## compare_amflow_json.py: what is compared

Every coefficient is read as an arb ball `[mid +/- rad]` (a plain number or an
empty ball has radius 0); the digit rule runs on the midpoints as before. Two
rules decide which components enter the digit count; everything else is named
and counted, never failed:

- **NUMERICAL ZERO.** A pair whose two members are both below the zero threshold
  `T = max(rad_ref, rad_test, 10^-(G-2))`, G the lower of the two goals, is
  printed by name (`int0 eps^-4 re: NUMERICAL ZERO (|ref| 2.7e-100, |test| 1.8e-142
  below T 1.0e-48) -> not compared`), excluded from the digit count and counted on
  the line `numerical zeros: N (int0 eps^-4 re, ...)`. A pair with exactly one
  member above T is zero-vs-nonzero and FAILS with its digit count printed (as
  before); a pair whose two midpoints are exactly 0 is `both exactly 0 -> PASS`
  (as before).
- **THE GOAL** has three sources, tried in this order; the `goal:` line names the
  one used. (1) `--goal G` (one value for both files) or `--goal GREF,GTEST`,
  which wins over every other source. (2) `--goal-from-config REF_CFG[,TEST_CFG]`:
  the amflow_cli input config(s) (one path applies to both files); the goal is
  the config's top-level `goal_digits` (a positive number; the input form of
  record), or `goal` when a config carries that key instead (the smoke fixtures'
  job configs); the key read is named. A config path that does not exist exits 4
  by name (`no config: <path> (REF, --goal-from-config) does not exist`); a config
  without either key, or with a non-numeric one, exits 2 by name (`no goal: config
  <path> carries none of goal_digits, goal at its top level`). (3) With neither
  flag: a file carrying `goal_digits` or `goal` at its top level or under
  `options` is read from there (the key named); the `solve_integrals` output form
  carries none of these keys (its `options` carry `working_pre` / `chop_pre` /
  `rationalize_pre`, working precisions, never read as the goal), so for such a
  file the config is looked up BY NAME beside it — `out/<name>.json` <->
  `configs/<name>.json`, i.e. `../configs/<name>.json` relative to the out file's
  directory, read as in (2) — and every location tried is named on the `goal:`
  line (`tried by name: <path> (TEST, absent)`). When no goal is derivable on
  either side the run REFUSES by name and compares nothing: exit 2,
  `no goal: pass --goal G|GREF,GTEST, --goal-from-config REF_CFG[,TEST_CFG], or
  --radii-only`. `--radii-only` is the explicit opt-in to the threshold from the
  two radii alone (the previous behaviour without a goal: the `goal:` line reads
  `goal: none (zero threshold = the two radii alone)` and a finite integral's
  pole-order zeros read as `0.0 digits -> FAIL`, as before); it cannot be combined
  with `--goal` or `--goal-from-config`.
- **MISSING BY WINDOW.** An order present in one file only is printed
  `int0 eps^-6: not compared (present at goal 50 (REF) only; outside the other
  member's window)` and counted on `orders not compared: N (...)`, never a FAIL.
  `--strict-window` restores the previous rule (`only one side, but zero -> PASS`
  / `MISSING on one side, nonzero -> FAIL`).
- **OVERALL** reads PASS / FAIL over the compared components (a component = the
  re or the im part of one order) with the counts beside it:
  `OVERALL: PASS (10 compared, 4 numerical zeros excluded, 2 orders not compared)`.
  Exit 0 iff every compared component is >= `--min-digits` (and, under
  `--strict-window`, no order is missing with a nonzero member); 1 otherwise; 2
  for a refused option, a malformed `--goal`, a config without a usable goal key
  or no derivable goal; 4 for a `--goal-from-config` path that does not exist.
  `--min-digits` and `--n-integrals` are unchanged.

Battery: `python3 selftest.py && python3 -m pytest tests -q -rs -p no:cacheprovider`
from this directory tests the whole kit (GUIDE.md, Battery). The comparator's legs
are `tests/test_compare_amflow_json.py` (pytest, stdlib + mpmath; fixtures
sha256-pinned in `tests/fixtures/PINS.json`).
The two amflow_cli configs of the reference pair ship beside it as
`pair_goal50.config.json` / `pair_goal70.config.json` (goal_digits 50 / 70 kept;
their machine paths and working names re-cut to `<re-cut>`, the re-cut fields
listed in PINS.json): the pair reads the same OVERALL PASS via
`--goal-from-config`, and by name with no flag on a copy laid out as
`out/<name>.json` + `configs/<name>.json`; with no flag and no config beside the
fixtures it refuses by name (rc 2, nothing compared); `--radii-only` reproduces
the previous no-goal reading `pair_radii_only_previous.out` byte for byte.
`selftest.py` compares the shipped vacuum fixture with itself under
`--goal-from-config` naming the fixture's own job config `in_vac2Bprobe.json`
(its key is `goal`).
The reference fixture pair is a finite eight-propagator integral solved at goals
50 and 70: its pole-order coefficients are numerical zeros and its eps window
differs by goal, so an earlier version of the tool without the two rules read it
OVERALL FAIL (that reading is the fixture `pair_compare_previous.out`) and the
current tool reads it OVERALL PASS with 4 numerical zeros and 2 orders not
compared. Set `COMPARE_AMFLOW_JSON_PREVIOUS=<path of the earlier tool>` to run
that positive control (skipped by name otherwise).
A second fixture pair, `pair_q1b_goal30.json` / `pair_q1b_goal60.json` (98
integrals of one family solved at goals 30 and 60), covers the imaginary parts:
46 im components of real coefficients are numerical zeros (|im| 1.7e-72 .. 9.2e-58 at
goal 30, 1.7e-132 .. 9.2e-108 at goal 60) that a radii-only comparison reads as `0.0
digits -> FAIL` (OVERALL FAIL over the 46); the tool with `--goal 30,60` reads
`OVERALL: PASS (1358 compared, 46 numerical zeros excluded, 0 orders not
compared)`.

Configuration and battery notes: `GUIDE.md` in this directory.

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).

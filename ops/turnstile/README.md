# Turnstile — admission control for long jobs on a shared machine

`ops/turnstile` is operations infrastructure, not a scientific instrument: the
one door every long BootLoops job walks through before it starts on a machine
that other jobs share. Behind the door sit three ledgers —

- a **priority token** (drop-in replacement for raw `flock $LOCK <cmd>`: the
  highest-priority live waiter takes the token, not whoever the kernel wakes
  first),
- a **RAM ledger** with a mechanized aggregate launch gate for the smaller
  "exempt" jobs that do not need the token (peak-RSS receipts become the
  pricing priors for the next launch), and
- a **CPU-width ledger** over the machine's reserved CPU span (declared worker
  widths gate co-fires; measured widths audit them).

Failure-mode classes it exists to prevent:

- **FIFO inversion.** flock's kernel queue is FIFO-ish; a multi-cycle
  top-priority chain releases at each cycle boundary and requeues at the BACK,
  losing every inter-cycle race to cheaper jobs queued ahead of it.
- **Simultaneous ramps.** Several sub-cap "exempt" ramps fired at once can
  crash available memory in under a minute, taking down an unrelated
  long-running job. Per-job exemption with no aggregate accounting.
- **TOCTOU.** `flock -n LOCK true` then launch-unheld: a probe is not a hold.
- Stale wrappers/entries from killed chains must never wedge the queue.

## Usage

```sh
# >300G-class big solve (takes the token, held across the ENTIRE run):
ops/turnstile/turnstile_run.sh --lane solveA_big --class critical --est-rss 400 --est-threads 16 --est-wall 12000 -- \
    nice -n 5 python3 -u solve_big.py ...

# sub-300G exempt job (no token; registers in the aggregate ledger, gated launch):
turnstile_run.sh --exempt --lane jobB_ramp --est-rss 60 --est-threads 8 -- python3 grow_table.py ...

# sub-1-thread helper (watcher/receipt writer): skips the width ledger entirely
turnstile_run.sh --exempt --lane jobB_watcher --est-rss 1 --est-threads 0 -- ./watch.sh

# numeric priority instead of a class name:
turnstile_run.sh --lane jobC --prio 40 --est-rss 520 -- ...

turnstile_status.sh       # read-only: holder, waiters in priority order, ledgers, gate arithmetic
bash selftest.sh          # scratch-lock battery (never touches the real lock)
```

`--lane <name>` names the job: the label under which it is queued, ledgered,
shown by `turnstile_status.sh` and matched by `PRIORITIES.conf` globs.
`--est-rss G` is **required** (caps-from-measurement: price from `LEDGER.log` priors or a measured
sibling datum). Exit 75 = refused/timed out (EX_TEMPFAIL); exit 69 = platform
refusal (below); otherwise the wrapped command's rc is propagated.
`TURNSTILE_HELD=1` (token) / `TURNSTILE_EXEMPT=1` (exempt), `TURNSTILE_LANE` and
`TURNSTILE_EST_THREADS` are set in the child's environment.

## Platform: Linux only

The wrapper is built on Linux-only tooling: `taskset` pin spans, `nproc`,
`/proc/meminfo` (the MemAvailable gate), `flock(1)`, and GNU `ps` per-thread
listing (`ps -eLo`). None of these exist on stock macOS. On darwin,
`turnstile_run.sh` and `turnstile_status.sh` refuse by name with exit 69
(EX_UNAVAILABLE) — `Linux-only: taskset unavailable on darwin ...` — rather
than mis-gate or half-run, and `selftest.sh` prints a named per-leg SKIP and
exits 0. There is no macOS port.

## Priority ladder (class names are CLI-stable)

| class | prio |
|---|---|
| critical | 10 |
| high | 20 |
| normal | 30 |
| low | 40 |
| misc (default) | 50 |

Lower number wins; ties by queue time (FIFO within a class). The **operator
overrides** priorities live via `<QDIR>/PRIORITIES.conf` (`<lane-glob> <prio>`,
first match wins, `default N` allowed) — re-read at every poll, so a re-ladder
takes effect without touching waiters.

## Mechanics

Queue dir `<lock>.q/` beside the lock (e.g. `/var/tmp/TURNSTILE.q/`):
`waiters/` (token queue), `exempt/` (aggregate ledger of running sub-300G jobs),
`width/` (pin-width ledger of FIRED jobs' est_threads), `HOLDER` (current
token holder record), `LEDGER.log` (append-only receipts: done/stale_pruned/
width_reserve/width_drift events with est vs measured peak — the pricing prior
for future `--est-rss`/`--est-threads`), plus operator-owned
`PRIORITIES.conf`, `HOLD`, `COTENANT_RESERVE_G`.

1. **Priority, not FIFO.** Waiters do NOT block on the kernel flock. Each polls
   (default 15 s): under a short queue lock it prunes dead entries, computes the
   best live waiter (effective priority, then age), and ONLY the best one attempts
   `flock -n`. A lower-priority waiter never takes a free token while a
   higher-priority waiter is registered — regardless of arrival order or how many
   cycle boundaries the leader crosses.
2. **TOCTOU-safe.** The winning wrapper takes the real kernel flock and spawns the
   job with the lock fd inherited: the token is held by the job's own fd tree
   across the entire run (even if the wrapper is SIGKILLed, the lock survives with
   the job — never test-then-launch, never an unheld big solve).
3. **Exempt aggregate gate (mechanized).** Before an exempt job fires:
   `MemAvailable >= own_est + SUM(running exempts' remaining ramp: est - current
   pgid RSS) + token-holder remaining ramp + COTENANT_RESERVE_G + floor(150G)`,
   plus one-ramp-at-a-time stagger (a peer younger than `TURNSTILE_RAMP_SECS`=300 s
   and below 90% of its est blocks new fires; `--no-stagger` to disable the
   stagger clause only). Token jobs run the same sum gate after winning the token,
   holding it while exempts drain — the mechanized EXEMPT-HOLD drain.
4. **Crash-safe.** Entries carry the wrapper pid + job pgid; any wrapper prunes
   entries whose pid AND pgid are dead (prune is ledgered). A stale `HOLDER`
   record is reaped only when a `flock -n` probe proves the kernel lock is
   actually free. The token itself is a kernel flock — released by the kernel
   when the holder dies, no lease files to expire.
5. **Operator levers** (files in `<lock>.q/`, no code):
   - `PRIORITIES.conf` — live re-ladder (above).
   - `HOLD` — presence blocks ALL new launches, token and exempt
     (EXEMPT-HOLD; content = reason, echoed in logs). Running jobs unaffected.
   - `COTENANT_RESERVE_G` — one number: protected co-tenant phase reserve
     (e.g. a co-tenant job's projected memory swell). Counted in every gate.

## Width ledger

The pin (`taskset` span; default = the whole machine) is a physical resource
like RAM: ~140 runnable threads on a 64-core pin (~2.2x oversubscribed) slow
every job on it and force contention pricing into every declared cap. The
width ledger makes CPU-width a ledgered, gated resource:

- **Declare:** `--est-threads N` (measured worker width). Undeclared -> WARN +
  default to the FULL pin width (`nproc` of the taskset span): undeclared jobs
  SERIALIZE on the pin. Declare a measured width to co-fire.
- **Ledger:** `width/` beside `waiters/`/`exempt/` under the same `.qlock`;
  entries `{lane, pgid, est_threads, ts, peak_r}` written at FIRE time, removed
  at exit, pruned by the same pid+session-scoped-pgid liveness as the RAM
  ledger.
- **Gate:** a new fire (token AND exempt — width is physical, exempt legs
  block too) waits while `SUM(est_threads of ALIVE entries) + own > cap`.
  Reservation happens under the qlock that passed the gate — no fire-fire race.
- **Only skip class:** `--est-threads 0` (sub-1-thread helpers: watchers,
  receipt writers). `--no-width-gate` is the operator escape — it still
  REGISTERS width, it only skips the wait. `--no-gate` bypasses the RAM gate
  only, never width.
- **Measured, not declared:** the monitor samples the pgid's live R-thread
  count (session-scoped, like `pg_rss_gib`); `peak_r` lands in the entry and
  the `done` ledger line; live > 2x declared ledgers a throttled `width_drift`
  event — the signal to re-measure the job's declared width.
- `width_reserve` ledger lines carry `pin_r_at_fire` (live pin R census at
  fire time) for auditing co-fires against unmanaged jobs.

## Env knobs

`TURNSTILE_LOCK` (default `/var/tmp/TURNSTILE.lock`; every launcher on a machine
must agree on it), `TURNSTILE_POLL` (15 s, integer), `TURNSTILE_MON` (peak-RSS sample,
30 s), `TURNSTILE_FLOOR_G` (150), `TURNSTILE_RAMP_SECS` (300),
`TURNSTILE_EXEMPT_MAX_G` (300; `--exempt` with est >= this refuses unless
`TURNSTILE_FORCE_EXEMPT=1` with an explicit override), `TURNSTILE_PIN`
(taskset span; default = the whole machine), `TURNSTILE_WIDTH_CAP` (default `nproc`
of the span), `TURNSTILE_DRIFT_LOG_SECS` (width_drift throttle, 600).
Flags: `--no-wait` (refuse instead of poll, rc 75), `--timeout S`, `--no-gate`
(RAM launch gate off — operator use only; width still gates), `--est-threads N`
(0 = sub-1-thread helper, skips width), `--no-width-gate` (register width,
skip the wait — operator use only), `--lock PATH`, `--floor G`, `--poll S`.

**Compatibility (former name).** This package was previously `tools/bigram-queue`
with scripts `bigram_run.sh` / `bigram_status.sh` and `BIGRAM_*` knobs. The
scripts honor every `BIGRAM_<K>` as a fallback when `TURNSTILE_<K>` is unset, and
the lock, queue-dir, entry and ledger formats are unchanged, so old and new
launch lines on one machine interoperate; migrate each launch line at its next
launch.

## Coexistence with raw flock

The wrapper holds the SAME kernel flock as raw `flock` users, so wrapped and
legacy jobs mutually exclude correctly. Caveats until every launcher on the
machine uses the wrapper: a legacy blocked `flock` waiter is woken by the kernel at release and
beats the wrapper's poll, and legacy holders show as "LEGACY raw-flock job" in
`turnstile_status.sh` (no est/peak record). Full priority + ledger semantics arrive
as launchers migrate — **each launch line replaces its `flock $LOCK cmd` with
`turnstile_run.sh ... -- cmd` at its NEXT launch; no live-job churn.**

Known limits: poll granularity leaves the token free up to ~`TURNSTILE_POLL` seconds
at handoff (negligible vs multi-hour solves); `kill -9` of wrapper+job together
leaves the entry until the next prune (pid-liveness, first poll of any peer);
priorities arbitrate the queue, they never preempt a running holder (parking a
holder stays an operator judgment call).

Deployment note: where launch lines call a deployed COPY of `turnstile_run.sh`
(a regular file, not a symlink), patch it BY RENAME only — running wrappers
read the script lazily and hold the old inode. Selftest: 23/23 (RAM legs,
width battery, 3-launcher 30/30/30 pilot, daemon-reap leg — the battery
terminates and joins every process it spawned before exiting, so no stray
holds an inherited output pipe past exit); `TURNSTILE_TEST_DIR=<dir>` redirects
its scratch. Width
adoption is per job stream, by operator decision: existing launch lines keep
working (undeclared = WARN + serialize); a job co-fires only once it declares a
measured `--est-threads`.

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).

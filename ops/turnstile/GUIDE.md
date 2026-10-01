# Turnstile — admission control for long jobs on a shared machine

Turnstile — `ops/turnstile` — https://bootloops.ai/tools/turnstile.html

KIND: operations wrapper (bash, no deps beyond flock/ps/awk/taskset). Infrastructure
for running BootLoops on a shared machine, not a scientific instrument — which is
why it lives under `ops/`, beside `tools/`, and outside the thematic tool roster.

PURPOSE: One door every long job walks through before it starts. Three ledgers
behind it: a PRIORITY TOKEN (replacement for raw `flock LOCK <cmd>`; the
highest-priority live waiter takes the token, not the first in the kernel queue),
a shared RAM ledger with a mechanized exempt-aggregate launch gate (peak-RSS
receipts), and a pin-global CPU-WIDTH ledger (declared `est_threads` gating
co-fires on the machine's reserved CPU span). Prevents three measured
scheduling-failure classes on a shared machine: FIFO inversion (a multi-cycle
chain requeues at the back at every cycle boundary), aggregate RAM surges from
individually-"exempt" jobs, and TOCTOU probe-then-launch races (a `flock -n`
probe is not a hold).

USE-WHEN:

- Launching ANY large-RSS job of the token class on a shared machine — use instead of
  raw flock.
- Launching smaller exempt work that could stack into an aggregate surge.
- Operator needs to re-ladder the waiting queue (PRIORITIES.conf), declare
  EXEMPT-HOLD (HOLD file), or reserve a protected co-tenant phase
  (COTENANT_RESERVE_G) without touching running jobs.

NOT-FOR: preempting a RUNNING holder (priorities arbitrate the queue only; parking a
holder stays an operator call). Not a cross-machine scheduler — locks/ledger are
machine-local. Not a memory cap: it admits or defers a launch, it never limits a
running process (address-space caps are a different door, e.g.
`tools/dipstick/run_msolve_capped.sh` for msolve).

INVOKE:

- token: `ops/turnstile/turnstile_run.sh --lane <name> --class <class> --est-rss G
  --est-threads N [--est-wall S] -- <cmd>` (or `--prio N`, lower wins)
- exempt: `turnstile_run.sh --exempt --lane <name> --est-rss G --est-threads N -- <cmd>`
- helper (sub-1-thread watcher/receipt writer): `--est-threads 0` (only class that
  skips the width ledger)
- view: `turnstile_status.sh` (shows the width ledger + live pin R-count); test:
  `bash ops/turnstile/selftest.sh` (scratch lock, safe anywhere)

RULES:

- `--est-rss` REQUIRED (price from LEDGER.log priors / a measured sibling).
- Token held by the JOB's inherited fd across the entire run — TOCTOU-safe by
  construction; never `flock -n ... true` then launch.
- Exempt gate: avail >= own_est + SUM(running exempts' remaining ramp) + holder
  remaining ramp + COTENANT_RESERVE_G + floor, plus one-ramp-at-a-time stagger
  (TURNSTILE_RAMP_SECS=300). Token jobs run the same sum gate holding the token.
- Stale entries expire by pid+pgid liveness at every poll; a stale HOLDER is reaped
  only after a `flock -n` proof.
- WIDTH: the fire gate waits while SUM(est_threads of ALIVE width entries)
  + own > pin cap. Width is PHYSICAL — exempt legs block too; `--no-gate` bypasses
  RAM only. Undeclared `--est-threads` WARNs + defaults to full pin width
  (serializes). live R > 2x declared => width_drift ledger event.

OUTPUTS: rc = wrapped cmd's rc; 75 = refused/timeout (EX_TEMPFAIL); 69 = platform
refusal (darwin). LEDGER.log append-only receipts (est vs measured peak, wall, rc,
stale prunes). Child env: TURNSTILE_HELD=1 / TURNSTILE_EXEMPT=1, TURNSTILE_LANE,
TURNSTILE_EST_THREADS.

ENV: TURNSTILE_LOCK (default `/var/tmp/TURNSTILE.lock`; every launcher on a machine
must agree on it), TURNSTILE_POLL=15 (integer s), TURNSTILE_MON=30,
TURNSTILE_FLOOR_G=150, TURNSTILE_RAMP_SECS=300, TURNSTILE_EXEMPT_MAX_G=300
(+TURNSTILE_FORCE_EXEMPT=1 override, operator-overridden only), TURNSTILE_PIN
(taskset span; default = the whole machine — set it to the reserved span on a
shared box), TURNSTILE_WIDTH_CAP (default nproc of the span),
TURNSTILE_DRIFT_LOG_SECS=600. Compatibility: the wrapper's former name was
`bigram_run.sh` with `BIGRAM_*` knobs; each `BIGRAM_<K>` is honored as a fallback
when `TURNSTILE_<K>` is unset (same lock, queue-dir and ledger formats), so a
machine can migrate launch lines one at a time.

INPUTS: operator files in `<lock>.q/`: PRIORITIES.conf ("<lane-glob> <prio>",
first match; "default N"), HOLD (presence blocks all new launches),
COTENANT_RESERVE_G (single number, GiB).

GATES: selftest 23/23 PASS, rc=0, run from a scratch cwd (hold-across-run,
3-waiter priority order, PRIORITIES.conf override, stale dead-pid prune, floor
refusal, sum-of-ramps refusal, HOLD refusal, RAM legs, width battery, 3-launcher
30/30/30 pilot 11/11, daemon-reap leg). The width legs pin TURNSTILE_WIDTH_CAP=64
explicitly, so the selftest is machine-independent. `TURNSTILE_TEST_DIR=<dir>`
redirects the scratch (default a fresh `mktemp -d` under `$TMPDIR`).

FOOTGUNS:

- Mixed use with raw flock: legacy raw-flock waiters are woken by the kernel and beat the
  wrapper's poll; legacy holders show as "LEGACY raw-flock job" in status (mutual
  exclusion still correct — same kernel lock). Full semantics apply once every launcher on the machine uses the wrapper.
- Token sits free up to ~TURNSTILE_POLL s at handoff (negligible vs multi-hour solves).
- kill -9 of the wrapper alone does NOT drop the token (the job's fd tree holds it —
  intentional); kill -9 of both leaves an entry until any peer's next prune.
- TURNSTILE_POLL must be an integer (jitter is appended as a decimal).
- Where launch lines call a deployed COPY of `turnstile_run.sh` (a regular file),
  patch it BY RENAME only — running wrappers read the script lazily and hold the
  old inode.
- The width ledger only sees wrapper-managed launches: already-running unmanaged
  jobs are invisible to the width sum.

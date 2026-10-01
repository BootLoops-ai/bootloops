#!/usr/bin/env bash
# turnstile_run.sh — Turnstile: admission control for long jobs on a shared machine.
#                    PRIORITY TOKEN + exempt aggregate RAM ledger + PIN-GLOBAL CPU-WIDTH ledger.
# Replaces raw `flock $LOCK <cmd>` for large solves (the token) and mechanizes the
# exempt-aggregate launch gate. See README.md / GUIDE.md beside this file.
#
#   Token (>300G-class):  turnstile_run.sh --lane solveA_big --class critical --est-rss 400 --est-threads 16 -- <cmd>
#   Exempt (sub-300G):    turnstile_run.sh --exempt --lane jobB_ramp --est-rss 60 --est-threads 8 -- <cmd>
#   (--lane names the job: the label under which it is queued, ledgered and matched by PRIORITIES.conf)
#
# WIDTH LEDGER: the pin (taskset span, e.g. 0-63 = 64
# CPUs) is a physical resource like RAM. A width ledger ($QDIR/width/, same
# qlock) tracks {lane,pgid,est_threads,ts} of every FIRED job; a new fire WAITS
# while sum(est_threads of ALIVE entries) + own > cap (64). Liveness is the same
# session-scoped pg_alive; dead entries are pruned like the RAM ledger's.
# Width is physical: EXEMPT legs block too. Only --est-threads 0 (sub-1-thread
# helpers) skips ledger+gate. Undeclared --est-threads defaults to the FULL pin
# width (nproc of the taskset span) with a WARN — undeclared jobs serialize on
# the pin; declare a measured width to co-fire. Drift >2x (live R-count vs
# declared) is ledgered for cap re-typing; declared widths gate, measured
# widths audit.
#
# Failure-mode classes this exists to prevent:
#   FIFO inversion   — the top-priority job requeues at the BACK every cycle
#                      and loses every race to cheaper jobs.
#   Simultaneous ramps — several exempt jobs ramping at once crash available
#                      memory in seconds, taking down an unrelated
#                      long-running job.
#   TOCTOU           — `flock -n LOCK true` then launch-unheld: the lock is
#                      tested, dropped, and the job launches unprotected.
#
# Guarantees:
#   (1) PRIORITY, not FIFO: waiters register in a queue dir; when the token frees,
#       only the highest-priority live waiter takes it (operator override file).
#   (2) EXEMPT LEDGER + LAUNCH GATE: avail >= own_est + sum(running exempts'
#       remaining ramp) + token-holder remaining ramp + co-tenant reserve + floor,
#       plus one-ramp-at-a-time stagger. Refuse (wait) otherwise.
#   (3) TOCTOU-safe: the kernel flock is taken by THIS process and the fd is
#       inherited by the job — held across the entire run, never test-then-launch.
#   (4) Crash-safe: queue entries expire by pid-liveness; the token is a kernel
#       flock, auto-released when the holder's fd tree dies; stale HOLDER records
#       are reaped by a flock -n probe.
#
# Env knobs are TURNSTILE_* (README "Env knobs"). Compatibility: this wrapper was
# formerly deployed as `bigram_run.sh` with BIGRAM_* knobs; every BIGRAM_<K> is
# honored as a fallback when TURNSTILE_<K> is unset, so old and new launch lines on
# one machine keep agreeing on the lock, queue dir and ledger (same formats).
set -uo pipefail

# ------------------------------------------------ legacy env-name fallback ----
for _k in LOCK POLL MON FLOOR_G RAMP_SECS EXEMPT_MAX_G PIN WIDTH_CAP DRIFT_LOG_SECS FORCE_EXEMPT; do
  _t="TURNSTILE_$_k"; _b="BIGRAM_$_k"
  if [ -z "${!_t:-}" ] && [ -n "${!_b:-}" ]; then printf -v "$_t" '%s' "${!_b}"; fi
done
unset _k _t _b

# ---------------------------------------------------------- platform gate ----
# The wrapper is built on Linux-only tooling: taskset pin spans, nproc,
# /proc/meminfo (MemAvailable gate), flock(1), and GNU ps per-thread listing
# (ps -eLo pgid,sid,stat,psr). None of these exist on stock macOS, so on darwin
# we refuse by NAME rather than mis-gate or half-run. There is no macOS port.
if [ "$(uname -s)" = Darwin ]; then
  echo "[turnstile] Linux-only: taskset unavailable on darwin (also required: nproc, /proc/meminfo, flock(1), ps -eLo thread listing). Refusing rather than running ungated — use a Linux host." >&2
  exit 69   # EX_UNAVAILABLE: named platform refusal, distinct from 75 (wait-refused) and 70 (fatal)
fi

# ---------------------------------------------------------------- defaults ---
# One lock file per machine: every launcher on the machine must agree on the path
# (set TURNSTILE_LOCK to a shared, persistent location; the default suits a
# single-machine deployment).
DEF_LOCK=/var/tmp/TURNSTILE.lock
LOCK="${TURNSTILE_LOCK:-$DEF_LOCK}"
POLL="${TURNSTILE_POLL:-15}"          # seconds between queue polls
MON="${TURNSTILE_MON:-30}"            # seconds between peak-RSS ledger updates
FLOOR_G="${TURNSTILE_FLOOR_G:-150}"   # machine memory floor; selftest overrides
RAMP_SECS="${TURNSTILE_RAMP_SECS:-300}"   # a running exempt younger than this and <90% est = "ramping" (stagger clause)
EXEMPT_MAX_G="${TURNSTILE_EXEMPT_MAX_G:-300}"  # token-class boundary
PIN_SPAN="${TURNSTILE_PIN:-0-$(( $(nproc) - 1 ))}"   # taskset span all turnstile jobs are pinned to (set TURNSTILE_PIN to the reserved span; default = whole machine)

MODE=token LANE="" PRIO="" CLASS="" EST_RSS="" EST_WALL=0 NOWAIT=0 NOGATE=0 NOSTAGGER=0 TIMEOUT=0
EST_THREADS="" NOWIDTHGATE=0 WIDTH_SKIP=0 MY_WIDTH="" FIRED_TS="" PEAK_R=0
usage() { sed -n '2,47p' "$0"; exit 2; }
while [ $# -gt 0 ]; do
  case "$1" in
    --exempt)     MODE=exempt; shift ;;
    --lane)       LANE="$2"; shift 2 ;;
    --prio)       PRIO="$2"; shift 2 ;;
    --class)      CLASS="$2"; shift 2 ;;
    --est-rss)    EST_RSS="$2"; shift 2 ;;
    --est-threads) EST_THREADS="$2"; shift 2 ;;
    --est-wall)   EST_WALL="$2"; shift 2 ;;
    --no-width-gate) NOWIDTHGATE=1; shift ;;   # operator escape: still REGISTERS width, skips only the wait
    --lock)       LOCK="$2"; shift 2 ;;
    --floor)      FLOOR_G="$2"; shift 2 ;;
    --poll)       POLL="$2"; shift 2 ;;
    --timeout)    TIMEOUT="$2"; shift 2 ;;
    --no-wait)    NOWAIT=1; shift ;;
    --no-gate)    NOGATE=1; shift ;;
    --no-stagger) NOSTAGGER=1; shift ;;
    --) shift; break ;;
    -h|--help) usage ;;
    *) echo "[turnstile] unknown arg: $1" >&2; usage ;;
  esac
done
[ $# -gt 0 ] || { echo "[turnstile] no command after --" >&2; exit 2; }
CMD=( "$@" )

# Priority ladder (class names are CLI-stable): critical > high > normal > low > misc
class_prio() { case "$1" in
  critical) echo 10 ;; high) echo 20 ;; normal) echo 30 ;;
  low) echo 40 ;; misc) echo 50 ;; *) echo "" ;; esac; }
if [ -z "$PRIO" ] && [ -n "$CLASS" ]; then PRIO=$(class_prio "$CLASS")
  [ -n "$PRIO" ] || { echo "[turnstile] unknown --class '$CLASS' (critical|high|normal|low|misc)" >&2; exit 2; }
fi
[ -n "$PRIO" ] || PRIO=50
[ -n "$LANE" ] || LANE="lane.$$"
if [ -z "$EST_RSS" ]; then
  echo "[turnstile] --est-rss G is REQUIRED (price from LEDGER.log priors or a measured sibling datum; caps-from-measurement)" >&2; exit 2
fi
if [ "$MODE" = exempt ] && [ "$EST_RSS" -ge "$EXEMPT_MAX_G" ] && [ "${TURNSTILE_FORCE_EXEMPT:-0}" != 1 ]; then
  echo "[turnstile] est-rss ${EST_RSS}G >= ${EXEMPT_MAX_G}G is token-class — drop --exempt (or TURNSTILE_FORCE_EXEMPT=1 with an explicit override)" >&2; exit 2
fi

QDIR="${LOCK%.lock}.q"
mkdir -p "$QDIR/waiters" "$QDIR/exempt" "$QDIR/width"
LEDGER="$QDIR/LEDGER.log"
HOST="$(hostname)"
QTS="$(date +%s)"
say() { echo "[turnstile $(date -u +%H:%M:%SZ)] $*" >&2; }

# ------------------------------------------------------- width token setup ---
pin_ncpu() { taskset -c "$PIN_SPAN" nproc 2>/dev/null || nproc; }
WIDTH_CAP="${TURNSTILE_WIDTH_CAP:-$(pin_ncpu)}"
if [ -z "$EST_THREADS" ]; then
  EST_THREADS=$WIDTH_CAP
  say "WARN: --est-threads undeclared — defaulting to FULL pin width ${EST_THREADS}t (nproc of taskset span $PIN_SPAN; width token). Undeclared jobs SERIALIZE on the pin — declare a measured width to co-fire."
fi
case "$EST_THREADS" in
  ''|*[!0-9]*) echo "[turnstile] --est-threads must be a non-negative integer (0 = sub-1-thread helper, skips the width ledger)" >&2; exit 2 ;;
esac
if [ "$EST_THREADS" -eq 0 ]; then
  WIDTH_SKIP=1   # sub-1-thread helper: the ONLY class that skips ledger + gate (width is physical; exempt legs BLOCK)
elif [ "$EST_THREADS" -gt "$WIDTH_CAP" ]; then
  say "WARN: --est-threads ${EST_THREADS} exceeds pin cap ${WIDTH_CAP}t — clamped (a job cannot exceed its physical span)"
  EST_THREADS=$WIDTH_CAP
fi

# ------------------------------------------------------------------ helpers ---
avail_gib() { awk '/MemAvailable/{printf "%d\n", $2/1048576; exit}' /proc/meminfo; }
pg_rss_gib() { ps -eo pgid=,sid=,rss= 2>/dev/null | awk -v p="$1" '$1==p||$2==p{s+=$3} END{printf "%d\n", s/1048576}'; }  # session-scope fallback: payloads re-group out of the setsid pgid, so pgid-only accounting reads 0 (ramp=full est, peak_rss_g=0)
pg_alive() { ps -eo pgid=,sid= 2>/dev/null | awk -v p="$1" '$1==p||$2==p{f=1} END{exit !f}'; }
pg_rcount() { ps -eLo pgid=,sid=,stat= 2>/dev/null | awk -v p="$1" '($1==p||$2==p) && $3 ~ /^R/ {n++} END{print n+0}'; }  # live R-state threads, session-scoped like pg_rss_gib
pin_rcount() { # live R-state threads whose last CPU sits in the pin span (fire-time census, diagnostic only)
  case "$PIN_SPAN" in
    [0-9]*-[0-9]*) ps -eLo stat=,psr= 2>/dev/null | awk -v lo="${PIN_SPAN%%-*}" -v hi="${PIN_SPAN##*-}" '$1 ~ /^R/ && $2>=lo && $2<=hi {n++} END{print n+0}' ;;
    *) ps -eLo stat= 2>/dev/null | awk '$1 ~ /^R/ {n++} END{print n+0}' ;;
  esac
}
entry_get() { awk -F= -v k="$2" '$1==k{sub(/^[^=]*=/,""); print; exit}' "$1" 2>/dev/null; }
cot_gib() { local v; v=$(cat "$QDIR/COTENANT_RESERVE_G" 2>/dev/null | tr -dc '0-9'); echo "${v:-0}"; }
now_iso() { date -u +%Y-%m-%dT%H:%M:%SZ; }
CMD_STR="$(printf '%s ' "${CMD[@]}" | tr '\n' ' ')"

ledger_line() { # append one key=val line to the shared ledger (caller holds qlock or accepts >>-atomicity)
  echo "ts=$(now_iso) host=$HOST $*" >> "$LEDGER"
}

# effective priority: operator file overrides registered prio. First match wins.
# PRIORITIES.conf lines: "<lane-glob> <prio>"  (e.g. "solveA_* 10", "jobB_* 40", "default 50")
eff_prio() { # $1=lane $2=registered
  local lane="$1" reg="$2" g p
  if [ -r "$QDIR/PRIORITIES.conf" ]; then
    while read -r g p _; do
      case "$g" in ''|'#'*) continue ;; esac
      # shellcheck disable=SC2254
      case "$lane" in $g) echo "$p"; return ;; esac
      [ "$g" = default ] && { echo "$p"; return; }
    done < "$QDIR/PRIORITIES.conf"
  fi
  echo "$reg"
}

# qlock: short critical-section lock over the queue dir (fd 9)
exec 9>>"$QDIR/.qlock"
qlock()   { flock -w 30 9 || { say "FATAL: qlock timeout"; exit 70; }; }
qunlock() { flock -u 9; }

prune_stale() { # caller holds qlock. Stale = wrapper pid dead AND (pgid==0 or pgid empty).
  local f pid pgid
  for f in "$QDIR"/waiters/* "$QDIR"/exempt/* "$QDIR"/width/*; do   # width entries pruned like the RAM I3 pruner
    [ -e "$f" ] || continue
    pid=$(entry_get "$f" pid); pgid=$(entry_get "$f" pgid)
    if ! kill -0 "${pid:-0}" 2>/dev/null; then
      if [ -z "${pgid:-}" ] || [ "${pgid:-0}" = 0 ] || ! pg_alive "$pgid"; then
        ledger_line "event=stale_pruned entry=$(basename "$f") lane=$(entry_get "$f" lane) mode=$(entry_get "$f" mode) est_rss_g=$(entry_get "$f" est_rss_g) est_threads=$(entry_get "$f" est_threads)"
        rm -f "$f"
      fi
    fi
  done
  # stale HOLDER: both pids dead AND the kernel flock is actually free
  if [ -e "$QDIR/HOLDER" ]; then
    pid=$(entry_get "$QDIR/HOLDER" pid); pgid=$(entry_get "$QDIR/HOLDER" pgid)
    if ! kill -0 "${pid:-0}" 2>/dev/null && { [ "${pgid:-0}" = 0 ] || ! pg_alive "$pgid"; }; then
      local pfd
      if exec {pfd}>>"$LOCK" && flock -n "$pfd"; then
        ledger_line "event=stale_holder_reaped lane=$(entry_get "$QDIR/HOLDER" lane)"
        rm -f "$QDIR/HOLDER"; flock -u "$pfd"
      fi
      [ -n "${pfd:-}" ] && exec {pfd}>&-
    fi
  fi
}

# sum of remaining ramps: running exempts (est - current pgid RSS, floored at 0)
# + token holder's remaining ramp (skipped when we ARE the holder).
sum_reserve_gib() { # $1 = "with_holder" | "no_holder"
  local total=0 f est pgid cur r pid
  for f in "$QDIR"/exempt/*; do
    [ -e "$f" ] || continue
    [ "$(entry_get "$f" state)" = running ] || continue
    [ "$f" = "${MY_ENTRY:-}" ] && continue
    pid=$(entry_get "$f" pid); pgid=$(entry_get "$f" pgid)
    kill -0 "${pid:-0}" 2>/dev/null || pg_alive "${pgid:-0}" || continue
    est=$(entry_get "$f" est_rss_g); cur=$(pg_rss_gib "${pgid:-0}")
    r=$(( est - cur )); [ "$r" -lt 0 ] && r=0
    total=$(( total + r ))
  done
  if [ "$1" = with_holder ] && [ -e "$QDIR/HOLDER" ]; then
    pid=$(entry_get "$QDIR/HOLDER" pid); pgid=$(entry_get "$QDIR/HOLDER" pgid)
    if kill -0 "${pid:-0}" 2>/dev/null || pg_alive "${pgid:-0}"; then
      est=$(entry_get "$QDIR/HOLDER" est_rss_g)
      cur=0; [ "${pgid:-0}" != 0 ] && cur=$(pg_rss_gib "$pgid")
      r=$(( ${est:-0} - cur )); [ "$r" -lt 0 ] && r=0
      total=$(( total + r ))
    fi
  fi
  echo "$total"
}

any_ramping_peer() { # stagger clause: another running exempt, younger than RAMP_SECS, below 90% est
  local f pid pgid est cur age ts
  for f in "$QDIR"/exempt/*; do
    [ -e "$f" ] || continue
    [ "$f" = "${MY_ENTRY:-}" ] && continue
    [ "$(entry_get "$f" state)" = running ] || continue
    pid=$(entry_get "$f" pid); pgid=$(entry_get "$f" pgid)
    kill -0 "${pid:-0}" 2>/dev/null || pg_alive "${pgid:-0}" || continue
    ts=$(entry_get "$f" started_ts); [ -n "$ts" ] || ts=$QTS
    age=$(( $(date +%s) - ts ))
    [ "$age" -ge "$RAMP_SECS" ] && continue
    est=$(entry_get "$f" est_rss_g); cur=$(pg_rss_gib "${pgid:-0}")
    if [ "$(( cur * 10 ))" -lt "$(( est * 9 ))" ]; then return 0; fi
  done
  return 1
}

# gate: the exempt-aggregate law, mechanized. $1 = with_holder|no_holder. Sets GATE_WHY on failure.
gate_ok() {
  local avail res cot need
  avail=$(avail_gib); res=$(sum_reserve_gib "$1"); cot=$(cot_gib)
  need=$(( EST_RSS + res + cot + FLOOR_G ))
  if [ "$avail" -lt "$need" ]; then
    GATE_WHY="avail=${avail}G < own_est=${EST_RSS} + ramps=${res} + cotenant=${cot} + floor=${FLOOR_G} (=${need}G)"
    return 1
  fi
  if [ "$MODE" = exempt ] && [ "$NOSTAGGER" = 0 ] && any_ramping_peer; then
    GATE_WHY="stagger: another exempt is still ramping (one ramp at a time)"
    return 1
  fi
  return 0
}

# ------------------------------------------------------------- width token ---
width_sum() { # caller holds qlock: sum est_threads over ALIVE width entries (excluding own)
  local total=0 f pid pgid et
  for f in "$QDIR"/width/*; do
    [ -e "$f" ] || continue
    [ "$f" = "${MY_WIDTH:-}" ] && continue
    pid=$(entry_get "$f" pid); pgid=$(entry_get "$f" pgid)
    kill -0 "${pid:-0}" 2>/dev/null || pg_alive "${pgid:-0}" || continue
    et=$(entry_get "$f" est_threads); et=${et//[^0-9]/}
    total=$(( total + ${et:-0} ))
  done
  echo "$total"
}

width_ok() { # caller holds qlock. Sets GATE_WHY on failure. Width is PHYSICAL: exempt legs block too.
  [ "$WIDTH_SKIP" = 1 ] && return 0
  [ "$NOWIDTHGATE" = 1 ] && return 0
  local ws; ws=$(width_sum)
  if [ $(( ws + EST_THREADS )) -gt "$WIDTH_CAP" ]; then
    GATE_WHY="width: pin ledger ${ws}t + own ${EST_THREADS}t > cap ${WIDTH_CAP}t (width token, pin $PIN_SPAN)"
    return 1
  fi
  return 0
}

write_width_entry() { # atomic rewrite of own width entry
  local tmp="$MY_WIDTH.tmp.$$"
  { echo "lane=$LANE"; echo "mode=$MODE"; echo "est_threads=$EST_THREADS"
    echo "pid=$$"; echo "pgid=${CHILD_PGID:-0}"; echo "ts=${FIRED_TS:-$QTS}"
    echo "peak_r=${PEAK_R:-0}"; echo "host=$HOST"; echo "cmd=$CMD_STR"; } > "$tmp"
  mv -f "$tmp" "$MY_WIDTH"
}

width_reserve() { # caller holds qlock: reserve BEFORE releasing qlock (no fire-fire race)
  [ "$WIDTH_SKIP" = 1 ] && return 0
  local ws; ws=$(width_sum)
  FIRED_TS=$(date +%s)
  MY_WIDTH="$QDIR/width/$QTS.$$"
  write_width_entry
  ledger_line "event=width_reserve lane=$LANE mode=$MODE est_threads=$EST_THREADS ledger_before_t=$ws cap_t=$WIDTH_CAP pin_r_at_fire=$(pin_rcount) gated=$(( 1 - NOWIDTHGATE ))"
}

write_entry() { # $1=path; state/pgid/started from globals
  local tmp="$1.tmp.$$"
  { echo "lane=$LANE"; echo "prio=$PRIO"; echo "mode=$MODE"; echo "est_rss_g=$EST_RSS"
    echo "est_threads=$EST_THREADS"
    echo "est_wall_s=$EST_WALL"; echo "pid=$$"; echo "pgid=${CHILD_PGID:-0}"
    echo "state=${STATE:-waiting}"; echo "queued_ts=$QTS"; echo "started_ts=${STARTED_TS:-}"
    echo "peak_rss_g=${PEAK_G:-0}"; echo "host=$HOST"; echo "cmd=$CMD_STR"; } > "$tmp"
  mv -f "$tmp" "$1"
}

CHILD_PGID=0 STATE=waiting STARTED_TS="" PEAK_G=0 GATE_WHY="" GOT_TOKEN=0 MY_ENTRY="" MON_PID=0
DEADLINE=0; [ "$TIMEOUT" -gt 0 ] && DEADLINE=$(( QTS + TIMEOUT ))

cleanup() {
  local rc=$1
  [ "$MON_PID" != 0 ] && kill "$MON_PID" 2>/dev/null
  if [ "$CHILD_PGID" != 0 ] && pg_alive "$CHILD_PGID"; then
    say "forwarding TERM to pgid $CHILD_PGID"
    kill -TERM -- "-$CHILD_PGID" 2>/dev/null; sleep 2
    pg_alive "$CHILD_PGID" && kill -KILL -- "-$CHILD_PGID" 2>/dev/null
  fi
  qlock 2>/dev/null || true
  [ -n "$MY_ENTRY" ] && rm -f "$MY_ENTRY"
  [ -n "$MY_WIDTH" ] && rm -f "$MY_WIDTH"   # release width at exit; SIGKILL leftovers are pruned by liveness
  if [ "$GOT_TOKEN" = 1 ]; then rm -f "$QDIR/HOLDER"; fi
  local wall=$(( $(date +%s) - QTS ))
  ledger_line "event=done mode=$MODE lane=$LANE prio=$PRIO est_rss_g=$EST_RSS peak_rss_g=$PEAK_G est_threads=$EST_THREADS peak_r=$PEAK_R wall_s=$wall rc=$rc queued_ts=$QTS cmd=\"$CMD_STR\""
  qunlock 2>/dev/null || true
}
on_sig() { say "signal — shutting down"; cleanup 143; exit 143; }
trap on_sig INT TERM

wait_or_fail() { # returns 1 (caller should exit 75) on --no-wait or timeout, else sleeps one poll
  if [ "$NOWAIT" = 1 ]; then say "REFUSED (--no-wait): ${GATE_WHY:-not next in queue / token busy}"; return 1; fi
  if [ "$DEADLINE" != 0 ] && [ "$(date +%s)" -ge "$DEADLINE" ]; then say "TIMEOUT after ${TIMEOUT}s: ${GATE_WHY:-waiting}"; return 1; fi
  sleep "$POLL.$(( RANDOM % 10 ))"
  return 0
}

spawn_and_wait() { # spawn CMD in its own pgid; monitor peak; wait. Sets RC.
  STARTED_TS=$(date +%s); STATE=running
  if [ "$GOT_TOKEN" = 1 ]; then
    # fd 8 (the token flock) is INTENTIONALLY inherited: the lock stays held by
    # the job's fd tree even if this wrapper is SIGKILLed. TOCTOU-safe by construction.
    setsid env TURNSTILE_HELD=1 TURNSTILE_LANE="$LANE" TURNSTILE_EST_THREADS="$EST_THREADS" taskset -c "$PIN_SPAN" "${CMD[@]}" 9>&- &
  else
    setsid env TURNSTILE_EXEMPT=1 TURNSTILE_LANE="$LANE" TURNSTILE_EST_THREADS="$EST_THREADS" taskset -c "$PIN_SPAN" "${CMD[@]}" 9>&- &
  fi
  local child=$!; CHILD_PGID=$child
  qlock; write_entry "$MY_ENTRY"; [ "$GOT_TOKEN" = 1 ] && write_entry "$QDIR/HOLDER"
  [ -n "$MY_WIDTH" ] && write_width_entry   # stamp the real pgid into the width entry
  qunlock
  say "launched pgid $child (mode=$MODE lane=$LANE est=${EST_RSS}G width=${EST_THREADS}t)"
  ( DRIFT_LOGT=0
    while kill -0 "$child" 2>/dev/null; do
      c=$(pg_rss_gib "$child"); p=$(entry_get "$MY_ENTRY" peak_rss_g); p=${p:-0}
      if [ "$c" -gt "$p" ]; then
        awk -F= -v v="$c" 'BEGIN{OFS="="} $1=="peak_rss_g"{$0="peak_rss_g="v} {print}' \
          "$MY_ENTRY" > "$MY_ENTRY.tmp.$$" 2>/dev/null && mv -f "$MY_ENTRY.tmp.$$" "$MY_ENTRY"
        if [ "$GOT_TOKEN" = 1 ] && [ "$MY_ENTRY" != "$QDIR/HOLDER" ]; then cp -f "$MY_ENTRY" "$QDIR/HOLDER" 2>/dev/null; fi
      fi
      if [ -n "${MY_WIDTH:-}" ] && [ -e "$MY_WIDTH" ]; then   # measured-vs-declared width audit
        rl=$(pg_rcount "$child"); pr=$(entry_get "$MY_WIDTH" peak_r); pr=${pr:-0}
        if [ "$rl" -gt "$pr" ]; then
          awk -F= -v v="$rl" 'BEGIN{OFS="="} $1=="peak_r"{$0="peak_r="v} {print}' \
            "$MY_WIDTH" > "$MY_WIDTH.tmp.$$" 2>/dev/null && mv -f "$MY_WIDTH.tmp.$$" "$MY_WIDTH"
        fi
        if [ "$WIDTH_SKIP" = 0 ] && [ "$rl" -gt $(( EST_THREADS * 2 )) ] && \
           [ $(( $(date +%s) - DRIFT_LOGT )) -ge "${TURNSTILE_DRIFT_LOG_SECS:-600}" ]; then
          ledger_line "event=width_drift lane=$LANE pgid=$child est_threads=$EST_THREADS live_r=$rl note=live>2x-declared--retype-width-cap"
          DRIFT_LOGT=$(date +%s)
        fi
      fi
      sleep "$MON"
    done ) 8>&- 9>&- & MON_PID=$!   # monitor must NOT inherit the token fd — its stragglers would hold the lock past exit
  wait "$child"; RC=$?
  PEAK_G=$(entry_get "$MY_ENTRY" peak_rss_g); PEAK_G=${PEAK_G:-0}
  [ -n "$MY_WIDTH" ] && { PEAK_R=$(entry_get "$MY_WIDTH" peak_r); PEAK_R=${PEAK_R:-0}; }
  kill "$MON_PID" 2>/dev/null; wait "$MON_PID" 2>/dev/null; MON_PID=0
  CHILD_PGID=0
}

# ================================================================ TOKEN path ==
if [ "$MODE" = token ]; then
  MY_ENTRY="$QDIR/waiters/$QTS.$$"
  qlock; write_entry "$MY_ENTRY"; qunlock
  say "queued for token (lane=$LANE prio=$PRIO est=${EST_RSS}G lock=$LOCK)"
  exec 8>>"$LOCK"
  LOGT=0
  while :; do
    if [ -e "$QDIR/HOLD" ]; then
      GATE_WHY="operator HOLD present ($(head -c 120 "$QDIR/HOLD" 2>/dev/null | tr '\n' ' '))"
      wait_or_fail || { cleanup 75; exit 75; }; continue
    fi
    qlock; prune_stale
    # pick best live waiter: min effective prio, then earliest queued_ts, then pid
    best=""; bp=999999; bt=9999999999; bpid=9999999
    for f in "$QDIR"/waiters/*; do
      [ -e "$f" ] || continue
      l=$(entry_get "$f" lane); rp=$(entry_get "$f" prio); ep=$(eff_prio "$l" "${rp:-50}")
      ep=${ep//[^0-9]/}; ep=${ep:-50}
      t=$(entry_get "$f" queued_ts); t=${t:-9999999999}; pd=$(entry_get "$f" pid); pd=${pd:-9999999}
      if [ "$ep" -lt "$bp" ] || { [ "$ep" -eq "$bp" ] && [ "$t" -lt "$bt" ]; } || \
         { [ "$ep" -eq "$bp" ] && [ "$t" -eq "$bt" ] && [ "$pd" -lt "$bpid" ]; }; then
        best="$f"; bp=$ep; bt=$t; bpid=$pd
      fi
    done
    if [ "$best" = "$MY_ENTRY" ]; then
      if flock -n 8; then
        GOT_TOKEN=1
        [ -e "$QDIR/HOLDER" ] && say "note: reaping leftover HOLDER record (kernel lock was free)"
        write_entry "$QDIR/HOLDER"
        rm -f "$MY_ENTRY"; MY_ENTRY="$QDIR/HOLDER"   # holder record replaces the waiter entry
        qunlock
        say "TOKEN ACQUIRED (prio=$bp; held across entire run)"
        break
      fi
      GATE_WHY="token busy (held by current holder or a legacy raw-flock job)"
    else
      GATE_WHY="not highest-priority waiter (best=$(entry_get "${best:-/dev/null}" lane) prio=$bp)"
    fi
    qunlock
    [ $(( $(date +%s) - LOGT )) -ge 300 ] && { say "waiting: $GATE_WHY"; LOGT=$(date +%s); }
    wait_or_fail || { cleanup 75; exit 75; }
  done
  # launch gate while HOLDING the token (never releases; exempts drain to us)
  # the WIDTH gate joins the RAM gate here. --no-gate bypasses only the
  # RAM gate — width is physical (only --est-threads 0 helpers skip; --no-width-gate
  # is the operator escape and still registers).
  WIDTH_ACTIVE=1; { [ "$WIDTH_SKIP" = 1 ] || [ "$NOWIDTHGATE" = 1 ]; } && WIDTH_ACTIVE=0
  if [ "$NOGATE" = 0 ] || [ "$WIDTH_ACTIVE" = 1 ]; then
    LOGT=0
    while :; do
      if [ "$NOGATE" = 0 ] && [ -e "$QDIR/HOLD" ]; then
        GATE_WHY="operator HOLD present"
      else
        qlock; prune_stale
        PASS=1
        if [ "$NOGATE" = 0 ] && ! gate_ok no_holder; then PASS=0; fi
        if [ "$PASS" = 1 ] && ! width_ok; then PASS=0; fi
        if [ "$PASS" = 1 ]; then width_reserve; qunlock; break; fi   # reserve width BEFORE releasing qlock
        qunlock
      fi
      [ $(( $(date +%s) - LOGT )) -ge 300 ] && { say "token held, gating launch: $GATE_WHY"; LOGT=$(date +%s); }
      if [ "$NOWAIT" = 1 ] || { [ "$DEADLINE" != 0 ] && [ "$(date +%s)" -ge "$DEADLINE" ]; }; then
        say "giving up at launch gate: $GATE_WHY"; cleanup 75; exit 75
      fi
      sleep "$POLL.$(( RANDOM % 10 ))"
    done
  else
    qlock; width_reserve; qunlock   # --no-gate + width bypassed: still REGISTER width (no-op for --est-threads 0)
  fi
  spawn_and_wait
  say "done rc=$RC peak=${PEAK_G}G"
  cleanup "$RC"; exit "$RC"
fi

# =============================================================== EXEMPT path ==
MY_ENTRY="$QDIR/exempt/$QTS.$$"
qlock; write_entry "$MY_ENTRY"; qunlock
say "exempt registered (lane=$LANE est=${EST_RSS}G ledger=$QDIR/exempt)"
LOGT=0
while :; do
  if [ -e "$QDIR/HOLD" ]; then
    GATE_WHY="operator HOLD present ($(head -c 120 "$QDIR/HOLD" 2>/dev/null | tr '\n' ' '))"
    wait_or_fail || { cleanup 75; exit 75; }; continue
  fi
  qlock; prune_stale
  PASS=1
  if [ "$NOGATE" = 0 ] && ! gate_ok with_holder; then PASS=0; fi
  if [ "$PASS" = 1 ] && ! width_ok; then PASS=0; fi   # width is physical — exempt legs BLOCK too
  if [ "$PASS" = 1 ]; then
    STATE=running; STARTED_TS=$(date +%s); write_entry "$MY_ENTRY"   # reserve BEFORE releasing qlock
    width_reserve                                                    # width reserved under the same qlock
    qunlock
    break
  fi
  qunlock
  [ $(( $(date +%s) - LOGT )) -ge 300 ] && { say "exempt gate: $GATE_WHY"; LOGT=$(date +%s); }
  wait_or_fail || { cleanup 75; exit 75; }
done
say "exempt gate PASSED (avail=$(avail_gib)G)"
spawn_and_wait
say "done rc=$RC peak=${PEAK_G}G (ledgered)"
cleanup "$RC"; exit "$RC"

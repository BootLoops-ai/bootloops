#!/usr/bin/env bash
# selftest.sh — Turnstile (ops/turnstile) self-test. Runs entirely against a SCRATCH lock in a
# temp dir; NEVER touches a real production lock file or any live job.
# Tests: (1) token held across entire run (TOCTOU), (2) 3-waiter PRIORITY order,
# (3) operator PRIORITIES.conf override, (4) stale dead-pid entry pruned,
# (5) exempt launch-gate refusal + ledger peak line, (6) sum-of-ramps gate with a
# live fake peer, (7) operator HOLD refusal, (8) WIDTH token: 30/30/30
# on cap 64 (third queues, fires on free), --est-threads 0 helper skip,
# undeclared-width WARN+default, (9) daemon reap: the battery kills every
# process it spawned before exiting — no stray may outlive the battery holding
# an inherited stdout/stderr pipe (that stalls any harness reading the pipe).
# Linux-only (like the wrapper under test): on darwin every leg SKIPs by name.
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
RUN="$HERE/turnstile_run.sh"

# ---- platform gate: named per-leg skip on darwin ---------------------------
# Every leg drives turnstile_run.sh, which needs Linux-only tooling; leg (1) also
# calls flock(1) directly. Skip loudly BY NAME (exit 0) instead of failing on
# absent tools.
if [ "$(uname -s)" = Darwin ]; then
  echo "== turnstile selftest =="
  echo "SKIP (darwin) legs (1)-(7): token/ledger legs drive turnstile_run.sh — Linux-only: flock(1) kernel lock token, /proc/meminfo MemAvailable gate, GNU ps session columns (ps -eo pgid=,sid=)"
  echo "SKIP (darwin) legs (8)-(8g): width-token legs — Linux-only: taskset pin spans, nproc pin-width cap, GNU ps thread listing (ps -eLo)"
  echo "SKIP (darwin) leg (9): daemon-reap leg — Linux-only: /proc/<pid>/environ scan"
  echo "== selftest: SKIPPED on darwin (Linux-only: taskset, nproc, flock(1), /proc); 0 fail =="
  exit 0
fi

T="${TURNSTILE_TEST_DIR:-${BIGRAM_TEST_DIR:-$(mktemp -d "${TMPDIR:-/tmp}/turnstile_selftest.XXXXXX")}}"
mkdir -p "$T"
# scrub inherited knobs (new and legacy names) so the battery is hermetic
for _k in LOCK POLL MON FLOOR_G RAMP_SECS EXEMPT_MAX_G PIN WIDTH_CAP DRIFT_LOG_SECS FORCE_EXEMPT; do unset "TURNSTILE_$_k" "BIGRAM_$_k"; done; unset _k
export TURNSTILE_LOCK="$T/TEST_TURNSTILE.lock"
export TURNSTILE_POLL=1 TURNSTILE_MON=1 TURNSTILE_FLOOR_G=0 TURNSTILE_RAMP_SECS=2
Q="${TURNSTILE_LOCK%.lock}.q"
PASS=0; FAIL=0
ok()  { PASS=$((PASS+1)); echo "PASS  $*"; }
bad() { FAIL=$((FAIL+1)); echo "FAIL  $*"; }
cleanup_q() { rm -rf "$Q"; rm -f "$TURNSTILE_LOCK"; }

# ---- daemon reap: this battery kills its own strays before exit ------------
# Wrapper payloads are setsid'd into their own sessions, so shell job control
# cannot see them; a survivor keeps the battery's inherited stdout pipe open
# and stalls any harness reading it to EOF. Match is PRECISE — this run's
# unique scratch lock path in the process environment — never a name-based
# pkill. (Linux /proc scan is fine here: the darwin gate above exits first.)
list_survivors() { # every OTHER live pid still carrying this run's lock path in its env
  local f p
  [ -d "$T" ] || return 0
  grep -laFs "TURNSTILE_LOCK=$TURNSTILE_LOCK" /proc/[0-9]*/environ 2>/dev/null > "$T/.psscan" || true
  while IFS= read -r f; do
    p=${f#/proc/}; p=${p%/environ}
    case "$p" in "$$"|"$PPID"|"$BASHPID") continue ;; esac
    kill -0 "$p" 2>/dev/null && echo "$p"
  done < "$T/.psscan"
}
reap_daemons() { # TERM every stray, join with a bounded wait (<=5 s), then KILL
  local p i strays
  strays="$(list_survivors)"
  [ -n "$strays" ] || return 0
  for p in $strays; do kill -TERM "$p" 2>/dev/null; done
  for i in $(seq 1 25); do
    strays="$(list_survivors)"; [ -z "$strays" ] && return 0
    sleep 0.2
  done
  for p in $strays; do kill -KILL "$p" 2>/dev/null; done
  for i in 1 2 3 4 5; do
    [ -z "$(list_survivors)" ] && break
    sleep 0.2
  done
}
trap reap_daemons EXIT     # every exit path reaps; signals route through exit
trap 'exit 130' INT
trap 'exit 143' TERM

echo "== turnstile selftest (scratch: $T) =="

# ---- (1) TOCTOU: token held across the ENTIRE run --------------------------
cleanup_q
"$RUN" --lane t1 --prio 50 --est-rss 1 -- sleep 3 2>"$T/t1.log" & W=$!
sleep 1.4
if flock -n "$TURNSTILE_LOCK" true 2>/dev/null; then bad "(1) lock free mid-run — NOT held across run"; else ok "(1) token held across entire run (flock -n refused mid-run)"; fi
wait "$W"; rc=$?
[ "$rc" = 0 ] && ok "(1b) wrapped cmd rc=0 propagated" || bad "(1b) rc=$rc"
flock -n "$TURNSTILE_LOCK" true 2>/dev/null && ok "(1c) lock released after run" || bad "(1c) lock still held after exit"
grep -q "event=done mode=token lane=t1" "$Q/LEDGER.log" && ok "(1d) token ledger line written" || bad "(1d) no ledger line"

# ---- (2) 3-waiter PRIORITY order (arrival order scrambled on purpose) ------
cleanup_q
ORD="$T/order.txt"; : > "$ORD"
"$RUN" --lane holder --prio 50 --est-rss 1 -- sleep 6 2>"$T/h2.log" & H=$!
sleep 1
"$RUN" --lane econ  --class low            --est-rss 1 -- sh -c "echo econ  >> $ORD" 2>"$T/w_econ.log" & W1=$!
sleep 0.4
"$RUN" --lane table --class critical       --est-rss 1 -- sh -c "echo table >> $ORD" 2>"$T/w_table.log" & W2=$!
sleep 0.4
"$RUN" --lane oracle --class high          --est-rss 1 -- sh -c "echo oracle >> $ORD" 2>"$T/w_oracle.log" & W3=$!
wait "$H" "$W1" "$W2" "$W3"
got=$(tr '\n' ' ' < "$ORD" | sed 's/ *$//')
if [ "$got" = "table oracle econ" ]; then ok "(2) 3-waiter priority order: '$got' (arrival was econ,table,oracle)"; else bad "(2) order '$got' != 'table oracle econ'"; fi

# ---- (3) operator PRIORITIES.conf override --------------------------------
cleanup_q
mkdir -p "$Q"
printf 'lateA 5\n' > "$Q/PRIORITIES.conf"     # operator promotes lateA over everything
: > "$ORD"
"$RUN" --lane holder --prio 50 --est-rss 1 -- sleep 5 2>"$T/h3.log" & H=$!
sleep 1
"$RUN" --lane earlyB --class high           --est-rss 1 -- sh -c "echo earlyB >> $ORD" 2>"$T/w_eb.log" & W1=$!
sleep 0.4
"$RUN" --lane lateA  --class misc           --est-rss 1 -- sh -c "echo lateA >> $ORD" 2>"$T/w_la.log" & W2=$!
wait "$H" "$W1" "$W2"
got=$(tr '\n' ' ' < "$ORD" | sed 's/ *$//')
if [ "$got" = "lateA earlyB" ]; then ok "(3) PRIORITIES.conf override: misc job promoted to prio 5 ran first"; else bad "(3) order '$got' != 'lateA earlyB'"; fi

# ---- (4) stale dead-pid waiter entry: pruned, does not block ---------------
cleanup_q
mkdir -p "$Q/waiters"
sleep 300 & ZPID=$!; kill -9 "$ZPID" 2>/dev/null; wait "$ZPID" 2>/dev/null
cat > "$Q/waiters/1000000000.$ZPID" <<EOF
lane=ghost
prio=1
mode=token
est_rss_g=1
est_wall_s=0
pid=$ZPID
pgid=0
state=waiting
queued_ts=1000000000
started_ts=
peak_rss_g=0
host=$(hostname)
cmd=ghost
EOF
: > "$ORD"
"$RUN" --lane real --prio 50 --est-rss 1 --timeout 86400 -- sh -c "echo real >> $ORD" 2>"$T/w_real.log"; rc=$?
if [ "$rc" = 0 ] && grep -q real "$ORD"; then ok "(4) dead-pid prio-1 ghost pruned; live prio-50 waiter ran (rc=0)"; else bad "(4) rc=$rc, ghost blocked the queue"; fi
[ -e "$Q/waiters/1000000000.$ZPID" ] && bad "(4b) stale entry file still present" || ok "(4b) stale entry file removed"
grep -q "event=stale_pruned.*lane=ghost" "$Q/LEDGER.log" && ok "(4c) stale prune ledgered" || bad "(4c) no stale_pruned ledger line"

# ---- (5) exempt gate refusal (floor) + ledger peak on success --------------
cleanup_q
TURNSTILE_FLOOR_G=99999999 "$RUN" --exempt --lane e1 --est-rss 1 --no-wait -- true 2>"$T/e1.log"; rc=$?
[ "$rc" = 75 ] && ok "(5) exempt gate REFUSES when floor law violated (rc=75)" || bad "(5) rc=$rc expected 75"
"$RUN" --exempt --lane e2 --est-rss 1 -- sleep 2 2>"$T/e2.log"; rc=$?
if [ "$rc" = 0 ] && grep -q "event=done mode=exempt lane=e2.*peak_rss_g=" "$Q/LEDGER.log"; then ok "(5b) exempt ran and registered peak in shared ledger"; else bad "(5b) rc=$rc or no exempt ledger line"; fi

# ---- (6) sum-of-ramps gate: live fake peer with huge est blocks launch -----
cleanup_q
mkdir -p "$Q/exempt"
sleep 300 & PPID2=$!
cat > "$Q/exempt/1000000001.$PPID2" <<EOF
lane=fatpeer
prio=50
mode=exempt
est_rss_g=99999999
est_wall_s=0
pid=$PPID2
pgid=$PPID2
state=running
queued_ts=1000000001
started_ts=1000000001
peak_rss_g=0
host=$(hostname)
cmd=fake
EOF
"$RUN" --exempt --lane e3 --est-rss 1 --no-wait --no-stagger -- true 2>"$T/e3.log"; rc=$?
[ "$rc" = 75 ] && ok "(6) SUM(ramping exempts)+floor > avail => REFUSED (rc=75, exempt-aggregate law)" || bad "(6) rc=$rc expected 75"
kill -9 "$PPID2" 2>/dev/null; wait "$PPID2" 2>/dev/null
"$RUN" --exempt --lane e4 --est-rss 1 --timeout 86400 --no-stagger -- true 2>"$T/e4.log"; rc=$?
[ "$rc" = 0 ] && ok "(6b) dead peer pruned by liveness; gate passes (rc=0)" || bad "(6b) rc=$rc expected 0"

# ---- (7) operator HOLD blocks all new launches ----------------------------
cleanup_q
mkdir -p "$Q"; echo "eta solve swell — operator hold" > "$Q/HOLD"
"$RUN" --exempt --lane e5 --est-rss 1 --no-wait -- true 2>"$T/e5.log"; rc1=$?
"$RUN" --lane t5 --prio 10 --est-rss 1 --no-wait -- true 2>"$T/t5.log"; rc2=$?
if [ "$rc1" = 75 ] && [ "$rc2" = 75 ]; then ok "(7) HOLD refuses both exempt and token launches (rc=75,75)"; else bad "(7) rc=$rc1/$rc2 expected 75/75"; fi
rm -f "$Q/HOLD"

# ---- (8) WIDTH token: 30/30/30 on cap 64 — third queues, fires on free
# (cap pinned explicitly so the leg is machine-independent; the wrapper's default
# cap is the nproc of the taskset span)
export TURNSTILE_WIDTH_CAP=64
cleanup_q
: > "$ORD"
"$RUN" --exempt --lane w1 --est-rss 1 --est-threads 30 --no-stagger -- sleep 4 2>"$T/w1.log" & X1=$!
sleep 0.6
"$RUN" --exempt --lane w2 --est-rss 1 --est-threads 30 --no-stagger -- sleep 9 2>"$T/w2.log" & X2=$!
sleep 0.6
"$RUN" --exempt --lane w3 --est-rss 1 --est-threads 30 --no-stagger -- sh -c "echo w3 >> $ORD" 2>"$T/w3.log" & X3=$!
sleep 1.5
grep -q "launched" "$T/w3.log" && bad "(8) w3 fired while ledger at 60t (should queue)" || ok "(8) w3 queued at 60t+30t>64t"
grep -q "width: pin ledger" "$T/w3.log" && ok "(8b) width GATE_WHY logged" || bad "(8b) no width gate message in w3 log"
wait "$X1"
wait "$X3"; rc=$?
if [ "$rc" = 0 ] && grep -q w3 "$ORD"; then ok "(8c) w3 fired after w1 freed 30t (rc=0)"; else bad "(8c) rc=$rc"; fi
"$RUN" --exempt --lane w4 --est-rss 1 --est-threads 0 --no-wait --no-stagger -- true 2>"$T/w4.log"; rc=$?
[ "$rc" = 0 ] && ok "(8d) --est-threads 0 helper skips width ledger+gate (rc=0 while w2 holds 30t)" || bad "(8d) rc=$rc"
wait "$X2"
grep -q "event=width_reserve lane=w1" "$Q/LEDGER.log" && ok "(8e) width_reserve ledgered" || bad "(8e) no width_reserve ledger line"
"$RUN" --exempt --lane w5 --est-rss 1 --no-wait --no-stagger -- true 2>"$T/w5.log"; rc=$?
grep -q "WARN: --est-threads undeclared" "$T/w5.log" && ok "(8f) undeclared width WARNs + defaults to pin nproc" || bad "(8f) no undeclared-width WARN"
[ "$rc" = 0 ] && ok "(8g) undeclared job fires alone on empty ledger (rc=0)" || bad "(8g) rc=$rc"

# ---- (9) daemon reap: no spawned process outlives the battery --------------
# Plant the exact hazard class: a payload that daemonizes (double-setsid) out of
# its wrapper's pgid, inheriting the battery's stdout — the survivor class that
# holds a harness's pipe open past battery exit. Then prove the reap clears it.
cleanup_q
"$RUN" --exempt --lane d1 --est-rss 1 --est-threads 0 --no-wait --no-stagger -- sh -c 'setsid sleep 300 &' 2>"$T/d1.log"
sleep 0.5
[ -n "$(list_survivors)" ] && ok "(9) planted daemonized payload survives its wrapper (hazard reproduced)" || bad "(9) plant failed — nothing to reap (reap path untested)"
reap_daemons
left="$(list_survivors)"
[ -z "$left" ] && ok "(9b) reap_daemons cleared all strays (TERM -> bounded join -> KILL); no pipe holder outlives the battery" || bad "(9b) survivors after reap: $left"

trap - EXIT
reap_daemons   # final net for legs (1)-(8) stragglers (e.g. monitor sleeps), before scratch scrub
echo "== selftest: $PASS pass, $FAIL fail =="
[ -z "${TURNSTILE_TEST_DIR:-}${BIGRAM_TEST_DIR:-}" ] && rm -rf "$T"
[ "$FAIL" = 0 ]

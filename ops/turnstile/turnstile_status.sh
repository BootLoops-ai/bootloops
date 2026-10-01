#!/usr/bin/env bash
# turnstile_status.sh — read-only view of the Turnstile priority queue + exempt RAM ledger + width ledger.
# Usage: turnstile_status.sh [--lock PATH]     (default: $TURNSTILE_LOCK or /var/tmp/TURNSTILE.lock)
# (legacy BIGRAM_* env names are honored as fallbacks, as in turnstile_run.sh)
set -uo pipefail
for _k in LOCK FLOOR_G PIN WIDTH_CAP; do
  _t="TURNSTILE_$_k"; _b="BIGRAM_$_k"
  if [ -z "${!_t:-}" ] && [ -n "${!_b:-}" ]; then printf -v "$_t" '%s' "${!_b}"; fi
done
unset _k _t _b
# Platform gate: same Linux-only tooling as turnstile_run.sh (nproc, taskset,
# /proc/meminfo, flock(1), ps -eLo). Refuse by name on darwin — a status view
# built on absent tools would print garbage, not status.
if [ "$(uname -s)" = Darwin ]; then
  echo "[turnstile] Linux-only: taskset unavailable on darwin (also required: nproc, /proc/meminfo, flock(1), ps -eLo thread listing). No macOS port — use a Linux host." >&2
  exit 69   # EX_UNAVAILABLE
fi
DEF_LOCK=/var/tmp/TURNSTILE.lock
LOCK="${TURNSTILE_LOCK:-$DEF_LOCK}"
[ "${1:-}" = --lock ] && LOCK="$2"
QDIR="${LOCK%.lock}.q"
FLOOR_G="${TURNSTILE_FLOOR_G:-150}"

eg() { awk -F= -v k="$2" '$1==k{sub(/^[^=]*=/,""); print; exit}' "$1" 2>/dev/null; }
pg_rss_gib() { ps -eo pgid=,rss= 2>/dev/null | awk -v p="$1" '$1==p{s+=$2} END{printf "%d\n", s/1048576}'; }
alive() { kill -0 "${1:-0}" 2>/dev/null && echo LIVE || echo DEAD; }
eff_prio() { local lane="$1" reg="$2" g p
  if [ -r "$QDIR/PRIORITIES.conf" ]; then
    while read -r g p _; do
      case "$g" in ''|'#'*) continue ;; esac
      # shellcheck disable=SC2254
      case "$lane" in $g) echo "$p"; return ;; esac
      [ "$g" = default ] && { echo "$p"; return; }
    done < "$QDIR/PRIORITIES.conf"
  fi
  echo "$reg"; }

avail=$(awk '/MemAvailable/{printf "%d", $2/1048576; exit}' /proc/meminfo)
cot=$(cat "$QDIR/COTENANT_RESERVE_G" 2>/dev/null | tr -dc '0-9'); cot=${cot:-0}
echo "== TURNSTILE queue @ $LOCK  ($(date -u +%Y-%m-%dT%H:%M:%SZ)) =="
echo "avail=${avail}G  floor=${FLOOR_G}G  cotenant_reserve=${cot}G  HOLD=$([ -e "$QDIR/HOLD" ] && echo "YES ($(head -c 80 "$QDIR/HOLD" | tr '\n' ' '))" || echo no)"
[ -r "$QDIR/PRIORITIES.conf" ] && { echo "-- PRIORITIES.conf --"; sed 's/^/   /' "$QDIR/PRIORITIES.conf"; }

echo "-- TOKEN HOLDER --"
if [ -e "$QDIR/HOLDER" ]; then
  h="$QDIR/HOLDER"; pg=$(eg "$h" pgid)
  echo "   lane=$(eg "$h" lane) prio=$(eg "$h" prio) est=$(eg "$h" est_rss_g)G cur=$(pg_rss_gib "${pg:-0}")G peak=$(eg "$h" peak_rss_g)G wrapper=$(eg "$h" pid)($(alive "$(eg "$h" pid)")) pgid=$pg since=$(eg "$h" queued_ts)"
else
  flock -n "$LOCK" true 2>/dev/null && echo "   (free)" || echo "   held by a LEGACY raw-flock job (no HOLDER record)"
fi

echo "-- WAITERS (priority order; eff prio after PRIORITIES.conf) --"
tmp=$(mktemp)
for f in "$QDIR"/waiters/*; do
  [ -e "$f" ] || continue
  lane=$(eg "$f" lane); rp=$(eg "$f" prio); ep=$(eff_prio "$lane" "${rp:-50}"); ep=${ep//[^0-9]/}; ep=${ep:-50}
  echo "$ep $(eg "$f" queued_ts) $lane est=$(eg "$f" est_rss_g)G wrapper=$(eg "$f" pid)($(alive "$(eg "$f" pid)")) reg_prio=$rp" >> "$tmp"
done
[ -s "$tmp" ] && sort -n -k1,1 -k2,2 "$tmp" | awk '{printf "   #%d prio=%s queued=%s %s\n", NR, $1, $2, substr($0, index($0,$3))}' || echo "   (none)"
rm -f "$tmp"

echo "-- EXEMPT LEDGER (running sub-300G) --"
sum=0
for f in "$QDIR"/exempt/*; do
  [ -e "$f" ] || continue
  pg=$(eg "$f" pgid); est=$(eg "$f" est_rss_g); cur=$(pg_rss_gib "${pg:-0}")
  r=$(( ${est:-0} - cur )); [ "$r" -lt 0 ] && r=0
  st=$(eg "$f" state); [ "$st" = running ] && sum=$(( sum + r ))
  echo "   $(eg "$f" lane) state=$st est=${est}G cur=${cur}G peak=$(eg "$f" peak_rss_g)G reserve=${r}G wrapper=$(eg "$f" pid)($(alive "$(eg "$f" pid)"))"
done
echo "   SUM remaining ramps = ${sum}G ; exempt-aggregate headroom for a new EST-G exempt = avail - sum - cotenant - floor = $(( avail - sum - cot - FLOOR_G ))G"

PIN_SPAN="${TURNSTILE_PIN:-0-$(( $(nproc) - 1 ))}"
wcap="${TURNSTILE_WIDTH_CAP:-$(taskset -c "$PIN_SPAN" nproc 2>/dev/null || nproc)}"
pg_alive_st() { ps -eo pgid=,sid= 2>/dev/null | awk -v p="${1:-0}" '$1==p||$2==p{f=1} END{exit !f}'; }
echo "-- WIDTH LEDGER (pin $PIN_SPAN cap ${wcap}t; ALIVE est_threads gate new fires) --"
wsum=0
for f in "$QDIR"/width/*; do
  [ -e "$f" ] || continue
  pid=$(eg "$f" pid); pg=$(eg "$f" pgid); et=$(eg "$f" est_threads); et=${et//[^0-9]/}; et=${et:-0}
  st=DEAD
  if kill -0 "${pid:-0}" 2>/dev/null || pg_alive_st "${pg:-0}"; then st=LIVE; wsum=$(( wsum + et )); fi
  echo "   $(eg "$f" lane) mode=$(eg "$f" mode) est_threads=${et}t peak_r=$(eg "$f" peak_r) wrapper=${pid}($st) pgid=$pg since=$(eg "$f" ts)"
done
liver=$(case "$PIN_SPAN" in [0-9]*-[0-9]*) ps -eLo stat=,psr= 2>/dev/null | awk -v lo="${PIN_SPAN%%-*}" -v hi="${PIN_SPAN##*-}" '$1 ~ /^R/ && $2>=lo && $2<=hi {n++} END{print n+0}' ;; *) echo "?" ;; esac)
echo "   SUM alive est_threads = ${wsum}t / cap ${wcap}t ; headroom for a new fire = $(( wcap - wsum ))t ; live pin R-count (all lanes, incl. unmanaged) = ${liver}"
echo "-- last 5 ledger events --"
tail -5 "$QDIR/LEDGER.log" 2>/dev/null | sed 's/^/   /' || true

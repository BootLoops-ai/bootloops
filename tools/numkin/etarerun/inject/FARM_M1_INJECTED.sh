#!/bin/bash
[ "$(uname -s)" = Linux ] || { echo "FATAL: $(basename "$0") requires Linux — its liveness checks read /proc" >&2; exit 2; }
# FARM_M1_INJECTED.sh  [N_CONC] [PAR_PER_SLICE] [MEM_GB_PER_SLICE] [MODE]
#
# Concurrency-controlled M1 injected-slice farm over the 49-prime eps grid.
# For each eps=1/p:
#   1. wait for source slices (M89_eps_1_p, [M6_eps_1_p]) to have kira_target.m
#   2. build rebase table (rebase_invert.py, once per source per eps)
#   3. extract_shared_rows.py -> uds_<p>
#   4. build_M1_injected.sh -> M1inj_<p>  (MODE=reuse|fresh)
#   5. on completion: merge_M1_injected.sh
#
# When ALL 49 done: reconstruct the 2-var A(d,η) from the slices via
# Thiele/Newton (tools/ratfit/thiele_gate.py is the gate for such a
# reconstruction) and inject it into the ibp-cache with ../inject_ibpcache.sh.
#
# WORKED EXAMPLE: the ETANUMD_DRAINQ layout below (work_M1_sec1022,
# work_M89_sec1018, work_M6_sec1014, work_M5_sec510 — one amflow work dir per
# master sector) is that of one reduction; substitute your own sector dirs.
#
# Usage: setsid nohup bash FARM_M1_INJECTED.sh 3 24 60 reuse \
#          > farm_M1_injected.log 2>&1 </dev/null &
#
set -u
# Template: every input root comes from the ETANUMD_* env vars; no defaults shipped.
HERE=${ETANUMD_INJECT_DIR:?set ETANUMD_INJECT_DIR}
BASE=${ETANUMD_DIR:?set ETANUMD_DIR}
DRAINQ=${ETANUMD_DRAINQ:?set ETANUMD_DRAINQ}
FAM=${ETANUMD_FAM:?set ETANUMD_FAM (kira family name)}
M1SRC=$DRAINQ/work_M1_sec1022/part_0/amf/system_0_diffeq/target_reduce
M1MAST=$M1SRC/results/$FAM/masters
M1TGT=$M1SRC/target
M89MAST=$DRAINQ/work_M89_sec1018/part_0/amf/system_0_diffeq/target_reduce/results/$FAM/masters
M6MAST=$DRAINQ/work_M6_sec1014/part_0/amf/system_0_diffeq/target_reduce/results/$FAM/masters
M5MAST=$DRAINQ/work_M5_sec510/part_0/amf/system_0_diffeq/target_reduce/results/$FAM/masters
M5KTM=$DRAINQ/work_M5_sec510/part_0/amf/system_0_diffeq/target_reduce/results/$FAM/kira_target.m

N_CONC=${1:-4}
PAR=${2:-16}
MEM_GB=${3:-40}
MODE=${4:-freshtgt}   # freshtgt = the faster route (measured: back-substitution block 89.8 s -> 23.9 s)

PRIMES=(7 11 13 17 19 23 29 31 37 41 43 47 53 59 61 67 71 73 79 83 89 97 101 103 107 109 113 127 131 137 139 149 151 157 163 167 173 179 181 191 193 197 199 211 223 227 229 233 239)

running() { ls -l /proc/[0-9]*/cwd 2>/dev/null | grep -c 'etanumd/inject/M1inj_1_'; }

for p in "${PRIMES[@]}"; do
  d="$HERE/M1inj_1_$p"
  # already merged?
  [ -s "$d/results/$FAM/kira_target.m" ] && [ -f "$d/results/$FAM/kira_target.residual.m" ] && \
    { echo "[farm] skip $p (merged)"; continue; }
  # in-progress?
  if [ -d "$d" ] && ls -l /proc/[0-9]*/cwd 2>/dev/null | grep -q "M1inj_1_$p\$"; then
    echo "[farm] skip $p (running)"; continue
  fi

  # ---- gather sources for this eps ----
  SRCARGS=(); RBARGS=()
  # M89 numd slice (primary — always present via farm_M89_cc)
  M89_K="$BASE/M89_eps_1_$p/results/$FAM/kira_target.m"
  if [ -s "$M89_K" ]; then
    RB="$HERE/rebase_M89_1_$p.m"
    [ -s "$RB" ] || python3 "$HERE/rebase_invert.py" --m1-masters "$M1MAST" \
        --src-masters "$M89MAST" --src-ktm "$M89_K" --out "$RB" >/dev/null 2>&1
    SRCARGS+=(--source "M89:$M89_K" --diff-prop M89:3)
    RBARGS+=(--rebase "$RB")
  fi
  # M6 numd slice (if landed)
  M6_K="$BASE/M6_eps_1_$p/results/$FAM/kira_target.m"
  if [ -s "$M6_K" ]; then
    RB6="$HERE/rebase_M6_1_$p.m"
    [ -s "$RB6" ] || python3 "$HERE/rebase_invert.py" --m1-masters "$M1MAST" \
        --src-masters "$M6MAST" --src-ktm "$M6_K" --out "$RB6" >/dev/null 2>&1
    SRCARGS+=(--source "M6:$M6_K" --diff-prop M6:4)
    # merge rebase tables
    if [ ${#RBARGS[@]} -gt 0 ]; then
      cat "$RB6" >> "${RBARGS[1]}" 2>/dev/null || true
    else
      RBARGS+=(--rebase "$RB6")
    fi
  fi
  # M5 symbolic (d,eta) — d-substitute for this slice (if landed)
  if [ -s "$M5KTM" ]; then
    read NUM DEN < <(python3 -c "from fractions import Fraction as F; x=4-2*F(1,$p); print(x.numerator,x.denominator)")
    M5_Kp="$HERE/M5_ktm_1_$p.m"
    [ -s "$M5_Kp" ] || sed -E "s/\bd\b/($NUM\/$DEN)/g" "$M5KTM" > "$M5_Kp"
    RB5="$HERE/rebase_M5_1_$p.m"
    [ -s "$RB5" ] || python3 "$HERE/rebase_invert.py" --m1-masters "$M1MAST" \
        --src-masters "$M5MAST" --src-ktm "$M5_Kp" --out "$RB5" >/dev/null 2>&1
    SRCARGS+=(--source "M5:$M5_Kp")
    if [ ${#RBARGS[@]} -gt 0 ]; then
      cat "$RB5" >> "${RBARGS[1]}" 2>/dev/null || true
    else
      RBARGS+=(--rebase "$RB5")
    fi
  fi

  if [ ${#SRCARGS[@]} -eq 0 ]; then
    echo "[farm] $p: no source slices ready — skip for now"; continue
  fi

  UDS="$HERE/uds_1_$p"
  python3 "$HERE/extract_shared_rows.py" \
    --m1-masters "$M1MAST" --m1-target "$M1TGT" \
    "${SRCARGS[@]}" "${RBARGS[@]}" --out "$UDS" \
    > "$HERE/extract_1_$p.log" 2>&1

  # wait for slot
  while [ "$(running)" -ge "$N_CONC" ]; do sleep 60; done
  echo "[farm] fire $p ($(running) running; shared=$(wc -l < $UDS.lhs.txt))"
  bash "$HERE/build_M1_injected.sh" "$d" "1/$p" "$UDS" "$PAR" "$MEM_GB" "$MODE" \
    || { echo "[farm] build $p FAILED"; continue; }

  # background merge-on-completion watcher (per-slice)
  ( while ps -p "$(cat $d/kira.pid 2>/dev/null)" >/dev/null 2>&1; do sleep 60; done
    [ -s "$d/results/$FAM/kira_target.m" ] && bash "$HERE/merge_M1_injected.sh" "$d"
  ) &
done

echo "[farm] all primes dispatched; waiting for drain"
while [ "$(running)" -gt 0 ]; do sleep 120; done
NDONE=$(ls $HERE/M1inj_1_*/results/$FAM/kira_target.residual.m 2>/dev/null | wc -l)
echo "[farm] FARM COMPLETE @ $(date +%H:%M%Z): $NDONE/49 slices merged"
touch $HERE/M1_FARM_COMPLETE.marker

# ---- reconstruct + inject_ibpcache path -----------------------------------
# Downstream (not auto-run here): reconstruct one kira_target.m in (d, eta)
# from the merged slices $HERE/M1inj_1_*/results/$FAM/kira_target.m with your
# Thiele/Newton reconstruction, then
#   bash $BASE/inject_ibpcache.sh <reconstructed kira_target.m> \
#       <amflow log of the stalled run> <its in.json> <its out.json>

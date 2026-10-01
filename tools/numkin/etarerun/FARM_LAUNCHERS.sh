#!/bin/bash
# FARM_LAUNCHERS.sh — detached (setsid+nohup) launcher templates for the
# etarerun routes; each job writes DST/kira.pid. Set the ETANUMD_* env vars
# before use.
set -euo pipefail
# Launcher templates: every input dir comes from the ETANUMD_* env vars —
# set them for your own tree before use. No defaults are shipped.
# WORKED EXAMPLE: the master labels (M89, M6, M1), core counts and memory
# caps below are those of one three-master reduction; substitute your own.
# JEMALLOC and FERMATPATH override the allocator and Fermat binary the
# etanumd_*.sh verbs preload (Debian/Ubuntu paths are the fallback).
HERE=${ETANUMD_DIR:?set ETANUMD_DIR (etanumd work dir)}
CONV=$HERE/etanumd_convert.sh
FF2=$HERE/etanumd_2varff.sh
SRC89=${ETANUMD_SRC89:?set ETANUMD_SRC89 (M89 target_reduce dir)}
SRC6=${ETANUMD_SRC6:?set ETANUMD_SRC6 (M6 target_reduce dir)}
SRC1=${ETANUMD_SRC1:?set ETANUMD_SRC1 (M1 target_reduce dir)}
FAM=${ETANUMD_FAM:-}   # Kira family name — names the results/<FAM>/kira_target.m output

case "${1:-help}" in

  route_aprime_M89)
    # ONE job: 2-var (d,eta) FireFly → exact A(d,η). 12 cores, 30GB.
    "$FF2" "$SRC89" "$HERE/M89_2varFF" 12 30
    ;;

  route_aprime_M6)
    "$FF2" "$SRC6" "$HERE/M6_2varFF" 12 40
    ;;

  route_aprime_M1)
    "$FF2" "$SRC1" "$HERE/M1_2varFF" 12 100
    ;;

  route_a_M89_farm)
    # Extended Route A: 50 fit + 3 verify small-prime d-nodes (deg_d≈47).
    # ~23 min/slice × 53 / N_CONC. Fire in batches of N_CONC=8.
    # eps=1/p, p prime: 50 fit primes + 3 verify primes.
    PRIMES=(7 11 13 17 19 23 29 31 37 41 43 47 53 59 61 67 71 73 79 83 89 97 101 103 107 109 113 127 131 137 139 149 151 157 163 167 173 179 181 191 193 197 199 211 223 227 229 233 239 241 251 257 263)
    : "${FAM:?route_a_M89_farm needs ETANUMD_FAM (kira family name)}"
    N_CONC=${2:-8}; i=0
    for p in "${PRIMES[@]}"; do
      d="$HERE/M89_eps_1_$p"
      [ -f "$d/results/$FAM/kira_target.m" ] && { echo "skip $p (done)"; continue; }
      [ -d "$d/ff_save" ] && { echo "skip $p (in-progress/partial)"; continue; }
      [ -d "$d" ] && kill -0 "$(cat $d/kira.pid 2>/dev/null)" 2>/dev/null && { echo "skip $p (running)"; continue; }
      "$CONV" "$SRC89" "$d" "1/$p" 6 15 || { echo "convert $p FAILED (continuing)"; continue; }
      # one batch takes ~23 min on the example system; adjust the sleep to your slice wall
      i=$((i+1)); [ $((i % N_CONC)) = 0 ] && { echo "batch of $N_CONC fired; sleep 1400s"; sleep 1400; }
    done
    ;;

  route_b_M89_farm)
    # Route B: geometric small-eps grid for fixed-eps transport + Laurent.
    # eps = 2^k / 160000, k=0..14 (ratio 2, eps_max=2^14/160000≈0.1).
    # ~33 min/slice (large-d). NOTE: the boundary data must be built
    # separately (the amflow boundary stage).
    for k in $(seq 0 14); do
      num=$((1<<k)); "$CONV" "$SRC89" "$HERE/M89_epsB_${num}_160000" "$num/160000" 6 15
    done
    ;;

  M1_probe_80g)
    # M1 rate probe with an 80 GB cap (a 1.78M-equation system does not load under 30 GB)
    "$CONV" "$SRC1" "$HERE/M1_eps_1_7_80g" "1/7" 8 80
    ;;

  M6_probe)
    "$CONV" "$SRC6" "$HERE/M6_eps_1_7" "1/7" 8 25
    ;;

  status)
    set +e
    for d in "$HERE"/M89_eps_* "$HERE"/M89_2varFF "$HERE"/M6_* "$HERE"/M1_*; do
      [ -d "$d" ] || continue
      pid=$(cat "$d/kira.pid" 2>/dev/null)
      A=$(kill -0 "$pid" 2>/dev/null && echo RUN || echo -)
      NP=$(grep -c 'Promote to new prime' "$d"/kira_*.log 2>/dev/null | tail -1)
      DN=$(grep -aoE 'Done: [0-9]+ / [0-9]+' "$d"/kira_*.log 2>/dev/null | tail -1)
      OUT=$(stat -c %s "$d/results/$FAM/kira_target.m" 2>/dev/null)
      printf "%-30s %s pid=%-8s primes=%-3s %-24s out=%s\n" "$(basename $d)" "$A" "${pid:-?}" "${NP:-?}" "${DN:-}" "${OUT:-none}"
    done
    ;;

  *)
    echo "Usage: $0 {route_aprime_M89|route_aprime_M6|route_aprime_M1|route_a_M89_farm [N_CONC]|route_b_M89_farm|M1_probe_80g|M6_probe|status}"
    ;;
esac

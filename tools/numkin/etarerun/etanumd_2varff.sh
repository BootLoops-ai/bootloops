#!/bin/bash
[ "$(uname -s)" = Linux ] || { echo "FATAL: $(basename "$0") requires Linux — its kira-liveness checks read /proc" >&2; exit 2; }
# etanumd_2varff.sh SRC_TARGET_REDUCE DST [PARALLEL] [MEM_GB] [--nolaunch]
#
# Route A': re-run the amflow η-DE target_reduce with FireFly instead of Fermat,
# keeping BOTH (d, eta) symbolic. NO substitution — uses the original SYSTEM_*.gz.
# One job → exact A(d,η) → drop-in kira_target.m for amflow ibp-cache injection.
#
# This is the simplest bypass IF 2-var FF cost is acceptable. Measured M89 params:
# deg_d≈47, deg_eta≈56, 23671 fns, black-box≈7.3s. 2-var Zippel probe count is
# NOT predictable a priori — MUST rate-probe (rule#2).
#
set -euo pipefail
SRC=${1:?need SRC}; DST=${2:?need DST}; PAR=${3:-12}; MEM_GB=${4:-30}
NOLAUNCH=0; [ "${5:-}" = "--nolaunch" ] && NOLAUNCH=1
[ -d "$SRC/tmp" ] || { echo "no $SRC/tmp"; exit 1; }
rm -rf "$DST"; mkdir -p "$DST"
for f in config sectormappings preferred preferred_masters target de_targets; do
  [ -e "$SRC/$f" ] && cp -a "$SRC/$f" "$DST/"
done
mkdir -p "$DST/results"
[ -f "$SRC/results/preferredMasters" ] && cp -a "$SRC/results/preferredMasters" "$DST/results/"
for fam in "$SRC"/results/*/; do
  fn=$(basename "$fam"); [ -f "$fam/masters" ] && { mkdir -p "$DST/results/$fn"; cp -a "$fam/masters" "$DST/results/$fn/"; }
done
# tmp: copy SYSTEM_*.gz + SYSTEMconfig + masters VERBATIM (no VER_*)
for fam in "$SRC"/tmp/*/; do
  fn=$(basename "$fam"); mkdir -p "$DST/tmp/$fn"
  for aux in SYSTEMconfig masters; do [ -f "$fam$aux" ] && cp -a "$fam$aux" "$DST/tmp/$fn/"; done
  for g in "$fam"SYSTEM_*.gz; do cp -a "$g" "$DST/tmp/$fn/"; done
done
# jobs.yaml: run_firefly:true
python3 - "$SRC/jobs.yaml" "$DST/jobs.yaml" <<'PY'
import sys,re
L=open(sys.argv[1]).read().splitlines(); O=[]; ff=False
for ln in L:
  s=ln.strip()
  if s.startswith('run_triangular') or s.startswith('run_back_substitution'): continue
  if s.startswith('run_firefly'): O.append(re.match(r'^(\s*)',ln).group(1)+'run_firefly: true'); ff=True; continue
  O.append(ln)
  if s.startswith('run_initiate') and not ff:
    O.append(re.match(r'^(\s*)',ln).group(1)+'run_firefly: true'); ff=True
open(sys.argv[2],'w').write('\n'.join(O)+'\n')
PY
echo "[2varff] built $DST (2-var d,eta FireFly)"
[ "$NOLAUNCH" = 1 ] && exit 0
JEMALLOC=${JEMALLOC:-/usr/lib/x86_64-linux-gnu/libjemalloc.so.2}  # env override; Debian/Ubuntu path is the fallback
[ -f "$JEMALLOC" ] || JEMALLOC=/lib/x86_64-linux-gnu/libjemalloc.so.2
FERMAT=${FERMATPATH:-fer64}  # fer64 (Fermat) on PATH or via FERMATPATH
cd "$DST"
env LD_PRELOAD="$JEMALLOC" FERMATPATH="$FERMAT" OPENBLAS_NUM_THREADS=2 \
  setsid nohup nice -n 5 prlimit --as=$((MEM_GB*1024*1024*1024)) \
  kira jobs.yaml --parallel="$PAR" > kira_2varff.log 2>&1 < /dev/null &
sleep 2
REAL=$(for p in $(pgrep -x kira); do [ "$(readlink -f /proc/$p/cwd 2>/dev/null)" = "$(readlink -f "$DST")" ] && echo $p; done | head -1)
echo "${REAL:-?}" > kira.pid
echo "[2varff] launched kira pid=${REAL:-?} in $DST (parallel=$PAR mem=${MEM_GB}GB)"

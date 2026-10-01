#!/bin/bash
# S1+S3 worker: for each "channel sigma pair" line: stages 1-2 (slice_run.py, resume-safe), semantics gate (slice_gate12.py), stage 3 + piece-oracle gate (stage3.py)
set -euo pipefail
# Relocatable: the code tree is resolved from this script's own location; receipts, ckpt/ and logs/ go under
# WORK = $SURD_WORK (default ./surd_work under the current directory; never inside the package tree).
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd); WORK=$(readlink -f "${SURD_WORK:-$PWD/surd_work}"); export SURD_WORK="${WORK:?}"
case "${WORK:?}/" in "$(dirname "$(dirname "${HERE:?}")")"/*) echo "s13_worker.sh: SURD_WORK=${WORK} is inside the package tree; use a scratch directory" >&2; exit 2;; esac
LIST=$(readlink -f "${1:?list file: one 'channel sigma pair' per line}"); NAME=${2:?worker name}; TS3=${3:-1/2,4/5}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
ulimit -v "${SURD_VMEM_KB:-8388608}" 2>/dev/null || true          # address-space cap per process (8 GB default; every measured group needs < 1 GB)
if [ -x /usr/bin/time ] && /usr/bin/time -v true >/dev/null 2>&1; then TIMEV="/usr/bin/time -v -o"; else TIMEV=""; fi   # GNU time resource log when available
tv() { if [ -n "$TIMEV" ]; then $TIMEV "$1" "${@:2}"; else "${@:2}"; fi; }
mkdir -p "${WORK:?}/logs"; cd "${HERE:?}"
echo "$(date -u +%Y-%m-%dT%H:%M:%SZ) pid=$$ tag=s13_worker:$NAME argv=$LIST" >> "${WORK:?}/logs/pids.log"
while read -r CH SG PR; do
  [ -n "${CH:-}" ] || continue
  G=$(python3 -c "import slice_common as sc; print(sc.CHART['$PR'])")
  TAG="${CH}_s${SG}_p${PR}_g${G}"
  if [ -f "${WORK:?}/SLICE_S3_${TAG}.json" ]; then echo "skip $TAG (S3 done)"; continue; fi
  if [ ! -f "${WORK:?}/SLICE_S12_${TAG}.json" ]; then
    echo "$(date -u +%H:%M:%S) S1 start $TAG"
    tv "${WORK:?}/logs/s12_${TAG}_time.txt" nice -n 19 python3 slice_run.py --channel "$CH" --sigma "$SG" --pair "$PR" >> "${WORK:?}/logs/s12_${TAG}.out" 2>&1 || { echo "FAIL S1 $TAG rc=$?"; continue; }
  fi
  if ! grep -q '"PASS(>=30)": true' "${WORK:?}/SLICE_G12_${TAG}.json" 2>/dev/null; then
    nice -n 19 python3 slice_gate12.py --tag "$TAG" >> "${WORK:?}/logs/g12_${TAG}.out" 2>&1 || { echo "FAIL gate $TAG rc=$?"; continue; }
    echo "$(date -u +%H:%M:%S) gate $TAG: $(tail -1 ${WORK:?}/logs/g12_${TAG}.out | grep -o 'min digits.*')"
  fi
  echo "$(date -u +%H:%M:%S) S3 start $TAG"
  tv "${WORK:?}/logs/s3_${TAG}_time.txt" nice -n 19 python3 stage3.py --tag "$TAG" --t "$TS3" >> "${WORK:?}/logs/s3_${TAG}.out" 2>&1 || { echo "FAIL S3 $TAG rc=$?"; continue; }
  echo "$(date -u +%H:%M:%S) S3 done $TAG: $(grep -o 'min digits.*' ${WORK:?}/logs/s3_${TAG}.out | tail -1)"
done < "$LIST"
echo "$(date -u +%H:%M:%S) worker $NAME finished"

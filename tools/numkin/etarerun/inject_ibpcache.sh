#!/bin/bash
# inject_ibpcache.sh KIRA_TARGET_M AMFLOW_LOG IN_JSON OUT_JSON
# Fabricate an ibp-cache entry from a Route-A'/A reconstructed kira_target.m so
# that a fresh amflow run HITs on the stalled target_reduce step.
#
# Method (checkpoint-resume recipe): find the MISS key in the amflow
# log, locate the ibp-cache dir (AMFLOW_IBP_CACHE default), write kira_target.m
# under <cache>/<key>/results/<fam>/kira_target.m + minimal masters.
#
set -euo pipefail
KTM=${1:?kira_target.m}; ALOG=${2:?amflow log}; IN=${3:?in.json}; OUT=${4:?out.json}
KEY=$(grep -oE 'MISS key=[0-9a-f]+' "$ALOG" | tail -1 | cut -d= -f2)
[ -n "$KEY" ] || { echo "no MISS key in $ALOG"; exit 1; }
CACHE=${AMFLOW_IBP_CACHE_DIR:?set AMFLOW_IBP_CACHE_DIR (the amflow ibp-cache root)}
FAM=$(python3 -c "import json;print(json.load(open('$IN'))['family'])" 2>/dev/null || echo "${ETANUMD_FAM:?input json lacks 'family' — set ETANUMD_FAM}")
echo "[inject] key=$KEY fam=$FAM cache=$CACHE"
mkdir -p "$CACHE/$KEY/results/$FAM"
cp -a "$KTM" "$CACHE/$KEY/results/$FAM/kira_target.m"
# also stash masters if adjacent
MDIR=$(dirname "$KTM")
[ -f "$MDIR/masters" ] && cp -a "$MDIR/masters" "$CACHE/$KEY/results/$FAM/"
echo "[inject] wrote $CACHE/$KEY/results/$FAM/kira_target.m ($(stat -c %s "$KTM") bytes)"
echo "[inject] relaunch amflow with the SAME env block as the stalled run (see the header of $ALOG):"
echo "  <your amflow launcher> amflow_cli $IN $OUT"

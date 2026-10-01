#!/bin/bash
[ "$(uname -s)" = Linux ] || { echo "FATAL: $(basename "$0") requires Linux — its kira-liveness checks read /proc" >&2; exit 2; }
# build_M1_injected.sh  DST  EPS_PQ  UDS_PREFIX  [PARALLEL] [MEM_GB] [MODE] [--nolaunch]
#
# Build an M1 numeric-d η-DE slice at eps=P/Q with subsector-injection:
#   - target restricted to the residual (M1_target \ shared_LHS)
#   - shared rows (already known from M5/M6/M89 slices) spliced post-hoc into
#     the final kira_target.m by merge_M1_injected.sh
#
# UDS_PREFIX is the --out prefix from extract_shared_rows.py, i.e.
#   ${UDS_PREFIX}.target_residual  ${UDS_PREFIX}.rows.m  ${UDS_PREFIX}.uds.kira
#
# MODE:
#   reuse    — reuse M1's existing tmp/SYSTEM_*.gz (etanumd_convert route).
#              bb time unchanged; win = fewer functions (~6% with M89-only).
#   fresh    — regenerate SYSTEM with `extra_relations: uds.kira`.  BLOCKED:
#              requires SYMBOLIC-(d,η) coeffs (numeric-d ⇒ over-determined,
#              uds_fmt_test).  Use only with M5's symbolic kira_target.m.
#   freshtgt — (WINNING ROUTE, probe-validated) regenerate SYSTEM with target
#              = residual only, NO extra_relations.  pyred select() prunes
#              n_eqns 1.78M→1.01M (hard-43 probe: bb 89.8→23.9s, nfunc
#              142696→24648).  ~12min init once → reusable SYSTEM for all
#              49 numeric-d slices via etanumd_convert.
#
# Copy-only (never touches SRC or running kiras).
#
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
# Template: input roots come from the ETANUMD_* env vars; no defaults shipped.
SRC=${ETANUMD_M1_SRC:?set ETANUMD_M1_SRC (M1 target_reduce dir)}
FAM=${ETANUMD_FAM:?set ETANUMD_FAM (kira family name)}
CONV=${ETANUMD_CONVERT:-"$(dirname "$0")/../etanumd_convert.sh"}

DST=${1:?need DST dir}
EPS_PQ=${2:?need eps as P/Q}
UDS=${3:?need UDS_PREFIX (from extract_shared_rows.py --out)}
PAR=${4:-16}
MEM_GB=${5:-60}
MODE=${6:-reuse}
NOLAUNCH=0; [ "${7:-}" = "--nolaunch" ] && NOLAUNCH=1

[ -f "$UDS.target_residual" ] || { echo "[inject] ERR: $UDS.target_residual missing"; exit 1; }
[ -f "$UDS.rows.m" ]          || { echo "[inject] ERR: $UDS.rows.m missing"; exit 1; }

NRES=$(grep -c '$FAM\[' "$UDS.target_residual")
NSHR=$(grep -c '^$FAM\[' "$UDS.lhs.txt")
echo "[inject] mode=$MODE  eps=$EPS_PQ  shared_rows=$NSHR  residual_targets=$NRES"

# ---- route 1: reuse existing SYSTEM (fast setup) ---------------------------
if [ "$MODE" = "reuse" ]; then
  # build the numd slice via etanumd_convert (SYSTEM d-substitute), no launch
  "$CONV" "$SRC" "$DST" "$EPS_PQ" "$PAR" "$MEM_GB" --nolaunch
  # overwrite target with residual list
  cp -a "$UDS.target_residual" "$DST/target"
  cp -a "$UDS.rows.m"          "$DST/shared_rows.m"
  cp -a "$UDS.lhs.txt"         "$DST/shared_lhs.txt"
  echo "reuse" > "$DST/INJECT_MODE"

# ---- route 2: fresh init with extra_relations ------------------------------
elif [ "$MODE" = "fresh" ]; then
  rm -rf "$DST"; mkdir -p "$DST"
  cp -a "$SRC/config" "$DST/"
  cp -a "$SRC/preferred" "$DST/"
  cp -a "$UDS.target_residual" "$DST/target"
  cp -a "$UDS.uds.kira"        "$DST/uds_shared.kira"
  cp -a "$UDS.rows.m"          "$DST/shared_rows.m"
  cp -a "$UDS.lhs.txt"         "$DST/shared_lhs.txt"
  # d-substitute inside uds file (numeric-d slice ⇒ coeffs must be eta-only)
  P=${EPS_PQ%%/*}; Q=${EPS_PQ##*/}
  read NUM DEN < <(python3 -c "from fractions import Fraction as F; d=4-2*F($P,$Q); print(d.numerator, d.denominator)")
  sed -i -E "s/\bd\b/($NUM\/$DEN)/g" "$DST/uds_shared.kira"
  # jobs.yaml: fresh init, run_firefly, extra_relations
  cat > "$DST/jobs.yaml" <<EOF
jobs:
 - reduce_sectors:
    reduce:
     - {topologies: [$FAM], sectors: [1022], r: 15, s: 0, d: 6}
    select_integrals:
      select_mandatory_list:
        - [$FAM, target]
    preferred_masters: preferred
    extra_relations: uds_shared.kira
    integral_ordering: 5
    run_initiate: true
    run_firefly: true
 - kira2math:
    target:
     - [$FAM, target]
EOF
  echo "fresh" > "$DST/INJECT_MODE"
  # provenance
  cat > "$DST/ETANUMD_SLICE.json" <<EOF
{"src":"$SRC","eps":"$EPS_PQ","d_num":$NUM,"d_den":$DEN,"mode":"fresh","uds":"$UDS","par":$PAR,"mem_gb":$MEM_GB,"built_utc":"$(date -u +%FT%TZ)"}
EOF

# ---- route 1b: fresh init, restricted target only (WINNING ROUTE) ---------
elif [ "$MODE" = "freshtgt" ]; then
  # Step 1: build (or reuse) the residual-SYSTEM template dir (symbolic d,η).
  # This is EPS-INDEPENDENT — do it once, then etanumd_convert per slice.
  TMPL="$HERE/M1_residual_SYSTEM"
  UDSHASH=$(md5sum "$UDS.target_residual" | cut -c1-8)
  if [ ! -f "$TMPL/tmp/$FAM/SYSTEMconfig" ] || \
     [ "$(cat $TMPL/UDSHASH 2>/dev/null)" != "$UDSHASH" ]; then
    echo "[inject] building residual-SYSTEM template (fresh init, ~12min once)"
    rm -rf "$TMPL"; mkdir -p "$TMPL"
    cp -a "$SRC/config" "$TMPL/"
    cp -a "$SRC/preferred" "$TMPL/"
    cp -a "$UDS.target_residual" "$TMPL/target"
    echo "$UDSHASH" > "$TMPL/UDSHASH"
    cat > "$TMPL/jobs.yaml" <<JEOF
jobs:
 - reduce_sectors:
    reduce:
     - {topologies: [$FAM], sectors: [1022], r: 15, s: 0, d: 6}
    select_integrals:
      select_mandatory_list:
        - [$FAM, target]
    preferred_masters: preferred
    integral_ordering: 5
    run_initiate: true
JEOF
    ( cd "$TMPL" && env LD_PRELOAD="${JEMALLOC:-/usr/lib/x86_64-linux-gnu/libjemalloc.so.2}" \
        FERMATPATH=${FERMATPATH:-fer64} \
        kira jobs.yaml --parallel="$PAR" > kira_init.log 2>&1 )
    [ -s "$TMPL/tmp/$FAM/SYSTEMconfig" ] || { echo "[inject] ERR: template init failed"; tail -20 "$TMPL/kira_init.log"; exit 2; }
    echo "[inject] template SYSTEM: $(cat $TMPL/tmp/$FAM/SYSTEMconfig) eqns, $(wc -l < $TMPL/results/$FAM/masters) masters"
  fi
  # Step 2: numeric-d slice from the template
  "$CONV" "$TMPL" "$DST" "$EPS_PQ" "$PAR" "$MEM_GB" --nolaunch
  cp -a "$UDS.rows.m"  "$DST/shared_rows.m"
  cp -a "$UDS.lhs.txt" "$DST/shared_lhs.txt"
  echo "freshtgt" > "$DST/INJECT_MODE"

else
  echo "[inject] ERR: MODE must be reuse|fresh|freshtgt"; exit 1
fi

# stamp injection metadata
cat > "$DST/INJECT_META.json" <<EOF
{"mode":"$MODE","eps":"$EPS_PQ","uds_prefix":"$UDS","n_shared_rows":$NSHR,"n_residual_targets":$NRES,"par":$PAR,"mem_gb":$MEM_GB,"built_utc":"$(date -u +%FT%TZ)"}
EOF

if [ "$NOLAUNCH" = 1 ]; then
  echo "[inject] --nolaunch: DST prepared at $DST"
  exit 0
fi

# ---- launch (setsid+jemalloc) ----------------------------------------------
JEMALLOC=${JEMALLOC:-/usr/lib/x86_64-linux-gnu/libjemalloc.so.2}  # env override; Debian/Ubuntu path is the fallback
[ -f "$JEMALLOC" ] || JEMALLOC=/lib/x86_64-linux-gnu/libjemalloc.so.2
FERMAT=${FERMATPATH:-fer64}
MEM_KB=$((MEM_GB*1024*1024))
cd "$DST"
env LD_PRELOAD="$JEMALLOC" FERMATPATH="$FERMAT" OPENBLAS_NUM_THREADS=2 \
  setsid nohup nice -n 5 prlimit --as=$((MEM_KB*1024)) \
  kira jobs.yaml --parallel="$PAR" > kira_numd.log 2>&1 < /dev/null &
sleep 2
REAL=$(for p in $(pgrep -x kira); do [ "$(readlink -f /proc/$p/cwd 2>/dev/null)" = "$(readlink -f "$DST")" ] && echo $p; done | head -1)
echo "[inject] launched kira in $DST (real pid=${REAL:-?}, parallel=$PAR, mem_cap=${MEM_GB}GB, mode=$MODE)"
echo "${REAL:-}" > "$DST/kira.pid"

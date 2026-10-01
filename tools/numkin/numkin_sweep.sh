#!/bin/bash
# numkin_sweep.sh — generalized numeric-kinematic point sweep for stuck multi-var
# symbolic FireFly reductions; works for any family/variable.
#
# Method: sed-substitute a WORD-BOUNDARY numeric value for one symbolic
# kinematic invariant directly into ALREADY-GENERATED SYSTEM_<FAM>_*.gz files
# (never regenerate — reuses the expensive Kira gen phase's sunk cost), leaving
# the remaining variable(s) (typically `d`) symbolic. Then runs ordinary
# run_initiate:false FireFly on the now-lower-variable-count system: turns an
# intractable N-var symbolic solve into many cheap (N-1)-var solves. Reconstruct
# the substituted variable's dependence afterward via exact Thiele/Padé fitting
# with held-out points (numkin_harvest.py in this package is the
# reconstruction half).
#
# When to reach for this: multi-var symbolic FireFly is "hours out" / stuck on
# a subset of high-degree entries, but single-var (d-only) FireFly on the same
# system is fast. See GUIDE.md.
#
# Usage:
#   numkin_sweep.sh stage <FAM> <SYS_DIR> <VAR> <VALUE> <OUT_DIR> <PT_ID> \
#                    <PREFERRED_MASTERS> <TARGETS_FILE> <REDUCE_YAML_FRAGMENT>
#   numkin_sweep.sh run   <OUT_DIR> <PT_ID> <NTHREADS> [KIRA_BIN] [FERMATPATH]
#
#   FAM                 topology name (as in the SYSTEM/results filenames)
#   SYS_DIR             dir containing tmp/<FAM>/SYSTEM_<FAM>_*.gz (+ SYSTEMconfig,
#                       masters) from a completed `run_initiate:true` gen pass,
#                       and config/ + sectormappings/ siblings
#   VAR                 symbol to substitute, e.g. m2 (NOT d -- d stays symbolic)
#   VALUE               numeric value (rational as "p/q" or integer), inserted
#                       as a parenthesized literal: sed s|\bVAR\b|(VALUE)|g
#   OUT_DIR             where per-point pt_<PT_ID>/ dirs are staged
#   PT_ID               point identifier (any string, used as dir suffix)
#   PREFERRED_MASTERS   path to a preferred_masters file to copy in
#   TARGETS_FILE        path to a target-list file to copy in (also becomes the
#                       select_mandatory_list + kira2math target name)
#   REDUCE_YAML_FRAGMENT path to a YAML fragment for the `reduce:` block (one
#                       or more `{topologies: [...], sectors: [...], r: .., s: ..}`
#                       entries, same syntax as a jobs.yaml reduce_sectors.reduce list)
#
# Example (single reduce-spec inline):
#   cat > reduce.yaml <<'EOF'
#   - {topologies: [box3L], sectors: [1023], r: 13, s: 1, d: 3}
#   - {topologies: [box3L], sectors: [170], r: 10, s: 1}
#   EOF
#   numkin_sweep.sh stage box3L /scratch/gen_pass_dir \
#     m2 101 /scratch/sec1023_numkin 0 \
#     preferred_masters de_targets_top30 reduce.yaml
#   numkin_sweep.sh run /scratch/sec1023_numkin 0 8
set -uo pipefail

stage(){
  local FAM=$1 SRC=$2 VAR=$3 VALUE=$4 OUT=$5 PT=$6 PM=$7 TGT=$8 REDUCE_YAML=$9
  local D="$OUT/pt_$PT"
  local TGTNAME
  TGTNAME=$(basename "$TGT")
  rm -rf "$D"; mkdir -p "$D/tmp/$FAM" "$D/results/$FAM"
  # COPY (not symlink) config so we can strip the substituted variable from
  # kinematic_invariants -- leaving it in place makes FireFly waste effort
  # certifying that a value we already eliminated by sed is "trivial" (the
  # tell in the log: "Reconstructing in: d, m2" instead of just "d" after
  # m2 was already substituted out of every equation).
  cp -r "$SRC/config" "$D/config"
  python3 -c "
import re, sys
p = '$D/config/kinematics.yaml'
txt = open(p).read()
# drop any 'kinematic_invariants' list entry whose first element is the substituted var
txt2 = re.sub(r'\n\s*-\s*\[\s*${VAR}\s*,[^\]]*\]', '', txt)
open(p, 'w').write(txt2)
" 2>/dev/null || echo "[numkin] WARNING: could not strip $VAR from kinematics.yaml (config format may differ; check manually)"
  # ALSO substitute VAR in integralfamilies.yaml (propagator masses). Kira builds
  # the FireFly variable registry from config symbols, NOT just kinematic_invariants:
  # leaving symbolic $VAR in propagators makes FireFly run MULTI-VAR interpolation
  # ("Reconstructing in: d, m2") on equations that no longer contain $VAR at all --
  # different (slower) interpolation strategy + wasted triviality probes.
  # (The tell: 2-var mode despite an m2-free system.)
  sed -i -E "s|\\b${VAR}\\b|(${VALUE})|g" "$D/config/integralfamilies.yaml" 2>/dev/null || true
  ln -s "$SRC/sectormappings" "$D/sectormappings"
  cp "$PM" "$D/preferred_masters"
  cp "$TGT" "$D/$TGTNAME"
  {
    echo "jobs:"
    echo "  - reduce_sectors:"
    echo "      reduce:"
    sed 's/^/        /' "$REDUCE_YAML"
    echo "      select_integrals:"
    echo "        select_mandatory_list:"
    echo "          - [$FAM, $TGTNAME]"
    echo "      preferred_masters: preferred_masters"
    echo "      integral_ordering: 8"
    echo "      run_initiate: false"
    echo "      run_triangular: false"
    echo "      run_firefly: true"
    echo "  - kira2math:"
    echo "      target:"
    echo "        - [$FAM, $TGTNAME]"
  } > "$D/jobs.yaml"
  # Parallel sed-substitution: files are independent; a serial loop costs
  # 15-20 min of pure staging overhead PER POINT (measured).
  # NPAR kept modest -- staging shares the machine with other jobs.
  local NPAR=${NUMKIN_STAGE_PAR:-12}
  ls "$SRC/tmp/$FAM"/SYSTEM_${FAM}_*.gz | \
    xargs -P "$NPAR" -I{} bash -c \
      'zcat "$1" | sed -E "s|\\b'"${VAR}"'\\b|('"${VALUE}"')|g" | gzip > "'"$D/tmp/$FAM"'/$(basename "$1")"' _ {}
  local n
  n=$(ls "$D/tmp/$FAM"/SYSTEM_${FAM}_*.gz 2>/dev/null | wc -l)
  cp "$SRC/tmp/$FAM/SYSTEMconfig" "$SRC/tmp/$FAM/masters" "$D/tmp/$FAM/" 2>/dev/null
  echo "$VALUE" > "$D/${VAR}val"
  echo "[numkin] staged pt_$PT ${VAR}=${VALUE} ($n SYSTEM files) -> $D"
}

run_pt(){
  local OUT=$1 PT=$2 NT=$3
  local KIRA=${4:-/usr/local/bin/kira}
  local FERMAT=${5:-${FERMATPATH:-fer64}}
  local D="$OUT/pt_$PT"
  date +%s > "$D/t0"
  env FERMATPATH="$FERMAT" \
    bash -c "ulimit -v unlimited; cd '$D' && exec $KIRA -p$NT jobs.yaml" \
    > "$D/run.log" 2>&1
  date +%s > "$D/t1"
  echo "[numkin] pt_$PT done in $(( $(cat "$D/t1") - $(cat "$D/t0") ))s -- $(ls "$D"/results/*/*.m 2>/dev/null | head -1)"
}

case "${1:-}" in
  stage) shift; stage "$@" ;;
  run)   shift; run_pt "$@" ;;
  *) sed -n '1,45p' "$0" | grep '^#' | sed 's/^# \?//'; exit 1 ;;
esac

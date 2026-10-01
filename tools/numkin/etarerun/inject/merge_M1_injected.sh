#!/bin/bash
# merge_M1_injected.sh  SLICE_DIR
#
# After an injected M1 slice's kira finishes, splice the shared rows
# ($SLICE_DIR/shared_rows.m, in kira_target.m fragment format) into the kira
# output → full 503-row kira_target.m at $SLICE_DIR/results/$FAM/.
#
# Idempotent: writes kira_target.m from kira_target.residual.m + shared_rows.m.
#
set -euo pipefail
D=${1:?need SLICE_DIR}
FAM=${ETANUMD_FAM:?set ETANUMD_FAM (kira family name)}
K=$D/results/$FAM/kira_target.m
S=$D/shared_rows.m
[ -s "$K" ] || { echo "[merge] ERR: $K missing/empty"; exit 1; }
[ -s "$S" ] || { echo "[merge] ERR: $S missing"; exit 1; }
# stash kira's own output as .residual.m (once)
[ -f "$D/results/$FAM/kira_target.residual.m" ] || cp -a "$K" "$D/results/$FAM/kira_target.residual.m"
R=$D/results/$FAM/kira_target.residual.m
# assemble: { residual_body , shared_body }  (both bodies end with ',')
{
  echo "{"
  # residual body: strip leading '{' and trailing '}' (kira2math wraps)
  sed -e '1s/^{//' -e '$s/}$//' "$R"
  # shared body already fragment (ends each row with ',')
  cat "$S"
  echo "}"
} > "$K"
NR=$(grep -c ' -> ' "$R"); NS=$(grep -c ' -> ' "$S"); NM=$(grep -c ' -> ' "$K")
echo "[merge] $D: residual=$NR + shared=$NS = merged=$NM rows"
[ "$NM" = "$((NR+NS))" ] || { echo "[merge] WARN: row-count mismatch"; }

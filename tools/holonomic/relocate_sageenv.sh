#!/bin/bash
# relocate_sageenv.sh — RELOCATED-SAGEENV member (holonomic family):
# relocate the pinned SageMath/ore_algebra environment into a destination
# directory and GATE the relocation before any science runs on it.
#
# THE PATTERN:
#   1. rsync the pinned env INTO the destination dir (source READ-ONLY,
#      untouched).
#   2. STALE-SHEBANG CAVEAT: the copy's bin/* scripts carry stale absolute
#      shebangs (#!/old-prefix/bin/python...) — sed-patch them IN THE
#      RELOCATED COPY ONLY.
#   3. Invoke the copy's own versioned bin python DIRECTLY, bypassing the
#      bin/sage wrapper (RPATH $ORIGIN/../lib; sys.prefix resolves to the
#      relocated copy).
#   4. SMOKE GATE: import sage.all + ore_algebra, engine versions equal the
#      pin (10.7 / 0.5), IC-convention probe, toy transition matrix.
#      Expected pinned warnings: the 'naive is the only working path'
#      Cython-fallback warning (reproduced verbatim = GOOD); a Singular
#      data-path warning appears at sage.all import — ore_algebra
#      numerical_transition_matrix does not route through Singular.
#   5. VERIFY GATE (relocation-correctness): re-run a transport on the
#      relocated env and require BIT-IDENTITY against the origin env —
#      deterministic 'naive' summation makes bit-identity, not mere
#      closeness, the right gate grade. This helper runs the toy transport
#      in verify_bit_identity.py on BOTH envs and compares the
#      full-precision serializations bit for bit. Before trusting deep
#      production runs on a relocated env, ALSO replicate a previously
#      recorded production receipt on it.
#
# Usage:
#   relocate_sageenv.sh SRC_ENV DEST_DIR HOLONOMIC_DIR [RECEIPT_OUT]
#     SRC_ENV       pinned env root (READ-ONLY; never pip-install into it)
#     DEST_DIR      destination dir; env lands at DEST_DIR/sageenv
#     HOLONOMIC_DIR dir containing holonomic_transport.py (this dir)
#     RECEIPT_OUT   receipt JSON path (default DEST_DIR/RECEIPT_relocation.json)
#   SAGEENV_PYBIN overrides the venv python basename (default: the first
#   versioned bin/python3.X found in SRC_ENV).
# All verdict fields in the receipt are SCRIPT-EMITTED by this helper.
set -e
SRC="${1:?SRC_ENV}"; DEST="${2:?DEST_DIR}"; HDIR="${3:?HOLONOMIC_DIR}"
OUT="${4:-$DEST/RECEIPT_relocation.json}"
PYBIN="${SAGEENV_PYBIN:-}"
if [ -z "$PYBIN" ]; then
  for c in "$SRC"/bin/python3.[0-9] "$SRC"/bin/python3.[0-9][0-9]; do
    [ -x "$c" ] || continue
    PYBIN=$(basename "$c"); break
  done
fi
[ -n "$PYBIN" ] && [ -x "$SRC/bin/$PYBIN" ] || { echo "FATAL: no versioned bin/python3.X in $SRC (set SAGEENV_PYBIN)"; exit 2; }
mkdir -p "$DEST"
STAMP0=$(date -u)

echo "[1/5] rsync $SRC -> $DEST/sageenv ($(date -u))"
T0=$SECONDS
rsync -a "$SRC/" "$DEST/sageenv/"
RSYNC_S=$((SECONDS - T0))

echo "[2/5] stale-shebang patch (relocated copy only)"
NPATCH=0
for f in "$DEST/sageenv/bin/"*; do
  [ -f "$f" ] || continue
  head1=$(head -c 200 "$f" 2>/dev/null | head -1) || continue
  case "$head1" in
    "#!"*python*)
      case "$head1" in
        "#!$DEST/sageenv/bin/"*|"#!/usr/bin/env"*) : ;;  # already fine
        *)
          sed -i "1s|^#!.*python[0-9.]*|#!$DEST/sageenv/bin/$PYBIN|" "$f"
          NPATCH=$((NPATCH+1)) ;;
      esac ;;
  esac
done
echo "  patched $NPATCH stale shebangs"

PY="$DEST/sageenv/bin/$PYBIN"
export HOLONOMIC_SAGEENV="$DEST/sageenv" HOLONOMIC_DIR="$HDIR"

echo "[3/5] smoke gate on the relocated env (direct bin/$PYBIN)"
SMOKE_LOG="$DEST/smoke_relocated.log"
if "$PY" "$HDIR/smoke_relocated_env.py" > "$SMOKE_LOG" 2>&1 \
   && grep -q "SMOKE_OK" "$SMOKE_LOG"; then SMOKE=PASS; else SMOKE=FAIL; fi
echo "  smoke: $SMOKE (log: $SMOKE_LOG)"

echo "[4/5] verify gate: toy transport on ORIGIN env"
ORIG_OUT="$DEST/verify_origin.txt"
HOLONOMIC_SAGEENV="$SRC" "$SRC/bin/$PYBIN" \
  "$HDIR/verify_bit_identity.py" > "$ORIG_OUT" 2> "$DEST/verify_origin.err" \
  || { echo "FATAL: origin verify run failed"; exit 3; }

echo "[5/5] verify gate: same toy transport on RELOCATED env + bit compare"
RELO_OUT="$DEST/verify_relocated.txt"
"$PY" "$HDIR/verify_bit_identity.py" > "$RELO_OUT" 2> "$DEST/verify_relocated.err" \
  || { echo "FATAL: relocated verify run failed"; exit 3; }
SHA_O=$(grep '^VERIFY_SHA256 ' "$ORIG_OUT" | awk '{print $2}')
SHA_R=$(grep '^VERIFY_SHA256 ' "$RELO_OUT" | awk '{print $2}')
if [ -n "$SHA_O" ] && [ "$SHA_O" = "$SHA_R" ]; then BITID=PASS; else BITID=FAIL; fi
echo "  bit-identity: $BITID (origin $SHA_O vs relocated $SHA_R)"

VERDICT=FAIL
[ "$SMOKE" = PASS ] && [ "$BITID" = PASS ] && VERDICT=PASS
cat > "$OUT" <<EOF
{
 "PRODUCER": {
  "script": "relocate_sageenv.sh (RELOCATED-SAGEENV member, holonomic family)",
  "stamp_utc_start": "$STAMP0",
  "stamp_utc_end": "$(date -u)"
 },
 "src_env": "$SRC",
 "relocated_env": "$DEST/sageenv",
 "rsync_wall_s": $RSYNC_S,
 "stale_shebangs_patched": $NPATCH,
 "smoke_gate": "$SMOKE",
 "verify_sha256_origin": "$SHA_O",
 "verify_sha256_relocated": "$SHA_R",
 "bit_identity_gate": "$BITID",
 "verdict": "$VERDICT",
 "note": "verdict fields script-emitted; bit-identity = full-precision serialization of the toy certified transport (verify_bit_identity.py) equal on origin and relocated envs. Before trusting deep production runs on a relocated env, ALSO replicate a previously recorded production receipt on it."
}
EOF
echo "receipt: $OUT"
echo "RELOCATION_VERDICT $VERDICT"
[ "$VERDICT" = PASS ]

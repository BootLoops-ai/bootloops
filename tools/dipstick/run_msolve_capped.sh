#!/usr/bin/env bash
# run_msolve_capped.sh -- run msolve under a MANDATORY virtual-memory cap.
#
# WHY: msolve's F4 matrix allocation can balloon to hundreds of GB with no
# internal limiter; an address-space cap makes it die alone (malloc failure)
# instead of taking the machine down. Mirrors the run_amflow_jemalloc.sh pattern
# (env-explicit wrapper, exec, stderr provenance; see the AMFlow.cpp fork,
# the sibling repository amflow-cpp), but
# REFUSES to run uncapped: there is deliberately no default-cap path.
#
# CHAR-P COEFFICIENT FOOTGUN: msolve (0.6.5) silently garbles char-p input
# coefficients > 2^31 — wrong parse, rc 0, no warning; the solved system is
# then NOT the emitted system. EVERY .ms emitter feeding this wrapper MUST
# reduce coefficients mod p to [0, p) at emission and refuse out-of-range
# bytes. This wrapper does NOT re-check the bytes.
#
# CAP UNITS (known footgun): `ulimit -v` takes KILOBYTES. A literal like
# "31G" is INVALID.
# We take an integer GB count in MSOLVE_CAP_GB and convert: KB = GB*1024*1024.
#
# Usage:   MSOLVE_CAP_GB=<int> run_msolve_capped.sh <msolve args ...>
# Example: MSOLVE_CAP_GB=31 run_msolve_capped.sh -f sys.ms -o sys.out -t 8
# Optional: MSOLVE_BIN=/path/to/msolve (default: msolve on PATH).
#
# rc: 2 on refusal (no/invalid cap, missing binary); otherwise msolve's own rc
# (exec => rc propagates verbatim; killed-at-cap typically shows as nonzero /
# malloc-failure exit from msolve, not a machine-wide OOM).
set -euo pipefail

if [ -z "${MSOLVE_CAP_GB:-}" ]; then
  echo "[msolve-cap] REFUSED: MSOLVE_CAP_GB is not set." >&2
  echo "[msolve-cap] Set an integer GB cap, e.g.: MSOLVE_CAP_GB=31 $0 <msolve args ...>" >&2
  exit 2
fi

case "$MSOLVE_CAP_GB" in
  ''|*[!0-9]*)
    echo "[msolve-cap] REFUSED: MSOLVE_CAP_GB='$MSOLVE_CAP_GB' is not a plain positive integer (GB)." >&2
    echo "[msolve-cap] Suffixed literals like '31G' are the known ulimit footgun; use MSOLVE_CAP_GB=31." >&2
    exit 2
    ;;
esac
if [ "$MSOLVE_CAP_GB" -lt 1 ]; then
  echo "[msolve-cap] REFUSED: MSOLVE_CAP_GB must be >= 1 (got '$MSOLVE_CAP_GB')." >&2
  exit 2
fi

MSOLVE="${MSOLVE_BIN:-msolve}"
if ! command -v "$MSOLVE" >/dev/null 2>&1; then
  echo "[msolve-cap] REFUSED: msolve binary '$MSOLVE' not found (set MSOLVE_BIN?)." >&2
  exit 2
fi

CAP_KB=$((MSOLVE_CAP_GB * 1024 * 1024))

# Cap and exec in THIS shell so the limit binds the msolve process itself.
ulimit -v "$CAP_KB"

echo "[msolve-cap] ulimit -v ${CAP_KB} KB (= ${MSOLVE_CAP_GB} GB address-space cap)" >&2
echo "[msolve-cap] exec: $MSOLVE $*" >&2
exec "$MSOLVE" "$@"

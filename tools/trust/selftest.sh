#!/bin/sh
# selftest.sh — one command for the whole trust package (public legs).
#   leg 1  trust package: import gates + fail-closed vendor pins + linked
#          member homes (python3 -E -B -c "import trust; trust.verify()")
#   leg 2  receipt member: sh receipt/tests/run_tests.sh (test_core +
#          mutation_test self-contained; test_banked SKIPS by name unless
#          RECEIPT_BANKED_ROOT points at the archived reference artifacts)
#   leg 3  strata member: python3 -m pytest strata/tests/test_certify.py
#          (self-contained synthetic system; SKIPS by name if pytest or
#          numpy is not installed)
# The adversarial battery (battery/battery.py <scratch>) is the separate,
# reference-data leg documented in GUIDE.md / MANUAL.md.
# Usage: sh selftest.sh [SCRATCH_DIR]   (default: $TMPDIR/trust_selftest)
# Runs from any cwd; writes only under SCRATCH_DIR. Exit 0 = all legs PASS.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
OUT="${1:-${TMPDIR:-/tmp}/trust_selftest}"
mkdir -p "$OUT"
export PYTHONDONTWRITEBYTECODE=1
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}" OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}" MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
rc=0

echo "== leg 1: trust package — import gates + vendor pins + linked member homes =="
( cd "$HERE" && env -u PYTHONPATH -u PYTHONSTARTUP -u PYTHONHOME \
    python3 -E -B -c "
import sys, trust
r = trust.verify()
bad = r['vendor_bad'] or r['linked_missing']
print('vendor pins ok: %d/%d; linked homes: %s; missing: %s; advisory drift: %s'
      % (len(r['vendor_ok']), len(trust.PINS), sorted(trust.LINKED),
         r['linked_missing'] or 'none', r['linked_drift'] or 'none'))
sys.exit(1 if bad else 0)
" ) && echo "leg 1: PASS" || { echo "leg 1: FAIL"; rc=1; }

echo
echo "== leg 2: receipt member battery =="
sh "$HERE/receipt/tests/run_tests.sh" "$OUT/receipt" && echo "leg 2: PASS" || { echo "leg 2: FAIL"; rc=1; }

echo
echo "== leg 3: strata member certify gates =="
if python3 -c "import pytest, numpy" >/dev/null 2>&1; then
  ( cd "$OUT" && nice -n 10 python3 -m pytest -q -p no:cacheprovider \
      "$HERE/strata/tests/test_certify.py" ) && echo "leg 3: PASS" || { echo "leg 3: FAIL"; rc=1; }
else
  echo "SKIP: leg 3 (strata/tests/test_certify.py) needs pytest + numpy"
fi

echo
if [ "$rc" -eq 0 ]; then echo "TRUST SELFTEST: ALL LEGS PASS";
else echo "TRUST SELFTEST: FAILURES (see above)"; fi
exit "$rc"

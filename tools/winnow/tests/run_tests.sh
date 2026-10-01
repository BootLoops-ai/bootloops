#!/bin/bash
# ibplapper battery. Fast self-contained gates first; the archived-replay
# batteries follow and SKIP cleanly unless WINNOW_BANKED_ROOT /
# WINNOW_T1_ROOT point at the archived reference artifacts (not shipped
# with this tree). Any FAIL = nonzero exit.
# Usage: run_tests.sh [scratch_dir]   (default /tmp/ibplapper_tests)
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
SCRATCH=${1:-/tmp/ibplapper_tests}
mkdir -p "$SCRATCH"
RC=0
for t in test_vendor_integrity test_gates test_kira_adapter \
         test_user_system test_m2_coeffs test_tbanks test_dense_backend \
         test_witness84 test_native_witness test_sabotage; do
  echo "=== $t ==="
  python3 "$HERE/$t.py" "$SCRATCH/$t" || RC=1
done
echo
if [ $RC -eq 0 ]; then echo "ibplapper battery: ALL SUITES PASS";
else echo "ibplapper battery: FAILURES (see above)"; fi
exit $RC

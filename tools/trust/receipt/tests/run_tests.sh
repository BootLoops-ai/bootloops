#!/usr/bin/env bash
# run_tests.sh — full receipt test battery.
# Usage: run_tests.sh [OUT_DIR]   (default: $TMPDIR/receipt_tests)
# test_core + mutation_test are self-contained; test_banked exercises the
# archived internal reference artifacts and SKIPS unless RECEIPT_BANKED_ROOT
# points at them. Budget: ~5 min wall, single-core numpy (measured prior:
# 84-row retrofit ~103 s solve/prime). Runs nice 10, 1 BLAS thread.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
OUT="${1:-${TMPDIR:-/tmp}/receipt_tests}"
mkdir -p "$OUT"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
rc=0

echo "== test_core (synthetic + sabotage + rc discipline) =="
nice -n 10 python3 "$HERE/test_core.py" "$OUT/test_core" || rc=1

echo "== mutation_test (checker mutants on scratch copies) =="
nice -n 10 python3 "$HERE/mutation_test.py" "$OUT/mutation_scratch" || rc=1

echo "== test_banked (30 archived pairs + 84/84 retrofit + wrong-table MUST-FAIL; skips without RECEIPT_BANKED_ROOT) =="
nice -n 10 python3 "$HERE/test_banked.py" "$OUT/test_banked" || rc=1

echo
if [ "$rc" -eq 0 ]; then echo "RECEIPT TEST BATTERY: ALL PASS";
else echo "RECEIPT TEST BATTERY: FAILURES (see above)"; fi
exit "$rc"

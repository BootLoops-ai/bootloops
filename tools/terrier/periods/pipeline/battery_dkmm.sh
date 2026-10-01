#!/bin/bash
# battery_dkmm.sh — full DKMM regression battery: all run_pipe.py stages in
# order, 13 gate groups; every stage must print PASS (needs TERRIER_KKLT_BANK).
# Resource caps: single-core, ulimit -v 32505856, nice 5, stages sequential.
set -e
cd "$(dirname "$0")"
ulimit -v 32505856
PY="nice -n 5 python3 run_pipe.py"
for st in g2 g1 g3 stage1; do
  echo "=== $st ==="; $PY $st
done
echo "=== route R 150 ===";  $PY route R 150
echo "=== route R 60 ===";   $PY route R 60
echo "=== route C 150 ===";  $PY route C 150
echo "=== jets ===";         $PY jets
echo "=== vac 150 R ===";    $PY vac 150 R
echo "=== vac 60 R ===";     $PY vac 60 R
echo "=== vac 150 C ===";    $PY vac 150 C
echo "=== fdgate ===";       $PY fdgate
echo "=== gates ===";        $PY gates
echo "BATTERY COMPLETE: all stages green in periods/pipeline"

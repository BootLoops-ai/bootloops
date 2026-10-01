#!/usr/bin/env bash
# Build (and optionally run) the GiNaC cross-check oracle.
#   ./build.sh        -> compile only
#   ./build.sh run    -> compile and run (writes results.txt next to this script)
set -euo pipefail
cd "$(dirname "$0")"

CXX=${CXX:-g++}
CXXFLAGS="-O2 -std=c++17 $(pkg-config --cflags ginac 2>/dev/null || true)"
LIBS="$(pkg-config --libs ginac 2>/dev/null || echo '-lginac -lcln')"

$CXX $CXXFLAGS -o ginac_oracle ginac_oracle.cpp $LIBS
echo "built ./ginac_oracle"

if [[ "${1:-}" == "run" ]]; then
    ./ginac_oracle "$(pwd)/results.txt"
fi

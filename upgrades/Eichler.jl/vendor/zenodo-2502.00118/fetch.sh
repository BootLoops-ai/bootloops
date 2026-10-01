#!/usr/bin/env bash
# Fetch the Zenodo 10.5281/zenodo.14733100 ancillaries (arXiv:2502.00118, CC-BY-4.0).
#   ./fetch.sh        -> the small files the library consumes (~620 KB)
#   ./fetch.sh --all  -> the complete record (~211 MB; includes the two large
#                        helicity-coefficient files)
set -euo pipefail
cd "$(dirname "$0")"

BASE="https://zenodo.org/records/14733100/files"
SMALL=(DiffRelationsEllipticFunctions.m Roots.m CanonicalBasis.m CanonicalDiffEqs.m README.txt)
BIG=(HelCoeffs2lbareExpLM.m HelCoeffs2lbareExpTt.m)

files=("${SMALL[@]}")
[[ "${1:-}" == "--all" ]] && files+=("${BIG[@]}")

curl -fsSL "https://zenodo.org/api/records/14733100" -o record.json
for f in "${files[@]}"; do
    echo "fetching $f"
    curl -fL --retry 3 -o "$f" "$BASE/$f?download=1"
done
echo "done. License: CC-BY-4.0 (see record.json)."

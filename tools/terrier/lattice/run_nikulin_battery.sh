#!/bin/bash
# Nikulin battery: pins (3/3, gate in selftest_nikulin.jl) then validate03
# (sig(0,3) disc-form route vs primitive_embeddings count, 3 transversal classes;
# validate03 reads the transversal-lattice list named by TERRIER_TRANSVERSAL3,
# which is not included in the package — see the header of nikulin.jl).
set -e
cd "$(dirname "$0")"
: "${TERRIER_JULIA_PROJECT:?set TERRIER_JULIA_PROJECT to your Oscar/Hecke julia project}"
nice -n 5 julia +1.10 --project="$TERRIER_JULIA_PROJECT" -t 1 selftest_nikulin.jl
nice -n 5 julia +1.10 --project="$TERRIER_JULIA_PROJECT" -t 1 nikulin.jl validate03
echo "NIKULIN-BATTERY-DONE"

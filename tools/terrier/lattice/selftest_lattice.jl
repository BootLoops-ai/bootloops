# Lattice selftest wrapper — the battery IS validate.jl (18 gates; writes its report
# RESULTS_BUILD.md next to this file, exits 1 on any FAIL). Run: nice julia +1.10
# --project=<your Oscar env> -t 1 selftest_lattice.jl
include(joinpath(@__DIR__, "validate.jl"))

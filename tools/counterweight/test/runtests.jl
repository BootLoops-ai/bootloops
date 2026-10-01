#!/usr/bin/env julia
# runtests.jl — the registered battery entry (tools/BATTERIES.json "counterweight"):
# the test files in ONE Julia process, so startup + Nemo load + package JIT are paid
# once.  Each leg hard-fails (@assert / error -> nonzero exit) and prints its PASS
# banner.  The test files still run standalone — their package
# include is guarded so the second leg here reuses the already-loaded module.
#
# Run:  julia --project=. test/runtests.jl

include(joinpath(@__DIR__, "test_orbit_shear.jl"))
include(joinpath(@__DIR__, "validate_henn.jl"))
include(joinpath(@__DIR__, "test_sec53_qorbit.jl"))

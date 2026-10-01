# L1b selftest wrapper — nikulin.jl pins mode (3 pins) + receipt gate: RUN_pins.json
# pins.pass must be true (nikulin.jl main() itself always exits 0; the gate is here).
empty!(ARGS); push!(ARGS, "pins")
include(joinpath(@__DIR__, "nikulin.jl"))
import JSON
r = JSON.parsefile(joinpath(@__DIR__, "RUN_pins.json"))
ok = get(r["pins"], "pass", false) == true
println("selftest_nikulin: pins pass = ", ok)
exit(ok ? 0 : 1)

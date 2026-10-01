#!/usr/bin/env julia
# test_sec53_qorbit.jl -- a 2x2 block with an irreducible quadratic letter as a Counterweight control set (see test/fixtures/sec53_qorbit/README.md).
#   engine legs: b53_chart1_x -> ok = true; b53_raw_s / b53_landau_x -> STOP by name (the half-integer exponents at the quadratic
#   orbits are not removable by the implemented moves);  verifier legs (python3, the vendored point verifier): PASS on the fresh
#   output and on both pinned forms; FAIL rc 3 on the three planted controls and on the row-major reading of all three (12 legs).
# Run:  julia --project=. test/test_sec53_qorbit.jl        (rc 0 = every leg as expected; rc 1 = a leg differs, named)
isdefined(Main, :Counterweight) ||
    include(joinpath(@__DIR__, "..", "src", "Counterweight.jl"))
using .Counterweight
using .Counterweight.RatFunc, .Counterweight.Fuchsia, .Counterweight.Normalize
using .Counterweight.KiraDE: parse_kira_coeff
using Nemo
import JSON

const FIX = joinpath(@__DIR__, "fixtures", "sec53_qorbit")
const VERIFIER = joinpath(FIX, "verify_epsform_point.py")
const STOP_PREFIX = "no eps^0-reducing balance found and constant eps-decoupling failed"

function load_block(path)
    inp = JSON.parsefile(path); n = Int(inp["n"]); ctx = make_ctx()
    kin = Dict{Symbol,Any}(:x => ctx.Qex(ctx.x), :d => ctx.Qex(4 - 2*ctx.eps))
    A = zero_matrix(ctx.Qex, n, n)
    for i in 1:n, j in 1:n
        s = String(inp["A"][i][j]); s == "0" && continue
        A[i, j] = parse_kira_coeff(ctx, replace(s, "**" => "^"), kin)
    end
    return ctx, n, A, inp
end

function run_engine(name)
    ctx, n, A, inp = load_block(joinpath(FIX, name * ".json"))
    t = time(); r = try_epsfactor(ctx, A; max_balance_rounds=400); dt = time() - t
    stop = get(r.report, "stop", "")
    return r, n, inp, stop, dt
end

function write_out(path, r, n, inp)
    res = Dict{String,Any}("n" => n, "rows" => inp["rows"], "masters" => inp["masters"], "ok" => r.ok, "stop" => get(r.report, "stop", ""))
    if r.ok
        res["Atilde"] = [string(r.Atilde[i, j]) for i in 1:n, j in 1:n]
        res["T"] = [string(r.T[i, j]) for i in 1:n, j in 1:n]
    end
    open(path, "w") do io; JSON.print(io, res); end
end

function verifier_rc(block, out; notranspose=false)
    args = ["python3", VERIFIER, block, out]
    notranspose && push!(args, "--no-transpose")
    p = run(pipeline(`$args`, stdout=devnull, stderr=devnull); wait=false); wait(p)
    return p.exitcode
end

fails = String[]
function expect(cond, what)
    println((cond ? "PASS " : "FAIL ") * what); cond || push!(fails, what)
end

# --- engine legs ---
r1, n1, inp1, stop1, dt1 = run_engine("b53_chart1_x")
expect(r1.ok, "engine b53_chart1_x ok = true (t = $(round(dt1, digits=1)) s)")
outdir = mktempdir(prefix="sec53_qorbit_test_")
fresh = joinpath(outdir, "out_engine_fresh.json")
r1.ok && write_out(fresh, r1, n1, inp1)
for name in ("b53_raw_s", "b53_landau_x")
    r, _, _, stop, dt = run_engine(name)
    expect(!r.ok && startswith(stop, STOP_PREFIX), "engine $name STOP by name (t = $(round(dt, digits=1)) s): $(first(stop, 60))")
end
# --- verifier legs ---
blk = joinpath(FIX, "b53_chart1_x.json")
if r1.ok
    expect(verifier_rc(blk, fresh) == 0, "verifier PASS on the engine's fresh output")
    expect(verifier_rc(blk, fresh; notranspose=true) == 3, "verifier row-major CONTROL rc 3 on the fresh output")
end
expect(verifier_rc(blk, joinpath(FIX, "out_lee_b53_chart1_x.json")) == 0, "verifier PASS on the pinned Lee form")
expect(verifier_rc(blk, joinpath(FIX, "out_engine_b53_chart1_x.json")) == 0, "verifier PASS on the pinned engine form")
expect(verifier_rc(blk, joinpath(FIX, "out_lee_b53_chart1_x.json"); notranspose=true) == 3, "verifier row-major CONTROL rc 3 on the pinned Lee form")
expect(verifier_rc(blk, joinpath(FIX, "out_engine_b53_chart1_x.json"); notranspose=true) == 3, "verifier row-major CONTROL rc 3 on the pinned engine form")
for c in ("planted_drop_last_balance", "planted_mutated_Atilde", "planted_eps_in_Atilde")
    expect(verifier_rc(blk, joinpath(FIX, c * ".json")) == 3, "verifier FAIL rc 3 on $c")
end
println(isempty(fails) ? "SEC53_QORBIT: ALL LEGS AS EXPECTED" : "SEC53_QORBIT: $(length(fails)) leg(s) differ: " * join(fails, " | "))
exit(isempty(fails) ? 0 : 1)

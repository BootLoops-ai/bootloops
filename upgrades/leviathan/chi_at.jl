#!/usr/bin/env julia
# chi_at.jl — χ_regulated at SPECIFIC kinematic points (degenerate-slice preflight).
# Kept out of src/ deliberately; thin caller around Leviathan.chi_regulated.
# Usage:  julia +1.10 --project=<dir with Oscar+JSON, e.g. this one> chi_at.jl REQ.json
# REQ.json:
#   { "in_file": "path/to/family.in",
#     "points":  [ null, {"s12":"-3","s23":"-5",...}, ... ],   # null = generic random F_p
#     "seed":    20260630,                                     # optional
#     "mode":    "regulated" | "top" }                         # top = top-sector only
# Output: one JSON line per point on stdout:
#   {"idx":i,"chi":N,"mode":"regulated","wall_s":t,"p":prime,"generic":bool}

include(joinpath(@__DIR__, "src", "Leviathan.jl"))
using .Leviathan, Random, JSON

# Minimal .in reader (mirrors leviathan.jl CLI; duplicated to avoid touching it).
function _read_in(path)
    kv = Dict{String,String}()
    for line in eachline(path)
        line = strip(split(line, '#')[1]); isempty(line) && continue
        m = match(r"^(\w+)\s*=\s*(.*)$", line); m === nothing && error("bad line: $line")
        kv[lowercase(m[1])] = strip(m[2])
    end
    kin = [strip(t) for t in split(get(kv,"kin",""), ",") if !isempty(strip(t))]
    xs  = [strip(t) for t in split(kv["x"], ",") if !isempty(strip(t))]
    cut = haskey(kv,"cut") && !isempty(kv["cut"]) ?
          [parse(Int,strip(t)) for t in split(kv["cut"],",")] : Int[]
    G = haskey(kv,"g") ? kv["g"] : "(" * kv["u"] * ") + (" * kv["f"] * ")"
    family(G, String.(kin), String.(xs); cut=cut)
end

"Rational string \"p/q\" or \"p\" → residue in F_p (Int in [0,p-1])."
function _rat_to_fp(s::AbstractString, p::Int)
    s = strip(s)
    if occursin('/', s)
        num, den = split(s, '/'; limit=2)
        n, d = parse(BigInt, strip(num)), parse(BigInt, strip(den))
    else
        n, d = parse(BigInt, s), big(1)
    end
    mod(d, p) == 0 && error("denominator $d ≡ 0 mod p=$p")
    Int(mod(n * invmod(d, big(p)), big(p)))
end

function main()
    length(ARGS) >= 1 || (println(stderr, "usage: chi_at.jl REQ.json"); exit(1))
    req = JSON.parsefile(ARGS[1])
    fam = _read_in(req["in_file"])
    seed = get(req, "seed", 20260630)
    mode = get(req, "mode", "regulated")
    p = default_p31()
    rng = MersenneTwister(seed)
    Stop = collect(1:nedges(fam))
    for (i, pt) in enumerate(req["points"])
        if pt === nothing
            kv = random_kinvals(rng, fam.kinnames, p)
        else
            for n in fam.kinnames
                haskey(pt, n) || error("point $i missing kinematic '$n'")
            end
            kv = Dict{String,Int}(n => _rat_to_fp(string(pt[n]), p) for n in fam.kinnames)
        end
        t0 = time()
        chi = mode == "top" ?
              chi_regulated(fam, Stop; p=p, kinvals=kv, rng=MersenneTwister(seed+i)) :
              chi_regulated(fam;       p=p, kinvals=kv, rng=MersenneTwister(seed+i))
        dt = time() - t0
        println(JSON.json(Dict("idx"=>i-1, "chi"=>chi, "mode"=>mode,
                               "wall_s"=>round(dt,digits=3), "p"=>p,
                               "generic"=>(pt===nothing))))
        flush(stdout)
    end
end

main()

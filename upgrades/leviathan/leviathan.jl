#!/usr/bin/env julia
# leviathan CLI — chi-drop Landau singularities (arXiv:2606.29612 reimplementation).
# Run under:  julia +1.10 --project=<dir with Oscar+JSON, e.g. this one> leviathan.jl <cmd> <input>
# Commands:
#   chi <input>         per-sector chi (nu=0) + regulated chi + additivity
#   candidates <input>  candidate singularity factors (function-field elimination)
#   diagnose <input>    the four Appendix-A diagnostics
#   landau <input>      full pipeline: chi + candidates + Type-2.2 verification
# Input file (key = value, one per line; '#' comments):
#   G   = x1 + x2 + mm*(x1+x2)^2 - s*x1*x2     (or U = ..., F = ... on two lines)
#   kin = s, mm
#   x   = x1, x2
#   cut =                                       (optional, comma-separated edge indices)

include(joinpath(@__DIR__, "src", "Leviathan.jl"))
using .Leviathan, Random

function read_input(path)
    kv = Dict{String,String}()
    for line in eachline(path)
        line = strip(split(line, '#')[1])
        isempty(line) && continue
        m = match(r"^(\w+)\s*=\s*(.*)$", line)
        m === nothing && error("bad input line: $line")
        kv[lowercase(m[1])] = strip(m[2])
    end
    kin = [strip(t) for t in split(kv["kin"], ",") if !isempty(strip(t))]
    xs = [strip(t) for t in split(kv["x"], ",") if !isempty(strip(t))]
    cut = haskey(kv, "cut") && !isempty(kv["cut"]) ?
          [parse(Int, strip(t)) for t in split(kv["cut"], ",")] : Int[]
    G = haskey(kv, "g") ? kv["g"] : "(" * kv["u"] * ") + (" * kv["f"] * ")"
    family(G, String.(kin), String.(xs); cut=cut)
end

function main()
    length(ARGS) >= 2 || (println(read(@__FILE__, String) |> s -> split(s, "\n\n")[1]); exit(1))
    cmd, path = ARGS[1], ARGS[2]
    rng = MersenneTwister(length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 20260630)
    fam = read_input(path)
    if cmd == "chi"
        kv = random_kinvals(rng, fam.kinnames, default_p31())
        tot = 0; ok = true
        for S in admissible_sectors(fam)
            c = chi_sector(fam, S; kinvals=kv, rng=rng)
            println("chi_S ", S, " = ", c === nothing ? "Indeterminate" : c)
            c === nothing ? (ok = false) : (tot += c)
        end
        cr = chi_regulated(fam; kinvals=kv, rng=rng)
        println("chi_regulated = ", cr)
        println("sum_sectors   = ", ok ? tot : "Indeterminate",
                ok && tot == cr ? "  [additivity OK]" : "  [MISMATCH/INDET]")
    elseif cmd == "candidates"
        cands, bysec = all_candidates(fam; rng=rng)
        for (S, f) in sort(collect(bysec), by=p -> (length(p[1]), p[1]))
            println("sector ", S, ": ", sort(collect(f)))
        end
        println("ALL CANDIDATES: ", sort(collect(cands)))
    elseif cmd == "diagnose"
        r = landau(fam; rng=rng, verify=false)
        d = diagnose(fam; rng=rng, candidates=r.candidates)
        for (k, v) in sort(collect(d), by=first)
            println(k, " => ", v[1] ? "PASS" : "FAIL", "\n    evidence: ", v[2])
        end
    elseif cmd == "landau"
        lane = get(Dict(split(a, "=") for a in ARGS[3:end] if occursin("=", a)), "--lane", "exact")
        if lane == "reconstruct"
            np = parse(Int, get(ENV, "LEV_NPRIMES", "3"))
            nk = parse(Int, get(ENV, "LEV_NPOINTS", "60"))
            dc = parse(Int, get(ENV, "LEV_DEGCAP", "12"))
            print(report(landau_reconstruct(fam; nprimes=np, npoints=nk, deg_cap=dc,
                                            savedir=get(ENV, "LEV_SAVEDIR", nothing))))
        else
            print(report(landau(fam; rng=rng)))
        end
    else
        error("unknown command $cmd")
    end
end

main()

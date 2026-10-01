# Interface to the amflow-cpp engine (the AMFlow.cpp fork, which lives in the
# sibling repository amflow-cpp, cloned beside this one as ../amflow-cpp)
# via its JSON CLI, implementing the mandatory staggered-goal consistency
# protocol: never trust printed radii — a digit is accepted only if runs at
# goal G and G+30 (different eps_order) agree on it.

export AmflowFamily, amflow_run, amflow_sample, parse_ball, MASSIVE_FAMILY_OPTS

# Path to the amflow_cli binary (build it from
# the sibling repository amflow-cpp); resolved from PATH unless AMFLOW_CLI
# points at the binary directly.
const AMFLOW_CLI = get(ENV, "AMFLOW_CLI", "amflow_cli")

"amf_options required for massive multi-loop families:
top-level ending_schemes + blackbox-level ibp_rank/ibp_dot"
const MASSIVE_FAMILY_OPTS = (
    top = Dict("ending_schemes" => ["Tradition", "SingleMass"]),
    blackbox = Dict("ibp_rank" => 2, "ibp_dot" => 1))
const NO_EXTRA_OPTS = (top = Dict(), blackbox = Dict())

struct AmflowFamily
    name::String
    loops::Vector{String}
    legs::Vector{String}
    conservation::Dict{String,String}
    replacement::Dict{String,String}
    propagators::Vector{String}
end

"parse amflow ball string '[1.23 +/- 4.5e-30]' (or plain number) into Arb"
function parse_ball(s::AbstractString, prec::Int)
    s = strip(s)
    if startswith(s, "[")
        body = strip(s[2:end-1])
        parts = split(body, "+/-")
        midstr = strip(parts[1])
        mid = isempty(midstr) ? Arb(0, prec = prec) : Arb(midstr, prec = prec)
        rad = Arb(strip(parts[2]), prec = 64)
        Arblib.add_error!(mid, magof(rad))
        return mid
    end
    Arb(s, prec = prec)
end

"""
    amflow_run(fam, indices, numerics; goal, eps_order, d0=4, amf_extra, workdir)
        -> Dict{Int,Acb}

One CLI run; returns ε-order → coefficient (printed radii attached, NOT trusted).
"""
function amflow_run(fam::AmflowFamily, indices::Vector{Int},
                    numerics::Dict{String,String};
                    goal::Int, eps_order::Int, d0::Int = 4,
                    amf_extra = NO_EXTRA_OPTS, workdir::String, nthread::Int = 8,
                    timeout_s::Int = 14400,
                    ibp_cache_dir::Union{Nothing,String} = get(ENV, "AMFLOW_IBP_CACHE", nothing),
                    symbolic_ibp::Bool = get(ENV, "AMFLOW_SYMBOLIC_IBP", "0") == "1")
    mkpath(workdir)
    famdict = Dict("name" => fam.name, "loops" => fam.loops, "legs" => fam.legs,
                   "propagators" => fam.propagators, "replacement" => fam.replacement)
    isempty(fam.conservation) || (famdict["conservation"] = fam.conservation)
    bb = merge(Dict("numeric_values" => numerics, "n_thread" => nthread),
               amf_extra.blackbox)
    # IBP-reduction cache (Layer A: content-addressed Kira-results cache;
    # Layer B: symbolic-kinematics IBP so multiple oracle points share one
    # reduction).  Both default off; set AMFLOW_IBP_CACHE / AMFLOW_SYMBOLIC_IBP
    # env or pass kwargs.  Cuts multi-hour-per-point oracle farms to one
    # symbolic reduction plus minutes per point.
    ibp_cache_dir === nothing || (bb["ibp_cache_dir"] = ibp_cache_dir)
    symbolic_ibp && (bb["symbolic_ibp"] = true)
    cfg = Dict(
        "mode" => "solve_integrals",
        "options" => Dict("silent_mode" => true, "d0" => string(d0)),
        "family" => famdict,
        "integrals" => [Dict("indices" => indices)],
        "goal_digits" => goal,
        "eps_order" => eps_order,
        "work_dir" => joinpath(workdir, "work"),
        "amf_options" => merge(Dict("blackbox" => bb), amf_extra.top))
    cfgpath = joinpath(workdir, "config.json")
    outpath = joinpath(workdir, "out.json")
    open(cfgpath, "w") do io
        JSON3.write(io, cfg)
    end
    logpath = joinpath(workdir, "run.log")
    # wall clamp: timeout(1) is absent on stock macOS (coreutils installs it
    # as gtimeout); without either, run unclamped — loudly
    tbin = Sys.which("timeout")
    tbin === nothing && (tbin = Sys.which("gtimeout"))
    if tbin === nothing
        @warn "no timeout(1)/gtimeout(1) on PATH — amflow_cli runs without the $(timeout_s)s wall clamp (macOS: brew install coreutils)"
    end
    argv = tbin === nothing ? String[] : String[tbin, string(timeout_s)]
    append!(argv, String[AMFLOW_CLI, cfgpath, outpath])
    ok = success(pipeline(Cmd(argv); stdout = logpath, stderr = logpath))
    ok || error("amflow_cli failed; see $logpath")
    out = JSON3.read(read(outpath, String))
    prec = ceil(Int, goal * 3.33) + 64
    res = Dict{Int,Acb}()
    for entry in out["result"][1]["coefficients"]
        re = parse_ball(string(entry["value"]["re"]), prec)
        im = parse_ball(string(entry["value"]["im"]), prec)
        res[Int(entry["order"])] = Acb(re, im, prec = prec)
    end
    res
end

"""
    amflow_sample(fam, indices, numerics; goal, eps_order, …) -> Dict{Int,Acb}

Staggered-goal consistency protocol: run at (goal, eps_order) and
(goal+30, eps_order+2); for each ε order return a ball centred at the
higher-goal midpoint whose radius is the observed disagreement (×8 safety),
floored at 10^(−goal). Orders that disagree at O(1) raise an error.
"""
function amflow_sample(fam::AmflowFamily, indices::Vector{Int},
                       numerics::Dict{String,String};
                       goal::Int, eps_order::Int, d0::Int = 4,
                       amf_extra = NO_EXTRA_OPTS, workdir::String, nthread::Int = 8)
    r1 = amflow_run(fam, indices, numerics; goal = goal, eps_order = eps_order,
                    d0 = d0, amf_extra = amf_extra,
                    workdir = joinpath(workdir, "g$(goal)"), nthread = nthread)
    r2 = amflow_run(fam, indices, numerics; goal = goal + 30, eps_order = eps_order + 2,
                    d0 = d0, amf_extra = amf_extra,
                    workdir = joinpath(workdir, "g$(goal+30)"), nthread = nthread)
    out = Dict{Int,Acb}()
    for (ord, v1) in r1
        haskey(r2, ord) || continue
        v2 = r2[ord]
        diff = abs(v1 - v2)
        scale = max(abs(v1), Arb(1))
        Float64(diff / scale) < 1e-3 ||
            error("staggered-goal protocol FAILED at eps^$ord: |Δ| = $(Float64(diff))")
        rad = max(ub(diff) * 8, Arb(BigFloat(10)^(-goal)))
        v = Acb(v2)
        Arblib.add_error!(v, magof(rad))
        out[ord] = v
    end
    out
end

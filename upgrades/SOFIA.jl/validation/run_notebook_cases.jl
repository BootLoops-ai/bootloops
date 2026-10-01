# Batch validation against the worked examples recorded in SOFIA_examples.nb.
# Usage: julia --project=. validation/run_notebook_cases.jl [case indices...]
#
# Each case: build the Diagram from the notebook input, replicate the SOFIA
# call options (SolverBound, IncludeSubtopologies, IncludeISPs, MaxCut,
# LoopEdges), run sofia_singularities, and compare the polynomial set with
# the recorded output up to proportionality and variable naming.

using SOFIA
import Nemo

function read_cases(path)
    cases = Vector{Dict{String,Any}}()
    for line in eachline(path)
        f = split(line, "\t")
        length(f) >= 4 || continue
        push!(cases, Dict("pos" => f[1], "diag" => f[2],
                          "input" => f[3], "output" => f[4],
                          "section" => getidx(f, 5, ""), "rules" => getidx(f, 6, "")))
    end
    return cases
end


getidx(v::AbstractVector, i::Int, d) = i <= length(v) ? v[i] : d

"""Apply a chain of label substitution rules like
`//.M3:>M//.M__:>0//.mi_:>m` to a diagram's edge/node mass labels.
Wildcards m__/m_/mi_ match any indexed m-label (likewise M). Unparseable
trailing segments (formatting junk) end the chain."""
function apply_label_rules(d::Diagram, rules::AbstractString)
    isempty(rules) && return d
    labelname(x) = x isa SOFIA.WLCall ? string(x.head) * join(string.(x.args)) : string(x)
    edges = [(e.ends, e.mass) for e in d.edges]
    nodes = [(n.vertex, n.mass) for n in d.nodes]
    for seg in split(rules, r"//?\.")
        isempty(seg) && continue
        m = match(r"^([\w|_]+):>([\w]+)", seg)
        m === nothing && break
        lhs, rhs = m.captures
        pats = split(lhs, "|")
        matches(lbl) = begin
            lbl == BigInt(0) && return false
            nm = labelname(lbl)
            any(pats) do p
                if p in ("m__", "m_", "mi_")
                    occursin(r"^m\d+$", nm)
                elseif p in ("M__", "M_", "Mi_")
                    occursin(r"^M\d+$", nm)
                else
                    nm == p
                end
            end
        end
        newlbl = rhs == "0" ? BigInt(0) : Symbol(rhs)
        edges = [(uv, matches(ms) ? newlbl : ms) for (uv, ms) in edges]
        nodes = [(v, matches(ms) ? newlbl : ms) for (v, ms) in nodes]
    end
    return Diagram([SOFIA.Edge(uv[1], uv[2], ms) for (uv, ms) in edges],
                   [SOFIA.Node(v, ms) for (v, ms) in nodes])
end

"""Spaces in box-linearized output are implicit multiplication; spaces next
to operators/separators are formatting noise."""
function clean_output(s)
    outs = replace(s, r"\s+" => " ")
    # implicit multiplication: space between an atom/closer and an atom/opener
    outs = replace(outs, r"([\w)\]]) (?=[\w({])" => s"\1*")
    return replace(outs, " " => "")
end

"""Parse an upstream diagram literal {{{{u,v},mass},...},{{vtx,Mass},...}}."""
function parse_diagram(s::AbstractString)
    e = wlparse(s)
    e isa WLCall && e.head == :List && length(e.args) == 2 ||
        error("bad diagram literal")
    eds, nds = e.args
    edges = Tuple{Tuple{Int,Int},Any}[]
    for it in eds.args
        uv, mass = it.args
        push!(edges, ((Int(uv.args[1]), Int(uv.args[2])), mass))
    end
    nodes = Tuple{Int,Any}[]
    for it in nds.args
        push!(nodes, (Int(it.args[1]), it.args[2]))
    end
    return Diagram(edges, nodes)
end

masssym(x::Symbol) = x
function normname(nm::AbstractString)
    return nm
end

const ALL_DIAGS = let
    out = Tuple{Int,String}[]
    for line in eachline(joinpath(@__DIR__, "notebook_diagrams.txt"))
        f = split(line, "\t")
        length(f) >= 3 && push!(out, (parse(Int, f[2]), String(f[3])))
    end
    sort!(out)
    out
end

"""Names a diagram's labels can produce: squared-mass names plus the
kinematic basis for its node count."""
function diagram_namespace(d::Diagram)
    names = Set{String}()
    for e in d.edges
        e.mass == BigInt(0) || push!(names, SOFIA.sqmass_name(e.mass))
    end
    for n in d.nodes
        n.mass == BigInt(0) || push!(names, SOFIA.sqmass_name(n.mass))
    end
    nn = length(d.nodes)
    if nn == 2
        push!(names, "s")
    elseif nn >= 4
        for w in SOFIA.cyclic_basis_indices(nn)
            push!(names, SOFIA.sname(w))
        end
    end
    return names
end

function run_case(c; verbose=true)
    haskey(c, "error") && return (:skip, "extractor error")
    input = get(c, "input", "")
    occursin("Join[", input) && return (:skip, "composite Join case")
    occursin("Substitutions", input) && return (:skip, "numeric Substitutions case")
    d = parse_diagram(c["diag"])
    d = apply_label_rules(d, get(c, "rules", ""))
    verbose && println("  diagram: ", length(d.edges), " edges, ",
                       length(d.nodes), " nodes, L=", nloops(d))
    limit = occursin("SolverBound->Infinity", input) ? typemax(Int) :
            begin
                m = match(r"SolverBound->(\d+)", input)
                m === nothing ? 100 : parse(Int, m.captures[1])
            end
    subt = !occursin("IncludeSubtopologies->False", input)
    maxcut = !occursin("MaxCut->False", input)
    # LoopEdges->{...} pins the loop-edge choice
    le = nothing
    mle = match(r"LoopEdges->\{([\d,]+)\}", input)
    mle !== nothing && (le = parse.(Int, split(mle.captures[1], ",")))
    # IncludeISPs
    isps = occursin("IncludeISPs->True", input)

    routes = parse(Int, get(ENV, "ROUTES", "1"))
    t0 = time()
    local sing, R
    budget = parse(Float64, get(ENV, "ROUTE_BUDGET", "Inf"))
    opb = parse(Int, get(ENV, "OPBOUND", string(typemax(Int))))
    sing, R = sofia_singularities(d; include_subtopologies=subt, maxcut, limit,
                                  include_isps=isps, routes, loopedges=le,
                                  time_budget=budget, maxopterms=opb)
    dt = round(time() - t0; digits=1)

    outs = clean_output(c["output"])
    isempty(outs) && return (:skip, "empty output")
    startswith(outs, "{") || return (:skip, "output not a list")
    exp_l = wlparse(outs)
    exp_l isa WLCall && exp_l.head == :List || return (:skip, "output not a list")
    keys_ = sort!(collect(SOFIA.wlvariables(exp_l)); by=SOFIA.varname)
    Re, gens_ = Nemo.polynomial_ring(Nemo.QQ, SOFIA.varname.(keys_))
    vmap = Dict{Any,Any}(zip(keys_, gens_))
    expected = Any[]
    for a in exp_l.args
        try
            push!(expected, SOFIA.tonemo(a, vmap))
        catch
            return (:skip, "non-polynomial expected entry")
        end
    end

    # common ring: union of names
    mynames = map(string, Nemo.symbols(R))
    exnames = map(string, Nemo.symbols(Re))
    common = unique(vcat(mynames, exnames))
    Rc, _ = Nemo.polynomial_ring(Nemo.QQ, common)
    idx = Dict(nm => k for (k, nm) in enumerate(common))
    lift(p, names) = Nemo.evaluate(p, [Nemo.gen(Rc, idx[nm]) for nm in names])
    mine = [lift(p, mynames) for p in sing]
    theirs = [lift(p, exnames) for p in expected]

    key(p) = SOFIA.normal_form_key(p)
    mineset = Set(key.(mine))
    theirset = Set(key.(theirs))
    missing_ = [p for p in theirs if !(key(p) in mineset)]
    extra_ = [p for p in mine if !(key(p) in theirset)]
    status = isempty(missing_) ? (isempty(extra_) ? :exact : :superset) : :fail
    if verbose
        println("  time $(dt)s  mine=$(length(mineset)) expected=$(length(theirset)) missing=$(length(missing_)) extra=$(length(extra_))")
        for p in missing_[1:min(end, 6)]
            println("    MISSING: ", p)
        end
        for p in extra_[1:min(end, 6)]
            println("    EXTRA:   ", p)
        end
    end
    return (status, "")
end

cases = read_cases(joinpath(@__DIR__, "notebook_cases.tsv"))
sel = isempty(ARGS) ? collect(eachindex(cases)) : [parse(Int, a) + 1 for a in ARGS]
results = Dict{Symbol,Int}()
for i in sel
    c = cases[i]
    print("case $(i-1): ")
    desc = first(get(c, "input", get(c, "error", "?")), 70)
    println(desc)
    st, why = try
        run_case(c)
    catch err
        (:error, sprint(showerror, err)[1:min(end, 120)])
    end
    println("  -> ", st, isempty(why) ? "" : "  ($why)")
    flush(stdout)
    results[st] = get(results, st, 0) + 1
end
println("\nSummary: ", results)

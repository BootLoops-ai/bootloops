# Top-level singularity pipeline.
#
# Port of `SOFIASingularitiesSEED` / `SOFIASingularities` (scr.m:1231-1396).
# The SymmetryQuotient acceleration IS ported (symmetry.jl: symmetry_map /
# symmetry_quotient; sofia_singularities defaults symmetries=true and analyzes
# one representative per symmetry class; symmetries=false is the
# every-subtopology path). Still unported: the δ-degeneration trick
# (degenerate-limit letters) — results otherwise agree; see task list.

"""
    singularities_seed(d; maxcut=true, include_isps=false, limit=100)

Landau singularities of a single diagram (no subtopologies): build the LBL
Baikov Gram determinants and eliminate the Baikov variables with FastFubini.
Returns polynomials in the diagram's kinematic ring.
Port of `SOFIASingularitiesSEED` (scr.m:1232).
"""
function singularities_seed(d::Diagram; maxcut::Bool=true,
                            include_isps::Bool=false, limit::Integer=100,
                            loopedges=nothing, dotperm::Int=0, rowperm::Int=0,
                            deadline::Float64=Inf, maxopterms::Integer=typemax(Int))
    gramlist, active, sys = prepare_landau_system(d; maxcut, loopedges, dotperm,
                                                  rowperm)
    R = sys.ring
    input = sort(gramlist; by=nterms)
    if maxcut
        include_isps && (input = vcat([Nemo.gen(R, i) for i in active], input))
    else
        extra = include_isps ? [Nemo.gen(R, i) for i in active] :
                               [Nemo.gen(R, i) for i in 1:sys.nprop]
        input = vcat(extra, input)
    end
    isempty(active) && return (dedup_proportional([p for p in input if !is_constant(p)]), sys)
    sols = fastfubini(input, active; limit, deadline, maxopterms)
    return (sols, sys)
end

"""
    singularities_routes(d; routes=1, kwargs...) -> Vector{(sols, sys)}

Run `singularities_seed` over up to `routes` distinct elimination routes
(loop-edge choices x dot-pivot orders). The candidate singularity set on the
maximal cut is route-dependent — both here and upstream, where it depends on
`FixLoopEdges` and Mathematica's `Solve` pivoting (the SOFIA notebook itself
shows automatic vs pinned LoopEdges giving different sets for the degenerate
acnode). The union over routes is a robust superset of any single route.
"""
function singularities_routes(d::Diagram; routes::Integer=1, maxcut::Bool=true,
                              include_isps::Bool=false, limit::Integer=100,
                              loopedges=nothing, time_budget::Real=Inf,
                              maxopterms::Integer=typemax(Int))
    combos = Tuple{Any,Int,Int}[]
    les = loopedges !== nothing ? [loopedges] : loop_edge_candidates(d)
    # loop-edge ORDER is itself a route choice (it fixes which momentum is
    # l_1 vs l_2, distinguishing edges in the candidate set)
    les = unique(vcat(les, [reverse(le) for le in les if length(le) > 1]))
    for rp in 0:2, dp in 0:3, le in les
        push!(combos, (le, dp, rp))
    end
    length(combos) > routes && resize!(combos, routes)
    out = Vector{Tuple{Any,Any}}()
    t0 = time()
    err = nothing
    for (le, dp, rp) in combos
        # a single bad route can be catastrophically slower than the rest;
        # stop adding routes once the budget is spent, and abandon a route
        # mid-elimination via the deadline (>=1 route always runs, but even
        # the first is bounded by 4x budget so a hang cannot block forever)
        !isempty(out) && time() - t0 > time_budget && break
        dl = isfinite(time_budget) ? t0 + (isempty(out) ? 4 : 1) * time_budget : Inf
        try
            push!(out, singularities_seed(d; maxcut, include_isps, limit,
                                          loopedges=le, dotperm=dp, rowperm=rp,
                                          deadline=dl, maxopterms))
        catch e
            err = e
            continue
        end
    end
    if isempty(out)
        if err isa RouteTimeout
            # over budget on every route: SKIP this diagram, like upstream's
            # UnclogTime mechanism (the skip is the caller's to report)
            return out
        end
        push!(out, singularities_seed(d; maxcut, include_isps, limit, maxopterms))
    end
    return out
end

"""
    sofia_singularities(d; include_subtopologies=true, maxcut=true,
                        limit=100) -> Vector

Candidate Landau singularities of `d` including all its subtopologies,
factored and deduplicated, expressed in a master kinematic ring
(generators: cyclic Mandelstam basis, external MM's, every internal mass
appearing in `d`). Port of `SOFIASingularities` with Symmetries->False.
"""
function sofia_singularities(d::Diagram; include_subtopologies::Bool=true,
                             maxcut::Bool=true, limit::Integer=100,
                             include_isps::Bool=false, routes::Integer=1,
                             loopedges=nothing, time_budget::Real=Inf,
                             symmetries::Bool=true, maxopterms::Integer=typemax(Int),
                             progress::Union{Nothing,Function}=nothing)
    tops = include_subtopologies ? subtopologies(d) : [d]
    # the symmetry quotient (upstream default Symmetries->True): compute each
    # equivalence class of subtopologies once and transport the result to the
    # members through the kinematic maps — semantically equivalent to the
    # direct per-subtopology run (docs/PORTING_NOTES.md), much cheaper on
    # symmetric diagrams
    local reps, classes
    if symmetries && length(tops) > 1
        reps, classes = symmetry_quotient(tops)
    else
        reps = tops
        classes = [Vector{Dict{String,Any}}() for _ in tops]
    end
    entries = Any[]            # (sols, seed ring, maps-for-this-class)
    for (i, t) in enumerate(reps)
        progress !== nothing && progress(i, length(reps))
        for (sols, sys) in singularities_routes(t; routes, maxcut, limit,
                                                include_isps, time_budget, maxopterms,
                                                loopedges=(t === d ? loopedges : nothing))
            push!(entries, (sols, sys.ring, classes[i]))
        end
    end
    # master ring: union of all generator names (x's excluded), including
    # names appearing in symmetry-map images (members may carry masses the
    # representative does not)
    names = String[]
    addname(nm) = (startswith(nm, "x") && occursin(r"^x\d+$", nm)) ||
                  (nm in names) || push!(names, nm)
    for (_, R, maps) in entries
        foreach(addname, map(string, Nemo.symbols(R)))
        for mp in maps, v in values(mp)
            foreach(addname, map(string, Nemo.symbols(Nemo.parent(v))))
        end
    end
    Rm, _ = Nemo.polynomial_ring(Nemo.QQ, names)
    gidx = Dict(nm => k for (k, nm) in enumerate(names))
    bynm(nm) = haskey(gidx, nm) ? Nemo.gen(Rm, gidx[nm]) : zero(Rm)
    out = Any[]
    for (sols, R, maps) in entries
        isempty(sols) && continue
        rnames = map(string, Nemo.symbols(R))
        identity_images = Any[bynm(nm) for nm in rnames]
        usemaps = isempty(maps) ? [nothing] : maps
        for mp in usemaps
            images = if mp === nothing || isempty(mp)
                identity_images
            else
                Rmap = Nemo.parent(first(values(mp)))
                mapimages = Any[bynm(string(s)) for s in Nemo.symbols(Rmap)]
                Any[haskey(mp, nm) ? Nemo.evaluate(mp[nm], mapimages) : bynm(nm)
                    for nm in rnames]
            end
            for p in sols
                push!(out, Nemo.evaluate(p, images))
            end
        end
    end
    final = [p for p in factor_list_unique(out) if !is_constant(p) && !iszero(p)]
    return final, Rm
end

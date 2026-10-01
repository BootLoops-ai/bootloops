# The FastFubini elimination engine — SOFIA's default Landau-singularity solver.
#
# Faithful port of scr.m lines 450-513 (FastFubini, Fubini, RunningFubini).
# Algorithm (arXiv:2503.16601 eq. 2.30): eliminate Baikov variables from the
# system of Gram determinants by iterated discriminants/resultants/leading
# coefficients, factoring at every step, and — the "Fast" part — taking
# intersections over all single-variable elimination orders at every level
# of the subset lattice of variables.
#
# Variables are addressed by generator index into the polynomials' common ring.

"""Thrown when a route exceeds its time budget mid-elimination."""
struct RouteTimeout <: Exception end

checkdeadline(deadline) = time() > deadline && throw(RouteTimeout())

"""
    fubini(polys, vars; limit=100, factorlast=true) -> Vector

Eliminate the variables `vars` (generator indices, in order) from `polys`.
For each variable: collect leading coefficients (which is the polynomial
itself when free of the variable), discriminants of every polynomial, and
pairwise resultants; then factor and dedup. Port of `Fubini` (scr.m:465).

`limit` skips discriminants/resultants of polynomials with more terms
(upstream `SolverBound`). `factorlast=false` reproduces upstream's option of
not factoring after the final elimination.
"""
function fubini(polys::AbstractVector, vars::AbstractVector{Int};
                limit::Integer=100, factorlast::Bool=true,
                deadline::Float64=Inf, maxopterms::Integer=typemax(Int))
    # `maxopterms` guards single resultant/discriminant calls whose input
    # term counts predict an intractable operation (FLINT native and PRS both
    # stall beyond a few thousand terms at high degree); skipped ops shrink
    # the candidate set exactly like upstream's SolverBound skip does.
    opok(f) = nterms(f) <= maxopterms
    Snew = collect(polys)
    for (k, var) in enumerate(vars)
        term = eltype(Snew)[]
        # term1: leading coefficients + discriminants
        for f in Snew
            push!(term, leading_coeff_in(f, var))
            if !free_of(f, var) && nterms(f) <= limit && opok(f)
                checkdeadline(deadline)
                push!(term, discriminant_in(f, var))
            end
        end
        # term2: pairwise resultants
        for i in eachindex(Snew), j in (i+1):lastindex(Snew)
            (free_of(Snew[i], var) || free_of(Snew[j], var)) && continue
            (nterms(Snew[i]) > limit || nterms(Snew[j]) > limit) && continue
            (opok(Snew[i]) && opok(Snew[j])) || continue
            checkdeadline(deadline)
            push!(term, resultant_in(Snew[i], Snew[j], var))
        end
        if !factorlast && k == length(vars)
            Snew = [f for f in term if !is_constant(f) && !iszero(f)]
        else
            Snew = eltype(term)[]
            for f in term
                append!(Snew, irreducible_factors(f))
            end
            Snew = [f for f in dedup_proportional(Snew) if !is_constant(f)]
        end
    end
    return Snew
end

"""
    intersect_proportional(lists) -> Vector

Elements of the first list having a proportional match in every other list
(upstream `Intersection[..., SameTest->ProportionalPolynomialsQ]`).
"""
function intersect_proportional(lists::AbstractVector)
    isempty(lists) && return Any[]
    length(lists) == 1 && return dedup_proportional(lists[1])
    keysets = [Set(normal_form_key(p) for p in l) for l in lists[2:end]]
    out = eltype(lists[1])[]
    for p in dedup_proportional(lists[1])
        k = normal_form_key(p)
        all(ks -> k in ks, keysets) && push!(out, p)
    end
    return out
end

"""
    fastfubini(polys, vars; limit=100, factorlast=true) -> Vector

The FastFubini solver (scr.m:450): walk the subset lattice of `vars` by size;
for each subset, eliminate its variables one at a time starting from the
already-computed result for each co-subset, and intersect the outcomes over
the choice of last-eliminated variable. Polynomials free of all remaining
target variables are set aside ("freePolynomials") and rejoined at the end.

`vars` are generator indices in the common parent ring of `polys`.
"""
function fastfubini(polys::AbstractVector, vars::AbstractVector{Int};
                    limit::Integer=100, factorlast::Bool=true,
                    deadline::Float64=Inf, maxopterms::Integer=typemax(Int),
                    progress::Union{Nothing,Function}=nothing)
    set = Dict{Vector{Int},Vector{eltype(polys)}}()
    set[Int[]] = collect(polys)
    free = eltype(polys)[]
    varset = sort(collect(vars))
    for size in 1:length(varset)
        for sub in sorted_subsets(varset, size)
            tables = Vector{eltype(polys)}[]
            for v in sub
                co = sort!(setdiff(sub, [v]))
                push!(tables, fubini(set[co], [v]; limit, factorlast, deadline,
                                     maxopterms))
            end
            pre = intersect_proportional(tables)
            progress !== nothing && progress(sub, length(pre))
            kept = eltype(polys)[]
            for f in pre
                if all(vi -> free_of(f, vi), varset)
                    push!(free, f)
                else
                    push!(kept, f)
                end
            end
            set[sub] = kept
        end
    end
    return dedup_proportional(vcat(set[varset], free))
end

"""All `size`-element subsets of `v` (sorted vectors, in lexicographic order
over positions — same enumeration as upstream `Subsets[v, {size}]`)."""
function sorted_subsets(v::Vector{Int}, size::Int)
    out = Vector{Int}[]
    n = length(v)
    idx = collect(1:size)
    size == 0 && return [Int[]]
    size > n && return out
    while true
        push!(out, v[idx])
        # advance combination
        i = size
        while i >= 1 && idx[i] == n - size + i
            i -= 1
        end
        i == 0 && break
        idx[i] += 1
        for j in (i+1):size
            idx[j] = idx[j-1] + 1
        end
    end
    return out
end

"""
    runningfubini(polys, vars; limit=8, factorlast=true, orders=nothing)

Intersection of `fastfubini` results over random permutations of `vars`
(upstream `RunningFubini`, scr.m:502). `orders` defaults to all permutations
in random order; pass a vector of orders for determinism.
"""
function runningfubini(polys::AbstractVector, vars::AbstractVector{Int};
                       limit::Integer=8, factorlast::Bool=true,
                       orders=nothing)
    perms = orders === nothing ? collect(permutations(collect(vars))) : orders
    running = nothing
    for perm in perms
        cur = fastfubini(polys, collect(perm); limit, factorlast)
        running = running === nothing ? cur : intersect_proportional([running, cur])
    end
    return running === nothing ? eltype(polys)[] : running
end

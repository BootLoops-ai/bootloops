# Symmetry maps between (sub)topologies and the symmetry quotient.
#
# Port of `Symmetry` (scr.m:1107) and `SymmetryQuotient0` / `SymmetryQuotient`
# (scr.m:1224-1226). Upstream finds graph isomorphisms with
# `FindGraphIsomorphism` on the simple internal graphs (scr.m:1123), filters
# them by parallel-edge multiplicities and mass-label availability
# (scr.m:1130-1149), and then recovers the kinematic replacement rules by
# solving a linear system equating the cyclic-basis invariants of the two
# diagrams built from vertex-total external momenta (scr.m:1162-1206).
#
# This port constructs the rules directly from the isomorphism data instead:
# a vertex map sigma, an external-leg (node-index) bijection pi, and a global
# mass-label map f.  The rule for a Mandelstam basis invariant s_W of A is
# the expansion of (sum_{j in W} p_{pi(j)})^2 in B's kinematics, and the rule
# for a squared-mass generator nm is the generator named f(nm) — equivalent
# to upstream's linear solve, but exact by construction.  Differences kept on
# purpose:
#   * `symmetry_map(A, B)` always maps A -> B; upstream's `Symmetry` may
#     return the rules in either direction, preferring the more generic
#     diagram (scr.m:1189-1210).
#   * the FIRST valid (pi, f) found in a deterministic search order is
#     returned; upstream picks a "best" solution by `SortBy[..., Length]`
#     (scr.m:1209).  Any valid map is correct, they differ only in shape.
#   * upstream memoizes `Symmetry`; the search here is cheap enough
#     (diagrams have <= ~10 internal vertices) that we do not.

export symmetry_map, symmetry_quotient

# ---------------------------------------------------------------- helpers

"""True for the massless label (`BigInt(0)` by convention; any numeric zero
is accepted)."""
_masszero(m) = m isa Number && iszero(m)

"""Parallel-edge bundles of the internal graph: `ends => [mass labels...]`
(the simple graph is the key set)."""
function _bundles(d::Diagram)
    b = Dict{Tuple{Int,Int},Vector{Any}}()
    for e in d.edges
        push!(get!(b, e.ends, Any[]), e.mass)
    end
    return b
end

"""All vertex labels of `d`: internal-edge endpoints plus external-leg
attachment points (vertex labels can be sparse after contractions)."""
_allvertices(d::Diagram) =
    sort!(unique!(vcat(Int[v for e in d.edges for v in e.ends],
                       Int[nd.vertex for nd in d.nodes])))

"""Distinct nonzero squared-mass generator names of `d` (edges first, then
nodes, in order of appearance).  Edge and node labels share one namespace
via `sqmass_name`, exactly as in `lbl` (src/baikov.jl)."""
function _massnames(d::Diagram)
    out = String[]
    for e in d.edges
        _masszero(e.mass) && continue
        nm = sqmass_name(e.mass)
        nm in out || push!(out, nm)
    end
    for nd in d.nodes
        _masszero(nd.mass) && continue
        nm = sqmass_name(nd.mass)
        nm in out || push!(out, nm)
    end
    return out
end

"""Isomorphism-invariant local signature of vertex `v`: sorted multiplicities
of the incident parallel-edge bundles, self-loop bundle multiplicity, number
of attached external legs, and how many of those legs are massless."""
function _vertex_signature(d::Diagram, bnd, v::Int)
    sizes = Int[]
    selfm = 0
    for (k, labels) in bnd
        if k == (v, v)
            selfm = length(labels)
        elseif k[1] == v || k[2] == v
            push!(sizes, length(labels))
        end
    end
    nleg = 0
    nleg0 = 0
    for nd in d.nodes
        nd.vertex == v || continue
        nleg += 1
        _masszero(nd.mass) && (nleg0 += 1)
    end
    return (sort!(sizes), selfm, nleg, nleg0)
end

# ------------------------------------------------------- rule construction

"""
    _build_rules(A, B, nmap, f) -> Dict{String,Any}

Kinematic replacement rules induced by the leg bijection `nmap` (node index
of A -> node index of B) and the mass-label name map `f`.  Keys are the
generator names of A's master kinematic ring (cyclic Mandelstam basis for
n = #legs, plus A's squared-mass names); values are polynomials in a fresh
ring `QQ[basis..., massnames(A) ∪ massnames(B)...]`.  Substituting these
rules (matching generators by name) into a singularity polynomial of A
yields the corresponding polynomial of B.  This replaces upstream's linear
solve between the two cyclic bases (scr.m:1170-1206).
"""
function _build_rules(A::Diagram, B::Diagram, nmap::Dict{Int,Int},
                      f::Dict{String,String})
    n = length(A.nodes)
    mA = _massnames(A)
    mB = _massnames(B)
    kin = n >= 2 ? generate_kinematics(n) : nothing
    names = kin === nothing ? String[] : copy(kin.basis)
    for nm in vcat(mA, mB)
        nm in names || push!(names, nm)
    end
    rules = Dict{String,Any}()
    isempty(names) && return rules      # fully massless few-point: nothing to map
    Rmap, g = Nemo.polynomial_ring(Nemo.QQ, names)
    gd = Dict(nm => g[k] for (k, nm) in enumerate(names))
    # mass rules: nm -> f(nm)
    for nm in mA
        rules[nm] = gd[f[nm]]
    end
    kin === nothing && return rules
    # Mandelstam rules: s_W -> (sum_{j in W} p_{nmap(j)})^2 expanded with B's
    # kinematics.  kin.ring generators are [basis..., MM_1..MM_n]; the generic
    # MM_j is sent to B's node-j labeled mass generator (zero if massless),
    # exactly like `lbl`'s kin_images (src/baikov.jl).
    images = Any[gd[nm] for nm in kin.basis]
    for j in 1:n
        m = B.nodes[j].mass
        push!(images, _masszero(m) ? zero(Rmap) : gd[sqmass_name(m)])
    end
    if n == 2
        # generate_kinematics(2) has the single invariant s = p1^2 outside the
        # cyclic-window scheme
        rules["s"] = Nemo.evaluate(kin.dots[(nmap[1], nmap[1])], images)
    else
        for w in cyclic_basis_indices(n)
            pw = [nmap[j] for j in w]
            poly = zero(kin.ring)
            for a in pw, b in pw
                poly += kin.dots[(min(a, b), max(a, b))]
            end
            rules[sname(w)] = Nemo.evaluate(poly, images)
        end
    end
    return rules
end

# ------------------------------------------- label / leg matching per sigma

"""
    _complete_map(A, B, sigma, bndA, bndB) -> Dict{String,Any} or nothing

Given a simple-graph isomorphism `sigma` (vertex of A -> vertex of B), try to
extend it to a full symmetry: choose a globally consistent mass-label map
`f` over the parallel-edge bundles (upstream's per-bundle `Permutations`,
scr.m:1134-1149) and a per-vertex bijection of external legs (upstream's
`vtomomenta`/`vtoM`, scr.m:1162-1168), then build the kinematic rules.
Backtracks over all bundle pairings and leg matchings; `nothing` if none is
consistent.  Consistency means: f is a well-defined function on label names
and massless maps to massless (upstream's `EnoughScalesQ`, scr.m:1142-1148).
"""
function _complete_map(A::Diagram, B::Diagram, sigma::Dict{Int,Int},
                       bndA, bndB)
    f = Dict{String,String}()      # A label name -> B label name
    nmap = Dict{Int,Int}()         # A node index -> B node index
    keysA = sort!(collect(keys(bndA)))
    nvertsA = sort!(unique!(Int[nd.vertex for nd in A.nodes]))

    # --- stage 2: per-vertex bijection on external legs ------------------
    function legs_rec(i)
        i > length(nvertsA) && return _build_rules(A, B, nmap, f)
        v = nvertsA[i]
        w = sigma[v]
        idxA = Int[j for (j, nd) in enumerate(A.nodes) if nd.vertex == v]
        idxB = Int[j for (j, nd) in enumerate(B.nodes) if nd.vertex == w]
        length(idxA) == length(idxB) || return nothing
        zA = Int[j for j in idxA if _masszero(A.nodes[j].mass)]
        zB = Int[j for j in idxB if _masszero(B.nodes[j].mass)]
        length(zA) == length(zB) || return nothing
        cA = Dict{String,Vector{Int}}()   # A name -> node indices at v
        for j in idxA
            _masszero(A.nodes[j].mass) && continue
            push!(get!(cA, sqmass_name(A.nodes[j].mass), Int[]), j)
        end
        rem = Dict{String,Vector{Int}}()  # B name -> still-free indices at w
        for j in idxB
            _masszero(B.nodes[j].mass) && continue
            push!(get!(rem, sqmass_name(B.nodes[j].mass), Int[]), j)
        end
        anames = sort!(collect(keys(cA)))
        function assign(ai)
            if ai > length(anames)
                # massless legs pair up in index order (any bijection works)
                for (a, b) in zip(zA, zB)
                    nmap[a] = b
                end
                r = legs_rec(i + 1)
                r !== nothing && return r
                for a in zA
                    delete!(nmap, a)
                end
                return nothing
            end
            a = anames[ai]
            need = cA[a]
            targets = haskey(f, a) ? [f[a]] : sort!(collect(keys(rem)))
            for b in targets
                avail = get(rem, b, nothing)
                avail === nothing && continue
                length(avail) >= length(need) || continue
                rem[b] = avail[length(need)+1:end]
                hadf = haskey(f, a)
                f[a] = b
                for (x, y) in zip(need, avail)
                    nmap[x] = y
                end
                r = assign(ai + 1)
                r !== nothing && return r
                for x in need
                    delete!(nmap, x)
                end
                hadf || delete!(f, a)
                rem[b] = avail
            end
            return nothing
        end
        return assign(1)
    end

    # --- stage 1: mass-label map over parallel-edge bundles --------------
    function bundles_rec(i)
        i > length(keysA) && return legs_rec(1)
        kA = keysA[i]
        kB = (min(sigma[kA[1]], sigma[kA[2]]), max(sigma[kA[1]], sigma[kA[2]]))
        LB = get(bndB, kB, nothing)
        LB === nothing && return nothing
        LA = bndA[kA]
        length(LA) == length(LB) || return nothing
        count(_masszero, LA) == count(_masszero, LB) || return nothing
        cA = Dict{String,Int}()
        for m in LA
            _masszero(m) && continue
            nm = sqmass_name(m)
            cA[nm] = get(cA, nm, 0) + 1
        end
        rem = Dict{String,Int}()
        for m in LB
            _masszero(m) && continue
            nm = sqmass_name(m)
            rem[nm] = get(rem, nm, 0) + 1
        end
        anames = sort!(collect(keys(cA)))
        function assign(ai)
            ai > length(anames) && return bundles_rec(i + 1)
            a = anames[ai]
            ca = cA[a]
            targets = haskey(f, a) ? [f[a]] : sort!(collect(keys(rem)))
            for b in targets
                get(rem, b, 0) >= ca || continue
                rem[b] -= ca
                hadf = haskey(f, a)
                f[a] = b
                r = assign(ai + 1)
                r !== nothing && return r
                hadf || delete!(f, a)
                rem[b] += ca
            end
            return nothing
        end
        return assign(1)
    end

    return bundles_rec(1)
end

# ----------------------------------------------------------------- public

"""
    symmetry_map(A::Diagram, B::Diagram) -> Dict{String,Any} or nothing

Port of `Symmetry` (scr.m:1107).  `A` and `B` must be subtopologies of the
same parent diagram (same number of external legs; masses drawn from the
same label set).  Returns `nothing` if they are not isomorphic as
mass-labeled multigraphs-with-legs; otherwise a `Dict` mapping generator
NAMES of A's master kinematic ring (Mandelstam basis names "s12", "s123",
... from `cyclic_basis_indices(n)`, and squared-mass names "MM1"/"mm1" via
`sqmass_name`) to polynomials over `QQ[basis..., masses of A and B...]`,
such that substituting the rules into a singularity polynomial of A yields
the corresponding polynomial of B.  An identity/trivial map (exact duplicate
or reflection) is still returned as a Dict — only non-isomorphism gives
`nothing`.

Search order is deterministic (sorted vertices/labels/indices) and the first
valid map is returned; upstream instead ranks candidate solutions by length
(scr.m:1209).
"""
function symmetry_map(A::Diagram, B::Diagram)
    # quick rejects (upstream scr.m:1109 edge count, scr.m:1115 multi-edge
    # multiplicity distribution; the rest are free additional invariants)
    nedges(A) == nedges(B) || return nothing
    length(A.nodes) == length(B.nodes) || return nothing
    nloops(A) == nloops(B) || return nothing
    vsA = _allvertices(A)
    vsB = _allvertices(B)
    length(vsA) == length(vsB) || return nothing
    bndA = _bundles(A)
    bndB = _bundles(B)
    sort!(Int[length(l) for l in values(bndA)]) ==
        sort!(Int[length(l) for l in values(bndB)]) || return nothing
    count(e -> _masszero(e.mass), A.edges) ==
        count(e -> _masszero(e.mass), B.edges) || return nothing
    count(nd -> _masszero(nd.mass), A.nodes) ==
        count(nd -> _masszero(nd.mass), B.nodes) || return nothing
    # f maps A's label set onto B's, so A needs at least as many scales
    length(_massnames(A)) >= length(_massnames(B)) || return nothing

    sigA = Dict(v => _vertex_signature(A, bndA, v) for v in vsA)
    sigB = Dict(w => _vertex_signature(B, bndB, w) for w in vsB)
    multA = Dict(k => length(l) for (k, l) in bndA)
    multB = Dict(k => length(l) for (k, l) in bndB)

    # plain backtracking enumeration of simple-graph isomorphisms with
    # degree/signature pruning (replaces upstream's FindGraphIsomorphism,
    # scr.m:1123; internal graphs have <= ~10 vertices)
    sigma = Dict{Int,Int}()
    used = Set{Int}()
    nv = length(vsA)
    function extend(k)
        k > nv && return _complete_map(A, B, sigma, bndA, bndB)
        v = vsA[k]
        for w in vsB
            w in used && continue
            sigA[v] == sigB[w] || continue
            ok = true
            for (u, su) in sigma
                ka = (min(u, v), max(u, v))
                kb = (min(su, w), max(su, w))
                if get(multA, ka, 0) != get(multB, kb, 0)
                    ok = false
                    break
                end
            end
            ok || continue
            sigma[v] = w
            push!(used, w)
            r = extend(k + 1)
            r !== nothing && return r
            delete!(sigma, v)
            delete!(used, w)
        end
        return nothing
    end
    return extend(1)
end

"""
    symmetry_quotient(tops::Vector{Diagram})
        -> (reps::Vector{Diagram}, classes::Vector{Vector{Dict{String,Any}}})

Port of `SymmetryQuotient0` / `SymmetryQuotient` (scr.m:1224-1226), with the
upstream `BestDiagram` fixed-point representative selection
(scr.m:1212-1223) simplified to a greedy first-seen partition: for each
diagram try `symmetry_map(rep, d)` against the existing class
representatives in order; on the first hit the map is recorded in that
class, otherwise `d` starts a new class.  `classes[i]` lists the maps from
`reps[i]` to every member of its class, INCLUDING the identity map for the
representative itself.
"""
function symmetry_quotient(tops::AbstractVector{Diagram})
    reps = Diagram[]
    classes = Vector{Vector{Dict{String,Any}}}()
    for d in tops
        placed = false
        for (i, r) in enumerate(reps)
            m = symmetry_map(r, d)
            if m !== nothing
                push!(classes[i], m)
                placed = true
                break
            end
        end
        if !placed
            m = symmetry_map(d, d)
            m === nothing &&
                error("symmetry_quotient: identity map not found (unreachable)")
            push!(reps, d)
            push!(classes, Dict{String,Any}[m])
        end
    end
    return reps, classes
end

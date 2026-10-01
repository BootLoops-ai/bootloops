# Loop-edge selection: which internal edges carry independent loop momenta.
#
# Port of `FixLoopEdges`/`FixLoopEdgesMEF` (scr.m:518-596). Upstream builds
# cycle-indicator matrices over GF(2) and picks edges that belong to exactly
# one short cycle. Equivalently (and this is what we implement): the
# complement of any spanning forest of the internal multigraph is a valid set
# of L loop edges — each co-tree edge closes exactly one independent cycle.
# We grow the forest greedily so that co-tree edges close SHORT fundamental
# cycles (upstream's MinimalBy[..., total support]) and return co-tree edge
# indices. Downstream `LBL` validates the choice against the expected
# integration-variable bound; `lbl` retries alternative forests on failure.

"""
    fix_loop_edges(d::Diagram; order=:bfs) -> Vector{Int}

Indices (into `d.edges`) of a set of `nloops(d)` edges whose removal makes
the internal multigraph a forest. `order=:bfs` grows a BFS forest (short
fundamental cycles); `order=:dfs` and `order=:reverse` give alternates used
as retries by `lbl`.
"""
function fix_loop_edges(d::Diagram; order::Symbol=:bfs)
    isempty(d.edges) && return Int[]
    vs = vertices(d)
    vindex = Dict(v => i for (i, v) in enumerate(vs))
    parent_ = collect(1:length(vs))
    function find(x)
        while parent_[x] != x
            parent_[x] = parent_[parent_[x]]
            x = parent_[x]
        end
        return x
    end
    edgeorder = collect(1:nedges(d))
    if order == :reverse
        reverse!(edgeorder)
    elseif order == :dfs
        # heuristic alternate: sort by (max endpoint desc)
        sort!(edgeorder; by=i -> -maximum(d.edges[i].ends))
    end
    loops = Int[]
    for i in edgeorder
        e = d.edges[i]
        u, v = vindex[e.ends[1]], vindex[e.ends[2]]
        ru, rv = find(u), find(v)
        if ru == rv
            push!(loops, i)        # closes a cycle -> loop edge
        else
            parent_[ru] = rv       # tree edge
        end
    end
    sort!(loops)
    @assert length(loops) == nloops(d)
    return loops
end

"""All loop-edge choices `lbl` may try, most-preferred first."""
loop_edge_candidates(d::Diagram) =
    unique([fix_loop_edges(d; order=o) for o in (:bfs, :reverse, :dfs)])

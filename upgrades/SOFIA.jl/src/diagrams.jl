# Feynman diagram representation and subtopology generation.
#
# Port of scr.m lines 208-314 (ContractEdge, Contractions, DeleteTadpoles,
# OneVertexIrreducibleQ, NotContactQ, Subtopologies) with the upstream's
# implicit edge-node lists replaced by an explicit typed Diagram.
#
# Upstream represents a diagram as {edges, nodes} with
#   edges = {{{u,v}, mass}, ...}   (internal lines, undirected)
#   nodes = {{vertex, Mass}, ...}  (external legs attached at vertex)
# Mass labels are symbolic atoms (0, m, m[1], M[2], ...) — we keep them as
# parser ASTs (Any) compared by structural equality.

struct Edge
    ends::Tuple{Int,Int}   # stored sorted, u <= v
    mass::Any
end
Edge(u::Int, v::Int, mass) = Edge((min(u, v), max(u, v)), mass)
Base.:(==)(a::Edge, b::Edge) = a.ends == b.ends && a.mass == b.mass
Base.hash(e::Edge, h::UInt) = hash(e.mass, hash(e.ends, hash(:Edge, h)))

struct Node
    vertex::Int
    mass::Any
end
Base.:(==)(a::Node, b::Node) = a.vertex == b.vertex && a.mass == b.mass
Base.hash(n::Node, h::UInt) = hash(n.mass, hash(n.vertex, hash(:Node, h)))

struct Diagram
    edges::Vector{Edge}
    nodes::Vector{Node}
end
Base.:(==)(a::Diagram, b::Diagram) = a.edges == b.edges && a.nodes == b.nodes
Base.hash(d::Diagram, h::UInt) = hash(d.nodes, hash(d.edges, hash(:Diagram, h)))

"""Construct a Diagram from upstream-style nested lists:
`edges = [((u,v), mass), ...]`, `nodes = [(vertex, Mass), ...]`."""
function Diagram(edges::AbstractVector{<:Tuple}, nodes::AbstractVector{<:Tuple})
    Diagram([Edge(e[1][1], e[1][2], e[2]) for e in edges],
            [Node(n[1], n[2]) for n in nodes])
end

nedges(d::Diagram) = length(d.edges)
vertices(d::Diagram) = sort!(unique!([v for e in d.edges for v in e.ends]))

"""Number of loops L = 1 + E - V over internal edges (scr.m:278 `NLoops`)."""
function nloops(d::Diagram)
    isempty(d.edges) && return 0
    return 1 + nedges(d) - length(vertices(d))
end

"""Contract internal edge at index `i`: identify its endpoints, keeping the
smaller vertex label (scr.m:208 `ContractEdge`)."""
function contract_edge(d::Diagram, i::Int)
    e = d.edges[i]
    u, v = e.ends                       # u <= v: v is replaced by u
    relabel(x) = x == v ? u : x
    newedges = Edge[]
    for (j, f) in enumerate(d.edges)
        j == i && continue
        push!(newedges, Edge(relabel(f.ends[1]), relabel(f.ends[2]), f.mass))
    end
    newnodes = [Node(relabel(n.vertex), n.mass) for n in d.nodes]
    return Diagram(newedges, newnodes)
end

"""All single-edge contractions (scr.m:218 `ContractSingle`)."""
contract_single(d::Diagram) = [contract_edge(d, i) for i in 1:nedges(d)]

"""Remove self-loop (tadpole) edges (scr.m:294 `DeleteTadpoles`)."""
delete_tadpoles(d::Diagram) =
    Diagram([e for e in d.edges if e.ends[1] != e.ends[2]], d.nodes)

"""True unless the diagram has no internal edges (scr.m:293 `NotContactQ`)."""
not_contact(d::Diagram) = !isempty(d.edges)

"""
One-vertex-irreducibility test (scr.m:284 `OneVertexIrreducibleQ`):
no single internal vertex disconnects the (simple) internal graph.
Upstream logic: if a minimum vertex cut has size > 1 → irreducible; else if
more than 2 vertices → reducible; for <= 2 vertices, irreducible iff the
external legs do not all attach at a single vertex.
"""
function one_vertex_irreducible(d::Diagram)
    isempty(d.edges) && return false
    vs = vertices(d)
    g = Graphs.SimpleGraph(maximum(vs))
    for e in d.edges
        e.ends[1] != e.ends[2] && Graphs.add_edge!(g, e.ends[1], e.ends[2])
    end
    sub = g[vs]   # induced subgraph on the actual vertices
    if length(vs) > 2
        Graphs.is_connected(sub) || return false
        return isempty(Graphs.articulation(sub))
    else
        # 1- or 2-vertex diagram (banana-like): irreducible iff external legs
        # attach at more than one vertex
        return length(unique([n.vertex for n in d.nodes])) > 1
    end
end

"""Canonical form used by upstream `Subtopologies` for dedup: stable-sort the
edge list (upstream `SortBy[edges, Length]` is stable and a no-op on the key,
so order of edges is preserved; we sort by (ends, mass-print) for a genuine
canonical order while remaining deterministic)."""
function canonical(d::Diagram)
    es = sort(d.edges; by=e -> (e.ends, repr(e.mass)))
    return Diagram(es, d.nodes)
end

"""
    subtopologies(d) -> Vector{Diagram}

All non-trivial subtopologies of `d`: iterated single-edge contractions,
with tadpoles deleted, keeping only one-vertex-irreducible non-contact
diagrams, deduplicated (scr.m:219 `Contractions` + scr.m:309
`Subtopologies`). Includes `d` itself first.
"""
function subtopologies(d::Diagram)
    level = [d]
    out = Diagram[]
    seen = Set{Diagram}()
    while !isempty(level)
        for x in level
            c = canonical(x)
            if !(c in seen)
                push!(seen, c)
                push!(out, x)
            end
        end
        next = Diagram[]
        for x in level
            for y in contract_single(x)
                y2 = delete_tadpoles(y)
                if one_vertex_irreducible(y2) && not_contact(y2)
                    c = canonical(y2)
                    c in seen && continue
                    push!(next, y2)
                end
            end
        end
        # dedup within the level
        lvlseen = Set{Diagram}()
        level = Diagram[]
        for y in next
            c = canonical(y)
            if !(c in lvlseen)
                push!(lvlseen, c)
                push!(level, y)
            end
        end
    end
    return out
end

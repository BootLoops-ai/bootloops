"""
SOFIA.jl — Singularities of Feynman Integrals Automatized, in Julia.

A Wolfram-free port of the Mathematica package SOFIA by Miguel Correia,
Mathieu Giroux, and Sebastian Mizera (arXiv:2503.16601, CPC 320 (2026) 109970).
All algorithmic credit belongs to the SOFIA authors; this package contributes
only the Julia implementation, built on Nemo/FLINT and Graphs.jl.

Pinned upstream sources used as the porting reference are vendored under
`reference/` (commit recorded in reference/UPSTREAM_COMMIT.txt).
"""
module SOFIA

import Nemo
import Graphs
using Combinatorics: permutations

export wlparse, wlparsefile, loadsingularities, WLCall
export isproportional, dedup_proportional, factor_list_unique, irreducible_factors
export fubini, fastfubini, runningfubini, intersect_proportional
export Diagram, Edge, Node, nloops, subtopologies, contract_edge,
       delete_tadpoles, one_vertex_irreducible
export Kinematics, generate_kinematics, fix_loop_edges
export lbl, BaikovSystem, prepare_landau_system, singularities_seed,
       sofia_singularities
export OddLetter, find_odd_letters, effortless_odd_letters, independent_letters
export PLDSystem, pld_delta, pld_system, pld_discriminants, pld_singularities

include("wlparser.jl")
include("polyutils.jl")
include("fastfubini.jl")
include("diagrams.jl")
include("kinematics.jl")
include("loopedges.jl")
include("baikov.jl")
include("singularities.jl")
include("symmetry.jl")
include("effortless.jl")
include("pld.jl")

end # module

"""cosmoflow -- FRW/dS wavefunction and correlator integrands and their
letter alphabets for an ARBITRARY site graph (chains, rings, stars,
multi-loop polygons): connected-subgraph facet forms q_g, folded and disc
letter classes in ansatzer schema (`alphabet_graph`, CLI kind `graph`), and
the Cayley--Menger Baikov polynomial of the 1-loop n-gon (`baikov_B`).
The triangle (1-loop 3-site) numerics -- oracle, maxcut periods -- are the
worked example under `cosmoflow.examples.triangle`.  See README.md / GUIDE.md.
"""
from .polytope import baikov_B, q_subgraph, baikov_B3, q_G12, B_z, B_num
from .alphabet import (alphabet_graph, alphabet_tree_chain,
                       alphabet_loop_ngon, connected_subgraphs,
                       q_of_subgraph, disc_kallen_block, parse_edges,
                       loop_ngon_edges, tree_chain_edges)

__all__ = [
    'alphabet_graph', 'alphabet_tree_chain', 'alphabet_loop_ngon',
    'connected_subgraphs', 'q_of_subgraph', 'disc_kallen_block',
    'parse_edges', 'loop_ngon_edges', 'tree_chain_edges',
    'baikov_B', 'q_subgraph', 'baikov_B3', 'q_G12', 'B_z', 'B_num',
]

"""Cell selection for a PIECE: first-admissible-pair group (hf_piece.group_cells) optionally refined to one denominator-SUPPORT class
(cells whose atom-name set equals `support_atoms`, powers ignored) -- the sub-split used for g_gggg, whose group-level common
denominators make the single-rational-function assembly (InstanceRF) too large (several GB)."""
import hf_piece
def support_classes(cellres, pair):
    sub = hf_piece.group_cells(cellres, pair)
    classes = sorted(set(frozenset(nm for nm, pw in dk) for dk, conv in sub["cells"]), key=lambda s: (len(s), sorted(s)))
    return classes
def select_cells(cellres, pair, support_atoms=None):
    sub = hf_piece.group_cells(cellres, pair)
    if support_atoms is None: return sub
    S = frozenset(support_atoms)
    out = dict(sub); out["cells"] = [(dk, conv) for dk, conv in sub["cells"] if frozenset(nm for nm, pw in dk) == S]
    return out

"""Shared plumbing: the engine-facing Partition shim and helpers."""
import numpy as np
from flint import arb, arb_mat, ctx


class Part:
    """Minimal partition object satisfying clinch.engine's contract:
    blocks (list of flat-index arrays), block_names, border_idx."""

    def __init__(self, blocks, block_names, border_idx):
        self.blocks = [np.asarray(b, dtype=int) for b in blocks]
        self.block_names = list(block_names)
        self.border_idx = np.asarray(border_idx, dtype=int)


def cur_prec():
    """The ambient flint precision (the engine's ctx_guard owns it)."""
    return ctx.prec


def neg(x):
    return -x


def mat1(v):
    return arb_mat([[v]])

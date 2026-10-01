"""clinch.mask — MaskedOracle: fix a set of coordinates (the ACTIVE SET of
a corner candidate) and present the REDUCED block-arrow system.

The A2 corner semantics: a published
optimum sitting on a parameter bound must certify AS a KKT/corner optimum.
The reduced system = the free coordinates' stationarity with the active
coordinates pinned at their bounds (exact arb constants); the KKT sign leg
is the ENGINE's business (full-oracle gradient over the certified box).

The verdict register follows the AT-CORNER audit convention (an AT-CORNER
flag alone leaves KKT adjudication out of scope; CLINCH adjudicates)."""
import numpy as np
from flint import arb

__all__ = ["MaskedOracle"]


class MaskedOracle:
    """Wrap a block-arrow oracle with a set of FIXED coordinates.

    fixed: dict {flat_theta_index: exact value (float)}. Fixed coordinates
    are removed from their blocks/border; empty blocks are dropped. The
    wrapped oracle's F/H are evaluated on the FULL system with the fixed
    coordinates as exact arb constants, then rows/columns are deleted —
    the reduced system is exactly the free-coordinate stationarity."""

    def __init__(self, base, fixed):
        self.base = base
        self.fixed = {int(k): float(v) for k, v in fixed.items()}
        bp = base.part
        self._base_blocks = [np.asarray(idx) for idx in bp.blocks]
        self._base_border = np.asarray(bp.border_idx)
        # free structure
        self.blocks_map = []      # per kept block: (base_block_i, keep_pos)
        blocks, names = [], []
        base_names = list(bp.block_names)
        for bi, idx in enumerate(self._base_blocks):
            keep = [j for j, t in enumerate(idx) if int(t) not in self.fixed]
            if not keep:
                continue
            blocks.append(idx[keep])
            names.append(base_names[bi])
            self.blocks_map.append((bi, keep))
        self.border_keep = [j for j, t in enumerate(self._base_border)
                            if int(t) not in self.fixed]
        self.border_idx = self._base_border[self.border_keep]
        if len(self.border_idx) == 0:
            raise ValueError("mask removed the whole border — unsupported")
        self._blocks = blocks
        self._names = names

    # engine-facing partition shim
    @property
    def part(self):
        m = self

        class _P:
            blocks = m._blocks
            block_names = m._names
            border_idx = m.border_idx
        return _P()

    @property
    def dims(self):
        return ([len(b) for b in self._blocks], len(self.border_idx))

    def _full_x(self, x):
        """Reduced (zs, g) -> full (zs, g) with fixed coords as exact arb."""
        zs, g = x
        keep_of = {bi: keep for (bi, keep), _ in zip(self.blocks_map, zs)}
        red_of = {bi: blk for (bi, _), blk in zip(self.blocks_map, zs)}
        full_z = []
        for bi, idx in enumerate(self._base_blocks):
            keep_map = {}
            if bi in keep_of:
                keep_map = {kp: red_of[bi][i]
                            for i, kp in enumerate(keep_of[bi])}
            full_z.append([keep_map[j] if j in keep_map
                           else arb(self.fixed[int(t)])
                           for j, t in enumerate(idx)])
        gm = {kp: g[i] for i, kp in enumerate(self.border_keep)}
        full_g = [gm[j] if j in gm else arb(self.fixed[int(t)])
                  for j, t in enumerate(self._base_border)]
        return full_z, full_g

    def F(self, x):
        Fz, Fg = self.base.F(self._full_x(x))
        out_z = [[Fz[bi][j] for j in keep] for bi, keep in self.blocks_map]
        out_g = [Fg[j] for j in self.border_keep]
        return out_z, out_g

    def F_full(self, x):
        """Full-system gradient at a reduced point (the KKT sign leg)."""
        return self.base.F(self._full_x(x))

    def H(self, x):
        from flint import arb_mat
        D, B, G = self.base.H(self._full_x(x))
        Dr, Br = [], []
        for bi, keep in self.blocks_map:
            n = len(keep)
            Dn = arb_mat(n, n)
            for a in range(n):
                for b in range(n):
                    Dn[a, b] = D[bi][keep[a], keep[b]]
            Bn = arb_mat(n, len(self.border_keep))
            for a in range(n):
                for b, gj in enumerate(self.border_keep):
                    Bn[a, b] = B[bi][keep[a], gj]
            Dr.append(Dn)
            Br.append(Bn)
        ngr = len(self.border_keep)
        Gr = arb_mat(ngr, ngr)
        for a in range(ngr):
            for b in range(ngr):
                Gr[a, b] = G[self.border_keep[a], self.border_keep[b]]
        return Dr, Br, Gr

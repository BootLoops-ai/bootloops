"""CondOracle — ball-rigorous oracle for the R4A variant-B fit objective
of record (step4_sealed_eval.fit_cond, frozen text):

  f(b) = sum_d w_d [ -sum_c y_dc z_dc + Y_d LSE_c(z_dc) ]
         + ridge ||b||^2,      z_dc = o_dc + Xs[d,c]·b,
  w_d = 1/max(Y_d,1) (weight='equal'),  o = log(max(lam_ref,1e-9)).

Days with Y_d = 0 contribute EXACTLY zero (y row zero and w_d Y_d = 0)
and are skipped identically in every implementation. The objective is
smooth everywhere (no rails; ZCLIP is emission-only). Derivatives are the
standard multinomial forms, assembled with arb_mat contractions:

  grad = gconst + sum_d wY_d (Xd^T s_d) + 2 r b
  H    = sum_d wY_d [ Xd^T diag(s_d) Xd - (Xd^T s_d)(Xd^T s_d)^T ] + 2 r I

with s_d = softmax(z_d) computed in balls via the exact-constant shift
M_d (any constant shift is exact): s_c = exp(z_c - M - log S), S =
sum_c exp(z_c - M). All outputs are rigorous enclosures over the box.
"""
import numpy as np
from flint import arb, arb_mat, ctx

from clinch.oracle_v31 import OracleRefusal

PREC_INIT = 192

__all__ = ["CondOracle", "Part", "OracleRefusal"]


class Part:
    def __init__(self, blocks, block_names, border_idx):
        self.blocks = [np.asarray(b, dtype=int) for b in blocks]
        self.block_names = list(block_names)
        self.border_idx = np.asarray(border_idx, dtype=int)


class CondOracle:
    def __init__(self, data, part, tag=""):
        old_prec = ctx.prec
        ctx.prec = PREC_INIT
        try:
            self._init(data, part, tag)
        finally:
            ctx.prec = old_prec

    def _init(self, data, part, tag):
        self.part = part
        self.tag = tag
        Xs = np.asarray(data["Xs"], float)      # (nd, k, p)
        o = np.asarray(data["o"], float)        # (nd, k)
        y = np.asarray(data["y"], float)        # (nd, k)
        self.ridge = arb(float(data["ridge"]))
        Yd = y.sum(1)
        act = np.where(Yd > 0)[0]
        self.act = act
        self.nd, self.k, self.p = Xs.shape
        self.n_active = len(act)
        wd = 1.0 / np.maximum(Yd, 1.0)
        self.wY = [arb(float(wd[d] * Yd[d])) for d in act]
        # cached arb structures per active day
        self._Xd = [arb_mat([[arb(float(v)) for v in Xs[d, c]]
                             for c in range(self.k)]) for d in act]
        self._Xd_rows = [[[arb(float(v)) for v in Xs[d, c]]
                          for c in range(self.k)] for d in act]
        self._o = [[arb(float(v)) for v in o[d]] for d in act]
        # exact constant parts: C0 = -sum_d wd sum_c y*o (active days),
        # gconst = -sum_d wd sum_c y*X
        C0 = arb(0)
        gconst = [arb(0) for _ in range(self.p)]
        for i, d in enumerate(act):
            wda = arb(float(wd[d]))
            for c in range(self.k):
                yv = y[d, c]
                if yv == 0.0:
                    continue
                wy = wda * arb(float(yv))
                C0 -= wy * self._o[i][c]
                row = self._Xd_rows[i][c]
                for a in range(self.p):
                    gconst[a] -= wy * row[a]
        self._C0 = C0
        self._gconst = gconst

    @property
    def dims(self):
        return ([len(b) for b in self.part.blocks],
                len(self.part.border_idx))

    def _b_from_box(self, x):
        zs, g = x
        th = [None] * self.p
        for bi, blk in enumerate(self.part.blocks):
            for pos, t in enumerate(blk):
                th[int(t)] = zs[bi][pos]
        for pos, t in enumerate(self.part.border_idx):
            th[int(t)] = g[pos]
        if any(v is None for v in th):
            raise OracleRefusal("partition", "part does not cover b")
        return th

    def _day_z(self, i, bcol):
        zc = self._Xd[i] * bcol
        return [zc[c, 0] + self._o[i][c] for c in range(self.k)]

    def _softmax(self, z):
        """(s list, lse) rigorous; shift by exact midpoint constant."""
        M = arb(max(float(v.mid()) for v in z))
        ex = [(v - M).exp() for v in z]
        S = arb(0)
        for e in ex:
            S += e
        if not (S > 0):
            raise OracleRefusal("softmax", "sum exp not provably > 0")
        lse = M + S.log()
        iS = 1 / S
        s = [e * iS for e in ex]
        return s, lse

    def _val(self, th, order, want_H=False):
        bcol = arb_mat([[v] for v in th])
        val = self._C0
        for a in range(self.p):
            val += self._gconst[a] * th[a]
        grad = ([self._gconst[a] for a in range(self.p)]
                if order >= 1 else None)
        H = ([[arb(0)] * self.p for _ in range(self.p)]
             if want_H else None)
        for i in range(self.n_active):
            z = self._day_z(i, bcol)
            s, lse = self._softmax(z)
            wY = self.wY[i]
            val += wY * lse
            if order >= 1:
                srow = arb_mat([s])          # 1 x k
                xs = srow * self._Xd[i]      # 1 x p
                for a in range(self.p):
                    grad[a] += wY * xs[0, a]
            if want_H:
                sc = arb_mat([[s[c] * v for v in self._Xd_rows[i][c]]
                              for c in range(self.k)])
                A = self._Xd[i].transpose() * sc          # p x p
                for a in range(self.p):
                    xa = xs[0, a]
                    for b in range(self.p):
                        H[a][b] += wY * (A[a, b] - xa * xs[0, b])
        # ridge
        rid = arb(0)
        for a in range(self.p):
            rid += th[a] * th[a]
        val += self.ridge * rid
        if order >= 1:
            for a in range(self.p):
                grad[a] += 2 * self.ridge * th[a]
        if want_H:
            for a in range(self.p):
                H[a][a] += 2 * self.ridge
            for a in range(self.p):
                for b in range(self.p):
                    if not H[a][b].is_finite():
                        raise OracleRefusal("nonfinite_oracle",
                                            f"H ({a},{b}) not finite")
        if not val.is_finite():
            raise OracleRefusal("nonfinite_oracle", "value not finite")
        return val, grad, H

    # ---------------- engine-facing contract -------------------------
    def nll(self, x):
        th = self._b_from_box(x)
        v, _, _ = self._val(th, 0)
        return v

    def _grad_full(self, x):
        th = self._b_from_box(x)
        v, g, _ = self._val(th, 1)
        return g, v

    def F(self, x):
        grad, _ = self._grad_full(x)
        Fz = [[grad[int(t)] for t in blk] for blk in self.part.blocks]
        Fg = [grad[int(t)] for t in self.part.border_idx]
        return Fz, Fg

    def F_full(self, x):
        return self.F(x)

    def H(self, x):
        th = self._b_from_box(x)
        _, _, Hfull = self._val(th, 1, want_H=True)
        bi_ = self.part.border_idx
        D = []
        B = []
        for blk in self.part.blocks:
            D.append(arb_mat([[Hfull[int(s)][int(t)] for t in blk]
                              for s in blk]))
            B.append(arb_mat([[Hfull[int(s)][int(t)] for t in bi_]
                              for s in blk]))
        G = arb_mat([[Hfull[int(s)][int(t)] for t in bi_] for s in bi_])
        return D, B, G

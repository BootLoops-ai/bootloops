"""CodedOracle — ball-rigorous oracle for the R4B fit objective of record
(fit_multipliers.coded_nll_factory):

    f(x) = (1/n) sum_r phi_r(eta_r) + l2 * ||x[1:]||^2,
    eta_r = x0 + Z[r]·x[1:],   lam_r = exp(o_r + eta_r),
    phi_r = lam                      (X_r = 0)
          = lam - (o+eta)            (X_r = 1)
          = lam + log 2 - 2(o+eta)   (X_r = 2)
          = lam - log S3(lam)        (X_r = 3), S3 = e^l - 1 - l - l^2/2

X = min(y, 3) (capped coded symbol). REFERENCE RAILS (certified target =
the analytic objective; the rails are censused per box, fail-closed):
  clip |eta| < 4  — every row's eta ball must satisfy |eta| < 4 provably;
  p3 floor 1e-9   — every X=3 row's p3 = e^{-lam} S3 must be > 1e-9
                    provably.
A box where either rail cannot be excluded raises OracleRefusal (named).

S3 wide-ball safety: for small lam the direct e^l-1-l-l^2/2 cancels; the
oracle uses the series sum_{k>=3} l^k/k! with a certified geometric tail
bound when lam is provably < 8, else the direct form (cancellation there
is benign at 192 bits). Branch choice by PROVABLE bound, never midpoint.

Derivatives: per-row second-order jets in eta (clinch.jets, N=1, by
identity); contractions as arb_mat products (C speed). Gram assembly:
  grad = A^T u / n + 2 l2 (0, beta),  H = A^T diag(w) A / n + 2 l2 D,
A = [1 | Z]. All entries arb balls — rigorous enclosures over the box.
"""
import numpy as np
from flint import arb, arb_mat, ctx

from clinch.jets import Jet, jvar
from clinch.oracle_v31 import OracleRefusal

__all__ = ["CodedOracle", "Part", "OracleRefusal"]

CLIP_ETA = 4.0
P3_FLOOR = 1e-9
LOG2 = None


class Part:
    def __init__(self, blocks, block_names, border_idx):
        self.blocks = [np.asarray(b, dtype=int) for b in blocks]
        self.block_names = list(block_names)
        self.border_idx = np.asarray(border_idx, dtype=int)


def _fact_inv(kmax):
    out = [arb(1)]
    f = arb(1)
    for k in range(1, kmax + 1):
        f *= k
        out.append(1 / f)
    return out


class CodedOracle:
    """data: dict(y int array, logref float array, Zg (n,q) float array,
    l2 float). part: Part over x indices (p = q+1, x0 = intercept)."""

    def __init__(self, data, part, tag=""):
        old_prec = ctx.prec
        ctx.prec = 192
        try:
            self._init(data, part, tag)
        finally:
            ctx.prec = old_prec

    def _init(self, data, part, tag):
        global LOG2
        self.part = part
        self.tag = tag
        y = np.asarray(data["y"])
        self.X = np.minimum(y, 3).astype(np.int8)
        self.o = [arb(float(v)) for v in data["logref"]]
        Z = np.asarray(data["Zg"], float)
        self.n, self.q = Z.shape
        self.p = self.q + 1
        self.l2 = arb(float(data["l2"]))
        # A = [1 | Z] cached as arb rows (floats exact in arb)
        self._A = arb_mat([[arb(1)] + [arb(float(v)) for v in Z[r]]
                           for r in range(self.n)])
        self._Arows = None   # lazily built python rows for scaled mats
        self._Zf = Z
        LOG2 = arb(2).log()
        self._finv = _fact_inv(26)

    @property
    def dims(self):
        return ([len(b) for b in self.part.blocks],
                len(self.part.border_idx))

    def _x_from_box(self, x):
        zs, g = x
        th = [None] * self.p
        for bi, blk in enumerate(self.part.blocks):
            for pos, t in enumerate(blk):
                th[int(t)] = zs[bi][pos]
        for pos, t in enumerate(self.part.border_idx):
            th[int(t)] = g[pos]
        if any(v is None for v in th):
            raise OracleRefusal("partition", "part does not cover x")
        return th

    def _etas(self, th):
        col = arb_mat([[v] for v in th])
        e = self._A * col
        return [e[r, 0] for r in range(self.n)]

    def _S123(self, lv):
        """(S1, S2, S3) BALLS at the lam ball lv: S_m = sum_{k>=m}
        lv^k/k!. lv provably < 1 -> positive-term series with certified
        tail pads (no cancellation: widths stay ~relwidth(lv)); else the
        direct exp forms (cancellation <= ~1 digit for lv >= 1).
        Design note: a direct-form JET for S3
        passed the center gates but blew the Krawczyk margins at box
        radii (jet-component subtraction turned relative widths into
        1/lam^2-amplified ones, H-enclosure width ~1e5 x Lipschitz);
        the hand-ratio assembly below is the cure."""
        if lv < 1:
            t = lv
            S1 = lv
            S2 = arb(0)
            S3 = arb(0)
            k = 1
            while k < 140:
                k += 1
                t = t * lv / k
                S1 += t
                S2 += t
                if k >= 3:
                    S3 += t
                u = lv.abs_upper() / (k + 1)
                if k >= 4 and u < 0.5:
                    tail = t.v.abs_upper() if hasattr(t, "v") else \
                        t.abs_upper()
                    tail = tail * u / (1 - u)
                    if tail < arb(2) ** (-200):
                        pad = arb(0, tail)
                        return S1 + pad, S2 + pad, S3 + pad
            raise OracleRefusal("s123_tail", "series tail not certified")
        e = lv.exp()
        S1 = e - 1
        S2 = S1 - lv
        S3 = S2 - lv * lv * arb(0.5)
        return S1, S2, S3

    def _row_jets(self, etas, order2):
        """Per-row (phi value sum, u=dphi/deta list, w=d2phi/deta2 list).
        Also runs both rail censuses fail-closed."""
        n = self.n
        u = [None] * n
        w = [None] * n if order2 else None
        val = arb(0)
        four = arb(CLIP_ETA)
        floor = arb(P3_FLOOR)
        for r in range(n):
            ev = etas[r]
            if not ((ev < four) and (ev > -four)):
                raise OracleRefusal(
                    "clip_rail", f"|eta| < 4 not provable at row {r}")
            ej = jvar(ev, 0, 1, order2)
            lam = (ej + self.o[r]).exp()
            X = int(self.X[r])
            if X == 0:
                phi = lam
            elif X == 1:
                phi = lam - (ej + self.o[r])
            elif X == 2:
                phi = lam + LOG2 - (ej + self.o[r]) * 2
            else:
                # X = 3: phi = lam - log S3, hand-assembled derivatives
                # (exact identities):
                #   phi'  = lam (1 - g),        g  = S2/S3
                #   phi'' = lam (1 - g) - g' lam^2,
                #   g' = (S1 S3 - S2^2) / S3^2
                lv = lam.v
                if lv < 3:
                    # small-lam regime: positive-term series ratios (the
                    # e^l - poly subtraction cancels catastrophically
                    # below l ~ 1; ratio cancellation here is <= ~5x)
                    S1, S2, S3 = self._S123(lv)
                    if not (S3 > 0):
                        raise OracleRefusal(
                            "s3_pos", f"S3 not provably > 0 at row {r}")
                    p3v = (-lv).exp() * S3
                    if not (p3v > floor):
                        raise OracleRefusal(
                            "floor_rail",
                            f"p3 > 1e-9 not provable at row {r}")
                    iS3 = 1 / S3
                    gg = S2 * iS3
                    one_m_g = 1 - gg
                    phi_v = lv - S3.log()
                    phi_g = lv * one_m_g
                    if order2:
                        gp = (S1 * S3 - S2 * S2) * iS3 * iS3
                        phi_h = phi_g - gp * lv * lv
                        phi = Jet(phi_v, [phi_g], [phi_h], 1)
                    else:
                        phi = Jet(phi_v, [phi_g], None, 1)
                else:
                    # large-lam regime: direct jet chain (S3 ~ e^l
                    # dominates — subtraction benign, ratios balanced;
                    # the RATIO form here cancels e^{2l}-vs-l^2 e^l/2 —
                    # the measured w-rad 0.08 rows, JOURNAL ledger)
                    S3j = lam.exp() - 1 - lam - lam * lam * arb(0.5)
                    if not (S3j.v > 0):
                        raise OracleRefusal(
                            "s3_pos", f"S3 not provably > 0 at row {r}")
                    p3v = (-lv).exp() * S3j.v
                    if not (p3v > floor):
                        raise OracleRefusal(
                            "floor_rail",
                            f"p3 > 1e-9 not provable at row {r}")
                    phi = lam - S3j.log()
            if not phi.v.is_finite():
                raise OracleRefusal("nonfinite_oracle",
                                    f"phi not finite at row {r}")
            val += phi.v
            u[r] = phi.g[0]
            if order2:
                w[r] = phi.h[0]
        return val, u, w

    # ---------------- engine-facing contract -------------------------
    def nll(self, x):
        th = self._x_from_box(x)
        etas = self._etas(th)
        val, _, _ = self._row_jets(etas, False)
        rid = arb(0)
        for t in range(1, self.p):
            rid += th[t] * th[t]
        return val / self.n + self.l2 * rid

    def _grad_full(self, x):
        th = self._x_from_box(x)
        etas = self._etas(th)
        val, u, _ = self._row_jets(etas, False)
        g = self._A.transpose() * arb_mat([[v] for v in u])
        inv_n = 1 / arb(self.n)
        grad = [g[a, 0] * inv_n for a in range(self.p)]
        for t in range(1, self.p):
            grad[t] += 2 * self.l2 * th[t]
        return grad, val * inv_n

    def F(self, x):
        grad, _ = self._grad_full(x)
        Fz = [[grad[int(t)] for t in blk] for blk in self.part.blocks]
        Fg = [grad[int(t)] for t in self.part.border_idx]
        return Fz, Fg

    def F_full(self, x):
        return self.F(x)

    def _gram(self, w, rows_idx=None):
        """A^T diag(w) A (optionally row-restricted) as arb_mat."""
        if self._Arows is None:
            self._Arows = [[self._A[r, c] for c in range(self.p)]
                           for r in range(self.n)]
        if rows_idx is None:
            Aw = arb_mat([[v * w[r] for v in self._Arows[r]]
                          for r in range(self.n)])
            return self._A.transpose() * Aw
        rows = [self._Arows[r] for r in rows_idx]
        Asub = arb_mat(rows)
        Aw = arb_mat([[v * w[j] for v in rows[j]]
                      for j in range(len(rows))])
        return Asub.transpose() * Aw

    def H(self, x):
        """Hessian enclosure over the box. Radius-0 box: plain interval
        Gram. Positive radii: MEAN-VALUE FORM for the m0/m1/m2 rows
        (phi'' = phi''' = lam there):
          H(X) in H(mid) + sum_k [-1,1] |T_k(X)| r_k,
          T_k = (1/n) A^T diag(lam(X) A_:k) A   (SIGNED interval sums —
        the abs-row-sum interval hull loses the cross-row cancellation
        and, against the F1-collinear Gram (|Y|~1.7e3), blocked the
        polished certificate at the pinned radius floor; JOURNAL gate
        catch #4). m3 rows keep their (tiny, post-#2) interval
        contribution. Rails are censused over the FULL box either way."""
        th = self._x_from_box(x)
        rads = [float(v.rad()) for v in th]
        etas = self._etas(th)
        _, u, w = self._row_jets(etas, True)   # rails censused on box
        inv_n = 1 / arb(self.n)
        m3_idx = np.where(self.X == 3)[0]
        m012 = np.where(self.X != 3)[0]
        if max(rads) == 0.0:
            Hm = self._gram(w)
            Hfull = [[Hm[a, b] * inv_n for b in range(self.p)]
                     for a in range(self.p)]
        else:
            th_mid = [arb(float(v.mid())) for v in th]
            etas_mid = self._etas(th_mid)
            _, _, w_mid = self._row_jets(etas_mid, True)
            Hc = self._gram(w_mid)
            # m3 over the box replaces m3 at mid
            H3_box = self._gram(w, list(m3_idx))
            H3_mid = self._gram(w_mid, list(m3_idx))
            # width: sum_k |T_k| r_k over m012 rows (lam = w there)
            lam_rho = [None] * self.n
            width = [[arb(0)] * self.p for _ in range(self.p)]
            old_prec = ctx.prec
            ctx.prec = 80          # 2nd-order term; rigorous at any prec
            try:
                for k in range(self.p):
                    if rads[k] == 0.0:
                        continue
                    wk = [w[r] * self._Arows[r][k] for r in m012]
                    Tk = self._gram(wk, list(m012))
                    rk = arb(rads[k])
                    for a in range(self.p):
                        for b in range(self.p):
                            width[a][b] += Tk[a, b].abs_upper() * rk
            finally:
                ctx.prec = old_prec
            Hfull = [[None] * self.p for _ in range(self.p)]
            for a in range(self.p):
                for b in range(self.p):
                    Hfull[a][b] = (Hc[a, b] - H3_mid[a, b] + H3_box[a, b]
                                   + arb(0, width[a][b].abs_upper())
                                   ) * inv_n
        for t in range(1, self.p):
            Hfull[t][t] = Hfull[t][t] + 2 * self.l2
        for a in range(self.p):
            for b in range(self.p):
                if not Hfull[a][b].is_finite():
                    raise OracleRefusal("nonfinite_oracle",
                                        f"H ({a},{b}) not finite")
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

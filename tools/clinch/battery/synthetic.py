"""CLINCH battery synthetic fixture — a small HIERARCHICAL planted-optimum
problem with an ANALYTICALLY KNOWN unique optimum (exact rational solve).

Model (Gaussian, so the MAP is an exact linear solve): cells c with
  y_c ~ N(g . x_c + z_{b(c)}, 1),   priors z_b ~ N(0,1), g_j ~ N(0,1).
Penalized objective F(theta) = 0.5 sum_c (g.x_c + z_{b(c)} - y_c)^2
                             + 0.5 |z|^2 + 0.5 |g|^2  — strictly convex
quadratic; the unique optimum solves (A^T A + I) theta = A^T y EXACTLY
over the rationals (Fraction Gauss elimination — no float in the truth).

The oracle implements the baller block-arrow contract via clinch.jets —
the same jet machinery the v3.1 oracle uses — so the battery exercises
CLINCH's own evaluation path, not a shortcut. saddle=True flips the sign
of one latent's PRIOR curvature (+1 -> -1, applied consistently in F and
H): the stationary point remains, PD must refuse naming that block.
flip_math=True plants a sign error in one cell's covariate (a DIFFERENT
objective; the optimum moves; Krawczyk must refuse).

CENTERED VARIANT (battery leg L8). CenteredSynthOracle is the same
fixture with HARD MEAN-ZERO CENTERING: the cell mean uses z_{b(c)} -
mean(z), so every latent couples to every other through the mean and the
plain block-arrow form is lost — the exact problem shape the centered
model spaces (clinch.oracle_v36) carry. The oracle presents the
AUXILIARY-BORDER KKT EXTENSION at rank 2: border = globals + (m, mu)
with L = Phi(z, g, m) + mu (m - mean(z)), declares kkt_multiplier_idx,
and the composed optimum stays an exact rational solve
(exact_centered_optimum: Fraction elimination on the composed quadratic,
then m = mean(z), mu = sum of residuals — exact multiplier resolution).
gflip scales the g_0 prior curvature: gflip=-6 makes the REDUCED Hessian
indefinite in a border direction while the KKT point stays unique — the
kkt-signature control. saddle_block plants latent-block indefiniteness —
the midpoint-Cholesky control.
"""
from fractions import Fraction

import numpy as np
from flint import arb, arb_mat

from clinch.jets import Jet

# fixture: 2 latent blocks (families) x sizes (3, 2), border 3 globals
NZ = [3, 2]
NG = 3
CELLS = [
    # (block, latent-in-block, x (3,), y)
    (0, 0, (1.0, 0.5, 0.0), 1.2),
    (0, 0, (0.5, -0.5, 1.0), 0.3),
    (0, 1, (1.0, 0.0, 0.5), -0.4),
    (0, 2, (0.0, 1.0, 1.0), 0.9),
    (0, 2, (1.0, 1.0, 0.0), 1.5),
    (1, 0, (0.5, 0.5, 0.5), -0.2),
    (1, 0, (1.0, -1.0, 0.0), 0.4),
    (1, 1, (0.0, 0.5, 1.0), -0.8),
    (1, 1, (0.5, 0.0, -1.0), 0.6),
]


def exact_optimum():
    """The unique optimum, exact rationals via Fraction elimination."""
    n = sum(NZ) + NG                      # latent-major then globals
    A = [[Fraction(0)] * n for _ in range(n)]
    b = [Fraction(0)] * n
    for i in range(n):
        A[i][i] = Fraction(1)             # the prior
    zoff = [0, NZ[0]]
    for blk, li, x, y in CELLS:
        idx = [zoff[blk] + li] + [sum(NZ) + j for j in range(NG)]
        co = [Fraction(1)] + [Fraction(x[j]).limit_denominator(10**6)
                              for j in range(NG)]
        yf = Fraction(y).limit_denominator(10**6)
        for a, ca in zip(idx, co):
            b[a] += ca * yf
            for bb, cb in zip(idx, co):
                A[a][bb] += ca * cb
    # Gaussian elimination over Fraction
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for col in range(n):
        piv = next(r for r in range(col, n) if M[r][col] != 0)
        M[col], M[piv] = M[piv], M[col]
        pv = M[col][col]
        M[col] = [v / pv for v in M[col]]
        for r in range(n):
            if r != col and M[r][col] != 0:
                f = M[r][col]
                M[r] = [vr - f * vc for vr, vc in zip(M[r], M[col])]
    sol = [M[i][n] for i in range(n)]
    zs = [sol[0:NZ[0]], sol[NZ[0]:NZ[0] + NZ[1]]]
    g = sol[sum(NZ):]
    return zs, g


class SynthOracle:
    """baller block-arrow oracle for the fixture, via clinch jets."""
    dims = (NZ, NG)

    def __init__(self, saddle_block=None, flip_math=False):
        self.saddle_block = saddle_block
        self.flip_math = flip_math

    def _fgh(self, x, order2):
        zs, g = x
        n = sum(NZ) + NG
        grad = [arb(0)] * n
        acc = {} if order2 else None
        zoff = [0, NZ[0]]

        def gidx(j):
            return sum(NZ) + j

        def add(a, b_, v):
            key = (a, b_) if a <= b_ else (b_, a)
            acc[key] = acc.get(key, arb(0)) + v

        cells = list(CELLS)
        for ci, (blk, li, x_, y) in enumerate(cells):
            xv = list(x_)
            if self.flip_math and ci == 0:
                xv[0] = -xv[0]              # the planted sign error
            # r = g.x + z - y as a 1-var jet in r itself is overkill;
            # quadratic: grad contribution r * coef, hess coef_a*coef_b
            r = zs[blk][li] - y
            for j in range(NG):
                r = r + xv[j] * g[j]
            coords = [(zoff[blk] + li, 1.0)] + \
                     [(gidx(j), xv[j]) for j in range(NG)]
            # use a jet anyway so the battery exercises the jet path
            rj = Jet.var(r, 0, 1, order2)
            fj = rj * rj * 0.5
            for a, ca in coords:
                grad[a] += fj.g[0] * ca
            if order2:
                for i2 in range(len(coords)):
                    a, ca = coords[i2]
                    for j2 in range(i2, len(coords)):
                        b_, cb = coords[j2]
                        add(a, b_, fj.h[0] * (ca * cb))
        for t in range(n):
            sgn = 1.0
            if self.saddle_block is not None:
                lo = zoff[self.saddle_block]
                hi = lo + NZ[self.saddle_block]
                if lo <= t < hi:
                    # strong flipped curvature: -6 beats the +2 data
                    # curvature -> the block is truly indefinite (a saddle)
                    sgn = -6.0
            th_t = (zs[0][t] if t < NZ[0] else
                    zs[1][t - NZ[0]] if t < sum(NZ) else g[t - sum(NZ)])
            grad[t] += sgn * th_t
            if order2:
                add(t, t, arb(sgn))
        return grad, acc

    def F(self, x):
        grad, _ = self._fgh(x, False)
        return ([grad[0:NZ[0]], grad[NZ[0]:sum(NZ)]], grad[sum(NZ):])

    def H(self, x):
        _, acc = self._fgh(x, True)
        n0, n1 = NZ
        D = [arb_mat(n0, n0), arb_mat(n1, n1)]
        B = [arb_mat(n0, NG), arb_mat(n1, NG)]
        G = arb_mat(NG, NG)

        def region(t):
            if t < n0:
                return (0, t)
            if t < n0 + n1:
                return (1, t - n0)
            return (-1, t - n0 - n1)
        for (a, b_), v in acc.items():
            ra, oa = region(a)
            rb, ob = region(b_)
            if ra < 0 and rb < 0:
                G[oa, ob] = G[oa, ob] + v
                if oa != ob:
                    G[ob, oa] = G[ob, oa] + v
            elif ra >= 0 and rb >= 0:
                assert ra == rb
                D[ra][oa, ob] = D[ra][oa, ob] + v
                if oa != ob:
                    D[ra][ob, oa] = D[ra][ob, oa] + v
            elif ra >= 0:
                B[ra][oa, ob] = B[ra][oa, ob] + v
            else:
                B[rb][ob, oa] = B[rb][ob, oa] + v
        return D, B, G

    # a Partition-like shim so clinch.engine can name blocks
    class _Part:
        blocks = [np.arange(NZ[0]), NZ[0] + np.arange(NZ[1])]
        block_names = ["proc:famA", "proc:famB"]
        border_idx = sum(NZ) + np.arange(NG)
        dims = (NZ, NG)

    part = _Part()


# ---------------------------------------------------------------------------
# the CENTERED fixture (extended-KKT, rank-2 auxiliary border) — leg L8
# ---------------------------------------------------------------------------
NZT = sum(NZ)
NAUX = 2                              # m (auxiliary mean) + mu (multiplier)
_IM = NZT + NG                        # flat index of m
_IMU = NZT + NG + 1                   # flat index of mu


def exact_centered_optimum(gflip=1, saddle_block=None):
    """Exact rational stationary point of the CENTERED composed problem
    F(z, g) = 0.5 sum_c (g.x_c + z_{b(c)} - mean(z) - y_c)^2
            + 0.5 sum_t p_t z_t^2 + 0.5 sum_j q_j g_j^2
    (p, q = 1 except the planted controls), lifted to the extended KKT
    variables: returns (zs, g, m, mu) with m = mean(z*) and
    mu = sum_c r_c — the exact multiplier resolution. Fraction end to
    end; no float in the truth."""
    n = NZT + NG
    A = [[Fraction(0)] * n for _ in range(n)]
    b = [Fraction(0)] * n
    zoff = [0, NZ[0]]
    for t in range(NZT):
        p = Fraction(1)
        if saddle_block is not None and \
                zoff[saddle_block] <= t < zoff[saddle_block] + \
                NZ[saddle_block]:
            p = Fraction(-6)
        A[t][t] += p
    for j in range(NG):
        A[NZT + j][NZT + j] += (Fraction(gflip).limit_denominator(10**6)
                                if j == 0 else Fraction(1))
    inv_nz = Fraction(1, NZT)
    for blk, li, x, y in CELLS:
        co = [-inv_nz] * NZT + [Fraction(x[j]).limit_denominator(10**6)
                                for j in range(NG)]
        co[zoff[blk] + li] += Fraction(1)
        yf = Fraction(y).limit_denominator(10**6)
        for a in range(n):
            if co[a]:
                b[a] += co[a] * yf
                for bb in range(n):
                    if co[bb]:
                        A[a][bb] += co[a] * co[bb]
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for col in range(n):
        piv = next(r for r in range(col, n) if M[r][col] != 0)
        M[col], M[piv] = M[piv], M[col]
        pv = M[col][col]
        M[col] = [v / pv for v in M[col]]
        for r in range(n):
            if r != col and M[r][col] != 0:
                f = M[r][col]
                M[r] = [vr - f * vc for vr, vc in zip(M[r], M[col])]
    sol = [M[i][n] for i in range(n)]
    z = sol[:NZT]
    g = sol[NZT:]
    m = sum(z) * inv_nz
    mu = Fraction(0)
    for blk, li, x, y in CELLS:
        r = z[zoff[blk] + li] - m - Fraction(y).limit_denominator(10**6)
        for j in range(NG):
            r += Fraction(x[j]).limit_denominator(10**6) * g[j]
        mu += r
    return [z[:NZ[0]], z[NZ[0]:]], g, m, mu


class CenteredSynthOracle:
    """Extended-KKT block-arrow oracle for the centered fixture: border =
    [g, m, mu], blocks = the two latent blocks (block-arrow again — the
    latent-latent mean coupling is carried by the auxiliary border).
    Declares kkt_multiplier_idx, so the engine's PD leg routes to
    clinch.kkt_pd.pd_reduced_congruence."""
    dims = (NZ, NG + NAUX)
    kkt_multiplier_idx = (NG + 1,)    # border offset of mu

    def __init__(self, saddle_block=None, gflip=1.0):
        self.saddle_block = saddle_block
        self.gflip = float(gflip)

    def _fgh(self, x, order2):
        zs, gext = x
        n = NZT + NG + NAUX
        grad = [arb(0)] * n
        acc = {} if order2 else None
        zoff = [0, NZ[0]]
        m_v, mu_v = gext[NG], gext[NG + 1]

        def add(a, b_, v):
            key = (a, b_) if a <= b_ else (b_, a)
            acc[key] = acc.get(key, arb(0)) + v

        for blk, li, x_, y in CELLS:
            r = zs[blk][li] - m_v - y
            for j in range(NG):
                r = r + x_[j] * gext[j]
            coords = ([(zoff[blk] + li, 1.0)]
                      + [(NZT + j, x_[j]) for j in range(NG)]
                      + [(_IM, -1.0)])
            rj = Jet.var(r, 0, 1, order2)
            fj = rj * rj * 0.5
            for a, ca in coords:
                grad[a] += fj.g[0] * ca
            if order2:
                for i2 in range(len(coords)):
                    a, ca = coords[i2]
                    for j2 in range(i2, len(coords)):
                        b_, cb = coords[j2]
                        add(a, b_, fj.h[0] * (ca * cb))
        for t in range(NZT + NG):
            sgn = 1.0
            if t < NZT and self.saddle_block is not None:
                lo = zoff[self.saddle_block]
                if lo <= t < lo + NZ[self.saddle_block]:
                    sgn = -6.0
            if t == NZT:
                sgn = self.gflip
            th_t = (zs[0][t] if t < NZ[0] else
                    zs[1][t - NZ[0]] if t < NZT else gext[t - NZT])
            grad[t] += sgn * th_t
            if order2:
                add(t, t, arb(sgn))
        # the centering constraint c = m - mean(z) and its couplings
        zbar = arb(0)
        for t in range(NZT):
            zbar += (zs[0][t] if t < NZ[0] else zs[1][t - NZ[0]])
        zbar = zbar / NZT
        grad[_IMU] += m_v - zbar
        grad[_IM] += mu_v
        w = -1.0 / NZT
        for t in range(NZT):
            grad[t] += mu_v * w
        if order2:
            add(_IM, _IMU, arb(1))
            for t in range(NZT):
                add(t, _IMU, arb(w))
        return grad, acc

    def F(self, x):
        grad, _ = self._fgh(x, False)
        return ([grad[0:NZ[0]], grad[NZ[0]:NZT]], grad[NZT:])

    def H(self, x):
        _, acc = self._fgh(x, True)
        n0, n1 = NZ
        ngb = NG + NAUX
        D = [arb_mat(n0, n0), arb_mat(n1, n1)]
        B = [arb_mat(n0, ngb), arb_mat(n1, ngb)]
        G = arb_mat(ngb, ngb)

        def region(t):
            if t < n0:
                return (0, t)
            if t < NZT:
                return (1, t - n0)
            return (-1, t - NZT)
        for (a, b_), v in acc.items():
            ra, oa = region(a)
            rb, ob = region(b_)
            if ra < 0 and rb < 0:
                G[oa, ob] = G[oa, ob] + v
                if oa != ob:
                    G[ob, oa] = G[ob, oa] + v
            elif ra >= 0 and rb >= 0:
                assert ra == rb
                D[ra][oa, ob] = D[ra][oa, ob] + v
                if oa != ob:
                    D[ra][ob, oa] = D[ra][ob, oa] + v
            elif ra >= 0:
                B[ra][oa, ob] = B[ra][oa, ob] + v
            else:
                B[rb][ob, oa] = B[rb][ob, oa] + v
        return D, B, G

    class _Part:
        blocks = [np.arange(NZ[0]), NZ[0] + np.arange(NZ[1])]
        block_names = ["proc:famA", "proc:famB"]
        border_idx = NZT + np.arange(NG + NAUX)
        dims = (NZ, NG + NAUX)

    part = _Part()

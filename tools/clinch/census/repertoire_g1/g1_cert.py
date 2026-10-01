"""Gate-1 complete-census certified layer (arb ball arithmetic end-to-end).

Certificate classes (every support of every cell is disposed by exactly one;
the disposition array + payloads ARE the completeness certificate):
  V    dead clone: g_i = sum_j b_ij - d_i < 0 certified => G_i(x) <= g_i < 0
       for all x>=0 (monotone), no equilibrium with i active (inherits upward).
  INV  conflict pair (i,l): certified M_il > 0 with
       M_il = min_{phi in [phimin,1]^m} sum_j (b_lj-b_ij) phi_j - (d_l-d_i);
       any equilibrium with i in S, l notin S has transversal eigenvalue
       G_l(x*) >= M_il > 0 => not stable. (Stable-census-complete pruning.)
  EXCL |S| = m+1 over live clones: certified det[B_S | d_S] != 0 => the active
       balance B_S p = d_S is unsolvable; infeasibility inherits to all
       supersets => no equilibria with |S| > m.
  BB   Krawczyk branch-and-bound over [0, Xmax_S]: every leaf excluded
       (exact-range corner test on the componentwise-decreasing G, or empty
       Krawczyk intersection) or absorbed by a verified Krawczyk root box.
       Verified roots get certified positivity + stability sign
       (block-triangular full Jacobian: Gershgorin on V^-1 J V in acb for the
       internal block; ball signs of G_l for the transversal eigenvalues).
UNKNOWN/UNCERTIFIED are recorded verbatim, never forced (D0 law).
Precision pinned: 128 bits (the measured precision; no escalation ladder).
"""
import itertools
import time

import numpy as np
from flint import arb, arb_mat, acb, acb_mat, ctx

import g1_model as M

PREC = 128
BOX_CAP = 20000          # B&B effort cap per support -> UNCERTIFIED
MINREL = 1e-10           # unresolvable-feature floor -> UNCERTIFIED
ctx.prec = PREC

D_V, D_INV, D_BB_EMPTY, D_ROOTS, D_UNCERT = 1, 2, 4, 5, 6


def ball(lo, hi):
    a, b = arb(float(lo)), arb(float(hi))
    return (a + b) / 2 + ((b - a) / 2) * arb("0 +/- 1")


def _lo(x):   # outward float lower bound of a ball
    v = float(x.mid()) - float(x.rad())
    return np.nextafter(np.nextafter(v, -np.inf), -np.inf)


def _hi(x):
    v = float(x.mid()) + float(x.rad())
    return np.nextafter(np.nextafter(v, np.inf), np.inf)


def amat(rows):
    return arb_mat([[arb(float(v)) for v in row] for row in rows])


def cabs_ub(z):
    """Certified UPPER bound of |z| over an acb ball (Gershgorin radii).
    (abs of a straddling arb ball keeps the straddle in python-flint, and
    squaring it straddles 0 => sqrt = nan; build from mid/rad instead.)"""
    a = abs(z.real.mid()) + z.real.rad()
    b = abs(z.imag.mid()) + z.imag.rad()
    return (a * a + b * b).sqrt()


def cabs_lb(z):
    """Certified LOWER bound of |z| over an acb ball (disc disjointness)."""
    a, b = z.real.mid(), z.imag.mid()
    return (a * a + b * b).sqrt() - (z.real.rad() + z.imag.rad())


class CellCert:
    def __init__(self, cell):
        self.cell = cell
        self.n, self.m = cell["n"], cell["m"]
        self.K = arb(float(cell["K"]))
        self.h = cell["h"]
        self.B = amat(cell["b"])                     # n x m
        self.C = amat(cell["c"])                     # n x m
        self.dv = [arb(float(v)) for v in cell["d"]]

    # ---- phi (h in {1,2} => rational; rigorous on balls) ----
    def phi(self, y):
        t = y / self.K
        return 1 / (1 + t) if self.h == 1 else 1 / (1 + t * t)

    def dphi(self, y):
        t = y / self.K
        if self.h == 1:
            f = 1 / (1 + t)
            return -(f * f) / self.K
        f = 1 / (1 + t * t)
        return -(2 * t / self.K) * (f * f)

    # ---- clone viability ----
    def viability(self):
        g = []
        for i in range(self.n):
            s = arb(0)
            for j in range(self.m):
                s += self.B[i, j]
            g.append(s - self.dv[i])
        dead = [i for i in range(self.n) if g[i] < 0]
        viable = [i for i in range(self.n) if g[i] > 0]
        vunk = [i for i in range(self.n) if i not in dead and i not in viable]
        return g, dead, viable, vunk

    # ---- certified per-clone equilibrium bound ----
    def xmax(self, i, fhint):
        def h_i(t):
            s = arb(0)
            for j in range(self.m):
                s += self.B[i, j] * self.phi(self.C[i, j] * t)
            return s - self.dv[i]
        hi = max(fhint * 1.5, 1e-6) if fhint else 1.0
        for _ in range(200):
            if h_i(arb(hi)) < 0:
                break
            hi *= 2.0
        else:
            return None
        lo = 0.0
        for _ in range(50):
            mid = 0.5 * (lo + hi)
            if h_i(arb(mid)) < 0:
                hi = mid
            else:
                lo = mid
        return hi

    # ---- support-restricted evaluation ----
    def sub(self, S):
        Bs = arb_mat([[self.B[i, j] for j in range(self.m)] for i in S])
        Cs = arb_mat([[self.C[i, j] for j in range(self.m)] for i in S])
        ds = arb_mat([[self.dv[i]] for i in S])
        return Bs, Cs, ds

    def G_at(self, Bs, Cs, ds, xs):
        k = len(xs)
        y = [arb(0) for _ in range(self.m)]
        for j in range(self.m):
            s = arb(0)
            for i in range(k):
                s += Cs[i, j] * xs[i]
            y[j] = s
        ph = arb_mat([[self.phi(y[j])] for j in range(self.m)])
        return Bs * ph - ds                                    # k x 1 balls

    def corner_excludes(self, Bs, Cs, ds, lo, hi):
        """Exact-range test: G componentwise non-increasing (b,c>=0, phi dec):
        max over box at lo-corner, min at hi-corner."""
        gmax = self.G_at(Bs, Cs, ds, [arb(float(v)) for v in lo])
        for i in range(len(lo)):
            if gmax[i, 0] < 0:
                return True
        gmin = self.G_at(Bs, Cs, ds, [arb(float(v)) for v in hi])
        for i in range(len(lo)):
            if gmin[i, 0] > 0:
                return True
        return False

    def ybox(self, Cs, lo, hi):
        """Rigorous per-niche crowding balls over the x-box (all in arb)."""
        k = len(lo)
        xb = [ball(lo[i], hi[i]) for i in range(k)]
        ys = []
        for j in range(self.m):
            s = arb(0)
            for i in range(k):
                s += Cs[i, j] * xb[i]
            ys.append(s)
        return ys

    def krawczyk(self, cell, S, Bs, Cs, ds, lo, hi):
        """Returns ('empty'|'contained'|'split', tightened (lo,hi) or None)."""
        k = len(S)
        x0 = 0.5 * (np.asarray(lo) + np.asarray(hi))
        G0 = self.G_at(Bs, Cs, ds, [arb(float(v)) for v in x0])
        Jf = M.f_DG(cell, list(S), x0)
        try:
            Yf = np.linalg.inv(Jf)
        except np.linalg.LinAlgError:
            return "split", None
        # interval Jacobian: DG = Bs * diag(dphi(Y_j)) * Cs^T on y-range balls
        ys = self.ybox(Cs, lo, hi)
        Dd = [[arb(0)] * self.m for _ in range(self.m)]
        for j in range(self.m):
            Dd[j][j] = self.dphi(ys[j])
        DGX = Bs * arb_mat(Dd) * Cs.transpose()
        Yb = amat(Yf)
        Ik = arb_mat([[arb(1) if i == j else arb(0) for j in range(k)]
                      for i in range(k)])
        dx = arb_mat([[ball(lo[i], hi[i]) - arb(float(x0[i]))]
                      for i in range(k)])
        x0m = arb_mat([[arb(float(v))] for v in x0])
        Kx = x0m - Yb * G0 + (Ik - Yb * DGX) * dx
        contained, empty = True, False
        nlo, nhi = list(lo), list(hi)
        for i in range(k):
            ki = Kx[i, 0]
            if ki < arb(float(lo[i])) or ki > arb(float(hi[i])):
                empty = True
                break
            if not (ki > arb(float(lo[i])) and ki < arb(float(hi[i]))):
                contained = False
            nlo[i] = max(lo[i], _lo(ki))
            nhi[i] = min(hi[i], _hi(ki))
        if empty:
            return "empty", None
        if contained:
            return "contained", (nlo, nhi)
        return "split", (nlo, nhi)

    # ---- k = m exact linear route (certificate class LIN) ----
    def lin_route(self, S, Bs, Cs, ds, xmax):
        """For |S| = m the active balance B_S p = d_S is square: certified
        solve => unique p*; p* outside (0, 1] or below the phi floor => no
        equilibrium; else y* = phi^-1(p*), x* = solve(C_S^T, y*) unique.
        Returns ('empty', None) | ('root', x*_ball_list) | ('fall', None).
        Nonsingularity is certified by the ball solves themselves."""
        k = len(S)
        if k != self.m:
            return "fall", None
        try:
            pstar = Bs.solve(ds)                      # m x 1 balls
        except Exception:
            return "fall", None
        # phi floor at this support's own crowding ceiling
        for j in range(self.m):
            ym = arb(0)
            for i in range(k):
                ym += Cs[i, j] * arb(xmax[i])
            if pstar[j, 0] < self.phi(ym) or pstar[j, 0] > 1:
                return "empty", None
            if not (pstar[j, 0] > 0 and pstar[j, 0] < 1):
                return "fall", None                   # can't invert safely
        ys = []
        for j in range(self.m):
            q = 1 / pstar[j, 0] - 1
            if not (q > 0):
                return "fall", None
            ys.append(self.K * q if self.h == 1 else self.K * q.sqrt())
        try:
            xstar = Cs.transpose().solve(arb_mat([[y] for y in ys]))
        except Exception:
            return "fall", None
        for i in range(k):
            if xstar[i, 0] < 0:
                return "empty", None                  # certified infeasible
        if all(xstar[i, 0] > 0 for i in range(k)):
            return "root", [xstar[i, 0] for i in range(k)]
        return "fall", None

    # ---- stability of a verified root ----
    def stability(self, cell, S, encl_lo, encl_hi, live_others):
        k = len(S)
        Bs, Cs, ds = self.sub(S)
        yb = self.ybox(Cs, encl_lo, encl_hi)
        ph = [self.phi(y) for y in yb]
        inv_signs = {}
        any_pos = False
        any_unk = False
        for l in live_others:
            s = arb(0)
            for j in range(self.m):
                s += self.B[l, j] * ph[j]
            gl = s - self.dv[l]
            if gl < 0:
                inv_signs[l] = -1
            elif gl > 0:
                inv_signs[l] = 1
                any_pos = True
            else:
                inv_signs[l] = 0
                any_unk = True
        # internal block J = diag(x) * Bs * diag(dphi) * Cs^T
        Dd = [[arb(0)] * self.m for _ in range(self.m)]
        for j in range(self.m):
            Dd[j][j] = self.dphi(yb[j])
        Jin = Bs * arb_mat(Dd) * Cs.transpose()
        Xd = [[arb(0)] * k for _ in range(k)]
        for i in range(k):
            Xd[i][i] = ball(encl_lo[i], encl_hi[i])
        Jb = arb_mat(Xd) * Jin
        Jf = np.array([[float(Jb[i, j].mid()) for j in range(k)]
                       for i in range(k)])
        internal = "UNKNOWN"
        try:
            w, V = np.linalg.eig(Jf)
            Vi = np.linalg.inv(V)
            Va = acb_mat([[acb(V[i, j].real, V[i, j].imag) for j in range(k)]
                          for i in range(k)])
            Via = acb_mat([[acb(Vi[i, j].real, Vi[i, j].imag)
                            for j in range(k)] for i in range(k)])
            Ja = acb_mat([[acb(Jb[i, j]) for j in range(k)] for i in range(k)])
            T = Via * Ja * Va
            centers = [T[i, i] for i in range(k)]
            radii = []
            for i in range(k):
                r = arb(0)
                for j in range(k):
                    if j != i:
                        r += cabs_ub(T[i, j])
                radii.append(r)
            if all((centers[i].real + radii[i]) < 0 for i in range(k)):
                internal = "STABLE"
            else:
                for i in range(k):
                    if (centers[i].real - radii[i]) > 0:
                        disjoint = all(
                            cabs_lb(centers[i] - centers[j])
                            > radii[i] + radii[j]
                            for j in range(k) if j != i)
                        if disjoint:
                            internal = "UNSTABLE"
                            break
        except np.linalg.LinAlgError:
            pass
        if any_pos or internal == "UNSTABLE":
            verdict = "UNSTABLE"
        elif internal == "STABLE" and not any_unk:
            verdict = "STABLE"
        else:
            verdict = "UNKNOWN"
        return verdict, inv_signs, internal

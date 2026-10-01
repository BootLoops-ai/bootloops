"""Different-I multi-sample Etienne oracle — the reference Panama-65 fit
(66 parameters: lntheta + lnI_p for 65 plots; the study's
named hole: "extending interval Newton to the [different-I] kernel
derivative is routine but not done").

STRUCTURE NOTE (soundness): the I_p x I_q Hessian couplings are NONZERO
(plots couple through shared species AND through the metacommunity
(theta)_A sum), so a "65 independent plot blocks" sketch is
NOT a valid block-arrow presentation of the exact Hessian. The sound
presentation used here: ONE dense 65-dim lnI block + the 1-dim lntheta
border. baller.certify.block_krawczyk consumes it unchanged.

MATH. With per-plot factors f_sp(z) = sum_{a=1}^{n_sp} sbar(n_sp,a) I_p^a
z^a, per-species u_s = prod_{p in O_s} f_sp, P_s = Borel(u_s) (coeff b
scaled by (b-1)!), M = prod_s P_s, and w_A = 1/(theta)_A:

  lnP = statics + S lntheta + ln F - sum_p [lgam(I_p+J_p) - lgam(I_p)]
  F   = <M, w> = sum_A M_A w_A

All lnI-derivatives of M are leaf-linear, so every contraction reduces to
(i) LEAF CONTEXTS c_s (one product-tree down-sweep per weight vector:
c_child = corr(sibling value, c_parent)) for the single-species terms, and
(ii) NODE-SPLIT PAIR TERMS for the cross-species terms:
  <d2_pq M, w> = sum_s <D2_s[pq], c_s>
               + sum_{nodes n} [B_{c_n}(g_p(L), g_q(R)) + (p<->q)]
with g_p(node) the subtree gradient polys (product rule up the tree) and
B_c(X,Y) = sum_{ij} X_i Y_j c_{i+j} (computed as corr + arb_mat product).
theta-rows use the psiD/psiT running recurrences (phase2_certify pattern)
against the root polynomial; the theta-I mixed row uses a second context
sweep with weight psiD_A w_A. Per-plot prefactor digamma/trigamma
differences are exact finite recurrences sum_k 1/(I_p+k) (no ball-digamma
anywhere). Every operation is arb ball arithmetic over the input box:
enclosures are rigorous by construction. NLL = -lnP sign convention.
"""
import math

from flint import arb, arb_poly, arb_mat, fmpz_poly

from . import ident
from .common import Part, cur_prec

from clinch.oracle_v31 import OracleRefusal


def _corr(P, c, need):
    """corr(P, c)[j] = sum_i P_i c_{i+j}, j = 0..need-1 (one poly mult)."""
    n = P.length()
    rev = arb_poly([P[n - 1 - i] for i in range(n)])
    prod = rev * c
    return arb_poly([prod[n - 1 + j] for j in range(need)])


def _dot(P, c, upto=None):
    """<P, c> = sum_i P_i c_i (c an arb_poly context)."""
    s = arb(0)
    n = P.length() if upto is None else min(P.length(), upto)
    for i in range(n):
        s += P[i] * c[i]
    return s


class _Node:
    __slots__ = ("left", "right", "leaf_idx", "poly", "supp")

    def __init__(self, leaf_idx=None, left=None, right=None):
        self.leaf_idx = leaf_idx
        self.left = left
        self.right = right
        self.poly = None
        self.supp = None


def _build_tree(nleaf):
    nodes = [_Node(leaf_idx=i) for i in range(nleaf)]
    while len(nodes) > 1:
        nxt = []
        for i in range(0, len(nodes) - 1, 2):
            nxt.append(_Node(left=nodes[i], right=nodes[i + 1]))
        if len(nodes) % 2:
            nxt.append(nodes[-1])
        nodes = nxt
    return nodes[0]


class EtienneDiffIOracle:
    """rows: S x N abundance matrix (list of tuples)."""

    def __init__(self, rows, tag=""):
        ident.pilot()  # pin verification (reference identity)
        self.rows = [tuple(int(x) for x in r) for r in rows]
        self.N = len(self.rows[0])
        self.S = len(self.rows)
        self.Js = [sum(r[i] for r in self.rows) for i in range(self.N)]
        self.J = sum(self.Js)
        self.tag = tag
        self.occ = [[p for p in range(self.N) if r[p] > 0]
                    for r in self.rows]
        self._rf = {}
        # flat layout: [0] = lntheta (border), [1..N] = lnI_p (one block)
        self.part = Part(blocks=[list(range(1, self.N + 1))],
                         block_names=["plots:lnI"], border_idx=[0])
        self._fact = None

    @property
    def dims(self):
        return ([self.N], 1)

    # ---------------- statics ---------------------------------------
    def _statics(self):
        from collections import Counter
        st = arb(0)
        for cnt in Counter(self.rows).values():
            st -= arb(cnt + 1).lgamma()
        for i in range(self.N):
            st += arb(self.Js[i] + 1).lgamma()
        for r in self.rows:
            for n in r:
                if n > 1:
                    st -= arb(n + 1).lgamma()
        return st

    def _rfints(self, n):
        if n not in self._rf:
            p = ident.pilot()["multisample_ball"].rf_poly_int(n)
            self._rf[n] = [int(p[a]) for a in range(n + 1)]
        return self._rf[n]

    def _factorials(self, upto):
        if self._fact is None or len(self._fact) <= upto:
            f = [arb(1)] * (upto + 1)
            for b in range(2, upto + 1):
                f[b] = f[b - 1] * (b - 1)
            self._fact = f          # f[b] = (b-1)! for b >= 1
        return self._fact

    # ---------------- per-species polys ------------------------------
    def _species_polys(self, Ivec, order):
        """Value polys P_s; if order >= 1 also D1[s][p]; if order >= 2 also
        D2[s][(p,q)] (p <= q). Borel twist applied to all."""
        fact = self._factorials(self.J + 2)

        def twist(u):
            return arb_poly([u[b] * fact[b] if b >= 1 else arb(0)
                             for b in range(u.length())])
        P, D1, D2 = [], [], []
        one = arb_poly([arb(1)])
        for s, row in enumerate(self.rows):
            occ = self.occ[s]
            fs, f1s, f2s = [], [], []
            for p in occ:
                n = row[p]
                rf = self._rfints(n)
                Ia = arb(1)
                c0, c1, c2 = [arb(0)], [arb(0)], [arb(0)]
                for a in range(1, n + 1):
                    Ia = Ia * Ivec[p] if a > 1 else Ivec[p]
                    base = arb(rf[a]) * Ia
                    c0.append(base)
                    c1.append(base * a)
                    c2.append(base * (a * a))
                fs.append(arb_poly(c0))
                f1s.append(arb_poly(c1))
                f2s.append(arb_poly(c2))
            k = len(occ)
            pre = [one] * (k + 1)
            for i in range(k):
                pre[i + 1] = pre[i] * fs[i]
            suf = [one] * (k + 1)
            for i in range(k - 1, -1, -1):
                suf[i] = suf[i + 1] * fs[i]
            P.append(twist(pre[k]))
            if order >= 1:
                d1 = {}
                for i, p in enumerate(occ):
                    d1[p] = twist(f1s[i] * (pre[i] * suf[i + 1]))
                D1.append(d1)
            if order >= 2:
                d2 = {}
                for i, p in enumerate(occ):
                    d2[(p, p)] = twist(f2s[i] * (pre[i] * suf[i + 1]))
                    mid = one
                    for jj in range(i + 1, k):
                        q = occ[jj]
                        d2[(p, q)] = twist(
                            (f1s[i] * f1s[jj]) * (pre[i] * mid * suf[jj + 1]))
                        mid = mid * fs[jj]
                D2.append(d2)
        return P, D1, D2

    # ---------------- tree passes ------------------------------------
    def _up_values(self, node, P):
        if node.leaf_idx is not None:
            node.poly = P[node.leaf_idx]
            node.supp = set(self.occ[node.leaf_idx])
            return
        self._up_values(node.left, P)
        self._up_values(node.right, P)
        node.poly = node.left.poly * node.right.poly
        node.supp = node.left.supp | node.right.supp

    def _down_contexts(self, node, ctxs, leaf_ctxs):
        """ctxs: list of context polys for this node (one per weight)."""
        if node.leaf_idx is not None:
            leaf_ctxs[node.leaf_idx] = ctxs
            return
        Lp, Rp = node.left.poly, node.right.poly
        cl = [_corr(Rp, c, Lp.length()) for c in ctxs]
        cr = [_corr(Lp, c, Rp.length()) for c in ctxs]
        self._down_contexts(node.left, cl, leaf_ctxs)
        self._down_contexts(node.right, cr, leaf_ctxs)

    def _pair_pass(self, node, ctx_of, D1, T2, top=False):
        """Post-order: returns grad dict {p: poly} of the subtree; adds
        this node's cross-child pair contributions into T2 (N x N python
        list of arb). ctx_of maps id(node)->context poly (weight w).
        top=True skips the (unused) root gradient combine."""
        if node.leaf_idx is not None:
            return dict(D1[node.leaf_idx])
        gL = self._pair_pass(node.left, ctx_of, D1, T2)
        gR = self._pair_pass(node.right, ctx_of, D1, T2)
        c = ctx_of[id(node)]
        supL = sorted(gL.keys())
        supR = sorted(gR.keys())
        if supL and supR:
            lenR = node.right.poly.length()
            # U[p] = corr(gL[p], c) truncated to lenR
            Umat = arb_mat(len(supL), lenR)
            for a, p in enumerate(supL):
                u = _corr(gL[p], c, lenR)
                for j in range(lenR):
                    Umat[a, j] = u[j]
            Gmat = arb_mat(lenR, len(supR))
            for b, q in enumerate(supR):
                gq = gR[q]
                ln = gq.length()
                for j in range(lenR):
                    Gmat[j, b] = gq[j] if j < ln else arb(0)
            M2 = Umat * Gmat
            for a, p in enumerate(supL):
                for b, q in enumerate(supR):
                    v = M2[a, b]
                    T2[p][q] += v
                    T2[q][p] += v
        # combine grads up
        if top:
            return {}
        out = {}
        Lp, Rp = node.left.poly, node.right.poly
        for p in set(supL) | set(supR):
            if p in gL and p in gR:
                out[p] = gL[p] * Rp + Lp * gR[p]
            elif p in gL:
                out[p] = gL[p] * Rp
            else:
                out[p] = Lp * gR[p]
        return out

    def _ctx_index(self, node, ctxs, ctx_of, which):
        """Record per-node context (weight `which`) during down-sweep."""
        ctx_of[id(node)] = ctxs[which]
        if node.leaf_idx is not None:
            return
        Lp, Rp = node.left.poly, node.right.poly
        cl = [_corr(Rp, c, Lp.length()) for c in ctxs]
        cr = [_corr(Lp, c, Rp.length()) for c in ctxs]
        self._ctx_index(node.left, cl, ctx_of, which)
        self._ctx_index(node.right, cr, ctx_of, which)

    # ---------------- weights ----------------------------------------
    def _weights(self, th, Amax, with_psi):
        """w_A = 1/(theta)_A and (optionally) wp_A = psiD_A w_A, plus the
        psi moment coefficient lists for the theta rows."""
        w = [arb(0)] * (Amax + 1)
        wp = [arb(0)] * (Amax + 1) if with_psi else None
        psiD = arb(0)
        psiT = arb(0)
        t = arb(1)
        psiDs = [arb(0)] * (Amax + 1)
        psiTs = [arb(0)] * (Amax + 1)
        for A in range(Amax + 1):
            w[A] = t
            psiDs[A] = psiD
            psiTs[A] = psiT
            if with_psi:
                wp[A] = psiD * t
            t = t / (th + A)
            psiD = psiD + 1 / (th + A)
            psiT = psiT - 1 / (th + A) ** 2
        return (arb_poly(w), arb_poly(wp) if with_psi else None,
                psiDs, psiTs)

    def _plot_digamma(self, Ivec):
        """Exact finite recurrences: (psiD_p, psiT_p) for each plot,
        psiD_p = psi(I_p+J_p)-psi(I_p) = sum_{k<J_p} 1/(I_p+k)."""
        out = []
        for p in range(self.N):
            I = Ivec[p]
            s1 = arb(0)
            s2 = arb(0)
            for k in range(self.Js[p]):
                d = 1 / (I + k)
                s1 += d
                s2 -= d * d
            out.append((s1, s2))
        return out

    # ---------------- main evaluation --------------------------------
    def _unpack(self, x):
        zs, g = x
        th = g[0].exp()
        Ivec = [z.exp() for z in zs[0]]
        if not th.is_finite() or not all(v.is_finite() for v in Ivec):
            raise OracleRefusal("domain", "exp(u) not finite over the box")
        return th, Ivec

    def evaluate(self, x, order=1):
        """Returns dict with lnP (order 0), grad (order 1), hess (order 2)
        in u = (lntheta, lnI_1..lnI_N); ball-rigorous over the box."""
        th, Ivec = self._unpack(x)
        N = self.N
        P, D1, D2 = self._species_polys(Ivec, order)
        root = _build_tree(self.S)
        self._up_values(root, P)
        M = root.poly
        Amax = M.length() - 1
        w, wp, psiDs, psiTs = self._weights(th, Amax, with_psi=(order >= 1))
        F = _dot(M, w)
        # Fail-closed register: over a fat box the contraction F may not
        # be provably positive; every downstream enclosure is then
        # nonfinite and block_krawczyk refuses the rung as
        # nonfinite_oracle (a NAMED refusal). Raising instead would abort
        # the whole ladder (the engine's best-rung report needs the
        # returned-dict path). SHORT-CIRCUIT: once F is not provably
        # positive no downstream quantity can be finite, so return the
        # NaN enclosures immediately instead of grinding the full
        # context/pair machinery through NaN balls (measured: 30 min ->
        # ~1 min per doomed rung; identical verdict semantics).
        out = {}
        if not (F > 0):
            nan = arb("nan")
            out["lnP"] = nan
            if order >= 1:
                out["grad"] = [nan] * (N + 1)
            if order >= 2:
                out["hess"] = [[nan] * (N + 1) for _ in range(N + 1)]
            return out
        statics = self._statics()
        lnP = statics + arb(self.S) * th.log() + F.log()
        pdg = self._plot_digamma(Ivec)
        for p in range(N):
            I = Ivec[p]
            lnP -= (I + self.Js[p]).lgamma() - I.lgamma()
        out["lnP"] = lnP
        if order < 1:
            return out
        # theta moments from the root
        Spsi = arb(0)
        Spsi2 = arb(0)
        SpsiT = arb(0)
        for A in range(self.S, Amax + 1):
            T = M[A] * w[A]
            Spsi += psiDs[A] * T
            Spsi2 += psiDs[A] * psiDs[A] * T
            SpsiT += psiTs[A] * T
        EY = Spsi / F
        # leaf contexts for w (and wp)
        leaf_ctxs = [None] * self.S
        ctxs0 = [w] + ([wp] if order >= 1 else [])
        self._down_contexts(root, ctxs0, leaf_ctxs)
        G = [arb(0)] * N
        Gpsi = [arb(0)] * N
        for s in range(self.S):
            cs = leaf_ctxs[s][0]
            cps = leaf_ctxs[s][1]
            for p, d1 in D1[s].items():
                G[p] += _dot(d1, cs)
                Gpsi[p] += _dot(d1, cps)
        grad = [arb(0)] * (N + 1)
        grad[0] = arb(self.S) - th * EY
        for p in range(N):
            grad[1 + p] = G[p] / F - Ivec[p] * pdg[p][0]
        out["grad"] = grad
        if order < 2:
            return out
        # Hessian
        H = [[arb(0)] * (N + 1) for _ in range(N + 1)]
        H[0][0] = (-th * EY
                   + th * th * (Spsi2 / F - EY * EY - SpsiT / F))
        for p in range(N):
            v = -th * (Gpsi[p] / F - EY * (G[p] / F))
            H[0][1 + p] = v
            H[1 + p][0] = v
        # single-species second-derivative terms
        C = [[arb(0)] * N for _ in range(N)]
        for s in range(self.S):
            cs = leaf_ctxs[s][0]
            for (p, q), d2 in D2[s].items():
                v = _dot(d2, cs)
                C[p][q] += v
                if p != q:
                    C[q][p] += v
        # cross-species pair terms (node-split scheme)
        ctx_of = {}
        self._ctx_index(root, [w], ctx_of, 0)
        T2 = [[arb(0)] * N for _ in range(N)]
        self._pair_pass(root, ctx_of, D1, T2, top=True)
        for p in range(N):
            for q in range(N):
                C[p][q] += T2[p][q]
        for p in range(N):
            for q in range(N):
                Hpq = C[p][q] / F - (G[p] / F) * (G[q] / F)
                if p == q:
                    I = Ivec[p]
                    Hpq += -I * pdg[p][0] - I * I * pdg[p][1]
                H[1 + p][1 + q] = Hpq
        out["hess"] = H
        return out

    # ---------------- engine-facing contract -------------------------
    def nll(self, x):
        return -self.evaluate(x, order=0)["lnP"]

    def F(self, x):
        g = self.evaluate(x, order=1)["grad"]
        return [[-g[1 + p] for p in range(self.N)]], [-g[0]]

    def H(self, x):
        r = self.evaluate(x, order=2)
        Hm = r["hess"]
        N = self.N
        D = arb_mat(N, N)
        for p in range(N):
            for q in range(N):
                D[p, q] = -Hm[1 + p][1 + q]
        B = arb_mat(N, 1)
        for p in range(N):
            B[p, 0] = -Hm[0][1 + p]
        G = arb_mat([[-Hm[0][0]]])
        self._last = r
        return [D], [B], G

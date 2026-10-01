"""Corner oracle for the m=1 (Ewens) boundary fits — the 16 microbiome
seawater per-sample fits of the reference study ("16/16 certified
boundary per-sample").

COORDINATES. u = (lntheta, d, q) with q = 1/I. m -> 1 corresponds to
q -> 0+, so the published boundary optimum is the DECLARED LOWER BOUND
q = 0 (active, exactly). In q the likelihood is analytic through the
boundary:

  lnP(theta, q) = statics(D) + S lntheta
                  + ln( sum_A K_A q^{J-A} / (theta)_A )
                  - sum_{k=0}^{J-1} ln(1 + k q)

(the identity lgamma(I+J)-lgamma(I) = J ln I + sum_k ln(1+kq) removes the
I -> inf divergence exactly). At q=0 the sum collapses to K_J/(theta)_J
(the Ewens likelihood; K_J = 1) and the q-gradient needs only
K_{J-1} = sum_{n_i>=2} n_i/2 — closed forms, gated against K_ball in the
battery.

SYSTEM PRESENTED (block-arrow contract needs >=1 block after masking, and
masking q empties its block, so an auxiliary identity coordinate d is
carried: F_d = d, H_dd = 1, cross terms exactly 0 — the augmented system's
stationary points are (d=0, x*) with block-diagonal Hessian, so existence/
uniqueness/PD statements transfer verbatim to the real system; the aux
block is named in the certificate):
  blocks = [[1]] (aux:identity, d), [[2]] (local:q)   border = [0] lntheta
NLL sign convention throughout. The oracle only evaluates at q pinned to
EXACT 0 (the MaskedOracle path); q != 0 raises a typed refusal — the
interior register for these fits would need the full kernel and is not
this oracle's scope. Dropped H rows/cols for q are NaN by construction
(fail-closed if ever consumed).
"""
from flint import arb, arb_mat, acb

from .common import Part, cur_prec

from clinch.oracle_v31 import OracleRefusal


class EtienneCornerOracle:
    def __init__(self, D, tag=""):
        self.D = sorted(int(x) for x in D)
        if min(self.D) < 1:
            raise OracleRefusal("data", "abundances must be >= 1")
        self.J, self.S = sum(self.D), len(self.D)
        self.tag = tag
        # closed-form top kernel coefficients (exact in binary arb)
        self.K_J = arb(1)
        s = arb(0)
        for n in self.D:
            if n >= 2:
                s += arb(n) / 2
        self.K_Jm1 = s
        self.part = Part(blocks=[[1], [2]],
                         block_names=["aux:identity", "local:q"],
                         border_idx=[0])

    @property
    def dims(self):
        return ([1, 1], 1)

    def _pieces(self, x):
        zs, g = x
        u_th = g[0]
        d = zs[0][0]
        q = zs[1][0]
        if not (q.is_zero()):
            raise OracleRefusal(
                "q_nonzero",
                "corner oracle evaluates the boundary register only "
                "(q pinned exactly 0 by the mask); interior-q evaluation "
                "is out of this oracle's scope")
        th = u_th.exp()
        if not th.is_finite():
            raise OracleRefusal("domain", "exp(lntheta) not finite")
        J = self.J
        psiD = (th + J).digamma() - th.digamma()
        psiT = (acb(th + J).polygamma(1).real
                - acb(th).polygamma(1).real)
        return th, d, psiD, psiT

    def F(self, x):
        th, d, psiD, _ = self._pieces(x)
        J = self.J
        # NLL = ln(theta)_J - S lntheta - statics + (q terms, 0 at q=0)
        g_th = th * psiD - arb(self.S)
        g_d = d
        # dNLL/dq at q=0 = J(J-1)/2 - (theta+J-1) K_{J-1}/K_J
        g_q = arb(J) * (J - 1) / 2 - (th + (J - 1)) * self.K_Jm1 / self.K_J
        return [[g_d], [g_q]], [g_th]

    def H(self, x):
        th, d, psiD, psiT = self._pieces(x)
        nan = arb("nan")
        D0 = arb_mat([[arb(1)]])           # aux block
        D1 = arb_mat([[nan]])              # q row (masked out; fail-closed)
        B0 = arb_mat([[arb(0)]])
        B1 = arb_mat([[nan]])
        G = arb_mat([[th * psiD + th * th * psiT]])
        return [D0, D1], [B0, B1], G

    def nll_ewens(self, x):
        """Ewens NLL at q=0 (value gate; statics included)."""
        from collections import Counter
        zs, g = x
        th = g[0].exp()
        J, S = self.J, self.S
        statics = arb(J + 1).lgamma()
        for n in self.D:
            statics -= arb(n).log()
        for c in Counter(self.D).values():
            statics -= arb(c + 1).lgamma()
        lnp = (statics + S * th.log()
               - ((th + J).lgamma() - th.lgamma()))
        return -lnp

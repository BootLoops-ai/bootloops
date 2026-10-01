"""Single-sample Etienne oracle (2-dim): u = (lntheta, lnI).

The system is the NLL gradient (NLL = -lnP), so a CLINCH certificate is a
certified strict local MINIMUM of the negative log-likelihood == strict
local maximum of the likelihood. Derivatives come from the study's own
phase2_certify.lnP_derivs2 (exact term-wise digamma/trigamma recurrences,
centered accumulators) consumed BY IDENTITY; every arb op is ball-rigorous
over the input box, so oracle.H(box) is a rigorous enclosure.

Layout: flat theta = [lntheta, lnI]; border = [0] (lntheta, the market
doc's border dim 1), one latent block [1] (lnI).
"""
from flint import arb, acb

from . import ident
from .common import Part, cur_prec, mat1

from clinch.oracle_v31 import OracleRefusal


def lnP_derivs2_anchored(D, KP, th, I, prec, centers=None):
    """Adapter-owned variant of phase2_certify.lnP_derivs2 with the
    running term RE-ANCHORED at the A-distribution's center a0.

    WHY (measured on the J=235,360 fit): the reference loop anchors the
    running product t_A = I^A/(theta)_A at A=S, so over a BOX every ball
    picks up relative width ~ (A-S) * boxwidth by the time A ~ E[A]
    (2.3e5 multiplies); with H_lnI2 a 3-digit cancellation
    (VarA ~ 38,300 vs prefactor sum ~ 38,280, H ~ 20) the Hessian
    enclosure blows to O(1) relative width at radius 1e-9 sigma and
    Krawczyk cannot contract. Anchoring t'_A = I^(A-A0) (theta)_A0 /
    (theta)_A at A0 = round(a0) makes the accumulated width
    ~ |A-a0| * boxwidth ~ 1e3 * boxwidth on all significant terms (the
    common factor I^A0/(theta)_A0 cancels in every normalized moment).
    Identical mathematics, tighter enclosures; gated against the
    reference route at tight centers (agreement to the ball radius) and
    by the per-fit FD/value gates. Same signature/returns as the
    reference function."""
    from collections import Counter
    from flint import ctx
    ctx.prec = prec
    D = sorted(D)
    J, S = sum(D), len(D)
    S0, P = KP
    p2 = ident.pilot()["phase2_certify"]
    if centers is None:
        _, (a0, p0) = p2._prepass_centers(D, P, th, I, S, J)
    else:
        a0, p0 = centers
    A0 = min(max(int(round(a0)), S), J)
    # anchored psiD/psiT at A0 (exact finite sums; tight)
    psiD0 = arb(0)
    psiT0 = arb(0)
    for k in range(A0):
        d = 1 / (th + k)
        psiD0 += d
        psiT0 -= d * d
    F = arb(0)
    MX = arb(0); MXX = arb(0)
    MY = arb(0); MYY = arb(0); MXY = arb(0)
    MP = arb(0)

    def acc(A, t, psiD, psiT):
        nonlocal F, MX, MXX, MY, MYY, MXY, MP
        T = arb(P[A]) * t
        X = A - a0
        Y = psiD - p0
        F += T
        MX += X * T
        MXX += X * X * T
        MY += Y * T
        MYY += Y * Y * T
        MXY += X * Y * T
        MP += psiT * T

    # upward sweep A0 -> J
    t = arb(1)
    psiD = psiD0
    psiT = psiT0
    A = A0
    while True:
        acc(A, t, psiD, psiT)
        if A == J:
            break
        t *= I / (th + A)
        d = 1 / (th + A)
        psiD += d
        psiT -= d * d
        A += 1
    # downward sweep A0-1 -> S
    t = arb(1)
    psiD = psiD0
    psiT = psiT0
    A = A0
    while A > S:
        t *= (th + A - 1) / I
        d = 1 / (th + A - 1)
        psiD -= d
        psiT += d * d
        A -= 1
        acc(A, t, psiD, psiT)
    EX = MX / F; EY = MY / F
    VarA = MXX / F - EX * EX
    VarY = MYY / F - EY * EY
    Cov = MXY / F - EX * EY
    EA = EX + a0
    EpsiD = EY + p0
    EpsiT = MP / F
    pref = arb(J + 1).lgamma()
    for n in D:
        pref -= arb(n).log()
    for c in Counter(D).values():
        pref -= arb(c + 1).lgamma()
    pref += S * th.log()
    pref -= (I + J).lgamma() - I.lgamma()
    dI_pg = (I + J).digamma() - I.digamma()
    tri = (acb(I + J).polygamma(1).real - acb(I).polygamma(1).real)
    # lnF = ln(sum T') + A0 ln I - ln(theta)_A0; the last two pieces are
    # the common factor (needed only for the VALUE, not the moments)
    lnC = A0 * I.log()
    lnC -= (th + A0).lgamma() - th.lgamma()
    lnP = pref + F.log() + lnC
    g = [S - th * EpsiD, -I * dI_pg + EA]
    H = [[-th * EpsiD + th * th * (VarY - EpsiT), -th * Cov],
         [-th * Cov, VarA + (-I * dI_pg - I * I * tri)]]
    return lnP, g, H


class EtienneSSOracle:
    """oracle for one abundance vector D (single sample)."""

    def __init__(self, D, kernel_prec=208, tag="", anchored=False):
        p = ident.pilot()
        self._derivs = (lnP_derivs2_anchored if anchored
                        else p["phase2_certify"].lnP_derivs2)
        self.anchored = bool(anchored)
        self.D = sorted(int(x) for x in D)
        if min(self.D) < 1:
            raise OracleRefusal("data", "abundances must be >= 1")
        self.J, self.S = sum(self.D), len(self.D)
        self.tag = tag
        self.kernel_prec = int(kernel_prec)
        self.KP = p["ball_engine"].K_ball(self.D, self.kernel_prec)
        self.part = Part(blocks=[[1]], block_names=["local:lnI"],
                         border_idx=[0])
        self._centers = None

    @property
    def dims(self):
        return ([1], 1)

    def _eval(self, x, order2=True):
        zs, g = x
        u_th, u_I = g[0], zs[0][0]
        th, I = u_th.exp(), u_I.exp()
        if not (th.is_finite() and I.is_finite()):
            raise OracleRefusal("domain", "exp(u) not finite over the box")
        lnP, gr, H = self._derivs(self.D, self.KP, th, I, cur_prec(),
                                  centers=self._centers)
        return lnP, gr, H

    def nll(self, x):
        lnP, _, _ = self._eval(x)
        return -lnP

    def F(self, x):
        _, gr, _ = self._eval(x)
        # NLL gradient: blocks first (lnI), then border (lntheta)
        return [[-gr[1]]], [-gr[0]]

    def H(self, x):
        _, _, Hm = self._eval(x)
        D = [mat1(-Hm[1][1])]
        B = [mat1(-Hm[1][0])]
        G = mat1(-Hm[0][0])
        return D, B, G

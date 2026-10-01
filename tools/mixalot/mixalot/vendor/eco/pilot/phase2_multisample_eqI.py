"""Certified equal-I multi-sample ML (Panama 3-plot).

For equal I the kernel factorizes: M(D,A,I) = I^A * Mtilde(A) with Mtilde
parameter-independent (built once as a certified arb_poly with I=1).  The
likelihood then has single-sample structure
  lnP = pref(th,I) + ln sum_A Mtilde_A I^A/(th)_A
  pref = -ln prod Phi_vec! + sum_i [ln J_i! - sum_k ln n_ik!]
         - sum_i [lnGamma(I+J_i) - lnGamma(I)] + S ln th
so lnP_derivs2's machinery applies with the prefactor I-terms summed over
samples.  Interval Newton as in phase2_certify.
"""
import json
import math
import time
from collections import Counter

from flint import arb, acb, ctx

from multisample_ball import species_poly_ball, _tree
from phase2_certify import _prepass_centers

import phase2_certify as p2


def build_Mtilde(Dmat, prec):
    ctx.prec = prec
    ones = [arb(1)] * len(Dmat[0])
    polys = [species_poly_ball(row, ones, prec) for row in Dmat]
    return _tree(polys)


def lnP_derivs_ms(Dmat, MT, th, I, prec, centers=None):
    """Equal-I multi-sample lnP, grad, Hess in u=(ln th, ln I)."""
    ctx.prec = prec
    N = len(Dmat[0])
    S = len(Dmat)
    Js = [sum(r[i] for r in Dmat) for i in range(N)]
    Jtot = sum(Js)
    P = MT  # arb_poly, coefficient of z^A = Mtilde_A (A = S..Jtot)
    if centers is None:
        _, centers = _prepass_centers(None, P, th, I, S, Jtot)
    a0, p0 = centers
    F = arb(0)
    MX = arb(0); MXX = arb(0)
    MY = arb(0); MYY = arb(0); MXY = arb(0)
    MP = arb(0)
    t = I ** S
    for k in range(S):
        t /= th + k
    psiD = arb(0); psiT = arb(0)
    for k in range(S):
        psiD += 1 / (th + k)
        psiT -= 1 / (th + k) ** 2
    A = S
    while True:
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
        if A == Jtot:
            break
        t *= I / (th + A)
        psiD += 1 / (th + A)
        psiT -= 1 / (th + A) ** 2
        A += 1
    EX = MX / F; EY = MY / F
    VarA = MXX / F - EX * EX
    VarY = MYY / F - EY * EY
    Cov = MXY / F - EX * EY
    EA = EX + a0
    EpsiD = EY + p0
    EpsiT = MP / F
    # prefactor
    pref = arb(0)
    for cnt in Counter(tuple(r) for r in Dmat).values():
        pref -= arb(cnt + 1).lgamma()
    for i in range(N):
        pref += arb(Js[i] + 1).lgamma()
        pref -= (I + Js[i]).lgamma() - I.lgamma()
    for row in Dmat:
        for n in row:
            if n > 1:
                pref -= arb(n + 1).lgamma()
    pref += S * th.log()
    dI_pg = arb(0)
    tri = arb(0)
    for i in range(N):
        dI_pg += (I + Js[i]).digamma() - I.digamma()
        tri += (acb(I + Js[i]).polygamma(1).real - acb(I).polygamma(1).real)
    lnP = pref + F.log()
    g = [S - th * EpsiD, -I * dI_pg + EA]
    H = [[-th * EpsiD + th * th * (VarY - EpsiT), -th * Cov],
         [-th * Cov, VarA + (-I * dI_pg - I * I * tri)]]
    return lnP, g, H


def interval_newton_ms(Dmat, MT, th0, I0, prec=128, polish=8,
                       rads=(1e-9, 1e-8, 1e-7, 1e-6)):
    u = [math.log(th0), math.log(I0)]
    laststep = 1.0
    for it in range(polish):
        thM, IM = arb(u[0]).exp(), arb(u[1]).exp()
        _, g, H = lnP_derivs_ms(Dmat, MT, thM, IM, prec)
        det = H[0][0] * H[1][1] - H[0][1] * H[1][0]
        d0 = (H[1][1] * g[0] - H[0][1] * g[1]) / det
        d1 = (H[0][0] * g[1] - H[1][0] * g[0]) / det
        u = [u[0] - float(d0.mid()), u[1] - float(d1.mid())]
        laststep = max(abs(float(d0.mid())), abs(float(d1.mid())))
        if laststep < 1e-14:
            break
    for rad in rads:
        if rad < 4 * laststep:
            continue
        thB = (arb(u[0]) + arb(0, rad)).exp()
        IB = (arb(u[1]) + arb(0, rad)).exp()
        _, _, HB = lnP_derivs_ms(Dmat, MT, thB, IB, prec)
        detB = HB[0][0] * HB[1][1] - HB[0][1] * HB[1][0]
        if detB.contains(arb(0)):
            continue
        thM, IM = arb(u[0]).exp(), arb(u[1]).exp()
        _, gM, _ = lnP_derivs_ms(Dmat, MT, thM, IM, prec)
        d0 = (HB[1][1] * gM[0] - HB[0][1] * gM[1]) / detB
        d1 = (HB[0][0] * gM[1] - HB[1][0] * gM[0]) / detB
        n0 = arb(u[0]) - d0
        n1 = arb(u[1]) - d1
        inside = (abs(float(n0.mid()) - u[0]) + float(n0.rad()) < rad and
                  abs(float(n1.mid()) - u[1]) + float(n1.rad()) < rad)
        if not inside:
            continue
        is_max = (float(HB[0][0].mid()) + float(HB[0][0].rad()) < 0 and
                  float(detB.mid()) - float(detB.rad()) > 0)
        thI, II = n0.exp(), n1.exp()
        lnPI, _, _ = lnP_derivs_ms(Dmat, MT, thI, II, prec)
        return {"ok": True, "certified_max": bool(is_max), "box_rad": rad,
                "theta": str(thI), "I": str(II), "lnP": str(lnPI),
                "theta_mid": float(thI.mid()), "I_mid": float(II.mid()),
                "theta_rad": float(thI.rad()), "I_rad": float(II.rad())}
    return {"ok": False, "why": "no containment", "u": u, "laststep": laststep}


if __name__ == "__main__":
    rows = [tuple(int(x) for x in line.split())
            for line in open("../data/packages/panama3_full.txt")]
    t0 = time.time()
    prec = 160
    MT = build_Mtilde(rows, prec)
    res = interval_newton_ms(rows, MT, 259.3, 44.24, prec=prec)
    res["secs"] = round(time.time() - t0, 1)
    res["S"] = len(rows)
    res["Js"] = [sum(r[i] for r in rows) for i in range(3)]
    print(json.dumps(res, indent=1))
    with open("PHASE2_MS_EQI_PANAMA.json", "w") as f:
        json.dump(res, f, indent=1)

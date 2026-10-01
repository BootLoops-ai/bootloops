"""Phase 2: certified ML for the single-sample Etienne likelihood.

Rigorous enclosure of (theta-hat, I-hat) by interval Newton on the exact
gradient in u = (ln theta, ln I) coordinates:
  lnP = pref(theta,I) + ln F,  F = sum_A K_A I^A/(theta)_A  (K_A certified balls)
Term-wise derivatives with running digamma/trigamma updates (one O(J) pass):
  psiD_A = psi(theta+A)-psi(theta):   psiD_{A+1} = psiD_A + 1/(theta+A)
  psiT_A = psi'(theta+A)-psi'(theta): psiT_{A+1} = psiT_A - 1/(theta+A)^2
Interval Newton: N(X) = mid(X) - H(X)^{-1} g(mid(X)); N(X) inside X proves a
unique critical point in X (enclosed by N(X)); Hessian ball-negative-definite
proves it is a maximum.
"""
import json
import math
import time
from collections import Counter
from math import factorial

from flint import arb, arb_mat, ctx

from ball_engine import K_ball, species_poly_int  # K_ball -> (S, arb_poly)


def lnP_derivs2(D, KP, th, I, prec, centers=None):
    """lnP, grad, Hess in u=(ln th, ln I), centered accumulators.

    T-weighted moments of X=(A - a0) and Y=(psiD - p0) with float centers
    (a0, p0) avoid the E[A^2]-E[A]^2 cancellation catastrophe at J~2e4.
      d_lnth lnF   = -th E[psiD]
      d2_lnth lnF  = -th E[psiD] + th^2 (Var(psiD) - E[psiT])
      d_lnI lnF    = E[A];  d2_lnI lnF = Var(A)
      cross        = -th Cov(A, psiD)
    """
    ctx.prec = prec
    D = sorted(D)
    J, S = sum(D), len(D)
    S0, P = KP
    if centers is None:
        # cheap float prepass at midpoints for (a0, p0)
        _, (a0, p0) = _prepass_centers(D, P, th, I, S, J)
    else:
        a0, p0 = centers
    F = arb(0)
    MX = arb(0); MXX = arb(0)          # E-sums of X, X^2
    MY = arb(0); MYY = arb(0); MXY = arb(0)
    MP = arb(0)                         # E-sum of psiT
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
        if A == J:
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
    pref = arb(J + 1).lgamma()
    for n in D:
        pref -= arb(n).log()
    for c in Counter(D).values():
        pref -= arb(c + 1).lgamma()
    pref += S * th.log()
    pref -= (I + J).lgamma() - I.lgamma()
    dI_pg = (I + J).digamma() - I.digamma()
    from flint import acb
    tri = (acb(I + J).polygamma(1).real - acb(I).polygamma(1).real)
    lnP = pref + F.log()
    g = [S - th * EpsiD, -I * dI_pg + EA]
    H = [[-th * EpsiD + th * th * (VarY - EpsiT), -th * Cov],
         [-th * Cov, VarA + (-I * dI_pg - I * I * tri)]]
    return lnP, g, H


def _prepass_centers(D, P, th, I, S, J):
    """Float centers (a0, p0) = (E[A], E[psiD]) at ball midpoints, low prec."""
    import math
    thf = float(th.mid()) if hasattr(th, "mid") else float(th)
    If = float(I.mid()) if hasattr(I, "mid") else float(I)
    # log-space float pass
    lt = S * math.lgamma(1)  # 0
    lt = S * math.log(If)
    for k in range(S):
        lt -= math.log(thf + k)
    psiD = sum(1.0 / (thf + k) for k in range(S))
    mx = -1e300; num_a = 0.0; num_p = 0.0; den = 0.0
    logs = []
    A = S
    vals = []
    while True:
        # ln T_A = ln K_A + lt   (K_A from ball poly midpoint)
        lk = float(arb(P[A]).log().mid()) if float(arb(P[A]).mid()) > 0 else -1e300
        vals.append((A, psiD, lk + lt))
        if A == J:
            break
        lt += math.log(If) - math.log(thf + A)
        psiD += 1.0 / (thf + A)
        A += 1
    mx = max(v[2] for v in vals)
    for A, p, lv in vals:
        w = math.exp(lv - mx)
        den += w; num_a += A * w; num_p += p * w
    return None, (num_a / den, num_p / den)


def interval_newton(D, KP, th0, I0, prec=128, polish=6, rads=(1e-9, 1e-8, 1e-7, 1e-6)):
    """Certified enclosure of the critical point near (th0, I0).

    Phase 1: point-Newton polish (exact-midpoint evals) until |step| tiny.
    Phase 2: containment test on a rad ladder: N(X) = mid - H(X)^{-1} g(mid)
             inside X proves unique critical point; H(X) neg-def proves max.
    """
    ctx.prec = prec
    u = [math.log(th0), math.log(I0)]
    J = sum(D)
    laststep = 1.0
    for it in range(polish):
        thM, IM = arb(u[0]).exp(), arb(u[1]).exp()
        _, g, H = lnP_derivs2(D, KP, thM, IM, prec)
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
        _, _, HB = lnP_derivs2(D, KP, thB, IB, prec)
        detB = HB[0][0] * HB[1][1] - HB[0][1] * HB[1][0]
        if detB.contains(arb(0)):
            continue
        thM, IM = arb(u[0]).exp(), arb(u[1]).exp()
        _, gM, _ = lnP_derivs2(D, KP, thM, IM, prec)
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
        thI = n0.exp()
        II = n1.exp()
        mI = II / (II + (J - 1))
        lnPI, _, _ = lnP_derivs2(D, KP, thI, II, prec)
        return {"ok": True, "certified_max": bool(is_max), "box_rad": rad,
                "theta": str(thI), "I": str(II), "m": str(mI),
                "lnP": str(lnPI),
                "theta_mid": float(thI.mid()), "m_mid": float(mI.mid()),
                "theta_rad": float(thI.rad()), "m_rad": float(mI.rad())}
    return {"ok": False, "why": "no containment on rad ladder",
            "u": u, "laststep": laststep}


def float_opt(D, KP, prec=96, th0=50.0, I0=2000.0):
    """Crude NM in u-space on ball midpoints to seed the certifier."""
    def f(u):
        th = arb(u[0]).exp()
        I = arb(u[1]).exp()
        lnP, _, _ = lnP_derivs2(D, KP, th, I, prec)
        return -float(lnP.mid())
    import itertools
    u = [math.log(th0), math.log(I0)]
    step = 0.3
    fu = f(u)
    for _ in range(200):
        improved = False
        for dim, sgn in itertools.product((0, 1), (1, -1)):
            v = list(u)
            v[dim] += sgn * step
            fv = f(v)
            if fv < fu:
                u, fu = v, fv
                improved = True
        if not improved:
            step *= 0.35
            if step < 1e-9:
                break
    return math.exp(u[0]), math.exp(u[1]), -fu


if __name__ == "__main__":
    import sys
    path = sys.argv[1]
    D = sorted(int(x) for x in open(path).read().split())
    J, S = sum(D), len(D)
    t0 = time.time()
    KP = K_ball(D, 192)
    th0, I0, ln0 = float_opt(D, KP)
    res = interval_newton(D, KP, th0, I0, prec=192)
    res.update(J=J, S=S, file=path.split("/")[-1],
               secs=round(time.time() - t0, 1))
    print(json.dumps(res, indent=1))

"""Random-fission (EH2011 eq 42) certified likelihood engine.

P_rf[D] = (J;n) / ((I)_J prod_k S_k!) * sum_A K_rf(D,A) * w_A
  K_rf(D,A) = [z^A] prod_i (sum_a sbar(n_i,a) a! z^a)     (a! twist, NOT (a-1)!)
  w_A = I^A * theta_rf^{S-A} * BesselI_{A+S-1}(2 theta_rf) / BesselI_1(2 theta_rf)

Bessel ratios rho_nu = I_{nu+1}(x)/I_nu(x) by certified backward continued
fraction: rho_{nu-1} = 1/(2 nu/x + rho_nu), seeded with the rigorous interval
0 < rho_nu < min(1, x/(2(nu+1))) far above the needed range (contraction makes
the seed width irrelevant). All weights positive; running-product A-sum.

Point-mutation twin evaluated with the same code path shape via ball_engine.
Gates in gate_rf.py (fail-loud): mpmath ratio cross-check + EH2011 Table 3 printed LLs, all six forests.
"""
from math import factorial
from collections import Counter

from flint import arb, arb_poly, ctx

from ball_engine import species_poly_int, _product_tree as _tree


def species_poly_rf_ball(n, prec):
    """arb_poly with coefficients sbar(n,a) * a! (rf twist)."""
    ctx.prec = prec
    vi = species_poly_int(n)   # coeffs sbar(n,a)*(a-1)! for a>=1
    # convert: multiply coefficient a by a  (a! = a*(a-1)!)
    return arb_poly([arb(0)] + [arb(int(vi[a])) * a for a in range(1, n + 1)])


def K_rf_ball(D, prec):
    ctx.prec = prec
    return _tree([species_poly_rf_ball(n, prec) for n in D])


def bessel_ratio_ladder(x, nu_max, prec, buffer=400):
    """Certified rho_nu = I_{nu+1}(x)/I_nu(x) for nu = 1..nu_max (arb balls).

    Backward CF from nu_max+buffer with rigorous seed interval.
    """
    ctx.prec = prec
    x = arb(x) if not isinstance(x, arb) else x
    start = nu_max + buffer
    seed_hi = x / (2 * (start + 1))
    if float(seed_hi.mid()) > 1:
        seed_hi = arb(1)
    # rho as ball [0, seed_hi] -> arb(mid, rad)
    r = arb(float(seed_hi.mid()) / 2, float(seed_hi.mid()) / 2 + 1e-30)
    out = [None] * (nu_max + 1)
    for nu in range(start, 0, -1):
        # identity: rho_{nu-1} = 1/(2 nu/x + rho_nu)
        r = 1 / (2 * nu / x + r)
        if 1 <= nu - 1 <= nu_max:
            out[nu - 1] = r
    return out  # out[nu] = I_{nu+1}(x)/I_nu(x)


def logP_rf(D, theta_rf_str, m_str, prec, KP=None, ladder=None):
    """Certified ln P_rf[D | theta_rf, m]. Returns arb."""
    ctx.prec = prec
    D = sorted(D)
    J, S = sum(D), len(D)
    th = arb(theta_rf_str)
    mm = arb(m_str)
    I = mm * (J - 1) / (1 - mm)
    if KP is None:
        KP = K_rf_ball(D, prec)
    x = 2 * th
    if ladder is None:
        ladder = bessel_ratio_ladder(x, J + S, prec)
    # weights: w_A/w_{A-1} = (I/th) * rho_{A+S-2};  w_S = I^S * BesselI_{2S-1}/BesselI_1
    # build w_S/th^0: from nu=1: I_{2S-1}/I_1 = prod_{nu=1}^{2S-2} rho_nu
    w = I ** S
    for nu in range(1, 2 * S - 1):
        w *= ladder[nu]
    tot = arb(0)
    A = S
    while True:
        tot += arb(KP[A]) * w
        if A == J:
            break
        w *= (I / th) * ladder[A + S - 1]
        A += 1
    # prefactor: J!/prod n_i! / (I)_J / prod_k S_k!
    lp = arb(J + 1).lgamma()
    for n in D:
        lp -= arb(n + 1).lgamma()
    for c in Counter(D).values():
        lp -= arb(c + 1).lgamma()
    lp -= (I + J).lgamma() - I.lgamma()
    return lp + tot.log()

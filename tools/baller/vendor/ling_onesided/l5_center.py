# baller engine — first-order centered cell enclosures (the L5 lever).
#!/usr/bin/env python3
"""l5_center.py — L5 lever measurement: FIRST-ORDER CENTERED cell enclosures
(forward s1-duals through the double-collapse circuit + MVT), the SD side's
documented lever, measured on representative cells of the live-K32 and
wing-K24 panels.

Centered form for the panel integrand g(s1) = prior(s1) * sum_j Wj S(p_j,s1):
  int_cell g  in  h*g(m)  +/-  (h^2/4) * sup_cell |g'|
  g' = prior * (F' - lam*F),  F = sum Wj S,  F' via forward duals:
  only y1 = exp(-beta*s1) carries s1-dependence; dy1 = -beta*y1; the
  (A,B,C,D) closed forms differentiate by the product rule; the pow-tree
  and contraction propagate (value, derivative) pairs (3 poly-mults per
  pair-mult, ~3x one eval).
vs the 0th-order HULL form (production evalcell): 1 eval on the cell ball.
Measured: enclosure width (nats, relative to the point value) + secs, both
forms, at 3 cell widths; the scaling exponent decides the lever.
"""
import json, math, os, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from flint import arb, arb_poly, ctx
import double_unc as du
import gate_price as gp
from taylor_p import GROUPS, LAM


# ------------------------------------------------- dual (value, d/ds1) algebra
def _pmul(a, b):
    return (a[0] * b[0], a[1] * b[0] + a[0] * b[1])


def abcd_cat_dual(p, s1, groups):
    """du.abcd_cat closed forms with forward s1-duals. p: arb (no p-dual).
    Returns list of (n, (A,dA), (B,dB), (C,dC), (D,dD))."""
    one = p * 0 + 1
    pi0 = one - p
    beta = one / (2 * pi0 * p)
    y1 = (-beta * s1).exp()
    dy1 = -beta * y1
    # P1 entries and their s1-derivatives (linear in y1)
    P1 = [[(pi0 + p * y1, p * dy1), (p - p * y1, -p * dy1)],
          [(pi0 - pi0 * y1, -pi0 * dy1), (p + pi0 * y1, pi0 * dy1)]]
    pip = [pi0, p]
    q = pi0 * p
    out = []
    for (ab, xC, xD), n in groups:
        xA, xB = ab
        G0 = _pmul(P1[0][xA], P1[0][xB])
        G1 = _pmul(P1[1][xA], P1[1][xB])
        g = (pi0 * G0[0] + p * G1[0], pi0 * G0[1] + p * G1[1])
        dG = (G0[0] - G1[0], G0[1] - G1[1])
        sC = 1 if xC == 0 else -1
        sD = 1 if xD == 0 else -1
        A = (pip[xC] * pip[xD] * g[0], pip[xC] * pip[xD] * g[1])
        y1dG = _pmul((y1, dy1), dG)
        B = (q * y1dG[0] * sC * pip[xD], q * y1dG[1] * sC * pip[xD])
        y1g = _pmul((y1, dy1), g)
        inner = (pip[xC] * dG[0] + (y1g[0] if sC == 1 else -y1g[0]),
                 pip[xC] * dG[1] + (y1g[1] if sC == 1 else -y1g[1]))
        qy1 = (q * y1, q * dy1)
        Cd = _pmul(qy1, inner)
        C = (Cd[0] * sD, Cd[1] * sD)
        y2dG = _pmul(_pmul((y1, dy1), (y1, dy1)), dG)
        D = (q * (p - pi0) * y2dG[0] * (sC * sD),
             q * (p - pi0) * y2dG[1] * (sC * sD))
        out.append((n, A, B, C, D))
    return out


def _ppow_dual(q, n):
    """(q, dq)^n by binary powering with the product rule (poly pairs)."""
    r = None
    while n:
        if n & 1:
            r = q if r is None else _dmul(r, q)
        n >>= 1
        if n:
            q = _dmul(q, q)
    return r


def _dmul(a, b):
    return (a[0] * b[0], a[1] * b[0] + a[0] * b[1])


def _dtree(facs):
    while len(facs) > 1:
        facs.sort(key=lambda x: x[0].length())
        facs.append(_dmul(facs.pop(0), facs.pop(0)))
    return facs[0]


def S_double_dual(pn, s1):
    """(S, dS/ds1) through the phi-weighted accumulation (arb path)."""
    old = ctx.prec
    ctx.prec = pn.prec
    try:
        abcd = abcd_cat_dual(pn.p, s1, pn.groups)
        z = arb(0)
        R = None
        for tk, wk in zip(pn.taus, pn.phis):
            facs = []
            for (n, A, B, C, D) in abcd:
                cub = arb_poly([A[0], z, B[0] + C[0] * tk, D[0] * tk])
                dcub = arb_poly([A[1], z, B[1] + C[1] * tk, D[1] * tk])
                facs.append(_ppow_dual((cub, dcub), n))
            Q = _dtree(facs)
            Qw = (arb_poly([wk]) * Q[0], arb_poly([wk]) * Q[1])
            R = Qw if R is None else (R[0] + Qw[0], R[1] + Qw[1])
        return ((R[0] * pn.revW2)[pn.M], (R[1] * pn.revW2)[pn.M])
    finally:
        ctx.prec = old


# ------------------------------------------------------------ cell measurers
def hull_cell(pns, Wjs, aa, bb):
    ctx.prec = du.DC_QUAD
    t0 = time.time()
    sball = gp.hull(aa, bb)
    tot = arb(0)
    for pn, Wj in zip(pns, Wjs):
        tot += Wj * du.S_double(pn, sball)
    g = tot * du.prior(sball)
    h = arb(bb) - arb(aa)
    lo = gp.lower(g)
    if not bool(lo > 0):
        lo = arb(0)
    return lo * h, gp.upper(g) * h, time.time() - t0


def centered_cell(pns, Wjs, aa, bb):
    """h*g(m) +/- h^2/4 * sup|g'|; g' = prior*(F' - lam F) on the cell hull."""
    ctx.prec = du.DC_QUAD
    t0 = time.time()
    m = arb(aa) / 2 + arb(bb) / 2
    F = arb(0)
    for pn, Wj in zip(pns, Wjs):
        F += Wj * du.S_double(pn, m)
    gm = F * du.prior(m)                      # tight point value
    sball = gp.hull(aa, bb)
    Fh, dFh = arb(0), arb(0)
    for pn, Wj in zip(pns, Wjs):
        v, d = S_double_dual(pn, sball)
        Fh += Wj * v
        dFh += Wj * d
    gp_ball = du.prior(sball) * (dFh - LAM * Fh)
    supd = abs(arb(gp_ball.mid())) + arb(gp_ball.rad())
    h = arb(bb) - arb(aa)
    slop = h * h / 4 * supd
    lo = gp.lower(gm * h - slop)
    if not bool(lo > 0):
        lo = arb(0)
    return lo, gp.upper(gm * h + slop), time.time() - t0, gm * h


def nats(lo, hi, ref):
    """Enclosure width in nats relative to ref (arb, >0)."""
    if not bool(ref > 0):
        return None
    w = hi - lo
    if not bool(w > 0):
        return 0.0
    return float(arb((w / ref).log().mid()))


def main():
    rows = {}
    for tag, pm, pr, K in (("live_peak_K32", 0.346, 0.002, 32),
                           ("wing_lo_K24", 0.1675, 0.00295, 24)):
        ctx.prec = du.DC_PREC
        pjs, Wjs, Mfac = gp.panel_setup(pm, pr, K)
        pns = [du.PNode(pj, GROUPS, "cat", 3.5) for pj in pjs]
        # dual-vs-value agreement control at a point (FD check, coarse):
        eps = arb("1e-8")
        v0, d0 = S_double_dual(pns[0], arb("0.3"))
        vp = du.S_double(pns[0], arb("0.3") + eps)
        fd = (vp - v0) / eps
        agree = float(arb((abs(fd - d0) / abs(d0)).mid()))
        cells = []
        for (aa, bb) in ((0.25, 0.30), (0.25, 0.2625), (0.25, 0.253125)):
            hlo, hhi, ht = hull_cell(pns, Wjs, aa, bb)
            clo, chi, ct, cref = centered_cell(pns, Wjs, aa, bb)
            ref = cref  # tight midpoint-value ball as the scale
            cells.append(dict(
                cell=[aa, bb], h=round(bb - aa, 6),
                hull_relwidth_nats=(round(nats(hlo, hhi, ref), 3)
                                    if nats(hlo, hhi, ref) is not None else None),
                centered_relwidth_nats=(round(nats(clo, chi, ref), 3)
                                        if nats(clo, chi, ref) is not None else None),
                hull_lower_positive=bool(hlo > 0),
                centered_lower_positive=bool(clo > 0),
                hull_secs=round(ht, 1), centered_secs=round(ct, 1)))
            print(tag, cells[-1], flush=True)
        rows[tag] = dict(pm=pm, pr=pr, K=K, dual_fd_reldiff=agree,
                         cells=cells)
        res = os.environ.get("WPG_RESULTS",
                             os.path.join(os.getcwd(), "onesided_results"))
        os.makedirs(res, exist_ok=True)
        name = os.environ.get("L5_OUT", "L5_CENTERED.json")
        json.dump(rows, open(os.path.join(res, name), "w"), indent=1)
    print("saved", name)


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""Reusable certified 1-D quadrature: composite Gauss-Legendre (OWN summation
— mp.quad is banned here: its error estimate is absolute, and it silently
underconverges on exponential-ramp panels) with log-range panel refinement and
a two-degree self-check.

integrate(f, lnf, cuts, dps, band=12, degrees=(6,7)) -> (I, selfcons_digits)
  f    : integrand (mpf -> mpf), evaluated at working precision dps+PAD
  lnf  : cheap log|f| estimator for panel refinement (mpf -> mpf; return
         <= -1e9 where f ~ 0); a wrong lnf costs panels, not correctness —
         correctness comes from the degree pair + caller's oracle gates
  cuts : initial panel boundaries (must include domain endpoints)
"""
from mpmath import mp, mpf, fabs, log10

PAD = 15
_NEG = mpf('-1e9')


def gl_nodes(degree):
    from mpmath.calculus.quadrature import GaussLegendre
    return GaussLegendre(mp).calc_nodes(degree, mp.prec + 20)


def refine(lnf, cuts, dps, band=12):
    pts = sorted(set(cuts))
    sup = max([lnf((a + b) / 2) for a, b in zip(pts[:-1], pts[1:])] +
              [lnf(p) for p in pts])
    floor = sup - (dps + PAD + 25) * mp.log(10)
    out, stack = [], list(zip(pts[:-1], pts[1:]))
    while stack:
        a, b = stack.pop()
        m = (a + b) / 2
        vals = [v for v in (lnf(a), lnf(m), lnf(b)) if v > mpf('-0.9e9')]
        relevant = bool(vals) and max(vals) > floor
        wide = len(vals) < 3 or (max(vals) - min(vals)) > band
        if relevant and wide and (b - a) > mpf(10) ** (-(dps + PAD)):
            stack.append((a, m))
            stack.append((m, b))
        else:
            out.append((a, b))
    out.sort()
    return out, floor


def integrate(f, lnf, cuts, dps, band=12, degrees=(6, 7)):
    with mp.workdps(dps + PAD):
        panels, floor = refine(lnf, [mpf(c) for c in cuts], dps, band)
        vals = []
        for deg in degrees:
            nodes = gl_nodes(deg)
            tot = mpf(0)
            for a, b in panels:
                if lnf(a) < floor and lnf((a + b) / 2) < floor and lnf(b) < floor:
                    continue
                h = (b - a) / 2
                mid = (a + b) / 2
                s = mpf(0)
                for t, w in nodes:
                    s += w * f(mid + h * t)
                tot += h * s
            vals.append(tot)
        if vals[-1] == 0 or vals[0] == vals[-1]:
            sc = 9999.0
        else:
            sc = float(-log10(fabs((vals[0] - vals[-1]) / vals[-1])))
        return vals[-1], sc

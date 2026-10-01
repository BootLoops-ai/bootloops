# baller engine — one-sided evidence route (fixed needle grid + hull).
#!/usr/bin/env python3
# sd_onesided.py -- ONE-SIDED SD-evidence route probe (topology-0).
# One-sided pattern: the integrand is
# POSITIVE, so
#   LOWER  = sum over disjoint interior cells of  inf_cell(f) * vol   (certified,
#            inf via arb interval extension of the SAME pinned peeling circuit;
#            NO adaptivity, NO acceptance geometry), and
#   UPPER  = sum over the same cells of sup_cell(f) * vol  (hull enclosure)
#            + the already-cheap certified complement bounds REUSED VERBATIM from
#            sd_engine.py: per-axis tail sups (axis_tail_sup), corner cube
#            (corner_bound).
# Cells tile [0,S1]x[0,S2]x[0,S3] minus the corner cube [0,delta]^3 (grid lines at
# delta on every axis => every kept cell has one coordinate >= delta => U >= delta,
# K_2(2 sqrt(10U)) branch-safe, exactly the L-shape of the axis-aligned route).
# Grid is FIXED at design time: needle-concentrated bands from the (float,
# design-only) scan peak+widths, edge-refined cells toward u2,u3 -> 0,
# geometric wings out to the tail cuts S_k.
#
# Controls law: the 45-digit anchor + both plants MUST pass through THIS evaluation
# path (subcommand gate, prec 256) before any number counts. Gate output goes to
# ONESIDED_GATE.txt under the output root (SD_ONESIDED_OUT, default: the CURRENT
# working directory — never the vendor tree, whose unpinned files trip
# baller.verify()).
import json, time, sys, math, argparse, resource, os
from flint import arb, ctx

# Import root: this module's own directory (the vendored sd_engine sibling —
# inserting a foreign dir here would bind a DIFFERENT sd_engine, the recorded
# footgun). Output root: SD_ONESIDED_OUT, default the current working directory.
SDDIR = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.environ.get("SD_ONESIDED_OUT", os.getcwd())
sys.path.insert(0, SDDIR)
from sd_engine import (TOPOLOGIES, COUNTS, ALLPAT, PAT_1100, topo_name, G_cat,
                       G_bal, G_bal_z, lcore_int, ball, ub_of, axis_scan_peak,
                       axis_tail_sup, corner_bound)

ONE = arb(1)

def lb_of(x):
    l = x.mid() - x.rad()
    return l - l.rad()

def f_dn(x):
    """Directed float export of a certified LOWER endpoint (F-WIDTH cure
    2026-07-19): arb's directed lower bound, converted to double and
    nudged DOWN until certified <= the bound. Plain float() rounds to
    NEAREST and can overclaim."""
    l = x.lower()
    f = float(l)
    while arb(f) > l:
        f = math.nextafter(f, -math.inf)
    assert not (arb(f) > l), "directed lower export violated"
    return f

def f_up(x):
    """Directed float export of a certified UPPER endpoint (mirror of
    f_dn)."""
    u = x.upper()
    f = float(u)
    while arb(f) < u:
        f = math.nextafter(f, math.inf)
    assert not (arb(f) < u), "directed upper export violated"
    return f

def w_up(lo_f, hi_f):
    """Double upper bound on hi_f - lo_f (exact by Sterbenz for our
    nearly-equal same-sign endpoints; guarded outward regardless)."""
    import fractions
    w = hi_f - lo_f
    if fractions.Fraction(w) < (fractions.Fraction(hi_f)
                                - fractions.Fraction(lo_f)):
        w = math.nextafter(w, math.inf)
    return w

# ---------------- the NEW evaluation path (gated before numbers count) ----------------
def G_eval(topo, u1, u2, u3, one, region=None, plant_a=False):
    """kind-dispatched peeling circuit (banked sd_engine G_cat / G_bal)."""
    kind, perm = topo
    if kind == "cat":
        return G_cat(perm, u1, u2, u3, one, plant_a=plant_a)
    return G_bal(perm, u1, u2, u3, one, region)

def lcore_ball(topo, counts, N, u1, u2, u3, region=None, plant_a=False,
               plant_b=False):
    """arb ball of Lcore at u-balls (point use); counts int OR arb (N-scaled
    real-count rows, registration s4 precedent); plant_b drops 1100 from Gtot."""
    G = G_eval(topo, u1, u2, u3, ONE, region, plant_a)
    Gtot = None
    for D in ALLPAT:
        g = G.get(D)
        if g is None or (plant_b and D == PAT_1100):
            continue
        Gtot = g if Gtot is None else Gtot + g
    num = None
    for D, n in counts.items():
        f = G[D] ** n
        num = f if num is None else num * f
    return num / Gtot ** N

def lcore_cell(topo, counts, N, u1, u2, u3, region=None, plant_a=False,
               plant_b=False):
    """(lb, ub) certified bounds on Lcore over the u-ball cell, via POSITIVE-interval
    endpoint arithmetic: ub = prod ub(G_D)^n / lb(Gtot)^N, lb = prod lb(G_D)^n /
    ub(Gtot)^N. Avoids arb centered-ball powering, whose wide-cell denominator ball
    contains 0 (measured: nan on a 0.1-wide cell); endpoint chains carry their own
    rounding radii, closed by the final lb_of/ub_of. lb is clipped to 0 when any
    G_D endpoint touches 0 (true inf may be 0 there; one-sided validity kept).
    For 'bal' with a region tag: the REGION FORMULA's interval extension over the
    whole box (used by the straddle hull: each formula's extension contains its
    own region's values; the union hull covers the box)."""
    G = G_eval(topo, u1, u2, u3, ONE, region, plant_a)
    Gtot = None
    for D in ALLPAT:
        g = G.get(D)
        if g is None or (plant_b and D == PAT_1100):
            continue
        Gtot = g if Gtot is None else Gtot + g
    num_hi = None; num_lo = None; dead = False
    for D, n in counts.items():
        g = G[D]
        hb = ub_of(g)
        if not (hb > 0):
            return arb(0), arb(0)          # G_D <= 0 on the whole cell -> Lcore = 0
        hi = hb ** n
        num_hi = hi if num_hi is None else num_hi * hi
        if not dead:
            lo = lb_of(g)
            if lo > 0:
                lo = lo ** n
                num_lo = lo if num_lo is None else num_lo * lo
            else:
                dead = True
    gtlb = lb_of(Gtot)
    if not (gtlb > 0):
        gtlb = arb(1)                       # structural G_tot >= 1 (SD_FEASIBILITY 1.2)
    ub = ub_of(num_hi / gtlb ** N)
    if dead or num_lo is None:
        return arb(0), ub
    lb = lb_of(num_lo / ub_of(Gtot) ** N)
    if not (lb > 0):
        lb = arb(0)
    return lb, ub

def w_exact(U):
    c = 10 * U
    return arb(2000) / c * (2 * c.sqrt()).bessel_k(2)

# w(U) = 200/U * K_2(2 sqrt(10U)) is STRICTLY DECREASING on (0, inf): both factors
# are positive and strictly decreasing (K_nu(x) decreasing in x > 0, DLMF 10.37;
# 1/U decreasing). So w over [Ulo, Uhi] lies in [w(Uhi), w(Ulo)], bounded via a
# cached POINT lattice (smoke-2 measured defect: interval-argument bessel calls
# missed their cache on almost every cell -> ~2.1 ms/cell; point lattice caps the
# total bessel cost at #lattice-points regardless of cell count).
WQ = 0.0005
_WPT = {}

def w_pt(j):
    v = _WPT.get(j)
    if v is None:
        v = w_exact(arb(j) * arb(WQ))       # exact arb product: lattice point j*WQ
        _WPT[j] = v
    return v

def w_bounds(Ulo, Uhi):
    """(lb, ub) of w over [Ulo, Uhi], 0 < Ulo <= Uhi, by monotone lattice pinch.
    jlo overshoots floor by 1 and jhi overshoots ceil by 1 so float rounding in
    the index arithmetic can never break jlo*WQ <= Ulo <= Uhi <= jhi*WQ
    (penalty <= ~2*WQ*|dlnw/dU| ~ 0.006 nats, carried in the hull)."""
    jlo = int(Ulo / WQ) - 1
    jhi = int(Uhi / WQ) + 2
    if jlo < 1:
        jlo = 1
    lb = lb_of(w_pt(jhi)); ub = ub_of(w_pt(jlo))
    if not (lb > 0):
        lb = arb(0)
    return lb, ub

# ---------------- first-order (centered) cell enclosures via dual numbers ----------
# 0th-order endpoint slop measured ~1270 nats/unit on axis 1 (pure interval
# dependency at the peak, N-scale) => width <= 1 nat needs ~1e8 cells. The centered
# form kills it: ln Lcore(u) in ln Lcore(m) + sum_k g_k * (u_k - m_k) with g_k =
# d ln Lcore / du_k enclosed over the WHOLE cell (mean value theorem, cell convex,
# differentiability certified by G_D > 0 balls over the cell). g_k via forward-mode
# dual arithmetic through the SAME banked peeling circuit (G_cat is generic in its
# scalar type). The weight w never needs a derivative: monotone endpoints above.
class DU:
    __slots__ = ("v", "d")
    def __init__(s, v, d):
        s.v = v; s.d = d
    def exp(s):
        e = s.v.exp()
        return DU(e, (e * s.d[0], e * s.d[1], e * s.d[2]))
    def expm1(s):
        e = s.v.exp()
        return DU(s.v.expm1(), (e * s.d[0], e * s.d[1], e * s.d[2]))
    def __neg__(s):
        return DU(-s.v, (-s.d[0], -s.d[1], -s.d[2]))
    def __mul__(s, o):
        if isinstance(o, DU):
            return DU(s.v * o.v, (s.d[0] * o.v + s.v * o.d[0],
                                  s.d[1] * o.v + s.v * o.d[1],
                                  s.d[2] * o.v + s.v * o.d[2]))
        return DU(s.v * o, (s.d[0] * o, s.d[1] * o, s.d[2] * o))
    __rmul__ = __mul__
    def __add__(s, o):
        if isinstance(o, DU):
            return DU(s.v + o.v, (s.d[0] + o.d[0], s.d[1] + o.d[1], s.d[2] + o.d[2]))
        return DU(s.v + o, s.d)
    __radd__ = __add__
    def __sub__(s, o):
        if isinstance(o, DU):
            return DU(s.v - o.v, (s.d[0] - o.d[0], s.d[1] - o.d[1], s.d[2] - o.d[2]))
        return DU(s.v - o, s.d)
    def __rsub__(s, o):
        return DU(o - s.v, (-s.d[0], -s.d[1], -s.d[2]))

A0 = arb(0); A1 = arb(1)
DONE = DU(arb(1), (A0, A0, A0))
FBC = [0]   # 0th-order fallback cell counter (reset per run)
STR = [0]   # balanced straddle (v1=v2 kink) cell counter (reset per run)
PRS = [0]   # straddle cells resolved by the PRISM evaluator (reset per run)

def _prism_side(perm, counts, N, q1, q2, l3, h3, side):
    """centered (lb, ub) of the region formula over its prism half of the diagonal
    square [q1,q2]^2 x [l3,h3], via the kink-adapted variables (v, d, v3):
      side A (v1 >= v2): v = v2, d = v1-v2 in [0,h]; e12 = v3, e34 = d+v3
      side B (v2 >= v1): v = v1, d = v2-v1 in [0,h]; e12 = d+v3, e34 = v3
    The formula is ANALYTIC on the (v,d,v3) BOX (e-arguments stay >= 0), the box's
    affine image contains the prism, so the MVT centered enclosure bounds f over
    the prism with TRUE-gradient slop -- no 0th diagonal blowup (measured +52 nats
    on the balanced smoke). Returns None if positivity certification fails."""
    h = q2 - q1
    V = DU(ball(q1, q2), (A1, A0, A0))
    Dd = DU(ball(0.0, h), (A0, A1, A0))
    V3 = DU(ball(l3, h3), (A0, A0, A1))
    VD = V + Dd
    if side == "A":
        zc1 = (-VD).exp(); omc1 = -(-VD).expm1()
        zc2 = (-V).exp(); omc2 = -(-V).expm1()
        e12 = V3; e34 = Dd + V3
    else:
        zc1 = (-V).exp(); omc1 = -(-V).expm1()
        zc2 = (-VD).exp(); omc2 = -(-VD).expm1()
        e12 = Dd + V3; e34 = V3
    z12 = (-e12).exp(); om12 = -(-e12).expm1()
    z34 = (-e34).exp(); om34 = -(-e34).expm1()
    G = G_bal_z(perm, zc1, omc1, zc2, omc2, z12, om12, z34, om34, DONE)
    Gtot = None
    for D in ALLPAT:
        g = G.get(D)
        if g is None:
            continue
        Gtot = g if Gtot is None else Gtot + g
    if not (Gtot.v > 0):
        return None
    for D in counts:
        if not (G[D].v > 0):
            return None
    # midpoint INSIDE the prism: (v, d) = (q1 + 3h/8, h/4) -> v + d = q1 + 5h/8 <= q2
    mv = q1 + 0.375 * h; md = 0.25 * h; m3 = 0.5 * (l3 + h3)
    va = arb(mv); da = arb(md); v3a = arb(m3)
    if side == "A":
        pz1 = (-(va + da)).exp(); po1 = -(-(va + da)).expm1()
        pz2 = (-va).exp(); po2 = -(-va).expm1()
        pe12 = v3a; pe34 = da + v3a
    else:
        pz1 = (-va).exp(); po1 = -(-va).expm1()
        pz2 = (-(va + da)).exp(); po2 = -(-(va + da)).expm1()
        pe12 = da + v3a; pe34 = v3a
    Gp = G_bal_z(perm, pz1, po1, pz2, po2, (-pe12).exp(), -(-pe12).expm1(),
                 (-pe34).exp(), -(-pe34).expm1(), ONE)
    Gtp = None
    for D in ALLPAT:
        g = Gp.get(D)
        if g is None:
            continue
        Gtp = g if Gtp is None else Gtp + g
    num = None
    for D, n in counts.items():
        f = Gp[D] ** n
        num = f if num is None else num * f
    Lm = num / Gtp ** N
    if not (Lm > 0):
        return None
    inv_gt = A1 / Gtot.v
    g0 = (-N) * (Gtot.d[0] * inv_gt)
    g1 = (-N) * (Gtot.d[1] * inv_gt)
    g2 = (-N) * (Gtot.d[2] * inv_gt)
    for D, n in counts.items():
        gd = G[D]
        inv = A1 / gd.v
        g0 = g0 + n * (gd.d[0] * inv)
        g1 = g1 + n * (gd.d[1] * inv)
        g2 = g2 + n * (gd.d[2] * inv)
    lnr = Lm.log() + g0 * ball(q1 - mv, q2 - mv) + g1 * ball(-md, h - md) \
        + g2 * ball(l3 - m3, h3 - m3)
    lo = lb_of(lnr); hi = ub_of(lnr)
    return lb_of(lo.exp()), ub_of(hi.exp())

def cell_bounds(topo, counts, N, c1, c2, c3):
    """(f_lb, f_ub, tag) over the cell (f = Lcore * w): centered form when every
    needed G_D is certified positive on the cell, else 0th-order endpoint pairs.
    Cells are (lo, hi, uball, vol, mid, deltaball). Balanced topologies: cells are
    classified against the v1=v2 region kink (shared axis-1/axis-2 breakpoints make
    straddlers exactly the diagonal cells); pure-region cells go centered on that
    region's analytic formula, straddlers take the certified two-region 0th hull
    (each formula's extension over the box contains its region's values)."""
    lo1, hi1, u1, v1, m1, d1 = c1
    lo2, hi2, u2, v2, m2, d2 = c2
    lo3, hi3, u3, v3, m3, d3 = c3
    wlb, wub = w_bounds(lo1 + lo2 + lo3, hi1 + hi2 + hi3)
    region = None
    if topo[0] == "bal":
        if lo1 >= hi2:
            region = "A"
        elif hi1 <= lo2:
            region = "B"
        else:
            STR[0] += 1
            la, ua = lcore_cell(topo, counts, N, u1, u2, u3, "A")
            lb2, ub2 = lcore_cell(topo, counts, N, u1, u2, u3, "B")
            clo = la if la < lb2 else lb2
            chi = ua if ua > ub2 else ub2
            f_lb = lb_of(clo * wlb)
            if not (f_lb > 0):
                f_lb = arb(0)
            f_ub = ub_of(chi * wub)
            tag = "s"
            if abs(lo2 - lo1) < 1e-14 and abs(hi2 - hi1) < 1e-14:
                rA = _prism_side(topo[1], counts, N, lo1, hi1, lo3, hi3, "A")
                rB = _prism_side(topo[1], counts, N, lo1, hi1, lo3, hi3, "B")
                if rA is not None and rB is not None:
                    PRS[0] += 1
                    # per-unit-volume bounds: each prism is exactly HALF the box;
                    # elementwise BEST of prism and two-region 0th hull (both valid)
                    p_lb = lb_of((rA[0] + rB[0]) * wlb / 2)
                    p_ub = ub_of((rA[1] + rB[1]) * wub / 2)
                    if p_lb > f_lb:
                        f_lb = p_lb
                    if p_ub < f_ub:
                        f_ub = p_ub
                    tag = "p"
            return f_lb, f_ub, tag
    G = G_eval(topo, DU(u1, (A1, A0, A0)), DU(u2, (A0, A1, A0)),
               DU(u3, (A0, A0, A1)), DONE, region)
    Gtot = None
    for D in ALLPAT:
        g = G.get(D)
        if g is None:
            continue
        Gtot = g if Gtot is None else Gtot + g
    ok = Gtot.v > 0
    if ok:
        for D in counts:
            if not (G[D].v > 0):
                ok = False
                break
    if ok:
        Lm = lcore_ball(topo, counts, N, arb(m1), arb(m2), arb(m3), region)
        ok = Lm > 0
    if ok:
        inv_gt = A1 / Gtot.v
        g0 = (-N) * (Gtot.d[0] * inv_gt)
        g1 = (-N) * (Gtot.d[1] * inv_gt)
        g2 = (-N) * (Gtot.d[2] * inv_gt)
        for D, n in counts.items():
            gd = G[D]
            inv = A1 / gd.v
            g0 = g0 + n * (gd.d[0] * inv)
            g1 = g1 + n * (gd.d[1] * inv)
            g2 = g2 + n * (gd.d[2] * inv)
        lnr = Lm.log() + g0 * d1 + g1 * d2 + g2 * d3
        f_lb = lb_of(lb_of(lnr).exp() * wlb)
        f_ub = ub_of(ub_of(lnr).exp() * wub)
        if not (f_lb > 0):
            f_lb = arb(0)
        if ub_of(lnr) - lb_of(lnr) < 8:
            return f_lb, f_ub, "c"
        # wide cell: the centered spread is valid but useless (g-ball x wide delta
        # blew the smoke-3 hull to e^+6e5); the 0th-order endpoint pair is tighter
        # out here -- take the elementwise BEST of the two valid enclosures
        clo, chi = lcore_cell(topo, counts, N, u1, u2, u3, region)
        f_lb0 = lb_of(clo * wlb)
        f_ub0 = ub_of(chi * wub)
        if f_lb0 > f_lb:
            f_lb = f_lb0
        if f_ub0 < f_ub:
            f_ub = f_ub0
        return f_lb, f_ub, "cb"
    FBC[0] += 1
    clo, chi = lcore_cell(topo, counts, N, u1, u2, u3, region)
    f_lb = lb_of(clo * wlb)
    if not (f_lb > 0):
        f_lb = arb(0)
    return f_lb, ub_of(chi * wub), "f"

# ---------------- gate: 45-digit anchor + plants through THIS path ----------------
def gate():
    ctx.prec = 256
    lines = ["", "== ONESIDED GATE (sd_onesided.py NEW evaluation path: lcore_cell + "
             "cell_bounds(centered duals) + w_bounds; prec 256; run: "
             + time.strftime("%Y-%m-%d %H:%M:%S %Z") + ") =="]
    ok_all = True
    def check(name, ok, got):
        nonlocal ok_all
        lines.append(f"{'PASS' if ok else 'FAIL'}  {name}   got: {got}")
        ok_all = ok_all and ok
    T0 = TOPOLOGIES[0]
    u1 = (arb(10) / 9).log(); u2 = (arb(5) / 4).log(); u3 = (arb(10) / 7).log()
    lnC = arb(318).lgamma()
    for D, n in COUNTS.items():
        lnC = lnC - arb(n + 1).lgamma()
    ANCHOR = arb("-264.252750460743437038228550850745165801524406")
    va = lnC + lcore_ball(T0, COUNTS, 318, u1, u2, u3).log()
    check("banked circuit reproduces 45-digit anchor", abs(va - ANCHOR) < arb("1e-41"),
          va.str(45))
    # the PRODUCTION cell path (lcore_cell endpoint pairs) at the anchor point:
    # both endpoints must pin the anchor, plants must FIRE through the same path
    def pair_ln(pa=False, pb=False):
        lo, hi = lcore_cell(T0, COUNTS, 318, u1, u2, u3, plant_a=pa, plant_b=pb)
        return lnC + lo.log(), lnC + hi.log()
    plo, phi = pair_ln()
    check("cell path (lb,ub) both reproduce 45-digit anchor",
          abs(plo - ANCHOR) < arb("1e-41") and abs(phi - ANCHOR) < arb("1e-41"),
          f"[{plo.str(45)}, {phi.str(45)}]")
    alo, ahi = pair_ln(pa=True)
    check("cell path plant a FIRES (wrong C-edge -> -311.493...)",
          abs(alo - arb("-311.49316474472421895")) < arb("1e-9") and abs(alo - ANCHOR) > 1
          and abs(ahi - alo) < arb("1e-9"), alo.str(20))
    blo, bhi = pair_ln(pb=True)
    check("cell path plant b FIRES (1100 dropped -> -222.822...)",
          abs(blo - arb("-222.8224339223726654")) < arb("1e-9") and abs(blo - ANCHOR) > 1
          and abs(bhi - blo) < arb("1e-9"), blo.str(20))
    # weight path: monotone-lattice w_bounds must CONTAIN direct values at interval
    # interior points; direct value must match the banked acb route (sd_engine
    # integrand formula, IDENTICAL exact input c = 10*U in acb) and scipy float
    from flint import acb
    from scipy.special import kv
    okc = True; okx = True; oks = True
    for Uf in (0.06, 0.3, 0.6863, 1.7, 4.0, 12.0):
        wd = w_exact(arb(Uf))
        wlbp, wubp = w_bounds(Uf, Uf)
        okc = okc and (wlbp <= lb_of(wd)) and (wubp >= ub_of(wd))
        ca = 10 * acb(Uf)
        wa = (acb(2000) / ca * (2 * ca.sqrt()).bessel_k(2)).real
        okx = okx and bool(abs(wd - wa) < arb("1e-60"))
        wf = 2000.0 / (10 * Uf) * float(kv(2, 2 * math.sqrt(10 * Uf)))
        oks = oks and abs(float(wd.mid()) / wf - 1) < 1e-12
    for (Ua, Ub) in ((0.051, 0.09), (0.3, 0.72), (0.69, 0.701), (2.0, 6.5)):
        wlbp, wubp = w_bounds(Ua, Ub)
        for t in (0.0, 0.37, 0.71, 1.0):
            wd = w_exact(arb(Ua + (Ub - Ua) * t))
            okc = okc and (wlbp <= lb_of(wd)) and (wubp >= ub_of(wd))
    check("w_bounds contains direct values (6 pts + 4 intervals x 4)", okc, "containment")
    check("w arb == banked acb formula < 1e-60 (6 pts)", okx, "match")
    check("w vs scipy float < 1e-12 rel (6 pts)", oks, "match")
    # 0th-order cell-pair containment: interior point value inside the fat-cell (lb,ub)
    c1, c2, c3 = (0.60, 0.70), (0.010, 0.020), (0.0010, 0.0020)
    clo, chi = lcore_cell(T0, COUNTS, 318, ball(*c1), ball(*c2), ball(*c3))
    wlbp, wubp = w_bounds(c1[0] + c2[0] + c3[0], c1[1] + c2[1] + c3[1])
    flo = clo * wlbp; fhi = chi * wubp
    fpt = lcore_ball(T0, COUNTS, 318, arb(0.65), arb(0.015), arb(0.0015)) \
        * w_exact(arb(0.65 + 0.015 + 0.0015))
    check("0th-order cell bounds contain interior point value",
          (lb_of(flo) <= lb_of(fpt)) and (ub_of(fhi) >= ub_of(fpt)),
          f"cell ln[{float(flo.log().mid()) if flo > 0 else float('-inf'):.2f},"
          f"{float(fhi.log().mid()):.2f}] pt {float(fpt.mid().log().mid()):.2f}")
    # CENTERED path containment: certified point values inside cell_bounds on 3
    # cells (needle-core, mid-band, off-needle) x 9 deterministic interior points
    okcc = True; wids = []
    for (a1, b1, a2, b2, a3, b3) in (
            (0.680, 0.690, 0.0004, 0.0022, 0.00005, 0.00125),
            (0.640, 0.655, 0.010, 0.014, 0.004, 0.007),
            (0.60, 0.70, 0.010, 0.020, 0.0010, 0.0020)):
        cc1, cc2, cc3 = mkcell(a1, b1), mkcell(a2, b2), mkcell(a3, b3)
        flb, fub, _tag = cell_bounds(T0, COUNTS, 318, cc1, cc2, cc3)
        wids.append(round(float((fub / flb).log().mid()), 3) if flb > 0 else None)
        for (t1, t2, t3) in ((0.0, 0.5, 0.9), (0.5, 0.0, 0.5), (0.9, 0.9, 0.0),
                             (1.0, 0.3, 0.7), (0.3, 1.0, 1.0), (0.13, 0.55, 0.91),
                             (0.5, 0.5, 0.5), (1.0, 1.0, 1.0), (0.0, 0.0, 0.0)):
            p1 = a1 + (b1 - a1) * t1; p2 = a2 + (b2 - a2) * t2; p3 = a3 + (b3 - a3) * t3
            fp = lcore_ball(T0, COUNTS, 318, arb(p1), arb(p2), arb(p3)) \
                * w_exact(arb(p1) + arb(p2) + arb(p3))
            okcc = okcc and (flb <= lb_of(fp)) and (fub >= ub_of(fp))
    check("centered cell_bounds contain 27 certified interior values", okcc,
          f"cell widths {wids} nats")
    # dual-derivative wiring sanity: g_k on a near-point cell vs central differences
    # of the certified point evaluator (loose float tolerance; rigor lives in the
    # containment checks above)
    okg = True; eps = 1e-6
    pt = (0.65, 0.015, 0.0015)
    tiny = [mkcell(pt[k], pt[k] + 1e-9) for k in range(3)]
    Gd = G_cat(T0[1], DU(tiny[0][2], (A1, A0, A0)), DU(tiny[1][2], (A0, A1, A0)),
               DU(tiny[2][2], (A0, A0, A1)), DONE)
    Gtd = None
    for D in ALLPAT:
        gg = Gd.get(D)
        if gg is not None:
            Gtd = gg if Gtd is None else Gtd + gg
    for k in range(3):
        gk = (-318) * (Gtd.d[k] / Gtd.v)
        for D, n in COUNTS.items():
            gk = gk + n * (Gd[D].d[k] / Gd[D].v)
        pp = list(pt); pm = list(pt)
        pp[k] += eps; pm[k] -= eps
        fd = float(((lcore_ball(T0, COUNTS, 318, arb(pp[0]), arb(pp[1]), arb(pp[2])).log()
                     - lcore_ball(T0, COUNTS, 318, arb(pm[0]), arb(pm[1]), arb(pm[2])).log())
                    / (2 * eps)).mid())
        gm = float(gk.mid())
        okg = okg and abs(gm - fd) <= 1e-4 * (abs(fd) + 1)
    check("dual g_k match central differences (3 axes, 1e-4 rel)", okg,
          f"g=({float(gk.mid()):.4f} last axis)")
    # ---------------- BALANCED-PATH GATE (jg0-pattern degeneration + symmetry) ----
    from sd_engine import brute_G, _frac_to_arb
    from fractions import Fraction
    v1 = arb(2).log(); v2 = (arb(3) / 2).log(); v3g = (arb(5) / 4).log()
    permB = ("A", "B", "C", "D")
    # (b1) balanced DUAL-path values vs independent exact-Fraction brute force,
    #      region A and region B points, both balanced label sets used in prod
    okbf = True
    for idx in (12, 13, 14):
        tb = TOPOLOGIES[idx]
        for reg, vv, params in (
                ("A", (v1, v2, v3g),
                 (Fraction(1, 2), Fraction(2, 3), Fraction(4, 5), Fraction(3, 5))),
                ("B", (v2, v1, v3g),
                 (Fraction(2, 3), Fraction(1, 2), Fraction(3, 5), Fraction(4, 5)))):
            Gdual = G_eval(tb, DU(vv[0], (A1, A0, A0)), DU(vv[1], (A0, A1, A0)),
                           DU(vv[2], (A0, A0, A1)), DONE, reg)
            Gf = brute_G(tb, params)
            okbf = okbf and set(Gf) == set(Gdual)
            for D in Gf:
                okbf = okbf and bool(abs(Gdual[D].v - _frac_to_arb(Gf[D])) < arb("1e-70"))
    check("bal dual path == exact brute force (3 topos x A,B, 15 pats)", okbf, "match")
    # (b2) STAR degeneration (jg0 star/epsilon0 pattern): bal(h,h,0) == cat(h,0,0)
    h = arb(2).log()
    Gb = G_eval(TOPOLOGIES[12], h, h, arb(0), ONE, "A")
    Gc = G_eval(T0, h, arb(0), arb(0), ONE)
    okst = set(Gb) == set(Gc)
    for D in Gb:
        okst = okst and bool(abs(Gb[D] - Gc[D]) < arb("1e-60"))
    check("bal star limit == cat star limit (v3=0,u2=u3=0; 15 pats)", okst, "match")
    # (b3) v1=v2 KINK CONTINUITY: region formulas agree on the kink
    okk = True
    for tb in (TOPOLOGIES[12], TOPOLOGIES[13], TOPOLOGIES[14]):
        GA = G_eval(tb, h, h, v3g, ONE, "A")
        GB = G_eval(tb, h, h, v3g, ONE, "B")
        for D in GA:
            okk = okk and bool(abs(GA[D] - GB[D]) < arb("1e-60"))
    check("bal kink continuity A==B at v1=v2 (3 topos, 15 pats)", okk, "match")
    # (b4) cherry-exchange relabel symmetry + WRONG-REGION PLANT (jg0 m3 pattern)
    G1 = G_eval(("bal", ("A", "B", "C", "D")), v1, v2, v3g, ONE, "A")
    G2 = G_eval(("bal", ("C", "D", "A", "B")), v2, v1, v3g, ONE, "B")
    oksy = all(bool(abs(G1[D] - G2[D]) < arb("1e-60")) for D in G1)
    check("bal cherry-exchange symmetry (perm+v1<->v2+region swap)", oksy, "match")
    G2w = G_eval(("bal", ("C", "D", "A", "B")), v2, v1, v3g, ONE, "A")
    nfire = sum(1 for D in G1 if bool(abs(G1[D] - G2w[D]) > arb("1e-6")))
    check("bal wrong-region PLANT FIRES (>=1 pattern differs)", nfire >= 1,
          f"{nfire}/15 patterns differ")
    # (b5) balanced centered containment: A-box, B-box, straddle box x 9 points
    tb = TOPOLOGIES[12]
    okbc = True; bwids = []
    for (a1, b1, a2, b2, a3, b3) in ((0.50, 0.52, 0.30, 0.32, 0.10, 0.12),
                                     (0.30, 0.32, 0.50, 0.52, 0.10, 0.12),
                                     (0.40, 0.50, 0.40, 0.50, 0.10, 0.12)):
        cc1, cc2, cc3 = mkcell(a1, b1), mkcell(a2, b2), mkcell(a3, b3)
        flb, fub, _t = cell_bounds(tb, COUNTS, 318, cc1, cc2, cc3)
        bwids.append(round(float((fub / flb).log().mid()), 3) if flb > 0 else None)
        for (t1, t2, t3) in ((0.0, 0.5, 0.9), (0.5, 0.0, 0.5), (0.9, 0.9, 0.0),
                             (1.0, 0.3, 0.7), (0.3, 1.0, 1.0), (0.13, 0.55, 0.91),
                             (0.5, 0.5, 0.5), (1.0, 1.0, 1.0), (0.0, 0.0, 0.0)):
            p1 = a1 + (b1 - a1) * t1; p2 = a2 + (b2 - a2) * t2; p3 = a3 + (b3 - a3) * t3
            regp = "A" if p1 >= p2 else "B"
            fp = lcore_ball(tb, COUNTS, 318, arb(p1), arb(p2), arb(p3), regp) \
                * w_exact(arb(p1) + arb(p2) + arb(p3))
            okbc = okbc and (flb <= lb_of(fp)) and (fub >= ub_of(fp))
    check("bal cell_bounds contain 27 certified interior values (A/B/straddle)",
          okbc, f"cell widths {bwids} nats")
    # (b6) balanced dual g_k vs central differences (region-A interior point)
    okbg = True
    ptb = (0.50, 0.30, 0.11)
    tinyb = [mkcell(ptb[k], ptb[k] + 1e-9) for k in range(3)]
    Gdb = G_eval(tb, DU(tinyb[0][2], (A1, A0, A0)), DU(tinyb[1][2], (A0, A1, A0)),
                 DU(tinyb[2][2], (A0, A0, A1)), DONE, "A")
    Gtb = None
    for D in ALLPAT:
        gg = Gdb.get(D)
        if gg is not None:
            Gtb = gg if Gtb is None else Gtb + gg
    for k in range(3):
        gk = (-318) * (Gtb.d[k] / Gtb.v)
        for D, n in COUNTS.items():
            gk = gk + n * (Gdb[D].d[k] / Gdb[D].v)
        pp = list(ptb); pm = list(ptb)
        pp[k] += eps; pm[k] -= eps
        fd = float(((lcore_ball(tb, COUNTS, 318, arb(pp[0]), arb(pp[1]), arb(pp[2]), "A").log()
                     - lcore_ball(tb, COUNTS, 318, arb(pm[0]), arb(pm[1]), arb(pm[2]), "A").log())
                    / (2 * eps)).mid())
        okbg = okbg and abs(float(gk.mid()) - fd) <= 1e-4 * (abs(fd) + 1)
    check("bal dual g_k match central differences (3 axes, 1e-4 rel)", okbg, "match")
    # (b6b) PRISM straddle evaluator: per-side containment of certified point
    # values on 2 diagonal cells (needle-adjacent + generic), and the prism width
    # must beat the 0th two-region hull it replaces
    okpc = True; okpw = True
    for (q1, q2, l3, h3) in ((0.68, 0.69, 0.0001, 0.0021), (0.44, 0.46, 0.10, 0.12)):
        w0lb, w0ub = None, None
        la, ua = lcore_cell(tb, COUNTS, 318, ball(q1, q2), ball(q1, q2),
                            ball(l3, h3), "A")
        lb2b, ub2b = lcore_cell(tb, COUNTS, 318, ball(q1, q2), ball(q1, q2),
                                ball(l3, h3), "B")
        oth_lo = la if la < lb2b else lb2b
        oth_hi = ua if ua > ub2b else ub2b
        for side in ("A", "B"):
            r = _prism_side(tb[1], COUNTS, 318, q1, q2, l3, h3, side)
            okpc = okpc and (r is not None)
            if r is None:
                continue
            plb, pub = r
            for (tv, td, t3) in ((0.0, 0.0, 0.5), (0.5, 0.5, 0.0), (0.9, 0.3, 1.0),
                                 (0.2, 0.9, 0.5), (1.0, 0.0, 0.2)):
                v = q1 + tv * (q2 - q1); d = td * (q2 - v); v3p = l3 + t3 * (h3 - l3)
                v1p, v2p = (v + d, v) if side == "A" else (v, v + d)
                fp = lcore_ball(tb, COUNTS, 318, arb(v1p), arb(v2p), arb(v3p), side)
                okpc = okpc and (plb <= lb_of(fp)) and (pub >= ub_of(fp))
            wid_p = float((pub / plb).log().mid()) if plb > 0 else 1e9
            wid_0 = float((oth_hi / oth_lo).log().mid()) if oth_lo > 0 else 1e9
            okpw = okpw and (wid_p < wid_0 or wid_0 > 1e8)
    check("prism sides contain certified prism-interior values (2 cells x 2 x 5)",
          okpc, "containment")
    check("prism width beats 0th two-region hull on diagonal cells", okpw, "tighter")
    # (b7) REAL-COUNTS path at fscale=1 must reproduce the 45-digit anchor + plants
    fs1 = arb(318) / 318
    counts1 = {D: arb(n) * fs1 for D, n in COUNTS.items()}
    var = lnC + lcore_ball(T0, counts1, 318, u1, u2, u3).log()
    vra = lnC + lcore_ball(T0, counts1, 318, u1, u2, u3, plant_a=True).log()
    vrb = lnC + lcore_ball(T0, counts1, 318, u1, u2, u3, plant_b=True).log()
    check("real-counts path reproduces 45-digit anchor (fscale=1)",
          abs(var - ANCHOR) < arb("1e-41"), var.str(45))
    check("real-counts path plants a+b FIRE",
          abs(vra - arb("-311.49316474472421895")) < arb("1e-9")
          and abs(vrb - arb("-222.8224339223726654")) < arb("1e-9"),
          f"[{vra.str(16)}, {vrb.str(16)}]")
    lines.append(f"ONESIDED GATE OVERALL: {'PASS' if ok_all else 'FAIL - NOTHING BANKS'}")
    out = "\n".join(lines)
    print(out, flush=True)
    with open(f"{OUTDIR}/ONESIDED_GATE.txt", "a") as f:
        f.write(out + "\n")
    return 0 if ok_all else 1

# ---------------- fixed grid design ----------------
def make_lnF(topo, counts, N):
    """float log-integrand for GRID DESIGN ONLY (never banked; same law as the
    engine's axis_scan_peak)."""
    from scipy.special import kv
    kind, perm = topo
    from sd_engine import F
    def lnF(u1, u2, u3):
        if u1 <= 0 or u2 <= 0 or u3 <= 0:
            return -1e300
        U = u1 + u2 + u3
        c = 10.0 * U
        k = kv(2, 2.0 * math.sqrt(c))
        if not (k > 0):
            return -1e300
        s = math.log(2000.0) - math.log(c) + math.log(k)
        if kind == "cat":
            G = G_cat(perm, F(u1), F(u2), F(u3), F(1.0))
        else:
            G = G_bal(perm, F(u1), F(u2), F(u3), F(1.0), "A" if u1 >= u2 else "B")
        gt = sum(g.v for g in G.values())
        s += -N * math.log(gt)
        for D, n in counts.items():
            gv = G[D].v
            if gv <= 0:
                return -1e300
            s += n * math.log(gv)
        return s
    return lnF

DGAP = 28.0          # wing-cell hull margin (nats below peak per wing cell, by design)

def ridge_drop(lnF, vpk, k, x, oth):
    """design-only: drop of the RIDGE max of lnF below the peak on the slab u_k = x,
    maximizing over the other two axes by warm-started multiplicative descent.
    Returns (drop, updated warm-start)."""
    p = list(oth); p[k] = x
    best = lnF(*p)
    for _ in range(3):
        moved = False
        for a in range(3):
            if a == k:
                continue
            for f in (0.5, 0.7, 1.4, 2.0):
                q = list(p); q[a] = max(p[a] * f, 1e-7)
                v = lnF(*q)
                if v > best:
                    best = v; p = q; moved = True
        if not moved:
            break
    return max(vpk - best, 0.0), p

def axis_breaks(k, upk, ws, S, delta, h0, ck, lnF, vpk, probe_ub, bal_diag=False,
                m_mid=2.0, m_out=3.5):
    """FIXED per-axis breakpoints: inner band p+-1.5w step h0; mid p+-3.5w step 2.5h0;
    outer p+-5.0w step 6h0; PROFILE-SIZED wings to 0 and S -- wing cell width
    (ridge drop - DGAP)/ck so every wing cell's endpoint-slop hull stays ~e^-DGAP
    below peak (measured smoke defect: geometric wings x slope ~1300/unit put the
    hull +170 nats over); delta inserted; edge refinement toward 0."""
    p = upk[k]; w = ws[k]
    pts = set([0.0, S, delta])
    lo_in, hi_in = p - 1.5 * w, p + 1.5 * w
    lo_mid, hi_mid = p - 3.5 * w, p + 3.5 * w
    lo_out, hi_out = p - 5.0 * w, p + 5.0 * w
    def fill(a, b, h):
        a = max(a, 0.0); b = min(b, S)
        if b <= a:
            return
        n = max(1, int(math.ceil((b - a) / h)))
        for i in range(n + 1):
            pts.add(a + (b - a) * i / n)
    fill(lo_in, hi_in, h0)
    mb = min(m_mid, 2.0) * h0 if bal_diag else m_mid * h0
    ob = min(m_out, 2.0) * h0 if bal_diag else m_out * h0
    fill(lo_mid, lo_in, mb); fill(hi_in, hi_mid, mb)
    fill(lo_out, lo_mid, ob); fill(hi_mid, hi_out, ob)
    # wings: candidate width from the ridge-drop rule, then a design-time CERTIFIED
    # shrink loop -- the smoke-4 hot cells showed the peak slop constant ck is wrong
    # near u1 -> 0 (local slop ~1/u1 from the (1-y1) patterns): shrink wc until the
    # probe cell's certified 0th-order ub sits DGAP below the peak
    def sized(a_from, direction, room):
        drop, oth_l = ridge_drop(lnF, vpk, k, a_from, sized.oth)
        sized.oth = oth_l
        wc = min(max((drop - DGAP) / ck, 3.5 * h0), 0.6, room)
        for _ in range(12):
            if wc <= 0.5 * h0 * 1.01:
                break
            a, b = (a_from, a_from + wc) if direction > 0 else (a_from - wc, a_from)
            pu = probe_ub(k, a, b, oth_l)
            if pu is not None and pu <= vpk - DGAP:
                break
            wc = max(wc * 0.5, 0.5 * h0)   # floor 0.5*h0 (3.5*h0 left balanced
        return wc                          # diagonal squares stuck-hot, smoke-5)
    sized.oth = list(upk)
    x = min(hi_out, S)                                   # upper wing
    while x < S - 1e-12:
        wc = sized(x, +1, S - x)
        x = x + max(wc, 1e-5); pts.add(min(x, S))
    sized.oth = list(upk)
    x = max(lo_out, 0.0)                                 # lower wing
    while x > 1e-12:
        wc = sized(x, -1, x)
        x = x - max(wc, 1e-5); pts.add(max(x, 0.0))
    if lo_out <= 0.0:                                    # edge refinement toward 0
        for f in (2, 4, 8):
            pts.add(min(h0, S) / f)
    srt = sorted(q for q in pts if 0.0 <= q <= S)
    out = [srt[0]]                                       # dedupe (min gap h0/16)
    for q in srt[1:]:
        if q - out[-1] > h0 / 16:
            out.append(q)
    if out[-1] < S:
        out.append(S)
    return out

def dball(a, b, m):
    """arb ball containing [a-m, b-m] built in arb arithmetic (float pre-subtraction
    could be 1 ulp narrow; arb ops carry rounding radii)."""
    dlo = arb(a) - arb(m); dhi = arb(b) - arb(m)
    return (dlo + dhi) / 2 + ((dhi - dlo) / 2) * arb(0, 1)

def mkcell(a, b):
    m = 0.5 * (a + b)
    return (a, b, ball(a, b), arb(b) - arb(a), m, dball(a, b, m))

def build_axis(breaks):
    return [mkcell(a, b) for a, b in zip(breaks[:-1], breaks[1:])]

# ---------------- the probe run ----------------
def run(topo_i, prec, hfrac, budget_s, out_path, scale_fix=0.0, chunk="0/1",
        Nt=318):
    ctx.prec = prec
    topo = TOPOLOGIES[topo_i]
    kind, perm = topo
    if Nt == 318:
        counts = dict(COUNTS)               # float/int shadow: scan + sup bounds
        counts_eval = counts                # certified evaluation counts
        N = 318
    else:
        # N-scaled denominator rows: pinned caterpillar only, real-valued counts
        # rescaled to the same shape (registration + engine convention)
        assert topo == TOPOLOGIES[0], "N-scaled rows registered on pinned caterpillar only"
        counts = {D: n * Nt / 318.0 for D, n in COUNTS.items()}
        fs = arb(Nt) / 318
        counts_eval = {D: arb(n) * fs for D, n in COUNTS.items()}
        N = Nt
    delta = 0.05
    t0 = time.time(); c0 = resource.getrusage(resource.RUSAGE_SELF)
    def cpu():
        r = resource.getrusage(resource.RUSAGE_SELF)
        return r.ru_utime + r.ru_stime - c0.ru_utime - c0.ru_stime

    # (a) float needle scan -- never banked, grid design only (same law as engine)
    vpk, upk, ws, lnZ_pre = axis_scan_peak(topo, counts, N)
    t_scan = cpu()
    print(f"[scan] peak lnF={vpk:.2f} u*=({upk[0]:.4f},{upk[1]:.4f},{upk[2]:.4f}) "
          f"widths=({ws[0]:.4f},{ws[1]:.4f},{ws[2]:.4f}) lnZ_pre~{lnZ_pre:.1f} "
          f"cpu={t_scan:.1f}s", flush=True)

    # (b) certified complement bounds -- REUSED sd_engine functions
    tgt_ln = lnZ_pre - 12.0
    S = []; tail_ubs = []; tail_lns = []
    for k in range(3):
        T = max(1.2, upk[k] * 3.0 + 0.5)
        while True:
            sb = axis_tail_sup(topo, k, T, counts)
            sl = float(sb.log().mid()) if sb > 0 else -1e300
            if sl < tgt_ln or T > 700:
                break
            T *= 1.6
        if sl >= tgt_ln:
            raise RuntimeError(f"axis-{k} tail will not close")
        S.append(T); tail_ubs.append(sb); tail_lns.append(sl)
    while True:
        cb = corner_bound(topo, delta, counts)
        cb_ln = float(cb.log().mid()) if cb > 0 else -1e300
        if cb_ln < tgt_ln or delta <= 0.0035:
            break
        delta /= 2                # t04/t05: corner sup at delta=0.05 sat ABOVE the
    if cb_ln >= tgt_ln:           # integral (needle nearer the origin) -- shrink
        raise RuntimeError(f"corner bound will not close (delta={delta}, ln={cb_ln:.1f})")
    t_bounds = cpu() - t_scan
    print(f"[bounds] S=({S[0]:.2f},{S[1]:.2f},{S[2]:.2f}) tail_lns="
          f"{[round(x,1) for x in tail_lns]} corner_ln<={cb_ln:.1f} "
          f"(vs lnZ_pre {lnZ_pre:.1f}) cpu={t_bounds:.1f}s", flush=True)

    def reg_of(a1, b1, a2, b2):
        """design-only region pick for thin probe cells (bal): pure-region tag if
        certain, else by midpoints (mis-pick costs grid quality, never rigor)."""
        if kind == "cat":
            return None
        if a1 >= b2:
            return "A"
        if b1 <= a2:
            return "B"
        return "A" if (a1 + b1) >= (a2 + b2) else "B"

    # (c) inline calibration: 0th-order slop slopes c_k (wing sizing) + CPU-measured
    # per-cell cost of the centered evaluator
    cks = []
    for k in range(3):
        h = ws[k] / 8
        ck = None
        while True:
            cs = []
            for ax in range(3):
                if ax == k:
                    lo = max(0.0, upk[ax] - h / 2); cs.append((lo, lo + h))
                else:
                    cs.append((upk[ax], upk[ax] + 1e-9))
            lb, ub = lcore_cell(topo, counts_eval, N,
                                ball(*cs[0]), ball(*cs[1]), ball(*cs[2]),
                                reg_of(cs[0][0], cs[0][1], cs[1][0], cs[1][1]))
            if lb > 0:
                wid = float((ub / lb).log().mid())
                if wid < 1.5:
                    ck = wid / h
                    break
            h /= 2
            if h < ws[k] / 4096:
                raise RuntimeError(f"calibration failed on axis {k}")
        cks.append(max(ck, 1e-3))
    cal_cells = [mkcell(upk[0], upk[0] + ws[0] / 7),
                 mkcell(upk[1], upk[1] + ws[1] / 7),
                 mkcell(max(0.0, upk[2] - ws[2] / 14), upk[2] + ws[2] / 14)]
    cc0 = resource.getrusage(resource.RUSAGE_SELF)
    ncal = 300
    for i in range(ncal):
        flb, fub, _tag = cell_bounds(topo, counts_eval, N, *cal_cells)
    cc1 = resource.getrusage(resource.RUSAGE_SELF)
    t_eval = (cc1.ru_utime + cc1.ru_stime - cc0.ru_utime - cc0.ru_stime) / ncal
    wid0 = float((fub / flb).log().mid()) if flb > 0 else float("inf")
    # inner-band recalibration (t10: large-absolute-width needles leave the INNER
    # band h2-slop-dominated -- band mults only ever tightened mid/outer): shrink
    # h0 so the peak cell certifies ~0.8 nats at scale 1 (h^2 regime => sqrt)
    hadj = min(1.0, math.sqrt(0.8 / wid0)) if wid0 > 0.8 else 1.0
    print(f"[calib] 0th slopes c=({cks[0]:.0f},{cks[1]:.0f},{cks[2]:.0f})/unit "
          f"centered t_eval={t_eval*1000:.3f}ms (CPU) peak-cell wid={wid0:.4f} nats "
          f"hadj={hadj:.3f}", flush=True)

    # (d) fixed grid design under the CPU budget (design-time choice, then frozen)
    lnF = make_lnF(topo, counts, N)
    def probe_ub(k, a, b, oth):
        """float ln of the certified 0th-order ub of Lcore*w on the probe cell
        [a,b] along axis k (design-time wing sizing only; the certified pass
        re-bounds every real cell). For BALANCED axes 0/1 the probe is the
        DIAGONAL SQUARE [a,b]^2, worst case over both region formulas -- the
        smoke-3 balanced hot cells were diagonal wing squares whose slop budget
        is spent on BOTH axes while the ridge drop is shared (non-additive)."""
        square = (kind == "bal" and k in (0, 1))
        oax = [ax for ax in range(3) if ax != k and not (square and ax in (0, 1))]
        # positions on the other axes: ridge point AND core-cell-sized boxes at the
        # needle / off-center -- t06 hot cells were wing-1 x off-center-core-u2
        # crosses the thin ridge-pinned probe never saw (local slop position-varies)
        combos = [tuple(("pt", max(oth[ax], 1e-9)) for ax in oax)]
        combos.append(tuple(("cell", upk[ax]) for ax in oax))
        combos.append(tuple(("cell", max(upk[ax] - ws[ax], 0.25 * ws[ax]))
                            for ax in oax))
        worst = None
        for combo in combos:
            cs = [None] * 3; slo = 0.0; shi = 0.0
            for ax in range(3):
                if ax == k or (square and ax in (0, 1)):
                    cs[ax] = ball(max(a, 0.0), b); slo += max(a, 0.0); shi += b
            for (mode, v), ax in zip(combo, oax):
                if mode == "pt":
                    cs[ax] = ball(v, v * (1 + 1e-9))
                    slo += v; shi += v * (1 + 1e-9)
                else:
                    h2 = ws[ax] * hfrac / 2
                    lo2p = max(v - h2, 0.0); hi2p = v + h2
                    cs[ax] = ball(lo2p, hi2p); slo += lo2p; shi += hi2p
            wlbp, wubp = w_bounds(max(slo, 1e-6), shi)
            if square:
                ua = lcore_cell(topo, counts_eval, N, *cs, "A")[1]
                ub2 = lcore_cell(topo, counts_eval, N, *cs, "B")[1]
                ub = ua if ua > ub2 else ub2
            else:
                lb, ub = lcore_cell(topo, counts_eval, N, *cs,
                                    reg_of(float(cs[0].mid() - cs[0].rad()),
                                           float(cs[0].mid() + cs[0].rad()),
                                           float(cs[1].mid() - cs[1].rad()),
                                           float(cs[1].mid() + cs[1].rad())))
            fu = ub * wubp
            if not (fu > 0):
                continue
            v = float(fu.log().mid())
            if worst is None or v > worst:
                worst = v
        return worst
    def band_mult(k, pos_sig, cands, tgt, h0k):
        """largest band step multiplier whose test cell (at pos_sig needle-widths
        from the peak, other axes = peak core cells) certifies width <= tgt nats
        through the PRODUCTION cell_bounds path; both sides tested."""
        best = cands[-1]
        for m in cands:
            ok = True
            for sgn in (+1, -1):
                pos = upk[k] + sgn * pos_sig * ws[k]
                h = m * h0k
                lo = max(0.0, min(pos, S[k] - h))
                cs = []
                for ax in range(3):
                    if ax == k:
                        cs.append(mkcell(lo, lo + h))
                    else:
                        c0 = max(upk[ax] - ws[ax] * hfrac / 2, 0.0)
                        cs.append(mkcell(c0, c0 + ws[ax] * hfrac))
                flb, fub, _t = cell_bounds(topo, counts_eval, N, *cs)
                if not (flb > 0) or float((fub / flb).log().mid()) > tgt:
                    ok = False
                    break
            if ok:
                return m
        return best

    def design(sc):
        h0s = [ws[k] * hfrac * sc * hadj for k in range(3)]
        mms = [band_mult(k, 2.0, (2.0, 1.4, 1.0, 0.7, 0.5), 4.5, h0s[k])
               for k in range(3)]
        mos = [band_mult(k, 4.2, (3.5, 2.4, 1.6, 1.0, 0.7, 0.5), 7.0, h0s[k])
               for k in range(3)]
        brs = [axis_breaks(k, upk, ws, S[k], delta, h0s[k], cks[k],
                           lnF, vpk, probe_ub, kind == "bal" and k in (0, 1),
                           mms[k], mos[k]) for k in range(3)]
        if kind == "bal":
            # v1=v2 kink alignment: axis-1 and axis-2 SHARE breakpoints on the
            # overlap, so the only straddling cells are the diagonal ones
            mg = min(h0s[0], h0s[1]) / 16
            m = sorted(set(brs[0]) | set(brs[1]))
            ded = [m[0]]
            for q in m[1:]:
                if q - ded[-1] > mg:
                    ded.append(q)
            for k in (0, 1):
                bk = [q for q in ded if q < S[k] - 1e-12]
                bk.append(S[k])
                brs[k] = bk
        axes = [build_axis(b) for b in brs]
        ns = tuple(len(a) for a in axes)
        ncells = ns[0] * ns[1] * ns[2]
        return h0s, axes, ns, ncells, ncells * t_eval * 1.2, mms, mos
    if scale_fix > 0:
        # fixed scale: grid fully deterministic given (topo, prec, hfrac, scale) --
        # REQUIRED for chunked runs (every chunk must build the identical grid)
        scale = scale_fix
        h0s, axes, ns, ncells, proj, mms, mos = design(scale)
    else:
        scale = 1.0
        h0s, axes, ns, ncells, proj, mms, mos = design(scale)
        while proj > budget_s and scale <= 64:
            scale *= 1.3
            h0s, axes, ns, ncells, proj, mms, mos = design(scale)
        while proj < budget_s / 2.6 and scale > 0.3:
            cand = design(scale / 1.25)
            if cand[4] > budget_s:
                break
            scale /= 1.25
            h0s, axes, ns, ncells, proj, mms, mos = cand
    n1, n2, n3 = ns
    ci, cn = (int(t) for t in chunk.split("/"))
    assert 0 <= ci < cn
    gridsig = [n1, n2, n3, round(sum(c[0] for c in axes[0]), 12),
               round(sum(c[0] for c in axes[1]), 12), round(sum(c[0] for c in axes[2]), 12)]
    if cn > 1:
        i0 = ci * n1 // cn; i1 = (ci + 1) * n1 // cn
        axes = [axes[0][i0:i1], axes[1], axes[2]]
        ncells = len(axes[0]) * n2 * n3
        print(f"[chunk] {ci}/{cn}: axis-1 cells [{i0},{i1}) of {n1} -> {ncells} cells",
              flush=True)
    # measured sample-cell widths at the FROZEN h0 (design receipt, not a bound)
    sw = []
    for k, off in ((0, 0.0), (0, 2.0), (0, 3.9)):
        lo = upk[0] + off * ws[0]
        cs = [mkcell(lo, lo + h0s[0] * (1 if off < 1.5 else (2 if off < 3.5 else 3.5))),
              mkcell(upk[1], upk[1] + h0s[1]),
              mkcell(max(0.0, upk[2] - h0s[2] / 2), upk[2] + h0s[2] / 2)]
        flb, fub, _tag = cell_bounds(topo, counts_eval, N, *cs)
        sw.append(round(float((fub / flb).log().mid()), 3) if flb > 0 else None)
    print(f"[grid] h0=({h0s[0]:.2e},{h0s[1]:.2e},{h0s[2]:.2e}) scale={scale:.3f} "
          f"cells={n1}x{n2}x{n3}={ncells} proj={proj:.0f}s "
          f"band mults mid={mms} out={mos} "
          f"sample cell widths (peak,+2w,+3.9w axis1)={sw} nats", flush=True)

    # (e) THE PASS: one ball evaluation per cell, lb and ub sums together
    S_lo = arb(0); S_hi = arb(0)
    nz_lb = 0; ndone = 0; nskip = 0
    diag_thresh = arb(lnZ_pre + 4.0).exp()   # cells whose ub*vol exceeds ~e^4 x Z/e^..
    diag_cells = []
    def diag_add(rec):
        diag_cells.append(rec)
        if len(diag_cells) >= 400:
            diag_cells.sort(reverse=True)
            del diag_cells[150:]
    hull_edge = arb(0)                    # ub mass in the u2/u3 edge-refined first cells
    e2 = axes[1][0][1] * 4; e3 = axes[2][0][1] * 4
    FBC[0] = 0; STR[0] = 0; PRS[0] = 0
    trep = time.time()
    for c1 in axes[0]:
        corner1 = c1[1] <= delta + 1e-15
        for c2 in axes[1]:
            lo2 = c2[0]
            c12 = corner1 and c2[1] <= delta + 1e-15
            v12 = c1[3] * c2[3]
            for c3 in axes[2]:
                if c12 and c3[1] <= delta + 1e-15:
                    nskip += 1
                    continue
                flb, fub, tag = cell_bounds(topo, counts_eval, N, c1, c2, c3)
                vol = v12 * c3[3]
                uv = fub * vol
                S_hi = S_hi + uv
                if uv > diag_thresh:
                    diag_add([float(uv.log().mid()), tag,
                              c1[0], c1[1], c2[0], c2[1], c3[0], c3[1]])
                if flb > 0:
                    S_lo = S_lo + flb * vol
                    nz_lb += 1
                if lo2 <= e2 or c3[0] <= e3:
                    hull_edge = hull_edge + fub * vol
                ndone += 1
                if ndone % 20000 == 0 and time.time() - trep > 30:
                    trep = time.time()
                    print(f"[pass] {ndone}/{ncells} cells cpu={cpu():.0f}s "
                          f"lnS_hi~{float(S_hi.log().mid()) if S_hi > 0 else float('-inf'):.2f} "
                          f"lnS_lo~{float(S_lo.log().mid()) if S_lo > 0 else float('-inf'):.2f}",
                          flush=True)
    t_pass = cpu() - t_bounds - t_scan
    print(f"[pass] done {ndone} cells ({nskip} corner-skipped) cpu_pass={t_pass:.1f}s "
          f"lb>0 on {nz_lb}", flush=True)

    # (f) assemble one-sided bounds (chunked runs bank partial sums; merge adds
    # them over the disjoint axis-1 stripes -- certified sums are additive)
    if cn > 1:
        cpu_total = cpu()
        part = {
            "probe_chunk": f"{ci}/{cn}", "grid_sig": gridsig,
            "topology": topo_name(topo), "topo_index": topo_i, "Nt": Nt,
            "prec": prec, "hfrac": hfrac, "scale": scale,
            "ln_Slo_lb": (f_dn(S_lo.log()) - 1e-9) if S_lo > 0 else None,
            "ln_Shi_ub": (f_up(S_hi.log()) + 1e-9) if S_hi > 0 else None,
            "complement_ln_ubs": {"corner": cb_ln, "tails": tail_lns,
                                  "S": S, "delta": delta},
            "cells_evaluated": ndone, "cells_corner_skipped": nskip,
            "cells_lb_positive": nz_lb, "fallback_cells": FBC[0],
            "straddle_cells": STR[0], "prism_cells": PRS[0],
            "hot_cells_ln_ubvol": diag_cells,
            "cost": {"cpu_s_total": round(cpu_total, 1),
                     "cpu_s_grid_pass": round(t_pass, 1),
                     "t_eval_ms": round(t_eval * 1000, 4)},
            "maxrss_mb": round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 1),
            "banked": time.strftime("%Y-%m-%d %H:%M:%S %Z"),
        }
        json.dump(part, open(out_path, "w"), indent=1)
        print(json.dumps(part, indent=1), flush=True)
        return 0
    Z_lo = lb_of(S_lo)
    Z_hi_ball = S_hi + cb + tail_ubs[0] + tail_ubs[1] + tail_ubs[2]
    Z_hi = ub_of(Z_hi_ball)
    if not (Z_lo > 0):
        raise RuntimeError("lower bound not positive -- grid too coarse")
    ln_lo = f_dn(arb(Z_lo).log())
    ln_hi = f_up(arb(Z_hi).log())
    width = ln_hi - ln_lo
    if Nt == 318:
        lnC = arb(318).lgamma()
        for D, n in COUNTS.items():
            lnC = lnC - arb(n + 1).lgamma()
        lnC_lo = f_dn(lnC); lnC_hi = f_up(lnC)
    else:
        lnC_lo = lnC_hi = None    # real-count rows: lnZcore only (engine convention)
    cpu_total = cpu()
    row = {
        "probe": "one-sided SD evidence, fixed needle grid + hull (topology-0)",
        "pattern_port": "positive integrand => interior "
                        "cell sums = rigorous LOWER; hull + one-sided tails = UPPER",
        "topology": topo_name(topo), "topo_index": topo_i, "Nt": Nt, "prec": prec,
        "lower_lnZcore": ln_lo, "upper_lnZcore": ln_hi, "width_nats": width,
        "lower_lnZ": (ln_lo + lnC_lo) if lnC_lo is not None else None,
        "upper_lnZ": (ln_hi + lnC_hi) if lnC_hi is not None else None,
        "complement_bounds": {"corner_ln_ub": cb_ln, "tail_ln_ubs": tail_lns,
                              "S": S, "delta": delta,
                              "reused_from": "sd_engine.axis_tail_sup/corner_bound (banked closers)"},
        "grid": {"u_peak": list(upk), "widths": ws, "h0": h0s,
                 "slopes_per_unit": cks, "hfrac": hfrac, "budget_scale": scale,
                 "hadj": hadj,
                 "sample_cell_widths_nats": sw, "fallback_cells": FBC[0],
                 "straddle_cells": STR[0], "prism_cells": PRS[0],
                 "band_mults_mid": mms, "band_mults_out": mos, "prism_cells": PRS[0],
                 "axis_cells": [n1, n2, n3], "cells_evaluated": ndone,
                 "cells_corner_skipped": nskip, "cells_lb_positive": nz_lb,
                 "edge_refined": "u2,u3 first cells split /2,/4,/8 toward 0",
                 "wq": WQ},
        "hull_edge_mass_ln": (float(hull_edge.log().mid()) if hull_edge > 0 else None),
        "hot_cells_ln_ubvol": diag_cells,
        "cost": {"cpu_s_total": round(cpu_total, 1), "cpu_s_scan": round(t_scan, 1),
                 "cpu_s_bounds": round(t_bounds, 1), "cpu_s_grid_pass": round(t_pass, 1),
                 "t_eval_ms": round(t_eval * 1000, 4),
                 "cpu_h_total": round(cpu_total / 3600, 3)},
        "maxrss_mb": round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 1),
        "wall_s": round(time.time() - t0, 1),
        "banked": time.strftime("%Y-%m-%d %H:%M:%S %Z"),
    }
    json.dump(row, open(out_path, "w"), indent=1)
    print(json.dumps(row, indent=1), flush=True)
    return 0

def merge(chunk_paths, out_path):
    """assemble the final probe row from chunked partial runs (disjoint axis-1
    stripes of ONE deterministic grid): lower = sum exp(ln_Slo_lb) rounded down;
    upper = sum exp(ln_Shi_ub) + complement bounds, rounded up."""
    ctx.prec = 128
    parts = [json.load(open(q)) for q in chunk_paths]
    sig = parts[0]["grid_sig"]
    for q in parts[1:]:
        assert q["grid_sig"] == sig, "grid signature mismatch -- chunks not identical"
    assert sorted(p["probe_chunk"] for p in parts) ==         sorted(f"{i}/{len(parts)}" for i in range(len(parts))), "chunk set incomplete"
    Zlo = arb(0); Zhi = arb(0)
    for q in parts:
        Zlo = Zlo + arb(q["ln_Slo_lb"]).exp()
        Zhi = Zhi + arb(q["ln_Shi_ub"]).exp()
    comp = parts[0]["complement_ln_ubs"]
    for lnu in [comp["corner"]] + comp["tails"]:
        Zhi = Zhi + arb(lnu + 1e-9).exp()
    Z_lo = lb_of(Zlo); Z_hi = ub_of(Zhi)
    ln_lo = f_dn(arb(Z_lo).log()); ln_hi = f_up(arb(Z_hi).log())
    NtM = parts[0].get("Nt", 318)
    if NtM == 318:
        lnC = arb(318).lgamma()
        for D, n in COUNTS.items():
            lnC = lnC - arb(n + 1).lgamma()
        lnZ_lo = ln_lo + f_dn(lnC); lnZ_hi = ln_hi + f_up(lnC)
    else:
        lnZ_lo = lnZ_hi = None
    cpu_tot = sum(q["cost"]["cpu_s_total"] for q in parts)
    row = {
        "probe": "one-sided SD evidence, fixed needle grid + centered duals + hull "
                 "(topology-0), merged from chunked capped runs",
        "topology": parts[0]["topology"], "topo_index": parts[0]["topo_index"],
        "Nt": NtM, "prec": parts[0]["prec"], "hfrac": parts[0]["hfrac"],
        "scale": parts[0]["scale"], "grid_sig": sig,
        "lower_lnZcore": ln_lo, "upper_lnZcore": ln_hi,
        "width_nats": ln_hi - ln_lo,
        "lower_lnZ": lnZ_lo, "upper_lnZ": lnZ_hi,
        "complement_ln_ubs": comp,
        "cells_evaluated": sum(q["cells_evaluated"] for q in parts),
        "fallback_cells": sum(q["fallback_cells"] for q in parts),
        "chunks": [{"chunk": q["probe_chunk"], "cpu_s": q["cost"]["cpu_s_total"],
                    "ln_Slo_lb": q["ln_Slo_lb"], "ln_Shi_ub": q["ln_Shi_ub"],
                    "maxrss_mb": q["maxrss_mb"], "banked": q["banked"]} for q in parts],
        "cost": {"cpu_s_total": round(cpu_tot, 1),
                 "cpu_h_total": round(cpu_tot / 3600, 3)},
        "banked": time.strftime("%Y-%m-%d %H:%M:%S %Z"),
    }
    json.dump(row, open(out_path, "w"), indent=1)
    print(json.dumps(row, indent=1), flush=True)
    return 0

def assemble(out_path):
    """SD_EVIDENCE.json: 15 certified per-topology lnZ enclosures, the registered
    uniform-1/15 evidence-weighted total, certified posterior-weight intervals,
    N=200/870 denominator rows, provenance."""
    ctx.prec = 128
    rows = []
    for i in range(15):
        f = (f"{OUTDIR}/rows/onesided_t{i:02d}_N318.json")
        rows.append(json.load(open(f)))
    Zlo = arb(0); Zhi = arb(0)
    for q in rows:
        Zlo = Zlo + arb(q["lower_lnZ"]).exp()
        Zhi = Zhi + arb(q["upper_lnZ"]).exp()
    lnZbar_lo = f_dn((Zlo / 15).log()) - 1e-12
    lnZbar_hi = f_up((Zhi / 15).log()) + 1e-12
    table = []
    for q in rows:
        plo = arb(q["lower_lnZ"]).exp() / Zhi
        phi = arb(q["upper_lnZ"]).exp() / Zlo
        table.append({
            "topo_index": q["topo_index"], "topology": q["topology"],
            "lnZ": [q["lower_lnZ"], q["upper_lnZ"]],
            "lnZcore": [q["lower_lnZcore"], q["upper_lnZcore"]],
            "width_nats": q["width_nats"],
            "posterior_weight": [max(0.0, f_dn(plo) - 1e-15),
                                 min(1.0, f_up(phi) + 1e-15)],
            "cpu_h": q["cost"]["cpu_h_total"],
            "cells": q.get("cells_evaluated"),
            "straddle_cells": q.get("straddle_cells"),
        })
    nrows = {}
    for nt in (200, 870):
        f = f"{OUTDIR}/rows/onesided_t00_N{nt}.json"
        if os.path.exists(f):
            q = json.load(open(f))
            nrows[str(nt)] = {"lnZcore": [q["lower_lnZcore"], q["upper_lnZcore"]],
                              "width_nats": q["width_nats"],
                              "cpu_h": q["cost"]["cpu_h_total"]}
    best = max(table, key=lambda r: r["lnZ"][0])
    seconds = sorted((r for r in table if r["topo_index"] != best["topo_index"]),
                     key=lambda r: -r["lnZ"][1])
    sep = (best["lnZ"][0] - seconds[0]["lnZ"][1]) if seconds else None
    out = {
        "artifact": "SD exact evidence, one-sided centered-grid route -- 15 rooted "
                    "topologies, certified lnZ enclosures (heirloom-C)",
        "registration": "priors s_k iid Exp(10), mu~Exp(1),"
                        " topology uniform 1/15, Omega=15 nonempty patterns; UNCHANGED",
        "data": "QUARTET_COLLAPSE.json, N=318, gate quartet A=Betta_Kurumba, B=Brahui,"
                " C=Gondi, D=Kannada",
        "per_topology": table,
        "evidence_weighted_total": {"lnZbar": [lnZbar_lo, lnZbar_hi],
                                    "prior": "uniform 1/15 (registered)"},
        "map_topology": {"topo_index": best["topo_index"], "topology": best["topology"],
                         "certified_lnBF_over_runner_up_at_least": sep},
        "N_scale_denominator_rows": nrows,
        "provenance": {
            "engine": "sd_onesided.py (centered duals + one-sided hull) on the pinned "
                      "sd_engine.py peeling circuit; complement bounds = sd_engine "
                      "axis_tail_sup/corner_bound verbatim",
            "gates": "ONESIDED_GATE.txt (anchor + plants + balanced degeneration/"
                     "symmetry/kink + real-counts anchor, all runs PASS)",
            "controls": "45-digit two-route anchor + plant a + plant b; "
                        "engine controls and registration untouched",
            "seeds": "none -- fully deterministic pipeline (no RNG anywhere)",
            "caps": "cap every run externally (kernel ulimit -t + timeout)",
        },
        "banked": time.strftime("%Y-%m-%d %H:%M:%S %Z"),
    }
    json.dump(out, open(out_path, "w"), indent=1)
    print(json.dumps({"lnZbar": out["evidence_weighted_total"]["lnZbar"],
                      "map": out["map_topology"],
                      "rows": len(table), "nrows": list(nrows)}, indent=1), flush=True)
    return 0

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["gate", "run", "merge", "assemble"])
    ap.add_argument("--topo", type=int, default=0)
    ap.add_argument("--prec", type=int, default=128)
    ap.add_argument("--hfrac", type=float, default=0.14)    # inner cell width / needle width
    ap.add_argument("--budget", type=float, default=1500.0) # grid-pass CPU budget (s)
    ap.add_argument("--scale", type=float, default=0.0)     # fixed grid scale (chunks)
    ap.add_argument("--chunk", type=str, default="0/1")
    ap.add_argument("--chunks", type=str, nargs="*", default=None)
    ap.add_argument("--out", type=str,
                    default=f"{OUTDIR}/ONESIDED_SD_PROBE.json")
    ap.add_argument("--nt", type=int, default=318)
    a = ap.parse_args()
    if a.cmd == "gate":
        sys.exit(gate())
    if a.cmd == "merge":
        sys.exit(merge(a.chunks, a.out))
    if a.cmd == "assemble":
        sys.exit(assemble(a.out))
    sys.exit(run(a.topo, a.prec, a.hfrac, a.budget, a.out, a.scale, a.chunk, a.nt))

if __name__ == "__main__":
    main()

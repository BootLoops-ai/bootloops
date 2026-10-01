# baller engine — certified evidence engine (peeling circuit over quartet pattern counts).
#!/usr/bin/env python3
# sd_engine.py -- certified SD evidence engine over quartet pattern counts
# (spec fixed BEFORE any evidence value; priors s_k iid Exp(10), mu ~ Exp(1), Omega = 15
# nonempty patterns; certified route = (U, a1, a2) decomposition with corner/tail
# one-sided ball bounds; nested flint acb.integral).
#
# Subcommands:
#   gate   -- engine-level 45-digit anchor + both planted mutations + independent
#             brute-force (2^k edge-config enumeration, exact Fractions) cross-checks
#             of the peeling circuit on non-pinned caterpillar and both balanced
#             regions. MUST pass before any row banks.
#   run    -- one certified topology row: --topo I --nt {318,200,870} [--rel] [--prec]
#
# Evaluation is the peeling circuit (never expanded monomials), 1-z via expm1
# (SD_FEASIBILITY implementation pin). Integer-count path uses integer powers only
# (entire integrand); real-count N-rows (pinned caterpillar only) use exp/log with an
# a1-edge strip trim entering the enclosure as a certified one-sided bound.
import json, time, sys, math, argparse, resource, os, itertools
from fractions import Fraction
from flint import arb, acb, ctx

# Fixture read root (QUARTET_COLLAPSE.json): defaults to this module's own
# directory (the vendored fixture). Output root: defaults to the CURRENT
# working directory — never this module's directory, whose unpinned files
# would trip baller.verify().
_HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.environ.get("WPG_BASE", _HERE)
OUT = os.environ.get("SD_ENGINE_OUT", os.getcwd())
QC = json.load(open(f"{BASE}/QUARTET_COLLAPSE.json"))
LEAVES = "ABCD"
BIT = {"A": 8, "B": 4, "C": 2, "D": 1}          # pattern string '1000' = A present = 8
COUNTS = {sum(BIT[LEAVES[i]] for i, ch in enumerate(p) if ch == "1"): int(n)
          for p, n in QC["pattern_counts"].items()}
N0 = sum(COUNTS.values())
assert N0 == 318
ALLPAT = list(range(1, 16))
LAM = 10                                          # increment prior rate (PINNED_SPEC pin 1)
PAT_1100 = BIT["A"] | BIT["B"]

TOPOLOGIES = []                                   # pinned order: 12 caterpillars, then 3 balanced
for c1, c2 in itertools.combinations(LEAVES, 2):
    rest = [x for x in LEAVES if x not in (c1, c2)]
    for k in (0, 1):
        TOPOLOGIES.append(("cat", (c1, c2, rest[k], rest[1 - k])))
TOPOLOGIES += [("bal", ("A", "B", "C", "D")), ("bal", ("A", "C", "B", "D")),
               ("bal", ("A", "D", "B", "C"))]
assert len(TOPOLOGIES) == 15

def topo_name(t):
    kind, (a, b, c, d) = t
    return f"((({a},{b}),{c}),{d})" if kind == "cat" else f"(({a},{b}),({c},{d}))"

EVALS = [0]

# ---------------- peeling circuit (generic in the scalar type) ----------------
def _join(qa, za, omza, qb, zb, omzb):
    out = {}
    for Da, pa in qa.items():
        fa = za * pa
        if Da == 0:
            fa = fa + omza
        for Db, pb in qb.items():
            fb = zb * pb
            if Db == 0:
                fb = fb + omzb
            D = Da | Db
            if D in out:
                out[D] = out[D] + fa * fb
            else:
                out[D] = fa * fb
    return out

def _assemble(qr, pieces, zero):
    """G_D = q_D(root) + sum over non-root nodes (1 - z_above) * q_D(node)."""
    G = dict(qr)
    for omz, q in pieces:
        for D, p in q.items():
            G[D] = G.get(D, zero) + omz * p
    G.pop(0, None)
    return G

def G_cat_z(perm, y1, om1, y2, om2, y3, om3, z3, om3e, z4, om4e, one):
    """caterpillar (((L1,L2),L3),L4) from edge survivals: leaf1/2 edges y1, edge v->w y2,
    leaf3 edge z3 (= y1*y2 unmutated), edge w->root y3, leaf4 edge z4 (= y1*y2*y3)."""
    b1, b2, b3, b4 = (BIT[l] for l in perm)
    zero = one * 0
    qv = {b1 | b2: y1 * y1, b1: y1 * om1, b2: om1 * y1, 0: om1 * om1}
    qw = _join(qv, y2, om2, {b3: one, 0: zero}, z3, om3e)
    qr = _join(qw, y3, om3, {b4: one, 0: zero}, z4, om4e)
    return _assemble(qr, [(om2, qv), (om3, qw), (om1, {b1: one}), (om1, {b2: one}),
                          (om3e, {b3: one}), (om4e, {b4: one})], zero)

def G_cat(perm, u1, u2, u3, one, plant_a=False):
    y1 = (-u1).exp(); y2 = (-u2).exp(); y3 = (-u3).exp()
    om1 = -(-u1).expm1(); om2 = -(-u2).expm1(); om3 = -(-u3).expm1()
    if plant_a:                                 # PLANT a: L3 edge length s1 instead of s1+s2
        z3, om3e = y1, om1
    else:
        z3 = y1 * y2; om3e = -(-(u1 + u2)).expm1()
    z4 = y1 * y2 * y3; om4e = -(-(u1 + u2 + u3)).expm1()
    return G_cat_z(perm, y1, om1, y2, om2, y3, om3, z3, om3e, z4, om4e, one)

def G_bal_z(perm, zc1, omc1, zc2, omc2, z12, om12, z34, om34, one):
    b1, b2, b3, b4 = (BIT[l] for l in perm)
    zero = one * 0
    q12 = {b1 | b2: zc1 * zc1, b1: zc1 * omc1, b2: omc1 * zc1, 0: omc1 * omc1}
    q34 = {b3 | b4: zc2 * zc2, b3: zc2 * omc2, b4: omc2 * zc2, 0: omc2 * omc2}
    qr = _join(q12, z12, om12, q34, z34, om34)
    return _assemble(qr, [(om12, q12), (om34, q34), (omc1, {b1: one}), (omc1, {b2: one}),
                          (omc2, {b3: one}), (omc2, {b4: one})], zero)

def G_bal(perm, v1, v2, v3, one, region):
    """balanced ((L1,L2),(L3,L4)); v1,v2 = mu-scaled cherry heights, v3 = mu-scaled root
    increment above the OLDER cherry (PINNED_SPEC pin 3). region 'A': v1 >= v2;
    'B': v2 >= v1. Formulas = analytic continuation of the respective region."""
    zc1 = (-v1).exp(); omc1 = -(-v1).expm1()
    zc2 = (-v2).exp(); omc2 = -(-v2).expm1()
    e12 = v3 if region == "A" else v2 - v1 + v3
    e34 = v1 - v2 + v3 if region == "A" else v3
    z12 = (-e12).exp(); om12 = -(-e12).expm1()
    z34 = (-e34).exp(); om34 = -(-e34).expm1()
    return G_bal_z(perm, zc1, omc1, zc2, omc2, z12, om12, z34, om34, one)

def lcore_int(G, counts, N, plant_b=False):
    """prod G_D^{n_D} / G_tot^N, integer counts (entire in the u's)."""
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

BADVAL = [None]   # set to acb nan at runtime

def lcore_real(G, counts_arb, N):
    """exp(sum n ln G - N ln Gtot), real (arb ball) counts; None if a log is branch-unsafe."""
    Gtot = None
    for D in ALLPAT:
        g = G.get(D)
        if g is None:
            continue
        Gtot = g if Gtot is None else Gtot + g
    lgt = _safe_log(Gtot)
    if lgt is None:
        return None
    s = -N * lgt
    for D, n in counts_arb.items():
        lg = _safe_log(G[D])
        if lg is None:
            return None
        s = s + n * lg
    return s.exp()

def _safe_log(g):
    if isinstance(g, acb):
        if not (g.real > 0):
            return None
        return g.log()
    if not (g > 0):
        return None
    return g.log()

# ---------------- certified sup bounds (arb interval balls) ----------------
def ball(lo, hi):
    lo = arb(lo); hi = arb(hi)
    return (lo + hi) / 2 + ((hi - lo) / 2) * arb(0, 1)

def hull01(b):
    """ball containing [0, ub(b)] for b >= 0 (one-sided bound folded into an enclosure)"""
    return b * arb(0.5, 0.5)

def ub_of(x):
    u = x.mid() + x.rad()
    return u + u.rad()

def sup_lcore(G, counts):
    """certified upper bound on Lcore = prod (G_D/G_tot)^{n_D} over the box the G-balls
    cover, using structural G_tot >= 1 and p_D <= 1 (sum_D n_D = N)."""
    out = arb(1)
    one = arb(1)
    for D, n in counts.items():
        g = G.get(D)
        p_ub = ub_of(g) if g is not None else arb(0)   # / G_tot lower bound 1
        if not (p_ub > 0):
            return arb(0)
        if not (p_ub < 1):
            p_ub = one
        out = out * (p_ub ** n if isinstance(n, int) else p_ub ** arb(n))
    return ub_of(out)

def _sup_over_zballs(topo, zdesc, counts):
    """sup Lcore where zdesc gives each edge-survival ball. For 'bal', evaluate both
    region formulas' enclosures and take the max (their union covers the region)."""
    kind, perm = topo
    one = arb(1)
    def om(z):
        return arb(1) - z
    if kind == "cat":
        y1, y2, y3, z3, z4 = zdesc
        G = G_cat_z(perm, y1, om(y1), y2, om(y2), y3, om(y3), z3, om(z3), z4, om(z4), one)
        return sup_lcore(G, counts)
    zc1, zc2, z12, z34 = zdesc
    G = G_bal_z(perm, zc1, om(zc1), zc2, om(zc2), z12, om(z12), z34, om(z34), one)
    return sup_lcore(G, counts)

def tail_bound(topo, Umax, counts):
    """one-sided bound: contribution of U > Umax <= 2 * max_k sup Lcore over
    {var_k >= Umax/3} (total fiber mass int w U^2 dU = 2). Edge survivals: any edge whose
    time-length is >= var_k lies in [0, e^-T]; every other edge survival lies in [0,1]
    (valid for ALL var values, so the unbounded region is fully covered)."""
    kind, perm = topo
    T = Umax / 3.0
    eT = ball(0, math.exp(-T))
    full = ball(0, 1.0)
    sups = []
    if kind == "cat":
        # (y1,y2,y3,z3=y1y2,z4=y1y2y3); var1>=T: y1,z3,z4 small; var2>=T: y2,z3,z4; var3>=T: y3,z4
        sups.append(_sup_over_zballs(topo, (eT, full, full, eT, eT), counts))
        sups.append(_sup_over_zballs(topo, (full, eT, full, eT, eT), counts))
        sups.append(_sup_over_zballs(topo, (full, full, eT, full, eT), counts))
    else:
        # (zc1,zc2,z12,z34); v1>=T: zc1 small (e12,e34 >= 0 in both regions -> [0,1]);
        # v2>=T: zc2 small; v3>=T: e12 >= v3 and e34 >= v3 in the applicable region union
        sups.append(_sup_over_zballs(topo, (eT, full, full, full), counts))
        sups.append(_sup_over_zballs(topo, (full, eT, full, full), counts))
        sups.append(_sup_over_zballs(topo, (full, full, eT, eT), counts))
    m = sups[0]
    for s in sups[1:]:
        if not (s < m):
            m = s
    return ub_of(arb(2) * m)

def corner_bound(topo, delta, counts):
    """one-sided bound: contribution of U <= delta <= sup Lcore over the corner
    (all mu-scaled vars in [0, delta]) times total fiber mass (<= 2)... mass of
    {U <= delta} is < 1, use 1."""
    kind, perm = topo
    lo = math.exp(-2 * delta)
    near1 = ball(lo, 1.0)          # any edge survival in the corner lies in [e^-2delta, 1]
    if kind == "cat":
        s = _sup_over_zballs(topo, (near1, near1, near1, near1, near1), counts)
    else:
        s = _sup_over_zballs(topo, (near1, near1, near1, near1), counts)
    return ub_of(s)

# ---------------- float scoping scan (never banked; sets tolerance floors only) ----------------
class F:
    __slots__ = ("v",)
    def __init__(s, v): s.v = v if not isinstance(v, F) else v.v
    def exp(s): return F(math.exp(s.v))
    def expm1(s): return F(math.expm1(s.v))
    def __mul__(s, o): return F(s.v * (o.v if isinstance(o, F) else o))
    __rmul__ = __mul__
    def __add__(s, o): return F(s.v + (o.v if isinstance(o, F) else o))
    __radd__ = __add__
    def __sub__(s, o): return F(s.v - (o.v if isinstance(o, F) else o))
    def __rsub__(s, o): return F((o.v if isinstance(o, F) else o) - s.v)
    def __neg__(s): return F(-s.v)

def float_lnzcore_prelim(topo, counts, N):
    from scipy.special import kv
    kind, perm = topo
    def lnL(u1, u2, u3, region):
        G = (G_cat(perm, F(u1), F(u2), F(u3), F(1.0)) if kind == "cat"
             else G_bal(perm, F(u1), F(u2), F(u3), F(1.0), region))
        gt = sum(g.v for g in G.values())
        s = -N * math.log(gt)
        for D, n in counts.items():
            gv = G[D].v
            if gv <= 0:
                return -1e300
            s += n * math.log(gv)
        return s
    best = -1e300; tot = -1e300
    for i in range(36):
        U = 0.15 * 1.32 ** i
        c = 10.0 * U
        kvv = kv(2, 2.0 * math.sqrt(c))
        if not (kvv > 0):
            continue
        lw = math.log(2000.0) - math.log(c) + math.log(kvv) + 2.0 * math.log(U)
        m = 12
        for a in range(1, m):
            for b in range(1, m - a):
                a1, a2 = a / m, b / m
                reg = "A" if a1 >= a2 else "B"
                v = lnL(U * a1, U * a2, U * (1 - a1 - a2), reg) + lw
                best = max(best, v)
                hi, lo2 = max(tot, v), min(tot, v)
                tot = hi + math.log1p(math.exp(lo2 - hi))
    return best, tot + math.log(0.28 / (12 * 12))   # crude cell-volume factor

# ---------------- certified evidence for one topology ----------------
# AXIS-ALIGNED ROUTE (v2):
# integrate directly in (u1,u2,u3) over [0,S1]x[0,S2]x[0,S3] minus the corner cube
# [0,delta]^3, decomposed into 3 axis-aligned boxes (each keeps one coordinate >= delta,
# so U = u1+u2+u3 >= delta on every box and K_2(2 sqrt(10 U)) stays branch-safe).
# w(U) moves into the innermost integrand. Corner cube and per-axis tails enter as
# one-sided certified bounds (fiber mass int w du = 1). Keying fix:
# peak-refined prelim, flat rel_goal/3 at all levels, abs floors spread over box
# lengths only (no stacked discounts).

def axis_scan_peak(topo, counts, N):
    """float u-space peak scan (never banked; keying only): coarse log-grid + refine,
    per-axis half-nat widths."""
    from scipy.special import kv
    kind, perm = topo
    def lnF(u1, u2, u3):
        if u1 <= 0 or u2 <= 0 or u3 <= 0:
            return -1e300
        U = u1 + u2 + u3
        c = 10.0 * U
        k = kv(2, 2.0 * math.sqrt(c))
        if not (k > 0):
            return -1e300
        s = math.log(2000.0) - math.log(c) + math.log(k)
        reg = "A" if u1 >= u2 else "B"
        G = (G_cat(perm, F(u1), F(u2), F(u3), F(1.0)) if kind == "cat"
             else G_bal(perm, F(u1), F(u2), F(u3), F(1.0), reg))
        gt = sum(g.v for g in G.values())
        s += -N * math.log(gt)
        for D, n in counts.items():
            gv = G[D].v
            if gv <= 0:
                return -1e300
            s += n * math.log(gv)
        return s
    grid = [5e-4 * (6.0 / 5e-4) ** (i / 15.0) for i in range(16)]
    best = (-1e300, 0.3, 0.3, 0.3)
    for u1 in grid:
        for u2 in grid:
            for u3 in grid:
                v = lnF(u1, u2, u3)
                if v > best[0]:
                    best = (v, u1, u2, u3)
    for _ in range(4):
        v0, x, y, z = best
        for f1 in (0.7, 0.85, 1.0, 1.18, 1.4):
            for f2 in (0.7, 0.85, 1.0, 1.18, 1.4):
                for f3 in (0.7, 0.85, 1.0, 1.18, 1.4):
                    v = lnF(x * f1, y * f2, z * f3)
                    if v > best[0]:
                        best = (v, x * f1, y * f2, z * f3)
    vpk, u1p, u2p, u3p = best
    upk = (u1p, u2p, u3p)
    ws = []
    for ax in range(3):
        x0 = upk[ax]
        h = x0 * 0.05 + 1e-5
        def at(t):
            c = list(upk); c[ax] = t
            return lnF(*c)
        while at(x0 + h) > vpk - 0.5 and h < 6.0:
            h *= 1.5
        ws.append(h)
    lnZ_pre = vpk + sum(math.log(min(2.0 * w, 6.0)) for w in ws)
    return vpk, upk, ws, lnZ_pre

def axis_tail_sup(topo, k, T, counts):
    """certified sup of Lcore over {var_k >= T}; contribution bound = sup * mass(<=1)."""
    kind, perm = topo
    eT = ball(0, math.exp(-T))
    full = ball(0, 1.0)
    if kind == "cat":
        zd = [(eT, full, full, eT, eT), (full, eT, full, eT, eT),
              (full, full, eT, full, eT)][k]
    else:
        zd = [(eT, full, full, full), (full, eT, full, full),
              (full, full, eT, eT)][k]
    return ub_of(_sup_over_zballs(topo, zd, counts))

def _segments(lo, hi, cuts):
    """segment boundaries: [lo..hi] split at interior cut points (floats, sorted)."""
    pts = [lo] + [c for c in sorted(cuts) if lo + 1e-12 < c < hi - 1e-12] + [hi]
    return list(zip(pts[:-1], pts[1:]))

def _seg_integral(f, lo, hi, cuts, rel_t, abs_t):
    """sum of certified integrals over needle-informed segments; per-segment abs
    budget = abs_t / nseg so the total stays <= abs_t."""
    segs = _segments(lo, hi, cuts)
    at = abs_t / len(segs)
    tot = None
    for a, b in segs:
        v = acb.integral(f, a, b, rel_tol=rel_t, abs_tol=at)
        tot = v if tot is None else tot + v
    return tot

def certified_topology(topo, Nt, rel_goal, prec, delta=0.05, verbose=True):
    ctx.prec = prec
    kind, perm = topo
    BADVAL[0] = acb(arb("nan"))
    if Nt == 318:
        counts = dict(COUNTS); N = 318; real = False; counts_arb = None
    else:
        assert topo == TOPOLOGIES[0], "N-scaled rows are registered on the pinned caterpillar only"
        fscale = arb(Nt) / 318
        counts = {D: n * Nt / 318.0 for D, n in COUNTS.items()}
        counts_arb = {D: arb(n) * fscale for D, n in COUNTS.items()}
        N = Nt; real = True
    t0 = time.time(); c0 = resource.getrusage(resource.RUSAGE_SELF)
    EVALS[0] = 0

    vpk, upk, ws, lnZ_pre = axis_scan_peak(topo, counts, N)
    if verbose:
        print(f"[scan-ax] peak lnF={vpk:.2f} at u*=({upk[0]:.4f},{upk[1]:.4f},{upk[2]:.4f}) "
              f"widths=({ws[0]:.4f},{ws[1]:.4f},{ws[2]:.4f}) lnZ_pre~{lnZ_pre:.1f}", flush=True)
    tgt_ln = lnZ_pre + math.log(rel_goal) - 1.0   # keying fix: was -3.0 (20x over-tight)

    S = []; tail_sum = arb(0); tl_lns = []
    for k in range(3):
        T = max(1.2, upk[k] * 3.0 + 0.5)
        while True:
            sb = axis_tail_sup(topo, k, T, counts)
            sl = float(sb.log().mid()) if sb > 0 else -1e300
            if sl < tgt_ln - 1.1 or T > 700:
                break
            T *= 1.6
        if sl >= tgt_ln - 1.1:
            raise RuntimeError(f"axis-{k} tail bound will not close (T={T}, ln={sl:.1f})")
        S.append(T); tail_sum = tail_sum + sb; tl_lns.append(sl)
    cb = corner_bound(topo, delta, counts)
    cb_ln = float(cb.log().mid()) if cb > 0 else -1e300
    if cb_ln >= tgt_ln:
        raise RuntimeError("corner bound too weak -- shrink delta")
    if verbose:
        print(f"[bounds-ax] S=({S[0]:.1f},{S[1]:.1f},{S[2]:.1f}) tail lns={[round(x,1) for x in tl_lns]} "
              f"corner ln<={cb_ln:.1f} target ln {tgt_ln:.1f}", flush=True)

    ETA = 1e-5 if real else 0.0
    strip = arb(0)
    if real:   # {u1 <= ETA}: A|B-splitting patterns vanish; sup * mass(<=1)
        strip = ub_of(_sup_over_zballs(topo, (ball(math.exp(-ETA), 1.0), ball(0, 1),
                                              ball(0, 1), ball(0, 1), ball(0, 1)), counts))
        st_ln = float(strip.log().mid()) if strip > 0 else -1e300
        if st_ln >= tgt_ln:
            raise RuntimeError("u1-strip bound too weak")
        if verbose:
            print(f"[bounds-ax] u1-strip eta={ETA} ln<={st_ln:.1f}", flush=True)

    Zp = arb(lnZ_pre).exp()
    rel3 = arb(rel_goal / 3.0)
    # needle-informed pre-splits per axis (validity unaffected: certified segments sum;
    # avoids per-node deep bisection cascades from needle scales ~1e-3 in unit intervals)
    CUTS = []
    for k in range(3):
        u, w = upk[k], ws[k]
        CUTS.append([u - 8 * w, u + 8 * w, u + 64 * w])

    def lcore_at(u1, u2, u3, region):
        EVALS[0] += 1
        G = (G_cat(perm, u1, u2, u3, one_a) if kind == "cat"
             else G_bal(perm, u1, u2, u3, one_a, region))
        if real:
            v = lcore_real(G, counts_arb, N)
            return BADVAL[0] if v is None else v
        return lcore_int(G, counts, N)
    one_a = acb(1)

    def integrand(u1, u2, u3, region, an):
        U = u1 + u2 + u3
        if an and not (U.real > 0):
            return BADVAL[0]
        c = 10 * U
        w = acb(2000) / c * (2 * c.sqrt()).bessel_k(2)
        return lcore_at(u1, u2, u3, region) * w

    def box_integral(b1, b2, b3, region_split):
        """nested certified integral over b1 x b2 x b3 (outer u1, mid u2, inner u3),
        needle pre-splits at every level. region_split: None (cat) or region for bal."""
        L1 = b1[1] - b1[0]; L2 = b2[1] - b2[0]
        abs_box = Zp * arb(rel_goal * 0.3)
        abs_mid = abs_box * arb(0.5 / max(L1, 1e-3))
        abs_in = abs_mid * arb(0.5 / max(L2, 1e-3))
        def f1(u1, an1):
            def f2(u2, an2):
                if an2 and not ((u1 + u2).real > 0):
                    return BADVAL[0]
                return _seg_integral(lambda u3, an3: integrand(u1, u2, u3, region_split, an3),
                                     b3[0], b3[1], CUTS[2], rel3, abs_in)
            return _seg_integral(f2, b2[0], b2[1], CUTS[1], rel3, abs_mid)
        return _seg_integral(f1, b1[0], b1[1], CUTS[0], rel3, abs_box)

    def box_integral_bal(b1, b2, b3):
        """balanced: mid (v2) integral splits at v2 = v1 (region kink); outer pieces are
        arranged by the caller so the split point is analytic in v1 on each piece.
        Needle pre-splits on outer and inner levels (mid keeps the region split)."""
        L1 = b1[1] - b1[0]; L2 = b2[1] - b2[0]
        abs_box = Zp * arb(rel_goal * 0.3)
        abs_mid = abs_box * arb(0.5 / max(L1, 1e-3))
        abs_in = abs_mid * arb(0.5 / max(L2, 1e-3))
        def inner(u1, u2, reg):
            return _seg_integral(lambda u3, an3: integrand(u1, u2, u3, reg, an3),
                                 b3[0], b3[1], CUTS[2], rel3, abs_in)
        def f1(u1, an1):
            def f2A(u2, an2):
                if an2 and not ((u1 + u2).real > 0):
                    return BADVAL[0]
                return inner(u1, u2, "A")
            def f2B(u2, an2):
                if an2 and not ((u1 + u2).real > 0):
                    return BADVAL[0]
                return inner(u1, u2, "B")
            lo, hi = b2
            if b1[1] <= lo:                      # v1 <= v2 everywhere: all region B
                return acb.integral(f2B, lo, hi, rel_tol=rel3, abs_tol=abs_mid)
            if b1[0] >= hi:                      # v1 >= v2 everywhere: all region A
                return acb.integral(f2A, lo, hi, rel_tol=rel3, abs_tol=abs_mid)
            return (acb.integral(f2A, lo, u1, rel_tol=rel3, abs_tol=abs_mid)
                    + acb.integral(f2B, u1, hi, rel_tol=rel3, abs_tol=abs_mid))
        return _seg_integral(f1, b1[0], b1[1], CUTS[0], rel3, abs_box)

    lo1 = ETA if real else 0.0
    boxes = [((delta, S[0]), (0.0, S[1]), (0.0, S[2])),
             ((lo1, delta), (delta, S[1]), (0.0, S[2])),
             ((lo1, delta), (0.0, delta), (delta, S[2]))]
    box_vals = []; box_evals = []; box_secs = []
    for i, (b1, b2, b3) in enumerate(boxes):
        e0, s0 = EVALS[0], time.time()
        if kind == "cat":
            v = box_integral(b1, b2, b3, None)
        else:
            # outer split at S2 boundary for box 1 (kink-split analyticity)
            if i == 0 and b1[1] > b2[1]:
                v = (box_integral_bal((b1[0], b2[1]), b2, b3)
                     + box_integral((b2[1], b1[1]), b2, b3, "A"))
            else:
                v = box_integral_bal(b1, b2, b3)
        box_vals.append(v)
        box_evals.append(EVALS[0] - e0); box_secs.append(round(time.time() - s0, 1))
        if verbose:
            print(f"[box{i+1}] evals={box_evals[-1]} {box_secs[-1]}s "
                  f"val~{v.real.str(6)}", flush=True)

    Z = box_vals[0].real + box_vals[1].real + box_vals[2].real \
        + hull01(cb) + hull01(tail_sum) + hull01(strip)
    if not (Z > 0):
        raise RuntimeError(f"Z ball not positive: {Z}")
    lnZcore = Z.log()
    if not real:
        lnC = arb(N).lgamma()
        for D, n in COUNTS.items():
            lnC = lnC - arb(n + 1).lgamma()
        lnZ = lnZcore + lnC
    else:
        lnZ = None
    c1 = resource.getrusage(resource.RUSAGE_SELF)
    return {
        "topology": topo_name(topo), "kind": kind, "Nt": Nt, "route": "axis-aligned-v2",
        "lnZcore_mid": float(lnZcore.mid()), "lnZcore_rad": float(lnZcore.rad()),
        "lnZcore_enclosure": lnZcore.str(30),
        "lnZ_mid": float(lnZ.mid()) if lnZ is not None else None,
        "lnZ_rad": float(lnZ.rad()) if lnZ is not None else None,
        "lnZ_enclosure": lnZ.str(30) if lnZ is not None else None,
        "corner_ln_ub": cb_ln, "tail_ln_ubs": tl_lns,
        "strip_ln_ub": (float(strip.log().mid()) if real and strip > 0 else None),
        "delta": delta, "S": S, "u_peak": list(upk), "widths": ws,
        "rel_goal": rel_goal, "prec": prec,
        "evals": EVALS[0], "box_evals": box_evals, "box_secs": box_secs,
        "wall_s": round(time.time() - t0, 2),
        "cpu_s": round(c1.ru_utime + c1.ru_stime - c0.ru_utime - c0.ru_stime, 2),
        "maxrss_mb": round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 1),
    }

# --- old sheared route kept for the diff record (structurally hostile) ---
def certified_topology_sheared_v1(topo, Nt, rel_goal, prec, delta=0.05, verbose=True):
    ctx.prec = prec
    kind, perm = topo
    BADVAL[0] = acb(arb("nan"))
    if Nt == 318:
        counts = dict(COUNTS); N = 318; real = False; counts_arb = None
    else:
        assert topo == TOPOLOGIES[0], "N-scaled rows are registered on the pinned caterpillar only"
        fscale = arb(Nt) / 318
        counts = {D: n * Nt / 318.0 for D, n in COUNTS.items()}   # float shadow (scan/sup)
        counts_arb = {D: arb(n) * fscale for D, n in COUNTS.items()}
        N = Nt; real = True
    t0 = time.time(); c0 = resource.getrusage(resource.RUSAGE_SELF)
    EVALS[0] = 0

    best, lnZ_pre = float_lnzcore_prelim(topo, counts, N)
    if verbose:
        print(f"[scan] float prelim: peak lnF ~ {best:.1f}, lnZcore ~ {lnZ_pre:.1f}", flush=True)

    tgt_ln = lnZ_pre + math.log(rel_goal) - 3.0
    Umax = 30.0
    while True:
        tb = tail_bound(topo, Umax, counts)
        tb_ln = float(tb.log().mid()) if tb > 0 else -1e300
        if tb_ln < tgt_ln or Umax > 4000:
            break
        Umax *= 1.6
    if tb_ln >= tgt_ln:
        raise RuntimeError("tail bound will not close")
    cb = corner_bound(topo, delta, counts)
    cb_ln = float(cb.log().mid()) if cb > 0 else -1e300
    if verbose:
        print(f"[bounds] delta={delta} corner ln<={cb_ln:.1f}; Umax={Umax:.1f} tail ln<={tb_ln:.1f}; abs target ln {tgt_ln:.1f}", flush=True)
    if cb_ln >= tgt_ln:
        raise RuntimeError("corner bound too weak -- shrink delta")

    ETA = 1e-6 if real else 0.0
    strip = arb(0)
    if real:  # a1-edge strip trim bound: {a1 <= ETA}: y1 in [e^-ETA*Umax, 1], rest free
        s_sup = _sup_over_zballs(topo, (ball(math.exp(-ETA * Umax), 1.0), ball(0, 1), ball(0, 1),
                                        ball(0, 1), ball(0, 1)), counts)
        strip = ub_of(arb(2) * arb(ETA) * s_sup)
        st_ln = float(strip.log().mid()) if strip > 0 else -1e300
        if verbose:
            print(f"[bounds] a1-strip eta={ETA} ln<={st_ln:.1f}", flush=True)
        if st_ln >= tgt_ln:
            raise RuntimeError("strip bound too weak -- shrink ETA")

    from scipy.special import kv as besselkv
    one_a = acb(1)
    abs_outer = (arb(lnZ_pre) + arb(math.log(rel_goal)) + arb(math.log(0.2))).exp()
    rel_outer = arb(rel_goal)
    rel_mid = arb(rel_goal * 0.3)
    rel_in = arb(rel_goal * 0.1)
    span = Umax - delta

    def lcore_at(u1, u2, u3, region):
        EVALS[0] += 1
        G = (G_cat(perm, u1, u2, u3, one_a) if kind == "cat"
             else G_bal(perm, u1, u2, u3, one_a, region))
        if real:
            v = lcore_real(G, counts_arb, N)
            return BADVAL[0] if v is None else v
        return lcore_int(G, counts, N)

    def T_of_U(U, abs_in):
        if kind == "cat":
            def inner(a1, an1):
                return acb.integral(lambda a2, an2: lcore_at(U * a1, U * a2, U * (1 - a1 - a2), None),
                                    0, 1 - a1, rel_tol=rel_in, abs_tol=abs_in)
            return acb.integral(inner, ETA, 1, rel_tol=rel_mid, abs_tol=abs_in)
        def inner_lo(a1, an1):   # a1 in [0,1/2]: a2 in [0,a1] reg A, [a1,1-a1] reg B
            pa = acb.integral(lambda a2, an2: lcore_at(U * a1, U * a2, U * (1 - a1 - a2), "A"),
                              0, a1, rel_tol=rel_in, abs_tol=abs_in)
            pb = acb.integral(lambda a2, an2: lcore_at(U * a1, U * a2, U * (1 - a1 - a2), "B"),
                              a1, 1 - a1, rel_tol=rel_in, abs_tol=abs_in)
            return pa + pb
        def inner_hi(a1, an1):   # a1 in [1/2,1]: a2 in [0,1-a1], all reg A
            return acb.integral(lambda a2, an2: lcore_at(U * a1, U * a2, U * (1 - a1 - a2), "A"),
                                0, 1 - a1, rel_tol=rel_in, abs_tol=abs_in)
        h = arb(1) / 2
        return (acb.integral(inner_lo, 0, h, rel_tol=rel_mid, abs_tol=abs_in)
                + acb.integral(inner_hi, h, 1, rel_tol=rel_mid, abs_tol=abs_in))

    def outer(U, an):
        if an and not (U.real > 0):
            return BADVAL[0]
        c = 10 * U
        w = acb(2000) / c * (2 * c.sqrt()).bessel_k(2)
        Umid = float(U.real.mid())
        kvv = float(besselkv(2, 2 * math.sqrt(10 * Umid)))
        wU2f = max(2000.0 / (10 * Umid) * kvv * Umid ** 2, 1e-280)
        abs_in = abs_outer * arb(0.3) / (arb(wU2f) * arb(span))
        if not (abs_in < 1):
            abs_in = arb(1)
        return w * U * U * T_of_U(U, abs_in)

    main = acb.integral(outer, arb(delta), arb(Umax), rel_tol=rel_outer, abs_tol=abs_outer)
    Z = main.real + hull01(cb) + hull01(tb) + hull01(strip)
    if not (Z > 0):
        raise RuntimeError(f"Z ball not positive: {Z}")
    lnZcore = Z.log()
    if not real:
        lnC = arb(N).lgamma()
        for D, n in COUNTS.items():
            lnC = lnC - arb(n + 1).lgamma()
        lnZ = lnZcore + lnC
    else:
        lnZ = None
    c1 = resource.getrusage(resource.RUSAGE_SELF)
    return {
        "topology": topo_name(topo), "kind": kind, "Nt": Nt,
        "lnZcore_mid": float(lnZcore.mid()), "lnZcore_rad": float(lnZcore.rad()),
        "lnZcore_enclosure": lnZcore.str(30),
        "lnZ_mid": float(lnZ.mid()) if lnZ is not None else None,
        "lnZ_rad": float(lnZ.rad()) if lnZ is not None else None,
        "lnZ_enclosure": lnZ.str(30) if lnZ is not None else None,
        "corner_ln_ub": cb_ln, "tail_ln_ub": tb_ln,
        "strip_ln_ub": (float(strip.log().mid()) if real and strip > 0 else None),
        "delta": delta, "Umax": Umax, "rel_goal": rel_goal, "prec": prec,
        "evals": EVALS[0], "wall_s": round(time.time() - t0, 2),
        "cpu_s": round(c1.ru_utime + c1.ru_stime - c0.ru_utime - c0.ru_stime, 2),
        "maxrss_mb": round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 1),
    }

# ---------------- independent brute-force route (exact Fractions) ----------------
def brute_G(topo, params):
    """G_D by 2^k enumeration over edge-survival configs, exact Fractions.
    cat params: (y1, y2, y3) Fractions -> edge z's (y1,y1,y2,y1*y2,y3,y1*y2*y3).
    bal params: (zc1, zc2, z12, z34) Fractions."""
    kind, perm = topo
    b = {l: BIT[l] for l in LEAVES}
    if kind == "cat":
        L1, L2, L3, L4 = perm
        y1, y2, y3 = params
        edges = {"e1": (y1, b[L1]), "e2": (y1, b[L2]), "ev": (y2, 0),
                 "e3": (y1 * y2, b[L3]), "ew": (y3, 0), "e4": (y1 * y2 * y3, b[L4])}
        below = {"v": ["e1", "e2"], "w": ["e1", "e2", "ev", "e3"],
                 "r": ["e1", "e2", "ev", "e3", "ew", "e4"]}
        # leaf reachability from node: leaf -> edges on the node->leaf path
        reach = {"v": {b[L1]: ["e1"], b[L2]: ["e2"]},
                 "w": {b[L1]: ["ev", "e1"], b[L2]: ["ev", "e2"], b[L3]: ["e3"]},
                 "r": {b[L1]: ["ew", "ev", "e1"], b[L2]: ["ew", "ev", "e2"],
                       b[L3]: ["ew", "e3"], b[L4]: ["e4"]}}
        above = {"v": "ev", "w": "ew"}
        leaf_above = {b[L1]: "e1", b[L2]: "e2", b[L3]: "e3", b[L4]: "e4"}
    else:
        L1, L2, L3, L4 = perm
        zc1, zc2, z12, z34 = params
        edges = {"e1": (zc1, b[L1]), "e2": (zc1, b[L2]), "e12": (z12, 0),
                 "e3": (zc2, b[L3]), "e4": (zc2, b[L4]), "e34": (z34, 0)}
        below = {"v": ["e1", "e2"], "w": ["e3", "e4"],
                 "r": ["e1", "e2", "e12", "e3", "e4", "e34"]}
        reach = {"v": {b[L1]: ["e1"], b[L2]: ["e2"]}, "w": {b[L3]: ["e3"], b[L4]: ["e4"]},
                 "r": {b[L1]: ["e12", "e1"], b[L2]: ["e12", "e2"],
                       b[L3]: ["e34", "e3"], b[L4]: ["e34", "e4"]}}
        above = {"v": "e12", "w": "e34"}
        leaf_above = {b[L1]: "e1", b[L2]: "e2", b[L3]: "e3", b[L4]: "e4"}

    def qtable(node):
        ed = below[node]
        q = {}
        for mask in range(1 << len(ed)):
            pr = Fraction(1)
            alive = set()
            for j, e in enumerate(ed):
                if mask >> j & 1:
                    pr *= edges[e][0]; alive.add(e)
                else:
                    pr *= 1 - edges[e][0]
            pat = 0
            for leafbit, path in reach[node].items():
                if all(e in alive for e in path):
                    pat |= leafbit
            q[pat] = q.get(pat, Fraction(0)) + pr
        return q

    G = {}
    qr = qtable("r")
    for D, p in qr.items():
        G[D] = G.get(D, Fraction(0)) + p
    for node in ("v", "w"):
        omz = 1 - edges[above[node]][0]
        for D, p in qtable(node).items():
            G[D] = G.get(D, Fraction(0)) + omz * p
    for leafbit, e in leaf_above.items():
        G[leafbit] = G.get(leafbit, Fraction(0)) + (1 - edges[e][0])
    G.pop(0, None)
    return G

def _frac_to_arb(fr):
    return arb(fr.numerator) / arb(fr.denominator)

def gate():
    ctx.prec = 256
    one = arb(1)
    lines = ["", "== ENGINE GATE (sd_engine.py arb peeling circuit, prec 256; run: "
             + time.strftime("%Y-%m-%d %H:%M:%S %Z") + ") =="]
    ok_all = True
    def check(name, ok, got):
        nonlocal ok_all
        lines.append(f"{'PASS' if ok else 'FAIL'}  {name}   got: {got}")
        ok_all = ok_all and ok
    # 1) 45-digit anchor + plants through the production evaluator path
    u1 = (arb(10) / 9).log(); u2 = (arb(5) / 4).log(); u3 = (arb(10) / 7).log()
    lnC = arb(318).lgamma()
    for D, n in COUNTS.items():
        lnC = lnC - arb(n + 1).lgamma()
    def logL(pa=False, pb=False):
        G = G_cat(("A", "B", "C", "D"), u1, u2, u3, one, plant_a=pa)
        return lnC + lcore_int(G, COUNTS, 318, plant_b=pb).log()
    ANCHOR = arb("-264.252750460743437038228550850745165801524406")
    va, vpa, vpb = logL(), logL(pa=True), logL(pb=True)
    check("engine reproduces 45-digit anchor", abs(va - ANCHOR) < arb("1e-41"), va.str(45))
    check("engine plant a FIRES (wrong C-edge -> -311.493...)",
          abs(vpa - arb("-311.49316474472421895")) < arb("1e-9") and abs(vpa - ANCHOR) > 1,
          vpa.str(20))
    check("engine plant b FIRES (1100 dropped -> -222.822...)",
          abs(vpb - arb("-222.8224339223726654")) < arb("1e-9") and abs(vpb - ANCHOR) > 1,
          vpb.str(20))
    # 2) independent brute-force cross-checks (exact Fractions vs peeling balls)
    def xcheck(name, topo, params, engineG):
        Gf = brute_G(topo, params)
        ok = set(Gf) == set(engineG)
        worst = arb(0)
        for D in Gf:
            d = abs(engineG[D] - _frac_to_arb(Gf[D]))
            ok = ok and d < arb("1e-70")
            if not (d < worst):
                worst = d
        check(name, ok, f"max|dG| < {worst.str(3)} over {len(Gf)} patterns")
    y = (Fraction(9, 10), Fraction(4, 5), Fraction(7, 10))
    for idx, nm in ((0, "pinned cat (((A,B),C),D)"), (9, "relabeled cat " + topo_name(TOPOLOGIES[9]))):
        kind, perm = TOPOLOGIES[idx]
        uu = tuple(-_frac_to_arb(t).log() for t in y)
        xcheck(f"brute-force x-check, {nm}", TOPOLOGIES[idx], y,
               G_cat(perm, uu[0], uu[1], uu[2], one))
    # balanced, region A point: v1=ln2 >= v2=ln(3/2), v3=ln(5/4)
    for idx in (12, 14):
        kind, perm = TOPOLOGIES[idx]
        v1 = arb(2).log(); v2 = (arb(3) / 2).log(); v3 = (arb(5) / 4).log()
        params = (Fraction(1, 2), Fraction(2, 3), Fraction(4, 5), Fraction(3, 5))
        # zc1=1/2, zc2=2/3, z12=e^-v3=4/5, z34=e^-(v1-v2+v3)=e^-ln(5/3)=3/5
        xcheck(f"brute-force x-check, {topo_name(TOPOLOGIES[idx])} region A",
               TOPOLOGIES[idx], params, G_bal(perm, v1, v2, v3, one, "A"))
        # region B point: v1=ln(3/2) <= v2=ln2, v3=ln(5/4): z12=e^-(v2-v1+v3)=3/5, z34=4/5
        v1b = (arb(3) / 2).log(); v2b = arb(2).log()
        paramsB = (Fraction(2, 3), Fraction(1, 2), Fraction(3, 5), Fraction(4, 5))
        xcheck(f"brute-force x-check, {topo_name(TOPOLOGIES[idx])} region B",
               TOPOLOGIES[idx], paramsB, G_bal(perm, v1b, v2b, v3, one, "B"))
    # 3) registration invariants through the engine: G_tot in [1,4] on a narrow ball
    # (narrow because naive interval growth widens wide-ball enclosures; structural
    # bounds are proved in SD_FEASIBILITY sec 1.2, this is a circuit sanity row only)
    e = ball(0.495, 0.505)
    Gb = G_cat(("A", "B", "C", "D"), e, e, e, one)
    gt = None
    for D in ALLPAT:
        gt = Gb[D] if gt is None else gt + Gb[D]
    check("G_tot in [1,4] on sample ball", (ub_of(gt) < 4.0001) and not (gt < 1 - 1e-30), gt.str(5))
    lines.append(f"ENGINE GATE OVERALL: {'PASS' if ok_all else 'FAIL - NOTHING BANKS'}")
    out = "\n".join(lines)
    print(out, flush=True)
    with open(f"{OUT}/CONTROLS.txt", "a") as f:
        f.write(out + "\n")
    return 0 if ok_all else 1

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["gate", "run"])
    ap.add_argument("--topo", type=int, default=0)
    ap.add_argument("--nt", type=int, default=318)
    ap.add_argument("--rel", type=float, default=1e-8)
    ap.add_argument("--prec", type=int, default=128)
    ap.add_argument("--out", type=str, default=None)
    a = ap.parse_args()
    if a.cmd == "gate":
        sys.exit(gate())
    topo = TOPOLOGIES[a.topo]
    print(f"[run] topo {a.topo} = {topo_name(topo)}  Nt={a.nt} rel={a.rel} prec={a.prec}", flush=True)
    row = certified_topology(topo, a.nt, a.rel, a.prec)
    row["topo_index"] = a.topo
    out = a.out or f"{OUT}/rows/row_t{a.topo:02d}_N{a.nt}.json"
    os.makedirs(os.path.dirname(out), exist_ok=True)
    json.dump(row, open(out, "w"), indent=1)
    print(json.dumps(row, indent=1), flush=True)

if __name__ == "__main__":
    main()

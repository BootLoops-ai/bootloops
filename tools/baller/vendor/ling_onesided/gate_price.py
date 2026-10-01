# baller engine — priced-lever gate layer over taylor_p.
#!/usr/bin/env python3
"""gate_price.py — WP-G Stage-1: full gate-object pricing runs.

Architecture (both arms):
  Z = Z_core (paneled in p, quadrature in s) + tails (certified sup bounds).

  p-layer: tiling of [pl,ph] by panels (center m, half-width r, order K).
  Per panel, order-K Taylor-Lagrange in p makes the integrand a degree-(K-1)
  polynomial + ball remainder; an n-point Gauss-Legendre rule (n>=K/2, nodes
  and weights as certified balls) integrates the polynomial part EXACTLY, so

    int_panel F dp  in  sum_j W_j F(p_j)  +/-  Mfac * sup-ish |f_K(P)|
    Mfac = 2 r^{K+1}/(K+1) + sum_j W_j (r|x_j|)^K      (both terms exact)

  and the s-quadrature is pulled outside: panel main value =
  nested acb quadrature over s of prior * sum_j W_j F(p_j, s)   (point-p
  assemblies only -> no p dependency blowup); panel remainder =
  Mfac * upper-bound integral of |f_K(p-ball; s)| via an adaptive grid
  (f_K from the truncated-series assembly on the full p-ball).

  Tails (p<pl, p>ph, s outside [0,shi]^d): recursive box engine; per box the
  16 pattern polynomials and P0000 are evaluated on (p-box, y-boxes) directly
  (small expressions -> mild dependency), giving certified sup of the
  318-product (corr: times (1 - inf P0000)^{-318} with a pair-bound fallback
  1-P0000 >= p(1+(1-p)(1-(y1y2y3)^2))); box bound * exact prior mass, refined
  by a worst-first heap. y-parametrization covers s -> inf exactly (y in
  [0, ymax]) — the memo's y-ball sup tail bound.

  unc arm: s3 closed exactly inside the panels (full closure weights);
  tail boxes carry s3 as a free y3 in [0,1] with prior mass 1.

Controls: T1/T2 anchors through this module's own value paths (gate 1e-30);
GL-panel vs mpmath high-order reference containment; split-panel overlap;
planted-error 0001:+1 end-to-end through both panel paths.
"""
import argparse, heapq, json, math, os, resource, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from flint import arb, acb, arb_poly, ctx
import taylor_p as tp
from taylor_p import (D, S, frac, grouped, mutate_counts, binomrow, s_var,
                      GROUPS, COUNTS, NTOT, LAM)

RESULTS = os.environ.get("WPG_RESULTS", os.path.join(os.getcwd(), "results"))
QUAD_PREC = int(os.environ.get("QUAD_PREC", "48"))
COLLAPSE_PREC = int(os.environ.get("COLLAPSE_PREC", "192"))
CALLS = [0]

def prior(s, lam=LAM):
    return lam * (-(s * lam)).exp()

def rss_mb():
    return round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 1)

def bank(tag, row):
    os.makedirs(RESULTS, exist_ok=True)
    json.dump(row, open(os.path.join(RESULTS, f"{tag}.json"), "w"), indent=1)
    print(json.dumps(row), flush=True)


def progress_hook(tag):
    """Per-job progress receipts (long runs otherwise die at timeouts with
    no rows). Writes {tag}.progress every 60 s when WPG_PROGRESS set."""
    if not os.environ.get("WPG_PROGRESS"):
        return
    import threading, time as _t
    t0 = _t.time()
    os.makedirs(RESULTS, exist_ok=True)

    def loop():
        while True:
            _t.sleep(60)
            try:
                open(os.path.join(RESULTS, f"{tag}.progress"), "w").write(
                    json.dumps(dict(evals=CALLS[0],
                                    secs=round(_t.time() - t0, 1),
                                    maxrss_mb=rss_mb())))
            except Exception:
                pass
    threading.Thread(target=loop, daemon=True).start()

def relw(b):
    b = b.real if isinstance(b, acb) else b
    return float("inf") if b.mid() == 0 else float(arb(b.rad()) / abs(arb(b.mid())))

def hull(a, b):
    h = arb(a).union(arb(b))
    assert h.contains(arb(a)) and h.contains(arb(b))
    return h

def fatq(z, thr=arb("1e-6")):
    return bool(arb(z.real.rad()) > thr) or bool(arb(z.imag.rad()) > thr)

def upper(x):
    """Certified upper bound of an arb ball as arb."""
    return arb(x.mid()) + arb(x.rad())

def lower(x):
    return arb(x.mid()) - arb(x.rad())

# ------------------------------------------------- plain value-only integrands
def L_corr_val(p, s1, s2, s3, counts=None, corrected=True):
    """Plain arb/acb corrected 4-dim integrand (bench_stage pattern)."""
    CALLS[0] += 1
    counts = COUNTS if counts is None else counts
    one = p * 0 + 1
    pi0 = one - p
    beta = one / (2 * pi0 * p)
    y1, y2, y3 = (-beta * s1).exp(), (-beta * s2).exp(), (-beta * s3).exp()
    P1, P2, P3 = _Pm(pi0, p, y1), _Pm(pi0, p, y2), _Pm(pi0, p, y3)
    P12, P123 = _Pm(pi0, p, y1 * y2), _Pm(pi0, p, y1 * y2 * y3)
    pi = [pi0, p]
    def pat(x):
        t = p * 0
        for r in (0, 1):
            for w in (0, 1):
                for v in (0, 1):
                    t += (pi[r] * P3[r][w] * P2[w][v] * P1[v][x[0]] * P1[v][x[1]]
                          * P12[w][x[2]] * P123[r][x[3]])
        return t
    val = one
    ntot = 0
    for pstr, n in counts.items():
        if n == 0:
            continue
        ntot += n
        val *= pat(tuple(int(ch) for ch in pstr)) ** n
    if corrected:
        val *= (one - pat((0, 0, 0, 0))) ** (-ntot)
    return val

def _Pm(pi0, p, y):
    return [[pi0 + p * y, p - p * y], [pi0 - pi0 * y, p + pi0 * y]]

def G_unc_pt(p, s1, s2, wv, groups=GROUPS):
    """Point-p collapsed s3-closed value; wv precomputed at this p (192-bit)."""
    CALLS[0] += 1
    old = ctx.prec
    ctx.prec = COLLAPSE_PREC
    try:
        (V, _), _ = tp.upoly_pair(D(p), s1, s2, groups, deriv=False)
        tot = V[0] * wv[0]
        for j in range(1, V.length()):
            tot += V[j] * wv[j]
        return tot
    finally:
        ctx.prec = old

def wv_at(p, J=NTOT):
    old = ctx.prec
    ctx.prec = COLLAPSE_PREC
    try:
        one = p * 0 + 1
        beta = one / (2 * (one - p) * p)
        return [LAM / (LAM + (2 * j) * beta) for j in range(J + 1)]
    finally:
        ctx.prec = old

# ----------------------------------------------- fast series assembly (unc fK)
def G_unc_series_fast(pS, s1, s2, groups=GROUPS):
    """y3-collapsed s3-closed object as order-K series in p; the 318-fold
    u-polynomial product uses K+1 arb_poly components (flint fast poly mult)
    instead of scalar-series convolutions."""
    old = ctx.prec
    ctx.prec = COLLAPSE_PREC
    try:
        K = len(pS.c) - 1
        ac, beta = tp.pattern_ac_D(pS, s1, s2, groups)   # S-series a,c per group
        facs = []
        for (n, a, c) in ac:
            bn = binomrow(n)
            one = pS._cs(1)
            apow, cpow = [one], [one]
            for _ in range(n):
                apow.append(apow[-1] * a)
                cpow.append(cpow[-1] * c)
            vco = [bn[j] * apow[n - j] * cpow[j] for j in range(n + 1)]  # S list
            comps = [arb_poly([vco[j].c[m] for j in range(n + 1)])
                     for m in range(K + 1)]
            facs.append(comps)
        def upmul(A, B):
            z = arb_poly([])
            out = [z for _ in range(K + 1)]
            for i, Ai in enumerate(A):
                for j in range(K + 1 - i):
                    out[i + j] = out[i + j] + Ai * B[j]
            return out
        facs.sort(key=lambda q: q[0].length())
        while len(facs) > 1:
            facs.sort(key=lambda q: q[0].length())
            facs.append(upmul(facs.pop(0), facs.pop(0)))
        C = facs[0]                       # C[m] = arb_poly in u, m-th t-order
        J = C[0].length() - 1
        wv = [(pS._cs(LAM)) / (beta * (2 * j) + LAM) for j in range(J + 1)]
        # only coefficient K is needed downstream: direct contraction
        cK = C[0][0] * 0
        for j in range(J + 1):
            for m in range(K + 1):
                cm = C[m][j]
                if cm.mid() == 0 and cm.rad() == 0:
                    continue
                cK += cm * wv[j].c[K - m]
        return cK
    finally:
        ctx.prec = old

# ------------------------------------------------------ certified GL rules
_GLCACHE = {}
GL_GUARD = os.environ.get("GL_GUARD", "1e-30")
GL_DPS = int(os.environ.get("GL_DPS", "50"))


def gl_rule(n, prec=256):
    """n-point Gauss-Legendre nodes/weights on [-1,1] as certified arb balls.
    mpmath root guess +/- GL_GUARD guard ball; containment certified by
    interval Newton (P_n(ball) contains 0, P_n'(ball) excludes 0).
    GL_GUARD/GL_DPS env-tunable (default 1e-30/50 = the registered probe
    values): the ONESIDED V2 diagnosis measured the double-collapse assembly
    amplifying the node guard width to an ABSOLUTE noise floor ~1e-366 at
    live-region panels (taylor_p convolution blowup on any p-width) — the
    guard, not working precision, set the floor. Tighter guards need
    GL_DPS ~ guard_digits + 20 and prec >= 3.5 * guard_digits bits."""
    if n in _GLCACHE:
        return _GLCACHE[n]
    import mpmath as mp
    mp.mp.dps = GL_DPS
    old = ctx.prec
    ctx.prec = prec
    def Pn(x, n):
        p0, p1 = x * 0 + 1, x
        for k in range(1, n):
            p0, p1 = p1, ((2 * k + 1) * x * p1 - k * p0) / (k + 1)
        return p1 if n >= 1 else p0
    def dPn(x, n):
        # (1-x^2) P_n' = n (P_{n-1} - x P_n)
        return n * (Pn(x, n - 1) - x * Pn(x, n)) / (1 - x * x)
    nodes, weights = [], []
    for k in range(1, n + 1):
        g = mp.mpf(math.cos(math.pi * (k - 0.25) / (n + 0.5)))
        for _ in range(60):
            g = g - mp.legendre(n, g) / mp.diff(lambda t: mp.legendre(n, t), g)
        xb = arb(str(g), GL_GUARD)
        assert Pn(xb, n).contains(arb(0)), (n, k)
        assert not dPn(xb, n).contains(arb(0)), (n, k)
        w = 2 / ((1 - xb * xb) * dPn(xb, n) ** 2)
        nodes.append(xb)
        weights.append(w)
    ctx.prec = old
    _GLCACHE[n] = (nodes, weights)
    return nodes, weights

def panel_setup(pm, pr, K):
    """Panel nodes p_j (arb balls), weights W_j = r*w_j, and the remainder
    moment factor Mfac = 2 r^{K+1}/(K+1) + sum_j W_j (r|x_j|)^K."""
    n = max(2, (K + 1) // 2)
    xs, ws = gl_rule(n)
    m, r = arb(pm), arb(pr)
    pjs = [m + r * x for x in xs]
    Wjs = [r * w for w in ws]
    Mfac = 2 * r ** (K + 1) / (K + 1)
    for w, x in zip(ws, xs):
        Mfac += (r * w) * (r * abs(x)) ** K
    return pjs, Wjs, Mfac

# --------------------------------------------- adaptive upper-bound grid (fK)
def sup_integral(fmag, dims, budget, maxev, prec=COLLAPSE_PREC):
    """Upper bound of int_region prod(prior) * fmag(s...) via worst-first
    refinement. dims = [(lo,hi),...]; fmag(cells)->arb upper bound on the box.
    Returns (bound_arb, evals)."""
    old = ctx.prec
    ctx.prec = prec
    try:
        def mass(c):
            m = arb(1)
            for (a, b) in c:
                m *= (-(arb(a) * LAM)).exp() - (-(arb(b) * LAM)).exp()
            return m
        def score(c):
            f = fmag([(arb(a).union(arb(b))) for (a, b) in c])
            return upper(f) * mass(c)
        def key(s):
            return -float(arb(s.mid()).log()) if bool(s > 0) else 1e308
        heap = []
        cnt = [0]
        root = tuple((float(a), float(b)) for (a, b) in dims)
        s0 = score(list(root))
        heapq.heappush(heap, (key(s0), 0, root, s0))
        ev = 1
        while ev < maxev:
            tot = arb(0)
            for x in heap:
                tot += x[3]
            if bool(tot < budget):
                break
            negs, _, cell, sc = heapq.heappop(heap)
            widths = [(math.log(b + 1e-9) - math.log(a + 1e-9) if a > 0
                       else b - a) for (a, b) in cell]
            d = widths.index(max(widths))
            a, b = cell[d]
            mid = 0.5 * (a + b)
            for half in ((a, mid), (mid, b)):
                nc = list(cell)
                nc[d] = half
                s = score(nc)
                cnt[0] += 1
                heapq.heappush(heap, (key(s), cnt[0], tuple(nc), s))
                ev += 1
        tot = arb(0)
        for x in heap:
            tot += x[3]
        return tot, ev
    finally:
        ctx.prec = old

# ------------------------------------------------------------- panel runners
def run_uncpanel(a):
    ctx.prec = QUAD_PREC
    progress_hook(a.tag)
    t0 = time.time()
    K = a.K
    groups = grouped(mutate_counts(a.mutate)) if a.mutate else GROUPS
    pjs, Wjs, Mfac = panel_setup(a.pm, a.pr, K)
    Jg = sum(n for _, n in groups)
    wvs = [wv_at(pj, Jg) for pj in pjs]
    lo, hi = 0.0, a.shi
    s1h = hull(lo, hi)
    FLOOR = arb(a.floor)
    def inner_sum(s1, s2):
        t = acb(0)
        for pj, Wj, wv in zip(pjs, Wjs, wvs):
            t += Wj * G_unc_pt(pj, s1, s2, wv, groups)
        return t
    def f(x2, _):
        if fatq(x2):
            return (hi - lo) * inner_sum(s1h, x2) * prior(s1h) * prior(x2)
        return acb.integral(
            lambda x1, __: inner_sum(x1, x2) * prior(x1),
            arb(lo), arb(hi), rel_tol=arb(a.toli),
            abs_tol=FLOOR * arb("1e-2")) * prior(x2)
    J = acb.integral(f, arb(lo), arb(hi), rel_tol=arb(a.tolo), abs_tol=FLOOR)
    tmain = time.time() - t0
    # remainder
    ctx.prec = COLLAPSE_PREC
    pball = arb(arb(a.pm).mid(), a.pr)
    def fmag(cells):
        return abs(G_unc_series_fast(s_var(pball, K), cells[0], cells[1], groups))
    rbud = arb(a.floor)
    U, rev = sup_integral(fmag, [(lo, hi), (lo, hi)], rbud / max(Mfac, arb("1e-300")),
                          a.maxrev)
    R = upper(Mfac * U)
    row = dict(kind="uncpanel", pm=a.pm, pr=a.pr, K=K, shi=a.shi,
               J_mid=arb(J.real.mid()).str(20), J_rad=arb(J.real.rad()).str(5, radius=False),
               R=R.str(5, radius=False), Mfac=float(Mfac.mid()),
               relw_J=relw(J), floor=a.floor, mutate=a.mutate or None,
               evals=CALLS[0], rem_evals=rev,
               secs_main=round(tmain, 1), secs=round(time.time() - t0, 1),
               maxrss_mb=rss_mb())
    bank(a.tag, row)
    return 0

def run_corrpanel(a):
    ctx.prec = QUAD_PREC
    t0 = time.time()
    K = a.K
    counts = mutate_counts(a.mutate) if a.mutate else COUNTS
    pjs, Wjs, Mfac = panel_setup(a.pm, a.pr, K)
    lo, hi = 0.0, a.shi
    sh = hull(lo, hi)
    FLOOR = arb(a.floor)
    def inner_sum(s1, s2, s3):
        t = acb(0)
        for pj, Wj in zip(pjs, Wjs):
            t += Wj * L_corr_val(pj, s1, s2, s3, counts)
        return t
    def g(x3, _):
        if fatq(x3):
            return ((hi - lo) ** 2 * inner_sum(sh, sh, x3)
                    * prior(sh) ** 2 * prior(x3))
        def f(x2, __):
            if fatq(x2):
                return (hi - lo) * inner_sum(sh, x2, x3) * prior(sh) * prior(x2)
            return acb.integral(
                lambda x1, ___: inner_sum(x1, x2, x3) * prior(x1),
                arb(lo), arb(hi), rel_tol=arb(a.toli),
                abs_tol=FLOOR * arb("1e-4")) * prior(x2)
        return acb.integral(f, arb(lo), arb(hi), rel_tol=arb(a.tolm),
                            abs_tol=FLOOR * arb("1e-2")) * prior(x3)
    J = acb.integral(g, arb(lo), arb(hi), rel_tol=arb(a.tolo), abs_tol=FLOOR)
    tmain = time.time() - t0
    ctx.prec = QUAD_PREC
    pball = arb(arb(a.pm).mid(), a.pr)
    def fmag(cells):
        sb = tp.L_prod_S(s_var(pball, K), cells[0], cells[1], cells[2],
                         counts, corrected=True)
        return abs(sb.c[K])
    rbud = arb(a.floor)
    U, rev = sup_integral(fmag, [(lo, hi)] * 3, rbud / max(Mfac, arb("1e-300")),
                          a.maxrev, prec=QUAD_PREC)
    R = upper(Mfac * U)
    row = dict(kind="corrpanel", pm=a.pm, pr=a.pr, K=K, shi=a.shi,
               J_mid=arb(J.real.mid()).str(20), J_rad=arb(J.real.rad()).str(5, radius=False),
               R=R.str(5, radius=False), Mfac=float(Mfac.mid()),
               relw_J=relw(J), floor=a.floor, mutate=a.mutate or None,
               evals=CALLS[0], rem_evals=rev,
               secs_main=round(tmain, 1), secs=round(time.time() - t0, 1),
               maxrss_mb=rss_mb())
    bank(a.tag, row)
    return 0

# ------------------------------------------------------------------ controls
def run_controls(a):
    ok = True
    t00 = time.time()
    ctx.prec = 384
    print("== controls: anchors through gate_price value paths ==")
    for name, ((aa, bb, cc, pp), corr50, unc50) in tp.ANCHORS.items():
        s1, s2, s3, p = frac(aa), frac(bb), frac(cc), frac(pp)
        lc = L_corr_val(p, s1, s2, s3).log()
        d1 = abs(lc - arb(corr50))
        h1 = bool(d1 < arb("1e-30"))
        print(f"{'OK  ' if h1 else 'FAIL'} {name} corr anchor |dlog| < {d1.str(3)}")
        ok = ok and h1
        # unc s3-closed value: G_unc_pt (production path) vs taylor_p dual path
        wv = wv_at(p)
        g1 = G_unc_pt(p, s1, s2, wv)
        g2 = tp.G_unc_D(D(p, p * 0 + 1), s1, s2).v
        d2 = abs(g1.log() - g2.log())
        # T1 full-closure contraction carries ~1e-31 relrad at 192-bit (both
        # paths identically); gate midpoint agreement + containment of 0.
        h2 = bool(abs(arb(d2.mid())) < arb("1e-30")) and d2.contains(arb(0))
        print(f"{'OK  ' if h2 else 'FAIL'} {name} unc G_unc_pt vs G_unc_D |dlog| < {d2.str(3)}")
        ok = ok and h2
    # GL panel vs mpmath reference + split consistency, thin-s fiber, corr arm
    print("\n== GL panel machinery: containment + split overlap (corr fiber) ==")
    ctx.prec = QUAD_PREC
    s1, s2, s3 = arb("0.1"), arb("0.2"), arb("0.3")
    pm, pr, K = 0.345, 1e-3, 16
    pjs, Wjs, Mfac = panel_setup(pm, pr, K)
    JS = acb(0)
    for pj, Wj in zip(pjs, Wjs):
        JS += Wj * L_corr_val(pj, s1, s2, s3)
    ctx.prec = QUAD_PREC
    pball = arb(arb(pm).mid(), pr)
    sb = tp.L_prod_S(s_var(pball, K), s1, s2, s3, corrected=True)
    E = upper(Mfac * abs(sb.c[K]))
    encl = JS.real.union(JS.real + E).union(JS.real - E)
    import mpmath as mp
    mp.mp.dps = 40
    def fmp(x):
        v = L_corr_val(arb(str(x)), s1, s2, s3)
        return mp.mpf(arb(v.mid()).str(30, radius=False))
    ref = mp.quad(fmp, [pm - pr, pm + pr])
    inref = encl.contains(arb(str(ref)))
    print(f"{'OK  ' if inref else 'FAIL'} mpmath reference inside panel enclosure "
          f"(relw {relw(encl):.2e})")
    ok = ok and inref
    halves = []
    for (m2, r2) in ((pm - pr / 2, pr / 2), (pm + pr / 2, pr / 2)):
        pj2, Wj2, Mf2 = panel_setup(m2, r2, K)
        J2 = acb(0)
        for pj, Wj in zip(pj2, Wj2):
            J2 += Wj * L_corr_val(pj, s1, s2, s3)
        sb2 = tp.L_prod_S(s_var(arb(arb(m2).mid(), r2), K), s1, s2, s3,
                          corrected=True)
        E2 = upper(Mf2 * abs(sb2.c[K]))
        halves.append(J2.real.union(J2.real + E2).union(J2.real - E2))
    hsum = halves[0] + halves[1]
    over = bool(hsum.overlaps(encl))
    print(f"{'OK  ' if over else 'FAIL'} split-panel overlap")
    ok = ok and over
    print(f"\nCONTROLS {'PASS' if ok else 'FAIL'} ({time.time()-t00:.1f}s)")
    return 0 if ok else 1

# ---------------------------------------------------------------------- cli
def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("controls")
    for nm in ("uncpanel", "corrpanel"):
        q = sub.add_parser(nm)
        q.add_argument("--pm", type=float, required=True)
        q.add_argument("--pr", type=float, required=True)
        q.add_argument("--K", type=int, required=True)
        q.add_argument("--shi", type=float, default=2.0)
        q.add_argument("--floor", required=True)
        q.add_argument("--tolo", default="3e-2")
        q.add_argument("--tolm", default="1e-2")
        q.add_argument("--toli", default="3e-3")
        q.add_argument("--maxrev", type=int, default=240)
        q.add_argument("--mutate", default=None)
        q.add_argument("--tag", required=True)
    a = ap.parse_args()
    return dict(controls=run_controls, uncpanel=run_uncpanel,
                corrpanel=run_corrpanel)[a.cmd](a)

if __name__ == "__main__":
    sys.exit(main())

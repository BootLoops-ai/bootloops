# baller engine — certified sup-bound tail integrals.
#!/usr/bin/env python3
"""tail_engine.py — certified sup-bound tail integrals for both gate arms.

Bound over a box (p-cell x s-cells): the 14 observed pattern polynomials and
P0000 are interval-evaluated on (p-hull, y-hulls) — small fixed expressions,
mild dependency — giving

  unc:  sup Prod pat_x^{n_x}        <= Prod min(1, sup pat_x)^{n_x}
  corr: sup Prod (pat_x/(1-P0000))^{n_x}
        <= Prod min(1, sup pat_x / inf(1-P0000))^{n_x}
     (pointwise pat_x <= 1-P0000 for every observed pattern — the event is a
      subset of {>=1 present} — so each ratio is capped at 1: corner-safe)
  inf(1-P0000) = max(1 - sup P0000,  inf p*(1+(1-p)(1-(y1y2y3)^2)))
     (second term: P0000 <= P(x_A=0,x_D=0) = (1-p) - p(1-p)(1-y_AD),
      y_AD = (y1 y2 y3)^2 the A-D path).

y from s via endpoint monotonicity of beta(p)=1/(2p(1-p)) (no interval
division): s in [a,b] (b may be inf), p-cell [pa,pb] ->
  y_hi = exp(-beta_lo * a)  (=1 at a=0),  y_lo = exp(-beta_hi * b) (=0 at inf).

All cell bounds kept as arb (exponents ~ e^-800 are fine); heap ranked by
float(log(bound)); totals are arb sums; the target is arb.

unc s3: panels close s3 exactly, so unc tail boxes carry s3 = (0, inf)
(prior mass 1, y3 free in [0,1]).
"""
import heapq, json, math, os, time
from flint import arb, ctx

LAM = 10

def upper(x):
    return arb(x.mid()) + arb(x.rad())

def lower(x):
    return arb(x.mid()) - arb(x.rad())

def hull(a, b):
    h = arb(a).union(arb(b))
    return h

def _Pm(pi0, p, y):
    return [[pi0 + p * y, p - p * y], [pi0 - pi0 * y, p + pi0 * y]]

def beta_range(pa, pb):
    """certified [beta_lo, beta_hi] over p in [pa,pb] (0<=pa<pb<=1)."""
    g = lambda p: 2 * arb(p) * (1 - arb(p))
    gvals = [g(pa), g(pb)]
    gmax = upper(gvals[0]).max(upper(gvals[1]))
    if pa <= 0.5 <= pb:
        gmax = arb(0.5)
    gmin = lower(gvals[0]).min(lower(gvals[1]))
    blo = 1 / gmax
    bhi = None if bool(gmin <= arb(0)) else 1 / gmin   # None == +inf
    return blo, bhi

def ybox(pa, pb, a, b):
    blo, bhi = beta_range(pa, pb)
    yhi = arb(1) if a == 0 else upper((-(lower(blo) * a)).exp())
    if math.isinf(b) or (bhi is None and b > 0):
        ylo = arb(0)
    else:
        ylo = lower((-(upper(bhi) * b)).exp())
    return hull(max(0.0, float(ylo)), min(1.0, float(yhi.max(arb(0)))))

def cell_logbound(arm, pc, sc, counts, NT):
    """log of certified upper bound of the integrand sup on the cell,
    plus log prior mass. Returns arb log-bound (or None if bound is 0)."""
    pa, pb = max(pc[0], 0.0), min(pc[1], 1.0)
    P = hull(pa, pb)
    ys = [ybox(pa, pb, a, b) for (a, b) in sc]
    one = P * 0 + 1
    pi0 = one - P
    P1, P2, P3 = _Pm(pi0, P, ys[0]), _Pm(pi0, P, ys[1]), _Pm(pi0, P, ys[2])
    P12, P123 = _Pm(pi0, P, ys[0] * ys[1]), _Pm(pi0, P, ys[0] * ys[1] * ys[2])
    pi = [pi0, P]
    def pat(x):
        t = P * 0
        for r in (0, 1):
            for w in (0, 1):
                for v in (0, 1):
                    t += (pi[r] * P3[r][w] * P2[w][v] * P1[v][x[0]] * P1[v][x[1]]
                          * P12[w][x[2]] * P123[r][x[3]])
        return t
    if arm == "corr":
        p0000u = upper(pat((0, 0, 0, 0)))
        l1 = 1 - p0000u
        ypath = (ys[0] * ys[1] * ys[2]) ** 2
        l2 = lower(P) * (1 + (1 - upper(P)) * (1 - upper(ypath)))
        l2 = lower(l2)
        inf1mp = l1.max(l2)
    lg = arb(0)
    for ps, n in counts.items():
        if n == 0:
            continue
        u = upper(pat(tuple(int(ch) for ch in ps)))
        if bool(u <= arb(0)):
            return None
        if arm == "corr":
            if bool(inf1mp > arb(0)):
                ratio = u / inf1mp
                r1 = ratio.min(arb(1))
            else:
                r1 = arb(1)
            lg += n * r1.log()
        else:
            lg += n * u.min(arb(1)).log()
    # log prior mass
    lmass = arb(pb - pa).log() if pb > pa else arb("-1e9")
    for (a, b) in sc:
        if a == 0 and math.isinf(b):
            continue
        mb = arb(0) if math.isinf(b) else (-(arb(b) * LAM)).exp()
        mm = (-(arb(a) * LAM)).exp() - mb
        if bool(mm <= arb(0)):
            return None
        lmass += mm.log()
    return lg + lmass

def tail_bound(arm, pbox, sboxes, log_target, counts, maxit=200000, prec=96):
    """Refine boxes until arb-sum of bounds <= exp(log_target).
    Returns (log10_bound, iters, ok)."""
    old = ctx.prec
    ctx.prec = prec
    NT = sum(counts.values())
    target = arb(log_target).exp()
    try:
        heap = []   # (-float(log), counter, pc, sc, arb_bound)
        cnt = [0]
        def push(pc, sc):
            lb = cell_logbound(arm, pc, sc, counts, NT)
            b = arb(0) if lb is None else lb.exp()
            key = -1e308 if lb is None else -float(arb(lb.mid()))
            cnt[0] += 1
            heapq.heappush(heap, (key, cnt[0], tuple(pc), tuple(map(tuple, sc)), b))
        push(list(pbox), [list(x) for x in sboxes])
        it = 1
        while it < maxit:
            tot = arb(0)
            for x in heap:
                tot += x[4]
            if bool(tot <= target):
                break
            key, _, pc, sc, b = heapq.heappop(heap)
            pc, sc = list(pc), [list(x) for x in sc]
            # split widest-in-effect dim
            spans = []
            pa, pb = pc
            spans.append(math.log((pb + 1e-14) / (pa + 1e-14)))
            for (a, bb) in sc:
                spans.append(3.0 if math.isinf(bb) else
                             math.log((bb + 1e-2 / LAM) / (a + 1e-2 / LAM)))
            d = spans.index(max(spans))
            if d == 0:
                mid = math.sqrt(pa * pb) if pa > 0 else pb / 4
                if not (pa < mid < pb):
                    mid = 0.5 * (pa + pb)
                push([pa, mid], sc)
                push([mid, pb], sc)
            else:
                a, bb = sc[d - 1]
                mid = a + math.log(2) / LAM if math.isinf(bb) else 0.5 * (a + bb)
                for half in ((a, mid), (mid, bb)):
                    nsc = [list(x) for x in sc]
                    nsc[d - 1] = list(half)
                    push(pc, nsc)
            it += 2
        tot = arb(0)
        for x in heap:
            tot += x[4]
        ok = bool(tot <= target)
        l10 = float(upper(tot).log() / arb(10).log()) if bool(tot > 0) else -9e9
        return l10, it, ok
    finally:
        ctx.prec = old

def run(arm, pl, ph, shi, log10_target_per_region, counts, maxit, tag,
        results_dir):
    t0 = time.time()
    inf = math.inf
    sfull, swin, sin_ = (0.0, inf), (shi, inf), (0.0, shi)
    if arm == "unc":
        regions = [("p_lo", (0.0, pl), [sfull, sfull, sfull]),
                   ("p_hi", (ph, 1.0), [sfull, sfull, sfull]),
                   ("s1_out", (pl, ph), [swin, sfull, sfull]),
                   ("s2_out", (pl, ph), [sin_, swin, sfull])]
    else:
        regions = [("p_lo", (0.0, pl), [sfull, sfull, sfull]),
                   ("p_hi", (ph, 1.0), [sfull, sfull, sfull]),
                   ("s1_out", (pl, ph), [swin, sfull, sfull]),
                   ("s2_out", (pl, ph), [sin_, swin, sfull]),
                   ("s3_out", (pl, ph), [sin_, sin_, swin])]
    lt = float(log10_target_per_region) * math.log(10)
    rows, oks = {}, True
    tot_l10 = -9e9
    for name, pbox, sboxes in regions:
        l10, it, ok = tail_bound(arm, pbox, sboxes, lt, counts, maxit=maxit)
        rows[name] = dict(log10_bound=l10, iters=it, ok=ok)
        print(f"{name}: log10 bound {l10:.2f} (target {log10_target_per_region})"
              f" iters {it} ok={ok}", flush=True)
        tot_l10 = max(tot_l10, l10) + math.log10(1 + 10 ** (min(tot_l10, l10) - max(tot_l10, l10))) \
            if tot_l10 > -8e9 else l10
        oks = oks and ok
    row = dict(kind="tails", arm=arm, pl=pl, ph=ph, shi=shi,
               log10_target_per_region=log10_target_per_region,
               regions=rows, total_log10=tot_l10, all_ok=oks,
               secs=round(time.time() - t0, 1))
    os.makedirs(results_dir, exist_ok=True)
    json.dump(row, open(os.path.join(results_dir, f"{tag}.json"), "w"), indent=1)
    print(json.dumps(row), flush=True)
    return 0 if oks else 1

if __name__ == "__main__":
    import argparse, sys
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", required=True)
    ap.add_argument("--pl", type=float, required=True)
    ap.add_argument("--ph", type=float, required=True)
    ap.add_argument("--shi", type=float, default=2.0)
    ap.add_argument("--l10target", type=float, required=True)
    ap.add_argument("--maxit", type=int, default=200000)
    ap.add_argument("--tag", required=True)
    a = ap.parse_args()
    HERE = os.path.dirname(os.path.abspath(__file__))
    counts = json.load(open(os.path.join(
        os.environ.get("WPG_BASE", HERE), "QUARTET_COLLAPSE.json")))["pattern_counts"]
    sys.exit(run(a.arm, a.pl, a.ph, a.shi, a.l10target, counts, a.maxit,
                 a.tag, os.environ.get("WPG_RESULTS", ".")))

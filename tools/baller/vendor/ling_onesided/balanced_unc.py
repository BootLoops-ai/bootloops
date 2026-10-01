# baller engine — balanced 2-dim quadrature support (double_unc/l5_center dep).
#!/usr/bin/env python3
"""balanced_unc.py — WP-G balanced-shape ((A,B),(C,D)) UNCORRECTED collapse
integrand, product-domain convention (PINNED_SPEC referee pin 3).

CONVENTION (pinned): cherry heights u1 (AB), u2 (CD) iid Exp(lam); root
increment s3 ~ Exp(lam) measured from the OLDER cherry; root height
h = max(u1,u2) + s3. Ranked-increment convention NOT used.

DOMAIN SPLIT + EDGE STRUCTURE. Region 1 (u1 >= u2), coordinates
t = u2, d = u1 - u2, s3 (all in [0,inf), Jacobian 1, prior density
lam^2 e^{-lam(2t+d)} * lam e^{-lam s3}); with yt = exp(-beta t),
yd = exp(-beta d), y3 = exp(-beta s3):

    edge         time            y
    A, B leaf    u1 = t+d        yt*yd
    C, D leaf    u2 = t          yt
    root->AB     s3              y3          (AB cherry is the older)
    root->CD     s3 + d          y3*yd

    pat(x) = sum_{r,v,w} pi_r P_{y3}[r][v] P_{y3 yd}[r][w]
             * P_{yt yd}[v][xA] P_{yt yd}[v][xB] * P_{yt}[w][xC] P_{yt}[w][xD]

Summing r with stationarity (sum_r pi_r d_{rv} = 0,
sum_r pi_r d_{rv} d_{rw} = (-1)^{v+w} pi0 pi1) kills the odd-y3 term:

    pat = a + c * y3^2   EXACTLY, with
    a = (pi0 QA0 + p QA1)(pi0 QC0 + p QC1)
    c = pi0 p yd (QA0 - QA1)(QC0 - QC1)
    QA_v = P_{yt yd}[v][xA] P_{yt yd}[v][xB],  QC_w = P_{yt}[w][xC] P_{yt}[w][xD]

(verified symbolically in `selftest`). Region 2 (u2 > u1) is the mirror:
the same form with pattern slots (xA,xB) <-> (xC,xD). The 318-character
product per region is a degree-318 polynomial in u = y3^2 (319 coeffs,
same as caterpillar); s3 closes EXACTLY with the identical weights
w_j = lam/(lam + 2j*beta) (windowed variant identical to collapse_unc.py).
Grouping: AB-swap AND CD-swap invariance -> 7 distinct factors per region
(exponents sum to 318). Evidence at fixed p:

  Z_bal(p) = int_0^inf int_0^inf lam^2 e^{-2 lam t} e^{-lam d}
             [ contract(C^{(1)}(t,d), w) + contract(C^{(2)}(t,d), w) ] dt dd

Subcommands: selftest | timeeval | ladder | slice2d | balpanel.
balpanel = eras p-panel (order-K GL + Taylor-Lagrange remainder in p,
gate_price.py architecture) over the balanced 2-dim (t,d) quadrature.
"""
import argparse, importlib.util, json, math, os, resource, sys, time
from flint import arb, acb, arb_poly, acb_poly, ctx

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import taylor_p as tp
from taylor_p import S, frac, binomrow, s_var, mutate_counts

RESULTS = os.environ.get("WPG_RESULTS",
                         os.path.join(os.getcwd(), "probe_results"))
COUNTS = json.load(open(os.environ.get(
    "WPG_COUNTS", os.path.join(HERE, "QUARTET_COLLAPSE.json"))))["pattern_counts"]
LAM = 10
QUAD_PREC = int(os.environ.get("QUAD_PREC", "48"))
COLLAPSE_PREC = int(os.environ.get("COLLAPSE_PREC", "192"))
CALLS = [0]


def grouped_bal(counts, swap=False):
    """Merge patterns by AB-swap x CD-swap symmetry. swap=True builds the
    region-2 (CD-cherry-older) grouping: slots (xA,xB) <-> (xC,xD) first."""
    g = {}
    for pat, n in counts.items():
        if n == 0:
            continue
        x = tuple(int(ch) for ch in pat)
        if swap:
            x = (x[2], x[3], x[0], x[1])
        key = (tuple(sorted(x[:2])), tuple(sorted(x[2:])))
        g[key] = g.get(key, 0) + n
    return sorted(g.items())


G1 = grouped_bal(COUNTS, swap=False)   # region 1: AB cherry older
G2 = grouped_bal(COUNTS, swap=True)    # region 2: CD cherry older
NTOT = sum(n for _, n in G1)           # 318


def _Pmat(pi0, p, y):
    return [[pi0 + p * y, p - p * y], [pi0 - pi0 * y, p + pi0 * y]]


def _lift(p, x):
    """Coerce a scalar/str to p's arithmetic type (arb/acb/D/S)."""
    if hasattr(p, "_c"):
        return p._c(frac(x)) if isinstance(x, str) else p._c(x)
    return p * 0 + frac(x)


def pattern_ac_bal(p, t, d, groups):
    """Per grouped pattern: (n, a, c) with pat = a + c*y3^2 exactly.
    Generic: p may be arb/acb (plain), tp.D (dual) or tp.S (series)."""
    one = _lift(p, 1)
    pi0 = one - p
    beta = one / ((pi0 * p) * 2)
    tt, dd_ = _lift(p, t), _lift(p, d)
    yt = (-(beta * tt)).exp()
    yd = (-(beta * dd_)).exp()
    ytd = yt * yd
    PA, PC = _Pmat(pi0, p, ytd), _Pmat(pi0, p, yt)
    cpref = pi0 * p * yd
    out = []
    for (ab, cd), n in groups:
        xA, xB = ab
        xC, xD = cd
        QA0 = PA[0][xA] * PA[0][xB]
        QA1 = PA[1][xA] * PA[1][xB]
        QC0 = PC[0][xC] * PC[0][xD]
        QC1 = PC[1][xC] * PC[1][xD]
        a = (pi0 * QA0 + p * QA1) * (pi0 * QC0 + p * QC1)
        c = cpref * (QA0 - QA1) * (QC0 - QC1)
        out.append((n, a, c))
    return out, beta


def upoly_bal(p, t, d, groups):
    """319-coefficient polynomial in u = y3^2 for one region arm."""
    ac, beta = pattern_ac_bal(p, t, d, groups)
    Poly = acb_poly if isinstance(ac[0][1], acb) else arb_poly
    facs = []
    for (n, a, c) in ac:
        bn = binomrow(n)
        apow = [a * 0 + 1]
        cpow = [a * 0 + 1]
        for _ in range(n):
            apow.append(apow[-1] * a)
            cpow.append(cpow[-1] * c)
        facs.append(Poly([bn[j] * apow[n - j] * cpow[j] for j in range(n + 1)]))
    while len(facs) > 1:
        facs.sort(key=lambda q: q.length())
        facs.append(facs.pop(0) * facs.pop(0))
    return facs[0], beta


def wvec_full(beta, J, lam=LAM):
    return [lam / (lam + (2 * j) * beta) for j in range(J + 1)]


def wvec_win(beta, J, A, B, lam=LAM):
    A = frac(A)
    Bv = None if B in (None, "inf") else frac(B)
    out = []
    for j in range(J + 1):
        r = lam + (2 * j) * beta
        v = (-r * A).exp()
        if Bv is not None:
            v = v - (-r * Bv).exp()
        out.append(lam * v / r)
    return out


def F_bal(p, t, d, wv, g1=G1, g2=G2):
    """Both-region s3-closed integrand at fixed p (COLLAPSE_PREC assembly).
    Does NOT include the (t,d) prior lam^2 e^{-lam(2t+d)}."""
    CALLS[0] += 1
    old = ctx.prec
    ctx.prec = COLLAPSE_PREC
    try:
        tot = None
        for g in (g1, g2):
            poly, _ = upoly_bal(p, t, d, g)
            s = poly[0] * wv[0]
            for j in range(1, poly.length()):
                s += poly[j] * wv[j]
            tot = s if tot is None else tot + s
        return tot
    finally:
        ctx.prec = old


def prior_td(t, d, lam=LAM):
    return (lam * lam) * (-(2 * lam) * t).exp() * (-lam * d).exp()


# --------------------------- direct (uncollapsed) generic balanced evaluator
def bal_pattern_labeled(pi1, u1, u2, s3, x, leaves=(0, 1, 2, 3)):
    """P(x) on the labeled balanced topology ((L1,L2),(L3,L4)): cherry
    (L1,L2) at height u1, cherry (L3,L4) at u2, root at max(u1,u2)+s3.
    Region logic explicit — works for either ordering. Independent of the
    collapse path (plain peeling)."""
    one = pi1 * 0 + 1
    pi0 = one - pi1
    beta = one / (2 * pi0 * pi1)
    hi = u1 if bool(arb(u1) >= arb(u2)) else u2
    bA = hi + s3 - u1          # root -> (L1,L2)-node branch time
    bC = hi + s3 - u2          # root -> (L3,L4)-node branch time
    P_A = _Pmat(pi0, pi1, (-beta * u1).exp())
    P_C = _Pmat(pi0, pi1, (-beta * u2).exp())
    P_rA = _Pmat(pi0, pi1, (-beta * bA).exp())
    P_rC = _Pmat(pi0, pi1, (-beta * bC).exp())
    pi = [pi0, pi1]
    x1, x2, x3, x4 = (x[i] for i in leaves)
    tot = pi1 * 0
    for r in (0, 1):
        for v in (0, 1):
            for w in (0, 1):
                tot += (pi[r] * P_rA[r][v] * P_rC[r][w]
                        * P_A[v][x1] * P_A[v][x2] * P_C[w][x3] * P_C[w][x4])
    return tot


def L_bal_direct(p, u1, u2, s3, counts=None):
    """Direct 318-character product on ((A,B),(C,D)) — cross-check path."""
    counts = COUNTS if counts is None else counts
    val = p * 0 + 1
    cache = {}
    one = p * 0 + 1
    pi0 = one - p
    beta = one / (2 * pi0 * p)
    hi = u1 if bool(arb(u1) >= arb(u2)) else u2
    P_A = _Pmat(pi0, p, (-beta * u1).exp())
    P_C = _Pmat(pi0, p, (-beta * u2).exp())
    P_rA = _Pmat(pi0, p, (-beta * (hi + s3 - u1)).exp())
    P_rC = _Pmat(pi0, p, (-beta * (hi + s3 - u2)).exp())
    pi = [pi0, p]
    for ps, n in counts.items():
        if n == 0:
            continue
        x = tuple(int(ch) for ch in ps)
        t = p * 0
        for r in (0, 1):
            for v in (0, 1):
                for w in (0, 1):
                    t += (pi[r] * P_rA[r][v] * P_rC[r][w]
                          * P_A[v][x[0]] * P_A[v][x[1]] * P_C[w][x[2]] * P_C[w][x[3]])
        val *= t ** n
    return val


# ------------------------------------------- series-in-p (eras remainder)
def G_bal_series_fast(pS, t, d, g1=G1, g2=G2):
    """Balanced s3-closed both-region object as order-K series in p
    (gate_price.G_unc_series_fast architecture: K+1 arb_poly components).
    Returns coefficient K only (the Taylor-Lagrange remainder term)."""
    old = ctx.prec
    ctx.prec = COLLAPSE_PREC
    try:
        K = len(pS.c) - 1
        cK_tot = None
        for groups in (g1, g2):
            ac, beta = pattern_ac_bal(pS, t, d, groups)
            facs = []
            for (n, a, c) in ac:
                bn = binomrow(n)
                one = pS._cs(1)
                apow, cpow = [one], [one]
                for _ in range(n):
                    apow.append(apow[-1] * a)
                    cpow.append(cpow[-1] * c)
                vco = [bn[j] * apow[n - j] * cpow[j] for j in range(n + 1)]
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

            while len(facs) > 1:
                facs.sort(key=lambda q: q[0].length())
                facs.append(upmul(facs.pop(0), facs.pop(0)))
            C = facs[0]
            J = C[0].length() - 1
            wv = [(pS._cs(LAM)) / (beta * (2 * j) + LAM) for j in range(J + 1)]
            cK = C[0][0] * 0
            for j in range(J + 1):
                for m in range(K + 1):
                    cm = C[m][j]
                    if cm.mid() == 0 and cm.rad() == 0:
                        continue
                    cK += cm * wv[j].c[K - m]
            cK_tot = cK if cK_tot is None else cK_tot + cK
        return cK_tot
    finally:
        ctx.prec = old


# ------------------------------------------------------------------ helpers
def upper(x):
    return arb(x.mid()) + arb(x.rad())


def relw(b):
    b = b.real if isinstance(b, acb) else b
    return float("inf") if b.mid() == 0 else float(arb(b.rad()) / abs(arb(b.mid())))


def hull(a, b):
    h = arb(a).union(arb(b))
    assert h.contains(arb(a)) and h.contains(arb(b))
    return h


def fatq(z, thr=arb("1e-6")):
    return bool(arb(z.real.rad()) > thr) or bool(arb(z.imag.rad()) > thr)


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
    import threading
    t0 = time.time()
    os.makedirs(RESULTS, exist_ok=True)

    def loop():
        while True:
            time.sleep(60)
            try:
                open(os.path.join(RESULTS, f"{tag}.progress"), "w").write(
                    json.dumps(dict(evals=CALLS[0],
                                    secs=round(time.time() - t0, 1),
                                    maxrss_mb=rss_mb())))
            except Exception:
                pass
    threading.Thread(target=loop, daemon=True).start()


def _load(modname, path):
    spec = importlib.util.spec_from_file_location(modname, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


# ------------------------------------------------------------------ selftest
def contains_zero(b):
    return bool(b.rad() >= abs(arb(b.mid())))


PATTERNS = [(a, b, c, d) for a in (0, 1) for b in (0, 1)
            for c in (0, 1) for d in (0, 1)]


def selftest():
    ok = True
    t00 = time.time()
    print("== balanced_unc selftest ==", flush=True)

    # (0) exact sympy structural check over QQ(p): parity, closed forms,
    #     support, sum-to-one — region-1 tree
    import sympy as sp
    YT, YD, Y3, P = sp.symbols("yt yd y3 p")
    pi0s = 1 - P
    Ps = lambda y: [[pi0s + P * y, P - P * y], [pi0s - pi0s * y, P + pi0s * y]]
    PAs, PCs = Ps(YT * YD), Ps(YT)
    PrA, PrC = Ps(Y3), Ps(Y3 * YD)
    pis = [pi0s, P]
    support = set()
    tot16 = 0
    for xi in range(16):
        x = tuple(int(b) for b in format(xi, "04b"))
        pat = sum(pis[r] * PrA[r][v] * PrC[r][w]
                  * PAs[v][x[0]] * PAs[v][x[1]] * PCs[w][x[2]] * PCs[w][x[3]]
                  for r in (0, 1) for v in (0, 1) for w in (0, 1))
        pat = sp.expand(pat)
        tot16 += pat
        if sp.Poly(pat, Y3).coeff_monomial(Y3) != 0:
            print(f"FAIL parity: pattern {x} has nonzero y3^1 coefficient")
            ok = False
        QA = [PAs[v][x[0]] * PAs[v][x[1]] for v in (0, 1)]
        QC = [PCs[w][x[2]] * PCs[w][x[3]] for w in (0, 1)]
        a_cf = (pi0s * QA[0] + P * QA[1]) * (pi0s * QC[0] + P * QC[1])
        c_cf = pi0s * P * YD * (QA[0] - QA[1]) * (QC[0] - QC[1])
        if sp.expand(pat - a_cf - c_cf * Y3**2) != 0:
            print(f"FAIL closed-form (a,c): pattern {x}")
            ok = False
        for mono in sp.Poly(pat, YT, YD, Y3).monoms():
            support.add(mono)
    if sp.simplify(tot16 - 1) != 0:
        print("FAIL sum-to-one (symbolic)")
        ok = False
    print("OK   parity: all 16 patterns even in y3; pat = a + c*y3^2 exact")
    print("OK   sum-to-one (symbolic, all 16 patterns)")
    print(f"OK   support in (yt,yd,y3): {len(support)} monomials = {sorted(support)}")

    # (a) collapsed vs direct product at generic ball points, both regions
    ctx.prec = 384
    print("\n-- (a) collapsed vs direct, prec 384 --")
    pts = [("R1", "0.22", "0.55", "0.20", "0.41"),   # u1>u2
           ("R1b", "1/3", "0.31", "0.09", "0.13"),
           ("R2", "0.61", "0.12", "0.83", "0.27"),   # u2>u1
           ("R2b", "2/5", "0.05", "0.95", "0.55")]
    for name, pp, su1, su2, ss3 in pts:
        p, u1, u2, s3 = frac(pp), frac(su1), frac(su2), frac(ss3)
        if bool(u1 >= u2):
            t, d, g = u2, u1 - u2, G1
        else:
            t, d, g = u1, u2 - u1, G2
        poly, beta = upoly_bal(p, t, d, g)
        hit0 = poly.length() == NTOT + 1
        u = ((-beta * s3).exp()) ** 2
        lv = poly(u).log()
        ld = L_bal_direct(p, u1, u2, s3).log()
        dl = abs(lv - ld)
        hit = hit0 and bool(dl < arb("1e-70"))
        print(f"{'OK  ' if hit else 'FAIL'} {name}: |dlog collapsed-direct| < {dl.str(3)} "
              f"(deg {poly.length()-1})")
        ok = ok and hit
    # region-boundary consistency: arm1(t,0) == arm2(t,0)
    p, t = frac("0.3"), frac("0.4")
    pol1, beta = upoly_bal(p, t, arb(0), G1)
    pol2, _ = upoly_bal(p, t, arb(0), G2)
    u = ((-beta * frac("0.17")).exp()) ** 2
    db = abs(pol1(u).log() - pol2(u).log())
    hit = bool(db < arb("1e-70"))
    print(f"{'OK  ' if hit else 'FAIL'} region boundary d=0: |dlog arm1-arm2| < {db.str(3)}")
    ok = ok and hit

    # (b) degeneration controls vs jg0 (independent code path)
    print("\n-- (b) degeneration: star (u1=u2, s3=0) + caterpillar-labeled (s3=0) --")
    jg0 = _load("jg0", os.path.join(os.path.dirname(HERE), "2a_prep", "jg0_controls.py"))
    ctx.prec = 256
    for hstr, pstr in (("3/10", "1/3"), ("1/2", "2/5")):
        h, p = frac(hstr), frac(pstr)
        bad = sum(1 for x in PATTERNS if not contains_zero(
            bal_pattern_labeled(p, h, h, arb(0), x) - jg0.star_pattern_prob(p, h, x)))
        badc = sum(1 for x in PATTERNS if not contains_zero(
            bal_pattern_labeled(p, h, h, arb(0), x)
            - jg0.pattern_prob(p, h, arb(0), arb(0), x)))
        hit = bad == 0 and badc == 0
        print(f"{'OK  ' if hit else 'FAIL'} star h={hstr} p={pstr}: vs star_pattern_prob "
              f"{16-bad}/16, vs caterpillar(h,0,0) {16-badc}/16 contain 0")
        ok = ok and hit
    # balanced(u1,u2,s3=0), u1>u2  ==  labeled caterpillar (((C,D),A),B) at
    # (s1=u2, s2=u1-u2, s3=0) — jg0's independently written evaluator
    p, u1, u2 = frac("0.28"), frac("0.5"), frac("0.2")
    bad = sum(1 for x in PATTERNS if not contains_zero(
        bal_pattern_labeled(p, u1, u2, arb(0), x)
        - jg0.cat_pattern_prob_labeled(p, u2, u1 - u2, arb(0), x, (2, 3, 0, 1))))
    hit = bad == 0
    print(f"{'OK  ' if hit else 'FAIL'} s3=0 vs jg0 labeled caterpillar (C,D,A,B): "
          f"{16-bad}/16 contain 0")
    ok = ok and hit
    # sum-to-one, numeric, generic point both regions
    for (su1, su2, ss3) in (("0.5", "0.2", "0.3"), ("0.1", "0.7", "0.05")):
        tot = arb(0)
        for x in PATTERNS:
            tot += bal_pattern_labeled(frac("0.37"), frac(su1), frac(su2), frac(ss3), x)
        hit = contains_zero(tot - arb(1))
        print(f"{'OK  ' if hit else 'FAIL'} sum-to-one at (u1,u2,s3)=({su1},{su2},{ss3})")
        ok = ok and hit

    # s3-closure: exact weights vs direct 1-dim s3 quadrature on a fiber
    ctx.prec = 384
    print("\n-- s3-fiber closure: exact vs direct quadrature, prec 384 --")
    for name, pp, st, sd in (("F1", "1/3", "0.2", "0.15"), ("F2", "0.22", "0.05", "0.6")):
        p, t, d = frac(pp), frac(st), frac(sd)
        poly1, beta = upoly_bal(p, t, d, G1)
        poly2, _ = upoly_bal(p, t, d, G2)
        for (A, B) in [("0", "2"), ("0.28", "0.32")]:
            wv = wvec_win(beta, NTOT, A, B)
            closed = poly1[0] * wv[0] + poly2[0] * wv[0]
            for j in range(1, NTOT + 1):
                closed += (poly1[j] + poly2[j]) * wv[j]
            u1a, u2a = t + d, t          # region-1 point
            u1b, u2b = t, t + d          # region-2 mirror point
            direct = acb.integral(
                lambda x, f: LAM * (-LAM * x).exp()
                * (L_bal_direct(p, u1a, u2a, x) + L_bal_direct(p, u1b, u2b, x)),
                arb(A), arb(B), rel_tol=arb("1e-40"), abs_tol=arb("1e-800"))
            dl = abs(closed.log() - direct.real.log())
            hit = bool(dl < arb("1e-28"))
            print(f"{'OK  ' if hit else 'FAIL'} {name} s3 in [{A},{B}]: |dlog| < {dl.str(3)}")
            ok = ok and hit

    # (c) planted error 0001:+1 through the new path (fiber level)
    print("\n-- (c) planted error 0001:+1, prec 384 --")
    mut = mutate_counts("0001:+1")
    g1m, g2m = grouped_bal(mut, False), grouped_bal(mut, True)
    p, t, d, s3 = frac("1/3"), frac("0.2"), frac("0.15"), frac("0.3")
    pol0, beta = upoly_bal(p, t, d, G1)
    polm, _ = upoly_bal(p, t, d, g1m)
    u = ((-beta * s3).exp()) ** 2
    shc = polm(u).log() - pol0(u).log()
    shd = (L_bal_direct(p, t + d, t, s3, mut).log()
           - L_bal_direct(p, t + d, t, s3).log())
    eq = bool(abs(shc - shd) < arb("1e-40")) and bool(abs(shc) > arb("1e-2"))
    print(f"{'OK  ' if eq else 'FAIL'} fiber shift equality: collapsed {shc.str(20)} "
          f"vs direct {shd.str(20)}, detectable")
    ok = ok and eq
    # and through the closed+2-arm object (wv long enough for the +1 plant)
    wv = wvec_full(beta, NTOT + 1)
    F0 = F_bal(p, t, d, wv)
    Fm = F_bal(p, t, d, wv, g1m, g2m)
    sh = Fm.log() - F0.log()
    fired = bool(not sh.contains(arb(0))) and bool(abs(arb(sh.mid())) > arb("1e-2"))
    print(f"{'OK  ' if fired else 'FAIL'} s3-closed both-arm shift {sh.str(10)} excludes 0")
    ok = ok and fired

    # (d) leaf-relabel symmetry within the balanced class (2A lemma numeric)
    print("\n-- (d) balanced-class relabel symmetry: B1 -> B2 under sigma=(B C) --")
    ctx.prec = 256
    B1, B2 = (0, 1, 2, 3), (0, 2, 1, 3)     # ((A,B),(C,D)), ((A,C),(B,D))
    sigma = (0, 2, 1, 3)                    # data-slot transposition B<->C
    pt = (frac("2/5"), frac("1/2"), frac("1/4"), frac("1/8"))  # PINNED T2 point
    p, u1, u2, s3 = pt
    bad = sum(1 for x in PATTERNS if not contains_zero(
        bal_pattern_labeled(p, u1, u2, s3, x, B2)
        - bal_pattern_labeled(p, u1, u2, s3, tuple(x[i] for i in sigma), B1)))
    hit = bad == 0
    print(f"{'OK  ' if hit else 'FAIL'} group action B1->B2: {16-bad}/16 contain 0")
    ok = ok and hit
    # cherry invariances on the production arm (AB swap, CD swap)
    bad = sum(1 for x in PATTERNS if not contains_zero(
        bal_pattern_labeled(p, u1, u2, s3, x)
        - bal_pattern_labeled(p, u1, u2, s3, (x[1], x[0], x[3], x[2]))))
    print(f"{'OK  ' if bad == 0 else 'FAIL'} AB+CD cherry invariance: {16-bad}/16")
    ok = ok and bad == 0
    # region-swap identity: ((A,B),(C,D)) at (u1,u2) == ((C,D),(A,B)) at (u2,u1)
    bad = sum(1 for x in PATTERNS if not contains_zero(
        bal_pattern_labeled(p, u1, u2, s3, x)
        - bal_pattern_labeled(p, u2, u1, s3, x, (2, 3, 0, 1))))
    print(f"{'OK  ' if bad == 0 else 'FAIL'} region-swap identity: {16-bad}/16")
    ok = ok and bad == 0
    # m3-analog: WRONG sigma (A<->B, a B1 self-symmetry) must FAIL
    wrong = (1, 0, 2, 3)
    badw = sum(1 for x in PATTERNS if not contains_zero(
        bal_pattern_labeled(p, u1, u2, s3, x, B2)
        - bal_pattern_labeled(p, u1, u2, s3, tuple(x[i] for i in wrong), B1)))
    fired = badw > 0
    print(f"{'OK  ' if fired else 'FAIL'} broken-relabel control FIRES: "
          f"{badw}/16 balls exclude 0 (must be > 0)")
    ok = ok and fired

    print(f"\nSELFTEST {'PASS' if ok else 'FAIL'} "
          f"({time.time()-t00:.1f}s, maxrss {rss_mb()} MB)", flush=True)
    return 0 if ok else 1


# ------------------------------------------------------------------ timeeval
def timeeval():
    ctx.prec = COLLAPSE_PREC
    p = arb(1) / 3
    one = arb(1)
    beta = one / (2 * (one - p) * p)
    wv = wvec_full(beta, NTOT)
    t, d = arb("0.2"), arb("0.15")
    t0 = time.time()
    n = 0
    while time.time() - t0 < 3.0:
        F_bal(p, t, d, wv)
        n += 1
    print(f"F_bal[both arms] @ {COLLAPSE_PREC}bit: {(time.time()-t0)/n*1e3:.2f} ms/eval")
    for K in (8, 16, 24, 32):
        t0 = time.time()
        G_bal_series_fast(s_var(arb(p.mid(), 2e-3), K), t, d)
        print(f"G_bal_series_fast K={K}: {time.time()-t0:.2f} s/eval")
    return 0


# -------------------------------------------------------------------- ladder
def ladder(args):
    """relwidth-by-K for the balanced panel object (K/r law transfer check)."""
    ctx.prec = COLLAPSE_PREC
    t, d = arb("0.2"), arb("0.15")
    rows = []
    for spec in args.cells.split(";"):
        p0s, rs, Ks = spec.split(",")
        pball = arb(arb(frac(p0s)).mid(), float(rs))
        rw = {}
        for K in (int(k) for k in Ks.split("/")):
            t0 = time.time()
            m = arb(pball.mid())
            Sm = tp.taylor_ladder(lambda pS: _gser(pS, t, d), pball, [K])[K]
            rw[K] = relw(Sm)
            print(f"p0={p0s} r={rs} K={K}: relw {rw[K]:.3e} ({time.time()-t0:.1f}s)",
                  flush=True)
        rows.append(dict(p0=p0s, r=rs, rw={str(k): v for k, v in rw.items()}))
    bank(args.tag, dict(kind="bal-ladder", rows=rows))
    return 0


def _gser(pS, t, d):
    """Full series object (all K+1 coefficients) for taylor_ladder — slower
    scalar-series path, only used by the ladder scan."""
    ac1, beta = pattern_ac_bal(pS, t, d, G1)
    ac2, _ = pattern_ac_bal(pS, t, d, G2)
    one = pS._cs(1)
    tot = None
    for ac in (ac1, ac2):
        polys = []
        for (n, a, c) in ac:
            bn = binomrow(n)
            apow, cpow = [one], [one]
            for _ in range(n):
                apow.append(apow[-1] * a)
                cpow.append(cpow[-1] * c)
            polys.append([bn[j] * apow[n - j] * cpow[j] for j in range(n + 1)])

        def pmul(A, B):
            z = A[0]._cs(0)
            out = [z] * (len(A) + len(B) - 1)
            for i, ai in enumerate(A):
                for j, bj in enumerate(B):
                    out[i + j] = out[i + j] + ai * bj
            return out

        while len(polys) > 1:
            polys.sort(key=len)
            polys.append(pmul(polys.pop(0), polys.pop(0)))
        C = polys[0]
        J = len(C) - 1
        wv = [(pS._cs(LAM)) / (beta * (2 * j) + LAM) for j in range(J + 1)]
        s = C[0] * wv[0]
        for j in range(1, J + 1):
            s = s + C[j] * wv[j]
        tot = s if tot is None else tot + s
    return tot


# ------------------------------------------------------------------- slice2d
def run_slice2d(a):
    """Fixed-p 2-dim (t,d) balanced evidence (full s3 closure) — scale keying
    + planted-error-at-quadrature-level when --mutate is set."""
    ctx.prec = QUAD_PREC
    progress_hook(a.tag)
    p = frac(a.p)
    g1 = grouped_bal(mutate_counts(a.mutate), False) if a.mutate else G1
    g2 = grouped_bal(mutate_counts(a.mutate), True) if a.mutate else G2
    ctx.prec = COLLAPSE_PREC
    one = arb(1)
    beta = one / (2 * (one - p) * p)
    wv = wvec_full(beta, sum(n for _, n in g1))
    ctx.prec = QUAD_PREC
    tlo, thi = (frac(x) for x in a.twin.split(","))
    dlo, dhi = (frac(x) for x in a.dwin.split(","))
    FLOOR = arb(a.floor)
    th = hull(tlo, thi)
    t0 = time.time()

    def f(x2, _):    # x2 = d (outer), inner = t
        if fatq(x2):
            return ((thi - tlo) * F_bal(p, th, x2, wv, g1, g2)
                    * prior_td(th, x2))
        return acb.integral(
            lambda x1, __: F_bal(p, x1, x2, wv, g1, g2) * prior_td(x1, x2),
            tlo, thi, rel_tol=arb(a.toli), abs_tol=FLOOR * arb("1e-2"))
    v = acb.integral(f, dlo, dhi, rel_tol=arb(a.tolo), abs_tol=FLOOR)
    row = dict(kind="bal-slice2d", p=a.p, twin=a.twin, dwin=a.dwin,
               tolo=a.tolo, toli=a.toli, floor=a.floor, mutate=a.mutate or None,
               value_mid=arb(v.real.mid()).str(20),
               value_rad=arb(v.real.rad()).str(5, radius=False),
               relwidth=relw(v),
               logZ=(v.real.log().str(20) if bool(v.real > 0) else "n/a"),
               evals=CALLS[0], secs=round(time.time() - t0, 1),
               maxrss_mb=rss_mb())
    bank(a.tag, row)
    return 0


# ------------------------------------------------------------------ balpanel
def run_balpanel(a):
    """eras p-panel for the balanced shape: order-K GL nodes (point-p
    2-dim (t,d) quadrature of the collapsed object) + Taylor-Lagrange
    remainder via order-K series on the full p-ball (gate_price law)."""
    import gate_price as gp
    ctx.prec = QUAD_PREC
    progress_hook(a.tag)
    t0 = time.time()
    K = a.K
    if a.mutate:
        mc = mutate_counts(a.mutate)
        g1, g2 = grouped_bal(mc, False), grouped_bal(mc, True)
    else:
        g1, g2 = G1, G2
    pjs, Wjs, Mfac = gp.panel_setup(a.pm, a.pr, K)
    Jg = sum(n for _, n in g1)
    wvs = []
    old = ctx.prec
    ctx.prec = COLLAPSE_PREC
    for pj in pjs:
        one = pj * 0 + 1
        beta = one / (2 * (one - pj) * pj)
        wvs.append(wvec_full(beta, Jg))
    ctx.prec = old
    FLOOR = arb(a.floor)
    tlo, thi_ = 0.0, a.thi
    dlo, dhi_ = 0.0, a.dhi
    th = hull(tlo, thi_)

    def inner_sum(t, d):
        s = acb(0)
        for pj, Wj, wv in zip(pjs, Wjs, wvs):
            s += Wj * F_bal(pj, t, d, wv, g1, g2)
        return s

    def f(x2, _):    # x2 = d outer, t inner
        if fatq(x2):
            return ((thi_ - tlo) * inner_sum(th, x2) * prior_td(th, x2))
        return acb.integral(
            lambda x1, __: inner_sum(x1, x2) * prior_td(x1, x2),
            arb(tlo), arb(thi_), rel_tol=arb(a.toli),
            abs_tol=FLOOR * arb("1e-2"))
    J = acb.integral(f, arb(dlo), arb(dhi_), rel_tol=arb(a.tolo), abs_tol=FLOOR)
    tmain = time.time() - t0
    # remainder: Mfac * int prior * |f_K(p-ball; t,d)| via worst-first grid
    ctx.prec = COLLAPSE_PREC
    pball = arb(arb(a.pm).mid(), a.pr)

    def fmag(cells):
        return abs(G_bal_series_fast(s_var(pball, K), cells[0], cells[1], g1, g2))
    U, rev = _sup_integral_bal(fmag, [(tlo, thi_), (dlo, dhi_)],
                               FLOOR / max(Mfac, arb("1e-300")), a.maxrev)
    R = upper(Mfac * U)
    row = dict(kind="balpanel", pm=a.pm, pr=a.pr, K=K, thi=a.thi, dhi=a.dhi,
               J_mid=arb(J.real.mid()).str(20),
               J_rad=arb(J.real.rad()).str(5, radius=False),
               R=R.str(5, radius=False), Mfac=float(Mfac.mid()),
               relw_J=relw(J), floor=a.floor, mutate=a.mutate or None,
               evals=CALLS[0], rem_evals=rev,
               secs_main=round(tmain, 1), secs=round(time.time() - t0, 1),
               maxrss_mb=rss_mb())
    bank(a.tag, row)
    return 0


def _sup_integral_bal(fmag, dims, budget, maxev):
    """Worst-first upper bound of int prior_td * fmag over (t,d) cells;
    exact per-cell prior mass with rates (2*LAM, LAM)."""
    import heapq
    old = ctx.prec
    ctx.prec = COLLAPSE_PREC
    rates = [2 * LAM, LAM]
    try:
        def mass(c):
            m = arb(1)
            for (aa, bb), rt in zip(c, rates):
                m *= (-(arb(aa) * rt)).exp() - (-(arb(bb) * rt)).exp()
            return m

        def score(c):
            f = fmag([(arb(aa).union(arb(bb))) for (aa, bb) in c])
            return upper(f) * mass(c)

        def key(s):
            return -float(arb(s.mid()).log()) if bool(s > 0) else 1e308
        heap = []
        cnt = [0]
        root = tuple((float(aa), float(bb)) for (aa, bb) in dims)
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
            widths = [(math.log(bb + 1e-9) - math.log(aa + 1e-9) if aa > 0
                       else bb - aa) for (aa, bb) in cell]
            di = widths.index(max(widths))
            aa, bb = cell[di]
            mid = 0.5 * (aa + bb)
            for half in ((aa, mid), (mid, bb)):
                nc = list(cell)
                nc[di] = half
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


# ----------------------------------------------------------------------- cli
def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("selftest")
    sub.add_parser("timeeval")
    q = sub.add_parser("ladder")
    q.add_argument("--cells", required=True,
                   help="p0,r,K1/K2/..;p0,r,Ks e.g. 0.345,2e-3,16/24/32")
    q.add_argument("--tag", required=True)
    q = sub.add_parser("slice2d")
    q.add_argument("--p", default="1/3")
    q.add_argument("--twin", default="0,1.75")
    q.add_argument("--dwin", default="0,3.5")
    q.add_argument("--tolo", default="3e-2")
    q.add_argument("--toli", default="1e-2")
    q.add_argument("--floor", required=True)
    q.add_argument("--mutate", default=None)
    q.add_argument("--tag", required=True)
    q = sub.add_parser("balpanel")
    q.add_argument("--pm", type=float, required=True)
    q.add_argument("--pr", type=float, required=True)
    q.add_argument("--K", type=int, required=True)
    q.add_argument("--thi", type=float, default=1.75)
    q.add_argument("--dhi", type=float, default=3.5)
    q.add_argument("--floor", required=True)
    q.add_argument("--tolo", default="3e-2")
    q.add_argument("--toli", default="1e-2")
    q.add_argument("--maxrev", type=int, default=120)
    q.add_argument("--mutate", default=None)
    q.add_argument("--tag", required=True)
    a = ap.parse_args()
    if a.cmd == "selftest":
        return selftest()
    if a.cmd == "timeeval":
        return timeeval()
    return dict(ladder=ladder, slice2d=run_slice2d,
                balpanel=run_balpanel)[a.cmd](a)


if __name__ == "__main__":
    sys.exit(main())

# baller engine — collapse-uncertainty core with selftest.
#!/usr/bin/env python3
"""WP-G collapse route, UNCORRECTED arm — y3-collapsed evidence integrand.

Per COLLAPSE_MEMO.md section 6 (y3-expansion route). Each caterpillar pattern
probability, written in y_k = exp(-beta*s_k), is EXACTLY quadratic-even in y3:

    pat_x(y1,y2,y3;p) = a_x(y1,y2,p) + c_x(y1,y2,p) * y3^2

with the odd coefficient identically zero: grouping the peeling sum over
(root state r, node-w state w),

    pat = sum_{r,w} pi_r (pi_w + d_{rw} y3)(pi_{xD} + d_{r,xD} y1 y2 y3) Q_w,
    d = [[pi1,-pi1],[-pi0,pi0]],  Q_w = P12[w][xC] sum_v P2[w][v] P1[v][xA] P1[v][xB],

and the y3^1 coefficient carries the factors sum_r pi_r d_{r,x} =
pi0*pi1 - pi1*pi0 = 0 (stationarity). This is the memo's "y3-exponent always
even". Closed forms used below (e_{w,x} = sum_r pi_r d_{rw} d_{rx} =
(-1)^{w+x} pi0 pi1):

    a_x = pi_{xD} (pi0 Q_0 + pi1 Q_1)
    c_x = (-1)^{xD} pi0 pi1 y1 y2 (Q_0 - Q_1)

(verified symbolically in `selftest`). The 318-character product is then a
degree-318 polynomial in u = y3^2 (319 even coefficients, the memo's "319
values"), assembled as arb/acb balls at COLLAPSE_PREC (default 192) bits via
per-pattern binomial expansion (each coefficient a single product — no
cancellation inside a coefficient) and a size-sorted product tree. The s3
integral closes EXACTLY (memo): with y3 = exp(-beta*s3),

    int_0^inf lam e^{-lam s} y3^{2j} ds = lam / (lam + 2j*beta)          (full)
    int_A^B  lam e^{-lam s} y3^{2j} ds
        = lam (e^{-(lam+2j beta)A} - e^{-(lam+2j beta)B}) / (lam+2j beta) (windowed)

leaving a 3-dim (s1,s2,pi1) object. (s1,s2) are integrated by the windowed
nested acb quadrature pattern of bench_stage.py at QUAD_PREC (default 48, the
measured p48 lever); pi1 outermost. AB-swap symmetry (patterns with xA,xB
exchanged have identical (a_x,c_x)) merges the 14 observed patterns into 10
distinct factors, exponents summing to 318.

Subcommands: selftest | timeeval | slice2d | direct3d | evidence.
"""
import json, math, os, resource, sys, time
from flint import arb, acb, arb_poly, acb_poly, ctx

# Fixture read root defaults to this module's own directory (the vendored
# QUARTET_COLLAPSE.json); row output root defaults to the CURRENT working
# directory — never the vendor tree (unpinned files there trip baller.verify()).
BASE = os.environ.get("WPG_BASE", os.path.dirname(os.path.abspath(__file__)))
OUT = os.environ.get("COLLAPSE_OUT", os.getcwd())
COUNTS = json.load(open(f"{BASE}/QUARTET_COLLAPSE.json"))["pattern_counts"]
LAM = 10
QUAD_PREC = int(os.environ.get("QUAD_PREC", "48"))
COLLAPSE_PREC = int(os.environ.get("COLLAPSE_PREC", "192"))
CALLS = [0]          # collapsed-integrand evals
DCALLS = [0]         # direct-integrand evals
BINOM = {}

ANCHORS = {  # PINNED_SPEC 50-digit uncorrected logL anchors
    "T1": (("1/10", "2/10", "3/10", "1/3"),
           "-976.09285161889490586499770937640740253545699456512481"),
    "T2": (("1/2", "1/4", "1/8", "2/5"),
           "-842.55310669243368577614288803238346376995933055630806"),
}

def frac(s):
    if isinstance(s, arb):
        return s
    n, _, d = str(s).partition("/")
    return arb(int(n)) / int(d) if d else arb(s)

def binomrow(n):
    if n not in BINOM:
        BINOM[n] = [math.comb(n, j) for j in range(n + 1)]
    return BINOM[n]

def grouped(counts):
    """Merge patterns by AB-swap symmetry: key ((xA,xB) sorted, xC, xD)."""
    g = {}
    for pat, n in counts.items():
        if n == 0:
            continue
        x = tuple(int(ch) for ch in pat)
        key = (tuple(sorted(x[:2])), x[2], x[3])
        g[key] = g.get(key, 0) + n
    return sorted(g.items())

GROUPS = grouped(COUNTS)
NTOT = sum(n for _, n in GROUPS)   # 318

# ---------------------------------------------------------------- integrands
def _Pmat(pi0, p, y):
    return [[pi0 + p * y, p - p * y], [pi0 - pi0 * y, p + pi0 * y]]

def L_direct(p, s1, s2, s3):
    """Direct (uncollapsed) uncorrected likelihood — bench_stage.py pattern."""
    DCALLS[0] += 1
    one = p * 0 + 1
    pi0 = one - p
    beta = one / (2 * pi0 * p)
    y1, y2, y3 = (-beta * s1).exp(), (-beta * s2).exp(), (-beta * s3).exp()
    P1, P2, P3 = _Pmat(pi0, p, y1), _Pmat(pi0, p, y2), _Pmat(pi0, p, y3)
    P12, P123 = _Pmat(pi0, p, y1 * y2), _Pmat(pi0, p, y1 * y2 * y3)
    pi = [pi0, p]
    val = one
    for pat, n in COUNTS.items():
        x = tuple(int(ch) for ch in pat)
        t = p * 0
        for r in (0, 1):
            for w in (0, 1):
                for v in (0, 1):
                    t += (pi[r] * P3[r][w] * P2[w][v] * P1[v][x[0]] * P1[v][x[1]]
                          * P12[w][x[2]] * P123[r][x[3]])
        val *= t ** n
    return val

def pattern_ac(p, s1, s2, groups=GROUPS):
    """Per grouped pattern: (n, a, c) with pat = a + c*y3^2 exactly."""
    one = p * 0 + 1
    pi0 = one - p
    beta = one / (2 * pi0 * p)
    y1, y2 = (-beta * s1).exp(), (-beta * s2).exp()
    y12 = y1 * y2
    P1, P2, P12 = _Pmat(pi0, p, y1), _Pmat(pi0, p, y2), _Pmat(pi0, p, y12)
    cpref = pi0 * p * y12
    out = []
    for (ab, xC, xD), n in groups:
        xA, xB = ab
        Q0 = P12[0][xC] * (P2[0][0] * P1[0][xA] * P1[0][xB]
                           + P2[0][1] * P1[1][xA] * P1[1][xB])
        Q1 = P12[1][xC] * (P2[1][0] * P1[0][xA] * P1[0][xB]
                           + P2[1][1] * P1[1][xA] * P1[1][xB])
        a = (pi0 * Q0 + p * Q1) * (pi0 if xD == 0 else p)
        d = cpref * (Q0 - Q1)
        out.append((n, a, d if xD == 0 else -d))
    return out, beta

def upoly(p, s1, s2, groups=GROUPS):
    """The 319-coefficient polynomial P(u), u = y3^2: the y3-collapse."""
    CALLS[0] += 1
    ac, beta = pattern_ac(p, s1, s2, groups)
    complexq = any(isinstance(v, acb) for v in (p, s1, s2))
    Poly = acb_poly if complexq else arb_poly
    facs = []
    for (n, a, c) in ac:
        if complexq:
            a, c = acb(a), acb(c)
        bn = binomrow(n)
        apow = [a * 0 + 1]
        cpow = [a * 0 + 1]
        for _ in range(n):
            apow.append(apow[-1] * a)
            cpow.append(cpow[-1] * c)
        facs.append(Poly([bn[j] * apow[n - j] * cpow[j] for j in range(n + 1)]))
    while len(facs) > 1:                      # size-sorted product tree
        facs.sort(key=lambda q: q.length())
        facs.append(facs.pop(0) * facs.pop(0))
    return facs[0], beta

def wvec_full(beta, J, lam=LAM):
    """w_j = int_0^inf lam e^{-lam s} y3^{2j} ds = lam/(lam + 2j*beta)."""
    return [lam / (lam + (2 * j) * beta) for j in range(J + 1)]

def wvec_win(beta, J, A, B, lam=LAM):
    """Windowed exact closure over s3 in [A,B] (B=None -> infinity)."""
    A = frac(A)
    Bv = None if B in (None, "inf") else frac(B)
    out = []
    for j in range(J + 1):
        r = lam + (2 * j) * beta
        t = (-r * A).exp()
        if Bv is not None:
            t = t - (-r * Bv).exp()
        out.append(lam * t / r)
    return out

def F_collapsed(p, s1, s2, wv, groups=GROUPS):
    """s3-closed integrand: sum_j C_j(s1,s2,p) * wv[j].  192-bit assembly."""
    old = ctx.prec
    ctx.prec = COLLAPSE_PREC
    try:
        poly, _ = upoly(p, s1, s2, groups)
        tot = poly[0] * wv[0]
        for j in range(1, poly.length()):
            tot += poly[j] * wv[j]
        return tot
    finally:
        ctx.prec = old

# ------------------------------------------------------------------ helpers
def prior(s, lam=LAM):
    return lam * (-lam * s).exp()

def lref(args=None):
    """Magnitude normalizer for abs_tol floors. Default: |L| at the T1 point
    (bench TOL convention, ~1e-424 — ~68 orders below the true 2-dim value
    scale, so the abs floor never activates). Override with --tolref keyed to
    the measured value scale (~1e-361 for the fixed-pi1 collapsed 2-dim) to
    let negligible regions terminate at the floor; rigorous either way (the
    returned enclosure is the claim, tolerances are only targets)."""
    if args is not None and getattr(args, "tolref", None):
        return arb(args.tolref)
    old = ctx.prec
    ctx.prec = 128
    v = L_direct(arb(1) / 3, arb("0.1"), arb("0.2"), arb("0.3"))
    ctx.prec = old
    return arb(v.mid())

def enclosure_row(v, extra):
    m = v.real if isinstance(v, acb) else v
    rw = float(m.rad() / abs(arb(m.mid()))) if m.mid() != 0 else float("inf")
    row = dict(extra)
    row.update(
        value_mid=arb(m.mid()).str(25),
        value_rad=arb(m.rad()).str(5, radius=False),
        relwidth=rw,
        logZ_enclosure=m.log().str(25) if m > 0 else "n/a",
        imag_rad=(arb(v.imag.rad()).str(5, radius=False)
                  if isinstance(v, acb) else "0"),
        quad_prec=QUAD_PREC, collapse_prec=COLLAPSE_PREC,
        maxrss_mb=round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 1),
    )
    return row

def writerow(tag, row):
    os.makedirs(f"{OUT}/COLLAPSE_ROWS", exist_ok=True)
    json.dump(row, open(f"{OUT}/COLLAPSE_ROWS/{tag}.json", "w"), indent=1)
    print(json.dumps(row), flush=True)

def mutate_counts(spec):
    """spec like '0001:+1' -> mutated copy of COUNTS."""
    c = dict(COUNTS)
    pat, _, dv = spec.partition(":")
    c[pat] = c[pat] + int(dv)
    return c

# ---------------------------------------------------------------- selftest
def selftest():
    ok = True
    print("== selftest: y3-collapse machinery ==", flush=True)

    # (0) exact symbolic parity + closed-form (a,c) check, sympy over QQ(p)
    import sympy as sp
    Y1, Y2, Y3, P = sp.symbols("y1 y2 y3 p")
    pi0s = 1 - P
    Ps = lambda y: [[pi0s + P * y, P - P * y], [pi0s - pi0s * y, P + pi0s * y]]
    P1s, P2s, P3s = Ps(Y1), Ps(Y2), Ps(Y3)
    P12s, P123s = Ps(Y1 * Y2), Ps(Y1 * Y2 * Y3)
    pis = [pi0s, P]
    support = set()
    for xi in range(16):
        x = tuple(int(b) for b in format(xi, "04b"))
        pat = sum(pis[r] * P3s[r][w] * P2s[w][v] * P1s[v][x[0]] * P1s[v][x[1]]
                  * P12s[w][x[2]] * P123s[r][x[3]]
                  for r in (0, 1) for w in (0, 1) for v in (0, 1))
        pat = sp.expand(pat)
        podd = sp.Poly(pat, Y3)
        if podd.coeff_monomial(Y3) != 0 and sp.simplify(podd.coeff_monomial(Y3)) != 0:
            print(f"FAIL parity: pattern {x} has nonzero y3^1 coefficient"); ok = False
        # closed-form a,c
        Q = [P12s[w][x[2]] * sum(P2s[w][v] * P1s[v][x[0]] * P1s[v][x[1]]
                                 for v in (0, 1)) for w in (0, 1)]
        a_cf = pis[x[3]] * (pi0s * Q[0] + P * Q[1])
        c_cf = (-1) ** x[3] * pi0s * P * Y1 * Y2 * (Q[0] - Q[1])
        if sp.expand(pat - a_cf - c_cf * Y3**2) != 0:
            print(f"FAIL closed-form (a,c): pattern {x}"); ok = False
        for mono in sp.Poly(pat, Y1, Y2, Y3).monoms():
            support.add(mono)
    memo_support = {(0,0,0),(2,0,0),(2,2,0),(2,2,2),(3,2,0),(3,2,2),(3,3,2),(4,2,2),(4,3,2)}
    print(f"OK   parity: all 16 patterns even in y3; pat = a + c*y3^2 exact")
    print(f"OK   support: union over 16 patterns = {sorted(support)}")
    print(f"     memo 9-monomial support {'MATCHES' if support == memo_support else 'DIFFERS'}")
    if support != memo_support:
        ok = False

    # (a) anchors: collapsed integrand at anchor s3 vs 50-digit logL_unc pins
    ctx.prec = 384
    print("\n== (a) pointwise anchors, prec 384, target 30+ digits ==")
    for name, ((a, b, c, pp), anch) in ANCHORS.items():
        s1, s2, s3, p = frac(a), frac(b), frac(c), frac(pp)
        poly, beta = upoly(p, s1, s2)
        assert poly.length() == NTOT + 1, poly.length()
        u = ((-beta * s3).exp()) ** 2
        lv = poly(u).log()
        diff = abs(lv - arb(anch))
        hit = bool(diff < arb("1e-30")) and bool(arb(lv.rad()) < arb("1e-40"))
        print(f"{'OK  ' if hit else 'FAIL'} {name}: |logL_collapsed - anchor| < {diff.str(3)}, "
              f"ball rad {arb(lv.rad()).str(3, radius=False)} (319 coeffs, deg {poly.length()-1} in u=y3^2)")
        ok = ok and hit

    # direct cross-check of (a,c) closed forms at a non-anchor ball point
    p, s1, s2, s3 = arb("0.22"), arb("0.07"), arb("0.9"), arb("0.41")
    poly, beta = upoly(p, s1, s2)
    u = ((-beta * s3).exp()) ** 2
    d1, d2 = poly(u), L_direct(p, s1, s2, s3)
    agree = bool(abs(d1.log() - d2.log()) < arb("1e-70"))
    print(f"{'OK  ' if agree else 'FAIL'} off-anchor point: collapsed vs direct product, "
          f"|dlog| < {abs(d1.log()-d2.log()).str(3)}")
    ok = ok and agree

    # (closure) exact s3-closure vs direct 1-dim s3 quadrature
    print("\n== s3-fiber closure: exact vs direct quadrature, prec 384 ==")
    pts = [("T1", "1/10", "2/10", "1/3"), ("T2", "1/2", "1/4", "2/5"),
           ("X1", "0.07", "0.9", "0.22"), ("X2", "0.33", "0.05", "0.61")]
    for name, a, b, pp in pts:
        s1, s2, p = frac(a), frac(b), frac(pp)
        poly, beta = upoly(p, s1, s2)
        for (A, B) in [("0", "2"), ("0.28", "0.32")]:
            wv = wvec_win(beta, poly.length() - 1, A, B)
            closed = poly[0] * wv[0]
            for j in range(1, poly.length()):
                closed += poly[j] * wv[j]
            direct = acb.integral(
                lambda x, f: prior(x) * L_direct(p, s1, s2, x),
                arb(A), arb(B), rel_tol=arb("1e-45"), abs_tol=arb("1e-800"))
            dl = abs(closed.log() - direct.real.log())
            hit = bool(dl < arb("1e-30"))
            print(f"{'OK  ' if hit else 'FAIL'} {name} s3 in [{A},{B}]: |dlog| < {dl.str(3)}")
            ok = ok and hit
        # full closure vs [0,20] + rigorous tail bound
        wvf = wvec_full(beta, poly.length() - 1)
        closedf = poly[0] * wvf[0]
        sabs = abs(poly[0])
        for j in range(1, poly.length()):
            closedf += poly[j] * wvf[j]
            sabs += abs(poly[j])
        wv20 = wvec_win(beta, poly.length() - 1, "0", "20")
        c20 = poly[0] * wv20[0]
        for j in range(1, poly.length()):
            c20 += poly[j] * wv20[j]
        tail = sabs * arb(-LAM * 20).exp()   # sum_j |C_j| lam e^{-r_j*20}/r_j <= sum|C_j| e^{-200}
        gap = abs(closedf - c20)
        hit = bool(gap <= tail * arb("1.0000001")) and bool(tail / closedf < arb("1e-30"))
        print(f"{'OK  ' if hit else 'FAIL'} {name} full-closure tail: |full-[0,20]| = {gap.str(3)} "
              f"<= bound {tail.str(3)}, negligible vs value")
        ok = ok and hit

    # (c-lite) planted error at anchor + fiber level (both routes move together)
    print("\n== (c) planted error 0001: 55->56, prec 384 ==")
    mut = mutate_counts("0001:+1")
    gmut = grouped(mut)
    for name, ((a, b, c, pp), anch) in ANCHORS.items():
        s1, s2, s3, p = frac(a), frac(b), frac(c), frac(pp)
        poly, beta = upoly(p, s1, s2, gmut)
        u = ((-beta * s3).exp()) ** 2
        lv = poly(u).log()
        moved = bool(abs(lv - arb(anch)) > arb("1e-3"))
        print(f"{'OK  ' if moved else 'FAIL'} {name}: mutated anchor check moves by "
              f"{(lv - arb(anch)).str(6)} (must be detected)")
        ok = ok and moved
    # fiber, both routes, mutated
    s1, s2, p = frac("1/10"), frac("2/10"), frac("1/3")
    OLD = dict(COUNTS)
    poly0, beta = upoly(p, s1, s2)
    polym, _ = upoly(p, s1, s2, gmut)
    wv = wvec_win(beta, polym.length() - 1, "0", "2")
    wv0 = wvec_win(beta, poly0.length() - 1, "0", "2")
    cm = polym[0] * wv[0]
    for j in range(1, polym.length()):
        cm += polym[j] * wv[j]
    c0 = poly0[0] * wv0[0]
    for j in range(1, poly0.length()):
        c0 += poly0[j] * wv0[j]
    COUNTS.clear(); COUNTS.update(mut)
    dm = acb.integral(lambda x, f: prior(x) * L_direct(p, s1, s2, x),
                      arb(0), arb(2), rel_tol=arb("1e-45"), abs_tol=arb("1e-800"))
    COUNTS.clear(); COUNTS.update(OLD)
    d0 = acb.integral(lambda x, f: prior(x) * L_direct(p, s1, s2, x),
                      arb(0), arb(2), rel_tol=arb("1e-45"), abs_tol=arb("1e-800"))
    shift_c = (cm.log() - c0.log())
    shift_d = (dm.real.log() - d0.real.log())
    agree = bool(abs(shift_c - shift_d) < arb("1e-30"))
    det = bool(abs(shift_c) > arb("1e-2"))
    print(f"{'OK  ' if agree and det else 'FAIL'} fiber shift: collapsed {shift_c.str(20)}, "
          f"direct {shift_d.str(20)} -> move together ({abs(shift_c-shift_d).str(3)}), detectable")
    ok = ok and agree and det

    print(f"\nSELFTEST {'PASS' if ok else 'FAIL'}", flush=True)
    return 0 if ok else 1

# ---------------------------------------------------------------- timeeval
def timeeval():
    ctx.prec = QUAD_PREC
    p = arb(1) / 3
    ctx.prec = COLLAPSE_PREC
    _, beta = pattern_ac(p, arb("0.1"), arb("0.2"))
    wv = wvec_full(beta, NTOT)
    ctx.prec = QUAD_PREC
    for tag, s1, s2 in [("arb", arb("0.1"), arb("0.2")),
                        ("acb-real", acb("0.1"), acb("0.2")),
                        ("acb-complex", acb("0.1", "0.01"), acb("0.2", "0.01"))]:
        t0 = time.time()
        n = 0
        while time.time() - t0 < 2.0:
            F_collapsed(p, s1, s2, wv)
            n += 1
        print(f"F_collapsed[{tag}] @ {COLLAPSE_PREC}bit: {(time.time()-t0)/n*1e3:.2f} ms/eval")
    t0 = time.time(); n = 0
    while time.time() - t0 < 1.0:
        L_direct(p, arb("0.1"), arb("0.2"), arb("0.3")); n += 1
    print(f"L_direct[arb] @ {QUAD_PREC}bit: {(time.time()-t0)/n*1e3:.3f} ms/eval")
    return 0

# ------------------------------------------------------- quadrature runners
def fatq(z, thr):
    """True if the acb ball z is a wide probe ball (acb_calc error-bound eval),
    not a quadrature node (node balls have radius ~2^-prec)."""
    return bool(arb(z.real.rad()) > thr) or bool(arb(z.imag.rad()) > thr)

def hull_ball(a, b):
    """acb ball covering the real interval [a,b] (verified containment)."""
    h = arb(a).union(arb(b))
    assert h.contains(arb(a)) and h.contains(arb(b))
    return acb(h)

def run_slice2d(args):
    """Collapsed 2-dim (s1,s2) at fixed pi1; s3 closed exactly (opt windowed).

    Fat-ball shortcut (--fat, default 1e-6; 0 disables): when the outer layer
    probes the integrand on a wide s2 ball (acb_calc's per-subinterval and
    Bernstein-ellipse error evaluations), do NOT run the nested inner
    quadrature (it grinds to eval_limit on the fat ball); return the rigorous
    hull enclosure  (b1-a1) * hull(F(s1_hull, s2_ball) * prior)  instead —
    one integrand eval, still a valid enclosure of the inner integral for
    every s2 in the ball, so the outer's subdivision logic is unchanged."""
    ctx.prec = QUAD_PREC
    p = frac(args.p)
    groups = grouped(mutate_counts(args.mutate)) if args.mutate else GROUPS
    J = sum(n for _, n in groups)
    ctx.prec = COLLAPSE_PREC
    one = arb(1)
    beta = one / (2 * (one - p) * p)
    if args.s3win:
        A, B = args.s3win.split(",")
        wv = wvec_win(beta, J, A, B)
    else:
        wv = wvec_full(beta, J)
    ctx.prec = QUAD_PREC
    TOL = lref(args)
    a1, b1 = (frac(x) for x in args.win.split(","))
    fat = arb(args.fat)
    s1h = hull_ball(a1, b1)
    t0 = time.time()
    def f(x2, _):
        if fat > 0 and fatq(x2, fat):
            return (b1 - a1) * F_collapsed(p, s1h, x2, wv, groups) \
                   * prior(s1h) * prior(x2)
        return acb.integral(
            lambda x1, __: F_collapsed(p, x1, x2, wv, groups) * prior(x1),
            a1, b1, rel_tol=arb(args.toli), abs_tol=TOL * arb(args.absi)) * prior(x2)
    v = acb.integral(f, frac(args.o1) if args.o1 else a1,
                     frac(args.o2) if args.o2 else b1,
                     rel_tol=arb(args.tolo), abs_tol=TOL * arb(args.abso))
    row = enclosure_row(v, dict(
        kind="slice2d-collapsed", p=args.p, win=args.win, s3win=args.s3win or "full",
        outer_panel=[args.o1 or args.win.split(",")[0], args.o2 or args.win.split(",")[1]],
        tolo=args.tolo, toli=args.toli, tolref=args.tolref, fat=args.fat,
        mutate=args.mutate or None,
        evals=CALLS[0], secs=round(time.time() - t0, 2)))
    writerow(args.tag, row)
    return 0

def run_direct3d(args):
    """Direct 3-dim (s3 outer / s2 mid / s1 inner) — bench_stage nesting."""
    ctx.prec = QUAD_PREC
    p = frac(args.p)
    TOL = lref(args)      # normalizer from UNMUTATED counts, both routes
    if args.mutate:
        mc = mutate_counts(args.mutate)
        COUNTS.clear(); COUNTS.update(mc)
    a1, b1 = (frac(x) for x in args.win.split(","))
    A, B = (frac(x) for x in args.s3win.split(","))
    fat = arb(args.fat)
    s1h = hull_ball(a1, b1)
    s3h = hull_ball(A, B)
    t0 = time.time()
    def g(x3, _):
        if fat > 0 and fatq(x3, fat):   # fat s3 probe: hull over (s1,s2) box
            return ((b1 - a1) ** 2 * L_direct(p, s1h, s1h, x3)
                    * prior(s1h) ** 2 * prior(x3))
        def f(x2, __):
            if fat > 0 and fatq(x2, fat):   # fat s2 probe: hull over s1
                return ((b1 - a1) * L_direct(p, s1h, x2, x3)
                        * prior(s1h) * prior(x2))
            return acb.integral(
                lambda x1, ___: L_direct(p, x1, x2, x3) * prior(x1),
                a1, b1, rel_tol=arb(args.toli), abs_tol=TOL * arb(args.absi)) * prior(x2)
        return acb.integral(f, a1, b1, rel_tol=arb(args.tolm),
                            abs_tol=TOL * arb(args.absm)) * prior(x3)
    v = acb.integral(g, A, B, rel_tol=arb(args.tolo), abs_tol=TOL * arb(args.abso))
    row = enclosure_row(v, dict(
        kind="direct3d", p=args.p, win=args.win, s3win=args.s3win,
        tolo=args.tolo, tolm=args.tolm, toli=args.toli, tolref=args.tolref,
        fat=args.fat, mutate=args.mutate or None,
        evals=DCALLS[0], secs=round(time.time() - t0, 2)))
    writerow(args.tag, row)
    return 0

def run_evidence(args):
    """pi1 outermost over a panel; per pi1 node: exact-wv + 2-dim (s2,s1)."""
    ctx.prec = QUAD_PREC
    TOL = lref(args)
    a1, b1 = (frac(x) for x in args.win.split(","))
    pa, pb = (frac(x) for x in args.ppanel.split(","))
    fat = arb(args.fat)
    s1h = hull_ball(a1, b1)
    t0 = time.time()
    def outer(pz, _):
        old = ctx.prec
        ctx.prec = COLLAPSE_PREC
        one = pz * 0 + 1
        beta = one / (2 * (one - pz) * pz)
        wv = wvec_full(beta, NTOT)
        ctx.prec = old
        if fat > 0 and fatq(pz, fat):   # fat pi1 probe: hull over (s1,s2) box
            return ((b1 - a1) ** 2 * F_collapsed(pz, s1h, s1h, wv)
                    * prior(s1h) ** 2)
        def f(x2, __):
            if fat > 0 and fatq(x2, fat):
                return ((b1 - a1) * F_collapsed(pz, s1h, x2, wv)
                        * prior(s1h) * prior(x2))
            return acb.integral(
                lambda x1, ___: F_collapsed(pz, x1, x2, wv) * prior(x1),
                a1, b1, rel_tol=arb(args.toli), abs_tol=TOL * arb(args.absi)) * prior(x2)
        return acb.integral(f, a1, b1, rel_tol=arb(args.tolm),
                            abs_tol=TOL * arb(args.absm))
    v = acb.integral(outer, pa, pb, rel_tol=arb(args.tolo),
                     abs_tol=TOL * arb(args.abso))
    row = enclosure_row(v, dict(
        kind="evidence-pi1-panel", ppanel=args.ppanel, win=args.win, s3win="full",
        tolo=args.tolo, tolm=args.tolm, toli=args.toli, tolref=args.tolref,
        fat=args.fat,
        evals=CALLS[0], secs=round(time.time() - t0, 2)))
    writerow(args.tag, row)
    return 0

# ---------------------------------------------------------------------- cli
def main():
    import argparse
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("selftest")
    sub.add_parser("timeeval")
    def common(q):
        q.add_argument("--win", default="0.01,1.5")
        q.add_argument("--mutate", default=None)
        q.add_argument("--tolref", default=None)
        q.add_argument("--fat", default="1e-6")
        q.add_argument("--tag", required=True)
    q = sub.add_parser("slice2d"); common(q)
    q.add_argument("--p", default="1/3"); q.add_argument("--s3win", default=None)
    q.add_argument("--tolo", default="3e-2"); q.add_argument("--toli", default="1e-4")
    q.add_argument("--abso", default="1e-3"); q.add_argument("--absi", default="1e-6")
    q.add_argument("--o1", default=None); q.add_argument("--o2", default=None)
    q = sub.add_parser("direct3d"); common(q)
    q.add_argument("--p", default="1/3"); q.add_argument("--s3win", required=True)
    q.add_argument("--tolo", default="3e-2"); q.add_argument("--tolm", default="1e-3")
    q.add_argument("--toli", default="3e-5")
    q.add_argument("--abso", default="1e-3"); q.add_argument("--absm", default="1e-5")
    q.add_argument("--absi", default="1e-7")
    q = sub.add_parser("evidence"); common(q)
    q.add_argument("--ppanel", required=True)
    q.add_argument("--tolo", default="3e-2"); q.add_argument("--tolm", default="1e-3")
    q.add_argument("--toli", default="3e-5")
    q.add_argument("--abso", default="1e-4"); q.add_argument("--absm", default="1e-6")
    q.add_argument("--absi", default="1e-8")
    args = ap.parse_args()
    if args.cmd == "selftest":
        return selftest()
    if args.cmd == "timeeval":
        return timeeval()
    return dict(slice2d=run_slice2d, direct3d=run_direct3d,
                evidence=run_evidence)[args.cmd](args)

if __name__ == "__main__":
    sys.exit(main())

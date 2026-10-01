# baller engine — two-sided enclosure supports (uncertainty panels).
#!/usr/bin/env python3
"""double_unc.py — WP-G (y2,y3) DOUBLE collapse, UNCORRECTED arm, both shapes.

COLLAPSE_MEMO.md section 6, second partial-collapse route. Per pattern the
peeling sum groups EXACTLY (stationarity kills all odd terms) as

  caterpillar:  pat = A(y1) + B(y1) y2^2 + [C(y1) y2^2 + D(y1) y2^3] y3^2
  balanced:     pat = A(yt) + B(yt) yd^2 + [C(yt) yd^2 + D(yt) yd^3] y3^2

(the memo's 9-monomial support, (e2,e3) in {(0,0),(2,0),(2,2),(3,2)}; closed
forms below, sympy-verified in `selftest` for all 16 patterns, both shapes).
The 318-character product is a bivariate polynomial P(w,u), w = y2 (or yd),
u = y3^2, degrees (954, 318). BOTH remaining exponential integrals close
exactly:

  s3 (full):    w3_j(p) = lam/(lam + 2j beta)  -> EXACT interpolatory rule on
                Jdeg+1 FIXED EXACT-DYADIC Chebyshev-distributed nodes tau_k in
                (0,1):  phi = V^{-T} w3  (V_{jk} = tau_k^j) gives
                sum_k phi_k q(tau_k) = sum_j w3_j [u^j]q  for EVERY polynomial
                q of degree <= Jdeg — a linear-algebra identity, no quadrature
                error term.  phi is computed by dual Bjorck-Pereyra (pvand) at
                1152 bits (measured: |phi| in [7e-11, 7e-3], ~all positive,
                radii ~1e-106; cond(V) ~ 6^n ~ 1e244 absorbed by precision).
                The p-node BALL enters only through w3(a), a = lam p(1-p) - 1;
                its width is carried by an order-T Taylor split in
                delta = a - mid(a) with the EXACT geometric remainder
                  w3_j(a) = 1 - j/(c0+delta),  c0 = 1 + mid(a) + j,
                  j/(c0+delta) = sum_{t<=T} j (-delta)^t / c0^{t+1}
                                 + j (-delta)^{T+1} / (c0^{T+1} (c0+delta)),
                one pvand solve per Taylor order (point RHS -> tight phi_t),
                phi(a) = sum_t delta^t phi_t + phi_rem — valid for every a in
                the ball. (A parametrized certified Gauss-Jacobi rule was
                measured IMPOSSIBLE here: the N-step recurrence amplifies the
                a-ball ~1e15, emptying the interval-Newton window.)
  s2/d (win):   w2_m = lam (e^{-rA}-e^{-rB})/r, r = lam + m beta  -> exact
                coefficient dot (reversal trick: one flint mult per eval)

leaving per p-panel ONE 1-dim certified quadrature (s1 for caterpillar,
t for balanced) plus the eras Taylor-Lagrange remainder in p, which is
UNCHANGED (identical F(p), so gate_price.G_unc_series_fast /
balanced_unc.G_bal_series_fast + sup_integral are reused verbatim).

Per tau node the product is a tree of flint mults of per-group cubics
  (A + (B + C tau_k) w^2 + D tau_k w^3)^n_g
(binary powering); the phi-weighted node polys are accumulated and the s2
dot is coeff M of (sum_k phi_k Q_k) * rev(w2) — ONE flint mult per eval.
No division anywhere in the eval path (weights precomputed per p node) ->
safe on the wide/complex balls acb_calc probes.

REF path (verification): materialize the full bivariate coefficient array
(list of u-blocks, arb/acb_polys in w) by the same binomial expansion +
block-convolution tree; supports point evaluation (anchors), arbitrary
(windowed) w2/w3 contraction, and coefficient-level comparison against the
single y3-collapse (collapse_unc.upoly): blocks[j](y2) == C_j(s1,s2,p).

Subcommands: selftest | rulecheck | timeeval | fiber | uncpanel2 | balpanel2.
"""
import argparse, importlib.util, json, math, os, resource, sys, time
from flint import arb, acb, arb_poly, acb_poly, ctx

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import taylor_p as tp
from taylor_p import frac, binomrow, grouped, mutate_counts, s_var
import gate_price as gp
import balanced_unc as bu

RESULTS = os.environ.get("WPG_RESULTS", os.path.join(os.getcwd(), "results"))
LAM = 10
DC_PREC = int(os.environ.get("DC_PREC", "384"))        # assembly precision
DC_QUAD = int(os.environ.get("DC_QUAD_PREC", "192"))   # 1-dim quadrature layer
COUNTS = tp.COUNTS
GROUPS = tp.GROUPS
NTOT = tp.NTOT
CALLS = [0]     # S_double evals (integrand calls x p-nodes)
QCALLS = [0]    # 1-dim integrand calls


def _load(modname, path):
    spec = importlib.util.spec_from_file_location(modname, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def load_cu():
    """collapse_unc.py: the copy vendored beside this module."""
    p = os.path.join(HERE, "collapse_unc.py")
    if os.path.exists(p):
        return _load("collapse_unc", p)
    raise FileNotFoundError("collapse_unc.py")


def _Pmat(pi0, p, y):
    return [[pi0 + p * y, p - p * y], [pi0 - pi0 * y, p + pi0 * y]]


# ------------------------------------------------------- per-group (A,B,C,D)
def abcd_cat(p, s1, groups=GROUPS):
    """Caterpillar: per grouped pattern (n, A, B, C, D), scalars in the
    arithmetic of (p, s1) (arb or acb).  pat = A + B y2^2 + (C y2^2 + D y2^3) y3^2.
      G_v = P1[v][xA] P1[v][xB];  g = pi0 G0 + p G1;  dG = G0 - G1
      A = pi_xC pi_xD g
      B = (-1)^xC pi_xD pi0 p y1 dG
      C = (-1)^xD pi0 p y1 (pi_xC dG + (-1)^xC y1 g)
      D = (-1)^(xC+xD) pi0 p (p - pi0) y1^2 dG"""
    one = p * 0 + 1
    pi0 = one - p
    beta = one / (2 * pi0 * p)
    y1 = (-beta * s1).exp()
    P1 = _Pmat(pi0, p, y1)
    pip = [pi0, p]
    q = pi0 * p
    out = []
    for (ab, xC, xD), n in groups:
        xA, xB = ab
        G0 = P1[0][xA] * P1[0][xB]
        G1 = P1[1][xA] * P1[1][xB]
        g = pi0 * G0 + p * G1
        dG = G0 - G1
        sC = 1 if xC == 0 else -1
        sD = 1 if xD == 0 else -1
        A = pip[xC] * pip[xD] * g
        B = (q * y1 * dG) * sC * pip[xD]
        C = (q * y1) * (pip[xC] * dG + (y1 * g if sC == 1 else -(y1 * g))) * sD
        D = (q * (p - pi0)) * (y1 * y1 * dG) * (sC * sD)
        out.append((n, A, B, C, D))
    return beta, out


def abcd_bal(p, t, groups):
    """Balanced (one region arm): per grouped pattern (n, A, B, C, D) with
    pat = A + B yd^2 + (C yd^2 + D yd^3) y3^2, functions of yt:
      QC_w = PC[w][xC] PC[w][xD];  h = pi0 QC0 + p QC1;  dQC = QC0 - QC1
      A = pi_xA pi_xB h
      B = (-1)^(xA+xB) pi0 p yt^2 h
      C = pi0 p [(-1)^xA pi_xB + (-1)^xB pi_xA] yt dQC
      D = (-1)^(xA+xB) pi0 p (p - pi0) yt^2 dQC"""
    one = p * 0 + 1
    pi0 = one - p
    beta = one / (2 * pi0 * p)
    yt = (-beta * t).exp()
    PC = _Pmat(pi0, p, yt)
    pip = [pi0, p]
    q = pi0 * p
    out = []
    for (ab, cd), n in groups:
        xA, xB = ab
        xC, xD = cd
        QC0 = PC[0][xC] * PC[0][xD]
        QC1 = PC[1][xC] * PC[1][xD]
        h = pi0 * QC0 + p * QC1
        dQC = QC0 - QC1
        sAB = 1 if (xA + xB) % 2 == 0 else -1
        dA1 = (pip[xB] if xA == 0 else -pip[xB]) + (pip[xA] if xB == 0 else -pip[xA])
        A = pip[xA] * pip[xB] * h
        B = (q * (yt * yt) * h) * sAB
        C = q * dA1 * yt * dQC
        D = (q * (p - pi0)) * (yt * yt) * dQC * sAB
        out.append((n, A, B, C, D))
    return beta, out


# ---------------------- exact-node interpolatory u-contraction (phi weights)
PHI_PREC = int(os.environ.get("DC_PHI_PREC", "1792"))   # phi abs radii ~1e-300:
# the u-contraction sum_k phi_k v_k carries up to ~100 digits of cancellation
# across nodes (measured, bal small-t fiber), so phi radii must sit far below
# 1e-(106+100); pvand is O(n^2) — the bump costs ~0.1 s per p node.
PHI_T = int(os.environ.get("DC_PHI_T", "12"))


def cheb_nodes01(N):
    """N distinct EXACT DYADIC Chebyshev-distributed nodes in (0,1),
    increasing (26-bit dyadic rounding of the first-kind points)."""
    nodes = []
    for k in range(N):
        v = 0.5 * (1.0 - math.cos(math.pi * (2 * k + 1) / (2 * N)))
        nodes.append(arb(round(v * (1 << 26))) / (1 << 26))
    assert len({float(x) for x in nodes}) == N
    return nodes


def pvand(al, b):
    """Solve V x = b with V_{jk} = al_k^j (dual Vandermonde / quadrature-
    weight system) by Bjorck-Pereyra, O(n^2), all ops ball arithmetic ->
    the returned balls rigorously contain the exact solution.
    (Verified against dense solve at n=6; at n=319/prec 1152 the measured
    weights are ~all positive in [7e-11, 7e-3] with radii ~1e-106.)"""
    n = len(b)
    x = list(b)
    for k in range(n - 1):
        for j in range(n - 1, k, -1):
            x[j] = x[j] - al[k] * x[j - 1]
    for k in range(n - 2, -1, -1):
        for j in range(k + 1, n):
            x[j] = x[j] / (al[j] - al[j - k - 1])
        for j in range(k, n - 1):
            x[j] = x[j] - x[j + 1]
    return x


def phi_weights(Jdeg, beta_ball, lam=LAM, T=None, prec=None):
    """Interpolatory u-contraction weights phi_k (arb balls) on the exact
    dyadic nodes, valid for EVERY p in the p-ball: with a+1 = lam/(2 beta),
    w3_j(a) = lam/(lam+2j beta) = 1 - j/(1+a+j).  Split in
    delta = (a - a0), a0 = exact midpoint:
      j/(c0+delta) = sum_{t=0..T} j (-delta)^t / c0^{t+1}
                     + j (-delta)^{T+1} / (c0^{T+1} (c0+delta))   [EXACT]
    One pvand solve per order (point RHS) + one for the remainder ball;
    phi(a) = sum_t delta^t phi_t + phi_rem contains V^{-T} w3(a) for every
    a in the ball. Returns (taus, phis)."""
    old = ctx.prec
    ctx.prec = prec or PHI_PREC
    T = T or PHI_T
    try:
        n = Jdeg + 1
        taus = cheb_nodes01(n)
        aball = lam / (2 * beta_ball) - 1          # a = lam p(1-p) - 1
        a0 = arb(aball.mid())
        delta = aball - a0
        # RHS vectors: t=0: w3_j(a0) = 1 - j/c0; t>=1: -(-1)^t j/c0^{t+1}
        rhs = []
        c0s = [1 + a0 + j for j in range(n)]
        rhs.append([1 - j / c0s[j] for j in range(n)])
        for t in range(1, T + 1):
            sgn = -1 if t % 2 == 0 else 1          # -(-1)^t
            rhs.append([sgn * j / c0s[j] ** (t + 1) for j in range(n)])
        dpow = 1 + delta * 0                       # (-delta)^(T+1) by repeated
        for _ in range(T + 1):                     # mult (arb ** on a ball
            dpow = dpow * (-delta)                 # containing 0 -> NaN)
        rem = [-(j * dpow) / (c0s[j] ** (T + 1) * (c0s[j] + delta))
               for j in range(n)]
        rhs.append(rem)
        sols = [pvand(taus, r) for r in rhs]
        phis = []
        for k in range(n):
            acc = sols[T][k]
            for t in range(T - 1, 0, -1):          # Horner in delta, t=T..1
                acc = acc * delta + sols[t][k]
            phis.append(sols[0][k] + acc * delta + sols[T + 1][k])
        # guard: vacuous weights (p-ball too wide for order T) must not
        # silently burn a panel run. GL-node balls (~1e-33) give ~1e-100s.
        mr = max(float(arb(x.rad())) for x in phis)
        if mr > 1e-20:
            raise ValueError(
                f"phi_weights: max phi radius {mr:.2e} — p-ball too wide for "
                f"Taylor order T={T} (supported input: GL-node-width balls)")
        return taus, phis
    finally:
        ctx.prec = old


# ------------------------------------------------------------ weight vectors
def w2vec(beta, M, B2, lam=LAM):
    """s2/d closure weights, m = 0..M: int_0^B2 lam e^{-lam s} w^m ds
    (B2 None/'inf' -> full).  w = exp(-beta s)."""
    Bv = None if B2 in (None, "inf") else frac(B2)
    out = []
    for m in range(M + 1):
        r = lam + m * beta
        if Bv is None:
            out.append(lam / r)
        else:
            out.append(lam * (1 - (-r * Bv).exp()) / r)
    return out


def w3vec(beta, J, lam=LAM, win=None):
    """s3 closure weights (u = y3^2): lam/(lam+2j beta), optional window."""
    if win is None:
        return [lam / (lam + (2 * j) * beta) for j in range(J + 1)]
    A, B = frac(win[0]), (None if win[1] in (None, "inf") else frac(win[1]))
    out = []
    for j in range(J + 1):
        r = lam + (2 * j) * beta
        v = (-r * A).exp()
        if B is not None:
            v = v - (-r * B).exp()
        out.append(lam * v / r)
    return out


# -------------------------------------------------- production eval (phi+dot)
class PNode:
    """Everything precomputed per p value: exact tau nodes, phi(a-ball)
    contraction weights, w2 reversal polys, sizes."""

    def __init__(self, p, groups, shape, B2, prec=None):
        old = ctx.prec
        self.prec = prec or DC_PREC
        ctx.prec = self.prec
        try:
            self.p = p
            self.shape = shape          # 'cat' | 'bal'
            self.groups = groups        # cat: groups; bal: (g1, g2)
            gg = groups if shape == "cat" else groups[0]
            self.ntot = sum(n for _, n in gg)
            self.M = 3 * self.ntot
            self.Jdeg = self.ntot
            one = p * 0 + 1
            self.beta = one / (2 * (one - p) * p)
            self.taus, self.phis = phi_weights(self.Jdeg, self.beta)
            wv = w2vec(self.beta, self.M, B2)
            self.revW2 = arb_poly([wv[self.M - i] for i in range(self.M + 1)])
            self.revW2c = acb_poly([acb(wv[self.M - i]) for i in range(self.M + 1)])
            self.tausc = [acb(t) for t in self.taus]
            self.phisc = [acb(w) for w in self.phis]
        finally:
            ctx.prec = old


def _ppow(q, n):
    """q**n by binary powering (flint poly mults)."""
    r = None
    while n:
        if n & 1:
            r = q if r is None else r * q
        n >>= 1
        if n:
            q = q * q
    return r


def _tree(facs):
    while len(facs) > 1:
        facs.sort(key=lambda x: x.length())
        facs.append(facs.pop(0) * facs.pop(0))
    return facs[0]


def S_double(pn, sx):
    """The double-closed object at (p, s1) [cat] or (p, t) [bal]:
    S = sum_{m,j} C_{mj} w2_m w3_j, exact in both closures.
    sx may be arb or acb (acb path used by acb_calc)."""
    CALLS[0] += 1
    old = ctx.prec
    ctx.prec = pn.prec
    try:
        cx = isinstance(sx, acb) or isinstance(pn.p, acb)
        Poly = acb_poly if cx else arb_poly
        z = (acb(0) if cx else arb(0))
        p = acb(pn.p) if cx and not isinstance(pn.p, acb) else pn.p
        sxx = acb(sx) if cx and not isinstance(sx, acb) else sx
        if pn.shape == "cat":
            armsets = [abcd_cat(p, sxx, pn.groups)[1]]
        else:
            armsets = [abcd_bal(p, sxx, pn.groups[0])[1],
                       abcd_bal(p, sxx, pn.groups[1])[1]]
        taus = pn.tausc if cx else pn.taus
        phis = pn.phisc if cx else pn.phis
        revW = pn.revW2c if cx else pn.revW2
        M = pn.M
        R = None                        # accumulated sum_k phi_k Q_k
        for tk, wk in zip(taus, phis):
            Q = None
            for abcd in armsets:
                facs = []
                for (n, A, B, C, D) in abcd:
                    cub = Poly([A, z, B + C * tk, D * tk])
                    facs.append(_ppow(cub, n))
                Qa = _tree(facs)
                Q = Qa if Q is None else Q + Qa
            Qw = Poly([wk]) * Q
            R = Qw if R is None else R + Qw
        return (R * revW)[M]            # single w2 dot
    finally:
        ctx.prec = old


# ------------------------------------------------- REF: materialized bivariate
def bivar_blocks(p, sx, groups, shape, prec=None):
    """Full bivariate product as u-blocks: blocks[j] = poly in w (the u^j
    coefficient), assembled by binomial expansion + block-convolution tree."""
    old = ctx.prec
    ctx.prec = prec or DC_PREC
    try:
        cx = isinstance(sx, acb) or isinstance(p, acb)
        Poly = acb_poly if cx else arb_poly
        if shape == "cat":
            armsets = [abcd_cat(p, sx, groups)[1]]
        else:
            armsets = [abcd_bal(p, sx, groups[0])[1],
                       abcd_bal(p, sx, groups[1])[1]]
        btot = None
        for abcd in armsets:
            gpolys = []
            for (n, A, B, C, D) in abcd:
                al = Poly([A, A * 0, B])
                be = Poly([A * 0, A * 0, C, D])
                apow = [Poly([A * 0 + 1])]
                bpow = [Poly([A * 0 + 1])]
                for _ in range(n):
                    apow.append(apow[-1] * al)
                    bpow.append(bpow[-1] * be)
                bn = binomrow(n)
                gpolys.append([Poly([A * 0 + bn[j]]) * (apow[n - j] * bpow[j])
                               for j in range(n + 1)])
            while len(gpolys) > 1:              # merge smallest u-degrees first
                gpolys.sort(key=len)
                Ab, Bb = gpolys.pop(0), gpolys.pop(0)
                out = [None] * (len(Ab) + len(Bb) - 1)
                for i, Ai in enumerate(Ab):
                    for jj, Bj in enumerate(Bb):
                        t = Ai * Bj
                        out[i + jj] = t if out[i + jj] is None else out[i + jj] + t
                gpolys.append(out)
            blocks = gpolys[0]
            if btot is None:
                btot = blocks
            else:
                btot = [x + y for x, y in zip(btot, blocks)]
        return btot
    finally:
        ctx.prec = old


def ref_point(blocks, w0, u0):
    """P(w0, u0) by Horner over u-blocks."""
    acc = blocks[-1](w0)
    for j in range(len(blocks) - 2, -1, -1):
        acc = acc * u0 + blocks[j](w0)
    return acc


def ref_contract(blocks, w2v, w3v):
    """sum_j w3_j * dot(blocks[j], w2) via reversal mults."""
    M = len(w2v) - 1
    rev = arb_poly([w2v[M - i] for i in range(M + 1)])
    tot = None
    for j, bl in enumerate(blocks):
        v = (bl * rev)[M] * w3v[j]
        tot = v if tot is None else tot + v
    return tot


# ------------------------------------------------------------------ helpers
def prior(s, lam=LAM):
    return lam * (-(s * lam)).exp()


def rss_mb():
    return round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 1)


def relw(b):
    b = b.real if isinstance(b, acb) else b
    return float("inf") if b.mid() == 0 else float(arb(b.rad()) / abs(arb(b.mid())))


def bank(tag, row):
    os.makedirs(RESULTS, exist_ok=True)
    json.dump(row, open(os.path.join(RESULTS, f"{tag}.json"), "w"), indent=1)
    print(json.dumps(row), flush=True)


def progress_hook(tag):
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
                    json.dumps(dict(evals=CALLS[0], qcalls=QCALLS[0],
                                    secs=round(time.time() - t0, 1),
                                    maxrss_mb=rss_mb())))
            except Exception:
                pass
    threading.Thread(target=loop, daemon=True).start()


def balgroups(counts=None):
    c = COUNTS if counts is None else counts
    return (bu.grouped_bal(c, False), bu.grouped_bal(c, True))


# ------------------------------------------------------------- panel runners
def run_uncpanel2(a):
    """Caterpillar eras p-panel through the DOUBLE collapse: order-K GL nodes,
    per node the 1-dim certified s1 quadrature of prior(s1) * S_double
    (s2 closed EXACTLY over [0,shi], s3 closed EXACTLY, full); remainder =
    UNCHANGED gate_price path (same F(p), same Mfac * sup|f_K|)."""
    progress_hook(a.tag)
    t0 = time.time()
    ctx.prec = DC_PREC          # panel nodes must carry GL-guard-width balls
    groups = grouped(mutate_counts(a.mutate)) if a.mutate else GROUPS
    pjs, Wjs, Mfac = gp.panel_setup(a.pm, a.pr, a.K)
    tR0 = time.time()
    pns = [PNode(pj, groups, "cat", a.shi) for pj in pjs]
    secs_rule = time.time() - tR0
    FLOOR = arb(a.floor)
    ctx.prec = DC_QUAD

    def f(x1, _):
        QCALLS[0] += 1
        tot = acb(0)
        for pn, Wj in zip(pns, Wjs):
            tot += acb(Wj) * S_double(pn, x1)
        return tot * prior(x1)

    J = acb.integral(f, arb(0), arb(a.shi), rel_tol=arb(a.tolo), abs_tol=FLOOR)
    tmain = time.time() - t0
    # remainder — identical to gate_price.run_uncpanel (same F(p))
    ctx.prec = gp.COLLAPSE_PREC
    pball = arb(arb(a.pm).mid(), a.pr)

    def fmag(cells):
        return abs(gp.G_unc_series_fast(s_var(pball, a.K), cells[0], cells[1],
                                        groups))
    U, rev = gp.sup_integral(fmag, [(0.0, a.shi)] * 2,
                             FLOOR / max(Mfac, arb("1e-300")), a.maxrev)
    R = gp.upper(Mfac * U)
    row = dict(kind="uncpanel2-double", pm=a.pm, pr=a.pr, K=a.K, shi=a.shi,
               J_mid=arb(J.real.mid()).str(20),
               J_rad=arb(J.real.rad()).str(5, radius=False),
               R=R.str(5, radius=False), Mfac=float(Mfac.mid()),
               relw_J=relw(J), floor=a.floor, tolo=a.tolo,
               dc_prec=DC_PREC, quad_prec=DC_QUAD,
               mutate=a.mutate or None,
               evals=CALLS[0], qcalls=QCALLS[0], rem_evals=rev,
               secs_rule=round(secs_rule, 1), secs_main=round(tmain, 1),
               secs=round(time.time() - t0, 1), maxrss_mb=rss_mb())
    bank(a.tag, row)
    return 0


def run_balpanel2(a):
    """Balanced eras p-panel, double collapse: d closed EXACTLY over [0,dhi],
    s3 closed EXACTLY (full), both region arms; 1-dim certified t quadrature
    with prior lam e^{-2 lam t}; remainder = balanced_unc path unchanged."""
    progress_hook(a.tag)
    t0 = time.time()
    ctx.prec = DC_PREC          # panel nodes must carry GL-guard-width balls
    if a.mutate:
        mc = mutate_counts(a.mutate)
        g12 = (bu.grouped_bal(mc, False), bu.grouped_bal(mc, True))
    else:
        g12 = balgroups()
    pjs, Wjs, Mfac = gp.panel_setup(a.pm, a.pr, a.K)
    tR0 = time.time()
    pns = [PNode(pj, g12, "bal", a.dhi) for pj in pjs]
    secs_rule = time.time() - tR0
    FLOOR = arb(a.floor)
    ctx.prec = DC_QUAD

    def f(x1, _):
        QCALLS[0] += 1
        tot = acb(0)
        for pn, Wj in zip(pns, Wjs):
            tot += acb(Wj) * S_double(pn, x1)
        return tot * LAM * (-(x1 * (2 * LAM))).exp()

    J = acb.integral(f, arb(0), arb(a.thi), rel_tol=arb(a.tolo), abs_tol=FLOOR)
    tmain = time.time() - t0
    ctx.prec = bu.COLLAPSE_PREC
    pball = arb(arb(a.pm).mid(), a.pr)

    def fmag(cells):
        return abs(bu.G_bal_series_fast(s_var(pball, a.K), cells[0], cells[1],
                                        g12[0], g12[1]))
    U, rev = bu._sup_integral_bal(fmag, [(0.0, a.thi), (0.0, a.dhi)],
                                  FLOOR / max(Mfac, arb("1e-300")), a.maxrev)
    R = gp.upper(Mfac * U)
    row = dict(kind="balpanel2-double", pm=a.pm, pr=a.pr, K=a.K, thi=a.thi,
               dhi=a.dhi, J_mid=arb(J.real.mid()).str(20),
               J_rad=arb(J.real.rad()).str(5, radius=False),
               R=R.str(5, radius=False), Mfac=float(Mfac.mid()),
               relw_J=relw(J), floor=a.floor, tolo=a.tolo,
               dc_prec=DC_PREC, quad_prec=DC_QUAD,
               mutate=a.mutate or None,
               evals=CALLS[0], qcalls=QCALLS[0], rem_evals=rev,
               secs_rule=round(secs_rule, 1), secs_main=round(tmain, 1),
               secs=round(time.time() - t0, 1), maxrss_mb=rss_mb())
    bank(a.tag, row)
    return 0


# ----------------------------------------------------------------- rulecheck
def rulecheck():
    """phi(a-ball) weights vs exact moments w3_j(p) for POINT p values inside
    the p-ball: the assembled phi must satisfy sum_k phi_k tau_k^j inside a
    ball containing w3_j(p_point) for every j <= Jdeg — the containment form
    of the exactness identity. Panel-width balls (pr 0.002/0.003, GL-node
    guard 1e-30) and a wide control ball are tested."""
    ok = True
    # radii = the ACTUAL PNode regime: GL panel nodes have rad ~ pr*1e-30
    # (gl_rule guard balls); 1e-20 as stress margin. Wide panel balls are NOT
    # a supported input (phi_weights guards against vacuous output).
    for pm, pr in ((0.072955, 3e-33), (0.346, 2e-33), (0.94, 3e-33),
                   (0.346, 1e-20)):
        ctx.prec = DC_PREC
        pball = arb(arb(pm).mid(), pr)
        one = arb(1)
        beta = one / (2 * (one - pball) * pball)
        t0 = time.time()
        taus, phis = phi_weights(NTOT, beta)
        secs = time.time() - t0
        maxrad = max(float(arb(x.rad())) for x in phis)
        worst = arb(0)
        nfail = 0
        ctx.prec = PHI_PREC
        for off in (-0.999, -0.5, 0.0, 0.5, 0.999):   # points inside the ball
            ppt = arb(arb(pm).mid()) + arb(str(off)) * pr
            assert pball.contains(ppt)
            bpt = 1 / (2 * (1 - ppt) * ppt)
            for j in (0, 1, 7, 50, 159, NTOT):
                s = arb(0)
                for tk, wk in zip(taus, phis):
                    s += wk * tk ** j
                ex = LAM / (LAM + 2 * j * bpt)
                if not s.contains(ex) and not s.overlaps(ex):
                    print(f"FAIL p={pm}+/-{pr} off={off} j={j}: "
                          f"{s.str(10)} vs {ex.str(10)}")
                    nfail += 1
                d = abs(s - ex)
                worst = max(worst, arb(d.mid()) + arb(d.rad()),
                            key=lambda x: float(x))
        ok = ok and nfail == 0
        print(f"{'OK  ' if nfail == 0 else 'FAIL'} p={pm}+/-{pr}: phi built "
              f"{secs:.1f}s (n={len(taus)}, max phi rad {maxrad:.2e}), "
              f"5x6 point-moments contained, worst |diff| bound {worst.str(3)}",
              flush=True)
    print(f"RULECHECK {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


# ------------------------------------------------------------------ selftest
def selftest():
    ok = True
    t00 = time.time()
    cu = load_cu()
    print("== double_unc selftest: (y2,y3) double collapse ==", flush=True)

    # (0) sympy: exact (A,B,C,D) closed forms, all 16 patterns, both shapes
    import sympy as sp
    P = sp.symbols("p")
    pi0s = 1 - P
    Ps = lambda y: [[pi0s + P * y, P - P * y], [pi0s - pi0s * y, P + pi0s * y]]
    pis = [pi0s, P]
    # caterpillar in (y1, y2, y3)
    Y1, Y2, Y3 = sp.symbols("y1 y2 y3")
    P1s, P2s, P3s = Ps(Y1), Ps(Y2), Ps(Y3)
    P12s, P123s = Ps(Y1 * Y2), Ps(Y1 * Y2 * Y3)
    bad = 0
    for xi in range(16):
        x = tuple(int(b) for b in format(xi, "04b"))
        pat = sp.expand(sum(
            pis[r] * P3s[r][w] * P2s[w][v] * P1s[v][x[0]] * P1s[v][x[1]]
            * P12s[w][x[2]] * P123s[r][x[3]]
            for r in (0, 1) for w in (0, 1) for v in (0, 1)))
        G0 = P1s[0][x[0]] * P1s[0][x[1]]
        G1 = P1s[1][x[0]] * P1s[1][x[1]]
        g = pi0s * G0 + P * G1
        dG = G0 - G1
        sC, sD = (-1) ** x[2], (-1) ** x[3]
        A = pis[x[2]] * pis[x[3]] * g
        B = sC * pis[x[3]] * pi0s * P * Y1 * dG
        C = sD * pi0s * P * Y1 * (pis[x[2]] * dG + sC * Y1 * g)
        D = sC * sD * pi0s * P * (P - pi0s) * Y1 ** 2 * dG
        if sp.expand(pat - (A + B * Y2 ** 2 + (C * Y2 ** 2 + D * Y2 ** 3) * Y3 ** 2)) != 0:
            print(f"FAIL cat closed form: pattern {x}")
            bad += 1
    print(f"{'OK  ' if bad == 0 else 'FAIL'} cat: pat == A + B y2^2 + "
          f"(C y2^2 + D y2^3) y3^2 exact, 16/16 patterns", flush=True)
    ok = ok and bad == 0
    # balanced in (yt, yd, y3), region-1 tree
    YT, YD = sp.symbols("yt yd")
    PAs, PCs = Ps(YT * YD), Ps(YT)
    PrA, PrC = Ps(Y3), Ps(Y3 * YD)
    bad = 0
    for xi in range(16):
        x = tuple(int(b) for b in format(xi, "04b"))
        pat = sp.expand(sum(
            pis[r] * PrA[r][v] * PrC[r][w]
            * PAs[v][x[0]] * PAs[v][x[1]] * PCs[w][x[2]] * PCs[w][x[3]]
            for r in (0, 1) for v in (0, 1) for w in (0, 1)))
        QC0 = PCs[0][x[2]] * PCs[0][x[3]]
        QC1 = PCs[1][x[2]] * PCs[1][x[3]]
        h = pi0s * QC0 + P * QC1
        dQC = QC0 - QC1
        sAB = (-1) ** (x[0] + x[1])
        dA1 = (-1) ** x[0] * pis[x[1]] + (-1) ** x[1] * pis[x[0]]
        A = pis[x[0]] * pis[x[1]] * h
        B = sAB * pi0s * P * YT ** 2 * h
        C = pi0s * P * dA1 * YT * dQC
        D = sAB * pi0s * P * (P - pi0s) * YT ** 2 * dQC
        if sp.expand(pat - (A + B * YD ** 2 + (C * YD ** 2 + D * YD ** 3) * Y3 ** 2)) != 0:
            print(f"FAIL bal closed form: pattern {x}")
            bad += 1
    print(f"{'OK  ' if bad == 0 else 'FAIL'} bal: pat == A + B yd^2 + "
          f"(C yd^2 + D yd^3) y3^2 exact, 16/16 patterns", flush=True)
    ok = ok and bad == 0

    # numeric (A,B,C,D) vs pattern_ac / pattern_ac_bal at a generic point
    ctx.prec = DC_PREC
    p, s1, s2 = arb("0.22"), arb("0.07"), arb("0.9")
    one = arb(1)
    beta = one / (2 * (one - p) * p)
    y2 = (-beta * s2).exp()
    ac, _ = cu.pattern_ac(p, s1, s2)
    _, abcd = abcd_cat(p, s1)
    worst = arb(0)
    for (n1, av, cv), (n2, A, B, C, D) in zip(ac, abcd):
        assert n1 == n2
        d1 = abs(av - (A + B * y2 ** 2))
        d2 = abs(cv - (C * y2 ** 2 + D * y2 ** 3))
        worst = max(worst, arb(d1.mid()) + arb(d1.rad()),
                    arb(d2.mid()) + arb(d2.rad()), key=float)
    hit = bool(worst < arb("1e-100"))
    print(f"{'OK  ' if hit else 'FAIL'} cat ABCD vs pattern_ac numeric: "
          f"worst |diff| {worst.str(3)}")
    ok = ok and hit
    t, d = arb("0.2"), arb("0.15")
    yd = (-beta * d).exp()
    for gi, g in enumerate((bu.G1, bu.G2)):
        acb_, _ = bu.pattern_ac_bal(p, t, d, g)
        _, abcdb = abcd_bal(p, t, g)
        worst = arb(0)
        for (n1, av, cv), (n2, A, B, C, D) in zip(acb_, abcdb):
            assert n1 == n2
            d1 = abs(av - (A + B * yd ** 2))
            d2 = abs(cv - (C * yd ** 2 + D * yd ** 3))
            worst = max(worst, arb(d1.mid()) + arb(d1.rad()),
                        arb(d2.mid()) + arb(d2.rad()), key=float)
        hit = bool(worst < arb("1e-100"))
        print(f"{'OK  ' if hit else 'FAIL'} bal arm{gi+1} ABCD vs "
              f"pattern_ac_bal: worst |diff| {worst.str(3)}")
        ok = ok and hit

    # (a) T1/T2 anchors through the ASSEMBLED bivariate (REF blocks), 30+ digits
    print("\n-- (a) anchors through the double-collapse bivariate, prec 384 --")
    for name, ((sa, sb, sc, pp), unc50) in cu.ANCHORS.items():
        s1, s2, s3, p = frac(sa), frac(sb), frac(sc), frac(pp)
        one = p * 0 + 1
        beta = one / (2 * (one - p) * p)
        blocks = bivar_blocks(p, s1, GROUPS, "cat")
        assert len(blocks) == NTOT + 1, len(blocks)
        wdeg = max(b.length() - 1 for b in blocks)
        w0 = (-beta * s2).exp()
        u0 = ((-beta * s3).exp()) ** 2
        lv = ref_point(blocks, w0, u0).log()
        dl = abs(lv - arb(unc50))
        hit = bool(dl < arb("1e-30")) and bool(arb(lv.rad()) < arb("1e-40"))
        print(f"{'OK  ' if hit else 'FAIL'} {name}: |dlog anchor| < {dl.str(3)}, "
              f"rad {arb(lv.rad()).str(3, radius=False)} "
              f"(u-deg {len(blocks)-1}, w-deg {wdeg})")
        ok = ok and hit

    # (b) coefficient-level agreement vs the y3-only collapse at 3 fibers
    print("\n-- (b) vs single collapse: 319 coefficients + closed values, 3 fibers --")
    fibers = [("T1", "1/3", "1/10", "2/10"), ("T2", "2/5", "1/2", "1/4"),
              ("LIVE", "0.346", "0.3", "0.7")]
    for name, pp, sa, sb in fibers:
        p, s1, s2 = frac(pp), frac(sa), frac(sb)
        one = p * 0 + 1
        beta = one / (2 * (one - p) * p)
        y2v = (-beta * s2).exp()
        poly1, _ = cu.upoly(p, s1, s2)          # y3-only collapse, 384-bit ctx
        blocks = bivar_blocks(p, s1, GROUPS, "cat")
        worst, novr = arb(0), 0
        for j in range(NTOT + 1):
            cj = blocks[j](y2v)
            c1 = poly1[j]
            if not cj.overlaps(c1):
                novr += 1
            dd = abs(cj - c1)
            worst = max(worst, arb(dd.mid()) + arb(dd.rad()), key=float)
        hit = novr == 0
        print(f"{'OK  ' if hit else 'FAIL'} {name}: all 319 u-coeffs overlap "
              f"single-collapse ({novr} fail), worst |diff| bound {worst.str(3)}")
        ok = ok and hit
        # closed fiber values: full + windowed s3
        for win in (None, ("0", "2"), ("0.28", "0.32")):
            wv3 = w3vec(beta, NTOT, win=win)
            fs = poly1[0] * wv3[0]
            for j in range(1, NTOT + 1):
                fs += poly1[j] * wv3[j]
            w2id = [y2v ** m for m in range(3 * NTOT + 1)]   # identity dot = eval
            fd = ref_contract(blocks, w2id, wv3)
            dl = abs(fd.log() - fs.log())
            hit = fd.overlaps(fs) and bool(dl < arb("1e-40"))
            print(f"{'OK  ' if hit else 'FAIL'}   s3win={win or 'full'}: "
                  f"overlap, |dlog| < {dl.str(3)}")
            ok = ok and hit

    # (b2) s2-closure: GJ production value vs certified s2 quadrature, 3 fibers
    print("\n-- (b2) exact s2 closure vs certified 1-dim s2 quadrature --")
    B2 = "3.5"
    for name, pp, sa in (("T1", "1/3", "1/10"), ("LIVE", "0.346", "0.3"),
                         ("WING", "0.072955", "0.5")):
        p, s1 = frac(pp), frac(sa)
        pn = PNode(p, GROUPS, "cat", B2)
        Sgj = S_double(pn, s1)
        one = p * 0 + 1
        beta = one / (2 * (one - p) * p)
        wv3 = w3vec(beta, NTOT)
        ctx.prec = DC_PREC
        direct = acb.integral(
            lambda x, _: prior(x) * cu.F_collapsed(p, s1, x, wv3),
            arb(0), arb(frac(B2)), rel_tol=arb("1e-30"),
            abs_tol=abs(arb(Sgj.mid())) * arb("1e-31"))
        dl = abs(Sgj.log() - direct.real.log())
        hit = Sgj.overlaps(direct.real) and bool(dl < arb("1e-25"))
        print(f"{'OK  ' if hit else 'FAIL'} {name}: phi vs quadrature overlap, "
              f"|dlog| < {dl.str(3)} (GJ rad {arb(Sgj.rad()).str(3, radius=False)})")
        ok = ok and hit
        # and GJ vs REF window-contraction (both exact routes)
        blocks = bivar_blocks(p, s1, GROUPS, "cat")
        w2v = w2vec(beta, 3 * NTOT, B2)
        Sref = ref_contract(blocks, w2v, wv3)
        dlr = abs(Sgj.log() - Sref.log())
        hit = Sgj.overlaps(Sref) and bool(dlr < arb("1e-25"))
        print(f"{'OK  ' if hit else 'FAIL'}   phi vs REF contraction overlap, "
              f"|dlog| < {dlr.str(3)}")
        ok = ok and hit

    # (b3) balanced: closure + agreement vs F_bal at 2 fibers, both arms
    print("\n-- (b3) balanced double collapse vs F_bal + d-quadrature --")
    # bal fibers at prec 576: at small t the (yd,y3) assembly cancellation is
    # MEASURED ~1.2 bits/char (F2 relwidth 7.5e-6 at 384 bits — still >>
    # tighter than the 3e-2 panel register, but too wide for the 1e-25
    # exact-agreement gates; 576 restores ~60 digits of margin).
    g12 = balgroups()
    BPREC = 576
    for name, pp, st in (("F1", "1/3", "0.2"), ("F2", "0.22", "0.05")):
        p, t = frac(pp), frac(st)
        pn = PNode(p, g12, "bal", "3.5", prec=BPREC)
        Sgj = S_double(pn, t)
        one = p * 0 + 1
        beta = one / (2 * (one - p) * p)
        wv3 = bu.wvec_full(beta, NTOT)
        ctx.prec = BPREC
        direct = acb.integral(
            lambda x, _: prior(x) * bu.F_bal(p, t, x, wv3),
            arb(0), arb("3.5"), rel_tol=arb("1e-30"),
            abs_tol=abs(arb(Sgj.mid())) * arb("1e-31"), eval_limit=2000000)
        dl = abs(Sgj.log() - direct.real.log())
        refw = relw(direct.real)
        hit = Sgj.overlaps(direct.real) and \
            bool(dl < max(arb("1e-25"), arb(refw) * 30))
        print(f"{'OK  ' if hit else 'FAIL'} {name}: bal phi vs d-quadrature, "
              f"|dlog| < {dl.str(3)} (reference relwidth {refw:.2e})")
        ok = ok and hit
        # exact-vs-exact: phi route vs REF materialized contraction
        blocks = bivar_blocks(p, t, g12, "bal", prec=BPREC)
        w2v = w2vec(beta, 3 * NTOT, "3.5")
        Sref = ref_contract(blocks, w2v, wv3)
        dlr = abs(Sgj.log() - Sref.log())
        hit = Sgj.overlaps(Sref) and bool(dlr < arb("1e-25"))
        print(f"{'OK  ' if hit else 'FAIL'}   bal phi vs REF contraction, "
              f"|dlog| < {dlr.str(3)}")
        ok = ok and hit

    # (d) planted error 0001:+1 through the new path (fiber level)
    print("\n-- (d) planted error 0001:+1 through the double collapse --")
    mut = mutate_counts("0001:+1")
    gmut = grouped(mut)
    p, s1 = frac("1/3"), frac("1/10")
    pn0 = PNode(p, GROUPS, "cat", B2)
    pnm = PNode(p, gmut, "cat", B2)
    S0 = S_double(pn0, s1)
    Sm = S_double(pnm, s1)
    sh = Sm.log() - S0.log()
    # same shift through the single collapse + s2 quadrature
    one = p * 0 + 1
    beta = one / (2 * (one - p) * p)
    ctx.prec = DC_PREC
    ATOL = abs(arb(S0.mid())) * arb("1e-31")
    d0 = acb.integral(lambda x, _: prior(x) * cu.F_collapsed(
        p, s1, x, w3vec(beta, NTOT)), arb(0), arb(frac(B2)),
        rel_tol=arb("1e-30"), abs_tol=ATOL)
    dm = acb.integral(lambda x, _: prior(x) * cu.F_collapsed(
        p, s1, x, w3vec(beta, NTOT + 1), gmut), arb(0), arb(frac(B2)),
        rel_tol=arb("1e-30"), abs_tol=ATOL)
    shq = dm.real.log() - d0.real.log()
    fired = bool(not sh.contains(arb(0))) and bool(abs(arb(sh.mid())) > arb("1e-2"))
    agree = bool(abs(sh - shq) < arb("1e-20"))
    print(f"{'OK  ' if fired and agree else 'FAIL'} plant shift {sh.str(10)} "
          f"excludes 0; matches single-collapse route |diff| < "
          f"{abs(sh - shq).str(3)}")
    ok = ok and fired and agree

    print(f"\nSELFTEST {'PASS' if ok else 'FAIL'} "
          f"({time.time()-t00:.1f}s, maxrss {rss_mb()} MB)", flush=True)
    return 0 if ok else 1


# ------------------------------------------------------------------ timeeval
def timeeval():
    ctx.prec = DC_PREC
    p = arb("0.346")
    t0 = time.time()
    pn = PNode(p, GROUPS, "cat", "3.5")
    print(f"PNode(cat) setup (phi n={len(pn.taus)} + w2): "
          f"{time.time()-t0:.1f}s")
    for tag, sx in (("arb", arb("0.3")), ("acb", acb("0.3")),
                    ("acb-wide", acb(arb("0.3", "0.05"), arb("0", "0.05")))):
        t0 = time.time()
        n = 0
        while time.time() - t0 < 5.0 and n < 50:
            S_double(pn, sx)
            n += 1
        print(f"S_double(cat)[{tag}] @ {DC_PREC}bit: "
              f"{(time.time()-t0)/n*1e3:.0f} ms/eval", flush=True)
    g12 = balgroups()
    t0 = time.time()
    pnb = PNode(p, g12, "bal", "3.5")
    print(f"PNode(bal) setup: {time.time()-t0:.1f}s")
    t0 = time.time()
    n = 0
    while time.time() - t0 < 5.0 and n < 50:
        S_double(pnb, acb("0.2"))
        n += 1
    print(f"S_double(bal)[acb] @ {DC_PREC}bit: {(time.time()-t0)/n*1e3:.0f} ms/eval")
    t0 = time.time()
    bivar_blocks(p, arb("0.3"), GROUPS, "cat")
    print(f"bivar_blocks(cat) REF materialization: {time.time()-t0:.1f}s/eval")
    return 0


# --------------------------------------------------------------------- fiber
def run_fiber(a):
    """One S_double value + optional cross-checks at (p, s1) — probe tool."""
    ctx.prec = DC_PREC
    p, s1 = frac(a.p), frac(a.s1)
    groups = grouped(mutate_counts(a.mutate)) if a.mutate else GROUPS
    t0 = time.time()
    pn = PNode(p, groups, "cat", a.B2)
    t1 = time.time()
    S = S_double(pn, s1)
    row = dict(kind="fiber-double", p=a.p, s1=a.s1, B2=a.B2,
               mutate=a.mutate or None,
               S_mid=arb(S.mid()).str(25), S_rad=arb(S.rad()).str(5, radius=False),
               logS=(S.log().str(25) if bool(S > 0) else "n/a"),
               secs_rule=round(t1 - t0, 1), secs_eval=round(time.time() - t1, 2),
               dc_prec=DC_PREC, maxrss_mb=rss_mb())
    bank(a.tag, row)
    return 0


# ----------------------------------------------------------------------- cli
def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("selftest")
    sub.add_parser("rulecheck")
    sub.add_parser("timeeval")
    q = sub.add_parser("fiber")
    q.add_argument("--p", required=True)
    q.add_argument("--s1", required=True)
    q.add_argument("--B2", default="3.5")
    q.add_argument("--mutate", default=None)
    q.add_argument("--tag", required=True)
    for nm in ("uncpanel2", "balpanel2"):
        q = sub.add_parser(nm)
        q.add_argument("--pm", type=float, required=True)
        q.add_argument("--pr", type=float, required=True)
        q.add_argument("--K", type=int, required=True)
        if nm == "uncpanel2":
            q.add_argument("--shi", type=float, default=3.5)
            q.add_argument("--maxrev", type=int, default=240)
        else:
            q.add_argument("--thi", type=float, default=1.75)
            q.add_argument("--dhi", type=float, default=3.5)
            q.add_argument("--maxrev", type=int, default=120)
        q.add_argument("--floor", required=True)
        q.add_argument("--tolo", default="3e-2")
        q.add_argument("--mutate", default=None)
        q.add_argument("--tag", required=True)
    a = ap.parse_args()
    if a.cmd in ("selftest", "rulecheck", "timeeval"):
        return dict(selftest=selftest, rulecheck=rulecheck,
                    timeeval=timeeval)[a.cmd]()
    return dict(fiber=run_fiber, uncpanel2=run_uncpanel2,
                balpanel2=run_balpanel2)[a.cmd](a)


if __name__ == "__main__":
    sys.exit(main())

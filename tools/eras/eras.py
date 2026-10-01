#!/usr/bin/env python3
"""eras.py — ERAS (Exact Remainder-Aware Series): mean-value / Taylor-Lagrange
enclosure forms in a shared parameter p.

The failure mode this tool exists to tame is interval dependency: when one
parameter enters every factor of a large product, naive ball evaluation
correlates nothing and the enclosure explodes. Measured on the reference
integrand shipped here (the collapsed arm at fiber s=(0.1,0.2,0.3), p0=1/3):
naive-ball enclosure relwidth 8.84e22 at p-ball radius 1e-4, on an object whose
true sensitivity dlogF/dp is O(1e2-1e3). (The blowup at fixed radius is steeply
fiber-dependent — measured span at r=1e-4 from 4.4e6 to 3.0e41 across thin
fibers — so 8.84e22 at this fiber is the reference number, and the verifier
regression-gates it on every run.) Both reference arms are provided in centered
form in p:

  (1) G_unc(p; s1,s2)  — the y3-collapsed uncorrected assembly with exact s3
      closure (full or windowed): a 319-coefficient polynomial in u = y3^2
      contracted against exact closure weights;
  (2) L_corr(p; s1,s2,s3) — the corrected 4-dim integrand in product form,
      prod_x pat_x^{n_x} * (1 - P0000)^{-N}.

Mean-value form: for a p-ball P with exact dyadic center m = mid(P),

      F(P)  ⊆  F(m) + F'(P) * (P - m)

(mean value theorem for real P; for complex/acb P the segment-integral form
F(p)-F(m) = (p-m) * int_0^1 F'(m+t(p-m)) dt with the convex ball enclosing the
segment gives the same enclosure for holomorphic F). F' is computed by
FORWARD-MODE differentiation of the SAME assembly: every scalar carries its
d/dp partner (class D below); the degree-318 polynomial layer carries (V, V')
poly pairs multiplied by the product rule ((V1,D1)*(V2,D2) -> (V1V2,
V1D2+D1V2)), so flint's fast polynomial multiplication is retained. The
derivative evaluated on the FULL ball only enters multiplied by the small
radius (P - m). K=1 mean-value is rigorous but can miss practical width
targets (measured); the order-K Taylor-Lagrange forms (class S below,
taylor_encl/taylor_ladder) carry the p-correlation to order K — order-K with
panel half-width <= ~3e-3 is the measured working configuration.

Subcommands:
  selftest — (a) T1/T2 anchors at point-p through the dual path, BOTH arms
             (unc via the collapsed 319-coefficient polynomial; corr via the
             4-dim product), gate 1e-30 (reference precision ~5e-51);
             forward-mode derivative vs central finite difference (h = 2^-64,
             prec 512), gate 1e-20; off-anchor collapsed-vs-product
             cross-check; (c) planted-error 0001:55->56 at anchor +
             fiber-shift-equality level through the dual path.
  blowup   — (b) the blowup-fix table: for p-ball radii {1e-4,1e-3,1e-2,5e-2}
             centered at p0 = 1/3, fixed thin (s1,s2)=(0.1,0.2) (T1 fiber),
             record naive-ball relwidth vs mean-value relwidth, CONTAINMENT of
             50 point evaluations across each ball, and the measured dlogF/dp;
             plus the planted-error shift through the mean-value enclosure at
             radius 1e-4. Both arms; corr also at 48-bit.
  taylork  — order-K Taylor-Lagrange relwidth vs (radius, K), both arms, with
             containment against the best-K enclosure.
  timeeval — dual-assembly cost vs plain value-only assembly (the ~2x claim).

Input constants must be constructed at working precision: a 53-bit
default-precision constant degrades a 4.7e-3 enclosure to ~2e9 (measured; the
verifier's precision probe gates this pathology). p0 must be exactly
representable or its representation error absorbed outward. Non-finite
enclosures (pole-adjacent arms at fat radii) are no-information, never a pass.

Data file QUARTET_COLLAPSE.json (the pinned 318-count pattern data shared with
tools/posq) is looked up next to this script; ERAS_BASE overrides the
directory, and a missing file raises at import — there is no fallback. Results
JSON go to $TMPDIR/eras_results (override: ERAS_RESULTS).

Any adaptation of these enclosure forms must be run through the shipped
battery (verify_eras_adversarial.py, beside this file) before it is used —
see that file's docstring for the consumer protocol.
"""
import json, math, os, resource, sys, time
from flint import arb, acb, arb_poly, acb_poly, ctx

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.environ.get("ERAS_BASE", HERE)
RESULTS = os.environ.get("ERAS_RESULTS", os.path.join(
    os.environ.get("TMPDIR", "/tmp"), "eras_results"))
COUNTS = json.load(open(f"{BASE}/QUARTET_COLLAPSE.json"))["pattern_counts"]
LAM = 10
COLLAPSE_PREC = int(os.environ.get("COLLAPSE_PREC", "192"))
BINOM = {}

ANCHORS = {  # PINNED_SPEC 50-digit anchors: name -> ((s1,s2,s3,pi1), corr, unc)
    "T1": (("1/10", "2/10", "3/10", "1/3"),
           "-854.30882122149070633100757217796207967803476410107557",
           "-976.09285161889490586499770937640740253545699456512481"),
    "T2": (("1/2", "1/4", "1/8", "2/5"),
           "-789.41183749194623446733341102350522257245211548545322",
           "-842.55310669243368577614288803238346376995933055630806"),
}

def frac(s):
    if isinstance(s, (arb, acb)):
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

def mutate_counts(spec):
    c = dict(COUNTS)
    pat, _, dv = spec.partition(":")
    c[pat] = c[pat] + int(dv)
    return c

# ------------------------------------------------ forward-mode dual numbers
class D:
    """Dual number (v, d) = (value, d/dp), components arb or acb balls.
    All mixed ops coerce scalars via v*0+x (never rely on flint __rmul__)."""
    __slots__ = ("v", "d")

    def __init__(self, v, d=None):
        self.v = v
        self.d = v * 0 if d is None else d

    def _c(self, x):
        if isinstance(x, D):
            return x
        z = self.v * 0
        return D(z + x, z)

    def __add__(self, o):
        o = self._c(o); return D(self.v + o.v, self.d + o.d)
    __radd__ = __add__

    def __sub__(self, o):
        o = self._c(o); return D(self.v - o.v, self.d - o.d)

    def __rsub__(self, o):
        o = self._c(o); return D(o.v - self.v, o.d - self.d)

    def __mul__(self, o):
        o = self._c(o); return D(self.v * o.v, self.d * o.v + self.v * o.d)
    __rmul__ = __mul__

    def __truediv__(self, o):
        o = self._c(o)
        q = self.v / o.v
        return D(q, (self.d - q * o.d) / o.v)

    def __rtruediv__(self, o):
        return self._c(o).__truediv__(self)

    def __neg__(self):
        return D(-self.v, -self.d)

    def exp(self):
        e = self.v.exp(); return D(e, e * self.d)

    def log(self):
        return D(self.v.log(), self.d / self.v)

    def __pow__(self, n):
        n = int(n)
        return D(self.v ** n, n * (self.v ** (n - 1)) * self.d)

def dconst(x, like):
    z = like * 0
    return D(z + x, z)

# --------------------------------------------------------- shared machinery
def _PmatD(pi0, p, y):
    """CTMC transition matrix in y = exp(-beta*t); entries Dual."""
    return [[pi0 + p * y, p - p * y], [pi0 - pi0 * y, p + pi0 * y]]

def _core(pD, s1, s2):
    one = pD._c(1)
    pi0 = one - pD
    beta = one / ((pi0 * pD) * 2)
    s1D, s2D = pD._c(s1), pD._c(s2)
    y1 = (-(beta * s1D)).exp()
    y2 = (-(beta * s2D)).exp()
    return one, pi0, beta, y1, y2

# ---------------------------------------- (1) y3-collapsed uncorrected arm
def pattern_ac_D(pD, s1, s2, groups=GROUPS):
    """Per grouped pattern: (n, a, c) duals with pat = a + c*y3^2 exactly."""
    one, pi0, beta, y1, y2 = _core(pD, s1, s2)
    y12 = y1 * y2
    P1, P2, P12 = _PmatD(pi0, pD, y1), _PmatD(pi0, pD, y2), _PmatD(pi0, pD, y12)
    cpref = pi0 * pD * y12
    out = []
    for (ab, xC, xD_), n in groups:
        xA, xB = ab
        Q0 = P12[0][xC] * (P2[0][0] * P1[0][xA] * P1[0][xB]
                           + P2[0][1] * P1[1][xA] * P1[1][xB])
        Q1 = P12[1][xC] * (P2[1][0] * P1[0][xA] * P1[0][xB]
                           + P2[1][1] * P1[1][xA] * P1[1][xB])
        a = (pi0 * Q0 + pD * Q1) * (pi0 if xD_ == 0 else pD)
        c = cpref * (Q0 - Q1)
        out.append((n, a, c if xD_ == 0 else -c))
    return out, beta

def upoly_pair(pD, s1, s2, groups=GROUPS, deriv=True):
    """The 319-coefficient polynomial in u = y3^2 as a (value, d/dp) poly
    pair; product rule at the poly level keeps flint fast multiplication."""
    ac, beta = pattern_ac_D(pD, s1, s2, groups)
    complexq = isinstance(pD.v, acb) or isinstance(frac(s1), acb) \
        or isinstance(frac(s2), acb)
    Poly = acb_poly if complexq else arb_poly
    facs = []
    for (n, a, c) in ac:
        bn = binomrow(n)
        onev = a.v * 0 + 1
        apow, cpow = [onev], [onev]
        for _ in range(n):
            apow.append(apow[-1] * a.v)
            cpow.append(cpow[-1] * c.v)
        vco, dco = [], []
        for j in range(n + 1):
            vco.append(bn[j] * apow[n - j] * cpow[j])
            if deriv:
                dv = a.v * 0
                if j < n:
                    dv += (bn[j] * (n - j)) * apow[n - j - 1] * cpow[j] * a.d
                if j > 0:
                    dv += (bn[j] * j) * apow[n - j] * cpow[j - 1] * c.d
                dco.append(dv)
        facs.append((Poly(vco), Poly(dco) if deriv else None))
    while len(facs) > 1:                      # size-sorted product tree
        facs.sort(key=lambda q: q[0].length())
        V1, D1 = facs.pop(0)
        V2, D2 = facs.pop(0)
        facs.append((V1 * V2, (V1 * D2 + D1 * V2) if deriv else None))
    return facs[0], beta

def wvecD_full(beta, J, lam=LAM):
    """w_j = int_0^inf lam e^{-lam s} y3^{2j} ds = lam/(lam + 2j*beta), dual."""
    return [lam / (beta * (2 * j) + lam) for j in range(J + 1)]

def wvecD_win(beta, J, A, B, lam=LAM):
    """Windowed exact closure over s3 in [A,B] (B None/'inf' -> infinity)."""
    Aa = frac(A)
    Bv = None if B in (None, "inf") else frac(B)
    out = []
    for j in range(J + 1):
        r = beta * (2 * j) + lam
        t = (-(r * Aa)).exp()
        if Bv is not None:
            t = t - (-(r * Bv)).exp()
        out.append((t * lam) / r)
    return out

def G_unc_D(pD, s1, s2, groups=GROUPS, s3win=None, prec=None):
    """s3-closed collapsed uncorrected object, dual: sum_j C_j(p,s1,s2) w_j(p).
    Assembly at `prec` (default COLLAPSE_PREC)."""
    old = ctx.prec
    ctx.prec = prec or COLLAPSE_PREC
    try:
        (V, Dp), beta = upoly_pair(pD, s1, s2, groups)
        J = V.length() - 1
        wv = wvecD_full(beta, J) if s3win is None else \
            wvecD_win(beta, J, s3win[0], s3win[1])
        totv = V[0] * wv[0].v
        totd = Dp[0] * wv[0].v + V[0] * wv[0].d
        for j in range(1, J + 1):
            totv += V[j] * wv[j].v
            totd += Dp[j] * wv[j].v + V[j] * wv[j].d
        return D(totv, totd)
    finally:
        ctx.prec = old

def G_unc_val(p, s1, s2, groups=GROUPS, s3win=None, prec=None):
    """Plain value-only assembly for cost comparison."""
    old = ctx.prec
    ctx.prec = prec or COLLAPSE_PREC
    try:
        (V, _), beta = upoly_pair(D(p), s1, s2, groups, deriv=False)
        J = V.length() - 1
        wv = wvecD_full(beta, J) if s3win is None else \
            wvecD_win(beta, J, s3win[0], s3win[1])
        tot = V[0] * wv[0].v
        for j in range(1, J + 1):
            tot += V[j] * wv[j].v
        return tot
    finally:
        ctx.prec = old

# --------------------------------------- (2) corrected 4-dim integrand arm
def L_prod_D(pD, s1, s2, s3, counts=None, corrected=True):
    """Product-form 4-dim integrand, dual in p.
    corrected=True multiplies by (1 - P0000)^{-N}, N = sum of counts."""
    counts = COUNTS if counts is None else counts
    one, pi0, beta, y1, y2 = _core(pD, s1, s2)
    y3 = (-(beta * pD._c(s3))).exp()
    P1, P2, P3 = _PmatD(pi0, pD, y1), _PmatD(pi0, pD, y2), _PmatD(pi0, pD, y3)
    P12 = _PmatD(pi0, pD, y1 * y2)
    P123 = _PmatD(pi0, pD, y1 * y2 * y3)
    pi = [pi0, pD]

    def pat(x):
        t = pD * 0
        for r in (0, 1):
            for w in (0, 1):
                for v in (0, 1):
                    t = t + (pi[r] * P3[r][w] * P2[w][v] * P1[v][x[0]]
                             * P1[v][x[1]] * P12[w][x[2]] * P123[r][x[3]])
        return t

    val = one
    ntot = 0
    for pstr, n in counts.items():
        if n == 0:
            continue
        ntot += n
        val = val * (pat(tuple(int(ch) for ch in pstr)) ** n)
    if corrected:
        val = val * ((one - pat((0, 0, 0, 0))) ** (-ntot))
    return val

# --------------------------------------- order-K Taylor form (series in p)
class S:
    """Truncated power series in t = p - center, coefficients arb/acb balls.
    With a BALL constant term, coefficient k encloses F^{(k)}(xi)/k! for every
    xi in the ball (all ops are the exact Taylor-coefficient recurrences), so
    the order-K Taylor-Lagrange form is rigorous. K=1 == mean-value form."""
    __slots__ = ("c",)

    def __init__(self, c):
        self.c = c

    def _z(self):
        return self.c[0] * 0

    def _cs(self, x):
        z = self._z()
        return S([z + x] + [z] * (len(self.c) - 1))

    _c = _cs   # duck-type D's coercion name so _core/_PmatD/pattern_ac_D work

    def __add__(self, o):
        o = o if isinstance(o, S) else self._cs(o)
        return S([a + b for a, b in zip(self.c, o.c)])
    __radd__ = __add__

    def __sub__(self, o):
        o = o if isinstance(o, S) else self._cs(o)
        return S([a - b for a, b in zip(self.c, o.c)])

    def __rsub__(self, o):
        return self._cs(o).__sub__(self)

    def __neg__(self):
        return S([-a for a in self.c])

    def __mul__(self, o):
        o = o if isinstance(o, S) else self._cs(o)
        K = len(self.c) - 1
        z = self._z()
        out = [z for _ in range(K + 1)]
        for i, a in enumerate(self.c):
            for j in range(K + 1 - i):
                out[i + j] += a * o.c[j]
        return S(out)
    __rmul__ = __mul__

    def recip(self):
        K = len(self.c) - 1
        b0 = 1 / self.c[0]
        out = [b0]
        for n in range(1, K + 1):
            s = self._z()
            for k in range(1, n + 1):
                s += self.c[k] * out[n - k]
            out.append(-b0 * s)
        return S(out)

    def __truediv__(self, o):
        o = o if isinstance(o, S) else self._cs(o)
        return self * o.recip()

    def __rtruediv__(self, o):
        return self._cs(o).__truediv__(self)

    def exp(self):
        K = len(self.c) - 1
        out = [self.c[0].exp()]
        for n in range(1, K + 1):
            s = self._z()
            for k in range(1, n + 1):
                s += (k * self.c[k]) * out[n - k]
            out.append(s / n)
        return S(out)

    def powi(self, n):
        if n < 0:
            # recip FIRST (constant term still tight), then positive power:
            # powi-then-recip can hit 1/(ball containing 0) = NaN once the
            # powered constant-term relwidth exceeds 1.
            return self.recip().powi(-n)
        r, b = self._cs(1), self
        while n:
            if n & 1:
                r = r * b
            b = b * b
            n >>= 1
        return r

def s_var(x0, K):
    """The series of p itself around the center: x0 + 1*t."""
    z = x0 * 0
    return S([x0, z + 1] + [z] * (K - 1))

def L_prod_S(pS, s1, s2, s3, counts=None, corrected=True):
    """Product-form 4-dim integrand as an order-K series in p."""
    counts = COUNTS if counts is None else counts
    one = pS._cs(1)
    pi0 = one - pS
    beta = one / ((pi0 * pS) * 2)
    ys = [(-(beta * t)).exp() for t in (frac(s1), frac(s2), frac(s3))]
    y1, y2, y3 = ys
    P1, P2, P3 = (_PmatD(pi0, pS, y) for y in ys)
    P12 = _PmatD(pi0, pS, y1 * y2)
    P123 = _PmatD(pi0, pS, y1 * y2 * y3)
    pi = [pi0, pS]

    def pat(x):
        t = pS._cs(0)
        for r in (0, 1):
            for w in (0, 1):
                for v in (0, 1):
                    t = t + (pi[r] * P3[r][w] * P2[w][v] * P1[v][x[0]]
                             * P1[v][x[1]] * P12[w][x[2]] * P123[r][x[3]])
        return t

    val = one
    ntot = 0
    for pstr, n in counts.items():
        if n == 0:
            continue
        ntot += n
        val = val * pat(tuple(int(ch) for ch in pstr)).powi(n)
    if corrected:
        val = val * (one - pat((0, 0, 0, 0))).powi(-ntot)
    return val

def G_unc_series(pS, s1, s2, groups=GROUPS, s3win=None):
    """y3-collapsed uncorrected s3-closed object as an order-K series in p:
    the 319 u-coefficients carried as series (dense convolution over S),
    contracted against the series closure weights w_j(p)."""
    ac, beta = pattern_ac_D(pS, s1, s2, groups)   # generic: works on S
    one = pS._cs(1)
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

    while len(polys) > 1:                      # size-sorted product tree
        polys.sort(key=len)
        polys.append(pmul(polys.pop(0), polys.pop(0)))
    C = polys[0]
    J = len(C) - 1
    if s3win is None:
        wv = [LAM / (beta * (2 * j) + LAM) for j in range(J + 1)]
    else:
        wv = []
        Aa, Bv = frac(s3win[0]), (None if s3win[1] in (None, "inf")
                                  else frac(s3win[1]))
        for j in range(J + 1):
            rr = beta * (2 * j) + LAM
            t = (-(rr * Aa)).exp()
            if Bv is not None:
                t = t - (-(rr * Bv)).exp()
            wv.append((t * LAM) / rr)
    tot = C[0] * wv[0]
    for j in range(1, J + 1):
        tot = tot + C[j] * wv[j]
    return tot

def taylor_ladder(fn_series, pball, Ks):
    """One point series + one ball series at max(Ks); the order-K enclosure
    for every K in Ks reuses them (coefficient k of a truncated series does
    not depend on the truncation order beyond k)."""
    Kmax = max(Ks)
    m = arb(pball.mid())
    Sm = fn_series(s_var(m, Kmax))
    SB = fn_series(s_var(pball, Kmax))
    dm = pball - m
    encls = {}
    for K in Ks:
        acc = SB.c[K]
        for i in range(K - 1, -1, -1):
            acc = acc * dm + Sm.c[i]
        encls[K] = acc
    return encls

def taylor_encl(fn_series, pball, K):
    """Order-K Taylor-Lagrange enclosure: sum_{i<K} f_i(m) dm^i + f_K(P) dm^K,
    f_i from the point series at m, f_K from the ball series (encloses
    F^{(K)}(xi)/K! for all xi in P). K=1 == the pinned mean-value form."""
    m = arb(pball.mid())
    Sm = fn_series(s_var(m, K))
    SB = fn_series(s_var(pball, K))
    dm = pball - m
    acc = SB.c[K]
    for i in range(K - 1, -1, -1):
        acc = acc * dm + Sm.c[i]
    return acc, Sm, SB

# ------------------------------------------------------- mean-value wrapper
def mv_form(fn, pball):
    """fn: Dual -> Dual. Returns (mv_enclosure, naive_ball, F0_dual, Fball_dual):
    mv = F(m) + F'(pball)*(pball - m), m = exact mid(pball). The naive ball
    F(pball) comes free from the same dual call (its value component)."""
    m = arb(pball.mid()) if isinstance(pball, arb) else \
        acb(arb(pball.real.mid()), arb(pball.imag.mid()))
    one = m * 0 + 1
    F0 = fn(D(m, one))
    FB = fn(D(pball, pball * 0 + 1))
    mv = F0.v + FB.d * (pball - m)
    return mv, FB.v, F0, FB

def relwidth(b):
    b = b.real if isinstance(b, acb) else b
    if b.mid() == 0:
        return float("inf")
    return float(arb(b.rad()) / abs(arb(b.mid())))

# ------------------------------------------------------------------ output
def bankrow(tag, row):
    os.makedirs(RESULTS, exist_ok=True)
    json.dump(row, open(os.path.join(RESULTS, f"{tag}.json"), "w"), indent=1)
    print(json.dumps(row), flush=True)

def rss_mb():
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform == "darwin":
        rss /= 1024  # ru_maxrss is bytes on macOS, KiB on Linux
    return round(rss / 1024, 1)

# ---------------------------------------------------------------- selftest
def selftest():
    ok = True
    t00 = time.time()
    print("== eras selftest: anchors + derivative + plants through the "
          "dual path ==", flush=True)

    # (a) T1/T2 anchors at point-p, BOTH arms, prec 384
    ctx.prec = 384
    print("\n-- (a) anchors, prec 384, gate 1e-30 (banked ~5e-51) --")
    for name, ((a, b, c, pp), corr50, unc50) in ANCHORS.items():
        s1, s2, s3, p = frac(a), frac(b), frac(c), frac(pp)
        pD = D(p, p * 0 + 1)
        # unc arm via collapsed 319-coefficient polynomial (dual value layer)
        (V, Dp), beta = upoly_pair(pD, s1, s2)
        assert V.length() == NTOT + 1 and Dp.length() >= NTOT, \
            (V.length(), Dp.length())
        u = ((-(beta.v * s3)).exp()) ** 2
        lv = V(u).log()
        d1 = abs(lv - arb(unc50))
        hit1 = bool(d1 < arb("1e-30")) and bool(arb(lv.rad()) < arb("1e-40"))
        print(f"{'OK  ' if hit1 else 'FAIL'} {name} unc/collapsed: "
              f"|dlog anchor| < {d1.str(3)}, rad {arb(lv.rad()).str(3, radius=False)}")
        # corr arm via 4-dim product form (dual)
        lc = L_prod_D(pD, s1, s2, s3, corrected=True).v.log()
        d2 = abs(lc - arb(corr50))
        hit2 = bool(d2 < arb("1e-30")) and bool(arb(lc.rad()) < arb("1e-40"))
        print(f"{'OK  ' if hit2 else 'FAIL'} {name} corr/product : "
              f"|dlog anchor| < {d2.str(3)}, rad {arb(lc.rad()).str(3, radius=False)}")
        # bonus: unc via product form (same dual scalar path as corr)
        lu = L_prod_D(pD, s1, s2, s3, corrected=False).v.log()
        d3 = abs(lu - arb(unc50))
        hit3 = bool(d3 < arb("1e-30"))
        print(f"{'OK  ' if hit3 else 'FAIL'} {name} unc/product  : "
              f"|dlog anchor| < {d3.str(3)}")
        ok = ok and hit1 and hit2 and hit3

    # off-anchor cross-check: collapsed dual value vs product dual value
    p, s1, s2, s3 = arb("0.22"), arb("0.07"), arb("0.9"), arb("0.41")
    pD = D(p, p * 0 + 1)
    (V, Dp), beta = upoly_pair(pD, s1, s2)
    u = ((-(beta.v * s3)).exp()) ** 2
    dd = abs(V(u).log() - L_prod_D(pD, s1, s2, s3, corrected=False).v.log())
    hit = bool(dd < arb("1e-70"))
    print(f"{'OK  ' if hit else 'FAIL'} off-anchor collapsed vs product: "
          f"|dlog| < {dd.str(3)}")
    ok = ok and hit

    # derivative check: forward-mode vs central finite difference, prec 512
    print("\n-- forward-mode derivative vs central FD (h=2^-64), prec 512, "
          "gate 1e-20 --")
    ctx.prec = 512
    h = arb(2) ** -64
    p0 = arb(1) / 3
    s1, s2, s3 = arb("0.1"), arb("0.2"), arb("0.3")
    cases = [
        ("G_unc (collapsed, full s3 closure)",
         lambda q: G_unc_D(D(q, q * 0 + 1), s1, s2, prec=512)),
        ("L_corr (4-dim product)",
         lambda q: L_prod_D(D(q, q * 0 + 1), s1, s2, s3, corrected=True)),
    ]
    for label, fn in cases:
        dual = fn(p0)
        fd = (fn(p0 + h).v - fn(p0 - h).v) / (2 * h)
        # gate: FD-vs-dual relative gap inside the h^2 truncation envelope
        # (~1e-33; balls themselves are far tighter than that, so ball
        # overlap is NOT expected) and the dual-derivative ball tight.
        rd = abs(arb(fd.mid()) - arb(dual.d.mid())) / abs(arb(dual.d.mid()))
        tight = bool(arb(dual.d.rad()) / abs(arb(dual.d.mid())) < arb("1e-30"))
        hit = bool(rd < arb("1e-20")) and tight
        sens = arb(dual.d.mid()) / arb(dual.v.mid())
        print(f"{'OK  ' if hit else 'FAIL'} {label}: rel FD-dual diff "
              f"{rd.str(3)} (gate 1e-20), dual-d relrad "
              f"{(arb(dual.d.rad()) / abs(arb(dual.d.mid()))).str(3)}; "
              f"dlogF/dp = {sens.str(6)} (true sensitivity)")
        ok = ok and hit

    # (c) planted error 0001: 55->56 through the new path
    print("\n-- (c) planted error 0001:+1, prec 384 --")
    ctx.prec = 384
    mut = mutate_counts("0001:+1")
    gmut = grouped(mut)
    for name, ((a, b, c, pp), corr50, unc50) in ANCHORS.items():
        s1, s2, s3, p = frac(a), frac(b), frac(c), frac(pp)
        pD = D(p, p * 0 + 1)
        (Vm, _), beta = upoly_pair(pD, s1, s2, gmut)
        u = ((-(beta.v * s3)).exp()) ** 2
        sh_u = Vm(u).log() - arb(unc50)
        sh_c = L_prod_D(pD, s1, s2, s3, mut, corrected=True).v.log() - arb(corr50)
        m1 = bool(abs(sh_u) > arb("1e-3"))
        m2 = bool(abs(sh_c) > arb("1e-3"))
        print(f"{'OK  ' if m1 else 'FAIL'} {name} unc/collapsed moves {sh_u.str(6)}")
        print(f"{'OK  ' if m2 else 'FAIL'} {name} corr/product  moves {sh_c.str(6)}")
        ok = ok and m1 and m2
    # fiber-shift equality: collapsed vs product carry the SAME plant shift
    s1, s2, s3, p = frac("1/10"), frac("2/10"), frac("3/10"), frac("1/3")
    pD = D(p, p * 0 + 1)
    (V0, _), beta = upoly_pair(pD, s1, s2)
    (Vm, _), _ = upoly_pair(pD, s1, s2, gmut)
    u = ((-(beta.v * s3)).exp()) ** 2
    shc = Vm(u).log() - V0(u).log()
    shd = (L_prod_D(pD, s1, s2, s3, mut, corrected=False).v.log()
           - L_prod_D(pD, s1, s2, s3, corrected=False).v.log())
    eq = bool(abs(shc - shd) < arb("1e-40")) and bool(abs(shc) > arb("1e-2"))
    print(f"{'OK  ' if eq else 'FAIL'} fiber shift equality: collapsed "
          f"{shc.str(20)} vs product {shd.str(20)}, |diff| < "
          f"{abs(shc - shd).str(3)}")
    ok = ok and eq

    print(f"\nSELFTEST {'PASS' if ok else 'FAIL'} "
          f"({time.time() - t00:.1f}s, maxrss {rss_mb()} MB)", flush=True)
    return 0 if ok else 1

# ------------------------------------------------------------------ blowup
def blowup():
    """(b) containment + blowup-fix table, both arms."""
    t00 = time.time()
    radii = ["1e-4", "1e-3", "1e-2", "5e-2"]
    s1, s2, s3 = arb("0.1"), arb("0.2"), arb("0.3")   # T1 fiber, thin
    arms = [
        ("unc-collapsed", COLLAPSE_PREC,
         lambda pD: G_unc_D(pD, arb("0.1"), arb("0.2"))),
        ("corr-4dim", COLLAPSE_PREC,
         lambda pD: L_prod_D(pD, arb("0.1"), arb("0.2"), arb("0.3"),
                             corrected=True)),
        ("corr-4dim", 48,
         lambda pD: L_prod_D(pD, arb("0.1"), arb("0.2"), arb("0.3"),
                             corrected=True)),
    ]
    table = []
    allok = True
    hdr = (f"{'arm':<14} {'prec':>4} {'radius':>7} {'naive relw':>11} "
           f"{'mv relw':>10} {'target':>8} {'relw<=tgt':>9} {'contain':>7} "
           f"{'dlogF/dp':>10} {'Fp relw':>9} {'secs':>6}")
    print("== blowup-fix table: naive ball vs mean-value form in p ==")
    print(f"   center p0 = 1/3, fixed thin (s1,s2,s3) = (0.1,0.2,0.3) "
          f"[T1 fiber]; 50 containment points per radius")
    print(hdr, flush=True)
    for arm, prec, fn in arms:
        ctx.prec = prec
        p0 = arb(1) / 3
        for rs in radii:
            t0 = time.time()
            r = float(rs)
            ctx.prec = prec
            pball = arb(p0.mid(), r)
            mv, naive, F0, FB = mv_form(fn, pball)
            nrw, mrw = relwidth(naive), relwidth(mv)
            target = 1e3 * r
            # containment: 50 point evals across the ball (offsets at 99.9%
            # of the nominal radius; the stored mag radius rounds UP, so all
            # sample points are certainly inside pball — asserted)
            m = arb(pball.mid())
            cont = 0
            for i in range(50):
                off = (arb(rs) * arb("0.999")) * (2 * i - 49) / 49
                pi_ = m + off
                assert pball.contains(pi_)
                vi = fn(D(pi_, pi_ * 0 + 1)).v
                if mv.contains(vi):
                    cont += 1
            sens = float(arb(FB.d.mid()) / arb(F0.v.mid())) \
                if F0.v.mid() != 0 else float("nan")
            fprw = relwidth(FB.d)   # dependency width of F' on the full ball
            hit = (mrw <= target) and (cont == 50)
            allok = allok and hit
            secs = time.time() - t0
            print(f"{arm:<14} {prec:>4} {rs:>7} {nrw:>11.3e} {mrw:>10.3e} "
                  f"{target:>8.1e} {str(mrw <= target):>9} {cont:>4}/50 "
                  f"{sens:>10.3e} {fprw:>9.2e} {secs:>6.2f}", flush=True)
            table.append(dict(
                arm=arm, prec=prec, radius=r, naive_relwidth=nrw,
                mv_relwidth=mrw, target_relwidth=target,
                relw_le_target=bool(mrw <= target), contained=f"{cont}/50",
                dlogF_dp_mid=sens, Fprime_ball_relwidth=fprw,
                secs=round(secs, 2)))

    # planted-error shift THROUGH the mean-value enclosure (radius 1e-4)
    print("\n== planted error 0001:+1 through the mean-value path, "
          "radius 1e-4, both arms ==")
    mutg = grouped(mutate_counts("0001:+1"))
    mutc = mutate_counts("0001:+1")
    plant = []
    for arm, prec, fn0, fnm in [
        ("unc-collapsed", COLLAPSE_PREC,
         lambda pD: G_unc_D(pD, s1, s2),
         lambda pD: G_unc_D(pD, s1, s2, groups=mutg)),
        ("corr-4dim", COLLAPSE_PREC,
         lambda pD: L_prod_D(pD, s1, s2, s3, corrected=True),
         lambda pD: L_prod_D(pD, s1, s2, s3, mutc, corrected=True)),
    ]:
        ctx.prec = prec
        p0 = arb(1) / 3
        pball = arb(p0.mid(), 1e-4)
        mv0, _, _, _ = mv_form(fn0, pball)
        mvm, _, _, _ = mv_form(fnm, pball)
        if mv0.contains(arb(0)) or mvm.contains(arb(0)):
            shift = arb("nan")
            fires = False
        else:
            shift = mvm.log() - mv0.log()
            fires = bool(not shift.contains(arb(0))) and \
                bool(abs(arb(shift.mid())) > arb("0.5"))
        allok = allok and fires
        print(f"{'OK  ' if fires else 'FAIL'} {arm}: mv log-shift = "
              f"{shift.str(10)} (excludes 0, |mid|>0.5)", flush=True)
        plant.append(dict(arm=arm, prec=prec, radius=1e-4,
                          mv_log_shift=shift.str(15), fires=fires))

    row = dict(
        kind="taylor-p-blowup-fix", date_utc=time.strftime("%Y-%m-%d"),
        center_p0="1/3", fiber=dict(s1="0.1", s2="0.2", s3="0.3"),
        note=("naive_relwidth = relwidth of F evaluated with p = full ball "
              "(the banked dependency blowup); mv_relwidth = mean-value form "
              "F(p0) + F'(ball)*(ball-p0); target relwidth <= 1e3*radius = "
              "the true-sensitivity scaling; containment = 50 point evals "
              "per ball inside the mv enclosure"),
        table=table, planted=plant, all_pass=bool(allok),
        secs_total=round(time.time() - t00, 2), maxrss_mb=rss_mb())
    bankrow("TAYLOR_P_BLOWUP", row)
    print(f"\nBLOWUP {'PASS' if allok else 'FAIL'} "
          f"({time.time() - t00:.1f}s, maxrss {rss_mb()} MB)", flush=True)
    return 0 if allok else 1

# ----------------------------------------------------------------- taylork
def taylork():
    """Order-K Taylor-Lagrange forms on the corr 4-dim integrand: measure
    relwidth vs (radius, K). K=1 is the pinned mean-value form; higher K
    tests whether carrying the p-correlation to order K kills the
    dependency explosion (and how K must scale with radius)."""
    t00 = time.time()
    ctx.prec = COLLAPSE_PREC
    s1, s2, s3 = arb("0.1"), arb("0.2"), arb("0.3")
    fnS = lambda pS: L_prod_S(pS, s1, s2, s3, corrected=True)
    fnD = lambda pD: L_prod_D(pD, s1, s2, s3, corrected=True)
    ok = True

    # sanity 1: series constant term reproduces the T1 corr anchor (prec 384)
    ctx.prec = 384
    p = frac("1/3")
    c0 = L_prod_S(s_var(p, 2), s1, s2, s3, corrected=True).c[0]
    d0 = abs(c0.log() - arb(ANCHORS["T1"][1]))
    hit = bool(d0 < arb("1e-30"))
    print(f"{'OK  ' if hit else 'FAIL'} series c0 vs T1 corr anchor: "
          f"|dlog| < {d0.str(3)}")
    ok = ok and hit
    # sanity 2: K=1 Taylor form == mean-value form (same enclosure class)
    ctx.prec = COLLAPSE_PREC
    pball = arb((arb(1) / 3).mid(), 1e-4)
    e1, _, SB = taylor_encl(fnS, pball, 1)
    mv, _, _, FB = mv_form(fnD, pball)
    same = bool(e1.overlaps(mv)) and bool(SB.c[1].overlaps(FB.d))
    print(f"{'OK  ' if same else 'FAIL'} K=1 == mean-value form: relw "
          f"{relwidth(e1):.3e} vs {relwidth(mv):.3e}; f1(P) vs F'(P) overlap")
    ok = ok and same

    radii = ["1e-4", "1e-3", "1e-2", "5e-2"]
    Ks = [1, 2, 4, 8, 16]
    table = []
    arms = [
        ("corr-4dim", fnS,
         lambda pt: L_prod_D(D(pt, pt * 0 + 1), s1, s2, s3, corrected=True).v),
        ("unc-collapsed", lambda pS: G_unc_series(pS, s1, s2),
         lambda pt: G_unc_D(D(pt, pt * 0 + 1), s1, s2).v),
    ]
    for arm, fnSer, fnPt in arms:
        print(f"\n== order-K Taylor form, {arm} arm, prec {COLLAPSE_PREC}, "
              f"center 1/3, T1 fiber ==")
        print(f"{'radius':>7} " + " ".join(f"{'K='+str(K):>10}" for K in Ks)
              + f" {'target':>8}   contain(bestK)")
        for rs in radii:
            t0 = time.time()
            r = float(rs)
            ctx.prec = COLLAPSE_PREC
            pball = arb((arb(1) / 3).mid(), r)
            encls = taylor_ladder(fnSer, pball, Ks)
            rws = [relwidth(encls[K]) for K in Ks]
            bestK = Ks[min(range(len(Ks)),
                           key=lambda i: (math.isnan(rws[i]), rws[i]))]
            encl_best = encls[bestK]
            # containment: 50 point evals across the ball vs best enclosure
            m = arb(pball.mid())
            cont = 0
            for i in range(50):
                off = (arb(rs) * arb("0.999")) * (2 * i - 49) / 49
                pi_ = m + off
                assert pball.contains(pi_)
                if encl_best.contains(fnPt(pi_)):
                    cont += 1
            target = 1e3 * r
            best = relwidth(encl_best)
            hit = (not math.isnan(best)) and best <= target and cont == 50
            ok = ok and hit
            print(f"{rs:>7} " + " ".join(f"{w:>10.2e}" for w in rws)
                  + f" {target:>8.1e}   {cont}/50 {'PASS' if hit else 'FAIL'}"
                  f" ({time.time()-t0:.1f}s)", flush=True)
            table.append(dict(
                arm=arm, radius=r, relwidth_by_K=dict(zip(map(str, Ks), rws)),
                target_relwidth=target, best_K=bestK, best_relwidth=best,
                contained=f"{cont}/50", passes=bool(hit),
                secs=round(time.time() - t0, 2)))

    row = dict(kind="taylor-p-orderK", date_utc=time.strftime("%Y-%m-%d"),
               arms=[a for a, _, _ in arms], prec=COLLAPSE_PREC,
               center_p0="1/3",
               fiber=dict(s1="0.1", s2="0.2", s3="0.3"), Ks=Ks,
               note=("order-K Taylor-Lagrange in p: sum_{i<K} f_i(m) dm^i + "
                     "f_K(ball) dm^K via truncated series arithmetic; K=1 == "
                     "the pinned mean-value form; containment tested against "
                     "the best-K enclosure"),
               table=table, all_pass=bool(ok),
               secs_total=round(time.time() - t00, 2), maxrss_mb=rss_mb())
    bankrow("TAYLOR_P_ORDERK", row)
    print(f"\nTAYLORK {'PASS' if ok else 'FAIL'} "
          f"({time.time() - t00:.1f}s, maxrss {rss_mb()} MB)", flush=True)
    return 0 if ok else 1

# ---------------------------------------------------------------- timeeval
def timeeval():
    p0 = arb(1) / 3
    s1, s2 = arb("0.1"), arb("0.2")
    out = {}
    for label, call in [
        ("plain-value", lambda: G_unc_val(p0, s1, s2)),
        ("dual",        lambda: G_unc_D(D(p0, p0 * 0 + 1), s1, s2)),
    ]:
        t0 = time.time(); n = 0
        while time.time() - t0 < 3.0:
            call(); n += 1
        ms = (time.time() - t0) / n * 1e3
        out[label] = round(ms, 3)
        print(f"G_unc[{label}] @ {COLLAPSE_PREC}bit: {ms:.2f} ms/eval")
    print(f"dual/plain cost ratio: {out['dual'] / out['plain-value']:.2f}x")
    t0 = time.time(); n = 0
    pD = D(p0, p0 * 0 + 1)
    while time.time() - t0 < 2.0:
        L_prod_D(pD, s1, s2, arb("0.3"), corrected=True); n += 1
    print(f"L_corr[dual] @ {ctx.prec}bit(ctx): {(time.time()-t0)/n*1e3:.3f} ms/eval")
    return 0

# --------------------------------------------------------------------- cli
def main():
    cmd = sys.argv[1] if len(sys.argv) > 1 else "selftest"
    return dict(selftest=selftest, blowup=blowup, taylork=taylork,
                timeeval=timeeval)[cmd]()

if __name__ == "__main__":
    sys.exit(main())

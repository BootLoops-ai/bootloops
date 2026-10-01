r"""box1_dilog.py — nestor example / dispersion KERNEL: the LBL3SE one-loop box
(3m+1M corner) in 1D-reduced dilogarithmic form, Vieta-stabilized small roots.

PROVENANCE: the math shipped here is the Vieta-stabilized branch of a
reference implementation that is not distributed (src-sha256
a4435b0841c5308eff1a2c53cc9b7b3656a9c6698036d6927f3a01e128c3dd36) —
small roots extracted as mal = r1/(M5*malp), mbe = r2/(M5*mbep) instead of
(p_k - sd_k)/(2 M5). The stabilization is essential: the direct
(p_k - sd_k)/(2 M5) extraction is a catastrophic cancellation losing
~log10(M5/(2 m2)) digits as M5 grows — invisible inside the dps+18 guard
band for M5 <= 5000, breaching a >=48d bar at dps=50 for M5 >~ 1e24
(measured, with a mutation control). See ../GUIDE.md, "examples/box1_dilog"
section.

SAME integral as box1_closed.box_closed (pySecDec normalization, factor 1):

    Box = int_simplex 1/F^2 ,  F = (M5 x0 + m2 x1 + m2 x2 + m2 x3) U
                                   - s x0 x2 - u x1 x3 ,
    U = x0+x1+x2+x3 ,  u = -s-t ,  p_i^2 = 0 .

DERIVATION (2D -> 1D):
  Cheng-Wu x3=1; lines 1,2,3 share mass m2 so the (x1,x2) quadratic part of F
  is the DEGENERATE form m2(x1+x2)^2. With v = x1+x2, w = x1, F is linear in
  w; the w- then v-integrals close, leaving
      Box = int_0^inf dx0 int_0^inf dv  v / (Q1 Q2),
      Q_k = M5 x0^2 + p_k x0 + r_k  (as quadratics in x0),
      p1 = a v + M5+m2,  r1 = m2 (v+1)^2,  a = M5+m2-s,
      p2 = (M5+m2)(v+1), r2 = r1 - u v.
  Doing dx0 analytically (4-root partial fractions; the linear pole x0=u/s
  cancels since Q1(u/s)=Q2(u/s)) leaves a v-integrand whose complex
  singularities sit at v=-1 and on |v|=1 (roots of r2) — M5-INDEPENDENT — so
  one fixed Gauss-Legendre rule in v converges uniformly in M5. (The dv-first
  variant is singular at x0=-m2/M5, collapsing onto the endpoint as M5->inf:
  6d at M5=500. Do not resurrect it.)
  In the deep-Euclidean region both discriminants are sums of non-negative
  terms, so all four roots are real negative: pure-real mpf arithmetic, four
  mp.log per node.

VIETA STABILIZATION: the SMALL
root of Q_k suffers cancellation in (p_k - sd_k)/(2 M5) once 4 M5 r_k <<
p_k^2 (large M5). Vieta gives it exactly as r_k / (M5 * large_root) — no
subtraction — so the rule is accurate for ARBITRARILY large M5.

REGION: s<0, t<0, m2>0, M5>0, u=-s-t < 4 m2 (sharp; ValueError past the
pseudo-threshold, same guard as box1_closed).

API: box_closed_dilog(s, t, m2, M5, dps=40) — drop-in for
box1_closed.box_closed, ~48x faster at equal precision; alias box_closed;
make_K(s, t, m2, dps) — memoized K(w')=Box(s,t;m2,M5=w') kernel factory for
nestor.dispersion panels (disp_subtracted). Restores mp.mp.dps on exit including error paths.

`python3 tools/nestor/examples/box1_dilog.py` reproduces the validation table
(5-point M5 gate vs the reference 2D oracle) + the M5-uniformity probe;
without BOX1_DILOG_ORACLE_DIR it prints SKIP and exits 0.
"""
import mpmath as mp

_GL_CACHE = {}


def _gl(n, work):
    """Cached Gauss-Legendre nodes/weights on [-1,1] at working precision."""
    key = (n, work)
    xw = _GL_CACHE.get(key)
    if xw is None:
        old = mp.mp.dps
        mp.mp.dps = work
        xw = mp.gauss_quadrature(n, 'legendre')
        mp.mp.dps = old
        _GL_CACHE[key] = xw
    return xw


def _g_of_v(v, s, u, m2, M5, a):
    """g(v) = int_0^inf dx0 (1/Q2-1/Q1)/(u-sx0) via 4-root partial fractions.

    All four roots real negative in the deep-Euclidean region; -rho > 0 is
    used directly so every log is real. Small roots by Vieta (r_k/(M5*big)),
    NOT by (p_k-sd_k)/(2M5): that difference IS the canonical-vs-reference drift
    (see the provenance header)."""
    vp1 = v + 1
    p1 = a * v + (M5 + m2)
    r1 = m2 * vp1 * vp1
    p2 = (M5 + m2) * vp1
    r2 = r1 - u * v
    sd1 = mp.sqrt(p1 * p1 - 4 * M5 * r1)
    sd2 = mp.sqrt(p2 * p2 - 4 * M5 * r2)
    twoM5 = 2 * M5
    # large roots (addition, stable); small roots by Vieta: alpha*alpha' = r1/M5
    malp = (p1 + sd1) / twoM5                                  # = -alpha'
    mal = r1 / (M5 * malp) if malp != 0 else (p1 - sd1) / twoM5    # = -alpha
    mbep = (p2 + sd2) / twoM5                                  # = -beta'
    mbe = r2 / (M5 * mbep) if mbep != 0 else (p2 - sd2) / twoM5    # = -beta
    # residues (linear pole at x0=u/s cancels; sum of residues = 0):
    return -((-1 / ((u + s * mal) * sd1)) * mp.log(mal)
             + (1 / ((u + s * malp) * sd1)) * mp.log(malp)
             + (1 / ((u + s * mbe) * sd2)) * mp.log(mbe)
             + (-1 / ((u + s * mbep) * sd2)) * mp.log(mbep))


def box_closed_dilog(s, t, m2, M5, dps=40):
    """eps^0 of the pySecDec box (factor 1), >= dps digits, 1D GL quadrature,
    convergence AND floating-point accuracy uniform in M5."""
    work = dps + 18
    old = mp.mp.dps
    mp.mp.dps = work
    try:
        s = mp.mpf(s); t = mp.mpf(t); m2 = mp.mpf(m2); M5 = mp.mpf(M5)
        u = -s - t
        a = M5 + m2 - s
        if not (s < 0 and t < 0 and m2 > 0 and M5 > 0) or u >= 4 * m2:
            raise ValueError(
                "box_closed_dilog: need deep-Euclidean s<0, t<0, m2>0, M5>0 "
                "with u=-s-t < 4 m2 (pseudo-threshold not supported)")
        n = work + 12
        nodes, weights = _gl(n, work)
        L = mp.mpf(2)
        one = mp.mpf(1)
        tot = mp.mpf(0)
        for tk, wk in zip(nodes, weights):
            v = L * (one + tk) / (one - tk)
            dv = L * 2 / (one - tk) ** 2
            tot += wk * _g_of_v(v, s, u, m2, M5, a) * dv
        res = +tot
    finally:
        mp.mp.dps = old
    return res


box_closed = box_closed_dilog          # drop-in alias


def make_K(s, t, m2, dps):
    """Memoized kernel factory K(w') = Box(s,t;m2,M5=w') for nestor.dispersion
    panels (cache key mp.nstr(w', dps+8); canonical-source pattern)."""
    cache = {}

    def K(w):
        key = mp.nstr(w, dps + 8)
        r = cache.get(key)
        if r is None:
            r = box_closed_dilog(s, t, m2, w, dps)
            cache[key] = r
        return r
    K.cache = cache
    return K


if __name__ == '__main__':
    import os, sys, time
    # Reference leg: the 2D validation oracle (box1_closed.py) is not
    # distributed with this package. Point BOX1_DILOG_ORACLE_DIR at a
    # directory holding box1_closed.py to reproduce the validation table.
    _od = os.environ.get('BOX1_DILOG_ORACLE_DIR', '')
    if not _od or not os.path.isfile(os.path.join(_od, 'box1_closed.py')):
        print('SKIP: reference 2D oracle not available '
              '(set BOX1_DILOG_ORACLE_DIR to run the validation table)')
        sys.exit(0)
    sys.path.insert(0, _od)
    import box1_closed as BC

    s, m2 = mp.mpf(-1), mp.mpf(1)
    mp.mp.dps = 80; t = mp.mpf(-1) / 3
    print("=== box_closed_dilog (member) vs reference box1_closed (2D oracle), dps=50 ===")
    _ = box_closed_dilog(s, t, m2, mp.mpf(2), dps=50)   # warm GL caches
    _ = BC.box_closed(s, t, m2, mp.mpf(2), dps=50)
    for M5 in ['2', '5', '50', '500', '0.1']:
        M5f = mp.mpf(M5)
        t0 = time.time(); v1 = box_closed_dilog(s, t, m2, M5f, dps=50); dt1 = time.time() - t0
        t0 = time.time(); v2 = BC.box_closed(s, t, m2, M5f, dps=50);    dt2 = time.time() - t0
        mp.mp.dps = 60
        d = -mp.log10(abs((v1 - v2) / v2)) if v1 != v2 else mp.mpf(99)
        print(f"M5={M5:>6}  {mp.nstr(v1, 50)}  agree {mp.nstr(d, 4)}d  "
              f"{dt1*1e3:6.1f}ms vs {dt2*1e3:6.1f}ms ({dt2/max(dt1,1e-9):4.1f}x)")
    print("(M5-uniformity probe:)")
    for M5 in ['5000', '0.01', '1']:
        v1 = box_closed_dilog(s, t, m2, mp.mpf(M5), dps=50)
        v2 = BC.box_closed(s, t, m2, mp.mpf(M5), dps=50)
        mp.mp.dps = 60
        d = -mp.log10(abs((v1 - v2) / v2)) if v1 != v2 else mp.mpf(99)
        print(f"  M5={M5:>6}: agree {mp.nstr(d, 4)}d")

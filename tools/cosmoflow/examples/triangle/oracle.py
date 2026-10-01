#!/usr/bin/env python3
r"""oracle.py -- WORKED EXAMPLE (cosmoflow.examples.triangle): INDEPENDENT
numerical oracle for the FRW 1-loop TRIANGLE elliptic sub-sector and the
two-site tree chain (arXiv:2408.16386, 2312.05303).  Problem-specific by
design; the general front door is cosmoflow.alphabet_graph / polytope.baikov_B.

Evaluates
    I_{00tau}[N] = \int_Gamma  B(y;X)^{-1/2+eps} * N(y) / q_G^tau  d^{n_s}y,
    Gamma = { y_e >= 0,  B >= 0 },   q_G linear in one y_e,
    B = polytope.B_z (Cayley--Menger).

n_s = 3 (triangle) is the ONLY fully-implemented case:

DOMAIN (n_s=3, derived):
    B is quadratic in each z_e = y_e^2 with leading coeff  -2*X_opp^2 < 0.
    disc_{z1}(B) = 4 * lambda(X1^2,X2^2,X3^2) * lambda(z2, z3, X3^2)
    => for lambda_X < 0  (X-triangle-ineq HOLDS):  inner z1-interval nonempty
       <=>  |y23 - y31| <= X3 <= y23 + y31.
    Gamma is a NON-compact tube along the z-diagonal; convergence at eps=0
    needs deg(N) < tau.

METHOD (n_s=3, eps=0): nested quadrature, exact-weight inner layer.
    inner y12  -- K(m)-SPLIT closed form (log-carrying piece = complete-K,
                  bounded remainder via periodic-trapezoid on phi in [0,pi/2])
    middle y31 -- tanh-sinh, split at y31=X1 (isolated interior log-sing at
                  (y23,y31)=(X2,X1) pushed to a ts corner)
    outer y23  -- tanh-sinh on [0,inf), split at {X2,X3}
    l-space cross-check: I = 1/(2*sqrt(8)) * int_{R^3} N/(q^tau*prod y) d^3 l
    (valid for lambda_X<0).

METHOD (eps!=0): inner + middle via tanh-sinh (algebraic-endpoint safe).

n_s = 2 (two-site chain, TREE): twisted 1D-reduced oracle F(X1,X2,Y;eps).
  F(eps) = -2Y * int_0^inf int_0^inf dx1 dx2 (x1 x2)^eps
             / ((x1+x2+X1+X2)(x1+X1+Y)(x2+X2+Y)),   -1 < eps <= 0.
  Inner x2 closed form: int_0^inf x2^eps dx2/((x2+A)(x2+B))
      = pi*csc(pi*eps)*(A^eps-B^eps)/(A-B),  A = x1+X1+X2, B = X2+Y
  (eps->0 recovers the log-rational reduction log(A/B)/(A-B)).
  Gates: 28/28 eps=0 receipts byte-exact + 3pts x 3eps
  vs print 2312.05303 eq.(3.42) >=40d.
  NOTE this is the TREE-chain wavefunction coefficient (poles never carry
  c: at tree level a connected subgraph cannot omit an edge without
  disconnecting), NOT the one-LOOP bubble masters.

n_s >= 4: NotImplementedError stub (see notes at bottom).
"""
import time
import mpmath as mp
import sympy as sp

from ...polytope import B_z, z1, z2, z3, X1, X2, X3


# ---- fixed spectral rules (cached) ---------------------------------------
_CHEB = {}
def _cheb(N):
    """Chebyshev-1 nodes cos(theta_j) and theta_j, j=0..N-1."""
    key = (N, mp.mp.dps)
    if key not in _CHEB:
        th = [mp.pi*(2*j+1)/(2*N) for j in range(N)]
        _CHEB[key] = (th, [mp.cos(t) for t in th], [mp.sin(t) for t in th])
    return _CHEB[key]

_TS = {}
def _ts_nodes(a, b, degree=6):
    """Tanh-sinh nodes+weights on [a,b] (hand-rolled, cached per dps).  Used by
    the corner fallback in inner_vec."""
    key = (degree, mp.mp.dps)
    if key not in _TS:
        h = mp.mpf(1)/mp.mpf(2)**(degree-1)
        eps = mp.mpf(10)**(-mp.mp.dps)
        pts = []; k = 0
        while True:
            t = k*h
            s = mp.pi/2 * mp.sinh(t)
            x = mp.tanh(s)
            w = h*(mp.pi/2)*mp.cosh(t)/mp.cosh(s)**2
            pts.append((x, w))
            if k > 0: pts.append((-x, w))
            if w < eps and k > 4: break
            k += 1
        _TS[key] = pts
    ctr, R = (a+b)/2, (b-a)/2
    return [ctr+R*x for x, w in _TS[key]], [R*w for x, w in _TS[key]]


_PHI = {}
def _phi_nodes(N):
    """Midpoint nodes on [0,pi/2]: return (h, [cos^2 phi_j], [sin^2 phi_j])."""
    key = (N, mp.mp.dps)
    if key not in _PHI:
        h = (mp.pi/2)/N
        ph = [h*(j+mp.mpf(1)/2) for j in range(N)]
        _PHI[key] = (h, [mp.cos(p)**2 for p in ph], [mp.sin(p)**2 for p in ph])
    return _PHI[key]


# ---- z1-quadratic coefficients (mpmath-callable) -------------------------
_p1 = sp.Poly(B_z, z1)
_A1e, _B1e, _C1e = _p1.all_coeffs()
_A1 = sp.lambdify((z2, z3, X1, X2, X3), _A1e, 'mpmath')
_B1 = sp.lambdify((z2, z3, X1, X2, X3), _B1e, 'mpmath')
_C1 = sp.lambdify((z2, z3, X1, X2, X3), _C1e, 'mpmath')

# C1 is quadratic in z3; its complex-conjugate
# roots are the C1=0 branches.  For X1>X3 kinematics they run INTERIOR to the
# middle y31-contour at small Im (p1-type points: 0.03-0.28), stalling the
# middle ts-quad at ~1e-22 for ALL y23.  Fix: split at Re(sqrt(z3root)).
def _pinch_y31(y23, Xv):
    X1v, X2v, X3v = Xv
    z2v = y23*y23
    a = -2*X2v**2
    b = 2*(X1v**2*X2v**2 + X1v**2*z2v - X2v**4 + X2v**2*X3v**2
           + X2v**2*z2v - X3v**2*z2v)
    c = -2*X1v**2*(X1v**2*z2v + X2v**2*X3v**2 - X2v**2*z2v - X3v**2*z2v
                   + z2v*z2v)
    disc = b*b - 4*a*c
    out = []
    if disc < 0:
        sq = mp.sqrt(mp.mpc(disc))
        for z3r in ((-b + sq)/(2*a), (-b - sq)/(2*a)):
            y31r = mp.sqrt(z3r)
            re, im = mp.re(y31r), abs(mp.im(y31r))
            if im < mp.mpf('0.5'):
                out.append((re, im))
    return out

# ---- 9 masters (numerator N, tau, degN)  per 2408.16386 elliptic_sector --
MASTERS = {
    'e1': (lambda y23,y31: y23*y31,   1, 2),
    'e2': (lambda y23,y31: y23,       1, 1),
    'e3': (lambda y23,y31: y23,       2, 1),
    'e4': (lambda y23,y31: y31,       1, 1),
    'e5': (lambda y23,y31: y31,       2, 1),
    'e6': (lambda y23,y31: mp.mpf(1), 2, 0),
    'e7': (lambda y23,y31: mp.mpf(1), 1, 0),
    'e8': (lambda y23,y31: y23**2,    1, 2),
    'e9': (lambda y23,y31: y31**2,    1, 2),
}
# UV-finite at eps=0  <=>  deg(N) < tau
FINITE_EPS0 = [k for k, (N, tau, d) in MASTERS.items() if d < tau]


# ---- inner y12 layer: K(m)-split closed form + bounded remainder ----------
def inner_vec(y23, y31, Xv, eps, Nch, taus=(1, 2), cs=(1, 2)):
    """Return {(tau,c): int_{y12 in Gamma} B^{-1/2+eps}/q^tau dy12} in one pass."""
    X1v, X2v, X3v = Xv
    z2v, z3v = y23*y23, y31*y31
    A  = _A1(z2v, z3v, *Xv); Bc = _B1(z2v, z3v, *Xv); Cc = _C1(z2v, z3v, *Xv)
    disc = Bc*Bc - 4*A*Cc
    zero = {(t, c): mp.mpf(0) for t in taus for c in cs}
    if disc <= 0:
        return zero
    sq = mp.sqrt(disc)
    z1m, z1p = sorted([(-Bc + sq)/(2*A), (-Bc - sq)/(2*A)])
    if z1p <= 0:
        return zero
    negA = -A; S = X1v + X2v + X3v
    if z1m <= 0 or eps != 0:    # rare corner / eps!=0: per-(tau,c) tanh-sinh
        out = {}
        lo = mp.sqrt(max(0, z1m)) if z1m >= 0 else mp.mpf(0)
        hi = mp.sqrt(z1p)
        for t in taus:
            for c in cs:
                def f(y, t=t, c=c):
                    Bv = negA*(y*y - z1m)*(z1p - y*y)
                    return mp.power(Bv, mp.mpf(-1)/2 + eps)/mp.power(S + c*y, t)
                out[(t, c)] = mp.quad(f, [lo, hi])
        return out
    b = mp.sqrt(z1p)
    ratio = z1m/z1p
    if ratio < mp.mpf(10)**(-mp.mp.dps):
        ratio = mp.mpf(10)**(-mp.mp.dps)   # clamp; K ~ 1/2 log(16/ratio)
    # midpoint-in-phi remainder converges as
    # exp(-4*Nch*sqrt(ratio)) -> ALGEBRAIC (not spectral) near the (X2,X1)
    # corner where ratio->0 — an unsplit rule biases e6,e7 by ~1e-8
    # (Nch-dependent).
    # Route alpha: threshold 0.20 and ADAPTIVE
    # ts degree (7 below ratio=1e-3, else 6): deg=6 self-converges to only ~21d
    # at ratio<1e-4 (measured); deg=7 gives >=33d everywhere.
    if ratio < mp.mpf('0.20'):
        a = mp.sqrt(z1m)
        nodes, wts = _ts_nodes(a, b, degree=(7 if ratio < mp.mpf('1e-3') else 6))
        out = {(t, c): mp.mpf(0) for t in taus for c in cs}
        for x, w in zip(nodes, wts):
            Bv = negA*(x - a)*(x + a)*(b - x)*(b + x)
            if Bv <= 0:
                continue
            iw = w/mp.sqrt(Bv)
            for c in cs:
                iq = 1/(S + c*x); p = iw
                for t in range(1, max(taus)+1):
                    p *= iq
                    if t in taus:
                        out[(t, c)] += p
        return out
    # phi-substitution:  y(phi)=sqrt(z1m*cos^2+z1p*sin^2),  B^{-1/2}dy = dphi/(sqrt(negA)*y).
    # Split off the log-carrying piece in CLOSED FORM via K(m):
    #   inner = (1/sqrt(negA)) * [ K(1 - z1m/z1p)/(b*S^tau)  +  int_0^{pi/2} g1(phi) dphi ],
    #   g1 = (1/q^tau - 1/S^tau)/y   -- BOUNDED uniformly in z1m -> midpoint-in-phi spectral.
    m = 1 - ratio
    K = mp.ellipk(m)
    invSnA = 1/mp.sqrt(negA)
    Kpc = K/b
    N = Nch
    hph, cs2, sn2 = _phi_nodes(N)
    acc = {(t, c): mp.mpf(0) for t in taus for c in cs}
    Spow = {t: S**t for t in taus}
    for cp2, sp2 in zip(cs2, sn2):
        y2 = z1m*cp2 + z1p*sp2
        y  = mp.sqrt(y2)
        for c in cs:
            q = S + c*y
            for t in taus:
                acc[(t, c)] += (1/mp.power(q, t) - 1/Spow[t])/y
    out = {}
    for t in taus:
        for c in cs:
            out[(t, c)] = invSnA*(Kpc/Spow[t] + hph*acc[(t, c)])
    return out


# ---- middle + outer: mp.quad tanh-sinh with corner splits -----------------
def middle_vec(y23, Xv, eps, Nch, mdeg, keys, cs, disk_R=None):
    """Return {(key,c): int dy31 N*inner_{tau,c}}.  tanh-sinh, split at X1.
    If disk_R is set and |y23-X2|<disk_R, EXCLUDE [X1-h, X1+h],
    h=sqrt(R^2-(y23-X2)^2) (polar disk carve-out; see disk_corner)."""
    X1v, X2v, X3v = Xv
    lo, hi = abs(y23 - X3v), y23 + X3v
    if hi <= lo:
        return {(k, c): mp.mpf(0) for k in keys for c in cs}
    taus = tuple(sorted({MASTERS[k][1] for k in keys}))
    cache = {}
    def inner_cached(y31):
        if y31 not in cache:
            cache[y31] = inner_vec(y23, y31, Xv, eps, Nch, taus, cs)
        return cache[y31]
    segs = None
    if disk_R is not None:
        dy = y23 - X2v
        if abs(dy) < disk_R:
            h = mp.sqrt(disk_R*disk_R - dy*dy)
            segs = [(lo, X1v - h), (X1v + h, hi)]
    if segs is None:
        pts = [lo] + ([X1v] if lo < X1v < hi else []) + [hi]
        segs = [(pts[i], pts[i+1]) for i in range(len(pts)-1)]
    # interior C1-pinch splits (X1>X3 kinematics)
    pinches = [re for (re, im) in _pinch_y31(y23, Xv)]
    if pinches:
        segs2 = []
        for (aa, bb) in segs:
            cuts = sorted({p for p in pinches if aa < p < bb})
            edges = [aa] + cuts + [bb]
            segs2 += [(edges[i], edges[i+1]) for i in range(len(edges)-1)]
        segs = segs2
    out = {}
    for k in keys:
        Nfun, tau, degN = MASTERS[k]
        for c in cs:
            out[(k, c)] = sum(mp.quad(
                lambda y31: Nfun(y23, y31)*inner_cached(y31)[(tau, c)],
                [aa, bb], maxdegree=mdeg) for aa, bb in segs)
    return out


# ---- Route beta: polar carve-out at the (X2,X1) codim-2 log --------------
# The C1=0 curve has complex-conjugate branches meeting at exactly ONE real
# point (y23,y31)=(X2,X1) => isolated 2D log singularity.  Fixed 1D ts splits
# are only first-order-adapted to it (caps mm-cross-check at 16-22d).  Fix:
# R-disk in polar coords (r ts-spectral incl. r*log r endpoint; theta periodic
# trapezoidal, nested => free Nth/2 self-cert) + clip disk out of middle/outer.
def disk_corner(Xv, eps, Nch, keys, cs, R, mdeg_r, Nth):
    """{(key,c): int_{disk R @ (X2,X1)} N*inner dy23 dy31} via polar coords.
    Returns (out_Nth, out_Nth2) for theta self-certification.  Nth even.
    theta-rate is point-dependent (measured 0.23-0.51 d/node); check the pair."""
    X1v, X2v, X3v = Xv
    taus = tuple(sorted({MASTERS[k][1] for k in keys}))
    assert Nth % 2 == 0
    hth = 2*mp.pi/Nth
    acc  = {(k, c): mp.mpf(0) for k in keys for c in cs}
    acc2 = {(k, c): mp.mpf(0) for k in keys for c in cs}
    for j in range(Nth):
        th = hth*j
        cth, sth = mp.cos(th), mp.sin(th)
        icache = {}
        def ic(r):
            if r not in icache:
                icache[r] = inner_vec(X2v+r*cth, X1v+r*sth, Xv, eps, Nch, taus, cs)
            return icache[r]
        for k in keys:
            Nfun, tau, _ = MASTERS[k]
            for c in cs:
                v = mp.quad(lambda r: r*Nfun(X2v+r*cth, X1v+r*sth)*ic(r)[(tau, c)],
                            [0, R], maxdegree=mdeg_r)
                acc[(k, c)] += v
                if j % 2 == 0:
                    acc2[(k, c)] += v
    return ({kc: hth*acc[kc] for kc in acc},
            {kc: 2*hth*acc2[kc] for kc in acc2})


def eval_all_triangle(a, lam, eps=0, dps=50, Nch=48, mdeg_mid=None,
                      mdeg_out=None, keys=None, cs=(1, 2), verbose=False,
                      disk_R='auto', Nth=128):
    """Evaluate all requested masters x c-values at X=(a*lam,lam,1).
    Returns ({(key,c): (val, digits_est)}, wall_s).
    disk_R: polar carve-out radius at (X2,X1); 'auto' -> 0.1*min(X1,X2)
    clipped inside the domain; None -> disable (caps ~16-22d)."""
    if keys is None:
        keys = FINITE_EPS0 if eps == 0 else list(MASTERS)
    mp.mp.dps = dps + 15
    if mdeg_mid is None: mdeg_mid = max(5, int(mp.log(dps, 2)) + 2)
    if mdeg_out is None: mdeg_out = max(5, int(mp.log(dps, 2)) + 2)
    Xv = (mp.mpf(a)*mp.mpf(lam), mp.mpf(lam), mp.mpf(1))
    X1v, X2v, X3v = Xv
    if disk_R == 'auto':
        Rmax = min(X1v - abs(X2v - X3v), (X2v + X3v) - X1v,
                   X2v, abs(X3v - X2v)) / 2
        disk_R = min(mp.mpf('0.1')*min(X1v, X2v), Rmax)
        if disk_R <= 0:
            disk_R = None
    if disk_R is not None and disk_R <= 0:
        disk_R = None
    splits = sorted({mp.mpf(0), X2v, X3v} |
                    ({X2v - disk_R, X2v + disk_R} if disk_R else set())) + [mp.inf]
    t0 = time.time()
    if disk_R is None:
        disk = {(k, c): mp.mpf(0) for k in keys for c in cs}
    else:
        disk, disk_half = disk_corner(Xv, eps, Nch, keys, cs, disk_R,
                                      min(mdeg_mid, 5), Nth)
        if verbose:
            dmin = min((float(-mp.log10(abs(disk[kc]-disk_half[kc])/abs(disk[kc])))
                        if disk[kc] != disk_half[kc] and disk[kc] != 0 else 999.0)
                       for kc in disk)
            print(f"  [disk_corner R={mp.nstr(disk_R,6)} Nth={Nth} "
                  f"theta_selfcert>={dmin:.1f}d] wall={time.time()-t0:.1f}s",
                  flush=True)
    cache = {}
    def mid_cached(y23):
        if y23 not in cache:
            cache[y23] = middle_vec(y23, Xv, eps, Nch, mdeg_mid, keys, cs,
                                    disk_R=disk_R)
        return cache[y23]
    out = {}
    for k in keys:
        for c in cs:
            v, e = disk[(k, c)], mp.mpf(0)
            for aa, bb in zip(splits[:-1], splits[1:]):
                vv, ee = mp.quad(lambda y23: mid_cached(y23)[(k, c)],
                                 [aa, bb], maxdegree=mdeg_out, error=True)
                v += vv; e += ee
            dig = float(-mp.log10(e/abs(v))) if (e > 0 and v != 0) else dps
            out[(k, c)] = (v, dig)
            if verbose:
                print(f"  {k} c={c}: {mp.nstr(v, min(dps, int(dig)+2))}  [~{dig:.1f}d]",
                      flush=True)
    wall = time.time() - t0
    return out, wall


# ---- l-space independent cross-check (n_s=3, eps=0, lambda_X<0 only) ------
def eval_master_lspace(key, a, lam, c, dps=30):
    """I = 1/(2*sqrt(8)) * int_{R^3} N(|k|)/(q^tau * |k1||k2||k3|) d^3 l.
    Cylindrical: l=(lx,ly,lz), P1=(X1,0,0), P2=(p2x,p2y,0).  Symmetric in lz."""
    Nfun, tau, degN = MASTERS[key]
    with mp.workdps(dps + 10):
        X1v, X2v, X3v = mp.mpf(a)*mp.mpf(lam), mp.mpf(lam), mp.mpf(1)
        S = X1v + X2v + X3v
        P1P2 = (X3v**2 - X1v**2 - X2v**2)/2
        gramP = X1v**2*X2v**2 - P1P2**2
        assert gramP > 0, "l-space rep needs X-triangle-ineq (lambda_X<0)"
        p2x = P1P2/X1v
        p2y = mp.sqrt(gramP)/X1v
        Q1 = (mp.mpf(0), mp.mpf(0))
        Q2 = (X1v, mp.mpf(0))
        Q3 = (X1v + p2x, p2y)
        def integ(lx, ly, lz):
            lz2 = lz*lz
            y31 = mp.sqrt((lx + Q1[0])**2 + (ly + Q1[1])**2 + lz2)
            y12 = mp.sqrt((lx + Q2[0])**2 + (ly + Q2[1])**2 + lz2)
            y23 = mp.sqrt((lx + Q3[0])**2 + (ly + Q3[1])**2 + lz2)
            q = S + c*y12
            return Nfun(y23, y31) / (mp.power(q, tau) * y12*y23*y31)
        def flz(lx, ly):
            return 2*mp.quad(lambda lz: integ(lx, ly, lz), [0, mp.inf])
        def fly(lx):
            return mp.quad(lambda ly: flz(lx, ly), [-mp.inf, mp.inf])
        val = mp.quad(fly, [-mp.inf, mp.inf])
        return val / (2*mp.sqrt(8))


# ===========================================================================
# n_s = 2: two-site chain twisted oracle
#
# F_twosite eps!=0 branch + F1_twosite + print form: gated vs the 2D-direct
# arbiter; F_twosite eps=0 branch: the log-rational reduction, gated vs the
# 28-pt receipt set; certification pattern: two-precision x two-level,
# quoted = min.
# Physical region: X1 > Y > 0, X2 > 0 (reduced integrand pole-free); the
# object is SYMMETRIC under (X1<->X2), so swap arguments if X2 > Y > X1.
# Twist domain: -1 < eps <= 0 (eps=0 = dS log branch).
# ===========================================================================
def _exact_mpf(v):
    """Exact-input conversion (Fraction/int/str/mpf); floats REFUSED (purity
    bar -- certified chains take exact rationals).  Call INSIDE workdps."""
    if isinstance(v, float):
        raise ValueError("float input refused: pass Fraction/int/str (exact)")
    if hasattr(v, 'numerator') and hasattr(v, 'denominator'):
        return mp.mpf(v.numerator) / v.denominator
    return mp.mpf(v)


def F_twosite(X1, X2, Y, eps=0, dps=50, level=10):
    """Two-site chain FRW coefficient F(X1,X2,Y;eps), -1 < eps <= 0.
    1D-reduced twisted oracle (closed-form inner x2 integral); exact inputs.
    Returns an mpf at working precision dps+15 (quote via F_twosite_certified)."""
    with mp.workdps(dps + 15):
        x1v, x2v, yv = _exact_mpf(X1), _exact_mpf(X2), _exact_mpf(Y)
        ev = _exact_mpf(eps)
        if not (-1 < ev <= 0):
            raise ValueError(f"eps={ev} outside the twist domain -1 < eps <= 0")
        if not (x1v > yv > 0 and x2v > 0):
            raise ValueError("physical region X1 > Y > 0, X2 > 0 required "
                             "(F is (X1<->X2)-symmetric: swap if X2 > Y > X1)")
        B = x2v + yv
        if ev == 0:
            # eps=0 branch: the log-rational reduction.
            def g(t):
                x = t / (1 - t)
                jac = 1 / (1 - t) ** 2
                A = x + x1v + x2v
                return jac * mp.log(A / B) / ((x + x1v + yv) * (x + x1v - yv))
            v = mp.quad(g, [0, 1], maxdegree=level)
            return -2 * yv * v
        pref = -2 * yv * mp.pi * mp.csc(mp.pi * ev)

        def g(t):
            x = t / (1 - t)
            jac = 1 / (1 - t) ** 2
            A = x + x1v + x2v
            # the denominator includes the x2-independent propagator
            # (x1+X1+Y): (x1+X1+Y)*(A-B), with A-B = x1+X1-Y; gated by the
            # 2D-direct arbiter.
            return jac * x ** ev * (A ** ev - B ** ev) / ((x + x1v + yv) * (x + x1v - yv))

        return pref * mp.quad(g, [0, 1], maxdegree=level)


def F1_twosite(X1, X2, Y, dps=50, level=10):
    """eps^1 Taylor coefficient of F(X1,X2,Y;eps) at eps=0 (analytic leg,
    log^2-rational 1D integrand).
    The eps^0 coefficient is F_twosite(..., eps=0)."""
    with mp.workdps(dps + 15):
        x1v, x2v, yv = _exact_mpf(X1), _exact_mpf(X2), _exact_mpf(Y)
        if not (x1v > yv > 0 and x2v > 0):
            raise ValueError("physical region X1 > Y > 0, X2 > 0 required")
        B = x2v + yv
        lB = mp.log(B)

        def g(t):
            x = t / (1 - t)
            jac = 1 / (1 - t) ** 2
            A = x + x1v + x2v
            lA = mp.log(A)
            # full denominator (x1+X1+Y)(x1+X1-Y), as in F_twosite
            return jac * ((lA ** 2 - lB ** 2) / 2 + mp.log(x) * (lA - lB)) / ((x + x1v + yv) * (x + x1v - yv))

        return -2 * yv * mp.quad(g, [0, 1], maxdegree=level)


def F_twosite_certified(X1, X2, Y, eps=0, dps_lo=50, dps_hi=60,
                        lev_lo=8, lev_hi=10, which='F'):
    """Certified two-site value: two-precision x two-level, quoted digits =
    min(level-pair agreement at each dps, cross-dps agreement).
    which: 'F' (value at eps; eps=0 = eps^0 coefficient) or 'F1' (eps^1 leg).
    Returns {'value': str@dps_hi, 'certified_digits', 'level_agree',
    'cross_dps_agree', 'wall_s'}."""
    fn = {'F': (lambda d, l: F_twosite(X1, X2, Y, eps, d, l)),
          'F1': (lambda d, l: F1_twosite(X1, X2, Y, d, l))}[which]
    t0 = time.time()
    vals = {}
    for dps in (dps_lo, dps_hi):
        v1 = fn(dps, lev_lo)
        v2 = fn(dps, lev_hi)
        with mp.workdps(dps_hi + 15):
            d = abs(v1 - v2)
            la = int(-mp.log10(d / abs(v2))) if d != 0 else dps
        vals[dps] = (v2, la)
    with mp.workdps(dps_hi + 15):
        d = abs(vals[dps_lo][0] - vals[dps_hi][0])
        cross = int(-mp.log10(d / abs(vals[dps_hi][0]))) if d != 0 else dps_lo
    cert = min(vals[dps_lo][1], vals[dps_hi][1], cross)
    with mp.workdps(dps_hi):
        s = mp.nstr(vals[dps_hi][0], dps_hi)
    return {'value': s, 'certified_digits': int(cert),
            'level_agree': [int(vals[dps_lo][1]), int(vals[dps_hi][1])],
            'cross_dps_agree': int(cross),
            'wall_s': round(time.time() - t0, 3)}


def F_twosite_print342(X1, X2, Y, eps, dps):
    """Pinned print form 2312.05303 eq.(3.42) with (3.37) constants.
    GATE-ONLY: NEVER used by the oracles (blind-side discipline)."""
    with mp.workdps(dps + 20):
        x1v, x2v, yv = _exact_mpf(X1), _exact_mpf(X2), _exact_mpf(Y)
        e = _exact_mpf(eps)
        cZ = 4 ** (-e) * mp.sqrt(mp.pi) * mp.csc(mp.pi * e) * mp.gamma(e) * mp.gamma(mp.mpf(1) / 2 - e)
        cD = mp.pi ** 2 * mp.csc(mp.pi * e) ** 2
        h1 = mp.hyp2f1(1, e, 1 - e, (yv - x2v) / (yv + x1v))
        h2 = mp.hyp2f1(1, e, 1 - e, (yv - x1v) / (yv + x2v))
        v = (cD * (x1v + yv) ** e * (x2v + yv) ** e
             + cZ * (x1v + x2v) ** (2 * e) * (1 - h1 - h2))
    with mp.workdps(dps):
        return +v


# ---- public dispatch ------------------------------------------------------
def eval_master(n_s, key, a, lam, c, eps=0, dps=50, Nch=48, verbose=False):
    """Evaluate a single master of the n_s-site 1-loop elliptic sub-sector.
    Returns (value, digits_est, wall_s).

    n_s=3: (key, a, lam, c) as documented above (triangle masters e1..e9).
    n_s=2: the three kinematic slots are read as (X1, X2, Y) of the TREE
           two-site chain, key in {'F','F0','F1'} ('F' = value at eps;
           'F0' = eps^0 coefficient = F at eps=0; 'F1' = eps^1 coefficient);
           digits_est = certified digits (two-precision x two-level).
           The one-LOOP bubble masters (gamma=0) are NOT this evaluator --
           they run through polytope + direct quadrature (positive control:
           51.33d vs the closed Li2 form)."""
    if n_s == 3:
        r, w = eval_all_triangle(a, lam, eps, dps, Nch,
                                 keys=[key], cs=(c,), verbose=verbose)
        v, d = r[(key, c)]
        return v, d, w
    if n_s == 2:
        X1_, X2_, Y_ = a, lam, c
        if key == 'F0':
            key, eps = 'F', 0
        if key not in ('F', 'F1'):
            raise ValueError(f"n_s=2 key must be 'F', 'F0' or 'F1', got {key!r}")
        r = F_twosite_certified(X1_, X2_, Y_, eps=eps, dps_lo=dps,
                                dps_hi=dps + 10, which=key)
        with mp.workdps(dps + 10):
            v = mp.mpf(r['value'])
        return v, float(r['certified_digits']), r['wall_s']
    if n_s >= 4:
        raise NotImplementedError(
            f"n_s={n_s}: 4-site and higher need an (n_s-1)-nested outer layer "
            "and the domain analysis (which y_e-interval is compact, where the "
            "log-singular corners sit) has not been done.  Build polytope.B "
            "with baikov_B(n_s,...) and extend inner_vec/middle_vec by analogy "
            "-- the K(m)-split inner layer generalises verbatim (B quadratic "
            "in each z_e).")
    raise ValueError(f"n_s={n_s} out of range")

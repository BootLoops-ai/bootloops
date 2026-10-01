r"""
nestor.dispersion.disp_sub -- singularity-subtracted dispersion quadrature (exponentially convergent)
=====================================================================================================
REUSABLE tool for any 2-particle-reducible self-energy-insertion dispersion integral

        I = (1/pi) \int_{wmin}^{Wcut..inf} rho(w') K(w') dw'

where
  * rho(w')  = the spectral density Im Sigma(w')  -- EXPENSIVE, evaluated by quadrature; carries
               the THRESHOLD non-analyticities (turn-on at the lowest threshold, cusp/jump at each
               higher normal threshold).  These non-analyticities make a plain Gauss-Legendre or
               tanh-sinh rule converge only POLYNOMIALLY in the node count.
  * K(w')    = a SMOOTH analytic kernel (here the one-loop box Box1(s,t;M^2=w')) -- cheap, real-
               analytic on each open panel, so it has a fast-converging local Taylor series.

THE TECHNIQUE  (Kahaner/Monegato/Davis-Rabinowitz singularity subtraction, specialised to the
dispersion + analytic-kernel setting):

  On a panel [a,b] whose LEFT endpoint a = w* is a threshold, rho has a known leading expansion

        rho(w') = sum_j c_j * (w'-w*)^{alpha_j} * [log(w'-w*)]^{m_j}  +  rho_smooth(w'),     (*)

  the exponents (alpha_j, m_j) being KNOWN analytically (turn-on power, log multiplicity, jump),
  and the coefficients c_j either KNOWN in closed form or fitted from a few high-precision rho
  evaluations near w*.  Define the SUBTRACTED integrand

        g_smooth(w') = [ rho(w') - S(w') ] * K(w'),     S(w') := sum_j c_j (w'-w*)^{alpha_j} [log]^{m_j}.

  rho - S is C^infty / analytic at w* (the leading singular structure is removed), so
  \int_a^b g_smooth converges EXPONENTIALLY under a tanh-sinh (double-exponential) rule.

  The subtracted piece is ADDED BACK in closed form.  Because K is analytic we expand it about w*,
        K(w') = sum_{n>=0} K_n (w'-w*)^n      (K_n = Taylor coeffs, cheap and geometrically small),
  and then every add-back term is an EXACT MOMENT of a power-(times-log) over the panel:

        \int_a^b (w'-w*)^{alpha} [log(w'-w*)]^m K(w') dw'
              = sum_n K_n * M_{alpha+n, m}(b-a),     M_{p,m}(L) := \int_0^L u^p (log u)^m du,

  with the closed forms
        M_{p,0}(L) = L^{p+1}/(p+1),
        M_{p,1}(L) = L^{p+1}/(p+1) * ( log L - 1/(p+1) ),
        M_{p,2}(L) = L^{p+1}/(p+1) * ( (log L)^2 - 2 log L/(p+1) + 2/(p+1)^2 ),    (p != -1).
  The K_n series is truncated when K_n L^{n} drops below tolerance -- EXPONENTIALLY few terms
  because K is analytic (radius of convergence = distance to the nearest kernel singularity).

  A RIGHT-endpoint threshold (a higher cut opening from above, with the panel to its LEFT having
  rho analytic up to a finite jump at b=w*) is handled by the symmetric expansion in (w*-w').
  A pure finite JUMP at an interior threshold is handled simply by BREAKING the panel there
  (each side analytic), which the assembler does automatically.

RESULT on the LBL3SE benchmark (Im Sigma_kite, thresholds w'=1 and w'=9; see report):
  the w'=1 turn-on is  rho ~ (w'-1)(A log(w'-1) + B)  [A,B fitted / closed]; subtracting the
  A(w'-1)log(w'-1) term removes the only non-analyticity at the dominant threshold and the panel
  [1, .] integral jumps from polynomial (~6 d at 180 nodes) to exponential (tanh-sinh saturates the
  rho-node precision in O(few-x-ten) nodes).  The w'=9 elliptic cut is a finite JUMP (the equal-mass
  sunrise discontinuity is analytic-with-jump at threshold, disc_direct.py), removed by a panel break.

API  (from nestor.dispersion import ...)
  moment_powerlog(p, m, L)                -> \int_0^L u^p (log u)^m du   (closed form)
  kernel_taylor(K, wstar, nmax, h, dps)   -> [K_0..K_nmax] local Taylor coeffs of analytic kernel K
  addback_endpoint(coeffs, Kc, L, side)   -> closed-form \int of S(w')K(w') over a one-threshold panel
  tanhsinh_panel(f, a, b, dps, level)     -> tanh-sinh (double-exponential) panel quadrature
  disp_subtracted(rho, K, panels, dps)    -> the full exponentially-convergent dispersion assembly

This file is self-contained (mpmath only).  A worked analytic kernel (the one-loop box in 1D
dilogarithmic form) is ../examples/box1_dilog.py.  The application driver that wired in
imsigma_1d/imsigma_derive (rho) and box1_closed (K) with the measured w'=1, w'=9 singular data
is not included in this package.
"""
import mpmath as mp


# ----------------------------------------------------------------------------- closed-form moments
def moment_powerlog(p, m, L):
    r"""\int_0^L u^p (log u)^m du, closed form.  p > -1 (so the endpoint u=0 is integrable),
    m in {0,1,2}.  Returns an mpf at the current mp.mp.dps."""
    p = mp.mpf(p); L = mp.mpf(L)
    pp = p + 1
    Lp = L**pp
    lnL = mp.log(L)
    if m == 0:
        return Lp / pp
    if m == 1:
        return Lp / pp * (lnL - 1/pp)
    if m == 2:
        return Lp / pp * (lnL*lnL - 2*lnL/pp + 2/pp**2)
    raise ValueError("moment_powerlog supports m in {0,1,2}")


# ----------------------------------------------------------------------- analytic-kernel Taylor exp
def kernel_taylor(K, wstar, nmax, radius=None, dps=40, method='cheb'):
    r"""Local Taylor coefficients [K_0, ..., K_nmax] of the analytic kernel K about w' = wstar,
    K(wstar + u) = sum_n K_n u^n.

    method='cheb' (default, REAL-ONLY kernels):  least-squares / collocation fit of the degree-nmax
        polynomial to K sampled at nmax+1 Chebyshev points on the REAL interval [wstar-radius,
        wstar+radius] (Chebyshev-Vandermonde solve at working precision dps+extra).  Robust and uses
        only real evaluations -- the one-loop box (box1_closed) is real-only, so this is the right
        default.  Accurate to ~working precision for nmax up to ~1.5*dps when radius < analyticity rad.
    method='quad' (COMPLEX-capable kernels): Cauchy-integral (mpmath) -- more accurate at very high
        order but requires K to accept complex arguments on a contour of the given radius.

    K must be real-analytic in |w'-wstar| <= radius (radius < distance to the nearest box pseudo-
    threshold).  `radius` default 0.4."""
    wstar = mp.mpf(wstar)
    if radius is None:
        radius = mp.mpf('0.4')
    radius = mp.mpf(radius)
    if method == 'quad':
        return mp.taylor(K, wstar, nmax, method='quad', radius=radius)
    # Chebyshev collocation: sample at x_j = wstar + radius*cos(pi (j+1/2)/(N)), fit poly in u=w'-wstar
    old = mp.mp.dps
    mp.mp.dps = dps + 25
    N = nmax + 1
    us = [radius*mp.cos(mp.pi*(j + mp.mpf('0.5'))/N) for j in range(N)]
    fs = [mp.mpf(K(wstar + u)) for u in us]
    # Vandermonde V[j,n] = us[j]^n  ; solve V a = f  for Taylor coeffs a[n]
    V = mp.matrix(N, N)
    for j in range(N):
        p = mp.mpf(1)
        for n in range(N):
            V[j, n] = p
            p *= us[j]
    a = mp.lu_solve(V, mp.matrix(fs))
    mp.mp.dps = old
    return [a[n] for n in range(N)]


def addback_endpoint(coeffs, Kc, L, side='left'):
    r"""Closed-form  \int_panel S(w') K(w') dw'  for a panel of length L with the threshold at one
    endpoint, where
        S(w') = sum over terms  c * (s)^{alpha} (log s)^m ,   s = (w'-w*) [left] or (w*-w') [right],
        K(w') = sum_n Kc[n] (w'-w*)^n .
    coeffs : list of (c, alpha, m) singular-term descriptors (alpha > -1).
    Kc     : kernel Taylor coefficients about w* (from kernel_taylor); Kc[n] multiplies (w'-w*)^n.
    L      : panel length (>0).
    side   : 'left'  -> threshold at the lower endpoint a, panel = [a, a+L], s = w'-a >= 0.
             'right' -> threshold at the upper endpoint b, panel = [b-L, b], s = b-w' >= 0;
                        then (w'-w*) = -s so (w'-w*)^n contributes (-1)^n s^n.
    Returns the exact value (sum of closed-form moments) at the current precision."""
    L = mp.mpf(L)
    tot = mp.mpf(0)
    sgn = mp.mpf(1) if side == 'left' else mp.mpf(-1)
    for (c, alpha, m) in coeffs:
        c = mp.mpf(c); alpha = mp.mpf(alpha)
        for n, Kn in enumerate(Kc):
            term = c * Kn * (sgn**n) * moment_powerlog(alpha + n, m, L)
            tot += term
    return tot


# ------------------------------------------------------------------------- tanh-sinh panel (D-E)
def tanhsinh_panel(f, a, b, dps=40, level=None):
    r"""Double-exponential (tanh-sinh) quadrature of a SMOOTH (analytic-on-[a,b]) integrand f over
    [a,b].  Exponentially convergent.  Thin wrapper on mpmath.quad with method='tanh-sinh', which
    picks/grows the level adaptively to the working precision; `level` forces a fixed level if given
    (useful for convergence studies / fixed-node precompute)."""
    a = mp.mpf(a); b = mp.mpf(b)
    if level is None:
        return mp.quad(f, [a, b], method='tanh-sinh')
    return mp.quad(f, [a, b], method='tanh-sinh', maxdegree=level)


# ---------------------------------------------------------------------- full subtracted assembly
def _Smodel(coeffs, wstar, side, wp):
    s = (wp - wstar) if side == 'left' else (wstar - wp)
    if s <= 0:
        return mp.mpf(0)
    val = mp.mpf(0)
    for (c, alpha, m) in coeffs:
        t = mp.mpf(c) * s**mp.mpf(alpha)
        if m:
            t *= mp.log(s)**int(m)
        val += t
    return val


def disp_subtracted(rho, K, panels, dps=40, kernel_nmax=None, maxdegree=None):
    r"""Exponentially-convergent dispersion  (1/pi) \int rho K dw'  via singularity subtraction.

    rho(w')  : callable spectral density (expensive ok -- evaluated only at tanh-sinh abscissae).
    K(w')    : callable analytic kernel.
    panels   : list of dicts, one per panel, each:
        { 'a': lo, 'b': hi,
          'thresh':  None | ('left', w*, [(c,alpha,m),...]) | ('right', w*, [(c,alpha,m),...]),
          'sub_width': delta,   # subtraction sub-panel width (MUST be < kernel analyticity radius)
          'radius':    r,       # kernel Taylor radius for the Cauchy-integral coeffs (< analyt. rad.)
          'map':       'tail'   # optional: map [a, inf) (b ignored) via w'=a+scale*(1+t)/(1-t) }

      THRESHOLD PANEL.  rho is non-analytic at the threshold w* (left or right endpoint).  We split
      the panel into the SUBTRACTION sub-panel of width `delta` adjacent to w* and the analytic
      REMAINDER.  On the sub-panel we subtract S (the known singular model) so rho-S is analytic and
      tanh-sinh is exponential, and ADD BACK \int S K in closed form (kernel Taylor about w*, valid
      because delta < kernel radius).  The remainder is plain tanh-sinh (rho analytic there).  A pure
      finite JUMP at an interior threshold needs NO subtraction -- just give it its own panel boundary
      and tanh-sinh each side (set thresh=None on both).
    kernel_nmax : kernel Taylor truncation for the add-back (auto ~1.7*dps if None).
    maxdegree   : force a fixed tanh-sinh level (for convergence studies / fixed precompute).
    Returns I (mpf)."""
    work = dps + 20
    old = mp.mp.dps
    mp.mp.dps = work
    if kernel_nmax is None:
        kernel_nmax = int(1.7 * dps) + 10

    def _ts(f, a, b):
        if maxdegree is None:
            return mp.quad(f, [a, b], method='tanh-sinh')
        return mp.quad(f, [a, b], method='tanh-sinh', maxdegree=maxdegree)

    tot = mp.mpf(0)
    for P in panels:
        a = mp.mpf(P['a'])
        thr = P.get('thresh', None)
        if P.get('map') == 'tail':
            scale = mp.mpf(P.get('scale', a if a > 0 else mp.mpf(1)))
            def ft(t, a=a, scale=scale):
                wp = a + scale*(1+t)/(1-t)
                return rho(wp)*K(wp)*scale*2/(1-t)**2
            tot += _ts(ft, mp.mpf(-1), mp.mpf(1))
            continue
        b = mp.mpf(P['b'])
        if thr is None:
            tot += _ts(lambda wp: rho(wp)*K(wp), a, b)
            continue
        side, wstar, coeffs = thr[0], mp.mpf(thr[1]), thr[2]
        delta = mp.mpf(P.get('sub_width', min(mp.mpf('0.5'), (b - a)/2)))
        radius = P.get('radius', None)
        # subtraction sub-panel adjacent to the threshold
        if side == 'left':
            sa, sb = a, a + delta          # [w*, w*+delta]
            ra, rb = a + delta, b          # analytic remainder
        else:
            sa, sb = b - delta, b          # [w*-delta, w*]
            ra, rb = a, b - delta
        gsm = lambda wp: (rho(wp) - _Smodel(coeffs, wstar, side, wp)) * K(wp)
        tot += _ts(gsm, sa, sb)
        Kc = kernel_taylor(K, wstar, kernel_nmax, radius=radius, dps=work)
        tot += addback_endpoint(coeffs, Kc, delta, side=side)
        if rb > ra:
            tot += _ts(lambda wp: rho(wp)*K(wp), ra, rb)

    res = tot / mp.pi
    mp.mp.dps = old
    return +res

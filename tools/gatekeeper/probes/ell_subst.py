#!/usr/bin/env python3
r"""ell_subst.py - weight-matched 1D substitution selector for sqrt-product chains.

Input: sympy polynomial factors P_i (one variable), each entering the
integrand to power -1/2, and an interval [a,b] whose endpoints are simple
roots of the product.  Returns a SubstPlan with nodes K_j / weights w_j:
    int_a^b h(K) prod_i P_i(K)^(-1/2) dK  ~=  sum_j w_j h(K_j),  h analytic.

Selection (roots of prod P_i on/near [a,b]):
 * exactly the 2 endpoint roots relevant -> GAUSS-CHEBYSHEV, affine tau map
   K = c + R tau, tau_j = cos(pi(2j+1)/2N): 1/sqrt((K-a)(b-K)) IS the
   Chebyshev weight; remaining factors fold into the weights (smooth).
 * 4 relevant = 2 endpoints + one real root on EACH side (the
   pinching case: 1/sqrt-factors vanishing just outside the interval)
   -> JACOBI-SN: tau = sn(u,m), m fixed by the
   cross-ratio of the 4 roots (Moebius-symmetrized to +-1, +-sigma;
   m = 1/sigma^2), midpoint-equispaced u on (-K(m),K(m)) = trapezoid on the
   period, rate exp(-pi N K'/K).  The sn rule matches an independent
   elliptic K-line reference implementation exactly for a symmetric outer
   pair (the optional T4 acceptance leg; supply it via ELL_KLINE_DIR).
 * 3 relevant = 2 endpoints + ONE real pinching root o on a single side
   -> JACOBI-SN2 (the G-R 3.131 map; cubic_under_sqrt_K below is its h=1,
   no-smooth-factor special case): K = a + (b-a) sn^2(u,m) with
   m = (b-a)/(o-a) for a right outer root (mirrored, K = b - (b-a) sn^2,
   m = (b-a)/(b-o), for a left one), so
   dK / sqrt((K-a)(b-K)|K-o|) = 2 du / sqrt(|o - far endpoint|)
   and all three sqrt-roots cancel exactly; midpoint-equispaced u on
   (0,K(m)) - the even, 2K-periodic extension of the integrand makes this
   a 2N-point trapezoid on the full period - rate exp(-2 pi N K'/K).
 * anything else: kind='unsupported' with an honest reason - do NOT
   silently tanh-sinh it (that is the ~3 digits/level trap).  Interior
   simple root: the product changes sign, the integral is not real as
   posed.  Endpoint multiplicity > 1: the integral diverges.  Two or more
   pinching roots on one side: no rule built.

Log-endpoint c-table hook: eps-tower moments
carry log(P_endpoint) = [smooth part] + log(1-tau_f^2) with tau_f the
Gauss-Chebyshev frame variable (cheb: tau itself; sn: tau_f = sin(pi u/2K),
equispaced u <-> Chebyshev nodes).  plan.log_frame() returns per-node
(theta_j, tau_f_j, log(1-tau_f_j^2)) for a log-moment table/weights split.
"""
import mpmath as mp
import sympy as sp


class SubstPlan:
    def __init__(self, **kw):
        self.__dict__.update(kw)

    def integrate(self, h=None):
        return mp.fsum(w if h is None else w*h(K)
                       for K, w in zip(self.nodes, self.weights))

    def log_frame(self):
        """(theta_j, tau_frame_j, log(1-tau_frame_j^2)) aligned with nodes.
        Defined for the chebyshev and jacobi-sn kinds only (the jacobi-sn2
        nodes live on a half-period; no c-table frame is defined there)."""
        if self.kind not in ('chebyshev', 'jacobi-sn'):
            raise ValueError(
                f'log_frame: no c-table frame for kind={self.kind!r}')
        th = [mp.pi*(2*j + 1)/(2*self.N) for j in range(self.N)]
        sgn = 1 if self.kind == 'chebyshev' else -1
        return (th, [sgn*mp.cos(t) for t in th],
                [2*mp.log(mp.sin(t)) for t in th])


def _num(c, dps):
    if getattr(c, 'is_Rational', False):
        return mp.mpf(c.p)/mp.mpf(c.q)
    return mp.mpmathify(str(sp.N(c, dps + 15)))


def select(factors, a, b, N=64, dps=50, x=None, pinch=10.0):
    """Build the weight-matched rule; kind in {chebyshev, jacobi-sn,
    jacobi-sn2, unsupported}.  pinch: real roots within pinch*(b-a) of
    [a,b] are 'relevant' (must be absorbed by the substitution, not the
    weights)."""
    with mp.workdps(dps + 15):
        a, b = mp.mpf(a), mp.mpf(b)
        assert b > a, 'need a < b'
        factors = [sp.sympify(f) for f in factors]
        if x is None:
            xs = set().union(*[f.free_symbols for f in factors])
            assert len(xs) == 1, 'one-variable chains only'
            x = xs.pop()
        roots, lctot = [], mp.mpf(1)
        for f in factors:
            cs = [_num(c, dps) for c in sp.Poly(f, x).all_coeffs()]
            lctot *= cs[0]
            if len(cs) > 1:
                roots += mp.polyroots(cs, maxsteps=400, extraprec=120)
        scale, tol = b - a, (b - a)*mp.mpf(10)**(-dps//2)
        na = nb = 0
        outs, rest = [], []
        for r in roots:
            re, im = mp.re(r), mp.im(r)
            if abs(im) > tol:
                rest.append(r)
            elif abs(re - a) < tol:
                na += 1
            elif abs(re - b) < tol:
                nb += 1
            elif a < re < b:
                return SubstPlan(kind='unsupported', report=dict(
                    reason=f'interior root at {mp.nstr(re, 8)}'))
            elif min(abs(re - a), abs(re - b)) < pinch*scale:
                outs.append(re)
            else:
                rest.append(re)
        if (na, nb) != (1, 1):
            return SubstPlan(kind='unsupported', report=dict(
                reason=f'endpoint root multiplicities ({na},{nb}) != (1,1)'))
        c, R = (a + b)/2, (b - a)/2

        def phi(K):  # smooth co-factor, cancellation-free (root-product form)
            v = lctot
            for r in rest:
                v *= (K - r)
            return mp.re(v)
        sgn = 1 if phi(c) > 0 else -1
        if not outs:                                  # -- Gauss-Chebyshev --
            nodes = [c + R*mp.cos(mp.pi*(2*j + 1)/(2*N)) for j in range(N)]
            wts = [mp.pi/N/mp.sqrt(sgn*phi(K)) for K in nodes]
            return SubstPlan(kind='chebyshev', N=N, nodes=nodes, weights=wts,
                             a=a, b=b, m=None, report=dict(n_smooth=len(rest)))
        if len(outs) == 2 and min(outs) < a and max(outs) > b:  # -- sn --
            ol, orr = (min(outs) - c)/R, (max(outs) - c)/R
            S, P = ol + orr, ol*orr
            g = mp.mpf(0) if abs(S) < tol/scale else \
                ((1 + P) - mp.sqrt((1 + P)**2 - S*S))/S   # Moebius |g|<1
            sig = (orr - g)/(1 - g*orr)                   # outer pair -> +-sig
            m = 1/sig**2
            Kell = mp.ellipk(m)
            pref = (2*Kell/N)*mp.sqrt(1 - g*g) \
                / (R*sig*mp.sqrt((1 - g*ol)*(1 - g*orr)))
            nodes, wts = [], []
            for j in range(N):
                tp = mp.ellipfun('sn', Kell*(mp.mpf(2*j + 1)/N - 1), m)
                K = c + R*(tp + g)/(1 + g*tp)
                nodes.append(K)
                wts.append(pref/mp.sqrt(sgn*phi(K)))
            dpn = float(mp.pi*mp.ellipk(1 - m)/(Kell*mp.log(10)))
            return SubstPlan(
                kind='jacobi-sn', N=N, nodes=nodes, weights=wts, a=a, b=b,
                m=m, Kell=Kell, gamma=g, sigma=sig,
                report=dict(n_smooth=len(rest), digits_per_node=dpn,
                            logmoment_digits_per_node=dpn/2,  # c-table rate
                            credit='jacobi-sn trapezoid-on-period rule'))
        if len(outs) == 1:                            # -- one-sided sn^2 --
            o = outs[0]
            right = o > b
            m = (b - a)/((o - a) if right else (b - o))
            Kell = mp.ellipk(m)
            pref = 2*Kell/(N*mp.sqrt((o - a) if right else (b - o)))
            nodes, wts = [], []
            for j in range(N):
                sj = mp.ellipfun('sn', Kell*mp.mpf(2*j + 1)/(2*N), m)
                K = a + (b - a)*sj*sj if right else b - (b - a)*sj*sj
                nodes.append(K)
                wts.append(pref/mp.sqrt(sgn*phi(K)))
            dpn = float(2*mp.pi*mp.ellipk(1 - m)/(Kell*mp.log(10)))
            return SubstPlan(
                kind='jacobi-sn2', N=N, nodes=nodes, weights=wts, a=a, b=b,
                m=m, Kell=Kell, outer=o, side='right' if right else 'left',
                report=dict(n_smooth=len(rest), digits_per_node=dpn,
                            credit='one-sided sn^2 half-period midpoint '
                                   'rule (three roots absorbed exactly)'))
        return SubstPlan(kind='unsupported', report=dict(
            reason=f'{len(outs)} pinching outer real roots at '
                   f'{[mp.nstr(o, 8) for o in outs]} (need 0, 1 total, or '
                   '1 per side); extend the selector before farming - '
                   'do not tanh-sinh'))


def cubic_under_sqrt_K(a3, a2, a1, a0, dps=50):
    r"""int_{r1}^{r2} dz / sqrt(a3 z^3 + a2 z^2 + a1 z + a0)  in closed form
    via complete K(m).  Gradshteyn--Ryzhik 3.131:
        int_{r1}^{r2} dz / sqrt((z-r1)(r2-z)(r3-z)) = 2/sqrt(r3-r1) * K(m),
        m = (r2 - r1)/(r3 - r1),   r1 < r2 < r3 real.
    Returns (value, {'roots': (r1,r2,r3), 'm': m}) or raises if the cubic
    doesn't have 3 real roots.  a3 > 0 required (else sqrt is imaginary on
    [r1,r2]).  Origin: cosmoflow inner-y12 layer -> lifted to a generic case.
    """
    with mp.workdps(dps + 20):
        a3 = mp.mpf(a3); a2 = mp.mpf(a2); a1 = mp.mpf(a1); a0 = mp.mpf(a0)
        assert a3 > 0, "need a3 > 0 (sign flip: negate all coeffs & endpoints)"
        rs = mp.polyroots([a3, a2, a1, a0], maxsteps=400, extraprec=120)
        tol = mp.mpf(10)**(-dps//2)
        if any(abs(mp.im(r)) > tol for r in rs):
            raise ValueError("cubic_under_sqrt_K: needs 3 real roots "
                             f"(got {[mp.nstr(r, 6) for r in rs]})")
        r1, r2, r3 = sorted(mp.re(r) for r in rs)
        m = (r2 - r1)/(r3 - r1)
        val = 2/mp.sqrt(a3*(r3 - r1)) * mp.ellipk(m)
    return val, {'roots': (r1, r2, r3), 'm': m}

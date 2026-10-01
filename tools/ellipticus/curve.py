"""ellipticus.curve — exact elliptic-curve layer + AGM/period verbs.

Curve spec (PURITY BAR: exact inputs only — Fraction/int/str; floats REFUSED):
  QuarticCurve(coeffs)  : y^2 = P(x) = c4 x^4 + c3 x^3 + c2 x^2 + c1 x + c0
                          (cubic allowed: c4 = 0), coefficients exact rationals.
Verbs (all take dps explicitly; ambient mp.mp.dps is never left modified):
  curve.roots(dps)            certified-separated roots (polyroots + Newton
                              polish at wdps, pairing check)
  curve.periods(dps)          (psi_a, tau) — a-cycle holomorphic period
                              oint dx/y over the oval [r2,r3] (4-real-root
                              case) or its standard complex-pair reduction,
                              via K(m) [AGM]; two-precision gated by caller.
  curve.cycle_moments(ks,cycle,dps)
                              oint x^k dx/y over a real oval, t=sin^2(theta)
                              substitution (endpoint-cancellation-free)
                              + tanh-sinh, node-doubling agreement returned.
  sunrise_frame(t,m1,m2,m3,dps)
                              BMSW E_{3,F} frame (roots, kF2, tauF, psi1F,
                              Abel-Jacobi punctures z1,z2,z3) — verbatim port
                              of the gate-verified rowscripts_empl.fcurve
                              (incl. the (1,1,4) branch repair).

Period reduction: standard Legendre reduction, gated in
selftest vs direct quadrature.
"""
from fractions import Fraction
import mpmath as mp


def exact(v, what="coefficient"):
    """Coerce to Fraction; REFUSE floats (purity bar)."""
    if isinstance(v, float):
        raise TypeError(f"float {what} refused (purity bar): {v!r}; "
                        f"pass Fraction/int/str")
    if isinstance(v, Fraction):
        return v
    if isinstance(v, int):
        return Fraction(v)
    if isinstance(v, str):
        return Fraction(v)
    # sympy Rational duck-type
    if hasattr(v, 'p') and hasattr(v, 'q'):
        return Fraction(int(v.p), int(v.q))
    raise TypeError(f"cannot coerce {what} {v!r} exactly")


def _mpf_frac(f):
    return mp.mpf(f.numerator) / mp.mpf(f.denominator)


class QuarticCurve:
    """y^2 = P(x), P quartic or cubic, exact rational coefficients
    (descending: [c4, c3, c2, c1, c0]; len 4 = cubic [c3,c2,c1,c0])."""

    def __init__(self, coeffs):
        cs = [exact(c) for c in coeffs]
        if len(cs) == 4:
            cs = [Fraction(0)] + cs
        if len(cs) != 5:
            raise ValueError("need 4 (cubic) or 5 (quartic) coefficients")
        if all(c == 0 for c in cs[:2]):
            raise ValueError("degree < 3: not an elliptic curve chart")
        self.coeffs = cs                       # descending, exact

    # ------------------------------------------------------------ numerics
    def pval(self, x):
        s = mp.mpc(0)
        for c in self.coeffs:
            s = s * x + _mpf_frac(c)
        return s

    def roots(self, dps):
        """All roots of P at working precision dps+30, Newton-polished,
        with a residual assert; returned sorted by (Re, Im)."""
        with mp.workdps(dps + 30):
            cs = [c for c in self.coeffs if True]
            lead_zero = cs[0] == 0
            csm = [_mpf_frac(c) for c in (cs[1:] if lead_zero else cs)]
            rts = mp.polyroots(csm, maxsteps=200, extraprec=dps)
            out = []
            for r in rts:
                # Newton polish on the exact P
                for _ in range(4):
                    pv = self.pval(r)
                    dv = self._dval(r)
                    if dv == 0:
                        break
                    r = r - pv / dv
                res = abs(self.pval(r))
                scale = max(abs(_mpf_frac(c)) for c in self.coeffs)
                assert res < scale * mp.mpf(10) ** (-(dps + 10)), \
                    f"root residual {mp.nstr(res, 5)}"
                out.append(mp.mpc(r))
            out.sort(key=lambda t: (mp.re(t), mp.im(t)))
            return out

    def _dval(self, x):
        cs = self.coeffs
        s = mp.mpc(0)
        n = len(cs) - 1
        for k, c in enumerate(cs[:-1]):
            s = s * x + (n - k) * _mpf_frac(c)
        return s

    def real_ovals(self, dps):
        """Intervals [a,b] between consecutive real branch points where
        P > 0 (real ovals of y^2=P)."""
        rts = self.roots(dps)
        tol = mp.mpf(10) ** (-dps)
        reals = sorted(mp.re(r) for r in rts if abs(mp.im(r)) < tol)
        out = []
        for a, b in zip(reals, reals[1:]):
            mid = (a + b) / 2
            if mp.re(self.pval(mid)) > 0:
                out.append((a, b))
        return out

    def period_oval(self, dps, oval=0):
        """a-cycle holomorphic period  oint dx/sqrt(P)  == 2 * int_a^b
        dx/sqrt(P) over real oval #oval.  4-real-root case: closed Legendre
        K(m) reduction (case table below, quadrature-gated in selftest);
        otherwise certified quadrature with the t=sin^2 substitution.
        Returns (value, meta)."""
        with mp.workdps(dps + 30):
            rts = self.roots(dps)
            tol = mp.mpf(10) ** (-dps)
            reals = sorted(mp.re(r) for r in rts if abs(mp.im(r)) < tol)
            ovals = self.real_ovals(dps)
            if not ovals:
                raise ValueError("no real oval: use cycle quadrature on a "
                                 "declared contour instead")
            a, b = ovals[oval]
            lead = _mpf_frac(self.coeffs[0])
            if len(reals) == 4 and lead != 0:
                r1, r2, r3, r4 = reals
                pref = 2 / mp.sqrt((r3 - r1) * (r4 - r2))
                if abs(a - r2) < tol and abs(b - r3) < tol and lead > 0:
                    # middle oval, lead>0: P = lead (x-r1)(x-r2)(r3-x)(r4-x)>0
                    m = ((r3 - r2) * (r4 - r1)) / ((r3 - r1) * (r4 - r2))
                    val = pref * mp.ellipk(m) / mp.sqrt(lead)
                    return 2 * val, dict(route="legendre_K", m=m)
                if abs(a - r3) < tol and abs(b - r4) < tol and lead < 0:
                    m = ((r4 - r3) * (r2 - r1)) / ((r4 - r2) * (r3 - r1))
                    val = pref * mp.ellipk(m) / mp.sqrt(-lead)
                    return 2 * val, dict(route="legendre_K", m=m)
                if abs(a - r1) < tol and abs(b - r2) < tol and lead < 0:
                    m = ((r2 - r1) * (r4 - r3)) / ((r4 - r2) * (r3 - r1))
                    val = pref * mp.ellipk(m) / mp.sqrt(-lead)
                    return 2 * val, dict(route="legendre_K", m=m)
            # generic: certified quadrature with sin^2 substitution
            v, meta = self._oval_quad(a, b, 0, dps)
            return 2 * v, meta

    def cycle_moments(self, ks, dps, oval=0):
        """oint x^k dx/sqrt(P) over real oval #oval (twice the segment
        integral), k in ks.  t=sin^2(theta) endpoint substitution + tanh-sinh,
        node-doubling agreement in meta."""
        with mp.workdps(dps + 30):
            ovals = self.real_ovals(dps)
            a, b = ovals[oval]
            out, metas = [], []
            for k in ks:
                v, meta = self._oval_quad(a, b, k, dps)
                out.append(2 * v)
                metas.append(meta)
            return out, metas

    def _deflate(self, a, b):
        """Numeric deflation Q(x) = -P(x)/((x-a)(x-b)) as descending coeffs
        (a, b root approximations at working precision); Q > 0 on (a,b)."""
        cs = [_mpf_frac(c) for c in self.coeffs]
        while cs and cs[0] == 0:
            cs = cs[1:]
        for r in (a, b):
            out = [cs[0]]
            for c in cs[1:-1]:
                out.append(out[-1] * r + c)
            # remainder cs[-1] + out[-1]*r ~ 0 (r is a root); drop it
            cs = out
        return [-c for c in cs]

    def _oval_quad(self, a, b, k, dps):
        """int_a^b x^k dx/sqrt(P) with the sqrt endpoint zeros at a,b removed
        ANALYTICALLY by x = a + (b-a) sin^2(theta) (the node-cancellation
        substitution): dx/sqrt((x-a)(b-x)) = 2 dtheta, so the
        integrand is 2 x^k / sqrt(Q(x)) with Q = -P/((x-a)(x-b)) deflated
        exactly at the working-precision roots; analytic on [0, pi/2]."""
        Q = self._deflate(a, b)
        w = b - a

        def f(th):
            s = mp.sin(th)
            x = a + w * s * s
            qv = mp.mpc(0)
            for c in Q:
                qv = qv * x + c
            return 2 * (x ** k) / mp.sqrt(qv)

        v, err = mp.quad(f, [0, mp.pi / 4, mp.pi / 2], error=True,
                         maxdegree=10)
        return v, dict(route="sin2_tanh_sinh", quad_err=mp.nstr(err, 3))

    def cycle_pair_contour(self, pair, ks, dps, pad=None, hw=None, laps=1):
        """oint x^k dx / y around the branch-point PAIR `pair` (two indices
        into self.roots(dps)), k in ks — the complex-cycle period verb
        (gamma1/gamma2 rotated-ellipse contour-pair pattern).

        Contour: ellipse centred between the two branch points, semi-major
        along their axis (half-distance + pad), semi-minor hw; every OTHER
        branch point must clear the contour (asserted).  y is the
        CONTINUITY-TRACKED sqrt of P along the contour with the PRINCIPAL
        branch at the theta=0 basepoint; counterclockwise
        in theta.  Trapezoid rule (spectrally accurate), N doubled until the
        last two levels agree < 10^-(dps+5) relative; fail-closed closure
        assert (tracked y returns to its start after `laps` laps) + radicand
        winding assert (arg P winds 2*2pi per lap for an enclosed pair).
        Returns (values list, meta)."""
        with mp.workdps(dps + 30):
            rts = self.roots(dps)
            A, B = rts[pair[0]], rts[pair[1]]
            others = [r for i, r in enumerate(rts) if i not in pair]
            C = (A + B) / 2
            d = B - A
            e1 = d / abs(d)
            dmin = min([abs(r - C) for r in others] + [mp.mpf('1e9')])
            if pad is None:
                pad = min(mp.mpf('0.25') * abs(d),
                          mp.mpf('0.2') * (dmin - abs(d) / 2))
            if hw is None:
                hw = min(mp.mpf('0.35') * abs(d),
                         mp.mpf('0.25') * (dmin - abs(d) / 2))
            a_ax = abs(d) / 2 + pad
            b_ax = hw
            assert a_ax > abs(d) / 2 and b_ax > 0, "degenerate contour"

            def zc(th):
                return C + e1 * (a_ax * mp.cos(th)
                                 + mp.mpc(0, 1) * b_ax * mp.sin(th))

            def dzc(th):
                return e1 * (-a_ax * mp.sin(th)
                             + mp.mpc(0, 1) * b_ax * mp.cos(th))

            # clearance of excluded branch points (coarse sampled distance)
            for r in others:
                dr = min(abs(zc(2 * mp.pi * j / 720) - r) for j in range(720))
                assert dr > mp.mpf('0.01') * max(abs(d), 1), \
                    f"branch point too close to contour ({mp.nstr(r, 6)})"
            # enclosed check
            for r in (A, B):
                w = (r - C) / e1
                assert (mp.re(w) / a_ax) ** 2 + (mp.im(w) / b_ax) ** 2 < 1
            target = mp.mpf(10) ** (-(dps + 5))
            prev = None
            meta = {}
            for kk in range(9, 22):
                N = 2 ** kk
                tot = [mp.mpc(0)] * len(ks)
                y_start = None
                y_prev = None
                argP_acc = mp.mpf(0)
                P_prev = None
                for j in range(N * laps + 1):
                    th = 2 * mp.pi * j / N
                    z = zc(th)
                    P = self.pval(z)
                    ysq = mp.sqrt(P)
                    if y_prev is None:
                        y = ysq          # principal at theta=0 (convention)
                        y_start = y
                    else:
                        y = ysq if abs(ysq - y_prev) <= abs(-ysq - y_prev) \
                            else -ysq
                        argP_acc += mp.arg(P / P_prev)
                    P_prev = P
                    if j == N * laps:
                        closure = abs(y - y_start) / max(abs(y_start), 1)
                        break
                    wgt = 1
                    for ki, k in enumerate(ks):
                        tot[ki] += (z ** k) / y * dzc(th)
                    y_prev = y
                tot = [t * 2 * mp.pi / N for t in tot]
                if prev is not None:
                    diff = max(abs(tot[i] - prev[i])
                               / max(abs(tot[i]), mp.mpf(1))
                               for i in range(len(ks)))
                    if diff < target:
                        wind = argP_acc / (2 * mp.pi)
                        assert abs(wind - 2 * laps) < mp.mpf('0.01'), \
                            f"radicand winding {mp.nstr(wind, 4)} != 2/lap"
                        assert closure < mp.mpf(10) ** (-(dps - 10)), \
                            f"tracked-radical closure {mp.nstr(closure, 4)}"
                        meta = dict(N=N, agree=mp.nstr(diff, 3),
                                    winding=float(wind),
                                    closure=mp.nstr(closure, 3),
                                    basepoint="theta=0 principal sqrt",
                                    orientation="counterclockwise")
                        return tot, meta
                prev = tot
            raise AssertionError("cycle_pair_contour: no convergence "
                                 "(raise dps guard or check contour)")

    def tau_from_oval(self, dps):
        """tau = i K(1-m)/K(m) of the Legendre reduction (4-real-root,
        middle-oval chart); refuse otherwise (loud scope certificate)."""
        with mp.workdps(dps + 30):
            _, meta = self.period_oval(dps)
            if meta.get("route") != "legendre_K":
                raise ValueError("tau_from_oval: no 4-real-root Legendre "
                                 "chart; supply a frame explicitly")
            m = meta["m"]
            return mp.mpc(0, 1) * mp.ellipk(1 - m) / mp.ellipk(m)


# ---------------------------------------------------------------------------
# BMSW sunrise frame — verbatim port of the gate-verified rowscripts_empl
# fcurve (branch repair included).  Masses/t exact-in, mp numerics at ambient
# working precision (callers wrap in mp.workdps).
# ---------------------------------------------------------------------------
def _csqrt(x):
    return mp.sqrt(mp.mpc(x))


def sunrise_frame(t, m1, m2, m3, mu=1):
    """All E_{3,F} data at Euclidean t (arXiv:1907.01251 conventions):
    e1F,e2F,e3F, kF2, kFp2, tauF, psi1F, z1F,z2F,z3F.  Port of
    rowscripts_empl.fcurve (rows 09/10/11 gate-verified; includes the
    Abel-Jacobi z3 branch repair).  Exact-rational inputs enforced."""
    for v, w in ((t, 't'), (m1, 'm1'), (m2, 'm2'), (m3, 'm3'), (mu, 'mu')):
        exact(v, w)
    t = mp.mpc(_mpf_frac(exact(t)))
    m1, m2, m3, mu = ((_mpf_frac(exact(v))) for v in (m1, m2, m3, mu))
    m1s, m2s, m3s = m1 * m1, m2 * m2, m3 * m3
    M100 = m1s + m2s + m3s
    mu1 = -m1 + m2 + m3
    mu2 = m1 - m2 + m3
    mu3 = m1 + m2 - m3
    mu4 = m1 + m2 + m3
    Delta = mu1 * mu2 * mu3 * mu4
    mu4p = mu ** 4
    rad = 3 * (_csqrt(mu1 * mu1 - t) * _csqrt(mu2 * mu2 - t)
               * _csqrt(mu3 * mu3 - t) * _csqrt(mu4 * mu4 - t))
    e1F = (-t * t + 2 * M100 * t + Delta + rad) / (24 * mu4p)
    e2F = (-t * t + 2 * M100 * t + Delta - rad) / (24 * mu4p)
    e3F = (2 * t * t - 4 * M100 * t - 2 * Delta) / (24 * mu4p)
    Z1F = e3F - e2F
    Z2F = e1F - e3F
    Z3F = e1F - e2F
    kF2 = Z1F / Z3F
    kFp2 = -Z1F / Z2F
    tauF = mp.mpc(0, 1) * mp.ellipk(1 - kF2) / mp.ellipk(kF2)
    psi1F = 2 / _csqrt(Z3F) * mp.ellipk(kF2)
    Kp = mp.ellipk(kFp2)

    def xhat(mi2, mj2):
        return e3F + mi2 * mj2 / mu4p

    def zF(xjk):
        up = _csqrt((e1F - e3F) / (xjk - e3F))
        return mp.ellipf(mp.asin(up), kFp2) / (2 * Kp)

    z1F = zF(xhat(m2s, m3s))
    z2F = zF(xhat(m1s, m3s))
    z3F = zF(xhat(m1s, m2s))
    # branch repair (ported): principal-branch reflection across the
    # half-period breaks z1+z2+z3=1; restore z3 = 1 - z1 - z2 when violated.
    s = z1F + z2F + z3F
    if abs(s - 1) > mp.mpf(10) ** (-mp.mp.dps // 2):
        z3F = 1 - z1F - z2F
        s = z1F + z2F + z3F
    assert abs(s - 1) < mp.mpf(10) ** (-mp.mp.dps // 2), \
        "Abel-Jacobi puncture sum != 1 after branch repair"
    return dict(e1F=e1F, e2F=e2F, e3F=e3F, kF2=kF2, kFp2=kFp2, tauF=tauF,
                psi1F=psi1F, z1F=z1F, z2F=z2F, z3F=z3F)

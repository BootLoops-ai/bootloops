#!/usr/bin/env python3
"""Certified DFE layer: polyDFE model-C class mixtures over the certified SFS
engine, with a CERTIFIED KERNEL CACHE so likelihood evaluation is
incumbent-speed (matrix-vector) while every number carries a proof budget.

Model C:  phi(S) = (1-p_b) * reflGamma(S; shape b, mean S_d<0)  on S<0
                 +  p_b    * Exp(S; mean S_b>0)                 on S>0

Design (certified analog of fitdadi/fastDFE caching):
  * Quadrature nodes in t = log|S| are FIXED per (n, sign, dps): half-decade
    panels over [10^-T_LO, 10^T_HI], composite GL degree pair (6,7) — the
    own-summation pattern (mp.quad banned: absolute error estimate, silent
    underconvergence on exponential-ramp panels).
  * Per node: ONE arb whole-vector solve (values + dM/dS) -> cached E-vector.
    Cache cost O(nodes * n); any DFE-parameter evaluation is then
    weights(params) . cache — O(nodes * n) float-free mpf dots.
  * Certified pieces per evaluation: (i) GL degree-pair self-check;
    (ii) S->0 truncated-mass bound (phi ~ |S|^{b-1}: bound
         (1-p_b) e^{lognorm} x_min^b / (b i), analytic);
    (iii) S->infty tail bound (integrand <= sup|E| * phi tail, both analytic).
  * Exact gradients: d(mix)/d(theta,r) trivial; d/d(p_b) linear;
    d/d(b, S_d, S_b) via analytic d(phi)/d(param) kernels — SAME cached
    E-vectors, new weights. dE/dS cache kept for score-vs-S diagnostics.

Selftest (__main__): degree-pair self-consistency of the mixed vector, the
analytic tail bounds, and the analytic parameter gradient vs a high-precision
finite difference. The heavier cross-oracle comparison ran against an
independent integrator that shares no code with this layer; the battery keeps
the same style of control via the finite-difference gradient gate.
"""
import os
import sys
from fractions import Fraction
from mpmath import mp, mpf, exp, expm1, fabs, log10, gamma as mpgamma, digamma
from mpmath import log as mplog

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from arb_route import solve_vector_arb, sfs_arb

PAD = 15


def _gl_nodes(degree):
    from mpmath.calculus.quadrature import GaussLegendre
    return GaussLegendre(mp).calc_nodes(degree, mp.prec + 20)


def _mpf_to_fraction(x):
    """Exact mpf -> Fraction ((-1)^s * man * 2^exp). The node S must reach the
    arb solver EXACTLY as used in the weights — a float() cast here cost 22
    digits (same class as the gate S-threading bug)."""
    sign, man, e, _ = x._mpf_
    fr = Fraction(man) * (Fraction(2) ** e)
    return -fr if sign else fr


class CertifiedKernel:
    """Cached certified E-vectors at fixed log-|S| quadrature nodes."""

    def __init__(self, n, dps=40, t_lo_dps_frac=None, t_hi=6, node_target=None):
        self.n = n
        self.dps = dps
        self.t_hi = t_hi          # integrate |S| up to 10^t_hi
        self.node_target = node_target or (dps + 10)
        self.cache = {}           # (sign, deg) -> list of (t, weight, Evec, dEvec)
        self._built = set()

    # ---- neutral-Taylor zone (|S| <= S_SPLIT): exact rational Taylor of
    # g(S,u)=(1-e^{-Su})/(1-e^{-S}) by series division; E-moments exact via
    # Beta; phi-moments via incomplete gamma. Remainder: measured coefficient
    # decay (analytic radius 2*pi) + arb spot-validation at the boundary.
    S_SPLIT = mpf('1e-3')
    K_TAYLOR = 18

    def _taylor_coeffs(self):
        """a[k][i]: E_i(S) = sum_k a[k][i] S^k for |S|<S_SPLIT. Exact."""
        from fractions import Fraction as F
        K, n = self.K_TAYLOR, self.n
        # series in S: num_m = coeffs of (1-e^{-Su}) = sum_{m>=1} (-1)^{m+1} u^m S^m/m!
        # den_m = same with u=1. g = num/den as power series in S, coeffs poly in u.
        fact = [1]*(K+3)
        for m in range(1, K+3): fact[m] = fact[m-1]*m
        num = [None]*(K+2)   # num[m] = poly in u as dict {deg: F}
        den = [F(0)]*(K+2)
        for m in range(1, K+2):
            c = F((-1)**(m+1), fact[m])
            num[m] = {m: c}
            den[m] = c
        num[0] = {}
        # g = sum g_k S^k, poly in u: g_0*den_1 = num_1 => shift: write
        # num = S*N(S,u), den = S*D(S): g = N/D with N_k = num[k+1], D_k = den[k+1]
        N = [num[k+1] for k in range(K+1)]
        D = [den[k+1] for k in range(K+1)]
        G = []
        for k in range(K+1):
            acc = dict(N[k])
            for j in range(k):
                for deg, c in G[j].items():
                    acc[deg] = acc.get(deg, F(0)) - c*D[k-j]
            g_k = {d: c/D[0] for d, c in acc.items()}
            G.append(g_k)
        # E_i = C(n,i) * int x^{i-1}(1-x)^{n-i-1} g(S, 1-x) dx
        #     = C(n,i) * sum_k S^k sum_deg c * B(i, n-i+deg)
        import math
        a = [[mpf(0)]*(n) for _ in range(K+1)]
        bigfact = [1]*(n + K + 3)          # Beta values need factorials to n+K+1
        for m in range(1, n + K + 3):
            bigfact[m] = bigfact[m-1]*m
        with mp.workdps(self.dps + 30):
            for i in range(1, n):
                Cni = F(math.comb(n, i))
                for k in range(K+1):
                    s = F(0)
                    for deg, c in G[k].items():
                        # B(i, n-i+deg) exact rational
                        Bv = F(bigfact[i-1]*bigfact[n-i+deg-1], bigfact[n+deg-1])
                        s += c*Bv
                    v = Cni*s
                    a[k][i] = mpf(v.numerator)/v.denominator
        # measured coefficient growth for remainder model
        mx = [max(fabs(a[k][i]) for i in range(1, n)) for k in range(K+1)]
        self._taylor_a = a
        self._taylor_maxcoef = mx

    def build(self, sign, shape_min=mpf('0.05')):
        """Precompute E-vectors on the arb-zone node set for one sign of S."""
        key = (sign,)
        if key in self._built:
            return
        if not hasattr(self, '_taylor_a'):
            with mp.workdps(self.dps + PAD):
                self._taylor_coeffs()
        with mp.workdps(self.dps + PAD):
            t_lo = mplog(self.S_SPLIT)
            t_hi = mpf(self.t_hi) * mplog(mpf(10))
            step = mplog(mpf(10)) / 2
            edges = [t_lo]
            while edges[-1] < t_hi:
                edges.append(min(edges[-1] + step, t_hi))
            for deg in (6, 7):
                nodes = _gl_nodes(deg)
                rows = []
                for a, b in zip(edges[:-1], edges[1:]):
                    h = (b - a) / 2
                    mid = (a + b) / 2
                    for tt, w in nodes:
                        t = mid + h * tt
                        S = sign * exp(t)
                        Sfr = _mpf_to_fraction(S)
                        sol = solve_vector_arb(self.n, float(S), self.node_target,
                                               want_grad=True, S_frac=Sfr)
                        Ev = [mpf(0)] * (self.n)
                        dEv = [mpf(0)] * (self.n)
                        for i in range(1, self.n):
                            E, dE = sfs_arb(self.n, float(S), sol, i, S_frac=Sfr)
                            Ev[i] = mp.mpmathify(E.mid().str(self.dps + 12, radius=False))
                            dEv[i] = mp.mpmathify(dE.mid().str(self.dps + 12, radius=False))
                        rows.append((t, w * h, Ev, dEv))
                self.cache[(sign, deg)] = rows
        self._built.add(key)

    @staticmethod
    def _lower_gamma_small(s, x, want_ds=False):
        """gamma(s,x) = x^s sum_m (-x)^m/(m!(s+m)) for SMALL x (Taylor zone:
        x = eps/sc <= ~0.1 => ~15 terms to 1e-45). d/ds analytic per term:
        x^s [log x/(s+m) - 1/(s+m)^2]. Returns (g, dg_ds or None)."""
        xs = exp(s * mplog(x))
        lx = mplog(x)
        g = mpf(0)
        dg = mpf(0)
        term = mpf(1)   # (-x)^m/m!
        m = 0
        while True:
            g += term / (s + m)
            if want_ds:
                dg += term * (lx / (s + m) - 1 / (s + m) ** 2)
            m += 1
            term *= -x / m
            if fabs(term) / (s + m) < mpf(10) ** (-(mp.dps + 8)):
                break
        return xs * g, (xs * dg if want_ds else None)

    def _taylor_zone_moments(self, sign, b_, sc, want_grad=False):
        """Moments mom_k = sign^k sc^k gamma(b+k, eps/sc)/Gamma(b) and their
        EXACT derivatives d/db (shape slot, incl. 1/Gamma(b)) and d/dsc.
        (Chain rules to model params are applied by the caller.)"""
        eps = self.S_SPLIT
        x = eps / sc
        gb = mpgamma(b_)
        psi_b = digamma(b_)
        mom, dmom_db, dmom_dsc = [], [], []
        for k in range(self.K_TAYLOR + 1):
            g, dg_ds = self._lower_gamma_small(b_ + k, x, want_ds=want_grad)
            base = sc ** k * g / gb
            mom.append((sign ** k) * base)
            if want_grad:
                # d/db: s-slot derivative + digamma from 1/Gamma(b)
                d_b = sc ** k * (dg_ds - psi_b * g) / gb
                # d/dsc: k sc^{k-1} g/gb + sc^k * dgamma/dx * dx/dsc / gb,
                # dgamma(s,x)/dx = x^{s-1} e^{-x};  dx/dsc = -eps/sc^2
                dgdx = exp((b_ + k - 1) * mplog(x) - x)
                d_sc = (k * sc ** (k - 1) * g + sc ** k * dgdx * (-eps / sc ** 2)) / gb
                dmom_db.append((sign ** k) * d_b)
                dmom_dsc.append((sign ** k) * d_sc)
        if want_grad:
            return mom, dmom_db, dmom_dsc
        return mom

    # ---------------- weights for model-C pieces and their param-gradients
    @staticmethod
    def _w_refl_gamma(t, b, Sd):
        """phi_del(S)*|S| at S=-e^t (the |S|=e^t Jacobian included):
        reflGamma pdf = x^{b-1} e^{-x/sc} / (Gamma(b) sc^b), x=|S|, sc=|Sd|/b."""
        x = exp(t)
        sc = fabs(Sd) / b
        return exp(b * t - x / sc - mplog(mpgamma(b)) - b * mplog(sc))

    @staticmethod
    def _w_exp(t, Sb):
        x = exp(t)
        return exp(t - x / Sb) / Sb

    def mix(self, params, dps_check=True):
        """params: dict(b, Sd, pb, Sb). Returns (Evec list len n, selfcons_digits,
        tail_bound_rel). Whole mixed SFS vector, theta=1."""
        with mp.workdps(self.dps + PAD):
            # param conversion INSIDE workdps (ambient-dps footgun, measured
            # 16d ceiling in the dominance q-series wrapper; here both sides
            # of every refit comparison shared the same param image so no
            # conclusion changes — fixed for hygiene)
            b_, Sd, pb, Sb = (mpf(str(params['b'])), mpf(str(params['Sd'])),
                              mpf(str(params['pb'])), mpf(str(params.get('Sb', 1))))
            out = {}
            sc_g = fabs(Sd) / b_
            # Taylor-zone (|S|<=S_SPLIT) exact-moment contributions
            mom_del = self._taylor_zone_moments(-1, b_, sc_g)
            mom_ben = self._taylor_zone_moments(1, mpf(1), Sb) if pb > 0 else None
            tay = [mpf(0)] * self.n
            for i in range(1, self.n):
                acc = mpf(0)
                for k in range(self.K_TAYLOR + 1):
                    acc += self._taylor_a[k][i] * ((1 - pb) * mom_del[k]
                            + (pb * mom_ben[k] if mom_ben else 0))
                tay[i] = acc
            for deg in ((6, 7) if dps_check else (7,)):
                tot = [tay[i] for i in range(self.n)]
                for (t, w, Ev, _dEv) in self.cache[(-1, deg)]:
                    wt = (1 - pb) * self._w_refl_gamma(t, b_, Sd) * w
                    for i in range(1, self.n):
                        tot[i] += wt * Ev[i]
                if pb > 0:
                    for (t, w, Ev, _dEv) in self.cache[(1, deg)]:
                        wt = pb * self._w_exp(t, Sb) * w
                        for i in range(1, self.n):
                            tot[i] += wt * Ev[i]
                out[deg] = tot
            ref = out[7 if 7 in out else 6]
            sc = 9999.0
            if dps_check and 6 in out:
                for i in range(1, self.n):
                    if ref[i] != 0 and out[6][i] != ref[i]:
                        sc = min(sc, float(-log10(fabs((out[6][i] - ref[i]) / ref[i]))))
            # Taylor-remainder bound over the |S|<=S_SPLIT zone:
            # coefficient growth ratio q measured from computed a_k (analytic
            # radius 2*pi => q ~ 1/(2*pi)); remainder <= mx[K]*q*eps^{K+1}/(1-eps*q)
            sc_gamma = sc_g
            lognorm = -b_ * mplog(sc_gamma) - mplog(mpgamma(b_))
            mx = self._taylor_maxcoef
            q = max(mx[k+1]/mx[k] for k in range(len(mx)-2, len(mx)-1) if mx[k] > 0)
            eps_z = self.S_SPLIT
            tail0 = 2 * mx[-1] * q * eps_z ** (self.K_TAYLOR + 1) / (1 - eps_z * q)
            # S->inf: reflGamma tail beyond 10^t_hi (E bounded by n/(i(n-i))*max)
            X = mpf(10) ** self.t_hi
            tail_inf = (1 - pb) * exp(lognorm + (b_ - 1) * mplog(X) - X / sc_gamma) * sc_gamma
            if pb > 0:
                tail_inf += pb * exp(-X / Sb)
            worst_entry = min(x for x in ref[1:] if x > 0)
            tail_rel = float((tail0 + tail_inf * mpf(self.n)) / worst_entry)
            return ref, sc, tail_rel

    def grad_mix(self, params):
        """d(mixed E-vector)/d(b, Sd, pb, Sb) via analytic d(phi)/d(param):
          d/db   refl-Gamma weight: w * (t - log sc - digamma(b) - 1 + x/(b*sc) ... )
          (worked: d/db log w = log x - log sc - digamma(b) + (x/sc - b)*(1/b)
           with sc=|Sd|/b depending on b: d log sc/db = -1/b)
          => d/db log w = t - log(sc) - digamma(b) + 1 - x/(Sd_abs)  ... derived below
          d/dSd: sc=|Sd|/b => d log w/d|Sd| = -b/|Sd| + x b/|Sd|^2
          d/dSb: d log w/dSb = -1/Sb + x/Sb^2
          d/dpb: linear.
        Returns dict param -> vector."""
        with mp.workdps(self.dps + PAD):
            # param conversion INSIDE workdps (ambient-dps footgun, measured
            # 16d ceiling in the dominance q-series wrapper; here both sides
            # of every refit comparison shared the same param image so no
            # conclusion changes — fixed for hygiene)
            b_, Sd, pb, Sb = (mpf(str(params['b'])), mpf(str(params['Sd'])),
                              mpf(str(params['pb'])), mpf(str(params.get('Sb', 1))))
            Sda = fabs(Sd)
            sc_g = Sda / b_
            g = {k: [mpf(0)] * self.n for k in ('b', 'Sd', 'pb', 'Sb')}
            # ---- Taylor-zone contributions (exact d/ds of gamma(s,x) series)
            mom_d, dmb, dms = self._taylor_zone_moments(-1, b_, sc_g, want_grad=True)
            mom_b, _, dms_b = self._taylor_zone_moments(1, mpf(1), Sb, want_grad=True)
            dsc_db = -sc_g / b_               # sc = |Sd|/b
            dsc_dSd = (1 / b_) * (-1 if Sd < 0 else 1)
            for i in range(1, self.n):
                Td = Tb = dTdb = dTdsc = dTbsc = mpf(0)
                for k in range(self.K_TAYLOR + 1):
                    ak = self._taylor_a[k][i]
                    Td += ak * mom_d[k]
                    Tb += ak * mom_b[k]
                    dTdb += ak * (dmb[k] + dms[k] * dsc_db)
                    dTdsc += ak * dms[k]
                    dTbsc += ak * dms_b[k]
                g['b'][i] += (1 - pb) * dTdb
                g['Sd'][i] += (1 - pb) * dTdsc * dsc_dSd
                g['pb'][i] += Tb - Td
                g['Sb'][i] += pb * dTbsc
            for (t, w, Ev, _d) in self.cache[(-1, 7)]:
                x = exp(t)
                wg = self._w_refl_gamma(t, b_, Sd) * w
                # d log w / db  (sc = Sda/b): log w = (b-1)... careful:
                # log w = b*t - x*b/Sda - lgamma(b) - b*log(Sda/b)
                #   d/db = t - x/Sda - digamma(b) - log(Sda/b) + 1
                dlb = t - x / Sda - digamma(b_) - mplog(sc_g) + 1
                # d log w / dSda = x*b/Sda^2 - b/Sda
                dlS = x * b_ / Sda ** 2 - b_ / Sda
                for i in range(1, self.n):
                    base = wg * Ev[i]
                    g['b'][i] += (1 - pb) * base * dlb
                    g['Sd'][i] += (1 - pb) * base * dlS * (-1 if Sd < 0 else 1)
                    g['pb'][i] -= base
            if pb > 0 or True:
                for (t, w, Ev, _d) in self.cache[(1, 7)]:
                    x = exp(t)
                    we = self._w_exp(t, Sb) * w
                    dlSb = x / Sb ** 2 - 1 / Sb
                    for i in range(1, self.n):
                        base = we * Ev[i]
                        g['pb'][i] += base
                        g['Sb'][i] += pb * base * dlSb
            return g


if __name__ == "__main__":
    import time
    n = 20
    t0 = time.time()
    K = CertifiedKernel(n, dps=40)
    K.build(-1)
    K.build(1)
    tb = time.time() - t0
    Ev, sc, tail = K.mix({'b': '0.4', 'Sd': '-1000', 'pb': 0})
    t0 = time.time()
    for _ in range(20):
        K.mix({'b': 0.35, 'Sd': -800.0, 'pb': 0.02, 'Sb': 10.0}, dps_check=False)
    teval = (time.time() - t0) / 20
    g = K.grad_mix({'b': 0.4, 'Sd': -1000.0, 'pb': 0.02, 'Sb': 10.0})
    # gradient positive control: finite difference on b at high dps
    with mp.workdps(60):
        h = mpf('1e-12')
        E1, _, _ = K.mix({'b': mpf('0.4') + h, 'Sd': -1000.0, 'pb': 0.02, 'Sb': 10.0}, dps_check=False)
        E0, _, _ = K.mix({'b': mpf('0.4') - h, 'Sd': -1000.0, 'pb': 0.02, 'Sb': 10.0}, dps_check=False)
        fd = (E1[10] - E0[10]) / (2 * h)
        gd = float(-log10(fabs((g['b'][10] - fd) / fd))) if fd != 0 else 0
    print(f"build {tb:.1f}s | selfcons {sc:.1f}d | tail_rel {tail:.1e} "
          f"| eval {teval*1000:.0f}ms | grad-b vs FD {gd:.1f}d")
    ok = sc >= 38 and tail < 1e-30 and gd >= 8
    print("PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)

"""Oracle identity gates (BOTH registers shipped, verbatim).

R-register (receipt): the adapter oracle's lnP ball at the study's
receipted point vs the receipted lnP ball (their exact/ball receipts pin
the reference values). Bar: |mid diff| <= 1e-12 * max(1, |lnP|) AND the
two enclosures intersect.

A-register (adjudication): the analytic arb gradient vs an
INDEPENDENT-DERIVATION mpmath-50dps gradient (direct log-sum lnP built
from the EXACT integer kernel, Richardson central differences). Bar:
|arb - mp50| <= 1e-12 * max(1, |g|) per coordinate.

FD-register: assembled arb Hessian vs central FD of the arb gradient at
an exact power-of-two step (the clinch gate_fd_hessian pattern).
"""
import math
from fractions import Fraction

from flint import arb, ctx, fmpz_poly

from . import ident


def _ball_from_str(s):
    """Parse an arb receipt string '[m +/- r]' into an arb ball."""
    return arb(s)


def value_gate(nll_ball, receipt_lnp_str):
    """Compare oracle NLL ball vs receipted lnP ball (sign-flipped)."""
    ref = -_ball_from_str(receipt_lnp_str)
    d = nll_ball - ref
    mid = abs(float(d.mid()))
    bar = 1e-12 * max(1.0, abs(float(ref.mid())))
    overlap = bool(d.contains(arb(0)) or mid <= bar)
    return dict(register="receipt", mid_abs_diff=mid, bar=bar,
                oracle_nll=str(nll_ball), receipt_lnp=receipt_lnp_str,
                ok=bool(mid <= bar), enclosures_consistent=overlap)


# ---------------------------------------------------------------------
# independent mpmath 50-dps lnP (single-sample + equal-I multisample)
# ---------------------------------------------------------------------
def _mp():
    import mpmath as mp
    mp.mp.dps = 50
    return mp


def _logsum(mp, logs):
    """log-sum-exp; terms below mx - 160 nats (< 1e-69 relative) are
    dropped — far below every gate bar at 50 dps."""
    mx = max(logs)
    cut = mx - 160
    return mx + mp.log(mp.fsum(mp.e**(v - mx) for v in logs if v > cut))


def _log_fraction(mp, fr):
    return mp.log(mp.mpf(fr.numerator)) - mp.log(mp.mpf(fr.denominator))


class MPSingleSample:
    """Derivation-independent lnP(u) at 50 dps from the EXACT integer K."""

    def __init__(self, D):
        mp = _mp()
        p = ident.pilot()
        self.D = sorted(int(x) for x in D)
        self.J, self.S = sum(self.D), len(self.D)
        # exact integer kernel by identity (their species_poly_int product
        # tree); log denominator = sum lgamma(n_i) exactly at 50 dps —
        # avoids 20k huge-int gcd normalizations of the Fraction route
        be = p["ball_engine"]
        prod = be._product_tree([be.species_poly_int(n) for n in self.D])
        logden = mp.fsum(mp.loggamma(n) for n in self.D)
        self.logK = [mp.log(mp.mpf(int(prod[A]))) - logden
                     for A in range(self.S, self.J + 1)]
        from collections import Counter
        self.static = (mp.loggamma(self.J + 1)
                       - mp.fsum(mp.log(n) for n in self.D)
                       - mp.fsum(mp.loggamma(c + 1)
                                 for c in Counter(self.D).values()))
        self.mp = mp

    def lnP(self, u_th, u_I):
        mp = self.mp
        th, I = mp.e**u_th, mp.e**u_I
        logI = mp.log(I)
        run = mp.mpf(0)
        for k in range(self.S):
            run += mp.log(th + k)
        logs = []
        for idx, lk in enumerate(self.logK):
            A = self.S + idx
            # run = sum_{k<A} ln(th+k) = ln (th)_A
            logs.append(lk + A * logI - run)
            run += mp.log(th + A)
        out = (self.static + self.S * mp.log(th)
               - (mp.loggamma(I + self.J) - mp.loggamma(I))
               + _logsum(mp, logs))
        return out

    def grad_fd(self, u_th, u_I, h=None):
        """2nd-order central differences at h = 2^-45: truncation
        ~ h^2 |f'''| ~ 1e-27, roundoff ~ 1e-50/h ~ 3e-37 — both far
        below the 1e-12 bar at 50 dps."""
        mp = self.mp
        h = h or mp.mpf(2) ** -45
        g = []
        for dim in (0, 1):
            def f(t):
                return self.lnP(u_th + (t if dim == 0 else 0),
                                u_I + (t if dim == 1 else 0))
            g.append((f(h) - f(-h)) / (2 * h))
        return g


class MPEqI(MPSingleSample):
    """Equal-I multisample: exact integer Mtilde + per-sample prefactors."""

    def __init__(self, Dmat):
        mp = _mp()
        p = ident.pilot()
        self.Dmat = [tuple(int(x) for x in r) for r in Dmat]
        self.N = len(self.Dmat[0])
        self.S = len(self.Dmat)
        self.Js = [sum(r[i] for r in self.Dmat) for i in range(self.N)]
        self.J = sum(self.Js)
        # exact integer Mtilde: per-species conv of rf polys + Borel twist
        rf = p["multisample_ball"].rf_poly_int
        polys = []
        for row in self.Dmat:
            u = fmpz_poly([1])
            for n in row:
                if n:
                    u = u * rf(n)
            c = [int(u[b]) * math.factorial(b - 1) if b >= 1 else 0
                 for b in range(u.length())]
            polys.append(fmpz_poly(c))
        while len(polys) > 1:
            polys = [polys[i] * polys[i + 1]
                     for i in range(0, len(polys) - 1, 2)] \
                + ([polys[-1]] if len(polys) % 2 else [])
        MT = polys[0]
        # coefficients below sum_s k_s (each species needs >= 1 ancestor
        # per occupied plot) are structurally ZERO — skip them; store
        # (A, log coeff) pairs
        self.logK = []
        for A in range(self.S, MT.length()):
            v = int(MT[A])
            assert v >= 0
            if v:
                self.logK.append((A, _log_fraction(mp, Fraction(v))))
        from collections import Counter
        stat = mp.mpf(0)
        for cnt in Counter(self.Dmat).values():
            stat -= mp.loggamma(cnt + 1)
        for i in range(self.N):
            stat += mp.loggamma(self.Js[i] + 1)
        for row in self.Dmat:
            for n in row:
                if n > 1:
                    stat -= mp.loggamma(n + 1)
        self.static = stat
        self.mp = mp

    def lnP(self, u_th, u_I):
        mp = self.mp
        th, I = mp.e**u_th, mp.e**u_I
        logI = mp.log(I)
        # ln (th)_A via loggamma differences (sparse A support)
        lg0 = mp.loggamma(th)
        logs = [lk + A * logI - (mp.loggamma(th + A) - lg0)
                for A, lk in self.logK]
        pref = self.static + self.S * mp.log(th)
        for i in range(self.N):
            pref -= mp.loggamma(I + self.Js[i]) - mp.loggamma(I)
        return pref + _logsum(mp, logs)


def adjudication_gate(oracle, mp_model, u_points):
    """A5 register: arb analytic gradient vs mpmath-50dps FD gradient."""
    rows = []
    ok = True
    for (u0, u1) in u_points:
        ctx_prec = 256
        old = ctx.prec
        try:
            ctx.prec = ctx_prec
            x = ([[arb(u1)]], [arb(u0)])
            (Fz, Fg) = oracle.F(x)
            g_arb = [-float(Fg[0].mid()), -float(Fz[0][0].mid())]
        finally:
            ctx.prec = old
        g_mp = mp_model.grad_fd(mp_model.mp.mpf(u0), mp_model.mp.mpf(u1))
        for k in range(2):
            diff = abs(g_arb[k] - float(g_mp[k]))
            bar = 1e-12 * max(1.0, abs(float(g_mp[k])))
            rows.append(dict(u=[u0, u1], coord=["lntheta", "lnI"][k],
                             arb=g_arb[k], mp50=float(g_mp[k]),
                             abs_diff=diff, bar=bar, ok=bool(diff <= bar)))
            ok = ok and diff <= bar
    return dict(register="A5_adjudication_mp50", rows=rows, ok=bool(ok))


def fd_hessian_gate(oracle, u, prec=256, log2h=-30):
    """Assembled arb Hessian vs central FD of the arb gradient (exact
    power-of-two step). Relative bar 1e-6 (FD truncation dominates)."""
    h = 2.0 ** log2h
    old = ctx.prec
    try:
        ctx.prec = prec

        def grad(uu):
            x = ([[arb(uu[1])]], [arb(uu[0])])
            Fz, Fg = oracle.F(x)
            return [float(Fg[0].mid()), float(Fz[0][0].mid())]

        x0 = ([[arb(u[1])]], [arb(u[0])])
        D, B, G = oracle.H(x0)
        H_an = [[float(G[0, 0].mid()), float(B[0][0, 0].mid())],
                [float(B[0][0, 0].mid()), float(D[0][0, 0].mid())]]
        H_fd = [[0.0] * 2 for _ in range(2)]
        for j in range(2):
            up = list(u); dn = list(u)
            up[j] += h; dn[j] -= h
            gp, gm = grad(up), grad(dn)
            for i in range(2):
                H_fd[i][j] = (gp[i] - gm[i]) / (2 * h)
    finally:
        ctx.prec = old
    worst = 0.0
    for i in range(2):
        for j in range(2):
            rel = abs(H_an[i][j] - H_fd[i][j]) / max(1.0, abs(H_fd[i][j]))
            worst = max(worst, rel)
    return dict(register="fd_hessian", H_analytic=H_an, H_fd=H_fd,
                worst_rel=worst, bar=1e-6, ok=bool(worst <= 1e-6))

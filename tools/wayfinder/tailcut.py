#!/usr/bin/env python3
r"""
tailcut.py — WP2: certified exponential tail bound for frozen-upper-cutoff
integrals (the x=50 class).

One shared helper — consumers import
besselk_tail_bound / choose_bessel_cutoff; nothing here is wired into any
row script by this module itself.

WHY: a frozen upper cutoff drops a dps-independent correlated tail that is
invisible to two-precision self-checks. Type specimen: the B-Bessel moment
B(eps) = 2^{d+1}/Gamma(d/2) * INT_0^inf x^{3-d} K_{d/2-1}(x)^4 dx, where a
frozen upper cutoff x=50 drops a tail of 2.7e-91 — a hidden ~91-digit
ceiling.
* Envelope (textbook estimate): for x >= max(1, |nu|^2),
      |K_nu(x)| <= sqrt(pi/(2x)) * e^{-x} * (1 + C_nu/x),
      C_nu = |4 nu^2 - 1| / 8   (explicit; |.| = complex modulus so the
      same expression serves the eps-fan's complex orders nu = 1 - eps).
  The validity domain x >= max(1, |nu|^2) is ASSERTED at runtime
  (TailDomainError, fail-closed); on request the envelope is additionally
  CHECKED numerically at the cutoff (verify=True evaluates |K_nu(X)| once
  against the envelope — an independent runtime control, not a proof).
* Tail bound: with m envelope factors and integrand
  prefactor x^sigma,
      |INT_X^inf x^sigma K_nu(x)^m dx|
        <= (pi/2)^{m/2} (1 + C_nu/X)^m * INT_X^inf x^a e^{-m x} dx,
      a = Re(sigma) - m/2 .
  The remaining integral is EXACT:
      a <= 0 :  <= X^a e^{-mX} / m                  (monotone x^a)
      a >  0 :  =  Gamma(a+1, mX) / m^{a+1}          (upper incomplete
                Gamma, mpmath.gammainc(a+1, m*X, mp.inf) — the spec's
                "Re(1-d) > 0" branch, closed form, still exact).
  Using Re(sigma) (|x^sigma| = x^{Re sigma} on x > 0) keeps the bound a
  genuine bound for complex sigma (the eps-fan case).
* Cutoff selection: X = f(dps) chosen BY THE BOUND at
  runtime — the current/legacy X is the SEED only; X grows until
  bound < tol = 10^-(dps+guard). The accepted bound VALUE is returned so
  callers propagate it into their printed certified error line. Growth is
  additive in the exponent (e^{-mX} decade steps), fail-closed at a cap
  (TailNonConvergence naming bound, tol, X, cap — the ctrl_raise message
  discipline of wiring-log items 12-14).
* eps-fan discipline: the bound must hold at the WORST node
  of the sampling circle — choose_bessel_cutoff takes the whole fan
  ([(nu, sigma), ...]) and maximizes the bound over nodes at each X.
  One shared X for every node keeps the per-node systematics identical.
* Import-dps footgun: all
  arithmetic runs inside mp.workprec(...); global mp.dps is NEVER touched.

WHAT IT DOES
============
besselk_tail_bound(nu, sigma, X, m=4, prefactor=1) -> mpf
    Certified bound on |prefactor * INT_X^inf x^sigma K_nu(x)^m dx|.
    Raises TailDomainError if X < max(1, |nu|^2) (validity domain).

choose_bessel_cutoff(dps, nodes, x_seed, m=4, guard=10, x_cap=100000,
                     verify=False) -> TailCert
    nodes = [(nu, sigma), ...] or [(nu, sigma, prefactor), ...] — the full
    eps-fan. Starts at X = max(x_seed, validity minimum over nodes) and
    increases X until the WORST-node bound < 10^-(dps+guard). Fail-closed
    TailNonConvergence at x_cap. TailCert carries X, bound, tol, the worst
    node index, and bound_line() for the printed certified error line.

FAILURE MODES (fail-closed, named)
==================================
* TailDomainError — envelope validity domain X >= max(1, |nu|^2) violated
  (only possible via besselk_tail_bound directly; choose_bessel_cutoff
  lifts the seed to the validity minimum and records it).
* TailNonConvergence — bound never beat tol by X = x_cap; message names
  the measured bound, tol, X, cap, worst node.
* verify=True envelope check failing raises TailEnvelopeError with the
  measured |K_nu(X)| and the envelope value.
"""

from mpmath import mp, mpf, mpc

__all__ = [
    "besselk_env_C", "besselk_tail_bound", "choose_bessel_cutoff",
    "TailCert", "TailDomainError", "TailNonConvergence",
    "TailEnvelopeError",
]


class TailDomainError(ValueError):
    """Envelope validity domain x >= max(1, |nu|^2) violated (fail-closed)."""


class TailEnvelopeError(RuntimeError):
    """Runtime envelope check |K_nu(X)| <= env(X) failed (fail-closed)."""


class TailNonConvergence(RuntimeError):
    """Fail-closed refusal: tail bound never beat tol by the X cap."""

    def __init__(self, *, bound, tol, X, cap, node_index):
        self.bound = bound
        self.tol = tol
        self.X = X
        self.cap = cap
        self.node_index = node_index
        with mp.workdps(6):
            msg = ("choose_bessel_cutoff FAILED (fail-closed): certified "
                   "tail bound %s >= tol %s at X=%s (cap %s), worst node %d"
                   % (mp.nstr(bound, 3), mp.nstr(tol, 3), mp.nstr(X, 6),
                      mp.nstr(mpf(cap), 6), node_index))
        super().__init__(msg)


class TailCert:
    """Certificate from choose_bessel_cutoff.

    X          — accepted cutoff (bound < tol there)
    bound      — the accepted WORST-node bound value (this is the number
                 callers add to their printed certified error line)
    tol        — 10^-(dps+guard)
    dps, guard, m — as run
    node_index — index (into nodes) of the worst node at acceptance
    x_seed     — the seed the caller passed (legacy/frozen X)
    x_start    — seed after the validity-domain lift
    verified   — True if the one-shot numeric envelope check ran and passed
    """

    __slots__ = ("X", "bound", "tol", "dps", "guard", "m", "node_index",
                 "x_seed", "x_start", "verified")

    def __init__(self, **kw):
        for k in self.__slots__:
            setattr(self, k, kw[k])

    def bound_line(self, label="tail"):
        with mp.workdps(6):
            return ("[certified] %s: cutoff X=%s (seed %s), worst-node "
                    "exp-tail bound %s <= tol %s (dps=%d guard=%d m=%d"
                    "%s)" % (label, mp.nstr(self.X, 6),
                             mp.nstr(mpf(self.x_seed), 6),
                             mp.nstr(self.bound, 3), mp.nstr(self.tol, 3),
                             self.dps, self.guard, self.m,
                             " env-checked" if self.verified else ""))

    def __repr__(self):
        with mp.workdps(6):
            return "TailCert(X=%s, bound=%s, tol=%s)" % (
                mp.nstr(self.X, 6), mp.nstr(self.bound, 3),
                mp.nstr(self.tol, 3))


def besselk_env_C(nu):
    """C_nu = |4 nu^2 - 1|/8 — the explicit envelope constant (spec)."""
    nu = mpc(nu)
    return abs(4 * nu * nu - 1) / 8


def _env_value(nu, x):
    """The envelope sqrt(pi/(2x)) e^{-x} (1 + C_nu/x) at x (validity NOT
    checked here — callers assert first)."""
    C = besselk_env_C(nu)
    return mp.sqrt(mp.pi / (2 * x)) * mp.e ** (-x) * (1 + C / x)


def _validity_min(nu):
    """Envelope validity threshold max(1, |nu|^2)."""
    nu = mpc(nu)
    return max(mpf(1), abs(nu) ** 2)


def besselk_tail_bound(nu, sigma, X, m=4, prefactor=1):
    """Certified bound on |prefactor * INT_X^inf x^sigma K_nu(x)^m dx|
    (module docstring has the derivation). Raises TailDomainError outside
    the envelope validity domain X >= max(1, |nu|^2)."""
    if m < 1 or int(m) != m:
        raise ValueError("besselk_tail_bound: m must be a positive integer "
                         "(got %r)" % (m,))
    X = mpf(X)
    vmin = _validity_min(nu)
    if not (X >= vmin):
        raise TailDomainError(
            "besselk_tail_bound: envelope validity domain violated: "
            "X=%s < max(1, |nu|^2)=%s (nu=%s) — raise the cutoff seed"
            % (mp.nstr(X, 6), mp.nstr(vmin, 6), mp.nstr(mpc(nu), 6)))
    C = besselk_env_C(nu)
    a = mp.re(mpc(sigma)) - mpf(m) / 2
    env_pref = (mp.pi / 2) ** (mpf(m) / 2) * (1 + C / X) ** m
    if a <= 0:
        core = X ** a * mp.e ** (-m * X) / m
    else:
        # exact: INT_X^inf x^a e^{-mx} dx = Gamma(a+1, mX) / m^(a+1)
        core = mp.gammainc(a + 1, m * X, mp.inf) / mpf(m) ** (a + 1)
    return abs(mpc(prefactor)) * env_pref * core


def choose_bessel_cutoff(dps, nodes, x_seed, m=4, guard=10, x_cap=100000,
                         verify=False):
    """Pick the cutoff X = f(dps) BY THE BOUND at runtime.

    dps    — target decimal digits; tol = 10^-(dps+guard)
    nodes  — [(nu, sigma), ...] or [(nu, sigma, prefactor), ...]: the FULL
             eps-fan; the bound is enforced at the WORST node
    x_seed — the current/legacy cutoff: the starting seed ONLY (it is
             lifted to the validity minimum if below it; it is never the
             acceptance criterion)
    verify — one-shot numeric envelope check |K_nu(X)| <= env(X) at the
             worst node at acceptance (independent runtime control)

    Returns TailCert. Fail-closed TailNonConvergence at x_cap.
    """
    if dps < 1:
        raise ValueError("choose_bessel_cutoff: dps must be >= 1")
    if guard < 1:
        raise ValueError("choose_bessel_cutoff: guard must be >= 1")
    ns = []
    for nd in nodes:
        if len(nd) == 2:
            nu, sigma = nd
            pref = 1
        else:
            nu, sigma, pref = nd
        ns.append((nu, sigma, pref))
    if not ns:
        raise ValueError("choose_bessel_cutoff: empty node list")

    wp = dps + guard + 20
    with mp.workdps(wp):
        tol = mpf(10) ** (-(dps + guard))
        vmin = max(_validity_min(nu) for nu, _, _ in ns)
        X = max(mpf(x_seed), vmin)
        # decade step in the dominating factor e^{-mX}: X += ln(10)/m per
        # missing digit, re-measured each round (never a formula accept).
        while True:
            worst = None
            iworst = -1
            for i, (nu, sigma, pref) in enumerate(ns):
                b = besselk_tail_bound(nu, sigma, X, m=m, prefactor=pref)
                if worst is None or b > worst:
                    worst, iworst = b, i
            if worst < tol:
                break
            if X >= x_cap:
                raise TailNonConvergence(bound=worst, tol=tol, X=X,
                                         cap=x_cap, node_index=iworst)
            missing = mp.log10(worst / tol)
            step = mp.log(10) / m * (missing + 1)
            X = min(mpf(x_cap), X + max(step, mpf(1)))
        verified = False
        if verify:
            nu, _, _ = ns[iworst]
            kv = abs(mp.besselk(mpc(nu), X)) if mpc(nu).imag != 0 or \
                mpc(nu).real != int(mpc(nu).real) else \
                abs(mp.besselk(mpc(nu).real, X))
            env = _env_value(nu, X)
            if not (kv <= env):
                raise TailEnvelopeError(
                    "envelope check FAILED at X=%s: |K_nu(X)|=%s > env=%s "
                    "(nu=%s)" % (mp.nstr(X, 6), mp.nstr(kv, 3),
                                 mp.nstr(env, 3), mp.nstr(mpc(nu), 6)))
            verified = True
        return TailCert(X=X, bound=worst, tol=tol, dps=dps, guard=guard,
                        m=int(m), node_index=iworst, x_seed=x_seed,
                        x_start=max(mpf(x_seed), vmin), verified=verified)

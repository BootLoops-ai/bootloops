#!/usr/bin/env python3
r"""nestor.singcheck -- load-time singularity verification of an integrand
spec, fail-closed.

Tanh-sinh absorbs integrable endpoint singularities (weights fall
double-exponentially), but ONLY if they are the kind the caller thinks they
are: a non-integrable pole, or a singularity on the WRONG side of what was
declared, produces confidently wrong numbers instead of a refusal.  This
module verifies the DECLARED endpoint behavior against measured samples
BEFORE any quadrature runs (the "load-time" leg):

  spec = {"a": {"type": "none"|"algebraic"|"log", "exponent": p},
          "b": {...}}          # either endpoint may be omitted => "none"

  * declaration sanity: algebraic exponent p <= -1 is non-integrable =>
    SpecRefusal at declaration time, before any sampling.
  * measured local slope: |f(a+eps)| ~ C eps^p gives
      slope = dlog|f| / dlog eps -> p
    sampled at eps = 10^-2..10^-6 (geometric); the DEEPEST pair is the
    measurement (transients at large eps are ignored).
  * "none": |f| must not GROW into the endpoint with a negative slope
    (growing samples AND slope < -0.05 = undeclared singularity =>
    refusal; the -0.05 gate catches the LOG class too, whose algebraic
    slope is only ~ -0.08, so the gate sits at -0.05).
  * "algebraic": |slope - p| <= tol (default 0.15) => pass.
  * "log": |f| must grow while the algebraic slope tends to 0
    (|slope| <= 0.2) -- the log signature.
  * any non-finite sample => refusal, always.

This is a SPEC verifier, not a certifier of the integral: it catches
mis-declared integrands at load time; the value-level certificate is the
ladder/farm agreement machinery.  No spec => nothing to verify; the caller's
receipt records "singcheck: skipped (no spec declared)" honestly.

"""
from __future__ import annotations

import mpmath as mp

from .ladder import SpecRefusal, prec_for

_TYPES = ("none", "algebraic", "log")
DEFAULT_TOL = 0.15
_EPS_EXPS = (2, 3, 4, 5, 6)      # eps = 10^-2 .. 10^-6


def _check_decl(side, d):
    typ = d.get("type", "none")
    if typ not in _TYPES:
        raise SpecRefusal("singcheck: endpoint %s declares unknown type %r "
                          "(known: %s)" % (side, typ, "/".join(_TYPES)))
    if typ == "algebraic":
        p = d.get("exponent")
        if p is None:
            raise SpecRefusal("singcheck: endpoint %s declares 'algebraic' "
                              "without an exponent" % side)
        if float(p) <= -1:
            raise SpecRefusal(
                "singcheck REFUSED (fail-closed): endpoint %s declares "
                "algebraic exponent %s <= -1 -- NON-INTEGRABLE; no "
                "quadrature will be attempted." % (side, p))
    return typ


def _slope(f, x0, sign, wp):
    """Deepest-pair local slope of log|f| vs log eps approaching x0, plus
    the |f| samples (eps descending)."""
    with mp.workprec(prec_for(wp)):
        absf = []
        for e in _EPS_EXPS:
            eps = mp.mpf(10) ** (-e)
            v = f(x0 + sign * eps)
            if v is None:
                absf.append(mp.mpf(0))
                continue
            if not (mp.isfinite(mp.re(v)) and mp.isfinite(mp.im(v))):
                raise SpecRefusal(
                    "singcheck REFUSED (fail-closed): non-finite integrand "
                    "sample at distance 1e-%d from endpoint %s"
                    % (e, mp.nstr(mp.mpf(x0), 8)))
            absf.append(abs(v))
        a1, a2 = absf[-2], absf[-1]            # eps 1e-5 -> 1e-6
        if a1 == 0 or a2 == 0:
            return 0.0, [float(x) for x in absf]
        # dlog|f|/dlog eps over the deepest decade:
        #   (log10 a2 - log10 a1) / (log10 1e-6 - log10 1e-5) = -log10(a2/a1)
        sl = -float(mp.log10(a2 / a1))
        return sl, [float(x) for x in absf]


def verify(f, a, b, spec, wp=30, tol=DEFAULT_TOL):
    """Verify the integrand spec at load time.  Returns a receipt fragment
    {side: {type, exponent?, slope, samples}}; raises SpecRefusal on any
    mismatch (fail-closed: the caller must not integrate)."""
    if spec is None:
        return {"skipped": "no spec declared"}
    rep = {}
    for side, x0, sign in (("a", mp.mpf(a), 1), ("b", mp.mpf(b), -1)):
        d = spec.get(side, {"type": "none"})
        typ = _check_decl(side, d)
        sl, samples = _slope(f, x0, sign, wp)
        growing = samples[-1] > samples[0] * 1.5 if samples[0] else \
            samples[-1] > 0
        if typ == "none":
            # growing |f| with ANY sustained negative slope is an undeclared
            # singularity; -0.05 (not -0.2) so the log class (slope ~ -0.08)
            # is caught too
            if sl < -0.05 and growing:
                raise SpecRefusal(
                    "singcheck REFUSED (fail-closed): endpoint %s declared "
                    "regular but measured slope %.3f (|f| grows toward the "
                    "endpoint) -- undeclared %s singularity."
                    % (side, sl,
                       "algebraic-type" if sl < -0.2 else "log/weak-type"))
        elif typ == "algebraic":
            p = float(d["exponent"])
            if abs(sl - p) > tol:
                raise SpecRefusal(
                    "singcheck REFUSED (fail-closed): endpoint %s declared "
                    "algebraic exponent %s but measured slope %.3f "
                    "(tol %.2f)." % (side, p, sl, tol))
        elif typ == "log":
            if not growing or abs(sl) > 0.2:
                raise SpecRefusal(
                    "singcheck REFUSED (fail-closed): endpoint %s declared "
                    "log-type but measured slope %.3f / growth %s does not "
                    "match the log signature (|slope|<=0.2, |f| growing)."
                    % (side, sl, growing))
        rep[side] = dict(type=typ, exponent=d.get("exponent"),
                         slope=round(sl, 4), samples=samples)
    return rep

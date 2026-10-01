#!/usr/bin/env python3
r"""
quad.py — THE standardized adaptive-quadrature refine-until-gate helper.
One shared helper — consumers import quad_refine; nothing here is wired into
any row script by this module itself.

THE CONTRACT (this module's own law):
    closed-form depth formula = STARTING seed only
    -> refine, measuring the engine's own error estimate at every depth
    -> accept only when the estimate beats 10^-(dps+guard)
    -> escalate by EXACT continuation of the same rule (nested levels;
       never recompute-differently)
    -> fail-closed RuntimeError at the cap NAMING the measured estimate,
       tol, depth, cap (e.g.
       "certified tail bound 1.39e-24 >= tol 1.0e-72 ... N=84 (cap 84)").
* Error estimate: double-refinement agreement — the integral is evaluated at
  nested rule depths d and d+1 and the ACCEPTED estimate is |I_{d+1} - I_d|
  (plain successive-depth agreement; deliberately MORE conservative than
  mpmath's Borwein-Bailey-Girgensohn extrapolation D1^2/D2, which is
  "not very conservative" by its own docstring). The returned value is the
  deeper evaluation I_{d+1}; for tanh-sinh / Gauss-Legendre each depth
  increment roughly doubles the accurate digits, so the true error of
  I_{d+1} sits far below the reported agreement — the guard digits absorb
  the heuristic gap (same accepted-pattern caveat as the pilot's
  trailing-window envelope, AXIS-3 PILOT SUMMARY).
* Node engine: mpmath.calculus.quadrature TanhSinh / GaussLegendre
  (mpmath 1.3.0) — nested tanh-sinh levels give the d -> d+1 escalation as
  EXACT continuation (sum_next reuses the previous level's sum; the new
  depth adds interleaved nodes, it does not re-derive anything). Tanh-sinh
  is the default: it handles endpoint singularities (the log(x) control)
  and infinite intervals; Gauss-Legendre is opt-in for smooth integrands
  (its node computation is expensive at high dps and its levels are NOT
  nested — each depth re-evaluates all nodes).
* Import-dps footgun: all
  arithmetic runs inside mp.workprec(...); global mp.dps is NEVER touched.
  Endpoint strings are converted at working precision.
* Zeno/near-pole discipline (see transport.py): a singularity ON or NEAR
  the integration path
  shows up here as non-convergence at the cap -> the fail-closed raise
  names the segment and suggests the fix (waypoint split), it never hangs
  or silently truncates.

WHAT IT DOES
============
quad_refine(f, interval, dps, *, guard=10, method="tanh-sinh", depth0=None,
            max_depth=None, wp_extra=20, scale=1, full_output=False)
    -> mpf/mpc value   (or QuadCert when full_output=True)

Integrate f along [interval[0], ..., interval[-1]] (>= 2 points; consecutive
pairs are segments — pass explicit waypoints to detour around/split at bad
points) and return the value certified so that the summed double-refinement
agreement over all segments is

    err_bound <= tol = 10^-(dps+guard) * scale.

Semantics (charter v2, all five points):
  1. NO frozen depths on the value path. Every truncation (rule depth) sits
     inside the refine-until-bound loop; acceptance is by the measured
     agreement, never by a formula. Non-convergence at the cap RAISES
     QuadNonConvergence with named diagnostics — a value the loop did not
     certify is never returned.
  2. The knobs are STARTING seeds / resource caps only:
       depth0    — first depth at which acceptance is TESTED (default:
                   mpmath guess_degree(prec), i.e. it scales with dps).
                   Levels below depth0 are accumulated (nested rule) but
                   never trusted. A different depth0 changes speed, not the
                   answer (tests pin this).
       max_depth — escalation cap (default depth0 + 6, i.e. up to 64x the
                   seeded node count). Hitting it raises; it cannot change
                   a returned value.
       wp_extra  — working-precision headroom in digits above dps+guard
                   (default 20). Too small -> the agreement plateaus at the
                   roundoff floor and the loop RAISES (fail-closed), it
                   does not return a degraded value.
       scale     — magnitude scale for the tolerance (default 1: absolute
                   tol 10^-(dps+guard), the right call for the O(1)
                   normalized constants these pipelines gate; pass the known
                   magnitude for very large/small integrals).
  3. Crank-safe: dps and dps+40 runs take genuinely different depths
     (depth0 scales with dps) and agree to dps (test_quad.py pins this).

Infinite endpoints: mp.inf / mp.ninf supported (mpmath's standard variable
transforms; (-inf, inf) is folded to [0, inf) of f(x)+f(-x) — mpmath
summation()'s own trick). Complex finite endpoints supported (straight
segments); f may return complex values.

QuadCert (full_output=True) carries the printable certificate:
value, agreement (raw gated estimate), err_bound (covers the RETURNED
dps-rounded value vs exact: agreement + measured dps rounding + wp roundoff
floor), tol, dps, guard, wp, method, depths (per segment), evals (f-call
count), segments (per-segment (a, b, depth, agreement)). Consumers print
cert.bound_line() next to values (pilot uniformity follow-up, wiring-log
item 14).

FAILURE MODES (all fail-closed, all named)
==========================================
* QuadNonConvergence(kind="cap") — agreement never beat tol by max_depth.
  Message names: measured agreement, tol, depth, cap, wp, method, segment.
* QuadNonConvergence(kind="nonfinite") — a level sum went nan/inf (true
  divergence, pole on the path, or f blowing up): named immediately, no
  further escalation.
* f raising is propagated untouched (the caller's exception is the
  diagnostic).
"""

import math

from mpmath import mp, mpf, mpc
from mpmath.calculus.quadrature import TanhSinh, GaussLegendre

__all__ = ["quad_refine", "QuadCert", "QuadNonConvergence"]

_METHODS = {
    "tanh-sinh": TanhSinh,
    "ts": TanhSinh,
    "gauss-legendre": GaussLegendre,
    "gl": GaussLegendre,
}

# rule instances cache their nodes; keep one per class (mpmath treats them
# as singletons). Node caches key on (a, b, degree, prec) so cross-call
# reuse is exact, never approximate.
_RULES = {}


def _get_rule(method):
    try:
        cls = _METHODS[method]
    except KeyError:
        raise ValueError(
            "quad_refine: unknown method %r (choose from %s)"
            % (method, sorted(set(_METHODS)))) from None
    if cls not in _RULES:
        _RULES[cls] = cls(mp)
    return _RULES[cls]


class QuadNonConvergence(RuntimeError):
    """Fail-closed refusal of quad_refine. Named diagnostics as attributes:
    kind ('cap' | 'nonfinite'), agreement, tol, depth, cap, wp, dps, guard,
    method, segment_index, a, b. The message carries all of them."""

    def __init__(self, kind, *, agreement, tol, depth, cap, wp, dps, guard,
                 method, segment_index, a, b):
        self.kind = kind
        self.agreement = agreement
        self.tol = tol
        self.depth = depth
        self.cap = cap
        self.wp = wp
        self.dps = dps
        self.guard = guard
        self.method = method
        self.segment_index = segment_index
        self.a = a
        self.b = b
        with mp.workdps(6):
            agr = mp.nstr(agreement, 3) if agreement is not None else "n/a"
            tols = mp.nstr(tol, 3)
            seg = "[%s, %s]" % (mp.nstr(a, 8), mp.nstr(b, 8))
        if kind == "nonfinite":
            head = ("quad_refine FAILED (fail-closed): non-finite level sum "
                    "(nan/inf) at depth %d" % depth)
            hint = ("integrand diverges or has a pole on/near the path — "
                    "check integrability, split at the bad point with a "
                    "waypoint, or detour into the complex plane")
        else:
            head = ("quad_refine FAILED (fail-closed): double-refinement "
                    "agreement %s >= tol %s at depth %d (cap %d)"
                    % (agr, tols, depth, cap))
            hint = ("raise max_depth / wp_extra seeds, or split the "
                    "interval at the difficult point with a waypoint; "
                    "knobs are speed seeds — the loop never returns an "
                    "uncertified value")
        msg = ("%s; method=%s, segment %d %s, dps=%d, guard=%d, wp=%d "
               "digits. %s." % (head, method, segment_index, seg, dps,
                                guard, wp, hint))
        super().__init__(msg)


class QuadCert:
    """Certificate returned by quad_refine(full_output=True).

    value      — the integral (deepest evaluation), rounded to dps
                 (the E1 transport_fixed_eps return convention)
    agreement  — summed raw double-refinement agreement over segments:
                 the estimate the acceptance GATED on (agreement <= tol)
    err_bound  — bound covering the RETURNED value vs the true integral:
                 agreement + measured dps-representation rounding + the
                 wp roundoff floor 10^-(wp-6)*max(1,|I|). This is the
                 number "|value - exact| <= err_bound" is asserted
                 against in the unit tests.
    tol        — the acceptance tolerance 10^-(dps+guard)*scale
    dps, guard, wp, method — as run
    depths     — accepted depth per segment (the returned value's depth)
    evals      — total integrand evaluations
    segments   — list of (a, b, depth, agreement) per segment
    """

    __slots__ = ("value", "agreement", "err_bound", "tol", "dps", "guard",
                 "wp", "method", "depths", "evals", "segments")

    def __init__(self, **kw):
        for k in self.__slots__:
            setattr(self, k, kw[k])

    def bound_line(self, label="quad"):
        """One-line printable certified-bound string (pilot uniformity)."""
        with mp.workdps(6):
            return ("[certified] %s: double-refinement agreement %s <= tol "
                    "%s; |err| <= %s incl. dps-rounding (dps=%d guard=%d "
                    "depth=%s method=%s evals=%d)"
                    % (label, mp.nstr(self.agreement, 3),
                       mp.nstr(self.tol, 3), mp.nstr(self.err_bound, 3),
                       self.dps, self.guard,
                       "/".join(str(d) for d in self.depths), self.method,
                       self.evals))

    def __repr__(self):
        with mp.workdps(min(self.dps, 30)):
            v = mp.nstr(self.value, min(self.dps, 30))
        with mp.workdps(6):
            e = mp.nstr(self.err_bound, 3)
        return "QuadCert(value=%s, err_bound=%s, depths=%s)" % (
            v, e, self.depths)


def _isfinite(z):
    z = mpc(z)
    return mp.isfinite(z.real) and mp.isfinite(z.imag)


def _segment(rule, f, a, b, prec, tol_seg, depth0, cap, wp, dps, guard,
             method, iseg, counter):
    """Refine ONE segment until |I_{d} - I_{d-1}| <= tol_seg.

    Returns (value_at_deepest, agreement, depth). Escalation is exact
    continuation: nested tanh-sinh levels reuse the previous sum
    (sum_next(previous=results)); Gauss-Legendre re-evaluates its full
    (non-nested) node set — same rule, deeper level, never a different
    formula."""
    if a == b:
        return mp.zero, mp.zero, 0
    # mpmath summation()'s own (-inf, inf) fold: better accuracy with 0 as
    # an endpoint.
    if a == mp.ninf and b == mp.inf:
        f_orig = f
        f = lambda x: f_orig(x) + f_orig(-x)
        a, b = mp.zero, mp.inf

    def fc(x):
        counter[0] += 1
        return f(x)

    results = []
    agreement = None
    for depth in range(1, cap + 1):
        nodes = rule.get_nodes(a, b, depth, prec)
        I = rule.sum_next(fc, nodes, depth, prec, results)
        results.append(I)
        if not _isfinite(I):
            raise QuadNonConvergence(
                "nonfinite", agreement=agreement, tol=tol_seg, depth=depth,
                cap=cap, wp=wp, dps=dps, guard=guard, method=method,
                segment_index=iseg, a=a, b=b)
        # acceptance is TESTED from the seeded depth on: the pair
        # (depth-1, depth) is the double refinement d -> d+1. Levels below
        # depth0 are accumulated (nesting) but never trusted.
        if depth >= max(depth0, 2):
            agreement = abs(results[-1] - results[-2])
            if agreement <= tol_seg:
                return results[-1], agreement, depth
    raise QuadNonConvergence(
        "cap", agreement=agreement, tol=tol_seg, depth=cap, cap=cap, wp=wp,
        dps=dps, guard=guard, method=method, segment_index=iseg, a=a, b=b)


def quad_refine(f, interval, dps, *, guard=10, method="tanh-sinh",
                depth0=None, max_depth=None, wp_extra=20, scale=1,
                full_output=False):
    """Adaptive-quadrature refine-until-gate (module docstring has the full
    contract). Returns the certified integral of f over interval at dps;
    QuadCert when full_output=True; raises QuadNonConvergence fail-closed."""
    if dps < 1:
        raise ValueError("quad_refine: dps must be >= 1 (got %r)" % (dps,))
    if guard < 1:
        raise ValueError("quad_refine: guard must be >= 1 (got %r)" % (guard,))
    pts = list(interval)
    if len(pts) < 2:
        raise ValueError("quad_refine: interval needs >= 2 points "
                         "[a, (waypoints...,) b]")
    rule = _get_rule(method)
    method = "gauss-legendre" if isinstance(rule, GaussLegendre) \
        else "tanh-sinh"

    wp = dps + guard + wp_extra          # working digits
    prec = int(wp * 3.3333) + 10         # binary precision for the rule
    nseg = len(pts) - 1

    with mp.workprec(prec):
        pts = [mp.convert(p) for p in pts]
        tol = mp.mpf(10) ** (-(dps + guard)) * mp.convert(scale)
        if not (tol > 0):
            raise ValueError("quad_refine: scale must be > 0")
        tol_seg = tol / nseg
        if depth0 is None:
            d0 = rule.guess_degree(prec)
        else:
            d0 = int(depth0)
        if d0 < 2:
            d0 = 2
        cap = (d0 + 6) if max_depth is None else int(max_depth)
        if cap < d0:
            raise ValueError(
                "quad_refine: max_depth=%d below starting depth %d"
                % (cap, d0))

        counter = [0]
        total = mp.zero
        agreement = mp.zero
        depths = []
        segments = []
        for i in range(nseg):
            a, b = pts[i], pts[i + 1]
            val, agr, depth = _segment(
                rule, f, a, b, prec, tol_seg, d0, cap, wp, dps, guard,
                method, i, counter)
            total += val
            agreement += agr
            depths.append(depth)
            segments.append((a, b, depth, agr))

        # round the RETURNED value to dps (E1 convention), then extend the
        # certificate to cover the returned representation:
        #   |value - exact| <= agreement            (truncation, gated)
        #                    + |value - total|      (dps rounding, MEASURED)
        #                    + 10^-(wp-6)*max(1,|I|) (wp roundoff floor:
        #                      <= 2^d*20 ~ 10^5.5 fdot terms of ulp(wp) slop
        #                      for every depth this loop can reach)
        with mp.workdps(dps):
            value = +total
        repr_delta = abs(mpc(value) - mpc(total))
        floor = mp.mpf(10) ** (-(wp - 6)) * max(mp.mpf(1), abs(total))
        err_bound = agreement + repr_delta + floor

    # bounds/tol are returned at working precision, NOT rounded to dps —
    # rounding a bound can shrink it.
    tol_out = tol

    if not full_output:
        return value
    return QuadCert(value=value, agreement=agreement, err_bound=err_bound,
                    tol=tol_out, dps=dps, guard=guard, wp=wp, method=method,
                    depths=depths, evals=counter[0], segments=segments)

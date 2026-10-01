#!/usr/bin/env python3
r"""nestor.ladder -- per-level precision ladder: refine-until-agreement,
fail-closed.

THE LADDER (verbatim ported math; NOT rewritten)
================================================
  -- cert_quad (refine-until-bound certified quadrature),
          digits(), the working-precision policy (WP_EXTRA=12, DEPTH_EXTRA=5,
          INNER_DROP=2, prec = int(wp*3.3333)+10), the two-height escalation
          schedule [h0, ceil(1.3 h0), ceil(1.69 h0)] (word_eval / piece_eval),
          and the exit-code refusal contract EXIT_OK/PIN/GATE/MISSING=0/2/3/4.
  src: the originating certified-quadrature evaluator,
       sha256 4d85f01dd4b54419af02804acdd372bf15439e5c03a7cbae3dc0fb59de38860a
       -- the original cert_quad + CertFail design; a sibling production
          evaluator ports it intact and this module ports the same text.
  The shipped evaluators are NOT modified by this package; this module
  extracts their LADDER PATTERN as a reusable layer.

THE PATTERN
===========
  * cert_quad: one integral, nested tanh-sinh depths escalated until the
    double-refinement agreement |I_d - I_{d-1}| beats 10^-tol_exp.  Success
    is ONLY by the measured agreement; the depth cap or a non-finite level
    sum RAISES (CertFail) with agreement/tol/depth named.  No value escapes
    a failed loop -- fail-closed.
  * ladder: a value-level escalation over an arbitrary level knob
    (grid height, md, dps, ...): evaluate at each level of a schedule and
    return as soon as two SUCCESSIVE levels agree to >= the bar; exhausting
    the schedule RAISES.  The reported digits are the agreed ones.
  * crank_check: rerun the whole evaluation at dps D+step and require
    >= D-digit agreement (the --check crank of the evaluators).
  * Per-level precision policy: working dps = D + guard (+WP_EXTRA inside
    quadratures); inner (nested) tolerances run TIGHTER by INNER_DROP.

"""
from __future__ import annotations

import mpmath as mp
from mpmath.calculus.quadrature import TanhSinh

# ---------------------------------------------------------------------------
# refusal contract (ported: the originating evaluator's exit codes + exception classes)
# ---------------------------------------------------------------------------
EXIT_OK, EXIT_PIN, EXIT_GATE, EXIT_MISSING = 0, 2, 3, 4


class NestorRefusal(RuntimeError):
    """Base of every named nestor refusal (fail-closed: no value printed)."""
    exit_code = EXIT_GATE


class PinRefusal(NestorRefusal):
    """Byte-level mutation of pinned/vendored data (exit 2)."""
    exit_code = EXIT_PIN


class GateRefusal(NestorRefusal):
    """A numeric gate or fail-closed loop missed its bar (exit 3)."""
    exit_code = EXIT_GATE


class MissingRefusal(NestorRefusal):
    """Required data or dependency absent (exit 4)."""
    exit_code = EXIT_MISSING


class CertFail(GateRefusal):
    """Fail-closed refusal: a loop could not certify its bound (ported)."""


class SpecRefusal(GateRefusal):
    """Integrand spec failed load-time verification (singcheck)."""


class AbortRefusal(NestorRefusal):
    """External abort (SIGTERM) mid-run: recorded in the receipt, ALWAYS
    re-raised; partial rungs are never presented as complete (exit 3).
    (A SIGTERMed farm run must never die without writing its receipt.)"""


class EvalRefusal(GateRefusal):
    """Integrand evaluation crashed with a non-numeric exception (e.g. a
    ZeroDivisionError at an interior pole node).  The original exception is
    CHAINED (`raise ... from`) and recorded in the receipt -- the entry
    contract is value-or-NAMED-refusal, so raw integrand crashes must not
    escape unnamed: a pole inside the contour must refuse NAMED."""


# ---------------------------------------------------------------------------
# precision policy constants (ported VERBATIM from the originating evaluators)
# ---------------------------------------------------------------------------
WP_EXTRA = 12       # working-precision headroom above dps+guard
DEPTH_EXTRA = 5     # depth cap = guess_degree + DEPTH_EXTRA
INNER_DROP = 2      # inner (nested) tolerance exponent runs TIGHTER by this
DEFAULT_GUARD = 10  # default guard digits above the certified request


def digits(a, b):
    """-log10 |a-b|/|b|: measured agreement in decimal digits.  (ported)"""
    if a == b:
        return float('inf')
    return float(-mp.log10(abs((a - b) / b)))


def wp_for(D, guard=DEFAULT_GUARD):
    """Working dps for a certified-D request (policy of the evaluators)."""
    return D + guard + WP_EXTRA


def prec_for(wp):
    """Binary working precision for working dps wp (ported: the
    int(wp*3.3333)+10 rule used by every mp.workprec block in the source)."""
    return int(wp * 3.3333) + 10


def inner_tol_exp(D, guard=DEFAULT_GUARD):
    """Tolerance exponent for a NESTED inner quadrature under an outer
    certified-D request (ported: inner_exp = D + guard + INNER_DROP)."""
    return D + guard + INNER_DROP


_TS = None


def _rule():
    global _TS
    if _TS is None:
        _TS = TanhSinh(mp.mp)
    return _TS


def cert_quad(f, a, b, tol_exp, wp, tag=""):
    """Certified quadrature: refine-until-bound, fail-closed.

    Integrate f over [a,b]; return (value, agreement, depth) with the
    double-refinement agreement |I_d - I_{d-1}| <= 10^-tol_exp.  Success is
    ONLY by the measured agreement; hitting the depth cap or a non-finite
    level sum raises CertFail with everything named.
    (Ported VERBATIM from the originating evaluator's cert_quad;
    only the module layout differs.)"""
    prec = prec_for(wp)
    with mp.workprec(prec):
        a, b = mp.convert(a), mp.convert(b)
        if a == b:
            return mp.mpf(0), mp.mpf(0), 0
        tol = mp.mpf(10) ** (-tol_exp)
        rule = _rule()
        cap = max(rule.guess_degree(prec), 3) + DEPTH_EXTRA
        results, agr = [], None
        for depth in range(1, cap + 1):
            nodes = rule.get_nodes(a, b, depth, prec)
            I = rule.sum_next(f, nodes, depth, prec, results)
            results.append(I)
            if not (mp.isfinite(mp.re(I)) and mp.isfinite(mp.im(I))):
                raise CertFail(
                    "[%s] certified quadrature FAILED (fail-closed): "
                    "non-finite level sum at depth %d on [%s, %s]" %
                    (tag, depth, mp.nstr(a, 8), mp.nstr(b, 8)))
            if depth >= 2:
                agr = abs(results[-1] - results[-2])
                if agr <= tol:
                    return results[-1], agr, depth
        # needed-dps: what the measured agreement WOULD have certified
        ach = float(-mp.log10(agr)) if (agr is not None and agr > 0) \
            else float('-inf')
        deficit = max(1, int(mp.ceil(tol_exp - ach)))
        raise CertFail(
            "[%s] certified quadrature FAILED (fail-closed): agreement %s "
            ">= tol %s at depth cap %d (wp %d) on [%s, %s] -- needed-dps: "
            "measured achievable tol exponent ~%.1f vs required %d; a "
            "request >= %d digits smaller would certify (or fix the "
            "integrand / raise wp)" %
            (tag, mp.nstr(agr, 3), mp.nstr(tol, 3), cap, wp,
             mp.nstr(a, 8), mp.nstr(b, 8), ach, tol_exp, deficit))


def escalation_schedule(D, factor_pair=(1.3, 1.69), floor=10):
    """The evaluators' three-rung level schedule (ported VERBATIM):
    h0 = max(floor, ceil((D+4)/2)); rungs [h0, ceil(1.3 h0), ceil(1.69 h0)].
    Interpretation of the knob (grid height denominator, md, ...) is the
    caller's; the ESCALATION LAW is the original one."""
    import math
    h0 = max(floor, int(math.ceil((D + 4) / 2.0)))
    return [h0] + [int(math.ceil(fc * h0)) for fc in factor_pair]


def ladder(eval_at_level, levels, D, bar=None, tag="ladder", wp=None,
           on_rung=None):
    """Refine-until-agreement over an explicit level schedule, fail-closed.

    eval_at_level(level) -> mpf value at that level.  Returns
    (value, agreed_digits, rungs) as soon as two SUCCESSIVE levels agree to
    >= bar digits (default D); exhausting the schedule raises CertFail with
    the schedule and bar named.  rungs = list of per-level dicts
    {level, wall_s, value_str, agreement_digits} for the receipt.
    (Pattern ported from the originating evaluator's word_eval/piece_eval two-height
    escalation; generalized only in the knob, not in the law.)"""
    import time
    bar = D if bar is None else bar
    wp = wp_for(D) if wp is None else wp
    prev = None
    rungs = []
    with mp.workprec(prec_for(wp)):
        for li, lv in enumerate(levels):
            t0 = time.time()
            val = eval_at_level(lv)
            wall = time.time() - t0
            agr = None if prev is None else digits(val, prev)
            rungs.append(dict(level=lv, wall_s=round(wall, 2),
                              value_str=mp.nstr(val, min(int(wp), 40)),
                              agreement_digits=(None if agr is None
                                                else round(agr, 2))))
            if on_rung is not None:
                on_rung(rungs[-1])
            if agr is not None and agr >= bar:
                return val, agr, rungs
            prev = val
    # needed-dps: the best successive-rung agreement IS the achievable bar
    agrs = [r["agreement_digits"] for r in rungs
            if r["agreement_digits"] is not None
            and r["agreement_digits"] == r["agreement_digits"]]  # drop nan
    best = max(agrs) if agrs else float('-inf')
    raise CertFail(
        "[%s] ladder FAILED (fail-closed): levels %s exhausted without two "
        "successive rungs agreeing to %s d (measured rungs: %s) -- "
        "needed-dps: best successive agreement %.1f d; a request of "
        "dps <= %d would certify on this schedule"
        % (tag, list(levels), bar,
           [r["agreement_digits"] for r in rungs], best,
           int(best) if best > 0 else 0))


def crank_check(run, D, step=8, tag="crank"):
    """The --check crank (ported design): run(D) and run(D+step) must agree
    to >= D digits, fail-closed.  Returns the measured agreement."""
    v1 = run(D)
    v2 = run(D + step)
    agr = digits(v2, v1)
    if agr < D:
        raise CertFail(
            "[%s] crank check FAILED (fail-closed): dps %d vs %d values "
            "agree to only %.1f d (bar %d) -- needed-dps: certify at most "
            "dps <= %d on this evaluation" % (tag, D, D + step, agr, D,
                                              max(0, int(agr))))
    return agr

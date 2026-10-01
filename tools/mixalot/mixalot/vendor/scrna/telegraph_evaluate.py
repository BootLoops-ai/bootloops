#!/usr/bin/env python3
"""telegraph_evaluate.py -- certified telegraph-model PMF / log-likelihood evaluator.

SINGLE-FILE deliverable.  Pip-only dependencies:
numpy, mpmath, python-flint (Arb), sympy (sympy imported lazily, selftest only).
The tail bound of `hyp1f1_pos_arb` is certified: the double q enters only
the tail BOUND, never the tail VALUE.

================================================================================
CONVENTIONS (embedded dictionary; full provenance in pilot/core/CONVENTIONS.md)
================================================================================
PRIMARY convention: Beta-mixture (a, b, r), degradation-rate units delta = 1:

    a := k_on  / delta      (promoter OFF->ON rate; FIRST Beta shape)
    b := k_off / delta      (promoter ON->OFF rate; second Beta shape)
    r := k_syn / delta      (mRNA synthesis rate while ON; Poisson scale)

    P(n) = int_0^1 Pois(n; r u) Beta(u; a, b) du,
           Beta(u; a, b) = u^(a-1) (1-u)^(b-1) / B(a, b)

Closed form (verified: quadrature 40+ digits, C2; truncated master equation to
4e-19, C5/D4 of the adversarial verification suite -- the CME check pins
a = k_on as the FIRST Beta shape; the swapped convention fails by factor ~39):

    P(n) = (r^n / n!) * B(a+n, b)/B(a, b) * 1F1(a+n; a+b+n; -r)            (*)
         = e^{-r} (r^n/n!) ((a)_n/(a+b)_n) * 1F1(b; a+b+n; r)              (***)

with 1F1(alpha; gamma; z) = sum_k (alpha)_k/((gamma)_k k!) z^k (alpha in the
NUMERATOR, (x)_k RISING factorial, DLMF 13.2.2 "M"), and the Kummer transform
1F1(a+n; a+b+n; -r) = e^{-r} 1F1(b; a+b+n; r) (DLMF 13.2.39; verified to 40+
digits, C4, and re-derived to ~1e-60 adversarially).  Form (***) has
ALL-POSITIVE terms for a, b, r > 0 -- cancellation-free; it is the certified
workhorse everywhere in this file.

Convention dictionary (all mappings verified at build time):

  | source                          | map to primary (a, b, r)                  |
  |---------------------------------|-------------------------------------------|
  | Beta-mixture (PRIMARY)          | identity                                   |
  | rates (k_on, k_off, k_syn, delta)| a=k_on/delta, b=k_off/delta, r=k_syn/delta|
  | legacy convention "b_wo"         | b_wo = a + b (Kummer DENOMINATOR base,    |
  |                                 |   NOT the second Beta shape)               |
  | txburst (Larsson et al. 2019)   | kon=a, koff=b, ksyn=r (identity; read from |
  |                                 |   txburstML.py code, SHA-256 82260daa..5257)|
  | Peccoud-Ycart 1995              | a=lambda/delta, b=mu/delta, r=nu/delta     |

CAUTION: rates fitted in per-hour units MUST be divided by delta before entering
these formulas (the CLI does this when given --kon/--koff/--ksyn/--delta).

RECURRENCE (two independent symbolic derivations, residual 0; pilot RECURRENCE.md):

    n (n+1) P(n+1) = n (n + a + b + r - 1) P(n) - r (n + a - 1) P(n-1),  n >= 1.

STABILITY (measured + Perron/Miller theory; adversarially confirmed): the wanted
solution is MINIMAL under upward recursion at ALL n; BACKWARD recursion from two
seeds at nmax+1, nmax is the stable direction everywhere (~1 digit lost per 2000
steps); forward recursion is unstable at essentially all n, including n < r.

================================================================================
WHAT "CERTIFIED" MEANS
================================================================================
mode='certified' returns arb BALLS (midpoint +/- radius, both dyadic): the true
P(n) is mathematically guaranteed to lie inside the ball.  Sources of rigor:
  - the all-positive series is summed in ball arithmetic; the truncation tail is
    enclosed by a geometric bound computed ENTIRELY in ball arithmetic and added
    as T*[0,1] (pilot adversarial fix: no double-precision quantity enters the
    bound value);
  - the backward vector recurrence propagates balls; because interval RADII grow
    like the dominant solution even in the stable direction, the implementation
    pre-estimates the growth, adds working precision, then VERIFIES achieved
    relative widths against 2^-prec and escalates (fail-closed RuntimeError if
    it cannot certify);
  - log-likelihoods use arb log of positive balls (rigorous), integer-weighted
    ball sums; the reported radius is a certified error bound.

================================================================================
MODES
================================================================================
  'certified' : arb balls at `prec` bits (default 150).  Seconds-class to
                nmax=2000.  Output: decimal mid/rad strings at requested digits.
  'mp'        : mpmath backward recurrence evaluated at TWO precisions
                (dps+10 vs dps+25); returned only if they agree to >= dps
                digits everywhere, else precision is escalated (x2 steps, 4
                attempts) and finally fails closed (RuntimeError).
  'float'     : fast float64 log-space production route (log-sum-exp style
                rescaled all-positive Kummer series) with the a-priori
                escalation policy below.  << 1 s to nmax=2000.

================================================================================
FLOAT-LANE ESCALATION POLICY (a-priori, documented, fail-closed)
================================================================================
The float lane computes log P(n) = lpref(n) + log 1F1(b; a+b+n; r) with
  lpref(n) = -r + n log r - lgamma(n+1) + sum_{j<n} log((a+j)/(a+b+j))
(the Pochhammer ratio via a compensated cumsum of per-term logs, each computed
in its well-conditioned branch -- avoids both lgamma cancellation and the
log1p-near(-1) conditioning trap at a << b), and the
series summed with Kahan compensation + decimal rescaling (equivalent to
log-sum-exp over the all-positive terms; no overflow to r = 1e4 and beyond).

A-priori absolute-error model for log P(n) in float64 (eps = 2^-52):

  ERR(n) = eps * [ 2*( r + |n log r| + lgamma(n+1) + |Pochhammer cumsum| ) <- T1a
                   + 2*|log 1F1|                                           <- T1b
                   + 3*n + 4*sqrt(K) + 32 ]                                <- T1c

T1a/T1b: each log-space component carries a few-ulp error PROPORTIONAL TO ITS
MAGNITUDE (this is the intrinsic float64 limit: representing log P ~ -1e4 in
float64 already costs ~2e-12 relative in P; no float64 implementation can beat
it).  T1c: rounding accumulation over the K series terms (K ~ 1.12*r + 9*sqrt(r),
Kahan-compensated so sqrt(K), not K).  Since a relative deviation delta in P is
|exp(dlog)-1| ~ dlog, ERR(n) directly bounds the relative PMF error.

TRIGGERS (checked on every float-lane call):
  T1 (condition / magnitude): ERR(n) > FLOAT_TOL = 1e-12 for any requested n
      -> escalate.  (Fires automatically for large r, large n log r, i.e. the
      exact regime where the pilot measured >1e-12 deviations.)
  T2 (series length / non-convergence): terms exhausted before the tail test
      t/S/(1-rho) < 1e-17 passes (K > K_LIMIT) -> escalate.
  T3 (underflow proximity): log P(n) < -700 -> linear-space P is not
      representable in float64; the log-space output remains valid and such
      entries are MARKED (p_underflow); linear-P consumers must escalate.
  T4 (non-finite catch-all): any NaN/inf in the pipeline -> escalate.

ESCALATION ACTION (policy 'auto', default): the WHOLE vector is recomputed in
mode='mp' (dps 30, two-precision self-check) and the output is tagged
route='mp'.  Policy 'off': raw float values are returned WITH the trigger flags
set and escalated=False -- callers see exactly what fired (used for policy
validation); the JSON carries "unsafe_entries" so degradation is never silent.
Validated on the pilot 228-cell grid: see GATE_EVALUATOR.json
(float_policy_validation block) and README.md.

Log-likelihood in float mode: certified-model bound ERR_LL = sum_i ERR(n_i)
(+ eps*sum|log P|); escalates when ERR_LL > 1e-9 nats or any used entry flags.

================================================================================
CLI
================================================================================
  python3 telegraph_evaluate.py pmf    --a 1.3 --b 4.7 --r 1000 --nmax 50 \
                                       --mode certified --prec 150 --digits 30
  python3 telegraph_evaluate.py pmf    --kon 3.4 --koff 9.4 --ksyn 180 --delta 2.0 \
                                       --nmax 100 --mode float
  python3 telegraph_evaluate.py loglik --a 0.4 --b 9.5 --r 300 --counts 0,3,7,2,9 \
                                       --mode certified
  python3 telegraph_evaluate.py loglik --a 0.4 --b 9.5 --r 300 \
                                       --counts-file counts.txt --mode float
  python3 telegraph_evaluate.py --selftest        # exact-rational lane + mutation control
  python3 telegraph_evaluate.py --version-info

Outputs deterministic JSON on stdout (sorted keys, no timestamps).  No network.

Precision discipline: no module-level mpf/arb constants; every mpf/arb
construction sits inside mpmath.mp.workdps / flint.ctx.workprec.
"""

import argparse
import hashlib
import json
import math
import os
import sys
from fractions import Fraction

import mpmath
import numpy as np
from flint import arb, ctx as _flint_ctx

__version__ = "1.0.0"

# float64 machine epsilon (plain float constant; not an mp/arb object)
_EPS64 = 2.220446049250313e-16

# float-lane policy constants (documented in the module docstring)
FLOAT_TOL = 1e-12          # T1 trigger threshold on the a-priori log-error model
FLOAT_UNDERFLOW_LOG = -700.0   # T3 linear-space underflow proximity
LOGLIK_FLOAT_TOL = 1e-9    # nats; float-lane loglik escalation threshold


# ----------------------------------------------------------------------------
# input coercion (exact; no precision dependence)
# ----------------------------------------------------------------------------

def _as_arb(x):
    """Coerce to arb exactly (int/float/str/Fraction/arb). Fractions via exact division."""
    if isinstance(x, arb):
        return x
    if isinstance(x, Fraction):
        return arb(x.numerator) / arb(x.denominator)   # rounded at current prec, ball-correct
    return arb(x)


def _as_mpf(x):
    """Coerce to mpf at current working dps (exact for int/float/str; Fraction via division)."""
    if isinstance(x, Fraction):
        return mpmath.mpf(x.numerator) / mpmath.mpf(x.denominator)
    return mpmath.mpf(x)


def arb_mid_mpf(x):
    """Midpoint of an arb ball as an mpf, rounded at the CURRENT mpmath working
    precision (use inside mp.workdps with enough digits). Bridge for cross-route
    comparisons; the ball radius must be inspected separately via x.rad()."""
    mid = x.mid()
    if mid == 0:
        return mpmath.mpf(0)
    m, e = mid.man_exp()
    return mpmath.ldexp(mpmath.mpf(int(m)), int(e))


def abr_from_rates(kon, koff, ksyn, delta):
    """Convert rate-parameter inputs (k_on, k_off, k_syn, delta) to the PRIMARY
    Beta-mixture convention: a = k_on/delta, b = k_off/delta, r = k_syn/delta.
    Exact when the inputs are Fractions; float division otherwise."""
    if any(v is None for v in (kon, koff, ksyn)):
        raise ValueError("rate input requires kon, koff, ksyn (and optional delta)")
    if isinstance(kon, Fraction) or isinstance(koff, Fraction) or isinstance(ksyn, Fraction):
        d = Fraction(delta)
        return Fraction(kon) / d, Fraction(koff) / d, Fraction(ksyn) / d
    d = float(delta)
    return float(kon) / d, float(koff) / d, float(ksyn) / d


# ----------------------------------------------------------------------------
# certified workhorse: all-positive Kummer-transformed series in arb balls
# [certified tail bound -- do not modify without re-running the
#  adversarial verification]
# ----------------------------------------------------------------------------

def hyp1f1_pos_arb(b, c, r, prec, max_terms=5_000_000):
    """Certified 1F1(b; c; r) for 0 < b < c, r >= 0, as an arb ball at `prec` bits.

    Series sum_k (b)_k / ((c)_k k!) r^k.  ALL terms positive (no cancellation).
    Term ratio rho_k = (b+k) r / ((c+k)(k+1)) <= r/(k+1) because b < c.
    Truncating after term K-1 (first omitted term t_K), with r/(K+1) < 1,
    the tail is <= t_K / (1 - r/(K+1)) (geometric bound).  The bound is computed
    ENTIRELY in ball arithmetic (denominator 1 - r/(K+1) as a ball with
    ball-certain positivity check; a double-precision q is used for loop control
    only) and added as the ball T * [0,1], which contains [0, tail]; hence the
    returned ball is a rigorous enclosure.
    [The double q must enter only the tail BOUND, never the tail VALUE:
    through t/(1-arb(q)), nearest-rounding could put q below sup(r)/(K+1) and
    shave ~1e-16 of the tail bound, below the ball's own radius at prec+20.]
    """
    with _flint_ctx.workprec(prec + 20):
        b_ = _as_arb(b)
        c_ = _as_arb(c)
        r_ = _as_arb(r)
        if not (b_ > 0 and c_ > b_ and r_ >= 0):
            raise ValueError("hyp1f1_pos_arb requires 0 < b < c and r >= 0 (as verifiable balls)")
        # upper bound on r as a float for the q test (heuristic loop control only;
        # the rigor comes from ball arithmetic on T below)
        r_ub = float(r_.mid()) + float(r_.rad())
        eps2 = arb(2) ** (-(prec + 10))
        S = arb(1)
        t = arb(1)
        k = 0
        while True:
            # t is currently t_k; compute t_{k+1}
            t = t * (b_ + k) * r_ / ((c_ + k) * (k + 1))
            k += 1
            # now t = t_k with the new k; decide whether t_k is the first OMITTED term
            q = r_ub / (k + 1)
            if q < 0.5:
                # candidate stop: tail <= t_k / (1 - r/(k+1)), computed as balls;
                # the double q above is loop control ONLY (its nearest-rounding
                # may sit below sup(r)/(k+1) and must not enter the bound)
                denom = 1 - r_ / (k + 1)      # ball; division rigorous
                if denom > 0:                 # ball-certain positivity
                    T = t / denom             # ball containing the true geometric bound
                    if T.mid() < (eps2 * S).mid() or k >= max_terms:
                        S = S + T * arb("0.5", "0.5")   # + [0,1]*T  (contains [0, tail])
                        break
            S = S + t
            if k >= max_terms:
                raise RuntimeError("hyp1f1_pos_arb: max_terms exceeded without convergence")
        return arb(S)  # snapshot at prec+20; caller sees a valid ball at any prec


def hyp1f1_native_arb(b, c, r, prec):
    """Native flint/Arb 1F1(b; c; r) (independent cross-check route ONLY).

    KNOWN LIMIT: returns NaN for b ~ 0.05, c >~ 1050, r ~ 2000 at ANY precision
    (tested to 2000 bits in the pilot).  Never used by production routes; the
    all-positive series `hyp1f1_pos_arb` is the workhorse there."""
    with _flint_ctx.workprec(prec):
        return _as_arb(r).hypgeom_1f1(_as_arb(b), _as_arb(c))


def _prefactor_arb(a_, b_, r_, n):
    """e^{-r} r^n/n! (a)_n/(a+b)_n as an arb ball at the CURRENT working precision."""
    return (-r_).exp() * (r_ ** n) / arb.fac_ui(n) * a_.rising(n) / (a_ + b_).rising(n)


def pmf_certified_single(a, b, r, n, prec, native=False):
    """Certified arb ball for P(n), Beta-mixture convention.

    P(n) = e^{-r} r^n/n! (a)_n/(a+b)_n * 1F1(b; a+b+n; r).
    native=False -> all-positive series (certified workhorse);
    native=True  -> flint's built-in 1F1 (independent certified route; may NaN,
                    see hyp1f1_native_arb docstring).
    """
    with _flint_ctx.workprec(prec + 30):
        a_, b_, r_ = _as_arb(a), _as_arb(b), _as_arb(r)
        if native:
            M = r_.hypgeom_1f1(b_, a_ + b_ + n)
        else:
            M = hyp1f1_pos_arb(b_, a_ + b_ + n, r_, prec + 30)
        return _prefactor_arb(a_, b_, r_, n) * M


def _margin_bits_estimate(a, b, r, nmax):
    """Cheap float estimate (bits) of ball-RADIUS growth under the backward
    recurrence (radii propagate with |coefficients|, i.e. like the dominant
    solution, even though midpoints are backward-stable).  Heuristic only; the
    certified route verifies achieved widths and escalates."""
    af, bf, rf = float(_as_arb(a).mid()), float(_as_arb(b).mid()), float(_as_arb(r).mid())
    A = af + bf + rf - 1
    margin = 0.0
    for n in range(1, nmax + 1):
        fn = (n * (n + A) + n * (n + 1)) / (rf * (n + af - 1)) if n + af - 1 > 0 else 4.0
        margin += max(0.0, math.log2(abs(fn)))
    return margin


def pmf_certified(a, b, r, nmax, prec, method="auto"):
    """Certified PMF vector [P(0), ..., P(nmax)] as arb balls.

    method='series': the all-positive series at every n (O(nmax * nterms); each
        entry an independent certified evaluation; used when the radius-growth
        margin of the recurrence would be extreme).
    method='recurrence': BACKWARD ball three-term recurrence from two certified
        series seeds at n = nmax+1, nmax down to n = 0.  Backward is the stable
        direction at EVERY n (wanted solution minimal under upward recursion
        against Q_n = C_n U(b; a+b+n; r); pilot RECURRENCE.md).  Ball arithmetic
        keeps every entry a rigorous enclosure regardless; backward keeps the
        relative widths ~2^-prec.  O(nmax + 2 series seeds).
    method='auto' (default): 'recurrence' unless the estimated radius margin
        exceeds 200k bits (protective cap), then 'series'.

    Fail-closed: every returned ball is VERIFIED to satisfy mid > 0 and
    rad < mid * 2^-prec; working precision escalates (4 attempts) else
    RuntimeError.
    """
    if method == "auto":
        method = "recurrence" if _margin_bits_estimate(a, b, r, nmax) <= 200_000 else "series"
    if method == "series":
        wp = prec
        for _attempt in range(4):
            P = [pmf_certified_single(a, b, r, n, wp) for n in range(nmax + 1)]
            with _flint_ctx.workprec(wp + 30):
                ok = all((p.mid() > 0 and (p.rad() + arb(0)) < (p.mid() + arb(0)) * arb(2) ** (-prec))
                         for p in P)
            if ok:
                return P
            wp = wp + max(64, wp // 2)
        raise RuntimeError("pmf_certified(series): could not certify widths after escalation")
    if method != "recurrence":
        raise ValueError("method must be 'auto', 'series' or 'recurrence'")

    # Ball radii do NOT enjoy the backward stability of midpoint errors: interval
    # radius propagation uses |coefficients|, so radii grow like the dominant
    # solution of the |coefficient| recurrence (~prod (2n+A+r)/(rn) factors).
    # Estimate that growth in cheap float arithmetic and add it as extra working
    # precision; then verify achieved widths and escalate if the estimate was
    # short (at most 4 attempts; each attempt is only O(nmax) ball ops).
    wp = prec + 60 + int(_margin_bits_estimate(a, b, r, nmax))

    for _attempt in range(4):
        with _flint_ctx.workprec(wp):
            a_, b_, r_ = _as_arb(a), _as_arb(b), _as_arb(r)
            P = [None] * (nmax + 1)
            Pn1 = _prefactor_arb(a_, b_, r_, nmax + 1) * hyp1f1_pos_arb(b_, a_ + b_ + nmax + 1, r_, wp)
            Pn = _prefactor_arb(a_, b_, r_, nmax) * hyp1f1_pos_arb(b_, a_ + b_ + nmax, r_, wp)
            P[nmax] = Pn
            for n in range(nmax, 0, -1):
                # P_{n-1} = [ n (n+a+b+r-1) P_n - n (n+1) P_{n+1} ] / ( r (n+a-1) )
                Pm = (n * (n + a_ + b_ + r_ - 1) * Pn - n * (n + 1) * Pn1) / (r_ * (n + a_ - 1))
                P[n - 1] = Pm
                Pn1, Pn = Pn, Pm
            # verify achieved relative widths against the 2^-prec target
            ok = True
            for n in range(nmax + 1):
                m, rd = P[n].mid(), P[n].rad() + arb(0)
                if not (m > 0 and rd < m * arb(2) ** (-prec)):
                    ok = False
                    break
            if ok:
                return P
        wp = wp + max(64, wp // 2)   # escalate and retry
    raise RuntimeError("pmf_certified: could not certify widths after escalation")


def loglik_certified(counts, a, b, r, prec):
    """Certified log-likelihood ball for integer counts under P(.; a, b, r).

    Sum over unique counts n of mult(n) * log(P(n)), with P(n) certified balls
    and log taken in ball arithmetic (rigorous: P-balls are verified strictly
    positive by pmf_certified).  Returns (ball, nmax_used)."""
    counts = _check_counts(counts)
    nmax = max(counts)
    P = pmf_certified(a, b, r, nmax, prec)
    mult = {}
    for c in counts:
        mult[c] = mult.get(c, 0) + 1
    with _flint_ctx.workprec(prec + 30):
        total = arb(0)
        for n, m in sorted(mult.items()):
            if not (P[n] > 0):   # defense in depth: arb log of a 0-containing
                raise RuntimeError(  # ball is a QUIET nan; never accept one
                    "loglik_certified: P(%d) ball not certainly positive" % n)
            total = total + arb(m) * P[n].log()
        return arb(total), nmax


# ----------------------------------------------------------------------------
# mpmath route (independent implementation) + two-precision self-check
# ----------------------------------------------------------------------------

# mpmath's hyp1f1 defaults to maxterms=6000; the all-positive series at z = r
# needs ~ r + O(sqrt(r)*digits) terms, so r >~ 5e3 spuriously raises NoConvergence
# (found adversarially at (a,b,r) = (1e-3, 999, 1e4): the large-b asymptotic
# branch also fails there, and the mp lane crashed with a raw traceback instead
# of the documented fail-closed behavior).  2e6 caps r ~ 1.9e6, far beyond the
# certification box; convergence still stops at the tail test, so no extra cost.
_MP_HYP_MAXTERMS = 2_000_000


def pmf_mp_single(a, b, r, n, dps):
    """P(n) by direct per-n mpmath evaluation (Kummer-transformed positive form)."""
    with mpmath.mp.workdps(dps):
        a_, b_, r_ = _as_mpf(a), _as_mpf(b), _as_mpf(r)
        M = mpmath.hyp1f1(b_, a_ + b_ + n, r_, maxterms=_MP_HYP_MAXTERMS)
        pref = (mpmath.e ** (-r_)) * (r_ ** n) / mpmath.factorial(n) \
            * mpmath.rf(a_, n) / mpmath.rf(a_ + b_, n)
        return pref * M


def pmf_recurrence(a, b, r, nmax, dps, direction="auto"):
    """PMF vector [P(0..nmax)] via the three-term recurrence (mpmath floats).

        n (n+1) P_{n+1} = n (n+a+b+r-1) P_n - r (n+a-1) P_{n-1},  n >= 1

    direction='auto' == 'backward' (default, STABLE at every n): seeds P(nmax+1),
        P(nmax) from two direct 1F1 evaluations, sweep down to n=0.  The wanted
        Poisson-Beta solution is the MINIMAL solution of the recurrence under
        upward recursion at ALL n (not just n>r; pilot RECURRENCE.md), so downward
        recursion is the stable direction everywhere: measured loss ~1 digit
        over sweeps of 2000 steps across all tested regimes.
    direction='forward': seeds P(0), P(1), sweep up.  UNSTABLE (kept for
        experiments only): e.g. ~17 digits lost by n=10 at (1.3,4.7,500).
    """
    with mpmath.mp.workdps(dps):
        a_, b_, r_ = _as_mpf(a), _as_mpf(b), _as_mpf(r)

        def seed(n):
            return ((mpmath.e ** (-r_)) * (r_ ** n) / mpmath.factorial(n)
                    * mpmath.rf(a_, n) / mpmath.rf(a_ + b_, n)
                    * mpmath.hyp1f1(b_, a_ + b_ + n, r_, maxterms=_MP_HYP_MAXTERMS))

        P = [mpmath.mpf(0)] * (nmax + 1)
        if direction == "forward":
            P[0] = seed(0)
            if nmax >= 1:
                P[1] = seed(1)
            for n in range(1, nmax):
                P[n + 1] = ((n + a_ + b_ + r_ - 1) * P[n]
                            - r_ * (n + a_ - 1) / n * P[n - 1]) / (n + 1)
        elif direction in ("auto", "backward"):
            Pn1 = seed(nmax + 1)
            Pn = seed(nmax)
            P[nmax] = Pn
            for n in range(nmax, 0, -1):
                Pm = (n * (n + a_ + b_ + r_ - 1) * Pn - n * (n + 1) * Pn1) / (r_ * (n + a_ - 1))
                P[n - 1] = Pm
                Pn1, Pn = Pn, Pm
        else:
            raise ValueError("direction must be 'forward', 'backward', or 'auto'")
        return P


def pmf_mp_selfcheck(a, b, r, nmax, dps=30, max_attempts=4):
    """mp-mode PMF with TWO-PRECISION SELF-CHECK (fail-closed).

    Evaluates the backward recurrence at dps+10 and dps+25; accepts only if the
    two agree to >= dps digits at EVERY n; otherwise escalates the working
    precision (adds 30 digits per attempt, up to max_attempts) and finally
    raises RuntimeError.  Returns (values_at_higher_precision, report_dict)."""
    import mpmath.libmp.libhyper as _libhyper
    d1 = dps + 10
    for attempt in range(max_attempts):
        d2 = d1 + 15
        try:
            P1 = pmf_recurrence(a, b, r, nmax, d1)
            P2 = pmf_recurrence(a, b, r, nmax, d2)
        except _libhyper.NoConvergence as e:
            raise RuntimeError(
                "pmf_mp_selfcheck: mpmath series did not converge (%s); "
                "fail-closed -- use mode='certified'" % e) from e
        with mpmath.mp.workdps(d2 + 10):
            worst = mpmath.mpf('inf')
            for n in range(nmax + 1):
                if P2[n] == 0:
                    worst = -mpmath.mpf('inf') if P1[n] != 0 else worst
                    continue
                dv = abs(P1[n] - P2[n]) / abs(P2[n])
                dgt = mpmath.mpf(d1) if dv == 0 else -mpmath.log10(dv)
                worst = min(worst, dgt)
            worst_f = float(worst)
        if worst_f >= dps:
            return P2, {"dps_requested": dps, "dps_pair": [d1, d2],
                        "min_agreement_digits": round(worst_f, 1),
                        "attempts": attempt + 1, "passed": True}
        d1 = d2 + 15
    raise RuntimeError(
        "pmf_mp_selfcheck: two-precision agreement < %d digits after %d escalations "
        "(fail-closed; use mode='certified')" % (dps, max_attempts))


def loglik_mp(counts, a, b, r, dps=30):
    """mp-mode log-likelihood with two-precision self-check on the final value.

    Returns (loglik_mpf_str, report)."""
    counts = _check_counts(counts)
    nmax = max(counts)
    P, rep = pmf_mp_selfcheck(a, b, r, nmax, dps=dps)
    mult = {}
    for c in counts:
        mult[c] = mult.get(c, 0) + 1
    with mpmath.mp.workdps(rep["dps_pair"][1] + 10):
        total = mpmath.mpf(0)
        for n, m in sorted(mult.items()):
            total += m * mpmath.log(P[n])
        return total, rep


# ----------------------------------------------------------------------------
# float fast lane (log-space, float64) + a-priori escalation policy
# ----------------------------------------------------------------------------

class FloatLaneEscalation(RuntimeError):
    """Raised internally when a float-lane trigger fires under policy 'error'."""


def _log1f1_pos_vec(b, c_vec, r, tol=1e-17):
    """log 1F1(b; c; r) for a VECTOR of denominator parameters c (= a+b+n), r >= 0.

    All-positive Kummer series in float64, Kahan-compensated summation,
    per-element decimal rescaling (log-sum-exp equivalent; immune to overflow
    since 1F1 ~ e^r).  Returns (log_values ndarray, K terms used).
    Raises FloatLaneEscalation on non-convergence (trigger T2)."""
    c = np.asarray(c_vec, dtype=np.float64)
    S = np.ones_like(c)
    comp = np.zeros_like(c)          # Kahan compensation
    t = np.ones_like(c)
    off = np.zeros_like(c)           # accumulated log offsets from rescaling
    log_rescale = 200.0 * math.log(10.0)
    k = 0
    k_limit = int(1.2 * r + 60.0 * math.sqrt(r + 1.0) + 800)
    while True:
        for _ in range(32):
            t = t * ((b + k) * r) / ((c + k) * (k + 1.0))
            y = t - comp             # Kahan add
            snew = S + y
            comp = (snew - S) - y
            S = snew
            k += 1
        m = S > 1e200
        if m.any():
            S[m] *= 1e-200
            t[m] *= 1e-200
            comp[m] *= 1e-200
            off[m] += log_rescale
        rho = r / (k + 1.0)
        if rho < 0.9 and float(np.max(t / S)) / (1.0 - rho) < tol:
            break
        if k > k_limit:
            raise FloatLaneEscalation("T2: series not converged after %d terms (r=%g)" % (k, r))
    return np.log(S) + off, k


def logpmf_float_raw(a, b, r, nmax):
    """Raw float64 log-space PMF with the a-priori error model (NO escalation).

    Returns dict with:
      logP      : ndarray, log P(n) for n = 0..nmax
      err_log   : ndarray, a-priori absolute-error bound on logP (== rel err of P)
      flags_t1  : boolean ndarray, ERR(n) > FLOAT_TOL      (trigger T1)
      flags_t3  : boolean ndarray, logP < FLOAT_UNDERFLOW_LOG (marker T3)
      K         : series terms used
    Raises FloatLaneEscalation for T2/T4-class failures."""
    a = float(a); b = float(b); r = float(r)
    if not (a > 0 and b > 0 and r > 0):
        raise ValueError("float lane requires a, b, r > 0")
    n = np.arange(nmax + 1, dtype=np.float64)
    lr = math.log(r)
    # lgamma(n+1), few-ulp accurate (C library)
    lg = np.array([math.lgamma(i + 1.0) for i in range(nmax + 1)])
    # Pochhammer ratio log (a)_n/(a+b)_n = cumsum of log((a+j)/(a+b+j)), Kahan
    # cumsum.  CONDITIONING BRANCH: log1p(-b/(a+b+j)) is catastrophically
    # ill-conditioned when b/(a+b+j) ~ 1 (i.e. a+j << b; amplification
    # 1/(1 - b/(a+b+j))); there use log((a+j)/(a+b+j)) whose argument is small
    # and perfectly conditioned.  [Found by the 228-cell policy validation at
    # (a,b)=(0.03,999): 2.4e-12 uncaught deviation before this branch.]
    cum = np.zeros(nmax + 1)
    s = 0.0
    kc = 0.0
    for j in range(nmax):
        if b <= a + j:      # ratio b/(a+b+j) <= 1/2: log1p well-conditioned
            term = math.log1p(-b / (a + b + j))
        else:               # (a+j)/(a+b+j) < 1/2: plain log well-conditioned
            term = math.log((a + j) / (a + b + j))
        y = term - kc
        snew = s + y
        kc = (snew - s) - y
        s = snew
        cum[j + 1] = s
    logS, K = _log1f1_pos_vec(b, a + b + n, r)
    lpref = -r + n * lr - lg + cum
    logP = lpref + logS
    if not np.all(np.isfinite(logP)):
        raise FloatLaneEscalation("T4: non-finite value in float lane")
    # a-priori error model (module docstring): few-ulp per log-space component;
    # the 3*n term covers per-term rounding accumulation in the Pochhammer
    # cumsum (each of the n terms carries ~2 ulp of input rounding)
    err_log = _EPS64 * (2.0 * (r + np.abs(n * lr) + lg + np.abs(cum))
                        + 2.0 * np.abs(logS)
                        + 3.0 * n
                        + 4.0 * math.sqrt(K) + 32.0)
    return {"logP": logP, "err_log": err_log,
            "flags_t1": err_log > FLOAT_TOL,
            "flags_t3": logP < FLOAT_UNDERFLOW_LOG,
            "K": K}


def pmf_float(a, b, r, nmax, policy="auto"):
    """Float fast-lane PMF with the documented escalation policy.

    policy='auto'  : if any T1/T2/T4 trigger fires, recompute the WHOLE vector
                     in mode='mp' (dps 30, self-checked); route tag tells which.
    policy='off'   : return raw float values with flags (validation use; the
                     output marks unsafe entries -- degradation is explicit).
    policy='error' : raise FloatLaneEscalation instead of escalating.

    Returns dict: logP (list), P (list; 0.0 where underflowed), route
    ('float'|'mp'), escalated (bool), triggers (dict), err_log (list, float
    route only), p_underflow (list of n with T3 marker)."""
    try:
        raw = logpmf_float_raw(a, b, r, nmax)
        t2 = None
    except FloatLaneEscalation as e:
        raw, t2 = None, str(e)

    triggers = {"t1_condition": bool(raw is not None and raw["flags_t1"].any()),
                "t2_series": t2 is not None,
                "t4_nonfinite": t2 is not None and t2.startswith("T4"),
                "t1_count": int(raw["flags_t1"].sum()) if raw is not None else None,
                "float_tol": FLOAT_TOL}
    need_escalation = triggers["t1_condition"] or triggers["t2_series"]

    if need_escalation and policy == "error":
        raise FloatLaneEscalation("float-lane triggers fired: %s" % json.dumps(triggers))

    if need_escalation and policy == "auto":
        P, rep = pmf_mp_selfcheck(a, b, r, nmax, dps=30)
        with mpmath.mp.workdps(40):
            logP = [float(mpmath.log(p)) for p in P]
            Pf = [float(p) for p in P]
        return {"logP": logP, "P": Pf, "route": "mp", "escalated": True,
                "triggers": triggers, "mp_selfcheck": rep,
                "p_underflow": [i for i, lp in enumerate(logP) if lp < FLOAT_UNDERFLOW_LOG]}

    if raw is None:   # policy == 'off' but the series itself failed: fail closed
        raise FloatLaneEscalation(t2)

    logP = raw["logP"]
    return {"logP": logP.tolist(), "P": np.exp(logP).tolist(), "route": "float",
            "escalated": False, "triggers": triggers,
            "err_log": raw["err_log"].tolist(),
            "unsafe_entries": np.nonzero(raw["flags_t1"])[0].tolist(),
            "p_underflow": np.nonzero(raw["flags_t3"])[0].tolist()}


def loglik_float(counts, a, b, r, policy="auto"):
    """Float fast-lane log-likelihood with escalation policy.

    Error bound: ERR_LL = sum_i err_log(n_i) + eps * sum_i |log P(n_i)|.
    Escalates (policy 'auto' -> mp lane) if ERR_LL > LOGLIK_FLOAT_TOL or any
    used entry carries a T1 flag.  Returns dict."""
    counts = _check_counts(counts)
    nmax = max(counts)
    mult = {}
    for c in counts:
        mult[c] = mult.get(c, 0) + 1
    ns = np.array(sorted(mult), dtype=int)
    ms = np.array([mult[n] for n in ns], dtype=float)

    try:
        raw = logpmf_float_raw(a, b, r, nmax)
        lp = raw["logP"][ns]
        el = raw["err_log"][ns]
        err_ll = float(np.sum(ms * el) + _EPS64 * np.sum(ms * np.abs(lp)))
        flagged = bool(raw["flags_t1"][ns].any())
        ll = float(np.sum(ms * lp))
        need = flagged or err_ll > LOGLIK_FLOAT_TOL
        fail = None
    except FloatLaneEscalation as e:
        need, fail = True, str(e)
        ll = err_ll = None
        flagged = True

    if need and policy == "auto":
        total, rep = loglik_mp(counts, a, b, r, dps=30)
        with mpmath.mp.workdps(40):
            return {"loglik": float(total), "route": "mp", "escalated": True,
                    "mp_selfcheck": rep,
                    "trigger_reason": fail or ("t1_flag" if flagged else "err_ll>tol"),
                    "err_bound_nats": 10.0 ** (-rep["min_agreement_digits"]) * abs(float(total))}
    if need and policy == "error":
        raise FloatLaneEscalation(fail or "loglik float-lane trigger")
    if fail:
        raise FloatLaneEscalation(fail)
    return {"loglik": ll, "route": "float", "escalated": False,
            "err_bound_nats": err_ll, "any_t1_flag": flagged}


# ----------------------------------------------------------------------------
# exact-rational self-test lane (integer a, b; rational r) + mutation control
# ----------------------------------------------------------------------------

def _check_counts(counts):
    counts = [int(c) for c in counts]
    if any(c < 0 for c in counts) or not counts:
        raise ValueError("counts must be a non-empty list of non-negative integers")
    return counts


def _exact_AB(a, b, r_frac, n):
    """For INTEGER a, b >= 1 and rational r: exact rationals (A_P, B_P) with
    P(n) = A_P + B_P e^{-r}, via sympy on the Euler integral
    I_n = int_0^1 u^{a+n-1} (1-u)^{b-1} e^{-r u} du = A + B e^{-r} (A, B rational
    in r), and P(n) = r^n/(n! B(a,b)) I_n.  Returns sympy Rationals."""
    import sympy as sp
    rs = sp.Rational(r_frac.numerator, r_frac.denominator)
    u = sp.Symbol("u", positive=True)
    R = sp.Symbol("R", positive=True)      # keep r symbolic so exp(-R) stays atomic
    I_R = sp.integrate(u ** (a + n - 1) * (1 - u) ** (b - 1) * sp.exp(-R * u), (u, 0, 1))
    E_R = sp.exp(-R)
    p = sp.Poly(sp.expand(I_R), E_R)
    if p.degree() > 1:
        raise RuntimeError("selftest: I_n not linear in exp(-r)")
    A_ = sp.Rational(sp.nsimplify(sp.simplify(p.coeff_monomial(1)).subs(R, rs)))
    B_ = sp.Rational(sp.nsimplify(sp.simplify(p.coeff_monomial(E_R)).subs(R, rs)))
    Bab = sp.Rational(sp.gamma(a) * sp.gamma(b) / sp.gamma(a + b))
    pref = rs ** n / (sp.factorial(n) * Bab)
    return sp.Rational(pref * A_), sp.Rational(pref * B_)


def _selftest_case(a, b, r_frac, n, prec=220, mutate=False):
    """One exact-rational containment test.  If mutate=True, deliberately corrupt
    the dominant exact coefficient by a relative 1e-30 bump; the certified ball
    (rel width ~1e-60 at prec 220) MUST then fail containment -- proving the
    test has teeth (fail-closed mutation control)."""
    A_P, B_P = _exact_AB(a, b, r_frac, n)
    ball = pmf_certified_single(a, b, r_frac, n, prec)
    with _flint_ctx.workprec(3 * prec + 60):
        e_mr = (-(arb(r_frac.numerator) / arb(r_frac.denominator))).exp()
        A_arb = arb(A_P.p) / arb(A_P.q)
        B_arb = arb(B_P.p) / arb(B_P.q)
        if mutate:
            # corrupt the coefficient with the LARGER contribution to P(n)
            bump = 1 + arb(1) / arb(10) ** 30
            if abs(A_arb) >= abs(B_arb * e_mr):
                A_arb = A_arb * bump
            else:
                B_arb = B_arb * bump
        exact_ball = A_arb + B_arb * e_mr
        diff = ball - exact_ball
        contains = bool(abs(diff.mid()) <= diff.rad())
    # mp and float routes vs exact (50-digit gate on mp; float to its own tol)
    with mpmath.mp.workdps(200):
        exact_mp = (mpmath.mpf(A_P.p) / A_P.q
                    + mpmath.mpf(B_P.p) / B_P.q * mpmath.exp(
                        -mpmath.mpf(r_frac.numerator) / r_frac.denominator))
        mp_val = pmf_mp_single(a, b, r_frac, n, 60)
        mp_digits = float(-mpmath.log10(abs(mp_val - exact_mp) / exact_mp)) \
            if mp_val != exact_mp else 60.0
        fl = logpmf_float_raw(float(a), float(b),
                              r_frac.numerator / r_frac.denominator, n)
        fdev = float(abs(mpmath.exp(mpmath.mpf(fl["logP"][n])) - exact_mp) / exact_mp)
        f_ok = fdev <= max(FLOAT_TOL, float(fl["err_log"][n]))
    return {"a": a, "b": b, "r": str(r_frac), "n": n,
            "A_P": str(A_P), "B_P": str(B_P),
            "certified_ball_contains_exact": contains,
            "mp_digits_vs_exact": round(mp_digits, 1),
            "float_dev_vs_exact": fdev, "float_within_model": f_ok,
            "mutated": mutate}


def selftest():
    """--selftest: exact-rational lane + mutation control.

    PASS requires (i) every unmutated case: certified ball CONTAINS the exact
    value, mp route >= 50 digits, float route within its error model; and
    (ii) the MUTATED control case is CAUGHT (containment fails).  Exit code 0
    only on full pass."""
    cases = [(2, 3, Fraction(7, 2), 5), (1, 4, Fraction(3), 10),
             (3, 2, Fraction(5, 4), 7), (5, 1, Fraction(10), 12),
             (4, 6, Fraction(25, 2), 20)]
    out = {"cases": [], "statement": "P(n) = A_P + B_P e^{-r}, A_P,B_P rational "
                                     "(integer a,b; rational r); certified ball must contain exact"}
    all_ok = True
    for (a, b, r, n) in cases:
        c = _selftest_case(a, b, r, n, mutate=False)
        ok = c["certified_ball_contains_exact"] and c["mp_digits_vs_exact"] >= 50 \
            and c["float_within_model"]
        c["pass"] = bool(ok)
        all_ok &= ok
        out["cases"].append(c)
    # mutation control (fail-closed proof): corrupt one coefficient, must be CAUGHT
    m = _selftest_case(2, 3, Fraction(7, 2), 5, mutate=True)
    mutation_caught = not m["certified_ball_contains_exact"]
    out["mutation_control"] = {"case": m, "caught": bool(mutation_caught),
                               "bump": "dominant coefficient * (1 + 1e-30)"}
    out["pass"] = bool(all_ok and mutation_caught)
    return out


# ----------------------------------------------------------------------------
# version info / provenance
# ----------------------------------------------------------------------------

# Core source sha-pin; for an on-site cross-check point
# TELEGRAPH_PILOT_CORE at your copy and version_info() re-hashes it.
_PILOT_CORE = os.environ.get("TELEGRAPH_PILOT_CORE", "")
_PILOT_CORE_SHA256 = "8e77904db5ae4a81b0c7835c248e90749bcedf6627ed2020871fc26d7c42d78b"


def version_info():
    import platform
    import flint as _flint
    import sympy as _sympy
    info = {
        "tool": "telegraph_evaluate.py",
        "version": __version__,
        "convention": "Beta-mixture (a,b,r) = (k_on,k_off,k_syn)/delta; "
                      "P(n) = int_0^1 Pois(n;ru) Beta(u;a,b) du",
        "python": platform.python_version(),
        "numpy": np.__version__,
        "mpmath": mpmath.__version__,
        "sympy": _sympy.__version__,
        "python-flint": getattr(_flint, "__version__", "unknown"),
        "self_sha256": _sha256_of(os.path.abspath(__file__)),
        "core": "telegraph core",
        "pilot_core_sha256": (_sha256_of(_PILOT_CORE) if _PILOT_CORE
                              else _PILOT_CORE_SHA256),
        "float_policy": {"FLOAT_TOL": FLOAT_TOL, "LOGLIK_FLOAT_TOL": LOGLIK_FLOAT_TOL,
                         "FLOAT_UNDERFLOW_LOG": FLOAT_UNDERFLOW_LOG},
    }
    return info


def _sha256_of(path):
    try:
        with open(path, "rb") as f:
            return hashlib.sha256(f.read()).hexdigest()
    except OSError:
        return "unavailable"


# ----------------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------------

def _ball_strings(ball, digits):
    """Decimal mid/rad strings of an arb ball at the requested digits."""
    with mpmath.mp.workdps(digits + 15):
        mid = arb_mid_mpf(ball)
        rad = arb_mid_mpf(ball.rad() + arb(0))
        ms = mpmath.nstr(mid, digits, strip_zeros=False)
        rs = mpmath.nstr(rad, 3)
        rel = float(rad / mid) if mid != 0 else float("inf")
    return ms, rs, rel


def _parse_counts(args):
    if args.counts is not None:
        return [int(x) for x in args.counts.replace(",", " ").split()]
    if args.counts_file is not None:
        with open(args.counts_file) as f:
            # same tokenization as --counts: commas or whitespace
            return [int(x) for x in f.read().replace(",", " ").split()]
    raise SystemExit("loglik requires --counts or --counts-file")


def _resolve_params(args):
    """(a, b, r) directly, or (kon, koff, ksyn, delta) converted + documented."""
    direct = args.a is not None or args.b is not None or args.r is not None
    rates = args.kon is not None or args.koff is not None or args.ksyn is not None
    if direct and rates:
        raise SystemExit("give EITHER --a/--b/--r OR --kon/--koff/--ksyn[/--delta], not both")
    if direct:
        if None in (args.a, args.b, args.r):
            raise SystemExit("need all of --a --b --r")
        return float(args.a), float(args.b), float(args.r), {
            "input": "primary (a,b,r)", "a": args.a, "b": args.b, "r": args.r}
    if rates:
        if None in (args.kon, args.koff, args.ksyn):
            raise SystemExit("need all of --kon --koff --ksyn (and optional --delta)")
        delta = args.delta if args.delta is not None else 1.0
        a, b, r = abr_from_rates(args.kon, args.koff, args.ksyn, delta)
        return a, b, r, {
            "input": "rates (k_on,k_off,k_syn,delta)",
            "kon": args.kon, "koff": args.koff, "ksyn": args.ksyn, "delta": delta,
            "mapping": "a=k_on/delta, b=k_off/delta, r=k_syn/delta (verified dictionary, "
                       "see module docstring / pilot CONVENTIONS.md)",
            "a": a, "b": b, "r": r}
    raise SystemExit("give --a/--b/--r or --kon/--koff/--ksyn")


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    # top-level flags per spec
    if "--selftest" in argv:
        out = selftest()
        print(json.dumps(out, indent=1, sort_keys=True, default=str))
        return 0 if out["pass"] else 1
    if "--version-info" in argv:
        print(json.dumps(version_info(), indent=1, sort_keys=True))
        return 0

    ap = argparse.ArgumentParser(
        description="Certified telegraph-model PMF / log-likelihood evaluator "
                    "(Beta-mixture convention (a,b,r) = (k_on,k_off,k_syn)/delta).")
    ap.add_argument("task", choices=["pmf", "loglik"])
    ap.add_argument("--a", type=float); ap.add_argument("--b", type=float)
    ap.add_argument("--r", type=float)
    ap.add_argument("--kon", type=float); ap.add_argument("--koff", type=float)
    ap.add_argument("--ksyn", type=float); ap.add_argument("--delta", type=float)
    ap.add_argument("--nmax", type=int, default=None)
    ap.add_argument("--n", type=int, default=None, help="single-n PMF evaluation")
    ap.add_argument("--mode", choices=["certified", "mp", "float"], default="certified")
    ap.add_argument("--prec", type=int, default=150, help="certified mode: bits")
    ap.add_argument("--dps", type=int, default=30, help="mp mode: target digits")
    ap.add_argument("--digits", type=int, default=30, help="output decimal digits")
    ap.add_argument("--policy", choices=["auto", "off", "error"], default="auto",
                    help="float-lane escalation policy")
    ap.add_argument("--floats", action="store_true",
                    help="certified/mp modes: also emit float64 midpoints")
    ap.add_argument("--counts", type=str, default=None)
    ap.add_argument("--counts-file", type=str, default=None)
    args = ap.parse_args(argv)

    a, b, r, param_block = _resolve_params(args)
    out = {"convention": "Beta-mixture (a,b,r)=(k_on,k_off,k_syn)/delta; "
                         "P(n)=int_0^1 Pois(n;ru)Beta(u;a,b)du",
           "params": param_block, "mode": args.mode, "version": __version__}

    try:
        return _run_task(args, a, b, r, out)
    except FloatLaneEscalation as e:
        print(json.dumps({"error": "float_lane_escalation", "detail": str(e),
                          "params": param_block, "mode": args.mode},
                         indent=1, sort_keys=True))
        return 2


def _run_task(args, a, b, r, out):

    if args.task == "pmf":
        if args.nmax is None and args.n is None:
            raise SystemExit("pmf requires --nmax (vector) or --n (single)")
        if args.mode == "certified":
            if args.n is not None and args.nmax is None:
                ball = pmf_certified_single(a, b, r, args.n, args.prec)
                ms, rs, rel = _ball_strings(ball, args.digits)
                out["P"] = [{"n": args.n, "mid": ms, "rad": rs, "rel_rad": rel}]
            else:
                P = pmf_certified(a, b, r, args.nmax, args.prec)
                out["P"] = []
                for n, ball in enumerate(P):
                    ms, rs, rel = _ball_strings(ball, args.digits)
                    ent = {"n": n, "mid": ms, "rad": rs, "rel_rad": rel}
                    if args.floats:
                        ent["float"] = float(ball.mid())
                    out["P"].append(ent)
            out["prec_bits"] = args.prec
            out["certified"] = True
        elif args.mode == "mp":
            nmax = args.nmax if args.nmax is not None else args.n
            P, rep = pmf_mp_selfcheck(a, b, r, nmax, dps=args.dps)
            with mpmath.mp.workdps(args.dps + 10):
                sel = range(nmax + 1) if args.nmax is not None else [args.n]
                out["P"] = [{"n": n, "value": mpmath.nstr(P[n], args.digits)} for n in sel]
            out["selfcheck"] = rep
        else:
            nmax = args.nmax if args.nmax is not None else args.n
            res = pmf_float(a, b, r, nmax, policy=args.policy)
            if args.n is not None and args.nmax is None:
                res["logP"] = [res["logP"][args.n]]
                res["P"] = [res["P"][args.n]]
            out.update(res)
    else:  # loglik
        counts = _parse_counts(args)
        out["n_obs"] = len(counts)
        if args.mode == "certified":
            ball, nmax_used = loglik_certified(counts, a, b, r, args.prec)
            ms, rs, rel = _ball_strings(ball, args.digits)
            out.update({"loglik_mid": ms, "loglik_rad": rs,
                        "certified_error_bound_nats": rs, "nmax_used": nmax_used,
                        "prec_bits": args.prec, "certified": True})
        elif args.mode == "mp":
            total, rep = loglik_mp(counts, a, b, r, dps=args.dps)
            with mpmath.mp.workdps(args.dps + 10):
                out["loglik"] = mpmath.nstr(total, args.digits)
            out["selfcheck"] = rep
        else:
            out.update(loglik_float(counts, a, b, r, policy=args.policy))

    print(json.dumps(out, indent=1, sort_keys=True, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())

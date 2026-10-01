# baller engine — certified arb-ball enclosure primitives.
"""Certified enclosure primitives for the matPTF alert computation.

Ball-arithmetic (python-flint arb 0.9.0) enclosures of a float64 reference
implementation of matPTF's hazard/alert/conventions layers (matPTF commit
fe4dd95).

BINDING LAWS
------------
(1) All numeric constants are constructed AT WORKING PRECISION inside a prec
    context (``working_prec``); never through a 53-bit float when the intended
    real differs from the float64 value (e.g. arb('0.1'), parsed at working
    precision, when the convention is DECIMAL 0.1).  Where the matPTF
    convention IS the float64 double (the 65-threshold grid, the shipped
    scenario weights, and — per the Matlab literal comparison — the 0.10 alert
    boundary), inputs are float64-PINNED: ``arb(float(v))`` is exact, no
    rounding.  The alert boundary carries an explicit flag
    ``boundary_semantics='float64'|'decimal'`` (default ``'float64'``:
    PTF_AlertLevels.m:26/139 compares against the compiled literal 0.10, i.e.
    the double nearest decimal 0.1 = 0.1000000000000000055511151231257827...).
(2) Fail-closed: any non-finite / blown enclosure returns INDET — the whole
    extended real line as a ball for the arb-valued functions (a still-valid
    enclosure that can never certify a pass), the string markers
    'INDET' / status 'indet' for the discrete-valued ones.  A ball that makes
    a convention path (dedup / branch) ambiguous is never resolved silently:
    the result is the UNION over all convention-consistent paths with
    status 'ambiguous-path'.
(3) Every numeric result is an arb ball, never a float.

Comparison semantics used throughout (python-flint): ``a < b`` / ``a <= b``
on arb balls are TRUE only when the relation holds for every pair of points
in the two balls ("certainly"); "possibly a < b" is ``not (b <= a)``.
"""

from contextlib import contextmanager

from flint import arb, ctx

DEFAULT_PREC = 128
BOUNDARY_SEMANTICS_DEFAULT = "float64"

ALERT_STRADDLES = "STRADDLES"
ALERT_INDET = "INDET"

# statuses returned by cert_percentile_bracket
STATUS_OK = "ok"
STATUS_GT100M_BRANCH = "gt100m-branch"
STATUS_AMBIGUOUS = "ambiguous-path"
STATUS_INDET = "indet"


@contextmanager
def working_prec(prec=DEFAULT_PREC):
    """Set the global arb working precision (bits), restoring on exit.

    Law (1): every certified primitive wraps its whole body — including all
    constant construction (sqrt(2), parsed decimal strings, ...) — in this
    context, so nothing is ever built at an unintended precision.
    """
    old = ctx.prec
    ctx.prec = int(prec)
    try:
        yield
    finally:
        ctx.prec = old


def INDET():
    """The fail-closed indeterminate ball: the whole extended real line.

    Still a mathematically valid enclosure of any real value, but it can never
    certify anything (law 2).  Detect with ``is_indet``.
    """
    return arb(0, float("inf"))


def is_indet(x):
    """True iff x is not a usable enclosure (nan, or non-finite mid/radius)."""
    if not isinstance(x, arb):
        return True
    return x.is_nan() or (not x.is_finite())


def pin(v):
    """Convert an input to arb under the pinning conventions.

    - arb: passed through unchanged (upstream enclosures compose here);
    - str: parsed AT WORKING PRECISION (small ball around the decimal value)
      — use for decimal-convention constants, per law (1);
    - int/float/numpy scalar: float64-PINNED, i.e. the exact double value
      (arb from a float is exact) — use for shipped float64 data (threshold
      grid, weights, Matlab-literal boundaries under 'float64' semantics).
    """
    if isinstance(v, arb):
        return v
    if isinstance(v, str):
        return arb(v)
    if isinstance(v, int):
        # exact integer pin — float(v) would silently round ints > 2^53
        # (documented minor finding); arb handles big ints exactly
        return arb(v)
    return arb(float(v))


def _is_point_neg_inf(m):
    """Exactly the point -inf (the mu = log(0) = -inf exclusion convention)."""
    return (not m.is_nan()) and m.is_exact() and (not m.is_finite()) \
        and bool(m < arb(0))


def _union_all(balls):
    out = balls[0]
    for b in balls[1:]:
        out = out.union(b)
    return out


# ---------------------------------------------------------------------------
# 1. Lognormal survival
# ---------------------------------------------------------------------------

def cert_lognormal_survival(x, mu, sigma=1.0, prec=DEFAULT_PREC):
    """Certified 1 - logncdf(x; mu, sigma) = 0.5*erfc((ln x - mu)/(sigma*sqrt 2)).

    Encloses the survival used at PTF_1_0.m:782-793 (the reference
    implementation's lognormal_survival).  Returns an arb ball (law 3).

    Conventions:
    - mu = -inf (exact point) is the glVal == 0 exclusion convention
      (PTF_1_0.m:768,776-778 via log(0) = -inf): returns EXACTLY arb(0).
    - x is a 65-threshold-grid value: pass the float64 — it is pinned exactly.
    - sigma: LOGNORMAL_SIGMA = 1 hardcoded in matPTF (conventions.py:37);
      exposed for sensitivity studies.

    Fail-closed (law 2): nan anywhere, x not certainly finite-positive, sigma
    not certainly finite-positive, mu non-finite other than the exact -inf
    point, or a blown enclosure -> INDET().  The result is intersected with
    [0, 1] (survival is a probability) — a sound tightener.
    """
    with working_prec(prec):
        xb, mb, sb = pin(x), pin(mu), pin(sigma)
        if xb.is_nan() or mb.is_nan() or sb.is_nan():
            return INDET()
        if not (xb.is_finite() and arb(0) < xb):
            return INDET()
        if not (sb.is_finite() and arb(0) < sb):
            return INDET()
        if _is_point_neg_inf(mb):
            return arb(0)                      # exact: exclusion convention
        if not mb.is_finite():
            return INDET()
        t = (xb.log() - mb) / (sb * arb(2).sqrt())
        res = arb(0.5) * t.erfc()
        if res.is_nan() or (not res.is_finite()):
            return INDET()
        try:
            res = res.intersection(arb(0.5, 0.5))     # [0, 1]
        except ValueError:                            # disjoint: soundness lost
            return INDET()
        return res


# ---------------------------------------------------------------------------
# 2. Weighted survival sum (one hazard-curve entry)
# ---------------------------------------------------------------------------

def cert_weighted_survival_sum(probs, mus, x, sigma=1.0, prec=DEFAULT_PREC):
    """Certified sum_s p_s * survival(x; mu_s, sigma) — one hc(ip, k) entry
    (PTF_1_0.m:793).  Returns an arb ball.

    probs are the SHIPPED float64 weights: they are pinned exactly (arb from
    float is exact), so this enclosure covers ONLY the evaluation error of
    THIS sum — log/erfc/mul/add rounding at working precision.  Enclosures of
    the weights themselves (scenario-probability uncertainty, upstream
    normalization error) compose upstream later: pass arb balls as probs and
    they propagate through unchanged (pin() passes arb through).

    Fail-closed: any nan/non-finite weight, any INDET survival term, or a
    blown running sum -> INDET().
    """
    if len(probs) != len(mus):
        raise ValueError("probs and mus must have the same length")
    with working_prec(prec):
        total = arb(0)
        for p, mu in zip(probs, mus):
            pb = pin(p)
            if pb.is_nan() or (not pb.is_finite()):
                return INDET()
            s = cert_lognormal_survival(x, mu, sigma, prec=prec)
            if is_indet(s):
                return INDET()
            total = total + pb * s
        if total.is_nan() or (not total.is_finite()):
            return INDET()
        return total


# ---------------------------------------------------------------------------
# 3. Normalized weighted sum
# ---------------------------------------------------------------------------

def cert_normalized_weighted_sum(numer_terms, denom_terms, values=None,
                                 prec=DEFAULT_PREC):
    """Certified (sum_i a_i) / (sum_i b_i): single ball sums + one division.

    Encloses the weight normalizations and the Average
    method's mean (weighted mean intensity, PTF_AlertLevels.m:131 path).
    Returns an arb ball.

    Convex-combination tightener (applied when applicable): if the caller
    asserts a_i = w_i * f_i with w_i = denom_terms[i] >= 0 by passing the f_i
    as ``values``, and all w_i are certainly >= 0 with sum certainly > 0, then
    the true value is a convex combination of the f_i and hence lies in
    [min_i f_i, max_i f_i]; the division result is INTERSECTED with the hull
    of the f_i balls.  The assertion is checked: each numer_terms[i] must
    overlap w_i * f_i (recomputed in ball arithmetic); a violation raises
    ValueError (caller misuse, not a numeric INDET).  If applicability cannot
    be certified (some w_i possibly negative), the tightener is skipped —
    never applied on an uncertain premise.

    Fail-closed: nan/non-finite terms, a denominator ball containing 0, a
    blown quotient, or an empty tightener intersection -> INDET().
    """
    if len(numer_terms) != len(denom_terms):
        raise ValueError("numer_terms and denom_terms must have the same length")
    if len(denom_terms) == 0:
        raise ValueError("empty term lists")
    with working_prec(prec):
        a = [pin(t) for t in numer_terms]
        b = [pin(t) for t in denom_terms]
        for t in a + b:
            if t.is_nan() or (not t.is_finite()):
                return INDET()
        sa = arb(0)
        for t in a:
            sa = sa + t
        sb = arb(0)
        for t in b:
            sb = sb + t
        if not (arb(0) < sb or sb < arb(0)):   # 0 possibly in denominator
            return INDET()
        r = sa / sb
        if r.is_nan() or (not r.is_finite()):
            return INDET()
        if values is not None:
            if len(values) != len(b):
                raise ValueError("values must match denom_terms in length")
            f = [pin(v) for v in values]
            for v in f:
                if v.is_nan() or (not v.is_finite()):
                    return INDET()
            # consistency check of the caller's a_i = w_i * f_i assertion
            for ai, wi, fi in zip(a, b, f):
                if not ai.overlaps(wi * fi):
                    raise ValueError(
                        "numer_terms[i] inconsistent with denom_terms[i]*values[i]"
                        " — convex-combination tightener premise violated")
            applicable = all(bool(arb(0) <= w) for w in b) and bool(arb(0) < sb)
            if applicable:
                hull = _union_all(f)
                try:
                    r = r.intersection(hull)
                except ValueError:             # disjoint: soundness lost
                    return INDET()
        return r


# ---------------------------------------------------------------------------
# 4. Percentile bracket (shipped interp1 convention)
# ---------------------------------------------------------------------------

def _interp_path(xxu, yyu, prb):
    """Evaluate the shipped interp1 convention on ONE deduped path.

    Returns (value_ball, gt_possible, gt_certain).  Segment selection is by
    candidacy (a segment contributes iff prb possibly lies in it), each
    contribution intersected with the segment's y-hull (interp at an interior
    query is a convex combination of the endpoint intensities); the union of
    candidates is a sound bracket even when hc balls overlap.
    """
    ZERO = arb(0)
    if len(xxu) < 2:
        # dedup collapsed everything onto the anchor: Matlab interp1 errors
        # ("at least two sample points") — edge geometry of the reference
        # code (index past the last threshold).
        return INDET(), False, False
    last = xxu[-1]
    gt_possible = not (last <= prb)           # possibly prb < xx(end)
    gt_certain = bool(prb < last)             # certainly prb < xx(end)
    cands = []
    if gt_possible:
        # PTF_AlertLevels.m:118-121: when prThr < xx(end) the reference code
        # assigns xx(end) (an exceedance probability) to valTmp and prints a
        # warning; replicated as-is and flagged so callers can see when it fires.
        cands.append(last)
    if not gt_certain:
        if not (prb <= xxu[0]):
            # possibly prb > anchor 1: interp1 out of range -> NaN in Matlab
            cands.append(INDET())
        for i in range(len(xxu) - 1):
            hi, lo = xxu[i], xxu[i + 1]
            # candidate iff possibly lo <= prb <= hi
            if (not (hi < prb)) and (not (prb < lo)):
                yint = yyu[i].union(yyu[i + 1])
                denom = hi - lo
                val = yint
                if denom.is_finite() and (ZERO < denom or denom < ZERO):
                    cand = yyu[i] + (yyu[i + 1] - yyu[i]) * ((hi - prb) / denom)
                    if cand.is_finite() and not cand.is_nan():
                        try:
                            val = cand.intersection(yint)
                        except ValueError:
                            val = yint          # keep the sound superset
                cands.append(val)
    if not cands:
        return INDET(), gt_possible, gt_certain
    return _union_all(cands), gt_possible, gt_certain


def cert_percentile_bracket(hc_balls, thresholds, pr, prec=DEFAULT_PREC,
                            max_paths=64):
    """Certified bracket for the shipped percentile convention
    (PTF_AlertLevels.m:104-121, percentile_value).

    Convention evaluated in ball arithmetic:
      nonNull = indices with hc > 0; if zeros exist and nonNull is nonempty,
      append nonNull[-1]+1 (raise geometry -> INDET, see below); anchor (1, 0)
      prepended; CONSECUTIVE-duplicate dedup (element kept iff != its
      predecessor in the ORIGINAL sequence, first of each run survives);
      linear interp1 in (exceedance prob, intensity) at pr; if pr < xx(end)
      the shipped >100m branch assigns xx(end) — a probability — as the value
      (status 'gt100m-branch').

    Monotone-envelope / path-union bracket: every discrete decision the balls
    cannot certify (hc_i > 0 vs == 0; consecutive equality for dedup; the
    >100m branch) is BRANCHED, the value is the UNION over all
    convention-consistent paths, and status = 'ambiguous-path' — never a
    silent pick (law 2).  If the path count would exceed max_paths, a global
    sound hull [0, max(1, thresholds, hc)] is returned with the same status.

    Returns (arb_ball, status) with status in
      'ok' | 'gt100m-branch' | 'ambiguous-path' | 'indet'.
    'indet' covers nan/non-finite inputs, certainly-negative hc (invalid
    geometry), the index-past-end edge case of the reference code (interior-zero curve
    with hc[end] > 0, PTF_AlertLevels.m:108), and dedup collapse to a single
    point (interp1 error geometry).

    All-zero curve returns (arb(0), 'ok') — the shipped valTmp = 0 (lines
    110-111); the CALLER maps a zero CURVE to alert level 1 (line 143) before
    ever calling this (PTF convention, alert layer line 105).
    """
    with working_prec(prec):
        hc = [pin(h) for h in hc_balls]
        th = [pin(t) for t in thresholds]
        n = len(hc)
        if len(th) != n:
            raise ValueError("hc_balls and thresholds must have the same length")
        if n == 0:
            raise ValueError("empty grid")
        prb = pin(pr)
        bad = prb.is_nan() or (not prb.is_finite()) or (not (arb(0) < prb))
        for v in hc + th:
            if v.is_nan() or (not v.is_finite()):
                bad = True
        if bad:
            return INDET(), STATUS_INDET

        ONE = arb(1)
        ZERO = arb(0)

        # sign classification of each hc entry (the nonNull decision)
        classes = []
        for h in hc:
            if h.is_zero():
                classes.append("zero")
            elif ZERO < h:
                classes.append("pos")
            elif h < ZERO:
                return INDET(), STATUS_INDET   # certainly negative: invalid
            else:
                classes.append("amb")
        amb = [i for i, c in enumerate(classes) if c == "amb"]

        def _fallback():
            hull = ZERO.union(ONE)
            for t in th:
                hull = hull.union(t)
            for h in hc:
                hull = hull.union(h.nonnegative_part())
            return hull, STATUS_AMBIGUOUS

        if 2 ** len(amb) > max_paths:
            return _fallback()

        ambiguous_any = bool(amb)
        path_vals = []
        gt_any = False
        gt_certain_all = True
        n_sign_paths = 2 ** len(amb)

        for mask in range(n_sign_paths):
            assign = {i: bool((mask >> k) & 1) for k, i in enumerate(amb)}
            non_null = [i for i in range(n)
                        if classes[i] == "pos"
                        or (classes[i] == "amb" and assign[i])]
            if not non_null:
                path_vals.append(ZERO)         # shipped lines 110-111
                gt_certain_all = False
                continue

            # xx / yy sequence: anchor + positive entries (+ appended zero)
            seq = [(ONE, ZERO, "anchor")]
            for i in non_null:
                v = hc[i].nonnegative_part() if classes[i] == "amb" else hc[i]
                seq.append((v, th[i], "pos"))
            if len(non_null) < n:              # zeros exist -> append (107-109)
                nxt = non_null[-1] + 1
                if nxt >= n:
                    # index-past-end edge case of the reference code (interior-zero + hc[end]>0)
                    return INDET(), STATUS_INDET
                seq.append((ZERO, th[nxt], "zero"))

            # dedup decisions between consecutive ORIGINAL entries
            decs = []
            for j in range(1, len(seq)):
                a, ka = seq[j - 1][0], seq[j - 1][2]
                b, kb = seq[j][0], seq[j][2]
                if kb == "zero" and ka == "pos":
                    d = "ne"        # entry is > 0 on this path, appended is 0
                elif a == b:
                    d = "eq"        # certain (both exact, equal)
                elif (a < b) or (b < a):
                    d = "ne"        # certain (disjoint balls)
                else:
                    d = "amb"
                decs.append(d)
            amb_pairs = [j for j, d in enumerate(decs) if d == "amb"]
            if n_sign_paths * (2 ** len(amb_pairs)) > max_paths:
                return _fallback()
            if amb_pairs:
                ambiguous_any = True

            for dmask in range(2 ** len(amb_pairs)):
                dd = list(decs)
                for k, j in enumerate(amb_pairs):
                    dd[j] = "ne" if (dmask >> k) & 1 else "eq"
                keep = [0] + [j + 1 for j, d in enumerate(dd) if d == "ne"]
                xxu = [seq[j][0] for j in keep]
                yyu = [seq[j][1] for j in keep]
                val, gt_possible, gt_certain = _interp_path(xxu, yyu, prb)
                path_vals.append(val)
                gt_any = gt_any or gt_possible
                gt_certain_all = gt_certain_all and gt_certain

        overall = _union_all(path_vals)
        if overall.is_nan() or (not overall.is_finite()):
            return INDET(), STATUS_INDET
        if ambiguous_any or len(path_vals) > 1 or \
                (gt_any and not gt_certain_all):
            return overall, STATUS_AMBIGUOUS
        if gt_any and gt_certain_all:
            return overall, STATUS_GT100M_BRANCH
        return overall, STATUS_OK


# ---------------------------------------------------------------------------
# 5. Alert level
# ---------------------------------------------------------------------------

def _boundary_ball(b, semantics):
    """One alert boundary under the chosen semantics.

    'float64': the boundary is the exact double — matPTF's actual boundary,
      since Matlab compares valTmp <= 0.10 against the compiled float64
      literal (PTF_AlertLevels.m:26,139).  arb(float) is exact.
    'decimal': the boundary is the DECIMAL value, parsed at working precision
      (a ball of radius ~2^-prec around the exact decimal).  Numeric bounds
      are serialized via repr(float(b)) — the shortest decimal that round-trips
      to the double, i.e. '0.1' for 0.10; pass a str for full control.
    """
    if semantics == "float64":
        if isinstance(b, str):
            b = float(b)
        return arb(float(b))
    if semantics == "decimal":
        if isinstance(b, str):
            return arb(b)
        return arb(repr(float(b)))
    raise ValueError("boundary_semantics must be 'float64' or 'decimal'")


def cert_alert_level(val_ball, bounds=(0.0, 0.10, 0.5),
                     boundary_semantics=BOUNDARY_SEMANTICS_DEFAULT,
                     prec=DEFAULT_PREC):
    """Certified alert binning: find(val <= [bounds..., inf], 1) - 1
    (PTF_AlertLevels.m:139; conventions.ALERT_BOUNDS_M = (0, 0.10, 0.5)).

    Returns int level (0..len(bounds)), or:
      'STRADDLES' — the ball contains a boundary and cannot be certified to
        either side.  A DEGENERATE (exact) ball equal to a float64 boundary
        under 'float64' semantics takes the <= side (the shipped comparison
        is exact there — no straddle exists).
      'INDET'     — nan / non-finite input (law 2: never a pass).

    boundary_semantics (default 'float64'): see _boundary_ball.  Note that
    under 'decimal' semantics an exact float64 0.1 value is CERTAINLY ABOVE
    the decimal-0.1 boundary (the double exceeds decimal 0.1 by ~5.55e-18,
    far more than the boundary ball's ~2^-prec radius).
    """
    with working_prec(prec):
        v = pin(val_ball)
        if v.is_nan() or (not v.is_finite()):
            return ALERT_INDET
        for i, braw in enumerate(bounds):
            b = _boundary_ball(braw, boundary_semantics)
            if v <= b:                 # certainly on the <= side
                return i
            if not (b < v):            # not certainly above either
                return ALERT_STRADDLES
        return len(bounds)

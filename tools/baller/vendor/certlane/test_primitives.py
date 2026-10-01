# baller engine — containment tests for certlane/primitives.py.
"""Containment tests for certlane.primitives against an mpmath dps=60 oracle.

Run from the vendor directory (the parent of certlane/):
    python3 -m pytest certlane/test_primitives.py -q

Coverage per the certified-enclosure spec:
  - containment vs mpmath dps=60 oracle, >=200 random cases per primitive
    (log-uniform x, mu in [-8,3], probs on a random simplex);
  - planted-error controls (perturbed weight must escape the enclosure;
    last-bit-corrupted survival must fail the containment check);
  - exact-convention cases (mu=-inf; float64-0.1 boundary straddle rules);
  - INDET behavior (nan/inf inputs, zero denominators, crash geometries).
"""

import math
import random
import time

import mpmath as mp
import numpy as np
import pytest
from flint import arb

from certlane.primitives import (
    ALERT_INDET,
    ALERT_STRADDLES,
    STATUS_AMBIGUOUS,
    STATUS_GT100M_BRANCH,
    STATUS_INDET,
    STATUS_OK,
    INDET,
    cert_alert_level,
    cert_lognormal_survival,
    cert_normalized_weighted_sum,
    cert_percentile_bracket,
    cert_weighted_survival_sum,
    is_indet,
    pin,
    working_prec,
)

mp.mp.dps = 60

NEG_INF = float("-inf")
THRESHOLDS_65 = np.logspace(-5, 2, 65)          # 1e-5..100 m, float64-pinned grid


# ---------------------------------------------------------------------------
# oracle helpers
# ---------------------------------------------------------------------------

def _bound_to_mpf(b):
    """Exact conversion of an exact arb bound to mpf (via mantissa*2^exp)."""
    assert not b.is_nan()
    if not b.is_finite():
        return mp.inf if b > arb(0) else mp.ninf
    if b.is_zero():
        return mp.mpf(0)
    m, e = b.man_exp()
    return mp.mpf(int(m)) * mp.power(2, int(e))


def ball_bounds(ball):
    with working_prec(128):
        return _bound_to_mpf(ball.lower()), _bound_to_mpf(ball.upper())


def contains_mpf(ball, v):
    lo, hi = ball_bounds(ball)
    return lo <= v <= hi


def ball_width(ball):
    lo, hi = ball_bounds(ball)
    return hi - lo


def oracle_survival(x, mu, sigma=1.0):
    """0.5*erfc((ln x - mu)/(sigma*sqrt 2)) at dps=60; mu=-inf -> 0 exactly."""
    if mu == NEG_INF:
        return mp.mpf(0)
    t = (mp.log(mp.mpf(x)) - mp.mpf(mu)) / (mp.mpf(sigma) * mp.sqrt(2))
    return mp.erfc(t) / 2


def oracle_weighted_sum(probs, mus, x, sigma=1.0):
    return mp.fsum(mp.mpf(p) * oracle_survival(x, mu, sigma)
                   for p, mu in zip(probs, mus))


def oracle_percentile(hc_floats, thresholds, pr):
    """The shipped interp1 convention (matPTF percentile_value) at
    dps=60 on exact float64 inputs.  Returns (value, flag) with flag in
    {'', 'gt100m'}; returns (None, 'nonmonotone') when the deduped xx is not
    strictly decreasing (Matlab interp1 would error — skip such cases)."""
    n = len(hc_floats)
    nn = [i for i in range(n) if hc_floats[i] > 0.0]
    if n > len(nn) and nn:
        nxt = nn[-1] + 1
        assert nxt < n, "crash geometry in oracle input"
        nn.append(nxt)
    if not nn:
        return mp.mpf(0), ""
    xx = [mp.mpf(1)] + [mp.mpf(hc_floats[i]) for i in nn]
    yy = [mp.mpf(0)] + [mp.mpf(float(thresholds[i])) for i in nn]
    keep = [0] + [j for j in range(1, len(xx)) if xx[j] != xx[j - 1]]
    xxu = [xx[j] for j in keep]
    yyu = [yy[j] for j in keep]
    if any(xxu[j] <= xxu[j + 1] for j in range(len(xxu) - 1)):
        return None, "nonmonotone"
    prm = mp.mpf(pr)
    if prm < xxu[-1]:
        return xxu[-1], "gt100m"
    for i in range(len(xxu) - 1):
        if xxu[i] >= prm >= xxu[i + 1]:
            hi, lo = xxu[i], xxu[i + 1]
            return yyu[i] + (yyu[i + 1] - yyu[i]) * (hi - prm) / (hi - lo), ""
    raise AssertionError("no segment found (pr > 1?)")


def oracle_alert_level(v, bounds=(0.0, 0.10, 0.5)):
    """find(val <= [bounds..., inf], 1) - 1 on float64 (shipped convention)."""
    for i, b in enumerate(bounds):
        if v <= b:
            return i
    return len(bounds)


def float_survival(x, mu, sigma=1.0):
    """The float64 reference lognormal_survival, scalar."""
    from scipy.special import erfc as _erfc
    if mu == NEG_INF:
        return 0.0
    return 0.5 * float(_erfc((math.log(x) - mu) / (sigma * math.sqrt(2.0))))


def make_hc(rng, n_scen, thresholds, mu_lo=-7.0, mu_hi=2.0, neginf_frac=0.15):
    """Random simplex weights + mus -> float64 hazard curve on the grid."""
    g = np.array([rng.random() for _ in range(n_scen)])
    probs = (g / g.sum()).astype(np.float64)
    mus = np.array([NEG_INF if rng.random() < neginf_frac
                    else rng.uniform(mu_lo, mu_hi)
                    for _ in range(n_scen)])
    hc = np.zeros(len(thresholds))
    for k, x in enumerate(thresholds):
        hc[k] = sum(p * float_survival(float(x), m) for p, m in zip(probs, mus))
    return probs, mus, hc


# ---------------------------------------------------------------------------
# 1. cert_lognormal_survival
# ---------------------------------------------------------------------------

def test_survival_containment_random():
    rng = random.Random(20260717)
    n_checked = 0
    for _ in range(300):
        x = 10.0 ** rng.uniform(-5, 2)          # log-uniform over the grid range
        if rng.random() < 0.10:
            mu = NEG_INF
        else:
            mu = rng.uniform(-8.0, 3.0)
        sigma = 1.0 if rng.random() < 0.7 else rng.uniform(0.5, 2.0)
        ball = cert_lognormal_survival(x, mu, sigma)
        assert not is_indet(ball)
        if mu == NEG_INF:
            assert ball.is_zero()               # exact convention value
            continue
        v = oracle_survival(x, mu, sigma)
        assert contains_mpf(ball, v), (x, mu, sigma)
        assert ball_width(ball) < mp.mpf("1e-25")
        n_checked += 1
    assert n_checked >= 200


def test_survival_mu_neg_inf_exact():
    for x in (1e-5, 0.1, 1.0, 100.0):
        ball = cert_lognormal_survival(x, NEG_INF)
        assert ball.is_zero() and ball.is_exact()
        ball2 = cert_lognormal_survival(x, arb(NEG_INF))
        assert ball2.is_zero() and ball2.is_exact()


def test_survival_lastbit_corruption_caught():
    rng = random.Random(4242)
    n_checked = 0
    tries = 0
    while n_checked < 200 and tries < 2000:
        tries += 1
        x = 10.0 ** rng.uniform(-3, 1.5)
        mu = rng.uniform(-4.0, 2.0)
        v = oracle_survival(x, mu)
        if not (mp.mpf("1e-6") < v < 1 - mp.mpf("1e-6")):
            continue                             # keep ulp geometry sensible
        ball = cert_lognormal_survival(x, mu)
        assert contains_mpf(ball, v)
        vf = float(v)                            # correctly rounded double
        for corrupted in (np.nextafter(vf, np.inf), np.nextafter(vf, -np.inf)):
            # one-ulp corruption sits ~1e-17 from the true value; the certified
            # ball at prec=128 is ~1e-38 wide -> containment must fail.
            assert not contains_mpf(ball, mp.mpf(float(corrupted)))
        n_checked += 1
    assert n_checked >= 200


def test_survival_indet_inputs():
    assert is_indet(cert_lognormal_survival(float("nan"), 0.0))
    assert is_indet(cert_lognormal_survival(float("inf"), 0.0))
    assert is_indet(cert_lognormal_survival(1.0, float("nan")))
    assert is_indet(cert_lognormal_survival(1.0, float("inf")))   # +inf mu
    assert is_indet(cert_lognormal_survival(0.0, 0.0))            # log(0)
    assert is_indet(cert_lognormal_survival(-1.0, 0.0))
    assert is_indet(cert_lognormal_survival(1.0, 0.0, sigma=0.0))
    assert is_indet(cert_lognormal_survival(1.0, 0.0, sigma=float("nan")))
    # ball x straddling zero: not certainly positive -> INDET
    assert is_indet(cert_lognormal_survival(arb(0, 1), 0.0))


# ---------------------------------------------------------------------------
# 2. cert_weighted_survival_sum
# ---------------------------------------------------------------------------

def test_weighted_sum_containment_random():
    rng = random.Random(31337)
    for _ in range(250):
        n = rng.randint(2, 20)
        g = [rng.random() for _ in range(n)]
        s = sum(g)
        probs = [gi / s for gi in g]             # random simplex, float64
        mus = [NEG_INF if rng.random() < 0.15 else rng.uniform(-8.0, 3.0)
               for _ in range(n)]
        x = 10.0 ** rng.uniform(-5, 2)
        ball = cert_weighted_survival_sum(probs, mus, x)
        assert not is_indet(ball)
        v = oracle_weighted_sum(probs, mus, x)
        assert contains_mpf(ball, v)
        assert ball_width(ball) < mp.mpf("1e-25")


def test_weighted_sum_planted_weight_error():
    rng = random.Random(777)
    n_effective = 0
    for _ in range(50):
        n = 8
        g = [rng.random() for _ in range(n)]
        s = sum(g)
        probs = [gi / s for gi in g]
        mus = [rng.uniform(0.0, 3.0) for _ in range(n)]   # survivals near 1
        x = 1e-3
        ball = cert_weighted_survival_sum(probs, mus, x)
        width = ball_width(ball)
        # planted error: perturb one weight by 1e-12
        probs_pert = list(probs)
        probs_pert[0] = probs_pert[0] + 1e-12
        v_pert = oracle_weighted_sum(probs_pert, mus, x)
        if width < mp.mpf("1e-13"):
            # the perturbed truth must land OUTSIDE the unperturbed enclosure
            assert not contains_mpf(ball, v_pert)
            n_effective += 1
        # sanity: the unperturbed truth stays inside
        assert contains_mpf(ball, oracle_weighted_sum(probs, mus, x))
    assert n_effective >= 45      # prec=128 widths are ~1e-36: all should count


def test_weighted_sum_upstream_ball_composition():
    # law: arb-ball weights pass through and their radius composes
    probs = [arb(0.6, 1e-6), arb(0.4, 1e-6)]
    mus = [0.5, -1.0]
    x = 0.5
    ball = cert_weighted_survival_sum(probs, mus, x)
    assert not is_indet(ball)
    # width must reflect the upstream weight uncertainty (~2e-6 * survival)
    assert ball_width(ball) > mp.mpf("1e-7")
    # truth for any weights within the balls must be contained: check center
    v = oracle_weighted_sum([0.6, 0.4], mus, x)
    assert contains_mpf(ball, v)


def test_weighted_sum_indet_inputs():
    assert is_indet(cert_weighted_survival_sum([float("nan"), 0.5], [0.0, 1.0], 1.0))
    assert is_indet(cert_weighted_survival_sum([float("inf"), 0.5], [0.0, 1.0], 1.0))
    assert is_indet(cert_weighted_survival_sum([0.5, 0.5], [0.0, float("nan")], 1.0))
    assert is_indet(cert_weighted_survival_sum([0.5, 0.5], [0.0, 1.0], float("nan")))
    with pytest.raises(ValueError):
        cert_weighted_survival_sum([0.5], [0.0, 1.0], 1.0)


# ---------------------------------------------------------------------------
# 3. cert_normalized_weighted_sum
# ---------------------------------------------------------------------------

def test_normalized_containment_random_float_terms():
    rng = random.Random(90210)
    for _ in range(200):
        n = rng.randint(2, 15)
        w = [rng.random() for _ in range(n)]              # unnormalized weights
        f = [rng.uniform(0.0, 100.0) for _ in range(n)]   # intensities
        a = [wi * fi for wi, fi in zip(w, f)]             # float64 products
        ball = cert_normalized_weighted_sum(a, w)
        assert not is_indet(ball)
        v = mp.fsum(mp.mpf(ai) for ai in a) / mp.fsum(mp.mpf(wi) for wi in w)
        assert contains_mpf(ball, v)
        assert ball_width(ball) < mp.mpf("1e-25")


def test_normalized_containment_with_tightener():
    rng = random.Random(11235)
    for _ in range(200):
        n = rng.randint(2, 15)
        w = [rng.random() for _ in range(n)]
        f = [rng.uniform(0.0, 100.0) for _ in range(n)]
        with working_prec(128):
            a = [pin(wi) * pin(fi) for wi, fi in zip(w, f)]  # exact at 128 bits
        ball = cert_normalized_weighted_sum(a, w, values=f)
        assert not is_indet(ball)
        v = mp.fsum(mp.mpf(wi) * mp.mpf(fi) for wi, fi in zip(w, f)) \
            / mp.fsum(mp.mpf(wi) for wi in w)
        assert contains_mpf(ball, v)
        # convex-combination tightener: result within [min f, max f]
        lo, hi = ball_bounds(ball)
        assert lo >= min(mp.mpf(fi) for fi in f) - mp.mpf("1e-30")
        assert hi <= max(mp.mpf(fi) for fi in f) + mp.mpf("1e-30")


def test_normalized_tightener_consistency_check():
    # numer inconsistent with weights*values -> loud caller-misuse error
    with pytest.raises(ValueError):
        cert_normalized_weighted_sum([1.0, 1.0], [0.5, 0.5], values=[10.0, 20.0])


def test_normalized_indet_inputs():
    assert is_indet(cert_normalized_weighted_sum([1.0], [0.0]))          # /0
    assert is_indet(cert_normalized_weighted_sum([1.0, 1.0], [1.0, -1.0]))
    assert is_indet(cert_normalized_weighted_sum([1.0, arb(0, 1)], [1.0, arb(0, 2)]))
    assert is_indet(cert_normalized_weighted_sum([float("nan")], [1.0]))
    assert is_indet(cert_normalized_weighted_sum([1.0], [float("inf")]))
    with pytest.raises(ValueError):
        cert_normalized_weighted_sum([1.0], [1.0, 2.0])


# ---------------------------------------------------------------------------
# 4. cert_percentile_bracket
# ---------------------------------------------------------------------------

def test_percentile_containment_random():
    rng = random.Random(65065)
    prs = [0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5]
    n_checked = 0
    n_skipped = 0
    n_gt = 0
    tries = 0
    while n_checked < 220 and tries < 1000:
        tries += 1
        if tries % 3:      # ordinary regime (zeros in the tail common)
            _, _, hc = make_hc(rng, rng.randint(2, 12), THRESHOLDS_65)
            pr = rng.choice(prs)
        else:              # heavy-tail regime: all-positive curve, small pr
            _, _, hc = make_hc(rng, rng.randint(2, 12), THRESHOLDS_65,
                               mu_lo=1.0, mu_hi=6.0, neginf_frac=0.0)
            pr = rng.choice([0.01, 0.05])
        v, flag = oracle_percentile(hc, THRESHOLDS_65, pr)
        if flag == "nonmonotone":
            n_skipped += 1
            continue
        ball, status = cert_percentile_bracket(hc, THRESHOLDS_65, pr)
        # exact float64-pinned inputs: every path decision is certain
        assert status in (STATUS_OK, STATUS_GT100M_BRANCH), status
        assert (status == STATUS_GT100M_BRANCH) == (flag == "gt100m"), (status, flag)
        assert contains_mpf(ball, v), (pr, status)
        assert ball_width(ball) < mp.mpf("1e-9")
        n_gt += status == STATUS_GT100M_BRANCH
        n_checked += 1
    assert n_checked >= 220, (n_checked, n_skipped)
    assert n_gt >= 10     # the >100 m branch must be exercised randomly


def test_percentile_gt100m_branch_certain():
    # all-positive curve whose tail prob exceeds pr: shipped >100m branch
    hc = [0.9, 0.5, 0.1, 0.005]
    th = [0.1, 1.0, 10.0, 100.0]
    ball, status = cert_percentile_bracket(hc, th, 0.001)
    assert status == STATUS_GT100M_BRANCH
    # the shipped value is the PROBABILITY 0.005 assigned as intensity
    assert contains_mpf(ball, mp.mpf(0.005))
    assert ball_width(ball) < mp.mpf("1e-20")


def test_percentile_dedup_exact_duplicates():
    # consecutive exact duplicates: dedup keeps the FIRST of the run
    hc = [0.8, 0.5, 0.5, 0.3, 0.1, 0.0]
    th = [0.01, 0.1, 0.5, 1.0, 10.0, 100.0]
    ball, status = cert_percentile_bracket(hc, th, 0.5)
    assert status == STATUS_OK
    v, flag = oracle_percentile(hc, th, 0.5)
    assert flag == ""
    assert v == mp.mpf(0.1)            # knot value at the first 0.5
    assert contains_mpf(ball, v)
    # certain dedup must keep the bracket tight: t2=0.5 excluded
    assert not contains_mpf(ball, mp.mpf(0.5))


def test_percentile_ambiguous_path_union():
    # tail entry ball straddles 0: nonNull membership ambiguous -> the bracket
    # must cover BOTH convention-consistent resolutions and say so
    th = [0.1, 1.0, 10.0]
    amb_entry = arb(0.02, 0.03)                     # possibly 0, possibly ~0.05
    hc_balls = [arb(0.9), arb(0.5), amb_entry]
    pr = 0.01
    ball, status = cert_percentile_bracket(hc_balls, th, pr)
    assert status == STATUS_AMBIGUOUS
    # resolution A: entry treated as zero -> hc = [0.9, 0.5, 0]
    vA, fA = oracle_percentile([0.9, 0.5, 0.0], th, pr)
    assert fA == ""
    assert contains_mpf(ball, vA)
    # resolution B: entry positive, e.g. 0.02 (inside its ball)
    vB, fB = oracle_percentile([0.9, 0.5, 0.02], th, pr)
    assert contains_mpf(ball, vB)


def test_percentile_all_zero_curve():
    ball, status = cert_percentile_bracket([0.0] * 5, [0.1] * 5, 0.1)
    assert status == STATUS_OK
    assert ball.is_zero()               # shipped valTmp = 0 (caller maps level)


def test_percentile_indet_inputs():
    th = [0.1, 1.0, 10.0]
    ball, status = cert_percentile_bracket([0.5, float("nan"), 0.0], th, 0.1)
    assert status == STATUS_INDET and is_indet(ball)
    ball, status = cert_percentile_bracket([0.5, float("inf"), 0.0], th, 0.1)
    assert status == STATUS_INDET and is_indet(ball)
    ball, status = cert_percentile_bracket([0.5, 0.2, 0.1], th, float("nan"))
    assert status == STATUS_INDET and is_indet(ball)
    # certainly-negative hc entry: invalid geometry
    ball, status = cert_percentile_bracket([0.5, -0.1, 0.0], th, 0.1)
    assert status == STATUS_INDET and is_indet(ball)
    # index-past-end edge case of the reference code: interior zero with hc[end] > 0
    ball, status = cert_percentile_bracket([0.5, 0.0, 0.1], th, 0.1)
    assert status == STATUS_INDET and is_indet(ball)


def test_percentile_fallback_hull_sound():
    # many ambiguous entries -> path cap -> global hull, still fail-closed
    th = list(THRESHOLDS_65[:20])
    hc_balls = [arb(1e-6, 2e-6) for _ in range(20)]     # 20 straddling entries
    ball, status = cert_percentile_bracket(hc_balls, th, 0.01, max_paths=8)
    assert status == STATUS_AMBIGUOUS
    lo, hi = ball_bounds(ball)
    assert lo <= 0 and hi >= 1          # hull covers value range + >100 m branch


# ---------------------------------------------------------------------------
# 5. cert_alert_level
# ---------------------------------------------------------------------------

def test_alert_random_degenerate_balls():
    rng = random.Random(555)
    bounds = (0.0, 0.10, 0.5)
    n = 0
    for _ in range(300):
        pick = rng.random()
        if pick < 0.5:
            v = 10.0 ** rng.uniform(-6, 1)
        elif pick < 0.8:
            b = rng.choice([0.0, 0.10, 0.5])
            v = b + rng.uniform(-1e-3, 1e-3)
        else:
            v = rng.choice([0.0, 0.10, 0.5, -1.0, 0.6, 100.0])
        lvl = cert_alert_level(arb(float(v)), bounds)   # degenerate exact ball
        assert lvl == oracle_alert_level(v, bounds), v
        n += 1
    assert n >= 300


def test_alert_nondegenerate_balls_random():
    rng = random.Random(556)
    bounds = (0.0, 0.10, 0.5)
    for _ in range(200):
        v = 10.0 ** rng.uniform(-6, 1)
        rad = 10.0 ** rng.uniform(-12, -2)
        ball = arb(float(v), float(rad))
        out = cert_alert_level(ball, bounds)
        straddles = any(v - rad <= b <= v + rad for b in bounds)
        if straddles:
            assert out == ALERT_STRADDLES
        else:
            assert out == oracle_alert_level(v, bounds)


def test_alert_boundary_conventions_float64():
    bounds = (0.0, 0.10, 0.5)
    # degenerate exact ball ON the float64 boundary -> the <= side (level 1)
    assert cert_alert_level(arb(0.10), bounds, boundary_semantics="float64") == 1
    assert cert_alert_level(arb(0.5), bounds, boundary_semantics="float64") == 2
    assert cert_alert_level(arb(0.0), bounds, boundary_semantics="float64") == 0
    # nondegenerate ball containing the boundary -> STRADDLES
    assert cert_alert_level(arb(0.10, 1e-20), bounds,
                            boundary_semantics="float64") == ALERT_STRADDLES
    assert cert_alert_level(arb(0.5, 1e-30), bounds,
                            boundary_semantics="float64") == ALERT_STRADDLES
    # clear cases
    assert cert_alert_level(arb(0.05), bounds) == 1
    assert cert_alert_level(arb(0.3), bounds) == 2
    assert cert_alert_level(arb(0.6), bounds) == 3
    assert cert_alert_level(arb(-1.0), bounds) == 0


def test_alert_boundary_conventions_decimal():
    bounds = (0.0, 0.10, 0.5)
    # float64 0.1 is ~5.55e-18 ABOVE decimal 0.1: under decimal semantics the
    # exact double is certainly on the > side -> level 2
    assert cert_alert_level(arb(0.10), bounds, boundary_semantics="decimal") == 2
    # a decimal-0.1 value ball overlaps the decimal-0.1 boundary ball -> STRADDLES
    with working_prec(128):
        v = arb("0.1")
    assert cert_alert_level(v, bounds,
                            boundary_semantics="decimal") == ALERT_STRADDLES
    # away from 0.1 the semantics agree
    assert cert_alert_level(arb(0.05), bounds, boundary_semantics="decimal") == 1
    assert cert_alert_level(arb(0.3), bounds, boundary_semantics="decimal") == 2


def test_alert_indet_inputs():
    bounds = (0.0, 0.10, 0.5)
    assert cert_alert_level(float("nan"), bounds) == ALERT_INDET
    assert cert_alert_level(float("inf"), bounds) == ALERT_INDET
    assert cert_alert_level(float("-inf"), bounds) == ALERT_INDET
    assert cert_alert_level(INDET(), bounds) == ALERT_INDET


def test_alert_bad_semantics_flag():
    with pytest.raises(ValueError):
        cert_alert_level(arb(0.05), (0.0, 0.10, 0.5), boundary_semantics="matlab")


# ---------------------------------------------------------------------------
# timings (informational; loose sanity bound only)
# ---------------------------------------------------------------------------

def test_timings_prec128(capsys):
    rng = random.Random(1)
    timings = {}

    t0 = time.perf_counter()
    reps = 2000
    for _ in range(reps):
        cert_lognormal_survival(0.5, -1.0)
    timings["cert_lognormal_survival (1 call)"] = (time.perf_counter() - t0) / reps

    n = 100
    probs = [1.0 / n] * n
    mus = [rng.uniform(-8, 3) for _ in range(n)]
    t0 = time.perf_counter()
    reps = 50
    for _ in range(reps):
        cert_weighted_survival_sum(probs, mus, 0.5)
    timings["cert_weighted_survival_sum (n=100)"] = (time.perf_counter() - t0) / reps

    w = [rng.random() for _ in range(100)]
    f = [rng.uniform(0, 100) for _ in range(100)]
    with working_prec(128):
        a = [pin(wi) * pin(fi) for wi, fi in zip(w, f)]
    t0 = time.perf_counter()
    reps = 500
    for _ in range(reps):
        cert_normalized_weighted_sum(a, w, values=f)
    timings["cert_normalized_weighted_sum (n=100)"] = (time.perf_counter() - t0) / reps

    _, _, hc = make_hc(rng, 8, THRESHOLDS_65)
    t0 = time.perf_counter()
    reps = 200
    for _ in range(reps):
        cert_percentile_bracket(hc, THRESHOLDS_65, 0.1)
    timings["cert_percentile_bracket (65-grid)"] = (time.perf_counter() - t0) / reps

    t0 = time.perf_counter()
    reps = 5000
    for _ in range(reps):
        cert_alert_level(arb(0.3), (0.0, 0.10, 0.5))
    timings["cert_alert_level (1 call)"] = (time.perf_counter() - t0) / reps

    with capsys.disabled():
        print("\n--- certified-lane per-call timings at prec=128 ---")
        for k, v in timings.items():
            print(f"  {k}: {v * 1e6:.1f} us")
    for k, v in timings.items():
        assert v < 0.25, (k, v)         # loose sanity bound


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-q"]))

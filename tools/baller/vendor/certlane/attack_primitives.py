# baller engine — adversarial attack harness for certlane/primitives.py.
"""ADVERSARIAL ATTACK on certlane/primitives.py — fresh-context skeptic harness.

Run from the vendor directory (the parent of certlane/):
    python3 certlane/attack_primitives.py

Four attack surfaces:
  1. CONTAINMENT   — 2000+ randomized cases per primitive vs mpmath dps=80,
                     incl. adversarial regimes (mu ~ ln x, sigma == 1, grid
                     endpoints 1e-5/100, 1e-300 probs, 1e5-term positive sums,
                     ball-input inclusion fuzz).
  2. PLANTED ERRORS — corrupt internals in a subprocess (shifted erfc, dropped
                     sum term, keep-last dedup, deleted >100 m branch) and demand
                     the shipped test suite FAILS each time.
  3. SEMANTICS     — recorded survival-curve archives (ATTACK_EVENT_NPZ):
                     the float64 reference percentile_value must lie inside
                     the certified bracket for every POI x Pr level; widths
                     measured at prec=128.  Named SKIP when the reference
                     implementation or the archives are absent.
  4. PRECISION LAWS — runtime probes of the prec-context / pinning laws.

Exit code 0 iff no FATAL failure.  Machine-readable JSON summary on stdout
(last line, prefixed ATTACK_JSON:).
"""

import json
import math
import os
import random
import subprocess
import sys
import time
from fractions import Fraction

import mpmath as mp
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from flint import arb, ctx

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

mp.mp.dps = 80
NEG_INF = float("-inf")
THRESHOLDS_65 = np.logspace(-5, 2, 65)

CHECKS = 0
FAILURES = []       # list of dicts {surface, desc, repro}
NOTES = []


def bump(n=1):
    global CHECKS
    CHECKS += n


def fail(surface, desc, repro):
    FAILURES.append({"surface": surface, "desc": desc, "repro": repro})
    print(f"  FAIL [{surface}] {desc}")


# ---------------------------------------------------------------------------
# ball -> exact rational bounds
# ---------------------------------------------------------------------------

def _bound_to_mpf(b):
    assert not b.is_nan()
    if not b.is_finite():
        return mp.inf if b > arb(0) else mp.ninf
    if b.is_zero():
        return mp.mpf(0)
    m, e = b.man_exp()
    return mp.mpf(int(m)) * mp.power(2, int(e))


def ball_bounds(ball):
    with working_prec(256):
        return _bound_to_mpf(ball.lower()), _bound_to_mpf(ball.upper())


def contains_mpf(ball, v):
    lo, hi = ball_bounds(ball)
    return lo <= v <= hi


def ball_width(ball):
    lo, hi = ball_bounds(ball)
    return hi - lo


# ---------------------------------------------------------------------------
# dps=80 oracles
# ---------------------------------------------------------------------------

def oracle_survival(x, mu, sigma=1.0):
    if mu == NEG_INF:
        return mp.mpf(0)
    t = (mp.log(mp.mpf(x)) - mp.mpf(mu)) / (mp.mpf(sigma) * mp.sqrt(2))
    return mp.erfc(t) / 2


def oracle_percentile(hc_floats, thresholds, pr):
    """Shipped interp1 convention at dps=80 on exact float64 inputs.
    Returns (value, flag) flag in {'', 'gt100m'}; (None,'nonmonotone') when the
    deduped xx is not strictly decreasing; (None,'crash') for the latent
    index-past-end edge geometry; (None,'range') for pr outside (0, xx[0]]."""
    n = len(hc_floats)
    nn = [i for i in range(n) if hc_floats[i] > 0.0]
    if n > len(nn) and nn:
        nxt = nn[-1] + 1
        if nxt >= n:
            return None, "crash"
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
    if prm > xxu[0]:
        return None, "range"
    for i in range(len(xxu) - 1):
        if xxu[i] >= prm >= xxu[i + 1]:
            hi, lo = xxu[i], xxu[i + 1]
            return yyu[i] + (yyu[i + 1] - yyu[i]) * (hi - prm) / (hi - lo), ""
    raise AssertionError("unreachable")


def oracle_alert_level(v, bounds=(0.0, 0.10, 0.5)):
    for i, b in enumerate(bounds):
        if v <= b:
            return i
    return len(bounds)


def oracle_alert_level_decimal(v):
    """Exact-rational alert binning against DECIMAL boundaries 0, 1/10, 1/2."""
    fv = Fraction(v)
    for i, b in enumerate((Fraction(0), Fraction(1, 10), Fraction(1, 2))):
        if fv <= b:
            return i
    return 3


def float_survival(x, mu, sigma=1.0):
    from scipy.special import erfc as _erfc
    if mu == NEG_INF:
        return 0.0
    return 0.5 * float(_erfc((math.log(x) - mu) / (sigma * math.sqrt(2.0))))


def make_hc(rng, n_scen, thresholds, mu_lo=-7.0, mu_hi=2.0, neginf_frac=0.15):
    g = np.array([rng.random() for _ in range(n_scen)])
    probs = (g / g.sum()).astype(np.float64)
    mus = np.array([NEG_INF if rng.random() < neginf_frac
                    else rng.uniform(mu_lo, mu_hi) for _ in range(n_scen)])
    hc = np.zeros(len(thresholds))
    for k, x in enumerate(thresholds):
        hc[k] = sum(p * float_survival(float(x), m) for p, m in zip(probs, mus))
    return probs, mus, hc


# ===========================================================================
# 1. CONTAINMENT
# ===========================================================================

def attack_survival():
    print("== containment: cert_lognormal_survival ==")
    rng = random.Random(80017)
    t0 = time.perf_counter()
    n = 0
    cases = []
    # regime A: broad random, sigma exactly 1 (70%) or random
    for _ in range(1100):
        x = 10.0 ** rng.uniform(-5, 2)
        mu = NEG_INF if rng.random() < 0.08 else rng.uniform(-8.0, 3.0)
        sigma = 1.0 if rng.random() < 0.7 else rng.uniform(0.3, 2.5)
        cases.append((x, mu, sigma))
    # regime B: mu near ln(x) — erfc argument sign boundary
    for _ in range(500):
        x = 10.0 ** rng.uniform(-5, 2)
        lx = math.log(x)
        d = 0.0 if rng.random() < 0.2 else \
            (rng.choice([-1, 1]) * 10.0 ** rng.uniform(-18, 0))
        cases.append((x, lx + d, 1.0))
    # regime C: exact grid endpoints
    for _ in range(300):
        x = rng.choice([1e-5, 100.0, float(THRESHOLDS_65[0]),
                        float(THRESHOLDS_65[-1])])
        mu = NEG_INF if rng.random() < 0.1 else rng.uniform(-8.0, 3.0)
        cases.append((x, mu, 1.0))
    # regime D: deep tails and extreme sigma
    for _ in range(300):
        x = rng.choice([1e-5, 100.0, 10.0 ** rng.uniform(-5, 2)])
        mu = rng.uniform(-60.0, 40.0)
        sigma = rng.choice([1.0, 1e-3, 1e3, rng.uniform(0.1, 10.0)])
        cases.append((x, mu, sigma))
    # regime E: extreme underflow (t enormous)
    for mu in (-1e6, -1e10, 1e6, 1e10):
        cases.append((1e-5, mu, 1.0))
        cases.append((100.0, mu, 1.0))

    for (x, mu, sigma) in cases:
        ball = cert_lognormal_survival(x, mu, sigma)
        bump()
        if is_indet(ball):
            fail("containment", f"survival INDET on valid input {(x, mu, sigma)}",
                 f"cert_lognormal_survival({x!r},{mu!r},{sigma!r})")
            continue
        if mu == NEG_INF:
            if not (ball.is_zero() and ball.is_exact()):
                fail("containment", f"mu=-inf not exact zero, x={x}",
                     f"cert_lognormal_survival({x!r},float('-inf'))")
            continue
        v = oracle_survival(x, mu, sigma)
        if not contains_mpf(ball, v):
            fail("containment",
                 f"TRUE VALUE OUTSIDE BALL survival({x},{mu},{sigma})",
                 f"cert_lognormal_survival({x!r},{mu!r},{sigma!r}) vs mpmath dps=80")
        n += 1
    # ball-input inclusion fuzz: sampled truths inside input balls
    for _ in range(200):
        xm = 10.0 ** rng.uniform(-4, 1.5)
        xr = xm * 10.0 ** rng.uniform(-12, -2)
        mm = rng.uniform(-6, 2)
        mr = 10.0 ** rng.uniform(-12, -2)
        xb = arb(xm, xr)
        mb = arb(mm, mr)
        ball = cert_lognormal_survival(xb, mb, 1.0)
        bump()
        if is_indet(ball):
            continue
        for _ in range(3):
            xs = xm + rng.uniform(-1, 1) * xr
            ms = mm + rng.uniform(-1, 1) * mr
            v = oracle_survival(xs, ms)
            if not contains_mpf(ball, v):
                fail("containment",
                     f"ball-input survival: sampled truth escapes ({xm},{xr},{mm},{mr})",
                     f"cert_lognormal_survival(arb({xm},{xr}),arb({mm},{mr}))")
                break
    print(f"  {n} exact + 200 ball cases, {time.perf_counter()-t0:.1f}s")


def attack_weighted_sum():
    print("== containment: cert_weighted_survival_sum ==")
    rng = random.Random(80018)
    t0 = time.perf_counter()
    n_cases = 0
    for i in range(2000):
        n = rng.randint(1, 30)
        g = [rng.random() for _ in range(n)]
        s = sum(g)
        probs = [gi / s for gi in g]
        # underflow entries regime: sprinkle 1e-300-scale weights
        if i % 4 == 0:
            k = rng.randint(1, max(1, n // 2))
            for j in rng.sample(range(n), k):
                probs[j] = probs[j] * 1e-300
        mus = [NEG_INF if rng.random() < 0.15 else rng.uniform(-8.0, 3.0)
               for _ in range(n)]
        if i % 5 == 0:
            x = rng.choice([1e-5, 100.0])
        else:
            x = 10.0 ** rng.uniform(-5, 2)
        ball = cert_weighted_survival_sum(probs, mus, x)
        bump()
        if is_indet(ball):
            fail("containment", f"weighted sum INDET on valid input (case {i})",
                 f"seed 80018 case {i}")
            continue
        v = mp.fsum(mp.mpf(p) * oracle_survival(x, m, 1.0)
                    for p, m in zip(probs, mus))
        if not contains_mpf(ball, v):
            fail("containment",
                 f"TRUE VALUE OUTSIDE BALL weighted sum case {i} n={n} x={x}",
                 f"attack_primitives.attack_weighted_sum seed 80018 case {i}")
        n_cases += 1
    # 1e5-term cancellation-free positive sum
    n_terms = 100_000
    n_distinct = 2000
    mus_d = [rng.uniform(-8.0, 3.0) for _ in range(n_distinct)]
    probs = []
    mus = []
    for j in range(n_terms):
        probs.append(rng.random() / n_terms)
        mus.append(mus_d[j % n_distinct])
    x = 0.37
    tb = time.perf_counter()
    ball = cert_weighted_survival_sum(probs, mus, x)
    t_big = time.perf_counter() - tb
    bump()
    surv_cache = {m: oracle_survival(x, m) for m in mus_d}
    v = mp.fsum(mp.mpf(p) * surv_cache[m] for p, m in zip(probs, mus))
    if is_indet(ball) or not contains_mpf(ball, v):
        fail("containment", "1e5-term positive sum escapes ball",
             "attack_weighted_sum big-sum block, seed 80018")
    w_big = ball_width(ball)
    NOTES.append(f"1e5-term sum: cert {t_big:.2f}s, width {mp.nstr(w_big, 3)}")
    print(f"  {n_cases} cases + 1e5-term sum (width {mp.nstr(w_big, 3)}), "
          f"{time.perf_counter()-t0:.1f}s")


def attack_normalized():
    print("== containment: cert_normalized_weighted_sum ==")
    rng = random.Random(80019)
    t0 = time.perf_counter()
    for i in range(2000):
        n = rng.randint(1, 20)
        w = [rng.random() * (1e-300 if (i % 7 == 0 and rng.random() < 0.3) else 1.0)
             for _ in range(n)]
        f = [rng.uniform(0.0, 100.0) for _ in range(n)]
        use_tight = i % 2 == 0
        if use_tight:
            with working_prec(128):
                a = [pin(wi) * pin(fi) for wi, fi in zip(w, f)]
            ball = cert_normalized_weighted_sum(a, w, values=f)
            v = mp.fsum(mp.mpf(wi) * mp.mpf(fi) for wi, fi in zip(w, f)) \
                / mp.fsum(mp.mpf(wi) for wi in w)
        else:
            a = [wi * fi for wi, fi in zip(w, f)]     # float64 products
            ball = cert_normalized_weighted_sum(a, w)
            v = mp.fsum(mp.mpf(ai) for ai in a) / mp.fsum(mp.mpf(wi) for wi in w)
        bump()
        if is_indet(ball):
            if sum(w) > 0:
                fail("containment", f"normalized sum INDET, case {i}",
                     f"seed 80019 case {i}")
            continue
        if not contains_mpf(ball, v):
            fail("containment", f"TRUE VALUE OUTSIDE BALL normalized case {i}",
                 f"seed 80019 case {i} tight={use_tight}")
        if use_tight:
            lo, hi = ball_bounds(ball)
            fmin = min(mp.mpf(fi) for fi in f)
            fmax = max(mp.mpf(fi) for fi in f)
            if lo < fmin - mp.mpf("1e-28") or hi > fmax + mp.mpf("1e-28"):
                fail("containment", f"tightener hull violated case {i}",
                     f"seed 80019 case {i}")
    print(f"  2000 cases, {time.perf_counter()-t0:.1f}s")


def attack_percentile():
    print("== containment: cert_percentile_bracket ==")
    rng = random.Random(80020)
    t0 = time.perf_counter()
    prs = [0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5]
    n_checked = n_skip = n_gt = 0
    tries = 0
    while n_checked < 2000 and tries < 6000:
        tries += 1
        r = tries % 4
        if r in (0, 1):
            _, _, hc = make_hc(rng, rng.randint(1, 12), THRESHOLDS_65)
            pr = rng.choice(prs)
        elif r == 2:   # heavy tail: all-positive, tiny pr -> >100 m branch
            _, _, hc = make_hc(rng, rng.randint(2, 12), THRESHOLDS_65,
                               mu_lo=1.0, mu_hi=7.0, neginf_frac=0.0)
            pr = rng.choice([0.001, 0.01, 0.05])
        else:          # pr planted exactly on a knot value
            _, _, hc = make_hc(rng, rng.randint(2, 10), THRESHOLDS_65)
            pos = [h for h in hc if h > 0]
            if not pos:
                continue
            pr = float(rng.choice(pos))
            if not (0 < pr <= 1):
                continue
        hcl = [float(h) for h in hc]
        v, flag = oracle_percentile(hcl, THRESHOLDS_65, pr)
        if v is None:
            n_skip += 1
            continue
        ball, status = cert_percentile_bracket(hcl, THRESHOLDS_65, pr)
        bump()
        if status not in (STATUS_OK, STATUS_GT100M_BRANCH):
            fail("containment",
                 f"exact float64 inputs gave status {status} (try {tries})",
                 f"seed 80020 try {tries}")
            continue
        if (status == STATUS_GT100M_BRANCH) != (flag == "gt100m"):
            fail("containment",
                 f">100 m-branch status mismatch cert={status} oracle={flag}",
                 f"seed 80020 try {tries}")
        if not contains_mpf(ball, v):
            fail("containment",
                 f"TRUE VALUE OUTSIDE BRACKET percentile try {tries} pr={pr}",
                 f"seed 80020 try {tries}")
        n_gt += status == STATUS_GT100M_BRANCH
        n_checked += 1

    # hand-crafted geometries
    crafted = [
        # (hc, th, pr) exact duplicates mid-curve
        ([0.8, 0.5, 0.5, 0.5, 0.3, 0.0], [0.01, 0.1, 0.5, 1.0, 10.0, 100.0], 0.5),
        ([0.8, 0.5, 0.5, 0.3, 0.1, 0.0], [0.01, 0.1, 0.5, 1.0, 10.0, 100.0], 0.4),
        # survival exactly 1 at head (anchor duplicate)
        ([1.0, 1.0, 0.6, 0.2, 0.0], [0.01, 0.1, 1.0, 10.0, 100.0], 0.5),
        ([1.0, 0.6, 0.2, 0.0], [0.01, 1.0, 10.0, 100.0], 1.0),   # pr == anchor
        # >100 m branch exact
        ([0.9, 0.5, 0.1, 0.005], [0.1, 1.0, 10.0, 100.0], 0.001),
        ([0.9, 0.5, 0.1, 0.005], [0.1, 1.0, 10.0, 100.0], 0.005),  # pr == xx(end)
        # single positive entry + appended zero
        ([0.5, 0.0, 0.0], [0.1, 1.0, 10.0], 0.25),
        ([0.5, 0.0, 0.0], [0.1, 1.0, 10.0], 0.5),
        # all-zero
        ([0.0] * 5, list(THRESHOLDS_65[:5]), 0.1),
        # tiny probabilities (underflow scale)
        ([1e-300, 1e-305, 0.0], [0.1, 1.0, 10.0], 0.2),
        ([0.7, 1e-300, 0.0], [0.1, 1.0, 10.0], 0.3),
        ([0.7, 1e-300, 0.0], [0.1, 1.0, 10.0], 1e-301),  # below xx(end): >100 m branch
    ]
    for (hcl, th, pr) in crafted:
        v, flag = oracle_percentile(hcl, th, pr)
        ball, status = cert_percentile_bracket(hcl, th, pr)
        bump()
        if v is None:
            if flag == "crash" and status != STATUS_INDET:
                fail("containment", f"crash geometry not INDET: {hcl}",
                     f"cert_percentile_bracket({hcl},{th},{pr})")
            continue
        if status not in (STATUS_OK, STATUS_GT100M_BRANCH) or not contains_mpf(ball, v):
            fail("containment", f"crafted percentile case fails: {hcl} pr={pr} "
                 f"status={status}",
                 f"cert_percentile_bracket({hcl},{th},{pr})")
        if (status == STATUS_GT100M_BRANCH) != (flag == "gt100m"):
            fail("containment", f"crafted >100 m-branch status mismatch {hcl} pr={pr}",
                 f"cert_percentile_bracket({hcl},{th},{pr})")

    # ball-input path-union fuzz: sampled truths must stay inside the bracket
    n_ball = 0
    for i in range(300):
        nn = rng.randint(2, 6)
        th = sorted(10.0 ** rng.uniform(-2, 2) for _ in range(nn))
        mids, rads = [], []
        prev = 1.0
        for k in range(nn):
            prev = prev * rng.uniform(0.2, 0.95)
            mids.append(prev)
            rads.append(prev * 10.0 ** rng.uniform(-10, -1)
                        if rng.random() < 0.6 else 0.0)
        # sometimes plant a near-duplicate pair or a near-zero tail
        if rng.random() < 0.4 and nn >= 3:
            mids[2] = mids[1]
            rads[2] = mids[1] * 1e-9
        if rng.random() < 0.4:
            mids[-1] = 1e-9
            rads[-1] = 2e-9        # straddles 0
        hc_balls = [arb(m, r) for m, r in zip(mids, rads)]
        pr = rng.choice([0.05, 0.1, 0.3])
        ball, status = cert_percentile_bracket(hc_balls, th, pr)
        bump()
        if status == STATUS_INDET:
            continue
        lo, hi = ball_bounds(ball)
        for _ in range(5):
            samp = [max(0.0, m + rng.uniform(-1, 1) * r)
                    for m, r in zip(mids, rads)]
            v, flag = oracle_percentile(samp, th, pr)
            if v is None:
                continue
            if not (lo <= v <= hi):
                fail("containment",
                     f"path-union bracket misses sampled truth (fuzz {i})",
                     f"seed 80020 ball-fuzz case {i}")
                break
        n_ball += 1
    print(f"  {n_checked} random ({n_gt} >100 m-branch, {n_skip} skipped) "
          f"+ {len(crafted)} crafted + {n_ball} ball-fuzz, "
          f"{time.perf_counter()-t0:.1f}s")


def attack_alert():
    print("== containment: cert_alert_level ==")
    rng = random.Random(80021)
    t0 = time.perf_counter()
    bounds = (0.0, 0.10, 0.5)
    for i in range(2200):
        pick = rng.random()
        if pick < 0.4:
            v = 10.0 ** rng.uniform(-8, 1)
        elif pick < 0.7:
            b = rng.choice([0.0, 0.10, 0.5])
            v = b + rng.choice([-1, 1]) * 10.0 ** rng.uniform(-18, -2)
        elif pick < 0.85:
            b = rng.choice([0.0, 0.10, 0.5])
            v = float(np.nextafter(b, rng.choice([-np.inf, np.inf])))
        else:
            v = rng.choice([0.0, 0.10, 0.5, -1.0, 0.6, 100.0,
                            float(np.nextafter(0.1, 0)),
                            float(np.nextafter(0.1, 1))])
        lvl = cert_alert_level(arb(float(v)), bounds)
        bump()
        want = oracle_alert_level(v, bounds)
        if lvl != want:
            fail("containment", f"alert float64 semantics: v={v!r} cert={lvl} "
                 f"oracle={want}", f"cert_alert_level(arb({v!r}))")
    # nondegenerate balls
    for i in range(600):
        v = 10.0 ** rng.uniform(-6, 1)
        rad = 10.0 ** rng.uniform(-14, -2)
        out = cert_alert_level(arb(v, rad), bounds)
        bump()
        straddles = any(v - rad <= b <= v + rad for b in bounds)
        if straddles:
            if out != ALERT_STRADDLES:
                fail("containment", f"ball over boundary not STRADDLES v={v} r={rad}",
                     f"cert_alert_level(arb({v},{rad}))")
        elif out != oracle_alert_level(v, bounds):
            fail("containment", f"clear ball wrong level v={v} r={rad}",
                 f"cert_alert_level(arb({v},{rad}))")
    # decimal semantics vs exact rational oracle
    for i in range(300):
        if i % 3 == 0:
            v = rng.choice([0.10, 0.5, 0.0, float(np.nextafter(0.1, 0)),
                            float(np.nextafter(0.1, 1)),
                            float(np.nextafter(0.5, 0))])
        else:
            v = 10.0 ** rng.uniform(-6, 1)
        out = cert_alert_level(arb(float(v)), bounds, boundary_semantics="decimal")
        bump()
        want = oracle_alert_level_decimal(v)
        # STRADDLES acceptable only if v is within the boundary ball radius
        # (~2^-128); no double except an exact-binary boundary is that close.
        if out == ALERT_STRADDLES:
            near = any(abs(Fraction(v) - fb) < Fraction(1, 2 ** 100)
                       for fb in (Fraction(0), Fraction(1, 10), Fraction(1, 2)))
            if not near:
                fail("containment", f"decimal semantics spurious STRADDLES v={v!r}",
                     f"cert_alert_level(arb({v!r}),boundary_semantics='decimal')")
        elif out != want:
            fail("containment", f"decimal semantics wrong level v={v!r} "
                 f"cert={out} rational-oracle={want}",
                 f"cert_alert_level(arb({v!r}),boundary_semantics='decimal')")
    print(f"  3100 cases, {time.perf_counter()-t0:.1f}s")


# ===========================================================================
# 2. PLANTED ERRORS (subprocess: corrupt internals, suite must fail)
# ===========================================================================

CORRUPTIONS = ["none", "erfc_shift", "drop_term", "dedup_swap", "kill_gt100m_branch"]


def _apply_corruption(name):
    import certlane.primitives as P
    if name == "none":
        return
    if name == "erfc_shift":
        orig = arb.erfc
        arb.erfc = lambda self: orig(self) + arb("1e-30")
        return
    if name == "drop_term":
        orig = P.cert_weighted_survival_sum
        def bad(probs, mus, x, sigma=1.0, prec=P.DEFAULT_PREC):
            return orig(list(probs)[:-1], list(mus)[:-1], x, sigma, prec=prec)
        P.cert_weighted_survival_sum = bad
        return
    src_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            "primitives.py")
    src = open(src_path).read()
    if name == "dedup_swap":
        old = 'keep = [0] + [j + 1 for j, d in enumerate(dd) if d == "ne"]'
        new = 'keep = [j for j, d in enumerate(dd) if d == "ne"] + [len(seq) - 1]'
        assert old in src
        src = src.replace(old, new)
    elif name == "kill_gt100m_branch":
        old = "        cands.append(last)"
        assert old in src
        src = src.replace(old, "        pass", 1)
    else:
        raise ValueError(name)
    exec(compile(src, src_path, "exec"), P.__dict__)


def _corruption_child(name):
    _apply_corruption(name)
    import pytest
    code = pytest.main(["certlane/test_primitives.py", "-q", "--no-header", "-p",
                        "no:cacheprovider"])
    sys.exit(int(code))


def attack_planted():
    print("== planted errors: suite must catch each corruption ==")
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    py = sys.executable
    for name in CORRUPTIONS:
        t0 = time.perf_counter()
        r = subprocess.run(
            [py, os.path.abspath(__file__), "--corruption", name],
            cwd=here, capture_output=True, text=True, timeout=600)
        dt = time.perf_counter() - t0
        bump()
        caught = r.returncode != 0
        if name == "none":
            if caught:
                fail("planted", "BASELINE suite fails uncorrupted (harness broken)",
                     f"{py} certlane/attack_primitives.py --corruption none")
            else:
                print(f"  baseline: suite passes ({dt:.1f}s)  [harness valid]")
            continue
        tail = "\n".join(r.stdout.strip().splitlines()[-2:])
        if caught:
            print(f"  {name}: CAUGHT (exit {r.returncode}, {dt:.1f}s) {tail}")
        else:
            fail("planted", f"corruption '{name}' NOT caught by test suite",
                 f"{py} certlane/attack_primitives.py --corruption {name}")


# ===========================================================================
# 3. SEMANTICS on real curves
# ===========================================================================

def _ulp(v):
    av = abs(v)
    if av == 0.0:
        return float(np.nextafter(0.0, 1.0))
    return float(np.nextafter(av, np.inf) - av)


def attack_semantics():
    print("== semantics: recorded survival curves vs the float64 reference ==")
    try:
        from replay.alert import percentile_value, AL_TYPES
        from replay import conventions as C
    except ImportError:
        print("  SKIP (named): the float64 reference implementation (replay) "
              "is not vendored — put a checkout on PYTHONPATH to run this "
              "surface")
        return
    events = [p for p in
              os.environ.get("ATTACK_EVENT_NPZ", "").split(os.pathsep) if p]
    if not events:
        print("  SKIP (named): ATTACK_EVENT_NPZ unset (os.pathsep-separated "
              ".npz survival-curve archives)")
        return
    pr_methods = list(C.PR_METHODS)
    pr_cols = {m: AL_TYPES.index(m) for m in pr_methods}
    all_widths = []
    t0 = time.perf_counter()
    total = 0
    n_gt = n_zero = n_straddle = 0
    n_float_inside = 0
    n_float_outside = 0
    worst_ulp = (0.0, None)          # (ulp distance, (ev, ip, meth))
    ulp_dists = []
    n_not_nearest = 0
    n_oracle_outside = 0
    stat_counts = {}
    stored_mismatch = level_mismatch = 0
    for ev in events:
        d = np.load(ev)
        hc_poi = d["hc_poi"]
        th = d["thresholds"]
        alert_values = d["alert_values"]
        alert_levels = d["alert_levels"]
        # monotonicity audit of the real curves (np.interp validity)
        n_nonmono = sum(1 for ip in range(hc_poi.shape[0])
                        if np.any(np.diff(hc_poi[ip]) > 0))
        bump()
        if n_nonmono:
            NOTES.append(f"{ev}: {n_nonmono} POIs with non-monotone hc rows "
                         f"(np.interp misuse risk in float lane)")
        ev_widths = []
        for ip in range(hc_poi.shape[0]):
            hc = hc_poi[ip]
            hc_sum_pos = float(hc.sum()) > 0
            hcl = [float(h) for h in hc]
            for meth in pr_methods:
                pr = C.pr_threshold(meth)
                total += 1
                bump()
                try:
                    v_float, flag = percentile_value(hc, th, pr)
                except AssertionError:
                    ball, status = cert_percentile_bracket(hc, th, pr)
                    if status != STATUS_INDET:
                        fail("semantics",
                             f"{ev} POI {ip} {meth}: float lane crash geometry "
                             f"but cert status {status}",
                             f"cert_percentile_bracket(hc_poi[{ip}], th, {pr})")
                    continue
                ball, status = cert_percentile_bracket(hc, th, pr, prec=128)
                stat_counts[status] = stat_counts.get(status, 0) + 1
                if status == STATUS_INDET:
                    fail("semantics", f"{ev} POI {ip} {meth}: cert INDET on "
                         f"real curve",
                         f"cert_percentile_bracket(hc_poi[{ip}], th, {pr})")
                    continue
                if status == STATUS_AMBIGUOUS:
                    fail("semantics",
                         f"{ev} POI {ip} {meth}: ambiguous on degenerate inputs",
                         f"cert_percentile_bracket(hc_poi[{ip}], th, {pr})")
                lo, hi = ball_bounds(ball)
                # (A) the primitives' OWN claim: exact convention value (dps=80
                #     oracle on the same float64 inputs) inside the bracket
                v_true, tflag = oracle_percentile(hcl, th, pr)
                if v_true is None:
                    fail("semantics", f"{ev} POI {ip} {meth}: oracle degenerate "
                         f"({tflag}) on real curve", f"POI {ip} {meth} {ev}")
                elif not (lo <= v_true <= hi):
                    n_oracle_outside += 1
                    fail("semantics",
                         f"FATAL {ev} POI {ip} {meth}: TRUE convention value "
                         f"outside bracket",
                         f"cert_percentile_bracket(np.load('{ev}')['hc_poi']"
                         f"[{ip}], thresholds, {pr}) vs mpmath dps=80")
                # (B) the literal task-spec claim: float-lane value inside
                vf = mp.mpf(v_float)
                if lo <= vf <= hi:
                    n_float_inside += 1
                else:
                    n_float_outside += 1
                    gap = max(lo - vf, vf - hi)
                    du = float(gap) / _ulp(v_float)
                    ulp_dists.append(du)
                    if du > worst_ulp[0]:
                        worst_ulp = (du, (ev, ip, meth, v_float,
                                          mp.nstr(lo, 25), mp.nstr(hi, 25)))
                    # is the float lane at least the correctly-rounded double
                    # of the bracketed truth?
                    mid = (lo + hi) / 2
                    if float(mp.mpf(mid)) != v_float:
                        n_not_nearest += 1
                # external reference versions label this branch either way
                if (status == STATUS_GT100M_BRANCH) != (
                        flag in ("gt100m-branch", "val>100m-bug")):
                    fail("semantics",
                         f"{ev} POI {ip} {meth}: >100 m-branch flag mismatch "
                         f"cert={status} float={flag}",
                         f"POI {ip} {meth} in {ev}")
                n_gt += status == STATUS_GT100M_BRANCH
                ev_widths.append(hi - lo)
                # cross-check stored pipeline values (only stored when hc>0)
                if hc_sum_pos:
                    stored = float(alert_values[ip, pr_cols[meth]])
                    if stored != v_float:
                        stored_mismatch += 1
                        fail("semantics",
                             f"{ev} POI {ip} {meth}: replayed float {v_float!r}"
                             f" != stored alert_values {stored!r}",
                             f"POI {ip} {meth} in {ev}")
                    lvl = cert_alert_level(ball)
                    stored_lvl = int(alert_levels[ip, pr_cols[meth]])
                    if lvl == ALERT_STRADDLES:
                        n_straddle += 1
                    elif lvl != stored_lvl:
                        level_mismatch += 1
                        fail("semantics",
                             f"{ev} POI {ip} {meth}: cert level {lvl} != stored"
                             f" {stored_lvl}", f"POI {ip} {meth} in {ev}")
                else:
                    n_zero += 1
        all_widths.extend(ev_widths)
        wme = sorted(ev_widths)[len(ev_widths) // 2] if ev_widths else None
        wmx = max(ev_widths) if ev_widths else None
        print(f"  {os.path.basename(ev)}: {hc_poi.shape[0]} POIs x 8 levels, "
              f"median width {mp.nstr(wme, 3)}, max {mp.nstr(wmx, 3)}, "
              f"nonmono rows {n_nonmono}")
    ws = sorted(all_widths)
    med = ws[len(ws) // 2]
    mx = ws[-1]
    if n_float_outside:
        ud = sorted(ulp_dists)
        fail("semantics-spec",
             f"literal task-spec check: float-lane value outside certified "
             f"bracket in {n_float_outside}/{total} POIxPr cases across both "
             f"events ({n_float_inside} inside). The brackets (median width "
             f"{mp.nstr(med, 3)}) enclose the EXACT convention value at "
             f"prec=128, not the float64-evaluated np.interp output, which "
             f"sits O(1 ulp) away. ulp-distance from bracket: median "
             f"{ud[len(ud)//2]:.3f}, max {ud[-1]:.3f} at {worst_ulp[1]}; "
             f"float lane == nearest-double(bracket mid) in "
             f"{n_float_outside - n_not_nearest}/{n_float_outside} of the "
             f"outside cases ({n_not_nearest} not correctly rounded).",
             "see attack_semantics(); e.g. "
             "cert_percentile_bracket(np.load(<event npz>)['hc_poi'][1], "
             "thresholds, 0.01) vs the float64 reference percentile_value "
             "on the same row")
    NOTES.append(f"semantics: {total} brackets, statuses {stat_counts}, "
                 f">100 m branch {n_gt}, all-zero curves {n_zero}, "
                 f"alert straddles {n_straddle}, stored-val mismatches "
                 f"{stored_mismatch}, cert-vs-stored level mismatches "
                 f"{level_mismatch}, TRUE-value-outside {n_oracle_outside}; "
                 f"width median {mp.nstr(med, 3)}, max {mp.nstr(mx, 3)}; "
                 f"{time.perf_counter()-t0:.1f}s")
    print(f"  total {total} brackets; width median {mp.nstr(med, 3)} "
          f"max {mp.nstr(mx, 3)}; float-inside {n_float_inside}, "
          f"float-outside {n_float_outside}; TRUE-outside {n_oracle_outside}; "
          f"{time.perf_counter()-t0:.1f}s")


# ===========================================================================
# 4. PRECISION LAWS
# ===========================================================================

def attack_precision_laws():
    print("== precision laws ==")
    # (a) ctx.prec restored after normal use and after in-context exception
    before = ctx.prec
    cert_lognormal_survival(0.5, -1.0, prec=192)
    bump()
    if ctx.prec != before:
        fail("precision", "ctx.prec not restored after primitive call",
             "cert_lognormal_survival(0.5,-1.0,prec=192); check ctx.prec")
    try:
        cert_normalized_weighted_sum([1.0, 1.0], [0.5, 0.5], values=[10.0, 20.0],
                                     prec=192)
    except ValueError:
        pass
    bump()
    if ctx.prec != before:
        fail("precision", "ctx.prec leaked after in-context ValueError",
             "trigger tightener ValueError at prec=192; check ctx.prec")

    # (b) decimal parse honors working precision (radius shrinks with prec)
    with working_prec(64):
        r64 = arb("0.1").rad()
    with working_prec(256):
        r256 = arb("0.1").rad()
    bump()
    if not (float(r256) < float(r64)):
        fail("precision", "pin('0.1') radius does not scale with working prec",
             "compare arb('0.1').rad() at prec 64 vs 256")

    # (c) float64 pinning is exact for doubles
    with working_prec(128):
        b = pin(0.1)
    bump()
    if not b.is_exact():
        fail("precision", "pin(float) not exact", "pin(0.1).is_exact()")

    # (d) pin() routes ints through float(): silent rounding above 2^53
    big = 2 ** 53 + 1
    with working_prec(128):
        bb = pin(big)
    bump()
    if not contains_mpf(bb, mp.mpf(big)):
        fail("precision-law-minor",
             "pin(int) goes through float(): pin(2**53+1) is an exact ball that "
             "does NOT contain the integer 2**53+1 (law-1 tension: constructed "
             "through a 53-bit float while the intended real differs; arb can "
             "represent the int exactly). Not reachable from the shipped "
             "float64 pipeline data.",
             "python3 -c \"from certlane.primitives import pin, "
             "working_prec; import mpmath as mp;\\nwith working_prec(128): "
             "b=pin(2**53+1);\\nprint(b)\"  # equals 2**53, not 2**53+1")

    # (e) boundary_semantics flag: float64 pin exact, decimal ball tightens
    with working_prec(128):
        from certlane.primitives import _boundary_ball
        bf = _boundary_ball(0.10, "float64")
        bd = _boundary_ball(0.10, "decimal")
        bump()
        if not bf.is_exact():
            fail("precision", "'float64' boundary not exact double",
                 "_boundary_ball(0.10,'float64')")
        # decimal ball must contain the true rational 1/10 and exclude the double
        if not contains_mpf(bd, mp.mpf(1) / 10):
            fail("precision", "'decimal' boundary ball misses rational 1/10",
                 "_boundary_ball(0.10,'decimal')")
        if contains_mpf(bd, mp.mpf(0.1)):   # the double, 2^-128 ball excludes it
            fail("precision", "'decimal' boundary ball contains the double 0.1 "
                 "(radius too large: semantics indistinguishable)",
                 "_boundary_ball(0.10,'decimal')")
        bump()
    # (f) INDET is a valid enclosure of anything and never certifies
    ind = INDET()
    bump()
    if not is_indet(ind):
        fail("precision", "INDET() not detected by is_indet", "is_indet(INDET())")
    if cert_alert_level(ind) != ALERT_INDET:
        fail("precision", "INDET ball certifies an alert level",
             "cert_alert_level(INDET())")
    print("  probes done")


# ===========================================================================

def main():
    if len(sys.argv) >= 3 and sys.argv[1] == "--corruption":
        _corruption_child(sys.argv[2])
        return
    t0 = time.perf_counter()
    attack_survival()
    attack_weighted_sum()
    attack_normalized()
    attack_percentile()
    attack_alert()
    attack_planted()
    attack_semantics()
    attack_precision_laws()
    dt = time.perf_counter() - t0
    fatal = [f for f in FAILURES if f["surface"] in ("containment", "semantics",
                                                     "semantics-spec",
                                                     "planted", "precision")]
    print(f"\nTOTAL checks {CHECKS}, failures {len(FAILURES)} "
          f"(fatal-class {len(fatal)}), {dt:.1f}s")
    for f in FAILURES:
        print("  -", f["surface"], "|", f["desc"])
    for nline in NOTES:
        print("  note:", nline)
    print("ATTACK_JSON:" + json.dumps({
        "checks_run": CHECKS,
        "failures": FAILURES,
        "notes": NOTES,
        "elapsed_s": round(dt, 1),
    }))
    sys.exit(1 if fatal else 0)


if __name__ == "__main__":
    main()

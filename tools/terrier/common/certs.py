#!/usr/bin/env python3
"""certs.py — terrier common chassis: two-dps agreement, ball honesty,
Arb-ball trim discipline.

Non-asserting certificate gates (two-dps, ball_gate, route_gate) and the
suite ball law:
- NEVER reduce certified balls through float64 (subnormal cliff: midpoint
  keeps a few bits, radius underflows to 0.0 and gates pass VACUOUSLY);
- the radius is QUOTED with every verdict (ball-honesty law);
- published-value gate = published in [ball +/- their-last-digit rounding];
- two-precision agreement does NOT catch shared systematics — it is one gate
  among several, never the whole verdict.
All arithmetic here is mpmath at DPS >= 250; floats never touch a ball.

DELEGATION: the ball PRIMITIVES parse_ball and matched_digits are
verdict.py's — imported here, not re-implemented. Split of surfaces:
verdict.py = asserting print-scoreboard gates; certs.py = non-asserting
dict/receipt gates + the honesty/trim/published-value laws verdict.py does
not carry.
"""
import os, sys

import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from verdict import matched_digits, parse_ball  # noqa: E402,F401

DPS = 250


def is_ball_str(s):
    """Trim/serialization discipline: a certified value must be shipped as an
    explicit ball string ('[mid +/- rad]', or '[mid]' exact) — a bare decimal
    is NOT a certified quantity."""
    s = s.strip()
    return s.startswith("[") and s.endswith("]")


def float64_safe(ball_s):
    """Trim-discipline guard (arb-to-float64 pitfall): False if reducing this ball
    through float64 would lie — radius underflows to 0.0 while nonzero, or the
    midpoint lands in the subnormal band (< ~2.3e-308) keeping only a few bits.
    Consumers must stay in ball arithmetic whenever this is False."""
    m, r = parse_ball(ball_s)
    with mp.workdps(DPS):
        if r > 0 and float(r) == 0.0:
            return False
        if m != 0 and abs(m) < mp.mpf("2.3e-308"):
            return False
    return True


def two_dps_gate(low_s, high_s, bar_digits):
    """c1 pattern: low-dps vs high-dps midpoints agree to >= bar digits.
    Catches precision-dependent bugs ONLY — shared systematics pass both."""
    d = matched_digits(parse_ball(low_s)[0], parse_ball(high_s)[0])
    return {"gate": "two-dps", "matched_digits": d, "bar": bar_digits,
            "pass": d >= bar_digits}


def ball_gate(mine_s, ref_s, bar_digits, max_rad_ratio=16.0):
    """same-route regression vs a reference ball: overlap + matched midpoint
    digits + comparable radius (dict form of verdict.ball_gate)."""
    m1, r1 = parse_ball(mine_s)
    m2, r2 = parse_ball(ref_s)
    with mp.workdps(DPS):
        ov = bool(abs(m1 - m2) <= r1 + r2)
        rat = float(r1 / r2) if r2 > 0 else (1.0 if r1 == 0 else float("inf"))
    d = matched_digits(m1, m2)
    return {"gate": "ball", "overlap": ov, "matched_digits": d, "bar": bar_digits,
            "radius_ratio": rat, "pass": ov and d >= bar_digits
            and rat <= max_rad_ratio}


def route_gate(a_s, b_s, bar_digits):
    """cross-ROUTE agreement: overlap + matched digits,
    NO radius-comparability bar — an independent route legitimately carries a
    wider certified ball (a C/R radius ratio ~2e9 is legitimate); like-for-like radius
    regression is ball_gate vs the SAME-route bank."""
    m1, r1 = parse_ball(a_s)
    m2, r2 = parse_ball(b_s)
    with mp.workdps(DPS):
        ov = bool(abs(m1 - m2) <= r1 + r2)
    d = matched_digits(m1, m2)
    return {"gate": "route", "overlap": ov, "matched_digits": d,
            "bar": bar_digits, "pass": ov and d >= bar_digits}


def rel_radius(ball_s):
    """rad/|mid| as an mpf string — computed at DPS, never in float64."""
    m, r = parse_ball(ball_s)
    with mp.workdps(DPS):
        assert m != 0, "rel_radius undefined at mid=0 — quote absolute radius"
        return mp.nstr(r / abs(m), 6)


def ball_honesty(ball_s, rel_bar):
    """Radius-reporting law: the radius is QUOTED with every verdict,
    and the relative radius must clear the stated bar. rel_bar is a decimal
    string (e.g. '1e-140') so the bar itself never passes through float64."""
    assert is_ball_str(ball_s), "certified value must be a ball string: %r" % ball_s
    m, r = parse_ball(ball_s)
    with mp.workdps(DPS):
        rr = r / abs(m)
        ok = bool(rr < mp.mpf(rel_bar))
        return {"gate": "ball-honesty", "quoted_ball": ball_s,
                "rel_radius": mp.nstr(rr, 6), "bar": rel_bar,
                "float64_safe": float64_safe(ball_s), "pass": ok}


def last_digit_ulp(published_s):
    """decimal weight of the last printed digit of a published value
    (e.g. '2.037e-8' -> 1e-11; '6.45964e-62' -> 1e-67)."""
    s = published_s.strip().lower().lstrip("+-")
    mant, _, ex = s.partition("e")
    exp = int(ex) if ex else 0
    frac = len(mant.partition(".")[2])
    with mp.workdps(DPS):
        return mp.mpf(10) ** (exp - frac)


def published_value_gate(ball_s, published_s):
    """H4 comparison law: PASS iff the published value lies within
    [ball +/- half its own last-digit rounding]. Ball width is reported with
    the verdict; a float match is a convention LOCK, not a confirmation."""
    m, r = parse_ball(ball_s)
    with mp.workdps(DPS):
        pub = mp.mpf(published_s)
        band = last_digit_ulp(published_s) / 2
        dist = abs(pub - m)
        ok = bool(dist <= r + band)
        return {"gate": "published-value", "published": published_s,
                "quoted_ball": ball_s, "rounding_band": mp.nstr(band, 3),
                "distance": mp.nstr(dist, 4), "pass": ok}

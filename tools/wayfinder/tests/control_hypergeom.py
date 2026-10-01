#!/usr/bin/env python3
"""
control_hypergeom.py — END-TO-END POSITIVE CONTROL for wayfinder.

Runnable standalone:  python3 tests/control_hypergeom.py     (~6-7 min)
Exit code 0 = PASS, 1 = FAIL.  Every reported digit is MEASURED against an
independent mpmath oracle — nothing is fabricated, nothing is softened.

WHAT IT EXERCISES (all layers together, on a system with a known closed form)
=============================================================================
The Gauss hypergeometric function w = 2F1(a,b;c;x) with eps-dependent
parameters
        a = eps,   b = -2*eps,   c = 1 - 3*eps
satisfies  x(1-x) w'' + [c - (a+b+1)x] w' - ab w = 0,  i.e. the first-order
system  y' = A(x, eps) y  with  y = (w, w'):

        A = [[0, 1],
             [ab/(x(1-x)),  ((a+b+1)x - c)/(x(1-x))]]

regular-singular at x = 0, 1 (the classic test case for endpoint transport:
the straight path 1/2 -> -1 passes THROUGH the singular point x = 0, so a
complex detour waypoint is mandatory — the MUM-crossing detour pattern).

  stage A  de_load-style eval: DESystem built inline against the shared
           contract (.n, .var, .A(x, eps, dps)); plus the documented
           duck-typed .A_series fast path (exact partial-fraction Taylor),
           cross-checked order-0 against pointwise A.
  stage B  transport.transport_fixed_eps at FIXED eps0 = 1/10 from x0 = 1/2
           (seeded from mpmath.hyp2f1 at the anchor) to x1 = -1 via the
           detour waypoint -0.6i, at dps 100 AND dps 140; landed values
           gated (gate.gate_strings, two-precision rule) against 90 stored
           digits of mpmath.hyp2f1(0.1, -0.2; 0.7; -1) — an oracle value the
           transport never saw (only the x = +1/2 anchor was consumed).
  stage C  epsfan.cauchy_laurent over eps of the FULL anchor->transport->land
           pipeline (default calibrated R = 1e-3, 32 nodes, measured R2
           second-pass diagnostics), window k = -2..6; the k = -2, -1
           coefficients must be zero-consistent (the function is analytic in
           eps — a nonzero pole coefficient means the fan or transport lies).
  stage D  independent eps-expansion oracle: mpmath.taylor of
           g(eps) = hyp2f1(eps, -2*eps; 1-3*eps; -1) at two working
           precisions (130/170), with the ANALYTIC anchor
           c2 = -2*Li2(-1) = pi^2/6 (exact) proving the oracle itself.
           Fan coefficients c_0..c_6 must match the taylor oracle.

PASS bars (all must hold):
  * gate_strings overall PASS at the 90-stored-digit bar (both components,
    matched >= 89 vs the full stored string, pair agreement >= 90),
  * fast-path (A_series) transport matches generic-path transport >= 95d,
  * fan min_agree_digits >= 30 (measured (R, 2R) two-pass agreement),
  * |c_-2|, |c_-1| < 10^-40 * gscale (zero-consistent),
  * fan c_k vs mpmath.taylor matched digits >= 30 for every k = 0..6,
  * fan c_2 vs EXACT pi^2/6 matched digits >= 30,
  * global mp.dps untouched at exit (import-dps footgun tripwire).

The >= 30 digit bar is the standing solved-bar floor; the measured numbers
(reported below) should land far above it — treat a barely-passing run as a
regression even if this script says PASS.
"""

import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from mpmath import mp, mpf, mpc  # noqa: E402

from transport import transport_fixed_eps  # noqa: E402
from epsfan import cauchy_laurent          # noqa: E402
from gate import gate_strings              # noqa: E402


# ---------------------------------------------------------------------------
# inline DESystem (de_load contract: .n, .var, .A(x, eps, dps), .meta)
# ---------------------------------------------------------------------------
class Hyp2F1System:
    """2x2 companion system of the hypergeometric ODE, a=eps, b=-2eps,
    c=1-3eps.  Pointwise A only (exercises transport's generic Cauchy-circle
    local-series path)."""
    n = 2
    var = "x"
    meta = {"source": "inline control_hypergeom.py", "format": "inline-2f1",
            "singular_points": [0, 1]}

    def A(self, x, eps, dps):
        with mp.workdps(dps):
            xv = mpc(x)
            e = mpc(eps)
            # a+b+1 = 1 - eps,  c = 1 - 3*eps,  ab = -2*eps^2
            den = xv * (1 - xv)
            return [[mpc(0), mpc(1)],
                    [-2 * e * e / den, ((1 - e) * xv - (1 - 3 * e)) / den]]


class Hyp2F1SystemFast(Hyp2F1System):
    """Same system + the documented duck-typed A_series fast path.

    Exact partial fractions:  1/(x(1-x)) = 1/x + 1/(1-x),
        A[1][0] = -2 eps^2 (1/x + 1/(1-x))
        A[1][1] = -(1-3 eps)/x + 2 eps/(1-x)
    Taylor about z0:  1/x -> (-1)^k / z0^(k+1),  1/(1-x) -> 1/(1-z0)^(k+1).
    """

    def A_series(self, z0, eps, dps, M):
        with mp.workdps(dps + 10):
            z0v = mpc(z0)
            e = mpc(eps)
            iz = 1 / z0v
            iw = 1 / (1 - z0v)
            c10 = -2 * e * e
            aX = -(1 - 3 * e)
            bX = 2 * e
            out = []
            pz, pw, sgn = iz, iw, 1
            for k in range(M + 1):
                fx = sgn * pz          # k-th Taylor coeff of 1/x about z0
                fw = pw                # k-th Taylor coeff of 1/(1-x) about z0
                out.append([[mpc(0), mpc(1) if k == 0 else mpc(0)],
                            [c10 * (fx + fw), aX * fx + bX * fw]])
                pz *= iz
                pw *= iw
                sgn = -sgn
            return out


# path 1/2 -> -1 must NOT go through x=0: detour below the real axis
# (stays inside the principal-branch domain C \ [1, inf), which is simply
# connected, so mpmath.hyp2f1's principal branch is the right oracle).
_DETOUR = [mpc(0, "-0.6")]
X0, X1 = "0.5", "-1"


def _seed(e, wp):
    """Anchor vector (w, w') at x = 1/2 from mpmath.hyp2f1 (the ONLY oracle
    input the transport consumes; the x = -1 oracle is held out)."""
    with mp.workdps(wp):
        a, b, c = mpc(e), -2 * mpc(e), 1 - 3 * mpc(e)
        h = mp.mpf(1) / 2
        w = mp.hyp2f1(a, b, c, h)
        wp_ = a * b / c * mp.hyp2f1(a + 1, b + 1, c + 1, h)
        return [w, wp_]


def _agree(a, b, cap):
    """Matched significant digits (measured, capped)."""
    a, b = mpc(a), mpc(b)
    diff = abs(a - b)
    scale = max(abs(a), abs(b))
    if diff == 0:
        return float(cap)
    if scale == 0:
        return 0.0
    return max(0.0, min(float(-mp.log10(diff / scale)), float(cap)))


def main():
    t_start = time.time()
    dps_ambient_before = mp.dps
    failures = []
    print("=" * 74)
    print("wayfinder END-TO-END POSITIVE CONTROL: 2F1(eps,-2eps;1-3eps;x)")
    print("  transport x: 1/2 -> -0.6i -> -1   (x=0 singular, detour forced)")
    print("=" * 74)

    generic = Hyp2F1System()
    fast = Hyp2F1SystemFast()

    # ---- stage A: contract + fast-path consistency ------------------------
    with mp.workdps(80):
        z0 = mpc("0.3", "-0.2")
        ep = mpc("0.07", "0.011")
        A0 = generic.A(z0, ep, 70)
        S = fast.A_series(z0, ep, 70, 4)
        errA = max(abs(mpc(A0[i][j]) - mpc(S[0][i][j]))
                   for i in range(2) for j in range(2))
        okA = errA < mpf(10) ** -65
    print(f"[A] A_series order-0 vs pointwise A: max err = {mp.nstr(errA, 3)} "
          f"-> {'ok' if okA else 'FAIL'}")
    if not okA:
        failures.append("stage A: A_series[0] != A pointwise")

    # ---- stage B: fixed-eps transport at dps 100 and 140, gate at 90d -----
    eps0 = "0.1"
    landed = {}
    for dps in (100, 140):
        y0 = _seed(mpf(1) / 10, dps + 30)
        t0 = time.time()
        y = transport_fixed_eps(generic, eps0, X0, X1, y0, dps, path=_DETOUR)
        dt = time.time() - t0
        landed[dps] = y
        # imag parts must be pure roundoff (everything is real at real eps,
        # real endpoints; the detour is a contour choice, not a branch change)
        with mp.workdps(dps + 10):
            im_rel = max(abs(mpc(v).imag) / max(abs(mpc(v)), mpf(1))
                         for v in y)
        print(f"[B] transport dps={dps}: {dt:.1f}s, |Im|/|y| <= "
              f"{mp.nstr(im_rel, 3)}")
        if im_rel > mpf(10) ** (-(dps - 10)):
            failures.append(f"stage B: imaginary contamination at dps={dps}")

    # held-out oracle at x=-1 (never consumed by the transport)
    with mp.workdps(200):
        e = mpf(1) / 10
        a, b, c = e, -2 * e, 1 - 3 * e
        ref1 = mp.hyp2f1(a, b, c, -1)
        ref2 = a * b / c * mp.hyp2f1(a + 1, b + 1, c + 1, -1)
        vendored = {"y1_re": mp.nstr(ref1.real if hasattr(ref1, "real") else ref1, 90),
                    "y2_re": mp.nstr(ref2.real if hasattr(ref2, "real") else ref2, 90)}
    computed = {}
    for key, comp in (("y1_re", 0), ("y2_re", 1)):
        pair = []
        for dps in (100, 140):
            with mp.workdps(dps + 5):
                pair.append(mp.nstr(mpc(landed[dps][comp]).real, dps))
        computed[key] = tuple(pair)
    report = gate_strings(computed, vendored, (100, 140))
    for key, rec in report["per_key"].items():
        print(f"[B] gate {key}: stored={rec['stored_digits']}d  "
              f"matched={rec['matched_digits']}  "
              f"pair_agree={rec['pair_agree_digits']}  "
              f"{'PASS' if rec['PASS'] else 'FAIL: ' + str(rec['reason'])}")
    print(f"[B] gate overall vs held-out mpmath.hyp2f1(x=-1): "
          f"{report['overall']}  (bar: 90 stored digits, two-precision rule)")
    if report["overall"] != "PASS":
        failures.append("stage B: gate FAIL vs held-out hyp2f1 oracle")

    # fast path must reproduce the generic path
    y0 = _seed(mpf(1) / 10, 130)
    yf = transport_fixed_eps(fast, eps0, X0, X1, y0, 100, path=_DETOUR)
    with mp.workdps(120):
        d_fg = min(_agree(yf[i], landed[100][i], 110) for i in range(2))
    print(f"[B] A_series fast path vs generic path @dps100: {d_fg:.1f} "
          f"matched digits (bar >= 95)")
    if d_fg < 95:
        failures.append(f"stage B: fast path vs generic only {d_fg:.1f}d")

    # ---- stage C: eps-fan of the full pipeline ----------------------------
    DPS_FAN = 60
    KMIN, KMAX = -2, 6
    n_evals = [0]

    def f(eps):
        n_evals[0] += 1
        wp_here = mp.dps                # inside cauchy_laurent's workdps
        y0e = _seed(eps, wp_here + 20)
        ye = transport_fixed_eps(fast, eps, X0, X1, y0e, wp_here + 5,
                                 path=_DETOUR)
        return ye[0]                    # w(-1; eps)

    t0 = time.time()
    fan = cauchy_laurent(f, KMIN, KMAX, DPS_FAN, R="1e-3", nodes=32,
                         check="R2")
    dt = time.time() - t0
    diag = fan["diag"]
    print(f"[C] eps-fan: {n_evals[0]} transports in {dt:.0f}s "
          f"(R={diag['R']}, nodes={diag['nodes']}, wp={diag['wp']})")
    print(f"[C] fan measured min_agree_digits = "
          f"{diag['min_agree_digits']:.1f} (bar >= 30)")
    if diag["min_agree_digits"] < 30:
        failures.append(f"stage C: fan min_agree "
                        f"{diag['min_agree_digits']:.1f} < 30")
    # pole coefficients must be zero-consistent (analytic in eps)
    with mp.workdps(DPS_FAN):
        gscale = mp.mpmathify(diag["gscale"])
        for k in (-2, -1):
            ck = abs(mpc(fan["coeffs"][k]))
            rel = ck / gscale
            print(f"[C] |c_{k}| / gscale = {mp.nstr(rel, 3)} "
                  f"(zero-consistency bar < 1e-40)")
            if rel > mpf(10) ** -40:
                failures.append(f"stage C: spurious eps^{k} pole, "
                                f"|c|/gscale = {mp.nstr(rel, 3)}")

    # ---- stage D: independent taylor oracle + exact anchor ----------------
    def g(e):
        return mp.hyp2f1(e, -2 * e, 1 - 3 * e, -1)

    tay = {}
    for wdps in (130, 170):
        with mp.workdps(wdps):
            tay[wdps] = mp.taylor(g, 0, KMAX)
    with mp.workdps(180):
        # taylor's finite differences return NOISE (not exact 0) for the
        # analytically-zero c_1: treat |c| below 10^-40 of the oracle scale
        # as zero-consistent, both here and in the per-k comparison below.
        t_scale = max(abs(mpc(v)) for v in tay[170])
        zfloor_o = t_scale * mpf(10) ** -40
        oracle_digits = min(
            _agree(tay[130][k], tay[170][k], 170) for k in range(KMAX + 1)
            if abs(mpc(tay[170][k])) > zfloor_o)
        exact_c2 = mp.pi ** 2 / 6      # c2 = -2*Li2(-1), analytic anchor
        anchor_d = _agree(tay[170][2], exact_c2, 170)
    print(f"[D] taylor oracle self-agreement (130/170): "
          f"{oracle_digits:.1f}d; c2 vs EXACT pi^2/6: {anchor_d:.1f}d")
    if anchor_d < 100:
        failures.append(f"stage D: taylor oracle broken — c2 vs pi^2/6 "
                        f"only {anchor_d:.1f}d")

    with mp.workdps(DPS_FAN + 20):
        print(f"[D] fan c_k vs independent mpmath.taylor oracle "
              f"(cap {DPS_FAN}d):")
        worst = float("inf")
        for k in range(0, KMAX + 1):
            ck_abs = abs(mpc(fan["coeffs"][k]))
            ok_abs = abs(mpc(tay[170][k]))
            fan_zero = ck_abs < gscale * mpf(10) ** -40
            orc_zero = ok_abs < zfloor_o
            if fan_zero or orc_zero:
                # zero-consistency regime (analytically-zero c_1): both
                # sides must be zero-consistent, else 0 matched digits
                d = float(DPS_FAN) if (fan_zero and orc_zero) else 0.0
            else:
                d = _agree(fan["coeffs"][k], tay[170][k], DPS_FAN)
            print(f"      k={k}:  {d:6.1f} matched digits   "
                  f"c_k = {mp.nstr(mpc(fan['coeffs'][k]).real, 15)}")
            worst = min(worst, d)
            if d < 30:
                failures.append(f"stage D: c_{k} only {d:.1f}d vs oracle")
        d_c2 = _agree(fan["coeffs"][2], exact_c2, DPS_FAN)
        print(f"[D] fan c_2 vs EXACT pi^2/6: {d_c2:.1f} matched digits")
        if d_c2 < 30:
            failures.append(f"stage D: c_2 vs exact pi^2/6 only {d_c2:.1f}d")

    # ---- import-dps footgun tripwire ---------------------------------------
    if mp.dps != dps_ambient_before:
        failures.append(f"global mp.dps mutated: {dps_ambient_before} -> "
                        f"{mp.dps}")

    print("-" * 74)
    verdict = "PASS" if not failures else "FAIL"
    print(f"CONTROL VERDICT: {verdict}   ({time.time() - t_start:.0f}s total)")
    for fmsg in failures:
        print("  FAIL:", fmsg)
    return 0 if not failures else 1


if __name__ == "__main__":
    sys.exit(main())

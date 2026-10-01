#!/usr/bin/env python3
"""
Tests for tools/wayfinder/boundary.py.

B1  tadpole plumbing      — string/mpc arg paths reproduce -Gamma(1-d/2) at
                            complex eps; rounding at dps.
B2  tadpole recursion     — INDEPENDENT check: T(nu+1) = (1/nu) dT/dm^2 via
                            mpmath central difference (no shared formula).
B3  vacuum table          — (1,1) == tadpole(1,1,d); (2,3) == massless-bubble
                            x single_mass_prefactor (the AMFlow SingleMass
                            loop-promotion identity, derived independently).
B4  TS_CLOSED gate        — tadpole(1,1,d)^3 Laurent coefficients via
                            epsfan.cauchy_laurent gated against the vendored
                            60d strings of the reference TS_CLOSED.json
                            fixture with gate.gate_strings (two-
                            precision).  Skipped if the file is absent.
B5  vacuum_ending_seed    — Tradition leaf rules, eta-power premultiply,
                            SingleMass == Tradition on the pinned (2,3)
                            overlap, and every error/NotImplemented path.

References are derived from mpmath closed forms in-test — no fabricated
digits.  Global mp.dps is asserted untouched.
"""
import json
import os
import sys
import types

import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
sys.path.insert(0, PKG)

import boundary  # noqa: E402
import epsfan    # noqa: E402
import gate      # noqa: E402

# Optional reference fixture, not shipped (leg SKIPs when absent).
TS_CLOSED = os.path.join(os.environ.get("WAYFINDER_REFERENCE_FIXTURES", ""),
                         "TS_CLOSED.json")

FAILS = []


def check(name, cond, msg=""):
    tag = "PASS" if cond else "FAIL"
    print("  [%s] %s%s" % (tag, name, (": " + msg) if msg else ""))
    if not cond:
        FAILS.append(name)


def agree_d(a, b, wp):
    with mp.workdps(wp):
        a, b = mp.mpc(a), mp.mpc(b)
        diff = abs(a - b)
        if diff == 0:
            return float(wp)
        s = max(abs(a), abs(b))
        return float(-mp.log10(diff / s)) if s > 0 else float(wp)


DPS0 = mp.mp.dps  # global-dps guard


# ---------------------------------------------------------------------------
# B1  tadpole plumbing
# ---------------------------------------------------------------------------
def test_tadpole_plumbing():
    dps = 80
    with mp.workdps(dps + 20):
        eps = mp.mpf(1) / 7 + mp.mpc(0, 1) / 11
        d = 4 - 2 * eps
        ref = -mp.gamma(1 - d / 2)          # closed form, derived here
    t1 = boundary.tadpole("1", 1, d, dps)   # string m2, mpc d
    check("B1-value", agree_d(t1, ref, dps + 10) >= dps - 2,
          "agree=%.1f d" % agree_d(t1, ref, dps + 10))
    # nu validation
    for bad in (0, -1, 1.5):
        try:
            boundary.tadpole(1, bad, 4 - 2 * mp.mpf("0.1"), 40)
            check("B1-nu-reject-%r" % bad, False, "no exception")
        except ValueError:
            check("B1-nu-reject-%r" % bad, True)
    # eps=0 pole guard: nu - d/2 = 1 - 2 = -1 -> Gamma pole
    try:
        boundary.tadpole(1, 1, 4, 40)
        check("B1-pole-guard", False, "no exception at d=4 exactly")
    except ValueError as ex:
        check("B1-pole-guard", "nonzero-eps" in str(ex))


# ---------------------------------------------------------------------------
# B2  tadpole m^2-derivative recursion (independent of the Gamma formula)
# ---------------------------------------------------------------------------
def test_tadpole_recursion():
    # T(nu+1; m2) = (1/nu) * d/dm2 T(nu; m2)  — follows from differentiating
    # the INTEGRAND 1/(k^2-m^2)^nu, so it holds for any correct closed form.
    dps = 110
    nu = 2
    with mp.workdps(dps + 20):
        eps = mp.mpf(1) / 7
        d = 4 - 2 * eps
        m2 = mp.mpf(3) / 2
        h = mp.mpf(10) ** (-30)
        tp = boundary.tadpole(m2 + h, nu, d, dps + 10)
        tm = boundary.tadpole(m2 - h, nu, d, dps + 10)
        lhs = (tp - tm) / (2 * h) / nu
        rhs = boundary.tadpole(m2, nu + 1, d, dps + 10)
    ad = agree_d(lhs, rhs, dps)
    # central difference is O(h^2) accurate -> ~60 digits available
    check("B2-recursion", ad >= 45, "agree=%.1f d (need >=45)" % ad)


# ---------------------------------------------------------------------------
# B3  vacuum table vs tadpole and vs loop-promotion identity
# ---------------------------------------------------------------------------
def test_vacuum_table():
    dps = 90
    with mp.workdps(dps + 20):
        eps = mp.mpf(1) / 9 + mp.mpc(0, 1) / 13
    v11 = boundary.singlemass_vacuum(1, 1, eps, dps)
    with mp.workdps(dps + 20):
        d = 4 - 2 * eps
    t = boundary.tadpole(1, 1, d, dps)
    check("B3-11-eq-tadpole", agree_d(v11, t, dps + 10) >= dps - 2)

    # (2,3) via loop promotion: massless bubble at q^2 = -1, then the
    # single_mass_prefactor with (n=1, m0=0, m_eps_coef=1) exactly as
    # amfsystem.cpp derives it for total=3, n=1, loop_num=2:
    #   m0 = 3 - 1 - 2*(2-1) = 0,  m_eps_coef = 2-1 = 1.
    with mp.workdps(dps + 20):
        c_bub = (mp.gamma(eps) * mp.gamma(1 - eps) ** 2 / mp.gamma(2 - 2 * eps))
    pref = boundary.single_mass_prefactor(1, 0, 1, eps, dps)
    with mp.workdps(dps + 20):
        lhs = mp.mpc(c_bub) * mp.mpc(pref)
    v23 = boundary.singlemass_vacuum(2, 3, eps, dps)
    ad = agree_d(lhs, v23, dps + 10)
    check("B3-23-loop-promotion", ad >= dps - 3, "agree=%.1f d" % ad)

    check("B3-known", boundary.vacuum_known(3, 5) and
          not boundary.vacuum_known(2, 2))
    try:
        boundary.singlemass_vacuum(2, 2, eps, dps)
        check("B3-unknown-reject", False, "no exception")
    except ValueError:
        check("B3-unknown-reject", True)


# ---------------------------------------------------------------------------
# B4  TS_CLOSED.json cross-check (boundary + epsfan + gate together)
# ---------------------------------------------------------------------------
def test_ts_closed_gate():
    if not os.path.exists(TS_CLOSED):
        print("  [skip] B4: %s missing" % TS_CLOSED)
        return
    ts = json.load(open(TS_CLOSED))
    vendored = {k: ts["verify"]["T"][k]["closed"] for k in
                ("-3", "-2", "-1", "0")}

    def coeffs_at(dps):
        inner = dps + 40

        def f(e):
            return boundary.tadpole(1, 1, 4 - 2 * e, inner) ** 3

        # nodes=64: T = -Gamma(eps-1)^3 has a deep Laurent tail; at nodes=32
        # the diag MEASURES an alias floor of ~84 d (c_{k+32}*R^32 leakage,
        # 2R check pass) — the (r,delta) discipline working as designed.
        # 64 nodes push the fold-in to R^64 and restore the full dps.
        out = epsfan.cauchy_laurent(f, -3, 0, dps, nodes=64)
        check("B4-diag-honest-dps%d" % dps,
              out["diag"]["min_agree_digits"] >= dps - 10,
              "min_agree=%.1f" % out["diag"]["min_agree_digits"])
        with mp.workdps(dps):
            return {str(k): mp.nstr(out["coeffs"][k].real, dps - 5)
                    for k in out["coeffs"]}

    lo, hi = coeffs_at(90), coeffs_at(110)
    computed = {k: (lo[k], hi[k]) for k in vendored}
    rep = gate.gate_strings(computed, vendored, (90, 110))
    check("B4-gate-PASS", rep["overall"] == "PASS",
          json.dumps({k: v for k, v in rep["per_key"].items()
                      if not v["PASS"]}) if rep["overall"] != "PASS"
          else "all %d keys" % rep["n_keys"])


# ---------------------------------------------------------------------------
# B5  vacuum_ending_seed
# ---------------------------------------------------------------------------
def test_vacuum_ending_seed():
    dps = 60
    with mp.workdps(dps + 20):
        eps = mp.mpf(1) / 8
    fake = types.SimpleNamespace(n=3, var="eta")
    masters = [[-1, 0, 0], [0, 0, 0], [1, 1, 1]]
    seeds = boundary.vacuum_ending_seed(
        fake, "Tradition", {"eps": eps, "loops": 2, "masters": masters}, dps)
    check("B5-len", len(seeds) == 3)
    check("B5-unit-negidx", seeds[0] == 1 and seeds[1] == 1)
    ref23 = boundary.singlemass_vacuum(2, 3, eps, dps)
    check("B5-leaf-23", agree_d(seeds[2], ref23, dps) >= dps - 2)

    # eta-power premultiply: lambda = L*(2-eps) - sum(nu)
    with mp.workdps(dps + 20):
        eta = mp.mpc(0, -50)   # NegIm-side eta
    seeds_eta = boundary.vacuum_ending_seed(
        fake, "Tradition",
        {"eps": eps, "loops": 2, "masters": masters, "eta": eta}, dps)
    with mp.workdps(dps + 20):
        lam2 = 2 * (2 - eps) - 3
        ref = mp.mpc(ref23) * mp.power(eta, lam2)
        lam0 = 2 * (2 - eps) - (-1)
        ref0 = mp.power(eta, lam0)
    check("B5-eta-power", agree_d(seeds_eta[2], ref, dps) >= dps - 3 and
          agree_d(seeds_eta[0], ref0, dps) >= dps - 3)

    # SingleMass pinned overlap: prefactor(1,0,1) * bubble == Tradition (2,3)
    fake1 = types.SimpleNamespace(n=1, var="eta")
    with mp.workdps(dps + 20):
        c_bub = (mp.gamma(eps) * mp.gamma(1 - eps) ** 2
                 / mp.gamma(2 - 2 * eps))
    sm = boundary.vacuum_ending_seed(
        fake1, "SingleMass",
        {"eps": eps, "loops": 2, "masters": [[1, 1, 1]],
         "gamma_params": [(1, 0, 1)], "sub_values": [c_bub]}, dps)
    check("B5-singlemass-overlap", agree_d(sm[0], ref23, dps) >= dps - 3)

    # error / NotImplemented paths
    def expect(name, exc, fn):
        try:
            fn()
            check(name, False, "no exception")
        except exc:
            check(name, True)
        except Exception as ex:
            check(name, False, "wrong exception %r" % ex)

    P = {"eps": eps, "loops": 2, "masters": masters}
    expect("B5-bad-scheme", ValueError,
           lambda: boundary.vacuum_ending_seed(fake, "Trad", P, dps))
    expect("B5-len-mismatch", ValueError,
           lambda: boundary.vacuum_ending_seed(
               fake, "Tradition", {**P, "masters": [[1]]}, dps))
    expect("B5-modes-NIE", NotImplementedError,
           lambda: boundary.vacuum_ending_seed(
               fake, "Tradition", {**P, "modes": ["k1"]}, dps))
    expect("B5-offtable-NIE", NotImplementedError,
           lambda: boundary.vacuum_ending_seed(
               fake1, "Tradition",
               {"eps": eps, "loops": 2, "masters": [[1, 1, 0]]}, dps))
    expect("B5-dots-NIE", NotImplementedError,
           lambda: boundary.vacuum_ending_seed(
               fake1, "Tradition",
               {"eps": eps, "loops": 2, "masters": [[2, 1, 1]]}, dps))
    expect("B5-sm-missing-NIE", NotImplementedError,
           lambda: boundary.vacuum_ending_seed(
               fake1, "SingleMass",
               {"eps": eps, "loops": 2, "masters": [[1, 1, 1]]}, dps))
    expect("B5-sm-eta-reject", ValueError,
           lambda: boundary.vacuum_ending_seed(
               fake1, "SingleMass",
               {"eps": eps, "loops": 2, "masters": [[1, 1, 1]],
                "gamma_params": [(1, 0, 1)], "sub_values": [1],
                "eta": 10}, dps))


if __name__ == "__main__":
    print("== test_boundary ==")
    test_tadpole_plumbing()
    test_tadpole_recursion()
    test_vacuum_table()
    test_ts_closed_gate()
    test_vacuum_ending_seed()
    check("global-dps-untouched", mp.mp.dps == DPS0,
          "dps %d -> %d" % (DPS0, mp.mp.dps))
    print("== %s ==" % ("ALL PASS" if not FAILS else "FAILURES: " + ",".join(FAILS)))
    sys.exit(0 if not FAILS else 1)

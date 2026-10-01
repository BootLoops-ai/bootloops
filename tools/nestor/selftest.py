#!/usr/bin/env python3
r"""nestor.selftest -- positive + mutation controls for the NESTOR package.

Run:  python3 -m nestor.selftest [--receipt PATH] [--nproc N]

CONTROLS (every one must land; synthetic-truth controls are mandatory --
consistency alone is not a control):
  P1  closed form:      int_0^1 4/(1+x^2) dx = pi            (serial ladder)
  P2  endpoint sing.:   int_0^1 x^{-1/2} dx = 2, spec algebraic -1/2
                        (singcheck PASS + value)
  P3  both endpoints:   int_0^1 dx/sqrt(x(1-x)) = pi
  P4  nested two-fold:  int_0^1 dy int_0^1 dx 1/(1+xy) = pi^2/12
                        (inner nestor.integrate inside the outer integrand
                        -- the production nesting pattern)
  P5  farm smoke:       int_0^{1/2} e^u cos(3u) du farmed at md=6 vs the
                        serial ladder value AND the closed form
                        (the original farm driver's smoke integrand; exercises grid,
                        nesting selfcons, pool discipline)
  P6  log endpoint:     int_0^1 log(x) dx = -1, spec log-type at a
                        (the singcheck 'log' positive control)
  M1  singcheck mutation: declare x^{-1/2} endpoint as 'none'
                        -> MUST refuse (SpecRefusal)
  M2  non-integrable declaration: algebraic exponent -1
                        -> MUST refuse at declaration time
  M3  refusal receipt:  the M1 run's receipt must carry the named refusal
                        and NO value
  D1  dispersion moments: closed-form power-log moments vs adaptive
                        tanh-sinh, >= 30 digits          (nestor.dispersion)
  D2  dispersion add-back: closed-form int S K over a threshold sub-panel
                        vs direct quadrature, >= 25 digits
  D3  dispersion assembly: level-5 subtracted assembly on the two-threshold
                        surrogate vs a high-precision reference, >= 30 d
  D4  dispersion mutation: the add-back with the WRONG sign must destroy
                        the agreement (closed-form term is load-bearing)
Mutation controls dent the DECLARATION only (in-memory); disk is untouched.
The D legs are the functions in nestor/tests/test_dispersion_*.py, so the
pytest files and this command exercise the same assertions.

Exit codes follow the evaluators' contract: 0 pass, 3 gate-class failure.

"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import mpmath as mp

# package importable both as `python3 -m nestor.selftest` and as a script
if __package__ in (None, ""):                  # script mode
    sys.path.insert(0, os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))))
    import nestor
else:
    import importlib
    nestor = importlib.import_module(__package__)

from nestor.ladder import SpecRefusal, digits


# ---- farm-importable integrands (module-level: workers import by spec) ----
def pilot_integrand(u):
    """Analytic smoke integrand from the original farm driver: e^u cos(3u)."""
    return mp.exp(u) * mp.cos(3 * u)


def _check(name, cond, detail=""):
    status = "PASS" if cond else "FAIL"
    print("[selftest] %-4s %s %s" % (name, status, detail))
    return bool(cond)


def run(nproc=2, receipt_path=None):
    ok = []
    t00 = time.time()
    D = 30
    # reference comparisons need headroom above every control's dps
    mp.mp.dps = 70

    # P1: closed form pi
    v = nestor.integrate(lambda x: 4 / (1 + x * x), 0, 1, dps=D)
    d = digits(v, mp.pi)
    ok.append(_check("P1", d >= D, "4/(1+x^2) vs pi: %.1f d (bar %d)"
                     % (d, D)))

    # P2: endpoint singularity with correct spec
    spec = {"a": {"type": "algebraic", "exponent": -0.5}}
    v = nestor.integrate(lambda x: 1 / mp.sqrt(x), 0, 1, dps=D, spec=spec)
    d = digits(v, mp.mpf(2))
    ok.append(_check("P2", d >= D, "x^-1/2 vs 2: %.1f d (spec verified)"
                     % d))

    # P3: both endpoints singular -- complement-pair route: the integrand
    # takes the cancellation-free endpoint distances (d_a, d_b)
    spec = {"a": {"type": "algebraic", "exponent": -0.5},
            "b": {"type": "algebraic", "exponent": -0.5}}
    v = nestor.integrate(lambda u, da, db: 1 / mp.sqrt(da * db), 0, 1,
                         dps=D, spec=spec, wants_dists=True)
    d = digits(v, mp.pi)
    ok.append(_check("P3", d >= D, "1/sqrt(x(1-x)) vs pi: %.1f d" % d))

    # P4: nested two-fold (inner nestor call inside outer integrand).
    # wp-scaled ladder: inner runs at dps 26 => inner agreement bar 1e-36,
    # safely below the outer 10^-(20+10) refine tolerance.
    Din = 20

    def outer(y):
        return nestor.integrate(lambda x: 1 / (1 + x * y), 0, 1, dps=26)
    v = nestor.integrate(outer, 0, 1, dps=Din)
    d = digits(v, mp.pi ** 2 / 12)
    ok.append(_check("P4", d >= Din, "nested 1/(1+xy) vs pi^2/12: %.1f d "
                     "(bar %d)" % (d, Din)))

    # P5: farm smoke vs serial ladder + closed form
    ref = (mp.exp(mp.mpf(1) / 2) * (mp.cos(mp.mpf(3) / 2)
           + 3 * mp.sin(mp.mpf(3) / 2)) - 1) / 10   # int e^u cos3u, [0,1/2]
    pe = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    rec = nestor.oracle("nestor.selftest:pilot_integrand", 0, "0.5",
                        dps=D, md=6, nproc=nproc, f_bound=2,
                        path_entry=pe, log=lambda *a: None)
    vf = mp.mpf(rec["value"])
    vs = nestor.integrate(pilot_integrand, 0, mp.mpf(1) / 2, dps=D)
    df, ds = digits(vf, ref), digits(vf, vs)
    ok.append(_check("P5", df >= D and ds >= D,
                     "farm md6 vs closed form %.1f d, vs serial %.1f d, "
                     "selfcons %.1f d, %d nodes, wall %.1fs"
                     % (df, ds, rec["rungs"][0]["agreement_digits"],
                        rec["rungs"][0]["nodes"], rec["wall_s"])))

    # P6: log-type endpoint positive control:
    # int_0^1 log x dx = -1 with the spec declared log at a
    v = nestor.integrate(lambda x: mp.log(x), 0, 1, dps=D,
                         spec={"a": {"type": "log"}})
    d = digits(v, mp.mpf(-1))
    ok.append(_check("P6", d >= D, "log(x) vs -1: %.1f d (log spec "
                     "verified)" % d))

    # M1: singcheck mutation -- wrong declaration MUST be caught
    caught = False
    try:
        nestor.integrate(lambda x: 1 / mp.sqrt(x), 0, 1, dps=D,
                         spec={"a": {"type": "none"}})
    except SpecRefusal:
        caught = True
    ok.append(_check("M1", caught, "undeclared-singularity dent CAUGHT"))

    # M2: non-integrable declaration refused before any sampling
    caught2 = False
    try:
        nestor.integrate(lambda x: 1 / x, 0, 1, dps=D,
                         spec={"a": {"type": "algebraic", "exponent": -1}})
    except SpecRefusal:
        caught2 = True
    ok.append(_check("M2", caught2, "exponent -1 refused at declaration"))

    # M3: refusal receipt carries the named refusal and NO value
    rec3 = nestor.oracle(lambda x: 1 / mp.sqrt(x), 0, 1, dps=D,
                         spec={"a": {"type": "none"}},
                         raise_on_refusal=False)
    ok.append(_check("M3", rec3["refusal"] is not None
                     and rec3["refusal"]["class"] == "SpecRefusal"
                     and rec3["value"] is None,
                     "refusal receipt: %s, value=%s"
                     % (rec3["refusal"]["class"], rec3["value"])))

    # D1-D4: the dispersion member (nestor.dispersion) -- same assertions
    # as the pytest files in nestor/tests/, run in-process here so ONE
    # command covers the whole package.
    from nestor.tests import test_dispersion_moments as _TDM
    from nestor.tests import test_dispersion_surrogate as _TDS
    for name, fn, what in (
            ("D1", _TDM.test_moments, "power-log moments closed form"),
            ("D2", _TDM.test_addback, "closed-form add-back vs quad"),
            ("D3", _TDS.test_disp_subtracted_level5,
             "subtracted assembly level 5 vs reference"),
            ("D4", _TDS.test_addback_sign_mutation_control,
             "add-back sign mutation control")):
        saved = mp.mp.dps
        try:
            fn()
            good, detail = True, what
        except AssertionError as e:
            good, detail = False, "%s: %s" % (what, e)
        finally:
            mp.mp.dps = saved
        ok.append(_check(name, good, detail))

    wall = time.time() - t00
    passed = all(ok)
    print("[selftest] %d/%d controls passed, wall %.1f s -> %s"
          % (sum(ok), len(ok), wall, "PASS" if passed else "FAIL"))
    if receipt_path:
        from nestor import receipts as R
        rr = R.new_receipt("selftest",
                           config=dict(nproc=nproc, controls=len(ok)))
        rr["value"] = "PASS" if passed else "FAIL"
        rr["wall_s"] = round(wall, 1)
        rr["rungs"] = [dict(control=n, passed=o) for n, o in
                       zip(("P1", "P2", "P3", "P4", "P5", "P6",
                            "M1", "M2", "M3",
                            "D1", "D2", "D3", "D4"), ok)]
        R.write(rr, receipt_path)
        print("[selftest] receipt -> %s" % receipt_path)
    return 0 if passed else 3


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nproc", type=int, default=2)
    ap.add_argument("--receipt", default=None)
    a = ap.parse_args()
    return run(nproc=a.nproc, receipt_path=a.receipt)


if __name__ == "__main__":
    sys.exit(main())

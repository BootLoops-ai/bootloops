#!/usr/bin/env python3
"""Truth-checked front-door case as a standalone gate over an arbitrary
frontdoor.py file ({FILE} target for the mutation harness). rc=0 iff the
front-door spec's Muller reference behaviors hold against an INDEPENDENT
raw-flint truth (never computed through the target): the 50-digit run's
x_30 enclosure + x_100 [mid-stream +/- inf] deviation, fail-closed render
(certified digits only, 1-ulp claim audit), the adaptive solve certifying
x_100 = 6.0000000160995649, and the typed cap refusal naming achieved
digits. Used by battery leg L18 (pristine PASS + mutants caught):
  mut0 render certification check disabled -> prints uncertified digits
  mut1 solve success check forced true -> returns without the target
  mut2 deviation continuation dropped -> nan instead of the mid stream
  mut3 escalation dead -> the driver never raises precision
"""
import importlib.util
import math
import random
import sys

from flint import arb, ctx


def muller(xp, x):
    return 111 - 1130 / x + 3000 / (x * xp)


def truth(n, dps):
    """Independent oracle: raw flint at dps, never through the target."""
    old = ctx.prec
    try:
        ctx.prec = int(dps * 3.3219281) + 16
        b = [arb(2), arb(-4)]
        for _ in range(2, n + 1):
            b.append(muller(b[-2], b[-1]))
        return b
    finally:
        ctx.prec = old


def main(path):
    spec = importlib.util.spec_from_file_location("frontdoor_target", path)
    fd = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(fd)
    t = truth(100, 400)
    t30, t100 = t[30], t[100]

    # gate 1: the 50-digit run — x_30 certified enclosure, x_100 deviation
    r = fd.run(muller, x0=(2, -4), n=100, dps=50)
    x30, x100 = r.values[30], r.values[100]
    if fd.render(x30) != "6.0056":
        print(f"FAIL x30 render {fd.render(x30)!r} != '6.0056'"); return 1
    if x30.mid().str(14, radius=False, more=True) != "6.0056486887714":
        print("FAIL x30 mid != reference 6.0056486887714"); return 1
    rad30 = float(x30.rad())
    if not (1e-6 < rad30 < 1e-3):
        print(f"FAIL x30 rad {rad30} outside the ~1e-5 band"); return 1
    if not x30.contains(t30):
        print("FAIL x30 ball does not contain the independent truth"); return 1
    if float(x100.rad()) != math.inf:
        print(f"FAIL x100 rad {float(x100.rad())} != +inf"); return 1
    m100 = float(x100.mid())
    if not abs(m100 - 100.0) < 1e-6:
        print(f"FAIL x100 mid {m100} is not the fixed-precision stream "
              f"(the deviation: central value continues to the float-trap "
              f"100)"); return 1
    s = fd.render(x100)
    if not s.startswith("UNCERTIFIED"):
        print(f"FAIL x100 render not marked: {s!r}"); return 1
    try:
        fd.render(x100, strict=True)
        print("FAIL strict render of the inf ball did not refuse"); return 1
    except fd.RenderRefused:
        pass

    # gate 2: fail-closed printing on a fat ball — certified digits ONLY
    b = arb("1.234567890123") + arb(0, 1e-6)
    if fd.render(b) != "1.23457":
        print(f"FAIL fat-ball render {fd.render(b)!r} != '1.23457'"); return 1
    if fd.render(b, digits=12) != "1.23457":
        print(f"FAIL clamp: 12-digit request printed "
              f"{fd.render(b, digits=12)!r}, beyond the 6 certified"); return 1
    try:
        fd.render(b, digits=12, strict=True)
        print("FAIL strict over-request did not refuse"); return 1
    except fd.RenderRefused as e:
        if e.certified != 6:
            print(f"FAIL refusal names {e.certified} certified, not 6"); return 1
    # 1-ulp claim audit on random balls: every printed string must be
    # within one unit of its last printed place of EVERY ball point
    random.seed(20260821)
    for _ in range(100):
        m = random.uniform(-1e3, 1e3)
        x = arb(m) + arb(0, 10 ** random.uniform(-12, 2))
        s = fd.render(x)
        if s.startswith("UNCERTIFIED"):
            continue
        old = ctx.prec
        try:
            ctx.prec = 128
            p = arb(s)
            ulp = arb(10) ** fd._ulp_exponent(s)
            ok = bool(abs(p - x) <= ulp)
        except ValueError:
            ok = False
        finally:
            ctx.prec = old
        if not ok:
            print(f"FAIL claim audit: {s!r} not within 1 ulp of {x}"); return 1

    # gate 3: adaptive solve certifies the spec's x_100
    res = fd.solve(muller, 16, x0=(2, -4), n=100)
    if res.certified_digits < 16:
        print(f"FAIL solve returned {res.certified_digits} < target 16"); return 1
    if fd.render(res, digits=17) != "6.0000000160995649":
        print(f"FAIL solve render {fd.render(res, digits=17)!r} != "
              f"reference 6.0000000160995649"); return 1
    if not res.ball.contains(t100):
        print("FAIL solve ball does not contain the independent truth"); return 1
    if not (150 <= res.dps <= 500):
        print(f"FAIL solve landed at dps={res.dps}, outside the ~300-digit "
              f"escalation band"); return 1

    # gate 4: cap -> typed refusal naming achieved digits, never a number
    try:
        fd.solve(muller, 16, x0=(2, -4), n=100, max_dps=60)
        print("FAIL capped solve RETURNED instead of refusing"); return 1
    except fd.SolveRefused as e:
        if e.achieved_digits != 0 or e.target_digits != 16:
            print(f"FAIL refusal fields achieved={e.achieved_digits} "
                  f"target={e.target_digits}"); return 1

    print("PASS front door (Muller x30/x100 vs independent truth, deviation "
          "mid-stream, fail-closed render + 1-ulp audit, solve + typed refusal)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))

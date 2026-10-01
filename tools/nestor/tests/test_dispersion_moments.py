r"""nestor.dispersion battery: closed-form moments and the kernel-Taylor add-back (no expensive rho).
Run standalone (`python3 tools/nestor/tests/test_dispersion_moments.py`), under pytest, or as
legs D1/D2 of `python3 -m nestor.selftest`."""
import sys
import mpmath as mp
import os
# tools/ on the path -> `nestor.dispersion` importable (script, pytest or selftest)
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from nestor.dispersion import disp_sub as DS


def test_moments():
    mp.mp.dps = 40
    L = mp.mpf('0.7')
    for p in [mp.mpf('0.5'), mp.mpf(1), mp.mpf('2.5'), mp.mpf(3)]:
        for m in [0, 1, 2]:
            closed = DS.moment_powerlog(p, m, L)
            num = mp.quad(lambda u: u**p * (mp.log(u)**m if m else 1), [0, L], method='tanh-sinh')
            dd = abs(closed - num); d = 99 if dd == 0 else -int(mp.log10(dd / abs(num)))
            assert d >= 30, f"moment p={p} m={m}: only {d} d"
    print("moment_powerlog: all closed forms match adaptive quad >=30 d  OK")


def test_addback():
    mp.mp.dps = 40
    c = mp.mpf('0.7')
    K = lambda w: 1 / (w + c) ** 2          # analytic, pole at w=-0.7 (radius about w=1 is 1.7)
    A = -2 * mp.pi
    delta = mp.mpf('0.5')
    ref = mp.quad(lambda w: A * (w - 1) * mp.log(w - 1) * K(w), [1, 1 + delta], method='tanh-sinh',
                  maxdegree=12)
    Kc = DS.kernel_taylor(K, mp.mpf(1), 50, radius=mp.mpf('0.7'), dps=40)
    ab = DS.addback_endpoint([(A, mp.mpf(1), 1)], Kc, delta, side='left')
    dd = abs(ab - ref); d = 99 if dd == 0 else -int(mp.log10(dd / abs(ref)))
    assert d >= 25, f"add-back only {d} d"
    print(f"addback_endpoint (A (w-1)log(w-1) * 1/(w+c)^2): {d} d  OK")


if __name__ == '__main__':
    test_moments()
    test_addback()
    print("ALL TESTS PASS")

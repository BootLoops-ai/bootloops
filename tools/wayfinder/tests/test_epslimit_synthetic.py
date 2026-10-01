#!/usr/bin/env python3
"""Synthetic validation of wayfinder.epslimit against functions with KNOWN
exact eps->0 Laurent expansions.  Covers: clean finite part, genuine pole
(k_min<0), resonance/Jordan block (log(eps)), a circle/eps-FFT recovery, an
ill-conditioned narrow grid, and the honest truncation floor (more eps samples
-> more digits).  12 graded checks in 6 legs.

Collected by `pytest tools/wayfinder` (legs test_epslimit_*); also runnable as
a script.  Ambient precision is a per-leg workdps context, never global.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import mpmath as mp
from wayfinder.epslimit import extrapolate, _digits_agree

DPS = 220
PI = mp.pi          # mpmath lazy constants: evaluated at the precision in force
Z3_CONST = lambda: mp.zeta(3)


def _report(name, ok, got, exact, achieved=None):
    d = _digits_agree(got, exact)
    tag = "PASS" if ok else "FAIL"
    extra = f", xval~{achieved:.1f}d" if achieved is not None else ""
    print(f"[{tag}] {name}: recovered {d:.1f} digits{extra}")
    print(f"        got   = {mp.nstr(got, 30)}")
    print(f"        exact = {mp.nstr(exact, 30)}")
    assert ok, f"{name}: recovered only {d:.1f} digits"
    return d


def _make_grid(f, eps_list):
    return [(e, f(mp.mpc(e))) for e in eps_list]


# ---------------------------------------------------------------------------
# Case 1: clean analytic function, finite part = f(0).  Real geometric ladder.
#   f(eps) = exp(eps) * cos(pi*eps) / (1 + eps)   ; f(0) = 1
# ---------------------------------------------------------------------------
def test_epslimit_finite_ladder():
    with mp.workdps(DPS):
        f = lambda e: mp.e ** e * mp.cos(PI * e) / (1 + e)
        eps = [mp.mpf(1) / (1000 * 2 ** k) for k in range(8)]   # geometric ladder
        g = _make_grid(f, eps)
        r = extrapolate(g, max_log=2)
        exact = mp.mpc(1)
        _report("finite ladder f(0)=1", _digits_agree(r.finite_part, exact) > 30,
                r.finite_part, exact, r.achieved_digits)


# ---------------------------------------------------------------------------
# Case 2: genuine simple pole.  f(eps) = Z3/eps + PI**2 + (1+ i) eps + ...
#   tests auto k_min detection and pole-coefficient recovery.
# ---------------------------------------------------------------------------
def test_epslimit_pole():
    with mp.workdps(DPS):
        Z3 = Z3_CONST()
        f = lambda e: Z3 / e + PI ** 2 + (1 + 1j) * e - 0.5 * e ** 2 + mp.mpc('0.1') * e ** 3
        eps = [mp.mpf(1) / (500 * mp.mpf('1.7') ** k) for k in range(9)]
        g = _make_grid(f, eps)
        r = extrapolate(g, max_log=2)
        _report("pole c_{-1}=zeta3", _digits_agree(r.c(-1), Z3) > 25, r.c(-1), Z3)
        _report("pole c_0 = pi^2", _digits_agree(r.c(0), PI ** 2) > 25, r.c(0), PI ** 2)
        _report("pole c_1 = 1+i", _digits_agree(r.c(1), 1 + 1j) > 20, r.c(1), 1 + 1j)


# ---------------------------------------------------------------------------
# Case 3: RESONANCE / Jordan block: a log(eps) appears (collapsing indicial
#   exponents -> log).
#   f(eps) = 2 + 3 eps*log(eps) + 5 eps + 7 eps^2 log(eps)^2 + ...
#   finite part (eps->0) = 2.
# ---------------------------------------------------------------------------
def test_epslimit_resonance_log():
    with mp.workdps(DPS):
        A, B, C, D = mp.mpc(2), mp.mpc(3), mp.mpc(5), mp.mpc(7)
        f = lambda e: A + B * e * mp.log(e) + C * e + D * e ** 2 * mp.log(e) ** 2 + mp.mpc('1.5') * e ** 2
        eps = [mp.mpf(1) / (200 * mp.mpf('1.6') ** k) for k in range(14)]
        g = _make_grid(f, eps)
        r = extrapolate(g, max_log=3)
        print(f"        [detected log_order={r.log_order}, scoreboard={r.diagnostics.get('log_scoreboard')}]")
        _report("resonance finite part = 2", _digits_agree(r.finite_part, A) > 25,
                r.finite_part, A)
        _report("resonance c[1,1] = 3 (eps log eps)", _digits_agree(r.c(1, 1), B) > 20,
                r.c(1, 1), B)
        _report("resonance c[1,0] = 5 (eps)", _digits_agree(r.c(1, 0), C) > 18,
                r.c(1, 0), C)


# ---------------------------------------------------------------------------
# Case 4: eps-FFT on a circle.  f analytic with a known Laurent (pole) part.
#   f(eps) = 1/eps + 4 + 2 eps + 9 eps^2 ; sample on |eps|=r circle.
# ---------------------------------------------------------------------------
def test_epslimit_circle_fft():
    # run at dps 80 / N=16: the FFT path's cost grows steeply with precision
    # and node count, and 80 dps still clears the 30-digit bar by ~45 digits
    # (battery speed: this case was 137 s at dps 220 / N=24, ~3 s here).
    with mp.workdps(80):
        f = lambda e: 1 / e + 4 + 2 * e + 9 * e ** 2
        N = 16
        r0 = mp.mpf('0.01')
        eps = [r0 * mp.e ** (2j * PI * n / N) for n in range(N)]
        g = _make_grid(f, eps)
        r = extrapolate(g, max_log=0)
        print(f"        [geometry={r.diagnostics['geometry']}, method={r.method}]")
        _report("circle c_{-1}=1", _digits_agree(r.c(-1), mp.mpc(1)) > 30, r.c(-1), mp.mpc(1))
        _report("circle c_0=4", _digits_agree(r.c(0), mp.mpc(4)) > 30, r.c(0), mp.mpc(4))
        _report("circle c_2=9", _digits_agree(r.c(2), mp.mpc(9)) > 30, r.c(2), mp.mpc(9))


# ---------------------------------------------------------------------------
# Case 5: ill-conditioned narrow real grid (a sparse engine run: 1/1000..1/8000
#   only, 3 close points) but with a KNOWN function so we can grade.
#   EXACT degree-2 polynomial: 3 points fully determine it, so the finite part
#   is recovered to full working precision.  Isolates the conditioning of the
#   close/narrow grid from the truncation question.
# ---------------------------------------------------------------------------
def test_epslimit_narrow_grid():
    with mp.workdps(DPS):
        f = lambda e: mp.mpf('-1.5') + 12 * e - 80 * e ** 2
        eps = [mp.mpf(1) / 1000, mp.mpf(1) / 4000, mp.mpf(1) / 8000]
        g = _make_grid(f, eps)
        r = extrapolate(g, max_log=0)
        print(f"        [narrow grid: working_dps={r.diagnostics['working_dps']}, "
              f"method={r.method}, n_powers={r.diagnostics['n_powers']}]")
        _report("narrow exact deg-2 finite part = -1.5",
                _digits_agree(r.finite_part, mp.mpf('-1.5')) > 60,
                r.finite_part, mp.mpf('-1.5'))


# ---------------------------------------------------------------------------
# Case 6: HONEST truncation floor.  3 points on a function WITH higher-order
#   terms cannot beat the truncation error c_n * prod(eps_i): the cap is
#   information-theoretic (too few eps samples), not a tool defect.  Adding
#   points lifts it; with 8 points the degree-4 function is over-resolved.
# ---------------------------------------------------------------------------
def test_epslimit_truncation_floor():
    with mp.workdps(DPS):
        f = lambda e: mp.mpf('-1.5') + 12 * e - 80 * e ** 2 + 500 * e ** 3 - 9000 * e ** 4
        d = 0.0
        for npts, ratio in [(3, 4), (5, 2.5), (8, 1.8)]:
            eps = [mp.mpf(1) / 1000 / ratio ** k for k in range(npts)]
            g = _make_grid(f, eps)
            r = extrapolate(g, max_log=0)
            d = _digits_agree(r.finite_part, mp.mpf('-1.5'))
            print(f"        npts={npts}: finite-part recovered {d:.1f} digits "
                  f"(working_dps={r.diagnostics['working_dps']})")
        ok = d > 60
        print(f"[{'PASS' if ok else 'FAIL'}] truncation floor: 8 points recover {d:.1f} digits "
              f"(demonstrates more eps samples -> more digits)")
        assert ok, f"8-point over-resolved grid recovered only {d:.1f} digits"


if __name__ == '__main__':
    legs = (test_epslimit_finite_ladder, test_epslimit_pole,
            test_epslimit_resonance_log, test_epslimit_circle_fft,
            test_epslimit_narrow_grid, test_epslimit_truncation_floor)
    n_ok = 0
    for fn in legs:
        print(f"\n=== {fn.__name__} ===")
        try:
            fn()
            n_ok += 1
        except AssertionError as e:
            print(f"[FAIL] {fn.__name__}: {e}")
    print(f"\n{n_ok}/{len(legs)} epslimit synthetic legs passed")
    sys.exit(0 if n_ok == len(legs) else 1)

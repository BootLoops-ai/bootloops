#!/usr/bin/env python3
"""epslimit fixture-grid leg: run the FILE pipeline (load_amflow_grid ->
extrapolate) on the committed synthetic eps-grid fixtures, whose eps->0
limits are known in closed form, and GRADE the recovered coefficients against
those limits.

Fixtures (tests/fixtures_epslimit/, engine-output JSON schema; each file's
'comment' states the generating closed form):
  series_quad_eps_1_*.json    one sample per file (multi-file loading) of the
                              exact quadratic 7/2 - 3 eps + 11 eps^2
  series_pole_ladder_8pt.json 8-point geometric ladder of
                              zeta(3)/eps + pi^2 + eps/(1+eps), plus a planted
                              engine-failure marker sample (bare value '0')
                              that the loader must drop

Collected by `pytest tools/wayfinder` (legs test_epslimit_*); also runnable as
a script (exit 0 iff every check passes).  The ambient mpmath precision is set
with a workdps context per leg, never globally.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import mpmath as mp
from wayfinder.epslimit import extrapolate, load_amflow_grid, _digits_agree

FIX = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fixtures_epslimit')
DPS = 120


def _check(name, got, exact, bar):
    d = _digits_agree(got, exact)
    ok = d > bar
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {d:.1f} digits (bar {bar})")
    print(f"        got   = {mp.nstr(got, 30)}")
    print(f"        exact = {mp.nstr(exact, 30)}")
    assert ok, f"{name}: {d:.1f} digits <= bar {bar}"


def test_epslimit_quad_multifile_grid():
    """sparse multi-file grid: exact quadratic, closed-form coefficients."""
    with mp.workdps(DPS):
        paths = [os.path.join(FIX, f"series_quad_eps_1_{d}.json")
                 for d in (1000, 4000, 8000)]
        samples = load_amflow_grid(paths, part='re')
        assert len(samples) == 3, f"expected 3 samples from 3 files, got {len(samples)}"
        samples.sort(key=lambda ev: -abs(ev[0]))
        r = extrapolate(samples, max_log=0)
        print(f"--- sparse 3-file quadratic grid "
              f"(method={r.method}, working_dps={r.diagnostics['working_dps']}) ---")
        _check("quad c_0 = 7/2", r.finite_part, mp.mpf(7) / 2, 40)
        _check("quad c_1 = -3", r.c(1), mp.mpf(-3), 35)
        _check("quad c_2 = 11", r.c(2), mp.mpf(11), 30)


def test_epslimit_pole_ladder_with_failure_marker():
    """dense ladder: genuine pole, zeta(3) and pi^2 in closed form; the planted
    failed-eps sample (bare '0') must be dropped by the loader."""
    with mp.workdps(DPS):
        p = os.path.join(FIX, 'series_pole_ladder_8pt.json')
        samples = load_amflow_grid([p], part='re')
        assert len(samples) == 8, \
            f"failure-marker sample not dropped: got {len(samples)} samples, want 8"
        print("\n--- dense 8-point pole ladder (planted failure marker dropped) ---")
        samples.sort(key=lambda ev: -abs(ev[0]))
        r = extrapolate(samples, max_log=0)
        print(f"(k_min={r.k_min}, method={r.method})")
        assert r.k_min == -1, f"pole not detected: k_min={r.k_min}"
        _check("ladder c_-1 = zeta(3)", r.c(-1), mp.zeta(3), 25)
        _check("ladder c_0 = pi^2", r.c(0), mp.pi ** 2, 20)


if __name__ == '__main__':
    n_ok = n = 0
    for fn in (test_epslimit_quad_multifile_grid,
               test_epslimit_pole_ladder_with_failure_marker):
        n += 1
        try:
            fn()
            n_ok += 1
        except AssertionError as e:
            print(f"[FAIL] {fn.__name__}: {e}")
    print(f"\n{n_ok}/{n} epslimit fixture legs passed")
    sys.exit(0 if n_ok == n else 1)

#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision.
"""selftest_parametric.py -- validation gate for the generic parametric pathway
(parametric.py).

Gate: one-loop box (E=4 -> 3-dim, Cheng-Wu chart) through spec_factory + the
CBC-QMC engine must agree with the cross-validated bench_box reference to
>=4 digits and within 5 sigma of its own error bar.
Also exercises check_spec homogeneity/positivity gates on a deliberately
broken spec (must raise).  Exits nonzero on failure.
"""
import sys, os, json, time, tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gmpy2
from gmpy2 import mpfr
import parametric
import hiprec


def main():
    # The QMC pool must only start from a real entry point: on platforms whose
    # multiprocessing start method is spawn (macOS default), workers re-import
    # this module, and unguarded top-level pool code would re-launch itself.
    FAIL = 0

    # (1) box positive control through the generic pathway
    spec = parametric.box_spec(s=-1, t='-1/3', m2=1, M5=2)
    with tempfile.NamedTemporaryFile('w', suffix='.json', delete=False) as fh:
        json.dump({'box': spec}, fh)
        bpath = fh.name
    t0 = time.time()
    res = hiprec.integrate_qmc(parametric.spec_factory, {'json': bpath, 'family': 'box'},
                               3, N=8009, n_shifts=8, korobov_p=3, dps=25,
                               nproc=8, seed=20260806)
    os.unlink(bpath)
    REF = mpfr('0.1245570572327092697104024441545542480922')
    dev = abs(res['value'] - REF)
    sig = float(dev / res['error']) if res['error'] > 0 else float('inf')
    agree = -int(gmpy2.log10(dev)) if dev > 0 else 30
    print(f"[parametric box] I={str(res['value'])[:16]} agree {agree}d pull {sig:.2f}sig "
          f"wall {time.time()-t0:.1f}s")
    if agree < 4 or sig > 5:
        print("  FAIL: generic-pathway box control below gate (>=4d, <=5sig)"); FAIL += 1
    else:
        print("  PASS")

    # (2) negative control: inhomogeneous spec must be rejected
    bad = {'name': 'bad', 'n_den': 3,
           'polys': {'F0': [[[1, 0, 0], [1, 1]], [[1, 1, 0], [1, 1]]]},
           'numerator': {'gamma': 1, 'terms': [[[1, 1], [['F0', -2]]]]}}
    try:
        parametric.check_spec(bad)
        print("[check_spec neg] FAIL: inhomogeneous poly accepted"); FAIL += 1
    except AssertionError:
        print("[check_spec neg] PASS (inhomogeneous poly rejected)")

    sys.exit(1 if FAIL else 0)


if __name__ == '__main__':
    main()

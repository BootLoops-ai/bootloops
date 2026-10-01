#!/usr/bin/env python3
"""
control_tadpole.py — POSITIVE CONTROL: boundary.tadpole vs a direct mpmath
Gamma-expression oracle at 150 digits.

Runnable standalone:  python3 tests/control_tadpole.py    (seconds)
Exit code 0 = PASS, 1 = FAIL.

For three (m2, nu, d) triples the AMFlow-normalization tadpole

    T(nu; m^2, d) = (-1)^nu * Gamma(nu - d/2) / Gamma(nu) * (m^2)^(d/2 - nu)

is computed two ways:
  * boundary.tadpole(m2, nu, d, 150)          (the ported implementation),
  * the SAME closed form typed out directly against mpmath.gamma/mp.power
    at workdps 175 (independent evaluation path: different parsing, different
    guard/rounding chain — this catches port typos, workdps bugs and the
    import-dps footgun, not just luck).

PASS bars per triple (all must hold):
  * matched digits vs the direct oracle >= 148 (150-dps result, 2d slack for
    the final rounding),
  * two-precision stability: tadpole @150 vs tadpole @165 agree >= 148
    (a value that moves with working precision is noise, not a result),
  * global mp.dps untouched at exit.

The triples deliberately include a non-dyadic string mass ('2.37') so a
truncated-parse (import-dps) regression shows up as a hard digit loss.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from mpmath import mp, mpf, mpc  # noqa: E402

import boundary  # noqa: E402

DPS = 150
TRIPLES = [
    # (m2, nu, d)   d = 4-2*eps: eps = 0.2, 0.1, 0.75
    (1, 1, "3.6"),
    ("2.37", 2, "3.8"),
    (4, 3, "2.5"),
]


def _direct(m2, nu, d, wdps):
    """The tadpole closed form, typed out directly (independent path)."""
    with mp.workdps(wdps):
        m2v = mp.mpmathify(m2) if isinstance(m2, str) else mp.mpf(m2)
        dv = mp.mpmathify(d) if isinstance(d, str) else mp.mpf(d)
        val = (mpc((-1) ** nu) * mp.gamma(nu - dv / 2) / mp.gamma(nu)
               * mp.power(mpc(m2v), dv / 2 - nu))
        return +mpc(val)


def _agree(a, b, cap):
    a, b = mpc(a), mpc(b)
    diff = abs(a - b)
    scale = max(abs(a), abs(b))
    if diff == 0:
        return float(cap)
    if scale == 0:
        return 0.0
    return max(0.0, min(float(-mp.log10(diff / scale)), float(cap)))


def main():
    dps_before = mp.dps
    failures = []
    print("=" * 74)
    print(f"TADPOLE CONTROL: boundary.tadpole vs direct mpmath Gamma @ {DPS}d")
    print("=" * 74)
    for (m2, nu, d) in TRIPLES:
        t_lo = boundary.tadpole(m2, nu, d, DPS)
        t_hi = boundary.tadpole(m2, nu, d, DPS + 15)
        ref = _direct(m2, nu, d, DPS + 25)
        with mp.workdps(DPS + 30):
            d_ref = _agree(t_lo, ref, DPS + 10)
            d_pair = _agree(t_lo, t_hi, DPS + 10)
            val_str = mp.nstr(mpc(t_lo).real, 25)
        ok = d_ref >= DPS - 2 and d_pair >= DPS - 2
        print(f"  (m2={m2!r:8} nu={nu} d={d!r:6}): T = {val_str}...")
        print(f"      vs direct oracle: {d_ref:6.1f} matched digits "
              f"(bar >= {DPS - 2});  two-precision {DPS}/{DPS + 15}: "
              f"{d_pair:6.1f}  -> {'PASS' if ok else 'FAIL'}")
        if not ok:
            failures.append(
                f"(m2={m2!r}, nu={nu}, d={d!r}): oracle {d_ref:.1f}d, "
                f"pair {d_pair:.1f}d")
    if mp.dps != dps_before:
        failures.append(f"global mp.dps mutated: {dps_before} -> {mp.dps}")
    verdict = "PASS" if not failures else "FAIL"
    print("-" * 74)
    print(f"TADPOLE CONTROL VERDICT: {verdict}")
    for fmsg in failures:
        print("  FAIL:", fmsg)
    return 0 if not failures else 1


if __name__ == "__main__":
    sys.exit(main())

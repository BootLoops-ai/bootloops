#!/usr/bin/env python3
"""
test_frob_scalar.py — unit self-checks + FAST positive controls for
wayfinder/frob_scalar.py (full validation record =
the Ising-integral constant C of arXiv:1110.1705 eq.(100) at 260d
two-precision + 2 oracle controls).

These are fast ports of the two control CLASSES — NOT the
full 260d production run (minutes-long).
Both recurrences here are HAND-DERIVED (derivations inline) and verified
EXACTLY on every series term before transport — no fitted input, no
fabricated digits. dps <= 120 throughout (low-CPU test discipline).

Control A (half-integer indicial class, boundary inhomogeneity VANISHES):
    a_k = C(2k,k)^3/64^k = ((1/2)_k)^3/(k!)^3,
    F(1) = 3F2([1/2,1/2,1/2];[1,1];1) = pi/Gamma(3/4)^4   (closed form).
    a_{k+1}/a_k = ((k+1/2)/(k+1))^3
        =>  8(n+1)^3 a_{n+1} - (2n+1)^3 a_n = 0,   n >= 0.

Control B (log-resonant integer class, boundary inhomogeneity NONZERO —
lesson 1: transport MUST compose D^{deg q + 1} o L):
    h_0 = 0,  h_k = C(2k,k)^2/(16^k k)  (k >= 1),
    H(1) = int_0^1 ((2/pi) K(m) - 1)/m dm     (independent tanh-sinh oracle;
    (2/pi)K(m) = 2F1(1/2,1/2;1;m), mpmath ellipk takes the PARAMETER m).
    h_{k+1}/h_k = k (2k+1)^2 / (4 (k+1)^3)  holds for k >= 1 only; the
    shift n -> n+1 gives a recurrence valid for ALL n >= 0:
        4(n+2)^3 h_{n+2} - (n+1)(2n+3)^2 h_{n+1} = 0.
    Boundary check (hand): q_0 = c_2(-2) h_0 = 0, q_1 = c_2(-1) h_1
    = 4*(1)^3*(1/4) = 1  =>  q(x) = x, deg 1  =>  compose D^2 o L
    (asserted on the run report below — the lesson-1 tripwire).

Runnable directly (python3 tests/test_frob_scalar.py) or under pytest.
"""

import os
import sys
import time
from fractions import Fraction as F
from math import comb

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import mpmath as mp  # noqa: E402
from frob_scalar import (transport_value, ptaylor_shift, subst_1_minus_s,  # noqa: E402
                         theta_to_D, peval_frac)

NSER = 700            # series length: tail at x0=1/2 ~ 2^-700 << 1e-140
DPS_PAIR = (70, 100)  # two-precision discipline, both <= 120
ORACLE_DPS = 150


def _check_recurrence_exact(seq, r, coeffs):
    """Certificate: sum_j c_j(n) a_{n+j} = 0 EXACTLY for all n >= 0."""
    for n in range(len(seq) - r):
        val = F(0)
        for (j, k), c in coeffs.items():
            val += c * seq[n + j] * F(n) ** k
        assert val == 0, f"recurrence FAILS exactly at n={n}: {val}"


def _rec_json(r, s, coeffs):
    return {'r': r, 's': s,
            'coeffs': {f"{j},{k}": [str(c.numerator), str(c.denominator)]
                       for (j, k), c in coeffs.items()}}


def _agree_digits(a_str, b_str, dps_hi):
    with mp.workdps(dps_hi + 20):
        a, b = mp.mpf(a_str), mp.mpf(b_str)
        d = abs(a - b)
        if d == 0:
            return dps_hi
        return int(-mp.log10(d / max(abs(a), abs(b))))


def test_exact_helpers():
    """Hand-checkable identities of the exact polynomial layer."""
    # p(x) = 1 + 2x + 3x^2; p(1-s) = 6 - 8s + 3s^2
    p = [F(1), F(2), F(3)]
    assert subst_1_minus_s(p) == [F(6), F(-8), F(3)]
    # p(c+e) at c=2: p(2+e) = 17 + 14e + 3e^2
    assert ptaylor_shift(p, F(2)) == [F(17), F(14), F(3)]
    assert peval_frac(p, F(1, 2)) == F(11, 4)
    # theta^2 = x^2 D^2 + x D  (Stirling S2(2,1)=1, S2(2,2)=1)
    pd = theta_to_D({(0, 2): F(1)}, 2)
    assert pd[0] == [] and pd[1] == [F(0), F(1)] and pd[2] == [F(0), F(0), F(1)]
    print("  exact polynomial helpers OK")


def _series_A():
    return [F(comb(2 * k, k)) ** 3 / F(64) ** k for k in range(NSER + 1)]


def _series_B():
    return [F(0)] + [F(comb(2 * k, k)) ** 2 / (F(16) ** k * k)
                     for k in range(1, NSER + 1)]


def test_control_A_halfint_3f2():
    """Half-integer class vs closed form pi/Gamma(3/4)^4; two-precision
    stability; boundary inhomogeeity vanishes for this control (the class
    that would NOT catch lesson 1 on its own)."""
    dps_saved = mp.mp.dps
    t0 = time.time()
    seq = _series_A()
    # c_1(n) = 8(n+1)^3,  c_0(n) = -(2n+1)^3
    coeffs = {(1, 3): F(8), (1, 2): F(24), (1, 1): F(24), (1, 0): F(8),
              (0, 3): F(-8), (0, 2): F(-12), (0, 1): F(-6), (0, 0): F(-1)}
    _check_recurrence_exact(seq, 1, coeffs)
    rec = _rec_json(1, 3, coeffs)
    with mp.workdps(ORACLE_DPS):
        oracle = mp.nstr(mp.pi / mp.gamma(mp.mpf(3) / 4) ** 4, ORACLE_DPS - 10,
                         strip_zeros=False)
    vals = {}
    for dps in DPS_PAIR:
        out = transport_value(rec, seq, dps, verbose=False)
        vals[dps] = out['value']
        assert any("inhomogeneity: none" in ln for ln in out['report']), \
            "control A must have a vanishing boundary term"
        agree = _agree_digits(out['value'], oracle, ORACLE_DPS)
        assert agree >= dps - 8, f"dps={dps}: only {agree} digits vs closed form"
        with mp.workdps(dps):
            assert mp.mpf(out['div_gate']) < mp.mpf('1e-30'), \
                f"divergent-branch gate {out['div_gate']}"
            assert mp.mpf(out['resid_check'][0]) < mp.mpf(10) ** (-(dps + 10)), \
                f"check-point residual {out['resid_check']}"
    stab = _agree_digits(vals[DPS_PAIR[0]], vals[DPS_PAIR[1]], DPS_PAIR[1])
    assert stab >= DPS_PAIR[0] - 5, f"two-precision agreement only {stab} digits"
    assert mp.mp.dps == dps_saved, "global mp.dps was mutated (import-dps footgun)"
    print(f"  control A (3F2 half-integer): {stab} stable digits, "
          f"oracle-confirmed at both dps  [{time.time()-t0:.1f}s]")


def test_control_B_log_resonant_boundary():
    """Log-resonant integer class WITH nonzero boundary inhomogeneity:
    asserts the D^2 o L composition fired (lesson 1 tripwire) and
    gates the value against an independent tanh-sinh quadrature oracle."""
    dps_saved = mp.mp.dps
    t0 = time.time()
    seq = _series_B()
    # c_2(n) = 4(n+2)^3,  c_1(n) = -(n+1)(2n+3)^2
    coeffs = {(2, 3): F(4), (2, 2): F(24), (2, 1): F(48), (2, 0): F(32),
              (1, 3): F(-4), (1, 2): F(-16), (1, 1): F(-21), (1, 0): F(-9)}
    _check_recurrence_exact(seq, 2, coeffs)
    rec = _rec_json(2, 3, coeffs)

    # independent oracle, itself two-precision checked (never trust one quad)
    def _quad_oracle(dps):
        with mp.workdps(dps):
            f = lambda m: (2 / mp.pi * mp.ellipk(m) - 1) / m
            return mp.nstr(mp.quad(f, [0, 1]), dps - 10, strip_zeros=False)
    o_lo = _quad_oracle(115)
    oracle = _quad_oracle(ORACLE_DPS)
    q_agree = _agree_digits(o_lo, oracle, ORACLE_DPS)
    assert q_agree >= 100, f"quad oracle unstable: {q_agree} digits"

    vals = {}
    for dps in DPS_PAIR:
        out = transport_value(rec, seq, dps, verbose=False)
        vals[dps] = out['value']
        assert any("composing D^2 o L" in ln for ln in out['report']), \
            "lesson-1 tripwire: D^{deg q+1} o L composition did not fire"
        agree = _agree_digits(out['value'], oracle, ORACLE_DPS)
        assert agree >= dps - 8, f"dps={dps}: only {agree} digits vs quad oracle"
        with mp.workdps(dps):
            assert mp.mpf(out['div_gate']) < mp.mpf('1e-30'), \
                f"divergent-branch gate {out['div_gate']}"
    stab = _agree_digits(vals[DPS_PAIR[0]], vals[DPS_PAIR[1]], DPS_PAIR[1])
    assert stab >= DPS_PAIR[0] - 5, f"two-precision agreement only {stab} digits"
    assert mp.mp.dps == dps_saved, "global mp.dps was mutated (import-dps footgun)"
    print(f"  control B (log-resonant + D^2 o L boundary): {stab} stable digits, "
          f"quad-oracle-confirmed  [{time.time()-t0:.1f}s]")


def main():
    t0 = time.time()
    print("test_frob_scalar.py self-checks:")
    test_exact_helpers()
    test_control_A_halfint_3f2()
    test_control_B_log_resonant_boundary()
    print(f"ALL PASS ({time.time() - t0:.1f}s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())

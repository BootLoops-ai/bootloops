#!/usr/bin/env python3
"""acb_fast_attach.py — the ACB FAST-ATTACH member (wayfinder family).

WHAT IT IS
  Given any scalar-operator companion system sysm satisfying the
  wayfinder DESystem duck-type contract with EXACT rational
  coefficients exposed as
      sysm.n                    order (companion dimension)
      sysm.coeffs               list of n+1 Fraction lists, coeffs[j][k]
                                = coeff of x^k in p_j; operator
                                sum_j p_j(x) (d/dx)^j
      sysm.A(x, eps, dps)       companion matrix evaluator
      sysm.meta['singular_points']  ALL lc roots (physical + apparent)
  these two attachers enable the transport_fixed_eps fast paths:

    attach_fast_path(sysm)      duck-typed A_series (mpmath): Taylor
                                shift of the exact p_i + truncated series
                                division, per the transport module's
                                _local_a_series contract.
    attach_acb_fast(sysm, wp)   the package's flint acb step-kernel
                                table (acbfast.AcbFastTable) built
                                DIRECTLY from the companion entries
                                -p_j/p_n (shared-denominator dedup
                                inside the table), enabling
                                transport_fixed_eps(..., backend="acb").
                                Same chassis, same certificates; state
                                mid-trimmed per accepted step (same
                                honesty class as mpmath).  Fixed at
                                eps=0 (pure-Q operator); refuses
                                dps > wp.

BACKEND LAW (measured): backend="acb" for every large-order high-dps
march.  Same short march, same steps, same truncation certificate:
measured on an order-17 march at dps 430, mpmath 73 s -> acb 2.2 s (33x).
See GUIDE.md for the F1/F2 precision laws that bind callers.
"""
import types
from fractions import Fraction as F  # noqa: F401 (contract documentation)

import mpmath as mp


def attach_fast_path(sysm):
    """Duck-typed A_series fast path for the scalar companion system:
    Taylor coefficients [A_0..A_M] of the companion matrix about z0, via
    numeric Taylor shift of the exact p_i and truncated series division.
    Contract per the transport module's _local_a_series."""
    n = sysm.n
    coeffs_f = sysm.coeffs   # Fraction lists low->high

    def A_series(self, z0, eps_v, wp, M):
        ev = mp.mpc(eps_v)
        if ev != 0:
            raise ValueError("pure-Q operator: eps must be 0")
        with mp.workdps(wp):
            z = mp.mpc(z0)
            # numeric coefficients once per call (cached per wp on self)
            key = ("pnum", wp)
            cache = getattr(self, "_fp_cache", None)
            if cache is None or cache.get("key") != key:
                pnum = []
                for row in coeffs_f:
                    pnum.append([mp.mpf(c.numerator) / mp.mpf(c.denominator)
                                 for c in row])
                self._fp_cache = {"key": key, "pnum": pnum}
            pnum = self._fp_cache["pnum"]
            # Taylor shift each p_i about z0 (Horner shift, O(deg^2))
            shifted = []
            for row in pnum:
                out = [mp.mpc(0)]
                for c in reversed(row):
                    # out = out*(z + h) + c  ==  scale + shift
                    new = [mp.mpc(0)] * (len(out) + 1)
                    for i, v in enumerate(out):
                        if v != 0:
                            new[i] += v * z
                            new[i + 1] += v
                    new[0] += c
                    # trim trailing zeros cheaply
                    while len(new) > 1 and new[-1] == 0:
                        new.pop()
                    out = new
                shifted.append(out)
            L = M + 1
            pn = shifted[n] + [mp.mpc(0)] * max(0, L - len(shifted[n]))
            pn = pn[:L]
            if pn[0] == 0:
                raise ZeroDivisionError("A_series: lc vanishes at z0")
            # inv = 1/pn as series to L terms
            inv = [mp.mpc(0)] * L
            inv[0] = 1 / pn[0]
            for k in range(1, L):
                acc = mp.mpc(0)
                for j in range(1, min(k, len(pn) - 1) + 1):
                    if pn[j] != 0:
                        acc += pn[j] * inv[k - j]
                inv[k] = -acc * inv[0]
            # rows: q_j = p_j * inv  (to L terms)
            qser = []
            for j in range(n):
                pj = shifted[j]
                out = [mp.mpc(0)] * L
                for a, va in enumerate(pj):
                    if va != 0 and a < L:
                        for b in range(L - a):
                            vb = inv[b]
                            if vb != 0:
                                out[a + b] += va * vb
                qser.append(out)
            # assemble A_k matrices
            mats = []
            zero = mp.mpc(0)
            one = mp.mpc(1)
            for k in range(L):
                rows = [[zero] * n for _ in range(n)]
                if k == 0:
                    for i in range(n - 1):
                        rows[i][i + 1] = one
                for j in range(n):
                    v = qser[j][k]
                    if v != 0:
                        rows[n - 1][j] = -v
                mats.append(rows)
            return mats

    sysm.A_series = types.MethodType(A_series, sysm)
    return sysm


def attach_acb_fast(sysm, wp):
    """Attach the package's flint acb step-kernel table (acbfast.AcbFastTable)
    to the scalar companion system, enabling
    transport_fixed_eps(..., backend="acb").  Entries: A[16][j] = -p_j/p_n
    (shared denominator dedup inside the table), superdiagonal ones.
    Semantics per README.md: same chassis, same certificates;
    state mid-trimmed per accepted step (same honesty class as mpmath)."""
    try:
        from . import acbfast
        from .ratfun import RF
    except ImportError:  # flat mode (tests import modules directly)
        import acbfast
        from ratfun import RF
    if not acbfast.HAVE_FLINT:
        raise ImportError("python-flint required for the acb backend")
    n = sysm.n
    with mp.workdps(wp):
        one = [mp.mpc(1)]
        ents = {}
        for i in range(n - 1):
            ents[(i, i + 1)] = RF(list(one), list(one))
        den = [mp.mpc(mp.mpf(c.numerator)) / mp.mpf(c.denominator)
               for c in sysm.coeffs[n]]
        for j in range(n):
            num = [-(mp.mpc(mp.mpf(c.numerator)) / mp.mpf(c.denominator))
                   for c in sysm.coeffs[j]]
            if any(v != 0 for v in num):
                ents[(n - 1, j)] = RF(num, list(den))
        sysm._acb_fast = acbfast.AcbFastTable(ents, wp)

        def _check(eps0, dps):
            e = mp.mpc(eps0)
            if e != 0:
                raise ValueError("pure-Q operator: acb backend fixed at eps=0")
            if dps > wp:
                raise ValueError(f"acb table wp={wp} < requested dps={dps}")
        sysm._acb_check_eps = _check
    return sysm

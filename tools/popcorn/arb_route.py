#!/usr/bin/env python3
"""Production (biobank-n) route: stable-direction contiguity recursion in ARB
ball arithmetic. Certification is BY CONSTRUCTION: ball radii track the
directional contamination exactly (no pad guessing); if the final worst
relative radius misses target, double precision and rerun.

Direction (measured): S>0 forward from
(M_0, M_1); S<0 backward from (M_n, M_{n-1}). Seeds are elementary or a
single cheap positive-term 1F1(1;n;+z) / 1F1(2;n+1;+z) series (certified
geometric tail folded into the ball).

Gradient vector: d/dS of the recurrence adds source +/-M_i with the SAME
homogeneous part => same stable direction, marched alongside.

  M_{i+1} = [(2i-n+S) M_i + (n-i) M_{i-1}] / i          (forward)
  M'_{i+1} = [(2i-n+S) M'_i + M_i + (n-i) M'_{i-1}] / i
  M_{i-1} = [i M_{i+1} + (n-2i-S) M_i] / (n-i)          (backward)
  M'_{i-1} = [i M'_{i+1} + (n-2i-S) M'_i - M_i] / (n-i)

E_i = n/(i(n-i)) (1 - M_i)/(1 - e^{-S}), dE_i/dS from (M_i, M'_i) by the
quotient rule — all in arb.
"""
from flint import arb, ctx


def _hyp1f1_pos(a, b, z_arb, prec):
    """1F1(a;b;z) for z>0 via positive-term series in arb; certified tail
    (once ratio<1/2, tail < 2*term) added to the ball. O(|z|) terms — use only
    when no elementary form applies."""
    old = ctx.prec
    ctx.prec = prec
    try:
        t = arb(1)
        s = arb(1)
        k = 0
        tol = arb(2) ** (-(prec + 10))
        while True:
            t = t * (a + k) * z_arb / ((b + k) * (k + 1))
            s = s + t
            k += 1
            ratio_num = (a + k) * z_arb
            ratio_den = arb(b + k) * (k + 1)
            if ratio_num < ratio_den / 2 and t < s * tol:
                tail = t * 2
                return s + arb(0, tail.upper() if hasattr(tail, 'upper') else float(tail))
            if k > 10 ** 7:
                raise RuntimeError("series stall")
    finally:
        ctx.prec = old


def _hyp1f1_1n_elem(n, z_arb, prec):
    """EXACT elementary forms (O(n) at ANY |z|; ball tracks the small-z
    cancellation, precision-doubling loop upstream handles it):
      F1 = 1F1(1;n;z)   = (n-1)! z^{1-n} (e^z - T_{n-2}(z)),
      F2 = 1F1(2;n+1;z) = n * dF1/dz
         = n (n-1)! [ (1-n) z^{-n} (e^z - T_{n-2}) + z^{1-n} (e^z - T_{n-3}) ],
    T_k(z) = sum_{m=0}^k z^m/m!.  Verified: 1F1(1;2;z)=(e^z-1)/z.
    Returns (F1, F2)."""
    old = ctx.prec
    ctx.prec = prec
    try:
        ez = z_arb.exp()
        # T_{n-2} and T_{n-3} in one pass
        term = arb(1)
        T2 = arb(1)          # T_{n-2}
        T3 = arb(1) if n >= 3 else arb(0)   # T_{n-3}
        for m in range(1, n - 1):
            term = term * z_arb / m
            T2 += term
            if m <= n - 3:
                T3 += term
        fac = arb(1)
        for m in range(2, n):
            fac *= m         # (n-1)!
        zi = z_arb ** (1 - n)
        F1 = fac * zi * (ez - T2)
        dF1 = fac * ((1 - n) * zi / z_arb * (ez - T2) + zi * (ez - T3))
        F2 = n * dF1
        return F1, F2
    finally:
        ctx.prec = old


def _arb_S(S, S_frac):
    # int() casts: mpmath's gmpy2 backend leaks mpz into Fraction components
    return (arb(int(S_frac.numerator)) / int(S_frac.denominator)) if S_frac is not None else arb(S)


def _seeds(n, S, prec, want_grad, S_frac=None):
    """Return (i0, M pair, Mp pair) at the seed end for the stable direction."""
    old = ctx.prec
    ctx.prec = prec
    try:
        # Seed function 1F1(1;n;z), z=|S|>0: TWO regimes (measured):
        #  z >= 2n: elementary form, no cancellation (e^z dominates T_{n-2});
        #  z <  2n: positive series, geometric ratio z/(n+k) — O(n+prec) terms.
        # Elementary at z<<n cancels by the Poisson-tail factor (measured
        # -890d at n=1000, z=37.5) — the threshold keeps both branches cheap.
        def _F12(nn, z):
            if float(z) >= 2 * nn:
                return _hyp1f1_1n_elem(nn, z, prec)
            return (_hyp1f1_pos(1, nn, z, prec), _hyp1f1_pos(2, nn + 1, z, prec))
        Sa = _arb_S(S, S_frac)
        if S > 0:
            # forward: M_0 = e^{-S}; M_1 = e^{-S} * 1F1(1;n;S)
            eS = (-Sa).exp()
            F1, F2 = _F12(n, Sa)
            M0, M1 = eS, eS * F1
            if not want_grad:
                return (M0, M1), None
            # M'_0 = -e^{-S}; M'_1 = e^{-S} [ (1/n) 1F1(2;n+1;S) - 1F1(1;n;S) ]
            Mp0 = -eS
            Mp1 = eS * (F2 / n - F1)
            return (M0, M1), (Mp0, Mp1)
        else:
            # backward: M_n = 1; M_{n-1} = 1F1(1;n;-S), -S>0
            z = -Sa
            F1, F2 = _F12(n, z)
            Mn, Mn1 = arb(1), F1
            if not want_grad:
                return (Mn, Mn1), None
            # M_{n-1} = 1F1(1;n;z), z=-S => d/dS = -(1/n) 1F1(2;n+1;z)
            Mpn = arb(0)
            Mpn1 = -F2 / n
            return (Mn, Mn1), (Mpn, Mpn1)
    finally:
        ctx.prec = old


def solve_vector_arb(n, S, target_dps, want_grad=True, max_prec=1 << 22, S_frac=None):
    """Whole M vector (and dM/dS) by stable-direction arb recursion with
    precision doubling until worst relative radius <= 10^-target_dps.
    Returns dict: M (list of arb, len n+1), Mp (or None), achieved_dps, prec."""
    assert S != 0, "S=0 is the neutral limit — handle upstream (E_i = 1/i)"
    prec = int(target_dps * 3.33) + 96
    while True:
        old = ctx.prec
        ctx.prec = prec
        try:
            Sa = _arb_S(S, S_frac)
            M = [arb(0)] * (n + 1)
            Mp = [arb(0)] * (n + 1) if want_grad else None
            (A0, A1), G = _seeds(n, S, prec, want_grad, S_frac)
            if S > 0:
                M[0], M[1] = A0, A1
                if want_grad:
                    Mp[0], Mp[1] = G
                for i in range(1, n):
                    M[i + 1] = ((2 * i - n + Sa) * M[i] + (n - i) * M[i - 1]) / i
                    if want_grad:
                        Mp[i + 1] = ((2 * i - n + Sa) * Mp[i] + M[i]
                                     + (n - i) * Mp[i - 1]) / i
            else:
                M[n], M[n - 1] = A0, A1
                if want_grad:
                    Mp[n], Mp[n - 1] = G
                for i in range(n - 1, 0, -1):
                    M[i - 1] = (i * M[i + 1] + (n - 2 * i - Sa) * M[i]) / (n - i)
                    if want_grad:
                        Mp[i - 1] = (i * Mp[i + 1] + (n - 2 * i - Sa) * Mp[i]
                                     - M[i]) / (n - i)
            # worst relative radius over interior entries
            worst = 0.0
            for i in range(1, n):
                m = M[i]
                if m.mid() != 0:
                    rr = float(m.rad() / abs(m.mid()))
                    worst = max(worst, rr)
            import math
            achieved = -math.log10(worst) if worst > 0 else float(prec) * 0.301
            if achieved >= target_dps:
                return {"M": M, "Mp": Mp, "achieved_dps": achieved, "prec": prec}
        finally:
            ctx.prec = old
        prec *= 2
        if prec > max_prec:
            raise RuntimeError(f"precision cap {max_prec} hit (n={n}, S={S})")


def sfs_arb(n, S, sol, i, S_frac=None):
    """E_i and dE_i/dS as arb from a solve_vector_arb result (at its prec)."""
    old = ctx.prec
    ctx.prec = sol["prec"]
    try:
        Sa = _arb_S(S, S_frac)
        den = 1 - (-Sa).exp()
        pref = arb(n) / (arb(i) * (n - i))
        E = pref * (1 - sol["M"][i]) / den
        dE = None
        if sol["Mp"] is not None:
            # d/dS [(1-M)/(1-e^{-S})] = [-M' (1-e^{-S}) - (1-M) e^{-S}] / (1-e^{-S})^2
            eS = (-Sa).exp()
            dE = pref * (-sol["Mp"][i] * den - (1 - sol["M"][i]) * eS) / (den * den)
        return E, dE
    finally:
        ctx.prec = old

#!/usr/bin/env python3
"""Dominance SFS via the q-DEFORMATION series over certified Kummer/Beta
objects, STABLE arrangement v2.

  E_k(S,h) = C(n,k) [ A_k - B_k / Z ],   p = 2hS, qt = (1-2h)S.

  A_k = K(k, n-k; p, qt) = sum_j qt^j/j! B(k+2j, n-k) 1F1(k+2j; n+2j; p)
        (1F1 arranged positive-term via Kummer flip when p<0).

  L(x) = INT_0^x e^{-(pq+qt q^2)} dq
       = sum_l (-qt)^l/l! G_{2l}(x),
    G_m(x) = x^{m+1}/(m+1) * 1F1(m+1; m+2; -px)      [stable: NO (2l)!/p^m]
    => e^{px} G_{2l}(x) = x^{2l+1}/(2l+1) * 1F1(1; 2l+2; px)   [Kummer]
    => B_k = INT x^{k-1}(1-x)^{c-1} e^{px+qt x^2} L(x) dx
           = sum_l (-qt)^l/(l!(2l+1)) sum_r p^r/(2l+2)_r sum_j qt^j/j!
                 * B(k+2l+1+r+2j, c)
  Z = L(1).

All coefficient magnitudes are l!/j!-damped and (2l+2)_r-damped — no
factorial NUMERATORS anywhere. Exponential-scale cancellation (alternating
sums at large |p|,|qt|) is bought with a COMPUTED precision pad
(~0.4343*(|p|+|qt|)) and certified by a two-dps agreement gate upstream.
v1 (elementary G_m) diverged at |p| < ~2|qt| — measured.
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mpmath import mp, mpf, exp, fabs, log10, binomial, beta as mpbeta

BASE_PAD = 25


def _hyp1f1_pos_mpf(a, b, z, dps):
    """1F1(a;b;z), z>=0, positive-term series, certified tail."""
    t = mpf(1); s = mpf(1); k = 0
    tol = mpf(10) ** (-(dps + 8))
    while True:
        t *= (a + k) * z / ((b + k) * (k + 1))
        s += t
        k += 1
        if (a + k) * z < (b + k) * (k + 1) / 2 and t < s * tol:
            return s
        if k > 10 ** 7:
            raise RuntimeError("1F1 stall")


def E_dominant_qseries(n, k, S, h, dps=40):
    """Route-A evaluator, stable arrangement. Returns value (mpf).
    Certify by calling twice (dps, dps+15) and comparing — the smoke/selftest
    and the validation harness do this."""
    # AMBIENT-DPS FOOTGUN (measured: uniform 16-digit ceiling vs oracle across
    # the whole validation grid): converting S,h OUTSIDE the workdps block
    # rounds them at ambient precision (15!). Convert INSIDE.
    p_f = 2 * float(h) * float(S)
    qt_f = (1 - 2 * float(h)) * float(S)
    pad = BASE_PAD + int(0.4343 * (abs(p_f) + 2 * abs(qt_f))) + n // 4
    with mp.workdps(dps + pad):
        S = mpf(str(S)); h = mpf(str(h))
        p = 2 * h * S
        qt = (1 - 2 * h) * S
        c = n - k
        if qt == 0:
            raise ValueError("h=1/2: use the exact engine")
        tol = mpf(10) ** (-(dps + pad - 10))
        # ---------- A_k
        Ak = mpf(0)
        tj = mpf(1)                      # qt^j / j!
        j = 0
        while True:
            A = k + 2 * j
            if p >= 0:
                F = _hyp1f1_pos_mpf(A, A + c, p, dps + pad)
            else:
                F = exp(p) * _hyp1f1_pos_mpf(c, A + c, -p, dps + pad)
            term = tj * mpbeta(A, c) * F
            Ak += term
            j += 1
            tj *= qt / j
            if fabs(tj) * mpbeta(k, c) * max(mpf(1), exp(p)) < tol and j > 4:
                break
            if j > 100000:
                raise RuntimeError("A_k stall")
        # ---------- Z = L(1) = sum_l (-qt)^l/(l!(2l+1)) 1F1(2l+1;2l+2;-p)
        Z = mpf(0)
        tl = mpf(1)                      # (-qt)^l / l!
        l = 0
        while True:
            if -p >= 0:
                F = _hyp1f1_pos_mpf(2 * l + 1, 2 * l + 2, -p, dps + pad)
            else:
                F = exp(-p) * _hyp1f1_pos_mpf(1, 2 * l + 2, p, dps + pad)
            Z += tl / (2 * l + 1) * F
            l += 1
            tl *= -qt / l
            if fabs(tl) * max(mpf(1), exp(-p)) < tol and l > 4:
                break
            if l > 100000:
                raise RuntimeError("Z stall")
        Lmax_used = l
        # ---------- B_k triple sum with incremental Beta
        # order loops l (outer), r, j; Beta(k+2l+1+r+2j, c) updated by ratio
        Bk = mpf(0)
        tl = mpf(1)
        # per-l bound: sup_x |1F1(1;2l+2;px)| <= e^{max(p,0)}; j-sum <= B*e^{|qt|}
        exp_pq_bound = exp(max(p, mpf(0))) * exp(fabs(qt))
        for l in range(Lmax_used + dps // 2 + 8):
            cl = tl / (2 * l + 1)
            tl *= -qt / (l + 1)
            base_A = k + 2 * l + 1
            # bound for early exit on l: |cl| * B(base_A,c) * e^{|p|+|qt|}
            Bb = mpbeta(base_A, c)
            if fabs(cl) * Bb * exp_pq_bound < tol * max(fabs(Bk), mpf(1)) and l > 4:
                break
            pr = mpf(1)                  # p^r / (2l+2)_r
            r = 0
            while True:
                # inner j-sum: sum_j qt^j/j! B(base_A + r + 2j, c)
                BA = mpbeta(base_A + r, c)
                tj = mpf(1)
                inner = mpf(0)
                jj = 0
                Bcur = BA
                while True:
                    inner += tj * Bcur
                    jj += 1
                    tj *= qt / jj
                    # B(A+2,c)/B(A,c) = (A)(A+1)/((A+c)(A+c+1))
                    Aax = base_A + r + 2 * (jj - 1)
                    Bcur = Bcur * (Aax) * (Aax + 1) / ((Aax + c) * (Aax + c + 1))
                    if fabs(tj) * Bcur < tol and jj > 3:
                        break
                Bk += cl * pr * inner
                r += 1
                pr *= p / (2 * l + 1 + r)
                if fabs(pr) * BA < tol and r > 4:
                    break
                if r > 200000:
                    raise RuntimeError("r stall")
        E = binomial(n, k) * (Ak - Bk / Z)
        # round to output precision
        with mp.workdps(dps + 10):
            return +E


if __name__ == "__main__":
    # Demo/self-check on a bounded point set (extreme-h points at large |S|
    # cost the full computed pad and can run tens of minutes — deliberate).
    import time
    for (n, k, S, h) in [(20, 3, -100.0, 0.4), (20, 10, -50.0, 0.3),
                         (20, 5, 30.0, 0.7), (20, 3, -5.0, 0.05),
                         (20, 7, -20.0, 0.6), (20, 15, 60.0, 0.35)]:
        t0 = time.time()
        E1 = E_dominant_qseries(n, k, S, h, dps=30)
        E2 = E_dominant_qseries(n, k, S, h, dps=45)
        with mp.workdps(55):
            d = float(-log10(fabs((E1 - E2) / E2))) if E1 != E2 else 999
        print(f"n={n} k={k:2d} S={S:+7.1f} h={h}: {mp.nstr(E2, 12)}  "
              f"dps30-vs-45 {d:5.1f}d  {time.time()-t0:.1f}s")

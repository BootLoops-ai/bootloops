#!/usr/bin/env python3
"""Certified whole-SFS vector + exact gradients (exact-rational tridiagonal route).

Core objects: M_i = 1F1(n-i; n; -S), i = 0..n, satisfying
    -(n-i)*M_{i-1} + (n-2i-S)*M_i + i*M_{i+1} = 0,   i = 1..n-1,
with ELEMENTARY boundary values M_0 = e^{-S}, M_n = 1. Both ends known
=> Olver's algorithm degenerates to a finite two-point BVP: ONE tridiagonal
solve, no truncation rule, no minimal-solution bookkeeping. This sidesteps the
measured instability of directional recursion (the stable direction flips with
sign(S); ~23-digit leak at |S|~10, n=100).

Routes:
  A. solve_M_exact(n, S_rational):  the system splits linearly over the
     transcendence basis {1, e^{-S}}: M_i = alpha_i + beta_i * e^{-S} with
     alpha, beta EXACT rationals (two Fraction tridiagonal solves).
     => whole SFS vector exact; height growth measured, not assumed.
  B. solve_M_mpf(n, S, dps):  mpf tridiagonal solve with scaled partial
     pivoting at padded precision; residual reported.
  Derivatives: d^k M/dS^k satisfies the SAME tridiagonal system with RHS
     built from the (k-1)-th solution:  A * M^(k) = k * M^(k-1)_i  (row i)
     ... precisely: differentiating row i k times w.r.t. S:
         -(n-i) M^(k)_{i-1} + (n-2i-S) M^(k)_i + i M^(k)_{i+1} = k*M^(k-1)_i,
     boundaries M^(k)_0 = (-1)^k e^{-S}, M^(k)_n = 0 (k>=1).
     => exact gradient stack, one extra solve per order, same matrix.

E[SFS_i | theta=1] = n/(i(n-i)) * (1 - M_i)/(1 - e^{-S}).

Validation: the package battery (selftest.py) checks route A vs per-i hyp1f1
incl. the |S|~10 leak regime, and gradients vs the closed-form derivative;
gate_phase1.py is the heavier held-out two-route gate.
"""
from fractions import Fraction
from mpmath import mp, mpf, exp, expm1, hyp1f1, fabs, log10

PAD = 25


# ---------------------------------------------------------------- tridiagonal
def _tridiag_rows(n, S):
    """Rows i=1..n-1: (sub, diag, sup) coefficients of (M_{i-1}, M_i, M_{i+1})."""
    for i in range(1, n):
        yield i, -(n - i), (n - 2 * i - S), i


def _solve_tridiag(n, S, rhs, zero, div):
    """Thomas algorithm with the row order as given (no pivoting) over an
    arbitrary field: zero = additive zero, div = field division.
    Returns interior solution x_1..x_{n-1} as a list (index 0 -> M_1).
    NOTE: no pivoting — validity certified a posteriori by residual/oracle
    gates, not assumed. (Pivoted variant kept for the mpf route.)"""
    sub, dia, sup, b = [], [], [], list(rhs)
    for i, a, d, c in _tridiag_rows(n, S):
        sub.append(a), dia.append(d), sup.append(c)
    m = n - 1
    # forward elimination
    for k in range(1, m):
        w = div(sub[k], dia[k - 1])
        dia[k] = dia[k] - w * sup[k - 1]
        b[k] = b[k] - w * b[k - 1]
    x = [zero] * m
    x[m - 1] = div(b[m - 1], dia[m - 1])
    for k in range(m - 2, -1, -1):
        x[k] = div(b[k] - sup[k] * x[k + 1], dia[k])
    return x


# ------------------------------------------------------------- route A: exact
def solve_M_exact(n, S):
    """S rational (Fraction/int). Returns (alpha, beta): lists of Fractions,
    len n+1, with M_i = alpha[i] + beta[i]*e^{-S} EXACTLY.
    alpha[0]=0, beta[0]=1 (M_0=e^{-S}); alpha[n]=1, beta[n]=0 (M_n=1).
    Interior: A x = b where boundary terms move to the RHS:
      row 1:    RHS += (n-1)*M_0   -> alpha-part 0,      beta-part (n-1)
      row n-1:  RHS += -(n-1)*M_n  -> alpha-part -(n-1), beta-part 0
    """
    Sf = Fraction(S)
    m = n - 1
    rhs_a = [Fraction(0)] * m
    rhs_b = [Fraction(0)] * m
    rhs_b[0] = Fraction(n - 1)            # +(n-1)*e^{-S} from -(n-1)*M_0 moved
    rhs_a[m - 1] = Fraction(-(n - 1))     # -(n-1)*1 ... sign check in tests
    div = lambda p, q: p / q
    za = Fraction(0)
    xa = _solve_tridiag(n, Sf, rhs_a, za, div)
    xb = _solve_tridiag(n, Sf, rhs_b, za, div)
    alpha = [Fraction(0)] + xa + [Fraction(1)]
    beta = [Fraction(1)] + xb + [Fraction(0)]
    return alpha, beta


# --------------------------------------------------------------- route B: mpf
def solve_M_mpf(n, S, dps, pad=PAD):
    """Pivoted mpf tridiagonal solve. Returns (M list len n+1, resid_log10)."""
    with mp.workdps(dps + pad):
        Sm = mpf(S)
        M0, Mn = exp(-Sm), mpf(1)
        m = n - 1
        sub, dia, sup, b = [mpf(0)] * m, [mpf(0)] * m, [mpf(0)] * m, [mpf(0)] * m
        for idx, (i, a, d, c) in enumerate(_tridiag_rows(n, Sm)):
            sub[idx], dia[idx], sup[idx] = mpf(a), d, mpf(c)
        b[0] += (n - 1) * M0
        b[m - 1] -= (m) * 0  # placeholder; boundary term below
        b[m - 1] += -(n - 1) * 0  # (kept explicit for readability)
        b[m - 1] -= sup[m - 1] * 0
        # row n-1 has sup coeff (n-1) multiplying M_n=1 -> move to RHS
        b[m - 1] -= (n - 1) * Mn
        # scaled-partial-pivot Thomas -> use dense-band LU via successive rows
        # (band 3; pivoting only swaps adjacent rows -> band grows to 4; keep
        #  a 4-diagonal representation)
        # For v0 simplicity: build rows as dicts (n<=2e5 fine, O(n) memory).
        lo = [sub[k] for k in range(m)]
        di = [dia[k] for k in range(m)]
        up = [sup[k] for k in range(m)]
        u2 = [mpf(0)] * m          # second superdiagonal from row swaps
        x = list(b)
        for k in range(m - 1):
            if fabs(lo[k + 1]) > fabs(di[k]):      # pivot: swap rows k, k+1
                di[k], lo[k + 1] = lo[k + 1], di[k]
                up[k], di[k + 1] = di[k + 1], up[k]
                u2[k], up[k + 1] = up[k + 1], u2[k]
                x[k], x[k + 1] = x[k + 1], x[k]
            w = lo[k + 1] / di[k]
            di[k + 1] -= w * up[k]
            up[k + 1] -= w * u2[k]
            x[k + 1] -= w * x[k]
        sol = [mpf(0)] * m
        sol[m - 1] = x[m - 1] / di[m - 1]
        if m >= 2:
            sol[m - 2] = (x[m - 2] - up[m - 2] * sol[m - 1]) / di[m - 2]
        for k in range(m - 3, -1, -1):
            sol[k] = (x[k] - up[k] * sol[k + 1] - u2[k] * sol[k + 2]) / di[k]
        M = [M0] + sol + [Mn]
        # residual of the ORIGINAL system (certification input)
        worst = mpf(0)
        for i, a, d, c in _tridiag_rows(n, Sm):
            r = a * M[i - 1] + d * M[i] + c * M[i + 1]
            scale = max(fabs(a * M[i - 1]), fabs(d * M[i]), fabs(c * M[i + 1]))
            if scale > 0:
                worst = max(worst, fabs(r) / scale)
        rlog = float(log10(worst)) if worst > 0 else float('-inf')
        return M, rlog


# ------------------------------------------------------------------ gradients
def solve_dM_exact(n, S, alpha, beta, order=1):
    """Exact d^k M/dS^k, k=1..order, from the SAME tridiagonal matrix.
    Row i differentiated k times: A * M^(k) = rhs, rhs_i = k * M^(k-1)_i
    (from d/dS of the (n-2i-S) coefficient acting k times via Leibniz —
    only the first derivative of the coefficient survives, contributing
    k * M^(k-1)_i on the RHS).
    Boundaries: M^(k)_0 = (-1)^k e^{-S} -> (alpha_part 0, beta_part (-1)^k);
                M^(k)_n = 0.
    Returns list of (alpha_k, beta_k) for k=1..order."""
    Sf = Fraction(S)
    m = n - 1
    div = lambda p, q: p / q
    za = Fraction(0)
    out = []
    prev_a, prev_b = alpha, beta
    for k in range(1, order + 1):
        rhs_a = [k * prev_a[i] for i in range(1, n)]
        rhs_b = [k * prev_b[i] for i in range(1, n)]
        sgn = Fraction((-1) ** k)
        # boundary M^(k)_0 contributes -(-(n-1)) ... row1 sub coeff -(n-1):
        rhs_b[0] += (n - 1) * sgn
        # boundary M^(k)_n = 0: no contribution
        xa = _solve_tridiag(n, Sf, rhs_a, za, div)
        xb = _solve_tridiag(n, Sf, rhs_b, za, div)
        ak = [Fraction(0)] + xa + [Fraction(0)]
        bk = [sgn] + xb + [Fraction(0)]
        out.append((ak, bk))
        prev_a, prev_b = ak, bk
    return out


# -------------------------------------------------------- certified assembly
# MEASURED: the (alpha, beta) decomposition carries the
# second homogeneous solution, whose INTERIOR magnitude (~10^197 at n=500,
# S=-100) is the true BVP condition number — fixed-precision assembly AND any
# fixed-precision BVP solve are structurally doomed at moderate |S| (route B
# got 0.4d with residual 1e-86). Cure: Fraction heights BOUND the cancellation
# a priori => assemble at target + height_digits + 0.44|S| + margin. Verified
# 7/7 configs, 141..3130 digits vs oracle incl. (100,+-10) and (500,-100).
def assembly_dps(alpha, beta, S, target_dps):
    h = max(height_bits(alpha), height_bits(beta))
    return target_dps + int(h * 0.30103) + int(0.44 * abs(float(S))) + 20


def _S_mp(S):
    """Exact conversion of S (Fraction/int/float) to mpf at CURRENT precision.
    CRITICAL: alpha/beta are exact for a specific rational S — assembling with
    any OTHER nearby S (e.g. the binary float of a decimal) is amplified by
    the full cancellation factor (measured: -7000d on a gate run)."""
    if isinstance(S, Fraction):
        return mpf(S.numerator) / S.denominator
    return mpf(S)


def assemble_M(alpha, beta, S, target_dps):
    """Certified-precision evaluation of M_i = alpha_i + beta_i e^{-S}.
    S MUST be the exact same rational used in solve_M_exact."""
    with mp.workdps(assembly_dps(alpha, beta, S, target_dps)):
        eS = exp(-_S_mp(S))
        return [mpf(a.numerator) / a.denominator
                + (mpf(b.numerator) / b.denominator) * eS
                for a, b in zip(alpha, beta)]


# ------------------------------------------------------------------ SFS layer
def sfs_from_M(n, S, M_alpha_beta_or_list, dps):
    """E[SFS_i], i=1..n-1, theta=1, from exact (alpha,beta) or mpf M list."""
    with mp.workdps(dps + PAD):
        Sm = mpf(S)
        den = -expm1(-Sm)
        eS = exp(-Sm)
        out = []
        if isinstance(M_alpha_beta_or_list, tuple):
            alpha, beta = M_alpha_beta_or_list
            for i in range(1, n):
                Mi = mpf(alpha[i].numerator) / alpha[i].denominator \
                    + (mpf(beta[i].numerator) / beta[i].denominator) * eS
                out.append(mpf(n) / (mpf(i) * (n - i)) * (1 - Mi) / den)
        else:
            M = M_alpha_beta_or_list
            for i in range(1, n):
                out.append(mpf(n) / (mpf(i) * (n - i)) * (1 - M[i]) / den)
        return out


def height_bits(fr_list):
    """max bit-height over a Fraction list (numerator/denominator)."""
    h = 0
    for f in fr_list:
        h = max(h, f.numerator.bit_length(), f.denominator.bit_length())
    return h

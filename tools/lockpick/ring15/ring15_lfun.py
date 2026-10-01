# lockpick ring15 member — L(f,s) at s=1,2,3 for both newforms.
"""ring15_lfun: L(f,s), s=1,2,3, for the two weight-3 level-15 CM newforms.

Method (standard Hecke integral, derived for F_f(y) := f(iy/sqrt(15))):
  Lambda_f(s) := (sqrt(15)/2pi)^s Gamma(s) L(f,s) = int_0^inf F_f(y) y^{s-1} dy.
Fricke W_15 acts on the 2-dim newspace; numerically we FIT the matrix C:
  F_f(1/y) = y^3 * sum_g C[f][g] F_g(y)        (fit at 2 y's, validated at a 3rd;
  entries must come out as exact +-1/0 -- recorded, not assumed).
Split the integral at y=t and substitute y->1/y in the lower piece:
  Lambda_f(s;t) = sum_n a^f_n (cn)^{-s} G(s, cnt)
               + sum_g C[f][g] sum_n a^g_n (cn)^{s-3} G(3-s, cn/t),  c = 2pi/sqrt(15)
with G(a,x) = upper incomplete Gamma(a,x):
  G(1,x)=e^{-x}, G(2,x)=e^{-x}(1+x), G(3,x)=e^{-x}(x^2+2x+2), G(0,x)=E1(x).
t-INDEPENDENCE (t=1 vs t=1.3) is the functional-equation check: it fails unless
C is right.  Terms with cnx > (dps+15)ln10 are dropped (rigorously negligible).
"""
import mpmath as mp
from .ring15_hecke import a_n_list

C_CONST = lambda: 2 * mp.pi / mp.sqrt(15)


def _G(a, x):
    if a == 0:
        return mp.e1(x)
    if a == 1:
        return mp.exp(-x)
    if a == 2:
        return mp.exp(-x) * (1 + x)
    if a == 3:
        return mp.exp(-x) * (x * x + 2 * x + 2)
    raise ValueError(a)


def F_val(an, y):
    """F(y) = sum a_n exp(-c n y), truncated at negligible terms."""
    c = C_CONST()
    cut = (mp.mp.dps + 15) * mp.log(10)
    s = mp.mpf(0)
    for n in range(1, len(an)):
        if an[n] == 0:
            continue
        x = c * n * y
        if x > cut:
            break
        s += an[n] * mp.exp(-x)
    return s


def fit_fricke(an_p, an_m):
    """Fit F_f(1/y) y^-3 = C[f][p] F_p(y) + C[f][m] F_m(y); validate at y3."""
    y1, y2, y3 = mp.mpf('1.15'), mp.mpf('1.35'), mp.mpf('1.6')
    rows = {}
    resid = {}
    for tag, an in (('p', an_p), ('m', an_m)):
        M = mp.matrix([[F_val(an_p, y1), F_val(an_m, y1)],
                       [F_val(an_p, y2), F_val(an_m, y2)]])
        b = mp.matrix([F_val(an, 1 / y1) / y1 ** 3,
                       F_val(an, 1 / y2) / y2 ** 3])
        sol = mp.lu_solve(M, b)
        rows[tag] = (sol[0], sol[1])
        resid[tag] = F_val(an, 1 / y3) / y3 ** 3 - \
            (sol[0] * F_val(an_p, y3) + sol[1] * F_val(an_m, y3))
    return rows, resid


def Lambda(an_f, C_row, an_p, an_m, s, t):
    """Lambda_f(s) with split point t (exact for any t>0)."""
    c = C_CONST()
    cut = (mp.mp.dps + 15) * mp.log(10)
    tot = mp.mpf(0)
    for n in range(1, len(an_f)):
        if an_f[n] == 0:
            continue
        x = c * n * t
        if x > cut:
            break
        tot += an_f[n] * (c * n) ** (-s) * _G(s, x)
    for cg, an_g in zip(C_row, (an_p, an_m)):
        if cg == 0:
            continue
        for n in range(1, len(an_g)):
            if an_g[n] == 0:
                continue
            x = c * n / t
            if x > cut:
                break
            tot += cg * an_g[n] * (c * n) ** (s - 3) * _G(3 - s, x)
    return tot


def L_values(NL=1400):
    """Returns dict: L for both forms at s=1,2,3 + t-invariance + Fricke data."""
    an_p = a_n_list(NL, +1)
    an_m = a_n_list(NL, -1)
    rows, resid = fit_fricke(an_p, an_m)
    out = {'fricke_rows': rows, 'fricke_resid_y3': resid}
    c = C_CONST()
    for tag, an in (('p', an_p), ('m', an_m)):
        C_row = rows[tag]
        for s in (1, 2, 3):
            lam1 = Lambda(an, C_row, an_p, an_m, s, mp.mpf(1))
            lam2 = Lambda(an, C_row, an_p, an_m, s, mp.mpf('1.3'))
            L = lam1 / ((1 / c) ** s * mp.gamma(s))
            out[f'L_{tag}_{s}'] = L
            out[f'tinv_{tag}_{s}'] = lam1 - lam2
    return out


def direct_sum_check(Lval_s3, sign, N=2 * 10 ** 6):
    """Float partial sum of sum a_n / n^3 -- ~5-digit sanity vs Lval_s3."""
    an = a_n_list(N, sign)
    s = 0.0
    for n in range(1, N + 1):
        if an[n]:
            s += an[n] / n ** 3
    return s - float(Lval_s3)

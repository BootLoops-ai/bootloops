#!/usr/bin/env python3
"""
ladder_reference.py — closed-form reference plugin for the 4L ladder family.

Conventions pinned from arXiv:1303.6909 (the conventions of the 100-digit
reference values phi_values.json, a record that is not shipped here):

    u = z*zbar,  v = (1-z)(1-zbar)
    Phi^(L)(u,v) = -f^(L)(w, wbar) / (z - zbar),   w = z/(z-1), wbar = zbar/(zbar-1)
    f^(L)(x,xbar) = sum_{r=0}^{L} c_r log^r(x*xbar) (Li_{2L-r}(x) - Li_{2L-r}(xbar)),
    c_r = (-1)^r (2L-r)! / (r! (L-r)! L!)   (exact Fractions).

CHART: the harness chart variables are (a, b) := (w, wbar) treated as two
INDEPENDENT complex variables. Inverse relations (exact rational):
    z = a/(a-1), zbar = b/(b-1),  z - zbar = (b-a)/((a-1)(b-1))
    prefactor  -1/(z-zbar) = (a-1)(b-1)/(a-b)
    u = a*b/((a-1)(b-1)),   v = 1/((a-1)(b-1))
On the Euclidean slice b = conj(a) with u,v in the Euclidean region, Phi is
real (the reference values).

MASTER BASIS (weight-graded, closes under d):
    ell := log(a*b)  (as a DE solution; anchored principal at real 0<a0,b0<1)
    masters: {'lA',r,s} = ell^r * Li_s(a);  {'lB',r,s} = ell^r * Li_s(b);
             {'l',r}   = ell^r   (r=0 -> the constant 1)
    for r=0..L, s=1..2L, r+s <= 2L. n = 65 at L=4.

THE SHEET OF ell (the sheet-selected default). With no explicit `ell_k`,
ell is the ANCHORED PRINCIPAL branch of the master, i.e. the principal
log(a*b), realized as log(a) + log(b) + 2*pi*i*k with the integer k chosen
so that the sum lands on the principal branch of the product
(`default_ell_k`: k = round((Im log(a*b) - Im(log a + log b)) / 2 pi)).
    k = 0  on the Euclidean slice b = conj(a) and at the anchor 0 < a0, b0 < 1
           (there log(a)+log(b) already is the principal log of the product);
    k = -1 on the real lambda > 0 sheet (z, zbar real, u + v < 1), where the
           chart point a = w, b = wbar has BOTH coordinates negative real and
           the principal sum log(a)+log(b) = ln(u/v) + 2*pi*i overshoots the
           real master ell = ln(u/v) by one turn.
An EXPLICIT cont['ell_k'] keeps its meaning RELATIVE TO log(a)+log(b) (below),
so on the real lambda > 0 sheet cont={'ell_k': -1} and the default agree, and
cont={'ell_k': 0} is the naive sum with the extra 2*pi*i (wrong there).
`ell_sheet(a, b, dps, cont)` reports the k in force and whether it was
selected by the default rule or given explicitly.

CONTINUATION-SAFE reference (the bug class: NAIVE principal-branch polylog
is WRONG on non-principal sheets — mpmath plateau law). The `cont` dict
applies EXPLICIT continuation data:
    cont = {'ell_k': integer  -> ell = log(a)+log(b) + 2*pi*i*ell_k
            'liA_disc': (m, sigma) -> Li_s(a) += sigma * m * 2*pi*i * log^{s-1}(a)/(s-1)!}
liA_disc encodes m crossings of the a in (1,inf) cut; sigma=+1 for a crossing
from Im a > 0 to Im a < 0 (derivation in GATES receipt; independently pinned
in gate 5 by the ELEMENTARY s=1 case Li_1 = -log(1-a), whose continuation is
computed here by dense_winding on 1-a — no shared code with the harness
tracker). Valid for paths that stay in Re a > 0 (principal log(a) continuous),
which covers every gate path.

Everything numeric is mpmath at caller dps inside mp.workdps; coefficients
are exact Fractions.
"""
from fractions import Fraction
from math import factorial
from mpmath import mp, mpf, mpc, polylog, log, pi, binomial

TWO_PI_I = None  # built per-dps


def masters_spec(L=4):
    """Ordered master list. Index order is THE basis order everywhere."""
    ms = []
    for r in range(L + 1):
        for s in range(1, 2 * L + 1):
            if r + s <= 2 * L:
                ms.append({"type": "lA", "r": r, "s": s})
    for r in range(L + 1):
        for s in range(1, 2 * L + 1):
            if r + s <= 2 * L:
                ms.append({"type": "lB", "r": r, "s": s})
    for r in range(L + 1):
        ms.append({"type": "l", "r": r})
    return ms


def f_coeffs(L=4):
    """c_r of eq (2.3), exact Fractions, r=0..L."""
    return [Fraction((-1) ** r * factorial(2 * L - r),
                     factorial(r) * factorial(L - r) * factorial(L))
            for r in range(L + 1)]


def _li_cont(s, a, cont, dps):
    """Li_s(a) with explicit continuation data (a-side only)."""
    with mp.workdps(dps):
        v = polylog(s, a)
        if cont and cont.get("liA_disc"):
            m, sigma = cont["liA_disc"]
            if m:
                v += sigma * m * 2 * pi * mpc(0, 1) * log(a) ** (s - 1) / factorial(s - 1)
        return v


def default_ell_k(a, b, dps):
    """The sheet-selected k: the integer with log(a) + log(b) + 2*pi*i*k equal
    to the PRINCIPAL log(a*b) (the anchored principal branch of the master
    ell). 0 on the Euclidean slice b = conj(a) and at the anchor 0<a0,b0<1;
    -1 on the real lambda > 0 sheet (a, b both negative real)."""
    with mp.workdps(dps):
        a = mpc(a); b = mpc(b)
        gap = mp.im(log(a * b)) - mp.im(log(a) + log(b))
        return int(mp.nint(gap / (2 * pi)))


def _ell_cont(a, b, cont, dps):
    """ell at (a, b). An explicit cont['ell_k'] (any integer, 0 included) is
    applied RELATIVE TO log(a) + log(b); with no explicit ell_k the
    sheet-selected default `default_ell_k` puts ell on the principal branch
    of the product."""
    with mp.workdps(dps):
        e = log(a) + log(b)
        if cont is not None and cont.get("ell_k") is not None:
            e += cont["ell_k"] * 2 * pi * mpc(0, 1)
        else:
            e += default_ell_k(a, b, dps) * 2 * pi * mpc(0, 1)
        return e


def ell_sheet(a, b, dps, cont=None):
    """Report the sheet of ell at (a, b): {'ell_k': k in force, 'selected_by':
    'default' | 'explicit', 'ell': the value}. Same rule as `_ell_cont`."""
    with mp.workdps(dps):
        a = mpc(a); b = mpc(b)
        explicit = cont is not None and cont.get("ell_k") is not None
        k = int(cont["ell_k"]) if explicit else default_ell_k(a, b, dps)
        return {"ell_k": k, "selected_by": "explicit" if explicit else "default",
                "ell": _ell_cont(a, b, cont, dps)}


def basis_values(a, b, dps, masters, cont=None):
    """J at (a,b): list of mpc in master order. cont: continuation spec."""
    with mp.workdps(dps + 10):
        a = mpc(a); b = mpc(b)
        ell = _ell_cont(a, b, cont, dps + 10)
        smax = max((m["s"] for m in masters if m["type"] != "l"), default=0)
        liA = {s: _li_cont(s, a, cont, dps + 10) for s in range(1, smax + 1)}
        liB = {s: polylog(s, b) for s in range(1, smax + 1)}  # b never continued in gates
        out = []
        for m in masters:
            if m["type"] == "l":
                out.append(ell ** m["r"])
            elif m["type"] == "lA":
                out.append(ell ** m["r"] * liA[m["s"]])
            else:
                out.append(ell ** m["r"] * liB[m["s"]])
        return out


def anchor_values(point, dps, masters):
    """J at the anchor (principal branch). point: dict {'a': mpf-able, 'b': ...}."""
    return basis_values(point["a"], point["b"], dps, masters, cont=None)


def prefactor(a, b):
    return (a - 1) * (b - 1) / (a - b)


def f_L(L, a, b, dps, cont=None):
    with mp.workdps(dps + 10):
        a = mpc(a); b = mpc(b)
        ell = _ell_cont(a, b, cont, dps + 10)
        cs = f_coeffs(L)
        s_ = mpc(0)
        for r in range(L + 1):
            c = mpf(cs[r].numerator) / cs[r].denominator
            s_ += c * ell ** r * (_li_cont(2 * L - r, a, cont, dps + 10) - polylog(2 * L - r, b))
        return s_


def phi_value(L, a, b, dps, cont=None):
    """Phi^(L) at chart point (a,b), optionally continued."""
    with mp.workdps(dps + 10):
        a = mpc(a); b = mpc(b)
        return prefactor(a, b) * f_L(L, a, b, dps, cont)


def phi_derivs_chart(L, a, b, dps):
    """(dPhi/da, dPhi/db) by TERMWISE analytic differentiation (principal
    polylogs; ell on the same sheet-selected default as `phi_value` with
    cont=None, so on the real lambda > 0 sheet the derivatives belong to the
    real master). Independent of the harness transport: direct mpmath
    polylog calls at the point. NOTE d Li_s(a)/da = Li_{s-1}(a)/a (s>=2),
    = 1/(1-a) for s=1 — all gate uses have s=2L-r >= L >= 2."""
    with mp.workdps(dps + 10):
        a = mpc(a); b = mpc(b)
        ell = _ell_cont(a, b, None, dps + 10)
        cs = f_coeffs(L)
        f = mpc(0); dfa = mpc(0); dfb = mpc(0)
        for r in range(L + 1):
            c = mpf(cs[r].numerator) / cs[r].denominator
            s = 2 * L - r
            dA = polylog(s, a) - polylog(s, b)
            f += c * ell ** r * dA
            t1a = (r * ell ** (r - 1) / a * dA) if r >= 1 else mpc(0)
            t1b = (r * ell ** (r - 1) / b * dA) if r >= 1 else mpc(0)
            dfa += c * (t1a + ell ** r * (polylog(s - 1, a) / a if s >= 2 else 1 / (1 - a)))
            dfb += c * (t1b - ell ** r * (polylog(s - 1, b) / b if s >= 2 else 1 / (1 - b)))
        pref = prefactor(a, b)
        dpref_a = -(b - 1) ** 2 / (a - b) ** 2
        dpref_b = (a - 1) ** 2 / (a - b) ** 2
        return dpref_a * f + pref * dfa, dpref_b * f + pref * dfb


def kin_jacobian(a, b, dps):
    """[[du/da, dv/da],[du/db, dv/db]] for u=ab/((a-1)(b-1)), v=1/((a-1)(b-1))."""
    with mp.workdps(dps + 10):
        a = mpc(a); b = mpc(b)
        du_da = -b / ((a - 1) ** 2 * (b - 1))
        du_db = -a / ((a - 1) * (b - 1) ** 2)
        dv_da = -1 / ((a - 1) ** 2 * (b - 1))
        dv_db = -1 / ((a - 1) * (b - 1) ** 2)
        return [[du_da, dv_da], [du_db, dv_db]]


def phi_derivs_uv(L, a, b, dps):
    """(dPhi/du, dPhi/dv) closed form: chart derivs + 2x2 Jacobian solve."""
    with mp.workdps(dps + 10):
        da_, db_ = phi_derivs_chart(L, a, b, dps)
        Jk = kin_jacobian(a, b, dps)
        det = Jk[0][0] * Jk[1][1] - Jk[0][1] * Jk[1][0]
        du = (da_ * Jk[1][1] - db_ * Jk[0][1]) / det
        dv = (db_ * Jk[0][0] - da_ * Jk[1][0]) / det
        return du, dv


def monodromy_a0_matrix(masters, dps, k=1):
    """EXACT reference monodromy for the branch-locus loop 'a winds 0 by k
    (ccw)' based at |a|<1, |b|<1: Li_s regular in the unit disk (unchanged);
    ell -> ell + 2*pi*i*k. Derived from the FUNCTION representation, i.e.
    independently of the connection/transport. Dense n x n mpc matrix M with
    J_after = M . J_before.
    NOTE (receipt): exp(2*pi*i*Res_{a=0}A) is NOT this matrix — the residue
    mixes Li-weights and only path-ordering corrections from the (1-a) letter
    cancel that mixing. The functional derivation is the honest reference."""
    with mp.workdps(dps + 10):
        idx = {}
        for i, m in enumerate(masters):
            key = (m["type"], m["r"], m.get("s"))
            idx[key] = i
        n = len(masters)
        M = [[mpc(0)] * n for _ in range(n)]
        shift = 2 * pi * mpc(0, 1) * k
        for i, m in enumerate(masters):
            r = m["r"]
            for j in range(r + 1):
                coef = binomial(r, j) * shift ** (r - j)
                key = (m["type"], j, m.get("s"))
                M[i][idx[key]] += coef
        return M


def dense_winding(w_of_t, n0=64, tol=None, dps=30, max_ref=14):
    """Independent continuous-arg tracker: total d(arg w) over t in [0,1] by
    plain principal-arg increment accumulation on a refining grid (NO shared
    code with the harness bisection tracker). Refine x2 until two successive
    totals agree to tol. w_of_t: callable t->mpc."""
    with mp.workdps(dps):
        tol = tol or mpf(10) ** (-(dps - 8))
        prev = None
        n = n0
        for _ in range(max_ref):
            tot = mpf(0)
            wprev = w_of_t(mpf(0))
            for i in range(1, n + 1):
                w = w_of_t(mpf(i) / n)
                if w == 0 or wprev == 0:
                    raise ValueError("dense_winding: path through zero")
                tot += mp.arg(w / wprev)
                wprev = w
            if prev is not None and abs(tot - prev) < tol:
                return tot
            prev = tot
            n *= 2
        raise ValueError("dense_winding: no convergence (path too close to zero?)")

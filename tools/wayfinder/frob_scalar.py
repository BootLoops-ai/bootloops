#!/usr/bin/env python3
r"""
frob_scalar.py — scalar Frobenius transport-to-singular-point library:
value of a holonomic power series AT a regular-singular point of its ODE,
from an exact P-recurrence for its Taylor coefficients.

VALIDATION (the engine's own record):
  * The Ising-integral constant C = 3I1^> - 4I2^> (eq. (100) of
    arXiv:1110.1705): 260 digits two-precision
    certified (260/320 dps pair; 229 at 230/260), check-point residual
    ~7.7e-299, divergent-branch coefficients exactly 0, first 200 digits
    confirm (and are confirmed by) arXiv:1110.1705.
  * Oracle control 1: sum_k C(2k,k)^2/(16^k k) — log-resonant integer class
    with a NONZERO boundary inhomogeneity — 260d agreement vs an independent
    tanh-sinh quadrature oracle.
  * Oracle control 2: 3F2([1/2]^3;1,1;1) — half-integer class — 259d
    agreement vs the closed form pi/Gamma(3/4)^4.
Fast ports of both controls (hand-derived recurrences, dps <= 120) live in
tests/test_frob_scalar.py.

Adaptations in this copy relative to the validated original (ONLY these;
everything else verbatim):
  * the tools/ directory for the annihilator import is resolved relative to
    this file instead of a hard-coded absolute path;
  * this docstring (provenance + design warnings) and __all__.

DESIGN WARNINGS (measured failure modes — keep them)
====================================================
1. SEQUENCE-RECURRENCE FITS DO NOT ANNIHILATE THE GENERATING FUNCTION.
   A P-recurrence fitted to the coefficients holds for n >= 0 only, so the
   corresponding operator gives L F = q(x) with q a polynomial of degree
   < r (boundary inhomogeneity), NOT L F = 0. You MUST build the exact q
   and compose L' = D^{deg q + 1} o L, which annihilates F exactly.
   transport_value() does this automatically — if you lift pieces of this
   module elsewhere, keep that step. History: control 1 caught the omission
   only via the independent-check-point gate (residual 3e-4); control 2
   alone would NOT have caught it (its boundary term vanishes). Two controls
   of different classes earned their keep.
2. THE FROBENIUS EPS-TRICK NEEDS THE PER-ROOT DERIVATIVE WINDOW.
   For each indicial root rho, seed a_0 = eps^{P_rho} with P_rho = total
   multiplicity of strictly-HIGHER congruent roots, and take solutions
   y_j = (1/j!) d^j/deps^j [...] for j in [P_rho, P_rho + mult) — per ROOT,
   not per congruence class from the smallest root. Degenerate NON-LOG
   resonances (the RHS vanishing exactly at a resonant level) otherwise
   silently drop basis elements; the resonance-drop residual gate reports
   but cannot repair a mis-windowed basis.

WHICH FROBENIUS TOOL WHEN (see also boundary_branches.py)
=========================================================
  * wayfinder.frobenius   — MATRIX first-order systems landing at a FINITE
    regular-singular point (Sylvester recursion; integer indicial gaps
    removed by shearing; real-eps0 matching policy).
  * wayfinder.frob_scalar — THIS module: SCALAR operators / P-recurrences
    with full resonant integer towers (e.g. indicial roots {0..6} u
    {3/2, 3/2}), transport of the series value to the singular point.
  * tools/frobenius_boundary (thin wrapper wayfinder.boundary_branches) —
    branch CLASSIFICATION at var = INFINITY for cut-block DEs
    (own selftest; eig(D1) caveat flagged there).

Dependencies: mpmath + stdlib at import; sympy and tools/annihilator.py are
imported LAZILY inside the functions that need them (rational_roots,
transport_value) — package import stays sympy-free.

PIPELINE (all exact until the final mpf evaluation)
===================================================
Given an exact P-recurrence  sum_{j=0}^r c_j(n) a_{n+j} = 0  for the Taylor
coefficients of F(x) = sum a_k x^k (radius 1, singular point at x=1), compute
F(1) = lim_{x->1^-} F(x) to arbitrary dps:

  1. rec -> theta-form ODE  L = sum c_{ij} x^i theta^j   (exact, annihilator.rec_to_theta)
  2. theta-form -> D-form   L = sum p_i(x) D^i           (exact, Stirling)
  3. shift s = 1-x          Lt = sum pt_i(s) Ds^i,  pt_i(s) = (-1)^i p_i(1-s)
  4. delta-form at s=0      s^v Lt = sum_m s^m f_m(delta),  delta = s d/ds
  5. Frobenius basis per congruence class of indicial roots (eps-jet method,
     a_0 = eps^P, P = sum of multiplicities of resonant higher roots;
     solutions y_j = (1/j!) d^j/deps^j [ s^{rho+eps} sum a_N(eps) s^N ] at eps=0)
  6. match sum_b c_b y_b(s) against exact-series values F^{(i)}(x0) at x0 = 1-s0
  7. gate: residual at independent check point s1; divergent-branch coeffs ~ 0
  8. F(1) = sum_b c_b * value(y_b),  value = coeff of s^0 log^0 branch.

All indicial data exact (Fraction); series recurrence in mpf eps-jets.
mp.workdps discipline: every entry point takes dps and wraps.
"""
from fractions import Fraction as F
import mpmath as mp

__all__ = ["transport_value", "load_series", "theta_to_D", "delta_form_at_s0",
           "build_basis", "frobenius_root", "eval_sol_derivs", "sol_value_at_0",
           "series_derivs_mpf", "rational_roots"]


# ---------- exact polynomial helpers (Fraction coefficient lists, low->high) ----------

def ptrim(p):
    while p and p[-1] == 0:
        p.pop()
    return p

def padd(a, b):
    n = max(len(a), len(b))
    return ptrim([ (a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)
                   for i in range(n) ])

def pmul(a, b):
    if not a or not b:
        return []
    out = [F(0)] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                out[i + j] += ai * bj
    return ptrim(out)

def pscal(c, a):
    return ptrim([c * x for x in a])

def pshift_scale(p):  # multiply by s (raise powers by 1)
    return [F(0)] + list(p)

def peval_frac(p, c):
    """exact evaluation at rational c"""
    v = F(0)
    for coef in reversed(p):
        v = v * c + coef
    return v

def ptaylor_shift(p, c):
    """coefficients of p(c + e) in e, exact"""
    # Horner in (c+e): repeatedly v = v*(c+e) + coef
    out = [F(0)]
    for coef in reversed(p):
        # out = out*(c+e): (c*out) + shift(out)
        out = padd(pscal(c, out), pshift_scale(out))
        out = padd(out, [coef])
    return out if out else [F(0)]

def subst_1_minus_s(p):
    """p(1-s) as poly in s, exact"""
    out = [F(0)]
    for coef in reversed(p):
        out = padd(pscal(F(-1), pshift_scale(out)), out)  # out*(1-s)
        out = padd(out, [coef])
    return out if out else [F(0)]


# ---------- step 1-4: exact operator manipulation ----------

def theta_to_D(cij, order):
    """cij: dict {(xpow,thetapow): Fraction}. Return list p[i] = poly in x (Fraction list),
    L = sum_i p[i](x) D^i, using theta^j = sum_k S2(j,k) x^k D^k."""
    # Stirling numbers of the second kind
    S2 = [[F(1)]]
    for j in range(1, order + 1):
        prev = S2[-1]
        row = [F(0)] * (j + 1)
        for k in range(j + 1):
            row[k] = (prev[k] if k < len(prev) else F(0)) * k + \
                     (prev[k - 1] if 1 <= k <= len(prev) else F(0))
        S2.append(row)
    p = [[] for _ in range(order + 1)]
    for (i, j), c in cij.items():
        for k in range(j + 1):
            s2 = S2[j][k] if k < len(S2[j]) else F(0)
            if s2:
                mono = [F(0)] * (i + k) + [c * s2]
                p[k] = padd(p[k], mono)
    return [ptrim(x) for x in p]


def falling_fact_poly(i):
    """delta(delta-1)...(delta-i+1) as Fraction list in delta"""
    out = [F(1)]
    for m in range(i):
        out = pmul(out, [F(-m), F(1)])
    return out


def delta_form_at_s0(pt):
    """pt[i] = poly in s (exact), operator sum pt_i(s) Ds^i.
    Return list fm: fm[m] = poly in delta (exact), with s^v L = sum_m s^m f_m(delta).
    """
    order = len(pt) - 1
    vals = []
    for i, poly in enumerate(pt):
        if not poly:
            vals.append(None)
            continue
        v0 = next(q for q, c in enumerate(poly) if c != 0)
        vals.append(v0)
    v = max(i - vals[i] for i in range(order + 1) if vals[i] is not None)
    fm = []
    for i, poly in enumerate(pt):
        if not poly:
            continue
        ff = falling_fact_poly(i)
        # s^{v-i} * poly: powers m = q + v - i
        for q, c in enumerate(poly):
            if c == 0:
                continue
            m = q + v - i
            assert m >= 0, "delta-form: negative power (irregular?)"
            while len(fm) <= m:
                fm.append([])
            fm[m] = padd(fm[m], pscal(c, ff))
    return fm


def rational_roots(poly):
    """all rational roots (with multiplicity) of exact Fraction poly; assert fully factored."""
    from sympy import Poly, Rational, symbols, roots
    x = symbols('x')
    P = Poly(sum(Rational(c.numerator, c.denominator) * x**i for i, c in enumerate(poly)), x)
    rts = roots(P)
    out = []
    tot = 0
    for r, m in rts.items():
        if r.is_rational:
            out.extend([F(int(r.p), int(r.q))] * m)
            tot += m
    deg = P.degree()
    if tot != deg:
        raise ValueError(f"indicial poly has non-rational roots: {rts} (deg {deg}, rational {tot})")
    return sorted(out)


# ---------- eps-jets (mpf coefficient lists, fixed truncation J) ----------

class Jet:
    __slots__ = ('c',)
    def __init__(self, c):
        self.c = c  # list of mpf/mpc, len J

def jet_from_fracs(fr, J):
    return Jet([mp.mpf(f.numerator) / mp.mpf(f.denominator) if i < len(fr) else mp.mpf(0)
                for i, f in enumerate(list(fr) + [F(0)] * J)][:J])

def jmul(a, b, J):
    out = [mp.mpf(0)] * J
    for i, ai in enumerate(a.c):
        if ai:
            for j in range(J - i):
                bj = b.c[j]
                if bj:
                    out[i + j] += ai * bj
    return Jet(out)

def jadd(a, b, J):
    return Jet([a.c[i] + b.c[i] for i in range(J)])

def jscale(s, a, J):
    return Jet([s * x for x in a.c])

def jdiv_unit(r, u, J):
    """r / u with u.c[0] != 0"""
    inv0 = 1 / u.c[0]
    out = [mp.mpf(0)] * J
    for i in range(J):
        acc = r.c[i]
        for k in range(i):
            acc -= out[k] * u.c[i - k]
        out[i] = acc * inv0
    return Jet(out)


# ---------- step 5: Frobenius class solutions ----------

def frobenius_root(fm_exact, rho, mult, Prho, Nmax, dps_work, report):
    """Frobenius solutions attached to ONE indicial root rho (Fraction) of
    multiplicity mult, with Prho = total multiplicity of strictly-higher
    congruent roots.  Classical recipe: lambda = rho + eps, a_0 = eps^Prho,
    solutions y_j = (1/j!) d^j/deps^j [ sum_N a_N(eps) s^{lambda+N} ]|_0
    for j in [Prho, Prho+mult).
    Each solution: dict terms[(N,p)] = coeff of s^{rho+N} log(s)^p."""
    J = 2 * Prho + mult + 2  # jet truncation with division-loss guard
    M = len(fm_exact) - 1

    a = [None] * (Nmax + 1)
    a0 = [mp.mpf(0)] * J
    a0[Prho] = mp.mpf(1)
    a[0] = Jet(a0)

    def fjet(midx, c):
        return jet_from_fracs(ptaylor_shift(fm_exact[midx], c), J)

    resid_gate = mp.mpf(10) ** (-(dps_work - 25))
    max_drop = mp.mpf(0)
    for N in range(1, Nmax + 1):
        rhs = Jet([mp.mpf(0)] * J)
        any_rhs = False
        for m in range(1, min(N, M) + 1):
            if not fm_exact[m]:
                continue
            fj = fjet(m, rho + N - m)
            if any(x != 0 for x in fj.c):
                rhs = jadd(rhs, jmul(fj, a[N - m], J), J)
                any_rhs = True
        rhs = Jet([-x for x in rhs.c])
        sh0 = ptaylor_shift(fm_exact[0], rho + N)
        mu = 0
        while mu < len(sh0) and sh0[mu] == 0:
            mu += 1
        u = jet_from_fracs(sh0[mu:], J)
        if mu > 0:
            scale = max([abs(x) for x in rhs.c] + [mp.mpf('1e-300')])
            for q in range(mu):
                drop = abs(rhs.c[q]) / scale
                if drop > max_drop:
                    max_drop = drop
            rhs = Jet(rhs.c[mu:] + [mp.mpf(0)] * mu)
        a[N] = jdiv_unit(rhs, u, J) if any_rhs else Jet([mp.mpf(0)] * J)

    if max_drop > resid_gate:
        report.append(f"WARN root rho={rho}: resonance drop residual {mp.nstr(max_drop,3)}")
    report.append(f"root rho={rho}: mult={mult} Prho={Prho} J={J} "
                  f"max resonance-drop {mp.nstr(max_drop,3)}")

    facs = [mp.factorial(p) for p in range(J + 1)]
    sols = []
    for j in range(Prho, Prho + mult):
        terms = {}
        for N in range(Nmax + 1):
            for p in range(j + 1):
                q = j - p
                if q < J:
                    v = a[N].c[q] / facs[p]
                    if v != 0:
                        terms[(N, p)] = v
        if not terms:
            report.append(f"WARN root rho={rho} j={j}: identically zero solution")
        sols.append({'rho': rho, 'terms': terms, 'j': j, 'P': Prho})
    return sols


def build_basis(fm_exact, Nmax, dps_work, report):
    """All roots, grouped by congruence class. Returns list of solution dicts
    (len = ODE order at s=0 counted with multiplicity)."""
    ind = fm_exact[0]
    roots = rational_roots(ind)
    report.append(f"indicial roots at s=0: {roots}")
    classes = {}
    for r in roots:
        classes.setdefault(r % 1, []).append(r)
    sols = []
    for key, rts in sorted(classes.items()):
        # distinct roots with multiplicities
        dist = {}
        for r in rts:
            dist[r] = dist.get(r, 0) + 1
        for rho in sorted(dist):
            Prho = sum(m for r2, m in dist.items() if r2 > rho)
            sols.extend(frobenius_root(fm_exact, rho, dist[rho], Prho,
                                       Nmax, dps_work, report))
    return sols


# ---------- step 6: evaluation of solutions and derivatives ----------

def eval_sol_derivs(sol, s0, nder, dps_work):
    """[d^i/ds^i y](s0) for i = 0..nder-1.

    y = s^rho * sum_p log^p(s) * g_p(s),  g_p(s) = sum_N terms[(N,p)] s^N.
    D^i y = sum_p sum_{q<=i} C(i,q) * D^q[s^rho log^p](s0) * g_p^{(i-q)}(s0).
    """
    rho = sol['rho']
    rhof = mp.mpf(rho.numerator) / mp.mpf(rho.denominator)
    ln_s0 = mp.log(s0)
    # group terms by log power p
    byp = {}
    Nmax = 0
    for (N, p), coef in sol['terms'].items():
        byp.setdefault(p, {})[N] = coef
        Nmax = max(Nmax, N)
    pmax = max(byp.keys()) if byp else 0
    # powers of s0
    pw = [mp.mpf(1)] * (Nmax + 1)
    for N in range(1, Nmax + 1):
        pw[N] = pw[N - 1] * s0
    inv_s0 = 1 / s0
    # g_p^{(d)}(s0) for d = 0..nder-1
    geval = {}
    for p, terms in byp.items():
        row = []
        for d in range(nder):
            acc = mp.mpf(0)
            for N, coef in terms.items():
                if N < d:
                    continue
                ff = mp.mpf(1)
                for q in range(d):
                    ff *= (N - q)
                acc += coef * ff * pw[N - d]
            row.append(acc)
        geval[p] = row
    # base[q][p] = D^q [s^rho log^p s](s0), q = 0..nder-1, p = 0..pmax
    s0_rho = s0 ** rhof
    lnp = [mp.mpf(1)]
    for p in range(1, pmax + 1):
        lnp.append(lnp[-1] * ln_s0)
    base = [[s0_rho * lnp[p] for p in range(pmax + 1)]]
    # maintain expansion D^q[s^rho L^p] = s^{rho-q} sum_pp c_{q,p,pp} L^pp
    coefs = {p: {p: mp.mpf(1)} for p in range(pmax + 1)}  # p -> {pp: c}
    for q in range(1, nder):
        newc = {}
        for p, d in coefs.items():
            nd = {}
            for pp, c in d.items():
                b = rhof - (q - 1)
                nd[pp] = nd.get(pp, mp.mpf(0)) + c * b
                if pp >= 1:
                    nd[pp - 1] = nd.get(pp - 1, mp.mpf(0)) + c * pp
            newc[p] = nd
        coefs = newc
        s0_pow = s0 ** (rhof - q)
        base.append([sum(c * lnp[pp] for pp, c in coefs[p].items()) * s0_pow
                     for p in range(pmax + 1)])
    # binomials
    out = []
    for i in range(nder):
        acc = mp.mpf(0)
        for p in byp:
            for q in range(i + 1):
                acc += mp.binomial(i, q) * base[q][p] * geval[p][i - q]
        out.append(acc)
    return out


def sol_value_at_0(sol):
    """(value, divergent_flag): value = coeff of s^0 log^0; divergent if any term
    with beta<0 or (beta==0 and p>=1) has nonzero coefficient."""
    rho = sol['rho']
    val = mp.mpf(0)
    div = mp.mpf(0)  # magnitude of divergent content
    for (N, p), coef in sol['terms'].items():
        beta = rho + N  # Fraction + int = Fraction
        if beta < 0:
            div += abs(coef)
        elif beta == 0:
            if p == 0:
                val += coef
            else:
                div += abs(coef)
    return val, div


# ---------- step 6b: target-side exact-series values ----------

def series_derivs_mpf(wfr, x0, nder, dps_work):
    """F^{(i)}(x0), i=0..nder-1, from exact Fraction coefficients wfr[k]."""
    out = []
    x0 = mp.mpf(x0)
    K = len(wfr)
    wm = [mp.mpf(f.numerator) / mp.mpf(f.denominator) for f in wfr]
    for i in range(nder):
        acc = mp.mpf(0)
        # sum_k w_k k(k-1)...(k-i+1) x0^{k-i}, backwards for stability
        pw = x0 ** (K - 1 - i)
        inv = 1 / x0
        for k in range(K - 1, i - 1, -1):
            ff = mp.mpf(1)
            for q in range(i):
                ff *= (k - q)
            acc += wm[k] * ff * pw
            pw *= inv
        out.append(acc)
    return out


# ---------- top level ----------

def transport_value(rec_json, wfr, dps, s0=F(1, 2), s_check=F(2, 5), verbose=True):
    """Full pipeline. rec_json: dict {'r', 's', 'coeffs': {'j,k': [num,den]}}
    (fit_ode.py / annihilator.pf_from_series format). wfr: exact Fraction coeffs.
    Returns dict with value C = F(1), diagnostics, gates.

    WARNING (lesson 1, module docstring): the recurrence is assumed to
    hold for n >= 0 ONLY; the exact boundary inhomogeneity q(x) is built here
    and L' = D^{deg q + 1} o L is composed automatically. Do not bypass."""
    import os
    import sys
    # annihilator lives as a module inside tools/annihilator/ in this repo
    _ann = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "annihilator")
    if _ann not in sys.path:
        sys.path.insert(0, _ann)
    from annihilator import rec_to_theta

    report = []
    guard = 40
    dps_work = dps + guard
    with mp.workdps(dps_work):
        r, s = rec_json['r'], rec_json['s']
        coeffs = {tuple(map(int, k.split(','))): F(int(v[0]), int(v[1]))
                  for k, v in rec_json['coeffs'].items()}
        ode = rec_to_theta(r, s, coeffs)   # {(xpow,thetapow): F}
        order = max(j for (_, j) in ode)
        p = theta_to_D(ode, order)
        # --- boundary inhomogeneity: the recurrence holds for n>=0 only, so
        # L F = q(x) with q_m = sum_j c_j(m-r) a_{m-r+j} (m=0..r-1, indices>=0).
        # Compose L' = D^{deg q + 1} o L, which annihilates F exactly.
        qpoly = []
        for m in range(r):
            n = m - r
            qm = F(0)
            for (j, k2), c in coeffs.items():
                if n + j >= 0:
                    qm += c * wfr[n + j] * F(n) ** k2
            qpoly.append(qm)
        qpoly = ptrim(qpoly)
        if qpoly:
            ncomp = len(qpoly)  # D^{deg+1}
            report.append(f"boundary inhomogeneity deg {len(qpoly)-1}: "
                          f"composing D^{ncomp} o L")
            for _ in range(ncomp):
                newp = [[] for _ in range(len(p) + 1)]
                for i, pi in enumerate(p):
                    if not pi:
                        continue
                    dpi = ptrim([pi[q2] * q2 for q2 in range(1, len(pi))])
                    newp[i] = padd(newp[i], dpi)
                    newp[i + 1] = padd(newp[i + 1], list(pi))
                p = newp
            order += ncomp
        else:
            report.append("boundary inhomogeneity: none (L F = 0 already)")
        # leading coefficient structure
        from sympy import Poly, Rational, symbols, factor
        x = symbols('x')
        plead = sum(Rational(c.numerator, c.denominator) * x**i for i, c in enumerate(p[order]))
        report.append(f"ODE order {order}; leading coeff factors: {factor(plead)}")
        # shift to s = 1-x
        pt = [pscal(F((-1) ** i), subst_1_minus_s(p[i])) for i in range(order + 1)]
        fm = delta_form_at_s0(pt)
        # Frobenius basis
        Nmax = int(dps_work * mp.log(10) / mp.log(1 / float(s0))) + 120
        sols = build_basis(fm, Nmax, dps_work, report)
        nb = len(sols)
        if nb != order:
            report.append(f"WARN: basis size {nb} != order {order}")
        # matching at s0  (x0 = 1 - s0)
        x0 = F(1) - s0
        s0m = mp.mpf(s0.numerator) / mp.mpf(s0.denominator)
        s1m = mp.mpf(s_check.numerator) / mp.mpf(s_check.denominator)
        A = mp.matrix(nb, nb)
        for b, sol in enumerate(sols):
            vals = eval_sol_derivs(sol, s0m, nb, dps_work)
            for i in range(nb):
                A[i, b] = vals[i]
        rhs_x = series_derivs_mpf(wfr, mp.mpf(x0.numerator) / mp.mpf(x0.denominator),
                                  nb, dps_work)
        rhs = mp.matrix([((-1) ** i) * rhs_x[i] for i in range(nb)])
        cvec = mp.lu_solve(A, rhs)
        # truncation tail diagnostic
        tail = mp.mpf(0)
        for sol in sols:
            NN = max(N for (N, _p) in sol['terms']) if sol['terms'] else 0
            mx = max((abs(c) for (N, _p), c in sol['terms'].items() if N >= NN - 5),
                     default=mp.mpf(0))
            tail = max(tail, mx * s0m ** NN)
        report.append(f"series tail estimate at s0: {mp.nstr(tail, 3)}")
        # gate: independent check point
        x1 = F(1) - s_check
        chk_x = series_derivs_mpf(wfr, mp.mpf(x1.numerator) / mp.mpf(x1.denominator),
                                  2, dps_work)
        resid = []
        solvals1 = [eval_sol_derivs(sol, s1m, 2, dps_work) for sol in sols]
        for i in range(2):
            tv = mp.mpf(0)
            for b in range(nb):
                tv += cvec[b] * solvals1[b][i]
            resid.append(abs(tv - ((-1) ** i) * chk_x[i]) / max(abs(chk_x[i]), mp.mpf(1)))
        report.append(f"check-point rel residuals (val, D): "
                      f"{mp.nstr(resid[0], 3)}, {mp.nstr(resid[1], 3)}")
        # value + divergence gates
        Cval = mp.mpf(0)
        divmax = mp.mpf(0)
        scale = max(abs(c) for c in cvec)
        for b, sol in enumerate(sols):
            v, d = sol_value_at_0(sol)
            Cval += cvec[b] * v
            if d > 0:
                divmax = max(divmax, abs(cvec[b]) * d / scale)
                report.append(f"  divergent branch b={b} (rho={sol['rho']},j={sol['j']}): "
                              f"|c_b|*div/scale = {mp.nstr(abs(cvec[b]) * d / scale, 3)}")
        report.append(f"divergent-branch max rel coeff: {mp.nstr(divmax, 3)}")
        Cstr = mp.nstr(Cval, dps, strip_zeros=False)
    if verbose:
        for line in report:
            print("   ", line)
    return {'value': Cstr, 'report': report,
            'resid_check': [mp.nstr(rv, 3) for rv in resid],
            'div_gate': mp.nstr(divmax, 3), 'order': order, 'nsols': nb}


def load_series(path):
    out = []
    with open(path) as f:
        for line in f:
            k, nu, de = line.split()
            out.append(F(int(nu), int(de)))
    return out

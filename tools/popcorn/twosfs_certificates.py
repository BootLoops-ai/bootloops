#!/usr/bin/env python3
"""twosfs_certificates — exact certificates that place an exact 2-SFS inside or
outside classes of variable-population-size Kingman coalescents.

Register: EXACT for every decision. Certificates, margins, projections and
verdicts are fractions.Fraction / sympy.Rational end to end. Floats appear in
exactly one role: inside class_verdict/hunt a float LP (scipy HiGHS on a numpy
grid) PROPOSES a candidate witness or support; the candidate is rationalized
(Fraction(x).limit_denominator(10**6)) and only exact arithmetic DECIDES.
A wrong or suboptimal LP answer can cost a verdict (UNDECIDED), never produce
a false certificate. Requires sympy at import; numpy + scipy only when
class_verdict / hunt is called.

Classes and what each certificate means
---------------------------------------
Let M = (E[L_i L_j])_{i<=j} be the exact second-moment matrix of a coalescent
model at sample size n (twosfs_engine.lambda_moments) and q its
normalization. On the Kingman side (twosfs_engine docstring): every neutral
single-population Kingman coalescent with a deterministic size history has
E[L_i L_j] = int int Csym_ij(x,y) nu(dx) nu(dy) for a measure nu on (0,1],
and pooling loci across histories adds such terms. Three nested statements,
strongest first:

  FULL CLASS (verdict 'OUT_FULLCLASS'). A rational vector W = (W_ij)_{i<=j}
    with   K_W(x, y) := sum_{i<=j} W_ij Csym_ij(x, y) >= 0 on [0,1]^2
    (certified exactly by tensor-Bernstein expansion with quadtree
    subdivision on the (s, u) chart x = s, y = s u: bern_2d_nonneg) gives
    <W, M'> = int int K_W dnu dnu >= 0 for EVERY history and every pooling;
    if also <W, M> < 0 exactly, then M (equivalently q: the statement is
    scale free) lies outside the expected-2-SFS set of every pooled
    variable-size Kingman coalescent. No condition on the 1-SFS is involved.
  ATOM POOLINGS (verdict 'OUT_ATOMPOOLINGS'). The same with the diagonal
    only: p_W(t) := K_W(t, t) = sum W_ij m2_ij(t) >= 0 on [0,1]
    (bern_1d_nonneg on the exact polynomial from m2_polys) and <W, M> < 0
    excludes every pooling of single-atom limit histories (nu a finite
    positive combination of point masses, any support), but NOT genuinely
    time-varying histories: the constant-size Kingman coalescent itself is
    OUT_ATOMPOOLINGS at n = 4, 5, 6 and is never OUT_FULLCLASS (it is in the
    class) -- the shipped control rows.
  SUPPORT HULL (twosfs_hull.project_hull): exact projection of q onto the
    hull of normalized atoms at a FIXED finite support; weakest register.
  INSIDE (verdict 'INSIDE_2SFS'): an exact convex combination of normalized
    atom second moments over a finite dyadic t-grid equal to q (LP proposes
    the support, twosfs_hull.solve_exact reconstructs exact weights, every
    equation re-verified exactly): q is realized by an atom pooling, so no
    2-SFS certificate of either kind can exist.
  Otherwise 'UNDECIDED'. Verdicts are per target and never averaged. On the
  way (only when no OUT certificate was found), hunt reports whether M has
  an exact conic representation by kernel values Csym(x_k, y_k) on the
  dyadic grid of spacing 1/16 ('outer_membership'): if so, no pointwise-
  nonnegative K_W can separate M (a certified obstruction to FULL-CLASS
  linear certificates, not a membership proof).

Soundness scope (do not soften): K_W >= 0 pointwise certifies nonnegativity
on the cone spanned by kernel values, which contains the class cone; failing
to find such a W does not prove membership (copositive-type witnesses could
exist). Margins are certificate specific, not maximal: norm_margin =
-<W, M> / E[L_tot^2] and, divided by max|W_ij|, a certified lower bound on
the l1 distance (over ordered pairs) between q and any class member's q.
bern_1d_nonneg / bern_2d_nonneg are THREE-valued: True = certified >= 0,
False = a definite negative value was found at a dyadic probe point,
None = undecided at max depth (a polynomial that touches zero in the
interior cannot be certified by subdivision alone; see
popcorn.certificates.positivity route B for touch-root deflation).

Functions
  poly_coeffs_2d(expr, s, u) -> {(a, b): Fraction}   sympy poly -> exact
      power-basis dict (also the input format of
      popcorn.certificates.region.to_bernstein).
  bern_1d_nonneg(coeffs, max_depth=12)   {deg: Fraction} on [0, 1].
  bern_2d_nonneg(expr, s, u, max_depth=7)   sympy polynomial on [0, 1]^2.
  m2_polys(n, s) -> {(i,j): sympy poly in s}   unnormalized single-time second
      moments m2_ij(t) of twosfs_engine.component_moments as polynomials
      (the u = 1 diagonal of the two-time kernel; hunt uses them for the
      atom-pooling certificate).
  hunt(n, M, kern, m2p, s, u, grid_ns=33, grid_nu=17, tag='') -> dict with
      'verdict' in {OUT_FULLCLASS, OUT_ATOMPOOLINGS, INSIDE_2SFS, UNDECIDED}
      and, for OUT verdicts, 'witness_W' (rational strings keyed 'i,j'),
      'G_M_lambda_exact' (= <W, M>, < 0), 'norm_margin_exact',
      'delta_l1_lb_exact', display floats, 'eta' (interior slack asked of
      the LP), 'bernstein': 'certified'; for INSIDE, 'inside_support' and
      'inside_weights' (exact strings); diagnostic keys 'lp_full_delta',
      'lp_atom_delta', 'full_class_search', 'atom_poolings_search',
      'outer_lp', 'outer_membership', 'inside'. LP: maximize delta subject
      to K_W >= eta * (sum of kernels) on the grid, <W, q> <= -delta,
      |W_ij| <= 1, for eta in {0, 1e-3, 3e-3, 1e-2} until a rationalized W
      certifies.
  kingman_kernel(n) -> (kern, m2p, s, u)   cached twotime_kernel_poly(n) +
      m2_polys, with the diagonal identity kern(s, 1) == m2p asserted.
  class_verdict(n, lam, **kw) -> hunt(...) for the coalescent with rates
      lam (e.g. lambda_exact.beta_rate(Fraction(3, 2))).
  certify_witness(n, lam, W) -> dict   re-derive, exactly, what a given
      rational witness W proves for the model lam: <W, M>, the margins, and
      fresh bern_2d (full class) and bern_1d (atom poolings) certificates.
  load_reference() -> the shipped reference verdicts
      (reference/twosfs/class_certificates.json).

Reference data (reference/twosfs/class_certificates.json): 19 exact rows for
n = 3..6 -- Beta(2-alpha, alpha) and Dirac coalescents whose normalized
expected 1-SFS lies exactly INSIDE the variable-size Kingman class (explicit
convex certificates in reference/twosfs/support_hull.json), plus the
constant-size Kingman control at n = 4, 5, 6. Every Lambda row at n >= 4 is
OUT_FULLCLASS with a rational witness and margin (e.g. n = 4, Beta alpha =
1/2: margin 11272245653/852514432000 ~ 1.3e-2); the n = 3 rows and the
Kingman controls come out OUT_ATOMPOOLINGS only: no full-class witness is
found for them (for the controls none can exist, constant size being in the
class). Reading: coalescent models whose expected 1-SFS lies exactly inside
the variable-size Kingman class are separated from that whole class by
their exact linked 2-SFS from four samples on.

Cost: kingman_kernel(n) 0.03-0.2 s for n <= 6 (~1 s at n = 8); one
class_verdict at n = 4..6 typically 0.02-0.3 s; certifying a full-class
witness by bern_2d_nonneg is milliseconds at n <= 5 and up to a few seconds
at n = 6 (deeper subdivision). Memory tens of MB plus sympy/scipy.

CLI:  python3 twosfs_certificates.py [--n N] [--family kingman|beta|dirac|
        msprime_dirac] [--param P] [--out FILE]
  builds the kernel at n (with its exact self-checks), runs class_verdict for
  the requested model, prints the verdict, the exact margin and the witness,
  and writes the result JSON only if --out is given.

References: the Bernstein / de Casteljau subdivision test for polynomial
positivity on a box is classical; the LP-proposes / exact-decides discipline
and the region-certificate twin (to_bernstein / bern_nonneg on power-basis
dicts) are in popcorn.certificates. For the coalescent objects see the
references in twosfs_engine.
"""
import argparse
import functools
import json
import os
import sys
from fractions import Fraction as F
from math import comb
import sympy as sp
from twosfs_engine import (lambda_moments, component_moments, jump_levels,
                           death_semigroup_coeffs, twotime_kernel_poly, _N,
                           rate_function)
from twosfs_hull import solve_exact, pairs

HERE = os.path.dirname(os.path.abspath(__file__))
REFERENCE_DIR = os.path.join(HERE, 'reference', 'twosfs')
CERTIFICATES_REF = os.path.join(REFERENCE_DIR, 'class_certificates.json')

PAIRS = pairs

# ------------------------------------------------ exact Bernstein certificates

def poly_coeffs_2d(expr, s_, u_):
    p = sp.Poly(expr, s_, u_)
    d = {}
    for (a, b), c in p.terms():
        d[(a, b)] = F(int(sp.numer(c)), int(sp.denom(c)))
    return d

def bern_1d_nonneg(coeffs, depth=0, max_depth=12, lo=F(0), hi=F(1)):
    """coeffs: {deg: Fraction} of poly on [lo,hi] mapped to [0,1]. True if
    certified >= 0 on the interval (exact), False if a definite negative
    value found, None if inconclusive at max depth."""
    m = max(coeffs) if coeffs else 0
    c = [coeffs.get(a, F(0)) for a in range(m + 1)]
    b = [sum(F(comb(k, a), comb(m, a)) * c[a] for a in range(k + 1))
         for k in range(m + 1)]
    if all(x >= 0 for x in b):
        return True
    # definite negativity test at endpoints / midpoint
    def ev(x):
        return sum(cc * x ** a for a, cc in enumerate(c))
    for x in (F(0), F(1, 2), F(1)):
        if ev(x) < 0:
            return False
    if depth >= max_depth:
        return None
    # subdivide: left half p(x/2), right half p(1/2 + x/2)
    for half in (0, 1):
        sub = {}
        for a, cc in enumerate(c):
            if cc == 0:
                continue
            # substitute x = (half + y)/2, expand binomially
            for k in range(a + 1):
                sub[k] = sub.get(k, F(0)) + cc * F(comb(a, k), 2 ** a) * \
                    (F(half) ** (a - k) if a - k >= 0 else F(0))
        r = bern_1d_nonneg(sub, depth + 1, max_depth)
        if r is not True:
            return r
    return True

def bern_2d_nonneg(expr, s_, u_, max_depth=7):
    """Exact certificate expr >= 0 on [0,1]^2 via tensor Bernstein +
    quadtree subdivision (subdivision by exact sympy substitution)."""
    def rec(e, depth):
        cf = poly_coeffs_2d(e, s_, u_)
        if not cf:
            return True
        ms = max(a for a, _ in cf)
        mu = max(b for _, b in cf)
        C = [[cf.get((a, b), F(0)) for b in range(mu + 1)]
             for a in range(ms + 1)]
        ok = True
        for k in range(ms + 1):
            for l in range(mu + 1):
                bkl = sum(F(comb(k, a), comb(ms, a)) * F(comb(l, b), comb(mu, b))
                          * C[a][b] for a in range(k + 1) for b in range(l + 1))
                if bkl < 0:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            return True
        # definite negativity probes
        for xs in (F(0), F(1, 2), F(1)):
            for xu in (F(0), F(1, 2), F(1)):
                v = sum(c * xs ** a * xu ** b for (a, b), c in cf.items())
                if v < 0:
                    return False
        if depth >= max_depth:
            return None
        for hs in (0, 1):
            for hu in (0, 1):
                e2 = sp.expand(e.subs({s_: (sp.Integer(hs) + s_) / 2,
                                       u_: (sp.Integer(hu) + u_) / 2}))
                r = rec(e2, depth + 1)
                if r is not True:
                    return r
        return True
    return rec(sp.expand(expr), 0)

# ------------------------------------------------------- symbolic ingredients

def m2_polys(n, s_):
    """m2_ij(t) as sympy polys in s_ (unnormalized component second moments)."""
    co = death_semigroup_coeffs(n)
    levels = jump_levels(n)
    out = {}
    for i, j in PAIRS(n):
        acc = sp.Integer(0)
        for k in range(2, n + 1):
            m = sum(p * _N(a, i) * _N(a, j) for a, p in levels[k].items())
            if m == 0:
                continue
            pk = sum(sp.Rational(c) * s_ ** comb(jj, 2)
                     for jj, c in co[n].get(k, {}).items())
            acc += sp.Rational(m) * pk
        out[(i, j)] = sp.expand(acc)
    return out

# --------------------------------------------------------------- verdict per M

def hunt(n, MLam, kern, m2p, s_, u_, grid_ns=33, grid_nu=17, tag=''):
    """Full-class then atom-pooling witness search + exact verification, then
    the outer-cone obstruction and the exact INSIDE reconstruction, for the
    target second moments MLam (dict pairs -> Fraction); kern, m2p, s_, u_
    from kingman_kernel(n). Returns the result dict described in the module
    docstring. numpy + scipy (HiGHS) propose; exact arithmetic decides."""
    import numpy as np
    from scipy.optimize import linprog
    prs = PAIRS(n)
    d = len(prs)
    Ks = {ij: sp.lambdify((s_, u_), kern[ij], 'numpy') for ij in prs}
    sg = np.linspace(0, 1, grid_ns)
    ug = np.linspace(0, 1, grid_nu)
    SS, UU = np.meshgrid(sg, ug, indexing='ij')
    phi = np.zeros(SS.shape)
    vals = {}
    for ij in prs:
        v = np.asarray(Ks[ij](SS, UU), dtype=float) + 0 * SS
        vals[ij] = v
        phi += v * (1 if ij[0] == ij[1] else 2)
    totL = sum((1 if i == j else 2) * MLam[(i, j)] for i, j in prs)
    mlam = [float(MLam[ij] / totL) for ij in prs]     # normalized q_Lambda
    out = {'tag': tag}

    def try_out(pts_vals, pts_phi, eta, kernel_mode):
        """max delta st K_W >= eta*phi on grid, <W, q_Lambda> <= -delta,
        |W|<=1; exact-verify winner. Returns result dict or None."""
        A_ub, b_ub = [], []
        for kv, kp in zip(pts_vals, pts_phi):
            A_ub.append([-x for x in kv] + [0.0])
            b_ub.append(-eta * kp)
        A_ub.append(mlam + [1.0])
        b_ub.append(0.0)
        res = linprog(c=[0.0] * d + [-1.0], A_ub=A_ub, b_ub=b_ub,
                      bounds=[(-1, 1)] * d + [(0, None)],
                      method='highs')
        if res.status != 0 or res.x[-1] <= 1e-8:
            return None, (float(res.x[-1]) if res.status == 0 else None)
        W = {ij: F(float(res.x[k])).limit_denominator(10 ** 6)
             for k, ij in enumerate(prs)}
        GM = sum(W[ij] * MLam[ij] for ij in prs)
        if GM >= 0:
            return None, float(res.x[-1])
        if kernel_mode == '2d':
            expr = sum(sp.Rational(W[ij]) * kern[ij] for ij in prs)
            cert = bern_2d_nonneg(expr, s_, u_)
        else:
            pexpr = sp.expand(sum(sp.Rational(W[ij]) * m2p[ij] for ij in prs))
            cf = {}
            for (a,), c in sp.Poly(pexpr, s_).terms():
                cf[a] = F(int(sp.numer(c)), int(sp.denom(c)))
            cert = bern_1d_nonneg(cf)
        if cert is not True:
            return None, float(res.x[-1])
        winf = max(abs(x) for x in W.values())
        return ({'witness_W': {f'{i},{j}': str(W[(i, j)]) for i, j in prs},
                 'G_M_lambda_exact': str(GM),
                 'norm_margin_exact': str(-GM / totL),
                 'norm_margin_float': float(-GM / totL),
                 'delta_l1_lb_exact': str(-GM / totL / winf),
                 'delta_l1_lb_float': float(-GM / totL / winf),
                 'eta': eta, 'bernstein': 'certified'},
                float(res.x[-1]))

    # ---- full class: kernel on the square (every pooled Kingman history)
    pts_vals = [[vals[ij][a, b] for ij in prs]
                for a in range(SS.shape[0]) for b in range(SS.shape[1])]
    pts_phi = [phi[a, b] for a in range(SS.shape[0])
               for b in range(SS.shape[1])]
    for eta in (0.0, 1e-3, 3e-3, 1e-2):
        win, z = try_out(pts_vals, pts_phi, eta, '2d')
        if eta == 0.0:
            out['lp_full_delta'] = z
        if win:
            out.update(verdict='OUT_FULLCLASS', **win)
            return out
        if z is None or z <= 1e-8:
            break
    out['full_class_search'] = 'no_certified_witness'
    # ---- atom poolings: diagonal only (any support)
    dvals = {}
    phid = np.zeros(sg.shape)
    for ij in prs:
        v = np.asarray(Ks[ij](sg, np.ones_like(sg)), dtype=float) + 0 * sg
        dvals[ij] = v
        phid += v * (1 if ij[0] == ij[1] else 2)
    pts_vals = [[dvals[ij][a] for ij in prs] for a in range(len(sg))]
    pts_phi = [phid[a] for a in range(len(sg))]
    for eta in (0.0, 1e-3, 3e-3, 1e-2):
        win, z = try_out(pts_vals, pts_phi, eta, '1d')
        if eta == 0.0:
            out['lp_atom_delta'] = z
        if win:
            out.update(verdict='OUT_ATOMPOOLINGS', **win)
            return out
        if z is None or z <= 1e-8:
            break
    out['atom_poolings_search'] = 'no_certified_witness'
    # ---- outer-relaxation membership: M_Lambda in cone{Csym(x,y)}?
    # If YES (exact conic witness), then NO pointwise-nonneg kernel
    # functional can separate: a full-class linear certificate is impossible
    # for this target, not merely unfound -- a certified obstruction.
    # (Membership of the OUTER cone does not imply class membership; it
    # rules out linear certificates.)
    import itertools as _it
    kdict = []          # (label, {pairs: Fraction})
    for a in range(0, 17):
        for b in range(0, a + 1):
            x, y = F(a, 16), F(b, 16)
            col = {}
            for ij in prs:
                e = kern[ij]
                val = e.subs({s_: sp.Rational(x),
                              u_: sp.Rational(y / x) if x != 0 else 0})
                col[ij] = F(int(sp.numer(val)), int(sp.denom(val)))
            kdict.append((f'({a}/16,{b}/16)', col))
    cols = [c for _, c in kdict]
    Aeq = np.array([[float(c[ij]) for c in cols] for ij in prs])
    beq = np.array([float(MLam[ij]) for ij in prs])
    res = linprog(c=np.zeros(len(cols)), A_eq=Aeq, b_eq=beq,
                  bounds=[(0, None)] * len(cols), method='highs')
    out['outer_lp'] = int(res.status)
    if res.status == 0:
        sup = sorted(range(len(cols)), key=lambda k: -res.x[k])[:2 * d]
        found = None
        tries = 0
        for S in _it.combinations(sup, d):
            tries += 1
            if tries > 3000:
                break
            A = [[cols[k][ij] for k in S] for ij in prs]
            b = [MLam[ij] for ij in prs]
            s0 = solve_exact(A, b)
            if s0 is None or any(x < 0 for x in s0):
                continue
            found = (S, s0)
            break
        if found:
            S, s0 = found
            out['outer_membership'] = {
                'generators': [kdict[k][0] for k in S],
                'weights': [str(x) for x in s0],
                'meaning': 'M_Lambda in cone{Csym points} EXACTLY: no '
                           'pointwise-nonneg kernel witness exists (no '
                           'full-class linear certificate exists)'}
        else:
            out['outer_membership'] = 'float_feasible_exact_open'
    # ---- INSIDE hunt: exact convex combination over a fine t-grid
    tgrid = [F(a, 64) for a in range(0, 65)] + \
            [F(1) - F(1, 1 << k) for k in range(7, 16)] + \
            [F(1, 1 << k) for k in range(7, 16)]
    tgrid = sorted(set(tgrid))
    qL = {}
    totL = sum((1 if i == j else 2) * MLam[(i, j)] for i, j in prs)
    for ij in prs:
        qL[ij] = MLam[ij] / totL
    comps = []
    for t in tgrid:
        _, m2 = component_moments(n, t)
        tot = sum((1 if i == j else 2) * m2[(i, j)] for i, j in prs)
        comps.append({ij: m2[ij] / tot for ij in prs})
    Aeq = np.array([[float(c[ij]) for c in comps] for ij in prs]
                   + [[1.0] * len(comps)])
    beq = np.array([float(qL[ij]) for ij in prs] + [1.0])
    res = linprog(c=np.zeros(len(comps)), A_eq=Aeq, b_eq=beq,
                  bounds=[(0, None)] * len(comps), method='highs')
    if res.status == 0:
        sup = [k for k in range(len(comps)) if res.x[k] > 1e-12]
        # exact solve on proposed support
        k_ = len(sup)
        rowsA = [[comps[k][ij] for k in sup] for ij in prs] + [[F(1)] * k_]
        rhs = [qL[ij] for ij in prs] + [F(1)]
        # least-squares-free exact: solve square subsystem, verify all rows
        import itertools
        sol = None
        for rows_idx in itertools.combinations(range(len(rowsA)), k_):
            A = [rowsA[r] for r in rows_idx]
            b = [rhs[r] for r in rows_idx]
            s0 = solve_exact(A, b)
            if s0 is None:
                continue
            if any(x < 0 for x in s0):
                continue
            if all(sum(rowsA[r][c] * s0[c] for c in range(k_)) == rhs[r]
                   for r in range(len(rowsA))):
                sol = s0
                break
        if sol is not None:
            out.update({
                'verdict': 'INSIDE_2SFS',
                'inside_support': [str(tgrid[k]) for k in sup],
                'inside_weights': [str(x) for x in sol]})
            return out
        out['inside'] = 'float_feasible_exact_reconstruction_failed'
    else:
        out['inside'] = 'lp_infeasible'
    out['verdict'] = 'UNDECIDED'
    return out

# ------------------------------------------------ convenience front doors

@functools.lru_cache(maxsize=None)
def kingman_kernel(n):
    """(kern, m2p, s, u): the two-time kernel polynomials of
    twosfs_engine.twotime_kernel_poly(n), the single-time second-moment
    polynomials m2_polys(n, s), and the sympy symbols; built once per n.
    Asserts the exact diagonal identity kern_ij(s, u=1) == m2p_ij(s)."""
    kern, (s_, u_) = twotime_kernel_poly(n)
    m2p = m2_polys(n, s_)
    for ij in PAIRS(n):
        assert sp.expand(kern[ij].subs(u_, 1) - m2p[ij]) == 0, (n, ij)
    return kern, m2p, s_, u_

def kernel_constant_size_check(n):
    """Exact identity: integrating the kernel against nu(dx) = dx/x (constant
    population size) reproduces the Kingman E[L_i L_j] of lambda_moments
    entrywise, 2 * sum_{a,b} c_ab / (a b) over the monomials c_ab s^a u^b of
    Csym_ij (every monomial must have a, b >= 1). Returns True/False."""
    from twosfs_engine import kingman_rate
    kern, _m2p, s_, u_ = kingman_kernel(n)
    _u, v = lambda_moments(n, kingman_rate())
    ok = True
    for ij in PAIRS(n):
        acc = F(0)
        for (a, b), c in poly_coeffs_2d(sp.expand(kern[ij]), s_, u_).items():
            if c == 0:
                continue
            if a < 1 or b < 1:
                return False
            acc += 2 * c * F(1, a * b)
        ok &= (acc == v[ij])
    return bool(ok)

def class_verdict(n, lam, **kw):
    """hunt() for the coalescent with merger rates lam at sample size n
    (kernel cached per n). Needs numpy + scipy for the proposing LPs."""
    kern, m2p, s_, u_ = kingman_kernel(n)
    _u, v = lambda_moments(n, lam)
    return hunt(n, v, kern, m2p, s_, u_, **kw)

def certify_witness(n, lam, W):
    """What the rational witness W (dict (i,j) -> Fraction or 'i,j' -> str,
    i <= j) proves, re-derived exactly for the model lam at sample size n:
      {'G_M_exact': <W, M>, 'norm_margin_exact', 'delta_l1_lb_exact'
       (Fractions; margins are None unless <W, M> < 0),
       'full_class': bern_2d_nonneg(K_W)      (True / False / None),
       'atom_poolings': bern_1d_nonneg(p_W)   (True / False / None),
       'verdict': 'OUT_FULLCLASS' | 'OUT_ATOMPOOLINGS' | 'NOT_A_CERTIFICATE'}.
    No LP, no floats: sympy + Fractions only."""
    prs = PAIRS(n)
    Wf = {}
    for k, val in W.items():
        ij = tuple(int(x) for x in k.split(',')) if isinstance(k, str) else tuple(k)
        Wf[ij] = F(val)
    assert set(Wf) == set(prs), f"witness keys {sorted(Wf)} != pairs({n})"
    kern, m2p, s_, u_ = kingman_kernel(n)
    _u, v = lambda_moments(n, lam)
    GM = sum(Wf[ij] * v[ij] for ij in prs)
    tot = sum((1 if i == j else 2) * v[(i, j)] for i, j in prs)
    winf = max(abs(x) for x in Wf.values())
    expr = sum(sp.Rational(Wf[ij].numerator, Wf[ij].denominator) * kern[ij]
               for ij in prs)
    full = bern_2d_nonneg(expr, s_, u_)
    pexpr = sp.expand(sum(sp.Rational(Wf[ij].numerator, Wf[ij].denominator)
                          * m2p[ij] for ij in prs))
    cf = {a: F(int(sp.numer(c)), int(sp.denom(c)))
          for (a,), c in sp.Poly(pexpr, s_).terms()}
    atom = bern_1d_nonneg(cf) if cf else True
    neg = GM < 0
    out = {'G_M_exact': GM,
           'norm_margin_exact': (-GM / tot) if neg else None,
           'delta_l1_lb_exact': (-GM / tot / winf) if (neg and winf) else None,
           'full_class': full, 'atom_poolings': atom}
    out['verdict'] = ('OUT_FULLCLASS' if (neg and full is True) else
                      'OUT_ATOMPOOLINGS' if (neg and atom is True) else
                      'NOT_A_CERTIFICATE')
    return out

def load_reference(path=None):
    """The shipped reference verdicts (dict with 'about' and 'rows')."""
    with open(path or CERTIFICATES_REF) as fh:
        return json.load(fh)

# ---------------------------------------------------------------- main

def main(argv=None):
    ap = argparse.ArgumentParser(
        description="exact 2-SFS class verdict for a Lambda-coalescent against "
                    "the variable-size Kingman class")
    ap.add_argument('--n', type=int, default=4, help="sample size (default 4)")
    ap.add_argument('--family', default='beta',
                    choices=['kingman', 'beta', 'dirac', 'msprime_dirac'])
    ap.add_argument('--param', default='1/2',
                    help="rational parameter (alpha | psi | 'psi,c'); "
                         "ignored for kingman (default 1/2)")
    ap.add_argument('--out', default=None, help="write the result JSON here")
    a = ap.parse_args(argv)
    n = a.n
    param = None if a.family == 'kingman' else a.param
    kingman_kernel(n)
    print(f"kernel n={n}: diagonal identity holds; constant-size "
          f"reconstruction {kernel_constant_size_check(n)}")
    r = class_verdict(n, rate_function(a.family, param))
    r.update({'n': n, 'family': a.family, 'param': param})
    tag = a.family + (f"({param})" if param is not None else '')
    print(f"n = {n}, {tag}: verdict {r['verdict']}")
    if r['verdict'].startswith('OUT'):
        print(f"  <W, M> = {r['G_M_lambda_exact']}  (< 0)")
        print(f"  normalized margin = {r['norm_margin_exact']} "
              f"~ {r['norm_margin_float']:.4e};  l1 lower bound "
              f"{r['delta_l1_lb_float']:.4e}")
        print(f"  witness W = {r['witness_W']}")
    elif r['verdict'] == 'INSIDE_2SFS':
        print(f"  exact atom pooling: support t = {r['inside_support']}, "
              f"weights {r['inside_weights']}")
    if a.out:
        with open(a.out, 'w') as fh:
            json.dump(r, fh, indent=1, default=str)
        print('wrote', a.out)
    return 0

if __name__ == '__main__':
    sys.exit(main())

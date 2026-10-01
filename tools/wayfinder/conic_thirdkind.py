#!/usr/bin/env python3
r"""
conic_thirdkind.py — the conic third-kind reduction core (wayfinder
member): builder (_pf_w1 / _herm / build_row_system) + the ExtSys numeric
wrapper + the three c2 engines polinv / odd_deriv / hermite_conic_sqf, in
one importable home.

Validated: 16/16 (row,pin) gates 35-37d (per-state march-vs-
quadrature 41-63d), mutation control fails as required, detour certificate
green.
Member battery: tests/test_conic_thirdkind.py — toy truth case (a small
synthetic fenced row rebuilt from tests/fixtures_conic_thirdkind/
toy_family.json and compared entry-exact against toy_banked.json, then
cross-checked row-by-row against independent quadrature at two
precisions), comparator + certificate mutation controls, deg-2 refusal
leg, ExtSys numeric leg. All legs run in seconds.

Adaptations in this copy relative to the validated original (ONLY these;
the algebra lines are verbatim):
  * symbols PARAMETERIZED — the fiber variable x and the
    connection variable t are explicit arguments, not module
    globals;
  * panel endpoints (A, B) are build_row_system arguments (originally
    hard-coded);
  * connection-input schema is an explicit adapter, rows_from_connection(),
    over the producing driver's JSON class, with a `subs` dict for pinning
    spectator variables (pin FIRST, then cancel — a measured lesson);
  * SCOPE NAMED, fail-closed: hermite_conic_sqf REFUSES deg(q, x) != 2
    with ConicScopeError (see SCOPE below) instead of the original's
    silent lc = None + downstream TypeError.

SCOPE — THE CONIC (deg-2) ASSUMPTION IS LOAD-BEARING
====================================================
The Hermite tail (polynomial-part absorption) eliminates the leading term
of the remaining polynomial P via d/dx(x^j Y) = (j x^{j-1} q + x^j q'/2)/Y,
whose leading coefficient is (j+1)*lc ONLY when deg(q, x) == 2. For
deg(q) > 2 (elliptic / hyperelliptic Y) the polynomial part is NOT
reducible to {d(H*Y), c0/Y, third-kind T/Y}: genuinely new second-kind
states x^i/Y (i = 1 .. deg(q)-2) survive, so the (H, c, T) return contract
and the builder's single-M connection column both change shape. That is a
basis EXTENSION, not a parameterization — deliberately NOT lifted here
(an "exactly this
atom class" claim would gloss this restriction). The refusal leg in the member
battery pins the fence. Second scope fence, also named: w1-dependent
denominator factors of T must have multiplicity 1 (squarefree third-kind
census; multiplicity >= 2 would need the Hermite loop extension — the
builder's assert fires loudly, as in the original driver).

WHAT IT DOES (the reduction stage, symbols generalized)
=======================================================
Input: a THIRD-KIND FENCED conic flat-connection row
    dL/dt = [H_L * Y]_A^B + c_L * M + INT_A^B T_L / Y dx ,   Y^2 = q(x, t)
with q quadratic in x, T_L != 0, T_L's x-dependent denominator factors
squarefree, and the t-poles MOVING over the family. It builds the EXACT
extended t-system over states [L, M, BcA, BcB, J_0.., F]: one third-kind
atom J = INT_A^B x^k/(f*Y) dx per (irreducible x-dependent factor f of
den T_L, k < deg f), endpoint carriers BcX = sqrt(q(X, t)), accumulator
F' = L. Row closures are exact: T_L partial-fractioned in x (residual==0
cert), polynomial parts + J-derivatives Hermite-reduced (cancel==0 certs;
J-row input sig = -x^k f_t/f^2 - x^k q_t/(2 f q)). The result is an exact
rational t-connection A(t) ready for wayfinder.transport_fixed_eps
(declared singular points = denominator roots = the moving poles; ExtSys
below is the duck-typed DESystem). Anchors/gates: quad.quad_refine.
The DE-marcher half of the original driver (window pick + march + gates,
across-pole detour certificates) is deliberately NOT included
— this member is the reduction core only.

Import: `from wayfinder import conic_thirdkind` (lazy — package import
stays sympy-free; sympy loads when THIS module does, frob_scalar pattern).
"""
import sys
import time

import sympy as sp

if sys.get_int_max_str_digits() < 4000000:
    sys.set_int_max_str_digits(4000000)

__all__ = [
    "ConicScopeError", "polinv", "odd_deriv", "hermite_conic_sqf",
    "partial_fractions", "hermite_certified", "load_expr",
    "rows_from_connection", "build_row_system", "ExtSys",
]


class ConicScopeError(NotImplementedError):
    """An input left the certified conic (deg-2) scope — see module SCOPE."""


def _require_conic(q, x):
    d = sp.degree(sp.expand(q), x)
    if d != 2:
        raise ConicScopeError(
            "hermite_conic_sqf: deg(q, %s) = %s; certified scope is the "
            "CONIC deg-2 class only (Y^2 = quadratic in the fiber "
            "variable). deg > 2 needs a second-kind basis extension "
            "(module docstring, SCOPE)." % (x, d))


# ---- c2 engines (verbatim; provenance twocurve_w2 via c2_parity.py) ---------
def polinv(a, f, x):
    """Inverse of a mod f in K[x]/(f), K the coefficient field."""
    a = sp.rem(sp.expand(a), f, x)
    s_, t_, g = sp.gcdex(sp.Poly(a, x), sp.Poly(f, x))
    assert g.degree() == 0, 'polinv: not coprime'
    return sp.expand((s_ / g.all_coeffs()[0]).as_expr())


def odd_deriv(gam, D, x):
    """d/dx (gam * Y) collapsed onto 1/Y, for Y^2 = D: returns the
    numerator N with d(gam*Y)/dx = N/Y."""
    return sp.cancel(sp.diff(gam, x) * D + gam * sp.diff(D, x) / 2)


def hermite_conic_sqf(sig, q, x):
    """sig/Y = d(H*Y) + c0/Y + T/Y over the quadratic conic Y^2 = q; exact,
    factorization-free (squarefree + branch-gcd only). SCOPE: deg(q,x)==2
    enforced (ConicScopeError otherwise)."""
    _require_conic(q, x)
    qp = sp.diff(q, x)
    Hnum = sp.Integer(0)
    sig = sp.cancel(sp.together(sig))
    nsteps = 0
    while True:
        nsteps += 1
        assert nsteps < 400, 'hermite_sqf: no convergence'
        N, Dn = sp.fraction(sp.cancel(sig))
        N, Dn = sp.expand(N), sp.expand(Dn)
        cont, sq = sp.sqf_list(sp.Poly(Dn, x))
        parts = []
        for f, mm in sq:
            f = f.as_expr()
            if sp.degree(f, x) == 0:
                continue
            fb = sp.gcd(sp.Poly(f, x), sp.Poly(q, x)).as_expr()
            if sp.degree(fb, x) > 0 and sp.degree(fb, x) < sp.degree(f, x):
                parts += [(sp.cancel(fb), mm), (sp.cancel(f / fb), mm)]
            else:
                parts.append((f, mm))
        todo = None
        for f, mm in parts:
            branch = sp.degree(sp.gcd(sp.Poly(f, x), sp.Poly(q, x)).as_expr(),
                               x) > 0
            if branch or mm >= 2:
                todo = (f, mm, branch)
                break
        if todo is None:
            break
        f, mm, branch = todo
        fp = sp.diff(f, x)
        rest = sp.cancel(Dn / f**mm)
        if not branch:
            S = sp.rem(sp.expand(-N * polinv(sp.expand(rest * (mm - 1) * fp * q),
                                             f, x)), f, x)
            dG = sp.cancel((sp.diff(S, x) * q * f + S * qp * f / 2
                            - (mm - 1) * S * fp * q) / f**mm)
            Hnum += sp.cancel(S / f**(mm - 1))
        else:
            if mm >= 2:
                S = sp.rem(sp.expand(-N * polinv(
                    sp.expand(rest * (sp.Rational(2 * mm - 1, 2)) * fp),
                    f, x)), f, x)
                dG = sp.cancel((sp.diff(S, x) * f * q - (mm - 1) * S * fp * q
                                - S * qp * f / 2) / (f**mm * q))
                Hnum += sp.cancel(S / (f**(mm - 1) * q))
            else:
                S = sp.rem(sp.expand(-2 * N * polinv(sp.expand(rest * fp),
                                                     f, x)), f, x)
                dG = sp.cancel((sp.diff(S, x) * q - S * qp / 2) / q)
                Hnum += sp.cancel(S / q)
        sig = sp.cancel(sig - dG)
    N, Dn = sp.fraction(sp.cancel(sig))
    P, Rr = sp.div(sp.expand(N), sp.expand(Dn), x)
    T = sp.cancel(Rr / Dn) if Rr != 0 else sp.Integer(0)
    P = sp.expand(P)
    lc = sp.LC(sp.Poly(q, x))          # deg-2 guaranteed by _require_conic
    while sp.degree(P, x) >= 1:
        d_ = sp.degree(P, x)
        j = d_ - 1
        lead = sp.LC(sp.Poly(P, x))
        cfac = sp.cancel(lead / ((j + 1) * lc))
        P = sp.expand(P - cfac * (j * x**(j - 1) * q + x**j * qp / 2))
        Hnum += cfac * x**j
    return sp.cancel(Hnum), sp.cancel(P), T


# ---- exact certified layers (verbatim from mpe_lib) -------------------------
def partial_fractions(T, factors, x):
    """Exact squarefree partial fractions of T over the given x-dependent
    irreducible factors (each multiplicity 1): returns
    (poly_part, {f_srepr: numer_poly (deg < deg f)}), exact residual cert."""
    n_, d_ = sp.fraction(sp.cancel(sp.together(T)))
    res = {}
    acc = sp.Integer(0)
    for f in factors:
        quo, rem_ = sp.div(sp.expand(d_), sp.expand(f), x)
        if sp.cancel(rem_) != 0:          # f does not divide den(T)
            res[sp.srepr(f)] = sp.Integer(0)
            continue
        rest = sp.cancel(quo)
        # R_f = n_/rest mod f  (inverse in K[x]/(f))
        Rf = sp.rem(sp.expand(sp.rem(sp.expand(n_), f, x)
                              * polinv(rest, f, x)), f, x)
        Rf = sp.cancel(Rf)
        res[sp.srepr(f)] = Rf
        acc += Rf / f
    poly = sp.cancel(T - acc)
    pn, pd = sp.fraction(sp.together(poly))
    assert not sp.expand(pd).has(x), \
        'partial_fractions: leftover %s denominator %s' % (x, str(pd)[:120])
    # exact residual certificate
    cert = sp.cancel(T - poly - acc)
    assert cert == 0, 'partial_fractions: residual nonzero'
    return poly, res


def hermite_certified(sig, Q, x):
    """Hermite-reduce sig/Y with exact certificate; returns (H, c, T)."""
    H_, c_, T_ = hermite_conic_sqf(sp.cancel(sp.together(sig)), Q, x)
    cert = sp.cancel(odd_deriv(H_, Q, x) + c_ + T_ - sig)
    assert cert == 0, 'hermite certificate nonzero'
    return H_, c_, T_


# ---- connection-input schema adapter ----------------------------------------
def load_expr(srep, symbols):
    """sympify a string/srepr and rebind its free symbols BY NAME to the
    caller's symbol table (assumption-safe cross-load).
    `symbols` = {name: Symbol}."""
    e = sp.sympify(srep)
    return e.xreplace({s_: symbols[s_.name] for s_ in e.free_symbols
                       if s_.name in symbols})


def rows_from_connection(conn, tag, mtag, cap, var_key, symbols, subs=None):
    """Adapter for the saved connection JSON class (schema
    keys parameterized):

        conn = {'rows': {tag: {var_key: {'H','c','T', ...}}},
                'Q': {cap: srep}, 'numerators': {tag: srep} (optional)}

    `var_key` names the connection variable's row block (e.g. 'z24');
    `subs` is an exact substitution dict applied to every expression
    (spectator pins, e.g. {u: Rational(1)}) — pin FIRST, then cancel
    (a measured ordering lesson). Returns (row, mrow, Q, numer) as exact sympy exprs;
    asserts the M row is T-free."""
    subs = subs or {}

    def _pin(v):
        return sp.cancel(load_expr(v, symbols).subs(subs))

    rowL = conn['rows'][tag][var_key]
    rowM = conn['rows'][mtag][var_key]
    row = {k: _pin(rowL[k]) for k in ('H', 'c', 'T')}
    mrow = {k: _pin(rowM[k]) for k in ('H', 'c')}
    assert load_expr(rowM['T'], symbols) == 0, \
        'rows_from_connection: M row carries a T block'
    Q = _pin(conn['Q'][cap])
    numer = (_pin(conn['numerators'][tag])
             if tag in conn.get('numerators', {}) else None)
    return row, mrow, Q, numer


# ---- the builder (verbatim from mpe_lib.build_row_system, parameterized) ----
def build_row_system(row, mrow, Q, x, t, endpoints, tag='row', numer=None,
                     verbose=True):
    """Build the extended exact t-system for a third-kind fenced conic row.

    row  = {'H','c','T'} exact sympy exprs in (x, t) (spectators pinned);
    mrow = {'H','c'} for the M row (its T must be 0 — adapter asserts);
    Q    = the conic, quadratic in x; endpoints = (A, B) exact panel ends;
    x, t = fiber / connection symbols.

    Returns dict: states (names), Amat (sympy Matrix, exact in t),
    atoms [(f_expr, k)], Q, x, t, endpoints, walls.
    State order: [L, M, BcA, BcB, J_0..J_{m-1}, F].
    """
    t0 = time.time()
    A_END, B_END = sp.sympify(endpoints[0]), sp.sympify(endpoints[1])
    HL, cL, TL = (sp.cancel(row[k]) for k in ('H', 'c', 'T'))
    HM, cM = sp.cancel(mrow['H']), sp.cancel(mrow['c'])
    _require_conic(Q, x)

    # ---- atom factor census -------------------------------------------------
    n_, d_ = sp.fraction(sp.cancel(sp.together(TL)))
    cont, fl = sp.factor_list(sp.expand(d_))
    fset = []
    for f, m in fl:
        fe = sp.expand(f.as_expr())
        if fe.has(x):
            assert m == 1, ('%s-dependent factor with mult %d — squarefree '
                            'third-kind scope only (Hermite loop extension '
                            'not built)' % (x, m))
            fset.append(fe)
    atoms = []
    for f in fset:
        for k in range(sp.degree(sp.Poly(f, x), x)):
            atoms.append((f, k))
    nJ = len(atoms)
    states = ['L', 'M', 'BcA', 'BcB'] + \
             ['J%d' % i for i in range(nJ)] + ['F']
    n = len(states)
    idx = {s_: i for i, s_ in enumerate(states)}
    Amat = sp.zeros(n, n)

    def add_boundary(i, H):
        """row i += [H*Y]_A^B  = H(B)*BcB - H(A)*BcA."""
        if H == 0:
            return
        Amat[i, idx['BcB']] += sp.cancel(H.subs(x, B_END))
        Amat[i, idx['BcA']] -= sp.cancel(H.subs(x, A_END))

    def add_T_as_atoms(i, T, allow_herm_loop=True):
        """row i += INT T/Y: PF into atoms; poly part Hermite-reduced."""
        if T == 0:
            return
        # only factors of den(T) that are in fset may appear
        nT, dT = sp.fraction(sp.cancel(sp.together(T)))
        _, flT = sp.factor_list(sp.expand(dT))
        for f, m in flT:
            fe = sp.expand(f.as_expr())
            if fe.has(x):
                match = [g for g in fset if sp.cancel(
                    sp.div(sp.expand(g), fe, x)[1]) == 0
                    and sp.degree(sp.Poly(g, x), x)
                    == sp.degree(sp.Poly(fe, x), x)]
                assert match, ('T-block factor outside atom set: %s'
                               % str(fe)[:120])
        poly, res = partial_fractions(T, fset, x)
        for f in fset:
            Rf = res[sp.srepr(f)]
            if Rf == 0:
                continue
            dfw = sp.degree(sp.Poly(f, x), x)
            Rp = sp.Poly(sp.expand(Rf), x)
            for k in range(dfw):
                ck = sp.cancel(Rp.coeff_monomial(x ** k) if k else
                               Rp.coeff_monomial(sp.Integer(1)))
                if ck != 0:
                    Amat[i, idx['J%d' % atoms.index((f, k))]] += ck
        if poly != 0:
            H2, c2, T2 = hermite_certified(poly, Q, x)
            assert T2 == 0 or allow_herm_loop, 'poly Hermite left remainder'
            add_boundary(i, H2)
            Amat[i, idx['M']] += c2
            if T2 != 0:
                add_T_as_atoms(i, T2, allow_herm_loop=False)

    # ---- L row --------------------------------------------------------------
    iL = idx['L']
    Amat[iL, idx['M']] += cL
    add_boundary(iL, HL)
    add_T_as_atoms(iL, TL)

    # ---- M row (stored, T==0) ----------------------------------------------
    iM = idx['M']
    Amat[iM, idx['M']] += cM
    add_boundary(iM, HM)

    # ---- Bc rows ------------------------------------------------------------
    for name, wv in (('BcA', A_END), ('BcB', B_END)):
        qX = sp.cancel(Q.subs(x, wv))
        Amat[idx[name], idx[name]] = sp.cancel(
            sp.diff(qX, t) / (2 * qX))

    # ---- J rows -------------------------------------------------------------
    Qz = sp.diff(Q, t)
    for a_i, (f, k) in enumerate(atoms):
        i = idx['J%d' % a_i]
        fz = sp.diff(f, t)
        sig = sp.cancel(-x ** k * fz / f ** 2
                        - x ** k * Qz / (2 * f * Q))
        Hjk, cjk, Tjk = hermite_certified(sig, Q, x)
        add_boundary(i, Hjk)
        Amat[i, idx['M']] += cjk
        add_T_as_atoms(i, Tjk)
        if verbose:
            print('  J row %d/%d (deg f=%d, k=%d) closed'
                  % (a_i + 1, nJ, sp.degree(sp.Poly(f, x), x), k),
                  flush=True)

    # ---- F row --------------------------------------------------------------
    Amat[idx['F'], idx['L']] = sp.Integer(1)

    Amat = Amat.applyfunc(sp.cancel)
    return {'tag': tag, 'states': states, 'atoms': atoms, 'Q': Q,
            'Amat': Amat, 'x': x, 't': t, 'endpoints': (A_END, B_END),
            'HL': HL, 'cL': cL, 'TL': TL,
            'numer': sp.cancel(numer) if numer is not None else None,
            'wall_s': round(time.time() - t0, 1)}


# ---- numeric wrapper for wayfinder (verbatim, symbol-parameterized) -------
class ExtSys:
    """Duck-typed DESystem: exact-coefficient univariate rational entries,
    numeric Horner at working precision (famhar _RatSegDE pattern).
    Feed to wayfinder.transport_fixed_eps; .singular_points = declared
    denominator roots (the MOVING POLES of the pinned family)."""

    def __init__(self, sysd, dps_build=40):
        import mpmath as mp
        self.mp = mp
        t = sysd['t']
        self.n = len(sysd['states'])
        self.var = str(t)
        self.meta = {'tag': sysd.get('tag', 'row')}
        self._ent = {}
        dens = set()
        for i in range(self.n):
            for j in range(self.n):
                e = sysd['Amat'][i, j]
                if e == 0:
                    continue
                nu, de = sp.fraction(sp.together(e))
                np_ = sp.Poly(sp.expand(nu), t)
                dp_ = sp.Poly(sp.expand(de), t)
                self._ent[(i, j)] = (
                    [sp.Rational(c) for c in np_.all_coeffs()],
                    [sp.Rational(c) for c in dp_.all_coeffs()])
                if dp_.degree() > 0:
                    dens.add(dp_)
        # singular points: roots of every distinct irreducible den factor
        sings = set()
        with mp.workdps(dps_build):
            facs = set()
            for dp_ in dens:
                for f, _m in sp.factor_list(dp_.as_expr())[1]:
                    if f.as_expr().has(t):
                        facs.add(sp.Poly(sp.expand(f.as_expr()), t))
            for f in facs:
                cs = [mp.mpf(int(sp.numer(c))) / mp.mpf(int(sp.denom(c)))
                      for c in f.all_coeffs()]
                try:
                    for r in mp.polyroots(cs, maxsteps=200,
                                          extraprec=120):
                        sings.add(complex(r))
                except Exception:
                    pass
        self.singular_points = [mp.mpc(s_) for s_ in sorted(
            sings, key=lambda c: (c.real, c.imag))]

    def A(self, x, eps, dps):
        mp = self.mp
        with mp.workdps(dps + 15):
            xv = mp.mpc(x)
            zero = mp.mpc(0)
            rows = [[zero] * self.n for _ in range(self.n)]
            cache = {}
            for (i, j), (nc, dc) in self._ent.items():
                key_n = id(nc)
                if key_n not in cache:
                    v = mp.mpc(0)
                    for c in nc:
                        v = v * xv + mp.mpf(int(sp.numer(c))) \
                            / mp.mpf(int(sp.denom(c)))
                    cache[key_n] = v
                key_d = id(dc)
                if key_d not in cache:
                    v = mp.mpc(0)
                    for c in dc:
                        v = v * xv + mp.mpf(int(sp.numer(c))) \
                            / mp.mpf(int(sp.denom(c)))
                    cache[key_d] = v
                rows[i][j] = cache[key_n] / cache[key_d]
        return rows

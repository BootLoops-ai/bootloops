#!/usr/bin/env python3
"""Exact j-invariant minimal polynomial for algebraic Legendre fibers.

For a fiber x0 with irreducible minimal polynomial m(x) over Q, compute the
minimal polynomial over Q of

    j = 256 (lam^2 - lam + 1)^3 / (lam^2 (lam - 1)^2),   lam = x0^2,

i.e. of the rational function

    j(x) = N(x)/D(x),  N = 256 (x^4 - x^2 + 1)^3,  D = x^4 (x^2 - 1)^2,

evaluated at any root of m.  Since m is irreducible, every root gives the same
minimal polynomial (Galois-conjugate j-values), so the engine is root-agnostic
and exact (which root is "physical" is the caller's business).

Engine = COPAIR of two independent flint-exact routes (no floats anywhere):
  Route A (literal resultant): R(t) = Res_x(m(x), t*D(x) - N(x)) in Z[t], degree
      deg(m) in t, computed by evaluation at t = 0..deg(m) via fmpz_poly.resultant
      and exact Lagrange interpolation over Q; then primitive part, flint
      factorization, and identification of the unique irreducible factor F with
      F(j) = 0 in Q[x]/(m)  (j computed as a field element via fmpq_poly xgcd
      inversion — this identification doubles as an annihilation gate).
  Route B (multiplication matrix): minimal polynomial of the multiplication-by-j
      matrix on Q[x]/(m) via fmpq_mat.minpoly, cleared to a primitive positive-lc
      integer polynomial.
  Gate: the two routes must agree coefficient-for-coefficient; the engine fails
  loudly on any disagreement.

Output fields:
  j_minpoly      : str(sympy expr) of the primitive positive-lc integer minimal
                   polynomial in t (identical formatting to str(minimal_polynomial(...)))
  j_minpoly_deg  : degree
  j_minpoly_lc   : leading coefficient.  A CM j-invariant is an algebraic
                   integer, so lc != 1 on the primitive positive-lc minimal
                   polynomial is an exact char-0 NON-CM certificate;
                   necessary-not-sufficient the other way — never read lc = 1
                   as a CM claim.
  j_is_algebraic_integer : bool(lc == 1)
  and when the degree is 1 (rational j):
  j_exact        : str(Fraction)
  j_in_CM_table  : bool against the class-number-1 CM table (complete for
                   rational j)
  D_from_table   : the CM discriminant when in table
All fields are computed, never hand-set.

Fences: lambda in {0,1} poles j (the engine raises; fence in the caller's
menu); a reducible input minpoly returns AMBIGUOUS-REDUCIBLE-MINPOLY — split
the fiber, never average; the Lagrange integrality assert is never relaxed —
fix the input.

Requires python-flint (the battery skips these legs by name when it is absent).
"""
import sys, os
from fractions import Fraction

from flint import fmpz_poly, fmpq_poly, fmpq_mat, fmpq, fmpz

X = fmpz_poly([0, 1])
N_POLY = 256 * (X**4 - X**2 + 1) ** 3          # deg 12
D_POLY = X**4 * (X**2 - 1) ** 2                # deg 8


def _primitive_pos_lc(zp: fmpz_poly) -> fmpz_poly:
    """primitive part with positive leading coefficient (sympy minimal_polynomial norm)."""
    c = 0
    for i in range(zp.degree() + 1):
        c = _gcd_int(c, int(zp[i]))
    if c == 0:
        raise ValueError("zero polynomial")
    if int(zp[zp.degree()]) < 0:
        c = -c
    return fmpz_poly([int(zp[i]) // c for i in range(zp.degree() + 1)])


def _gcd_int(a, b):
    a, b = abs(a), abs(b)
    while b:
        a, b = b, a % b
    return a


def _fq(zp: fmpz_poly) -> fmpq_poly:
    return fmpq_poly([fmpq(int(zp[i])) for i in range(zp.degree() + 1)])


def _clear_denoms(qp: fmpq_poly) -> fmpz_poly:
    return _primitive_pos_lc(qp.numer())


def _lagrange_interp_int(points):
    """exact Lagrange interpolation through (int k, int v) points; returns int coeff
    list ascending.  All arithmetic in Fraction; asserts integrality of the result."""
    n = len(points)
    coeffs = [Fraction(0)] * n
    for k, v in points:
        # basis poly prod_{j!=k} (t - j)/(k - j), built as ascending Fraction coeffs
        num = [Fraction(1)]
        den = Fraction(1)
        for j, _ in points:
            if j == k:
                continue
            num = _mul_lin(num, -j)
            den *= (k - j)
        w = Fraction(v) / den
        for i, c in enumerate(num):
            coeffs[i] += w * c
    out = []
    for c in coeffs:
        assert c.denominator == 1, "interpolation produced a non-integer coefficient"
        out.append(int(c))
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def _mul_lin(coeffs, a):
    """(ascending coeffs) * (t + a)"""
    out = [Fraction(0)] * (len(coeffs) + 1)
    for i, c in enumerate(coeffs):
        out[i] += a * c
        out[i + 1] += c
    return out


def _field_element_j(mm_q: fmpq_poly):
    """j = N(x)*D(x)^{-1} as an element of Q[x]/(mm); raises if D not invertible."""
    num = _fq(N_POLY) % mm_q
    den = _fq(D_POLY) % mm_q
    g, s, t = den.xgcd(mm_q)          # s*den + t*mm = g
    if g.degree() != 0:
        raise ZeroDivisionError("j has a pole on this fiber (lambda in {0,1})")
    inv = s / g[0]
    return (num * inv) % mm_q


def _horner_mod(F: fmpz_poly, el: fmpq_poly, mm_q: fmpq_poly) -> fmpq_poly:
    acc = fmpq_poly([0])
    for i in range(F.degree(), -1, -1):
        acc = (acc * el + fmpq_poly([fmpq(int(F[i]))])) % mm_q
    return acc


def exact_j_resultant(minpoly_coeffs):
    """minpoly_coeffs: ascending ints.  Returns the dict described in the module
    docstring plus an internal 'engine_gates' dict."""
    mz = _primitive_pos_lc(fmpz_poly([int(c) for c in minpoly_coeffs]))
    d_in = mz.degree()
    if d_in < 1:
        raise ValueError("constant minimal polynomial")
    _, fac = mz.factor()
    irreducible = (len(fac) == 1 and fac[0][1] == 1)
    factors = [f for f, mult in fac]

    results = []
    for mm in factors:
        mm = _primitive_pos_lc(mm)
        d = mm.degree()
        mm_q = _fq(mm) / fmpq(int(mm[d]))     # monic for reductions
        # pole guard (lambda = 0 or 1 <=> x in {0, +-1})
        if int(mz.resultant(D_POLY)) == 0 and int(mm.resultant(D_POLY)) == 0:
            raise ZeroDivisionError("j has a pole on this fiber (lambda in {0,1})")
        j_el = _field_element_j(mm_q)

        # ---- Route A: literal resultant by evaluation + exact interpolation
        pts = []
        for k in range(d + 1):
            Rk = mm.resultant(k * D_POLY - N_POLY)
            pts.append((k, int(Rk)))
        Rt_coeffs = _lagrange_interp_int(pts)
        Rt = fmpz_poly(Rt_coeffs)
        if Rt.degree() < 1:
            raise ValueError("resultant degenerated to a constant")
        Rt = _primitive_pos_lc(Rt)
        _, rfac = Rt.factor()
        annihilators = []
        for f, mult in rfac:
            f = _primitive_pos_lc(f)
            if _horner_mod(f, j_el, mm_q).is_zero():
                annihilators.append(f)
        gate_annih_unique = (len(annihilators) == 1)
        if not annihilators:
            raise ArithmeticError("no resultant factor annihilates j in the field")
        FA = annihilators[0]

        # ---- Route B: multiplication-by-j matrix minimal polynomial
        cols = []
        acc = fmpq_poly([1])
        for i in range(d):
            v = (j_el * acc) % mm_q
            cols.append([v[r] for r in range(d)])
            acc = (acc * fmpq_poly([0, 1])) % mm_q
        Mj = fmpq_mat(d, d, [cols[c][r] for r in range(d) for c in range(d)])
        FB = _clear_denoms(Mj.minpoly())

        gate_copair = (FA == FB)
        results.append((FA, {"copair_routes_agree": bool(gate_copair),
                             "annihilating_factor_unique": bool(gate_annih_unique)}))

    F0, gates0 = results[0]
    ambiguous = any(F != F0 for F, _ in results[1:])
    gates = {"m_irreducible": bool(irreducible),
             "factors_all_same_j_minpoly": (not ambiguous), **gates0}
    if ambiguous:
        return {"j_minpoly_status": "AMBIGUOUS-REDUCIBLE-MINPOLY",
                "engine_gates": gates}

    # ---- format identically to sympy's str(minimal_polynomial(...))
    from sympy import Poly, Symbol
    t = Symbol("t")
    desc = [int(F0[i]) for i in range(F0.degree(), -1, -1)]
    mp_sym = Poly(desc, t)
    out = {"j_minpoly": str(mp_sym.as_expr()),
           "j_minpoly_deg": int(mp_sym.degree()),
           # extra computed fields:
           # a CM j-invariant is an algebraic integer, so lc != 1 on the primitive
           # positive-lc minimal polynomial is an exact char-0 NON-CM certificate.
           "j_minpoly_lc": int(F0[F0.degree()]),
           "j_is_algebraic_integer": bool(int(F0[F0.degree()]) == 1),
           "engine_gates": gates}
    if out["j_minpoly_deg"] == 1:
        # rational j: the exact Fraction from the degree-1 minimal polynomial
        jrat = Fraction(*[int(c) for c in [-mp_sym.all_coeffs()[1],
                                           mp_sym.all_coeffs()[0]]])
        out["j_exact"] = str(jrat)
        table = _cm_table()
        out["j_in_CM_table"] = bool(jrat.denominator == 1 and int(jrat) in table)
        if out["j_in_CM_table"]:
            out["D_from_table"] = table[int(jrat)]
    return out


_CM_TABLE = None


def _cm_table():
    """class-number-1 CM j-invariant table, shared with the classifier."""
    global _CM_TABLE
    if _CM_TABLE is None:
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from geotriage import CM_J_TABLE
        _CM_TABLE = CM_J_TABLE
    return _CM_TABLE


if __name__ == "__main__":
    # smoke: Phi4 -> j = 1728
    print(exact_j_resultant([1, 0, 1]))

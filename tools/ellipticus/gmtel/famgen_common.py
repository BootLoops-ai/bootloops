#!/usr/bin/env python3
"""famgen_common — shared exact utilities for the gmtel member.

Conventions (matching the packaged L5_theta.txt reference operator):
  moments p_{2m}(A,B,C) = C(2m,m)/4^m * sum_{i+j+k=m} [m!/(i!j!k!)]^2 A^i B^j C^k
       = <lambda^{2m}>, lambda = a cos k1 + b cos k2 + c cos k3, (A,B,C)=(a^2,b^2,c^2)
  f(x) = w*W = sum_m p_{2m} x^m, x = 1/w^2;  L5 = sum_{p=0}^6 x^p P_p(theta_x)
  annihilates f;  coefficient of x^n:  sum_p P_p(n-p) p_{2(n-p)} = 0 for all n.

Operator helper: an operator in canonical D-form is a dict {j: Poly-in-var} meaning
  sum_j  c_j(v) * d^j/dv^j   (polynomial coefficients, left of the derivations).
All exact (sympy Rational / Fraction).
"""
from fractions import Fraction as F
from math import comb, factorial
import sympy as sp

# ---------------------------------------------------------------- moments
def moments(A, B, C, mmax):
    """Exact LGF moments p_{2m} for m=0..mmax as Fractions. A,B,C ints/Fractions."""
    A, B, C = F(A), F(B), F(C)
    out = []
    for m in range(mmax + 1):
        s = F(0)
        for i in range(m + 1):
            for j in range(m - i + 1):
                k = m - i - j
                mult = F(factorial(m), factorial(i) * factorial(j) * factorial(k))
                s += mult * mult * (A ** i) * (B ** j) * (C ** k)
        out.append(F(comb(2 * m, m), 4 ** m) * s)
    return out

# ---------------------------------------------------------------- L5 loader
import os as _os
# Reference L5 operator: the packaged fixture ships beside the package
# (fixtures/gmtel/L5_theta.txt); FAMGEN_L5 overrides.
L5_PATH = _os.environ.get('FAMGEN_L5', _os.path.join(
    _os.path.dirname(_os.path.abspath(__file__)),
    '..', 'fixtures', 'gmtel', 'L5_theta.txt'))

def load_L5(path=L5_PATH):
    """L5_theta.txt -> {p: sympy expr in (theta, A, B, C)}, p=0..6."""
    th, A_, B_, C_ = sp.symbols('theta A B C')
    P = {}
    for line in open(path):
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        head, expr = line.split(':', 1)
        P[int(head.strip()[2:])] = sp.sympify(expr.strip())
    assert sorted(P) == list(range(7)), sorted(P)
    return P

def L5_moment_rec(P, Av, Bv, Cv):
    """Moment recurrence sum_{j=0}^{6} c_j(m) a_{m+j} = 0 implied by L5 at (A,B,C).
    From sum_p P_p(n-p) a_{n-p} = 0: set j = 6-p, m = n-6 -> c_j(m) = P_{6-j}(m+j).
    Returns {j: sympy Poly in m} (may share overall content; caller normalizes)."""
    th, A_, B_, C_ = sp.symbols('theta A B C')
    m = sp.symbols('m')
    out = {}
    for j in range(7):
        e = P[6 - j].subs({A_: Av, B_: Bv, C_: Cv, th: m + j})
        out[j] = sp.Poly(sp.expand(e), m)
    return out

# ---------------------------------------------------------------- Ore helper (1 var)
def op_mul_poly(op, poly, v):
    """Left-multiply operator by polynomial: poly(v) * op."""
    return {j: sp.expand(poly * c) for j, c in op.items()}

def op_mul_D(op, v):
    """Left-multiply by d/dv: D * (c_j D^j) = c_j' D^j + c_j D^{j+1}."""
    out = {}
    for j, c in op.items():
        out[j] = sp.expand(out.get(j, 0) + sp.diff(c, v))
        out[j + 1] = sp.expand(out.get(j + 1, 0) + c)
    return {j: c for j, c in out.items() if sp.expand(c) != 0}

def op_add(o1, o2):
    out = dict(o1)
    for j, c in o2.items():
        out[j] = sp.expand(out.get(j, 0) + c)
    return {j: c for j, c in out.items() if sp.expand(c) != 0}

def theta_power_op(k, v):
    """theta^k as canonical D-form dict, theta = v*d/dv."""
    op = {0: sp.Integer(1)}
    for _ in range(k):
        op = op_mul_poly(op_mul_D(op, v), v, v)   # theta o op = v * D * op
    return op

def theta_form_to_D(P_of_theta, v, theta_sym):
    """Given sympy expr polynomial in theta_sym with coeffs in v (and rationals),
    return canonical D-form dict."""
    poly = sp.Poly(sp.expand(P_of_theta), theta_sym)
    out = {}
    for k, coeff in enumerate(reversed(poly.all_coeffs())):
        if coeff == 0:
            continue
        tk = theta_power_op(k, v)
        out = op_add(out, op_mul_poly(tk, coeff, v))
    return out

def op_invert_var(op, v):
    """Exact change of variable v -> 1/v on a canonical D-form operator.
    If y(v) satisfies sum c_j(v) y^{(j)}(v) = 0 then Y(u) := y(1/u) satisfies
    the transformed operator in u.  d/dv = -u^2 d/du under v = 1/u.
    Implementation: convert to theta form first (theta_v = v d/dv = -theta_u),
    i.e. write op = sum_j c_j(v) v^{-j} * theta_v(theta_v-1)...(theta_v-j+1),
    then substitute v -> 1/u, theta_v -> -theta_u, and clear denominators.
    Returns canonical D-form dict in u (same symbol object v reused as u)."""
    u = v
    th = sp.Symbol('__theta_tmp__')
    # theta-form expression: op = sum_j c_j(v) * v^{-j} * falling(th, j)
    expr = 0
    for j, c in op.items():
        fall = sp.prod([th - i for i in range(j)]) if j > 0 else sp.Integer(1)
        expr += c * v ** (-j) * fall
    # substitute v -> 1/u, theta -> -theta
    expr = sp.together(sp.expand(expr.subs({th: -th}).subs({v: 1 / u})))
    num, den = sp.fraction(sp.cancel(expr))
    # den is a power of u (up to rational): the operator num (poly in u, th) is
    # the transformed operator up to the unit den — drop den (left unit).
    poly = sp.Poly(sp.expand(num), th)
    out = {}
    for k, coeff in enumerate(reversed(poly.all_coeffs())):
        if coeff == 0:
            continue
        tk = theta_power_op(k, u)
        out = op_add(out, op_mul_poly(tk, coeff, u))
    # clear u-power content and rational content -> caller normalizes fully
    return out

def op_primitive(op, v):
    """Normalize: clear denominators, remove integer + polynomial content,
    fix sign so the leading coefficient of the highest j has positive lead."""
    if not op:
        return op
    # common denominator
    dens = [sp.denom(sp.together(c)) for c in op.values()]
    den = sp.lcm(dens) if dens else 1
    op2 = {j: sp.expand(sp.cancel(c * den)) for j, c in op.items()}
    polys = {j: sp.Poly(c, v) for j, c in op2.items()}
    # rational content
    cont = None
    for pp in polys.values():
        cc = pp.content()
        cont = cc if cont is None else sp.gcd(cont, cc)
    g = None
    for pp in polys.values():
        g = pp if g is None else g.gcd(pp)
    out = {}
    for j, pp in polys.items():
        q = sp.expand(pp.as_expr() / cont)
        if g is not None and g.degree() > 0:
            q = sp.expand(sp.cancel(q / g.as_expr()))
        out[j] = q
    jmax = max(out)
    lead = sp.Poly(out[jmax], v).LC()
    if lead < 0:
        out = {j: sp.expand(-c) for j, c in out.items()}
    return out

def op_equal(o1, o2, v):
    """Exact equality of primitive canonical forms."""
    a, b = op_primitive(o1, v), op_primitive(o2, v)
    if set(a) != set(b):
        return False
    return all(sp.expand(a[j] - b[j]) == 0 for j in a)

# ---------------------------------------------------------------- rec <-> op
def rec_to_op_x(rec, v):
    """Moment recurrence {j: Poly c_j(m)} (sum_j c_j(m) a_{m+j} = 0) -> the
    x-space operator sum_p x^p Q_p(theta) with p = r-j, Q_p(theta)=c_{r-p}(theta-p),
    returned in canonical D-form in v=x.  (Inverse of the L5_moment_rec map.)"""
    th = sp.Symbol('__theta_r__')
    m = sp.symbols('m')
    r = max(rec)
    expr = 0
    for j, c in rec.items():
        p = r - j
        cc = c.as_expr() if isinstance(c, sp.Poly) else c
        expr += v ** p * cc.subs(m, th - j)     # Q_p(theta) = c_{r-p}(theta-(r-p))
    return theta_form_to_D(sp.expand(expr), v, th)

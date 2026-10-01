"""Certified ball-arithmetic + exact product-tree engine for the Etienne formula.

Two modes, both overflow-impossible by construction:
  - exact: K(D,A) integer numerators via fmpz_poly balanced product tree
           (quasi-linear #ops; exact rational likelihood at rational (theta,m))
  - ball:  K(D,A) via arb_poly product tree at chosen precision; positivity
           (all terms >= 0) makes enclosures provably tight — the certified
           fast path for large J and for ML fitting.

Species polynomial: v_n(z) = sum_a s(n,a)(a-1)! z^a  (integer coefficients),
built from the rising factorial z(z+1)...(z+n-1) [= sum_a s(n,a) z^a] via an
fmpz_poly product tree, then Borel-twisted by (a-1)!.
K(D,A) = [z^A] prod_i v_{n_i}(z) / prod_i (n_i-1)!.
"""
from math import factorial, lgamma, log
from fractions import Fraction
from collections import Counter

from flint import fmpz_poly, arb_poly, arb, ctx


# ---------- exact integer layer ----------

_rf_cache = {}


def rising_factorial_poly(n):
    """fmpz_poly z(z+1)...(z+n-1) via balanced product tree."""
    if n in _rf_cache:
        return _rf_cache[n]
    factors = [fmpz_poly([k, 1]) for k in range(n)]  # (z + k)
    while len(factors) > 1:
        nxt = []
        for i in range(0, len(factors) - 1, 2):
            nxt.append(factors[i] * factors[i + 1])
        if len(factors) % 2:
            nxt.append(factors[-1])
        factors = nxt
    _rf_cache[n] = factors[0]
    return factors[0]


_sp_cache = {}


def species_poly_int(n):
    """fmpz_poly v_n(z) = sum_a s(n,a)(a-1)! z^a (integer Borel twist)."""
    if n in _sp_cache:
        return _sp_cache[n]
    rf = rising_factorial_poly(n)
    coeffs = [0] + [int(rf[a]) * factorial(a - 1) for a in range(1, n + 1)]
    p = fmpz_poly(coeffs)
    _sp_cache[n] = p
    return p


def _product_tree(polys):
    while len(polys) > 1:
        nxt = []
        for i in range(0, len(polys) - 1, 2):
            nxt.append(polys[i] * polys[i + 1])
        if len(polys) % 2:
            nxt.append(polys[-1])
        polys = nxt
    return polys[0]


def K_exact(D):
    """(S, [K(D,A) for A=S..J] as Fractions) — exact, product-tree."""
    D = list(D)
    S = len(D)
    prod = _product_tree([species_poly_int(n) for n in D])
    denom = 1
    for n in D:
        denom *= factorial(n - 1)
    # coefficient of z^A for A=S..J (lower coeffs are exactly 0)
    return S, [Fraction(int(prod[A]), denom) for A in range(S, sum(D) + 1)]


def L_exact(D, theta, m=None, I=None):
    """Exact rational L(theta,I) = sum_A K_A I^A/(theta)_A."""
    D = list(D)
    J, S = sum(D), len(D)
    theta = Fraction(theta)
    if I is None:
        m = Fraction(m)
        I = m * (J - 1) / (1 - m)
    else:
        I = Fraction(I)
    _, K = K_exact(D)
    tot = Fraction(0)
    t = I ** S
    for k in range(S):
        t /= theta + k
    for idx, KA in enumerate(K):
        A = S + idx
        tot += KA * t
        t *= I / (theta + A)
    return tot


def P_exact(D, theta, m=None, I=None):
    """Exact rational P[D] via the product-tree engine (fast brute force)."""
    D = sorted(D)
    J, S = sum(D), len(D)
    theta = Fraction(theta)
    if I is None:
        m = Fraction(m)
        if m == 1:
            from etienne_oracle import ewens_P
            return ewens_P(D, theta)
        I = m * (J - 1) / (1 - m)
    else:
        I = Fraction(I)
    if I == 0:
        return Fraction(1) if S == 1 else Fraction(0)
    pref = Fraction(factorial(J))
    for n in D:
        pref /= n
    for cnt in Counter(D).values():
        pref /= factorial(cnt)
    IJ = Fraction(1)
    for k in range(J):
        IJ *= I + k
    return pref * theta ** S / IJ * L_exact(D, theta, I=I)


# ---------- certified ball layer ----------

def species_poly_ball(n, prec, big_n_threshold=3000):
    """arb_poly of c_n(a) = s(n,a)(a-1)!/(n-1)!, certified at prec bits.

    For n > big_n_threshold the exact fmpz rising-factorial poly is huge
    (coefficients ~n log n bits; n=4e4 needs ~3+ GB), so build directly in
    ball arithmetic: RF poly by arb product tree of (z+k), then multiply
    coefficient a by (a-1)!/(n-1)! = 1/prod_{k=a}^{n-1} k via a descending
    running product. All positive => certified tight.
    """
    ctx.prec = prec
    if n <= big_n_threshold:
        vi = species_poly_int(n)
        d = factorial(n - 1)
        return arb_poly([arb(int(vi[a])) / arb(d) for a in range(n + 1)])
    factors = [arb_poly([arb(k), arb(1)]) for k in range(n)]
    while len(factors) > 1:
        factors = [factors[i] * factors[i + 1]
                   for i in range(0, len(factors) - 1, 2)] \
                  + ([factors[-1]] if len(factors) % 2 else [])
    rf = factors[0]
    # descending running product R_a = prod_{k=a}^{n-1} k  (R_n = 1)
    R = [arb(1)] * (n + 1)
    for a in range(n - 1, 0, -1):
        R[a] = R[a + 1] * a
    return arb_poly([arb(0)] + [rf[a] / R[a] for a in range(1, n + 1)])


def K_ball(D, prec):
    """(S, arb_poly with K(D,A) as coeff of z^A) — certified enclosures."""
    ctx.prec = prec
    D = list(D)
    cnt = Counter(D)
    polys = []
    for n, mult in cnt.items():
        p = species_poly_ball(n, prec)
        polys.extend([p] * mult)  # product tree handles balance
    return len(D), _product_tree(polys)


def logP_ball(D, theta, m, prec, Kpoly=None):
    """Certified log P[D] (arb). theta, m floats/Fractions/strings.

    Kpoly: pass the (S, arb_poly) from K_ball to reuse across (theta, m).
    Returns arb ball containing the true log-likelihood.
    """
    ctx.prec = prec
    D = sorted(D)
    J, S = sum(D), len(D)
    th = arb(str(theta)) if not isinstance(theta, arb) else theta
    mm = arb(str(m)) if not isinstance(m, arb) else m
    I = mm * (J - 1) / (1 - mm)
    if Kpoly is None:
        Kpoly = K_ball(D, prec)
    S0, KP = Kpoly
    assert S0 == S
    # A-sum with running term t_A = I^A/(theta)_A  (all positive)
    tot = arb(0)
    t = I ** S
    for k in range(S):
        t /= th + k
    for A in range(S, J + 1):
        tot += KP[A] * t
        t *= I / (th + A)
    # log prefactor + log(theta^S/(I)_J) + log(sum)
    lp = arb(J + 1).lgamma()
    for n in D:
        lp -= arb(n).log()
    for c in Counter(D).values():
        lp -= arb(c + 1).lgamma()
    lp += S * th.log()
    lp -= (I + J).lgamma() - I.lgamma()
    lp += tot.log()
    return lp

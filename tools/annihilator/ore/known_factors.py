"""
Exact-arithmetic differential-operator factors for the diagonal Ising
susceptibility χ̃_d^(3) and χ̃_d^(5).

Operator coefficients transcribed from the published equations in:
  * arXiv:1110.1705, Appendix C (label appL11) — U^(5)_1, V^(5)_1, W^(5)_1,
    L^(5)_4, and eq.(72) for the L^(5)_11 direct-sum structure.
  * arXiv:1110.1705 §2 eq.(l33) and the companion χ̃_d^(3) literature —
    L^(3)_1, L^(3)_2, L^(3)_3.

CAVEAT: L5_4 below is the PRINTED arXiv:1110.1705 Appendix C operator,
transcribed as published; the printed L^(5)_4 has been flagged as the typo
source in the L4·V1·U1 direct summand.  The LCLM/right-division algebra
exercised by the selftests here is unaffected (right-divisibility of an LCLM
by each summand holds for whatever operators are supplied), but do not treat
L5_4 as physics-validated.

Conventions
-----------
Variable: x  (x = t^{1/2}, high-temperature).
D_x = d/dx.

Each operator L of order n is stored as a Python list

    L = [c0, c1, ..., cn]

meaning  L = c0(x) + c1(x)·D_x + ... + cn(x)·D_x^n,  with each c_j a
sympy.Poly in x with EXACT INTEGER coefficients (a common polynomial
denominator has been cleared; the operator is defined only up to an overall
rational-function left-multiple, which is irrelevant for its kernel).
"""

from fractions import Fraction
import sympy as sp

x = sp.symbols('x')
Dx = sp.Function('Dx')  # placeholder symbol only; operators are coefficient lists


def _P(expr):
    """Shorthand: expand expr and return as sympy.Poly in x over ZZ."""
    return sp.Poly(sp.expand(expr), x, domain='ZZ')


# ---------------------------------------------------------------------------
#  χ̃_d^(3) factors:  L^(3)_1 (order 1),  L^(3)_2 (order 2),  L^(3)_3 (order 3)
#  Source: the companion χ̃_d^(3) literature;  1110.1705 eq.(l33) for L^(3)_3.
# ---------------------------------------------------------------------------

# L^(3)_1 = D_x + 1/(x-1)           (annihilates 1/(1-x) = χ̃_d^(1))
# cleared: (x-1)·D_x + 1
L3_1 = [_P(1), _P(x - 1)]

# L^(3)_2 = D_x^2 + 2(1+2x)/((1+x)(x-1)) · D_x + (1+2x)/((1+x)(x-1)x)
# cleared by x(1+x)(x-1):
L3_2 = [
    _P(1 + 2*x),
    _P(2*x*(1 + 2*x)),
    _P(x*(1 + x)*(x - 1)),
]

# L^(3)_3 = D_x^3 + (3/2)·n2/d · D_x^2 + n1/((x+1)(x-1)·x·d) · D_x
#                 + n0/((x+1)(x-1)^2·x·d)
# with  d(x) = (x+2)(1+2x)(x+1)(x-1)(1+x+x^2)·x
# cleared by  2·(x+2)(1+2x)(x+1)^2(x-1)^3(1+x+x^2)·x^2 :
_n0 = 2*x**8 + 8*x**7 - 7*x**6 - 13*x**5 - 58*x**4 - 88*x**3 - 52*x**2 - 13*x + 5
_n1 = (14*x**8 + 71*x**7 + 146*x**6 + 170*x**5 + 38*x**4
       - 112*x**3 - 94*x**2 - 19*x + 2)
_n2 = 8*x**6 + 36*x**5 + 63*x**4 + 62*x**3 + 21*x**2 - 6*x - 4
_d  = (x + 2)*(1 + 2*x)*(x + 1)*(x - 1)*(1 + x + x**2)*x

L3_3 = [
    _P(2*_n0),
    _P(2*(x - 1)*_n1),
    _P(3*x*(x + 1)*(x - 1)**2 * _n2),
    _P(2*(x + 2)*(1 + 2*x)*(x + 1)**2*(x - 1)**3*(1 + x + x**2)*x**2),
]


# ---------------------------------------------------------------------------
#  χ̃_d^(5) order-1 factors  U^(5)_1, V^(5)_1, W^(5)_1
#  Source: 1110.1705 Appendix C (eqs. C.1–C.3).
#
#  Each is  D_x − (d/dx) ln R(x)   or  D_x − (1/2)(d/dx) ln R(x),
#  i.e. annihilates R(x) resp. R(x)^{1/2}.  We clear to integer-poly form
#  [c0, c1]  =  [ −R'·(denom)/R (or −R'/(2R)·denom),  denom ].
# ---------------------------------------------------------------------------

def _order1_from_dlog(R_num, R_den, half=False):
    """Return [c0, c1] for  D_x − k·(d/dx) ln(R_num/R_den),  k = 1/2 if half else 1,
    cleared to integer-polynomial coefficients."""
    R = sp.together(R_num / R_den)
    dlog = sp.cancel(sp.diff(sp.log(R), x))
    coeff = sp.Rational(1, 2) * dlog if half else dlog
    # operator = D_x − coeff  ;  write coeff = p/q, clear:
    p, q = sp.fraction(sp.cancel(coeff))
    # q·D_x − p ; ensure integer polys
    q = sp.expand(q)
    p = sp.expand(p)
    # If half introduced a 1/2 that survived in q, multiply through:
    lc = sp.lcm_list([t.q for t in sp.Poly(q, x).all_coeffs()] +
                     [t.q for t in sp.Poly(p, x).all_coeffs()])
    return [_P(-lc*p), _P(lc*q)]


# U^(5)_1 = D_x − d/dx ln( x / (1−x)^3 )
U5_1 = _order1_from_dlog(x, (1 - x)**3, half=False)

# V^(5)_1 = D_x − (1/2) d/dx ln( (1+x+x^2)^3 / ((1+x)^2 (1−x)^6 x^2) )
V5_1 = _order1_from_dlog((1 + x + x**2)**3,
                         (1 + x)**2 * (1 - x)**6 * x**2, half=True)

# W^(5)_1 = D_x − (1/2) d/dx ln( (x^2+1)^2 / ((1+x)^2 (1−x)^6 x) )
W5_1 = _order1_from_dlog((x**2 + 1)**2,
                         (1 + x)**2 * (1 - x)**6 * x, half=True)


# ---------------------------------------------------------------------------
#  L^(5)_4  (order 4)
#  Source: 1110.1705 Appendix C, eq.(C.4) with p0..p4 given explicitly.
#  L^(5)_4 = D_x^4 + (p3/p4) D_x^3 + (p2/p4) D_x^2 + (p1/p4) D_x + p0/p4
#  ⇒ cleared form  [p0, p1, p2, p3, p4].
# ---------------------------------------------------------------------------

_p4_inner = (160 + 3148*x + 24988*x**2 + 86008*x**3 + 141698*x**4
             + 69707*x**5 - 141750*x**6 - 358707*x**7 - 356606*x**8
             - 1071*x**9 + 347302*x**10 + 510214*x**11 + 347302*x**12
             - 1071*x**13 - 356606*x**14 - 358707*x**15 - 141750*x**16
             + 69707*x**17 + 141698*x**18 + 86008*x**19 + 24988*x**20
             + 3148*x**21 + 160*x**22)
p4 = _P(x**3 * (1 + x + x**2) * (x + 1)**3 * (x - 1)**4 * _p4_inner)

_p3_inner = (-880 - 16620*x - 126586*x**2 - 421558*x**3 - 520547*x**4
             + 733378*x**5 + 3794648*x**6 + 6252130*x**7 + 3922367*x**8
             - 4349032*x**9 - 12817741*x**10 - 12881692*x**11
             - 2612141*x**12 + 10986996*x**13 + 16830947*x**14
             + 12283572*x**15 + 729267*x**16 - 8919176*x**17
             - 10905121*x**18 - 5398478*x**19 + 866024*x**20
             + 3665682*x**21 + 3069821*x**22 + 1351818*x**23
             + 323590*x**24 + 36308*x**25 + 1680*x**26)
p3 = _P(2 * x**2 * (x + 1)**2 * (x - 1)**3 * _p3_inner)

_p2_inner = (2400 + 38692*x + 228422*x**2 + 366806*x**3 - 1591741*x**4
             - 8948446*x**5 - 18137183*x**6 - 10301088*x**7
             + 31576074*x**8 + 82978356*x**9 + 80098415*x**10
             - 8308172*x**11 - 123518048*x**12 - 158759046*x**13
             - 65285821*x**14 + 78248130*x**15 + 152708392*x**16
             + 124727752*x**17 + 26488355*x**18 - 65301174*x**19
             - 90679899*x**20 - 47527872*x**21 + 4032496*x**22
             + 27473954*x**23 + 23107094*x**24 + 9927812*x**25
             + 2288564*x**26 + 245416*x**27 + 10800*x**28)
p2 = _P(2 * x * (x - 1)**2 * _p2_inner)

_p1_inner = (-1440 - 15176*x - 3552*x**2 + 632252*x**3 + 3988986*x**4
             + 11012538*x**5 + 10122851*x**6 - 31358640*x**7
             - 125311964*x**8 - 166380144*x**9 + 20063039*x**10
             + 375202188*x**11 + 523233277*x**12 + 189830162*x**13
             - 422078559*x**14 - 747281488*x**15 - 440223099*x**16
             + 161161298*x**17 + 530901457*x**18 + 491902752*x**19
             + 168466049*x**20 - 168274188*x**21 - 282329480*x**22
             - 158906808*x**23 - 754525*x**24 + 72189798*x**25
             + 61435092*x**26 + 25677392*x**27 + 5672988*x**28
             + 577984*x**29 + 24000*x**30)
p1 = _P(2 * (x - 1) * _p1_inner)

p0 = _P(-3600 - 52880*x - 324108*x**2 - 1147996*x**3 - 1575180*x**4
        + 8228874*x**5 + 52977905*x**6 + 108476130*x**7 - 739178*x**8
        - 371064711*x**9 - 563202298*x**10 - 29824206*x**11
        + 842725375*x**12 + 1075242362*x**13 + 273493047*x**14
        - 909934423*x**15 - 1189246308*x**16 - 414338515*x**17
        + 420114304*x**18 + 702981552*x**19 + 447865799*x**20
        + 30467322*x**21 - 270639170*x**22 - 233990685*x**23
        - 67035676*x**24 + 45089100*x**25 + 61580064*x**26
        + 29851532*x**27 + 7030080*x**28 + 714400*x**29 + 28800*x**30)

L5_4 = [p0, p1, p2, p3, p4]


# ---------------------------------------------------------------------------
#  Operator algebra on power series (exact, over Fraction).
# ---------------------------------------------------------------------------

def _poly_coeffs(c):
    """Return list [a0, a1, ...] of Fraction coeffs of a sympy Poly / expr in x."""
    p = c if isinstance(c, sp.Poly) else sp.Poly(c, x)
    d = p.degree()
    if d < 0:
        return [Fraction(0)]
    out = [Fraction(0)] * (d + 1)
    for (k,), v in p.terms():
        out[k] = Fraction(int(v))
    return out


def _series_derivative(a):
    """Given [a0,a1,...,aN] for Σ a_k x^k, return derivative series (length N)."""
    return [Fraction(k) * a[k] for k in range(1, len(a))]


def _poly_times_series(pcoeffs, s, N):
    """Multiply polynomial (list of Fraction) by series (list of Fraction),
    return first N coefficients."""
    out = [Fraction(0)] * N
    for i, pi in enumerate(pcoeffs):
        if pi == 0:
            continue
        for k in range(min(len(s), N - i)):
            out[i + k] += pi * s[k]
    return out


def apply_op(op, series_coeffs):
    """Apply operator op = [c0,...,cn] (each a sympy Poly in x) to a power
    series given as a list of Fraction coefficients [a0,a1,...,a_{N-1}].
    Returns the resulting series truncated to the same length N.

    NOTE: the top `n` output coefficients (n = order of op) are UNRELIABLE
    since they depend on input coefficients a_N, a_{N+1}, ... which are not
    supplied; callers testing for annihilation should ignore the last
    (len(op)-1) entries.
    """
    N = len(series_coeffs)
    s = [Fraction(a) for a in series_coeffs]
    # precompute derivatives f, f', f'', ... up to order n, each padded to length N
    derivs = [s]
    for _ in range(1, len(op)):
        derivs.append(_series_derivative(derivs[-1]) + [Fraction(0)])
    # pad all to length N
    derivs = [d + [Fraction(0)] * (N - len(d)) for d in derivs]
    out = [Fraction(0)] * N
    for j, cj in enumerate(op):
        pc = _poly_coeffs(cj)
        term = _poly_times_series(pc, derivs[j], N)
        for k in range(N):
            out[k] += term[k]
    return out


def compose_ops(L, R):
    """Operator product L·R (apply R first, then L).
    L = Σ_i a_i D^i,  R = Σ_j b_j D^j.  Uses Leibniz:
        a_i D^i ∘ b_j D^j = a_i · Σ_{k=0}^{i} C(i,k) b_j^{(k)} D^{i-k+j}.
    Returns coefficient list of sympy.Poly's.
    """
    m = len(L) - 1
    n = len(R) - 1
    out = [sp.Integer(0)] * (m + n + 1)
    # Express each b_j as sympy expression for easy differentiation.
    b = [c.as_expr() if isinstance(c, sp.Poly) else sp.sympify(c) for c in R]
    a = [c.as_expr() if isinstance(c, sp.Poly) else sp.sympify(c) for c in L]
    for i in range(m + 1):
        ai = a[i]
        if ai == 0:
            continue
        for j in range(n + 1):
            bj = b[j]
            # D^i ∘ (bj D^j) = Σ_{k=0}^{i} C(i,k) bj^{(k)} D^{i-k+j}
            bder = bj
            for k in range(i + 1):
                if k > 0:
                    bder = sp.diff(bder, x)
                coef = sp.binomial(i, k) * ai * bder
                out[i - k + j] += sp.expand(coef)
    return [_P(c) for c in out]


# ---------------------------------------------------------------------------
#  L^(5)_11 direct-sum structure (1110.1705 eq.(72) / label factoL11):
#
#    L^(5)_11 = L^(3)_1  ⊕  L^(3)_3  ⊕  (W^(5)_1 · U^(5)_1)
#                        ⊕  (L^(5)_4 · V^(5)_1 · U^(5)_1).
#
#  We build each direct-summand product exactly.  (A full LCLM is not
#  constructed here; each summand suffices to test annihilation of the
#  corresponding solution component.)
# ---------------------------------------------------------------------------

WU_2  = compose_ops(W5_1, U5_1)                       # order 2
VU_2  = compose_ops(V5_1, U5_1)                       # order 2
LVU_6 = compose_ops(L5_4, VU_2)                       # order 6

L11_5_summands = {
    'L3_1':  L3_1,
    'L3_3':  L3_3,
    'W1.U1': WU_2,
    'L4.V1.U1': LVU_6,
}


# ---------------------------------------------------------------------------
#  Quick self-checks when run as a script.
# ---------------------------------------------------------------------------
if __name__ == '__main__':
    ok = True
    report = []

    def annihilates(op, ser, tag):
        r = apply_op(op, ser)
        n = len(op) - 1
        good = all(c == 0 for c in r[:len(r) - n])
        report.append((tag, good))
        return good

    # 1.  L^(3)_1 must annihilate 1/(1-x) = Σ x^k.
    N = 60
    ones = [Fraction(1)] * N
    pass1 = annihilates(L3_1, ones, "L3_1 annihilates 1/(1-x) [60 terms]")
    ok &= pass1

    # 2.  U^(5)_1 must annihilate x/(1-x)^3 = Σ_{k≥1} k(k+1)/2 · x^k.
    ser_U = [Fraction(k*(k+1), 2) for k in range(N)]
    pass2 = annihilates(U5_1, ser_U, "U5_1 annihilates x/(1-x)^3 [60 terms]")
    ok &= pass2

    # 2b.  W·U and L4·V·U must both annihilate x/(1-x)^3 (right-factor U5_1).
    pass2b = annihilates(WU_2, ser_U, "W1.U1 annihilates x/(1-x)^3")
    ok &= pass2b
    pass2c = annihilates(LVU_6, ser_U, "L4.V1.U1 annihilates x/(1-x)^3")
    ok &= pass2c

    # 2d.  L^(3)_2 must annihilate  χ̃^(3)_{d;2} = E(x²)/(1-x)² − K(x²)/(1-x)
    #      (1110.1705 eq.(chi2) and the companion χ̃_d^(3) literature).
    E = sp.hyperexpand(sp.hyper([sp.Rational(1,2), -sp.Rational(1,2)], [1], x**2))
    K = sp.hyperexpand(sp.hyper([sp.Rational(1,2),  sp.Rational(1,2)], [1], x**2))
    chi2_expr = E/(1-x)**2 - K/(1-x)
    chi2_ser = sp.series(chi2_expr, x, 0, N).removeO()
    chi2_coeffs = [Fraction(sp.Rational(chi2_ser.coeff(x, k))) for k in range(N)]
    pass2d = annihilates(L3_2, chi2_coeffs,
                         "L3_2 annihilates E(x^2)/(1-x)^2 - K(x^2)/(1-x)")
    ok &= pass2d

    # 3.  Palindrome check on p4 inner degree-22 factor (paper states palindromic).
    inner = sp.Poly(_p4_inner, x).all_coeffs()
    pass3 = inner == inner[::-1]
    report.append(("p4 inner deg-22 factor is palindromic", pass3))
    ok &= pass3

    # 4.  Degrees of p0..p4 as expected from Appendix C.
    degs = [c.degree() for c in L5_4]
    pass4 = degs == [30, 31, 31, 33, 34]
    report.append((f"deg(p0..p4) = {degs} == [30,31,31,33,34]", pass4))
    ok &= pass4

    # 5.  Orders of composed summands.
    pass5 = (len(WU_2) - 1 == 2) and (len(LVU_6) - 1 == 6)
    report.append(("orders: W·U=2, L4·V·U=6", pass5))
    ok &= pass5

    for msg, p in report:
        print(("PASS  " if p else "FAIL  ") + msg)
    print("\nOVERALL:", "PASS" if ok else "FAIL")

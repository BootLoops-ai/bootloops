#!/usr/bin/env python3
"""Exact nonnegativity of sparse rational polynomials on [0,1] — two routes.

Route S (Sturm-class, complete): sympy squarefree decomposition + exact real
  root counting of the odd-multiplicity part in the open interval.
Route B (Bernstein/de Casteljau, independent): exact rational Bernstein
  expansion; touch-root deflation by synthetic division; dyadic subdivision.
Both routes are pure rational arithmetic end-to-end.
"""
from gmpy2 import mpq, lcm
from math import comb
import flint
import sympy as sp

_t = sp.symbols('t')

def _ball_to_interval(x):
    mm, me = x.mid().man_exp()
    rm, re = x.rad().man_exp()
    mid = mpq(int(mm)) * mpq(2) ** int(me) if int(mm) else mpq(0)
    rad = mpq(int(rm)) * mpq(2) ** int(re) if int(rm) else mpq(0)
    return mid - rad, mid + rad

def _to_fmpz(p):
    den = 1
    for v in p.values():
        den = lcm(den, v.denominator)
    dmax = max(p)
    dense = [int(p.get(e, mpq(0)) * den) for e in range(dmax + 1)]
    return flint.fmpz_poly(dense)

def nonneg_route_R(coeffs, max_tries=5):
    """Complete exact nonnegativity on [0,1] (Sturm-class replacement; the
    sympy sqf/count_roots route measured >6 min per functional and was
    retired). Method: certified root isolation (flint complex_roots of the
    squarefree part, arb balls) + EXACT rational sign sampling of the original
    polynomial in the gaps between root regions. Sign logic: p >= 0 on [0,1]
    iff p(0),p(1) >= 0 and p > 0 at one exact rational point in every gap
    between consecutive distinct-root regions (single-root regions cannot hide
    a sign excursion; multi-root regions force a precision escalation).
    Returns (True/False/None, info, dips)."""
    p = {e: mpq(c) for e, c in coeffs.items() if c != 0}
    if not p:
        return True, 'zero-poly', []
    v0, v1 = poly_eval(p, mpq(0)), poly_eval(p, mpq(1))
    if v0 < 0 or v1 < 0:
        return False, 'endpoint negative', [t for t, v in ((mpq(0), v0), (mpq(1), v1)) if v < 0]
    Pz = _to_fmpz(p)
    sf = Pz // Pz.gcd(Pz.derivative())
    old = flint.ctx.prec
    try:
        for attempt in range(max_tries):
            flint.ctx.prec = 128 << attempt
            try:
                cr = sf.complex_roots()
            except Exception:
                continue
            iv = []
            for r, _m in cr:
                if not r.imag.contains(flint.arb(0)):
                    continue
                lo, hi = _ball_to_interval(r.real)
                if hi <= 0 or lo >= 1:
                    continue
                iv.append((max(lo, mpq(0)), min(hi, mpq(1)), 1))
            iv.sort()
            merged = []
            for a, b, c in iv:
                if merged and a <= merged[-1][1]:
                    pa, pb, pc = merged[-1]
                    merged[-1] = (pa, max(pb, b), pc + c)
                else:
                    merged.append((a, b, c))
            if any(c > 1 for _a, _b, c in merged):
                continue  # roots not separated at this precision
            pts = []
            prev = mpq(0)
            for a, b, _c in merged:
                if a > prev:
                    pts.append((prev + a) / 2)
                prev = max(prev, b)
            if prev < 1:
                pts.append((prev + 1) / 2)
            dips = []
            for s in pts:
                v = poly_eval(p, s)
                k = 1
                while v == 0 and k < 50:
                    s = s + mpq(k, 10 ** 12)
                    v = poly_eval(p, s)
                    k += 1
                if v < 0:
                    dips.append(s)
            if dips:
                return False, f'sign change ({len(dips)} dip pts, {len(merged)} root regions)', dips
            return True, f'{len(merged)} isolated root regions, all gaps positive', []
        return None, 'root separation failed at max precision', []
    finally:
        flint.ctx.prec = old

def poly_eval(coeffs, x):
    """coeffs: dict{exponent: mpq}; x: mpq."""
    x = mpq(x)
    return sum((c * x ** e for e, c in coeffs.items()), mpq(0))

def _to_sympy(coeffs):
    return sp.Poly({int(e): sp.Rational(int(c.numerator), int(c.denominator))
                    for e, c in coeffs.items() if c != 0}, _t, domain='QQ')

def nonneg_route_S(coeffs):
    """Returns (verdict_bool_or_None, info)."""
    P = _to_sympy(coeffs)
    if P.is_zero:
        return True, 'zero-poly'
    if P.eval(0) < 0 or P.eval(1) < 0:
        return False, 'endpoint negative'
    _, factors = P.sqf_list()
    odd = sp.Poly(1, _t, domain='QQ')
    for f, e in factors:
        if e % 2 == 1:
            odd = odd * f
    if odd.degree() == 0:
        # perfect even power times constant: sign = sign of leading const anywhere
        v = P.eval(sp.Rational(1, 2))
        if v > 0:
            return True, 'even-square structure'
        if v == 0:
            # find a non-root rational point
            for dq in range(3, 200):
                v = P.eval(sp.Rational(1, dq))
                if v != 0:
                    return (v > 0), f'sign at 1/{dq}'
        return (v > 0), 'sign probe'
    n_closed = odd.count_roots(0, 1)
    n_end = int(odd.eval(0) == 0) + int(odd.eval(1) == 0)
    n_open = n_closed - n_end
    if n_open == 0:
        # no interior sign change; interior sign decides
        for dq in [2, 3, 5, 7, 11, 13, 17, 19, 23]:
            v = P.eval(sp.Rational(1, dq))
            if v != 0:
                return (v > 0), 'no interior odd roots'
        return None, 'could not find non-root probe'
    return False, f'{n_open} interior odd-multiplicity roots'

def find_dips(coeffs, extra_points=()):
    """Rational t in (0,1) with p(t) < 0: fine dyadic grid + near-endpoint
    ladders + sympy isolating intervals of the squarefree part."""
    dips = []
    grid = [mpq(j, 256) for j in range(1, 256)]
    grid += [1 - mpq(1, 1 << k) for k in range(9, 26)]
    grid += [mpq(1, 1 << k) for k in range(9, 26)]
    grid += [mpq(x) for x in extra_points]
    for x in grid:
        if poly_eval(coeffs, x) < 0:
            dips.append(x)
    return sorted(set(dips))

def simplify_dip(coeffs, s):
    """Replace a huge-denominator dip point by a coarse dyadic that still
    certifies p < 0 exactly (keeps LP column bit-size bounded)."""
    s = mpq(s)
    for k in range(8, 100):
        num = (s.numerator * (1 << k)) // s.denominator
        for cand in (mpq(num, 1 << k), mpq(num + 1, 1 << k)):
            if 0 < cand < 1 and poly_eval(coeffs, cand) < 0:
                return cand
    return s

# ------------------------- Route B: Bernstein ------------------------------

def _deflate(coeffs, r):
    """Divide dict-poly by (t - r) exactly; returns (quotient dict, remainder)."""
    d = max(coeffs)
    dense = [coeffs.get(e, mpq(0)) for e in range(d + 1)]
    q = [mpq(0)] * d
    carry = mpq(0)
    for e in range(d, 0, -1):
        c = dense[e] + carry
        q[e - 1] = c
        carry = c * r
    rem = dense[0] + carry
    return {e: q[e] for e in range(d) if q[e] != 0}, rem

def _bernstein_coeffs(coeffs, d):
    dense = [coeffs.get(e, mpq(0)) for e in range(d + 1)]
    return [sum((dense[j] * mpq(comb(k, j), comb(d, j)) for j in range(0, k + 1)), mpq(0))
            for k in range(d + 1)]

def _decasteljau_split(b):
    d = len(b) - 1
    left, right = [b[0]], [b[-1]]
    cur = list(b)
    half = mpq(1, 2)
    for _ in range(d):
        cur = [(cur[i] + cur[i + 1]) * half for i in range(len(cur) - 1)]
        left.append(cur[0])
        right.append(cur[-1])
    right.reverse()
    return left, right

def nonneg_route_B(coeffs, touch_points=(), depth_cap=12, node_cap=4000):
    """touch_points: rational t in (0,1) where p is known/suspected to vanish
    (LP-binding curve samples). Deflates even-multiplicity rational roots first.
    Returns (True/False/None, info)."""
    p = {e: mpq(c) for e, c in coeffs.items() if c != 0}
    if not p:
        return True, 'zero-poly'
    # deflate interior rational touch roots
    for r in touch_points:
        r = mpq(r)
        if not (0 < r < 1):
            continue
        m = 0
        while poly_eval(p, r) == 0 and p:
            q, rem = _deflate(p, r)
            assert rem == 0
            p, m = q, m + 1
        if m % 2 == 1:
            if poly_eval(p, r) < 0:
                return False, f'odd-mult root at {r} with negative cofactor'
            # odd multiplicity with positive cofactor => sign change => negative side
            return False, f'odd-mult interior root at {r}'
    # deflate t=0 and t=1 roots (any multiplicity; both factors >=0 on [0,1])
    while p and min(p) > 0:
        p = {e - 1: c for e, c in p.items()}
    while p and poly_eval(p, mpq(1)) == 0:
        q, rem = _deflate(p, mpq(1))
        assert rem == 0
        p = {e: -c for e, c in q.items()}  # divided by (t-1); flip to (1-t) factor
    if not p:
        return True, 'monomial factors only'
    d = max(p)
    b = _bernstein_coeffs(p, d)
    stack = [b]
    nodes = 0
    while stack:
        cur = stack.pop()
        nodes += 1
        if nodes > node_cap:
            return None, 'node cap'
        if cur[0] < 0 or cur[-1] < 0:
            return False, 'negative endpoint value'
        if all(x >= 0 for x in cur):
            continue
        if nodes > (1 << depth_cap):
            return None, 'depth cap'
        L, R = _decasteljau_split(cur)
        stack.append(L)
        stack.append(R)
    return True, f'bernstein certified ({nodes} nodes)'

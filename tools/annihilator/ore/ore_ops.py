"""
ore_ops.py — minimal Ore-algebra operations over ℚ(x) for linear differential operators.

Operators are lists [c0, c1, ..., cn] of sympy rationals-in-x meaning
    L = c0 + c1·Dx + ... + cn·Dx^n.

Provides: mul (compose), apply_to_series, lclm, right_divide,
formal_solution (indicial-exponent series from a recurrence).
"""
import sympy as sp
from fractions import Fraction
from functools import reduce

x = sp.symbols('x')


def _norm(op):
    """Strip trailing zeros; ensure sympy exprs."""
    op = [sp.together(sp.sympify(c)) for c in op]
    while len(op) > 1 and op[-1] == 0:
        op.pop()
    return op


def order(op):
    return len(_norm(op)) - 1


def mul(L, R):
    """Operator product L·R (apply R first, then L). Leibniz on Dx·f = f·Dx + f'."""
    L, R = _norm(L), _norm(R)
    # Dx^i · (r_j Dx^j) = sum_k C(i,k) r_j^{(k)} Dx^{i-k+j}
    m, n = order(L), order(R)
    out = [sp.Integer(0)] * (m + n + 1)
    # precompute derivatives of R coeffs
    Rderiv = [[sp.diff(rj, x, k) for k in range(m + 1)] for rj in R]
    for i, li in enumerate(L):
        if li == 0:
            continue
        for j, rj in enumerate(R):
            if rj == 0:
                continue
            for k in range(i + 1):
                out[i - k + j] += li * sp.binomial(i, k) * Rderiv[j][k]
    return _norm([sp.together(sp.expand(c)) for c in out])


def apply_to_series(op, coeffs):
    """Apply op to a power series given as list of Fractions [a0,a1,...,aN].
    Returns list of Fractions of length N - deg_drop (safe truncation).
    op coeffs are sympy polys in x with integer coeffs (denominator-cleared)."""
    op = _norm(op)
    N = len(coeffs)
    n = order(op)
    # d^k/dx^k series: [a_m -> (m)_k a_{m}] shifted: coefficient of x^{m} in f^{(k)} is (m+k)!/m! * a_{m+k}
    # Represent op as poly in x with poly-in-Dx coeffs, act on series.
    # Easiest: expand each c_i(x) as poly, then (c_i * f^{(i)})_m = sum_j [x^j]c_i * [x^{m-j}]f^{(i)}
    out = [Fraction(0)] * N
    for i, ci in enumerate(op):
        ci_poly = sp.Poly(sp.expand(ci), x)
        ci_coeffs = ci_poly.all_coeffs()[::-1]  # [c0, c1, ...] low-to-high
        # f^{(i)}[m] = falling(m+i, i) * a_{m+i}
        for m in range(N):
            # sum over j: ci_coeffs[j] * fderiv_i[m-j]
            s = Fraction(0)
            for j, cj in enumerate(ci_coeffs):
                idx = m - j
                if idx < 0 or idx + i >= N:
                    continue
                fall = 1
                for t in range(i):
                    fall *= (idx + i - t)
                s += Fraction(int(cj)) * fall * coeffs[idx + i]
            out[m] += s
    # last (n + max_deg) entries unreliable
    maxdeg = max((sp.Poly(sp.expand(c), x).degree() for c in op if c != 0), default=0)
    return out[: max(0, N - n - maxdeg)]


def lclm2(A, B):
    """LCLM of two operators via linear algebra over ℚ(x):
    find lowest-order P,Q with P·A = Q·B; return P·A."""
    A, B = _norm(A), _norm(B)
    a, b = order(A), order(B)
    # Try target orders p+a = q+b = a+b down to max(a,b); generically LCLM has order a+b - ord(GCRD).
    for tgt in range(max(a, b), a + b + 1):
        p, q = tgt - a, tgt - b
        if p < 0 or q < 0:
            continue
        # unknowns: P coeffs (p+1) rat-funcs, Q coeffs (q+1). Equation: (P·A - Q·B) = 0 as operator of order tgt.
        # This is (tgt+1) equations over ℚ(x) in (p+q+2) unknowns. Nontrivial null iff rank deficient.
        Psyms = sp.symbols(f'P0:{p+1}', cls=sp.Dummy)
        Qsyms = sp.symbols(f'Q0:{q+1}', cls=sp.Dummy)
        PA = mul(list(Psyms), A)
        QB = mul(list(Qsyms), B)
        diff = [sp.expand((PA[k] if k < len(PA) else 0) - (QB[k] if k < len(QB) else 0))
                for k in range(tgt + 1)]
        # Solve linear system for Psyms, Qsyms over ℚ(x)
        sol = sp.solve(diff, list(Psyms) + list(Qsyms), dict=True)
        if sol:
            s = sol[0]
            # need nontrivial: check not all zero
            Pval = [s.get(Pi, Pi) for Pi in Psyms]
            # sp.solve leaves free params; substitute a generic value for any remaining symbol
            free = [v for v in (list(Psyms) + list(Qsyms)) if v not in s]
            if free:
                subs = {free[-1]: 1, **{f: 0 for f in free[:-1]}}
                Pval = [sp.simplify(pv.subs(subs)) if hasattr(pv, 'subs') else pv for pv in Pval]
            if any(pv != 0 for pv in Pval):
                L = mul(Pval, A)
                # normalize: clear denominators
                den = sp.lcm_list([sp.fraction(sp.together(c))[1] for c in L])
                L = _norm([sp.cancel(den * c) for c in L])
                return L
    raise RuntimeError("LCLM not found")


def lclm(*ops):
    return reduce(lclm2, ops)


def right_divide(N, R):
    """If N = Q·R exactly, return Q. Else raise."""
    N, R = _norm(N), _norm(R)
    n, r = order(N), order(R)
    q = n - r
    if q < 0:
        raise ValueError("order(N) < order(R)")
    Qsyms = sp.symbols(f'Q0:{q+1}', cls=sp.Dummy)
    QR = mul(list(Qsyms), R)
    eqs = [sp.expand((QR[k] if k < len(QR) else 0) - (N[k] if k < len(N) else 0))
           for k in range(n + 1)]
    sol = sp.solve(eqs, list(Qsyms), dict=True)
    if not sol:
        raise RuntimeError("R does not right-divide N")
    s = sol[0]
    Q = [sp.together(s.get(Qi, Qi)) for Qi in Qsyms]
    return _norm(Q)


def op_to_recurrence(op):
    """Convert L = Σ c_i(x) Dx^i to a recurrence on series coeffs a_m:
    Σ_k r_k(m) a_{m+k} = 0. Returns (offset_min, [r_0,...,r_K]) with r_j sympy polys in m."""
    m = sp.symbols('m', integer=True)
    # x^j Dx^i acting on Σ a_n x^n: coeff of x^m is (m-j+i)!/(m-j)! a_{m-j+i}
    # So shift index s = i - j.
    op = _norm(op)
    terms = {}  # shift s -> poly in m
    for i, ci in enumerate(op):
        for j, cj in enumerate(sp.Poly(sp.expand(ci), x).all_coeffs()[::-1]):
            if cj == 0:
                continue
            s = i - j
            fall = sp.prod([m - j + i - t for t in range(i)]) if i > 0 else sp.Integer(1)
            terms[s] = terms.get(s, sp.Integer(0)) + sp.Integer(int(cj)) * fall
    smin, smax = min(terms), max(terms)
    rec = [sp.expand(terms.get(smin + k, 0)) for k in range(smax - smin + 1)]
    return smin, rec, m


def formal_solution(op, lead_exp, N):
    """Compute the analytic formal series solution starting x^{lead_exp} to N terms,
    assuming lead_exp is a simple highest-integer indicial root (unique solution)."""
    smin, rec, m = op_to_recurrence(op)
    K = len(rec) - 1
    a = [Fraction(0)] * N
    a[lead_exp] = Fraction(1)
    # recurrence at index m: Σ_k rec[k](m) a_{m+smin+k} = 0.
    # Solve forward for a_{m+smin+K} given lower.
    # Rearranged: for each output index n ≥ lead_exp+1, use m such that m+smin+K = n, i.e. m = n - smin - K.
    for n in range(lead_exp + 1, N):
        mval = n - smin - K
        top = rec[K].subs(m, mval)
        if top == 0:
            # resonance — leave a[n] free (=0) if consistent
            # check consistency
            rhs = sum(Fraction(int(rec[k].subs(m, mval))) * a[mval + smin + k]
                      for k in range(K) if 0 <= mval + smin + k < N)
            if rhs != 0:
                raise RuntimeError(f"inconsistent at n={n}")
            continue
        rhs = sum(Fraction(int(rec[k].subs(m, mval))) * a[mval + smin + k]
                  for k in range(K) if 0 <= mval + smin + k < N)
        a[n] = -rhs / Fraction(int(top))
    return a


if __name__ == '__main__':
    # smoke test: L = (1-x)Dx - 1 annihilates 1/(1-x)
    L = [-1, 1 - x]
    ser = [Fraction(1)] * 20
    out = apply_to_series(L, ser)
    assert all(c == 0 for c in out), out
    print("smoke: apply_to_series OK")
    # LCLM(Dx-1, Dx-2) = Dx^2 - 3Dx + 2
    A = [-1, 1]
    B = [-2, 1]
    L = lclm2(A, B)
    print("LCLM(Dx-1,Dx-2) =", L)
    # right_divide
    Q = right_divide(L, A)
    print("L / A =", Q)
    print("smoke tests passed")

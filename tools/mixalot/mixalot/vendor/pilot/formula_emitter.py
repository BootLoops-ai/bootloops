#!/usr/bin/env python3
"""formula_emitter.py — THEOREM (constructive): for the 2-mixture of one
k-state variable with uniform priors, the LSX evidence integral is

    Z(U) = c_0 + sum_j c_j * H_{n_j}        (c's in Q, H = harmonic numbers)

with an explicit finite reduction. Derivation implemented here:
  Z = P * sum_m W(m) / [prod_{i=1}^{k-1}(m+i) * prod_{j=1}^{k-1}(N-m+j)],
  P = prod_v u_v! * (k-1)!^2 / (N+1)!,  W(m) = #{x <= U, |x| = m}.
Partial fractions (A_i = prod_{i'!=i} 1/(i'-i)):
  1/prod(m+i) = sum_i A_i/(m+i);  1/((m+i)(N-m+j)) = [1/(m+i)+1/(N-m+j)]/(N+i+j)
  and W(N-m)=W(m)  =>  Z = P * sum_{i,j} A_i A_j (S_i + S_j)/(N+i+j),
  S_i = sum_m W(m)/(m+i).
Inclusion-exclusion W(m) = sum_{T subset states} (-1)^{|T|} C(m-a_T+k-1, k-1),
a_T = sum_{v in T}(u_v+1) (C(n,r)=0 for n<r), and polynomial division
  q(m)/(m+i) = s(m) + q(-i)/(m+i)
turn each S_i into  rational + q(-i) * (H_{N+i} - H_{a_T+i-1}).

Emits {'const': Fraction, ('H', n): Fraction}; verified exactly against the
independent direct evaluator for k=2..5 on random count vectors, and against
the recorded k=2 closed form.
"""

import itertools
import random
import os
import sys
from fractions import Fraction
from math import comb, factorial

# sibling imports resolve to the vendored copies beside this file
# (MIXALOT_PILOT_DIR overrides)
sys.path.insert(0, os.environ.get("MIXALOT_PILOT_DIR",
                                  os.path.dirname(os.path.abspath(__file__))))
from zseries import z_1var_exact  # noqa: E402


def _binom_poly(k, shift):
    """Coefficients (ascending, Fractions) of q(m) = C(m - shift + k - 1, k - 1)
    as a polynomial in m: prod_{t=1}^{k-1} (m - shift + t) / (k-1)!."""
    coeffs = [Fraction(1)]
    for t in range(1, k):
        c = t - shift          # factor (m + c)
        new = [Fraction(0)] * (len(coeffs) + 1)
        for a, co in enumerate(coeffs):
            new[a + 1] += co
            new[a] += co * c
        coeffs = new
    inv = Fraction(1, factorial(k - 1))
    return [co * inv for co in coeffs]


def _poly_eval(coeffs, x):
    v = Fraction(0)
    for co in reversed(coeffs):
        v = v * x + co
    return v


def _poly_divide_linear(coeffs, i):
    """q(m) = (m+i) s(m) + R: synthetic division; returns (s_coeffs, R=q(-i))."""
    s = [Fraction(0)] * (len(coeffs) - 1)
    rem = Fraction(0)
    for a in range(len(coeffs) - 1, -1, -1):
        cur = coeffs[a] + rem
        if a == 0:
            return s, cur
        s[a - 1] = cur
        rem = -Fraction(i) * cur
    return s, rem


def _poly_sum_range(coeffs, lo, hi):
    """sum_{m=lo}^{hi} s(m), exact (evaluate term by term is O(N); fine here —
    use Faulhaber-free direct sum with Fractions for clarity/correctness)."""
    if hi < lo:
        return Fraction(0)
    # Hockey-stick style: expand in falling factorials would be faster; direct
    # summation is exact and adequate for N up to a few thousand.
    return sum(_poly_eval(coeffs, m) for m in range(lo, hi + 1))


def emit_formula(U):
    """Returns (formula_dict, P) with formula keys 'const' and ('H', n);
    Z = P * (formula evaluated with H_n = harmonic numbers)."""
    k = len(U)
    N = sum(U)
    A = {i: Fraction(1) for i in range(1, k)}
    for i in range(1, k):
        for ip in range(1, k):
            if ip != i:
                A[i] /= (ip - i)
    # S_i = const_i + sum coeff * H_n   for i = 1..k-1
    S = {}
    for i in range(1, k):
        const = Fraction(0)
        hterms = {}
        for tsize in range(k + 1):
            for T in itertools.combinations(range(k), tsize):
                aT = sum(U[v] + 1 for v in T)
                if aT > N:
                    continue
                sign = -1 if tsize % 2 else 1
                q = _binom_poly(k, aT)
                s, R = _poly_divide_linear(q, i)
                const += sign * _poly_sum_range(s, aT, N)
                if R:
                    # sum_{m=aT}^{N} 1/(m+i) = H_{N+i} - H_{aT+i-1}
                    hterms[N + i] = hterms.get(N + i, Fraction(0)) + sign * R
                    hterms[aT + i - 1] = hterms.get(aT + i - 1, Fraction(0)) - sign * R
        S[i] = (const, hterms)
    # assemble sum_{i,j} A_i A_j (S_i + S_j) / (N+i+j)
    out = {"const": Fraction(0)}
    for i in range(1, k):
        for j in range(1, k):
            w = A[i] * A[j] / (N + i + j)
            for src in (S[i], S[j]):
                out["const"] += w * src[0]
                for n_idx, co in src[1].items():
                    key = ("H", n_idx)
                    out[key] = out.get(key, Fraction(0)) + w * co
    out = {kk: vv for kk, vv in out.items() if kk == "const" or vv != 0}
    P = Fraction(1)
    for u in U:
        P *= factorial(u)
    P *= Fraction(factorial(k - 1) ** 2, factorial(N + 1))
    return out, P


def eval_formula(formula, P):
    def H(n):
        return sum(Fraction(1, t) for t in range(1, n + 1))
    tot = formula["const"]
    for key, co in formula.items():
        if key != "const":
            tot += co * H(key[1])
    return P * tot


def pretty(formula, P):
    parts = [f"{formula['const']}"]
    for key in sorted((kk for kk in formula if kk != "const"),
                      key=lambda t: t[1]):
        parts.append(f"({formula[key]})*H_{key[1]}")
    return f"({P}) * [ " + " + ".join(parts) + " ]"


def main():
    ok = True
    rng = random.Random(20260708)
    for k in (2, 3, 4, 5):
        for trial in range(6):
            U = [rng.randint(0, 9) for _ in range(k)]
            got = eval_formula(*emit_formula(U))
            want = z_1var_exact(U)
            good = got == want
            ok &= good
            if not good:
                print(f"  FAIL k={k} U={U}: {got} != {want}")
        print(f"k={k}: 6 random count vectors exact: {'PASS' if ok else 'FAIL'}")
    # cross-check vs the recorded k=2 closed form
    from closed_form_1var import Z_closed
    ok2 = all(eval_formula(*emit_formula([a, b])) == Z_closed(a, b)
              for a in range(8) for b in range(8))
    ok &= ok2
    print(f"k=2 emitter == recorded closed form [0..7]^2: {'PASS' if ok2 else 'FAIL'}")
    # display: the k=3 symmetric formula at (u,u,u), u=2 as an example
    f, P = emit_formula([2, 2, 2])
    print("example k=3, U=(2,2,2):")
    print(" ", pretty(f, P))
    print("  =", eval_formula(f, P), "== direct:", z_1var_exact([2, 2, 2]))
    print("FORMULA_EMITTER:", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""closed_form_1var.py — exact closed-form evidence formulas for the
2-mixture of one binary variable, uniform priors (LSX conventions, bare
integral eq (25)).

PROPOSITION (symmetric counts): Z(n,n) = (n!)²/(2n+1)! · [1/(n+1) + 2(H_{2n+1}−H_{n+1})]
PROPOSITION (general counts, u0<=u1, N=u0+u1):
  Z(u0,u1) = u0! u1!/(N+1)! · 2/(N+2) ·
             [ (u0+1) + (u0+1)(H_{u1+1}−H_{u0+1}) + (N+2)(H_{N+1}−H_{u1+1}) − (N−u1) ]
Derivation: the A=I collapse Z = u0!u1!/(N+1)! Σ_m W(m)/((m+1)(N−m+1)) with
W(m) = #{x<=U, |x|=m} (trapezoid), partial fractions 1/((m+1)(N−m+1)) =
[1/(m+1)+1/(N−m+1)]/(N+2), W symmetric ⇒ 2/(N+2)·Σ W(m)/(m+1), and the
trapezoid sum evaluates in harmonic numbers.

EXACT ASYMPTOTICS (from the operator's indicial equation at n=∞, roots
{−1/2,−3/2} after peeling 4^{−n}; one branch = n!/(4ⁿ(3/2)ₙ)):
  −log Z(n,n) = n log 4 + (1/2) log n − log(√π log 2) + O(1/n)
i.e. learning coefficient λ = 1/2 EXACT (indicial root), constant √π·log2 EXACT
(from H_{2n+1}−H_{n+1} → log 2 and B(n+1,n+1) ~ √(π/n)·4^{−n}).

Gates: byte-exact vs the independent direct evaluator on a grid; asymptotic
ratio → 1 numerically.
"""

import os
import sys
from fractions import Fraction
from math import factorial

# sibling imports resolve to the vendored copies beside this file
# (MIXALOT_PILOT_DIR overrides the default sibling path)
sys.path.insert(0, os.environ.get("MIXALOT_PILOT_DIR",
                                  os.path.dirname(os.path.abspath(__file__))))
from zseries import z_1var_exact  # noqa: E402


def H(k):
    return sum(Fraction(1, j) for j in range(1, k + 1))


def Z_closed_sym(n):
    return Fraction(factorial(n) ** 2, factorial(2 * n + 1)) * (
        Fraction(1, n + 1) + 2 * (H(2 * n + 1) - H(n + 1)))


def Z_closed(u0, u1):
    if u0 > u1:
        u0, u1 = u1, u0
    N = u0 + u1
    bracket = (Fraction(u0 + 1)
               + (u0 + 1) * (H(u1 + 1) - H(u0 + 1))
               + (N + 2) * (H(N + 1) - H(u1 + 1))
               - (N - u1))
    return Fraction(factorial(u0) * factorial(u1), factorial(N + 1)) * \
        Fraction(2, N + 2) * bracket


def main():
    ok = True
    # symmetric formula, n = 0..40
    for n in range(41):
        ok &= (Z_closed_sym(n) == z_1var_exact([n, n]))
    print(f"symmetric closed form n=0..40: {'PASS' if ok else 'FAIL'}")
    # general formula on the full grid [0..25]^2
    ok2 = True
    for u0 in range(26):
        for u1 in range(26):
            ok2 &= (Z_closed(u0, u1) == z_1var_exact([u0, u1]))
    print(f"general closed form [0..25]^2: {'PASS' if ok2 else 'FAIL'}")
    # consistency: general reduces to symmetric
    ok3 = all(Z_closed(n, n) == Z_closed_sym(n) for n in range(30))
    print(f"general == symmetric on diagonal: {'PASS' if ok3 else 'FAIL'}")
    # exact asymptotic constant check
    import mpmath as mp
    mp.mp.dps = 40
    n = 20000
    val = mp.beta(n + 1, n + 1) * (mp.mpf(1) / (n + 1)
                                   + 2 * (mp.harmonic(2 * n + 1) - mp.harmonic(n + 1)))
    ratio = val / (mp.sqrt(mp.pi) * mp.log(2) * mp.mpf(4) ** (-n) / mp.sqrt(n))
    print(f"asymptotic ratio z/(sqrt(pi)*log2*4^-n/sqrt(n)) at n=20000: "
          f"{mp.nstr(ratio, 12)} (→1)")
    allok = ok and ok2 and ok3
    print("CLOSED_FORM_1VAR:", "PASS" if allok else "FAIL")
    sys.exit(0 if allok else 1)


if __name__ == "__main__":
    main()

"""THE STEP verification: P-finite A-recurrence for K(D,A), exact vs oracle.

Derivation (NOTES_MATH.md): p_n(z) = P_n(z)/z (constant term 1),
R~(z) = prod_n p_n(z)^{Phi_n}  (so K(D, S+B) = [z^B] R~),
first-order ODE:  Q R~' = Q~ R~,
  Q  = prod_{distinct n} p_n            (deg D* = sum_{distinct}(n-1))
  Q~ = sum_{distinct n} Phi_n p_n' prod_{m != n} p_m
Coefficient recurrence (K~_B := K(D, S+B), K~_0 = 1):
  (B+1) K~_{B+1} = sum_j Q~_j K~_{B-j} - sum_{j>=1} Q_j (B-j+1) K~_{B-j+1}

Checks on 10 random datasets (exact Fractions):
  R1: ODE identity Q R~' - Q~ R~ == 0 as exact polynomial
  R2: forward recurrence from K~_0=1 reproduces ALL oracle K(D,A) exactly
  R3: contiguity P_{n+1}(z) = P_n(z) + z^2 P_n'(z)/n for n = 1..40
"""
import random
import sys
from fractions import Fraction
from math import factorial
from collections import Counter

from etienne_oracle import K_vector, stirling1_row
from gate_engines import random_dataset


def p_poly(n):
    """p_n(z) = P_n(z)/z as Fraction list, p[k] = coeff z^k, k=0..n-1."""
    row = stirling1_row(n)
    d = factorial(n - 1)
    return [Fraction(row[a] * factorial(a - 1), d) for a in range(1, n + 1)]


def poly_mul(a, b):
    out = [Fraction(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                out[i + j] += x * y
    return out


def poly_deriv(a):
    return [k * a[k] for k in range(1, len(a))] or [Fraction(0)]


random.seed(77)
fails = 0
datasets = [random_dataset(J) for J in (8, 12, 16, 20, 25, 30, 35, 40, 50, 60)]

# R3 contiguity
for n in range(1, 41):
    P = [Fraction(0)] + p_poly(n)          # P_n
    Pn1 = [Fraction(0)] + p_poly(n + 1)    # P_{n+1}
    dP = poly_deriv(P)
    rhs = P + [Fraction(0)]
    for k, c in enumerate(dP):             # + z^2 P'/n
        rhs[k + 2] += c / n
    if [x for x in rhs] != [x for x in Pn1] + [Fraction(0)] * (len(rhs) - len(Pn1)):
        fails += 1
        print(f"R3 FAIL n={n}")
print("R3: contiguity P_(n+1) = P_n + z^2 P_n'/n verified n=1..40")

for D in datasets:
    S, J = len(D), sum(D)
    cnt = Counter(D)
    dist = sorted(cnt)
    # Q, Q~
    Q = [Fraction(1)]
    for n in dist:
        Q = poly_mul(Q, p_poly(n))
    Qt = [Fraction(0)]
    for n in dist:
        term = poly_deriv(p_poly(n))
        term = [c * cnt[n] for c in term]
        for m in dist:
            if m != n:
                term = poly_mul(term, p_poly(m))
        L = max(len(Qt), len(term))
        Qt = [(Qt[k] if k < len(Qt) else 0) + (term[k] if k < len(term) else 0)
              for k in range(L)]
    Dstar = len(Q) - 1

    # oracle K
    S0, Kor = K_vector(tuple(D))
    Ktil = [Kor[B] for B in range(J - S + 1)]  # K~_B = K(D,S+B)

    # R1: ODE identity with the FULL R~ (from oracle coefficients)
    Rt = Ktil
    lhs = poly_mul(Q, poly_deriv(Rt))
    rhs = poly_mul(Qt, Rt)
    L = max(len(lhs), len(rhs))
    ok = all((lhs[k] if k < len(lhs) else 0) == (rhs[k] if k < len(rhs) else 0)
             for k in range(L))
    if not ok:
        fails += 1
        print(f"R1 FAIL D={D}")

    # R2: forward recurrence from K~_0 = 1
    K = [Fraction(1)]
    for B in range(J - S):
        acc = Fraction(0)
        for j in range(len(Qt)):
            if 0 <= B - j < len(K):
                acc += Qt[j] * K[B - j]
        for j in range(1, len(Q)):
            if 0 <= B - j + 1 < len(K):
                acc -= Q[j] * (B - j + 1) * K[B - j + 1]
        K.append(acc / (B + 1))
    if K != Ktil:
        fails += 1
        print(f"R2 FAIL D={D}")
    else:
        print(f"R2 OK  D(J={J},S={S},D*={Dstar}): {J-S+1} K-values reproduced exactly")

if fails:
    print(f"*** {fails} FAILURES ***")
    sys.exit(1)
print("RECURRENCE VERIFIED (R1 certificate identity, R2 forward generation, R3 contiguity)")

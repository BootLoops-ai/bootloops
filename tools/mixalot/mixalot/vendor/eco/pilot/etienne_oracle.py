"""Brute-force exact-rational oracle for the Etienne (2005) sampling formula.

TRANSCRIPTION STATUS: from the standard form of the formula; checked against
Etienne's own worked example in gate_worked_example.py.

Formula (Etienne 2005, Ecol. Lett. 8:253, eq. 6):

    P[D | theta, m, J] = J! / (prod_i n_i * prod_j Phi_j!)
                         * theta^S / (I)_J
                         * sum_{A=S}^{J} K(D,A) * I^A / (theta)_A

with
    I = m (J-1) / (1-m),          (x)_A = Pochhammer rising factorial,
    K(D,A) = sum_{{a_i}: sum a_i = A, 1<=a_i<=n_i}
             prod_i  s(n_i,a_i) s(a_i,1) / s(n_i,1)
           = sum_{{a_i}} prod_i  s(n_i,a_i) (a_i-1)! / (n_i-1)!

where s(n,a) are UNSIGNED Stirling numbers of the first kind, D is the
unlabelled abundance dataset {n_1..n_S}, Phi_j = #{i : n_i = j}.

Everything is exact rational (fractions.Fraction); no floats anywhere.
K(D,.) is computed by integer polynomial convolution of the per-species
vectors v_i(a) = s(n_i,a) (a-1)!  (integers), divided by prod_i (n_i-1)!.
Complexity O(J^2) big-int ops — the brute-force baseline the recurrence
route must beat.
"""
from fractions import Fraction
from functools import lru_cache
from math import factorial
from collections import Counter


@lru_cache(maxsize=None)
def stirling1_row(n):
    """Unsigned Stirling numbers of the first kind: row s(n, 0..n).

    x(x+1)...(x+n-1) = sum_a s(n,a) x^a.  s(n,0)=0 for n>=1; s(0,0)=1.
    Recurrence: s(n+1,a) = n*s(n,a) + s(n,a-1).
    """
    if n == 0:
        return (1,)
    prev = stirling1_row(n - 1)
    row = [0] * (n + 1)
    for a in range(n):
        row[a] += (n - 1) * prev[a]
        row[a + 1] += prev[a]
    return tuple(row)


@lru_cache(maxsize=None)
def species_vector(n):
    """Integer vector v(a) = s(n,a) * (a-1)!, a = 1..n (index 0 -> a=1)."""
    row = stirling1_row(n)
    return tuple(row[a] * factorial(a - 1) for a in range(1, n + 1))


def K_vector(D):
    """K(D,A) for A = S..J as a list of Fractions.

    Convolution of the per-species integer vectors, then division by
    prod_i (n_i-1)!.  Returns (A_min, [K(D,A_min), ..., K(D,J)]).
    """
    D = list(D)
    S = len(D)
    # integer convolution; poly[k] is coefficient of A = S + k offset handled below
    poly = [1]
    for n in D:
        v = species_vector(n)
        new = [0] * (len(poly) + len(v) - 1)
        for i, c in enumerate(poly):
            if c:
                for j, w in enumerate(v):
                    new[i + j] += c * w
        poly = new
    denom = 1
    for n in D:
        denom *= factorial(n - 1)
    return S, [Fraction(c, denom) for c in poly]


def pochhammer(x, k):
    """Rising factorial (x)_k, x Fraction/int, exact."""
    x = Fraction(x)
    out = Fraction(1)
    for i in range(k):
        out *= x + i
    return out


def etienne_P(D, theta, m=None, I=None):
    """Exact P[D | theta, m] (or pass I directly).  D = abundance list."""
    D = sorted(D)
    J = sum(D)
    S = len(D)
    theta = Fraction(theta)
    if I is None:
        m = Fraction(m)
        if m == 1:
            return ewens_P(D, theta)
        I = m * (J - 1) / (1 - m)
    else:
        I = Fraction(I)

    if I == 0:
        # m -> 0 limit: whole community descends from one immigrant
        # (lim I^A/(I)_J nonzero only for A=1, forcing S=1; limit value 1)
        return Fraction(1) if S == 1 else Fraction(0)

    # combinatorial prefactor J! / (prod n_i * prod Phi_j!)
    pref = Fraction(factorial(J))
    for n in D:
        pref /= n
    for cnt in Counter(D).values():
        pref /= factorial(cnt)

    S0, Kvec = K_vector(D)
    assert S0 == S and len(Kvec) == J - S + 1
    total = Fraction(0)
    IA = I ** S
    for idx, KA in enumerate(Kvec):
        A = S + idx
        total += KA * IA / pochhammer(theta, A)
        IA *= I
    return pref * theta ** S / pochhammer(I, J) * total


def ewens_P(D, theta):
    """Ewens sampling formula for the unlabelled abundance dataset D."""
    D = sorted(D)
    J = sum(D)
    S = len(D)
    theta = Fraction(theta)
    pref = Fraction(factorial(J))
    for n in D:
        pref /= n
    for cnt in Counter(D).values():
        pref /= factorial(cnt)
    return pref * theta ** S / pochhammer(theta, J)


def partitions(J, max_part=None):
    """All partitions of J (as sorted tuples), for unitarity sums."""
    if max_part is None:
        max_part = J
    if J == 0:
        yield ()
        return
    for k in range(min(J, max_part), 0, -1):
        for rest in partitions(J - k, k):
            yield tuple(sorted(rest + (k,)))


if __name__ == "__main__":
    # smoke: singleton dataset
    print(etienne_P([3, 2, 1], Fraction(3, 2), m=Fraction(1, 3)))

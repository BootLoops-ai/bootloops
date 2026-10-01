"""Multi-sample Etienne (2007/2009) formula + independent model DP oracle.

Formula (Etienne 2009 JTB eqs 1-6, transcription per LITREVIEW_MULTISAMPLE):
  P[D | I_1..I_N, theta] = (1/prod_vec Phi_vec!) * prod_i [ J_i! / ((I_i)_{J_i} prod_k n_ik!) ]
      * sum_A M(D,A,I) * theta^S / (theta)_A
  M(D,A,I) = sum_{{a_ik}: sum=A} prod_k [ (a_k-1)! prod_i s(n_ik,a_ik) I_i^{a_ik} ],
  a_k = sum_i a_ik;  a_ik in 1..n_ik if n_ik>0 else 0.

Engine (nested positive convolutions — NOTES_MULTISAMPLE):
  per species k:  u_k(z) = prod_i RF_{n_ik}(I_i z)      (convolve across samples)
                  S_k(z) = sum_b (b-1)! [z^b]u_k * z^b   (Borel twist)
  M(D,A,I) = [z^A] prod_k S_k(z)                         (convolve across species)
with RF_n(z) = sum_a s(n,a) z^a = z(z+1)...(z+n-1) coefficients, RF_0 = 1.

DP oracle: sequential two-level construction, M samples sharing one
metacommunity Hoppe urn. State per species: (n_1..n_N, a) with a = total
metacommunity lineages (only the per-species TOTAL a matters for dynamics —
the formula's (a_k-1)! structure reflects exactly this).

Everything exact (Fractions).
"""
from fractions import Fraction
from functools import lru_cache
from math import factorial
from collections import defaultdict, Counter
from itertools import product as iproduct


@lru_cache(maxsize=None)
def rf_row(n):
    """Unsigned Stirling row: coefficients of z^0..z^n in z(z+1)...(z+n-1)."""
    if n == 0:
        return (1,)
    prev = rf_row(n - 1)
    row = [0] * (n + 1)
    for a in range(n):
        row[a] += (n - 1) * prev[a]
        row[a + 1] += prev[a]
    return tuple(row)


def poly_mul(a, b):
    out = [Fraction(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                out[i + j] += x * y
    return out


def species_poly(nvec, Ivec):
    """S_k(z) coefficients for one species with per-sample counts nvec."""
    u = [Fraction(1)]
    for n, I in zip(nvec, Ivec):
        if n == 0:
            continue
        row = rf_row(n)
        u = poly_mul(u, [row[a] * I ** a for a in range(n + 1)])
    # Borel twist (b-1)! for b >= 1; b=0 coeff must vanish for occupied species
    return [u[b] * factorial(b - 1) if b >= 1 else Fraction(0)
            for b in range(len(u))]


def pochhammer(x, k):
    x = Fraction(x)
    out = Fraction(1)
    for i in range(k):
        out *= x + i
    return out


def multisample_P(Dmat, theta, Ivec, ij_norm="per-sample"):
    """Exact P[D] via the convolution engine.

    Dmat: list of species rows, each a tuple of per-sample counts n_ik
          (k = species index = row; i = sample index = column).
    ij_norm: 'per-sample'  -> prod_i (I_i)_{J_i}   [our reading of eq 5]
             'pooled'      -> (I)_{Jtot} with common I [literal print, equal-I only]
    """
    theta = Fraction(theta)
    Ivec = [Fraction(I) for I in Ivec]
    N = len(Ivec)
    S = len(Dmat)
    Js = [sum(row[i] for row in Dmat) for i in range(N)]

    # prefactor
    pref = Fraction(1)
    for cnt in Counter(tuple(r) for r in Dmat).values():
        pref /= factorial(cnt)
    for i in range(N):
        pref *= factorial(Js[i])
        for row in Dmat:
            pref /= factorial(row[i])
    if ij_norm == "per-sample":
        for i in range(N):
            pref /= pochhammer(Ivec[i], Js[i])
    elif ij_norm == "pooled":
        assert len(set(Ivec)) == 1
        pref /= pochhammer(Ivec[0], sum(Js))
    else:
        raise ValueError(ij_norm)

    # M-generating polynomial
    M = [Fraction(1)]
    for row in Dmat:
        M = poly_mul(M, species_poly(row, Ivec))
    total = Fraction(0)
    for A in range(S, len(M)):
        if M[A]:
            total += M[A] * theta ** S / pochhammer(theta, A)
    return pref * total


# ---------------- independent model DP oracle ----------------

def process_distribution_multi(Jvec, theta, Ivec):
    """Exact distribution over unlabelled abundance-VECTOR multisets.

    Builds samples sequentially (exchangeability makes order irrelevant);
    state = sorted tuple of (n_1..n_N, a) per species.
    Returns dict: multiset-of-vectors key -> Fraction.
    """
    theta = Fraction(theta)
    Ivec = [Fraction(I) for I in Ivec]
    N = len(Jvec)
    states = {(): Fraction(1)}
    for i in range(N):
        I = Ivec[i]
        for j in range(Jvec[i]):  # adding individual j+1 of sample i
            new = defaultdict(Fraction)
            for st, p in states.items():
                A = sum(sp[-1] for sp in st)
                sp_list = [list(s) for s in st]
                if j > 0:
                    # local copy: species c with weight n_ic
                    for c, sp in enumerate(sp_list):
                        if sp[i] == 0:
                            continue
                        q = p * Fraction(j, 1) / (I + j) * Fraction(sp[i], j)
                        ns = [tuple(x) for x in sp_list]
                        ns[c] = tuple(sp[:i] + [sp[i] + 1] + sp[i + 1:-1] + [sp[-1]])
                        new[tuple(sorted(ns))] += q
                # immigrant, existing metacommunity species c (weight a_c)
                for c, sp in enumerate(sp_list):
                    q = p * I / (I + j) * Fraction(sp[-1]) / (theta + A)
                    ns = [tuple(x) for x in sp_list]
                    ns[c] = tuple(sp[:i] + [sp[i] + 1] + sp[i + 1:-1] + [sp[-1] + 1])
                    new[tuple(sorted(ns))] += q
                # immigrant, new species
                q = p * I / (I + j) * theta / (theta + A)
                fresh = [0] * N + [1]
                fresh[i] = 1
                new[tuple(sorted(st + (tuple(fresh),)))] += q
            states = dict(new)
    out = defaultdict(Fraction)
    for st, p in states.items():
        key = tuple(sorted(tuple(sp[:-1]) for sp in st))
        out[key] += p
    return dict(out)

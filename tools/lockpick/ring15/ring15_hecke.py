# lockpick ring15 member — exact Hecke a_n of the two weight-3 level-15 CM newforms.
"""ring15_hecke: a_n of the two weight-3 level-15 CM newforms (exact integers).

Grossencharacter construction, K = Q(sqrt(-15)), h=2, w=2, weight k=3:
  level 15 = |D|*N(cond)^2 forces conductor (1); infinity type alpha^2.
  psi((alpha)) = alpha^2 on principal ideals (well defined since (-a)^2=a^2).
  p2 = (2,beta), beta = (1+sqrt(-15))/2, p2^2 = (beta);  psi(p2) = s*beta, s=+-1.
  The two signs give the two rational newforms (genus-character pair).
Element rep: (A,B) := (A + B*sqrt(-15))/2 with A=B mod 2; all arithmetic exact.
  split principal p (x^2+15y^2=4p):   a_p = (x^2-15y^2)/2
  split non-principal p:  alpha generates p*p2 (norm 2p, X^2+15Y^2=8p),
    filtered by  p2|(alpha) and p2bar∤(alpha);  psi(p) = alpha^2*conj(beta)/(4s);
    a_p = 2*Re(psi) = s*A(alpha^2 conj(beta))/4.
  ramified p=3,5: same route, psi real;  a_p = s*A/8.   inert: a_p = 0.
  a_{p^{r+1}} = a_p a_{p^r} - chi(p) p^2 a_{p^{r-1}};  a_n multiplicative.
Validation (run_validation): exact equality of n<=NCMP coefficients against the
eta products f_s = eta(3t)^3 eta(5t)^3 + s*eta(t)^3 eta(15t)^3 (Ligozat: both
factors live in S_3(Gamma_0(15), chi_{-15})), and against PARI/GP mfeigenbasis
(vectors recorded below verbatim from PARI/GP 2.15.4 mfeigenbasis).
"""
from math import isqrt
from .ring15_core import chi15

BETA = (1, 1)

PARI_FORM_A2_PLUS1 = [1, 1, -3, -3, 5, -3, 0, -7, 9, 5,
                      0, 9, 0, 0, -15, 5, -14, 9, -22, -15]   # n=1..20
PARI_FORM_A2_MINUS1 = [1, -1, 3, -3, -5, -3, 0, 7, 9, 5,
                       0, -9, 0, 0, -15, 5, 14, -9, -22, 15]  # n=1..20


def _mul(u, v):
    A1, B1 = u
    A2, B2 = v
    return ((A1 * A2 - 15 * B1 * B2) // 2, (A1 * B2 + A2 * B1) // 2)


def _conj(u):
    return (u[0], -u[1])


def _half_integral(u):
    """Is ((A+B sqrt(-15))/2)/2 in the ring of integers?"""
    A, B = u
    return A % 2 == 0 and B % 2 == 0 and ((A // 2) - (B // 2)) % 2 == 0


def _solve(N):
    """All (X,Y), X,Y >= 0, X^2 + 15 Y^2 = N, X = Y mod 2."""
    out, Y = [], 0
    while 15 * Y * Y <= N:
        r = N - 15 * Y * Y
        X = isqrt(r)
        if X * X == r and (X - Y) % 2 == 0:
            out.append((X, Y))
        Y += 1
    return out


def a_p(p, s):
    c = chi15(p)
    if c == -1:
        return 0
    if c == 1:
        prin = [v for v in _solve(4 * p)]
        if prin:
            X, Y = prin[0]
            return (X * X - 15 * Y * Y) // 2
    # non-principal split (incl. p=2) or ramified (p=3,5)
    vals = set()
    for X, Y in _solve(8 * p):
        for al in {(X, Y), (X, -Y)}:
            if not (_half_integral(_mul(al, _conj(BETA)))
                    and not _half_integral(_mul(al, BETA))):
                continue   # keep only p2 | (alpha), p2bar ∤ (alpha)
            A, B = _mul(_mul(al, al), _conj(BETA))
            if c == 0:                      # ramified: psi(p) real, a_p = psi
                assert B == 0, (p, al)
                assert A % 8 == 0
                vals.add(s * (A // 8))
            else:                           # split: a_p = 2 Re psi
                assert A % 4 == 0
                vals.add(s * (A // 4))
    assert len(vals) == 1, (p, vals)
    return vals.pop()


def a_n_list(N, s):
    """[a_0(=0), a_1, ..., a_N] via smallest-prime-factor sieve, exact ints."""
    spf = list(range(N + 1))
    for i in range(2, isqrt(N) + 1):
        if spf[i] == i:
            for j in range(i * i, N + 1, i):
                if spf[j] == j:
                    spf[j] = i
    a = [0] * (N + 1)
    a[1] = 1
    ap_cache = {}
    for n in range(2, N + 1):
        p = spf[n]
        m, e = n, 0
        while m % p == 0:
            m //= p
            e += 1
        if p not in ap_cache:
            ap_cache[p] = a_p(p, s)
        if m > 1:
            a[n] = a[m] * a[n // m]         # coprime factorization
        else:                               # n = p^e: Hecke recursion
            ap, ch = ap_cache[p], chi15(p)
            aprev, acur = 1, ap
            for _ in range(e - 1):
                aprev, acur = acur, ap * acur - ch * p * p * aprev
            a[n] = acur
    return a


def _eta3(N):
    c = [0] * (N + 1)
    k = 0
    while k * (k + 1) // 2 <= N:
        c[k * (k + 1) // 2] = (-1) ** k * (2 * k + 1)
        k += 1
    return c


def _scale(c, m, N):
    out = [0] * (N + 1)
    for i, v in enumerate(c):
        if i * m <= N:
            out[i * m] = v
    return out


def _mulpoly(a, b, N):
    out = [0] * (N + 1)
    for i, va in enumerate(a):
        if va:
            for j in range(0, min(N - i, N) + 1):
                if b[j]:
                    out[i + j] += va * b[j]
    return out


def _shift(c, k, N):
    return [0] * k + c[:N + 1 - k]


def eta_products(N):
    """E35 = q prod(1-q^{3n})^3(1-q^{5n})^3,  E115 = q^2 prod(1-q^n)^3(1-q^{15n})^3."""
    e3 = _eta3(N)
    E35 = _shift(_mulpoly(_scale(e3, 3, N), _scale(e3, 5, N), N), 1, N)
    E115 = _shift(_mulpoly(e3, _scale(e3, 15, N), N), 2, N)
    return E35, E115


def run_validation(NCMP=60):
    """Exact 3-way comparison: Grossencharacter vs eta products vs PARI."""
    E35, E115 = eta_products(NCMP)
    rep = {}
    for s in (+1, -1):
        an = a_n_list(NCMP, s)
        eta_f = [E35[n] + s * E115[n] for n in range(NCMP + 1)]
        rep[f's{s:+d}_eta_match'] = (an == eta_f)
        pari = PARI_FORM_A2_PLUS1 if s == 1 else PARI_FORM_A2_MINUS1
        rep[f's{s:+d}_pari_match'] = (an[1:21] == pari)
        rep[f's{s:+d}_an_1_20'] = an[1:21]
    rep['all_ok'] = all(v for k, v in rep.items() if k.endswith('_match'))
    return rep

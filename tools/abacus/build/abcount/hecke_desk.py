"""hecke_desk — in-house Hecke-character / Jacobi-sum evaluation engine
(M3T_REGISTRATION.md sec 1 Builds item 3a).

Exact desk arithmetic over cyclotomic residue rings; PARI-INDEPENDENT — shares
no code with `ellap` / `hyperellcharpoly` beyond integer arithmetic.  The named
route-D implementation for battery members B3 (this stage, S2), B4/B4s/B5 (S3).

S2 scope: Z[zeta_5] Jacobi sums for B3 = J(C5), C5: y^2 = x^5 + 1, under the
COMMITTED SIGN PIN (review item 3, battery row B3 — pinned in the spec
before any run):

    Frobenius eigenvalues  alpha_a = -J(phi, chi^a),  a = 1..4,

phi the quadratic character (Euler criterion), chi the quintic character with
chi(g) = zeta_5 for the canonical generator g (SEEDS generator_rule: smallest
g fully primitive mod p).  J(phi, chi^a) = sum_{x != 0,1} phi(x) chi^a(1-x),
an algebraic integer of Q(zeta_10) = Q(zeta_5), computed in the exact basis
(1, z, z^2, z^3), z^4 = -1 - z - z^2 - z^3.

The charpoly prod_a (T - alpha_a) is expanded EXACTLY in Z[zeta_5]; its
coefficients must land in Z (the alpha_a are a Galois orbit) — a non-rational
coefficient raises FAIL-HECKE-DESK (never silently truncated).  A different
primitive root permutes {alpha_a} (a -> a*k mod 5), so the charpoly is
generator-independent; the J-sum VECTORS are reported for the canonical g.

S3 extension point: O_K = Z[sqrt(-5)] Hecke characters for B4/B4s ride the
same exact-desk substrate (this module), per the registration.
"""


def factorize(n):
    fac = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            fac[d] = fac.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        fac[n] = fac.get(n, 0) + 1
    return fac


def canonical_generator(p):
    """SEEDS generator_rule: smallest g >= 2 with g^((p-1)/q) != 1 mod p for
    ALL primes q | p-1 (full primitivity — p = 31 needs q = 3 too)."""
    qs = list(factorize(p - 1))
    for g in range(2, p):
        if all(pow(g, (p - 1) // q, p) != 1 for q in qs):
            return g
    raise ValueError("FAIL-HECKE-DESK: no primitive root found (p not prime?)")


# ---------------------------------------------------------------------------
# exact Z[zeta_5] arithmetic, basis (1, z, z^2, z^3), z^4 = -(1 + z + z^2 + z^3)
# ---------------------------------------------------------------------------

def z5_mul(a, b):
    raw = [0] * 7
    for i in range(4):
        if a[i] == 0:
            continue
        for j in range(4):
            raw[i + j] += a[i] * b[j]
    out = raw[:4]
    for k, c in ((4, raw[4]), (5, raw[5]), (6, raw[6])):
        if c == 0:
            continue
        if k == 4:                       # z^4 = -1 - z - z^2 - z^3
            for i in range(4):
                out[i] -= c
        elif k == 5:                     # z^5 = 1
            out[0] += c
        else:                            # z^6 = z
            out[1] += c
    return out


def z5_add(a, b):
    return [x + y for x, y in zip(a, b)]


ZETA_POW = ([1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0],
            [0, 0, 0, 1], [-1, -1, -1, -1])          # zeta_5^0 .. zeta_5^4


def jacobi_sum_vectors(p):
    """J(phi, chi^a) for a = 1..4 as Z[zeta_5] vectors (canonical g).
    Precondition: p = 1 mod 10 (else UNDECIDED-PRECONDITION)."""
    if p % 10 != 1:
        raise ValueError("UNDECIDED-PRECONDITION: p = %d not 1 mod 10" % p)
    g = canonical_generator(p)
    dlog = {}
    v = 1
    for k in range(p - 1):
        if v in dlog:
            raise ValueError("FAIL-HECKE-DESK: generator not primitive")
        dlog[v] = k
        v = v * g % p
    J = {}
    for a in (1, 2, 3, 4):
        s = [0, 0, 0, 0]
        for x in range(2, p):            # x != 0, 1 (and 1 - x != 0)
            sgn = -1 if dlog[x] % 2 else 1
            zp = ZETA_POW[(a * dlog[(1 - x) % p]) % 5]
            s = z5_add(s, [sgn * c for c in zp])
        J[a] = s
    return J, g


def charpoly_from_eigen_sign(J, sign):
    """Expand prod_{a=1..4} (T - sign*J_a) exactly in Z[zeta_5]; return the
    monic integer coefficient list [1, c3, c2, c1, c0] (descending powers).
    Raises FAIL-HECKE-DESK if any coefficient is non-rational."""
    poly = [[1, 0, 0, 0]]
    for a in (1, 2, 3, 4):
        neg_alpha = [-sign * c for c in J[a]]
        newp = [[0, 0, 0, 0] for _ in range(len(poly) + 1)]
        for i, co in enumerate(poly):
            newp[i] = z5_add(newp[i], co)                 # co * T
            newp[i + 1] = z5_add(newp[i + 1], z5_mul(co, neg_alpha))
        poly = newp
    coeffs = []
    for co in poly:
        if co[1] != 0 or co[2] != 0 or co[3] != 0:
            raise ValueError("FAIL-HECKE-DESK: non-rational charpoly coefficient %r" % (co,))
        coeffs.append(co[0])
    return coeffs


def b3_charpoly(p):
    """The B3 route-D result under the committed pin alpha_a = -J(phi, chi^a),
    plus the +J discrimination variant (SEEDS pin_discrimination_receipt).
    Returns (charpoly_pinned, charpoly_plusJ, receipt)."""
    J, g = jacobi_sum_vectors(p)
    pinned = charpoly_from_eigen_sign(J, -1)
    plus = charpoly_from_eigen_sign(J, +1)
    odd_differ = (pinned[1] != plus[1]) or (pinned[3] != plus[3])
    even_equal = (pinned[0] == plus[0] and pinned[2] == plus[2]
                  and pinned[4] == plus[4])
    receipt = {
        "engine": "hecke_desk (in-house exact Z[zeta_5] desk arithmetic; NO PARI code path; Builds item 3a)",
        "sign_pin": "alpha_a = -J(phi, chi^a) (committed pin — review item 3)",
        "generator_canonical_g": g,
        "jsum_vectors_basis_1_z_z2_z3": {("a%d" % a): J[a] for a in J},
        "charpoly_pinned_minusJ": pinned,
        "charpoly_plusJ_variant": plus,
        "pin_discrimination": {"odd_coefficients_differ": odd_differ,
                               "even_coefficients_equal": even_equal},
    }
    return pinned, plus, receipt

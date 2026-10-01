# Tonelli-Shanks sqrt mod p and quadratic-point reduction.
def modsqrt(a, p):
    a %= p
    if a == 0:
        return 0
    if pow(a, (p - 1) // 2, p) != 1:
        return None
    if p % 4 == 3:
        return pow(a, (p + 1) // 4, p)
    q, s = p - 1, 0
    while q % 2 == 0:
        q //= 2; s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    m, c, t, r = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
    while t != 1:
        i, tt = 0, t
        while tt != 1:
            tt = tt * tt % p; i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c, t, r = i, b * b % p, t * b * b % p, r * b % p
    return r

def reduce_quad(u, v, d, p):
    """Reductions of u + v*sqrt(d) mod p (u,v rationals as Fractions or ints).
    Returns [] if d is a nonresidue (inert) else the two residues (may coincide)."""
    from fractions import Fraction
    r = modsqrt(d % p, p)
    if r is None:
        return []
    def red(x):
        f = Fraction(x)
        return f.numerator * pow(f.denominator, p - 2, p) % p
    return sorted({(red(u) + red(v) * r) % p, (red(u) - red(v) * r) % p})

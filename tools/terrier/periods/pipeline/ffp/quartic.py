# Frobenius quartic tools. R(T) = 1 + a T + b p T^2 + a p^3 T^3 + p^6 T^4.
# Rank-2 split: R = (1 - alpha p T + p^3 T^2)(1 - beta T + p^3 T^2)
#   <=> a = -(alpha p + beta), b = 2 p^2 + alpha beta.
import math

def isqrt_exact(n):
    if n < 0:
        return None
    r = math.isqrt(n)
    return r if r * r == n else None

def weil_ok(a, b, p):
    """Necessary bounds: |a| <= 4 p^{3/2} + slack; |b| <= 6 p^2 + slack (e2 = b p, |e2| <= 6 p^3)."""
    return abs(a) <= 4 * p ** 1.5 + 1 and abs(b) <= 6 * p ** 2 + 1

def split_12(a, b, p):
    """Return list of (alpha, beta) integer splits with weight-(1,3) pairing, or []."""
    out = []
    disc = a * a - 4 * p * (b - 2 * p * p)
    r = isqrt_exact(disc)
    if r is None:
        return out
    for sgn in ((1, -1) if r else (1,)):
        num = -a + sgn * r
        if num % (2 * p):
            continue
        alpha = num // (2 * p)
        beta = -a - alpha * p
        if alpha * alpha <= 4 * p and beta * beta <= 4 * p ** 3:
            out.append((alpha, beta))
    return out

if __name__ == '__main__':
    import groundtruth_p19 as gt
    for phi in sorted(gt.ROWS):
        v = gt.ab(phi)
        if v is None:
            continue
        a, b = v
        print(phi, (a, b), split_12(a, b, 19))

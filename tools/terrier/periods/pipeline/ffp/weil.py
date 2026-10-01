# Exact Weil-circle test for R = 1 + aT + bpT^2 + ap^3T^3 + p^6T^4 (weight-3, smooth fibre).
# w = p^{3/2} T: R = w^4 + at w^3 + (b/p^2) w^2 + at w + 1, at = a/p^{3/2};
# u = w + 1/w: g(u) = u^2 + at u + (b/p^2 - 2); need both roots real in [-2,2].
def _ge0(c_sqrtp, c_1, p):
    """Test c_sqrtp*sqrt(p) + c_1 >= 0 exactly (integer coefficients)."""
    if c_sqrtp >= 0 and c_1 >= 0:
        return True
    if c_sqrtp <= 0 and c_1 <= 0:
        return False
    if c_sqrtp > 0:
        return c_sqrtp * c_sqrtp * p >= c_1 * c_1
    return c_sqrtp * c_sqrtp * p <= c_1 * c_1

def weil_circle(a, b, p):
    """All four roots on |lambda| = p^{3/2}. Exact integer arithmetic."""
    if a * a - 4 * p * b + 8 * p ** 3 < 0:        # disc(g) >= 0
        return False
    if not _ge0(2 * p * p + b, 2 * a, p):         # g(2) >= 0
        return False
    if not _ge0(2 * p * p + b, -2 * a, p):        # g(-2) >= 0
        return False
    return a * a <= 16 * p ** 3                   # vertex |at| <= 4

if __name__ == '__main__':
    import groundtruth_p19 as gt
    ok = all(weil_circle(*gt.ab(f), 19) for f in gt.ROWS if gt.ab(f))
    print("p=19 ground truth all pass weil_circle:", ok)
    # perturbed b must fail for most rows
    bad = sum(weil_circle(a, b + 40, 19) for (a, b) in (gt.ab(f) for f in gt.ROWS if gt.ab(f)))
    print("rows passing with b+40 (want few):", bad, "/15")

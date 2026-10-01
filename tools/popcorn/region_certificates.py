#!/usr/bin/env python3
"""popcorn.certificates.region — exact-rational REGION certificates.

Certificates for REGION claims — "the data point is provably at least this
far from the WHOLE model region", not just from a grid of samples.

WHAT IT CERTIFIES (four legs, fully generic):
  1. BOX POSITIVITY  p(s,u) >= 0 on [0,1]^2 by exact-rational Bernstein
     coefficients + dyadic de-Casteljau subdivision. Pure Fractions end to end.
  2. LP-DUAL RATIONALIZATION  take the float dual functional phi of a distance/
     membership LP (a float LP solver is PROPOSER-ONLY), rationalize it, and
     certify c*T - Phi >= 0 on the box: every point of the region
     conv{v(s,u)/T(s,u)} then satisfies phi.q <= c — a certified bound on the
     TRUE convex hull, valid regardless of the float LP's optimality
     ("validity independent of optimality" is the load-bearing sentence).
  3. CAUCHY-SCHWARZ POWER BOUND  chi2(v) = d'S^{-1}d >= (phi.target - c)^2 /
     (phi'S phi) for every region point v — a certified lower bound on the
     continuum chi-square from ONE certified sup. phi'S phi enters exactly via
     Fraction(float) (exact binary value, no rounding).
  4. RATIONAL RECONSTRUCTION + STURM  reconstruct a rational function P/Q from
     exact samples with HELD-OUT degree certification, then count its real
     roots in an open interval by Sturm's theorem (sympy exact isolation over
     QQ) with a denominator no-pole check — turns "we scanned and saw one root"
     into "exactly this many roots, certified".

MEASURED CAVEATS (do not soften):
  - Hull relaxation cost: measured on an n=8 test problem, the hull retained
    52-70% of the grid separation — a region certificate captures most, not
    all, of a grid fit's power.
  - Edge-aware sampling: a first run failed to seed because a corner of the
    parameter box lay outside every interior sample set; charts sampled only
    at interior points inner-approximate MORE than stated. Sample edges and
    corners too.
  - bern_nonneg is subdivision-only: a polynomial that TOUCHES zero on the box
    can be certified-nonnegative only with root deflation — that route lives
    in the sibling module positivity.py (route B, touch-root deflation); this
    module's kernel is for strictly-positive margins.

USAGE
  from popcorn import certificates
  R = certificates.region
  b = R.to_bernstein(poly, ds, du); R.bern_nonneg(b, depth)
  c = R.certified_sup(Ppoly, Tpoly, lo, hi, depth=7)      # exact Fraction
  chi2 = R.cs_chi2_lower_bound(phi_dot_target, c, phi_S_phi)
  P, Q, d = R.reconstruct_rational(f, dmax=40)            # exact samples of f
  nr, info = R.count_roots(P, lo, hi)                     # Sturm, open interval
Polynomials are dicts {(es, eu): Fraction} in the power basis on [0,1]^2.
Selftest: `python3 region_certificates.py` — prints REGION CERTS PASS, rc=0
(any failure: rc=1). Includes planted MUST-FAIL negative controls. The package
battery (selftest.py) runs it as its region-certificates leg.
"""
import sys
from fractions import Fraction as F
from math import comb

__all__ = [
    "poly_add", "poly_eval", "to_bernstein", "decasteljau_split",
    "bern_nonneg", "certified_sup", "rationalize", "cs_chi2_lower_bound",
    "solve_exact", "reconstruct_rational", "count_roots",
]


# ---------------- leg 1: exact Bernstein machinery --------------------------
def poly_add(a, b, coef=F(1)):
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, F(0)) + coef * v
    return {k: v for k, v in out.items() if v != 0}


def poly_eval(p, s, u):
    return sum(c * s ** es * u ** eu for (es, eu), c in p.items())


def to_bernstein(p, ds, du):
    """Bernstein coefficients of poly p (power basis dict) on [0,1]^2, exact.
    b[i][j] = sum_{k<=i, l<=j} a_{kl} C(i,k)C(j,l)/(C(ds,k)C(du,l))."""
    a = [[F(0)] * (du + 1) for _ in range(ds + 1)]
    for (es, eu), c in p.items():
        a[es][eu] = c
    b = [[F(0)] * (du + 1) for _ in range(ds + 1)]
    for i in range(ds + 1):
        for j in range(du + 1):
            tot = F(0)
            for k in range(i + 1):
                for l in range(j + 1):
                    if a[k][l] != 0:
                        tot += (a[k][l] * comb(i, k) * comb(j, l)
                                / (comb(ds, k) * comb(du, l)))
            b[i][j] = tot
    return b


def decasteljau_split(b, axis):
    """Split Bernstein patch in half along axis (0=s, 1=u); exact."""
    if axis == 1:
        rows = [list(r) for r in b]
    else:
        rows = [list(r) for r in zip(*b)]
    left, right = [], []
    for row in rows:
        n = len(row) - 1
        cur = list(row)
        lo = [cur[0]]
        hi = [cur[-1]]
        for _ in range(n):
            cur = [(cur[i] + cur[i + 1]) / 2 for i in range(len(cur) - 1)]
            lo.append(cur[0])
            hi.append(cur[-1])
        left.append(lo)
        right.append(list(reversed(hi)))
    if axis == 1:
        return left, right
    return ([list(r) for r in zip(*left)], [list(r) for r in zip(*right)])


def bern_nonneg(b, depth):
    """Certify all-values >= 0 on the patch via Bernstein coeffs + subdivision.
    True = CERTIFIED nonnegative. False = NOT certified (negative corner value
    is an exact refutation; depth exhaustion is merely inconclusive)."""
    mn = min(min(r) for r in b)
    if mn >= 0:
        return True
    if b[0][0] < 0 or b[0][-1] < 0 or b[-1][0] < 0 or b[-1][-1] < 0:
        return False
    if depth == 0:
        return False
    l, r = decasteljau_split(b, 0)
    ll, lr = decasteljau_split(l, 1)
    rl, rr = decasteljau_split(r, 1)
    return all(bern_nonneg(x, depth - 1) for x in (ll, lr, rl, rr))


# ---------------- leg 2: certified sup of a rational chart ------------------
def certified_sup(Ppoly, Tpoly, lo, hi, depth=6, iters=40):
    """Smallest c (binary search, exact rationals) with certificate
    c*T - P >= 0 on [0,1]^2. lo must FAIL, hi must PASS. With T >= 0 certified
    separately, every region point q = v/T then satisfies phi.q <= c."""
    ds = max(max(es for (es, eu) in p) for p in (Ppoly, Tpoly))
    du = max(max(eu for (es, eu) in p) for p in (Ppoly, Tpoly))

    def passes(c):
        g = poly_add({k: c * v for k, v in Tpoly.items()}, Ppoly, F(-1))
        b = to_bernstein(g, ds, du)
        return bern_nonneg(b, depth)

    assert passes(hi), "upper seed does not certify; raise hi or depth"
    for _ in range(iters):
        mid = (lo + hi) / 2
        if passes(mid):
            hi = mid
        else:
            lo = mid
    return hi


def rationalize(x, max_den=10 ** 12):
    """Float -> Fraction for DUAL functionals (validity independent of
    optimality; the exactness of the downstream certificate never depends on
    how well this approximates the float)."""
    return F(x).limit_denominator(max_den)


# ---------------- leg 3: Cauchy-Schwarz certified power bound ---------------
def cs_chi2_lower_bound(phi_dot_target, cert_sup, phi_S_phi):
    """chi2_continuum >= (phi.target - certified_sup)^2 / (phi' S phi), all
    Fractions (convert float S-quadratic-forms via Fraction(float) — exact
    binary value, no rounding). Returns Fraction(0) when the margin is
    nonpositive (no certified separation)."""
    margin = phi_dot_target - cert_sup
    if margin <= 0:
        return F(0)
    return margin * margin / phi_S_phi


# ---------------- leg 4: rational reconstruction + Sturm --------------------
def solve_exact(A, b):
    """exact Fraction least-structure solve of A x = b (A: list of rows)."""
    m, n = len(A), len(A[0])
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    piv_cols, r = [], 0
    for c in range(n):
        p = next((k for k in range(r, m) if M[k][c] != 0), None)
        if p is None:
            continue
        M[r], M[p] = M[p], M[r]
        inv = F(1) / M[r][c]
        M[r] = [x * inv for x in M[r]]
        for k in range(m):
            if k != r and M[k][c] != 0:
                f = M[k][c]
                M[k] = [M[k][j] - f * M[r][j] for j in range(n + 1)]
        piv_cols.append(c)
        r += 1
        if r == m:
            break
    for k in range(r, m):
        if M[k][n] != 0:
            return None                     # inconsistent
    x = [F(0)] * n
    for i, c in enumerate(piv_cols):
        x[c] = M[i][n]
    return x


def reconstruct_rational(f, pts=None, dmax=40, n_held=6):
    """P, Q (Fraction coeff lists, low->high, Q monic at degree d) with
    f = P/Q, degree certified by EXACT agreement on held-out points. `f` is a
    callable returning exact Fractions at rational points (ZeroDivisionError
    at a pole is tolerated: the point is dropped). Returns (None, None, None)
    when no degree <= dmax certifies."""
    if pts is None:
        pts = [F(1) + F(t, 97) for t in range(1, 97)]
    vals = {}
    for a in pts:
        try:
            vals[a] = f(a)
        except ZeroDivisionError:
            pass
    pts = [a for a in pts if a in vals]
    for d in range(1, dmax + 1):
        need = 2 * d + 2
        if len(pts) < need + n_held:
            return None, None, None         # not enough samples
        fit, held = pts[:need], pts[need:need + n_held]
        A, b = [], []
        for a in fit:
            row = [a ** i for i in range(d + 1)]
            row += [-vals[a] * a ** i for i in range(d)]
            A.append(row)
            b.append(vals[a] * a ** d)
        sol = solve_exact(A, b)
        if sol is None:
            continue
        P = sol[:d + 1]
        Q = sol[d + 1:] + [F(1)]
        ok = True
        for a in held:
            qa = sum(Q[i] * a ** i for i in range(len(Q)))
            pa = sum(P[i] * a ** i for i in range(len(P)))
            if qa == 0 or pa / qa != vals[a]:
                ok = False
                break
        if ok:
            return P, Q, d
    return None, None, None


def count_roots(coeffs, lo=1, hi=2):
    """exact count of real roots in the OPEN interval (lo, hi) via sympy
    (Sturm / exact real-root isolation over QQ). Returns (count, info)."""
    import sympy as sp
    a = sp.Symbol('a')
    poly = sum(sp.Rational(c.numerator, c.denominator) * a ** i
               for i, c in enumerate(coeffs))
    if poly == 0:
        return None, "identically zero"
    P = sp.Poly(sp.expand(poly), a)
    total = P.count_roots(lo, hi)
    ends = sum(1 for e in (lo, hi) if P.eval(sp.Integer(e)) == 0)
    return int(total) - ends, f"deg {P.degree()}"


# ---------------- selftest (battery LEG 8) ----------------------------------
def _selftest():
    ok = True

    def check(name, cond):
        nonlocal ok
        print(f"{'PASS' if cond else 'FAIL'} {name}")
        ok = ok and cond

    # S1 positivity via subdivision: (s+u-1)^2 + 1/64 > 0 on the box, but its
    # depth-0 Bernstein patch has a negative interior coefficient (-31/64) —
    # the certificate must come from actual de-Casteljau subdivision.
    p1 = {(0, 0): F(1) + F(1, 64), (1, 0): F(-2), (0, 1): F(-2),
          (2, 0): F(1), (0, 2): F(1), (1, 1): F(2)}
    b1 = to_bernstein(p1, 2, 2)
    check("S1a depth-0 patch not already nonneg (test is real)",
          min(min(r) for r in b1) < 0)
    check("S1b bern_nonneg certifies (s+u-1)^2 + 1/64 on [0,1]^2",
          bern_nonneg(b1, 8))

    # S2 planted MUST-FAIL controls (synthetic-truth law): a certifier that
    # passes a negative polynomial is broken.
    p2 = {(1, 1): F(1), (0, 0): F(-1, 64)}          # s*u - 1/64 < 0 at (0,0)
    check("S2a planted negative-corner poly is REFUSED",
          not bern_nonneg(to_bernstein(p2, 1, 1), 8))
    p3 = dict(p1)
    p3[(0, 0)] = F(1) - F(1, 1024)                  # dips below 0 near s+u=1
    check("S2b planted thin-band negative poly is NOT certified",
          not bern_nonneg(to_bernstein(p3, 2, 2), 6))

    # S3 certified_sup on a known chart: P = s, T = 1 -> sup = 1 exactly.
    c = certified_sup({(1, 0): F(1)}, {(0, 0): F(1)}, F(1, 2), F(2), depth=4)
    check("S3 certified_sup(P=s, T=1) in [1, 1+2^-38]",
          F(1) <= c <= F(1) + F(1, 2 ** 38))

    # S4 C-S bound arithmetic + nonpositive-margin guard.
    check("S4a cs bound (3-1)^2/4 == 1",
          cs_chi2_lower_bound(F(3), F(1), F(4)) == F(1))
    check("S4b nonpositive margin -> 0",
          cs_chi2_lower_bound(F(1), F(2), F(4)) == F(0))

    # S5 rational reconstruction with held-out degree cert + Sturm count:
    # f = (a-5/4)(a-7/4) / ((a+3)(a+2)); numerator has EXACTLY 2 roots in
    # (1,2), denominator none (no pole masking).
    def f(a):
        return ((a - F(5, 4)) * (a - F(7, 4))) / ((a + 3) * (a + 2))

    P, Q, d = reconstruct_rational(f)
    check("S5a reconstruction certifies degree 2", d == 2 and P is not None)
    if P is not None:
        nr, _ = count_roots(P, 1, 2)
        nq, _ = count_roots(Q, 1, 2)
        check("S5b Sturm: numerator roots in (1,2) == 2", nr == 2)
        check("S5c Sturm: denominator has no root in (1,2)", nq == 0)

    if ok:
        print("REGION CERTS PASS")
        return 0
    print("REGION CERTS FAIL")
    return 1


if __name__ == "__main__":
    sys.exit(_selftest())

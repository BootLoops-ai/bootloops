#!/usr/bin/env python3
"""seg_blind.py -- BLIND BUILD 2 of the contiguous-segment author-mixture engine.

STATE DISCIPLINE: written from the SPEC ALONE (blind build). The author did
NOT read seg_v1.py, its derivation note, or any
other derivation/implementation file before or during this build. The route is
deliberately different from a direct segmentation sum: dynamic programming over
(position, segments-used) with Dirichlet-multinomial segment scores -- the
O(n^2 g) DP -- plus an independent brute-force segmentation enumerator used
only inside the selfcheck battery.

MODEL (from SPEC):
  Text = n ordered units, each a length-k count vector over channels.
  M_g: text is g contiguous segments; changepoints unknown, uniform prior over
  the C(n-1, g-1) segmentations. Each segment s has its own p_s ~ Dir(alpha,
  symmetric, k channels), INTEGRATED out; unit counts within a segment pool.
  Segment score (DirMult sequence likelihood; the multinomial coefficient is
  omitted because it depends only on unit-level counts, hence is identical
  across all g and all segmentations and cancels in every posterior):
      S(c) = [prod_ch rising(alpha, c_ch)] / rising(k*alpha, sum(c)),
      rising(x, m) = x (x+1) ... (x+m-1),  rising(x, 0) = 1,  S(empty) = 1.
  Raw sum    T_g = sum_{segmentations} prod_s S(counts_s)
  Evidence   Z_g = T_g / C(n-1, g-1)
  Posteriors: P(g | text) propto prior(g) * Z_g  (prior declared by caller,
  default uniform on g = 1..G); P(boundary after unit t | g, text) by
  forward/backward DP ratio.
  Free-assignment variant (g = 2, n <= 25): assignment a in {0,1}^n, uniform
  prior 2^-n,  Z2_free = 2^-n * sum_a S(pool_0(a)) * S(pool_1(a)).
All arithmetic is exact (fractions.Fraction); no floats in the evidence path.
"""

from fractions import Fraction
from math import comb
from itertools import combinations, permutations
import sys

# ---------------------------------------------------------------- core scores

def rising(x, m):
    """Rising factorial x(x+1)...(x+m-1), exact."""
    out = Fraction(1)
    for i in range(m):
        out *= (x + i)
    return out

def seg_score(counts, alpha):
    """DirMult sequence likelihood S(counts) for one segment, exact Fraction."""
    k = len(counts)
    tot = sum(counts)
    num = Fraction(1)
    for c in counts:
        if c:
            num *= rising(alpha, c)
    return num / rising(k * alpha, tot)

def prefix_sums(units):
    """P[i] = channelwise sum of units[0:i]; P[0] = zero vector."""
    k = len(units[0])
    P = [(0,) * k]
    for u in units:
        P.append(tuple(a + b for a, b in zip(P[-1], u)))
    return P

def all_block_scores(units, alpha):
    """S[i][j] = seg_score of pooled units i..j-1 (0-indexed), via prefix sums."""
    n = len(units)
    P = prefix_sums(units)
    S = [[None] * (n + 1) for _ in range(n + 1)]
    for i in range(n):
        for j in range(i + 1, n + 1):
            S[i][j] = seg_score(tuple(P[j][c] - P[i][c] for c in range(len(P[0]))),
                                alpha)
    return S

# ------------------------------------------------------------------- forward/backward DP

def forward_dp(units, alpha, gmax, S=None):
    """F[j][m] = sum over splittings of units 0..j-1 into m contiguous
    segments of the product of segment scores. F[0][0] = 1."""
    n = len(units)
    if S is None:
        S = all_block_scores(units, alpha)
    F = [[Fraction(0)] * (gmax + 1) for _ in range(n + 1)]
    F[0][0] = Fraction(1)
    for m in range(1, gmax + 1):
        for j in range(m, n + 1):
            acc = Fraction(0)
            for i in range(m - 1, j):
                if F[i][m - 1]:
                    acc += F[i][m - 1] * S[i][j]
            F[j][m] = acc
    return F, S

def backward_dp(units, alpha, gmax, S):
    """B[i][m] = sum over splittings of units i..n-1 into m contiguous
    segments of the product of segment scores. B[n][0] = 1."""
    n = len(units)
    B = [[Fraction(0)] * (gmax + 1) for _ in range(n + 1)]
    B[n][0] = Fraction(1)
    for m in range(1, gmax + 1):
        for i in range(n - m, -1, -1):
            acc = Fraction(0)
            for j in range(i + 1, n + 1):
                if B[j][m - 1]:
                    acc += S[i][j] * B[j][m - 1]
            B[i][m] = acc
    return B

# ------------------------------------------------------------------- evidences

def raw_sums(units, alpha, gmax):
    """[T_0, T_1, ..., T_gmax]; T_g = 0 automatically when g > n."""
    F, _ = forward_dp(units, alpha, gmax)
    n = len(units)
    return [F[n][m] for m in range(gmax + 1)]

def evidence(units, alpha, g):
    """Z_g = T_g / C(n-1, g-1); Fraction(0) when g > n or g < 1."""
    n = len(units)
    if g < 1 or g > n:
        return Fraction(0)
    T = raw_sums(units, alpha, g)
    return T[g] / comb(n - 1, g - 1)

def posterior_g(units, alpha, G, prior=None):
    """P(g | text) for g = 1..G. prior: list of G Fractions (default uniform).
    Returns list of Fractions summing to exactly 1."""
    n = len(units)
    gmax = min(G, n)
    T = raw_sums(units, alpha, gmax)
    if prior is None:
        prior = [Fraction(1, G)] * G
    w = []
    for g in range(1, G + 1):
        if g <= n:
            w.append(prior[g - 1] * T[g] / comb(n - 1, g - 1))
        else:
            w.append(Fraction(0))
    tot = sum(w)
    return [x / tot for x in w]

def boundary_posterior(units, alpha, g):
    """out[t-1] = P(changepoint after unit t | g, text), t = 1..n-1.
    P(bnd at t) = sum_m F[t][m] * B[t][g-m] / T_g   (m = segs used left of t)."""
    n = len(units)
    F, S = forward_dp(units, alpha, g)
    B = backward_dp(units, alpha, g, S)
    Tg = F[n][g]
    out = []
    for t in range(1, n):
        acc = Fraction(0)
        for m in range(1, g):
            if F[t][m] and B[t][g - m]:
                acc += F[t][m] * B[t][g - m]
        out.append(acc / Tg)
    return out

# ------------------------------------------------- independent cross-checkers

def brute_T(units, alpha, g):
    """T_g by direct enumeration of all C(n-1, g-1) segmentations.
    Independent of the DP route; used only in the selfcheck."""
    n = len(units)
    P = prefix_sums(units)
    k = len(units[0])
    tot = Fraction(0)
    for cuts in combinations(range(1, n), g - 1):
        bounds = (0,) + cuts + (n,)
        term = Fraction(1)
        for a, b in zip(bounds, bounds[1:]):
            term *= seg_score(tuple(P[b][c] - P[a][c] for c in range(k)), alpha)
        tot += term
    return tot

def brute_boundary(units, alpha, g):
    """Boundary posterior by direct enumeration; selfcheck only."""
    n = len(units)
    P = prefix_sums(units)
    k = len(units[0])
    tot = Fraction(0)
    hits = [Fraction(0)] * (n - 1)
    for cuts in combinations(range(1, n), g - 1):
        bounds = (0,) + cuts + (n,)
        term = Fraction(1)
        for a, b in zip(bounds, bounds[1:]):
            term *= seg_score(tuple(P[b][c] - P[a][c] for c in range(k)), alpha)
        tot += term
        for t in cuts:
            hits[t - 1] += term
    return [h / tot for h in hits]

def free_g2_raw(units, alpha):
    """sum_{a in {0,1}^n} S(pool_0(a)) S(pool_1(a)); Z2_free = this / 2^n.
    SPEC bound n <= 25 (2^n enumeration)."""
    n = len(units)
    k = len(units[0])
    assert n <= 25, "free-assignment enumeration limited to n <= 25 by SPEC"
    tot = Fraction(0)
    for mask in range(1 << n):
        c0 = [0] * k
        c1 = [0] * k
        for t in range(n):
            tgt = c1 if (mask >> t) & 1 else c0
            for c in range(k):
                tgt[c] += units[t][c]
        tot += seg_score(tuple(c0), alpha) * seg_score(tuple(c1), alpha)
    return tot

def free_g2_prefix_part(units, alpha):
    """Sum of free-assignment terms over the 2(n-1) proper prefix/suffix
    patterns (0^t 1^(n-t)) and (1^t 0^(n-t)), t = 1..n-1. Must equal
    exactly 2 * T_2 (each contiguous 2-split appears twice, by label swap)."""
    n = len(units)
    P = prefix_sums(units)
    k = len(units[0])
    tot = Fraction(0)
    for t in range(1, n):
        left = tuple(P[t][c] for c in range(k))
        right = tuple(P[n][c] - P[t][c] for c in range(k))
        tot += 2 * seg_score(left, alpha) * seg_score(right, alpha)
    return tot

# ------------------------------------------------------------ selfcheck battery

# Deterministic synthetic datasets (no RNG, so the battery is reproducible).
UNITS_PLANT = [  # n=8, k=3; planted changepoint after unit 4 (A: ch0-heavy, B: ch2-heavy)
    (9, 1, 0), (8, 2, 1), (9, 0, 1), (7, 2, 1),
    (1, 1, 8), (0, 2, 9), (1, 0, 9), (2, 1, 8),
]
UNITS_MISC = [  # n=6, k=3; no planted structure
    (3, 1, 2), (0, 4, 1), (2, 2, 2), (5, 0, 0), (1, 3, 1), (2, 0, 3),
]
UNITS_TINY = [(1, 0), (0, 1)]  # n=2, k=2; hand-computable analytic case

def run_selfcheck(outpath):
    lines = []
    results = []

    def check(name, ok, detail=""):
        results.append(ok)
        lines.append("%s  %s%s" % ("PASS" if ok else "FAIL", name,
                                   ("  [" + detail + "]") if detail else ""))

    lines.append("seg_blind.py selfcheck battery -- blind build 2 (DP route)")
    lines.append("state discipline: seg_v1.py / DERIVATION_SEGMENT.md NOT read")
    lines.append("all arithmetic exact Fractions; datasets deterministic")
    lines.append("")

    alphas3 = [Fraction(1, 3), Fraction(1, 2), Fraction(1)]  # k=3: {1/k, 1/2, 1}

    # 1. Analytic tiny case, alpha=1, k=2: hand-computed targets.
    a1 = Fraction(1)
    z1 = evidence(UNITS_TINY, a1, 1)
    z2 = evidence(UNITS_TINY, a1, 2)
    pg = posterior_g(UNITS_TINY, a1, 2)
    check("tiny analytic Z_1 == 1/6", z1 == Fraction(1, 6), str(z1))
    check("tiny analytic Z_2 == 1/4", z2 == Fraction(1, 4), str(z2))
    check("tiny analytic P(g=2) == 3/5", pg[1] == Fraction(3, 5), str(pg[1]))

    for ds_name, units in (("PLANT", UNITS_PLANT), ("MISC", UNITS_MISC)):
        n = len(units)
        for al in alphas3:
            tag = "%s a=%s" % (ds_name, al)
            # 2. Z_1 equals DirMult of the fully pooled counts (single segment).
            pooled = tuple(sum(u[c] for u in units) for c in range(3))
            check("Z_1 == S(pooled)  " + tag,
                  evidence(units, al, 1) == seg_score(pooled, al))
            # 3. DP == brute-force enumeration, g = 1..4, exact equality.
            T = raw_sums(units, al, 4)
            ok = all(T[g] == brute_T(units, al, g) for g in range(1, 5))
            check("DP T_g == brute enumeration, g=1..4  " + tag, ok)
            # 4. Boundary posterior: DP == brute; sums to exactly g-1; in [0,1].
            for g in (2, 3):
                bp = boundary_posterior(units, al, g)
                ok = (bp == brute_boundary(units, al, g)
                      and sum(bp) == g - 1
                      and all(0 <= p <= 1 for p in bp))
                check("boundary posterior g=%d: DP==brute, sum==g-1, in [0,1]  %s"
                      % (g, tag), ok)
            # 5. P(g|text) sums to exactly 1 (uniform prior, G=4).
            check("sum_g P(g|text) == 1  " + tag,
                  sum(posterior_g(units, al, 4)) == 1)
            # 6. Reversal symmetry: Z_g invariant; boundary posterior mirrors.
            rev = list(reversed(units))
            ok = all(evidence(units, al, g) == evidence(rev, al, g)
                     for g in range(1, 5))
            bp = boundary_posterior(units, al, 2)
            ok = ok and (list(reversed(boundary_posterior(rev, al, 2))) == bp)
            check("reversal: Z_g invariant, boundaries mirror  " + tag, ok)
            # 7. Channel permutation invariance (alpha symmetric).
            ok = True
            for perm in permutations(range(3)):
                pu = [tuple(u[p] for p in perm) for u in units]
                ok = ok and raw_sums(pu, al, 3) == raw_sums(units, al, 3)
            check("channel permutation invariance  " + tag, ok)
            # 8. Free-assignment g=2 (2^n enum): prefix/suffix part == 2*T_2.
            check("free-assign prefix patterns == 2*T_2  " + tag,
                  free_g2_prefix_part(units, al) == 2 * raw_sums(units, al, 2)[2])
            # 9. Type discipline: every evidence is a Fraction (no float leak).
            check("all outputs exact Fraction  " + tag,
                  all(isinstance(x, Fraction) for x in
                      raw_sums(units, al, 4) + boundary_posterior(units, al, 2)))

        # 10. g > n gives T_g == 0 (DP handles it without special-casing).
        check("T_g == 0 for g > n  " + ds_name,
              raw_sums(units[:3], Fraction(1), 5)[5] == Fraction(0)
              and evidence(units[:3], Fraction(1), 5) == Fraction(0))

    # 11. Planted changepoint: g=2 boundary posterior peaks at t=4 (all alphas).
    for al in alphas3:
        bp = boundary_posterior(UNITS_PLANT, al, 2)
        check("planted boundary argmax at t=4  a=%s" % al,
              bp.index(max(bp)) == 3, "P(t=4)=%s" % bp[3])

    # 12. Free vs contiguous g=2 on PLANT: report ratio (contiguity fence check;
    # reported, not gated -- the SPEC asks for the cross-check number).
    for al in alphas3:
        n = len(UNITS_PLANT)
        z2c = raw_sums(UNITS_PLANT, al, 2)[2] / comb(n - 1, 1)
        z2f = free_g2_raw(UNITS_PLANT, al) / Fraction(2 ** n)
        lines.append("INFO  free-vs-contig g=2 PLANT a=%s: Z2_contig=%s  Z2_free=%s"
                     % (al, z2c, z2f))
        check("free-assign Z2 is a valid Fraction > 0  a=%s" % al,
              isinstance(z2f, Fraction) and z2f > 0)

    lines.append("")
    verdict = "ALL PASS (%d checks)" % len(results) if all(results) else \
              "FAILURES: %d of %d" % (results.count(False), len(results))
    lines.append("VERDICT: " + verdict)
    with open(outpath, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    return all(results)

if __name__ == "__main__":
    if len(sys.argv) > 1:
        out = sys.argv[1]
    else:  # public port: no machine-local default output path
        raise SystemExit("usage: seg_blind.py <selfcheck-out-path>")
    ok = run_selfcheck(out)
    sys.exit(0 if ok else 1)

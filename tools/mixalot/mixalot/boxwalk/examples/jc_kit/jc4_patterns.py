"""
B1: 15 Jukes-Cantor site-pattern probabilities for the unrooted 4-taxon tree
    topology 12|34 (taxa 1,2 on one side of the internal edge; 3,4 on the other).

Edges: x1..x4 pendant (to taxa 1..4), x5 internal.
JC transition prob along edge e:  P_ii = 1/4 + 3/4 x_e,  P_{i!=j} = 1/4 - 1/4 x_e
where x_e = exp(-4/3 * mu * t_e) in [0,1].

Pattern classes = orbits of S_4 (nucleotide relabelling) on {A,C,G,T}^4
              == set partitions of {1,2,3,4}  (Bell(4) = 15).

Verifies sum_k m_k * p_k(x) == 1 identically.
"""
from fractions import Fraction
import itertools
import sympy as sp

x1, x2, x3, x4, x5 = sp.symbols('x1 x2 x3 x4 x5')
EDGES = [x1, x2, x3, x4, x5]
STATES = (0, 1, 2, 3)  # A,C,G,T


def Pij(i, j, x):
    """JC transition prob along an edge with parameter x."""
    if i == j:
        return sp.Rational(1, 4) + sp.Rational(3, 4) * x
    return sp.Rational(1, 4) - sp.Rational(1, 4) * x


def site_prob(s1, s2, s3, s4):
    """Felsenstein pruning on the 12|34 quartet.
    Internal nodes: u (adj to taxa 1,2 and edge5), v (adj to taxa 3,4 and edge5).
    Root at u with uniform stationary dist 1/4."""
    expr = sp.Integer(0)
    for a in STATES:
        for b in STATES:
            expr += (sp.Rational(1, 4)
                     * Pij(a, s1, x1) * Pij(a, s2, x2)
                     * Pij(a, b, x5)
                     * Pij(b, s3, x3) * Pij(b, s4, x4))
    return sp.expand(expr)


def canonical_partition(s):
    """Map a 4-tuple of states to its set-partition of {0,1,2,3} (sorted blocks)."""
    blocks = {}
    for i, si in enumerate(s):
        blocks.setdefault(si, []).append(i)
    return tuple(sorted(tuple(sorted(b)) for b in blocks.values()))


# Enumerate all 256 patterns, group by set-partition orbit
_orbit_rep = {}       # partition -> representative (s1,s2,s3,s4)
_orbit_mult = {}      # partition -> multiplicity
_orbit_poly = {}      # partition -> polynomial p_k(x)

for s in itertools.product(STATES, repeat=4):
    pi = canonical_partition(s)
    if pi not in _orbit_rep:
        _orbit_rep[pi] = s
        _orbit_mult[pi] = 0
        _orbit_poly[pi] = site_prob(*s)
    _orbit_mult[pi] += 1

# Human-readable labels (e.g. xxxx, xxyy, xyzw) using first-appearance lettering
def label(s):
    seen, out = {}, []
    letters = 'xyzw'
    for si in s:
        if si not in seen:
            seen[si] = letters[len(seen)]
        out.append(seen[si])
    return ''.join(out)

PATTERNS = []  # list of (label, partition, rep, mult, poly)
for pi in sorted(_orbit_rep.keys(), key=lambda p: (-max(len(b) for b in p), len(p), p)):
    rep = _orbit_rep[pi]
    PATTERNS.append((label(rep), pi, rep, _orbit_mult[pi], _orbit_poly[pi]))

assert len(PATTERNS) == 15, f"expected 15 pattern classes, got {len(PATTERNS)}"
assert sum(m for (_, _, _, m, _) in PATTERNS) == 256, "multiplicities must sum to 256"

# Critical check: sum_k m_k p_k == 1 identically
_total = sp.expand(sum(m * p for (_, _, _, m, p) in PATTERNS))
assert _total == 1, f"sum m_k p_k = {_total} != 1"

# Sanity: each p_k is multilinear (degree <=1 in each x_e)
for lab, _, _, _, p in PATTERNS:
    poly = sp.Poly(p, *EDGES)
    for d in poly.degree_list():
        assert d <= 1, f"{lab}: not multilinear, degree_list={poly.degree_list()}"


def pretty_table():
    lines = []
    for i, (lab, pi, rep, m, p) in enumerate(PATTERNS):
        lines.append(f"[{i:2d}] {lab:5s} mult={m:3d}  p = {sp.factor(256*p)}/256")
    return "\n".join(lines)


if __name__ == "__main__":
    print(f"15 JC4 pattern classes (topology 12|34), sum m_k p_k = {_total}")
    print(pretty_table())
    # Spot-check: at x=(1,1,1,1,1) (zero branch lengths), p_xxxx = 1/4, all others 0
    subs1 = {x1: 1, x2: 1, x3: 1, x4: 1, x5: 1}
    for lab, _, _, _, p in PATTERNS:
        v = p.subs(subs1)
        expected = sp.Rational(1, 4) if lab == 'xxxx' else 0
        assert v == expected, f"{lab}: p(1,1,1,1,1)={v}, expected {expected}"
    print("[OK] x=1 spot-check passed (p_xxxx=1/4, rest=0)")
    # Spot-check: at x=(0,...,0) (infinite branches), all p_k = m_k/256? No --
    # at x=0 every P_ij=1/4, so p(s)=1/256 for every s, hence p_k=1/256.
    subs0 = {x1: 0, x2: 0, x3: 0, x4: 0, x5: 0}
    for lab, _, _, _, p in PATTERNS:
        assert p.subs(subs0) == sp.Rational(1, 256), f"{lab} at x=0"
    print("[OK] x=0 spot-check passed (all p_k=1/256)")

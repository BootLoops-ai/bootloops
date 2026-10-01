"""No-exclusion-zone guarantee: structural + statistical.

Statistical target: density (2-x)(1+y) on [0,1]^2 (integral 9/4), two
synthetic special points with Chebyshev eps-tubes (eps = 1/16):
  s1 = (1/4, 1/4): tube mass = (7/32)(5/32) = 35/1024,
       fraction = (35/1024)/(9/4) = 35/2304
  s2 = (3/4, 1/2): tube mass = (5/32)(3/16) = 15/512,
       fraction = (15/512)/(9/4) = 5/384
Accepted samples must populate BOTH tubes at their closed-form fractions
(4 SE) -- the sampler provably does not avoid tube neighborhoods.
"""
import ast
import inspect
import os
import sys
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import rejection
import seeds as S
from ad_density import ADDensity
from cp_bound import cp_on_fctl
from domain import BoxDomain
from rejection import sample

SEEDS_JSON = os.environ.get("TERRIER_SEEDS_JSON", "")
if not SEEDS_JSON:
    raise SystemExit("REFUSE: env TERRIER_SEEDS_JSON unset — must point at "
                     "the seeds JSON file (reference data not included in "
                     "the package)")
EPS = Fraction(1, 16)
S_SYNTH = [(Fraction(1, 4), Fraction(1, 4)), (Fraction(3, 4), Fraction(1, 2))]
FRACS = [Fraction(35, 2304), Fraction(5, 384)]


def _in_tube(x, s):
    return all(abs(xi - si) <= EPS for xi, si in zip(x, s))


def test_structural_no_locus_identifiers():
    """The sampler module cannot even NAME locus data: no identifier in
    its code (args, names, attributes) matches the locus vocabulary."""
    tree = ast.parse(inspect.getsource(rejection))
    ids = {n.id for n in ast.walk(tree) if isinstance(n, ast.Name)}
    ids |= {a.arg for f in ast.walk(tree)
            if isinstance(f, (ast.FunctionDef, ast.AsyncFunctionDef))
            for a in list(f.args.args) + list(f.args.kwonlyargs)}
    ids |= {n.attr for n in ast.walk(tree) if isinstance(n, ast.Attribute)}
    bad = {"S", "S_frozen", "eps0", "eps", "locus", "exclusion", "mask",
           "special_points", "tube"}
    assert ids.isdisjoint(bad), ids & bad


def _R(z):
    x, y = z
    return [[x - 3, 0], [0, -2 - y]]


DENS = ADDensity(2, _R, lambda z: [[1, 0], [0, 1]])


def test_tubes_populated_at_closed_form_rates():
    m = S.load_master(SEEDS_JSON)
    gen = S.bin_generator(m, "SYNTH-NOEXCL", 1)
    pts, c = sample(BoxDomain([0, 0], [1, 1]), DENS, 4, gen, 4000)
    n = len(pts)
    for s, frac in zip(S_SYNTH, FRACS):
        k = sum(1 for x in pts if _in_tube(x, s))
        assert k > 0, "exclusion zone at %r" % (s,)
        mu = float(frac) * n
        se = (n * float(frac) * (1 - float(frac))) ** 0.5
        assert abs(k - mu) <= 4 * se, (s, k, mu, se)


def test_cp_on_fctl_covers_true_tube_fraction():
    """Integration with CP-on-F_ctl: the upper 95% CP bound on the
    sampled tube fraction sits at/above the closed-form fraction (this is
    a fixed-seed draw, so the check is deterministic)."""
    m = S.load_master(SEEDS_JSON)
    gen = S.bin_generator(m, "SYNTH-NOEXCL", 1)
    pts, _ = sample(BoxDomain([0, 0], [1, 1]), DENS, 4, gen, 4000)
    k = sum(1 for x in pts if _in_tube(x, S_SYNTH[0]))
    p_u = cp_on_fctl(k, len(pts))
    assert p_u >= FRACS[0], (float(p_u), float(FRACS[0]))

"""dist_cert.zchart -- G4/G5 certified distance in the pinned algebraic z chart.

Definition: dist_cert(v, s) = |z*(v) - z_s| in BALL arithmetic, z*(v) a
certified transported ball. Conifold-locus rule: conifolds are measured in
the SAME z chart, never WP. eps0 = 1e-3 (z chart); validity gate:
certified ball radius <= eps0/10 = 1e-4, else the candidate is
N_failed(reason=ball-blowup). Duality invariance: dist_cert on the
QUOTIENT = min over orbit representatives. The constants below are the
pinned strata only (conifolds; the G5 rank-2 class point), exposed for
tests -- a census run loads its full frozen special set S_frozen from a
JSON list via load_loci().
NO census data is touched by this module.
"""
from fractions import Fraction
import json

from ball import CBall, RIv, dist_iv

F = Fraction

EPS0 = {"G4": F(1, 1000), "G5": F(1, 1000)}  # pinned per geometry, a NUMBER


class Locus:
    """One special point s: exact rational (rad 0) or algebraic (minpoly +
    isolating ball certified to contain exactly that root)."""

    def __init__(self, name, ball, minpoly=None):
        self.name, self.ball, self.minpoly = name, ball, minpoly

    def __repr__(self):
        return f"Locus({self.name})"


# Pinned strata (cited special points of the two one-parameter families):
# G4 = AESZ34: conifolds {1/25, 1/9, 1} in z chart.
PINNED_LOCI_G4 = [
    Locus("conifold_1/25", CBall(F(1, 25))),
    Locus("conifold_1/9", CBall(F(1, 9))),
    Locus("conifold_1", CBall(F(1))),
]

# G5 = mirror Reye AESZ22/118: rank-2 class point z = -1 (CITED
# arXiv:2410.07107); conifolds 1/32 and (11 +- 5 sqrt(2))/2 (minpoly
# 4z^2-44z+71); apparent singularity 7/4 EXCLUDED by convention pin.
_SQRT2_LO = F(199032864766430, 2**47)   # certified: lo^2 < 2 < hi^2 (tested)
_SQRT2_HI = _SQRT2_LO + F(1, 2**47)
_C_PLUS = (11 + 5 * (_SQRT2_LO + _SQRT2_HI) / 2) / 2
_C_MINUS = (11 - 5 * (_SQRT2_LO + _SQRT2_HI) / 2) / 2
_C_RAD = 5 * (_SQRT2_HI - _SQRT2_LO) / 4
PINNED_LOCI_G5 = [
    Locus("rank2_T3_z=-1", CBall(F(-1))),
    Locus("conifold_1/32", CBall(F(1, 32))),
    Locus("conifold_(11+5r2)/2", CBall(_C_PLUS, 0, _C_RAD), (71, -44, 4)),
    Locus("conifold_(11-5r2)/2", CBall(_C_MINUS, 0, _C_RAD), (71, -44, 4)),
]


def load_loci(path):
    """Load S_frozen from a JSON list. Entries:
    {name, re:[num,den], im:[num,den], rad:[num,den], minpoly optional}."""
    with open(path) as fh:
        rows = json.load(fh)
    out = []
    for r in rows:
        b = CBall(F(*r["re"]), F(*r.get("im", [0, 1])), F(*r.get("rad", [0, 1])))
        out.append(Locus(r["name"], b, tuple(r["minpoly"]) if "minpoly" in r
                         else None))
    return out


def dist_to_locus(vball, locus, same_point=False, bits=64):
    """Certified RIv for |z*(v) - z_s|. same_point=True asserts the candidate
    IS the locus point (exact algebraic identity, e.g. sieve output equals a
    frozen algebraic point) => exact 0 (the algebraic-locus exactness path)."""
    if same_point:
        return RIv(0, 0)
    return dist_iv(vball, locus.ball, bits)


def dist_cert_quotient(vballs, loci, bits=64):
    """Duality-quotient distance: min over ORBIT REPRESENTATIVES
    of the candidate x min over the duality-closed S_frozen list.
    vballs: list of CBall orbit representatives. Returns (RIv, argmin locus)."""
    best, arg = None, None
    for vb in vballs:
        for s in loci:
            d = dist_to_locus(vb, s, bits=bits)
            if best is None or d.hi < best.hi:
                best, arg = d, s
    lo = min(min(dist_to_locus(vb, s, bits=bits).lo for s in loci)
             for vb in vballs)
    return RIv(lo, best.hi), arg


def classify(vballs, geometry, loci=None, bits=64):
    """Per-candidate verdict for geometry "G4"/"G5". Returns
    (verdict, RIv|reason, locus).
    ON  : certified min-dist upper bound <= eps0
    OFF : certified min-dist lower bound  > eps0
    FAILED(ball-blowup): any orbit-rep radius > eps0/10 (validity gate FIRST)
    FAILED(indeterminate): enclosure straddles eps0 (counts into N_failed)."""
    eps0 = EPS0[geometry]
    if loci is None:
        loci = PINNED_LOCI_G4 if geometry == "G4" else PINNED_LOCI_G5
    if any(vb.rad > eps0 / 10 for vb in vballs):
        return "FAILED", "ball-blowup", None
    d, s = dist_cert_quotient(vballs, loci, bits)
    if d.hi <= eps0:
        return "ON", d, s
    if d.lo > eps0:
        return "OFF", d, s
    return "FAILED", "indeterminate", s

r"""membound.assemble — Gamma-prefactor Laurent builders, closed-form Mellin
slots, IBP application, and the BoundaryVector output contract.

Layer basis: a memory-region boundary VECTOR at omega->0 is a list of layers
|omega|^{-k eps}; k = (sum of retarded-factor eps-weights) + (#K-lines)
[e.g. k=7 for the 2SF planar memory sector — matching the 2^{7eps}
prefactors of arXiv:2601.16256].  Each layer carries named slots; a slot is either
  gamma-closed : exact Gamma/hypergeometric Laurent (single-frequency cores
                 close via the Mellin law
                 Int_0^inf w^{s-1} K_eps(w) dw = 2^{s-2} G((s+eps)/2) G((s-eps)/2))
  numeric      : eps-layer moments from membound.core quadrature
  pslq-slot    : numeric + recognition hook (tools/lockpick/pslq_gate.py —
                 NEVER raw mp.pslq; two-precision + planted controls)
"""
import json
import mpmath as mp
from .laurent import L


# --- closed-form single-frequency slot ---------------------------------------
def mellin_K(s0, s1, N=6):
    """Laurent of Int_0^inf w^{(s0 + s1*eps) - 1} |w|^{-eps} K_eps(w) dw.

    = 2^{s-2} Gamma((s+eps)/2) Gamma((s-eps)/2) with s = s0 + (s1-1)*eps.
    Analytic for s0 > 0 (pole layers must be factored by the caller when
    s0 + (s1-1)*eps hits the Gamma poles; this slot covers the generic case).
    """
    def fn(e):
        s = s0 + (s1 - 1) * e
        return mp.mpf(2) ** (s - 2) * mp.gamma((s + e) / 2) * mp.gamma((s - e) / 2)
    return L.taylor(fn, N)


# --- c_M-sector exact prefactors (arXiv:2601.16256 eqs. 23-24, poles factored)
NTAY = 6


def prefactor_I1():
    """-2^{7e-1} G((2e-1)/2)^2 G(4e-1) / [G(2-5e) (4pi)^{5-4e}], 1/eps factored."""
    def num_an(e):
        return (-mp.mpf(2) ** (7 * e - 1) * mp.gamma((2 * e - 1) / 2) ** 2
                * mp.gamma(1 + 4 * e) / ((4 * e - 1) * 4)
                / (mp.gamma(2 - 5 * e) * (4 * mp.pi) ** (5 - 4 * e)))
    return L(-1, [1]) * L.taylor(num_an, NTAY)


def prefactor_I3():
    """+2^{7e-3} G((2e+1)/2)^2 G(4e-2) / [G(3-5e) (4pi)^{5-4e}], 1/eps factored."""
    def num_an(e):
        return (mp.mpf(2) ** (7 * e - 3) * mp.gamma((2 * e + 1) / 2) ** 2
                * mp.gamma(1 + 4 * e) / ((4 * e - 1) * (4 * e - 2) * 4)
                / (mp.gamma(3 - 5 * e) * (4 * mp.pi) ** (5 - 4 * e)))
    return L(-1, [1]) * L.taylor(num_an, NTAY)


def ibp_I2(I1, I3):
    """I2^(M) from eq:IBP of 2601.16256 supplement (1/eps^3 factored)."""
    def c3_hat(e):
        num = 28 * (5 * e - 2) * (14 * e - 9) * (1088 * e ** 4 - 898 * e ** 3
                                                 + 275 * e ** 2 - 38 * e + 2)
        den = 3 * (96 * e ** 3 - 116 * e ** 2 + 43 * e - 5)
        return num / den

    def c1_hat(e):
        num = (1 - 2 * e) ** 2 * (121024 * e ** 6 - 237340 * e ** 5
                                  + 181552 * e ** 4 - 70363 * e ** 3
                                  + 14767 * e ** 2 - 1612 * e + 72)
        den = 6 * (96 * e ** 3 - 116 * e ** 2 + 43 * e - 5)
        return num / den

    c3 = L(-3, [1]) * L.taylor(lambda e: c3_hat(mp.mpf(e)), NTAY)
    c1 = L(-3, [1]) * L.taylor(lambda e: c1_hat(mp.mpf(e)), NTAY)
    return c3 * I3 + c1 * I1


def core_to_L(j):
    """eps-layer moments {0: j0, 1: j1, 2: j2} -> Laurent L(0, [j0, j1, j2])."""
    return L(0, [j[0], j[1], j[2]])


# --- BoundaryVector contract --------------------------------------------------
def boundary_vector(family, region, sector, prescription, layers, meta=None):
    """Assemble the JSON-serializable boundary-vector record.

    layers: list of dicts {"k": int, "slots": [ {name, closure, eps_series
            (dict str(order)->str value), digits (measured; None = unknown -
            measuring), notes} ]}
    """
    return {
        "contract": "membound.boundary_vector.v1",
        "family": family,
        "region": region,
        "sector": sector,
        "prescription": prescription,
        "layer_basis": "|omega|^{-k*eps}",
        "layers": layers,
        "meta": meta or {},
    }


def emit(path, bv):
    with open(path, "w") as fh:
        json.dump(bv, fh, indent=1, default=str)

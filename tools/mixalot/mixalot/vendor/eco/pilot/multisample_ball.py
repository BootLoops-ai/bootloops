"""Certified ball-mode multi-sample Etienne likelihood (Panama-scale).

Same structure as multisample.multisample_P but in arb ball arithmetic with
fmpz/arb product trees; positivity (theta>0, I_i>0) keeps enclosures tight.
M(D,A,I) is parameter-dependent (I_i inside), so the species-poly build is
per-evaluation — still quasi-linear (the across-species product tree
dominates and is rebuilt each call).

Transcription certified vs model DP (run_multisample_controls.py, ALL PASS):
per-sample (I_i)_{J_i} normalization.
"""
from collections import Counter
from math import factorial

from flint import arb, arb_poly, fmpz_poly, ctx

_rf_cache = {}


def rf_poly_int(n):
    """fmpz_poly z(z+1)...(z+n-1) (product tree)."""
    if n in _rf_cache:
        return _rf_cache[n]
    if n == 0:
        p = fmpz_poly([1])
    else:
        fs = [fmpz_poly([k, 1]) for k in range(n)]
        while len(fs) > 1:
            fs = [fs[i] * fs[i + 1] for i in range(0, len(fs) - 1, 2)] \
                 + ([fs[-1]] if len(fs) % 2 else [])
        p = fs[0]
    _rf_cache[n] = p
    return p


def _tree(polys):
    while len(polys) > 1:
        polys = [polys[i] * polys[i + 1] for i in range(0, len(polys) - 1, 2)] \
                + ([polys[-1]] if len(polys) % 2 else [])
    return polys[0]


def species_poly_ball(nvec, Ivec_arb, prec):
    """arb_poly S_k(z): across-sample convolution + (b-1)! Borel twist.

    Coefficients scaled: we fold I_i^a into the polynomial via z -> composition
    with scaled coefficient lists (coeff_a * I_i^a).
    """
    ctx.prec = prec
    u = arb_poly([arb(1)])
    for n, I in zip(nvec, Ivec_arb):
        if n == 0:
            continue
        rf = rf_poly_int(n)
        Ia = arb(1)
        coeffs = []
        for a in range(n + 1):
            coeffs.append(arb(int(rf[a])) * Ia)
            Ia *= I
        u = u * arb_poly(coeffs)
    uc = [u[b] for b in range(u.length())]
    return arb_poly([uc[b] * arb(factorial(b - 1)) if b >= 1 else arb(0)
                     for b in range(len(uc))])


def logP_multisample_ball(Dmat, theta_str, I_strs, prec):
    """Certified log P[D | theta, I_1..I_N] (arb ball).

    Dmat: list of per-species tuples (n_i1..n_iN) — zeros allowed for absent.
    """
    ctx.prec = prec
    th = arb(theta_str)
    Ivec = [arb(s) for s in I_strs]
    N = len(Ivec)
    S = len(Dmat)
    Js = [sum(row[i] for row in Dmat) for i in range(N)]

    # M-generating polynomial: product tree over species
    polys = [species_poly_ball(row, Ivec, prec) for row in Dmat]
    M = _tree(polys)

    # sum_A M_A * theta^S / (theta)_A, running Pochhammer (positive terms)
    tot = arb(0)
    t = arb(1)
    for k in range(S):
        t /= th + k
    for A in range(S, M.length()):
        tot += M[A] * t
        t /= th + A
    # note: t entering loop is 1/(theta)_S; inside we divide progressively

    lp = tot.log() + S * th.log()
    # prefactor logs
    for cnt in Counter(tuple(r) for r in Dmat).values():
        lp -= arb(cnt + 1).lgamma()
    for i in range(N):
        lp += arb(Js[i] + 1).lgamma()
        lp -= (Ivec[i] + Js[i]).lgamma() - Ivec[i].lgamma()
    for row in Dmat:
        for n in row:
            if n > 1:
                lp -= arb(n + 1).lgamma()
    return lp

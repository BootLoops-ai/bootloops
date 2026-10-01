# FFP pilot driver. Blind a via E-law (a = M1 - E, E = p^3-p^2+p-3); k via J'-quantization;
# N2-side via class Weil-windows in r (step offset rel. T1 class); FIRE := exists r in window
# with integral b and split quartic. T1 anchor used ONLY to place D2 reference & validate.
import json, math, sys
import numpy as np
from field import inv_table_fp
from hv_count import torus_counts, twisted12_counts
from ext_field import mu_row
from anchors import ap_14a, load_nf
from weil import weil_circle
from quartic import split_12
from modsqrt import reduce_quad

def b_interval(a, p):
    """Closed-form Weil window in b for given a (region: both u-roots in [-2,2])."""
    if a * a > 16 * p ** 3:
        return None
    hi = (a * a + 8 * p ** 3) // (4 * p)                      # disc >= 0
    lo = -2 * p * p + math.isqrt(4 * a * a * p) // p          # (2p^2+b)sqrt(p) >= 2|a|
    while not weil_circle(a, lo, p):
        lo += 1
    return lo, hi

def compute(p):
    Kp, N1 = torus_counts(p, 1)
    W12 = twisted12_counts(p)
    W5 = mu_row(p, 5, 1)
    K2, N2 = torus_counts(p, 2)
    invp = inv_table_fp(p)
    E = p ** 3 - p * p + p - 3
    sing = {pow(c, p - 2, p) for c in (1, 9, 25)}
    smooth = [f for f in range(1, p) if f not in sing]
    psi = {f: int(invp[f]) for f in smooth}
    n1 = {f: int(N1[psi[f]]) for f in smooth}
    w12 = {f: int(W12[psi[f]]) for f in smooth}
    n2 = {f: int(N2[K2.embed_fp(psi[f])]) for f in smooth}
    m1r = {f: n1[f] + 4 * int(W5[psi[f]]) for f in smooth}
    assert all(v % 5 == 0 for v in m1r.values())
    a = {f: m1r[f] // 5 - E for f in smooth}
    Jp = {f: n1[f] + w12[f] for f in smooth}
    base = max(Jp.values())
    k = {}
    for f in smooth:
        d = base - Jp[f]
        assert d % (12 * p) == 0
        k[f] = d // (12 * p)
    # law checks (validates E-law + k at every phi)
    f0 = smooth[0]
    C1 = n1[f0] - a[f0] + 24 * p * k[f0]
    Cw = w12[f0] + a[f0] - 12 * p * k[f0]
    for f in smooth:
        assert n1[f] == C1 + a[f] - 24 * p * k[f], (p, f, 'N1law')
        assert w12[f] == Cw - a[f] + 12 * p * k[f], (p, f, 'W12law')
    # T1 anchor for the N2 reference constant
    f1 = (-pow(7, p - 2, p)) % p
    alpha1, beta1 = ap_14a(p), load_nf('nf_14_4_a_a.json')[p]
    a1 = -(p * alpha1 + beta1)
    assert a[f1] == a1, (p, 'T1 alpha/beta anchor vs blind a', a[f1], a1)
    b1 = 2 * p * p + alpha1 * beta1
    C2p = n2[f1] + a1 * a1 - 2 * b1 * p           # = D2 + shat(class(T1)) p^2
    # class r-windows (intersection of per-phi closed-form b windows)
    wins = {}
    for kv in sorted(set(k.values())):
        lo, hi = -10 ** 9, 10 ** 9
        for f in smooth:
            if k[f] != kv:
                continue
            iv = b_interval(a[f], p)
            assert iv is not None, (p, f, 'a out of Weil range')
            # b = (a^2 - s2)/(2p), s2 = C2p + r p^2 - n2[f]  =>  r = (a^2-2pb - C2p + n2)/p^2
            rlo = -(-(a[f] * a[f] - 2 * p * iv[1] - C2p + n2[f]) // (p * p))
            rhi = (a[f] * a[f] - 2 * p * iv[0] - C2p + n2[f]) // (p * p)
            lo, hi = max(lo, rlo), min(hi, rhi)
        assert lo <= hi, (p, kv, 'empty class window')
        wins[kv] = (lo, hi)
    return dict(p=p, E=E, a=a, k=k, n2=n2, C2p=C2p, wins=wins, smooth=smooth, sing=sing,
                f1=f1, alpha1=alpha1, beta1=beta1)

def test_point(ctx, f):
    """FIRE semantics for residue f at prime ctx['p']."""
    p, a, k, n2, C2p, wins = (ctx[x] for x in ('p', 'a', 'k', 'n2', 'C2p', 'wins'))
    if f not in a:
        return dict(status='singular-reduction')
    av, kv = a[f], k[f]
    lo, hi = wins[kv]
    hits = []
    for r in range(lo, hi + 1):
        s2 = C2p + r * p * p - n2[f]
        num = av * av - s2
        if num % (2 * p):
            continue
        b = num // (2 * p)
        if not weil_circle(av, b, p):
            continue
        for al, be in split_12(av, b, p):
            hits.append(dict(r=r, b=b, alpha=al, beta=be))
    return dict(status='ok', a=av, kclass=kv, window=[lo, hi], hits=hits, fired=bool(hits))

# eichler.genus2 — invariant harness: recognition-free candidate arbitration.
# True side: the oracle json's leg1 Rosenhain invariants (370d, VALID numerics) -> abs invariants.
# Candidate side: exact sextic coefficients -> exact specialization (A,B,C;w)=(2,3,5;6)
# -> 350d roots -> abs invariants (+ Richelot orbit) -> compare. No PSLQ anywhere.
import json, itertools
from mpmath import mp, mpf, mpc, polyroots, nstr
import sympy as sp

mp.dps = 340
import os
# oracle json: EICHLER_GENUS2_TRUE_JSON (package name) or G2KIT_TRUE_JSON
# (the original name, still honored) — either one, the first that is set wins.
_TRUE_JSON = (os.environ.get('EICHLER_GENUS2_TRUE_JSON')
              or os.environ.get('G2KIT_TRUE_JSON'))
if not _TRUE_JSON:
    raise SystemExit(
        "G2KIT_TRUE_JSON (alias EICHLER_GENUS2_TRUE_JSON) is not set. The "
        "harness arbitrates candidate curves against a high-dps numeric oracle "
        "json (reference oracle not distributed): "
        "{'leg1': {'rosenhain': [three ~370d complex strings]}}. "
        "Point G2KIT_TRUE_JSON at your own oracle file.")
TRUE = json.load(open(_TRUE_JSON))

def pc(s):
    s = s.replace(' ', '').strip('()')
    k = None
    for p in range(1, len(s)):
        if s[p] in '+-' and s[p-1] not in 'eE': k = p
    if s.endswith('j'):
        return mpc(mpf(s[:k]), mpf(s[k:-1])) if k else mpc(0, mpf(s[:-1]))
    return mpc(mpf(s), 0)

lamT, muT, nuT = [pc(t).real for t in TRUE['leg1']['rosenhain']]

def mob(t):  # (2t-1)/(t+1): {0,1,oo}->{-1,1/2,2}
    return (2*t-1)/(t+1)

def invariants_from_roots(r):
    """absolute (j1,j2,j3) from 6 finite branch points; Moebius/twist-invariant."""
    idx = range(6)
    def matchings(elems):
        if not elems: yield []
        else:
            a = elems[0]
            for k in range(1, len(elems)):
                rest = elems[1:k] + elems[k+1:]
                for m in matchings(rest):
                    yield [(a, elems[k])] + m
    A = mpc(0)
    for m in matchings(list(idx)):
        t = mpc(1)
        for (i, j) in m: t *= (r[i]-r[j])**2
        A += t
    trips = [(S, tuple(i for i in idx if i not in S)) for S in itertools.combinations(idx, 3) if 0 in S]
    def tri(T):
        (i, j, k) = T
        return (r[i]-r[j])*(r[j]-r[k])*(r[k]-r[i])
    B = sum(tri(S)**2*tri(Sc)**2 for (S, Sc) in trips)
    C = mpc(0)
    for (S, Sc) in trips:
        for perm in itertools.permutations(Sc):
            t = tri(S)**2*tri(Sc)**2
            for (i, j) in zip(S, perm): t *= (r[i]-r[j])**2
            C += t
    D = mpc(1)
    for (i, j) in itertools.combinations(idx, 2): D *= (r[i]-r[j])**2
    return (A**5/D, A**3*B/D, A**2*C/D)

def rel_dist(J1, J2):
    return max(abs(a-b)/max(abs(a), abs(b), mpf(1)) for a, b in zip(J1, J2))

# ---- true invariants (two Moebius frames: self-test) ----
rootsT_frame1 = [mpc(-1), mpc(mpf(1)/2), mpc(2), mob(lamT), mob(muT), mob(nuT)]
def mob2(t): return (t-3)/(t+2)     # different frame, {0,1,oo}->{-3/2,-2/3,1}
rootsT_frame2 = [mpc(-3)/2, mpc(-2)/3, mpc(1), mob2(lamT), mob2(muT), mob2(nuT)]
JT1 = invariants_from_roots(rootsT_frame1)
JT2 = invariants_from_roots(rootsT_frame2)
print('SELF-TEST (two frames):', nstr(rel_dist(JT1, JT2), 4))
# planted-mismatch control
rootsP = [mpc(-1), mpc(mpf(1)/2), mpc(2), mob(lamT + mpf(10)**-30), mob(muT), mob(nuT)]
JP = invariants_from_roots(rootsP)
print('MISMATCH CONTROL (1e-30 perturb):', nstr(rel_dist(JT1, JP), 4), '(must be ~1e-25, i.e. DETECTED)')

# ---- candidate machinery ----
def roots_of_sextic(coeffs_desc):
    """coeffs_desc: exact sympy Rationals, descending. If degree 5, add root at 'oo'
    handled by Moebius u = 1/(t-shift) trick: simpler — treat deg-5 as 6th root at oo:
    map all roots r -> mob(r), plus mob(oo)=2."""
    cs = [mpf(sp.Rational(c).p)/mpf(sp.Rational(c).q) for c in coeffs_desc]
    rts = polyroots(cs, maxsteps=400, extraprec=400)
    return list(rts)

def candidate_invariants(poly_sigma, sigma, is_quintic):
    p = sp.Poly(sp.expand(poly_sigma), sigma)
    rts = roots_of_sextic(p.all_coeffs())
    rr = [mob(r) for r in rts]
    if is_quintic: rr.append(mpc(2))     # mob(oo)
    assert len(rr) == 6
    return invariants_from_roots(rr), rr

def richelot_orbit_invs(rr):
    """15 Richelot partners from quadratic-factor partitions of the 6 roots."""
    out = []
    idx = list(range(6))
    seen_parts = set()
    for pairing in _pairings(idx):
        key = tuple(sorted(tuple(sorted(p)) for p in pairing))
        if key in seen_parts: continue
        seen_parts.add(key)
        Gs = []
        ok = True
        for (i, j) in pairing:
            # G(t) = (t - r_i)(t - r_j)
            Gs.append((mpf(1), -(rr[i]+rr[j]), rr[i]*rr[j]))
        # Richelot: H_i = [G_j, G_k] = G_j' G_k - G_j G_k' (deg <=2), new curve y^2 = H1 H2 H3 / delta
        def bracket(Ga, Gb):
            a2, a1, a0 = Ga; b2, b1, b0 = Gb
            # G' coefficients: (2a2, a1)
            # [Ga,Gb] = Ga' Gb - Ga Gb':
            # do polynomial arithmetic explicitly
            Gap = [2*a2, a1]; Gbp = [2*b2, b1]
            def pmul(u, v):
                out = [mpc(0)]*(len(u)+len(v)-1)
                for i, ui in enumerate(u):
                    for j, vj in enumerate(v):
                        out[i+j] += ui*vj
                return out
            t1 = pmul(Gap, [b2, b1, b0])
            t2 = pmul([a2, a1, a0], Gbp)
            L = max(len(t1), len(t2))
            t1 = [mpc(0)]*(L-len(t1)) + t1
            t2 = [mpc(0)]*(L-len(t2)) + t2
            return [u-v for u, v in zip(t1, t2)]
        H = [bracket(Gs[a], Gs[b]) for (a, b) in ((1, 2), (2, 0), (0, 1))]
        # product H1*H2*H3 -> degree 6 poly; roots -> invariants
        def pmul(u, v):
            out = [mpc(0)]*(len(u)+len(v)-1)
            for i, ui in enumerate(u):
                for j, vj in enumerate(v):
                    out[i+j] += ui*vj
            return out
        prod = pmul(pmul(H[0], H[1]), H[2])
        if abs(prod[0]) < mpf(10)**(-100):
            prod = prod[1:]
        try:
            rts = polyroots(prod, maxsteps=400, extraprec=400)
        except Exception:
            continue
        rr2 = [mob(r) for r in rts]
        while len(rr2) < 6: rr2.append(mpc(2))
        try:
            out.append(invariants_from_roots(rr2))
        except Exception:
            continue
    return out

def _pairings(elems):
    if not elems:
        yield []
        return
    a = elems[0]
    for k in range(1, len(elems)):
        b = elems[k]
        rest = elems[1:k] + elems[k+1:]
        for rest_pairing in _pairings(rest):
            yield [(a, b)] + rest_pairing

# ---- exact candidates at (A,B,C;w) = (2,3,5;6) ----
sigma = sp.Symbol('sigma')
Aq, Bq, Cq = sp.Integer(2), sp.Integer(3), sp.Integer(5)
wq = sp.Integer(6)
alpha_p_beta = 8*(Aq+Bq)            # α+β
alpha_t_beta = 16*(Aq-Bq)**2        # αβ
quadAB = sigma**2 - alpha_p_beta*sigma + alpha_t_beta
R = (4*Cq - 4*wq**2 - sigma)**2 - 16*wq**2*sigma

candidates = {
  'C1_prym: y2 = sigma*(sigma2-quadAB)*R  [deg5+oo]': (sigma*quadAB*R, True),
  'C3_naive_s-plane': (None, None),   # built from roots directly below
}
results = {}
for name, (poly, isq) in candidates.items():
    if poly is None: continue
    try:
        J, rr = candidate_invariants(poly, sigma, isq)
        d0 = rel_dist(JT1, J)
        dR = min([d0] + [rel_dist(JT1, Jr) for Jr in richelot_orbit_invs(rr)])
        results[name] = (nstr(d0, 4), nstr(dR, 4))
        print(f'{name}: direct {nstr(d0,4)}  best-in-Richelot-orbit {nstr(dR,4)}', flush=True)
    except Exception as e:
        results[name] = ('ERR', str(e)[:80]); print(name, 'ERR', str(e)[:80])

# naive s-plane sextic: roots {±(2a-2b), ±(2a+2b), 2w±2c} with a=sqrt2 etc — numeric roots
from mpmath import sqrt as msqrt
a_, b_, c_ = msqrt(2), msqrt(3), msqrt(5)
naive_roots = [2*a_-2*b_, -(2*a_-2*b_), 2*a_+2*b_, -(2*a_+2*b_), 2*mpf(6)-2*c_, 2*mpf(6)+2*c_]
rrN = [mob(mpc(r)) for r in naive_roots]
JN = invariants_from_roots(rrN)
d0 = rel_dist(JT1, JN)
dR = min([d0] + [rel_dist(JT1, Jr) for Jr in richelot_orbit_invs(rrN)])
print(f'C3_naive_s-plane: direct {nstr(d0,4)}  best-in-Richelot-orbit {nstr(dR,4)}', flush=True)
results['C3_naive_s-plane'] = (nstr(d0, 4), nstr(dR, 4))
json.dump({'selftest': nstr(rel_dist(JT1, JT2), 6), 'mismatch_ctrl': nstr(rel_dist(JT1, JP), 6),
           'true_invariants': [nstr(t, 60) for t in JT1], 'results': results},
          open('pilot_harness.json', 'w'), indent=1)
print('wrote pilot_harness.json')

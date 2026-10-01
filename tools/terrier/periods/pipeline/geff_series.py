#!/usr/bin/env python3
"""
geff_series.py — h_eff-parameter effective Gamma-series (general pipeline).

Effective Gamma-series for any h (the h=2 DKMM example is the regression case):
  A(nu) = prod_num Gamma(1+v.nu)^m / prod_den Gamma(1+v.nu)^m,  nu = n + rho,
computed as exact truncated rho-jets (total degree <= 3: CY3 towers stop at
log^3 regardless of h_eff) over Q[gamma, zeta2, zeta3]; gamma-cancellation
(CY condition) asserted per term.  Provides:
  * GammaJets(card, nmax)      — A(n) jets on a multi-index box
  * cfrac(card, n)             — exact rational w0 coefficient (0 out of cone)
  * w0_curve / jet_tower       — restriction to the monomial curve
                                 z~_a = s^{nu_a} (tilde = sigma-twisted frame)
  * gv_extract                 — h-dim HKTY double-log GV extraction (log-free
                                 identity P_a = k_abc/2 [S_bc/w0-(S_b/w0)(S_c/w0)])
                                 + sigma/eps sign pin against pinned GV leaders
  * mirror_curve               — restricted mirror map q_s(s), s(q_s), N_m peel
  * eff_route                  — SIMPLICIALITY CHECK FIRST (exact extreme-ray
                                 reduction + facet H-rep of the card's
                                 effective cone); routes the GV layer
NON-SIMPLICIAL EFFECTIVE CONES: when the CY-effective cone has MORE
extreme rays than h_eff, no effective-generator basis exists and the h-param
double-log identity AND the mirror-map instanton peel of gv_extract fail in
ANY containing basis (proved: exact kappa-fit rank-full inconsistent; the
plain semigroup-filtered lattice sum was ALSO measured to fail, already at
total degree 1 — the CY-phase double-log data needs sectors the h-param jets
do not carry).  Route: eff_route(card) decides FIRST from card
eff_rays; simplicial cards (incl. every card without eff_rays) take the
simplicial path byte-identically; non-simplicial cards FAIL-CLOSE gv_extract
and route pin_sign_frame through the provenance-gated card gv_table, gated by
cone containment + integrality + leader pins here and by the racetrack peel
(mirror_curve) downstream — exactly the layers measured to stay EXACT (w0,
single-log tower, mirror map, naive m-graded sums, N_m peel).  Also enforced
here: deg(-K.ell)>=0 basis rule surfaced by name (GammaJets/cfrac);
negative-arg factorial law in cfrac (dominated num-Gamma poles -> exact 0,
never a silent negative-index lookup); gkz_box_degree_factors (the GKZ
numerator-column shift) for card derivers.
GATE (run_pipe G1): from the DKMM card this module must reproduce the
reference restricted series, jet tower and GV table (TERRIER_KKLT_BANK
restrict/) EXACTLY.  Single-thread, exact fractions throughout.
"""
from fractions import Fraction as Fr
from math import factorial, gcd
from itertools import product, combinations

TR = 3  # rho-jet truncation (CY3: up to triple logs)

# ------------- symbol coefficients: {(g,e2,e3): Fr} = gamma^g zeta2^e2 zeta3^e3
def sadd(a, b):
    out = dict(a)
    for k, v in b.items():
        w = out.get(k, Fr(0)) + v
        if w: out[k] = w
        elif k in out: del out[k]
    return out

def smul(a, b):
    out = {}
    for ka, va in a.items():
        for kb, vb in b.items():
            k = (ka[0]+kb[0], ka[1]+kb[1], ka[2]+kb[2])
            w = out.get(k, Fr(0)) + va*vb
            if w: out[k] = w
            elif k in out: del out[k]
    return out

def sscale(a, r):
    return {k: v*r for k, v in a.items() if v*r != 0}

S_ONE = {(0,0,0): Fr(1)}
S_G   = {(1,0,0): Fr(1)}
S_Z2  = {(0,1,0): Fr(1)}
S_Z3  = {(0,0,1): Fr(1)}

def s_rat(a):
    for k in a:
        assert k[0] == 0, f"gamma survived: {a}"
    return a.get((0,0,0), Fr(0))

# ------------- rho-jets: {alpha (h-tuple): symdict}, |alpha| <= TR
def jmul(A, B):
    out = {}
    for ka, ca in A.items():
        for kb, cb in B.items():
            if sum(ka) + sum(kb) <= TR:
                key = tuple(x + y for x, y in zip(ka, kb))
                out[key] = sadd(out.get(key, {}), smul(ca, cb))
    return {k: v for k, v in out.items() if v}

def jadd(A, B):
    out = dict(A)
    for k, v in B.items():
        out[k] = sadd(out.get(k, {}), v)
    return {k: v for k, v in out.items() if v}

def jscale(A, r):
    return {k: sscale(v, r) for k, v in A.items() if sscale(v, r)}

def zvec(h):
    return (0,) * h

def jconst(r, h):
    return {zvec(h): {(0,0,0): Fr(r)}} if r != 0 else {}

def jlin(v, h, c=Fr(1)):
    """c * (v . rho)"""
    out = {}
    for a in range(h):
        if v[a] * c:
            key = tuple(1 if i == a else 0 for i in range(h))
            out[key] = {(0,0,0): Fr(v[a]) * c}
    return out

def jinv_unit(A, h):
    c0 = A.get(zvec(h), {}).get((0,0,0), Fr(0))
    assert c0 != 0 and list(A.get(zvec(h), {})) == [(0,0,0)], "non-unit jet"
    r = jscale(jadd(A, jconst(-c0, h)), Fr(1)/c0)
    out = jconst(1, h); pw = jconst(1, h)
    for k in range(1, TR + 1):
        pw = jmul(pw, r)
        out = jadd(out, jscale(pw, Fr((-1)**k)))
    return jscale(out, Fr(1)/c0)

def invgamma1(v, h):
    """jet of 1/Gamma(1+x), x = v.rho: exp(g x - z2 x^2/2 + z3 x^3/3 - ...)"""
    x = jlin(v, h); x2 = jmul(x, x); x3 = jmul(x2, x)
    out = jconst(1, h)
    out = jadd(out, {k: smul(c, S_G) for k, c in x.items()})
    t2 = sadd(sscale(smul(S_G, S_G), Fr(1,2)), sscale(S_Z2, Fr(-1,2)))
    out = jadd(out, {k: smul(c, t2) for k, c in x2.items()})
    t3 = sadd(sadd(sscale(smul(smul(S_G,S_G),S_G), Fr(1,6)),
                   sscale(smul(S_G,S_Z2), Fr(-1,2))), sscale(S_Z3, Fr(1,3)))
    out = jadd(out, {k: smul(c, t3) for k, c in x3.items()})
    return out

def recip_gamma_tables(v, h, amin, amax):
    """R[a] = jet of 1/Gamma(1+a+x), x = v.rho, for a in [amin, amax]."""
    R = {0: invgamma1(v, h)}
    x = jlin(v, h)
    for a in range(0, amax):
        inv = jconst(Fr(1, a+1), h)
        pw = jconst(1, h)
        for k in range(1, TR + 1):
            pw = jmul(pw, x)
            inv = jadd(inv, jscale(pw, Fr((-1)**k, (a+1)**(k+1))))
        R[a+1] = jmul(R[a], inv)
    for a in range(0, amin, -1):
        R[a-1] = jmul(R[a], jadd(x, jconst(a, h)))
    return R

class GammaJets:
    """A(n+rho) jets for a general card, exact, all sectors on a box n<=nmax."""
    def __init__(self, card, nmax):
        self.card = card; self.h = card.h
        self.T = []
        for (v, mult, is_num) in card.factors:
            lo = sum(min(0, v[a] * nmax[a]) for a in range(self.h))
            hi = sum(max(0, v[a] * nmax[a]) for a in range(self.h))
            self.T.append((v, mult, is_num,
                           recip_gamma_tables(v, self.h, lo, hi)))
    def A(self, n):
        J = jconst(1, self.h)
        for (v, mult, is_num, R) in self.T:
            a = sum(v[i] * n[i] for i in range(self.h))
            if is_num:
                # Pitfall guarded here: a numerator-Gamma pole means the
                # card basis violates deg(-K.ell^(a)) >= 0 — the Laurent
                # 1/(v.rho) jet is direction-dependent and MUST NOT enter
                # silently, so it is reported by name.
                assert a >= 0, (f"numerator Gamma pole at n={n}: card basis "
                                "violates the deg(-K.ell)>=0 rule "
                                "(see eff_route: non-simplicial effective cones)")
            f = jinv_unit(R[a], self.h) if is_num else R[a]
            for _ in range(mult):
                J = jmul(J, f)
        for k, c in J.items():
            for mono in c:
                assert mono[0] == 0, f"gamma survives at n={n}"
        return J

def cfrac(card, n):
    """cone-sector coefficient c(n) as exact Fraction (0 out of cone).
    Negative-argument rule: factorials are evaluated
    by SIGNED Gamma-argument with no table lookups (a python factorial-table
    indexed by a signed arg silently wraps and corrupts the series while it
    still looks support-consistent).  Numerator-Gamma poles (arg < 0) are
    admissible ONLY when denominator zeros strictly dominate (-> exact 0);
    a balanced pole/zero count is direction-dependent and fail-closes on the
    deg(-K.ell)>=0 basis rule."""
    num = 1; den = 1; npole = 0; nzero = 0
    for (v, mult, is_num) in card.factors:
        a = sum(v[i] * n[i] for i in range(card.h))
        if is_num:
            if a < 0:
                npole += mult
            else:
                num *= factorial(a) ** mult
        else:
            if a < 0:
                nzero += mult
            else:
                den *= factorial(a) ** mult
    if nzero > npole:
        return Fr(0)
    assert npole == 0, (f"numerator Gamma pole not dominated at n={n}: card "
                        "basis violates the deg(-K.ell)>=0 rule "
                        "(see eff_route: non-simplicial effective cones)")
    return Fr(num, den)

def alphas(h, K):
    """all h-tuples >= 0 with total degree <= K, sorted."""
    out = []
    def rec(a, rem, cur):
        if a == h:
            out.append(tuple(cur)); return
        for x in range(rem + 1):
            cur.append(x); rec(a + 1, rem - x, cur); cur.pop()
    rec(0, K, [])
    return sorted(out)

def curve_points(card, nu, M, cone=True):
    """all n with nu.n <= M (cone=True: only c(n) != 0 sector)."""
    h = card.h
    dens = [v for (v, m, isn) in card.factors if not isn]
    out = []
    def rec(a, rem, cur):
        if a == h:
            n = tuple(cur)
            if not cone or all(sum(v[i]*n[i] for i in range(h)) >= 0
                               for v in dens):
                out.append(n)
            return
        for na in range(rem // nu[a] + 1):
            cur.append(na); rec(a + 1, rem - nu[a]*na, cur); cur.pop()
    rec(0, M, [])
    return out

def sig_tw(sigma, n):
    tw = 1
    for i, x in enumerate(n):
        tw *= sigma[i] ** x
    return Fr(tw)

def w0_curve(card, nu, sigma, M):
    """restricted fundamental period sum a_m s^m, exact integers."""
    a = [Fr(0)] * (M + 1)
    for n in curve_points(card, nu, M, cone=True):
        m = sum(nu[i]*n[i] for i in range(card.h))
        a[m] += sig_tw(sigma, n) * cfrac(card, n)
    for x in a:
        assert x.denominator == 1, "non-integer w0 coefficient"
    return a

def jet_tower(card, nu, sigma, M):
    """C[beta][m] = {(e2,e3): Fr} — restricted rho-jet tower, |beta| <= TR.
    C_beta = beta! * (jet component beta) = del_rho^beta w-part (convention)."""
    h = card.h
    nmax = tuple(M // nu[i] for i in range(h))
    GJ = GammaJets(card, nmax)
    BET = alphas(h, TR)
    C = {b: [dict() for _ in range(M + 1)] for b in BET}
    for n in curve_points(card, nu, M, cone=False):
        m = sum(nu[i]*n[i] for i in range(h))
        J = GJ.A(n)
        tw = sig_tw(sigma, n)
        for b, sym in J.items():
            if b not in C: continue
            fb = 1
            for x in b: fb *= factorial(x)
            cur = C[b][m]
            for mono, v in sym.items():
                key = (mono[1], mono[2])
                w = cur.get(key, Fr(0)) + tw * v * fb
                if w: cur[key] = w
                elif key in cur: del cur[key]
    return C

# ------------- h-var truncated series on a box N (tuple), symdict coeffs
def tmul(A, B, N):
    h = len(N)
    out = {}
    for ka, ca in A.items():
        for kb, cb in B.items():
            key = tuple(ka[i] + kb[i] for i in range(h))
            if all(key[i] <= N[i] for i in range(h)):
                out[key] = sadd(out.get(key, {}), smul(ca, cb))
    return {k: v for k, v in out.items() if v}

def tadd(A, B):
    out = dict(A)
    for k, v in B.items():
        out[k] = sadd(out.get(k, {}), v)
    return {k: v for k, v in out.items() if v}

def tscale(A, r):
    return {k: sscale(v, r) for k, v in A.items() if sscale(v, r)}

def tinv_unit(A, N):
    h = len(N)
    z = zvec(h)
    assert A.get(z, {}).get((0,0,0), Fr(0)) == 1 and list(A.get(z, {})) == [(0,0,0)]
    r = dict(A); del r[z]
    out = {z: S_ONE}; pw = {z: S_ONE}
    for k in range(1, sum(N) + 1):
        pw = tmul(pw, r, N)
        if not pw: break
        out = tadd(out, tscale(pw, Fr((-1)**k)))
    return out

def texp(X, N):
    h = len(N)
    z = zvec(h)
    assert z not in X
    out = {z: S_ONE}; pw = {z: S_ONE}
    for k in range(1, sum(N) + 1):
        pw = tmul(pw, X, N)
        if not pw: break
        out = tadd(out, tscale(pw, Fr(1, factorial(k))))
    return out

def tcompose(F, Gs, N):
    """F(z) with z_a = q_a * G_a(q); Gs = list of the unit series G_a in q."""
    h = len(N)
    P = []
    for a in range(h):
        pw = [{zvec(h): S_ONE}]
        for i in range(1, N[a] + 1):
            pw.append(tmul(pw[-1], Gs[a], N))
        P.append(pw)
    out = {}
    for k, c in F.items():
        red = tuple(N[i] - k[i] for i in range(h))
        base = {zvec(h): S_ONE}
        for a in range(h):
            if k[a]:
                base = tmul(base, P[a][k[a]], red)
        for kk, s2 in base.items():
            key = tuple(k[i] + kk[i] for i in range(h))
            out[key] = sadd(out.get(key, {}), smul(s2, c))
    return {k: v for k, v in out.items() if v}

def e_(a, h):
    return tuple(1 if i == a else 0 for i in range(h))

def build_box(card, N, sig):
    """jet series W[alpha] (|alpha| <= 2) on box N, coeffs twisted by sigma^n;
    W[alpha] = del^alpha w (the alpha! * jet-component convention)."""
    h = card.h
    GJ = GammaJets(card, N)
    ALS = alphas(h, 2)
    W = {al: {} for al in ALS}
    for n in product(*[range(N[a] + 1) for a in range(h)]):
        J = GJ.A(n)
        tw = sig_tw(sig, n)
        for al in ALS:
            c = J.get(al)
            if c:
                fac = 1
                for x in al: fac *= factorial(x)
                W[al][n] = sscale(c, tw * fac)
    return W

# ------------- effective-cone route (non-simplicial cards)
def _prim(v):
    g = 0
    for x in v:
        g = gcd(g, abs(int(x)))
    return tuple(int(x) // g for x in v) if g else tuple(int(x) for x in v)

def _rref(rows, h):
    """exact row reduction; returns (rank, pivot->row map, matrix)."""
    M = [[Fr(x) for x in r] for r in rows]
    piv = {}; r = 0
    for c in range(h):
        p = next((i for i in range(r, len(M)) if M[i][c] != 0), None)
        if p is None:
            continue
        M[r], M[p] = M[p], M[r]
        pv = M[r][c]
        M[r] = [x / pv for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c] != 0:
                f = M[i][c]
                M[i] = [M[i][t] - f * M[r][t] for t in range(h)]
        piv[c] = r; r += 1
    return r, piv, M

def _nullvec(rows, h):
    """primitive integer null vector when the row space has corank 1."""
    r, piv, M = _rref(rows, h)
    free = [c for c in range(h) if c not in piv]
    if len(free) != 1:
        return None
    fc = free[0]
    v = [Fr(0)] * h; v[fc] = Fr(1)
    for c, rr in piv.items():
        v[c] = -M[rr][fc]
    den = 1
    for x in v:
        den = den * x.denominator // gcd(den, x.denominator)
    return _prim([int(x * den) for x in v])

def cone_facets(rays, h):
    """exact inward facet normals of the full-dimensional cone(rays); the
    cheap-LP layer of the simpliciality check (no floats anywhere)."""
    assert _rref(rays, h)[0] == h, "effective cone not full-dimensional"
    facets = set()
    for sub in combinations(range(len(rays)), h - 1):
        nv = _nullvec([rays[i] for i in sub], h)
        if nv is None:
            continue
        vals = [sum(nv[j] * r[j] for j in range(h)) for r in rays]
        if all(v >= 0 for v in vals):
            facets.add(nv)
        elif all(v <= 0 for v in vals):
            facets.add(tuple(-x for x in nv))
    facets = sorted(facets)
    for f in facets:   # H-rep gate: each facet supports h-1 independent rays
        on = [r for r in rays if sum(f[j] * r[j] for j in range(h)) == 0]
        assert _rref(on, h)[0] == h - 1, "facet gate failed"
    return facets

def cone_extreme(rays, h):
    """reduce a generating set to its extreme rays (exact, small ray counts:
    r is redundant iff r in cone(rays - {r}), tested on that cone's facets)."""
    rays = sorted(set(_prim(r) for r in rays))
    assert all(any(r) for r in rays), "zero ray"
    ext = []
    for i, r in enumerate(rays):
        others = rays[:i] + rays[i + 1:]
        if len(others) >= h and _rref(others, h)[0] == h:
            F = cone_facets(others, h)
            if all(sum(f[j] * r[j] for j in range(h)) >= 0 for f in F):
                continue                        # r redundant
        ext.append(r)
    return ext

class EffRoute:
    def __init__(self, simplicial, rays, facets):
        self.simplicial = simplicial
        self.rays = rays          # extreme rays (primitive, sorted)
        self.facets = facets      # inward facet normals (H-rep)
    def member(self, n):
        return all(sum(f[j] * n[j] for j in range(len(n))) >= 0
                   for f in self.facets)

def eff_route(card):
    """SIMPLICIALITY CHECK FIRST — run this before ANY GV work on a card
    (non-simplicial effective cones break gv_extract).  Cards without eff_rays
    declare the basis orthant as the effective cone (legacy simplicial
    contract): every downstream path is then byte-identical to the pre-route
    code.  Cards with eff_rays get an exact extreme-ray reduction; simplicial
    iff the extreme set is exactly the h unit vectors."""
    rt = getattr(card, "_eff_route", None)
    if rt is not None:
        return rt
    h = card.h
    raw = getattr(card, "eff_rays", None)
    if raw is None:
        raw = card.raw.get("eff_rays")
    if not raw:
        units = [e_(a, h) for a in range(h)]
        rt = EffRoute(True, units, list(units))
    else:
        rays = [tuple(int(x) for x in r) for r in raw]
        assert all(len(r) == h for r in rays), "eff_rays dimension mismatch"
        ext = cone_extreme(rays, h)
        simp = sorted(ext) == sorted(e_(a, h) for a in range(h))
        rt = EffRoute(simp, ext, cone_facets(ext, h))
        for d in card.gv_pinned:
            assert rt.member(d), f"pinned GV class {d} outside effective cone"
        tbl = getattr(card, "gv_table", None) or {}
        for d in tbl:
            assert rt.member(d), f"gv_table class {d} outside effective cone"
    card._eff_route = rt
    return rt

def gkz_box_degree_factors(l0):
    """Rule: the degree (numerator) column of a GKZ
    box operator contributes factors (D0 - k) for k = 1..|l0| with overall
    sign (-1)^{|l0|} — NOT k = 0..|l0|-1 (the naive transcription fails
    termwise; caught exactly on ads-5-113-4627-main, ratio (a0+8)/a0).
    Ordinary denominator columns keep k = 0..|lj|-1.  Returns (sign, [k...])
    for card derivers building box ideals."""
    m = abs(int(l0))
    return (-1) ** m, list(range(1, m + 1))

def gv_extract(card, N, sig, eps):
    """h-dim HKTY double-log GV extraction, holomorphic limit; returns
    (nGV dict d -> Fr, Gs mirror-map inverse units).  SIMPLICIAL ROUTE ONLY:
    on a non-simplicial effective cone the identity is provably inconsistent
    in any containing basis (exact kappa-fit on a 7-ray h=5 example) and
    semigroup-filtering the sum fails at total degree 1 (measured) — so this
    fail-closes; use pin_sign_frame's gv_table route instead."""
    assert eff_route(card).simplicial, (
        "STRUCTURAL (non-simplicial effective cone): effective "
        f"cone has {len(eff_route(card).rays)} > h extreme rays; double-log "
        "identity + instanton peel are inconsistent in ANY containing basis. "
        "Route: pin_sign_frame -> card gv_table + racetrack-peel validation.")
    h = card.h
    W = build_box(card, N, sig)
    z = zvec(h)
    w0 = W[z]
    w0inv = tinv_unit(w0, N)
    S = {a: W[e_(a, h)] for a in range(h)}
    Sn = {a: tmul(S[a], w0inv, N) for a in range(h)}
    def SS(b, c):
        al = tuple(e_(b, h)[i] + e_(c, h)[i] for i in range(h))
        return W[al]
    # P_a = (1/2) kappa_abc [ S_bc/w0 - (S_b/w0)(S_c/w0) ]  (log-free identity)
    P = {}
    for a in range(h):
        acc = {}
        for b in range(h):
            for c in range(h):
                kv = card.kap(a, b, c)
                if kv == 0: continue
                t = tadd(tmul(SS(b, c), w0inv, N),
                         tscale(tmul(Sn[b], Sn[c], N), Fr(-1)))
                acc = tadd(acc, tscale(t, kv / 2))
        P[a] = acc
    # mirror map: q_a = z_a exp(S_a/w0) -> invert z_a = q_a G_a(q)
    E = {a: texp(Sn[a], N) for a in range(h)}
    Gs = [{z: S_ONE} for _ in range(h)]
    for _ in range(sum(N) + 2):
        newGs = []
        for a in range(h):
            Eq = tcompose(E[a], Gs, N)
            newGs.append(tinv_unit(Eq, N))
        if newGs == Gs:
            Gs = newGs; break
        Gs = newGs
    Pq = {a: tcompose(P[a], Gs, N) for a in range(h)}
    for a in range(h):
        for d, c in Pq[a].items():
            if d != z:
                for mono, v in c.items():
                    assert mono == (0,0,0), f"zeta survives at q^{d}: {c}"
    # peel multicovers: coef_a(e) = eps * sum_{e=k d} n_d d_a / k^2
    nGV = {}
    degs = sorted(set(d for a in range(h) for d in Pq[a] if d != z),
                  key=lambda d: (sum(d), d))
    for e in degs:
        vals = {}
        for a in range(h):
            coef = s_rat(Pq[a].get(e, {}))
            rest = Fr(0)
            k = 2
            while k <= max(e):
                if all(x % k == 0 for x in e):
                    d = tuple(x // k for x in e)
                    if d in nGV:
                        rest += nGV[d] * Fr(d[a], k**2)
                k += 1
            if e[a] != 0:
                vals[a] = (Fr(eps) * coef - rest) / e[a]
            else:
                assert Fr(eps) * coef - rest == 0, f"deg-{e} a={a} inconsistency"
        vv = sorted(set(vals.values()))
        assert len(vv) == 1, f"axis mismatch at {e}: {vals}"
        nGV[e] = vv[0]
    return nGV, Gs

def pin_sign_frame(card, N=None):
    """search sigma in {+-1}^h, eps in {+-1} reproducing the pinned GV leaders;
    returns (sigma, eps, nGV).  If the card pins sigma/eps, verify only.
    NON-SIMPLICIAL ROUTE (eff_route): gv_extract is structurally unavailable;
    sigma/eps must be card-pinned and nGV comes from the provenance-gated
    card gv_table, gated here (cone containment via eff_route, integrality,
    leader pins) and downstream by the racetrack peel (mirror_curve) — the
    layers measured to stay exact on non-simplicial cones."""
    h = card.h
    rt = eff_route(card)
    if not rt.simplicial:
        assert card.sigma is not None and card.eps is not None, \
            "non-simplicial route: sigma/eps must be pinned on the card"
        tbl = getattr(card, "gv_table", None)
        assert tbl, ("non-simplicial route: card gv_table required "
                     "(provenance-gated spec GV classes in the card basis)")
        nGV = {tuple(d): Fr(v) for d, v in tbl.items()}
        assert all(v.denominator == 1 for v in nGV.values()), \
            "gv_table non-integer entry"
        for d, v in card.gv_pinned.items():
            assert nGV.get(d) == v, \
                f"gv_table vs pinned leader mismatch at {d}: {nGV.get(d)} != {v}"
        return tuple(card.sigma), card.eps, nGV
    N = N or card.gv_window
    leaders = {d: v for d, v in card.gv_pinned.items() if sum(d) <= 2}
    cands = ([ (card.sigma, card.eps) ] if card.sigma else
             [(s, e) for s in product(*[(1, -1)] * h) for e in (1, -1)])
    for sig, eps in cands:
        try:
            nGV, _ = gv_extract(card, N, sig, eps)
        except AssertionError:
            continue
        if all(nGV.get(d) == v for d, v in leaders.items()):
            return tuple(sig), eps, nGV
    raise AssertionError("no sign convention reproduces the pinned GV leaders")

# ------------- 1-d exact series helpers (lists of Fr on [0..M])
def smul1(x, y, M):
    out = [Fr(0)] * (M + 1)
    for i, xi in enumerate(x[:M + 1]):
        if xi:
            for j in range(0, M + 1 - i):
                if y[j]: out[i + j] += xi * y[j]
    return out

def sinv1(x, M):
    assert x[0] == 1
    out = [Fr(0)] * (M + 1); out[0] = Fr(1)
    for m in range(1, M + 1):
        out[m] = -sum(x[k] * out[m - k] for k in range(1, m + 1) if x[k])
    return out

def sexp1(x, M):
    assert x[0] == 0
    out = [Fr(0)] * (M + 1); out[0] = Fr(1)
    for m in range(1, M + 1):
        out[m] = sum((Fr(k) * x[k] * out[m - k] for k in range(1, m + 1)
                      if x[k]), Fr(0)) / m
    return out

def spow1(x, k, M):
    out = [Fr(0)] * (M + 1); out[0] = Fr(1)
    for _ in range(k):
        out = smul1(out, x, M)
    return out

def scomp1(ser, inner, M):
    """ser(inner(q)) — inner has zero constant term."""
    out = [Fr(0)] * (M + 1)
    pw = [Fr(0)] * (M + 1); pw[0] = Fr(1)
    for k in range(M + 1):
        if k < len(ser) and ser[k]:
            for i, v in enumerate(pw):
                if v: out[i] += ser[k] * v
        pw = smul1(pw, inner, M)
        if all(v == 0 for v in pw): break
    return out

def bezout_vec(nu):
    """integer c with c.nu = gcd(nu) (= 1 for a primitive curve)."""
    from math import gcd
    c = [0] * len(nu)
    g, c[0] = nu[0], 1
    for i in range(1, len(nu)):
        # extended gcd(g, nu[i])
        a, b = g, nu[i]
        x0, x1, y0, y1 = 1, 0, 0, 1
        while b:
            q, a, b = a // b, b, a % b
            x0, x1 = x1, x0 - q * x1
            y0, y1 = y1, y0 - q * y1
        for j in range(i):
            c[j] *= x0
        c[i] = y0
        g = a
    return c, g

def mirror_curve(card, nu, C, nGV, M, MN=20):
    """restricted mirror map + instanton peel on the curve.
    C = jet tower from jet_tower();  nGV = extracted GV table.
    Flat curve coordinate: q_s = prod_a q_a^{c_a} with c.nu = 1 (Bezout).
    Returns dict with q_s(s)/s, s(q_s)/q_s, E_a, N_m, naive sums, gates."""
    h = card.h
    z = zvec(h)
    def ser_rat(b):
        return [C[b][m].get((0, 0), Fr(0)) for m in range(M + 1)]
    w0 = ser_rat(z)
    w0i = sinv1(w0, M)
    Sg = [smul1(ser_rat(e_(a, h)), w0i, M) for a in range(h)]
    E = [sexp1(Sg[a], M) for a in range(h)]        # q_a = s^{nu_a} E_a(s)
    cvec, g = bezout_vec(list(nu))
    assert g == 1, "curve exponents not primitive"
    D = [sum(cvec[a] * Sg[a][m] for a in range(h)) for m in range(M + 1)]
    Es = sexp1(D, M)                               # q_s = s * Es(s)
    # invert q_s = s Es(s):  s = q_s * G(q_s)   (fixed-point iteration)
    G = [Fr(0)] * (M + 1); G[0] = Fr(1)
    for _ in range(M + 1):
        sq = [Fr(0)] * (M + 1)
        for m in range(M): sq[m + 1] = G[m]
        comp = scomp1(Es, sq, M)
        Gn = sinv1(comp, M)
        if Gn == G: break
        G = Gn
    sq = [Fr(0)] * (M + 1)
    for m in range(M): sq[m + 1] = G[m]
    ver = smul1(sq, scomp1(Es, sq, M), M)
    assert ver[1] == 1 and all(ver[m] == 0 for m in range(2, M // 2))
    # instanton sum in q_s to order MN, peel N_m
    sq_t = sq[:MN + 1]
    Eq = [scomp1(E[a], sq_t, MN) for a in range(h)]
    qq = [smul1(spow1(sq_t, nu[a], MN), Eq[a], MN) for a in range(h)]
    for a in range(h):
        if nu[a] <= MN:
            assert qq[a][nu[a]] == 1 and all(qq[a][i] == 0
                                             for i in range(nu[a]))
        else:
            # Rule: nu_a > MN — q_a = s^{nu_a}(1+...) vanishes identically
            # at this truncation; classes needing it have lead > MN and are
            # skipped below.
            assert all(x == 0 for x in qq[a])
    Fsum = [Fr(0)] * (MN + 1)
    for d, nd in sorted(nGV.items()):
        lead = sum(nu[a] * d[a] for a in range(h))
        if lead == 0 or lead > MN: continue
        qd = [Fr(1) if i == 0 else Fr(0) for i in range(MN + 1)]
        for a in range(h):
            qd = smul1(qd, spow1(qq[a], d[a], MN), MN)
        xk = qd[:]
        k = 1
        while k * lead <= MN:
            for i, v in enumerate(xk):
                if v: Fsum[i] += Fr(nd) * v / k**3
            xk = smul1(xk, qd, MN)
            k += 1
    Nm = {}
    for m in range(1, MN + 1):
        c = Fsum[m]
        for k in range(2, m + 1):
            if m % k == 0 and (m // k) in Nm:
                c -= Fr(Nm[m // k], k**3)
        if c != 0 or m in card.racetrack:
            Nm[m] = c
    naive = {}
    for d, nd in nGV.items():
        mm = sum(nu[a] * d[a] for a in range(h))
        if 0 < mm <= MN: naive[mm] = naive.get(mm, 0) + nd
    rt_pass = all(Nm.get(m) == v for m, v in card.racetrack.items())
    return dict(qs_of_s=Es, s_of_qs=G, E=E, cvec=cvec, Nm=Nm, naive=naive,
                racetrack_pass=rt_pass,
                integrality_pass=all(v.denominator == 1 for v in Nm.values()))


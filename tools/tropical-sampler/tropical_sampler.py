"""BSST tropical importance sampling (arXiv 2204.06414, Algorithm 1 reimplementation,
python, no Polymake) for the phylo evidence integral

  Z(T;u,alpha=1) = 4^{-(E+1)N} * int_{R+^E} prod_k Q_k(t)^{u_k} prod_e (1+t_e)^{-(N+2)} dt

with Q_k the subtraction-free (positivity-gated) pattern presentation.
Log coords y=log t. Tropical PL exponent:
  Phi(y) = sum_k u_k max_{a in supp Q_k} <a,y> - (N+2) sum_e relu(y_e) + sum_e y_e.
Fan = normal fan of M = sum_k conv(supp Q_k) + sum_e [0,e_e]  (Minkowski, exact ints).
Per max cone sigma (vertex of M): gradient g = sum_k u_k a_k(sigma) + 1 - (N+2)s(sigma);
sector integral c_sigma = |det R| / prod_i beta_i, beta_i = -<g, r_i> > 0 (exact Fraction).
I^tr = sum c_sigma. Sample: cone ~ c_sigma/I^tr, y = R @ lambda, lambda_i ~ Exp(beta_i).
Weight w = exp(sum_k u_k [log Q_k - Phi_k] - (N+2) sum_e log1p(exp(-|y_e|))), bounded.
Z_hat = prefac * I^tr * mean(w).
"""
import numpy as np
from fractions import Fraction
from itertools import combinations
from scipy.spatial import ConvexHull


# ---------- exact integer determinant (Bareiss) ----------
def idet(M):
    A = [[int(x) for x in row] for row in M]
    n = len(A); sign = 1; prev = 1
    for k in range(n - 1):
        if A[k][k] == 0:
            for i in range(k + 1, n):
                if A[i][k]:
                    A[k], A[i] = A[i], A[k]; sign = -sign; break
            else:
                return 0
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                A[i][j] = (A[i][j] * A[k][k] - A[i][k] * A[k][j]) // prev
        prev = A[k][k]
    return sign * A[n - 1][n - 1]


def primitive_normal(pts):
    """Exact primitive integer normal of the hyperplane through E integer points
    (rows). Returns int vector n via cofactor expansion; None if degenerate."""
    import math
    P = [ [int(x) for x in p] for p in pts ]
    E = len(P[0])
    M = [[P[i][j] - P[0][j] for j in range(E)] for i in range(1, E)]  # (E-1) x E
    n = []
    for j in range(E):
        sub = [[row[c] for c in range(E) if c != j] for row in M]
        n.append((-1) ** j * idet(sub))
    if all(v == 0 for v in n):
        return None
    g = 0
    for v in n:
        g = math.gcd(g, abs(v))
    return tuple(v // g for v in n)


def hull_vertices(pts):
    """Vertices of conv(pts) (int array N x E) via qhull. Handles lower-dimensional
    point sets by projecting onto independent coordinates of the affine span
    (coords are small ints; float rank is safe)."""
    pts = np.unique(np.asarray(pts, dtype=np.int64), axis=0)
    E = pts.shape[1]
    if len(pts) <= 2:
        return pts
    D = (pts - pts[0]).astype(float)
    # greedy pivot columns
    cols, r = [], 0
    for j in range(E):
        if np.linalg.matrix_rank(D[:, cols + [j]], tol=1e-9) > r:
            cols.append(j); r += 1
    if r == 0:
        return pts[:1]
    if len(pts) <= r + 1:
        return pts
    P = pts[:, cols].astype(float)
    if r == 1:
        keep = [int(np.argmin(P[:, 0])), int(np.argmax(P[:, 0]))]
        return pts[sorted(set(keep))]
    h = ConvexHull(P, qhull_options='Qx Q12')
    return pts[np.array(sorted(set(h.vertices)))]


def minkowski_vertices(point_sets, log=None):
    """Iterated Minkowski sum, pruning to hull vertices each step."""
    V = np.zeros((1, point_sets[0].shape[1]), dtype=np.int64)
    for i, P in enumerate(point_sets):
        cand = (V[:, None, :] + P[None, :, :]).reshape(-1, V.shape[1])
        V = hull_vertices(cand)
        if log:
            log(f"  minkowski step {i+1}/{len(point_sets)}: +{len(P)} pts -> {len(V)} vertices")
    return V


# ---------- fan construction ----------
def build_fan(groups, E, rng, log=None, refine_extra=None, skip_chain=0):
    """groups: list of (weight:int, supp:int array S x E). Includes pattern groups
    (weight u_k) and segment groups (weight -(N+2), supp {0,e_e}).
    refine_extra: optional int array of points added as weight-0 summand (fan refine).
    Returns dict with cones (R int rays ExE cols, g int vec, beta, csigma Fraction),
    Itr (Fraction)."""
    chain = groups[:len(groups) - skip_chain] if skip_chain else groups
    psets = [hull_vertices(s) for (_, s) in chain]
    if refine_extra is not None:
        psets = psets + [hull_vertices(refine_extra)]
    V = minkowski_vertices(psets, log=log)
    if log:
        log(f"  Minkowski polytope: {len(V)} vertices, dim {E}")
    hull = ConvexHull(V.astype(float), qhull_options='Qt Qx Q12')
    # exact primitive outward facet normals from simplicial facets
    normals = {}
    for simp, eq in zip(hull.simplices, hull.equations):
        key = tuple(np.round(eq[:E] / max(1e-300, np.linalg.norm(eq[:E])), 8))
        if key in normals:
            continue
        n = primitive_normal(V[simp])
        if n is None:
            continue
        nv = np.array(n, dtype=object)
        vals = V.astype(object) @ nv
        mx, mn = max(vals), min(vals)
        h0 = int(np.array([int(x) for x in V[simp][0]], dtype=object) @ nv)
        if mx == mn:
            continue                      # degenerate direction
        if h0 == mn:                      # inward -> flip
            n = tuple(-x for x in n)
        elif h0 != mx:
            raise RuntimeError("facet normal neither supports max nor min")
        normals[key] = n
    # exact dedupe + exact incidence (int64 with overflow guard)
    NN = sorted(set(normals.values()))
    NNa = np.array(NN, dtype=np.int64)                     # F x E
    maxprod = float(np.abs(NNa).max()) * float(np.abs(V).max()) * E
    assert maxprod < 2**62, f"incidence overflow risk: {maxprod:.2e}"
    VALS = V.astype(np.int64) @ NNa.T                      # V x F
    HMAX = VALS.max(axis=0)
    INC = (VALS == HMAX[None, :])                          # vertex-facet incidence
    if log:
        log(f"  {len(NN)} exact facet normals")
    # per-vertex normal cone -> simplicial subcones
    cones = []
    n_nonsimple = 0
    for vi in range(len(V)):
        rays = [NN[f] for f in range(len(NN)) if INC[vi, f]]
        if len(rays) < E:
            continue  # not a vertex / degenerate
        if len(rays) == E:
            tris = [tuple(range(E))]
        else:
            n_nonsimple += 1
            Rn = np.array(rays, dtype=float)
            Rn = Rn / np.linalg.norm(Rn, axis=1, keepdims=True)
            if np.linalg.matrix_rank(Rn, tol=1e-9) < E:
                continue          # lower-dim normal cone: not a vertex (measure zero)
            P0 = np.vstack([np.zeros(E), Rn])
            try:
                ch = ConvexHull(P0, qhull_options='Qt Qx Q12')
            except Exception:
                # near-degenerate ray configuration: joggle for the combinatorics
                # only (exact rays still used downstream; tiling probe guards)
                ch = ConvexHull(P0, qhull_options='QJ Qt Q12')
            tris = [tuple(s - 1 for s in simp if s != 0)
                    for simp in ch.simplices if 0 not in simp]
            tris = [t for t in tris if len(t) == E]
        for t in tris:
            R = np.array([rays[i] for i in t], dtype=object).T   # E x E, cols=rays
            d = abs(idet(R.T))
            if d == 0:
                continue
            cones.append(dict(vi=vi, R=R, det=d))
    if log:
        log(f"  {len(cones)} simplicial cones ({n_nonsimple} non-simple vertices triangulated)")
    # per-cone gradient via exact interior direction + argmax per group
    Itr = Fraction(0)
    ok_cones = []
    for c in cones:
        R = c['R']
        Ri = R.astype(np.int64)
        assert float(np.abs(Ri).max()) * 1000 * E < 2**62
        for attempt in range(60):
            coef = rng.integers(1, 1000, size=E)
            y0 = Ri @ coef                   # exact integer interior direction (int64)
            gs, good = [], True
            for (w, supp) in groups:
                sc = supp @ y0               # supp entries 0/1: no overflow beyond y0 sums
                mx = sc.max()
                arg = np.flatnonzero(sc == mx)
                if len(arg) != 1:
                    good = False; break
                gs.append((w, supp[arg[0]]))
            if good:
                break
        else:
            raise RuntimeError("no generic interior direction found (cone degenerate?)")
        g = np.array([1] * E, dtype=object)
        for (w, a) in gs:
            g = g + int(w) * a.astype(object)
        beta = [-int(g @ r) for r in R.T]
        assert all(b > 0 for b in beta), \
            f"tropical integral divergent on cone (beta={beta}, g={g.tolist()})"
        cs = Fraction(int(c['det']))
        for b in beta:
            cs /= b
        c.update(g=np.array([int(x) for x in g]), beta=np.array(beta, dtype=float),
                 csigma=cs, Rf=R.astype(float))
        Itr += cs
        ok_cones.append(c)
    if log:
        log(f"  I^tr = {float(Itr):.6e} (exact rational, {len(ok_cones)} cones)")
    return dict(cones=ok_cones, Itr=Itr, E=E)


# ---------- integrand data ----------
def make_eval_data(patterns, E):
    """patterns: list of (u_k, Dsets list-of-frozensets). Returns per-pattern
    (u_k, membership int8 array n_terms x E deduped with multiplicity, logmult)."""
    out = []
    for u, Dsets in patterns:
        cnt = {}
        for D in Dsets:
            cnt[D] = cnt.get(D, 0) + 1
        Ds = sorted(cnt.keys(), key=lambda s: (len(s), sorted(s)))
        Mb = np.zeros((len(Ds), E))
        for i, D in enumerate(Ds):
            for e in D:
                Mb[i, e] = 1.0
        lm = np.log(np.array([cnt[D] for D in Ds], dtype=float))
        out.append((u, Mb, lm))
    return out


def log_weight(Y, evdata, N, E):
    """Y: n x E sample block. Returns log w (n,)."""
    lw = np.zeros(len(Y))
    l4 = np.log(4.0)
    lf = np.logaddexp(0.0, l4 + Y)          # log(1+4 e^y)  per edge
    rl = np.maximum(Y, 0.0)                 # relu
    for (u, Mb, lm) in evdata:
        s = lf @ Mb.T + lm[None, :]         # n x terms : log latent products
        logQ = s.max(axis=1)
        logQ = logQ + np.log(np.exp(s - logQ[:, None]).sum(axis=1))
        phi = (rl @ Mb.T).max(axis=1)       # tropical max over supp (downward-closed)
        lw += u * (logQ - phi)
    lw -= (N + 2) * np.log1p(np.exp(-np.abs(Y))).sum(axis=1)
    return lw


def sample_tropical(fan, nsamples, rng, batch=500_000):
    """Yields blocks Y (m x E) sampled from the tropical density."""
    cones = fan['cones']
    p = np.array([float(c['csigma'] / fan['Itr']) for c in cones])
    p = p / p.sum()
    E = fan['E']
    left = nsamples
    while left > 0:
        m = min(batch, left)
        counts = rng.multinomial(m, p)
        Ys = []
        for ci in np.nonzero(counts)[0]:
            k = counts[ci]
            lam = rng.exponential(1.0 / cones[ci]['beta'], size=(k, E))
            Ys.append(lam @ cones[ci]['Rf'].T)
        Y = np.vstack(Ys)
        rng.shuffle(Y)
        yield Y
        left -= m


def ln_fraction(fr):
    import math
    return math.log(fr.numerator) - math.log(fr.denominator)


def tiling_probe(fan, rng, nprobe=200, tol=1e-9):
    """Completeness/tiling check: random directions must lie in EXACTLY one cone
    (up to boundary ties). Returns (ok, ncover_min, ncover_max)."""
    E = fan['E']
    Y = rng.normal(0, 2.0, size=(nprobe, E))
    cover = np.zeros(nprobe, dtype=int)
    for c in fan['cones']:
        lam = np.linalg.solve(c['Rf'], Y.T)          # E x nprobe
        cover += (lam.min(axis=0) >= -tol)
    return bool(np.all(cover == 1)), int(cover.min()), int(cover.max())


def _cone_sample_w(cone, k, rng, evdata, N, E):
    lam = rng.exponential(1.0 / cone['beta'], size=(k, E))
    Y = lam @ cone['Rf'].T
    return np.exp(log_weight(Y, evdata, N, E))


def estimate_logZ_tilt(fan, evdata, N, E, n_pilot, n_main, rng, ln_prefac,
                       batch=400_000, min_pilot=64):
    """Stratified tropical sampling with per-cone exponential tilting.
    Pilot (tropical Exp(beta_i) per cone) -> linear fit ln w ~ a + sum_i g_i lam_i
    -> tilted rates bt_i = clip(beta_i - g_i, 0.1 beta_i, beta_i)  (only fatten:
    density ratio q/qt = prod (b/bt) exp(-(b-bt)lam) stays bounded).
    Allocation: Neyman on the pilot-reweighted tilted second moment.
    Main run samples Exp(bt), weight w * q/qt. Unbiased for Z_int = sum c_s E[w]."""
    cones = fan['cones']
    ncone = len(cones)
    cfl = np.array([float(c['csigma']) for c in cones])
    npil = np.maximum(min_pilot, (n_pilot * cfl / cfl.sum()).astype(int))
    tilts, m_est, v_est = [], np.zeros(ncone), np.zeros(ncone)
    for i, c in enumerate(cones):
        k = int(npil[i])
        lam = rng.exponential(1.0 / c['beta'], size=(k, E))
        lw = log_weight(lam @ c['Rf'].T, evdata, N, E)
        A = np.hstack([np.ones((k, 1)), lam])
        coef, *_ = np.linalg.lstsq(A, lw, rcond=None)
        g = coef[1:]
        bt = np.clip(c['beta'] - g, 0.1 * c['beta'], c['beta'])
        tilts.append(bt)
        # pilot-reweighted moments under the tilted proposal
        lr = (np.log(c['beta'] / bt).sum() - lam @ (c['beta'] - bt))  # log q/qt
        w = np.exp(lw)
        m_est[i] = w.mean()
        v_est[i] = max(0.0, np.mean(w * w * np.exp(lr)) - m_est[i] ** 2)
    sd = np.sqrt(v_est)
    sd = np.maximum(sd, 1e-12 * np.maximum(m_est, 1e-300))
    alloc = cfl * sd
    alloc = np.maximum(20, (n_main * alloc / alloc.sum()).astype(int))
    s1 = np.zeros(ncone); s2 = np.zeros(ncone); cnt = np.zeros(ncone, dtype=np.int64)
    wmin, wmax = np.inf, -np.inf
    for i, c in enumerate(cones):
        bt = tilts[i]
        left = int(alloc[i])
        while left > 0:
            k = min(left, batch)
            lam = rng.exponential(1.0 / bt, size=(k, E))
            lw = log_weight(lam @ c['Rf'].T, evdata, N, E)
            lw = lw + np.log(c['beta'] / bt).sum() - lam @ (c['beta'] - bt)
            w = np.exp(lw)
            s1[i] += w.sum(); s2[i] += (w * w).sum(); cnt[i] += k
            wmin = min(wmin, w.min()); wmax = max(wmax, w.max())
            left -= k
    m = s1 / cnt
    v = np.maximum(0.0, s2 / cnt - m * m)
    Zint = float((cfl * m).sum())
    var = float((cfl * cfl * v / cnt).sum())
    return (np.log(Zint) + ln_prefac, np.sqrt(var) / Zint, Zint,
            int(cnt.sum() + npil.sum()), wmin, wmax)


def estimate_logZ_strat(fan, evdata, N, E, n_pilot, n_main, rng, ln_prefac,
                        batch=400_000):
    """Stratified tropical sampling: Z_int = sum_sigma c_sigma E_sigma[w].
    Pilot (prop. to c_sigma, floor 40/cone) -> Neyman allocation n_sigma prop. to
    c_sigma*std_sigma for the main run; pooled per-cone means. Unbiased; the
    proposal within each cone is exactly the tropical density restricted to it."""
    cones = fan['cones']
    ncone = len(cones)
    cfl = np.array([float(c['csigma']) for c in cones])
    # pilot
    np_alloc = np.maximum(40, (n_pilot * cfl / cfl.sum()).astype(int))
    s1 = np.zeros(ncone); s2 = np.zeros(ncone); cnt = np.zeros(ncone, dtype=np.int64)
    wmin, wmax = np.inf, -np.inf
    for i, c in enumerate(cones):
        w = _cone_sample_w(c, int(np_alloc[i]), rng, evdata, N, E)
        s1[i] += w.sum(); s2[i] += (w * w).sum(); cnt[i] += len(w)
        wmin = min(wmin, w.min()); wmax = max(wmax, w.max())
    m = s1 / cnt
    sd = np.sqrt(np.maximum(0.0, s2 / cnt - m * m))
    sd = np.maximum(sd, 1e-12 * np.maximum(m, 1e-300))     # keep every cone alive
    alloc = cfl * sd
    alloc = np.maximum(10, (n_main * alloc / alloc.sum()).astype(int))
    # main (estimate from main samples only -> allocation-independent, unbiased)
    s1 = np.zeros(ncone); s2 = np.zeros(ncone); cnt = np.zeros(ncone, dtype=np.int64)
    for i, c in enumerate(cones):
        left = int(alloc[i])
        while left > 0:
            k = min(left, batch)
            w = _cone_sample_w(c, k, rng, evdata, N, E)
            s1[i] += w.sum(); s2[i] += (w * w).sum(); cnt[i] += k
            wmin = min(wmin, w.min()); wmax = max(wmax, w.max())
            left -= k
    m = s1 / cnt
    v = np.maximum(0.0, s2 / cnt - m * m)
    Zint = float((cfl * m).sum())
    var = float((cfl * cfl * v / cnt).sum())
    se_rel = np.sqrt(var) / Zint
    logZ = np.log(Zint) + ln_prefac
    return logZ, se_rel, Zint, int(cnt.sum()), wmin, wmax


def estimate_logZ(fan, evdata, N, E, nsamples, rng, ln_prefac):
    """Returns (logZ_natural, stderr_logZ, mean_w, se_w, wmin, wmax)."""
    s1 = s2 = 0.0; n = 0; wmin = np.inf; wmax = -np.inf
    for Y in sample_tropical(fan, nsamples, rng):
        w = np.exp(log_weight(Y, evdata, N, E))
        s1 += w.sum(); s2 += (w * w).sum(); n += len(w)
        wmin = min(wmin, w.min()); wmax = max(wmax, w.max())
    mean = s1 / n
    var = max(0.0, s2 / n - mean * mean)
    se = np.sqrt(var / n)
    logZ = np.log(mean) + ln_fraction(fan['Itr']) + ln_prefac
    return logZ, se / mean, mean, se, wmin, wmax


# ---------- selftest ----------
def _selftest():
    """Import + fan-construction smoke (the GUIDE's documented battery):
    an E=2 three-group fan (4 cones, exact I^tr = 25/16), the tiling probe
    on 500 random directions (each must lie in exactly one cone), and a
    weight-0 refine_extra summand (7 cones) that must reproduce the SAME
    exact Fraction I^tr — the fan-refinement-invariance consistency check.
    Deterministic (seeded rng), seconds-scale, no external data.
    Returns the number of failed legs (0 = PASS)."""
    fails = 0
    E = 2
    # Three groups, all 0/1 supports (pattern-membership style):
    #   pattern group  weight +2, supp = the four 0/1 corners of [0,1]^2
    #   edge segments  weight -7, supp {0, e_i} for i = 1, 2
    # Minkowski polytope = [0,2] x [0,2]: 4 vertices -> 4 simplicial cones.
    # Per quadrant cone the gradient g = 1 + 2*argmax - 7*(segment argmaxes):
    #   (-,-): g=( 1, 1)  c=1      (+,-): g=(-4, 1)  c=1/4
    #   (-,+): g=( 1,-4)  c=1/4    (+,+): g=(-4,-4)  c=1/16
    # I^tr = 1 + 1/4 + 1/4 + 1/16 = 25/16 exactly.
    pat = np.array([[0, 0], [1, 0], [0, 1], [1, 1]], dtype=np.int64)
    seg1 = np.array([[0, 0], [1, 0]], dtype=np.int64)
    seg2 = np.array([[0, 0], [0, 1]], dtype=np.int64)
    groups = [(2, pat), (-7, seg1), (-7, seg2)]

    rng = np.random.default_rng(0)
    fan = build_fan(groups, E, rng)

    # 1: base fan cone count
    if len(fan['cones']) != 4:
        fails += 1
        print(f"FAIL 1: base fan has {len(fan['cones'])} cones, expected 4",
              flush=True)

    # 2: exact rational tropical integral
    if fan['Itr'] != Fraction(25, 16):
        fails += 1
        print(f"FAIL 2: I^tr = {fan['Itr']}, expected 25/16", flush=True)

    # 3: tiling probe — 500 directions, each covered exactly once
    ok, cmin, cmax = tiling_probe(fan, rng, nprobe=500)
    if not (ok and cmin == 1 and cmax == 1):
        fails += 1
        print(f"FAIL 3: tiling probe ok={ok} cover=[{cmin},{cmax}], "
              f"expected exactly-once coverage", flush=True)

    # 4: weight-0 refine_extra summand — triangle whose three edge normals
    # (1,-2), (1,1), (-2,1) each split one quadrant: 4 -> 7 cones.
    extra = np.array([[0, 0], [2, 1], [1, 2]], dtype=np.int64)
    fan2 = build_fan(groups, E, rng, refine_extra=extra)
    if len(fan2['cones']) != 7:
        fails += 1
        print(f"FAIL 4: refined fan has {len(fan2['cones'])} cones, expected 7",
              flush=True)

    # 5: fan-refinement invariance — the SAME exact Fraction
    if fan2['Itr'] != fan['Itr']:
        fails += 1
        print(f"FAIL 5: refined I^tr = {fan2['Itr']} != base {fan['Itr']}",
              flush=True)

    # 6: refined fan still tiles exactly
    ok2, cmin2, cmax2 = tiling_probe(fan2, rng, nprobe=500)
    if not (ok2 and cmin2 == 1 and cmax2 == 1):
        fails += 1
        print(f"FAIL 6: refined tiling probe ok={ok2} cover=[{cmin2},{cmax2}]",
              flush=True)

    print(f"tropical_sampler selftest: base {len(fan['cones'])} cones, "
          f"I^tr = {fan['Itr']}; refined {len(fan2['cones'])} cones, "
          f"I^tr = {fan2['Itr']}; tiling 500/500 exact x2 — "
          + ("PASS" if fails == 0 else f"{fails} FAILS"), flush=True)
    return fails


def main(argv):
    if argv == ['--selftest']:
        return _selftest()
    import sys
    print("usage: tropical_sampler.py --selftest", file=sys.stderr, flush=True)
    return 2


if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv[1:]))

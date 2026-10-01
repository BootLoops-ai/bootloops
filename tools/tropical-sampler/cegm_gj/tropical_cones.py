#!/usr/bin/env python3
"""Tropical-cone (normal-fan) decomposition for the X(3,6) string integral.

The corner-singularity failure of plain tensor quadrature (measured order
~n^-1.4) is cured by decomposing log-space R^4 into the maximal cones of the
common refinement of the 10 Newton-polytope normal fans (= normal fan of
their Minkowski sum; combinatorially the Trop+G(3,6) / P(3,6) fan of AHL
1912.08707). In each simplicial subcone with primitive integer rays r_1..r_4:
    x = exp(w), w = sum z_i r_i, y_i = exp(-z_i)
    I_cone = |det R| int_{(0,1)^4} prod y_i^{kappa_i-1} prod_t Q_t(y)^{s_t} dy
    Q_t(y) = sum_{m in N_t} prod_i y_i^{(v_t-m).r_i},  Q_t(0)=1, coeffs +1
so every Q_t is analytic and >=1 on the closed cube: tensor Gauss-Jacobi
(weight y^{kappa-1}) converges geometrically. kappa_i = -(A+sum_t s_t v_t).r_i.
"""
import numpy as np
import itertools, json, math, sys
from fractions import Fraction as Fr

VARS = ['a', 'b', 'c', 'd']
PHAT_MONOS = {
 '135': [(0,0,0,0),(1,0,0,0)],
 '136': [(0,0,0,0),(1,0,0,0),(1,1,0,0)],
 '146': [(0,0,0,0),(0,1,0,0)],
 '235': [(0,0,0,0),(1,0,0,0),(0,0,1,0)],
 '245': [(1,0,0,0),(0,0,1,0)],
 '236': [(0,0,0,0),(1,0,0,0),(0,0,1,0),(1,1,0,0),(0,1,1,0),(0,0,1,1)],
 '246': [(1,0,0,0),(0,0,1,0),(1,1,0,0),(0,1,1,0),(0,0,1,1)],
 '256': [(1,1,0,0),(0,1,1,0),(0,0,1,1)],
 '346': [(0,0,0,0),(0,1,0,0),(0,0,0,1)],
 '356': [(0,1,0,0),(0,0,0,1),(1,0,0,1)],
}
POLY = list(PHAT_MONOS)
VMAT = {t: np.array(PHAT_MONOS[t], dtype=np.int64) for t in POLY}

def vertex_tuples(nsamp=400000, seed=1):
    rng = np.random.default_rng(seed)
    W = rng.integers(-997, 998, size=(nsamp, 4))
    keys = {}
    valid = np.ones(nsamp, dtype=bool)
    idxs = {}
    for t in POLY:
        S = VMAT[t] @ W.T                      # (#verts, nsamp)
        mx = S.max(axis=0)
        ties = (S == mx).sum(axis=0) > 1
        valid &= ~ties
        idxs[t] = S.argmax(axis=0)
    for i in np.nonzero(valid)[0]:
        key = tuple(int(idxs[t][i]) for t in POLY)
        keys.setdefault(key, 0)
        keys[key] += 1
    return keys

def det3(M):
    (a,b,c),(d,e,f),(g,h,i) = M
    return a*(e*i-f*h) - b*(d*i-f*g) + c*(d*h-e*g)

def nullvec4(rows):
    """integer null vector of 3x4 integer matrix (generalized cross product)."""
    M = [list(map(int, r)) for r in rows]
    r = []
    for j in range(4):
        Mj = [[row[k] for k in range(4) if k != j] for row in M]
        r.append((-1)**j * det3(Mj))
    return r

def primitive(r):
    g = 0
    for v in r: g = math.gcd(g, abs(v))
    if g == 0: return None
    return tuple(v//g for v in r)

def cone_rays(key):
    """extreme rays of cone {w: (v_t - u).w >= 0 for all u in N_t}."""
    rows = []
    for ti, t in enumerate(POLY):
        v = VMAT[t][key[ti]]
        for u in VMAT[t]:
            d = v - u
            if d.any():
                rows.append(tuple(int(x) for x in d))
    rows = sorted(set(rows))
    R = np.array(rows, dtype=np.int64)
    rays = set()
    for sub in itertools.combinations(range(len(rows)), 3):
        r = nullvec4([rows[i] for i in sub])
        p = primitive(r)
        if p is None: continue
        pv = np.array(p, dtype=np.int64)
        for cand in (pv, -pv):
            if (R @ cand >= 0).all():
                tight = (R @ cand == 0)
                if np.linalg.matrix_rank(R[tight].astype(float)) == 3:
                    rays.add(tuple(int(x) for x in cand))
    return sorted(rays), rows

def det4(M):
    import numpy as _np
    from fractions import Fraction
    # exact integer 4x4 determinant via cofactor on first row
    M = [list(map(int, r)) for r in M]
    tot = 0
    for j in range(4):
        Mj = [[M[i][k] for k in range(4) if k != j] for i in range(1, 4)]
        tot += (-1)**j * M[0][j] * det3(Mj)
    return tot

def triangulate(rays):
    """triangulate pointed cone (rays as tuples) into simplicial subcones."""
    if len(rays) == 4:
        d = det4(rays)
        assert d != 0
        return [tuple(range(4))]
    Rm = np.array(rays, dtype=float)
    Ru = Rm / np.linalg.norm(Rm, axis=1)[:, None]
    h = Ru.sum(axis=0)                     # interior direction (pointed cone)
    den = Rm @ h
    assert (den > 1e-9).all(), f"cone not pointed w.r.t. h: {den}"
    pts = Rm / den[:, None]                # cross-section {h.w=1}
    ctr = pts.mean(axis=0)
    Q = pts - ctr
    _, _, Vt = np.linalg.svd(Q, full_matrices=False)
    coords = Q @ Vt[:3].T                  # 3-dim coords in the hyperplane
    from scipy.spatial import Delaunay
    tri = Delaunay(coords, qhull_options="QJ")
    out = []
    for simp in tri.simplices:
        sub = [rays[i] for i in simp]
        if det4(sub) != 0:
            out.append(tuple(int(i) for i in simp))
    return out

def build_fan(nsamp=400000):
    keys = vertex_tuples(nsamp)
    cones = []
    for key, hits in sorted(keys.items(), key=lambda kv: -kv[1]):
        rays, rows = cone_rays(key)
        assert len(rays) >= 4, (key, rays)
        tris = triangulate(rays)
        cones.append(dict(key=list(key), hits=int(hits),
                          rays=[list(r) for r in rays],
                          simplices=[list(s) for s in tris]))
    return cones

def coverage_check(cones, nsamp=20000, seed=5):
    """every random direction must land in exactly one simplicial subcone."""
    rng = np.random.default_rng(seed)
    W = rng.normal(size=(nsamp, 4))
    mats = []
    for ci, c in enumerate(cones):
        for s in c['simplices']:
            R = np.array([c['rays'][i] for i in s], dtype=float).T
            mats.append(np.linalg.inv(R))
    counts = np.zeros(nsamp, dtype=int)
    for Minv in mats:
        z = W @ Minv.T
        inside = (z > 1e-9).all(axis=1)
        counts += inside
    return dict(exactly1=int((counts == 1).sum()), zero=int((counts == 0).sum()),
                more=int((counts > 1).sum()), n=nsamp)

if __name__ == '__main__':
    cones = build_fan()
    nsub = sum(len(c['simplices']) for c in cones)
    print(f"maximal cones: {len(cones)}  simplicial subcones: {nsub}")
    cov = coverage_check(cones)
    print("coverage:", cov)
    dets = [abs(det4([c['rays'][i] for i in s]))
            for c in cones for s in c['simplices']]
    print("subcone |det| range:", min(dets), max(dets))
    with open('fan_x36.json', 'w') as f:
        json.dump(dict(cones=cones, coverage=cov), f)
    print("saved fan_x36.json")

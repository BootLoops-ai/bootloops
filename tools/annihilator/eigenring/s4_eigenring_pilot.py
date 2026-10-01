#!/usr/bin/env python3
# What this is: eigenring/hom-space right-factor splitter for linear ODE
# operators mod p (an open reimplementation of the van Hoeij / DFactor-style eigenring factorization).  Implements the
# eigenring E(L) and the hom-space Hom(D/D.A, D/D.B) over F_p(x) with the
# annihilator-kit primitives (ore_rdiv_modp, lclm_modp, flint nmod_mat
# nullspace), and extracts Jordan/log-entangled right factors that plain
# right-division cannot separate.
# Verification is COPAIR: every candidate operator is applied to independently
# generated series (never accepted on prime-agreement or division alone).
# Usage: s4_eigenring_pilot.py <p> <cache_L17_npz> <exact_inputs.json> <outdir> [--selftest] [--rung3]
# (`--selftest` is self-contained: toy LCLM, no input files needed.)
import sys, os, json, time, resource
import numpy as np, flint
from fractions import Fraction

# tools/ on sys.path so `annihilator` resolves when run as a script
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))))
from annihilator.ore.ore_rdiv_modp import ore_rightdiv_modp
from annihilator.ore.lclm_modp import op_to_nmod, op_from_nmod, op_mul_nm

T0 = time.time()
def log(m): print(f"[{time.time()-T0:8.1f}s] {m}", flush=True)
def rss_mb(): return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss // 1024

# ---------------- basic op helpers (int-list ops: op[i] = coeff list of D^i, low-to-high in x)
def npoly(c, p): return flint.nmod_poly([int(v) % p for v in c], p)
def op_nm(il, p): return [npoly(c, p) for c in il]
def op_il(nm):
    return [[int(c[i]) for i in range(c.degree() + 1)] if c.degree() >= 0 else [0] for c in nm]
def op_order(il):
    for i in range(len(il) - 1, -1, -1):
        if any(int(v) != 0 for v in il[i]): return i
    return -1
def op_deg(il):
    return max((len(c) - 1 for c in il if any(int(v) != 0 for v in c)), default=-1)
def rem_zero(L_il, R_il, p):
    Q, S, e, cont = ore_rightdiv_modp(L_il, R_il, p)
    return (not any(any(int(v) != 0 for v in s) for s in S)), Q
def op_content_normalize(il, p):
    """primitive + leading poly top-coefficient -> 1"""
    nm = op_nm(il, p)
    g = None
    for c in nm:
        if c != 0: g = c if g is None else g.gcd(c)
    if g is not None and g.degree() > 0:
        nm = [c // g for c in nm]
    k = op_order(op_il(nm))
    lead = nm[k]
    top = int(lead[lead.degree()])
    inv = pow(top, p - 2, p)
    nm = [c * inv for c in nm]
    return op_il(nm)

def op_shift(il, c, p):
    """coeff polys in x -> polys in t, x = t + c."""
    sub = npoly([c, 1], p)
    return [npoly(cf, p).compose(sub) for cf in il]
def op_unshift_il(il_t, c, p):
    sub = npoly([(-c) % p, 1], p)
    return op_il([npoly(cf, p).compose(sub) for cf in il_t])

# ---------------- series machinery at an ordinary point (t = x - c)
def ff_table(K, N, p):
    FF = np.zeros((K + 1, N), dtype=np.int64); FF[0, :] = 1
    j = np.arange(N, dtype=np.int64)
    for i in range(1, K + 1):
        FF[i] = FF[i - 1] * ((j - i + 1) % p) % p
    return FF

def sol_basis(op_t_nm, p, N):
    """Fundamental solutions of the operator (nmod_poly coeffs in t) at ordinary t=0.
    Returns C (K, N) int64 with y_j = t^j + O(t^K)."""
    K = len(op_t_nm) - 1
    while K >= 0 and op_t_nm[K].degree() < 0: K -= 1
    A = [np.array([int(op_t_nm[i][j]) for j in range(max(op_t_nm[i].degree() + 1, 1))],
                  dtype=np.int64) for i in range(K + 1)]
    assert int(A[K][0]) != 0, "not an ordinary point (lc vanishes at c)"
    FF = ff_table(K, N, p)
    C = np.zeros((K, N), dtype=np.int64)
    for j in range(K): C[j, j] = 1
    lc0 = int(A[K][0])
    for n in range(K, N):
        s = n - K
        acc = np.zeros(K, dtype=np.int64)
        for i in range(K + 1):
            di = len(A[i]) - 1
            jlo = max(0, s + i - di); jhi = s + i
            if i == K: jhi = min(jhi, n - 1)
            if jhi < jlo: continue
            w = A[i][s + i - jhi: s + i - jlo + 1][::-1].copy()
            v = (w * FF[i, jlo:jhi + 1]) % p
            acc = (acc + C[:, jlo:jhi + 1] @ v) % p
        piv = lc0 * int(FF[K, n]) % p
        C[:, n] = (-acc) % p * pow(piv, p - 2, p) % p
    return C

def series_derivs(C, p, kmax):
    """list [i] -> array (K, N) of i-th derivatives (tail zero-padded)."""
    K, N = C.shape
    out = [C.copy()]
    cur = C
    for i in range(1, kmax + 1):
        nxt = np.zeros_like(cur)
        mult = (np.arange(1, N, dtype=np.int64)) % p
        nxt[:, :N - 1] = cur[:, 1:] * mult[None, :] % p
        out.append(nxt); cur = nxt
    return out

def apply_op_series(op_t_nm, y, p, nvalid):
    """apply operator (nmod_poly coeffs in t) to series y (int64 array); returns np array,
    valid to nvalid - order - maxdeg margin handled by caller."""
    K = len(op_t_nm) - 1
    N = len(y)
    yp = flint.nmod_poly([int(v) % p for v in y], p)
    acc = flint.nmod_poly([0], p)
    cur = yp
    for i in range(K + 1):
        if op_t_nm[i].degree() >= 0:
            acc = acc + op_t_nm[i] * cur
        cur = cur.derivative()
    out = np.array([int(acc[i]) for i in range(N)], dtype=np.int64)
    return out

# ---------------- linear algebra
def nullspace_modp(M2d, p):
    nrow, ncol = M2d.shape
    ent = (M2d % p).flatten().tolist()
    Mf = flint.nmod_mat(nrow, ncol, ent, p)
    V, nullity = Mf.nullspace()
    del Mf, ent
    out = []
    for k in range(nullity):
        out.append(np.array([int(V[i, k]) % p for i in range(ncol)], dtype=np.int64))
    return out

def mat_nullspace_small(M, p):
    return nullspace_modp(np.array(M, dtype=np.int64) % p, p)

# ---------------- hom space / eigenring
def hom_space(src_derivs, tgt, Dpoly_t, d, p, margin=250):
    """Solve for R = sum_i n_i(t) D^i (i < Ks=len(src_derivs)) with deg n_i <= d and
    R(src_j) in span(tgt) after clearing denominator Dpoly_t:
        sum_i n_i * y_j^(i)  =  sum_k m_kj * (Dpoly * tgt_k)   (as series)
    Returns (kernel vectors, layout dict). Kernel vector layout:
    [ n_0 coeffs (d+1) | n_1 ... | n_{Ks-1} | m_{k,j} col-major j-blocks (Kt per j) ]"""
    Ks = len(src_derivs); Kt = tgt.shape[0]; N = tgt.shape[1]
    ncol_n = Ks * (d + 1); ncol = ncol_n + Kt * src_derivs[0].shape[0]
    nsrc = src_derivs[0].shape[0]
    ncol = ncol_n + Kt * nsrc
    Nr = (ncol + margin + nsrc - 1) // nsrc
    Nr = min(Nr, N - 5)
    assert nsrc * Nr > ncol + margin // 2, f"series too short: Nr={Nr} ncol={ncol}"
    Dt = np.array([int(Dpoly_t[i]) for i in range(Dpoly_t.degree() + 1)], dtype=np.int64)
    # premultiply targets by D
    Dtgt = np.zeros((Kt, Nr), dtype=np.int64)
    for k in range(Kt):
        prod = npoly(Dt.tolist(), p) * npoly(tgt[k].tolist(), p)
        Dtgt[k] = np.array([int(prod[i]) for i in range(Nr)], dtype=np.int64)
    M = np.zeros((nsrc * Nr, ncol), dtype=np.int64)
    for a in range(nsrc):
        r0 = a * Nr
        for i in range(Ks):
            ya = src_derivs[i][a]
            base = i * (d + 1)
            for mu in range(d + 1):
                ln = Nr - mu
                if ln <= 0: break
                M[r0 + mu: r0 + Nr, base + mu] = ya[:ln]
        for k in range(Kt):
            M[r0: r0 + Nr, ncol_n + a * Kt + k] = (p - Dtgt[k]) % p
    kern = nullspace_modp(M, p)
    return kern, {'Ks': Ks, 'Kt': Kt, 'd': d, 'Nr': Nr, 'ncol': ncol, 'nsrc': nsrc}

def hom_residual_check(vec, lay, src_derivs_long, tgt_long, Dpoly_t, p, Ncheck):
    """recompute the defining relation on longer series; return # nonzero coeffs in
    rows [Nr, Ncheck) summed over sources (0 = clean holdout)."""
    Ks, Kt, d, Nr, nsrc = lay['Ks'], lay['Kt'], lay['d'], lay['Nr'], lay['nsrc']
    ncol_n = Ks * (d + 1)
    bad = 0
    Dt = Dpoly_t
    for a in range(nsrc):
        acc = flint.nmod_poly([0], p)
        for i in range(Ks):
            ni = npoly(vec[i * (d + 1):(i + 1) * (d + 1)].tolist(), p)
            acc = acc + ni * npoly(src_derivs_long[i][a].tolist(), p)
        for k in range(Kt):
            mk = int(vec[ncol_n + a * Kt + k])
            if mk: acc = acc - npoly([mk], p) * Dt * npoly(tgt_long[k].tolist(), p)
        for r in range(Nr, Ncheck):
            if int(acc[r]) != 0: bad += 1
    return bad

def action_matrix(vec, lay, Dpoly0inv=None):
    """m_kj matrix (Kt x nsrc) from kernel vector (eigenring: Kt = nsrc = K)."""
    Ks, Kt, d, nsrc = lay['Ks'], lay['Kt'], lay['d'], lay['nsrc']
    ncol_n = Ks * (d + 1)
    M = np.zeros((Kt, nsrc), dtype=np.int64)
    for a in range(nsrc):
        for k in range(Kt):
            M[k, a] = int(vec[ncol_n + a * Kt + k])
    return M

# ---------------- multi-series operator fit
def fit_multi(series_list, order, deg, p, holdout=40):
    """Fit sum_i b_i(t) y^(i) = 0 (deg b_i <= deg) to ALL series simultaneously.
    Returns (op_il_t, nullity, ho_bad) or None."""
    ncol = (order + 1) * (deg + 1)
    Nmin = min(len(s) for s in series_list)
    nrows_per = (ncol + 200 + len(series_list) - 1) // len(series_list)
    Nr = min(nrows_per, Nmin - order - holdout - 2)
    if Nr * len(series_list) < ncol + 30: return None
    blocks = []
    for s in series_list:
        D = series_derivs(np.array([s], dtype=np.int64), p, order)
        Mb = np.zeros((Nr, ncol), dtype=np.int64)
        for i in range(order + 1):
            base = i * (deg + 1)
            ya = D[i][0]
            for mu in range(deg + 1):
                ln = Nr - mu
                if ln <= 0: break
                Mb[mu:Nr, base + mu] = ya[:ln]
        blocks.append((Mb, D))
    M = np.vstack([b[0] for b in blocks])
    kern = nullspace_modp(M, p)
    if not kern: return None
    vec = kern[0]
    opil = [vec[i * (deg + 1):(i + 1) * (deg + 1)].tolist() for i in range(order + 1)]
    # holdout: apply to each series on rows [Nr, Nmin-order-2)
    bad = 0
    op_nm_t = op_nm(opil, p)
    for s in series_list:
        res = apply_op_series(op_nm_t, np.array(s, dtype=np.int64), p, len(s))
        for r in range(Nr, Nmin - order - 8):
            if int(res[r]) != 0: bad += 1
    return opil, len(kern), bad

def fit_multi_ramp(series_list, order, degs, p):
    for dg in degs:
        r = fit_multi(series_list, order, dg, p)
        if r is not None and r[2] == 0 and r[1] >= 1:
            return r[0], dg, r[1]
    return None, None, None

def min_deg_fit(series_list, order, p, dmax=1100, coarse=120, verbose=None):
    """coarse ramp then bisect to the MINIMAL clean-fit degree (for CRT sizing)."""
    found = None; prev_fail = 0
    for dg in range(60, dmax + 1, coarse):
        r = fit_multi(series_list, order, dg, p)
        ok = (r is not None and r[2] == 0 and r[1] >= 1)
        if verbose: verbose(f"  min_deg_fit coarse d={dg}: {'OK nul=%d' % r[1] if ok else 'fail'}")
        if ok: found = (dg, r); break
        prev_fail = dg
    if found is None: return None, None, None
    lo, hi = prev_fail, found[0]; best = found
    while hi - lo > 1:
        mid = (lo + hi) // 2
        r = fit_multi(series_list, order, mid, p)
        ok = (r is not None and r[2] == 0 and r[1] >= 1)
        if verbose: verbose(f"  min_deg_fit bisect d={mid}: {'OK' if ok else 'fail'}")
        if ok: hi = mid; best = (mid, r)
        else: lo = mid
    return best[1][0], best[0], best[1][1]

# ---------------- reference hypergeometric solution series, pure flint/Fraction
def hyp_series_x2(a, b, c, p, N):
    coeffs = [0] * N
    t = 1; coeffs[0] = 1 % p; k = 1
    while 2 * k < N:
        r = Fraction((a.numerator + (k - 1) * a.denominator) * (b.numerator + (k - 1) * b.denominator) * a.denominator * b.denominator, 1)
        # careful exact: t_k = t_{k-1} * (a+k-1)(b+k-1)/((c+k-1) k)
        num = (a + k - 1) * (b + k - 1)
        den = (c + k - 1) * k
        rr = num / den
        rp = (rr.numerator % p) * pow(rr.denominator % p, p - 2, p) % p
        t = t * rp % p
        coeffs[2 * k] = t
        k += 1
    return flint.nmod_poly(coeffs, p)

def build_z(p, N):
    F0 = hyp_series_x2(Fraction(1, 2), Fraction(1, 2), Fraction(1), p, N)
    F1 = hyp_series_x2(Fraction(1, 2), Fraction(3, 2), Fraction(2), p, N)
    def trunc(g): return flint.nmod_poly([int(g[i]) for i in range(min(N, g.degree() + 1))], p)
    F0_2 = trunc(F0 * F0); F0_3 = trunc(F0_2 * F0)
    F1_2 = trunc(F1 * F1); F1_3 = trunc(F1_2 * F1)
    F1F0_2 = trunc(F1 * F0_2); F1_2F0 = trunc(F1_2 * F0)
    P = lambda c: npoly(c, p)
    inner = trunc(P([0,0,0,0,0,3]) * F1_3) + trunc(P([16,16,8,16,16]) * F0_3) \
          + trunc(P([0,20,4,-12,4,20]) * F1F0_2) + trunc(P([0,0,4,-20,-34,-20,4]) * F1_2F0)
    inner = trunc(inner)
    D = P([8]) * P([-1,1])**2 * P([1,1])**2 * P([1,1,1])
    z = trunc(trunc(inner) * D.inverse_series_trunc(N))
    return np.array([int(z[i]) for i in range(N)], dtype=np.int64)

# ---------------- invariant-subspace refinement over F_p
def subspace_rank(vecs, p):
    if not vecs: return 0
    M = np.array(vecs, dtype=np.int64) % p
    Mf = flint.nmod_mat(M.shape[0], M.shape[1], M.flatten().tolist(), p)
    return Mf.rank()

def rref_basis(vecs, p):
    """row-reduce; return (rows (k,n) int64 in rref, pivot column list)."""
    M = np.array(vecs, dtype=np.int64) % p
    Mf = flint.nmod_mat(M.shape[0], M.shape[1], M.flatten().tolist(), p)
    Rf, rank = Mf.rref()
    rows = []
    pivs = []
    for i in range(rank):
        row = np.array([int(Rf[i, j]) % p for j in range(M.shape[1])], dtype=np.int64)
        nz = np.nonzero(row)[0]
        rows.append(row); pivs.append(int(nz[0]))
    return np.array(rows, dtype=np.int64), pivs

def restrict_action(Bs, pivs, Lam, p):
    """Bs (k,n) rref rows spanning an invariant subspace; return (Rm (k,k), invariant?)."""
    k = Bs.shape[0]
    W = (Lam % p) @ Bs.T % p          # (n,k): images as columns
    Rm = W[pivs, :] % p               # coords via pivot entries
    recon = Bs.T @ Rm % p
    inv = np.array_equal(recon % p, W % p)
    return Rm.T % p, inv              # Rm.T rows: image coords? keep consistent below

def split_by_action(Bs, pivs, Lam, p):
    """Split invariant subspace (rref rows Bs) by generalized eigenspaces of the
    restricted action. Returns list of (rows, pivs) sub-blocks (possibly [original])."""
    k = Bs.shape[0]
    W = (Lam % p) @ Bs.T % p
    Rm = W[pivs, :] % p               # (k,k): column i = coords of Lam@Bs_i
    recon = Bs.T @ Rm % p
    if not np.array_equal(recon % p, W % p):
        return None                   # not invariant under this action
    Mf = flint.nmod_mat(k, k, (Rm % p).flatten().tolist(), p)
    cp = Mf.charpoly()
    fac = cp.factor()
    parts = []
    for (f, mult) in fac[1]:
        # generalized kernel of f(Rm)^mult
        fR = np.zeros((k, k), dtype=np.int64)
        # evaluate f at Rm: f = sum f_j x^j
        P = np.eye(k, dtype=np.int64)
        for j in range(f.degree() + 1):
            fj = int(f[j]) % p
            if fj: fR = (fR + fj * P) % p
            P = P @ Rm % p
        fRk = np.eye(k, dtype=np.int64)
        for _ in range(int(mult)):
            fRk = fRk @ fR % p
        kv = mat_nullspace_small(fRk, p)
        if kv:
            sub = np.array(kv, dtype=np.int64) % p @ Bs % p   # lift coords to ambient
            parts.append(rref_basis(sub, p))
    if sum(b.shape[0] for b, _ in parts) != k:
        return None
    if len(parts) <= 1:
        return [ (Bs, pivs) ]
    return parts

def refine_blocks(actions, n, p, max_rounds=6, seed=12345):
    """Finest joint invariant decomposition obtainable from the action list
    (basis elements + deterministic pseudo-random combos)."""
    import random
    rng = random.Random(seed)
    ops = [a % p for a in actions]
    for _ in range(8):
        combo = sum(rng.randrange(1, p) * a for a in actions) % p
        ops.append(combo)
    blocks = [rref_basis(np.eye(n, dtype=np.int64), p)]
    for _ in range(max_rounds):
        changed = False
        for Lam in ops:
            newblocks = []
            for (Bs, pivs) in blocks:
                if Bs.shape[0] == 1:
                    newblocks.append((Bs, pivs)); continue
                parts = split_by_action(Bs, pivs, Lam, p)
                if parts is None or len(parts) == 1:
                    newblocks.append((Bs, pivs))
                else:
                    newblocks.extend(parts); changed = True
            blocks = newblocks
        if not changed: break
    return blocks

# =====================================================================
def main():
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    flags = {a for a in sys.argv[1:] if a.startswith('--')}
    if '--selftest' in flags:
        selftest(); return
    p = int(args[0]); cachef = args[1]; inputsf = args[2]; outdir = args[3]
    os.makedirs(outdir, exist_ok=True)
    R = {'p': p, 'stamp_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
         'walls': {}, 'producer': 'annihilator.eigenring s4_eigenring_pilot.py'}
    tw = lambda k, t: R['walls'].__setitem__(k, round(time.time() - t, 2))

    # ---- stage 0: inputs
    t = time.time()
    ex = json.load(open(inputsf))
    cache = np.load(cachef)
    L17 = [cache['C17'][i].tolist() for i in range(cache['C17'].shape[0])]
    U1 = ex['U1']; Vt = ex['V_true']; Wt = ex['W_true']; Vp = ex['V_paper']
    L1 = ex['L1']; L3 = ex['L3']; L4p = ex['L4_paper']
    VU = op_from_nmod(op_mul_nm(op_to_nmod(Vt, p), op_to_nmod(U1, p), p))
    WU = op_from_nmod(op_mul_nm(op_to_nmod(Wt, p), op_to_nmod(U1, p), p))
    VpU = op_from_nmod(op_mul_nm(op_to_nmod(Vp, p), op_to_nmod(U1, p), p))
    z1, A15 = rem_zero(L17, VU, p)
    assert z1, "VU does not right-divide L17"
    R['A15'] = {'order': op_order(A15), 'deg': op_deg(A15)}
    log(f"A15: order {op_order(A15)} deg {op_deg(A15)}")
    tw('stage0_inputs', t)

    # ---- stage 1: reference series z, L4z fit, decisive gates
    t = time.time()
    Nz = 500
    z = build_z(p, Nz)
    L4z, dgz, nulz = fit_multi_ramp([z.tolist()], 4, list(range(20, 60, 2)), p)
    assert L4z is not None, "L4z fit failed"
    L4z = op_content_normalize(L4z, p)
    R['L4z'] = {'deg': op_deg(L4z), 'order': op_order(L4z)}
    g_a15_l4z, _ = rem_zero(A15, L4z, p)
    g_a15_l4p, _ = rem_zero(A15, L4p, p)
    a15z = apply_op_series(op_nm(A15, p), z, p, Nz)
    z_in_kerA15 = not any(int(v) != 0 for v in a15z[:Nz - 260])
    R['gates_stage1'] = {'rem(A15,L4z)==0': bool(g_a15_l4z),
                        'rem(A15,L4paper)==0': bool(g_a15_l4p),
                        'A15(z)==0': bool(z_in_kerA15)}
    log(f"stage1: L4z deg={op_deg(L4z)}; rem(A15,L4z)=0:{g_a15_l4z} "
        f"rem(A15,L4p)=0:{g_a15_l4p} A15(z)=0:{z_in_kerA15}")
    tw('stage1_L4z', t)

    # ---- stage 2: solution bases at ordinary point c
    t = time.time()
    c = 7
    lcA15 = npoly(A15[op_order(A15)], p)
    assert int(lcA15(flint.nmod(c, p))) != 0
    A15t = op_shift(A15, c, p)
    d_r1 = op_deg(A15) + 40
    Nlong = d_r1 + 40 + 700
    Y = sol_basis(A15t, p, Nlong)
    K = Y.shape[0]
    Yd = series_derivs(Y, p, K - 1)
    # internal consistency: A15 applied to first and last basis solution -> 0
    for jchk in (0, K - 1):
        r0 = apply_op_series(A15t, Y[jchk], p, Nlong)
        assert not any(int(v) != 0 for v in r0[:Nlong - K - 8]), "sol_basis inconsistent"
    tw('stage2_solbasis', t)
    log(f"sol basis A15 at c={c}: {K} sols x {Nlong} terms (consistency OK)")

    # known-block images: VU(ker L1), VU(ker L3), VU(ker Wt.U1)
    t = time.time()
    VUt = op_shift(VU, c, p)
    known_vecs = []; known_labels = []
    for (opk, lbl) in ((L1, 'L1'), (L3, 'L3'), (WU, 'WtU1')):
        opkt = op_shift(opk, c, p)
        Sk = sol_basis(opkt, p, Nlong)
        for r_ in range(Sk.shape[0]):
            w = apply_op_series(VUt, Sk[r_], p, Nlong)
            if not any(int(v) != 0 for v in w[:Nlong - 300]):
                continue  # killed by VU
            coords = w[:K].copy()
            resid = (w - coords @ Y % p) % p
            ok = not any(int(v) != 0 for v in resid[:Nlong - 320])
            known_vecs.append(coords.tolist()); known_labels.append((lbl, bool(ok)))
    R['known_block_images'] = {'labels': known_labels,
                               'rank': subspace_rank(known_vecs, p)}
    log(f"known-block image vectors: {known_labels} rank={R['known_block_images']['rank']}")
    tw('stage2b_knownblocks', t)

    # ---- stage 3: eigenring E(A15), rung 1 (denominator = lc)
    t = time.time()
    Dden = lcA15.compose(npoly([c, 1], p))
    kern, lay = hom_space(Yd, Y, Dden, d_r1, p)
    R['E_A15_rung1'] = {'d': d_r1, 'dim': len(kern), 'Nr': lay['Nr'], 'ncol': lay['ncol'],
                        'wall_s': round(time.time() - t, 2)}
    log(f"E(A15) rung1: dim={len(kern)} (d={d_r1})  [{R['E_A15_rung1']['wall_s']}s]")
    hol = [hom_residual_check(v, lay, Yd, Y, Dden, p, min(lay['Nr'] + 250, Nlong - K - 5))
           for v in kern]
    R['E_A15_rung1']['holdout_bad'] = hol
    tw('stage3_EA15_r1', t)

    # rung 2 if trivial
    if len([h for h in hol if h == 0]) <= 1:
        t = time.time()
        sing = npoly([0, 1], p) * npoly([-1, 1], p) * npoly([1, 1], p) * npoly([1, 1, 1], p) * npoly([1, 1, 1, 1, 1], p)
        Dden2 = (lcA15 * sing ** 2).compose(npoly([c, 1], p))
        d_r2 = op_deg(A15) + 20 + 40 + 60
        if d_r2 + 40 + 700 > Nlong:
            Y2 = sol_basis(A15t, p, d_r2 + 40 + 700); Yd2 = series_derivs(Y2, p, K - 1)
        else:
            Y2, Yd2 = Y, Yd
        kern2, lay2 = hom_space(Yd2, Y2, Dden2, d_r2, p)
        hol2 = [hom_residual_check(v, lay2, Yd2, Y2, Dden2, p,
                                   min(lay2['Nr'] + 250, Y2.shape[1] - K - 5)) for v in kern2]
        R['E_A15_rung2'] = {'d': d_r2, 'dim': len(kern2), 'holdout_bad': hol2,
                            'wall_s': round(time.time() - t, 2)}
        log(f"E(A15) rung2: dim={len(kern2)} holdout={hol2} [{R['E_A15_rung2']['wall_s']}s]")
        if len([h for h in hol2 if h == 0]) > 1:
            kern, lay, Dden, Y, Yd = kern2, lay2, Dden2, Y2, Yd2
            hol = hol2

    # ---- stage 4: known-block peel: A5 fit, Q10 = A15 // A5, E(Q10)
    t = time.time()
    kseries = []
    for (opk, lbl) in ((L1, 'L1'), (L3, 'L3'), (WU, 'WtU1')):
        opkt = op_shift(opk, c, p)
        Sk = sol_basis(opkt, p, Nlong)
        for r_ in range(Sk.shape[0]):
            w = apply_op_series(VUt, Sk[r_], p, Nlong - 60)
            if any(int(v) != 0 for v in w[:Nlong - 300]):
                kseries.append(w[:Nlong - 60].tolist())
    A5t, dgA5, nulA5 = fit_multi_ramp(kseries, 5, list(range(20, 320, 20)), p)
    R['A5'] = {'deg_t': None if A5t is None else op_deg(A5t), 'fit_deg': dgA5}
    Q10 = None
    if A5t is not None:
        A5 = op_content_normalize(op_unshift_il(A5t, c, p), p)
        gA5, Q10c = rem_zero(A15, A5, p)
        R['A5']['rem(A15,A5)==0'] = bool(gA5)
        R['A5']['deg_x'] = op_deg(A5)
        if gA5:
            Q10 = Q10c
            R['Q10'] = {'order': op_order(Q10), 'deg': op_deg(Q10)}
    log(f"A5 fit deg={R['A5']} Q10={R.get('Q10')}")
    tw('stage4_A5_Q10', t)

    if Q10 is not None:
        t = time.time()
        lcQ = npoly(Q10[op_order(Q10)], p)
        if op_deg(Q10) > 650:
            R['E_Q10'] = {'skipped': f'Q10 deg {op_deg(Q10)} > 650 cap (matrix too large for pilot)'}
        elif int(lcQ(flint.nmod(c, p))) != 0:
            Q10t = op_shift(Q10, c, p)
            dq = op_deg(Q10) + 40
            Nq = dq + 40 + 500
            YQ = sol_basis(Q10t, p, Nq)
            YQd = series_derivs(YQ, p, YQ.shape[0] - 1)
            DdenQ = lcQ.compose(npoly([c, 1], p))
            kq, layq = hom_space(YQd, YQ, DdenQ, dq, p)
            holq = [hom_residual_check(v, layq, YQd, YQ, DdenQ, p,
                                       min(layq['Nr'] + 200, Nq - 15)) for v in kq]
            R['E_Q10'] = {'d': dq, 'dim': len(kq), 'holdout_bad': holq,
                          'wall_s': round(time.time() - t, 2)}
            log(f"E(Q10): dim={len(kq)} holdout={holq} [{R['E_Q10']['wall_s']}s]")
        else:
            R['E_Q10'] = {'skipped': 'lcQ vanishes at c'}
        tw('stage4b_EQ10', t)

    # ---- stage 5: hom probe Hom(D/D.A15, D/D.L4z): R order<=3, R(ker L4z) subset ker A15
    t = time.time()
    lcz = npoly(L4z[op_order(L4z)], p)
    if int(lcz(flint.nmod(c, p))) == 0:
        R['hom_L4z'] = {'skipped': 'lc L4z vanishes at c'}
    else:
        L4zt = op_shift(L4z, c, p)
        Yz = sol_basis(L4zt, p, Nlong)
        Yzd = series_derivs(Yz, p, 3)
        Dh = (lcA15 * lcz).compose(npoly([c, 1], p))
        dh = op_deg(A15) + op_deg(L4z) + 60
        kh, layh = hom_space(Yzd, Y, Dh, dh, p)
        holh = [hom_residual_check(v, layh, Yzd, Y, Dh, p,
                                   min(layh['Nr'] + 250, Nlong - K - 5)) for v in kh]
        R['hom_L4z'] = {'d': dh, 'dim': len(kh), 'holdout_bad': holh,
                        'wall_s': round(time.time() - t, 2)}
        log(f"Hom(A15->L4z-type): dim={len(kh)} holdout={holh}")
    tw('stage5_hom', t)

    # ---- stage 6: T4 extraction from whichever structure appeared
    t = time.time()
    T4 = None; T4_route = None
    NL2 = 2600
    Ylong = None; U4_eigen = None; U4_hom = None
    kv = np.array(known_vecs, dtype=np.int64)
    kvrank = subspace_rank(known_vecs, p)
    clean = [i for i, h in enumerate(hol) if h == 0]
    if len(clean) > 1:
        # eigenring route on A15: refine ker A15 into joint invariant blocks
        acts = [action_matrix(kern[i], lay) for i in clean]
        blocks = refine_blocks(acts, K, p)
        R['eigen_A15_blocks'] = sorted(b.shape[0] for b, _ in blocks)
        overlaps = []
        comp = []
        for (Bs, _) in blocks:
            u = subspace_rank(np.vstack([kv, Bs]).tolist(), p)
            ov = kvrank + Bs.shape[0] - u
            overlaps.append({'dim': int(Bs.shape[0]), 'overlap_known5': int(ov)})
            if ov == 0:
                comp.extend(Bs.tolist())
        R['eigen_A15_overlaps'] = overlaps
        R['eigen_A15_compdim'] = subspace_rank(comp, p)
        log(f"eigenring blocks: {overlaps} compdim={R['eigen_A15_compdim']}")
        if R['eigen_A15_compdim'] == 4:
            U4_eigen = rref_basis(np.array(comp, dtype=np.int64), p)[0]
    # hom-route U4 image
    if 'hom_L4z' in R and R['hom_L4z'].get('dim', 0) > 0:
        cleanh = [i for i, h in enumerate(R['hom_L4z']['holdout_bad']) if h == 0]
        if cleanh:
            v = kh[cleanh[0]]
            # image coords in Y basis: R(y_{z,a}) = sum_k m_{k,a} y_k  (m-part of vector)
            Mimg = action_matrix(v, layh)          # (Kt=15, nsrc=4)
            imgs = [Mimg[:, a] % p for a in range(layh['nsrc'])]
            U4h = rref_basis(np.array(imgs, dtype=np.int64), p)[0]
            R['hom_L4z']['image_rank'] = int(U4h.shape[0])
            ovh = kvrank + U4h.shape[0] - subspace_rank(np.vstack([kv, U4h]).tolist(), p)
            R['hom_L4z']['image_overlap_known5'] = int(ovh)
            if U4h.shape[0] == 4:
                U4_hom = U4h
    if U4_eigen is not None and U4_hom is not None:
        R['U4_eigen_eq_hom'] = bool(subspace_rank(np.vstack([U4_eigen, U4_hom]).tolist(), p) == 4)
    # pick U4 and fit T4 at minimal degree on LONG series
    U4 = U4_hom if U4_hom is not None else U4_eigen
    if U4 is None and U4_eigen is None and 'eigen_A15_overlaps' in R:
        # fall back: any 4-dim block even if overlapping (diagnostic fit)
        pass
    if U4 is not None:
        log("generating long solution basis for minimal-degree T4 fit ...")
        Ylong = sol_basis(A15t, p, NL2)
        ws = (U4 % p) @ Ylong % p
        T4t, dgt, nult = min_deg_fit([w.tolist() for w in ws], 4, p, dmax=1100,
                                     coarse=120, verbose=log)
        R['T4_fit'] = {'min_deg_t': dgt, 'nullity': nult}
        if T4t is not None:
            T4 = op_content_normalize(op_unshift_il(T4t, c, p), p)
            T4_route = 'hom_L4z' if U4 is U4_hom else 'eigenring_A15'
    R['T4_route'] = T4_route
    tw('stage6_extract', t)

    # ---- stage 7: COPAIR verification of T4
    if T4 is not None:
        t = time.time()
        R['T4'] = {'order': op_order(T4), 'deg': op_deg(T4)}
        gate_div, B11 = rem_zero(A15, T4, p)
        R['T4']['rem(A15,T4)==0'] = bool(gate_div)
        R['T4']['B11_order'] = op_order(B11)
        lcT = npoly(T4[op_order(T4)], p)
        c2 = None
        for cc in (11, 13, 17, 19, 23):
            if int(lcT(flint.nmod(cc, p))) != 0 and int(lcA15(flint.nmod(cc, p))) != 0:
                c2 = cc; break
        if c2 is not None:
            T4c2 = op_shift(T4, c2, p)
            N3 = 400 + op_deg(A15)
            S4 = sol_basis(T4c2, p, N3)
            A15c2 = op_shift(A15, c2, p)
            bad = 0; checked = 0
            for a in range(S4.shape[0]):
                res = apply_op_series(A15c2, S4[a], p, N3)
                for r_ in range(N3 - op_deg(A15) - 20):
                    checked += 1
                    if int(res[r_]) != 0: bad += 1
            R['T4']['COPAIR_A15_kills_kerT4'] = {'bad': bad, 'checked': checked}
            # chain gate: P6 = T4.VU ; L17 applied to ker P6
            P6 = op_from_nmod(op_mul_nm(op_to_nmod(T4, p), op_to_nmod(VU, p), p))
            lcP = npoly(P6[op_order(P6)], p)
            if int(lcP(flint.nmod(c2, p))) != 0:
                S6 = sol_basis(op_shift(P6, c2, p), p, N3)
                L17c2 = op_shift(L17, c2, p)
                bad6 = 0; chk6 = 0
                for a in range(S6.shape[0]):
                    res = apply_op_series(L17c2, S6[a], p, N3)
                    for r_ in range(N3 - op_deg(A15) - 40):
                        chk6 += 1
                        if int(res[r_]) != 0: bad6 += 1
                R['T4']['COPAIR_L17_kills_kerT4VU'] = {'bad': bad6, 'checked': chk6}
            # reference-series relation
            t4z = apply_op_series(op_nm(T4, p), z, p, Nz)
            R['T4']['T4(z)==0'] = not any(int(v) != 0 for v in t4z[:Nz - op_deg(T4) - 20])
        json.dump({'p': p, 'T4': T4, 'route': T4_route},
                  open(os.path.join(outdir, f'T4_p{p}.json'), 'w'))
        tw('stage7_copair', t)
        log(f"T4: {R['T4']}")

    R['peak_rss_mb'] = rss_mb()
    R['total_wall_s'] = round(time.time() - T0, 2)
    json.dump(R, open(os.path.join(outdir, f'receipt_p{p}.json'), 'w'), indent=1)
    log(f"DONE p={p} wall={R['total_wall_s']}s rss={R['peak_rss_mb']}MB")

# =====================================================================
def selftest():
    """Toy: L = LCLM of (1-x)D-1 [sol 1/(1-x)] and (1-2x)D-2 [sol 1/(1-2x)].
    E(L) must be 2-dim; extraction must recover an order-1 right factor that
    COPAIR-verifies. Exercises sol_basis, hom_space, action, fits, unshift, gates."""
    p = 1048573
    A = [[-1], [1, -1]]           # (1-x)D - 1
    B = [[-2], [1, -2]]           # (1-2x)D - 2
    # LCLM via fit on the two sols
    N = 200
    c = 5
    At = op_shift(A, c, p); Bt = op_shift(B, c, p)
    sA = sol_basis(At, p, N); sB = sol_basis(Bt, p, N)
    L2t, dg, nul = fit_multi_ramp([sA[0].tolist(), sB[0].tolist()], 2, [1, 2, 3, 4], p)
    assert L2t is not None, "selftest LCLM fit failed"
    L2 = op_content_normalize(op_unshift_il(L2t, c, p), p)
    print("L2 (toy LCLM):", L2)
    L2t2 = op_shift(L2, c, p)
    Y = sol_basis(L2t2, p, N)
    Yd = series_derivs(Y, p, 1)
    lc = npoly(L2[op_order(L2)], p).compose(npoly([c, 1], p))
    kern, lay = hom_space(Yd, Y, lc, op_deg(L2) + 6, p, margin=60)
    hol = [hom_residual_check(v, lay, Yd, Y, lc, p, min(lay['Nr'] + 60, N - 10)) for v in kern]
    print("E(toy) dim:", len(kern), "holdout:", hol)
    # both toy modules are trivial-type (rational solutions), so E is the full
    # 2x2 matrix algebra: dim 4. Require >= 2 (non-scalar structure present).
    assert len([h for h in hol if h == 0]) >= 2, "toy eigenring should be non-scalar"
    clean = [i for i, h in enumerate(hol) if h == 0]
    acts = [action_matrix(kern[i], lay) for i in clean]
    blocks = refine_blocks(acts, 2, p)
    print("blocks:", [b.shape[0] for b, _ in blocks])
    got = False
    for (Bs, _) in blocks:
        if Bs.shape[0] == 1:
            w = (Bs[0] % p) @ Y % p
            T1t, dgt, _ = fit_multi_ramp([w.tolist()], 1, [1, 2, 3], p)
            if T1t is None: continue
            T1 = op_content_normalize(op_unshift_il(T1t, c, p), p)
            ok, _ = rem_zero(L2, T1, p)
            # COPAIR: sols of T1 at c2=3 killed by L2
            S1 = sol_basis(op_shift(T1, 3, p), p, 80)
            res = apply_op_series(op_shift(L2, 3, p), S1[0], p, 80)
            copair = not any(int(v) != 0 for v in res[:60])
            print("  factor:", T1, "rem0:", ok, "copair:", copair)
            got = got or (ok and copair)
    assert got, "toy extraction failed"
    print("SELFTEST PASS")

if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""branches_mod.py — Frobenius branch series via the D_p-block recursion.

Design rules:
 (1) sign convention (README.md, section 'Convention'): lambda_j = eig(D1),
     branch  N^(j)(u) = u^{-lambda_j} * sum_n N_n u^n  (u = 1/x, N ~ x^lambda),
     so the order-n matrix is  D1 + (n - lambda_j) I  (seed: (D1-lambda_j)v=0);
 (2) the n=0 consistency equation's RHS is -(D1 - lambda_j) N_0 (=0 up to seed
     residual), never -D1*N_0 (dropping the lambda term is WRONG).
Resonance handling: orders n where lambda_j - n collides
with another eigenvalue make D1+(n-lam)I singular; the pure-backward path
solves the compatible/min-norm system there (SVD projection of the RHS onto
the range; the numerically-asserted incompatibility is reported as
refine_rfinal and, if large, flags LOG_BRANCH_REQUIRED — the standard
Frobenius prescription: homogeneous freedom stays with the resonant
eigenvalue's own seeds, span of Phi unchanged).
Forward blocks: fit-noise blocks (~1e-120) are dropped; a
GENUINE nilpotent D_0 block (rank-1 irregular-looking N-frame)
is removed by an exact shear reduction (u^-1 scalings on range(B0)
subspaces, Moser-style) and the branch is built by the pure-backward
recursion in the Fuchsian reduced frame, then mapped back; a
forward-coupled dense lstsq is kept only as a fallback (diag['reduce_failed'])
for structures the reduction does not cover (fdepth > 1).
Injective seed->column map: callers pass the requested seed rank as
seed_index (injective by construction; cli.py enumerates the seed basis).
Argmax-overlap recovery against a freshly recomputed nullspace basis of
(D1 - lam I) is NOT injective on degenerate families (the two SVD bases
differ by an internal rotation and the argmax can collide, returning the
same orbit column twice and dropping the Phi rank); _seed_col survives only
as the legacy fallback for direct API callers that do not pass seed_index.
diag['seed_col_used'] records the column actually fetched.
"""
import mpmath as mp
mp.mp.dps = max(mp.mp.dps, 50)      # set FIRST
import numpy as np


def resonance_orders(lam, eigvals, n_ord, tol=1e-8):
    """Orders n>=1 with lambda - n ~ another eigenvalue (descending ladder)."""
    out = []
    for n in range(1, n_ord + 1):
        for z in eigvals:
            if abs(z - (lam - n)) < tol:
                out.append({'n': n, 'collides_with': mp.nstr(z, 20)})
    return out


def _matnorm(M):
    return max(abs(M[i, j]) for i in range(M.rows) for j in range(M.cols))


def _minnorm_solve(Mt, rhs, prec):
    """Compatible/min-norm solve of a (numerically) singular square system:
    project rhs onto range(Mt) via SVD, return the pseudo-inverse solution and
    the incompatibility |P_null rhs| (must be ~0 for a true resonant branch —
    the numerical assertion of Frobenius compatibility)."""
    n = Mt.rows
    U, S, V = mp.svd(mp.matrix(Mt))              # Mt = U * diag(S) * V
    smax = S[0] if S[0] > 0 else mp.mpf(1)
    tol = smax * mp.mpf(10) ** (-(prec // 2))
    x = mp.matrix(n, 1)
    incomp = mp.mpf(0)
    for i in range(n):
        ui = sum(mp.conj(U[k, i]) * rhs[k, 0] for k in range(n))
        if S[i] > tol:
            c = ui / S[i]
            for k in range(n):
                x[k, 0] += mp.conj(V[i, k]) * c
        else:
            incomp = max(incomp, abs(ui))
    return x, incomp


def _backward(Dp, N, lam, seed, n_ord, prec, diag, eigvals):
    """Pure backward recursion (no forward blocks): at resonant orders solve
    the compatible/min-norm system instead of lu_solve."""
    Nn = {0: mp.matrix(seed)}
    D1 = Dp[1]
    res_n = set()
    if eigvals is not None:
        for n in range(1, n_ord + 1):
            for z in eigvals:
                if abs(z - (lam - n)) < 1e-8:
                    res_n.add(n)
    rmax = mp.mpf(0)
    for n in range(1, n_ord + 1):
        rhs = mp.matrix(N, 1)
        for m in range(0, n):
            pp = n + 1 - m
            if pp in Dp:
                rhs -= Dp[pp] * Nn[m]
        Mt = D1 + (n - lam) * mp.eye(N)
        if n in res_n:
            Nn[n], inc = _minnorm_solve(Mt, rhs, prec)
        else:
            try:
                Nn[n] = mp.lu_solve(Mt, rhs)
                continue
            except ZeroDivisionError:      # unlisted resonance (eigvals=None)
                Nn[n], inc = _minnorm_solve(Mt, rhs, prec)
        rn = max(abs(rhs[i, 0]) for i in range(N))
        rmax = max(rmax, inc / rn if rn > 0 else mp.mpf(0))
    diag['refine_rfinal'] = float(rmax)
    diag['log_branch_required'] = rmax > mp.mpf(10) ** (-(prec // 3))
    return Nn, diag


# ---------------------------------------------------------------- reduction
_RED_CACHE = {}     # (id(Dp), fingerprint, prec) -> reduction data + orbits


def _reduce_fuchs(Dp, N, prec, maxsteps=8):
    """Shear-reduce  t N = sum_p D_p u^{p-1} N  with a genuine nilpotent D_0
    block to a Fuchsian frame (u^2 N' = B(u) N, B_k = -D_k; each step scales
    an orthonormal basis of range(B_0) by u^-1).  Returns (Dt {p>=1}, Tj) with
    N(u) = T(u) Ntilde(u), T(u) = sum_j Tj[j] u^-j, or None on failure."""
    mp.mp.dps = prec
    kmax = max(Dp.keys())
    B = {k: (-1) * Dp[k] if k in Dp else mp.matrix(N, N)
         for k in range(0, kmax + 1)}
    Tj = [mp.eye(N)]
    for step in range(maxsteps):
        n0 = _matnorm(B[0])
        cur = max(_matnorm(B[k]) for k in range(1, min(3, kmax) + 1))
        if n0 <= cur * mp.mpf(10) ** -35:
            Dt = {p: (-1) * B[p] for p in B
                  if p >= 1 and _matnorm(B[p]) > 0}
            return Dt, Tj
        U0, S0, V0 = mp.svd(mp.matrix(B[0]))
        tol = S0[0] * mp.mpf(10) ** -30
        r = sum(1 for i in range(N) if S0[i] > tol)
        if r == 0 or 2 * r > N:
            return None
        Q = mp.matrix(N, N)
        for c in range(r):                                   # range(B0)
            for k in range(N):
                Q[k, c] = U0[k, c]
        for c in range(r):                                   # rowspace(B0)
            for k in range(N):
                Q[k, N - r + c] = mp.conj(V0[c, k])
        mid = []                                             # ker \ range
        for i in range(r, N):
            v = mp.matrix([mp.conj(V0[i, k]) for k in range(N)])
            for c in range(r):
                ip = sum(mp.conj(Q[k, c]) * v[k, 0] for k in range(N))
                for k in range(N):
                    v[k, 0] -= ip * Q[k, c]
            for w in mid:
                ip = sum(mp.conj(w[k, 0]) * v[k, 0] for k in range(N))
                for k in range(N):
                    v[k, 0] -= ip * w[k, 0]
            nv = mp.sqrt(sum(abs(v[k, 0]) ** 2 for k in range(N)))
            if nv > mp.mpf(10) ** -20:
                for k in range(N):
                    v[k, 0] /= nv
                mid.append(v)
        if len(mid) != N - 2 * r:
            return None
        for c, v in enumerate(mid):
            for k in range(N):
                Q[k, r + c] = v[k, 0]
        if _matnorm(Q.H * Q - mp.eye(N)) > mp.mpf(10) ** -25:
            return None
        Bp = {k: Q.H * (B[k] * Q) for k in B}
        Bn = {}
        for k in range(0, kmax + 2):
            Mtk = mp.matrix(N, N)
            src, up, dn = Bp.get(k), Bp.get(k - 1), Bp.get(k + 1)
            for i in range(N):
                for j in range(N):
                    if i < r and j >= r:                     # block (1,2): *u
                        if up is not None:
                            Mtk[i, j] = up[i, j]
                    elif i >= r and j < r:                   # block (2,1): /u
                        if dn is not None:
                            Mtk[i, j] = dn[i, j]
                    else:
                        if src is not None:
                            Mtk[i, j] = src[i, j]
            if k == 1:                                       # -u^2 S^-1 S'
                for i in range(r):
                    Mtk[i, i] += 1
            Bn[k] = Mtk
        B = Bn
        kmax += 1
        A = mp.matrix(N, N); Bm = mp.matrix(N, N)            # Q*S(u)=A+Bm/u
        for i in range(N):
            for j in range(N):
                if j < r:
                    Bm[i, j] = Q[i, j]
                else:
                    A[i, j] = Q[i, j]
        newT = [mp.matrix(N, N) for _ in range(len(Tj) + 1)]
        for j, Tm in enumerate(Tj):
            newT[j] += Tm * A
            newT[j + 1] += Tm * Bm
        Tj = newT
    return None


def _orbit_key(Dp, n_ord, prec):
    D1 = Dp[1]
    return (id(Dp), len(Dp), mp.nstr(D1[0, 0], 25), n_ord, prec)


def _reduced_orbits(Dp, N, lam, n_ord, prec):
    """Reduction + per-Z-orbit branch construction, cached per Dp object.
    Returns dict: leading-lambda(str,25) -> list of (Nn_reindexed, meta)."""
    from . import infinity
    key = _orbit_key(Dp, n_ord, prec)
    if key not in _RED_CACHE:
        red = _reduce_fuchs(Dp, N, prec)
        if red is None:
            _RED_CACHE[key] = None
            return None
        Dt, Tj = red
        evr = infinity.eig_spectrum(Dt[1], min(prec, 80))
        _RED_CACHE[key] = {'Dt': Dt, 'Tj': Tj, 'evr': evr, 'orbits': {}}
    C = _RED_CACHE[key]
    if C is None:
        return None
    Dt, Tj, evr = C['Dt'], C['Tj'], C['evr']
    jmax = len(Tj) - 1
    # orbit members: reduced eigenvalues congruent to lam mod 1
    mems = []
    for z in evr:
        dz = z - lam
        if abs(mp.im(dz)) < 1e-10 and \
           abs(mp.re(dz) - mp.nint(mp.re(dz))) < 1e-10 and \
           abs(mp.re(dz)) < jmax + n_ord + 2:
            if not any(abs(z - w) < 1e-15 for w in mems):
                mems.append(z)
    if not mems:
        return None
    obase = mp.nstr(min(mp.re(z) for z in mems), 25)
    if obase in C['orbits']:
        return C['orbits'][obase]
    # build all orbit branch series in the reduced (Fuchsian) frame
    n_full = n_ord + 2 * jmax + 4
    lam_ref = max(mp.re(z) for z in mems) + jmax
    raw = []
    rmax_orbit = mp.mpf(0)
    for lt in sorted(mems, key=lambda z: -mp.re(z)):
        Vv, ndim, sres = infinity.seed_vectors(Dt[1], lt)
        for c in range(ndim):
            seed = mp.matrix([Vv[r, c] for r in range(N)])
            dg = {}
            Nn, dg = _backward(Dt, N, lt, seed, n_full, prec, dg, evr)
            rmax_orbit = max(rmax_orbit, mp.mpf(dg['refine_rfinal']))
            off = int(mp.nint(mp.re(lam_ref - (lt + jmax))))
            Mn = {}
            for m in range(0, n_full + 1 + off):
                acc = mp.matrix(N, 1)
                for j in range(jmax + 1):
                    idx = (m - off) - (jmax - j)
                    if idx in Nn:
                        acc += Tj[j] * Nn[idx]
                Mn[m] = acc
            raw.append(Mn)
    # echelon by leading index -> well-defined leading exponents
    out = {}
    tolr = mp.mpf(10) ** (-(prec // 2) + 8)
    work = raw
    while work:
        lis = []
        for Mn in work:
            mx = max(_matnorm(v) for v in Mn.values())
            li = min((m for m in Mn if _matnorm(Mn[m]) > mx * tolr),
                     default=None)
            lis.append(li)
        work = [w for w, li in zip(work, lis) if li is not None]
        lis = [li for li in lis if li is not None]
        if not work:
            break
        limin = min(lis)
        piv_i = max((i for i in range(len(work)) if lis[i] == limin),
                    key=lambda i: _matnorm(work[i][limin]))
        piv = work.pop(piv_i)
        lis.pop(piv_i)
        le = lam_ref - limin
        # normalize: unit max-norm leading coefficient
        sc = _matnorm(piv[limin])
        piv = {m: piv[m] / sc for m in piv}
        out.setdefault(mp.nstr(le, 25), []).append(
            {'Mn': piv, 'li': limin, 'lam_ref': lam_ref})
        pv = piv[limin]
        pn2 = sum(abs(pv[i, 0]) ** 2 for i in range(N))
        for i, w in enumerate(work):
            if lis[i] == limin:
                ip = sum(mp.conj(pv[k, 0]) * w[limin][k, 0]
                         for k in range(N)) / pn2
                for m in w:
                    if m in piv:
                        w[m] = w[m] - ip * piv[m]
    res = {'branches': out, 'rmax': rmax_orbit, 'jmax': jmax, 'evr': evr,
           'multiset': {k: len(v) for k, v in out.items()}}
    C['orbits'][obase] = res
    return res


def genuine_fwd_blocks(Dp):
    """Forward blocks p<1 above the fit-noise floor (rel 1e-45 vs |D_1|,|D_2|;
    noise blocks sit at ~1e-120 on the reference family).  Shared by branch_series and
    the cli spectrum stage (reduced-residue classification)."""
    scale = max(_matnorm(Dp[p]) for p in (1, 2) if p in Dp)
    return sorted(pp for pp in Dp if pp < 1
                  and _matnorm(Dp[pp]) > scale * mp.mpf(10) ** -45)


def _seed_col(Dp, N, lam, seed):
    """LEGACY fallback (non-injective on degenerate families — see module
    docstring; prefer seed_index):
    argmax_j |<seed, V_j>| over the nullspace basis of (D1 - lam I)
    (same construction as infinity.seed_vectors)."""
    from . import infinity
    Vv, ndim, _ = infinity.seed_vectors(Dp[1], lam)
    best, bc = -1, 0
    for c in range(ndim):
        ip = abs(sum(mp.conj(Vv[k, c]) * seed[k, 0] for k in range(N)))
        if ip > best:
            best, bc = ip, c
    return bc


def branch_series(Dp, N, lam, seed, n_ord, prec, eigvals=None, seed_index=None):
    """Dp: {p: mp.matrix} Laurent blocks (build_Dp); lam = eig(D1) of the
    branch; seed: mp.matrix Nx1 with (D1 - lam I) seed = 0.
    Returns (Nn dict {0..n_ord: mp.matrix Nx1}, diag).
    Forward-block systems (genuine D_0): the seed only labels the request —
    the true seeds live in the reduced frame.  seed_index (the rank of the
    requested seed in ITS caller's nullspace basis) selects the orbit column
    injectively; when None, the legacy _seed_col argmax
    fallback is used (non-injective on degenerate families — avoid)."""
    mp.mp.dps = prec
    diag = {'resonances': resonance_orders(lam, eigvals or [], n_ord)
            if eigvals is not None else None}
    # genuine vs fit-noise forward blocks (noise ~1e-120 on the reference family)
    fwd = genuine_fwd_blocks(Dp)
    if not fwd:
        # pure backward recursion: (D1 + (n-lam) I) N_n = -sum_{m<n} D_{n+1-m} N_m
        return _backward(Dp, N, lam, seed, n_ord, prec, diag, eigvals)
    if fwd == [0]:
        orb = _reduced_orbits(Dp, N, lam, n_ord, prec)
        if orb is not None:
            cands = None
            for k, v in orb['branches'].items():
                if abs(mp.mpf(k) - mp.re(lam)) < 1e-10:
                    cands = v
                    break
            col = seed_index if seed_index is not None \
                else _seed_col(Dp, N, lam, seed)
            if cands is not None and col < len(cands):
                b = cands[col]
                Nn = {n: mp.matrix(b['Mn'][b['li'] + n])
                      for n in range(0, n_ord + 1)}
                diag['fwd_reduced'] = True
                diag['fwd_blocks'] = fwd
                diag['seed_col_used'] = col
                diag['orbit_multiset'] = orb['multiset']
                diag['refine_rfinal'] = float(orb['rmax'])
                diag['log_branch_required'] = \
                    orb['rmax'] > mp.mpf(10) ** (-(prec // 3))
                return Nn, diag
            diag['orbit_shortfall'] = {
                'requested': mp.nstr(lam, 20), 'col': col,
                'available': orb['multiset']}
        diag['reduce_failed'] = True
    # forward-coupled dense solve (irregular blocks D_{p<1}): joint system for
    # N_1..N_{n_ord}; equations n = 1-fdepth .. n_ord.  Eq(n):
    #   (n - lam) N_n [if 1<=n] + sum_m D_{n+1-m} N_m = 0,  N_0 = seed.
    fdepth = 1 - min(fwd)
    nU = n_ord
    Nn = {0: mp.matrix(seed)}
    rows, rhsl = [], []
    for n in range(1 - fdepth, n_ord + 1):
        r0 = mp.matrix(N, 1)
        if (n + 1) in Dp:
            r0 = -Dp[n + 1] * Nn[0]
        if n == 0:
            r0 = r0 + lam * Nn[0]      # full RHS = -(D1 - lam) N_0 (fix (2))
        blk = {}
        for m in range(1, nU + 1):
            pp = n + 1 - m
            if pp in Dp:
                blk[m] = Dp[pp].copy()
        if 1 <= n <= nU:
            if n not in blk:
                blk[n] = mp.matrix(N, N)
            else:
                blk[n] = blk[n].copy()
            for i in range(N):
                blk[n][i, i] += (n - lam)
        for k in range(N):
            nz = any(any(abs(blk[m][k, l]) > mp.mpf(10) ** (-(prec - 5))
                         for l in range(N)) for m in blk)
            if nz:
                rows.append((n, k, blk)); rhsl.append(r0[k, 0])
    nR = len(rows)
    BIGnp = np.zeros((nR, nU * N), dtype=np.complex128)
    BIGmp = mp.matrix(nR, nU * N)
    for ri, (n, k, blk) in enumerate(rows):
        for m, C in blk.items():
            for l in range(N):
                v = C[k, l]
                if v != 0:
                    BIGnp[ri, (m - 1) * N + l] = complex(v)
                    BIGmp[ri, (m - 1) * N + l] = v
    RHSmp = mp.matrix(nR, 1)
    for ri in range(nR):
        RHSmp[ri, 0] = rhsl[ri]
    rsc = np.ones(nR); csc = np.ones(nU * N); Beq = BIGnp.copy()
    for _ in range(3):
        rs = np.sqrt(np.maximum(np.max(np.abs(Beq), axis=1), 1e-300))
        Beq /= rs[:, None]; rsc *= rs
        cm = np.max(np.abs(Beq), axis=0); cm[cm < 1e-290] = 1.0
        cs = np.sqrt(cm); Beq /= cs[None, :]; csc *= cs
    sol = mp.matrix(nU * N, 1); rhs = RHSmp.copy(); rhist = []
    for it in range(80):
        rnp = np.array([complex(rhs[i, 0]) for i in range(nR)])
        rn = float(np.max(np.abs(rnp))); rhist.append(rn)
        if it > 0 and (rn < 10 ** (-(prec - 15)) or (it > 2 and rn > 0.95 * rhist[-2])):
            break
        dnp, _, _, _ = np.linalg.lstsq(Beq, rnp / rsc, rcond=None)
        dnp /= csc
        for i in range(nU * N):
            if not (np.isnan(dnp[i].real) or np.isnan(dnp[i].imag)):
                sol[i, 0] += mp.mpc(float(dnp[i].real), float(dnp[i].imag))
        rhs = RHSmp - BIGmp * sol
    for ui in range(nU):
        v = mp.matrix(N, 1)
        for i in range(N):
            v[i, 0] = sol[ui * N + i, 0]
        Nn[ui + 1] = v
    diag['refine_iters'] = len(rhist) - 1
    diag['refine_rfinal'] = rhist[-1]
    diag['log_branch_required'] = rhist[-1] > 10 ** (-(prec // 3))
    return Nn, diag


def phi_matrix(Mstart, branch_list, ai, dd, n_loops=3):
    """Phi(Mstart) [N x J]:  Phi[k,j] = Mstart^(L*d/2 + ai_k + lambda_j)
    * sum_n N_n^(j)[k] * Mstart^(-n)   — the branch-solution matrix whose
    column span the physical cut solution lives in (c-vector solve:
    R_rows.Phi.c = O_rows).  branch_list: [(lambda_j, Nn_dict), ...]."""
    N = len(ai)
    td2 = mp.mpf(n_loops) * dd / 2
    J = len(branch_list)
    Phi = mp.matrix(N, J)
    Ms = mp.mpf(Mstart) if not isinstance(Mstart, (mp.mpf, mp.mpc)) else Mstart
    x = 1 / Ms
    for j, (lam, Nn) in enumerate(branch_list):
        n_ord = max(Nn.keys())
        ser = mp.matrix(N, 1)
        for n in range(n_ord, -1, -1):           # Horner in x
            ser = ser * x + Nn[n]
        pref = mp.power(Ms, lam)
        for k in range(N):
            Phi[k, j] = mp.power(Ms, td2 + ai[k]) * pref * ser[k, 0]
    return Phi

#!/usr/bin/env python3
"""infinity.py — exact Laurent of x*A at x=inf, the D_p blocks and D1.

Family-agnostic: DE variable x (any name), L loops (alpha_k =
L*d/2 + ai_k).  Convention (README.md, section 'Convention'):
  N_k = x^{-alpha_k} M_k;  C-Laurent blocks D_p in u=1/x;  D1 = R_Delta -
  diag(alpha);  branch solutions  M_k ~ x^{alpha_k + lambda_j},
  lambda_j = eig(D1)   (N ~ x^{lambda_j} = u^{-lambda_j}).
"""
import mpmath as mp
mp.mp.dps = max(mp.mp.dps, 50)      # set FIRST


def laurent_xA(Dc, Gc, N, n_ord, prec):
    """R_j = coef of (1/x)^j in x*A = x*G/D, j = -qdeg .. n_ord+6 (exact:
    polynomial long division of x*G by D + 1/x-series of 1/D)."""
    mp.mp.dps = prec
    p = len(Dc) - 1; q = len(Gc) - 1
    qdeg = max(0, q + 1 - p)
    Rj = {j: mp.matrix(N, N) for j in range(-qdeg, n_ord + 7)}
    # series of 1/D at inf: 1/D = (1/x)^p * sum_k S_k (1/x)^k
    Dp = Dc[p]
    nS = n_ord + 8 + qdeg
    S = [mp.mpc(0)] * nS
    S[0] = 1 / Dp
    for k in range(1, nS):
        s = mp.mpc(0)
        for a in range(1, min(k, p) + 1):
            s += Dc[p - a] * S[k - a]
        S[k] = -s / Dp
    # x*G has coefs H_b = G_{b-1} for b=1..q+1, H_0=0; R_j = sum_b H_b S_{j+b-p}
    Hb = [mp.matrix(N, N)] + [Gc[b] for b in range(q + 1)]
    for j in range(-qdeg, n_ord + 7):
        R = mp.matrix(N, N)
        for b in range(q + 2):
            kk = j + b - p
            if 0 <= kk < nS:
                R += S[kk] * Hb[b]
        Rj[j] = R
    return Rj


def build_Dp(Rj, ai, dd, Delta, N, n_ord, n_loops=3):
    """Laurent blocks D_p of the connection C for the rescaled N-system:
    D_p[k,l] = R_{p-1+Delta_kl}[k,l] - delta_kl*alpha_k*[p==1].
    Returns dict {p: mp.matrix}; p < 1 are the forward (irregular) blocks."""
    td2 = mp.mpf(n_loops) * dd / 2
    alpha = [td2 + aik for aik in ai]
    pmax = max(max(row) for row in Delta)
    jmin = min(Rj.keys())
    Dp = {}
    for pp in range(jmin + 1 - pmax, n_ord + 3):
        M = mp.matrix(N, N); nz = False
        for k in range(N):
            for l in range(N):
                j = pp - 1 + Delta[k][l]
                if j in Rj and Rj[j][k, l] != 0:
                    M[k, l] = Rj[j][k, l]; nz = True
            if pp == 1:
                M[k, k] -= alpha[k]; nz = True
        if nz:
            Dp[pp] = M
    return Dp


def D1_matrix(Rj, ai, dd, Delta, N, n_loops=3):
    """The p=1 block alone: D1 = R_Delta - diag(L*d/2 + ai)."""
    td2 = mp.mpf(n_loops) * dd / 2
    D1 = mp.matrix(N, N)
    for k in range(N):
        for l in range(N):
            j = Delta[k][l]
            if j in Rj:
                D1[k, l] = Rj[j][k, l]
        D1[k, k] -= (td2 + ai[k])
    return D1


def eig_spectrum(D1, dps_eig=60):
    """Sorted eigenvalues of D1 (= the lambda_j) via mpmath at dps_eig."""
    save = mp.mp.dps
    mp.mp.dps = dps_eig
    E = mp.eig(D1, left=False, right=False)
    out = sorted([mp.mpc(e) for e in E], key=lambda z: (mp.re(z), mp.im(z)))
    mp.mp.dps = save
    return out


def seed_vectors(D1, lam, tol_exp=15):
    """Orthonormal nullspace of (D1 - lam*I) via mpmath SVD: the Frobenius
    seeds v_j of the branch lambda=lam.  Returns (V [N x m], max residual)."""
    N = D1.rows
    A = D1 - lam * mp.eye(N)
    U, S, V = mp.svd(mp.matrix(A))
    tol = S[0] * mp.mpf(10) ** (-tol_exp) if S[0] > 0 else mp.mpf(10) ** (-tol_exp)
    null_idx = [i for i in range(N) if S[i] < tol]
    Vm = mp.matrix(N, max(len(null_idx), 1))
    for c, i in enumerate(null_idx):
        for r in range(N):
            Vm[r, c] = mp.conj(V[i, r])
    R = A * Vm
    resid = max(abs(R[i, j]) for i in range(N) for j in range(Vm.cols)) \
        if null_idx else mp.mpf('inf')
    return Vm, len(null_idx), float(resid)

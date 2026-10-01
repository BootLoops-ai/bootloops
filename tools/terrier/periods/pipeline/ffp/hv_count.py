# HV/AESZ34 torus fiber counts over F_q for ALL psi in one pass.
# E_k(A,B) = #{x in (F_q^*)^k : sum x = A, sum 1/x = B}; E_k(A,B) = F_k(A*B) for A != 0,
# g_k = E_k(0, B!=0) (const), e_k = E_k(0,0).
# Step (convolve with one F_q^* variable):
#   F'(t) = g + [t==1](e-g) + sum_{x != 0,1} F[(1-x)*t + (1 - 1/x)]
#   g' = sum_u F[u] - F[1];  e' = (q-1)*F[1]
# N_torus(psi) = E_5(1, psi) = F_5(psi), fiber (sum X)(sum 1/X) = psi, psi = 1/phi.
import numpy as np
from field import Fq

def step(K, F, g, e, chunk=256):
    q = K.q
    one = K.one()
    Fn = np.full(q, g, dtype=np.int64)
    Fn[one] += e - g
    xs = np.array([x for x in range(1, q) if x != one], dtype=np.int64)
    # for each x: idx(t) = s*t + c, s = 1-x, c = 1 - 1/x
    for lo in range(0, len(xs), chunk):
        xc = xs[lo:lo + chunk]
        acc = np.zeros(q, dtype=np.int64)
        for x in xc:
            s = K.sub(one, x)
            c = K.sub(one, K.inv[x])
            idx = K.add_scalar(K.mul_scalar_all(s), c)
            acc += F[idx]
        Fn += acc
    gn = int(F.sum() - F[one])
    en = int((q - 1) * F[one])
    return Fn, gn, en

def base_mu1(K):
    F = np.zeros(K.q, dtype=np.int64)
    F[K.one()] = 1
    return F, 0, 0

def base_mu2(K2sub_p):
    """mu_2 over the *prime* field plane: one F_{p^2} conjugate-orbit variable.
    Returns (F, g, e) as functions on F_p (t = A*B).  K2sub_p: Fq(p,2)."""
    K = K2sub_p
    p = K.p
    y = np.arange(1, K.q, dtype=np.int64)
    a, b = y // p, y % p
    tr = (2 * a) % p                      # Tr(a+bw) = 2a
    iv = K.inv[y]
    tri = (2 * (iv // p)) % p
    # bucket by (Tr y, Tr 1/y)
    cnt = np.bincount(tr * p + tri, minlength=p * p).reshape(p, p)
    Frow = cnt[1, :].copy()               # Tr y = 1
    F = np.zeros(p, dtype=np.int64)
    # F(t) = mu2(1, t); need scaling: mu2(1,t) row is exactly Tr=1 row indexed by t
    F[:] = Frow
    g_vals = cnt[0, 1:]
    assert g_vals.min() == g_vals.max(), "mu2(0,B) not constant?!"
    return F, int(g_vals[0]), int(cnt[0, 0])

def torus_counts(p, deg, nvar=5, verbose=False):
    """F_5 array over F_q: N_torus(psi) at index psi (q=p^deg)."""
    K = Fq(p, deg)
    F, g, e = base_mu1(K)
    for k in range(nvar - 1):
        F, g, e = step(K, F, g, e)
        mass = (K.q - 1) * int(F.sum()) + (K.q - 1) * g + e
        assert mass == (K.q - 1) ** (k + 2), (k, mass)
    return K, F

def twisted12_counts(p, nvar_fixed=3):
    """(12)-twisted count over F_p: one conjugate-pair + nvar_fixed F_p^* variables."""
    Kp = Fq(p, 1)
    K2 = Fq(p, 2)
    F, g, e = base_mu2(K2)
    total = p * p - 1
    for k in range(nvar_fixed):
        F, g, e = step(Kp, F, g, e)
        total *= (p - 1)
        mass = (p - 1) * int(F.sum()) + (p - 1) * g + e
        assert mass == total, (k, mass, total)
    return F

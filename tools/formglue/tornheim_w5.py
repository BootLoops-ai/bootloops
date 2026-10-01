#!/usr/bin/env python3
r"""tornheim_w5.py — 5-leg signed lattice sum W₅ + 4-fold Mordell–Tornheim MT4
(formglue member).

This member imports tornheim by
identity for the shared primitives (_S2_asym/_asym_mul/_tail_sum/_zetaHp1_asym/
_default_Mmax_Pmax/_S2_pf_coeffs) and adds on top:
  MT4(a,b,c,d;e) = Σ_{n₁..n₄≥1} 1/(n₁^a n₂^b n₃^c n₄^d (Σn)^e)  — recursion to
    min-index 0; boundary via exact S₃ arrays + Euler–Maclaurin asymptotic tails
    (S₃ asym by the exact 1/M recursion down to PS₂, constant matched at Mmax).
  W5(p;q) = Σ_{m∈(ℤ\0)⁵:Σm=0} ∏|m_r|^{-p_r}(Σ|m_r|)^{-q}  — sign sectors
    (npos=0,5 impossible): (1,4) → 2·Σ_i MT4(p_{S\i};p_i+q)/2^q; (2,3) →
    2·Σ_{pairs} conv23 (S₂·S₃ same-argument convolution, product tail).

Validation (records not distributed with this package): published
D₅ = C_{1,1,1,1,1} Laurent (1502.06698 eq powerd5) 78–79d all 7 coefficients at
dps80 through laurent_full_w5; MT4 vs certified brute 15–24d (trunc-limited);
W5 vs box-lattice brute 9–18d; cross-dps 59d. All caches dps-keyed.
"""
import mpmath as mp
from mpmath import mpf
from itertools import combinations
try:
    from . import tornheim as _tornheim
    from .tornheim import (_S2_pf_coeffs, _S2_asym, _asym_mul, _tail_sum,
                           _zetaHp1_asym, _default_Mmax_Pmax)
except ImportError:  # script-mode / flat-path fallback (house shim pattern)
    import tornheim as _tornheim
    from tornheim import (_S2_pf_coeffs, _S2_asym, _asym_mul, _tail_sum,
                          _zetaHp1_asym, _default_Mmax_Pmax)

# ---- term algebra on A={(l,p):c} meaning Σ c·(ln M)^l/M^p ----
def _t_add(A, B):
    C = dict(A)
    for k, v in B.items(): C[k] = C.get(k, mpf(0)) + v
    return C

def _t_deriv(A):
    """d/dM of Σ c ln^l M/M^p."""
    C = {}
    for (l, p), c in A.items():
        if l: C[(l-1, p+1)] = C.get((l-1, p+1), mpf(0)) + c*l
        if p: C[(l, p+1)] = C.get((l, p+1), mpf(0)) - c*p
    return C

def _t_antideriv_term(l, p, c):
    """∫ c ln^l t·t^{-p} dt as term dict. Requires p≥1."""
    if p == 1:
        return {(l+1, 0): c/mpf(l+1)}
    out = {}
    # ∫ln^l t·t^{-p} = ln^l M·M^{1-p}/(1-p) − (l/(1-p))∫ln^{l-1}t·t^{-p}
    coef = c
    for ll in range(l, -1, -1):
        out[(ll, p-1)] = out.get((ll, p-1), mpf(0)) + coef/mpf(1-p)
        coef = -coef*ll/mpf(1-p)
        if ll == 0: break
    return out

def _t_antideriv(A):
    C = {}
    for (l, p), c in A.items():
        for k, v in _t_antideriv_term(l, p, c).items():
            C[k] = C.get(k, mpf(0)) + v
    return C

def _psum_Mpart(A, Pmax):
    """M-dependent part of Σ_{N=N0}^{M-1} f(N), f~A: F(M) − f(M)/2 +
    Σ_k B_{2k}/(2k)!·f^{(2k-1)}(M).  Constant of integration NOT included."""
    A = {k: v for k, v in A.items() if k[1] <= Pmax}
    out = _t_antideriv(A)
    out = _t_add(out, {k: -v/2 for k, v in A.items()})
    D = _t_deriv(A)  # f'
    k = 1
    while D and min(p for (_, p) in D) <= Pmax:
        fac = mp.bernoulli(2*k)/mp.factorial(2*k)
        out = _t_add(out, {kk: fac*v for kk, v in D.items() if kk[1] <= Pmax})
        D = _t_deriv(_t_deriv(D))  # f^{(2k+1)}
        k += 1
    return {k2: v for k2, v in out.items() if k2[1] <= Pmax and v != 0}

def _t_eval(A, M):
    lM = mp.log(M)
    return sum(c*lM**l/mpf(M)**p for (l, p), c in A.items())

# ---- exact S2 array: S2(c,d;K) for K=0..Mmax via running H^{(i)} ----
def _S2_exact_arr(c, d, Mmax):
    e = _S2_pf_coeffs(c, d); imax = len(e)-1; cd = c+d
    arr = [mpf(0)]*(Mmax+1)
    H = [mpf(0)]*(imax+1)
    for K in range(2, Mmax+1):
        invKm1 = mpf(1)/(K-1); pw = mpf(1)
        for i in range(1, imax+1):
            pw *= invKm1; H[i] += pw
        Kf = mpf(K); s = mpf(0)
        for i in range(1, imax+1):
            if e[i]: s += e[i]*H[i]/Kf**(cd-i)
        arr[K] = s
    return arr

# ---- PS2(c,d;M)=Σ_{N=2}^{M-1}S2(c,d;N): asym with matched constant ----
_PS2a_cache = {}
def _PS2_asym(c, d, Pmax, Mmax):
    c, d = sorted((c, d))
    key = (c, d, Pmax, Mmax, mp.mp.dps)
    if key in _PS2a_cache: return _PS2a_cache[key]
    D = _psum_Mpart(_S2_asym(c, d, Pmax), Pmax)
    arr = _S2_exact_arr(c, d, Mmax)
    ps = mpf(0)
    for N in range(2, Mmax):  # PS2 at M=Mmax
        ps += arr[N]
    Cst = ps - _t_eval(D, Mmax)
    A = _t_add(D, {(0, 0): Cst})
    _PS2a_cache[key] = A
    return A

# ---- S3(b,c,d;M)=Σ_{n1+n2+n3=M,ni≥1}∏ni^{-·}: asym via exact 1/M recursion ----
_S3a_cache = {}
def S3_asym(b, c, d, Pmax, Mmax):
    b, c, d = sorted((b, c, d))
    key = (b, c, d, Pmax, Mmax, mp.mp.dps)
    if key in _S3a_cache: return _S3a_cache[key]
    if b == 0:
        A = _PS2_asym(c, d, Pmax, Mmax)   # S3(0,c,d;M)=PS2(c,d;M)
    else:
        S = _t_add(_t_add(S3_asym(b-1, c, d, Pmax, Mmax),
                          S3_asym(b, c-1, d, Pmax, Mmax)),
                   S3_asym(b, c, d-1, Pmax, Mmax))
        A = {(l, p+1): v for (l, p), v in S.items() if p+1 <= Pmax}  # ·1/M
    _S3a_cache[key] = A
    return A

# ---- exact S3 array: S3(b,c,d;M), M=0..Mmax, via conv of j^{-b} with S2 ----
_S3arr_cache = {}
def _S3_exact_arr(b, c, d, Mmax):
    b, c, d = sorted((b, c, d))
    key = (b, c, d, Mmax, mp.mp.dps)
    if key in _S3arr_cache: return _S3arr_cache[key]
    s2 = _S2_exact_arr(c, d, Mmax)
    jb = [mpf(0)]*(Mmax+1)
    for j in range(1, Mmax+1): jb[j] = mpf(1)/mpf(j)**b
    arr = [mpf(0)]*(Mmax+1)
    for M in range(3, Mmax+1):
        s = mpf(0)
        for j in range(1, M-1):
            s += jb[j]*s2[M-j]
        arr[M] = s
    _S3arr_cache[key] = arr
    return arr

# ---- MT4 ----
_MT4_cache = {}
def MT4(a, b, c, d, e):
    """MT4(a,b,c,d;e)=Σ_{n₁..n₄≥1} 1/(n₁^a n₂^b n₃^c n₄^d (Σn)^e)."""
    a, b, c, d = sorted((a, b, c, d))
    key = (a, b, c, d, e, mp.mp.dps)
    if key in _MT4_cache: return _MT4_cache[key]
    if a == 0:
        v = _MT4_boundary(b, c, d, e)
    else:
        v = (MT4(a-1, b, c, d, e+1) + MT4(a, b-1, c, d, e+1)
             + MT4(a, b, c-1, d, e+1) + MT4(a, b, c, d-1, e+1))
    _MT4_cache[key] = v
    return v

_MT4b_cache = {}
def _MT4_boundary(b, c, d, e):
    """MT4(0,b,c,d;e)=Σ_{M≥3}S₃(b,c,d;M)·ζ_H(e,M+1),  b,c,d≥1, e≥2."""
    b, c, d = sorted((b, c, d))
    key = (b, c, d, e, mp.mp.dps)
    if key in _MT4b_cache: return _MT4b_cache[key]
    Mmax, Pmax = _default_Mmax_Pmax()
    s3 = _S3_exact_arr(b, c, d, Mmax)
    zhp1 = mp.zeta(e, 4)          # ζ_H(e,M+1) at M=3
    head = mpf(0)
    for M in range(3, Mmax+1):
        head += s3[M]*zhp1
        zhp1 -= mpf(1)/mpf(M+1)**e
    A = _asym_mul(S3_asym(b, c, d, Pmax, Mmax), _zetaHp1_asym(e, Pmax), Pmax)
    v = head + _tail_sum(A, Mmax)
    _MT4b_cache[key] = v
    return v

def _clear_mt4_caches():
    _PS2a_cache.clear(); _S3a_cache.clear(); _S3arr_cache.clear()
    _MT4_cache.clear(); _MT4b_cache.clear()
    _tornheim._clear_caches()

_conv23_cache = {}
def conv23(al, be, ga, de, ep, q):
    """Σ_{K≥3} S₂(α,β;K)·S₃(γ,δ,ε;K)/(2K)^q,  all indices ≥1, q≥1."""
    al, be = sorted((al, be)); ga, de, ep = sorted((ga, de, ep))
    key = (al, be, ga, de, ep, q, mp.mp.dps)
    if key in _conv23_cache: return _conv23_cache[key]
    Mmax, Pmax = _default_Mmax_Pmax()
    s2 = _S2_exact_arr(al, be, Mmax)
    s3 = _S3_exact_arr(ga, de, ep, Mmax)
    tq = mpf(2)**q
    head = mpf(0)
    for K in range(3, Mmax+1):
        head += s2[K]*s3[K]/(tq*mpf(K)**q)
    A = _asym_mul(_S2_asym(al, be, Pmax), S3_asym(ga, de, ep, Pmax, Mmax), Pmax)
    A = {(l, p+q): c/tq for (l, p), c in A.items() if p+q <= Pmax}
    v = head + _tail_sum(A, Mmax)
    _conv23_cache[key] = v
    return v

def W5(p, q):
    """W₅(p;q), p=(p₁..p₅), p_r≥1, q≥1. Sign-sector decomposition."""
    assert len(p) == 5
    idx = (0, 1, 2, 3, 4)
    total = mpf(0)
    # (1,4) sectors: lone leg i, ×2 for the m→−m mirror.
    for i in idx:
        rest = [p[j] for j in idx if j != i]
        total += 2*MT4(*rest, p[i]+q)/mpf(2)**q
    # (2,3) sectors: positive pair {i,j}, ×2 for the mirror.
    for pair in combinations(idx, 2):
        al, be = p[pair[0]], p[pair[1]]
        ga, de, ep = [p[j] for j in idx if j not in pair]
        total += 2*conv23(al, be, ga, de, ep, q)
    return total

def _clear_caches():
    _conv23_cache.clear(); _clear_mt4_caches()

if __name__ == '__main__':
    # regression smoke vs the recorded reference value of W5((1,1,1,1,1);1)
    mp.mp.dps = 30
    v = W5((1, 1, 1, 1, 1), 1)
    ref = mpf('32.84136187715471988163475202')
    err = abs(v - ref)
    print(f"W5((1,1,1,1,1);1) = {mp.nstr(v, 28)}  |err vs reference| = {mp.nstr(err, 3)}")
    raise SystemExit(0 if err < mpf('1e-25') else 1)

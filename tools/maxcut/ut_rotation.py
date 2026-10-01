#!/usr/bin/env python3
r"""
maxcut_ut_rotation.py — UT canonical rotation for a 1-dim quartic Baikov maxcut.

THEORY
======
A 2-loop N-prop maxcut whose Baikov polynomial reduces (after two nested
Beta integrals — the Sylvester cancellation) to a 1-dim integral

    m_k(ε;P) = Nhat(ε;P) · ∫_C R_k(u;ε) · P4(u;P)^{-1-ε} du,
    P4(u) = lc · ∏_{b=1}^{4}(u - r_b),

with R_k polynomial (deg ≤ 2) numerator insertions, has an EXPLICIT canonical
ε-form basis given by the partial-fraction (twisted-cohomology) generators

    χ_a(ε;C) := ∫_C (u-r_a)^{-1} · W(u)^{-ε} du,      W := ∏_b(u-r_b),
    Σ_a χ_a = 0   (closed cycle).

These satisfy the Aomoto / KZ connection in PURE ε-form

    ∂_v χ_a = -ε · Σ_{c≠a} ∂_v log(r_a-r_c) · (χ_a - χ_c)          (★)

i.e. ∂χ = ε·B·χ with B the dlog connection on root differences.  A 3-dim
canonical basis (using Σχ=0) is g_can = S·χ; when two of the four roots are
rational (here r=s15 and r=s15+2(s45-mm4), a 5-point one-mass kinematic
configuration), take S so that two rows are those LINEAR roots and the third
is the SURD-antisymmetric χ_{r3}-χ_{r4} (the √Gram5-odd combination).

The UT rotation g_can = T · m_raw is then  T = D^{-1} / (Nhat · lc^{-1-ε})
where the 3×3 matrix D is the partial-fraction matrix C reduced by S⁺:

    m_k = Nhat·lc^{-1-ε}·Σ_a c_{k,a}·χ_a,   c_{k,a}=R_k(r_a)/∏_{b≠a}(r_a-r_b),
    D_{k,1}=c_{k,L1}-½(c_{k,S1}+c_{k,S2}),  D_{k,2}=c_{k,L2}-½(...),
    D_{k,3}=½(c_{k,S1}-c_{k,S2}).

D_{k,1},D_{k,2} are rational; D_{k,3} ∝ 1/√Gram5 (one (r3-r4) in the
denominator of c).  Hence  D = T_pre^{-1}·diag(1,1,1/√Gram5)  with T_pre
rational, i.e.  T = diag(1,1,√Gram5)·T_pre / (Nhat·lc^{-1-ε})  — a
"T_pre·diag(1,√G,√G)" structure (up to slot permutation).

VALIDATION
==========
(★) is verified by gauge-transforming an independently derived (Kira/IBP)
3×3 maxcut DE M_v(d) by the explicit T(ε;P) at several ε and checking the
result is ε·(ε-independent).

PUBLIC API
==========
    classify_roots(roots, r_L1, r_L2)        → (iL1,iL2,iS1,iS2) sorted-idx map
    pf_matrix(roots, Rk_funcs, eps)          → c[k,a] 3×4
    D_matrix(c, idx_map)                     → 3×3 D
    T_of_eps(roots, Rk_funcs, Nhat, lc, eps) → 3×3 numeric T
    chi_integral(roots, a, cyc, eps)         → χ_a on cycle cyc (direct quad)
    kz_connection_4(roots, droots)           → 4×4 B (dlog on root diffs)
    kz_connection_3(roots, droots, idx_map)  → 3×3 reduced B
    apply_T_laurent(Tcoeffs, m_laurent)      → g_can Laurent (series product)
"""
from __future__ import annotations
import mpmath as mp


# ---------------------------------------------------------------------------
def classify_roots(roots, r_L1, r_L2, tol=None):
    """Given the 4 sorted roots and the two known rational roots r_L1=s15,
    r_L2=s15+2(s45-mm4), return (iL1,iL2,iS1,iS2) sorted-indices, with
    (iS1,iS2) the surd pair, iS1<iS2."""
    if tol is None:
        tol = mp.mpf(10)**(-mp.mp.dps//2)
    idx = list(range(4))
    iL1 = min(idx, key=lambda i: abs(roots[i]-r_L1)); idx.remove(iL1)
    iL2 = min(idx, key=lambda i: abs(roots[i]-r_L2)); idx.remove(iL2)
    assert abs(roots[iL1]-r_L1) < tol and abs(roots[iL2]-r_L2) < tol, \
        f"linear roots not matched: {roots} vs {r_L1},{r_L2}"
    iS1, iS2 = sorted(idx)
    return iL1, iL2, iS1, iS2


def pf_matrix(roots, Rk_funcs, eps=None):
    """c[k][a] = R_k(r_a) / ∏_{b≠a}(r_a-r_b),  k=0..2, a=0..3."""
    c = [[None]*4 for _ in range(len(Rk_funcs))]
    for a in range(4):
        ra = roots[a]
        denom = mp.mpf(1)
        for b in range(4):
            if b != a:
                denom *= (ra - roots[b])
        for k, Rk in enumerate(Rk_funcs):
            num = Rk(ra) if eps is None else Rk(ra, eps)
            c[k][a] = num/denom
    return c


def D_matrix(c, idx_map):
    """3×3 D from 3×4 c via the (L1,L2,S-) reduction (Σχ=0)."""
    iL1, iL2, iS1, iS2 = idx_map
    D = mp.matrix(3, 3)
    for k in range(3):
        sym = (c[k][iS1]+c[k][iS2])/2
        D[k, 0] = c[k][iL1] - sym
        D[k, 1] = c[k][iL2] - sym
        D[k, 2] = (c[k][iS1]-c[k][iS2])/2
    return D


def T_of_eps(roots, Rk_funcs, Nhat_val, lc, eps, idx_map):
    """Full UT rotation T(ε;P) = D^{-1} / (Nhat · lc^{-1-ε}).
    Returns 3×3 mpmath matrix."""
    c = pf_matrix(roots, Rk_funcs, eps)
    D = D_matrix(c, idx_map)
    pre = Nhat_val * mp.power(lc, -1-eps)
    return D**-1 / pre


# ---------------------------------------------------------------------------
def chi_integral(roots, a, cyc, eps, maxd=None):
    """χ_a(ε) = ∫_{r_i}^{r_j} (u-r_a)^{-1} · ∏_b(u-r_b)^{-ε} du  on the
    OPEN interval (r_i,r_j) with i=cyc[0],j=cyc[1] sorted-ADJACENT (no
    interior root).  Endpoint-subtracted analytic continuation as in
    maxcut_close.cycle_integral_reg, but exponent is -ε (not -1-ε) on
    the non-a factors and -1-ε on (u-r_a) iff a∈cyc."""
    i, j = cyc
    ri, rj = roots[i], roots[j]
    dr = rj - ri
    oth = [roots[b] for b in range(4) if b not in (i, j)]
    # split the a-factor:  if a is endpoint, it carries -1-ε; else it is in oth
    a_is_i = (a == i); a_is_j = (a == j)
    # prefactor pulls out (u-r_i)^{..}(u-r_j)^{..} → t,1-t powers
    pw_i = -1-eps if a_is_i else -eps
    pw_j = -1-eps if a_is_j else -eps
    pre = mp.power(dr, pw_i)*mp.power(-dr, pw_j)*dr
    ra = roots[a]

    def h(t):
        ut = ri + dr*t
        val = mp.mpc(1)
        for ro in oth:
            val *= mp.power(ut-ro, -eps)
        if not (a_is_i or a_is_j):
            val *= 1/(ut-ra)
        return val
    half = mp.mpf('0.5')
    md = maxd or max(6, int(mp.mp.dps*0.08)+4)

    def piece(pw, h0, side):
        """∫_0^{1/2} s^{pw} · H(s) ds with single subtraction at s=0.
        Always subtracts (valid for all pw≠-1): the previous Re(pw)>-1
        direct-quad shortcut lost ~2d/|pw+1| precision when pw was close
        to -1 (the ε<0 quartic-Pu case)."""
        H0 = h0(mp.mpf(0))
        A1 = H0*mp.power(half, pw+1)/(pw+1)
        A2 = mp.quad(lambda s: mp.power(s, pw)*(h0(s)-H0), [0, half],
                     method='tanh-sinh', maxdegree=md)
        return A1+A2
    # left half t∈[0,1/2]:  t^{pw_i}·(1-t)^{pw_j}·h(t)
    L = piece(pw_i, lambda t: mp.power(1-t, pw_j)*h(t), 'L')
    # right half s=1-t∈[0,1/2]: s^{pw_j}·(1-s)^{pw_i}·h(1-s)
    R = piece(pw_j, lambda s: mp.power(1-s, pw_i)*h(1-s), 'R')
    return pre*(L+R)


# ---------------------------------------------------------------------------
def kz_connection_4(roots, droots):
    """4×4 B with ∂_v χ = ε·B·χ.  droots[a] = ∂_v r_a."""
    B = mp.matrix(4, 4)
    for a in range(4):
        diag = mp.mpf(0)
        for c in range(4):
            if c == a:
                continue
            dl = (droots[a]-droots[c])/(roots[a]-roots[c])
            B[a, c] = dl
            diag -= dl
        B[a, a] = diag
    return B


def kz_connection_3(roots, droots, idx_map):
    """3×3 reduced B in the (χ_L1, χ_L2, χ_S-) basis using Σχ=0."""
    iL1, iL2, iS1, iS2 = idx_map
    B4 = kz_connection_4(roots, droots)
    # χ = S⁺·g  with χ_L1=g1, χ_L2=g2, χ_S1=(-g1-g2+g3)/2, χ_S2=(-g1-g2-g3)/2
    Sp = mp.matrix(4, 3)
    Sp[iL1, 0] = 1; Sp[iL2, 1] = 1
    Sp[iS1, 0] = Sp[iS1, 1] = mp.mpf(-1)/2; Sp[iS1, 2] = mp.mpf(1)/2
    Sp[iS2, 0] = Sp[iS2, 1] = mp.mpf(-1)/2; Sp[iS2, 2] = mp.mpf(-1)/2
    # g = S·χ
    S = mp.matrix(3, 4)
    S[0, iL1] = 1; S[1, iL2] = 1; S[2, iS1] = 1; S[2, iS2] = -1
    return S*B4*Sp


# ---------------------------------------------------------------------------
def laurent_mul(A, B, ordA0, ordB0, n_keep):
    """Multiply two Laurent series (lists of coeffs starting at ε^{ordA0},
    ε^{ordB0}) → list starting at ε^{ordA0+ordB0}, length n_keep.
    Coeffs may be scalars or mp.matrix (compatible shapes)."""
    out = [0]*n_keep
    for i, a in enumerate(A):
        for j, b in enumerate(B):
            k = i+j
            if k < n_keep:
                out[k] = out[k] + a*b if out[k] != 0 else a*b
    return out, ordA0+ordB0


def matrix_laurent_inverse(Dser, n_keep):
    """Given D(ε)=Σ D_n ε^n (D_0 invertible), return D^{-1} to n_keep orders."""
    D0i = Dser[0]**-1
    out = [D0i]
    for n in range(1, n_keep):
        acc = mp.matrix(D0i.rows, D0i.cols)
        for k in range(1, n+1):
            if k < len(Dser):
                acc += Dser[k]*out[n-k]
        out.append(-D0i*acc)
    return out


# ---------------------------------------------------------------------------
if __name__ == '__main__':
    # smoke test on a toy quartic
    mp.mp.dps = 40
    rts = [mp.mpf(x) for x in (-3, -1, 1, 4)]
    eps = mp.mpf('0.1')
    # Σχ_a = 0 ?
    s = sum(chi_integral(rts, a, (1, 2), eps) for a in range(4))
    print("Σχ_a (should be 0):", mp.nstr(s, 20))
    # KZ: FD vs formula
    h = mp.mpf('1e-12')
    rtsP = [rts[0]+h]+rts[1:]
    chiP = [chi_integral(rtsP, a, (1, 2), eps) for a in range(4)]
    chiM = [chi_integral([rts[0]-h]+rts[1:], a, (1, 2), eps) for a in range(4)]
    chi0 = [chi_integral(rts, a, (1, 2), eps) for a in range(4)]
    dchi = [(chiP[a]-chiM[a])/(2*h) for a in range(4)]
    B = kz_connection_4(rts, [mp.mpf(1), 0, 0, 0])
    Bchi = [eps*sum(B[a, c]*chi0[c] for c in range(4)) for a in range(4)]
    for a in range(4):
        print(f"a={a}: |∂χ - εBχ|/|∂χ| =",
              mp.nstr(abs(dchi[a]-Bchi[a])/max(abs(dchi[a]), mp.mpf('1e-30')), 6))

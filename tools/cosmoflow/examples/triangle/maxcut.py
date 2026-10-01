#!/usr/bin/env python3
r"""maxcut.py -- WORKED EXAMPLE (cosmoflow.examples.triangle): maxcut
residue-period sampler at q_{G12}=0 + symbolic L_2 proof for the FRW 1-loop
3-site (triangle) elliptic sub-sector (arXiv:2408.16386).

  * period_residue(a, lam, cq): fast (~0.1s @100d) evaluation of
        I(lam;c) = (1/c) INT_{cycle} B^{-1/2} dy23 dy31   at y12^2 = (S/c)^2
    via inner ellipk closed form + sin^2-theta substitution.  Verified
    against direct 2D nested quad to 25d.  NOTE: this residue cycle carries
    the FULL 9-master lower-block (1+2+2) inhomogeneity, so I(lam) has NO
    order<=3 poly-coeff annihilator (numeric-annihilator route STOPPED).
  * L2_paper_coeffs, K2_modulus, varpi0, varpi1: the paper's 2nd-order PF
    factor L_2 and its two solutions
        varpi_{0,1}(lam;a) = K^{(')}(K^2) / sqrt(1 - (a+1)^2 lam^2),
        K^2 = [(a-1)^2 lam^2 - 1] / [(a+1)^2 lam^2 - 1].
  * prove_L2_varpi0_symbolic(): sympy chain-rule reduction of L_2[varpi0]
    to a multiple of the Legendre ODE; vanishes identically for generic a.

n_s = 3 SPECIFIC -- this is the worked example, not the general tool.
Import as `from cosmoflow.examples.triangle import maxcut` (unrelated to
the top-level tools/maxcut package).
"""
import sympy as sp
import mpmath as mp
from mpmath import mpf, sqrt, ellipk

from ...polytope import B_z, z1, z2, z3, X1, X2, X3


# ---------------------------------------------------------------------------
# Maxcut residue-period sampler
# ---------------------------------------------------------------------------
_P = sp.Poly(B_z, z2, z3)
_cf = {(i, j): sp.expand(_P.coeff_monomial(z2**i * z3**j))
       for i in range(3) for j in range(3)}
assert _cf[(1, 2)] == 0 and _cf[(2, 2)] == 0 and _cf[(2, 1)] == 0, "B not biquadratic"
_lam = {k: sp.lambdify((z1, X1, X2, X3), v, 'mpmath') for k, v in _cf.items()}


def period_residue(a, lam, cq, dps=80):
    """Residue period at q_{G12}=0 (see module docstring). ~0.1s @100d."""
    mp.mp.dps = dps + 20
    a = mpf(a); cq = mpf(cq)
    Xv = (a * lam, lam, mpf(1)); S = sum(Xv)
    z1v = (S / cq) ** 2
    A2  = _lam[(0, 2)](z1v, *Xv)
    a11 = _lam[(1, 1)](z1v, *Xv); a10 = _lam[(0, 1)](z1v, *Xv)
    a02 = _lam[(2, 0)](z1v, *Xv); a01 = _lam[(1, 0)](z1v, *Xv); a00 = _lam[(0, 0)](z1v, *Xv)
    D2 = a11**2 - 4*A2*a02; D1 = 2*a11*a10 - 4*A2*a01; D0 = a10**2 - 4*A2*a00
    dd = sqrt(D1**2 - 4*D2*D0)
    r1 = (-D1 - dd)/(2*D2); r2 = (-D1 + dd)/(2*D2)
    z2lo, z2hi = (r1, r2) if r1 < r2 else (r2, r1)
    assert D2 > 0 and z2lo > 0
    nA2 = -A2; sqlo = sqrt(z2lo)
    def Jth(th):
        s = mp.sin(th); c = mp.cos(th); z2v = z2lo*s**2
        A1 = a11*z2v + a10
        r = sqrt(D2 * z2lo * c**2 * (z2hi - z2v))
        z3a = (-A1 - r)/(2*A2); z3b = (-A1 + r)/(2*A2)
        z3m, z3p = (z3a, z3b) if z3a < z3b else (z3b, z3a)
        return (2/sqrt(z3p)) * ellipk(1 - z3m/z3p) / sqrt(nA2) * sqlo * c / 2
    I = mp.quad(Jth, [0, mp.pi/2]) / cq
    mp.mp.dps = dps
    return +I


# ---------------------------------------------------------------------------
# Paper's L_2 and its analytic periods
# ---------------------------------------------------------------------------
def L2_paper_coeffs(a, lam):
    D = (a**2-1)**2*lam**4 - 2*(a**2+1)*lam**2 + 1
    return ((5*(a**2-1)**2*lam**4 - 6*(a**2+1)*lam**2 + 1)/(lam*D),
            (3*(a**2-1)**2*lam**2 - 2*(a**2+1))/D)


def K2_modulus(a, lam):
    return ((a-1)**2*lam**2 - 1) / ((a+1)**2*lam**2 - 1)


def varpi0(a, lam):
    return ellipk(K2_modulus(a, lam)) / sqrt(1 - (a+1)**2 * lam**2)


def varpi1(a, lam):
    return ellipk(1 - K2_modulus(a, lam)) / sqrt(1 - (a+1)**2 * lam**2)


def prove_L2_varpi0_symbolic():
    """Return True iff L_2[K(K^2)/sqrt(1-(a+1)^2 lam^2)] == 0 identically
    (sympy, generic a), by eliminating K'' via the Legendre ODE."""
    lamS, aS, u, v = sp.symbols('lam a u v')   # u=K(m), v=K'(m)
    D4 = (aS**2-1)**2*lamS**4 - 2*(aS**2+1)*lamS**2 + 1
    p1 = (5*(aS**2-1)**2*lamS**4 - 6*(aS**2+1)*lamS**2 + 1)/(lamS*D4)
    q0 = (3*(aS**2-1)**2*lamS**2 - 2*(aS**2+1))/D4
    m  = ((aS-1)**2*lamS**2 - 1)/((aS+1)**2*lamS**2 - 1)
    f  = 1/sp.sqrt(1 - (aS+1)**2*lamS**2)
    m1 = sp.diff(m, lamS); m2 = sp.diff(m, lamS, 2)
    Kpp = (u/4 - (1-2*m)*v)/(m*(1-m))          # Legendre ODE
    psi0 = f*u
    psi1 = sp.diff(f, lamS)*u + f*v*m1
    psi2 = sp.diff(f, lamS, 2)*u + 2*sp.diff(f, lamS)*v*m1 + f*(Kpp*m1**2 + v*m2)
    return sp.simplify(sp.together(psi2 + p1*psi1 + q0*psi0)) == 0

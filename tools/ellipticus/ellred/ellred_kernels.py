"""ellred consolidated module: standard-elliptic (Legendre F/E/Pi) reduction of
three named elliptic kernel classes K_asin, K_W, K_alg.

CURVES (demonstration family, ell_normal.P4_exact):
  E      : y^2 = P4(x;tau,z23),  x=u^2   -- K_asin (784-core groups) + K_W.
           P4 = 16 tau^2 (z23-1)^4 x^4 + ... ; lead/const perfect squares.
           Legendre chart (ellred_engine.chart): 2 complex-conj root pairs;
           Mobius zeta=(x-t1)/(x-t2), phi=atan(sqrt(a1/b1) zeta), m=1-(a2 b1)/(a1 b2).
  E2     : the w2s-qb2-dL asin curve (j=-242.2): the 103849-core, cubic image
           y^2=-x(x^2+25950 x/103849+5625/103849) = 1 real + 2 complex roots.
  E_alg  : y^2=(1-x)(k(z) x^2 - k(z) x + c0(z)), x=u^2 -- K_alg. 3 real roots;
           j non-constant in z (non-isotrivial).

STANDARD ATOMS (all gated >=30d vs direct quad; ellred_engine / ellred_cubic):
  Quartic E (K_asin/K_W measure): F-atom 1/sqrt(P4), E-atom (a2 z^2+b2)/sqrt..,
    Pi-atom 1/((z^2-c) sqrt..). int x^k dx/sqrt(P4)=J_k (recurrence). All 40-50d.
  Cubic (K_alg E_alg; K_asin even-quartic image x*(784x^2-784x+125)):
    Jmoments_cubic J_k -> F/E (3-real-root or 1-real-2-complex region maps). 30-38d.

KERNEL REDUCTIONS:
  K_asin: by-parts in u; assembled EVEN cofactor is an EXACT POLYNOMIAL over
    sqrt(R), R=784u^4-784u^2+125 (cofactor probe resid 0.0 for the 784-core).
    => the genuine elliptic part of K_asin is PURE 1st/2nd-kind: F + E (no Pi,
    no weight-2 residue) on the cubic image y^2 = -x(784x^2-784x+125) [3 real
    roots, region (e2,e3)].  I_ev = (1/2) sum_k p_k J_k.  (One lone group on E2
    103849-core, 1r2c.)
  K_W: WEIGHT-2 elliptic dilogarithm (16 kL1 Li2/loglog groups). NOT F/E/Pi.
    By-parts lowers Li2/loglog to weight-1 letter x measure dx/sqrt(P4): an
    elliptic dilogarithm E_2 / BDGS eMPL Gamma-tilde of length 2.  Named kernel.
  K_alg: PURE measure (asin/dilog-free) = int (cofactor over Q(sqrt5,sqrt x,
    sqrt(1-x))) sqrt(E_alg cubic) dx.  Genuine elliptic part: 3rd-kind F/E/Pi on
    E_alg; coupled to genus-0 sqrt(1-x) letters.  No weight-2 part.

This module re-exports the engines; battery leg L5 gates the atoms against
direct quadrature.
"""
import os, sys, json
from . import ellred_engine as ENG      # quartic F/E/Pi atoms + Jmoments + reduce_pole
from . import ellred_cubic as CUB       # cubic Jmoments (3r + 1r2c)
from . import ell_normal as EN

# convenience re-exports
chart = ENG.chart
atomF = ENG.atomF
atomE_closed = ENG.atomE_closed
atomPi_closed = ENG.atomPi_closed
Jmoments_cubic = CUB.Jmoments_cubic
I0_cubic = CUB.I0_cubic

def P4_at(tau, z23):
    """exact P4 coeffs (ascending mpf) at rational (tau,z23) for curve E."""
    import sympy as sp
    from mpmath import mpf
    t, z, x = sp.symbols('tau z23 x', positive=True)
    core, _, _ = EN.P4_exact()
    P4c = sp.Poly(core.subs({t: sp.nsimplify(tau), z: sp.nsimplify(z23)}), x)
    return [mpf(str(sp.nsimplify(c))) for c in reversed(P4c.all_coeffs())]

if __name__ == '__main__':
    # smoke test
    from mpmath import mp, mpf, sqrt as msqrt, quad as mquad, nstr
    mp.dps = 40
    P4 = P4_at('3/8', '3/10')
    ch = chart(P4)
    print("E modulus m(A) =", nstr(ch['m'], 30))
    print("OK: engines importable; battery leg L5 carries the per-atom gates.")

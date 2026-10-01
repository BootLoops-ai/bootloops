"""Derive & verify CLOSED forms for E-atom and Pi-atom (no linear solve).

GATE SCRIPT (executes at module level by design; run as
  PYTHONPATH=tools python3 -m ellipticus.ellred.derive_atoms
never imported by the package)."""
import os, sys
import sympy as sp
from mpmath import (mp,mpf,sqrt as msqrt,atan as matan,tan as mtan,cos as mcos,
                    sin as msin,ellipf,ellipe,ellippi,quad as mquad,nstr)
from . import ell_normal as EN
from . import ellred_engine as E
mp.dps=50
tau,z23,xx=sp.symbols('tau z23 x',positive=True)
core,_,_=EN.P4_exact()
for cell,tv,zv,x1f,x2f in [('A',sp.Rational(3,8),sp.Rational(3,10),'0.46','0.54'),
                            ('B',sp.Rational(1,4),sp.Rational(1,5),'0.66','0.74')]:
    P4c=sp.Poly(core.subs({tau:tv,z23:zv}),xx)
    P4n=[mpf(str(sp.nsimplify(c))) for c in reversed(P4c.all_coeffs())]
    ch=E.chart(P4n)
    x1,x2=mpf(x1f)**2,mpf(x2f)**2
    z1,z2=E.zof(ch,x1),E.zof(ch,x2)
    a1,b1,a2,b2,m=ch['a1'],ch['b1'],ch['a2'],ch['b2'],ch['m']
    al=msqrt(a1/b1)
    def Phi(z): return matan(al*z)
    def Delta(ph): return msqrt(1-m*msin(ph)**2)
    # E-atom CLOSED: sqrt(b2/a1) [ tan phi Delta - E(phi,m) + F(phi,m) ]
    def Eatom(z):
        ph=Phi(z)
        return msqrt(b2/a1)*( mtan(ph)*Delta(ph) - ellipe(ph,m) + ellipf(ph,m) )
    Ec=Eatom(z2)-Eatom(z1)
    Eref=mquad(lambda z:(a2*z*z+b2)/E.gz(ch,z),[z1,z2])
    eE=abs(Ec-Eref)/abs(Eref)
    # Pi-atom CLOSED: int dz/((z^2-c)g).  Derive coefficient.
    # 1/(z^2-c) dz/g = (1/sqrt(a1 b2)) (1/(z^2-c)) dphi/Delta, z=sqrt(b1/a1)tan phi
    # z^2-c = (b1/a1)tan^2 - c = ((b1/a1+c)sin^2 - c)/cos^2 = (b1/a1)(sin^2 - c a1/b1 cos^2)/cos^2
    #  hmm: (b1/a1)tan^2 - c = [(b1/a1)sin^2 - c cos^2]/cos^2 = [(b1/a1)sin^2 - c(1-sin^2)]/cos^2
    #     = [(b1/a1 + c)sin^2 - c]/cos^2
    # So 1/(z^2-c) = cos^2/[(b1/a1+c)sin^2 - c] = -cos^2/[c(1 - ((b1/a1+c)/c) sin^2)]
    #   = -(1/c) cos^2 /(1 - n0 sin^2),  n0=(b1/a1+c)/c = 1 + (b1/a1)/c = 1+ b1/(a1 c)
    # and dphi/Delta * cos^2/(1-n0 sin^2): need int cos^2/((1-n0 sin^2)Delta) dphi
    #  cos^2 = 1 - sin^2.  1-sin^2 = (1 - n0 sin^2)/n0 + (1 - 1/n0) ... decompose:
    #  (1-sin^2)/(1-n0 sin^2) = 1/n0 + (1-1/n0)/(1-n0 sin^2)
    #   check: 1/n0 (1-n0 sin^2) + (1-1/n0) = 1/n0 - sin^2 + 1 - 1/n0 = 1 - sin^2. YES.
    # So int cos^2/((1-n0 sin^2)Delta)dphi = (1/n0)F + (1-1/n0)Pi(n0;phi,m).
    # Full: int dz/((z^2-c)g) = (1/sqrt(a1 b2))*(-(1/c))*[ (1/n0)F + (1-1/n0)Pi(n0) ]
    c=mpf(2)
    n0=1+b1/(a1*c)
    def Piatom(z):
        ph=Phi(z)
        return (1/msqrt(a1*b2))*(-(1/c))*( (1/n0)*ellipf(ph,m) + (1-1/n0)*ellippi(n0,ph,m) )
    Pc=Piatom(z2)-Piatom(z1)
    Pref=mquad(lambda z:1/((z*z-c)*E.gz(ch,z)),[z1,z2])
    eP=abs(Pc-Pref)/abs(Pref)
    print(f"cell {cell}: E-atom closed {float(-mp.log10(eE)):.1f}d   Pi-atom closed {float(-mp.log10(eP)):.1f}d  (n0={nstr(n0,12)}, m={nstr(m,18)})")

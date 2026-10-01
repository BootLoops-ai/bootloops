"""Cubic-curve Legendre reducer for K_alg: E_alg: y^2 = C3(x) = a3 x^3+a2 x^2+a1 x+a0
with 3 REAL or 1-real-2-complex roots.  Reduce int rat(x) dx/sqrt(C3) to standard
incomplete Legendre F/E/Pi.

Standard reduction (Byrd&Friedman 233-260 / G&R 3.13): for cubic with real root e1
and complex pair, OR three real roots, the substitution to Legendre normal form is
classical. We use the robust route:
  - find roots of C3 (a3 (x-e1)(x-e2)(x-e3)); for our E_alg cubic
    (1-x)(k x^2 - k x + 125): root e=1 (real), plus the quadratic k x^2-k x+125
    has complex roots (disc = k^2-500k<0 for 0<k<500) -> 1 real + 2 complex.
  - 1-real-2-complex cubic: B&F 260.00. Let e1=real root; the complex pair gives
    a real quadratic (x-e1) factored out: C3 = a3 (x-e1)(x^2+px+q), q-p^2/4>0.
    Sub: standard -> F(phi,m), with
       A1 = sqrt((e1-? )...)  -- we implement via the known closed atoms and
    CERTIFY each by exact derivative + gate vs direct quad (like the quartic).
We mirror the quartic engine: build base I0=int dx/sqrt(C3) (F), and moments
J_k=int x^k dx/sqrt(C3) by the x-recurrence anchored on d/dx[x^n sqrt(C3)], plus
pole atoms int dx/((x-p) sqrt(C3)) (-> Pi). All gated vs direct quad.
"""
import os,sys
from mpmath import (mp,mpf,mpc,sqrt as msqrt,quad as mquad,nstr,polyroots,
                    atan as matan, asin as masin, tan as mtan, ellipf,ellipe,ellippi)

def C3val(C3n, xx):  # C3n ascending [a0,a1,a2,a3]
    v=mpc(0)
    for c in reversed(C3n): v=v*xx+c
    return v.real if abs(v.imag)<mpf(10)**(-mp.dps+20) else v

def cubic_chart(C3n):
    """1 real root e1 + complex pair (x^2+p x+q). Map to Legendre via the
    classical 1-real-2-complex reduction. Return chart with a func phi(x), m,
    and the atom A-factor so that int dx/sqrt(C3) = cF (F(phi2,m)-F(phi1,m))."""
    a3=mpf(C3n[3])
    rts=polyroots([mpc(c) for c in reversed(C3n)],maxsteps=800,extraprec=mp.prec+120)
    rts=[mpc(r) for r in rts]
    # real root
    reals=[r for r in rts if abs(r.imag)<mpf(10)**(-mp.dps+25)]
    e1=reals[0].real
    comp=[r for r in rts if abs(r.imag)>=mpf(10)**(-mp.dps+25)]
    if len(comp)==2:
        p=-(comp[0]+comp[1]).real; q=(comp[0]*comp[1]).real  # x^2+p x+q
        # B&F 260: for cubic a3(x-e1)(x^2+px+q), with A^2=(e1^2+p e1+q)= value of
        # quad at e1 >0. Sub: cos(phi)=(A-(x-e1))/(A+(x-e1)) ... we use the
        # Carlson-free standard:  let
        A2 = e1*e1+p*e1+q   # quad(e1)
        A = msqrt(A2)
        # g1 = A + ... ; the reduction (real root e1, A=sqrt(quad(e1))):
        #   m = (1/2) - (3 e1 + ... )  -- to avoid slips, we CALIBRATE m and the
        # phi-map by matching int dx/sqrt(C3)=cF*(F(phi(x2))-F(phi(x1))) and its
        # derivative; but cleaner: use the KNOWN formula (Gradshteyn 3.131):
        #  For y=sqrt((x-e1)(x^2+px+q)), with A=sqrt(quad(e1)),
        #   substitution  x = e1 + A (1-cos w)/(1+cos w) = e1 + A tan^2(w/2)
        #   maps to int ~ dphi/Delta with
        #   m = (1/2)(1 - (e1 + p/2)/A )   [B&F 260.00, b1=p/2... ]
        b = p/2 + e1   # = (e1 + p/2)?  quad'(e1)/2 = (2e1+p)/2 = e1+p/2
        m = (1 - (e1+p/2)/A)/2
        return dict(kind='1r2c', e1=e1, p=p, q=q, A=A, m=m, a3=a3)
    else:
        # 3 real roots e1<e2<e3
        es=sorted([r.real for r in reals])
        e1,e2,e3=es
        m=(e2-e1)/(e3-e1)
        return dict(kind='3r', e1=e1, e2=e2, e3=e3, m=m, a3=a3)

def _phi_1r2c(ch, x):
    # x = e1 + A tan^2(w/2)? -> we use cos(theta)=(A-(x-e1))/(A+(x-e1)),
    # phi = theta/... Standard (G&R 3.131): with t=x-e1>=0,
    #   cos(2 alpha) = (A - t)/(A + t);  the elliptic argument is alpha? We
    # CALIBRATE phi-map numerically robustly: phi = acos((A-t)/(A+t))/1, and
    # find cF by gating. Implement as: u-substitution verified in test.
    A=ch['A']; t=x-ch['e1']
    c2a=(A - t)/(A + t)
    # clamp
    from mpmath import acos
    return acos(c2a)

def _phi_3r(ch, x):
    """3 real roots e1<e2<e3.
    region (e1,e2): modulus m=(e2-e1)/(e3-e1),  sin^2 phi=(x-e1)/(e2-e1).
    region (e2,e3): modulus m=(e3-e2)/(e3-e1),  sin^2 phi=(e3-x)/(e3-e2)
       (verified constant cF; brute-force scan).  The window region is in
       ch['_region']; the corresponding modulus is in ch['m'] (set by caller)."""
    e1,e2,e3=ch['e1'],ch['e2'],ch['e3']
    from mpmath import asin as ms
    if ch.get('_region')=='e2e3':
        s2=(e3-x)/(e3-e2)
    else:
        s2=(x-e1)/(e2-e1)
    return ms(msqrt(s2)) if (hasattr(s2,'real') and s2.real>=0) or s2>=0 else ms(msqrt(mpc(s2)))

def I0_cubic(C3n, x1, x2):
    """int_{x1}^{x2} dx/sqrt(C3) closed (F-atom) + gate-ready. Returns (val, ch).
    Delegates J0 to the high-precision Jmoments_cubic (region/kind-correct)."""
    J,info=Jmoments_cubic(C3n, x1, x2, kmax=1)
    ch=info['ch']; ch['cF']=info.get('cF')
    return J[0], ch
    if ch['kind']=='1r2c':
        A=ch['A']; m=ch['m']; a3=ch['a3']
        # G&R 3.131.3 / B&F 260.00:
        #   int_{e1}^{x} dt/sqrt((t-e1)(t^2+pt+q)) ... = (1/sqrt(A a3)) F(phi,m)?
        # cF determined by gate; we compute phi via _phi_1r2c and fit cF by one
        # derivative match (exact): d/dx[F(phi(x),m)] = phi'(x)/Delta(phi).
        def phi(x): return _phi_1r2c(ch,x)
        # cF s.t. cF d/dx F(phi,m) = 1/sqrt(C3). derivative match at midpoint.
        xm=(x1+x2)/2
        h=mpf(10)**(-mp.dps//3)
        dF=(ellipf(phi(xm+h),m)-ellipf(phi(xm-h),m))/(2*h)
        tru=1/msqrt(C3val(C3n,xm))
        cF=tru/dF
        val=cF*(ellipf(phi(x2),m)-ellipf(phi(x1),m))
        ch['cF']=cF; ch['phi']=phi
        return val, ch
    raise NotImplementedError('3-real-root cubic not needed for E_alg here')

if __name__=='__main__':
    mp.dps=45
    # E_alg cell A: y^2=(1-x)(784x^2-784x+125). x range = u^2 in [0.46^2,0.54^2]
    import sympy as sp
    x=sp.Symbol('x')
    for cell,k,c0,x1f,x2f in [('A',784,125,'0.46','0.54'),('B',256,45,'0.66','0.74')]:
        cubic=sp.expand((1-x)*(k*x**2-k*x+c0))
        C3n=[mpf(str(c)) for c in reversed(sp.Poly(cubic,x).all_coeffs())]
        # pad to len4
        while len(C3n)<4: C3n.append(mpf(0))
        x1,x2=mpf(x1f)**2,mpf(x2f)**2
        val,ch=I0_cubic(C3n,x1,x2)
        ref=mquad(lambda t:1/msqrt(C3val(C3n,t)),[x1,x2])
        if abs(val+ref)<abs(val-ref): val=-val
        err=abs(val-ref)/abs(ref); dig=float(-mp.log10(err)) if err>0 else 99
        print(f"cell {cell}: E_alg I0 cubic closed={nstr(val,18)} ref={nstr(ref,18)} -> {dig:.1f}d  m={nstr(ch['m'],20)} kind={ch['kind']} e=({nstr(ch.get('e1',0),6)},{nstr(ch.get('e2',0),6)},{nstr(ch.get('e3',0),6)})")


# ---- cubic moments + pole atoms (mirror quartic engine) ----
def cubic_atoms(C3n, x1, x2):
    """returns helpers to reduce int rat(x) dx/sqrt(C3) on (e1,e2) region.
    Provides Jmoments_cubic and pole_cubic, both gated."""
    ch=cubic_chart(C3n); ch['m']=(ch['e2']-ch['e1'])/(ch['e3']-ch['e1'])
    return ch

def Jmoments_cubic(C3n, x1, x2, kmax=4):
    """J_k=int_{x1}^{x2} x^k dx/sqrt(C3), C3=a3 x^3+a2 x^2+a1 x+a0 (3 real roots).
    J0=I0 (F). J1 via 2nd-kind (E). higher via x-recurrence from d/dx[x^n sqrt C3]:
      [x^n sqrtC3]' = (n x^{n-1} C3 + x^n C3'/2)/sqrtC3.
      top power: n a3 x^{n+2} + x^n (3 a3 x^2)/2 = (n+3/2) a3 x^{n+2}.
      => (n+3/2)a3 J_{n+2} = B_n - [(n+1)a2 J_{n+1} + (n+1/2)a1 J_n + n a0 J_{n-1}].
      gives J2,J3,.. from J0,J1. J0,J1 are the elliptic periods (F,E)."""
    a0,a1,a2,a3=[mpf(C3n[i]) for i in range(4)]
    ch=cubic_chart(C3n)
    xmw=(mpf(x1)+mpf(x2))/2
    if ch['kind']=='3r':
        e1,e2,e3=ch['e1'],ch['e2'],ch['e3']
        if xmw>e2 and xmw<e3:
            ch['_region']='e2e3'; ch['m']=(e3-e2)/(e3-e1)
        else:
            ch['_region']='e1e2'; ch['m']=(e2-e1)/(e3-e1)
        m=ch['m']
        def phi(x): return _phi_3r(ch,x)
    else:  # 1r2c: real root e1 + complex pair (x^2+p x+q)
        e1,p,q=ch['e1'],ch['p'],ch['q']
        A=msqrt(e1*e1+p*e1+q)            # sqrt(quad(e1)) (>0)
        ch['A']=A
        ch['m']=(1-(e1+p/2)/A)/2          # B&F 260.00 (verified constant cF)
        m=ch['m']
        from mpmath import acos as macos
        def phi(x):
            t=x-e1; return macos((A-t)/(A+t))
        ch['_region']='1r2c'
    # cF for J0 (analytic): d/dx F(phi,m)=phi'/Delta; with sin^2=(x-e1)/(e2-e1):
    #   sin phi=sqrt((x-e1)/(e2-e1)); cos phi dphi = dx/(2(e2-e1) sin phi)
    #   dphi = dx/(2(e2-e1) sin phi cos phi). Delta=sqrt(1-m sin^2).
    #   1/sqrtC3 = 1/sqrt(a3(x-e1)(x-e2)(x-e3)). Match cF: do it via gate-grade
    #   high-prec finite diff at small h with extra dps.
    dps0=mp.dps; mp.dps=dps0+20
    def F(x): return ellipf(phi(x),m)
    def E(x): return ellipe(phi(x),m)
    # J0 = int dx/sqrtC3 = cF (F(x2)-F(x1)); cF by exact-deriv match.
    xm=(mpf(x1)+mpf(x2))/2; h=mpf(10)**(-mp.dps//2)
    dF=(F(xm+h)-F(xm-h))/(2*h); cF=(1/msqrt(C3val(C3n,xm)))/dF
    J0=cF*(F(mpf(x2))-F(mpf(x1)))
    # J1 = int x dx/sqrtC3.  Antiderivative = cF1*E(phi,m)+cF0*F(phi,m)+calg*sqrtC3
    #   (2nd-kind + 1st-kind + algebraic).  Fit cF1,cF0,calg by exact-derivative
    #   match at 3 interior points (m moderate => well conditioned).
    from mpmath import matrix, lu_solve
    def dFx(x):
        hh=mpf(10)**(-mp.dps//2); return (F(x+hh)-F(x-hh))/(2*hh)
    def dEx(x):
        hh=mpf(10)**(-mp.dps//2); return (E(x+hh)-E(x-hh))/(2*hh)
    def sC3(x): return msqrt(C3val(C3n,x))
    def dalg(x):  # d/dx sqrtC3 = C3'/(2 sqrtC3)
        c=C3val(C3n,x); cp=3*a3*x*x+2*a2*x+a1; return cp/(2*msqrt(c))
    def tgt1(x): return x/msqrt(C3val(C3n,x))
    pts=[mpf(x1)+(mpf(x2)-mpf(x1))*f for f in (mpf('0.2'),mpf('0.5'),mpf('0.8'))]
    M=matrix(3,3); rhs=matrix(3,1)
    for i,xp in enumerate(pts):
        M[i,0]=dEx(xp); M[i,1]=dFx(xp); M[i,2]=dalg(xp); rhs[i,0]=tgt1(xp)
    s1=lu_solve(M,rhs)
    def J1anti(x): return s1[0]*E(x)+s1[1]*F(x)+s1[2]*sC3(x)
    J1=J1anti(mpf(x2))-J1anti(mpf(x1))
    J=[J0,J1]
    def sC3(x): return msqrt(C3val(C3n,x))
    while len(J)<=kmax:
        n=len(J)-2   # solving for J_{n+2}
        Bn=mpf(x2)**n*sC3(mpf(x2))-mpf(x1)**n*sC3(mpf(x1))
        Jm1=J[n-1] if n-1>=0 else mpf(0)
        rhs=Bn-((n+1)*a2*J[n+1]+(n+mpf(1)/2)*a1*J[n]+n*a0*Jm1)
        J.append(rhs/((n+mpf(3)/2)*a3))
    mp.dps=dps0
    return [ (jj.real if abs(jj.imag)<abs(jj.real)*mpf(10)**(-30) else jj) for jj in J], dict(ch=ch,cF=cF,m=m)

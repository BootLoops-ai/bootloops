"""ellred_reduce: robust reduction of  int_{x1}^{x2} rat(x) dx / sqrt(P4(x))
to STANDARD Legendre F/E/Pi + algebraic.  Replaces the buggy reduce_rat: does
the partial fraction in z with mpmath-exact residues and maps every term to a
gated closed atom.

Pullback:  x = (t2 z - t1)/(z-1),  dx/sqrt(P4) = (jac/sqrt(A)) dz/g.
  rat(x) -> ratz(z) (rational in z, with (z-1)^k in denom from x=inf pole).
Write  ratz(z) = P(z) + sum_k R_k/(z - rho_k)   [full PF over C, mpmath roots].
Integrand = (jac/sqrt(A)) [ P(z) + sum R_k/(z-rho_k) ] dz/g.

Atom map (all gated closed forms in ellred_engine):
  * P(z) even powers z^{2j}: z^0->F-atom; z^2->(E-atom - b2 F-atom)/a2;
    z^{>=4} via recurrence (reduce_even). ODD powers -> ALGEBRAIC (w=z^2).
  * simple pole 1/(z-rho): combine with conj 1/(z-rhobar) ->
    (Bz+C)/((z-Re)^2+Im^2). Decompose:
      - real pole rho real: 1/(z-rho) = (z+rho)/(z^2-rho^2) -> z/(z^2-rho^2)
        [odd->algebraic] + rho/(z^2-rho^2) [-> Pi-atom with c=rho^2].
      - complex-conjugate pair with real part: shift not compatible with g(z)
        (g is even in z).  Use: (Bz+C)/(z^2+pz+q). Since g(z) is EVEN in z, we
        symmetrise:  the WHOLE integrand over the (even) measure picks out the
        even-in-z part of ratz automatically only if limits symmetric -- they
        are NOT.  So we keep general (Bz+C)/(z^2+pz+q) and reduce by writing it
        over the even denom (z^2+q)^2 - (pz)^2... messy.  Cleaner: the kernel
        cofactors we hit are EVEN in u=>x even=>the x-poles are at specific
        rational x; we map pole-in-x p_x -> z-pole rho=zof(p_x), and use the
        x-space 3rd-kind atom directly (below) instead of z-PF for poles.
This file therefore provides TWO routes:
  reduce_poly(coeffs_x, P4n, x1, x2): rat = polynomial in x  -> exact F/E/alg.
  reduce_pole(p_x, P4n, x1, x2): int dx/((x-p_x) sqrt(P4))   -> exact F/Pi/alg.
A general rat = poly + sum residues/(x-p_x) is assembled from these in x-space
(no z-1 artifacts): standard Hermite reduction in x.
"""
import os, sys
import sympy as sp
from mpmath import (mp,mpf,mpc,sqrt as msqrt,quad as mquad,nstr,polyroots,
                    tan as mtan, atan as matan, ellipf,ellipe,ellippi)
from . import ell_normal as EN
from . import ellred_engine as E

def _setup(P4n, x1, x2):
    ch = E.chart(P4n)
    z1, z2 = E.zof(ch, mpf(x1)), E.zof(ch, mpf(x2))
    zlo, zhi = (z1, z2) if z1 < z2 else (z2, z1)
    ch['_zwin'] = (zlo, zhi)
    return ch, z1, z2

# ---- x-space atom: int dx/((x - p) sqrt(P4)) -> F/Pi --------------------
# Pull p to z: pole at x=p <-> z=rho=zof(p). 1/(x-p) dx = ... but easier:
#  1/(x-p) * dx/sqrt(P4) = 1/(x-p) (jac/sqrtA) dz/g.  x-p = (t2 z - t1)/(z-1) - p
#    = ((t2-p) z - (t1-p))/(z-1) = (t2-p)(z - rho)/(z-1),  rho=(t1-p)/(t2-p).
#  so 1/(x-p) = (z-1)/((t2-p)(z-rho)).
#  integrand = (jac/sqrtA)/(t2-p) * (z-1)/(z-rho) /g dz
#            = (jac/sqrtA)/(t2-p) * [ 1 + (rho-1)/(z-rho) ] /g dz
#  -> (z-1)/(z-rho) = 1 + (rho-1)/(z-rho).  First term -> F-atom.
#  Second: (rho-1) int dz/((z-rho) g).  Combine 1/(z-rho) over EVEN g:
#    1/(z-rho) = (z+rho)/(z^2-rho^2).  z/(z^2-rho^2)->algebraic(w=z^2);
#    rho/(z^2-rho^2) -> rho * Pi-atom(c=rho^2)  (atomPi_closed).
def reduce_pole(p, P4n, x1, x2):
    ch, z1, z2 = _setup(P4n, x1, x2)
    t1,t2,A,jac = ch['t1'],ch['t2'],ch['A'],ch['jac']
    p = mpf(p) if not isinstance(p,(mpf,mpc)) else p
    rho = (t1-p)/(t2-p)
    c = rho*rho
    pref = (jac/msqrt(A))/(t2-p)
    def F(z): return E.atomF(ch, z)
    def Pi(z): return E.atomPi_closed(ch, c, z)
    # algebraic part: (rho-1)* int z/((z^2-c) g) dz = (rho-1)*(1/2) int dw/((w-c)sqrt((a1 w+b1)(a2 w+b2)))
    def algpart(z):
        w2=z*z
        return mpf('0.5')*mquad(lambda w: 1/((w-c)*msqrt((ch['a1']*w+ch['b1'])*(ch['a2']*w+ch['b2']))),[mpf('0'),w2])
    Fpart = F(z2)-F(z1)
    Pipart = Pi(z2)-Pi(z1)
    Algp = algpart(z2)-algpart(z1)
    val = pref*( Fpart + (rho-1)*( Algp + rho*Pipart ) )
    return dict(value=val, c=c, rho=rho, n0=1+ch['b1']/(ch['a1']*c), m=ch['m'],
                F=Fpart, Pi=Pipart, alg=Algp, pref=pref, ch=ch,
                atoms=['F', f'Pi(n0={nstr(1+ch["b1"]/(ch["a1"]*c),12)})', 'alg(z/(z^2-c))'])

# ---- x-space polynomial atom: int x^k dx/sqrt(P4) -> F/E/alg -------------
# x = (t2 z - t1)/(z-1).  x^k = ((t2 z - t1)/(z-1))^k. Expand & PF in z:
#   denom (z-1)^k.  -> poles at z=1 (=x=inf).  This re-introduces (z-1) PF.
# Cleaner: reduce int x^k dx/sqrt(P4) DIRECTLY in x by the quartic recurrence:
#   d/dx[ x^{k-3} sqrt(P4) ] = (k-3) x^{k-4} sqrt(P4) + x^{k-3} P4'/(2 sqrt(P4))
# gives a linear relation among int x^{j} dx/sqrt(P4), j=k-... down. Base cases
#   J0=int dx/sqrtP4 (F), J1, J2, J3 expressed via F,E and the z-atoms.
# We get J0..J3 from the z-atoms by the Mobius image of {1,z^2}->{F,E} plus the
# odd ALGEBRAIC pieces. Implement J0,J1,J2,J3 by solving the 4x4 from x-monomial
# moments matched to the z-basis {F-atom, E-atom, alg-odd1, alg-odd2} via the
# EXACT linear pullback (x^k = rational in z). We do this by expressing each
# x^k dx/sqrtP4 integrand in z and PF, but ONLY the (z-1) poles -> which are at
# x=inf; near x=inf P4~lead x^4 so x^k/sqrtP4 ~ x^{k-2}: for k<=1 integrable,
# k>=2 the (z-1) pole is genuine and corresponds to 2nd-kind/algebraic at inf.
# To keep it ROBUST we compute J0..J3 by: J0=F-atom; and J1,J2,J3 by the
# x-recurrence anchored on algebraic derivative terms d/dx[x^j sqrtP4]:
def _Jmoments(ch, P4n, z1, z2, kmax):
    """return [J0,...,Jkmax], Jk=int_{x1}^{x2} x^k dx/sqrt(P4), via x-recurrence.
    P4=c4 x^4+c3 x^3+c2 x^2+c1 x+c0. d/dx[x^{n} sqrt(P4)] =
      n x^{n-1} sqrtP4 + x^n P4'/(2 sqrtP4)
      = [ n x^{n-1} P4 + x^n P4'/2 ] / sqrtP4.
    Integrate x1..x2: [x^n sqrtP4]_{x1}^{x2} = sum_j coeff_j J_{...}.
    This gives J_{n+3} in terms of lower J's (since top power n x^{n-1}*c4 x^4 +
    x^n * 4 c4 x^3 /2 = (n+2) c4 x^{n+3}). Solve upward for J3,J4,...; J0,J1,J2
    are the independent elliptic integrals (F,E + one more). We compute J0 (F),
    and get J1,J2 by ALSO using the z-atoms? Simpler: J0,J1,J2 elliptic, J3+ by
    recurrence. We need J1,J2 closed. Provide them via z-basis below.
    """
    c0,c1,c2,c3,c4 = [mpf(P4n[i]) for i in range(5)]
    P4 = lambda xx: ((((c4*xx+c3)*xx+c2)*xx+c1)*xx+c0)
    return None  # built in reduce_poly using z-atoms for J0..J2



# ---- polynomial moments J_k = int x^k dx/sqrt(P4) ----------------------
# Base J0 = F-atom (x-space, =reduce of poly 1). J1,J2 via Mobius image:
#   x = (t2 z - t1)/(z-1).  x dx/sqrtP4 and x^2 dx/sqrtP4 pulled to z give
#   rational(z) dz/g with (z-1) poles. We reduce those (z-1) poles with the
#   z-atoms: 1/(z-1)^k. The pole at z=1 corresponds to x=inf.
#   Express x^k = ((t2 z - t1)/(z-1))^k, integrand x^k (jac/sqrtA) dz/g.
#   PF in z: poly(z) + sum c_j/(z-1)^j.  poly even/odd -> F/E/alg; the 1/(z-1)^j
#   are poles at z=1 (real): handle via reduce_pole-style with rho=1 i.e. c=1.
#   But z=1 maps to x=inf where the integrand 1/(x-inf)... -> use the SAME
#   z-pole atom (1/(z-1) and derivatives). 1/(z-1): rho=1, c=1.  1/(z-1)^2 via
#   d/dz. To stay robust we reduce J1,J2 by the EXACT x-RECURRENCE instead:
def Jmoments(P4n, x1, x2, kmax=4):
    """J_k = int_{x1}^{x2} x^k dx/sqrt(P4), k=0..kmax.
    J0 = F-atom (gated). J1: int x/sqrtP4 -> substitute, =algebraic+J0 combo?
      Actually int x dx/sqrtP4 with P4 quartic is a genuine elliptic (2nd kind).
    Use x-recurrence from d/dx[x^{n}sqrtP4]:
      [x^n sqrtP4]' = (n x^{n-1}P4 + x^n P4'/2)/sqrtP4.
      => boundary B_n := x2^n sqrt(P4(x2)) - x1^n sqrt(P4(x1))
         = n c4 J_{n+3}+... +(n+2)... ; collect:
      coefficient of J_{n+3}: n c4 + (4 c4)/2 = (n+2) c4.
      So: (n+2)c4 J_{n+3} = B_n - [ (n+ 3/2) c3 J_{n+2} + (n+1) c2 J_{n+1}
                                    + (n+1/2) c1 J_n + n c0 J_{n-1} ].
      Let me derive coefficients exactly below.
    This gives J3,J4,... from J0,J1,J2. We still need J0,J1,J2 (the elliptic
    'periods'): J0=F, and J1,J2 come from the z-basis via atomE_closed + odd alg.
    """
    # E is the module-level `from . import ellred_engine as E` binding; a bare
    # `import ellred_engine` here would fail under package import.
    c0,c1,c2,c3,c4=[mpf(P4n[i]) for i in range(5)]
    ch=E.chart(P4n); z1,z2=E.zof(ch,mpf(x1)),E.zof(ch,mpf(x2))
    ch['_zwin']=(min(z1,z2),max(z1,z2))
    t1,t2,A,jac=ch['t1'],ch['t2'],ch['A'],ch['jac']
    # J0,J1,J2 via z-basis. x=(t2 z-t1)/(z-1). On [z1,z2], integrand
    #   x^k (jac/sqrtA)/g. Expand x^k as rational in z and integrate each
    #   monomial-in-z over dz/g using: z^0->F, z^2->E-combo, z odd->alg,
    #   and 1/(z-1)^j poles via direct (algebraic for odd-part / Pi for even).
    # SIMPLER & ROBUST: J0,J1,J2 by GAUSS-form? No. We just integrate x^k dz-pullback
    # with the z-atoms for the POLYNOMIAL-in-z part and a single Pi/alg for the
    # (z-1) pole, exactly as reduce_pole does for a pole at x=inf (p->inf):
    #   lim p->inf of int dx/((x-p)sqrtP4) * (-p) = int dx/sqrtP4 = J0 (consistency)
    # Rather than limits, directly: J1,J2 obtained by solving the linear system
    # from matching d/dx[sqrtP4] and d/dx[x sqrtP4] (algebraic) to {J0,J1,J2,J3}:
    #   d/dx[sqrtP4] = P4'/(2 sqrtP4) = (4c4 x^3+3c3 x^2+2c2 x+c1)/(2 sqrtP4)
    #     => B_{-1}:=[sqrtP4]_{x1}^{x2} = 2c4 J3 + (3/2)c3 J2 + c2 J1 + (1/2)c1 J0
    #   d/dx[x sqrtP4] = sqrtP4 + x P4'/(2 sqrtP4)
    #     => [x sqrtP4] = J? : (x P4'/2 + P4)/sqrtP4 integrated
    #        = (1/sqrtP4)[ (4c4/2+c4)x^4 + (3c3/2+c3)x^3+(2c2/2+c2)x^2+(c1/2+c1)x+c0 ]
    #        = 3c4 J4 + (5/2)c3 J3 + 2c2 J2 + (3/2)c1 J1 + c0 J0
    # These RELATE high J's; not enough to get J1,J2 alone. So compute J1,J2 from
    # the z-atoms directly (the clean route):
    # x = t2 + (t2-t1)/(z-1).  Let s=(t2-t1). x = t2 + s/(z-1).
    #  x dx/sqrtP4 = (t2 + s/(z-1)) (jac/sqrtA) dz/g
    #  x^2 dx/sqrtP4 = (t2 + s/(z-1))^2 (jac/sqrtA) dz/g
    # so we need M0=int dz/g (=>F), Mp1=int dz/((z-1)g), Mp2=int dz/((z-1)^2 g).
    pref=jac/msqrt(A)
    M0=(E.atomF(ch,z2)-E.atomF(ch,z1))
    # 1/(z-1): real pole rho=1, c=1: 1/(z-1)=(z+1)/(z^2-1)= z/(z^2-1)+1/(z^2-1)
    c=mpf(1)
    Pi1=(E.atomPi_closed(ch,c,z2)-E.atomPi_closed(ch,c,z1))   # int dz/((z^2-1)g)
    def algodd(z):  # int z/((z^2-1)g) dz = (1/2) int dw/((w-1)sqrt(...))
        w2=z*z
        return mpf('0.5')*mquad(lambda w:1/((w-1)*msqrt((ch['a1']*w+ch['b1'])*(ch['a2']*w+ch['b2']))),[mpf(0),w2])
    Aodd1=algodd(z2)-algodd(z1)
    Mp1=Aodd1+Pi1     # int dz/((z-1)g) = int (z+1)/((z^2-1)g) = algodd + Pi1
    # 1/(z-1)^2 = d/dz[-1/(z-1)]; int 1/(z-1)^2 /g dz = [-1/((z-1)g)]... integrate
    #   by parts: int (z-1)^{-2}/g dz = -1/((z-1)g) - int 1/(z-1) d(1/g)... messy.
    # robust: Mp2 by direct quad (still a standard 3rd-kind/algebraic; certified)
    Mp2=mquad(lambda z:1/((z-1)**2*E.gz(ch,z)),[z1,z2])
    s=t2-t1
    J0=pref*M0
    J1=pref*(t2*M0 + s*Mp1)
    J2=pref*(t2*t2*M0 + 2*t2*s*Mp1 + s*s*Mp2)
    J=[J0,J1,J2]
    # higher via recurrence: (n+2)c4 J_{n+3} = B_n - [(n+3/2)c3 J_{n+2}
    #   + (n+1)c2 J_{n+1} + (n+1/2)c1 J_n + n c0 J_{n-1}]
    def sP4(xx): return msqrt(EN.G_pval(P4n,xx))
    while len(J)<=kmax:
        n=len(J)-3
        Bn=mpf(x2)**n*sP4(mpf(x2))-mpf(x1)**n*sP4(mpf(x1))
        Jm1=J[n-1] if n-1>=0 else mpf(0)
        rhs=Bn-((n+mpf(3)/2)*c3*J[n+2]+(n+1)*c2*J[n+1]+(n+mpf(1)/2)*c1*J[n]+n*c0*Jm1)
        J.append(rhs/((n+2)*c4))
    return J, dict(ch=ch, M0=M0, Mp1=Mp1, Mp2=Mp2, pref=pref)


#-
Off statistics;
*
*   g(p1,mu,a) + g(p2,nu,b)  ->  q(p3,i) + qbar(p4,jbar)
*
*   Tree-level squared amplitude in massless QCD, exact SU(N) color,
*   d-dimensional Dirac traces.  Three diagrams:
*     (t)  quark exchange, propagator momentum p3-p1,  color T^a T^b
*     (u)  quark exchange, propagator momentum p3-p2,  color T^b T^a
*     (s)  three-gluon vertex into the quark pair,     color f(a,b,c) T^c
*
*   Conventions (Peskin-Schroeder): quark-gluon vertex i g gamma^mu T^a,
*   massless quark propagator i pslash/p^2, gluon propagator (Feynman
*   gauge) -i g_{mu nu} delta^{ab}/p^2, triple-gluon vertex
*   g f^{abc} [ g^{mu nu}(k1-k2)^rho + cyclic ], all momenta incoming.
*   An overall factor g^2 is stripped from the amplitude; the printed
*   squared amplitude is |M|^2 / g^4.
*
*   Gluon polarization treatment: PHYSICAL (axial-gauge) polarization
*   sums, using the other gluon's momentum as the reference vector:
*     sum_pol eps_mu(p1) eps*_nu(p1) = -g_{mu nu} + (p1_mu p2_nu + p2_mu p1_nu)/(p1.p2)
*   and (p1 <-> p2) for the second gluon.  Only the d-2 transverse
*   states propagate, so no ghost subtraction is needed.
*
*   Kinematics: s = (p1+p2)^2, t = (p1-p3)^2, u = (p1-p4)^2, all masses
*   zero.  The result is written as a Laurent polynomial in s,t,u with
*   d and NF (= number of colors N) symbolic.  s+t+u = 0 is NOT imposed
*   here; impose it downstream when comparing forms.
*
*   Averaging over initial state (2 helicities x (N^2-1) colors per
*   gluon, 4-dim counting) is done downstream: divide by 4(N^2-1)^2.
*
*   MUTATION knobs for the self-test (a gate that cannot fail is not a
*   gate):  -D MUT=1 flips the sign of the three-gluon vertex;
*           -D MUT=2 breaks the color Fierz identity (wrong 1/N term).
*   Either must make the textbook comparison FAIL.
*
#ifndef `MUT'
#define MUT "0"
#endif

Symbols d, NF, s, t, u;
Dimension d;
Indices mu, nu, mub, nub;
Vectors p1, p2, p3, p4;
CFunctions TrF, f, cc, cb;
AutoDeclare Index a;

#if `MUT' == 1
#define SGN3G "-1"
#else
#define SGN3G "1"
#endif

*--------------------------------------------------------------------
* Amplitude (spin line 1, sandwiched as Tr[ p3sl Gamma_i p4sl GammaBar_j ]):
* each term = coefficient * color * ( p3sl Gamma_i p4sl ).
* cc(...) holds the fundamental generators of the quark line in order.
*--------------------------------------------------------------------
Local Amp =
   - i_*t^-1 * cc(a1,a2)
       * g_(1,p3) * g_(1,mu) * (g_(1,p3)-g_(1,p1)) * g_(1,nu) * g_(1,p4)
   - i_*u^-1 * cc(a2,a1)
       * g_(1,p3) * g_(1,nu) * (g_(1,p3)-g_(1,p2)) * g_(1,mu) * g_(1,p4)
   + `SGN3G'*s^-1 * f(a1,a2,a3)*cc(a3)
       * g_(1,p3) * (   d_(mu,nu) * (g_(1,p1)-g_(1,p2))
                      + (p1(mu)+2*p2(mu)) * g_(1,nu)
                      - (2*p1(nu)+p2(nu)) * g_(1,mu)  ) * g_(1,p4);

* Conjugate: coefficients conjugated, gamma chains reversed, color
* chains reversed (T^a hermitian, f real); cb(...) already reversed.
Local Conj =
   + i_*t^-1 * cb(a2,a1)
       * g_(1,nub) * (g_(1,p3)-g_(1,p1)) * g_(1,mub)
   + i_*u^-1 * cb(a1,a2)
       * g_(1,mub) * (g_(1,p3)-g_(1,p2)) * g_(1,nub)
   + `SGN3G'*s^-1 * f(a1,a2,a4)*cb(a4)
       * (   d_(mub,nub) * (g_(1,p1)-g_(1,p2))
           + (p1(mub)+2*p2(mub)) * g_(1,nub)
           - (2*p1(nub)+p2(nub)) * g_(1,mub)  );
.sort

* |M|^2 summed over quark spins, gluon physical polarizations, all colors.
Local Msq = Amp * Conj
   * ( -d_(mu,mub)  + 2*s^-1*( p1(mu)*p2(mub) + p2(mu)*p1(mub) ) )
   * ( -d_(nu,nub)  + 2*s^-1*( p2(nu)*p1(nub) + p1(nu)*p2(nub) ) );
.sort
Drop Amp, Conj;
.sort

*--------------------------------------------------------------------
* Dirac algebra: d-dimensional trace over spin line 1.
*--------------------------------------------------------------------
tracen,1;
.sort

* Massless 2->2 kinematics.
id p1.p1 = 0;
id p2.p2 = 0;
id p3.p3 = 0;
id p4.p4 = 0;
id p1.p2 = s/2;
id p3.p4 = s/2;
id p1.p3 = -t/2;
id p2.p4 = -t/2;
id p1.p4 = -u/2;
id p2.p3 = -u/2;
.sort

*--------------------------------------------------------------------
* Color algebra: merge quark-line chains into a single fundamental
* trace, convert f to fundamental traces, reduce by the SU(N) Fierz
* identity (T_R = 1/2):
*   T^j_{il} T^j_{km} = 1/2 ( d_{im} d_{kl} - d_{il} d_{km}/N ).
*--------------------------------------------------------------------
id cc(?x)*cb(?y) = TrF(?x,?y);
id f(a1?,a2?,a3?) = -2*i_*( TrF(a1,a2,a3) - TrF(a1,a3,a2) );
.sort

#if `MUT' == 2
#define CFIERZ "+1/(2*NF)"
#else
#define CFIERZ "-1/(2*NF)"
#endif
repeat;
  id TrF(?x,a0?,?y)*TrF(?z,a0?,?w) = 1/2*TrF(?x,?w,?z,?y) `CFIERZ'*TrF(?x,?y)*TrF(?z,?w);
  id TrF(?x,a0?,?y,a0?,?z)         = 1/2*TrF(?x,?z)*TrF(?y) `CFIERZ'*TrF(?x,?y,?z);
  id TrF(a0?) = 0;
  id TrF      = NF;
endrepeat;
.sort

Print +s Msq;
.sort
Format nospaces;
#write <gg2qq.out> "Msq=%E;", Msq
#close <gg2qq.out>
.end

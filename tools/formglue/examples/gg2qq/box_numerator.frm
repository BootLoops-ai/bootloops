#-
Off statistics;
*
* 1-loop probe: full numerator algebra for one representative box of
* g(p1,mu,a) g(p2,nu,b) -> q(p3) qbar(p4), interfered with the full
* conjugate tree amplitude.
*
* Box (quark line p3 -> p4, internal gluon exchanged across the line):
*   ubar(p3) [ig ga^al T^e] iS(k+p3) [ig ga^mu T^a] iS(k+p3-p1)
*            [ig ga^nu T^b] iS(k-p4) [ig ga^be T^e'] v(p4)
*   x  (-i d_(al,be) delta_{ee'} / D4),
* loop denominators  D1=(k+p3)^2, D2=(k+p3-p1)^2, D3=(k-p4)^2, D4=k^2.
* Overall factor: i^4 (vertices) x i^3 (quark props) x (-i) (gluon
* prop) = -1; a factor g^4 is stripped.  Color: T^e T^a T^b T^e.
*
* Output: the numerator N(k) of  sum_{spins,pols,colors} M_box M_tree^*
* x (D1 D2 D3 D4) / g^6, as a polynomial in the loop scalar products
* kk=k.k, kp1=k.p1, kp2=k.p2, kp3=k.p3 (k.p4 eliminated by momentum
* conservation), coefficients rational in s,t,u and symbolic in d, NF.
* Tree denominators appear as monomials t^-1, u^-1, s^-1.
* Same conventions and physical polarization sums as gg2qq.frm.
*
Symbols d, NF, s, t, u, kk, kp1, kp2, kp3;
Dimension d;
Indices mu, nu, mub, nub, al;
Vectors p1, p2, p3, p4, k;
CFunctions TrF, f, cc, cb;
AutoDeclare Index a;

* Box amplitude, stripped:  - (spinor chain) * T^e T^a T^b T^e
Local Abox =
   - cc(a3,a1,a2,a3)
     * g_(1,p3) * g_(1,al)
     * (g_(1,k)+g_(1,p3))          * g_(1,mu)
     * (g_(1,k)+g_(1,p3)-g_(1,p1)) * g_(1,nu)
     * (g_(1,k)-g_(1,p4))          * g_(1,al)
     * g_(1,p4);

* Conjugate full tree (same object as in gg2qq.frm).
Local Conj =
   + i_*t^-1 * cb(a2,a1)
       * g_(1,nub) * (g_(1,p3)-g_(1,p1)) * g_(1,mub)
   + i_*u^-1 * cb(a1,a2)
       * g_(1,mub) * (g_(1,p3)-g_(1,p2)) * g_(1,nub)
   + s^-1 * f(a1,a2,a4)*cb(a4)
       * (   d_(mub,nub) * (g_(1,p1)-g_(1,p2))
           + (p1(mub)+2*p2(mub)) * g_(1,nub)
           - (2*p1(nub)+p2(nub)) * g_(1,mub)  );
.sort

Local Nbox = Abox * Conj
   * ( -d_(mu,mub)  + 2*s^-1*( p1(mu)*p2(mub) + p2(mu)*p1(mub) ) )
   * ( -d_(nu,nub)  + 2*s^-1*( p2(nu)*p1(nub) + p1(nu)*p2(nub) ) );
.sort
Drop Abox, Conj;
.sort

tracen,1;
.sort
#$nraw = termsin_(Nbox);
#write "TERMS after d-dim trace = `$nraw'"

* Kinematics: external on-shell dots and momentum conservation.
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
id k.p4 = k.p1 + k.p2 - k.p3;
id k.k  = kk;
id k.p1 = kp1;
id k.p2 = kp2;
id k.p3 = kp3;
.sort

* Color reduction (SU(N) Fierz, T_R = 1/2), as in gg2qq.frm.
id cc(?x)*cb(?y) = TrF(?x,?y);
id f(a1?,a2?,a3?) = -2*i_*( TrF(a1,a2,a3) - TrF(a1,a3,a2) );
.sort
repeat;
  id TrF(?x,a0?,?y)*TrF(?z,a0?,?w) = 1/2*TrF(?x,?w,?z,?y) - 1/(2*NF)*TrF(?x,?y)*TrF(?z,?w);
  id TrF(?x,a0?,?y,a0?,?z)         = 1/2*TrF(?x,?z)*TrF(?y) - 1/(2*NF)*TrF(?x,?y,?z);
  id TrF(a0?) = 0;
  id TrF      = NF;
endrepeat;
.sort
#$nfin = termsin_(Nbox);
#write "TERMS integrand-ready    = `$nfin'"

Format nospaces;
#write <box_numerator.out> "Nbox=%E;", Nbox
#close <box_numerator.out>
.end

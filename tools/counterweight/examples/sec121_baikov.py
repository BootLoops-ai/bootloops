#!/usr/bin/env python3
"""
Sector-121 (props D1,D4,D5,D6,D7; double-bubble) connection matrix A_s, A_t via the
2-loop Baikov representation, exact in sympy. m^2=1.
D1=k1^2, D4=k2^2-1, D5=(k2-q)^2-1, D6=(k1-k2)^2-1, D7=(k1-k2+p3)^2-1, q=p1+p2+p3.
For the maximal cut / DE we need the Baikov poly P = Gram(k1,k2 . k1,k2,p1,p2,p3) expressed
in the z_i = D_i. We then compute d/ds, d/dt of each block master by:
  dI/dx = gamma * <(dP/dx)/P> with the bracket expanded into z-shifts (dots), gamma=(d-L-E-1)/2.
This is standard; we implement P from scalar products.
"""
import sympy as sp
s,t,eps = sp.symbols('s t eps')
d = 4-2*eps
# scalar product basis. Loop k1,k2; ext p1,p2,p3 with p_i^2=0, p1.p2=s/2,p1.p3=t/2,p2.p3=-(s+t)/2.
# ISP/ propagator scalar products for sector 121 (5 props) - 2-loop has 2*3+3=... 
# independent scalar products among (k1,k2) with (k1,k2,p1,p2,p3): k1^2,k2^2,k1.k2, k1.pi(3), k2.pi(3) = 9 sps.
# 5 propagators -> 4 ISPs. Building full Baikov is involved; verify we can at least form P det.
k1k1,k2k2,k1k2 = sp.symbols('k1k1 k2k2 k1k2')
k1p1,k1p2,k1p3 = sp.symbols('k1p1 k1p2 k1p3')
k2p1,k2p2,k2p3 = sp.symbols('k2p1 k2p2 k2p3')
# Gram matrix of (k1,k2,p1,p2,p3)
p1p1=p2p2=p3p3=0
p1p2=s/2; p1p3=t/2; p2p3=-(s+t)/2
G = sp.Matrix([
 [k1k1,k1k2,k1p1,k1p2,k1p3],
 [k1k2,k2k2,k2p1,k2p2,k2p3],
 [k1p1,k2p1,p1p1,p1p2,p1p3],
 [k1p2,k2p2,p1p2,p2p2,p2p3],
 [k1p3,k2p3,p1p3,p2p3,p3p3]])
P = G.det()
print("Baikov P built; #terms:", len(sp.Poly(sp.expand(P), k1k1,k2k2,k1k2,k1p1,k1p2,k1p3,k2p1,k2p2,k2p3).terms()))
# d/ds and d/dt of P at fixed scalar products:
dPds = sp.diff(P,s); dPdt=sp.diff(P,t)
print("dP/ds nonzero:", dPds!=0, " dP/dt nonzero:", dPdt!=0)

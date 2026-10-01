#!/usr/bin/env python3
"""General exact audit for fiber .ms systems (2 roamer vars).
Parses rows with sympy; dummy rows (no x) = chart-infinity lines (excluded).
Reports: affine census, deg Res, arrangement valuations, crit count."""
import sys, itertools
import sympy as sp
def audit(fn, p=999999937):
    txt=[l.strip().rstrip(',') for l in open(fn).read().split('\n')]
    vars_=txt[0].split(',')
    x1,x2=sp.symbols('x1 x2')
    ms=sp.symbols(','.join(v for v in vars_ if v.startswith('m')))
    rows=[l for l in txt[2:] if l]
    nl=len(ms)
    L=[];U=[];dummy=[]
    for i in range(nl):
        e=sp.sympify(rows[i])
        mi=ms[i]
        A=e.coeff(mi*x1); B=e.coeff(mi*x2)
        C=sp.expand(e - A*mi*x1 - B*mi*x2)
        Cc=C.coeff(mi); u=-sp.expand(C - Cc*mi)
        A,B,Cc,u=int(A),int(B),int(Cc),int(u)
        if A==0 and B==0: dummy.append(i); L.append(None); U.append(u); continue
        L.append((A,B,Cc)); U.append(u)
    act=[i for i in range(nl) if L[i] is not None]
    x,y=sp.symbols('x y')
    ells={i:(L[i][0]*x+L[i][1]*y+L[i][2]) for i in act}
    P1=sp.expand(sum(U[i]*L[i][0]*sp.prod([ells[j] for j in act if j!=i]) for i in act))
    P2=sp.expand(sum(U[i]*L[i][1]*sp.prod([ells[j] for j in act if j!=i]) for i in act))
    P1=sp.Poly(P1,x,y,modulus=p); P2=sp.Poly(P2,x,y,modulus=p)
    R=sp.Poly(sp.resultant(P1,P2,y),x,modulus=p)
    pts={}; par=0
    for i,j in itertools.combinations(act,2):
        A1,B1,C1=L[i]; A2,B2,C2=L[j]
        den=(A1*B2-A2*B1)%p
        if den==0: par+=1; continue
        xi=((-C1*B2+C2*B1)*pow(den,-1,p))%p
        yi=((-A1*C2+A2*C1)*pow(den,-1,p))%p
        pts.setdefault((xi,yi),set()).update([i,j])
    from collections import Counter
    cen=Counter(len(v) for v in pts.values())
    acc=0
    for xi in {q[0] for q in pts}:
        Rr=R; v=0; div=sp.Poly([1,-xi],x,modulus=p)
        while True:
            q_,r_=sp.div(Rr,div)
            if r_.is_zero: v+=1; Rr=q_
            else: break
        acc+=v
    nlines=len(act)
    return dict(file=fn, lines=nlines, dummies=len(dummy), parallel_pairs=par,
                census=dict(sorted(cen.items())), degP=(P1.total_degree(),P2.total_degree()),
                degR=R.degree(), arr_val=acc, crit=R.degree()-acc)
if __name__=='__main__':
    for fn in sys.argv[1:]:
        print(audit(fn), flush=True)

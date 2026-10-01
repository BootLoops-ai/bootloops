#!/usr/bin/env python3
"""Fast exact audit via evaluate-interpolate resultants over F_p."""
import sys, itertools
import sympy as sp
import numpy as np
def audit(fn, p=999999937):
    txt=[l.strip().rstrip(',') for l in open(fn).read().split('\n')]
    vars_=txt[0].split(',')
    x1,x2=sp.symbols('x1 x2')
    ms=sp.symbols(','.join(v for v in vars_ if v.startswith('m')))
    rows=[l for l in txt[2:] if l]
    nl=len(ms)
    L=[];U=[]
    for i in range(nl):
        e=sp.sympify(rows[i]); mi=ms[i]
        A=e.coeff(mi*x1); B=e.coeff(mi*x2)
        C=sp.expand(e-A*mi*x1-B*mi*x2).coeff(mi)
        u=-sp.expand(e-A*mi*x1-B*mi*x2-C*mi)
        A,B,C,u=int(A),int(B),int(C),int(u)
        if A==0 and B==0: L.append(None); U.append(u); continue
        L.append((A%p,B%p,C%p)); U.append(u%p)
    act=[i for i in range(nl) if L[i] is not None]
    n=len(act)
    # P1(x,y), P2: for fixed x evaluate coefficients in y then resultant via
    # Euclidean algorithm on univariate polys over F_p. Simpler: R(x0) =
    # Res_y(P1(x0,y),P2(x0,y)) via poly gcd-free resultant (Euclid with tracking)
    def polymulmod(a,b):
        r=np.zeros(len(a)+len(b)-1,dtype=object)
        for i,ai in enumerate(a):
            if ai:
                for j,bj in enumerate(b):
                    r[i+j]=(r[i+j]+ai*bj)%p
        return list(r)
    def resultant_uni(f,g):
        # f,g lists of coeffs (ascending), over F_p; standard Euclid resultant
        f=[c%p for c in f]; g=[c%p for c in g]
        def deg(h):
            d=len(h)-1
            while d>=0 and h[d]==0: d-=1
            return d
        res=1
        df,dg=deg(f),deg(g)
        while True:
            if dg<0: return 0
            if dg==0: return res*pow(g[0],df,p)%p
            # f = q*g + r
            fdeg=deg(f)
            if fdeg<dg:
                f,g=g,f; df,dg=deg(f),deg(g)
                if df%2==1 and dg%2==1: res=(-res)%p
                continue
            lc=g[dg]; lcinv=pow(lc,-1,p)
            r=f[:]
            for k in range(fdeg-dg,-1,-1):
                c=r[dg+k]*lcinv%p
                if c:
                    for j in range(dg+1):
                        r[j+k]=(r[j+k]-c*g[j])%p
            dr=deg(r)
            res=res*pow(lc,fdeg-dr,p)%p
            if df%2==1 and dg%2==1: res=(-res)%p  # swap sign convention
            f,g=g,r[:dr+1] if dr>=0 else [0]
            df,dg=deg(f),deg(g)
            if dg<0: return 0
    def P_at(x0):
        # build P1(y), P2(y) coeff lists at x=x0
        ells=[]
        for i in act:
            A,B,C=L[i]
            ells.append([ (A*x0+C)%p, B%p ])  # a + b*y
        P1=[0]; P2=[0]
        for i,ii in enumerate(act):
            pr=[1]
            for j,jj in enumerate(act):
                if j!=i: pr=polymulmod(pr,ells[j])
            A,B,_=L[ii]
            t1=[(U[ii]*A%p)*c%p for c in pr]
            t2=[(U[ii]*B%p)*c%p for c in pr]
            if len(t1)>len(P1): P1=P1+[0]*(len(t1)-len(P1))
            if len(t2)>len(P2): P2=P2+[0]*(len(t2)-len(P2))
            for k,c in enumerate(t1): P1[k]=(P1[k]+c)%p
            for k,c in enumerate(t2): P2[k]=(P2[k]+c)%p
        return P1,P2
    # R(x) degree <= (n-1)*(n-1)*2? deg P in y: n-1; coeffs deg in x: n-1
    # deg Res <= 2*(n-1)^2 ~ safe bound
    D=2*(n-1)*(n-1)+10
    xs=list(range(1,D+2))
    Rv=[]
    for x0 in xs:
        f,g=P_at(x0)
        Rv.append(resultant_uni(f,g))
    # interpolate R over F_p -> get degree & valuations at arrangement x's.
    # Exact Newton interpolation in F_p (integer arithmetic throughout: a
    # rational-coefficient route via interpolation over Q does not coerce
    # into GF(p) and is far slower at D ~ hundreds of sample points).
    m=len(xs)
    dd=list(Rv)  # divided differences, in place
    for j in range(1,m):
        for i in range(m-1,j-1,-1):
            dd[i]=(dd[i]-dd[i-1])*pow((xs[i]-xs[i-j])%p,-1,p)%p
    coeffs=[0]*m  # ascending coefficient list of the Newton form, expanded
    cur=[1]       # running product (x-xs[0])...(x-xs[k-1])
    for kk in range(m):
        for idx,c in enumerate(cur):
            coeffs[idx]=(coeffs[idx]+dd[kk]*c)%p
        new=[0]*(len(cur)+1)  # cur *= (x - xs[kk])
        for idx,c in enumerate(cur):
            new[idx]=(new[idx]-c*xs[kk])%p
            new[idx+1]=(new[idx+1]+c)%p
        cur=new
    xg=sp.Symbol('x')
    P=sp.Poly(list(reversed(coeffs)), xg, modulus=p)
    degR=P.degree()
    # arrangement points
    ptsA={}
    for i,j in itertools.combinations(act,2):
        A1,B1,C1=L[i]; A2,B2,C2=L[j]
        den=(A1*B2-A2*B1)%p
        if den==0: continue
        xi=((-C1*B2+C2*B1)*pow(den,-1,p))%p
        yi=((-A1*C2+A2*C1)*pow(den,-1,p))%p
        ptsA.setdefault((xi,yi),set()).update([i,j])
    from collections import Counter
    cen=Counter(len(v) for v in ptsA.values())
    acc=0
    for xi in {q[0] for q in ptsA}:
        Rr=P; v=0; div=sp.Poly([1,-xi],xg,modulus=p)
        while True:
            q_,r_=sp.div(Rr,div)
            if r_.is_zero: v+=1; Rr=q_
            else: break
        acc+=v
    return dict(file=fn,lines=n,census=dict(sorted(cen.items())),degR=degR,
                arr_val=acc, crit=degR-acc)
if __name__=='__main__':
    for fn in sys.argv[1:]:
        print(audit(fn), flush=True)

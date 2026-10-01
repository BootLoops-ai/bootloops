#!/usr/bin/env python3
"""Affine flat semilattice + Whitney/Mobius for bounded regions.
b = |chi_A(1)| (Zaslavsky). Planes given as (a1..ak, c): a.x + c = 0."""
import sys
import sympy as sp

def affine_chi(planes, k):
    n=len(planes)
    def flat_of(idx):
        """solve intersection; return (closure set, rank) or None if empty."""
        A=sp.Matrix([list(planes[i][:k]) for i in idx])
        b=sp.Matrix([[-planes[i][k]] for i in idx])
        r=A.rank(); Ab=A.row_join(b)
        if Ab.rank()>r: return None
        # closure: all planes containing the flat
        # flat = particular solution + nullspace
        try:
            x0=A.pinv()*b  # rational pinv ok? use solve on independent rows
        except Exception:
            return None
        # verify
        if sp.simplify(A*x0-b)!=sp.zeros(len(idx),1): 
            # use generic solve
            xs=sp.symbols(f'x0:{k}')
            eqs=[sum(planes[i][j]*xs[j] for j in range(k))+planes[i][k] for i in idx]
            s=sp.solve(eqs,xs,dict=True)
            if not s: return None
            x0=sp.Matrix([[s[0].get(v,v)] for v in xs])
        N=A.nullspace()
        cl=set()
        for i in range(n):
            a=sp.Matrix(1,k,list(planes[i][:k])); c=planes[i][k]
            v=(a*x0)[0,0]+c
            if sp.simplify(v)!=0: continue
            if all(sp.simplify((a*w)[0,0])==0 for w in N): cl.add(i)
        return frozenset(cl), r
    allf={}
    for i in range(n):
        f=flat_of([i]); allf[f[0]]=1
    for r in range(2,k+1):
        new={}
        prev=[S for S,rr in allf.items() if rr==r-1]
        for S in prev:
            for j in range(n):
                if j in S: continue
                f=flat_of(list(S)+[j])
                if f is None: continue
                cl,rr=f
                if rr==r and cl not in allf and cl not in new: new[cl]=rr
        allf.update(new)
    order=sorted(allf.items(), key=lambda kv: kv[1])
    mu={}; rank={}
    for S,r in order: rank[S]=r
    mu0=1
    for S,r in order:
        s=mu0
        for T,rT in order:
            if rT<r and T<S: s+=mu[T]
        mu[S]=-s
    q=sp.Symbol('q')
    chi=q**k
    for S,r in order: chi+=mu[S]*q**(k-r)
    b=chi.subs(q,1)
    reg=((-1)**k)*chi.subs(q,-1)
    return sp.expand(chi), b, reg

if __name__=='__main__':
    fn=sys.argv[1]
    txt=[l.strip().rstrip(',') for l in open(fn).read().split('\n')]
    vars_=txt[0].split(',')
    xs=[v for v in vars_ if v.startswith('x')]
    k=len(xs)
    X=sp.symbols(' '.join(xs)); ms=sp.symbols(','.join(v for v in vars_ if v.startswith('m')))
    rows=[l for l in txt[2:] if l]
    planes=[]
    for i in range(len(ms)):
        e=sp.sympify(rows[i]); mi=ms[i]
        co=[e.coeff(mi*x) for x in X]
        rest=sp.expand(e-sum(c*mi*x for c,x in zip(co,X)))
        cc=rest.coeff(mi)
        if all(c==0 for c in co): continue
        planes.append(tuple([sp.Rational(c) for c in co]+[sp.Rational(cc)]))
    print(f'{fn}: {len(planes)} active planes, dim {k}')
    chi,b,reg=affine_chi(planes,k)
    print('chi_A(q) =', chi)
    print('bounded (|chi(1)|):', abs(b), ' regions:', reg)

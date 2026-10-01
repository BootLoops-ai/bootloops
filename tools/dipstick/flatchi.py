#!/usr/bin/env python3
"""Exact chi_top of P^k-complement of a hyperplane arrangement via the flat
lattice + Mobius + char poly. chi_top(P^k minus A) = d/dq [chihat(q)/(q-1)]
at q=1 = chihat'(1) (since chihat(1)=0)."""
import itertools, sys
from fractions import Fraction as Fr
import sympy as sp

def flats_and_chi(normals, k):
    """normals: list of (k+1)-tuples (central arrangement in C^{k+1})."""
    n=len(normals)
    # rank-2..k flats: iterate: start from hyperplane pairs, close under join;
    # represent a flat by the set of ALL hyperplanes containing it
    def closure(idx):
        M=sp.Matrix([list(normals[i]) for i in idx])
        r=M.rank()
        ns=M.nullspace()
        out=set()
        for i in range(n):
            v=sp.Matrix(1,k+1,list(normals[i]))
            if all((v*b)[0,0]==0 for b in ns):
                out.add(i)
        return frozenset(out), r
    allflats={}  # closure-set -> rank
    for i in range(n): allflats[frozenset(
        closure([i])[0])]=1
    # BFS joins
    for r in range(2, k+2):
        new={}
        prev=[S for S,rr in allflats.items() if rr==r-1]
        for S in prev:
            for j in range(n):
                if j in S: continue
                cl,rr=closure(list(S)+[j])
                if rr==r and cl not in allflats and cl not in new:
                    new[cl]=r
        allflats.update(new)
        if not new: break
    # Mobius
    order=sorted(allflats.items(), key=lambda kv: kv[1])
    mu={frozenset(): 1}
    rankof={frozenset(): 0}
    for S,r in order: rankof[S]=r
    for S,r in order:
        s=1  # mu of bottom
        for T,rT in order:
            if rT<r and T<S: s+=mu[T]
        mu[S]=-s
    q=sp.Symbol('q')
    chihat=q**(k+1)
    for S,r in order:
        chihat+=mu[S]*q**(k+1-r)
    chihat=sp.expand(chihat)
    chiproj=sp.simplify(sp.cancel(chihat/(q-1)))
    chitop=sp.diff(chihat,q).subs(q,1)
    return chihat, chiproj, chitop

def planes_from_points(C, k):
    """all hyperplanes through k of the given points in P^k (as normal vecs)"""
    out=[]
    for idx in itertools.combinations(range(len(C)), k):
        M=sp.Matrix([list(C[i]) for i in idx])
        ns=M.nullspace()
        assert len(ns)==1
        v=ns[0].T
        den=sp.lcm([sp.fraction(sp.Rational(x))[1] for x in v])
        out.append(tuple(sp.Rational(x)*den for x in v))
    return out

if __name__=='__main__':
    # CALIBRATION (mathematical cages, self-contained):
    #   15 lines from 6 generic pts in P^2 -> chi_top 42
    #   --quick stops there (seconds-scale; the battery leg uses it);
    #   default continues to the 35-planes-from-7-pts P^3 cage (~minutes).
    import random
    quick = '--quick' in sys.argv[1:]
    rng=random.Random(99)
    def rnd(): return Fr(rng.randint(-15,15), rng.randint(1,5))
    while True:
        C=[(1,0,0),(0,1,0),(0,0,1),(1,1,1),(1,rnd(),rnd()),(1,rnd(),rnd())]
        C=[tuple(sp.Rational(x) for x in p) for p in C]
        if all(sp.Matrix([C[i] for i in t]).det()!=0
               for t in itertools.combinations(range(6),3)): break
    N=planes_from_points(C,2)
    ch,cp,ct=flats_and_chi(N,2)
    print('CAL 15-lines-from-6-points cage: chi_top =', ct, '(expect 42)')
    if quick:
        sys.exit(0)
    # 35 planes from 7 generic pts in P^3
    while True:
        C=[(1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1),(1,1,1,1),
           (1,rnd(),rnd(),rnd()),(1,rnd(),rnd(),rnd())]
        C=[tuple(sp.Rational(x) for x in p) for p in C]
        if all(sp.Matrix([C[i] for i in t]).det()!=0
               for t in itertools.combinations(range(7),4)): break
    N=planes_from_points(C,3)
    ch,cp,ct=flats_and_chi(N,3)
    print('35-planes-from-7-points cage: chihat =', sp.factor(ch))
    print('35-planes-from-7-points cage: chi_top(P^3 complement) =', ct)

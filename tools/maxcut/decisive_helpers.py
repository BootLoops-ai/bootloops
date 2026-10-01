#!/usr/bin/env python3
r"""maxcut_decisive_helpers.py — decidable-algebra kernels (see maxcut_decisive.py)."""
import sympy as sp, random

def gram_baikov(loops, ext, ext_gram, props, prefer_free=None, kin_syms=()):
    """Return (B, free_ISPs).  props = [(mom_dict, m²), ...]."""
    bv={}
    for i,a in enumerate(loops):
        for b in loops[i:]: bv[(a,b)]=sp.symbols(f"{a}{b}")
    for a in loops:
        for b in ext: bv[(a,b)]=sp.symbols(f"{a}{b}")
    def dot(a,b):
        if a in ext and b in ext: return ext_gram.get((a,b), ext_gram.get((b,a)))
        return bv[(a,b)] if (a,b) in bv else bv[(b,a)]
    def sq(q):
        return sp.expand(sum(ca*cb*dot(a,b) for a,ca in q.items()
                                            for b,cb in q.items()))
    bvars=list(bv.values())
    eqs=[sp.expand(sq(q)-m) for q,m in props]
    if prefer_free:
        pf=[v for v in bvars if str(v) in prefer_free]
        sol=sp.solve(eqs,[v for v in bvars if v not in pf],dict=True)[0]
    else:
        sol=dict(zip(bvars, list(sp.linsolve(eqs,bvars))[0]))
    ALL=loops+ext; n=len(ALL); M=sp.zeros(n,n)
    for i,a in enumerate(ALL):
        for j,b in enumerate(ALL): M[i,j]=dot(a,b)
    B=sp.expand(M.det().subs(sol))
    free=sorted(set(B.free_symbols)-set(kin_syms),key=str)
    return B,free

def curve_j(P, x):
    """j-invariant of y²=P(x) via cubic depression / binary-quartic I,J."""
    p=sp.Poly(sp.expand(P),x); d=p.degree()
    if d==3:
        a,b,c,e=[p.nth(i) for i in (3,2,1,0)]
        pp=(3*a*c-b*b)/(3*a*a); qq=(2*b**3-9*a*b*c+27*a*a*e)/(27*a**3)
        D=-16*(4*pp**3+27*qq*qq)
        return sp.oo if D==0 else sp.nsimplify(1728*(-4*pp)**3/D)
    if d==4:
        a,b,c,dd,e=[p.nth(i) for i in (4,3,2,1,0)]
        I=12*a*e-3*b*dd+c*c; J=72*a*c*e+9*b*c*dd-27*a*dd*dd-27*b*b*e-2*c**3
        D=4*I**3-J*J
        return sp.oo if D==0 else sp.nsimplify(1728*4*I**3/D)
    return None

def carrier_slice_genus(B, isp, kin_sub, n_slices=4, seed=2026):
    """Highest-deg ISP kept symbolic; spectators sliced; sqfree deg → genus."""
    rng=random.Random(seed)
    Bn=sp.expand(B.subs(kin_sub))
    free=[v for v in isp if v in Bn.free_symbols]
    if not free: return dict(genus=0, degs={}, carrier=None, slices=[])
    degs={str(v):sp.Poly(Bn,v).degree() for v in free}
    carrier=max(free,key=lambda v:degs[str(v)])
    spect=[v for v in free if v is not carrier]
    slices=[]
    for _ in range(n_slices):
        sub={v:sp.Rational(rng.randint(-7,7),rng.randint(1,7)) for v in spect}
        Pn=sp.Poly(sp.expand(Bn.subs(sub)),carrier)
        c,fl=sp.factor_list(Pn.as_expr(),carrier)
        sqf=sum(sp.Poly(f,carrier).degree() for f,m in fl
                if m%2==1 and carrier in f.free_symbols)
        j=curve_j(Pn.as_expr(),carrier) if sqf in (3,4) else None
        slices.append(dict(deg=Pn.degree(),sqfree=sqf,j=str(j) if j is not None else None))
    g=max((0 if s['sqfree']<=2 else 1 if s['sqfree']<=4 else (s['sqfree']-1)//2)
          for s in slices)
    jset={s['j'] for s in slices if s['j'] not in (None,'oo')}
    return dict(genus=g, degs=degs, carrier=str(carrier), slices=slices,
                j_distinct=len(jset))

def disc_chain(B, isp, kin_sub, maxsteps=12):
    """Iterated discriminant on deg-2 ISPs (never touches deg-≥3 carriers)."""
    cur=sp.expand(B.subs(kin_sub)); chain=[]
    for _ in range(maxsteps):
        free=[v for v in isp if v in sp.sympify(cur).free_symbols]
        if not free: break
        d2=[v for v in free if sp.Poly(cur,v).degree()==2]
        if not d2: break
        v=d2[0]; A,b,C=[sp.Poly(cur,v).nth(i) for i in (2,1,0)]
        cur=sp.factor(sp.expand(b*b-4*A*C)); chain.append(str(v))
    rem=[v for v in isp if v in sp.sympify(cur).free_symbols]
    out=dict(drained=chain, residual_vars=[str(v) for v in rem])
    if len(rem)==1:
        c,fl=sp.factor_list(cur,rem[0])
        out['sqfree']=sum(sp.Poly(f,rem[0]).degree() for f,m in fl
                          if m%2==1 and rem[0] in f.free_symbols)
    elif not rem:
        out['sqfree']=0
    else:
        # multi-factor chain (BaikovNestedBeta mode='multi' greedy step):
        # eliminate any var in exactly one factor at deg≤2 via its disc.
        _,fl=sp.factor_list(cur)
        facs=[f for f,m in fl if m%2==1 and set(f.free_symbols)&set(rem)]
        for _ in range(maxsteps):
            live=[v for v in isp if any(v in f.free_symbols for f in facs)]
            if len(live)<=1: break
            pick=None
            for v in live:
                act=[(j,sp.Poly(f,v).degree()) for j,f in enumerate(facs)
                     if v in f.free_symbols]
                if len(act)==1 and act[0][1]<=2:
                    pick=(v,act[0][0],act[0][1]); break
            if pick is None: break
            v,j,dg=pick; f=facs.pop(j)
            if dg==2:
                A,b,C=[sp.Poly(f,v).nth(i) for i in (2,1,0)]
                _,nfl=sp.factor_list(sp.expand(b*b-4*A*C))
                facs.extend(g for g,m in nfl if m%2==1
                            and set(g.free_symbols)&set(isp))
            chain.append(str(v))
        live=[v for v in isp if any(v in f.free_symbols for f in facs)]
        out['drained']=chain; out['residual_vars']=[str(v) for v in live]
        out['n_factors']=len(facs)
        if len(live)<=1:
            sq=sum(sp.Poly(f,live[0]).degree() for f in facs) if live else 0
            out['sqfree']=sq
        else:
            out['max_factor_pervar_deg']=max(
                sp.Poly(f,v).degree() for f in facs for v in live
                if v in f.free_symbols)
    return cur,out

def find_theta_op(a, order, degz_max=8, extra=25):
    """Krylov θ-annihilator of Σaₙzⁿ (mod-p LP-poly analogue: tools/pf_rank.jl)."""
    z=sp.Symbol('z')
    for degz in range(1,degz_max+1):
        cols=[(j,k) for k in range(degz+1) for j in range(order+1)]
        N=len(cols)+order+extra
        if len(a)<=N: return None
        M=sp.zeros(N+1,len(cols))
        for ci,(j,k) in enumerate(cols):
            for m in range(N+1):
                mm=m-k
                if 0<=mm: M[m,ci]=a[mm]*(mm**j if mm or j==0 else 0)
        ns=M.nullspace()
        if ns:
            v=ns[0]; L=sp.ilcm(*[sp.Rational(c).q for c in v])
            iv=[int(c*L) for c in v]; g=sp.igcd(*iv) or 1
            iv=[c//g for c in iv]
            Pj={j:sp.Integer(0) for j in range(order+1)}
            for ci,(j,k) in enumerate(cols): Pj[j]+=iv[ci]*z**k
            if Pj[order]!=0: return Pj,degz
    return None

def minimal_pf_order(a, max_order=6, degz_max=8):
    for k in range(1,max_order+1):
        r=find_theta_op(a,k,degz_max)
        if r: return k,r
    return None,None

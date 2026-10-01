#!/usr/bin/env python3
r"""tornheim.py — fast Mordell–Tornheim–Witten sums + depth-2 MZVs (general tool).

The primitives are general:
depth-2 MZVs `mzv2(s,t)`, 2/3-fold Tornheim `T2(a,b,c)`/`MT3(a,b,c,d)`, Euler
harmonic sum `eulerH(s)`, and the signed lattice sums `W3`/`W4` show up in
massless multi-loop propagators, banana/sunrise ε-expansions, and lattice-sum
problems — not just MGFs.  All caches are dps-keyed (mp.mpf lru_cache footgun).


W_S(p;q) := Σ_{m∈(ℤ\0)^S:Σm=0} ∏|m_r|^{-p_r} · (Σ|m_r|)^{-q}

|S|=3: sign-sector → 2·Σ_{3 perms} T(α,β;γ+q)/2^q,  T(a,b;c)=Σ_{m,n≥1}1/(m^a n^b (m+n)^c).
       T reduced via T(a,b,c)=T(a-1,b,c+1)+T(a,b-1,c+1) to base cases:
         T(0,b,c)=ζ(b)ζ(c)-ζ(b+c)-ζ(b,c)    [b≥2]
         T(1,0,c)=T(0,1,c)=Σ_{k≥2}H_{k-1}/k^c = EulerH(c)-ζ(c+1)
         T(1,1,c)=2·T(1,0,c+1)               [via 1/(mn)=(1/(m+n))(1/m+1/n)]
       ζ(b,c) via Richardson-accelerated hurwitz sum (~4s @ 200d).
       EulerH(s):=Σ_{k≥1}H_k/k^s = (1+s/2)ζ(s+1)-(1/2)Σ_{i=1}^{s-2}ζ(i+1)ζ(s-i)  [Euler].

|S|=4: sign-sector patterns (1,3) and (2,2):
  (1,3): 8 patterns → Σ_{n∈ℤ_{>0}³} 1/(n₁^α n₂^β n₃^γ (Σn)^δ (2Σn)^q) = MT₃(α,β,γ;δ+q)/2^q.
  (2,2): 6 patterns → Σ_K S_{αβ}(K)·S_{γδ}(K)/(2K)^q,  S_{αβ}(K)=Σ_{a=1}^{K-1}1/(a^α(K-a)^β).
  MT₃ via recursion (a,b,c;d)→Σ(…-1…;d+1) to a=0, then MT₃(0,b,c;d)=Σ_M S₂(b,c;M)ζ_H(d,M+1).
  S₂ has EXACT partial-fraction S₂(b,c;M)=Σᵢeᵢ·H_{M-1}^{(i)}/M^{b+c-i}; head-sum to Mmax +
  asymptotic tail (H,ζ_H asymptotics ⇒ Σ(ln M)^l/M^p ⇒ (-1)^l ζ_H^{(l)}(p,Mmax+1)).
  Full working precision (verified: W₄((1,1,1,1);1)=30ζ₅−12ζ₂ζ₃ @150d; dps-stable).
"""
import mpmath as mp
from mpmath import mpf, zeta
from functools import lru_cache

# ---- double MZV ζ(s,t)=Σ_{m>n≥1}1/(m^s n^t), s≥2 ----
_MZV2={}
def mzv2(s,t):
    key=(s,t,mp.mp.dps)  # dps-keyed: a dps-free key cache-poisons across precision changes
    if key in _MZV2: return _MZV2[key]
    v=mp.nsum(lambda n: mp.hurwitz(s,int(n)+1)/n**t,[1,mp.inf],method='richardson')
    _MZV2[key]=v; return v

def eulerH(s):
    """Σ_{k≥1}H_k/k^s (Euler)."""
    v=(mpf(1)+mpf(s)/2)*zeta(s+1)
    for i in range(1,s-1): v-=zeta(i+1)*zeta(s-i)/2
    return v

# ---- 2-fold Tornheim T(a,b,c) ----
_TC={}
def T2(a,b,c):
    if a>b: a,b=b,a
    key=(a,b,c,mp.mp.dps)  # dps-keyed: same cache-poison class as _MZV2
    if key in _TC: return _TC[key]
    if a==0 and b==0:
        raise ValueError("T(0,0,c) diverges")
    if a==0 and b==1:
        v=eulerH(c)-zeta(c+1)  # T(0,1,c)=Σ_{k≥2}H_{k-1}/k^c
    elif a==0:
        v=zeta(b)*zeta(c)-zeta(b+c)-mzv2(b,c)
    elif a==1 and b==1:
        v=2*T2(0,1,c+1)
    else:
        v=T2(a-1,b,c+1)+T2(a,b-1,c+1)
    _TC[key]=v; return v

def W3(p,q):
    """W₃(p₁,p₂,p₃;q). Sign-sector: 2·Σ_{i=lone} T(p_a,p_b; p_i+q)/2^q."""
    p1,p2,p3=p
    perms=[(p2,p3,p1),(p1,p3,p2),(p1,p2,p3)]
    return 2*sum(T2(a,b,g+q) for a,b,g in perms)/mpf(2)**q

# ==== high-precision MT3 / conv22 via exact S₂ partial-fraction + asymptotic tail ====
# S₂(b,c;M)=Σ_{j=1}^{M-1}1/(j^b(M-j)^c) has EXACT partial-fraction decomposition:
#   S₂(b,c;M) = Σ_{i=1}^{max(b,c)} e_i(b,c)·H_{M-1}^{(i)}/M^{b+c-i},
#   e_i = C(b+c-i-1,c-1)[i≤b] + C(b+c-i-1,b-1)[i≤c].
# Each H_{M-1}^{(i)} has a known asymptotic in 1/M (with a single ln M for i=1).
# ⇒ tails Σ_{M>Mmax}(…) reduce to (-1)^m ζ_H^{(m)}(p,Mmax+1) exactly.
from math import comb as _comb

def _S2_pf_coeffs(b,c):
    m=max(b,c); e=[0]*(m+1)
    for i in range(1,m+1):
        v=0
        if i<=b: v+=_comb(b+c-i-1,c-1)
        if i<=c: v+=_comb(b+c-i-1,b-1)
        e[i]=v
    return e

def _H_asym(i,Pmax):
    """H_{M-1}^{(i)} ~ Σ c·(ln M)^l/M^p, return {(l,p):c} for p≤Pmax."""
    A={}
    if i==1:
        A[(1,0)]=mpf(1); A[(0,0)]=mp.euler
        if Pmax>=1: A[(0,1)]=mpf(-1)/2
        k=1
        while 2*k<=Pmax:
            A[(0,2*k)]=-mp.bernoulli(2*k)/(2*k); k+=1
    else:
        A[(0,0)]=zeta(i)
        if i-1<=Pmax: A[(0,i-1)]=mpf(-1)/(i-1)
        if i  <=Pmax: A[(0,i)]  =mpf(-1)/2
        poch=mpf(i); k=1
        while i+2*k-1<=Pmax:
            A[(0,i+2*k-1)]=-mp.bernoulli(2*k)/mp.factorial(2*k)*poch
            poch*=(i+2*k-1)*(i+2*k); k+=1
    return A

def _zetaHp1_asym(d,Pmax):
    """ζ_H(d,M+1)=ζ_H(d,M)-1/M^d ~ {(0,p):c}."""
    A={}
    if d-1<=Pmax: A[(0,d-1)]=mpf(1)/(d-1)
    if d  <=Pmax: A[(0,d)]  =mpf(-1)/2
    poch=mpf(d); k=1
    while d+2*k-1<=Pmax:
        A[(0,d+2*k-1)]=mp.bernoulli(2*k)/mp.factorial(2*k)*poch
        poch*=(d+2*k-1)*(d+2*k); k+=1
    return A

def _S2_asym(b,c,Pmax):
    e=_S2_pf_coeffs(b,c); A={}
    for i in range(1,len(e)):
        if e[i]==0: continue
        for (l,p),v in _H_asym(i,Pmax).items():
            pp=p+b+c-i
            if pp>Pmax: continue
            A[(l,pp)]=A.get((l,pp),mpf(0))+e[i]*v
    return A

def _asym_mul(A,B,Pmax):
    C={}
    for (la,pa),ca in A.items():
        for (lb,pb),cb in B.items():
            p=pa+pb
            if p>Pmax: continue
            k=(la+lb,p); C[k]=C.get(k,mpf(0))+ca*cb
    return C

_tailz_cache={}
def _tail_zeta(l,p,a,guard,dps0):
    """(-1)^l ζ_H^{(l)}(p,a) at dps0+guard (mpmath's ζ_H has ~10^{-dps} ABSOLUTE
    error for moderate p, so elevate). Cached across calls."""
    key=(l,p,a,dps0)
    if key in _tailz_cache and _tailz_cache[key][0]>=guard:
        return _tailz_cache[key][1]
    mp.mp.dps=dps0+guard
    z=mp.zeta(p,a) if l==0 else ((-1)**l)*mp.zeta(p,a,l)
    mp.mp.dps=dps0
    _tailz_cache[key]=(guard,z); return z

def _tail_sum(A,Mmax):
    """Σ_{M>Mmax} Σ c·(ln M)^l/M^p via Hurwitz-ζ derivatives."""
    dps0=mp.mp.dps; la=float(mp.log10(Mmax+1))
    s=mpf(0)
    for (l,p),c in A.items():
        if c==0: continue
        guard=max(0,int(float(mp.log10(abs(c))))+2)+int((p-1)*la)+10
        s+=c*_tail_zeta(l,p,Mmax+1,guard,dps0)
    return s

def _default_Mmax_Pmax():
    # Asymptotic coeffs grow ~10^{~1·p} (from (d)_{2k-1} Pochhammer), so contrib
    # decays only ~10^{-(log10 Mmax - 1)·p}. Need generous Pmax; keep Mmax large
    # enough that the divergent series is well inside its decreasing regime.
    dps=mp.mp.dps
    Mmax=max(500,3*dps)
    Pmax=2*int(dps/float(mp.log10(Mmax)))+20
    return Mmax,Pmax

# ---- MT₃(a,b,c;d) via recursion to a=0 then head+tail ----
_MT3_cache={}
def MT3(a,b,c,d):
    a,b,c=sorted((a,b,c))
    key=(a,b,c,d,mp.mp.dps)
    if key in _MT3_cache: return _MT3_cache[key]
    if a==0:
        v=_MT3_boundary(b,c,d)
    else:
        v=MT3(a-1,b,c,d+1)+MT3(a,b-1,c,d+1)+MT3(a,b,c-1,d+1)
    _MT3_cache[key]=v; return v

_MT3b_cache={}
def _MT3_boundary(b,c,d):
    """MT₃(0,b,c;d)=Σ_{M≥2}S₂(b,c;M)·ζ_H(d,M+1),  b,c≥1, d≥2."""
    key=(b,c,d,mp.mp.dps)
    if key in _MT3b_cache: return _MT3b_cache[key]
    Mmax,Pmax=_default_Mmax_Pmax()
    e=_S2_pf_coeffs(b,c); imax=len(e)-1; bc=b+c
    # head
    H=[mpf(0)]*(imax+1)
    zhp1=mp.zeta(d,3)           # ζ_H(d,3) at M=2
    head=mpf(0)
    for M in range(2,Mmax+1):
        invMm1=mpf(1)/(M-1); p=mpf(1)
        for i in range(1,imax+1):
            p*=invMm1; H[i]+=p
        Mf=mpf(M)
        S2M=mpf(0)
        for i in range(1,imax+1):
            if e[i]: S2M+=e[i]*H[i]/Mf**(bc-i)
        head+=S2M*zhp1
        zhp1-=mpf(1)/mpf(M+1)**d
    # tail
    A=_asym_mul(_S2_asym(b,c,Pmax),_zetaHp1_asym(d,Pmax),Pmax)
    v=head+_tail_sum(A,Mmax)
    _MT3b_cache[key]=v; return v

# ---- (2,2)-sector convolution sum ----
_conv22_cache={}
def conv22(al,be,ga,de,q):
    """Σ_{K≥2} S₂(α,β;K)·S₂(γ,δ;K)/(2K)^q,  α,β,γ,δ≥1, q≥1."""
    al,be=sorted((al,be)); ga,de=sorted((ga,de))
    if (al,be)>(ga,de): al,be,ga,de=ga,de,al,be
    key=(al,be,ga,de,q,mp.mp.dps)
    if key in _conv22_cache: return _conv22_cache[key]
    Mmax,Pmax=_default_Mmax_Pmax()
    e1=_S2_pf_coeffs(al,be); e2=_S2_pf_coeffs(ga,de)
    imax=max(len(e1),len(e2))-1; s1=al+be; s2=ga+de
    # head
    H=[mpf(0)]*(imax+1); head=mpf(0); tq=mpf(2)**q
    for K in range(2,Mmax+1):
        invKm1=mpf(1)/(K-1); p=mpf(1)
        for i in range(1,imax+1):
            p*=invKm1; H[i]+=p
        Kf=mpf(K)
        S2a=mpf(0)
        for i in range(1,len(e1)):
            if e1[i]: S2a+=e1[i]*H[i]/Kf**(s1-i)
        S2b=mpf(0)
        for i in range(1,len(e2)):
            if e2[i]: S2b+=e2[i]*H[i]/Kf**(s2-i)
        head+=S2a*S2b/(tq*Kf**q)
    # tail
    A=_asym_mul(_S2_asym(al,be,Pmax),_S2_asym(ga,de,Pmax),Pmax)
    A={(l,p+q):c/tq for (l,p),c in A.items() if p+q<=Pmax}
    v=head+_tail_sum(A,Mmax)
    _conv22_cache[key]=v; return v

def _clear_caches():
    _MT3_cache.clear(); _MT3b_cache.clear(); _conv22_cache.clear()
    _tailz_cache.clear(); _MZV2.clear(); _TC.clear()

def W4(p,q):
    """W₄(p;q), p=(p₁..p₄), p_r≥1. Sign-sector decomposition."""
    from itertools import combinations
    idx=(0,1,2,3)
    total=mpf(0)
    # npos=1: 4 choices; lone positive m_i = Σ_{j≠i}n_j, n_j>0. Σ|m|=2Σn.
    for i in idx:
        neg=[j for j in idx if j!=i]
        a,b,c=[p[j] for j in neg]; d=p[i]
        total += MT3(a,b,c,d+q)/mpf(2)**q
    # npos=3: symmetric to npos=1 (m→-m). Another 4.
    for i in idx:
        neg=[j for j in idx if j!=i]
        a,b,c=[p[j] for j in neg]; d=p[i]
        total += MT3(a,b,c,d+q)/mpf(2)**q
    # npos=2: 6 choices; Σ_{pos}=Σ_{neg}=K.
    for pos in combinations(idx,2):
        neg=tuple(j for j in idx if j not in pos)
        al,be=p[pos[0]],p[pos[1]]; ga,de=p[neg[0]],p[neg[1]]
        total += conv22(al,be,ga,de,q)
    return total

# ---- validation ----
if __name__=='__main__':
    import time
    mp.mp.dps=60
    print("W₃((1,1,1);1) [expect π⁴/60]:")
    t0=time.time(); w=W3((1,1,1),1); print(f"  {mp.nstr(w,30)} vs {mp.nstr(mp.pi**4/60,30)} in {time.time()-t0:.2f}s")
    print("W₃((1,1,1);2) [expect 6ζ₅−3ζ₂ζ₃]:")
    w=W3((1,1,1),2); ref=6*zeta(5)-3*zeta(2)*zeta(3)
    print(f"  {mp.nstr(w,30)} vs {mp.nstr(ref,30)} diff={mp.nstr(w-ref,3)}")

    # ---- W4 precision-stability + brute-force cross-check ----
    print("\n== W4 validation ==")
    def MT3_brute(a,b,c,d,N):
        s=mpf(0)
        for n1 in range(1,N):
            for n2 in range(1,N):
                for n3 in range(1,N):
                    s+=mpf(1)/(mpf(n1)**a*n2**b*n3**c*(n1+n2+n3)**d)
        return s
    def conv22_brute(al,be,ga,de,q,N):
        s=mpf(0)
        for K in range(2,N):
            S1=sum(mpf(1)/(mpf(j)**al*(K-j)**be) for j in range(1,K))
            S2=sum(mpf(1)/(mpf(j)**ga*(K-j)**de) for j in range(1,K))
            s+=S1*S2/(mpf(2)*K)**q
        return s
    mp.mp.dps=40; _clear_caches()
    print("MT3(2,2,2;7) brute vs new:",
          mp.nstr(abs(MT3(2,2,2,7)-MT3_brute(2,2,2,7,60)),3),"(expect ≲1e-14)")
    print("conv22(2,2,2,2;5) brute vs new:",
          mp.nstr(abs(conv22(2,2,2,2,5)-conv22_brute(2,2,2,2,5,200)),3),"(expect ≲1e-14)")
    print("W4((1,1,1,1);1) vs 30ζ5−12ζ2ζ3:")
    _clear_caches(); w=W4((1,1,1,1),1); ref=30*zeta(5)-12*zeta(2)*zeta(3)
    print(f"  diff={mp.nstr(abs(w-ref),3)}")
    print("\nW4((2,2,2,2);1) precision stability:")
    vals={}
    for dps in (60,100,150):
        mp.mp.dps=dps; _clear_caches()
        t0=time.time(); vals[dps]=W4((2,2,2,2),1)
        print(f"  dps={dps:3d}: {mp.nstr(vals[dps],dps-5)}  ({time.time()-t0:.2f}s)")
    mp.mp.dps=160
    print(f"  |v60−v150| ={mp.nstr(abs(vals[60]-vals[150]),3)}")
    print(f"  |v100−v150|={mp.nstr(abs(vals[100]-vals[150]),3)}")

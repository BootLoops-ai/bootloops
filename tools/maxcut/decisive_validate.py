#!/usr/bin/env python3
r"""
maxcut_decisive_validate.py — 5 KNOWN-answer cases (time budget: any >5min → skip).

  g0-3loop     genus 0 → polylog       Track A  (3-loop 10-prop family, inline spec)
  ell-3loop    elliptic, j non-const   Track A  (3-loop 10-prop family, inline spec)
  pent5        rational (Appell, 1-d)  Track A  (5-point leg; needs external data,
                                                 SKIPs loudly without it)
  K3-ord3      K3 (order 3)            Track B  (period series K(x²)² in z=x²)
  sunpair-red  order 4 REDUCIBLE 2×2   Track B  (sunrise(z) + sunrise(1/z) pieces)

Calibration: topbox → genus 1 (a case a naive saddle heuristic collapses to
genus 0).

Also runs the sparse_cascade GF(p)-matmul boundary probe (gfp-matmul) as a
gated leg outside the 5-target count.

Exit status (script mode): 0 only if every executed leg PASSes; any FAIL —
including the calibration and the gfp-matmul probe — exits 1.  The pent5
SKIP without external data is the documented state, not a failure.
"""
import sys, os, time, json, signal
# All intra-family imports are package-relative; the external baikov_ls
# dependency of the pent5 leg is loaded explicitly at its use site below.
import sympy as sp
from math import factorial, comb
from .decisive import classify

s,t,m2=sp.symbols('s t m2',real=True); u=-s-t
PP4={('p1','p1'):0,('p2','p2'):0,('p3','p3'):0,
     ('p1','p2'):s/2,('p2','p3'):t/2,('p1','p3'):-(s+t)/2}
KIN={s:7,t:sp.Rational(-3,2),m2:1}
def cf(*tm):
    d={}
    for c,v in tm: d[v]=d.get(v,0)+c
    return d

def timed(fn, limit=300):
    def _h(sig,frm): raise TimeoutError
    signal.signal(signal.SIGALRM,_h); signal.alarm(limit)
    t0=time.time()
    try: r=fn(); return r, time.time()-t0, None
    except TimeoutError: return None, time.time()-t0, "TIMEOUT>5min"
    except Exception as e: return None, time.time()-t0, f"{type(e).__name__}:{e}"
    finally: signal.alarm(0)

RESULTS=[]
def run(name, expect, spec):
    r,w,err=timed(lambda: classify(spec))
    got=r['type'] if r else err
    ok=(expect in str(got)) if r else False
    print(f"[{name:12s}] expect={expect:22s} got={got!s:34s} "
          f"{'PASS' if ok else 'FAIL':4s} {w:6.1f}s"
          + (f"  j_distinct={r.get('j_distinct')}" if r and r.get('j_distinct') else ""))
    RESULTS.append(dict(name=name,expect=expect,got=str(got),ok=ok,
                        wall_s=round(w,2),detail=r))

# ---- 0. gfp-matmul: sparse_cascade envelope probe (gated, not in 5-count) --
from .sparse_cascade import matmul_boundary_probe
def _gfp_probe():
    t0=time.time()
    try: ok,checks=matmul_boundary_probe(); got='exact' if ok else \
        'FAIL:'+','.join(n for n,o in checks if not o)
    except Exception as e: ok,checks,got=False,[],f"{type(e).__name__}:{e}"
    w=time.time()-t0
    print(f"[{'gfp-matmul':12s}] expect={'exact-at-boundary':22s} got={got!s:34s} "
          f"{'PASS' if ok else 'FAIL':4s} {w:6.1f}s")
    return dict(name='gfp-matmul',expect='exact-at-boundary',got=got,ok=ok,
                wall_s=round(w,2),detail=[list(c) for c in checks])
PROBE=_gfp_probe()

# ---- calibration: topbox (a saddle heuristic collapses this to genus 0) -----
tb=[cf((1,'k1')),cf((1,'k1'),(-1,'p1')),cf((1,'k1'),(-1,'p1'),(-1,'p2')),
    cf((1,'k2')),cf((1,'k2'),(-1,'p1'),(-1,'p2'),(-1,'p3')),
    cf((1,'k1'),(-1,'k2')),cf((1,'k1'),(-1,'k2'),(1,'p3'))]
run("topbox-cal","elliptic",{'baikov':dict(
    gram=dict(loops=['k1','k2'],ext=['p1','p2','p3'],ext_gram=PP4,
              props=list(zip(tb,[0,0,0,m2,m2,m2,m2])),kin_syms=(s,t,m2)),
    kin=KIN)})

# ---- 1. g0-3loop ------------------------------------------------------------
mopp=[cf((1,'l'),(1,'p1')),cf((1,'l'),(1,'p1'),(-1,'k1')),
      cf((1,'l'),(1,'p1'),(-1,'k1'),(-1,'k2')),cf((1,'l'),(1,'p1'),(-1,'k2')),
      cf((1,'l'),(1,'p1'),(1,'p2'),(-1,'k2')),
      cf((1,'l'),(1,'p1'),(1,'p2'),(1,'p3'),(-1,'k2')),
      cf((1,'l'),(1,'p1'),(1,'p2'),(1,'p3')),cf((1,'l')),
      cf((1,'k1')),cf((1,'k2'))]
run("g0-3loop","polylog",{'baikov':dict(
    gram=dict(loops=['l','k1','k2'],ext=['p1','p2','p3'],ext_gram=PP4,
              props=list(zip(mopp,[m2]*8+[0,0])),kin_syms=(s,t,m2)),
    kin=KIN)})

# ---- 2. ell-3loop -----------------------------------------------------------
lq=[cf((1,'k1')),cf((1,'k1'),(-1,'k2')),cf((1,'k1'),(-1,'k2'),(1,'p2')),
    cf((1,'k1'),(-1,'k2'),(-1,'k3'),(1,'p2')),
    cf((1,'k1'),(-1,'k2'),(-1,'k3'),(1,'p2'),(1,'p3')),
    cf((1,'k1'),(-1,'k3'),(1,'p2'),(1,'p3')),
    cf((1,'k1'),(-1,'k3'),(-1,'p1')),cf((1,'k1'),(-1,'p1')),
    cf((1,'k2')),cf((1,'k3'))]
run("ell-3loop","elliptic",{'baikov':dict(
    gram=dict(loops=['k1','k2','k3'],ext=['p1','p2','p3'],ext_gram=PP4,
              props=list(zip(lq,[m2]*8+[0,0])),
              prefer_free=['k2k3','k2p1','k2p2','k3p2','k3p3'],
              kin_syms=(s,t,m2)),
    kin=KIN)})

# ---- 3. pent5 (pre-built B0 at P0; external data) --------------------------
# The external dir is APPENDED to sys.path (never inserted at [0], which would
# shadow top-level names for the rest of the process), and baikov_ls is loaded
# EXPLICITLY by file path under a private name so sys.modules['baikov_ls'] is
# never aliased.
_R3C=os.environ.get('MAXCUT_R3C','')  # dir providing baikov_ls.py (not shipped)
if _R3C and os.path.isdir(_R3C):
    if _R3C not in sys.path: sys.path.append(_R3C)
    import importlib.util as _ilu
    _bl_spec=_ilu.spec_from_file_location('_maxcut_validate_baikov_ls',
                                          os.path.join(_R3C,'baikov_ls.py'))
    _baikov_ls=_ilu.module_from_spec(_bl_spec); _bl_spec.loader.exec_module(_baikov_ls)
    baikov_B0, POINTS = _baikov_ls.baikov_B0, _baikov_ls.POINTS
    B0,(z9,z10,z11),_=baikov_B0(POINTS['P0'],verbose=False)
    run("pent5","polylog",{'baikov':dict(B=B0,isp=[z9,z10,z11],kin={})})
else:
    print('SKIP pent5: external data dir not present '
          '(set MAXCUT_R3C to a dir containing baikov_ls.py)', flush=True)

# ---- 4. K3-ord3 (period series K(x²)² in z=x²) -----------------------------
def K2_series(N):
    c=[sp.Rational(comb(2*n,n)**2,16**n) for n in range(N+1)]
    return [sum(c[i]*c[n-i] for i in range(n+1)) for n in range(N+1)]
run("K3-ord3","K3",{'series':K2_series(70)})

# ---- 5. sunpair-red (holomorphic pieces = sunrise(z) + sunrise(1/z)) -------
def sun_series(N):
    a=[]
    for n in range(N+1):
        nf=factorial(n); tot=0
        for i in range(n+1):
            for j in range(n+1-i):
                k=n-i-j; m=nf//(factorial(i)*factorial(j)*factorial(k)); tot+=m*m
        a.append(tot)
    return a
# reciprocal-argument piece expanded at the SAME point z→0:
#   Φ(1/x) hol series in z=1/t at t→∞ is sun_series; at z→0 the second
#   Frobenius piece has a_n^rec = a_n·(shifted).  For the reducibility test we
#   supply BOTH holomorphic pieces at their own MUM points; distinct order-2
#   ops with distinct singular sets ⇒ LCLM order 4.
def sun_recip_series(N):
    # hol period of L_sun(1/x) at x→0: indicials {0,0} ⇒ series with same
    # recursion but from R_j(t)=(t-1)(9t-1)θ²+…; a_0=1 fixed by R_j(0)θ²=θ².
    # Solve the 3-term recursion (9n²)a_n = (10n²-10n+3)a_{n-1} - (n-1)²a_{n-2}.
    a=[1,sp.Rational(1,3)]
    for n in range(2,N+1):
        a.append(((10*n*n-10*n+3)*a[-1]-(n-1)**2*a[-2])/(9*n*n))
    return a
run("sunpair-red","2×elliptic",{'series':sun_series(80),
                         'series_aux':sun_recip_series(80)})

# ---- summary ---------------------------------------------------------------
np=sum(1 for r in RESULTS[1:] if r['ok'])   # exclude calibration from 5-count
print(f"\nSUMMARY: {np}/5 validation targets PASS "
      f"(calibration topbox: {'PASS' if RESULTS[0]['ok'] else 'FAIL'}; "
      f"gfp-matmul probe: {'PASS' if PROBE['ok'] else 'FAIL'})")
# receipt: written beside the package module — the canonical copy; each run
# rewrites it.
json.dump([PROBE]+RESULTS,
          open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
          'decisive_validate.json'),'w'),indent=2,default=str)
# battery honesty: a FAIL on any executed leg must be visible in the exit
# status (a SKIPped leg is absent from RESULTS and does not gate).
if __name__ == '__main__':
    sys.exit(0 if (PROBE['ok'] and all(r['ok'] for r in RESULTS)) else 1)

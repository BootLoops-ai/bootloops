#!/usr/bin/env python3
r"""
maxcut_decisive.py — maximal-cut geometry classifier (χ / GB-staircase).

A naive saddle heuristic (differentiating the carrier ISP)
mis-classifies the topbox as genus 0; this path uses
decidable algebraic invariants:

  Track A (Baikov B(z;kin) → curve invariants)
    1. per-ISP degree profile of B on the cut (exact ℚ).
    2. carrier-slice: highest-deg ISP kept symbolic, spectators sliced to
       generic rationals; squarefree degree → hyperelliptic genus.
    3. disc-chain: iterated discriminant on deg-2 ISPs only (never the
       carrier).  Drains-completely ⇒ rational; residual ≥2 vars ⇒ higher.
    4. genus-1 → j-invariant on ≥3 slices (constant ⇒ isotrivial).

  Track B (period series → PF order, Krylov mod ℚ)
    1. minimal-order θ-annihilator of the holomorphic-period series a_n.
    2. reducibility: order-⌊m/2⌋ right-factor annihilates a_n ⇒ reducible
       (the ice-cone LCLM signature).
    Cross-check via mod-p Griffiths–Dwork: tools/pf_rank.jl (LP-polynomial
    input, GB-staircase χ + rank(∂_s^k mod exact) over GF(p)).

Verdict map:  order 1 / genus 0 → polylog;  2 → elliptic;  3 → K3;
              4 irreducible → CY₃;  4 reducible → 2×elliptic (ice-cone).

Input spec (dict):
    {'baikov': {'B': expr|None, 'gram': {loops,ext,ext_gram,props,
                 prefer_free}, 'isp': [...], 'kin': {sym:val}}}
    {'series': list|callable, 'series_aux': list (2nd Frobenius, optional)}
    {'lp': {'G': str, 'kin': [...], 'x': [...], 'svar': str}}   # → pf_rank.jl
"""
import sympy as sp, subprocess, json, os, time
from .decisive_helpers import (gram_baikov, carrier_slice_genus,
    disc_chain, minimal_pf_order, find_theta_op)

HERE=os.path.dirname(os.path.abspath(__file__))
# pf_rank.jl lives in the FLAT tools/ dir (the package's parent), shared
# with annihilator.py — do NOT move it (a {HERE}-relative include would
# silently soft-fail Track B-jl to an error dict).
PF_RANK_JL=os.path.join(os.path.dirname(HERE),'pf_rank.jl')
_TYPE={0:'polylog',1:'elliptic',2:'hyperelliptic'}

def classify_baikov(spec):
    """Track A."""
    t0=time.time()
    B=spec.get('B')
    if B is None:
        g=spec['gram']
        B,isp=gram_baikov(g['loops'],g['ext'],g['ext_gram'],g['props'],
                          prefer_free=g.get('prefer_free'),
                          kin_syms=g.get('kin_syms',()))
    else:
        isp=list(spec['isp'])
    kin=spec.get('kin',{})
    cs=carrier_slice_genus(B,isp,kin)
    _,dc=disc_chain(B,isp,kin)
    # verdict.  Disc-chain (single- + multi-factor) is the primary decidable
    # invariant: it never differentiates a carrier (the naive saddle heuristic's bug) and
    # cannot spuriously LOWER genus.  Drains to ≤1 residual var ⇒ curve
    # (sqfree deg → genus).  Stalls at ≥2 vars ⇒ higher: if a genuine
    # deg-≥3 carrier is present read genus from carrier-slice; else K3/CY
    # candidate → Track B.
    caveats=[]
    nres=len(dc['residual_vars'])
    if nres==0 or dc.get('sqfree')==0:
        pf,typ,g=1,'polylog',0
    elif nres==1:
        sq=dc['sqfree']
        g=0 if sq<=2 else 1 if sq<=4 else (sq-1)//2
        pf=max(1,g+1); typ=_TYPE.get(g,f'genus-{g}')
    elif max(cs['degs'].values(),default=0)>=3:
        g=cs['genus']; pf=max(1,g+1); typ=_TYPE.get(g,f'genus-{g}')
    else:
        g=None; pf=None; typ='higher'
        caveats.append(f"multi-factor disc-chain stalls at {nres} vars "
                       f"({dc.get('n_factors')} factors, max per-var deg "
                       f"{dc.get('max_factor_pervar_deg')}) ⇒ K3/CY "
                       "candidate; supply LP-G or series for Track B")
    return dict(track='A', type=typ, genus=g, pf_order=pf,
                degs=cs['degs'], carrier=cs['carrier'],
                j_distinct=cs.get('j_distinct'),
                disc_chain=dc, method='baikov-carrier+disc',
                caveats=caveats, wall_s=round(time.time()-t0,2))

def classify_series(series, series_aux=None, max_order=6, degz_max=8):
    """Track B."""
    t0=time.time()
    a=series if isinstance(series,(list,tuple)) else series(80)
    m,op=minimal_pf_order(a,max_order,degz_max)
    m1=m; reducible=None; parts=[m]
    if series_aux is not None:
        a2=series_aux if isinstance(series_aux,(list,tuple)) else series_aux(80)
        m2,op2=minimal_pf_order(a2,max_order,degz_max)
        parts=[m,m2]
        # LCLM order = m+m2 unless the two ops coincide (same right-factor):
        same=(m==m2 and op and op2 and
              all(sp.expand(op[0][j]*op2[0][m]-op2[0][j]*op[0][m])==0
                  for j in range(m+1)))
        m=(m if same else (m or 0)+(m2 or 0))
        reducible=(not same) and m>max(m1 or 0,m2 or 0)
    elif m and m>=3:
        reducible=find_theta_op(a,m//2,degz_max=degz_max+2) is not None
    typ={1:'polylog',2:'elliptic',3:'K3',4:'CY3'}.get(m,f'order-{m}')
    if m==4 and reducible: typ='2×elliptic (reducible)'
    if m==3 and reducible: typ='K3 (Sym² check needed)'
    lead=str(sp.factor(op[0][m1])) if op and m1 in op[0] else None
    return dict(track='B', pf_order=m, pf_parts=parts, reducible=reducible,
                type=typ, leading_symbol=lead,
                method='series-krylov', wall_s=round(time.time()-t0,2))

def classify_lp(G, kin, x, svar, Lmax=6, timeout=300):
    """Track B via pf_rank.jl (mod-p Griffiths–Dwork + χ_regulated GB-staircase)."""
    t0=time.time()
    jl=(f'include("{PF_RANK_JL}");using JSON;'
        f'r=pf_probe("{G}",{json.dumps(kin)},{json.dumps(x)},"{svar}";Lmax={Lmax});'
        'print(JSON.json(Dict(pairs(r))))')
    out=None   # bound before try: TimeoutExpired/missing-julia raise BEFORE
               # assignment, and the handler must not UnboundLocalError
    try:
        jcmd=["julia","+1.10"]
        if os.environ.get("PLD_ENV"):  # a julia env with the pf_rank deps
            jcmd.append("--project="+os.environ["PLD_ENV"])
        out=subprocess.run(jcmd+["-e",jl],capture_output=True,text=True,timeout=timeout)
        r=json.loads(out.stdout.strip().splitlines()[-1])
    except Exception as e:
        return dict(track='B-jl', error=f"{type(e).__name__}: {e}",
                    stderr=(out.stderr if out is not None and out.stderr else '')[:400])
    m=r['pf_order']; hr=r['holonomic_rank']
    typ={1:'polylog',2:'elliptic',3:'K3',4:'CY3'}.get(m,f'order-{m}')
    if r.get('reducible') and m>=3: typ+=' (reducible sub-D-module)'
    return dict(track='B-jl', pf_order=m, holonomic_rank=hr,
                reducible=r.get('reducible'), type=typ,
                method='pf_rank.jl GD mod-p + χ_regulated',
                wall_s=round(time.time()-t0,2))

def classify(spec):
    if 'baikov' in spec:  return classify_baikov(spec['baikov'])
    if 'series' in spec:  return classify_series(spec['series'],
                                                 spec.get('series_aux'))
    if 'lp' in spec:      return classify_lp(**spec['lp'])
    raise ValueError("spec needs one of: baikov | series | lp")

if __name__=='__main__':
    import sys
    print(json.dumps(classify(json.load(open(sys.argv[1]))),indent=2,default=str))

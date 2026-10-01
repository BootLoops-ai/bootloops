#!/usr/bin/env python3
"""check_degd.py — estimate max deg_d of A(d,eta) from etanumd slice outputs.
Method: for each slice at d=p/q, max coefficient digit-length L satisfies
  L ≈ deg_d · log10(max(p,q)) + C   (C = eta-independent integer prefactor size)
Linear fit L vs log10(max(p,q)) over ≥2 slices → slope = deg_d, intercept = C.
This is the DIRECT route-decision measurement (Route A viable iff deg_d ≲ 24).
"""
import sys, re, glob, json, math
import os
here = os.environ.get('ETANUMD_DIR', '.')  # dir with the etanumd slice receipts
fam = os.environ.get('ETANUMD_FAM', '*')    # kira family name (glob-ok)
pts = []
for d in sorted(glob.glob(here + '/M89_eps_*/results/' + fam + '/kira_target.m')):
    slj = json.load(open(d.rsplit('/results/',1)[0] + '/ETANUMD_SLICE.json'))
    p, q = slj['d_num'], slj['d_den']
    mx = 0
    with open(d) as f:
        for ln in f:
            for m in re.finditer(r'\b\d{10,}\b', ln):
                if len(m.group(0)) > mx: mx = len(m.group(0))
    pts.append((max(abs(p), q), mx, slj['eps']))
    print(f"  eps={slj['eps']:>10s}  d={p}/{q}  max_digits={mx}")
if len(pts) < 2:
    print("need ≥2 finished slices"); sys.exit(1)
xs = [math.log10(m) for m,_,_ in pts]; ys = [L for _,L,_ in pts]
n=len(xs); sx=sum(xs); sy=sum(ys); sxx=sum(x*x for x in xs); sxy=sum(x*y for x,y in zip(xs,ys))
slope = (n*sxy - sx*sy)/(n*sxx - sx*sx); intercept = (sy - slope*sx)/n
print(f"\nlinear fit: max_digits ≈ {slope:.1f} · log10(max(p,q)) + {intercept:.1f}")
print(f"=> deg_d ≈ {slope:.0f}  (Route A needs ~{int(slope)+4} fit + 3 verify = {int(slope)+7} slices)")
route = 'A' if slope < 24 else ('A-extended' if slope < 50 else 'B')
print(f"=> ROUTE: {route}")
json.dump({'deg_d_est': slope, 'intercept': intercept, 'pts': [(e,m,L) for m,L,e in pts],
           'route': route, 'slices_needed': int(slope)+7}, open(here+'/DEGD.json','w'), indent=2)

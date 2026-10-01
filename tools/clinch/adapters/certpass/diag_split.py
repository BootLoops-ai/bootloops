import sys, json, numpy as np
sys.path.insert(0,'.')
import __init__ as pkg
from oracle_etas import EtasOracle, Part, ser, deser, _pairs_chunk, _sources_chunk, _init_worker
from flint import arb, ctx
tag = 'q3_hi'
z = dict(np.load(f'data/etas_{tag}.npz'))
c = json.load(open(f'receipts/etas_{tag}_corner/CERT.json'))
thc = np.asarray(c['polished_center'], float)
part = Part(blocks=[list(range(1,9))], block_names=['triggering'], border_idx=[0])
o = EtasOracle(z, part, 'd', nworkers=1)
ctx.prec = 192
# build theta balls with the corner-rung radii (from diag: scales*1e-9; ltau exact)
scales = [0.001275, 0.002993, 0.007901, 0.000675, 0.001208, 0.002875, 0.00276, 0.006085]
# order in diag: [lk0,a,lc,om,ld,gam,rho] + border lmu
rad9 = [scales[7]*1e-9, scales[0]*1e-9, scales[1]*1e-9, scales[2]*1e-9,
        scales[3]*1e-9, 0.0, scales[4]*1e-9, scales[5]*1e-9, scales[6]*1e-9]
th = [arb(float(thc[k]), rad9[k]) for k in range(9)]
th_ser = [ser(v) for v in th]
_init_worker(o._shared())
# sources H on a SLICE of sources (first 3000) + pairs H on first 3000 targets
outs = _sources_chunk((0, 3000, th_ser, 192, 2))
idx = 10
hs = {}
for k in range(9):
    for l in range(k,9):
        hs[(k,l)] = deser(outs[idx]); idx += 1
print('sources[0:3000] H(om,om): mid', float(hs[(4,4)].mid()), 'rad', float(hs[(4,4)].rad()))
print('sources[0:3000] H(lc,om): mid', float(hs[(3,4)].mid()), 'rad', float(hs[(3,4)].rad()))
outp = _pairs_chunk((0, 3000, th_ser, 192, 2))
idx = 10
hp = {}
for k in range(9):
    for l in range(k,9):
        hp[(k,l)] = deser(outp[idx]); idx += 1
print('pairs[t<3000] H(om,om): mid', float(hp[(4,4)].mid()), 'rad', float(hp[(4,4)].rad()))
print('pairs[t<3000] H(lk0,lk0): mid', float(hp[(1,1)].mid()), 'rad', float(hp[(1,1)].rad()))
o.close()

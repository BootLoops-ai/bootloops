import sys, numpy as np
sys.path.insert(0,'.')
import __init__ as pkg
from oracle_etas import EtasOracle, Part
from twin_etas import TwinEtas
from clinch.mask import MaskedOracle
from clinch import engine
from baller.hygiene import ctx_guard
import baller.certify as bc
bk = bc.block_krawczyk
tag = sys.argv[1]
z = dict(np.load(f'data/etas_{tag}.npz'))
import json
c = json.load(open(f'receipts/etas_{tag}_corner/CERT.json'))
thc = np.asarray(c['polished_center'], float)
part = Part(blocks=[list(range(1,9))], block_names=['triggering'], border_idx=[0])
o = EtasOracle(z, part, 'diag', nworkers=16)
m = MaskedOracle(o, {5: 12.26})
with ctx_guard(prec=192):
    box0 = bk.box_around(([[float(thc[t]) for t in m.part.blocks[0]]],[float(thc[0])]), 0.0)
    Ht = m.H(box0)
    sz, sg = engine._radius_scales(bk, Ht)
    ru = 1e-9
    r = ([[ru*float(v) for v in sz[0]]],[ru*float(v) for v in sg])
    box = bk.box_around(([[float(thc[t]) for t in m.part.blocks[0]]],[float(thc[0])]), r)
    D,B,G = m.H(box)
n = D[0].nrows()
names = ['lk0','a','lc','om','ld','gam','rho']
Dm = np.array([[float(D[0][i,j].mid()) for j in range(n)] for i in range(n)])
Dr = np.array([[float(D[0][i,j].rad()) for j in range(n)] for i in range(n)])
print('scales:', [round(float(v),6) for v in sz[0]], round(float(sg[0]),6), flush=True)
print('D mid diag:', np.diag(Dm).round(1))
print('D rad matrix max per row:', Dr.max(1).round(6))
i,j = np.unravel_index(Dr.argmax(), Dr.shape)
print('max rad', Dr.max(), 'at', names[i], names[j], 'mid', Dm[i,j])
Ym = np.linalg.inv(Dm)
print('|Y| max:', round(np.abs(Ym).max(),1), 'cond:', round(np.linalg.cond(Dm),1))
print('|Y|@Dr max:', (np.abs(Ym)@Dr).max())
# tight-center rads for comparison
D0m = np.array([[float(Ht[0][0][i,j].rad()) for j in range(n)] for i in range(n)])
print('tight-center rad max:', D0m.max())
o.close()

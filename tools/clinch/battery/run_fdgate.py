import sys, time, json
import os as _os
_TOOLS = _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, _os.path.join(_TOOLS, "clinch"))
sys.path.insert(0, _os.path.join(_TOOLS, "baller"))
import numpy as np
import clinch; clinch.verify(quiet=True)
from clinch import adapter_v31 as A
from clinch.oracle_v31 import ModelV31Oracle
from flint import arb
from baller.hygiene import ctx_guard

tag = sys.argv[1] if len(sys.argv)>1 else 'fit13'
arrs, consts, theta, prov = A.load(tag)
o = ModelV31Oracle(arrs, consts)
t0 = time.time()
with ctx_guard(prec=256):
    th = [arb(float(t)) for t in theta]
    _, _, acc = o.evaluate(th, order2=True)
    print('H(tight) pass wall:', round(time.time()-t0,1), 's', flush=True)
    t1 = time.time()
    D, B, G = acc.to_mats()
    print('to_mats wall:', round(time.time()-t1,1), 's', flush=True)
cols = [0, 34, 43, 44, 1710, 2345]
g = A.gate_fd_hessian(o, theta, cols, prec=256, hess_mats=(D, B, G))
print(json.dumps(g, indent=1, default=float))
json.dump(g, open(f'GATE_FDHESS_{tag}.json', 'w'), indent=1, default=float)

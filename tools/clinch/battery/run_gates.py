import sys, time, json
import os as _os
_TOOLS = _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, _os.path.join(_TOOLS, "clinch"))
sys.path.insert(0, _os.path.join(_TOOLS, "baller"))
import numpy as np
import clinch; clinch.verify(quiet=True)
from clinch import adapter_v31 as A
from clinch.oracle_v31 import ModelV31Oracle

tag = sys.argv[1] if len(sys.argv) > 1 else 'fit13'
arrs, consts, theta, prov = A.load(tag)
o = ModelV31Oracle(arrs, consts)
g = A.gate_gradient_adjudicated(o, theta, tag)
print(json.dumps(g, indent=1, default=float))
json.dump(g, open(f'GATE_GRADIENT_{tag}.json','w'), indent=1, default=float)

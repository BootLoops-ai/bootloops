import sys, time, json
import os as _os
_TOOLS = _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, _os.path.join(_TOOLS, "clinch"))
sys.path.insert(0, _os.path.join(_TOOLS, "baller"))
import numpy as np
import clinch; clinch.verify(quiet=True)
from clinch import adapter_v31 as A, engine, cert
from clinch.oracle_v31 import ModelV31Oracle

tag = sys.argv[1]
arrs, consts, theta, prov = A.load(tag)
o = ModelV31Oracle(arrs, consts)
gates = {}
for nm, key in ((f'GATE_GRADIENT_{tag}.json', 'gradient'), (f'GATE_FDHESS_{tag}.json', 'fdhess')):
    try: gates[key] = json.load(open(nm))
    except FileNotFoundError: pass

print('=== direct certificate (as-receipted theta)', flush=True)
res = engine.certify(o, theta, log=lambda s: print('  '+s, flush=True))
print('direct:', res['verdict'], res.get('failed'), flush=True)

print('=== newton polish + polished-center certificate', flush=True)
thp, prcpt = engine.newton_polish(o, theta, log=lambda s: print('  '+s, flush=True))
resp = engine.certify(o, thp, log=lambda s: print('  '+s, flush=True))
print('polished:', resp['verdict'], resp.get('failed'), 'radius', resp.get('radius'), flush=True)

out_dir = _os.path.join(A.REFERENCE_DIR, f'clinch_{tag}_null')
c = cert.write_cert(res, out_dir, f'{tag}_null (as-receipted)', gates=gates,
                    extra=dict(assembly_provenance={k: v for k, v in prov.items() if not k.startswith('_')},
                               polished_center_certificate=dict(
                                   verdict=resp['verdict'],
                                   radius=resp.get('radius'), radius_min=resp.get('radius_min'),
                                   prec_dps=resp.get('prec_dps'),
                                   theta_polished_sha16=resp.get('theta_sha16'),
                                   distance_inf_from_receipted=prcpt['distance_inf'],
                                   polish_receipt=prcpt,
                                   contraction_margins=(resp.get('krawczyk') or {}).get('margins_named'),
                                   max_margin=(resp.get('krawczyk') or {}).get('max_margin'),
                                   pd_margins=(resp.get('pd') or {}).get('margins_named') or (resp.get('pd') or {}).get('margins'),
                                   statement=cert.human_statement(resp, f'{tag}_null POLISHED-CENTER'))))
np.savez(out_dir + '/THETA_POLISHED.npz', theta=thp)
print(json.dumps(dict(direct=res['verdict'], polished=resp['verdict'],
                      pol_radius=resp.get('radius'), pol_radius_min=resp.get('radius_min'),
                      dist=prcpt['distance_inf'], wall_direct=res.get('wall_s_total'),
                      wall_pol=resp.get('wall_s_total')), indent=1, default=float))
if resp['verdict'].startswith('CERTIFIED'):
    print('POLISHED STATEMENT:', cert.human_statement(resp, f'{tag}_null POLISHED-CENTER')[:400])

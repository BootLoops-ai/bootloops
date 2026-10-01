"""Generate examples/jc_quartet.json from Jukes-Cantor quartet kernels.

Imports KERNEL (the 15 multilinear 2^5 kernels, = 256*p_k) from the
Jukes-Cantor kit shipped beside this script in examples/jc_kit/ (two modules
from the jackandjill phylogenetics package, the sibling repository jackandjill;
PHYLO_GKZ overrides the kit directory; an optional argv[1] overrides the
output path).  NOTE: with these kernels
Z_spec(u) = 256^N * Z_true(u), N = sum u — divide by 256^N for the
biological value.
"""
import os, sys, json, itertools
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
# kernel source: the JC quartet kit beside this script (needs numpy + sympy)
GKZ = os.environ.get('PHYLO_GKZ') or os.path.join(HERE, 'jc_kit')
sys.path.insert(0, GKZ)
from b3_recurrence import KERNEL  # noqa: E402

EPS = list(itertools.product((0, 1), repeat=5))

polys = {}
for lab, K in sorted(KERNEL.items()):
    pd = {}
    for e in EPS:
        c = int(np.asarray(K)[e])
        if c:
            pd[' '.join(map(str, e))] = c
    polys[lab] = pd

spec = {
    'name': 'jc_quartet',
    'comment': ('JC69 quartet site-pattern kernels (256*p_k), n=5 branch '
                'lengths; Z_spec = 256^N * Z_true. target_u below is the '
                'banked N=20 synthetic (b4_bayes_factors); replace with '
                'real pattern counts. order auto = ascending counts, '
                'largest class (xxxx) last as the pure final ray.'),
    'variables': ['x1', 'x2', 'x3', 'x4', 'x5'],
    'polynomials': polys,
    'target_u': {'xxxx': 10, 'xxyy': 4, 'xyxy': 2, 'xyyx': 2,
                 'xxxy': 1, 'xyzw': 1},
    'order': 'auto',
    'options': {},
}
out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, 'jc_quartet.json')
json.dump(spec, open(out, 'w'), indent=1)
print(f"wrote {out}: {len(polys)} kernels, "
      f"N={sum(spec['target_u'].values())}")

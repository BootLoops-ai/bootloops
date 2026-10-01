# eichler.genus2 — parametrized complex-basepoint transport:
#   python3 exact_transport_driver.py <re> <im> <tag> [dps] [column]
# PUBLIC-REPO NOTE: this is the exact-rational-basepoint transport PATTERN (the
# dyadic basepoint — a decimal basepoint poisons PSLQ). It imports `monodromy_transport.L5System`,
# the originating project's DE-system spec (not shipped) — swap in your
# own system class exposing the same interface. Wayfinder ships in
# this repo's tools/ (three levels up from tools/eichler/genus2/) and is found
# via the relative path insert below.
import sys, json, time, os
_TOOLS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, _TOOLS)
from mpmath import mp, mpf, mpc
from wayfinder.transport import transport_fixed_eps
from monodromy_transport import L5System

re_n, re_d = sys.argv[1].split('/')
im_n, im_d = sys.argv[2].split('/')
tag = sys.argv[3]
dps = int(sys.argv[4]) if len(sys.argv) > 4 else 120
xb = mpf(1)/100      # dyadic (frame lives here — do NOT change)
with mp.workdps(dps + 30):
    xc = mpc(mpf(int(re_n))/int(re_d), mpf(int(im_n))/int(im_d))
sysm = L5System()
cols = {}
t0 = time.time()
onecol = int(sys.argv[5]) if len(sys.argv) > 5 else None
rng = [onecol] if onecol is not None else range(5)
for j in rng:
    y0 = [mpf(1) if i == j else mpf(0) for i in range(5)]
    y1 = transport_fixed_eps(sysm, 0, xb, xc, y0, dps)
    cols[j] = [mp.nstr(v, dps) for v in y1]
    print(f'{tag} col {j}: {time.time()-t0:.1f}s', flush=True)
suffix = f'_col{onecol}' if onecol is not None else ''
json.dump({'dps': dps, 'xc': [sys.argv[1], sys.argv[2]], 'cols': {str(k): v for k, v in cols.items()}},
          open(f'xcf_{tag}{suffix}.json', 'w'))
print(f'wrote xcf_{tag}{suffix}.json')

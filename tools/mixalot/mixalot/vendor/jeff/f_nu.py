"""Thermodynamic limit of the exact 2-mixture evidence: f(nu) = lim_k [lnZ + N ln k]/k
(baseline N ln k = frozen/uniform likelihood scale), three data shapes at fixed filling nu."""
import sys, time
from math import log
from bigk2 import Z_fast2

def lnint(x):
    sh = max(0, x.bit_length()-500); return log(x >> sh) + sh*0.6931471805599453
def lnZ(U):
    Z = Z_fast2(list(U)); return lnint(Z.numerator) - lnint(Z.denominator)

def shape(name, k, nu):
    if name == 'uniform':  return [nu]*k
    if name == 'twoblock': return [2*nu]*(k//2) + [0]*(k//2)          # half the categories empty
    if name == 'tworate':  return [(3*nu)//2]*(k//2) + [nu - (3*nu)//2 + nu]*0 + [nu//2]*(k//2)  # 1.5nu / 0.5nu
    raise ValueError

KS = [25, 50, 100, 200, 400, 800]
print(f"{'shape':>9} {'nu':>4} | " + " ".join(f"k={k:<5}" for k in KS) + "   (Delta_k = (lnZ + N ln k)/k)")
for name in ['uniform', 'twoblock', 'tworate']:
    for nu in [2, 10, 50]:
        row = []
        for k in (KS if name=='uniform' else [x for x in KS if x%2==0]):
            U = shape(name, k, nu)
            N = sum(U)
            assert N == nu*k, (name, k, nu, N)
            v = (lnZ(U) + N*log(k))/k
            row.append(v)
        print(f"{name:>9} {nu:>4} | " + " ".join(f"{v:8.4f}" for v in row) + f"   conv: {row[-1]-row[-2]:+.4f}")

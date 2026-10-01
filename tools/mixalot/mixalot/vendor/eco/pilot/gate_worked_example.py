"""Gate vs Etienne's own worked example D=(1,1,2,3,5,8) (SA5/untb lineage).

Targets (from Etienne 2005 and the SA5/untb lineage, confirmed against the literature):
  K(D,6)=K(D,20)=1, K(D,7)=1507/210, K(D,13)=27751573/60480, K(D,14)=1127899/2520
  ln K(D,7) = 1.9707686679090268
  ln L(theta=7.10, m=0.225)            = -4.796695178274603
  ln L(theta=7.047958, m=0.22635923)   = -4.796686701878364  (float-noise +/- 4e-15)
  P[D] at optimum                      = 8.25706e-3
"""
import sys
from fractions import Fraction
from flint import arb, ctx

from ball_engine import K_exact, K_ball, logP_ball

D = [1, 1, 2, 3, 5, 8]
fails = 0

S, K = K_exact(D)
targets = {6: Fraction(1), 7: Fraction(1507, 210), 13: Fraction(27751573, 60480),
           14: Fraction(1127899, 2520), 20: Fraction(1)}
for A, t in targets.items():
    got = K[A - S]
    ok = got == t
    fails += not ok
    print(f"K(D,{A}) = {got} {'== ' + str(t) + ' OK' if ok else 'FAIL expected ' + str(t)}")

prec = 128
Kp = K_ball(D, prec)
for (th, m), tgt in [(("7.10", "0.225"), -4.796695178274603),
                     (("7.047958", "0.22635923"), -4.796686701878364)]:
    lp = logP_ball(D, th, m, prec, Kpoly=Kp)
    diff = abs(float(lp.mid()) - tgt)
    ok = diff < 5e-13
    fails += not ok
    print(f"lnL({th},{m}) = {lp} vs {tgt}  |diff|={diff:.2e} {'OK' if ok else 'FAIL'}")

lp = logP_ball(D, "7.047958", "0.22635923", prec, Kpoly=Kp)
ctx.prec = prec
P = lp.exp()
ok = abs(float(P.mid()) - 8.25706e-3) < 1e-8
fails += not ok
print(f"P[D] at optimum = {P} vs 8.25706e-3 {'OK' if ok else 'FAIL'}")

import math
lnK7 = math.log(1507 / 210)
ok = abs(lnK7 - 1.9707686679090268) < 1e-15
fails += not ok
print(f"ln K(D,7) = {lnK7!r} vs 1.9707686679090268 {'OK' if ok else 'FAIL'}")

sys.exit(1 if fails else 0)

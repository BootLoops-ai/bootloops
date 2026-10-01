#!/usr/bin/env python3
"""Positive-control oracle self-check.

Formula under test (transcribed from Gimenez Umbert--Sturmfels 2501.10805,
Simplex Amplitudes section, and GS split-kinematics section for X(3,6)):

  J(u0,u1,u2) := int_0^inf int_0^inf (dx/x)(dy/y) x^{u1} y^{u2} (1+x+y)^{-(u0+u1+u2)}
              =  Gamma(u0) Gamma(u1) Gamma(u2) / Gamma(u0+u1+u2),   u0,u1,u2 > 0.

At X(3,6) split kinematics the full 20-minor integrand collapses exactly to
J(u_{0,1},u_{1,1},u_{2,1}) * J(u_{0,2},u_{1,2},u_{2,2}).
Verify J's Gamma-form against direct 2-dim tanh-sinh quadrature at TWO
parameter choices before using it as the pilot oracle.
"""
import mpmath as mp
import time, json

mp.mp.dps = 30

def gamma_form(u0, u1, u2):
    return mp.gamma(u0)*mp.gamma(u1)*mp.gamma(u2)/mp.gamma(u0+u1+u2)

def direct(u0, u1, u2):
    # substitute x=t/(1-t), y=s/(1-s):  1+x+y = (1-ts)/((1-t)(1-s))
    # integrand -> t^{u1-1} s^{u2-1} (1-t)^{u0+u2-1} (1-s)^{u0+u1-1} (1-ts)^{-T}
    T = u0+u1+u2
    with mp.workdps(mp.mp.dps + 15):
        f = lambda t, s: (t**(u1-1) * s**(u2-1) * (1-t)**(u0+u2-1)
                          * (1-s)**(u0+u1-1) * (1-t*s)**(-T))
        v = mp.quad(lambda t: mp.quad(lambda s: f(t, s), [0, 1]), [0, 1])
    return +v

cases = [(mp.mpf(3)/4, mp.mpf(1)/3, mp.mpf(1)/2),
         (mp.mpf(2)/3, mp.mpf(2)/5, mp.mpf(9)/7)]
out = []
for (u0, u1, u2) in cases:
    t0 = time.time()
    g = gamma_form(u0, u1, u2)
    d = direct(u0, u1, u2)
    dig = int(-mp.log10(abs(g-d)/abs(g)))
    out.append(dict(u=(str(u0), str(u1), str(u2)), gamma=mp.nstr(g, 25),
                    direct=mp.nstr(d, 25), agree_digits=dig,
                    wall_s=round(time.time()-t0, 2)))
    print(f"u=({u0},{u1},{u2}): gamma={mp.nstr(g,20)} direct={mp.nstr(d,20)} "
          f"agree to {dig} digits [{out[-1]['wall_s']}s]")

with open('control_oracle_check.json', 'w') as f:
    json.dump(out, f, indent=1)
ok = all(c['agree_digits'] >= 25 for c in out)
print("ORACLE SELF-CHECK:", "PASS" if ok else "FAIL")

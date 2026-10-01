"""MANDATORY gate for the random-fission engine (fail-loud, rc!=0 on any fail).

G-A: certified Bessel-ratio ladder vs mpmath at moderate and large arguments.
G-B: EH2011 Table 3 (ALL SIX forests, printed values incl. Lambir; digits
     verified against the publisher PDF)
     — the engine evaluated AT the printed optima must reproduce the printed
     log-likelihoods to |diff| < 0.02 (their print rounding + parameter
     rounding; BCI rf measured 0.013 from their rounded params).
     Korup rf theta=inf: skipped here (ridge; see fit_sixforest assertion).
"""
import sys

from mpmath import mp, besseli

from rf_engine import bessel_ratio_ladder, logP_rf
from ball_engine import K_ball, logP_ball

mp.dps = 40
fails = 0

# G-A
for x, nus in [(7.0, [1, 5, 20, 80]), (1190.2, [1, 100, 1000, 5000])]:
    lad = bessel_ratio_ladder(x, max(nus) + 10, 128)
    for nu in nus:
        ref = besseli(nu + 1, x) / besseli(nu, x)
        d = abs(float(lad[nu].mid()) - float(ref))
        if d > 1e-25:
            fails += 1
            print(f"G-A FAIL x={x} nu={nu} diff={d:.2e}")
print("G-A Bessel ladder: PASS" if not fails else "G-A: FAILURES")

# G-B: printed Table 3 (six rows)
T3 = {
 'bci':       ((47.67, 0.093, -308.73),  (595.1, 0.0029, -311.92)),
 'korup':     ((52.73, 0.547, -317.04),  (None, 0.0020, -318.67)),
 'pasoh':     ((190.9, 0.093, -359.38),  (1528.0, 0.0098, -363.75)),
 'sinharaja': ((436.8, 0.0019, -252.93), (927.6, 0.0019, -252.88)),
 'yasuni':    ((204.2, 0.429, -297.15),  (10980.0, 0.0111, -306.75)),
 'lambir':    ((285.6, 0.115, -386.38),  (2500.0, 0.0111, -402.32)),
}
for forest, (pm, rf) in T3.items():
    D = sorted(int(x) for x in open(f'../data/volkov2005/derived/{forest}_abundvec.txt'))
    th, m, LL = pm
    v = float(logP_ball(D, repr(th), repr(m), 128, Kpoly=K_ball(D, 128)).mid())
    if abs(v - LL) > 0.02:
        fails += 1
        print(f"G-B FAIL {forest} pm: {v:.4f} vs printed {LL}")
    th_rf, m_rf, LL_rf = rf
    if th_rf is not None:
        vr = float(logP_rf(D, repr(th_rf), repr(m_rf), 128).mid())
        if abs(vr - LL_rf) > 0.02:
            fails += 1
            print(f"G-B FAIL {forest} rf: {vr:.4f} vs printed {LL_rf}")
    print(f"G-B {forest}: pm {v:.4f} (print {LL})"
          + (f" | rf {vr:.4f} (print {LL_rf})" if th_rf is not None else " | rf skipped (theta=inf)"))

if fails:
    print(f"*** {fails} FAILURES ***")
    sys.exit(1)
print("GATE_RF: ALL PASS")

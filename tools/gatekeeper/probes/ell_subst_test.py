#!/usr/bin/env python3
"""Acceptance tests for gatekeeper.probes.ell_subst (gates >=40d):  T1 plain elliptic K(m);  T2 sn-rule vs mp.quad with nontrivial h;
T3 Gauss-Chebyshev branch vs mp.quad;  T4 one elliptic K-line vs a reference
ell_kline.py (optional, not shipped);
T5/T5b honest-unsupported (same-side pair, interior root);
T6 cubic_under_sqrt_K closed form;
T7/T7b one-sided sn^2 branch, right outer root (closed form, nontrivial h,
smooth-factor folding);  T8/T8b/T8c one-sided sn^2, left mirror + hard
pinch (closed form) + moderate pinch with nontrivial h.
Run: python3 tools/gatekeeper/probes/ell_subst_test.py"""
import os, sys
import mpmath as mp
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
# T4 cross-checks against a reference ell_kline.py (not shipped in this repo);
# point ELL_KLINE_DIR at a dir containing it, else T4 SKIPs loudly.
KLINE_DIR = os.environ.get('ELL_KLINE_DIR', '')
sys.path.insert(0, HERE)
if KLINE_DIR: sys.path.insert(0, KLINE_DIR)
from ell_subst import select, cubic_under_sqrt_K

x = sp.symbols('x')
FAILED = []


def gate(name, got, ref, d=40):
    err = abs(got - ref)/max(abs(ref), mp.mpf('1e-30'))
    dig = float(-mp.log10(err)) if err > 0 else 999.0
    ok = dig >= d
    print(f"  {name}: {dig:.1f} digits vs ref (gate {d})  "
          f"{'PASS' if ok else 'FAIL'}")
    if not ok:
        FAILED.append(name)


# T1: int_{-1}^{1} dx/sqrt((1-x^2)(4-x^2)) = K(1/4)  (jacobi-sn branch, h=1)
mp.mp.dps = 80
plan = select([1 - x**2, 4 - x**2], -1, 1, N=64, dps=70)
assert plan.kind == 'jacobi-sn', plan.report
gate('T1 K(1/4) h=1      ', plan.integrate(), mp.ellipk(mp.mpf(1)/4), 60)
gate('T1 modulus m       ', plan.m, mp.mpf(1)/4, 60)

# T2: same chain, h = exp(x)*(2+x)  vs adaptive mp.quad at dps 110
h2 = lambda t: mp.exp(t)*(2 + t)
with mp.workdps(110):
    ref2 = mp.quad(lambda t: h2(t)/mp.sqrt((1 - t*t)*(4 - t*t)),
                   [-1, 0, 1], maxdegree=12)
gate('T2 sn nontrivial h ', plan.integrate(h2), ref2, 40)

# T3: chebyshev branch: int_{-1}^1 exp(x)/sqrt((1-x^2)(x^2+4)) dx
plan3 = select([1 - x**2, x**2 + 4], -1, 1, N=96, dps=70)
assert plan3.kind == 'chebyshev', plan3.report
with mp.workdps(110):
    ref3 = mp.quad(lambda t: mp.exp(t)/mp.sqrt((1 - t*t)*(t*t + 4)),
                   [-1, 0, 1], maxdegree=12)
gate('T3 chebyshev branch', plan3.integrate(mp.exp), ref3, 40)

# T4: elliptic K-line (s,t,q1,q2)=(-3,-5,2,3) vs ell_kline.kline S0[0]
try:
    import ell_kline
except ImportError:
    ell_kline = None
if ell_kline is None:
    print('  T4 SKIP: reference ell_kline.py not present (set ELL_KLINE_DIR;'
          ' optional cross-check leg)')
else:
    ell_kline.NPOOL = [96]                   # trim pool: probe-size run
    mp.mp.dps = 60
    s_v, t_v = mp.mpf(-3), mp.mpf(-5)
    q1, q2 = mp.mpf(2), mp.mpf(3)
    P5 = q1*q2 + s_v                             # = 3
    ftsut_xy = 4*t_v*(s_v + t_v)*q1*q2           # = 960
    ctx = ell_kline.build_pool(60)
    Nk, S0, SK, SK2 = ell_kline.kline(s_v, t_v, q1, q2, P5, ftsut_xy, ctx)
    with mp.workdps(90):
        sd = mp.sqrt(mp.mpf(13))
        a4, b4 = (1 - sd)/2, (1 + sd)/2          # roots of P2 = -K^2+K+3
    mp.mp.dps = 80
    plan4 = select([3 - x, x + 2, -x**2 + x + 3, (3*x + 5)**2 + 960],
                   a4, b4, N=160, dps=70)
    assert plan4.kind == 'jacobi-sn', plan4.report
    gate('T4 m = 13/25       ', plan4.m, mp.mpf(13)/25, 60)
    gate('T4 elliptic K-line ', plan4.integrate(), S0[0], 40)
    print(f"  T4 info: ell_kline used N={Nk}; plan gamma="
          f"{mp.nstr(plan4.gamma, 5)} (symmetric outer pair -> 0), "
          f"rate {plan4.report['digits_per_node']:.2f} d/node "
          f"(log-moments {plan4.report['logmoment_digits_per_node']:.2f})")

    # T4b: log_frame hook matches ell_kline's lt2 table (same c-table frame)
    th, tf, lt2 = plan4.log_frame()
    chk = max(abs(l - mp.log(1 - t*t)) for t, l in zip(tf, lt2))
    print(f"  T4b log_frame self-consistency max err {mp.nstr(chk, 3)} "
          f"{'PASS' if chk < mp.mpf('1e-70') else 'FAIL'}")
    if chk >= mp.mpf('1e-70'):
        FAILED.append('T4b')

# T5: honest unsupported: TWO pinching outer roots on the SAME side
plan5 = select([x*(1 - x)*(3 - x)*(sp.Rational(7, 2) - x)], 0, 1, N=32, dps=40)
ok5 = plan5.kind == 'unsupported' and 'outer real roots' in plan5.report['reason']
print(f"  T5 unsupported same-side pair: kind={plan5.kind}  "
      f"{'PASS' if ok5 else 'FAIL'}")
if not ok5:
    FAILED.append('T5')

# T5b: honest unsupported: interior root (the product changes sign)
plan5b = select([x*(1 - x)*(x - sp.Rational(1, 2))], 0, 1, N=32, dps=40)
ok5b = (plan5b.kind == 'unsupported'
        and 'interior root' in plan5b.report['reason'])
print(f"  T5b unsupported interior root: kind={plan5b.kind}  "
      f"{'PASS' if ok5b else 'FAIL'}")
if not ok5b:
    FAILED.append('T5b')

# T6: cubic_under_sqrt_K vs mp.quad  (G-R 3.131; roots 0,1,3)
mp.mp.dps = 80
val6, info6 = cubic_under_sqrt_K(2, -8, 6, 0, dps=70)   # 2z(z-1)(z-3)
with mp.workdps(110):
    ref6 = mp.quad(lambda z: 1/mp.sqrt(2*z*(1 - z)*(3 - z)),
                   [mp.mpf(0), mp.mpf(1)/2, mp.mpf(1)], maxdegree=12)
gate('T6 cubic->K(m)     ', val6, ref6, 40)
gate('T6 modulus m=1/3   ', info6['m'], mp.mpf(1)/3, 40)

# T7: one-sided sn^2, right outer root:
#     int_0^1 dx/sqrt(x(1-x)(3-x)) = (2/sqrt(3)) K(1/3)   (G-R 3.131)
mp.mp.dps = 80
plan7 = select([x*(1 - x)*(3 - x)], 0, 1, N=64, dps=70)
assert plan7.kind == 'jacobi-sn2' and plan7.side == 'right', plan7.report
ref7 = 2/mp.sqrt(3)*mp.ellipk(mp.mpf(1)/3)
gate('T7 sn2 right h=1   ', plan7.integrate(), ref7, 60)
gate('T7 modulus m=1/3   ', plan7.m, mp.mpf(1)/3, 60)
with mp.workdps(110):
    ref7h = mp.quad(lambda t: mp.exp(t)/mp.sqrt(t*(1 - t)*(3 - t)),
                    [0, mp.mpf(1)/2, 1], maxdegree=12)
gate('T7 sn2 nontrivial h', plan7.integrate(mp.exp), ref7h, 40)

# T7b: sn^2 with an extra smooth factor folded into the weights
plan7b = select([x, 1 - x, 3 - x, x**2 + 4], 0, 1, N=96, dps=70)
assert plan7b.kind == 'jacobi-sn2', plan7b.report
with mp.workdps(110):
    ref7b = mp.quad(lambda t: mp.exp(t)/mp.sqrt(t*(1 - t)*(3 - t)*(t*t + 4)),
                    [0, mp.mpf(1)/2, 1], maxdegree=12)
gate('T7b sn2 smooth fold', plan7b.integrate(mp.exp), ref7b, 40)

# T8: left mirror (x -> 1-x of T7): int_0^1 dx/sqrt(x(1-x)(x+2)) = same ref
plan8 = select([x*(1 - x)*(x + 2)], 0, 1, N=64, dps=70)
assert plan8.kind == 'jacobi-sn2' and plan8.side == 'left', plan8.report
gate('T8 sn2 left mirror ', plan8.integrate(), ref7, 60)
with mp.workdps(110):
    ref8h = mp.quad(lambda t: mp.exp(t)/mp.sqrt(t*(1 - t)*(t + 2)),
                    [0, mp.mpf(1)/2, 1], maxdegree=12)
gate('T8 sn2 left h      ', plan8.integrate(mp.exp), ref8h, 40)

# T8b: hard pinch, outer root at 101/100: 2/sqrt(101/100) K(100/101)
plan8b = select([x*(1 - x)*(sp.Rational(101, 100) - x)], 0, 1, N=96, dps=70)
assert plan8b.kind == 'jacobi-sn2', plan8b.report
ref8b = 2/mp.sqrt(mp.mpf(101)/100)*mp.ellipk(mp.mpf(100)/101)
gate('T8b sn2 hard pinch ', plan8b.integrate(), ref8b, 40)
print(f"  T8b info: m={mp.nstr(plan8b.m, 8)}, "
      f"rate {plan8b.report['digits_per_node']:.2f} d/node")

# T8c: moderate pinch (outer root 11/10), nontrivial h vs mp.quad
plan8c = select([x*(1 - x)*(sp.Rational(11, 10) - x)], 0, 1, N=96, dps=70)
assert plan8c.kind == 'jacobi-sn2', plan8c.report
with mp.workdps(110):
    ref8c = mp.quad(lambda t: mp.exp(t)/mp.sqrt(t*(1 - t)*(mp.mpf(11)/10 - t)),
                    [0, mp.mpf(1)/2, 1], maxdegree=12)
gate('T8c sn2 pinch + h  ', plan8c.integrate(mp.exp), ref8c, 40)

print('ELL_SUBST ACCEPTANCE:', 'PASS' if not FAILED else f'FAIL {FAILED}')
sys.exit(1 if FAILED else 0)

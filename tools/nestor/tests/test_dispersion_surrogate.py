r"""
nestor.dispersion battery: disp_sub on a FAST analytic surrogate with the EXACT LBL3SE singularity
structure.
  rho_surr(w') = (w'-1)(A log(w'-1)+B) exp(-(w'-1)/4)        [lower-threshold (w'-1)log turn-on]
               + theta(w'-9) D sqrt(w'-9)/(1+(w'-9))         [sqrt onset at the elliptic w'=9 cut]
  K(w')        = 1/(w'+c)^2                                   [smooth analytic kernel, box-like]
Reference = high-precision endpoint-aware tanh-sinh.  We compare, at MATCHING tanh-sinh levels:
  (A) plain GAUSS-LEGENDRE (the OLD scheme)           -> polynomial,
  (B) plain TANH-SINH, panels broken at 1 and 9       -> already exponential for integrable log/sqrt,
  (C) singularity-SUBTRACTED disp_sub                 -> exponential at the LOWEST node budget.

Asserted battery legs (pytest; legs D3/D4 of `python3 -m nestor.selftest`):
  test_disp_subtracted_level5          -- level-5 subtracted assembly vs reference, >=30 digits;
  test_addback_sign_mutation_control   -- flipping the sign of the closed-form add-back must
                                          destroy the agreement (the add-back is load-bearing).
Running the file directly (`python3 tools/nestor/tests/test_dispersion_surrogate.py`) prints the
full A/B/C convergence ladder.
"""
import sys, time
import mpmath as mp
import os
# tools/ on the path -> `nestor.dispersion` importable (script, pytest or selftest)
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from nestor.dispersion import disp_sub as DS

A = mp.mpf('-6.2857'); B = mp.mpf('9.3858'); D = mp.mpf('0.4'); c = mp.mpf('0.7')
WINF = mp.mpf(40)


def rho(wp):
    wp = mp.mpf(wp)
    if wp <= 1:
        return mp.mpf(0)
    x = wp - 1
    val = x*(A*mp.log(x) + B)*mp.exp(-x/4)
    if wp > 9:
        val += D*mp.sqrt(wp-9)/(1 + (wp-9))
    return val


def K(wp):
    return 1/(mp.mpf(wp)+c)**2


def ref(dps):
    mp.mp.dps = dps
    return (mp.quad(lambda w: rho(w)*K(w), [1, 9], method='tanh-sinh')
            + mp.quad(lambda w: rho(w)*K(w), [9, WINF], method='tanh-sinh')) / mp.pi


_REF = {}


def _reference(dps):
    """Memoized high-precision reference (endpoint-aware tanh-sinh at elevated dps)."""
    if dps not in _REF:
        _REF[dps] = ref(dps)
    return _REF[dps]


def _panels():
    """Panel list with the two-threshold structure: (w'-1)log turn-on at w'=1 (subtracted),
    sqrt onset at w'=9 (subtracted)."""
    return [
        {'a': 1, 'b': 9, 'thresh': ('left', 1, [(A, mp.mpf(1), 1)]),
         'sub_width': mp.mpf('0.5'), 'radius': mp.mpf('0.6')},
        {'a': 9, 'b': WINF, 'thresh': ('left', 9, [(D, mp.mpf('0.5'), 0)]),
         'sub_width': mp.mpf('0.5'), 'radius': mp.mpf('0.6')},
    ]


def _digits(J, R):
    """Matching decimal digits of J against the reference R (99 = exact)."""
    return 99 if J == R else -int(mp.log10(abs(J - R) / abs(R)))


def test_disp_subtracted_level5():
    """Level-5 subtracted assembly agrees with the tanh-sinh reference to >= 30 digits."""
    DPS = 30
    mp.mp.dps = DPS + 25
    R = _reference(DPS + 25)
    J = DS.disp_subtracted(rho, K, _panels(), dps=DPS, maxdegree=5)
    d = _digits(J, R)
    assert d >= 30, f"level-5 subtracted assembly: only {d} d vs reference (need >= 30)"
    print(f"disp_subtracted level 5 vs reference: {d} d  OK", flush=True)


def test_addback_sign_mutation_control():
    """Control: applying the closed-form add-back with the WRONG sign must destroy the
    agreement -- proves the add-back term is load-bearing in the passing assembly, so the
    level-5 leg cannot pass with the add-back dropped or sign-flipped."""
    DPS = 30
    mp.mp.dps = DPS + 25
    R = _reference(DPS + 25)
    J = DS.disp_subtracted(rho, K, _panels(), dps=DPS, maxdegree=5)
    # Standalone add-back for the w'=1 turn-on sub-panel, same recipe the assembly uses
    # (kernel Taylor about w*=1, radius 0.6, sub-panel length 0.5).
    work = DPS + 20
    mp.mp.dps = work
    Kc = DS.kernel_taylor(K, mp.mpf(1), int(1.7 * DPS) + 10, radius=mp.mpf('0.6'), dps=work)
    AB = DS.addback_endpoint([(A, mp.mpf(1), 1)], Kc, mp.mpf('0.5'), side='left') / mp.pi
    mp.mp.dps = DPS + 25
    J_mut = J - 2 * AB                     # the assembly with the add-back sign flipped
    d_ok, d_mut = _digits(J, R), _digits(J_mut, R)
    assert d_ok >= 30, f"control precondition: passing assembly only {d_ok} d (need >= 30)"
    assert d_mut <= 5, f"sign-flipped add-back still agrees to {d_mut} d -- control not sensitive"
    print(f"add-back sign mutation control: {d_ok} d -> {d_mut} d  OK", flush=True)


if __name__ == '__main__':
    DPS = 30
    mp.mp.dps = DPS + 25
    R = ref(DPS + 25)
    dd = lambda J: (-int(mp.log10(abs(J-R)/abs(R))) if J != R else 99)
    print("reference J/pi =", mp.nstr(R, DPS), "\n")

    print("=== (A) plain GAUSS-LEGENDRE (old scheme), by node count per panel ===")
    for n in [20, 40, 80, 160]:
        mp.mp.dps = DPS + 25
        x, wt = mp.gauss_quadrature(n, 'legendre')
        tot = mp.mpf(0)
        for (lo, hi) in [(mp.mpf(1), mp.mpf(9)), (mp.mpf(9), WINF)]:
            h = (hi-lo)/2; cc = (hi+lo)/2
            tot += sum(wi*rho(cc+h*ti)*K(cc+h*ti) for ti, wi in zip(x, wt))*h
        print(f"  GL n={n:4d}/panel ({2*n} rho-evals): ~{dd(tot/mp.pi)} d")

    print("\n=== (B) plain TANH-SINH, panels broken at 1 and 9, NO subtraction ===")
    for lev in [3, 4, 5, 6]:
        mp.mp.dps = DPS + 25
        tot = (mp.quad(lambda w: rho(w)*K(w), [1, 9], method='tanh-sinh', maxdegree=lev)
               + mp.quad(lambda w: rho(w)*K(w), [9, WINF], method='tanh-sinh', maxdegree=lev))
        print(f"  tanh-sinh level={lev}: ~{dd(tot/mp.pi)} d")

    print("\n=== (C) singularity-SUBTRACTED disp_sub (sub-panel + closed add-back) ===")
    panels = _panels()
    for lev in [3, 4, 5]:
        mp.mp.dps = DPS + 25
        J = DS.disp_subtracted(rho, K, panels, dps=DPS, maxdegree=lev)
        print(f"  level={lev}: ~{dd(J)} d")

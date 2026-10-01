#!/usr/bin/env python3
"""certsfs — certified selection-SFS evaluator (CLI front door).

The single-file front door to the certified selection-SFS engine. Anyone with
python3 + mpmath (+ python-flint for large n) can reproduce certified fixed-S
selection-SFS numbers.

LAYOUT: 'entry', 'check' and 'selftest' modes are self-contained (python3 +
mpmath only). 'vector' mode imports the engine modules shipped beside this
file: sfs_engine.py (exact-rational tridiagonal route, n<=2000) and
arb_route.py (ball recursion above; needs python-flint). Both are
SHA-256-pinned below and verified before import; so is REFERENCE_VALUES.json
(pinned reference values, shipped beside this file). Any mismatch is a named
CertSFSGateError refusal with rc=3 — never a silent fallback.

Usage:
  certsfs.py selftest [--dps D]            verify pinned files + reproduce
                                           two pinned reference values
                                           digit-for-digit + dps-consistency
  certsfs.py entry  N I S [--dps D]        one SFS entry E_i(S), theta=1
  certsfs.py vector N S [--dps D] [--grad] whole vector (i=1..N-1) [+ dE/dS]
  certsfs.py check  N I S [--dps D]        entry by TWO independent routes +
                                           agreement digits (the honesty gate)

Conventions: S = 4*Ne*s, S>0 advantageous, unfolded derived-allele SFS.
polyDFE/fastDFE S == this S; dadi/fitdadi gamma == S/2.

Routes: closed form theta*n/(i(n-i))*(1-1F1(n-i;n;-S))/(1-e^{-S}) via mpmath;
exact-rational certified series (positive-term Kummer arrangement, geometric
tail bound) as the independent check; exact QQ(e^{-S}) tridiagonal solve for
vectors at n<=2000, ball recursion above (needs python-flint).
"""
import argparse, hashlib, json, os, sys


class CertSFSGateError(Exception):
    """Named failure: a gate did not pass (route disagreement below target,
    pinned-file tamper/missing, or selftest receipt mismatch) — do not use
    any value from this run. Raised so callers cannot mistake a gate failure
    for an evaluation result. Process exit code is 3."""


HERE = os.path.dirname(os.path.abspath(__file__))
ENGINE = HERE                      # flat layout: engine files ship beside the CLI

# SHA-256 pins of the shipped engine files and reference values (regenerated
# whenever the shipped bytes change; see _pins.py for the package-wide set).
PINNED_SHA256 = {
    'sfs_engine.py':
        'cff99bdfd1f530a9ca9e4bd0328ca176a47adafb5465810d7998e09bf195a0df',
    'arb_route.py':
        '8ecd6a94cd2d898259229e85607533b871a99d623544075b7f962915517cea92',
    'REFERENCE_VALUES.json':
        '22c019a601a6efafbebad0df0b036bcf272690df7a38a7615d929a6b499794fc',
}


def _pin_dir(name):
    return HERE


# The two pinned reference values reproduced by selftest (h=1/2 genic points
# of REFERENCE_VALUES.json, where the dominance oracle collapses to the
# closed form; oracle strings carry ~60 certified digits, printed to 40/32).
SELFTEST_POINTS = [
    {'n': 20, 'k': 10, 'S': '500', 'h': 0.5},
    {'n': 20, 'k': 3, 'S': '-100', 'h': 0.5},
]


def vendored(name):
    """Return the path of pinned file `name` after SHA-256 verification.
    Missing file or hash mismatch -> named CertSFSGateError (fail closed)."""
    path = os.path.join(_pin_dir(name), name)
    if not os.path.isfile(path):
        raise CertSFSGateError(f"pinned file {name} MISSING at {_pin_dir(name)}")
    with open(path, 'rb') as f:
        h = hashlib.sha256(f.read()).hexdigest()
    if h != PINNED_SHA256[name]:
        raise CertSFSGateError(
            f"pinned file {name} SHA-256 mismatch (got {h[:16]}..., "
            f"pinned {PINNED_SHA256[name][:16]}...) — tampered or wrong copy; refusing")
    return path


from fractions import Fraction
from mpmath import mp, mpf, exp, expm1, hyp1f1, fabs, log10


def entry_closed(n, i, S, dps):
    with mp.workdps(dps + 20):
        Sm = mpf(S.numerator) / S.denominator
        if Sm == 0:
            return mpf(1) / i
        return mpf(n) / (mpf(i) * (n - i)) * (1 - hyp1f1(n - i, n, -Sm)) / (-expm1(-Sm))


def entry_series(n, i, S, dps):
    """Independent route: exact-rational positive-term series."""
    zfr = -S if S < 0 else S
    a, b = (n - i, n) if S < 0 else (i, n)
    t = Fraction(1); s = Fraction(1); k = 0
    tol = Fraction(1, 10 ** (dps + 30))
    while True:
        t *= Fraction(zfr.numerator * (a + k), zfr.denominator * (b + k) * (k + 1))
        s += t
        k += 1
        r = Fraction(zfr.numerator * (a + k), zfr.denominator * (b + k) * (k + 1))
        if r < Fraction(1, 2) and t < s * tol:
            break
    with mp.workdps(dps + 20):
        R = mpf(s.numerator) / mpf(s.denominator)
        Sm = mpf(S.numerator) / S.denominator
        M = R if S < 0 else exp(-Sm) * R
        return mpf(n) / (mpf(i) * (n - i)) * (1 - M) / (-expm1(-Sm))


def _sig_digits(dec_str):
    """Count significant digits of a plain decimal string."""
    ds = dec_str.lstrip('-').replace('.', '').lstrip('0')
    if 'e' in dec_str or 'E' in dec_str:
        ds = dec_str.lstrip('-').split('e')[0].split('E')[0].replace('.', '').lstrip('0')
    return len(ds)


def _agree_digits(v, w, dps):
    with mp.workdps(dps + 20):
        if v == w:
            return 9999.0
        return float(-log10(fabs((v - w) / w)))


def selftest(dps):
    """Verify pinned files, reproduce the two pinned reference values
    digit-for-digit, cross-check the two independent routes, and check
    dps vs 2*dps arbitrary-precision consistency. Any failure -> named
    CertSFSGateError, rc=3. Prints one PASS line per sub-check."""
    for name in ('sfs_engine.py', 'arb_route.py', 'REFERENCE_VALUES.json'):
        vendored(name)
        print(f"PASS sha256 {name}")
    with open(vendored('REFERENCE_VALUES.json')) as f:
        bank = json.load(f)
    for want in SELFTEST_POINTS:
        n, k, Sfr, h = want['n'], want['k'], Fraction(want['S']), want['h']
        rec = [p for p in bank['points']
               if p['n'] == n and p['k'] == k and p['h'] == h
               and Fraction(p['S']) == Sfr]
        if len(rec) != 1:
            raise CertSFSGateError(
                f"selftest point n={n} k={k} S={Sfr} h={h} not found uniquely "
                f"in REFERENCE_VALUES.json ({len(rec)} matches)")
        oracle = rec[0]['oracle']
        nd = _sig_digits(oracle)
        v = entry_closed(n, k, Sfr, dps)
        with mp.workdps(dps + 20):
            got = mp.nstr(v, nd)
        if got != oracle:
            raise CertSFSGateError(
                f"receipt MISMATCH at E_{k}(n={n}, S={Sfr}): computed {got} "
                f"!= pinned reference {oracle}")
        print(f"PASS receipt E_{k}(n={n},S={Sfr}) digit-for-digit "
              f"({nd} oracle digits)")
        w = entry_series(n, k, Sfr, dps)
        d = _agree_digits(v, w, dps)
        if d < dps - 5:
            raise CertSFSGateError(
                f"route cross-check at E_{k}(n={n},S={Sfr}): closed vs series "
                f"agree only {d:.1f} digits (target {dps})")
        print(f"PASS routes closed/series agree {min(d, 9999):.1f} digits "
              f"(bar {dps - 5})")
        v2 = entry_closed(n, k, Sfr, 2 * dps)
        d2 = _agree_digits(v, v2, 2 * dps)
        if d2 < dps - 2:
            raise CertSFSGateError(
                f"dps-consistency at E_{k}(n={n},S={Sfr}): dps {dps} vs "
                f"{2 * dps} agree only {d2:.1f} digits (bar {dps - 2})")
        print(f"PASS dps-consistency {dps} vs {2 * dps}: {min(d2, 9999):.1f} "
              f"digits (bar {dps - 2})")
    print(f"SELFTEST PASS (all checks, dps {dps})")
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('mode', choices=['entry', 'vector', 'check', 'selftest'])
    ap.add_argument('N', type=int, nargs='?')
    ap.add_argument('args', nargs='*')
    ap.add_argument('--dps', type=int, default=40)
    ap.add_argument('--grad', action='store_true')
    a = ap.parse_args()
    if a.mode == 'selftest':
        return selftest(a.dps)
    if a.N is None:
        ap.error('N required for entry/vector/check')
    n = a.N
    if a.mode in ('entry', 'check'):
        i = int(a.args[0]); S = Fraction(a.args[1])
        v = entry_closed(n, i, S, a.dps)
        if a.mode == 'entry':
            print(mp.nstr(v, a.dps))
            return 0
        w = entry_series(n, i, S, a.dps)
        d = _agree_digits(v, w, a.dps)
        print(f"closed : {mp.nstr(v, a.dps)}")
        print(f"series : {mp.nstr(w, a.dps)}")
        print(f"agree  : {d:.1f} digits (target {a.dps})")
        if d < a.dps - 5:
            raise CertSFSGateError(
                f"routes agree to only {d:.1f} digits (target {a.dps}); value NOT certified")
        return 0
    # vector: verify pinned engine bytes BEFORE import (fail closed)
    vendored('sfs_engine.py')
    vendored('arb_route.py')
    sys.path.insert(0, ENGINE)
    S = Fraction(a.args[0])
    if n <= 2000:
        from sfs_engine import solve_M_exact, solve_dM_exact, assemble_M
        alpha, beta = solve_M_exact(n, S)
        M = assemble_M(alpha, beta, S, a.dps + 10)
        G = None
        if a.grad:
            ders = solve_dM_exact(n, S, alpha, beta, 1)
            a1, b1 = ders[0]
            G = assemble_M(a1, b1, S, a.dps + 10)
        with mp.workdps(a.dps + 15):
            Sm = mpf(S.numerator) / S.denominator
            den = -expm1(-Sm)
            for i in range(1, n):
                pref = mpf(n) / (mpf(i) * (n - i))
                E = pref * (1 - M[i]) / den
                line = f"{i}\t{mp.nstr(E, a.dps)}"
                if G is not None:
                    eS = exp(-Sm)
                    dE = pref * (-G[i] * den - (1 - M[i]) * eS) / (den * den)
                    line += f"\t{mp.nstr(dE, a.dps)}"
                print(line)
    else:
        from arb_route import solve_vector_arb, sfs_arb
        sol = solve_vector_arb(n, float(S), a.dps, want_grad=a.grad, S_frac=S)
        for i in range(1, n):
            E, dE = sfs_arb(n, float(S), sol, i, S_frac=S)
            line = f"{i}\t{E.mid().str(a.dps, radius=False)}"
            if a.grad and dE is not None:
                line += f"\t{dE.mid().str(a.dps, radius=False)}"
            print(line)
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except CertSFSGateError as e:
        print(f"CertSFSGateError: {e}", file=sys.stderr)
        sys.exit(3)

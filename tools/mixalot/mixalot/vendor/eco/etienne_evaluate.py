#!/usr/bin/env python3
"""etienne_evaluate.py — exact / certified evaluation of the Etienne (2005)
neutral-biodiversity sampling formula.

    P[D|theta,m,J] = J!/(prod_i n_i * prod_j Phi_j!) * theta^S/(I)_J
                     * sum_{A=S}^J K(D,A) I^A/(theta)_A,   I = m(J-1)/(1-m)

Transcription certified against an independent exact implementation of the
neutral model itself (two-level Hoppe/CRP urn DP) and against Etienne's own
published worked-example values (see pilot/ gates).

Modes:
  --mode ball   (default) certified ball arithmetic (python-flint arb).
                Overflow-impossible; positivity makes enclosures tight.
                Output: log-likelihood ball [mid +/- rad] at --prec bits.
  --mode exact  exact rational P[D] at rational theta, m (integer common-
                denominator scheme, no per-term gcds). Slow-ish at J>10^4 but
                EXACT (the Phase-1 gate object).
  --lnk-table   print ln K(D,A) for A=S..J (16 digits, Etienne SA3 format)
  --gate-lnk F  compare our ln K(D,A) against a reference file (one value per
                line, A=S..J; e.g. Etienne's SA3 *.out files) and report
                agreement statistics.
  --selftest    run the worked-example gate D=(1,1,2,3,5,8).

Dataset input: --abund "1,1,2,3,5,8" or --file PATH (whitespace/newline-
separated integers; lines starting with # ignored).

Usage examples:
  python3 etienne_evaluate.py --abund "1,1,2,3,5,8" --theta 7.047958 --m 0.22635923
  python3 etienne_evaluate.py --file bci.txt --theta 48.4 --m 0.133 --prec 256
  python3 etienne_evaluate.py --file bci.txt --theta 47.94 --m 0.092 --mode exact
  python3 etienne_evaluate.py --file bci1982.txt --gate-lnk BCICensus1982-lnK.out
"""
import argparse
import json
import sys
import time
from collections import Counter
from fractions import Fraction
from math import factorial

from flint import fmpz_poly, arb_poly, arb, ctx

sys.set_int_max_str_digits(0) if hasattr(sys, "set_int_max_str_digits") else None


# ---------------- kernel (integer product tree) ----------------

_rf_cache, _sp_cache = {}, {}


def _rising_factorial_poly(n):
    if n in _rf_cache:
        return _rf_cache[n]
    factors = [fmpz_poly([k, 1]) for k in range(n)]
    while len(factors) > 1:
        factors = [factors[i] * factors[i + 1] for i in range(0, len(factors) - 1, 2)] \
                  + ([factors[-1]] if len(factors) % 2 else [])
    _rf_cache[n] = factors[0]
    return factors[0]


def species_poly_int(n):
    """fmpz_poly v_n(z) = sum_a s(n,a)(a-1)! z^a."""
    if n in _sp_cache:
        return _sp_cache[n]
    rf = _rising_factorial_poly(n)
    p = fmpz_poly([0] + [int(rf[a]) * factorial(a - 1) for a in range(1, n + 1)])
    _sp_cache[n] = p
    return p


def _tree(polys):
    while len(polys) > 1:
        polys = [polys[i] * polys[i + 1] for i in range(0, len(polys) - 1, 2)] \
                + ([polys[-1]] if len(polys) % 2 else [])
    return polys[0]


def K_int(D):
    """Integer K-numerators: K(D,A) = Kint[A-S] / prod_i (n_i-1)!."""
    prod = _tree([species_poly_int(n) for n in D])
    S, J = len(D), sum(D)
    F = 1
    for n in D:
        F *= factorial(n - 1)
    return S, J, [int(prod[A]) for A in range(S, J + 1)], F


def K_ball_poly(D, prec):
    ctx.prec = prec
    polys = []
    for n in D:
        vi = species_poly_int(n)
        d = factorial(n - 1)
        polys.append(arb_poly([arb(int(vi[a])) / arb(d) for a in range(n + 1)]))
    return _tree(polys)


# ---------------- evaluation ----------------

def logP_ball(D, theta_str, m_str, prec, KP=None):
    ctx.prec = prec
    D = sorted(D)
    J, S = sum(D), len(D)
    th, mm = arb(theta_str), arb(m_str)
    I = mm * (J - 1) / (1 - mm)
    if KP is None:
        KP = K_ball_poly(D, prec)
    tot = arb(0)
    t = I ** S
    for k in range(S):
        t /= th + k
    for A in range(S, J + 1):
        tot += KP[A] * t
        t *= I / (th + A)
    lp = arb(J + 1).lgamma()
    for n in D:
        lp -= arb(n).log()
    for c in Counter(D).values():
        lp -= arb(c + 1).lgamma()
    lp += S * th.log()
    lp -= (I + J).lgamma() - I.lgamma()
    lp += tot.log()
    return lp


def P_exact(D, theta, m):
    """Exact rational P[D]; theta, m Fractions. Integer common-denominator
    scheme: L = N / (F s^J T_0) with all-integer accumulation, gcd once."""
    D = sorted(D)
    J, S = sum(D), len(D)
    theta, m = Fraction(theta), Fraction(m)
    if m == 0:
        return Fraction(int(S == 1))
    I = m * (J - 1) / (1 - m)
    S_, J_, Kint, F = K_int(D)
    p, q = theta.numerator, theta.denominator
    r, s = I.numerator, I.denominator
    # suffix products T_A = prod_{k=A}^{J-1} (p + k q); T_J = 1
    T = [1] * (J + 1)
    for A in range(J - 1, -1, -1):
        T[A] = (p + A * q) * T[A + 1]
    N = 0
    rqA = (r * q) ** S
    sJA = s ** (J - S)
    for idx, KA in enumerate(Kint):
        A = S + idx
        N += KA * rqA * sJA * T[A]
        rqA *= r * q
        if A < J:
            sJA //= s
    L = Fraction(N, F * s ** J * T[0])
    # prefactor and theta^S/(I)_J
    pref = Fraction(factorial(J))
    for n in D:
        pref /= n
    for c in Counter(D).values():
        pref /= factorial(c)
    IJ = Fraction(1)
    for k in range(J):
        IJ *= I + k
    return pref * theta ** S / IJ * L


def lnK_table(D, prec=96):
    """[ln K(D,A) for A=S..J] as arb balls (>= ~25 certified digits)."""
    ctx.prec = prec
    S, J, Kint, F = K_int(D)
    lF = arb(F).log()
    return [(arb(k).log() - lF) for k in Kint]


# ---------------- CLI ----------------

def read_dataset(args):
    if args.abund:
        D = [int(x) for x in args.abund.replace(",", " ").split()]
    else:
        toks = []
        with open(args.file) as f:
            for line in f:
                line = line.split("#")[0].replace(",", " ")
                toks += line.split()
        D = [int(x) for x in toks]
    if any(n < 1 for n in D):
        sys.exit("abundances must be >= 1")
    return sorted(D)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--abund")
    ap.add_argument("--file")
    ap.add_argument("--theta")
    ap.add_argument("--m")
    ap.add_argument("--mode", choices=["ball", "exact"], default="ball")
    ap.add_argument("--prec", type=int, default=128)
    ap.add_argument("--lnk-table", action="store_true")
    ap.add_argument("--gate-lnk")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()

    if args.selftest:
        D = [1, 1, 2, 3, 5, 8]
        S, J, Kint, F = K_int(D)
        assert Fraction(Kint[1], F) == Fraction(1507, 210)
        assert Fraction(Kint[7], F) == Fraction(27751573, 60480)
        assert Fraction(Kint[8], F) == Fraction(1127899, 2520)
        assert Kint[0] == F and Kint[-1] == F
        lp = logP_ball(D, "7.047958", "0.22635923", 128)
        assert abs(float(lp.mid()) - (-4.796686701878364)) < 5e-13
        pe = P_exact(D, Fraction(71, 10), Fraction(225, 1000))
        import math
        assert abs(math.log(pe) - (-4.796695178274603)) < 5e-13
        print("SELFTEST PASS (worked example, exact + ball + cross)")
        return

    D = read_dataset(args)
    out = {"J": sum(D), "S": len(D)}

    if args.gate_lnk:
        t0 = time.time()
        ours = lnK_table(D)
        refs = [float(x) for x in open(args.gate_lnk).read().split()]
        if len(refs) != len(ours):
            sys.exit(f"GATE FAIL: length mismatch ours={len(ours)} ref={len(refs)} "
                     f"(J-S+1={out['J']-out['S']+1})")
        worst = 0.0
        for o, r in zip(ours, refs):
            d = abs(float(o.mid()) - r)
            rel = d / max(abs(r), 1e-12)
            worst = max(worst, min(d, rel))
        out.update(gate="PASS" if worst < 5e-14 else "CHECK",
                   worst_diff=worst, n=len(ours), secs=round(time.time() - t0, 2))
        print(json.dumps(out))
        return

    if args.lnk_table:
        for v in lnK_table(D):
            print(f"{float(v.mid()):.16E}")
        return

    if not (args.theta and args.m):
        sys.exit("need --theta and --m (or --lnk-table/--gate-lnk/--selftest)")

    t0 = time.time()
    if args.mode == "ball":
        lp = logP_ball(D, args.theta, args.m, args.prec)
        out.update(mode="ball", prec=args.prec, logP=str(lp),
                   logP_mid=float(lp.mid()), logP_rad=float(lp.rad()))
    else:
        P = P_exact(D, Fraction(args.theta), Fraction(args.m))
        out.update(mode="exact",
                   P_exact=f"{P.numerator}/{P.denominator}" if P.denominator.bit_length() < 4000
                           else "(huge rational; printing suppressed)",
                   P_float=float(P))
        ctx.prec = 128
        out["logP"] = str((arb(P.numerator).log() - arb(P.denominator).log()))
    out["secs"] = round(time.time() - t0, 3)
    print(json.dumps(out))


if __name__ == "__main__":
    main()

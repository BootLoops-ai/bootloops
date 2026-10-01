#!/usr/bin/env python3
"""
restrict_op.py — restricted PF operator from the exact curve series,
general h_eff / general order (production path: multi-prime CRT + rational
reconstruction + EXACT proof, never pure-Fraction elimination).

Pipeline (engine: the annihilator module, tools/annihilator):
 1. mod-p s-first (r,s)-pin (theta-order-first scan) at prime p1 on the full
    series;
 2. cross-truncation gate: SECOND prime AND SHORTER truncation must reproduce
    the pin (same-method-different-data robustness; the independent METHOD
    gate is the D-module route in pipe_vac.py);
 3. exact reconstruction: the (r,s) linear system solved mod n_primes 31-bit
    primes, CRT-combined, Wang-rational-reconstructed, integer-normalized;
 4. r2 PROOF gate (fully exact): the recurrence annihilates ALL windows of
    the exact series; the fit uses only a strict prefix -> true held-out tail;
 5. theta-form (rec_to_theta) + indicial factorization over integer roots
    (assert all roots integer >= 0; MUM block = exponent-0 multiplicity >= 4).
"""
import json, os, sys, time
from fractions import Fraction as Fr
from math import gcd, isqrt

sys.path.insert(0, os.path.abspath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "..", "..", "..", "annihilator")))
from annihilator import (pf_from_series, find_recurrence_fast,
                         nullspace_mod_fast, rec_to_theta, P31)
import numpy as np


def pin_rs(A, rmax, smax, p2, trunc2):
    """(r,s) pin at p=2^31-1 full series + robustness at (p2, trunc2)."""
    res = pf_from_series(A, rmax=rmax, smax=smax, mode='modp', p=P31)
    assert res, f"no operator found in (rmax,smax)=({rmax},{smax})"
    r1, s1 = res[0], res[1]
    ap = [x.numerator % p2 * pow(x.denominator % p2, p2 - 2, p2) % p2
          for x in A[:trunc2]]
    hit2 = find_recurrence_fast(ap, p2, rmax=rmax, dmax=smax, overdet=5,
                                prefer='ode_order')
    assert hit2 and hit2[0] == r1 and hit2[1] == s1, \
        f"cross-truncation pin mismatch: {hit2 and hit2[:2]} vs ({r1},{s1})"
    return r1, s1


def primes31(count, start=(1 << 31) - 1):
    import sympy
    out, c = [], start
    while len(out) < count:
        if sympy.isprime(c):
            out.append(c)
        c -= 2
    return out


def ratrec(a, m):
    """Wang rational reconstruction: a mod m -> p/q, |p|, q <= sqrt(m/2)."""
    bound = isqrt(m // 2)
    r0, r1 = m, a % m
    s0, s1 = 0, 1
    while r1 > bound:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        s0, s1 = s1, s0 - q * s1
    if s1 == 0 or abs(s1) > bound:
        return None
    return Fr(r1, s1) if s1 > 0 else Fr(-r1, -s1)


def reconstruct(A, R, S, nfit, n_primes):
    """multi-prime CRT exact reconstruction of the (R,S) recurrence."""
    nunk = (R + 1) * (S + 1)
    nrows = nfit - R
    assert nrows > nunk + 20, "fit window too short"
    resids, mods, freecol = [], [], None
    for p in primes31(n_primes):
        ap = [x.numerator % p * pow(x.denominator % p, p - 2, p) % p
              for x in A[:nfit]]
        rows = np.zeros((nrows, nunk), dtype=np.int64)
        for n in range(nrows):
            npows = [pow(n, k, p) for k in range(S + 1)]
            for j in range(R + 1):
                base = ap[n + j]
                for k in range(S + 1):
                    rows[n, j * (S + 1) + k] = base * npows[k] % p
        v = nullspace_mod_fast(rows, p)
        if v is None:
            continue
        piv = max(i for i, x in enumerate(v) if x != 0)
        if freecol is None:
            freecol = piv
        if piv != freecol:
            continue
        inv = pow(int(v[freecol]), p - 2, p)
        resids.append([int(x) * inv % p for x in v])
        mods.append(p)
    comb, Mtot = resids[0][:], mods[0]
    for pi in range(1, len(mods)):
        p = mods[pi]
        dinv = pow(Mtot, -1, p)
        for i in range(nunk):
            t = (resids[pi][i] - comb[i]) * dinv % p
            comb[i] = comb[i] + Mtot * t
        Mtot *= p
    coef = {}
    for i in range(nunk):
        f = ratrec(comb[i], Mtot)
        assert f is not None, f"rational reconstruction failed at unknown {i}"
        if f != 0:
            coef[(i // (S + 1), i % (S + 1))] = f
    den = 1
    for f in coef.values():
        den = den * f.denominator // gcd(den, f.denominator)
    g = 0
    for f in coef.values():
        g = gcd(g, abs(f.numerator * (den // f.denominator)))
    coeffs = {k: Fr(f.numerator * (den // f.denominator), g)
              for k, f in coef.items()}
    assert all(v.denominator == 1 for v in coeffs.values())
    return coeffs, Mtot.bit_length()

def prove_annihilation(A, coeffs, R, nfit):
    """r2 PROOF: exact recurrence kills ALL windows; report held-out count."""
    M = len(A) - 1
    bad = []
    for n in range(0, M - R + 1):
        tot = Fr(0)
        for (j, k), c in coeffs.items():
            if A[n + j]:
                tot += c * Fr(n) ** k * A[n + j]
        if tot != 0:
            bad.append(n)
            if len(bad) > 3: break
    assert not bad, f"recurrence fails at n = {bad}"
    return M - R + 1, M - R + 1 - (nfit - R)


def indicial_factor(theta_form):
    """indicial polynomial at s=0 from the theta-form; factor over integer
    roots.  Returns (lead, [(root, mult)...], mult0).  Asserts full integer
    factorization with roots >= 0 (MUM-type restricted operators)."""
    ind = {k: int(c) for (i, k), c in theta_form.items() if i == 0}
    r = max(ind)
    poly = [ind.get(k, 0) for k in range(r + 1)]          # ascending
    roots = []
    cur = poly[:]
    for e in range(0, 4 * r + 40):
        while sum(c * e ** k for k, c in enumerate(cur)) == 0 and len(cur) > 1:
            # synthetic division by (x - e)
            q = [0] * (len(cur) - 1)
            q[-1] = cur[-1]
            for k in range(len(cur) - 2, 0, -1):
                q[k - 1] = cur[k] + e * q[k]
            rem = cur[0] + e * q[0]
            assert rem == 0
            cur = q
            roots.append(e)
        if len(cur) == 1:
            break
    assert len(cur) == 1, f"indicial not fully integer-factorable: {cur}"
    lead = cur[0]
    fac = {}
    for e in roots:
        fac[e] = fac.get(e, 0) + 1
    mult0 = fac.get(0, 0)
    return lead, sorted(fac.items()), mult0


def find_operator(card, A, outfile):
    """full Route-A production run on exact series A (list of Fr)."""
    P = card.op_search
    M = len(A) - 1
    t0 = time.time()
    r1, s1 = pin_rs(A, P["rmax"], P["smax"], P["p2"], P["trunc2"])
    print(f"[pin] (r,s) = ({r1},{s1}) at p=2^31-1 (N={M+1}) AND "
          f"p={P['p2']} (N={P['trunc2']})   {time.time()-t0:.1f}s", flush=True)
    nfit = P.get("fit_window") or (r1 + (r1 + 1) * (s1 + 1) + 60)
    assert nfit + 60 <= M, "no held-out tail left"
    t0 = time.time()
    coeffs, crt_bits = reconstruct(A, r1, s1, nfit, P["n_primes"])
    print(f"[exact] {len(coeffs)} nonzero c_(j,k), CRT modulus {crt_bits} bits"
          f"   {time.time()-t0:.1f}s", flush=True)
    t0 = time.time()
    nwin, nheld = prove_annihilation(A, coeffs, r1, nfit)
    print(f"[r2 GATE] PASS  annihilates all {nwin} windows exactly; "
          f"{nheld} TRUE held-out beyond fit window m < {nfit}"
          f"   {time.time()-t0:.0f}s", flush=True)
    ode = rec_to_theta(r1, s1, coeffs)
    assert all(c.denominator == 1 for c in ode.values())
    lead, fac, mult0 = indicial_factor(ode)
    print(f"[indicial] lead {lead}, factors {fac}; exponent-0 multiplicity "
          f"{mult0} -> MUM block {'PRESENT' if mult0 >= 4 else 'ABSENT'}")
    assert mult0 >= 4, "weight-3 MUM block absent on the restricted operator"
    out = {"r_shift": r1, "s_theta": s1,
           "recurrence": {f"{j},{k}": str(c) for (j, k), c in coeffs.items()},
           "theta_form": {f"{i},{k}": str(c) for (i, k), c in ode.items()},
           "fit_window": nfit, "M": M, "n_primes": P["n_primes"],
           "indicial": {"lead": lead, "factors": fac, "mult0": mult0},
           "pins": {"p31_full": [r1, s1],
                    f"p{P['p2']}_N{P['trunc2']}": [r1, s1]}}
    with open(outfile, "w") as f:
        json.dump(out, f, indent=1)
    print(f"[wrote] {outfile}", flush=True)
    return out

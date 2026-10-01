#!/usr/bin/env python3
"""boxops_u1.py — B1: HV4 6-var box-operator supply on the orthant Z^6_{>=0}
(the conifold-frame box-operator generator pattern of the Pfaffian route re-instantiated on
the (5.18) coefficient law). c(k) = multinom(|k|;k)^2 = (|k|!)^2 / prod(k_I!)^2,
support = FULL orthant (orthant form; no curve grading, no negative rows).
Exact ratio law: c(k+e_I)/c(k) = (|k|+1)^2/(k_I+1)^2.
L_l = phi^{l-} DEN'(th-l) - phi^{l+} NUM'(th) after primitive cancellation;
MANDATORY termwise gate on exact c(k) (NBOX=4). Supply: 0<|l|_1<=2.
Engine reuse VERBATIM: reduce_sides/poly_from_factors from the ads-5-81 bank
(H monkeypatched 5->6; guarded by the D_I printed-form exact assert below)."""
import sys, time
from fractions import Fraction as Fr
from math import factorial, gcd
from itertools import product as iproduct

# engine bank path via env TERRIER_BOXOPS_BANK (required — the box-operator
# bank is not included in the package; fail-closed).
import os
BANK = os.environ.get("TERRIER_BOXOPS_BANK")
if not BANK:
    raise ImportError("TERRIER_BOXOPS_BANK is not set (boxops bank dir)")
sys.path.insert(0, BANK)
import boxops as _bx                       # generic engine (ads-5-81 bank)
_bx.H = 6                                  # poly_from_factors/theta_mul dim
reduce_sides, poly_from_factors = _bx.reduce_sides, _bx.poly_from_factors

H = 6
ONES = (1,) * H
E = [tuple(1 if j == i else 0 for j in range(H)) for i in range(H)]
# factor list: (vector, multiplicity, is_numerator); (5.18) law
FACTORS = [(ONES, 2, True)] + [(E[i], 2, False) for i in range(H)]


def cn(k):
    """exact c(k) (Fraction); orthant support."""
    if min(k) < 0:
        return Fr(0)
    den = 1
    for x in k:
        den *= factorial(x) ** 2
    return Fr(factorial(sum(k)) ** 2, den)


def ratio_sides(l):
    """(num_factors, den_factors) of c(k+l)/c(k); factor = (vec6, const)."""
    num, den = [], []
    for (v, m, isn) in FACTORS:
        kk = sum(v[a] * l[a] for a in range(H))
        rising = num if isn else den
        falling = den if isn else num
        for _ in range(m):
            if kk > 0:
                for j in range(1, kk + 1):
                    rising.append((v, j))
            elif kk < 0:
                for j in range(0, -kk):
                    falling.append((v, -j))
    return num, den


def box_op(l):
    """L_l as {(zmono, thmono): int}, factor-reduced, content-normalized."""
    num, den = ratio_sides(l)
    Sn, Sd, rn, rd = reduce_sides(num, den)
    lp = tuple(max(x, 0) for x in l)
    lm = tuple(max(-x, 0) for x in l)
    mshift = tuple(-x for x in l)
    D = poly_from_factors(rd, mshift)
    N = poly_from_factors(rn, None)
    q = Sn / Sd
    wD, wN = q.denominator, q.numerator
    op = {}
    for m, c in D.items():
        op[(lm, m)] = op.get((lm, m), 0) + wD * c
    for m, c in N.items():
        op[(lp, m)] = op.get((lp, m), 0) - wN * c
    g = 0
    for c in op.values():
        g = gcd(g, abs(c))
    return {k: c // g for k, c in op.items() if c}


def termwise_gate(op, C, NBOX=4):
    """L_l annihilates sum c(k) phi^k: every interior output monomial == 0."""
    checked = 0
    for m in iproduct(range(NBOX + 1), repeat=H):
        tot, ok = Fr(0), True
        for (zm, tm), c in op.items():
            src = tuple(m[a] - zm[a] for a in range(H))
            if min(src) < 0:
                continue
            if max(src) > NBOX:
                ok = False
                break
            cv = C[src]
            if cv:
                w = c
                for a in range(H):
                    w *= src[a] ** tm[a] if tm[a] else 1
                tot += w * cv
        if ok:
            assert tot == 0, f"termwise FAIL at {m}: {tot}"
            checked += 1
    return checked


def perm_op(op, s):
    """S6 action: permute lattice + theta indices by s (s[i] = image of i)."""
    def pm(t):
        out = [0] * H
        for i in range(H):
            out[s[i]] = t[i]
        return tuple(out)
    return {(pm(zm), pm(tm)): c for (zm, tm), c in op.items()}


def main():
    import random
    random.seed(20260717)
    t0 = time.time()
    # gate 0: exact ratio law c(k+e_I)/c(k) = (|k|+1)^2/(k_I+1)^2, 200 draws
    for _ in range(200):
        k = tuple(random.randrange(0, 9) for _ in range(H))
        i = random.randrange(H)
        ke = tuple(k[a] + (1 if a == i else 0) for a in range(H))
        assert cn(ke) / cn(k) == Fr((sum(k) + 1) ** 2, (k[i] + 1) ** 2)
    # gate 1: D_I printed form: th_I^2 - phi^I (sum th + 1)^2
    for i in range(H):
        want = {((0,) * H, tuple(2 * x for x in E[i])): 1}
        for a in range(H):
            for b in range(H):
                m = tuple(E[a][j] + E[b][j] for j in range(H))
                want[(E[i], m)] = want.get((E[i], m), 0) - 1
            want[(E[i], E[a])] = want.get((E[i], E[a]), 0) - 2
        want[(E[i], (0,) * H)] = -1
        got = box_op(E[i])
        assert got == {k: c for k, c in want.items() if c}, f"D_{i+1} form"
    print(f"[gate] ratio-law 200/200 + D_I printed-form 6/6 PASS "
          f"({time.time()-t0:.2f}s)", flush=True)
    # supply: all 0 < |l|_1 <= 2
    supply = [l for l in iproduct(range(-2, 3), repeat=H)
              if 0 < sum(abs(x) for x in l) <= 2]
    C = {}
    for n in iproduct(range(5), repeat=H):
        C[n] = cn(n)
    print(f"[c-table] 5^6 = {len(C)} exact c(k) cached "
          f"({time.time()-t0:.2f}s)", flush=True)
    ops = {}
    t1 = time.time()
    l0 = supply[0]
    ops[l0] = box_op(l0)
    n0 = termwise_gate(ops[l0], C)
    dt = time.time() - t1
    print(f"[pilot-probe] first op l={l0}: {n0} identities, {dt:.2f}s "
          f"-> {len(supply)} ops ~ {dt*len(supply)/60:.1f} min", flush=True)
    assert dt * len(supply) < 1800, "timed-pilot: projected past 30-min envelope"
    tot_checked = n0
    maxterms = len(ops[l0])
    for l in supply[1:]:
        ops[l] = box_op(l)
        tot_checked += termwise_gate(ops[l], C)
        maxterms = max(maxterms, len(ops[l]))
    # gate 3: S6 equivariance of the supply: op(s.l) == s.op(l), 3 perms x all l
    perms = [(1, 0, 2, 3, 4, 5), (5, 0, 1, 2, 3, 4), (2, 4, 0, 5, 1, 3)]
    neq = 0
    for s in perms:
        for l in supply:
            sl = [0] * H
            for i in range(H):
                sl[s[i]] = l[i]
            assert ops[tuple(sl)] == perm_op(ops[l], s), f"S6 fail {s} {l}"
            neq += 1
    # bank the supply (integer normal-ordered reps)
    import json
    out = {" ".join(map(str, l)): sorted(
        [[list(zm), list(tm), c] for (zm, tm), c in ops[l].items()])
        for l in supply}
    with open(sys.argv[1] if len(sys.argv) > 1 else
              "b1_boxops_supply.json", "w") as f:
        json.dump(out, f)
    print(f"[gate] {len(supply)} box ops (|l|_1<=2, orthant) termwise-verified"
          f" NBOX=4: {tot_checked} identities ALL PASS; S6-equivariance "
          f"{neq} checks PASS; max terms/op = {maxterms}", flush=True)
    print(f"[done] wall {time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()

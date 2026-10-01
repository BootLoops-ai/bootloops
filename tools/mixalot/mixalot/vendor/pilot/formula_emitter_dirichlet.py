#!/usr/bin/env python3
"""formula_emitter_dirichlet.py — Dirichlet extension of the harmonic-class
closed forms: 2-mixture of one k-state variable, INTEGER Dirichlet priors
Dir(beta) on the mixing weight, Dir(at) on theta, Dir(ar) on rho.

Structure: Z = C * sum_m W(m) * G(m),
  W(m) = sum_{|x|=m, x<=U} prod_v C(u_v,x_v) * (at_v)^rise_{x_v} * (ar_v)^rise_{u_v-x_v}
       — a convolution of k polynomial-times-indicator sequences, hence
       PIECEWISE POLYNOMIAL in m (breakpoints among subset sums of u_v+1);
  G(m) = Gamma(m+b0)Gamma(N-m+b1) / [Gamma(m+A)Gamma(N-m+B)],
       A=|at|, B=|ar| — ratio of rising factorials: polynomial numerator and/or
       SIMPLE-POLE denominator prod (m+i), prod (N-m+j)  (integer parameters
       => all poles simple => plain harmonic numbers suffice, no H^(2)).
Reduction: partial fractions + per-piece polynomial division (piece polynomials
recovered by exact interpolation, verified on EVERY point of the piece — a
per-instance proof), ranges of 1/(m+i) => H terms.  Emits
{'const': Fraction, ('H', n): Fraction};  Z = C * eval(formula).
Gates: == brute-force Dirichlet engine (lsx_direct.Z_xsum) on random small
instances; uniform (all-ones) case == uniform emitter.
"""

import itertools
import random
import os
import sys
from fractions import Fraction
from math import comb, factorial

# sibling imports resolve to the vendored copies beside this file
# (MIXALOT_PILOT_DIR overrides the default sibling path)
sys.path.insert(0, os.environ.get("MIXALOT_PILOT_DIR",
                                  os.path.dirname(os.path.abspath(__file__))))
from formula_emitter import (_poly_divide_linear, _poly_eval,  # noqa: E402
                             emit_formula, eval_formula)


def rising(a, n):
    r = 1
    for t in range(n):
        r *= a + t
    return r


def _poly_mul(p, q):
    """Exact product of two coefficient lists (ascending powers of m)."""
    out = [Fraction(0)] * (len(p) + len(q) - 1)
    for a, pa in enumerate(p):
        if pa:
            for b, qb in enumerate(q):
                out[a + b] += pa * qb
    return out


def W_sequence(U, at, ar):
    """Exact W(m) values, m=0..N (integer convolution)."""
    W = [1]
    for u, a_t, a_r in zip(U, at, ar):
        f = [comb(u, x) * rising(a_t, x) * rising(a_r, u - x)
             for x in range(u + 1)]
        new = [0] * (len(W) + u)
        for m, w in enumerate(W):
            if w:
                for x, fx in enumerate(f):
                    new[m + x] += w * fx
        W = new
    return W


def piecewise_fit(W, breakpoints, maxdeg):
    """Fit each piece [b_t, b_{t+1}) by exact polynomial interpolation of
    degree <= maxdeg; VERIFY on every integer point of the piece."""
    pieces = []
    bps = sorted(set(breakpoints))
    for lo, hi in zip(bps, bps[1:] + [len(W)]):
        if lo >= len(W):
            break
        hi = min(hi, len(W))
        pts = list(range(lo, hi))
        deg = min(maxdeg, len(pts) - 1)
        # Newton interpolation on first deg+1 points (exact Fractions)
        xs = pts[:deg + 1]
        coeffs_newton = []
        table = [Fraction(W[x]) for x in xs]
        for lvl in range(len(xs)):
            coeffs_newton.append(table[0])
            table = [(table[i + 1] - table[i]) / (xs[i + 1 + lvl] - xs[i])
                     for i in range(len(table) - 1)]
        # convert to monomial coefficients
        poly = [Fraction(0)] * (deg + 1)
        basis = [Fraction(1)]           # prod (m - xs[j]) ascending
        for lvl, cN in enumerate(coeffs_newton):
            for a, b in enumerate(basis):
                poly[a] += cN * b
            if lvl < deg:
                newb = [Fraction(0)] * (len(basis) + 1)
                for a, b in enumerate(basis):
                    newb[a + 1] += b
                    newb[a] -= b * xs[lvl]
                basis = newb
        # verify on ALL points of the piece (per-instance proof)
        for m in pts:
            assert _poly_eval(poly, m) == W[m], \
                f"piece [{lo},{hi}) fit fails at m={m} (need finer breakpoints)"
        pieces.append((lo, hi - 1, poly))
    return pieces


def emit_dirichlet(U, at, ar, beta):
    """Returns (formula, C) with Z = C * [formula evaluated with harmonic H]."""
    k = len(U)
    N = sum(U)
    A, B = sum(at), sum(ar)
    b0, b1 = beta
    # normalization: Gamma(|beta|)/prod Gamma(beta) * prod_i Gamma(|a_i|)/prod Gamma(a_ij)
    # (densities wrt Lebesgue PROBABILITY measure include the (t)!-style factors:
    #  Dir-vs-uniform ratio; matches lsx_direct.E_dir/D_dir conventions)
    C = Fraction(factorial(b0 + b1 - 1),
                 factorial(b0 - 1) * factorial(b1 - 1))
    # NOTE: W's rising(a,x) = Gamma(a+x)/Gamma(a) already carries the
    # 1/prod Gamma(alpha) of E_dir — do NOT divide again here (dividing
    # again is off by prod (alpha-1)! exactly).
    C *= Fraction(factorial(A - 1))
    C *= Fraction(factorial(B - 1))
    # binomials folded into W; leftover Gamma structure:
    # D_dir: Gamma(m+b0)Gamma(N-m+b1)/Gamma(N+b0+b1) ; E: 1/[Gamma(m+A)Gamma(N-m+B)]
    C /= factorial(N + b0 + b1 - 1)
    # G(m) = (m+b0)...(m+A-1)^{-1} type: represent numerator polys/denominators
    #   Gamma(m+b0)/Gamma(m+A) = 1/prod_{i=b0}^{A-1}(m+i)   if b0 <= A
    #   if b0 > A: polynomial prod_{i=A}^{b0-1}(m+i)   (numerator factor)
    # and mirrored in (N-m) for b1 vs B; each range below is empty on the
    # side that does not apply, so both orderings are handled per side.
    den_m = list(range(b0, A))          # pole factors (m+i)
    den_c = list(range(b1, B))          # pole factors (N-m+j)
    num_m = list(range(A, b0))          # polynomial factors (m+i)
    num_c = list(range(B, b1))          # polynomial factors (N-m+j)
    W = W_sequence(U, at, ar)
    maxdeg = sum(a_t + a_r - 2 for a_t, a_r in zip(at, ar)) + k - 1
    bps = {0}
    for T in itertools.chain.from_iterable(
            itertools.combinations(range(k), s) for s in range(1, k + 1)):
        bps.add(sum(U[v] + 1 for v in T))
    pieces = piecewise_fit(W, {b for b in bps if b <= N}, maxdeg)
    if num_m or num_c:
        # fold the polynomial part of G(m) into the pieces exactly; the
        # remaining den_m/den_c poles go through partial fractions as usual
        gnum = [Fraction(1)]
        for i in num_m:
            gnum = _poly_mul(gnum, [Fraction(i), Fraction(1)])
        for j in num_c:
            gnum = _poly_mul(gnum, [Fraction(N + j), Fraction(-1)])
        pieces = [(lo, hi, _poly_mul(poly, gnum)) for lo, hi, poly in pieces]

    # partial fractions over den_m x den_c
    Am = {i: Fraction(1) for i in den_m}
    for i in den_m:
        for ip in den_m:
            if ip != i:
                Am[i] /= (ip - i)
    Ac = {j: Fraction(1) for j in den_c}
    for j in den_c:
        for jp in den_c:
            if jp != j:
                Ac[j] /= (jp - j)

    out = {"const": Fraction(0)}

    def add_H(n_idx, co):
        if co:
            key = ("H", n_idx)
            out[key] = out.get(key, Fraction(0)) + co

    def S_of(i, mirror):
        """sum_m W(m)/(m+i) (mirror=False) or /(N-m+i) (mirror=True):
        rational part -> const; range parts -> H terms. Returns via out."""
        tot_const = Fraction(0)
        hlist = []
        for lo, hi, poly in pieces:
            if mirror:
                # substitute m -> N-m': piece [lo,hi] of m = [N-hi, N-lo] of m'
                p2 = poly[:]
                # poly(N - m') expand
                q = [Fraction(0)] * len(p2)
                for a, coef in enumerate(p2):
                    # (N - m')^a
                    for t in range(a + 1):
                        q[t] += coef * comb(a, t) * (Fraction(N) ** (a - t)) * ((-1) ** t)
                poly_use, lo_u, hi_u = q, N - hi, N - lo
            else:
                poly_use, lo_u, hi_u = poly, lo, hi
            s, R = _poly_divide_linear(poly_use, i)
            tot_const += sum(_poly_eval(s, m) for m in range(lo_u, hi_u + 1))
            if R:
                hlist.append((hi_u + i, R))
                hlist.append((lo_u + i - 1, -R))
        return tot_const, hlist

    # cases: both denominators nonempty (generic), or one/both empty
    if den_m and den_c:
        for i in den_m:
            for j in den_c:
                w = Am[i] * Ac[j] / (N + i + j)
                c1, h1 = S_of(i, mirror=False)
                c2, h2 = S_of(j, mirror=True)
                out["const"] += w * (c1 + c2)
                for n_idx, co in h1 + h2:
                    add_H(n_idx, w * co)
    elif den_m and not den_c:
        for i in den_m:
            c1, h1 = S_of(i, mirror=False)
            out["const"] += Am[i] * c1
            for n_idx, co in h1:
                add_H(n_idx, Am[i] * co)
    elif den_c and not den_m:
        for j in den_c:
            c2, h2 = S_of(j, mirror=True)
            out["const"] += Ac[j] * c2
            for n_idx, co in h2:
                add_H(n_idx, Ac[j] * co)
    else:
        # no poles at all: sum the (numerator-folded) pieces directly
        out["const"] += sum(_poly_eval(poly, m)
                            for lo, hi, poly in pieces
                            for m in range(lo, hi + 1))
    out = {kk: vv for kk, vv in out.items() if kk == "const" or vv != 0}
    return out, C


def main():
    from lsx_direct import Z_xsum, model_1var, Model
    rng = random.Random(20260709)
    ok = True

    def gate(U, at, ar, beta):
        nonlocal ok
        k = len(U)
        if k == 2:
            m = model_1var(tuple(U))
        else:
            cells = [(tuple(tuple(1 if j == v else 0 for j in range(k))
                            for _ in range(1)), U[v]) for v in range(k)]
            m = Model([1], [((tuple(1 if j == v else 0 for j in range(k)),), U[v])
                            for v in range(k)])
        want = Z_xsum(m, dirichlet=(tuple(beta), [tuple(at)], [tuple(ar)]))
        f, C = emit_dirichlet(U, at, ar, beta)
        got = C * eval_formula({"const": f["const"],
                                **{kk: vv for kk, vv in f.items() if kk != "const"}},
                               Fraction(1))
        good = got == want
        ok &= good
        print(f"  {'PASS' if good else 'FAIL'} U={U} at={at} ar={ar} beta={beta}"
              + ("" if good else f"  got {got} want {want}"))

    # uniform sanity: all-ones == uniform emitter
    for U in ([2, 3], [1, 4], [3, 2, 2]):
        k = len(U)
        f, C = emit_dirichlet(U, [1] * k, [1] * k, (1, 1))
        got = C * eval_formula(f, Fraction(1))
        want = eval_formula(*emit_formula(U))
        good = got == want
        ok &= good
        print(f"  {'PASS' if good else 'FAIL'} uniform-reduction U={U}")
    # random Dirichlet gates (k=2,3) — both orderings of beta vs (A,B)
    for k in (2, 3):
        for _ in range(4):
            U = [rng.randint(0, 6) for _ in range(k)]
            at = [rng.randint(1, 3) for _ in range(k)]
            ar = [rng.randint(1, 3) for _ in range(k)]
            beta = (rng.randint(1, 3), rng.randint(1, 3))
            gate(U, at, ar, beta)
    # polynomial-branch gates: b0 > A and/or b1 > B (numerator factors of G)
    for U, at, ar, beta in (
            ([3, 2], [1, 1], [1, 1], (5, 2)),           # b0 > A only
            ([2, 4], [1, 1], [1, 2], (4, 5)),           # both sides polynomial
            ([1, 3], [1, 2], [1, 1], (3, 3)),           # b0 = A edge, b1 > B
            ([4, 1], [2, 1], [1, 1], (1, 6)),           # b1 > B only
            ([2, 2, 1], [1, 1, 1], [1, 1, 1], (6, 2)),  # k=3, b0 > A
    ):
        gate(U, at, ar, beta)
    print("DIRICHLET EMITTER:", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()

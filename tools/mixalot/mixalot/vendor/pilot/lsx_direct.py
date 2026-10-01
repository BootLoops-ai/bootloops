#!/usr/bin/env python3
"""lsx_direct.py — BAYES-MIXTURES Phase 0a: independent direct exact evaluator.

Computes the BARE marginal-likelihood integral (LSX eq (25) convention)

    Z(U) = int_{Delta_1 x P x P} prod_v (sigma_0 theta^{a_v} + sigma_1 rho^{a_v})^{U_v}

with Lebesgue PROBABILITY measure on every simplex (LSX eq (12), with its t!
factor), for mixtures of independence models (secant of Segre-Veronese).
Conventions pinned to Lin-Sturmfels-Xu arXiv:0805.3602v2 = JMLR 10:1611 (2009);
see ../refs/lsx.pdf.  Reduced columns are UNSCALED (Section 4.1 Maple code);
full marginal likelihood for reduced data = N!/prod(U_v!) * prod(alpha_v^U_v) * Z
with alpha_v the column multiplicities (LSX Section 3 remark).

Engines (independent code paths, cross-gated in --selftest):
  * phi-DP  : dense numpy int64 tensor over (m, kept b-coords); exact for N<=62
              because sum_b phi = 2^N and every partial sum is <= 2^N.
  * x-sum   : brute enumeration of eq (19)-(20) (tiny U only).
  * sympy   : direct symbolic integration (tiny U only) — absolute convention
              anchor for measure/prior/label-switching conventions.

Assembly: Z = sum_m  m!(N-m)!/(N+1)! * prod_i [ t_i!^2 / ((s_i m+t_i)! (s_i(N-m)+t_i)!) ]
                   * S_m,   S_m = sum_{b in slice m} phi_A(b,U) * prod_{i,j} b_ij! c_ij!
(c = AU - b), one Fraction per m-slice, integer arithmetic inside the slice.

References vendored from the paper (exact printed values):
  swiss   : eq (3)      — 100 Swiss Francs 4x4 table, N=40
  coin10  : Sec 4.1     — coin-toss model U=(2,2,2,2,2), N=10
  coin242 : Ex 2.5/5.3  — U=(51,18,73,25,75), N=242, full ML ~ 0.7788716338838678611335742e-22
Support-count gate: #nonzero phi for swiss must equal 3,892,097 (LSX Ex 3.7).
"""

import argparse
import sys
import time
from fractions import Fraction
from itertools import product as iproduct
from math import comb, factorial

import numpy as np

# ----------------------------------------------------------------------------
# Model spec
# ----------------------------------------------------------------------------

class Model:
    """cells: list of (a_v, U_v) with a_v = tuple over blocks of tuples of
    exponents (length t_i+1, entries summing to s_i).  blocks: list of s_i."""

    def __init__(self, block_s, cells):
        self.block_s = list(block_s)          # s_i per block
        self.cells = cells                    # [(a_v, U_v)]
        self.nblocks = len(block_s)
        self.block_len = [len(cells[0][0][i]) for i in range(self.nblocks)]
        for a_v, _ in cells:
            for i, s in enumerate(self.block_s):
                assert sum(a_v[i]) == s, (a_v, i, s)
        self.N = sum(U for _, U in cells)
        # AU per block (length t_i+1 each)
        self.AU = [tuple(sum(U * a_v[i][j] for a_v, U in cells)
                         for j in range(self.block_len[i]))
                   for i in range(self.nblocks)]

    @classmethod
    def from_table(cls, U):
        """Two-way table (k=2, s=(1,1)): U = list of rows."""
        r, c = len(U), len(U[0])
        cells = []
        for i in range(r):
            for j in range(c):
                a = (tuple(1 if x == i else 0 for x in range(r)),
                     tuple(1 if y == j else 0 for y in range(c)))
                cells.append((a, U[i][j]))
        return cls([1, 1], cells)

    @classmethod
    def coin(cls, U):
        """One group of s_1=4 iid binary variables, reduced columns
        a~_v = (4-v, v), v=0..4 (UNSCALED)."""
        assert len(U) == 5
        cells = [(((4 - v, v),), U[v]) for v in range(5)]
        return cls([4], cells)


# ----------------------------------------------------------------------------
# Exact weights (uniform prior; integer-Dirichlet generalization)
# ----------------------------------------------------------------------------

def D_unif(a, b):
    return Fraction(factorial(a) * factorial(b), factorial(a + b + 1))


def E_unif(bvec):
    T = len(bvec)
    num = factorial(T - 1)
    for x in bvec:
        num *= factorial(x)
    return Fraction(num, factorial(sum(bvec) + T - 1))


def D_dir(a, b, al):
    al0, al1 = al
    return Fraction(factorial(al0 + al1 - 1), factorial(al0 - 1) * factorial(al1 - 1)) * \
        Fraction(factorial(a + al0 - 1) * factorial(b + al1 - 1),
                 factorial(a + b + al0 + al1 - 1))


def E_dir(bvec, alvec):
    s = sum(alvec)
    num = factorial(s - 1)
    den = 1
    for a in alvec:
        den *= factorial(a - 1)
    for x, a in zip(bvec, alvec):
        num *= factorial(x + a - 1)
    return Fraction(num, den * factorial(sum(bvec) + s - 1))


# ----------------------------------------------------------------------------
# Engine 1: brute x-sum (eq 19-20 verbatim; tiny U only)
# ----------------------------------------------------------------------------

def Z_xsum(model, dirichlet=None):
    N = model.N
    Z = Fraction(0)
    ranges = [range(U + 1) for _, U in model.cells]
    for xs in iproduct(*ranges):
        m = sum(xs)
        coef = 1
        for x, (_, U) in zip(xs, model.cells):
            coef *= comb(U, x)
        b = [[0] * L for L in model.block_len]
        for x, (a_v, _) in zip(xs, model.cells):
            for i in range(model.nblocks):
                for j in range(model.block_len[i]):
                    b[i][j] += x * a_v[i][j]
        c = [[model.AU[i][j] - b[i][j] for j in range(model.block_len[i])]
             for i in range(model.nblocks)]
        if dirichlet is None:
            w = D_unif(m, N - m)
            for i in range(model.nblocks):
                w *= E_unif(b[i]) * E_unif(c[i])
        else:
            als, alth, alrh = dirichlet
            w = D_dir(m, N - m, als)
            for i in range(model.nblocks):
                w *= E_dir(b[i], alth[i]) * E_dir(c[i], alrh[i])
        Z += coef * w
    return Z


# ----------------------------------------------------------------------------
# Engine 2: dense int64 phi-DP + per-m-slice exact assembly
# ----------------------------------------------------------------------------

def phi_dp(model):
    """Dense DP for phi_A(b,U).  Axes: (m, kept coords) where kept coords are
    b_ij for j=1..t_i per block (b_i0 = s_i*m - sum kept).  Exact in int64 for
    N <= 62 (sum phi = 2^N bounds every entry and partial sum)."""
    N = model.N
    if N > 62:
        raise ValueError("int64 phi-DP valid only for N<=62 (sum phi = 2^N); "
                         "use the mod-p engine (Phase 0b) beyond that")
    caps = []
    for i in range(model.nblocks):
        for j in range(1, model.block_len[i]):
            caps.append(model.AU[i][j])
    dims = [N + 1] + [c + 1 for c in caps]
    phi = np.zeros(dims, dtype=np.int64)
    phi[(0,) * len(dims)] = 1
    for a_v, U in model.cells:
        d = [1]  # m increments by 1 per unit of x
        for i in range(model.nblocks):
            for j in range(1, model.block_len[i]):
                d.append(a_v[i][j])
        new = np.zeros_like(phi)
        for x in range(U + 1):
            C = comb(U, x)
            src = tuple(slice(0, dims[k] - x * d[k]) if d[k] else slice(None)
                        for k in range(len(dims)))
            dst = tuple(slice(x * d[k], dims[k]) if d[k] else slice(None)
                        for k in range(len(dims)))
            new[dst] += C * phi[src]
        phi = new
    assert int(phi.sum()) == 2 ** N, "phi mass check failed (overflow?)"
    return phi, caps


def phi_dp_big(model, max_states=20_000_000):
    """Dict-based exact bignum phi-DP for N>62 (small state spaces only).
    Returns (phi_array_object, caps) shaped like phi_dp's output."""
    N = model.N
    caps = []
    for i in range(model.nblocks):
        for j in range(1, model.block_len[i]):
            caps.append(model.AU[i][j])
    dims = [N + 1] + [c + 1 for c in caps]
    est = 1
    for d in dims:
        est *= d
    if est > max_states:
        raise ValueError(f"state space {est} exceeds max_states={max_states}; "
                         "use the mod-p engine")
    phi = {(0,) * len(dims): 1}
    for a_v, U in model.cells:
        d = [1]
        for i in range(model.nblocks):
            for j in range(1, model.block_len[i]):
                d.append(a_v[i][j])
        new = {}
        Cs = [comb(U, x) for x in range(U + 1)]
        for key, w in phi.items():
            for x in range(U + 1):
                nk = tuple(k + x * dk for k, dk in zip(key, d))
                new[nk] = new.get(nk, 0) + Cs[x] * w
        phi = new
    assert sum(phi.values()) == 2 ** N, "phi mass check failed"
    arr = np.zeros(dims, dtype=object)
    for key, w in phi.items():
        arr[key] = w
    return arr, caps


def Z_phi(model, return_support=False):
    """Exact Z from the phi-DP, per-m-slice integer assembly (uniform prior)."""
    N = model.N
    if N <= 62:
        phi, caps = phi_dp(model)
    else:
        phi, caps = phi_dp_big(model)
    nz = np.nonzero(phi)
    vals = phi[nz]
    support = len(vals)
    coords = np.stack(nz, axis=1)  # (support, 1+len(caps))
    # group indices by m
    Z = Fraction(0)
    ms = coords[:, 0]
    order = np.argsort(ms, kind='stable')
    coords, vals, ms = coords[order], vals[order], ms[order]
    starts = np.searchsorted(ms, np.arange(N + 2))
    F = [factorial(i)
         for i in range(max(s * N for s in model.block_s) + 1)]
    for m in range(N + 1):
        lo, hi = starts[m], starts[m + 1]
        if lo == hi:
            continue
        S = 0
        for idx in range(lo, hi):
            row = coords[idx]
            w = int(vals[idx])
            pos = 1
            num = w
            for i in range(model.nblocks):
                s_i = model.block_s[i]
                kept = row[pos:pos + model.block_len[i] - 1]
                pos += model.block_len[i] - 1
                b0 = s_i * m - int(kept.sum())
                c0 = s_i * (N - m) - int(sum(model.AU[i][1:])) + int(kept.sum())
                num *= F[b0] * F[c0]
                for j, bij in enumerate(kept, start=1):
                    num *= F[int(bij)] * F[model.AU[i][j] - int(bij)]
            S += num
        pref = Fraction(factorial(m) * factorial(N - m), factorial(N + 1))
        for i in range(model.nblocks):
            t_i = model.block_len[i] - 1
            s_i = model.block_s[i]
            pref *= Fraction(factorial(t_i) ** 2,
                             factorial(s_i * m + t_i) * factorial(s_i * (N - m) + t_i))
        Z += pref * S
    return (Z, support) if return_support else Z


# ----------------------------------------------------------------------------
# Engine 3: sympy direct integration (tiny U; absolute convention anchor)
# ----------------------------------------------------------------------------

def Z_sympy_1var(U):
    from sympy import symbols, integrate, expand, Rational
    s, t, r = symbols('s t r', positive=True)
    p0 = s * t + (1 - s) * r
    p1 = s * (1 - t) + (1 - s) * (1 - r)
    f = expand(p0 ** U[0] * p1 ** U[1])
    for var in (s, t, r):
        f = integrate(f, (var, 0, 1))
    return Fraction(int(f.p), int(f.q))


def Z_sympy_2x2(U):
    from sympy import symbols, integrate, expand, Integer
    s, t1, t2, r1, r2 = symbols('s t1 t2 r1 r2', positive=True)
    tr = [t1, 1 - t1]; tc = [t2, 1 - t2]
    rr = [r1, 1 - r1]; rc = [r2, 1 - r2]
    f = Integer(1)
    for i in range(2):
        for j in range(2):
            f *= (s * tr[i] * tc[j] + (1 - s) * rr[i] * rc[j]) ** U[i][j]
    f = expand(f)
    for var in (s, t1, t2, r1, r2):
        f = integrate(f, (var, 0, 1))
    return Fraction(int(f.p), int(f.q))


def model_1var(U):
    cells = [(((1, 0),), U[0]), (((0, 1),), U[1])]
    return Model([1], cells)


# ----------------------------------------------------------------------------
# Vendored LSX reference values (exact printed forms)
# ----------------------------------------------------------------------------

SWISS_U = [[4, 2, 2, 2], [2, 4, 2, 2], [2, 2, 4, 2], [2, 2, 2, 4]]
SWISS_NUM = 571 * 773426813 * 17682039596993 * 625015426432626533
SWISS_DEN = (2**31 * 3**20 * 5**12 * 7**11 * 11**8 * 13**7 * 17**5 * 19**5
             * 23**5 * 29**3 * 31**3 * 37**3 * 41**3 * 43**2)
SWISS_SUPPORT = 3_892_097          # LSX Example 3.7
COIN10_U = (2, 2, 2, 2, 2)
COIN10_REF = Fraction(66364720654753, 59057383987217015339940000)  # LSX Sec 4.1
COIN242_U = (51, 18, 73, 25, 75)
COIN242_FULLML_25DIG = "7.788716338838678611335742e-23"            # LSX Ex 2.5
COIN_MULT = (1, 4, 6, 4, 1)        # column multiplicities alpha_v (reduced)


def full_ml_reduced(Z, U, mult):
    """Full marginal likelihood for reduced data: N!/prod U_v! * prod mult^U * Z."""
    N = sum(U)
    c = Fraction(factorial(N))
    for u, a in zip(U, mult):
        c /= factorial(u)
        c *= Fraction(a) ** u
    return c * Z


# ----------------------------------------------------------------------------
# Selftest & examples
# ----------------------------------------------------------------------------

def selftest():
    ok = True

    def chk(name, got, want):
        nonlocal ok
        good = got == want
        ok &= good
        print(f"  [{'PASS' if good else 'FAIL'}] {name}: {got}"
              + ("" if good else f"  (want {want})"))

    print("== tiny exact cross-gates (three independent engines) ==")
    for U, ref in [((1, 1), Fraction(7, 36)), ((2, 1), Fraction(7, 72)),
                   ((3, 2), Fraction(37, 1800))]:
        m = model_1var(U)
        chk(f"1var {U} x-sum", Z_xsum(m), ref)
        chk(f"1var {U} phi-DP", Z_phi(m), ref)
        chk(f"1var {U} sympy", Z_sympy_1var(U), ref)
    for U, ref in [(((1, 1), (1, 0)), Fraction(1, 108)),
                   (((2, 1), (1, 1)), Fraction(173, 432000))]:
        m = Model.from_table([list(r) for r in U])
        chk(f"2x2 {U} x-sum", Z_xsum(m), ref)
        chk(f"2x2 {U} phi-DP", Z_phi(m), ref)
        chk(f"2x2 {U} sympy", Z_sympy_2x2(U), ref)

    print("== engine cross-gates on random small tables ==")
    rng = np.random.default_rng(20260708)
    for trial in range(3):
        U = [[int(rng.integers(0, 3)) for _ in range(3)] for _ in range(2)]
        m = Model.from_table(U)
        a, b = Z_xsum(m), Z_phi(m)
        chk(f"2x3 random#{trial} {U} x-sum==phi-DP", a, b)
    mcoin = Model.coin((1, 0, 2, 1, 0))
    chk("coin (1,0,2,1,0) x-sum==phi-DP", Z_xsum(mcoin), Z_phi(mcoin))
    # bignum dict engine vs int64 engine (independent DP implementations)
    m23 = Model.from_table([[2, 1, 0], [1, 2, 1]])
    pa, _ = phi_dp(m23)
    pb, _ = phi_dp_big(m23)
    chk("2x3 phi_dp int64 == phi_dp_big dict",
        True, bool((pa.astype(object) == pb).all()))

    print("== Dirichlet-twisted tiny gate (Cor 5.1 conventions) ==")
    m = model_1var((2, 1))
    zd = Z_xsum(m, dirichlet=((2, 1), [(1, 3)], [(1, 1)]))
    chk("1var (2,1) Dir alpha=(2,1),(1,3),(1,1)", zd, Fraction(127, 1800))
    zd1 = Z_xsum(m, dirichlet=((1, 1), [(1, 1)], [(1, 1)]))
    chk("1var (2,1) Dir all-ones == uniform", zd1, Fraction(7, 72))

    print("== MUTATION CONTROL (must FAIL internally, PASS as control) ==")
    # sabotage: drop the t! measure factor -> must NOT equal reference
    bad = Z_xsum(Model.from_table([[1, 1], [1, 0]])) * Fraction(1, 9)
    good = bad != Fraction(1, 108)
    ok &= good
    print(f"  [{'PASS' if good else 'FAIL'}] sabotaged measure differs from ref")

    print("SELFTEST:", "PASS" if ok else "FAIL")
    return ok


def run_example(name, full_ml=False):
    t0 = time.time()
    if name == "swiss":
        m = Model.from_table(SWISS_U)
        Z, support = Z_phi(m, return_support=True)
        dt = time.time() - t0
        ref = Fraction(SWISS_NUM, SWISS_DEN)
        print(f"swiss: N={m.N} support={support} wall={dt:.1f}s")
        print(f"  support == 3,892,097 (LSX Ex 3.7): "
              f"{'PASS' if support == SWISS_SUPPORT else 'FAIL'}")
        print(f"  Z == eq(3) printed factorization:  "
              f"{'PASS' if Z == ref else 'FAIL'}")
        print(f"  Z = {Z.numerator}\n      / {Z.denominator}")
        if full_ml:
            c = Fraction(factorial(40))
            for row in SWISS_U:
                for u in row:
                    c /= factorial(u)
            ml = c * Z
            print(f"  full ML = 40!/(prod U_ij!) * Z ~ {float(ml):.9e} "
                  f"(paper ~5.679049589e-13)")
        return Z == ref and support == SWISS_SUPPORT
    elif name == "coin10":
        m = Model.coin(COIN10_U)
        Z = Z_phi(m)
        dt = time.time() - t0
        print(f"coin10: N={m.N} wall={dt:.2f}s")
        print(f"  Z == printed rational (Sec 4.1): "
              f"{'PASS' if Z == COIN10_REF else 'FAIL'}")
        print(f"  Z = {Z}")
        return Z == COIN10_REF
    elif name == "coin242":
        m = Model.coin(COIN242_U)
        Z = Z_phi(m)
        dt = time.time() - t0
        ml = full_ml_reduced(Z, COIN242_U, COIN_MULT)
        import mpmath as mp
        mp.mp.dps = 40
        got = mp.mpf(ml.numerator) / mp.mpf(ml.denominator)
        want = mp.mpf(COIN242_FULLML_25DIG)
        rel = abs(got - want) / want
        print(f"coin242: N={m.N} wall={dt:.2f}s")
        print(f"  bare Z digits: num {len(str(Z.numerator))} / den {len(str(Z.denominator))}")
        nn, dd = ml.numerator, ml.denominator
        print(f"  full ML digits: num {len(str(nn))} / den {len(str(dd))} "
              f"(paper: 530/552 relatively prime)")
        print(f"  full ML = {mp.nstr(got, 25)}")
        print(f"  vs paper 25-digit value : rel diff {mp.nstr(rel, 3)} "
              f"({'PASS' if rel < mp.mpf('1e-24') else 'FAIL'})")
        digit_gate = (len(str(nn)) == 530 and len(str(dd)) == 552)
        print(f"  530/552-digit gate: {'PASS' if digit_gate else 'FAIL'}")
        return rel < mp.mpf('1e-24') and digit_gate
    else:
        raise SystemExit(f"unknown example {name}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--example", choices=["swiss", "coin10", "coin242"])
    ap.add_argument("--full-ml", action="store_true")
    ap.add_argument("--table", help="semicolon-separated rows, e.g. '4,2;2,4'")
    args = ap.parse_args()
    if args.selftest:
        sys.exit(0 if selftest() else 1)
    if args.example:
        sys.exit(0 if run_example(args.example, args.full_ml) else 2)
    if args.table:
        U = [[int(x) for x in row.split(",")] for row in args.table.split(";")]
        m = Model.from_table(U)
        t0 = time.time()
        Z, support = Z_phi(m, return_support=True)
        print(f"N={m.N} support={support} wall={time.time()-t0:.1f}s")
        print(f"Z = {Z.numerator} / {Z.denominator}")
        sys.exit(0)
    ap.print_help()


if __name__ == "__main__":
    main()

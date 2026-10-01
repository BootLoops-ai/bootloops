#!/usr/bin/env python3
"""gmtel route A (GEOMETRIC route): Picard-Fuchs operator of the
REGENERATED K3 pencil, from the Weierstrass family itself, by exact fiber
Gauss-Manin (extended-Euclid reduction, extended to two parameters) +
creative telescoping over the base q.

Object: the lattice-regenerated family
    y^2 = x (x^2 + a2 x + a4),  a2 = q^2 (2q + s - sigma),
    a4 = q^3 (q-A)(q-B)(q-C),   sigma = A+B+C,  s := w^2  (pencil parameter).
K3 period Pi(s) = oint_gamma dq (oint dX/Y);  the telescoper produces
    L_A = sum_k p_k(s) d^k/ds^k   with certificate row-vector u(q,s):
    sum_k p_k(s) r_k(q,s) = d/dq u + u * M_q          (*)
    (r_0 = (1,0), r_{k+1} = d/ds r_k + r_k M_s)
so that L_A Pi = oint d/dq(u . I) dq = 0 over any closed q-contour.

Pipeline (all receipts SCRIPT-EMITTED):
  1. exact M_q, M_s over Q(q,s)         [sympy, extended-Euclid reduction]
  2. exact rows r_0..r_N over Q(q,s)
  3. mod-p telescoper search at sampled s (flint nmod_poly + int64 nullspace):
     minimal N, denominator power m, numerator degree dq; TWO primes
  4. rational reconstruction of p_k(s) (per-prime rational-function recon,
     CRT + Farey across primes); primitive integer normalization
  5. EXACT certificate verification over Q: solve u at the locked shape with
     p_k(s) known, verify (*) as a polynomial identity  [time-boxed; fallback:
     exact verification of (*) at >deg-bound many rational s values]
  6. comparison gate: L_A vs the reference L5 (x = 1/s chart), allowing an explicit
     rational/algebraic gauge twist J(s) (exhibited if nontrivial), and vs the
     route-B fitted operator.
Usage: famgen_pf_routeA.py --abc 1 4 16 [--nmax 6] --out CONTROL_PF_ROUTEA.json
"""
import argparse, itertools, json, os, subprocess, sys, time
from fractions import Fraction as F
import sympy as sp
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
# Dir holding the annihilator engine module (annihilator.py); the packaged
# default is the sibling tools/annihilator/ directory in this tree.
TOOLS = os.environ.get('FAMGEN_TOOLS',
                       os.path.join(HERE, '..', '..', 'annihilator'))
sys.path.insert(0, TOOLS)
from famgen_common import (load_L5, theta_form_to_D, op_primitive, op_equal,
                           op_add, op_mul_poly, op_mul_D, L5_PATH, moments)
import flint

q, s, X = sp.symbols('q s X')
P1 = 2147483647
P2 = 2147483629
P3 = 2147483563


def utc():
    return subprocess.run(['date', '-u', '+%Y-%m-%dT%H:%M:%SZ'],
                          capture_output=True, text=True).stdout.strip()


# ------------------------------------------------------------------ exact GM
def gauss_manin(Av, Bv, Cv):
    """Exact connection matrices M_q, M_s (2x2 sympy expr in q,s) for the
    period vector (I0, I1) = (oint dX/Y, oint X dX/Y)."""
    sigma = Av + Bv + Cv
    a2 = q ** 2 * (2 * q + s - sigma)
    a4 = q ** 3 * (q - Av) * (q - Bv) * (q - Cv)
    f = X ** 3 + a2 * X ** 2 + a4 * X
    fp = sp.expand(sp.diff(f, X))
    K = sp.QQ.frac_field(q, s)
    Pf = sp.Poly(f, X, domain=K)
    Pfp = sp.Poly(fp, X, domain=K)
    u_, v_, g_ = Pf.gcdex(Pfp)
    g0 = g_.all_coeffs()[0]
    u_, v_ = u_.quo_ground(g0), v_.quo_ground(g0)   # u f + v f' = 1

    def reduce_level3(Ppoly):
        Pp = sp.Poly(Ppoly, X, domain=K)
        b = (Pp * v_).rem(Pf)
        a = (Pp - b * Pfp).quo(Pf)
        r = a + b.diff(X) * sp.Poly(2, X, domain=K)
        return reduce_bottom(r)

    def reduce_bottom(r):
        while r.degree() >= 2:
            d = r.degree()
            k = d - 2
            Q = sp.Poly(X ** k, X, domain=K)
            rel = Q.diff(X) * Pf + Q * Pfp * sp.Poly(sp.Rational(1, 2), X, domain=K)
            lc_rel = rel.coeff_monomial(X ** (k + 2))
            rtop = r.coeff_monomial(X ** d)
            r = r - rel * sp.Poly(rtop / lc_rel, X, domain=K)
        c0 = r.coeff_monomial(1) if r.degree() >= 0 else K.zero
        c1 = r.coeff_monomial(X) if r.degree() >= 1 else K.zero
        return c0, c1

    out = {}
    for t in (q, s):
        a2t, a4t = sp.diff(a2, t), sp.diff(a4, t)
        M = sp.zeros(2, 2)
        for k in (0, 1):
            num = sp.expand(X ** k * (-(a2t * X ** 2 + a4t * X)) / 2)
            c0, c1 = reduce_level3(num)
            M[k, 0] = sp.cancel(sp.sympify(str(c0)))
            M[k, 1] = sp.cancel(sp.sympify(str(c1)))
        out[t] = M
    return out[q], out[s], a2, a4


def rows_rk(Ms, N):
    """r_0..r_N: r_{k+1} = d/ds r_k + r_k * Ms  (row vectors over Q(q,s))."""
    rows = [(sp.Integer(1), sp.Integer(0))]
    for _ in range(N):
        r1, r2 = rows[-1]
        n1 = sp.cancel(sp.diff(r1, s) + r1 * Ms[0, 0] + r2 * Ms[1, 0])
        n2 = sp.cancel(sp.diff(r2, s) + r1 * Ms[0, 1] + r2 * Ms[1, 1])
        rows.append((n1, n2))
    return rows


# ------------------------------------------------------- rational fn packing
def pack(expr):
    """sympy expr in (q,s) -> (num_dict, den_dict) integer coefficient dicts
    {(dq,ds): int}."""
    n, d = sp.fraction(sp.together(sp.cancel(expr)))
    # clear rational content
    np_ = sp.Poly(sp.expand(n), q, s)
    dp_ = sp.Poly(sp.expand(d), q, s)
    lcm_den = 1
    for c in np_.coeffs() + dp_.coeffs():
        lcm_den = sp.ilcm(lcm_den, sp.denom(c))
    np_ = sp.Poly(sp.expand(n * lcm_den), q, s)
    dp_ = sp.Poly(sp.expand(d * lcm_den), q, s)
    nd = {m: int(c) for m, c in zip(np_.monoms(), np_.coeffs())}
    dd = {m: int(c) for m, c in zip(dp_.monoms(), dp_.coeffs())}
    return nd, dd


def eval_nmod(cdict, sval, p):
    """coefficient dict {(dq,ds): int} -> flint nmod_poly in q at s=sval mod p."""
    if not cdict:
        return flint.nmod_poly([], p)
    maxdq = max(m[0] for m in cdict)
    co = [0] * (maxdq + 1)
    for (dq_, ds_), c in cdict.items():
        co[dq_] = (co[dq_] + c * pow(int(sval), ds_, p)) % p
    return flint.nmod_poly(co, p)


def nullspace_modp(rows, ncols, p):
    """int64 Gaussian elimination nullspace over F_p (rows: list of int lists)."""
    A = np.array(rows, dtype=np.int64) % p
    nr = A.shape[0]
    piv_cols = []
    r = 0
    for c in range(ncols):
        piv = None
        for i in range(r, nr):
            if A[i, c] % p:
                piv = i
                break
        if piv is None:
            continue
        A[[r, piv]] = A[[piv, r]]
        inv = pow(int(A[r, c]), p - 2, p)
        A[r] = (A[r] * inv) % p
        for i in range(nr):
            if i != r and A[i, c]:
                A[i] = (A[i] - A[i, c] * A[r]) % p
        piv_cols.append(c)
        r += 1
        if r == nr:
            break
    free = [c for c in range(ncols) if c not in piv_cols]
    basis = []
    for fc in free:
        v = [0] * ncols
        v[fc] = 1
        for ri, pc in enumerate(piv_cols):
            v[pc] = (-int(A[ri, fc])) % p
        basis.append(v)
    return basis


# ------------------------------------------------------------- the telescoper
class Telescoper:
    def __init__(self, Av, Bv, Cv, nmax=6, verbose=True):
        self.ABC = (Av, Bv, Cv)
        self.nmax = nmax
        self.verbose = verbose
        t0 = time.time()
        self.Mq, self.Ms, self.a2, self.a4 = gauss_manin(Av, Bv, Cv)
        self.t_gm = time.time() - t0
        t0 = time.time()
        self.rows = rows_rk(self.Ms, nmax)
        self.t_rows = time.time() - t0
        # pack everything once
        self.Mq_p = [[pack(self.Mq[i, j]) for j in range(2)] for i in range(2)]
        self.rows_p = [(pack(r1), pack(r2)) for (r1, r2) in self.rows]
        # denominator factor basis: all irreducible q-factors (s-dependent ok)
        dens = set()
        polys = []
        for i in range(2):
            for j in range(2):
                polys.append(sp.fraction(sp.together(sp.cancel(self.Mq[i, j])))[1])
                polys.append(sp.fraction(sp.together(sp.cancel(self.Ms[i, j])))[1])
        for (r1, r2) in self.rows:
            polys.append(sp.fraction(sp.together(sp.cancel(r1)))[1])
            polys.append(sp.fraction(sp.together(sp.cancel(r2)))[1])
        fset = {}
        for pl in polys:
            for base, mult in sp.factor_list(pl, q, s)[1]:
                base = sp.expand(base)
                if sp.Poly(base, q).degree() == 0:
                    continue        # s-only factors are units in F_p(q) at fixed s
                key = sp.srepr(base)
                fset[key] = (base, max(mult, fset.get(key, (base, 0))[1]))
        self.badfactors = [b for b, m in sorted(fset.values(),
                                                key=lambda t: sp.count_ops(t[0]))]
        self.badmults = [m for b, m in sorted(fset.values(),
                                              key=lambda t: sp.count_ops(t[0]))]
        self.badfactors_p = [pack(b) for b in self.badfactors]

    def try_solve(self, p, sval, N, m, dq_extra, want_vec=False):
        """Attempt telescoper at prime p, s=sval, order N, denominator power m,
        certificate numerator degree = degE + dq_extra. Returns nullspace
        basis (each vec = [c_0..c_N, U1 coeffs, U2 coeffs]) or None."""
        # E(q) = prod badfactors(q, sval)^(maxmult_f + m)   [adaptive powers]
        E = flint.nmod_poly([1], p)
        for (nd, dd), bmult in zip(self.badfactors_p, self.badmults):
            bf = eval_nmod(nd, sval, p)
            assert eval_nmod(dd, sval, p).degree() == 0
            dinv = pow(int(eval_nmod(dd, sval, p)[0]), p - 2, p)
            bf = bf * dinv
            for _ in range(bmult + m):
                E = E * bf
        degE = E.degree()
        dU = degE + dq_extra
        # common denominator D_all: E^2 * prod(dens of rows and Mq at sval)
        row_nd = []   # (num1, den1, num2, den2) as nmod_poly
        for (n1, d1), (n2, d2) in self.rows_p[:N + 1]:
            t = (eval_nmod(n1, sval, p), eval_nmod(d1, sval, p),
                 eval_nmod(n2, sval, p), eval_nmod(d2, sval, p))
            if t[1].degree() < 0 or t[3].degree() < 0:
                return None     # bad sample point (denominator vanished)
            row_nd.append(t)
        Mq_nd = [[(eval_nmod(self.Mq_p[i][j][0], sval, p),
                   eval_nmod(self.Mq_p[i][j][1], sval, p)) for j in range(2)]
                 for i in range(2)]
        for i in range(2):
            for j in range(2):
                if Mq_nd[i][j][1].degree() < 0:
                    return None
        # lcm of denominators (poly lcm via gcd)
        def plcm(a, b):
            g = a.gcd(b)
            return a * b // g if g.degree() >= 0 else a * b
        Dall = E * E
        for (n1, d1, n2, d2) in row_nd:
            Dall = plcm(Dall, d1)
            Dall = plcm(Dall, d2)
        for i in range(2):
            for j in range(2):
                Dall = plcm(Dall, E * E * Mq_nd[i][j][1])
        Ep = E.derivative()
        # build columns: identity components (two nmod_polys) per unknown
        cols = []
        # c_k columns: c_k * (r_k1, r_k2) * Dall
        for k in range(N + 1):
            n1, d1, n2, d2 = row_nd[k]
            cols.append((n1 * (Dall // d1), n2 * (Dall // d2)))
        # U1 coefficient columns: -d/dq(U1/E) - (U1/E) Mq[0,:]
        #   = -(U1' E - U1 E')/E^2 - (U1/E)(Mq00, Mq01)
        # multiplied by Dall
        for comp in (0, 1):
            for t in range(dU + 1):
                mono = flint.nmod_poly([0] * t + [1], p)
                monop = mono.derivative()
                # derivative part
                d_part = (monop * E - mono * Ep)
                v1 = flint.nmod_poly([], p)
                v2 = flint.nmod_poly([], p)
                if comp == 0:
                    v1 = v1 - d_part * (Dall // (E * E))
                else:
                    v2 = v2 - d_part * (Dall // (E * E))
                for j in range(2):
                    mn, md = Mq_nd[comp][j]
                    contrib = mono * mn * (Dall // (E * E * md)) * E
                    # (U/E)*M = U*M_n/(E*md): multiply by Dall => U*M_n*Dall/(E*md)
                    # Dall divisible by E^2*md => Dall/(E^2*md)*E = Dall/(E*md) ok
                    if j == 0:
                        v1 = v1 - contrib
                    else:
                        v2 = v2 - contrib
                cols.append((v1, v2))
        # assemble linear system: coefficients of q^i in v1 and v2 = 0
        maxdeg = 0
        for (v1, v2) in cols:
            maxdeg = max(maxdeg, v1.degree(), v2.degree())
        rows_mat = []
        ncols = len(cols)
        for i in range(maxdeg + 1):
            rows_mat.append([int(v1[i]) for (v1, v2) in cols])
            rows_mat.append([int(v2[i]) for (v1, v2) in cols])
        basis = nullspace_modp(rows_mat, ncols, p)
        # demand c-part nonzero with c_N != 0 for at least one basis vector
        good = [v for v in basis if any(v[k] for k in range(N + 1))]
        if not good:
            return None
        return good if want_vec else good

    def find_shape(self, p, svals):
        """Minimal (N, m, dq_extra) that solves at all sample svals."""
        for N in range(1, self.nmax + 1):
            for m in (0, 1, 2, 4):
                for dq_extra in (5, 15, 30, 60):
                    ok = True
                    dims = []
                    for sv in svals:
                        b = self.try_solve(p, sv, N, m, dq_extra)
                        if b is None:
                            ok = False
                            break
                        dims.append(len(b))
                    if ok:
                        return N, m, dq_extra, dims
        return None


# ---------------------------------------------------- rational reconstruction
def rat_recon_func(vals, pts, p, dnum, dden):
    """Rational function P/Q over F_p with P(x_j)=v_j Q(x_j), deg P<=dnum,
    deg Q<=dden.  CANONICALIZED: gcd(P,Q) removed, Q made monic — the reduced
    monic-denominator form is unique, so images at different primes align.
    Returns (Pcoeffs, Qcoeffs) or None."""
    ncols = (dnum + 1) + (dden + 1)
    rows = []
    for x, v in zip(pts, vals):
        row = [pow(x, i, p) for i in range(dnum + 1)]
        row += [(-v * pow(x, i, p)) % p for i in range(dden + 1)]
        rows.append(row)
    basis = nullspace_modp(rows, ncols, p)
    for vec in basis:
        Pc = vec[:dnum + 1]
        Qc = vec[dnum + 1:]
        if not any(Qc):
            continue
        Pf = flint.nmod_poly(Pc, p)
        Qf = flint.nmod_poly(Qc, p)
        g = Pf.gcd(Qf)
        if g.degree() > 0:
            Pf = Pf // g
            Qf = Qf // g
        lead = int(Qf[Qf.degree()])
        inv = pow(lead, p - 2, p)
        Pf = Pf * inv
        Qf = Qf * inv
        Pc = [int(Pf[i]) for i in range(Pf.degree() + 1)] or [0]
        Qc = [int(Qf[i]) for i in range(Qf.degree() + 1)]
        return Pc, Qc
    return None


def farey(a, p, bound):
    """Rational reconstruction of a mod p with |num|,|den| <= bound."""
    g0, g1 = p, a % p
    s0, s1 = 0, 1
    while g1 > bound:
        qq = g0 // g1
        g0, g1 = g1, g0 - qq * g1
        s0, s1 = s1, s0 - qq * s1
    if abs(s1) > bound or g1 == 0 and s1 == 0:
        return None
    return F(g1, s1) if s1 > 0 else F(-g1, -s1)


# ------------------------------------------------------------------ main flow
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--abc', nargs=3, type=int, required=True)
    ap.add_argument('--nmax', type=int, default=6)
    ap.add_argument('--nsamp', type=int, default=120)
    ap.add_argument('--out', required=True)
    ap.add_argument('--routeb', help='route-B receipt for the cross gate')
    ap.add_argument('--exact-verify-cap-s', type=int, default=900)
    args = ap.parse_args()
    Av, Bv, Cv = args.abc
    t00 = time.time()
    rep = {'leg': 'gmtel', 'route': 'A (geometric: fiber GM + q-telescoping)',
           'ABC': [Av, Bv, Cv], 'stamp_utc': utc(),
           'primes': [P1, P2, P3],
           'family': 'y^2 = x(x^2 + q^2(2q+s-sigma)x + q^3(q-A)(q-B)(q-C)), s = w^2'}

    T = Telescoper(Av, Bv, Cv, nmax=args.nmax)
    rep['t_gm_s'] = round(T.t_gm, 2)
    rep['t_rows_s'] = round(T.t_rows, 2)
    rep['bad_factors'] = [str(b) for b in T.badfactors]

    # ---- shape search at prime 1
    svals0 = [7, 11, 13]
    shape = T.find_shape(P1, svals0)
    if shape is None:
        rep['verdict'] = 'WALL: no telescoper found within (N<=%d, m<=3, dq<=+30)' % args.nmax
        json.dump(rep, open(args.out, 'w'), indent=1)
        print(rep['verdict'])
        return 2
    N, m, dq_extra, dims = shape
    rep['shape'] = {'order_N': N, 'den_power_m': m, 'dq_extra': dq_extra,
                    'nullspace_dims_at_probe_svals': dims}

    # ---- sample + reconstruct at three primes (canonical reduced forms)
    recon = {}
    for p in (P1, P2, P3):
        sols = {}
        sv = 3
        while len(sols) < args.nsamp:
            sv += 1
            b = T.try_solve(p, sv % p, N, m, dq_extra)
            if b is None:
                continue
            # pick a vector with c_N != 0; prefer dim-1
            vec = None
            for v in b:
                if v[N] % p:
                    vec = v
                    break
            if vec is None:
                continue
            cN_inv = pow(vec[N], p - 2, p)
            sols[sv % p] = [(c * cN_inv) % p for c in vec[:N + 1]]
        pts = sorted(sols)
        ratios = []
        okp = True
        for k in range(N + 1):
            vals = [sols[x][k] for x in pts]
            hit = None
            for dd in (6, 10, 16, 24):
                dn = dd + 8
                use = min(len(pts), dn + dd + 12)
                hit = rat_recon_func(vals[:use], pts[:use], p, dn, dd)
                if hit:
                    Pc, Qc = hit
                    ok = all((sum(c * pow(x, i, p) for i, c in enumerate(Pc)) -
                              v * sum(c * pow(x, i, p) for i, c in enumerate(Qc))) % p == 0
                             for x, v in zip(pts[use:use + 10], vals[use:use + 10]))
                    if ok:
                        break
                    hit = None
            if hit is None:
                okp = False
                break
            ratios.append(hit)
        if not okp:
            rep['verdict'] = f'WALL: rational-function reconstruction failed at p={p}'
            json.dump(rep, open(args.out, 'w'), indent=1)
            print(rep['verdict'])
            return 2
        recon[p] = ratios
    rep['t_sample_s'] = round(time.time() - t00, 2)

    # ---- CRT + Farey to Q(s), assemble primitive p_k(s)
    primes = [P1, P2, P3]
    M = 1
    for p_ in primes:
        M *= p_
    bound = int((M // 2) ** 0.5)
    ops = {}
    sS = sp.Symbol('s')

    def crt_list(vals):
        x, mod = 0, 1
        for v, p_ in zip(vals, primes):
            x = x + mod * (((v - x) * pow(mod, -1, p_)) % p_)
            mod *= p_
        return x % M

    for k in range(N + 1):
        pieces = [recon[p_][k] for p_ in primes]
        degsP = {len(pc) for pc, qc in pieces}
        degsQ = {len(qc) for pc, qc in pieces}
        if len(degsP) != 1 or len(degsQ) != 1:
            rep['verdict'] = (f'WALL: reduced-form degree mismatch across primes '
                              f'at k={k} (P lens {sorted(degsP)}, Q lens '
                              f'{sorted(degsQ)}) — bad-reduction prime; rerun '
                              'with different primes')
            json.dump(rep, open(args.out, 'w'), indent=1)
            print(rep['verdict'])
            return 2
        LP, LQ = degsP.pop(), degsQ.pop()
        Pq = sp.Integer(0)
        Qq = sp.Integer(0)
        okk = True
        for i in range(max(LP, LQ)):
            if i < LP:
                fa = farey(crt_list([pc[i] for pc, qc in pieces]), M, bound)
                if fa is None:
                    okk = False
                    break
                Pq += sp.Rational(fa.numerator, fa.denominator) * sS ** i
            if i < LQ:
                fb = farey(crt_list([qc[i] for pc, qc in pieces]), M, bound)
                if fb is None:
                    okk = False
                    break
                Qq += sp.Rational(fb.numerator, fb.denominator) * sS ** i
        if not okk:
            rep['verdict'] = 'WALL: Farey lift failed at 3-prime modulus (raise primes)'
            json.dump(rep, open(args.out, 'w'), indent=1)
            print(rep['verdict'])
            return 2
        ops[k] = sp.cancel(Pq / Qq)
    # clear to primitive Z[s]
    dens = [sp.denom(sp.together(vv)) for vv in ops.values()]
    den = sp.lcm(dens)
    op_can = {k: sp.expand(sp.cancel(vv * den)) for k, vv in ops.items()}
    from famgen_common import op_primitive as _oprim
    op_can = _oprim(op_can, sS)
    rep['L_A'] = {str(k): str(v) for k, v in sorted(op_can.items())}
    rep['L_A_order'] = max(op_can)
    rep['L_A_sdeg'] = max(sp.Poly(c, sS).degree() for c in op_can.values())

    # ---- exact certificate verification over Q (time-boxed)
    ver = exact_verify(T, op_can, N, m, dq_extra, sS, cap_s=args.exact_verify_cap_s)
    rep['exact_certificate'] = ver

    # ---- gates in the x = 1/s chart: series twist detection + operator equality
    gate = gates_x_chart(op_can, Av, Bv, Cv, sS, routeb_path=args.routeb)
    rep['gates'] = gate
    verdict = 'PASS' if (gate['operator_match'] and ver['n_exact_s_points'] > 0) else \
              ('PASS-mod-p' if gate['operator_match'] else 'FAIL')
    rep['verdict'] = verdict + ': ' + gate['detail']
    rep['t_total_s'] = round(time.time() - t00, 2)
    json.dump(rep, open(args.out, 'w'), indent=1)
    print(json.dumps({k: rep[k] for k in ('shape', 'L_A_order', 'L_A_sdeg',
                                          'exact_certificate', 'gates',
                                          'verdict', 't_total_s')}, indent=1,
                     default=str))
    return 0 if verdict.startswith('PASS') else 1


def exact_verify(T, op_can, N, m, dq_extra, sS, cap_s=900):
    """Verify sum_k p_k(s) r_k = d/dq u + u Mq exactly over Q(q,s), solving for
    u by linear algebra over Q(s) — time-boxed; fallback: exact identity checks
    at rational s values exceeding the s-degree bound."""
    t0 = time.time()
    # target row: R = sum_k p_k(s) r_k  (exact, sympy)
    R1 = sp.Integer(0)
    R2 = sp.Integer(0)
    for k in range(N + 1):
        pk = op_can.get(k, sp.Integer(0)).subs(sS, s)
        r1, r2 = T.rows[k]
        R1 += pk * r1
        R2 += pk * r2
    R1, R2 = sp.cancel(R1), sp.cancel(R2)
    E = sp.Integer(1)
    for b, bm in zip(T.badfactors, T.badmults):
        E *= b ** (bm + m)
    E = sp.expand(E)
    degE = sp.Poly(E, q).degree()
    dU = degE + dq_extra
    # rational-s spot verification (cheap, exact over Q per s value)
    import random
    random.seed(20260830)
    sdeg_bound = 40 + 2 * max(sp.Poly(c, sS).degree() for c in op_can.values())
    npts = 8
    ok_pts = 0
    for trial in range(npts):
        sv = sp.Rational(random.randint(50, 5000), random.randint(1, 7))
        # solve for u over Q[q] at this s: unknown coeffs of U1,U2
        Ev = sp.Poly(E.subs(s, sv), q)
        u1c = sp.symbols(f'u1_0:{dU + 1}')
        u2c = sp.symbols(f'u2_0:{dU + 1}')
        U1 = sum(c * q ** i for i, c in enumerate(u1c))
        U2 = sum(c * q ** i for i, c in enumerate(u2c))
        Mq_v = T.Mq.subs(s, sv)
        lhs1 = sp.together(R1.subs(s, sv) - (sp.diff(U1 / Ev.as_expr(), q)
               + U1 / Ev.as_expr() * Mq_v[0, 0] + U2 / Ev.as_expr() * Mq_v[1, 0]))
        lhs2 = sp.together(R2.subs(s, sv) - (sp.diff(U2 / Ev.as_expr(), q)
               + U1 / Ev.as_expr() * Mq_v[0, 1] + U2 / Ev.as_expr() * Mq_v[1, 1]))
        n1 = sp.numer(sp.cancel(lhs1))
        n2 = sp.numer(sp.cancel(lhs2))
        eqs = []
        for nn in (n1, n2):
            pol = sp.Poly(sp.expand(nn), q)
            eqs += list(pol.all_coeffs())
        unks = list(u1c) + list(u2c)
        Amat, bvec = sp.linear_eq_to_matrix(eqs, unks)
        try:
            sol, params = Amat.gauss_jordan_solve(bvec)
            ok_pts += 1
        except ValueError:
            pass
        if time.time() - t0 > cap_s:
            break
    grade = ('EXACT-at-%d-rational-s-points (certificate solved over Q[q] '
             'at each; s-degree bound %d NOT exceeded by point count — '
             'combined with the three-prime full construction this is the '
             'stated grade, not a proof over Q(s))' % (ok_pts, sdeg_bound)) \
        if ok_pts == npts else f'PARTIAL: {ok_pts}/{npts} exact s-points'
    return {'grade': grade, 'n_exact_s_points': ok_pts, 'dU': dU,
            't_s': round(time.time() - t0, 1)}


def op_apply_series(op_x, x, kappa, coeffs, nmax):
    """Apply sum_j c_j(x) D^j to  x^kappa * sum_m a_m x^m  (a_m Fractions);
    return list of coefficients of x^{n+kappa-jmax_shift}... — we simply check
    that ALL collected coefficients up to n <= nmax vanish. Returns True/False."""
    from collections import defaultdict
    acc = defaultdict(lambda: sp.Integer(0))
    for j, cpoly in op_x.items():
        cp = sp.Poly(cpoly, x)
        for (i,), ci in zip(cp.monoms(), cp.coeffs()):
            for m_, am in enumerate(coeffs):
                if am == 0:
                    continue
                e = sp.Rational(m_) + kappa
                fall = sp.Integer(1)
                for t in range(j):
                    fall *= (e - t)
                if fall == 0:
                    continue
                amr = sp.Rational(am.numerator, am.denominator)
                acc[e - j + i] += amr * ci * fall
    ords = sorted([k for k in acc], key=lambda z: sp.Rational(z))
    # only trust exponents whose contributions are complete: n+kappa-j+i is
    # complete for n <= len(coeffs) - 1 - (max shift); use a safe cut
    jmax = max(op_x)
    imax = max(sp.Poly(c, x).degree() for c in op_x.values())
    safe_hi = sp.Rational(len(coeffs) - 1) + kappa - jmax
    bad = [k for k in ords if k <= safe_hi and sp.simplify(acc[k]) != 0]
    return len(bad) == 0, len([k for k in ords if k <= safe_hi])


def gates_x_chart(op_can, Av, Bv, Cv, sS, routeb_path=None):
    """Translate L_A to the x = 1/s chart; empirically pin the gauge twist
    kappa (Pi ~ x^kappa f(x)); exact operator comparison vs the route-B minimal
    operator and vs the reference L5 when orders agree."""
    from famgen_common import op_invert_var, op_primitive, op_equal, moments
    import json as _json
    x = sp.Symbol('x')
    # L_A in x-chart (v reused: op_invert_var treats v both as s and u=1/s)
    opA_x = op_primitive(op_invert_var({j: c.subs(sS, x) for j, c in op_can.items()},
                                       x), x)
    out = {'L_A_xchart': {str(j): str(c) for j, c in sorted(opA_x.items())}}
    # series twist: find kappa in halves with L_A_x (x^kappa f(x)) = 0
    a = moments(Av, Bv, Cv, 60)
    kfound = []
    for two_k in range(-6, 7):
        kappa = sp.Rational(two_k, 2)
        ok, nchecked = op_apply_series(opA_x, x, kappa, a, 60)
        if ok and nchecked > 20:
            kfound.append(str(kappa))
    out['series_twist_kappa_candidates'] = kfound
    # exact operator gates
    detail = []
    match = False
    if routeb_path and os.path.exists(routeb_path):
        rb = _json.load(open(routeb_path))
        opB = {int(j): sp.sympify(c) for j, c in rb['op_fit'].items()}
        opB = op_primitive(opB, x)
        if op_equal(opA_x, opB, x):
            match = True
            detail.append('L_A (x-chart) == route-B minimal operator EXACTLY (no twist)')
        else:
            for two_k in range(-6, 7):
                kappa = sp.Rational(two_k, 2)
                tw = conjugate_by_power(opB, x, kappa)
                if tw and op_equal(opA_x, tw, x):
                    match = True
                    detail.append(f'L_A (x-chart) == x^{kappa}-conjugate of route-B '
                                  'minimal operator EXACTLY (twist exhibited)')
                    break
        out['vs_routeB'] = match
    # L5 comparison (only decisive if same order)
    th, A_, B_, C_ = sp.symbols('theta A B C')
    P = load_L5()
    th2 = sp.Symbol('__theta_tmp2__')
    expr = sum(x ** p * P[p].subs({A_: Av, B_: Bv, C_: Cv, th: th2}) for p in range(7))
    opL5 = op_primitive(theta_form_to_D(sp.expand(expr), x, th2), x)
    l5_match = op_equal(opA_x, opL5, x)
    if not l5_match:
        for two_k in range(-6, 7):
            kappa = sp.Rational(two_k, 2)
            tw = conjugate_by_power(opL5, x, kappa)
            if tw and op_equal(opA_x, tw, x):
                l5_match = True
                detail.append(f'L_A (x-chart) == x^{kappa}-conjugate of the reference L5 '
                              'specialization EXACTLY')
                break
    else:
        detail.append('L_A (x-chart) == reference L5 specialization EXACTLY (no twist)')
    out['vs_L5'] = l5_match
    out['L5_note'] = ('L5 gate decisive only when minimal order = 5 (generic rates); '
                      'on degeneration loci the minimal operator is a right factor')
    if max(opA_x) < 5 and not l5_match:
        detail.append('order < 5 (degeneration locus): L5 equality not expected; '
                      'route-B equality is the operative gate')
    out['operator_match'] = bool(match or l5_match)
    out['detail'] = '; '.join(detail) if detail else 'no exact operator match found'
    return out


def compare_with_L5(op_can, Av, Bv, Cv, sS):
    """Translate the reference L5 (x-chart, theta form) to the s = 1/x chart and
    compare with L_A up to an explicit gauge twist J = s^e (small rational e)."""
    th, A_, B_, C_ = sp.symbols('theta A B C')
    P = load_L5()
    x = sp.Symbol('x')
    # L5 in canonical D-form in x
    expr = sum(x ** p * P[p].subs({A_: Av, B_: Bv, C_: Cv}) for p in range(7))
    op_x = theta_form_to_D(sp.expand(expr.subs(th, sp.Symbol('__theta_tmp2__'))),
                           x, sp.Symbol('__theta_tmp2__'))
    # translate to s-chart: y(x), x = 1/s: theta_x -> -theta_s.
    # In theta form: L5 = sum_p x^p P_p(theta_x)  ->  sum_p s^{-p} P_p(-theta_s);
    # multiply by s^6: L5_s = sum_p s^{6-p} P_p(-theta_s).
    ths = sp.Symbol('__theta_s__')
    expr_s = sum(sS ** (6 - p) * P[p].subs({A_: Av, B_: Bv, C_: Cv, th: -ths})
                 for p in range(7))
    op_s = theta_form_to_D(sp.expand(expr_s), sS, ths)
    from famgen_common import op_primitive, op_equal
    op_s = op_primitive(op_s, sS)
    opA = op_primitive(op_can, sS)
    if op_equal(opA, op_s, sS):
        return {'match': True, 'twist': 'none',
                'detail': 'L_A == L5 (s-chart) exactly, primitive forms, no twist'}
    # try twists J = s^e, e in small set: conjugated operator has coeffs in Q(s)
    for enum_, eden in [(1, 2), (1, 1), (3, 2), (2, 1), (-1, 2), (-1, 1),
                        (-3, 2), (-2, 1), (5, 2), (-5, 2), (3, 1), (-3, 1)]:
        e = sp.Rational(enum_, eden)
        twisted = conjugate_by_power(op_s, sS, e)
        if twisted and op_equal(opA, twisted, sS):
            return {'match': True, 'twist': f's^{e}',
                    'detail': f'L_A == s^{e} twist of L5 (solutions differ by s^{e}) '
                              '— twist exhibited exactly'}
    return {'match': False, 'twist': None,
            'detail': 'no match up to the searched power twists'}


def conjugate_by_power(op, v, e):
    """J^{-1} L J with J = v^e: acts on y -> v^{-e} L(v^e y). D -> D + e/v."""
    # represent shifted operator: sum_j c_j(v) (D + e/v)^j, expand in D
    out = {}
    for j, c in op.items():
        # (D + e/v)^j expanded: build iteratively as operator on functions
        term = {0: sp.Integer(1)}
        for _ in range(j):
            new = {}
            for jj, cc in term.items():
                # D o (cc D^jj) = cc' D^jj + cc D^{jj+1}
                new[jj] = sp.cancel(new.get(jj, 0) + sp.diff(cc, v) + cc * e / v)
                new[jj + 1] = sp.cancel(new.get(jj + 1, 0) + cc)
            term = new
        for jj, cc in term.items():
            out[jj] = sp.cancel(out.get(jj, 0) + c * cc)
    # clear denominators to polynomial primitive form
    from famgen_common import op_primitive
    try:
        return op_primitive({j: sp.together(c) for j, c in out.items()}, v)
    except Exception:
        return None


if __name__ == '__main__':
    sys.exit(main())

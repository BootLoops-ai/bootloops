"""Exact Cheng-Wu cell assembly of the LO collinear E4C integrand, per channel and
chart (gauge variable x_g = 1), with the SHAPE KEPT SYMBOLIC through the six pair
distances d_ab (a<b).  sigma-instances and shapes are obtained later by numeric
substitution of d_ab -> d_{sigma(a) sigma(b)}(shape).

Object (the P files are the tree-level 1->4 splitting-function term lists, homogeneous in z of
degree h (=1 quark files, 0/.. gluon) so z -> x/Om is exact):
  I_sigma = int_{(0,inf)^3} d^3x  (x1 x2 x3 x4)^2|_{x_g=1} * Phat(x; d^sigma) / ( s1234^3 * Om^6 )
  Phat: z_i -> x_i/Om, z_S -> ell_S/Om, s_ij -> d_ij x_i x_j, s_S -> sum_{a<b in S} d_ab x_a x_b,
        CF=4/3, CA=3, D=4, gg -> -1, kikj -> (1/2) x_i x_j G_ij / Om^2 (x-version bookkeeping, see build_cells),
        G_ij = sum_{k,l} x_k x_l c_{ik,jl},  c_{ik,jl} = (w_i-w_k).(w_j-w_l) = (d_il + d_jk - d_ij - d_kl)/2.
  G4_channel = xL^3 * w_flavor(nf) * sfac * sum_sigma I_sigma   (prefactor applied by the evaluator).

Cell = denominator multiset over atoms {x_j, Om, ell_S, s_S(|S|>=3 incl. s1234)};
numerator = dict { x-exponent tuple (3 chart vars) : fmpq_mpoly in D=(d12,d13,d14,d23,d24,d34) }
with a global d-shift (all d exponents shifted by +DSHIFT so they are >= 0; the evaluator
divides by prod d^DSHIFT).

Cache: gate/cells/<channel>_g<g>.pkl (shipped for the channels listed in DATA.md); a cache miss builds the cells from the
term lists under $SURD_TERMS and writes the new pickle under $SURD_WORK/cells (never into the
package tree).  Gate (exact): cell sum == direct Fraction evaluation of the flat term list at random rational (x, d) points
(function gate_cells).
"""
import os, sys, time, pickle, itertools, random
from fractions import Fraction as F
from flint import fmpq, fmpq_mpoly_ctx

GATE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_ROOT_CH = os.path.dirname(GATE)
# term lists (<file>.terms.pkl, read ONLY on a cells-cache miss): $SURD_TERMS, else <root>/data/terms
CACHE_TERMS = os.environ.get("SURD_TERMS") or os.path.join(_ROOT_CH, "data", "terms")
CELLDIR = os.path.join(GATE, "cells")          # shipped cell pickles (read-only)
# new cell pickles built on a cache miss go under the output root ($SURD_WORK, default ./surd_work), never into the package tree
WORK_CELLDIR = os.path.join(os.path.abspath(os.environ.get("SURD_WORK") or os.path.join(os.getcwd(), "surd_work")), "cells")

CHFILES = {"q_qbpqpgq": "P_qbp_qp_g_q.m", "q_qbqgq": "P_qb_q_g_q.m", "q_gggq": "P_symm_g_g_g_q.m",
           "g_qbpqpqbq": "P_qbp_qp_qb_q.m", "g_qbqqbq": "P_qb_q_qb_q.m", "g_qbggq": "P_qb_g_g_q.m",
           "g_gggg": "P_symm_g_g_g_g.m",
           "n4": "P_N4_coll.m"}      # CONTROL P: N=4 1->4 collinear integrand from |F_5|^2, <P>_N4 = s1234^3 P_N4
# candidate exact particle-relabeling symmetries tau of each channel's squared splitting amplitude (I_{sigma tau} = I_sigma); consumers
# (slice/scripts/slice_classes.py, oracle.py) VERIFY each tau exactly on the cells before using it, so this table is an input, not a claim
RELABELING_SYMMETRIES = {"q_qbpqpgq": [], "q_qbqgq": [[1, 4, 3, 2]], "q_gggq": [],
                         "g_qbpqpqbq": [[2, 1, 4, 3], [3, 4, 1, 2], [4, 3, 2, 1]],
                         "g_qbqqbq": [[1, 4, 3, 2], [2, 1, 4, 3], [2, 3, 4, 1], [3, 2, 1, 4], [3, 4, 1, 2], [4, 1, 2, 3], [4, 3, 2, 1]],
                         "g_qbggq": [[1, 3, 2, 4], [4, 2, 3, 1], [4, 3, 2, 1]], "g_gggg": [[2, 1, 3, 4]], "n4": []}
# flavor weight at nf, identical-particle factor
def flavor_weight(ch, nf):
    nf = F(nf)
    return {"q_qbpqpgq": nf - 1, "q_qbqgq": F(1), "q_gggq": F(1), "g_qbpqpqbq": nf * (nf - 1) / 2,
            "g_qbqqbq": nf, "g_qbggq": nf, "g_gggg": F(1), "n4": F(1)}[ch]
SFAC = {"q_qbpqpgq": F(1), "q_qbqgq": F(1, 2), "q_gggq": F(1), "g_qbpqpqbq": F(1),
        "g_qbqqbq": F(1, 4), "g_qbggq": F(1, 2), "g_gggg": F(1, 2), "n4": F(1, 24)}
QUARK_JET = ["q_qbpqpgq", "q_qbqgq", "q_gggq"]
GLUON_JET = ["g_qbpqpqbq", "g_qbqqbq", "g_qbggq", "g_gggg"]

PAIRS = [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]
PIDX = {p: i for i, p in enumerate(PAIRS)}
def pidx(a, b):
    return PIDX[(a, b) if a < b else (b, a)]
DSHIFT = 4   # global shift of d exponents (max negative d power per term is <= 4: s_ij^-2 twice... checked at build)
CFv, CAv, Dv = F(4, 3), F(3), F(4)
DCTX = fmpq_mpoly_ctx.get(("d12", "d13", "d14", "d23", "d24", "d34"))


def load_terms(ch):
    p = os.path.join(CACHE_TERMS, CHFILES[ch] + ".terms.pkl")
    if not os.path.exists(p):
        raise FileNotFoundError("surd: no cell pickle for channel %r in gate/cells or %s, and the term list %s is not available: "
                                "set SURD_TERMS to a directory holding %s.terms.pkl (format in DATA.md)" % (ch, WORK_CELLDIR, p, CHFILES[ch]))
    with open(p, "rb") as f:
        return pickle.load(f)


class ChartAlgebra:
    """Polynomials in the 3 chart variables (x_j, j != g, increasing j) and the 6 d's,
    represented as dicts {(ex(3), ed(6)) : Fraction}."""
    def __init__(self, g):
        self.g = g
        self.xvars = [j for j in (1, 2, 3, 4) if j != g]
        self.xpos = {j: i for i, j in enumerate(self.xvars)}

    def mono_x(self, j, e=1):
        ex = [0, 0, 0]
        if j != self.g:
            ex[self.xpos[j]] = e
        return tuple(ex)

    def poly_ell(self, S):
        p = {}
        for j in S:
            k = (self.mono_x(j), (0,) * 6)
            p[k] = p.get(k, F(0)) + 1
        return p

    def poly_Om(self):
        return self.poly_ell((1, 2, 3, 4))

    def poly_s(self, S):
        p = {}
        for a, b in itertools.combinations(sorted(S), 2):
            ex = [0, 0, 0]
            for j in (a, b):
                if j != self.g:
                    ex[self.xpos[j]] += 1
            ed = [0] * 6; ed[pidx(a, b)] = 1
            p[(tuple(ex), tuple(ed))] = F(1)
        return p

    def poly_G(self, i, j):
        """G_ij = sum_{k,l} x_k x_l c_{ik,jl}, c = (d_il + d_jk - d_ij - d_kl)/2 (d_aa = 0)."""
        p = {}
        def addterm(k, l, pair, coef):
            if pair[0] == pair[1]:
                return
            ex = [0, 0, 0]
            for m in (k, l):
                if m != self.g:
                    ex[self.xpos[m]] += 1
            ed = [0] * 6; ed[pidx(*pair)] = 1
            key = (tuple(ex), tuple(ed))
            p[key] = p.get(key, F(0)) + coef
        for k in (1, 2, 3, 4):
            for l in (1, 2, 3, 4):
                addterm(k, l, (i, l), F(1, 2))
                addterm(k, l, (j, k), F(1, 2))
                addterm(k, l, (i, j), F(-1, 2))
                addterm(k, l, (k, l), F(-1, 2))
        return {k: c for k, c in p.items() if c != 0}

    @staticmethod
    def mul(p, q):
        out = {}
        for (ax, ad), ac in p.items():
            for (bx, bd), bc in q.items():
                kx = (ax[0] + bx[0], ax[1] + bx[1], ax[2] + bx[2])
                kd = tuple(ad[i] + bd[i] for i in range(6))
                k = (kx, kd)
                out[k] = out.get(k, F(0)) + ac * bc
        return {k: c for k, c in out.items() if c != 0}

    @staticmethod
    def pow(p, e):
        r = {((0, 0, 0), (0,) * 6): F(1)}
        for _ in range(e):
            r = ChartAlgebra.mul(r, p)
        return r


def parse_name(nm):
    """-> (kind, data): kind in {'const','z','s','gg','kk'}"""
    if nm in ("CF", "CA", "D"):
        return "const", {"CF": CFv, "CA": CAv, "D": Dv}[nm]
    if nm == "gg":
        return "gg", None
    if nm[0] == "z":
        return "z", tuple(int(c) for c in nm[1:])
    if nm[0] == "s":
        return "s", tuple(int(c) for c in nm[1:])
    if nm[0] == "k" and len(nm) == 4 and nm[2] == "k":
        return "kk", (int(nm[1]), int(nm[3]))
    raise ValueError(nm)


def build_cells(ch, g, verbose=True):
    for cdir in (CELLDIR, WORK_CELLDIR):
        cpath = os.path.join(cdir, "%s_g%d.pkl" % (ch, g))
        if os.path.exists(cpath):
            with open(cpath, "rb") as f:
                return _hydrate(pickle.load(f))
    cpath = os.path.join(WORK_CELLDIR, "%s_g%d.pkl" % (ch, g))
    t0 = time.time()
    terms, names = load_terms(ch)
    A = ChartAlgebra(g)
    kinds = [parse_name(nm) for nm in names]
    # per-name precomputed numerator polynomial (for composite s_S numerators, G_ij) and power caches
    polycache = {}
    def numpoly(ni, e):
        key = (ni, e)
        if key not in polycache:
            kind, data = kinds[ni]
            if kind == "s":
                base = A.poly_s(data)
            elif kind == "kk":
                base = A.poly_G(*data)
            elif kind == "z":
                base = A.poly_ell(data)
            else:
                raise ValueError
            polycache[key] = ChartAlgebra.pow(base, e)
        return polycache[key]
    cells = {}      # denkey -> {(ex,ed): Fraction}
    ONEKEY = ((0, 0, 0), (0,) * 6)
    nt = 0
    hz_check = set()
    for key, c in terms.items():
        nt += 1
        coeff = F(c)
        ex = [0, 0, 0]; ed = [0] * 6
        om = -6          # measure Om^-6
        den = {}          # atom name -> power (>0)
        polyf = []        # list of numerator poly dicts
        # measure (x1x2x3x4)^2 / s1234^3
        for j in A.xvars:
            ex[A.xpos[j]] += 2
        den["s1234"] = 3
        netz = 0
        for (ni, e) in key:
            kind, data = kinds[ni]
            if kind == "const":
                coeff *= data ** e
            elif kind == "gg":
                assert e == 1
                coeff = -coeff
            elif kind == "z":
                netz += e
                S = data
                if len(S) == 4:
                    continue            # z_1234 = (x1+x2+x3+x4)/Om = 1 identically in the chart (homogenizing factor of the N=4 term list)
                om -= e
                if len(S) == 1:
                    j = S[0]
                    if j != g:
                        ex[A.xpos[j]] += e      # may go negative -> denominator monomial
                else:
                    nm = "ell" + "".join(map(str, S))
                    if e > 0:
                        polyf.append(numpoly(ni, e))
                    else:
                        den[nm] = den.get(nm, 0) - e
            elif kind == "s":
                S = data
                if len(S) == 2:
                    a, b = S
                    ed[pidx(a, b)] += e
                    for j in (a, b):
                        if j != g:
                            ex[A.xpos[j]] += e
                else:
                    nm = "s" + "".join(map(str, S))
                    if e > 0:
                        polyf.append(numpoly(ni, e))
                    else:
                        den[nm] = den.get(nm, 0) - e
            elif kind == "kk":
                assert e == 1
                i, j = data
                coeff *= F(1, 2)
                for m in (i, j):
                    if m != g:
                        ex[A.xpos[m]] += 1
                # kikj^z = (1/2) z_i z_j G_ij(z) = Om^-4 * (1/2) x_i x_j G_ij(x); the s-like quantities are
                # all converted to their x-versions (s^x = Om^2 s^z) and P is degree 0 in {s, kk} jointly,
                # so the kk term needs Om^{-4} * Om^{+2} = Om^-2 relative to the x-version (exact gate).
                om -= 2
                polyf.append(numpoly(ni, 1))
        hz_check.add(netz)
        # Om: net power om (<0 always expected)
        if om < 0:
            den["Om"] = -om
        elif om > 0:
            polyf.append(ChartAlgebra.pow(A.poly_Om(), om))
        # negative x exponents -> denominator atoms
        for j in A.xvars:
            p = A.xpos[j]
            if ex[p] < 0:
                den["x%d" % j] = -ex[p]
                ex[p] = 0
        # d shift
        eds = tuple(ed[i] + DSHIFT for i in range(6))
        assert min(eds) >= 0, ("DSHIFT too small", ed)
        base = {(tuple(ex), eds): coeff}
        for pf in polyf:
            base = ChartAlgebra.mul(base, pf)
        dk = tuple(sorted(den.items()))
        cell = cells.get(dk)
        if cell is None:
            cells[dk] = dict(base)
        else:
            for k2, c2 in base.items():
                v = cell.get(k2)
                if v is None:
                    cell[k2] = c2
                else:
                    v = v + c2
                    if v:
                        cell[k2] = v
                    else:
                        del cell[k2]
        if verbose and nt % 200000 == 0:
            print("  %s g%d: %d/%d terms, %d cells, %.0fs" % (ch, g, nt, len(terms), len(cells), time.time() - t0), flush=True)
    assert len(hz_check) == 1, ("P not z-homogeneous", hz_check)
    hz = hz_check.pop()
    # convert: per cell, group by x-exponent -> fmpq_mpoly in d's
    out_cells = []
    nmon = 0
    for dk, poly in cells.items():
        byx = {}
        for (exx, edd), cc in poly.items():
            byx.setdefault(exx, {})[edd] = (cc.numerator, cc.denominator)
        if byx:
            out_cells.append((dk, byx))
            nmon += sum(len(p) for p in byx.values())
    res = {"channel": ch, "gauge": g, "xvars": A.xvars, "dshift": DSHIFT, "hz": hz,
           "n_terms": len(terms), "cells_raw": out_cells, "n_cells": len(out_cells), "n_monomials": nmon,
           "build_s": time.time() - t0}
    os.makedirs(WORK_CELLDIR, exist_ok=True)
    tmp = cpath + ".tmp%d" % os.getpid()
    with open(tmp, "wb") as f:
        pickle.dump(res, f, protocol=4)
    os.replace(tmp, cpath)
    res = _hydrate(res)
    if verbose:
        print("%s g%d: %d terms -> %d cells, %d (x,d)-monomials, hz=%d, %.0fs" % (ch, g, len(terms), len(out_cells), nmon, hz, time.time() - t0), flush=True)
    return res


def _hydrate(res):
    """raw {exx: {edd: (num,den)}} -> {exx: fmpq_mpoly over the d's}; drops zero polys/cells."""
    cells = []
    for dk, byx in res["cells_raw"]:
        conv = {}
        for exx, dd in byx.items():
            p = DCTX.from_dict({edd: fmpq(n, d) for edd, (n, d) in dd.items()})
            if not p.is_zero():
                conv[exx] = p
        if conv:
            cells.append((dk, conv))
    res = dict(res)
    del res["cells_raw"]
    res["cells"] = cells
    res["n_cells"] = len(cells)
    return res


# ---------------------------------------------------------------- exact gate
def eval_terms_direct(ch, g, xs, ds):
    """Direct Fraction evaluation of  (x1x2x3x4)^2 Phat / (s1234^3 Om^6)  at chart point xs (dict j->Fraction,
    x_g = 1) and pair distances ds (dict (a,b)->Fraction), from the flat term list."""
    terms, names = load_terms(ch)
    X = {j: (F(1) if j == g else F(xs[j])) for j in (1, 2, 3, 4)}
    Om = sum(X.values())
    z = {j: X[j] / Om for j in X}
    env = {}
    for nm in names:
        kind, data = parse_name(nm)
        if kind == "const":
            env[nm] = data
        elif kind == "gg":
            env[nm] = F(-1)
        elif kind == "z":
            env[nm] = sum(z[j] for j in data)
        elif kind == "s":
            S = data
            env[nm] = sum(ds[(a, b)] * z[a] * z[b] for a, b in itertools.combinations(sorted(S), 2))
        elif kind == "kk":
            i, j = data
            def c(i, k, j, l):
                def dd(a, b):
                    return F(0) if a == b else ds[(a, b) if a < b else (b, a)]
                return (dd(i, l) + dd(j, k) - dd(i, j) - dd(k, l)) / 2
            G = sum(z[k] * z[l] * c(i, k, j, l) for k in (1, 2, 3, 4) for l in (1, 2, 3, 4))
            env[nm] = F(1, 2) * z[i] * z[j] * G
    # note: with z (not x) in s and kk, P(z; s(z)) is the simplex-form value; Phat(x) = P(z) since P is
    # s-degree-0 homogeneous... careful: s-degree 0 => s(z) vs s(x) differ by Om^2 per s -> cancels. z-degree hz
    # is accounted by using z directly. So P evaluated here == Phat(x)/Om^0 with z=x/Om exactly.
    vals = [env[nm] for nm in names]
    P = F(0)
    for key, c in terms.items():
        t = F(c)
        for (ni, e) in key:
            t *= vals[ni] ** e
        P += t
    s1234x = sum(ds[(a, b)] * X[a] * X[b] for a, b in PAIRS)
    meas = F(1)
    for j in (1, 2, 3, 4):
        meas *= X[j] ** 2
    return meas * P / (s1234x ** 3 * Om ** 6)


def eval_cells_exact(cellres, xs, ds):
    """Exact Fraction evaluation of the cell sum at chart point xs, distances ds."""
    g = cellres["gauge"]; xv = cellres["xvars"]; sh = cellres["dshift"]
    X = {j: (F(1) if j == g else F(xs[j])) for j in (1, 2, 3, 4)}
    dvals = [fmpq(ds[p].numerator, ds[p].denominator) for p in PAIRS]
    dshiftinv = F(1)
    for p in PAIRS:
        dshiftinv /= ds[p] ** sh
    def atom(nm):
        if nm == "Om":
            return sum(X.values())
        if nm[0] == "x":
            return X[int(nm[1])]
        if nm.startswith("ell"):
            return sum(X[int(c)] for c in nm[3:])
        if nm[0] == "s":
            S = [int(c) for c in nm[1:]]
            return sum(ds[(a, b)] * X[a] * X[b] for a, b in itertools.combinations(S, 2))
        raise ValueError(nm)
    tot = F(0)
    for dk, conv in cellres["cells"]:
        den = F(1)
        for nm, pw in dk:
            den *= atom(nm) ** pw
        num = F(0)
        for exx, p in conv.items():
            v = p(*dvals)
            v = F(int(v.numer()), int(v.denom())) if hasattr(v, "numer") else F(str(v))
            mono = F(1)
            for jj, j in enumerate(xv):
                mono *= X[j] ** exx[jj]
            num += v * mono
        tot += num * dshiftinv / den
    return tot


def gate_cells(ch, g, npts=2, seed=1):
    cellres = build_cells(ch, g)
    rng = random.Random(seed)
    worst = []
    for _ in range(npts):
        xs = {j: F(rng.randint(1, 9), rng.randint(1, 9)) for j in (1, 2, 3, 4)}
        ds = {p: F(rng.randint(1, 9), rng.randint(1, 9)) for p in PAIRS}
        a = eval_cells_exact(cellres, xs, ds)
        b = eval_terms_direct(ch, g, xs, ds)
        worst.append(a == b)
        print("  gate %s g%d: cells == direct : %s   (value %s)" % (ch, g, a == b, float(a)), flush=True)
    return all(worst)


def gate_vs_e4c(ch, g, seed=7):
    """Optional float check of the RECIPE against an external float integrand engine (a module e4c.py providing get_channel/env_for
    and mparse.rpn_eval, found under $SURD_E4C_ENGINE; not part of this package): direct evaluator (z=x/Om, tensor
    contraction) vs that engine at one random interior simplex point, sigma = id, planar rational shape."""
    eng = os.environ.get("SURD_E4C_ENGINE")
    if not eng: raise RuntimeError("surd: gate_vs_e4c needs SURD_E4C_ENGINE (directory of an external float engine e4c.py); it is an optional check")
    sys.path.insert(0, eng)
    import e4c, mparse, numpy as np
    rng = random.Random(seed)
    pts = [(F(0), F(0)), (F(1), F(0)), (F(rng.randint(1, 9), 10), F(rng.randint(1, 9), 10)),
           (F(rng.randint(1, 9), 10), F(-rng.randint(1, 9), 10))]
    ds = {}
    for a, b in PAIRS:
        dx = pts[a - 1][0] - pts[b - 1][0]; dy = pts[a - 1][1] - pts[b - 1][1]; ds[(a, b)] = dx * dx + dy * dy
    xs = {j: F(rng.randint(1, 9), rng.randint(1, 9)) for j in (1, 2, 3, 4)}
    xs[g] = F(1)
    mine = eval_terms_direct(ch, g, xs, ds)
    Om = sum(xs.values())
    z = [float(xs[j] / Om) for j in (1, 2, 3, 4)]
    w = [complex(float(p[0]), float(p[1])) for p in pts]
    rpn, tensor, wfun, sfac = e4c.get_channel(ch)
    env = e4c.env_for((0, 1, 2, 3), w, tuple(np.array([zz]) for zz in z), tensor)
    P = mparse.rpn_eval(*rpn, env)[0]
    val = np.prod([zz ** 2 for zz in z]) / env["s1234"][0] ** 3 * P * float(Om) ** -4
    rel = abs(float(mine) - val) / abs(val)
    print("  gate-F %s g%d: direct(recipe) vs e4c.py integrand rel diff %.2e" % (ch, g, rel), flush=True)
    return rel


if __name__ == "__main__":
    import argparse, json, subprocess
    ap = argparse.ArgumentParser()
    ap.add_argument("--channels", default="q_qbpqpgq")
    ap.add_argument("--gauges", default="4")
    ap.add_argument("--gate", type=int, default=1, help="number of random exact points for the cells-vs-term-list gate (needs the term list); 0 = build only")
    ap.add_argument("--e4c", action="store_true", help="also run the optional float check against $SURD_E4C_ENGINE")
    a = ap.parse_args()
    rec = {"script": os.path.abspath(__file__), "builds": []}
    for ch in a.channels.split(","):
        for g in [int(x) for x in a.gauges.split(",")]:
            r = build_cells(ch, g)
            ok = gate_cells(ch, g, npts=a.gate) if a.gate else None
            relF = gate_vs_e4c(ch, g) if a.e4c else None
            rec["builds"].append({"channel": ch, "gauge": g, "n_terms": r["n_terms"], "n_cells": r["n_cells"],
                                  "n_monomials": r["n_monomials"], "hz": r["hz"], "build_s": round(r["build_s"], 1),
                                  "exact_gate_pass": ok, "recipe_vs_e4c_rel": relF})
    rec["stamp_utc"] = subprocess.check_output(["date", "-u", "+%Y-%m-%dT%H:%M:%SZ"]).decode().strip()
    os.makedirs(os.path.dirname(WORK_CELLDIR), exist_ok=True)
    p = os.path.join(os.path.dirname(WORK_CELLDIR), "CELLS_%s.json" % a.channels.replace(",", "_")[:60])
    json.dump(rec, open(p, "w"), indent=1)
    print("wrote", p)

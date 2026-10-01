"""Slice common: the dipole slice of the collinear four-point energy correlator: detector positions w1 = +(t/2) e^{i phi}, w2 = -(t/2) e^{i phi},
w3 = -1/2, w4 = +1/2, cos phi = 3/5, i.e. the labeled distance matrix (sigma = id: parton a at position a)
    d12 = t^2, d34 = 1, d13 = d24 = q_+(t), d14 = d23 = q_-(t),  q_+- = t^2/4 +- 3t/10 + 1/4  (all POLYNOMIAL in t).
The integrand depends on the shape only through the six d_ab (s_S = sum d_ab x_a x_b; gluon spin dot products are (d_il+d_jk-d_ij-d_kl)/2,
already inside gate/scripts/oracle_cells.py cells), so the exact group integrand on the slice is the d-symbolic cell object composed with
d_ab -> d_{sigma(a) sigma(b)}(t): a polynomial N(x_a, x_b, x_c, t) over prod atoms^M (external minimal d-powers divided out and recorded)."""
import os, sys, itertools, pickle
from fractions import Fraction as Fr
import slice_prov
from flint import fmpq, fmpq_poly, fmpq_mpoly_ctx
import alphabet, fibr
import oracle_cells as OC, hf_piece
ROOT = slice_prov.ROOT; GATE = os.path.join(ROOT, "gate")
CHART = {"12": 4, "13": 4, "14": 3, "23": 4, "24": 3, "34": 2}      # chart (gauge variable) used for each first pair = (gauge rest[1] | last rest[0])
PAIRS = OC.PAIRS
QP = (Fr(1, 4), Fr(3, 10), Fr(1, 4)); QM = (Fr(1, 4), Fr(-3, 10), Fr(1, 4))     # q_+-(t) coefficients of t^2, t, 1 ... stored as (c2, c1, c0)
def qplus(t): t = Fr(t); return t * t / 4 + Fr(3, 10) * t + Fr(1, 4)
def qminus(t): t = Fr(t); return t * t / 4 - Fr(3, 10) * t + Fr(1, 4)
def d_of_t(t):
    """labeled distance matrix of the slice at rational t (sigma = id)"""
    t = Fr(t); return {(1, 2): t * t, (3, 4): Fr(1), (1, 3): qplus(t), (2, 4): qplus(t), (1, 4): qminus(t), (2, 3): qminus(t)}
def d_polys_t():
    """the six d_ab as fmpq_poly in t"""
    return {(1, 2): fmpq_poly([0, 0, 1]), (3, 4): fmpq_poly([1]), (1, 3): fmpq_poly([fmpq(1, 4), fmpq(3, 10), fmpq(1, 4)]), (2, 4): fmpq_poly([fmpq(1, 4), fmpq(3, 10), fmpq(1, 4)]),
            (1, 4): fmpq_poly([fmpq(1, 4), fmpq(-3, 10), fmpq(1, 4)]), (2, 3): fmpq_poly([fmpq(1, 4), fmpq(-3, 10), fmpq(1, 4)])}
def sigma_d(d, sigma):
    out = {}
    for a, b in PAIRS:
        p, q = sigma[a - 1], sigma[b - 1]; out[(a, b)] = d[(p, q) if p < q else (q, p)]
    return out
def sigma_classes_slice():
    """the 24 assignments grouped by labeled d-matrix as POLYNOMIALS in t (exact): -> list of (sigma_rep, mult, members)"""
    dp = d_polys_t(); cl = {}
    for sig in itertools.permutations((1, 2, 3, 4)):
        ds = sigma_d(dp, sig); key = tuple(tuple(str(c) for c in ds[p].coeffs()) for p in PAIRS)
        cl.setdefault(key, []).append(sig)
    return [(v[0], len(v), v) for v in cl.values()]
def order_for(pair, gauge):
    k, j = sorted(int(c) for c in pair); l = [i for i in (1, 2, 3, 4) if i not in (k, j, gauge)][0]
    return [k, j, l]
def slice_ctx(xvars):
    """ctx (x_a, x_b, x_c, t) and the images of (x_a, x_b, x_c, d12, d13, d14, d23, d24, d34) for sigma: callable"""
    names = tuple("x%d" % j for j in xvars) + ("t",)
    C = fmpq_mpoly_ctx.get(names); g = C.gens(); n = len(names); tt = g[len(xvars)]
    def cst(f): f = Fr(f); return C.from_dict({(0,) * n: fmpq(f.numerator, f.denominator)})
    dpoly = {(1, 2): tt * tt, (3, 4): cst(1), (1, 3): cst(Fr(1, 4)) * tt * tt + cst(Fr(3, 10)) * tt + cst(Fr(1, 4)), (1, 4): cst(Fr(1, 4)) * tt * tt - cst(Fr(3, 10)) * tt + cst(Fr(1, 4))}
    dpoly[(2, 4)] = dpoly[(1, 3)]; dpoly[(2, 3)] = dpoly[(1, 4)]
    return C, list(g[:len(xvars)]), dpoly
def build_group(channel, gauge, pair, sigma=(1, 2, 3, 4), cell=None, support=None):
    """exact (channel, chart x_gauge=1, first-pair group [optionally one cell or one denominator-support class], sigma) integrand on the slice:
    {"N": fmpq_mpoly in (x_a,x_b,x_c,t), "atoms": {name: poly}, "M": powers, "external_d_powers": {dab: e}, "ext_factor_poly": fmpq_poly in t, ...};
    true group integrand = N / prod atoms^M / ext_factor(t),  ext_factor = prod_p d^sigma_p(t)^ext_p."""
    cellres = OC.build_cells(channel, gauge, verbose=False)
    pair = tuple(sorted(int(c) for c in pair))
    sub = hf_piece.group_cells(cellres, pair)
    if support is not None:
        S = frozenset(support); sub = dict(sub); sub["cells"] = [(dk, conv) for dk, conv in sub["cells"] if frozenset(nm for nm, pw in dk) == S]
    if cell is not None:
        cl = sorted(sub["cells"], key=lambda kv: (len(kv[0]), str(kv[0]))); sub = dict(sub); sub["cells"] = [cl[cell]]
    xv = cellres["xvars"]
    names9 = tuple("x%d" % j for j in xv) + ("d12", "d13", "d14", "d23", "d24", "d34")
    C9 = fmpq_mpoly_ctx.get(names9); g9 = {n: C9.from_dict({tuple(1 if i == k else 0 for i in range(9)): fmpq(1)}) for k, n in enumerate(names9)}
    one = C9.from_dict({(0,) * 9: fmpq(1)})
    X = {j: (one if j == gauge else g9["x%d" % j]) for j in (1, 2, 3, 4)}; D = {p: g9["d%d%d" % p] for p in PAIRS}
    def atom(nm):
        if nm == "Om": return X[1] + X[2] + X[3] + X[4]
        if nm[0] == "x": return X[int(nm[1])]
        if nm.startswith("ell"):
            Sx = [int(c) for c in nm[3:]]; r = X[Sx[0]]
            for j in Sx[1:]: r = r + X[j]
            return r
        Sx = [int(c) for c in nm[1:]]; r = None
        for a_, b_ in itertools.combinations(Sx, 2):
            tt = D[(a_, b_)] * X[a_] * X[b_]; r = tt if r is None else r + tt
        return r
    M = {}
    for dk, conv in sub["cells"]:
        for nm, pw in dk: M[nm] = max(M.get(nm, 0), pw)
    atoms = sorted(M); apoly = {nm: atom(nm) for nm in atoms}
    dgen = [g9["d%d%d" % p] for p in PAIRS]
    dmin = [10 ** 9] * 6
    for dk, conv in sub["cells"]:
        for exx, p in conv.items():
            for ed in p.to_dict().keys(): dmin = [min(a_, b_) for a_, b_ in zip(dmin, ed)]
    dmin = [0 if m == 10 ** 9 else m for m in dmin]
    N = C9.from_dict({(0,) * 9: fmpq(0)})
    for dk, conv in sub["cells"]:
        Nc = None
        for exx, p in conv.items():
            pdict = {tuple(e - m for e, m in zip(ed, dmin)): c for ed, c in p.to_dict().items()}
            p2 = p.context().from_dict(pdict); pd = p2.compose(*dgen, ctx=C9)
            mono = C9.from_dict({tuple(list(exx) + [0] * 6): fmpq(1)}); tt = pd * mono
            Nc = tt if Nc is None else Nc + tt
        have = dict(dk); cof = one
        for nm in atoms:
            e = M[nm] - have.get(nm, 0)
            if e: cof = cof * apoly[nm] ** e
        N = N + Nc * cof
    CS, gx, dpoly = slice_ctx(xv)
    dsig = sigma_d(dpoly, sigma)
    images = list(gx) + [dsig[p] for p in PAIRS]
    NS = N.compose(*images, ctx=CS); atomsS = {nm: apoly[nm].compose(*images, ctx=CS) for nm in atoms}
    shift = cellres["dshift"]
    ext = {("d%d%d" % p): int(shift - m) for p, m in zip(PAIRS, dmin)}
    dpt = sigma_d(d_polys_t(), sigma); extpoly = fmpq_poly([1])
    for p in PAIRS:
        e = ext["d%d%d" % p]
        if e > 0: extpoly = extpoly * dpt[p] ** e
        elif e < 0: raise ValueError("negative external d power")
    return {"N": NS, "atoms": atomsS, "M": dict(M), "xvars": list(xv), "n_cells": len(sub["cells"]), "external_d_powers": ext, "ext_factor_poly": extpoly,
            "dsigma_polys": {"d%d%d" % p: str(dpt[p]) for p in PAIRS}, "channel": channel, "gauge": gauge, "pair": pair, "sigma": tuple(sigma), "cell": cell, "support": support, "ctx": CS}
def ext_factor_value(G, tval):
    v = G["ext_factor_poly"](fmpq(Fr(tval).numerator, Fr(tval).denominator)); return Fr(int(v.p), int(v.q))

class SliceAlphabet(alphabet.Alphabet):
    """auto-extending alphabet on the slice (unknown irreducible factors are ADDED and recorded, never refused)"""
    def __init__(self, ctx):
        super().__init__(ctx); self.auto_added = []
    def factor(self, p, context=""):
        k = str(p)
        if k in self._fact_cache: return self._fact_cache[k]
        if p.is_zero(): raise ValueError("factor(0) " + context)
        cont, facs = p.factor(); c = Fr(int(cont.p), int(cont.q)); E = {}
        for f, e in facs:
            cf, q = alphabet.normalize(f); c *= cf ** int(e)
            if f.total_degree() == 0: continue
            kk = str(q)
            if kk not in self.key2idx:
                idx = self.add(q, "auto:" + context.split(":")[0]); self.auto_added.append((idx, kk, context[:60]))
            i = self.key2idx[kk]; E[i] = E.get(i, 0) + int(e)
        res = (c, E); self._fact_cache[k] = res
        return res

def ts_dump(ts, path, ctx_names, alph, extra=None):
    ser = [([list(x) for x in Dk], [[[list(h[0]), [list(y) for y in h[1]]] for h in wd] for wd in Wk], [(list(int(z) for z in e), int(cf.p), int(cf.q)) for e, cf in N.to_dict().items()]) for (Dk, Wk), N in ts.items()]
    tmp = path + ".tmp%d" % os.getpid()
    with open(tmp, "wb") as fh: pickle.dump({"ctx_names": list(ctx_names), "alphabet": [str(p) for p in alph.letters], "alphabet_dicts": [[(list(int(z) for z in e), int(c.p), int(c.q)) for e, c in p.to_dict().items()] for p in alph.letters], "ts": ser, "extra": extra}, fh, protocol=4)
    os.replace(tmp, path)
def ts_load(path):
    """-> (ctx, letters list (fmpq_mpoly), ts dict, extra) rebuilt in a fresh context (letters from exact dicts, no sympy)"""
    d = pickle.load(open(path, "rb"))
    C = fmpq_mpoly_ctx.get(tuple(d["ctx_names"]))
    letters = [C.from_dict({tuple(e): fmpq(p, q) for e, p, q in lst}) for lst in d["alphabet_dicts"]]
    ts = {}
    for Dk, Wk, lst in d["ts"]:
        key = (tuple(tuple(x) for x in Dk), tuple(tuple((tuple(h[0]), tuple(tuple(y) for y in h[1])) for h in wd) for wd in Wk))
        ts[key] = C.from_dict({tuple(e): fmpq(p, q) for e, p, q in lst})
    return C, letters, ts, d.get("extra")

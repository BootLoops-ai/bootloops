"""Arbitrary-precision evaluator core for the LO collinear E4C (per channel, per S4
instance sigma, per Cheng-Wu chart), certified by Arb ball arithmetic + level doubling.

Route = (gauge g, analytic variable a, outer o, inner i):
  * chart x_g = 1; the sigma-instance integrand is assembled EXACTLY as ONE rational
    function  N(x_o,x_i,x_a) / prod_k A_k^{M_k}  (cells combined over the common
    denominator with exact fmpq_mpoly arithmetic; every atom A_k is linear in x_a);
  * innermost x_a-integral over (0,inf) done in CLOSED FORM:
        int_0^inf R dx = - sum_poles Res[ R(x) log(-x) ]   (keyhole; poles at x<0, log real),
    residues of arbitrary order by truncated power series in Arb (ball arithmetic; the
    node value carries a rigorous radius; catastrophic cancellation only widens balls);
  * the remaining 2-D integral by the tensor exp-sinh (double-exponential) trapezoid rule
    x = exp(c sinh s), two nested levels (h and 2h) from one grid, adaptive precision per
    node (precision doubled until the weighted node radius is below eps_abs).
"""
import os, sys, time, math, itertools
from fractions import Fraction as F
from flint import arb, arb_poly, arb_series, fmpq, fmpq_poly, fmpq_mpoly_ctx, ctx

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import oracle_cells as OC

PAIRS = OC.PAIRS


def shape_points(u, v):
    """shape = the four points {0, 1, u, v}; u, v complex with Fraction real/imag parts -> list of (re,im)."""
    return [(F(0), F(0)), (F(1), F(0)), (F(u[0]), F(u[1])), (F(v[0]), F(v[1]))]


def dmatrix(pts):
    d = {}
    for a, b in PAIRS:
        dx = pts[a - 1][0] - pts[b - 1][0]; dy = pts[a - 1][1] - pts[b - 1][1]
        d[(a, b)] = dx * dx + dy * dy
    return d


def dmatrix_split(u, ub, v, vb):
    """split-signature (complexified) shape in the frame w=(u,1,0,v), wb=(ub,1,0,vb):
    d_ab = (w_a-w_b)(wb_a-wb_b) with u,ub,v,vb independent rationals (all d_ab must be > 0 for a
    convergent Euclidean-like integral)."""
    W = {1: (F(u), F(ub)), 2: (F(1), F(1)), 3: (F(0), F(0)), 4: (F(v), F(vb))}
    return {(a, b): (W[a][0] - W[b][0]) * (W[a][1] - W[b][1]) for a, b in PAIRS}


def sigma_dmatrix(d, sigma):
    """labeled d-matrix for instance sigma (particle i at point sigma[i-1], points 1..4)."""
    out = {}
    for a, b in PAIRS:
        p, q = sigma[a - 1], sigma[b - 1]
        out[(a, b)] = d[(p, q) if p < q else (q, p)]
    return out


def sigma_classes(d):
    """group the 24 assignments by labeled d-matrix -> list of (sigma, mult, dsig)."""
    cl = {}
    for sig in itertools.permutations((1, 2, 3, 4)):
        ds = sigma_dmatrix(d, sig)
        key = tuple(ds[p] for p in PAIRS)
        if key not in cl:
            cl[key] = [sig, 0, ds]
        cl[key][1] += 1
    return [tuple(v) for v in cl.values()]


# ------------------------------------------------------------------ exact assembly
class InstanceRF:
    """Exact rational function of one (channel, chart g, sigma): N / prod atoms^M in the
    3 chart variables (ascending particle index)."""
    def __init__(self, cellres, dsig):
        self.g = cellres["gauge"]; self.xvars = cellres["xvars"]
        C3 = fmpq_mpoly_ctx.get(tuple("x%d" % j for j in self.xvars))
        self.C3 = C3
        gens = {j: C3.from_dict({tuple(1 if jj == j else 0 for jj in self.xvars): fmpq(1)}) for j in self.xvars}
        one = C3.from_dict({(0, 0, 0): fmpq(1)})
        X = {j: (one if j == self.g else gens[j]) for j in (1, 2, 3, 4)}
        self.X = X
        dq = {p: fmpq(dsig[p].numerator, dsig[p].denominator) for p in PAIRS}
        def atom_poly(nm):
            if nm == "Om":
                return X[1] + X[2] + X[3] + X[4]
            if nm[0] == "x":
                return X[int(nm[1])]
            if nm.startswith("ell"):
                S = [int(c) for c in nm[3:]]
                r = X[S[0]]
                for j in S[1:]:
                    r = r + X[j]
                return r
            if nm[0] == "s":
                S = [int(c) for c in nm[1:]]
                r = None
                for a, b in itertools.combinations(S, 2):
                    t = dq[(a, b)] * X[a] * X[b]
                    r = t if r is None else r + t
                return r
            raise ValueError(nm)
        # atoms and max powers
        M = {}
        for dk, conv in cellres["cells"]:
            for nm, pw in dk:
                M[nm] = max(M.get(nm, 0), pw)
        self.atom_names = sorted(M)
        self.M = M
        self.atom_polys = {nm: atom_poly(nm) for nm in self.atom_names}
        # power cache
        pc = {}
        def apow(nm, e):
            if (nm, e) not in pc:
                pc[(nm, e)] = self.atom_polys[nm] ** e if e > 0 else one
            return pc[(nm, e)]
        # d substitution
        dvals = [dq[p] for p in PAIRS]
        shiftinv = fmpq(1)
        for p in PAIRS:
            shiftinv /= dq[p] ** cellres["dshift"]
        N = C3.from_dict({(0, 0, 0): fmpq(0)})
        t0 = time.time()
        for dk, conv in cellres["cells"]:
            nd = {}
            for exx, p in conv.items():
                v = p(*dvals)
                if v != 0:
                    nd[exx] = v * shiftinv
            if not nd:
                continue
            Nc = C3.from_dict(nd)
            have = dict(dk)
            cof = one
            for nm in self.atom_names:
                e = M[nm] - have.get(nm, 0)
                if e:
                    cof = cof * apow(nm, e)
            N = N + Nc * cof
        self.N = N
        self.build_s = time.time() - t0
        self.degN = N.degrees()
        self.nN = len(N)


# ------------------------------------------------------------------ route structures
class RouteEval:
    """Numerical structures for one InstanceRF and roles (o,i,a) = particle indices."""
    def __init__(self, rf, o, i, a):
        assert sorted([o, i, a]) == rf.xvars
        self.rf = rf; self.o, self.i, self.a = o, i, a
        pos = {j: k for k, j in enumerate(rf.xvars)}
        po, pi, pa = pos[o], pos[i], pos[a]
        # numerator table: T[e_a][e_i] = fmpq_poly in x_o
        dN = rf.N.to_dict()
        tab = {}
        for exx, c in dN.items():
            tab.setdefault(exx[pa], {}).setdefault(exx[pi], {})[exx[po]] = c
        self.dega = max(tab) if tab else 0
        self.T = []
        for ea in range(self.dega + 1):
            row = tab.get(ea, {})
            degi = max(row) if row else -1
            polys = []
            for ei in range(degi + 1):
                cd = row.get(ei, {})
                dego = max(cd) if cd else -1
                polys.append(fmpq_poly([cd.get(eo, fmpq(0)) for eo in range(dego + 1)]))
            self.T.append(polys)
        # atoms: A = alpha(x_o,x_i) * x_a + beta(x_o,x_i); alpha/beta as {(eo,ei): fmpq}
        self.atoms = []
        for nm in rf.atom_names:
            d = rf.atom_polys[nm].to_dict()
            al, be = {}, {}
            for exx, c in d.items():
                assert exx[pa] <= 1, ("atom not linear in analytic var", nm)
                (al if exx[pa] == 1 else be)[(exx[po], exx[pi])] = c
            self.atoms.append({"name": nm, "M": rf.M[nm], "alpha": al, "beta": be,
                               "afree": len(al) == 0, "is_xa": (len(be) == 0)})
        # pole atoms (alpha != 0, beta != 0), the pure x_a atom (beta == 0), x_a-free atoms
        self.M0 = sum(A["M"] for A in self.atoms if A["is_xa"])
        self.poles = [A for A in self.atoms if (not A["afree"]) and (not A["is_xa"])]
        self.free = [A for A in self.atoms if A["afree"]]
        self.degD_rest = sum(A["M"] for A in self.poles)
        self._arbcache = {}

    # -------- per precision arb structures
    def _arb_T(self):
        key = ("T", ctx.prec)
        if key not in self._arbcache:
            self._arbcache[key] = [[arb_poly(p) for p in row] for row in self.T]
        return self._arbcache[key]

    @staticmethod
    def _eval2(dct, xo, xi):
        s = arb(0)
        for (eo, ei), c in dct.items():
            s += arb(c) * xo ** eo * xi ** ei
        return s

    def outer_setup(self, xo):
        """substitute x_o: returns per-e_a arb_poly in x_i, and atom alpha/beta as callables of x_i."""
        T = self._arb_T()
        Pi = [arb_poly([p(xo) for p in row]) for row in T]
        ats = []
        for A in self.atoms:
            # alpha, beta as arb_poly in x_i
            def topoly(dct):
                if not dct:
                    return None
                deg = max(ei for (eo, ei) in dct)
                cs = [arb(0)] * (deg + 1)
                cs = list(cs)
                for (eo, ei), c in dct.items():
                    cs[ei] = cs[ei] + arb(c) * xo ** eo
                return arb_poly(cs)
            ats.append((topoly(A["alpha"]), topoly(A["beta"])))
        return Pi, ats

    def node(self, setup, xi):
        """g(x_o, x_i) = int_0^inf N/D dx_a  as an arb ball (may be wide)."""
        Pi, ats = setup
        ncoef = [p(xi) for p in Pi]                      # coefficients of N in x_a
        # constant factor from x_a-free atoms, pole data
        C = arb(1)
        pole_ab = []
        for A, (alp, bep) in zip(self.atoms, ats):
            if A["afree"]:
                C *= bep(xi) ** A["M"]
            elif A["is_xa"]:
                pass
            else:
                pole_ab.append((alp(xi), bep(xi), A["M"]))
        # remove x^{M0}: low coefficients must vanish
        M0 = self.M0
        for k in range(min(M0, len(ncoef))):
            if not ncoef[k].contains(0):
                raise RuntimeError("x_a^-%d pole does not cancel: coeff %s" % (M0 - k, ncoef[k]))
        nt = ncoef[M0:]
        # degree bound: deg <= degD_rest - 2
        top = self.degD_rest - 2
        for k in range(top + 1, len(nt)):
            if not nt[k].contains(0):
                raise RuntimeError("degree excess at x_a^inf: coeff[%d] %s" % (k, nt[k]))
        nt = nt[:top + 1]
        # drop exact-zero information: replace must-vanish coefficients by their (zero-containing) balls? No:
        # they ARE dropped (set to 0) -- justified by the exact identities (integrability at generic fiber).
        Npoly = arb_poly(nt) if nt else arb_poly([arb(0)])
        total = arb(0)
        npole = len(pole_ab)
        for k in range(npole):
            al, be, M = pole_ab[k]
            p = -be / al
            # A(t) = N(p+t) truncated to M terms
            sh = Npoly(arb_poly([p, arb(1)]))
            cs = sh.coeffs() if hasattr(sh, "coeffs") else [sh[i] for i in range(sh.degree() + 1)]
            cs = list(cs)[:M] + [arb(0)] * max(0, M - len(cs))
            ser = arb_series(cs, M)
            for j in range(npole):
                if j == k:
                    continue
                al2, be2, M2 = pole_ab[j]
                fac = arb_series([al2 * p + be2, al2], M)
                ser = ser * (fac.inv() ** M2 if M2 > 1 else fac.inv())
            if M0:
                ser = ser * (arb_series([p, arb(1)], M).inv() ** M0)
            L = arb_series([-p, arb(-1)], M).log()
            ser = ser * L
            co = ser.coeffs()
            co = list(co) + [arb(0)] * (M - len(co))
            res = co[M - 1] / al ** M
            total += res
        return -total / C


# ------------------------------------------------------------------ 2-D DE quadrature
def de_node(hh, c, k, offset):
    """s = (k+offset) h ; x = exp(c sinh s) ; w = h c cosh(s) x   (arb at ctx.prec)."""
    off = arb(F(offset).numerator) / arb(F(offset).denominator) if offset else arb(0)
    s = hh * (arb(k) + off)
    x = (c * s.sinh()).exp()
    return x, hh * c * s.cosh() * x


def integrate_instance(rev, ho, hi, S, digits, prec0=384, precmax=12288, center_inner=True, So=None, Si=None, log=None, center_pow=F(1, 2)):
    """Tensor exp-sinh rule on the (x_o, x_i) quadrant:
         x_o = exp((pi/2) sinh((k+1/3) ho)),      |k ho| <= So
         x_i = xc(x_o) * exp((pi/2) sinh((k+1/2) hi)),  |.| <= Si,  xc = sqrt(x_o) if center_inner else 1
       (inner map centered on the geometric middle of the inner branch points, which live at
       O(1) and O(x_o) times distance ratios; the half offset excludes exact pole coincidences).
       Node acceptance: rad(w_o w_i g) <= 10^-(digits+6), precision doubled until met.
       Returns {'value': arb ball (quadrature sum incl. arithmetic radius), 'stats':...}."""
    guard = 6
    So = S if So is None else So; Si = S if Si is None else Si
    hoq, hiq = F(ho), F(hi)
    Ko = int(math.ceil(So / float(hoq))); Ki = int(math.ceil(Si / float(hiq)))
    stats = {"nodes": 0, "prec_hist": {}, "max_prec": prec0, "t_setup": 0.0, "t_nodes": 0.0, "Ko": Ko, "Ki": Ki}
    oldprec = ctx.prec
    total = arb(0)
    try:
        ctx.prec = prec0
        eps = arb(10) ** (-(digits + guard))
        for ko in range(-Ko - 1, Ko + 1):
            prec = prec0
            setups = {}
            row = arb(0)
            for ki in range(-Ki - 1, Ki + 1):
                # precision demand varies smoothly along a row: keep the neighbor's precision and
                # probe a halving every 6th node (a failed low-precision attempt costs a full node)
                if ki % 6 == 0:
                    prec = max(prec0, prec // 2)
                while True:
                    ctx.prec = prec
                    if prec not in setups:
                        t0 = time.time()
                        hho = arb(hoq.numerator) / arb(hoq.denominator); hhi = arb(hiq.numerator) / arb(hiq.denominator)
                        c = arb.pi() / 2
                        xo, wo = de_node(hho, c, ko, F(1, 3))
                        xc = (xo.sqrt() if center_pow == F(1, 2) else xo ** (arb(center_pow.numerator) / arb(center_pow.denominator))) if center_inner else arb(1)
                        setups[prec] = (rev.outer_setup(xo), hhi, c, wo, xc)
                        stats["t_setup"] += time.time() - t0
                    setup, hhi, c, wo, xc = setups[prec]
                    xi, wi = de_node(hhi, c, ki, F(1, 2))
                    xi = xi * xc; wi = wi * xc
                    t0 = time.time()
                    gval = None
                    try:
                        gval = rev.node(setup, xi)
                        ok = gval.is_finite() and (gval * wo * wi).rad() <= eps
                    except ZeroDivisionError:
                        ok = False
                    stats["t_nodes"] += time.time() - t0
                    if ok:
                        break
                    if prec >= precmax:
                        raise RuntimeError("node (%d,%d) failed at precmax %d: %s" % (ko, ki, prec, gval))
                    prec *= 2
                stats["nodes"] += 1
                stats["prec_hist"][prec] = stats["prec_hist"].get(prec, 0) + 1
                stats["max_prec"] = max(stats["max_prec"], prec)
                row += gval * wi
            ctx.prec = prec0
            hho = arb(hoq.numerator) / arb(hoq.denominator); c = arb.pi() / 2
            xo, wo = de_node(hho, c, ko, F(1, 3))
            total += row * wo
            if log and ko % 16 == 0:
                log("  row ko=%d done, nodes=%d maxprec=%d" % (ko, stats["nodes"], stats["max_prec"]))
        res = {"value": total, "ho": str(hoq), "hi": str(hiq), "So": So, "Si": Si, "center_inner": center_inner, "stats": stats}
    finally:
        ctx.prec = oldprec
    return res

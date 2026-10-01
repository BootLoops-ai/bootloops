"""Numeric evaluation of SLICE term stores (stage-2 one-fold objects over Q[x,t,j]/(j^2+1)) at rational t and a point x:
exact specialization t -> rational, j -> i (pairs of fmpq_poly A + i B in the last variable), then mpmath evaluation of
N/D * prod ZIP[words] with gate/scripts/hlog_eval.ZipEvalAlg (complex algebraic letters off the contour)."""
import os, sys
import slice_prov
from fractions import Fraction as Fr
import mpmath as mp
from flint import fmpq, fmpq_poly, fmpq_mpoly_ctx
import hlog_eval, hpath
ZERO_HL = ((0, 1), ())
class Specialized:
    """term store at fixed rational t: numerators/letters as (A, B) fmpq_poly pairs in x (value A(x) + i B(x))"""
    def __init__(self, C, letters, ts, tval, xname):
        self.names = list(C.names()); self.jx = self.names.index(xname); self.jt = self.names.index("t"); self.jj = self.names.index("j")
        tq = fmpq(Fr(tval).numerator, Fr(tval).denominator)
        C2 = fmpq_mpoly_ctx.get((xname, "j")); gx, gj = C2.gens(); J2 = gj * gj + C2.from_dict({(0, 0): fmpq(1)})
        imgs = []
        for k, n in enumerate(self.names):
            if k == self.jx: imgs.append(gx)
            elif k == self.jj: imgs.append(gj)
            elif k == self.jt: imgs.append(C2.from_dict({(0, 0): tq}))
            else: imgs.append(C2.from_dict({(0, 0): fmpq(0)}))
        def spec(P):
            P2 = P.compose(*imgs, ctx=C2)
            if int(P2.degrees()[1]) >= 2: P2 = divmod(P2, J2)[1]
            d = P2.to_dict(); A = {}; B = {}
            for e, c in d.items(): (A if int(e[1]) == 0 else B)[int(e[0])] = c
            def up(dd):
                if not dd: return fmpq_poly([])
                m = max(dd); return fmpq_poly([dd.get(i, fmpq(0)) for i in range(m + 1)])
            return up(A), up(B)
        self.letters = [spec(L) for L in letters]
        self.terms = [(Dk, Wk, spec(N)) for (Dk, Wk), N in ts.items()]
    @staticmethod
    def pval(AB, x):
        A, B = AB
        def ev(p):
            cs = p.coeffs(); tot = mp.mpf(0); xp = mp.mpf(1) if not isinstance(x, mp.mpc) else mp.mpc(1)
            acc = 0
            for c in reversed(cs): acc = acc * x + mp.mpf(int(c.p)) / int(c.q)
            return acc
        return mp.mpc(ev(A), ev(B))
    def hl_value(self, h, x):
        """value of a hyperlog letter at x; exact rational constants stay Fractions"""
        if h == ZERO_HL: return Fr(0)
        (p, q), Eh = h
        if not Eh: return Fr(p, q)
        v = mp.mpc(mp.mpf(p) / q)
        for i, e in Eh: v = v * self.pval(self.letters[i], x) ** e
        return v
    def value(self, x, dps=60, return_terms=False, rule="below"):
        """F(x) at real rational/complex x; on-contour letters (exactly real positive) get the universal 'letter - i0' prescription (rule='below')."""
        mp.mp.dps = dps; H = hpath.HPath(dps); tot = mp.mpc(0); lc = {}; mx = mp.mpf(0); self.n_oncontour = 0
        x = mp.mpf(Fr(x).numerator) / Fr(x).denominator if isinstance(x, (Fr, int)) else x
        tol = mp.mpf(10) ** (-(dps - 8))
        for Dk, Wk, NAB in self.terms:
            v = self.pval(NAB, x)
            if v == 0: continue
            for i, e in Dk:
                if i not in lc: lc[i] = self.pval(self.letters[i], x)
                v = v / lc[i] ** e
            for w in Wk:
                lets = []
                for h in w:
                    lv = self.hl_value(h, x)
                    if isinstance(lv, Fr):
                        lets.append((lv, bool(lv > 0))); self.n_oncontour += int(lv > 0); continue
                    if abs(lv.imag) < tol * (1 + abs(lv.real)):
                        re = lv.real
                        if re > 0:
                            if rule != "below": raise AssertionError(("on-contour ZIP letter", lv))
                            lets.append((mp.mpc(re, 0), True)); self.n_oncontour += 1
                        else: lets.append((mp.mpc(re, 0), False))
                    else: lets.append((lv, False))
                v = v * H.zip_value(lets)
            tot += v; mx = max(mx, abs(v))
        return (tot, mx) if return_terms else tot

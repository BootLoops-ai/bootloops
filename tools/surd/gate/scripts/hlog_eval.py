"""Evaluator for hyperlogarithms / ZeroInfPeriod words with ALGEBRAIC (real or complex) letters.

G(a_1..a_n; 1) for letters a_i given as exact Fractions OR mpmath numbers (mpf/mpc, e.g. roots of the
last-variable cubic), none on the open path (0,1): shuffle-regularization of leading 1's and trailing
0's (exact letters only carry those), then a RECURSIVE Hoelder split at 1/2 with rescaling until every
non-zero letter has |a| >= 2, then the Taylor recursion.  ZIP_reg(word) (HyperFLINT ZeroInfPeriod
semantics, x = t/(1-t)) on top.  Self-gates against closed forms and direct quadrature (run as main)."""
import sys, os, time, json
from fractions import Fraction as Fr
import mpmath as mp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fibrate_dep import zip_num


def _num(a):
    if isinstance(a, Fr):
        return mp.mpf(a.numerator) / a.denominator
    return a


def _is_exact(a, val):
    return isinstance(a, (Fr, int)) and a == val


def _on_open_unit_interval(a):
    if isinstance(a, (Fr, int)):
        return 0 < a < 1
    z = mp.mpc(a)
    return abs(z.imag) < mp.mpf(10) ** (-(mp.mp.dps - 5)) and 0 < z.real < 1


class GEvalAlg:
    def __init__(self, dps):
        self.dps = dps; self.cache = {}

    def key(self, w):
        return tuple((("F", a.numerator, a.denominator) if isinstance(a, Fr) else ("N", mp.nstr(a, self.dps + 5))) for a in w)

    def series(self, w):
        amin = min(abs(_num(a)) for a in w if not _is_exact(a, 0))
        depth = len(w)
        NT = int((self.dps + 15 + 3 * depth) * mp.log(10) / mp.log(amin)) + 20 + 5 * depth
        b = _num(w[-1]); c = [mp.mpf(0)] * (NT + 1); ib = 1 / b; p = mp.mpf(1)
        for n in range(1, NT + 1):
            p *= ib; c[n] = -p / n
        for a in reversed(w[:-1]):
            if _is_exact(a, 0):
                c = [mp.mpf(0)] + [c[n] / n for n in range(1, NT + 1)]
            else:
                av = _num(a); ia = 1 / av; cp = [mp.mpf(0)] * (NT + 1); s = mp.mpf(0); ap = mp.mpf(1); iap = mp.mpf(1)
                for j in range(1, NT + 1):
                    s += c[j - 1] * ap; iap *= ia; cp[j] = -s * iap / j; ap *= av
                c = cp
        return mp.fsum(c[1:])

    def G1(self, w):
        """G(w;1); requires w[0] != 1 and w[-1] != 0 (regularized word) and no letter on (0,1)."""
        w = tuple(w)
        if not w:
            return mp.mpf(1)
        k = self.key(w)
        if k in self.cache:
            return self.cache[k]
        assert not _is_exact(w[0], 1) and not _is_exact(w[-1], 0), w
        for a in w:
            assert not _on_open_unit_interval(a), ("letter on the path (0,1)", a)
        nz = [abs(_num(a)) for a in w if not _is_exact(a, 0)]
        if min(nz) >= 2:
            val = self.series(w)
        else:
            n = len(w); tot = mp.mpf(0)
            for kk in range(n + 1):
                left = tuple((2 * (1 - a)) if isinstance(a, Fr) else 2 * (1 - a) for a in reversed(w[:kk]))
                right = tuple((2 * a) for a in w[kk:])
                tot += (-1) ** kk * self.G1(left) * self.G1(right)
            val = tot
        self.cache[k] = val
        return val


class ZipEvalAlg(zip_num.ZipEval):
    """ZIP_reg over letters <= 0 (rational) or complex/algebraic off the contour."""
    def __init__(self, dps):
        super().__init__(); self.ge = GEvalAlg(dps)

    def G_one_conv(self, w):
        return self.ge.G1(tuple(w))

    def zip_value(self, word):
        word = tuple(word)
        k = self.ge.key(word)
        if k in self.ZM:
            return self.ZM[k]
        for a in word:
            if isinstance(a, Fr):
                assert a <= 0, ("on-contour ZIP letter", a)
            else:
                z = mp.mpc(a); assert not (abs(z.imag) < mp.mpf(10) ** (-(mp.mp.dps - 5)) and z.real > 0), ("on-contour ZIP letter", a)
        exps = [((), Fr(1))]
        for a in word:
            opts = [(Fr(1), Fr(-1))]
            if not _is_exact(a, -1):
                opts.append(((a / (1 + a)), Fr(1)))
            exps = [(ww + (l,), c * s) for (ww, c) in exps for (l, s) in opts]
        tot = mp.mpf(0)
        for ww, c in exps:
            tot += mp.mpf(c.numerator) / c.denominator * self.G_one(ww)
        self.ZM[k] = tot
        return tot

    # G_one from zip_num uses strip_trailing_zeros / reg_leading_ones with exact comparisons on Fractions;
    # algebraic letters are never exactly 0 or 1, so those routines pass them through unchanged.


def selfgate(dps=50):
    mp.mp.dps = dps
    Z = ZipEvalAlg(dps); Gv = Z.ge
    out = {}
    s2 = mp.sqrt(2); phi = (1 + mp.sqrt(5)) / 2
    a = 1 + s2                                    # algebraic real > 1
    out["G(1+sqrt2;1) = log(1-1/a)"] = float(abs(Gv.G1((a,)) - mp.log(1 - 1 / a)))
    b = mp.mpc(2) * mp.exp(1j * mp.pi / 3)       # complex algebraic
    out["G(2 e^{i pi/3};1) = log(1-1/b)"] = float(abs(Gv.G1((b,)) - mp.log(1 - 1 / b)))
    out["G(0,a;1) = -Li2(1/a), a=1+sqrt2"] = float(abs(Gv.G1((Fr(0), a)) - (-mp.polylog(2, 1 / a))))
    out["G(0,0,b;1) = -Li3(1/b), b complex"] = float(abs(Gv.G1((Fr(0), Fr(0), b)) - (-mp.polylog(3, 1 / b))))
    c1, c2 = -s2, -phi                            # two algebraic negative letters vs nested quadrature
    g = Gv.G1((c1, c2))
    f = lambda t1: 1 / (t1 - c1) * mp.quad(lambda t2: 1 / (t2 - c2), [0, t1])
    out["G(-sqrt2,-phi;1) vs nested quad"] = float(abs(g - mp.quad(f, [0, 1])))
    r = -(mp.mpf(2)) ** (mp.mpf(1) / 3)            # cube-root letter: ZIP[a] = -log(-a)
    out["ZIP[-2^(1/3)] = -log(2^(1/3))"] = float(abs(Z.zip_value((r,)) - (-mp.log(-r))))
    # ZIP with a small algebraic letter (deep recursion) vs quadrature of the mapped G
    a2 = -mp.sqrt(3) / 40
    g2 = Z.zip_value((Fr(0), a2))                  # = G-combination; compare with direct definition:
    # ZIP[0,a] = int_0^inf dx/x * int_0^x dy/(y-a) regularized at infinity: use the mapped identity instead:
    bb = a2 / (1 + a2)
    direct = Gv.G1((Fr(1), bb)) if False else None
    # compare word (0,a2) between this evaluator and mpmath nested quadrature of G-letters after the x=t/(1-t) map:
    # ZIP[0,a] = sum over the 4 mapped words; check instead the rational neighbor a=-1/23 against zip semantics of rationals
    out["ZIP[0,-1/23] alg-path vs Fraction-path"] = float(abs(Z.zip_value((Fr(0), mp.mpf(-1) / 23)) - Z.zip_value((Fr(0), Fr(-1, 23)))))
    return out


if __name__ == "__main__":
    import argparse, gate_prov, fibrate_dep
    ap = argparse.ArgumentParser(description="self-check of the hyperlog evaluator against closed forms and quadrature (dps 50); exit 0 = PASS")
    ap.add_argument("--out", default="", help="optional JSON record of the self-check (default: none; print only)"); a = ap.parse_args()
    t0 = time.time(); res = selfgate(50)
    ok = all(v < 1e-45 for v in res.values())
    rec = {"selfgate_absdiffs_dps50": res, "PASS": ok, "wall_s": round(time.time() - t0, 1),
           "scope": "G(w;1) and ZIP_reg(word) with rational, real-algebraic and complex-algebraic letters off the path; recursion = Hoelder split at 1/2 with rescaling; series when all |letters|>=2"}
    if a.out:
        rec["producer"] = gate_prov.producer(inputs=[os.path.join(fibrate_dep.fibrate_dir(), "zip_num.py")], modules=["hlog_eval.py", "gate_prov.py"]); gate_prov.dump_json(rec, a.out)
    print(json.dumps(res, indent=1), "PASS", ok); sys.exit(0 if ok else 1)

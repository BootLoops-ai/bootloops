"""ellipticus.api — front door: MomentFamily (engine i), evaluate() with the
built-in two-precision self-gate.

The certified word classes (v1, scope certificates — anything else REFUSES
loudly rather than answering wrong):
  march engine (workhorse): words  [dlog a_w, ..., dlog a_2, mom k]  on a
    declared family y^2 = F(x;z) with exact rational window — the inner
    kernel is the weight-2 incomplete moment I_k(z), each outer dlog letter
    adds one weight.  a_j exact rational (or exact algebraic root supplied
    at working precision); a_j = 0 allowed for the OUTERMOST letter only
    (single-log tangential-base regularization; deeper log letters need the
    full shuffle L-algebra — use the tau-word engine for those).
  q engine: Kronecker tau-words at a declared torus frame (vendored itint
    shuffle-regularized engine) + the BMSW sunrise closed forms.
  AGM verb: weight-1/period layer (curve.periods / cycle_moments).
"""
import json
from fractions import Fraction
import sympy as sp
import mpmath as mp

from .curve import exact, QuarticCurve, sunrise_frame
from . import extend
from .march import (System, letter_base_log, letter_base_r, seed_from_series)


class MomentFamily:
    """y^2 = F(x;z) family with exact rational window [x1,x2]:
    the engine-(i) object.  coeffs_x_desc: descending-in-x coefficient list,
    each entry a string/sympy expression in the parameter (exact rationals
    only).  param: parameter symbol name."""

    def __init__(self, coeffs_x_desc, window, kmax=None, param='z'):
        self.z = sp.Symbol(param)
        x = sp.Symbol('_x_mf')
        self.x = x
        cs = [sp.sympify(c, locals={param: self.z}) for c in coeffs_x_desc]
        for c in cs:
            if c.free_symbols - {self.z}:
                raise ValueError(f"foreign symbol in coefficient {c}")
        deg = len(cs) - 1
        if deg not in (3, 4):
            raise ValueError("family must be cubic or quartic in x")
        self.F = sp.Poly(sum(c * x ** (deg - k) for k, c in enumerate(cs)), x)
        self.X1 = sp.Rational(exact(window[0]).numerator,
                              exact(window[0]).denominator)
        self.X2 = sp.Rational(exact(window[1]).numerator,
                              exact(window[1]).denominator)
        # dim H^1(affine curve) = 2g + (#points at infinity - 1)
        self.kmax = kmax if kmax is not None else (2 if deg == 4 else 1)
        self._pf = None

    def polyform(self):
        """exact extended system, common-denominator polynomial form
        (cached per instance)."""
        if self._pf is None:
            M = extend.moment_system(self.F, self.kmax, self.x, self.z,
                                     self.X1, self.X2)
            self._pf = extend.polyform(M, self.z)
        return self._pf

    def seed_series(self, zbase, N):
        """exact local z-series of (I_0..I_kmax) at an analytic (perfect-
        square or square-x-linear fibre) base point; ValueError otherwise.
        Disk-cached under $ELLIPTICUS_CACHE (former name $EMPL_EVAL_CACHE
        also read) if set (keyed by family sha,
        zbase, N, dps — serializer-width policy: strings at working dps)."""
        import os, json, hashlib
        from . import env
        cdir = env("ELLIPTICUS_CACHE", "EMPL_EVAL_CACHE")
        key = None
        if cdir:
            h = hashlib.sha256(
                (str(self.F.as_expr()) + f"|{self.X1}|{self.X2}|{zbase}|"
                 f"{N}|{mp.mp.dps}").encode()).hexdigest()[:16]
            key = os.path.join(cdir, f"empl_seed_{h}.json")
            if os.path.exists(key):
                raw = json.load(open(key))
                return {int(k): [mp.mpf(s) for s in raw[k]] for k in raw}
        out = extend.analytic_seed_series(self.F, self.x, self.z, zbase,
                                          Fraction(self.X1.p, self.X1.q),
                                          Fraction(self.X2.p, self.X2.q),
                                          self.kmax, N)
        if key:
            json.dump({str(k): [mp.nstr(v, mp.mp.dps) for v in out[k]]
                       for k in out}, open(key, "w"))
        return out

    def boundary_states(self, zv):
        """exact boundary components B_i = F(x_i;z)^(-1/2) at mp precision."""
        out = []
        for XI in (self.X1, self.X2):
            Fx = self.F.as_expr().subs({self.x: XI})
            num, den = sp.fraction(sp.together(sp.cancel(Fx)))
            pz = lambda p: [Fraction(sp.Rational(c).p, sp.Rational(c).q)
                            for c in sp.Poly(p, self.z).all_coeffs()]
            def ev(coefs):
                s = mp.mpc(0)
                for c in coefs:
                    s = s * zv + mp.mpf(c.numerator) / mp.mpf(c.denominator)
                return s
            out.append(1 / mp.sqrt(ev(pz(num)) / ev(pz(den))))
        return out

    # ------------------------------------------------------------- evaluate
    def eval_words(self, words, endpoints, dps, zbase=0, z1=None,
                   waypoints=None, wdps_extra=80):
        """March-engine evaluation of a list of words at each endpoint.
        words: [ [a_w, ..., a_2, ('mom', k)] ] — outer dlog letters a_j
        (exact rationals or mp values from exact data), innermost the moment
        kernel.  endpoints: exact rationals (transported in given order).
        Returns dict endpoint -> {'kernels': [...], 'words': [...]} at
        working precision wdps = dps + wdps_extra.  Certification: per-step
        proven-tail asserts (fail-closed) + the caller-level two-precision
        gate (ellipticus.evaluate)."""
        zbase = exact(zbase, 'zbase')
        if zbase != 0:
            raise NotImplementedError(
                "v1: analytic seed base at z=0 chart only; shift your "
                "family so the seed fibre sits at z=0")
        if z1 is None:
            z1 = Fraction(1, 50)
        z1 = exact(z1, 'z1')
        pf = self.polyform()
        wdps = dps + wdps_extra
        with mp.workdps(wdps):
            NB = int(wdps / 1.7) + 30
            S0 = self.seed_series(0, NB)
            z1m = mp.mpf(z1.numerator) / mp.mpf(z1.denominator)
            Iseed = seed_from_series(S0, z1m, NB)
            Bseed = self.boundary_states(z1m)
            # letters: build the System letter list + seed values
            letters = []
            lseed = []
            n = self.kmax + 3
            for w in words:
                inner = w[-1]
                assert isinstance(inner, tuple) and inner[0] == 'mom', \
                    "innermost entry must be ('mom', k)"
                k = int(inner[1])
                assert 0 <= k <= self.kmax
                src = k
                depth = len(w) - 1
                for pos in range(depth - 1, -1, -1):
                    a = w[pos]
                    a_ex = a if isinstance(a, (mp.mpf, mp.mpc)) \
                        else mp.mpf(exact(a, 'letter').numerator) \
                        / mp.mpf(exact(a, 'letter').denominator)
                    if a_ex == 0 and pos != 0:
                        raise NotImplementedError(
                            "v1: dlog(z) letter allowed OUTERMOST only "
                            "(single-log reg); use the tau-word engine")
                    letters.append((a_ex, src))
                    # seed value of this letter state at z1
                    if src < n:      # source is the kernel I_k
                        ser = S0[k]
                        if a_ex == 0:
                            lseed.append(letter_base_log(ser[0], ser, z1m,
                                                         NB))
                        else:
                            lseed.append(letter_base_r(ser, z1m, a_ex, NB))
                    else:
                        raise NotImplementedError(
                            "v1: nested letter seeds beyond depth-1 need "
                            "letter series at the base — not yet wired; "
                            "single outer letter per word for now")
                    src = n + len(letters) - 1
            sysm = System(pf, letters=letters, dps=wdps)
            Y = list(Iseed) + list(Bseed) + lseed
            out = {}
            zc = z1m
            for ep in endpoints:
                epf = exact(ep, 'endpoint')
                zt = mp.mpf(epf.numerator) / mp.mpf(epf.denominator)
                if waypoints:
                    Y, _ = sysm.transport_path(zc, Y, list(waypoints) + [zt])
                    waypoints = None
                else:
                    Y, _ = sysm.transport(zc, Y, zt)
                zc = zt
                out[str(epf)] = dict(
                    kernels=[mp.mpc(v) for v in Y[:self.kmax + 1]],
                    boundary=[mp.mpc(v) for v in
                              Y[self.kmax + 1:self.kmax + 3]],
                    words=[mp.mpc(v) for v in Y[self.kmax + 3:]])
            return out


# --------------------------------------------------------------------------
def shared_digits(a, b):
    """decimal digits shared by two mp values (relative)."""
    a, b = mp.mpc(a), mp.mpc(b)
    d = abs(a - b) / max(abs(b), mp.mpf(10) ** (-mp.mp.dps))
    return float(-mp.log10(d)) if d > 0 else float(mp.mp.dps)


def evaluate(verb, dps, extra=20, **kw):
    """Front door with the built-in TWO-PRECISION SELF-GATE: runs the verb
    at dps and dps+extra, returns dict(value=..., shared_digits=...,
    lo=..., hi=...).  verb: callable taking dps -> mp value (or list of).
    Values are compared entrywise; shared_digits = worst entry."""
    lo = verb(dps, **kw) if kw else verb(dps)
    hi = verb(dps + extra, **kw) if kw else verb(dps + extra)
    with mp.workdps(dps + extra + 20):
        if isinstance(lo, (list, tuple)):
            sd = min(shared_digits(a, b) for a, b in zip(lo, hi))
        else:
            sd = shared_digits(lo, hi)
    return dict(value=hi, shared_digits=sd, lo=lo, hi=hi)

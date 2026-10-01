#!/usr/bin/env python3
"""fibrate — SOUND fibration engine for ZIP_reg words with ze-dependent
letters -> G({0,1,...}; ze) x level-2 constants.

WHY THIS TOOL EXISTS: in our v1.2.8 build HyperFLINT's `fibration_basis`
returned incorrect values on regulator-divergent / on-contour ZIP words
(631/635 split-job fibrations in one large production run, while all 635
step-1 regulators were exactly right).  Upstream's own fibration route
shells out to HyperInt's fibrationBasis (reference/SubTropica.wl:11411,
:15565 STToFibrationBasis) without per-word validation.  This engine is an
independently validated alternative for that case.

METHOD (370/370 words fibrated + validated in the production
pilot): Goncharov total differential on the regularized iterated
integral I(0; b1..bn; oo) (b = reversed ZIP word) — differentiation peels
letters via dlog forms (paper eqs 3.40/3.41; total-differential remark
paper:1423-1426), recursive EXACT sympy integration in ze with pole
letters, integration constant fixed numerically at ze=1/2 and identified
by PSLQ over the level-2 weight<=5 constant basis.

HARDENING:
 (1) ZIP series lengths are ADAPTIVE per word (fixed-length dps-ceiling
     class closed in zip_num.py; --selfcheck verifies it).
 (2) The PSLQ integration-constant step carries its FULL control battery:
     dual-precision stability, height cap, POSITIVE control (run at every
     Engine construction; engine refuses to run if it fails), residual
     gate.  Provenance "pslq_constants" is STAMPED into the output for
     every word whose representation (transitively) contains a
     PSLQ-identified nonzero constant, so downstream LaurentSeries can
     set provenance honestly (never silently :analytic).
 (3) Per-word 2-point numeric validation against the TRUSTED step-1 ZIP
     semantics (zip_num.ZipEval) is BAKED IN — every fibrated word is
     checked at >=2 kinematic points distinct from the fit point; there
     is deliberately NO flag to disable it.

CLI:  python3 fibrate.py WORDLIST.json [--out OUT.json] [--dps N]
                         [--validation-points P1 P2 ...]
      python3 fibrate.py --selfcheck
Wordlist JSON = list of words; word = list of letter strings in `ze`
(e.g. ["-1", "-ze", "1/(ze - 1)", "0"]), same format as the pilot's
unres_words.json.  See README.md in this directory for the API.
"""
import argparse
import json
import os
import sys
from fractions import Fraction as Fr

import mpmath as mp
import sympy as sp

try:                                    # package import (scripts.fibrate.*)
    from .zip_num import ZipEval
    from .gpl_num import GEval
except ImportError:                     # flat script / sys.path import
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from zip_num import ZipEval
    from gpl_num import GEval

ze = sp.Symbol('ze')

TOOL_NAME = "subtropica/fibrate"
TOOL_VERSION = "1.0.0"
SOURCE_LINEAGE = ("F^red1-validated fibration engine "
                  "(370/370 words validated)")

HEIGHT_CAP = 10 ** 10      # PSLQ maxcoeff (sane height, integrity rule)
PSLQ_MAXSTEPS = 600000


class FibrateError(RuntimeError):
    """Loud typed failure: validation / PSLQ / input refusal."""


# ---------------------------------------------------------------------------
# Level-2 constant basis (weight <= 5).  Built at REQUESTED precision —
# NEVER at import-time ambient dps (mplll-ambient-dps footgun class).
# ---------------------------------------------------------------------------
def build_consts(dps):
    """[(name, mpf value at dps+30, weight)] — level-2 weight<=5 basis."""
    with mp.workdps(dps + 30):
        l2 = mp.log(2)
        z2 = mp.zeta(2)
        z3 = mp.zeta(3)
        z4 = mp.zeta(4)
        z5 = mp.zeta(5)
        li4 = mp.polylog(4, mp.mpf(1) / 2)
        li5 = mp.polylog(5, mp.mpf(1) / 2)
        return [
            ('one', mp.mpf(1), 0),
            ('l2', l2, 1),
            ('z2', z2, 2), ('l2p2', l2 ** 2, 2),
            ('z3', z3, 3), ('z2l2', z2 * l2, 3), ('l2p3', l2 ** 3, 3),
            ('z4', z4, 4), ('li4', li4, 4), ('z3l2', z3 * l2, 4),
            ('z2l2p2', z2 * l2 ** 2, 4), ('l2p4', l2 ** 4, 4),
            ('z5', z5, 5), ('z4l2', z4 * l2, 5), ('z2z3', z2 * z3, 5),
            ('li4l2', li4 * l2, 5), ('li5', li5, 5), ('z3l2p2', z3 * l2 ** 2, 5),
            ('z2l2p3', z2 * l2 ** 3, 5), ('l2p5', l2 ** 5, 5),
        ]


ONE = ()          # the empty constant multiset (rational coefficient)


def cmul(k1, k2):
    """Product of two constant multisets (sorted name tuples)."""
    return tuple(sorted(k1 + k2))


# ---------------------------------------------------------------------------
# Letter parsing / representation helpers (exact sympy)
# rep: dict {(gword tuple of sympy numbers, const name tuple): ratfn in ze}
# ---------------------------------------------------------------------------
def parse_letter(s):
    return sp.nsimplify(sp.sympify(s.replace('^', '**'), locals={'ze': ze}),
                        rational=False)


def rep_add(r, key, coef):
    if coef == 0:
        return
    r[key] = sp.cancel(r.get(key, sp.Integer(0)) + coef)
    if r[key] == 0:
        del r[key]


def rep_scale(r, rho):
    out = {}
    for k, c in r.items():
        cc = sp.cancel(rho * c)
        if cc != 0:
            out[k] = cc
    return out


# ---------------------------------------------------------------------------
# GEvalX — gpl_num.GEval generalized to arbitrary exact sympy letters with
# |letter| > |x| (validation-side evaluator for the target G(...; ze) basis).
# ---------------------------------------------------------------------------
class GEvalX:
    def __init__(self, x_fr):
        self.x = mp.mpf(x_fr.numerator) / x_fr.denominator
        self.logx = mp.log(self.x)
        self.N = int((mp.mp.dps + 12) * mp.log(10) / (-mp.log(self.x))) + 10
        self.memo = {}

    def series(self, w):
        if w in self.memo:
            return self.memo[w]
        N = self.N
        b = mp.mpf(str(sp.N(w[-1], mp.mp.dps + 10)))
        assert abs(b) > self.x, ('letter inside radius', w[-1], self.x)
        c = [mp.mpf(0)] * (N + 1)
        ib = 1 / b
        p = mp.mpf(1)
        for n in range(1, N + 1):
            p *= ib
            c[n] = -p / n
        for a in reversed(w[:-1]):
            if a == 0:
                c = [mp.mpf(0)] + [c[n] / n for n in range(1, N + 1)]
            else:
                av = mp.mpf(str(sp.N(a, mp.mp.dps + 10)))
                assert abs(av) > self.x, ('letter inside radius', a, self.x)
                ia = 1 / av
                cp = [mp.mpf(0)] * (N + 1)
                s = mp.mpf(0)
                ap = mp.mpf(1)
                iap = mp.mpf(1)
                for j in range(1, N + 1):
                    s += c[j - 1] * ap
                    iap *= ia
                    cp[j] = -s * iap / j
                    ap *= av
                c = cp
        self.memo[w] = c
        return c

    def G(self, w):
        """G(w; x), w = tuple of exact sympy numbers; trailing zeros via
        shuffle-reg strip (log(x)^k pieces)."""
        w = tuple(sp.nsimplify(a) for a in w)

        def strip(u):
            if not u or u[-1] != 0:
                return {(0, u): Fr(1)}
            if all(a == 0 for a in u):
                k = len(u)
                f = Fr(1)
                for j in range(1, k + 1):
                    f /= j
                return {(k, ()): f}
            v = u[:-1]
            ins = {}
            for i in range(len(v) + 1):
                cand = tuple(v[:i] + (sp.Integer(0),) + v[i:])
                ins[cand] = ins.get(cand, 0) + 1
            t = ins.pop(u)
            out = {}
            for (k, vv), c in strip(v).items():
                out[(k + 1, vv)] = out.get((k + 1, vv), Fr(0)) + c / t
            for cand, m in ins.items():
                for (k, vv), c in strip(cand).items():
                    out[(k, vv)] = out.get((k, vv), Fr(0)) - Fr(m) * c / t
            return out

        tot = mp.mpf(0)
        for (k, v), q in strip(w).items():
            if v:
                c = self.series(v)
                val = mp.mpf(0)
                xp = mp.mpf(1)
                for n in range(1, self.N + 1):
                    xp *= self.x
                    val += c[n] * xp
            else:
                val = mp.mpf(1)
            tot += (mp.mpf(q.numerator) / q.denominator) * self.logx ** k * val
        return tot


# ---------------------------------------------------------------------------
# Exact integration of ratfn(ze) * G(u; ze)  (recursive, sympy).
# Antiderivative construction = paper §3.3.2 step 2 ("integrating back"
# eq 3.41); pole terms prepend letters, polynomial terms integrate by parts.
# ---------------------------------------------------------------------------
def integrate_term(coef, u, out):
    """out += antiderivative of coef(ze)*G(u; ze).  Exact.  Refuses loudly
    on algebraic (non-linear irreducible) pole factors."""
    coef = sp.cancel(sp.together(coef))
    num, den = sp.fraction(coef)
    num = sp.Poly(num, ze)
    den = sp.Poly(den, ze)
    if den.degree() == 0:
        poly = (num / den.as_expr()).as_poly(ze)
        _int_poly_times_G(poly, u, out)
        return
    # loud refusal guard: only rational (linear-factor) poles are in scope
    for base, _mult in sp.factor_list(den.as_expr())[1]:
        pb = sp.Poly(base, ze)
        if pb.degree() > 1:
            raise FibrateError(
                f"nonlinear irreducible denominator factor (algebraic pole): "
                f"{base}")
    ap = sp.apart(coef, ze, full=False)
    for t in sp.Add.make_args(ap):
        n2, d2 = sp.fraction(sp.together(t))
        pd = sp.Poly(d2, ze)
        if pd.degree() == 0:
            _int_poly_times_G((n2 / d2).as_poly(ze), u, out)
            continue
        # match c/(ze-beta)^m
        fl = sp.factor_list(d2)
        assert len(fl[1]) == 1 and sp.Poly(fl[1][0][0], ze).degree() == 1, t
        base, m = fl[1][0]
        pb = sp.Poly(base, ze)
        a1, a0 = pb.all_coeffs()
        beta = sp.nsimplify(-a0 / a1)
        cnum = sp.cancel(n2 / (fl[0] * a1 ** m))
        pn = sp.Poly(cnum, ze)
        assert pn.degree() == 0, ('numerator not constant after apart', t)
        c = pn.all_coeffs()[0]
        if m == 1:
            rep_add(out, ((beta,) + u, ONE), c)     # prepend pole letter
        else:
            _int_pole_m(c, beta, m, u, out)


def _int_pole_m(c, beta, m, u, out):
    """int c*G(u;ze)/(ze-beta)^m by parts (m >= 2 reduces to m-1)."""
    if m == 1:
        rep_add(out, ((beta,) + u, ONE), c)
        return
    rep_add(out, (u, ONE), -c / ((m - 1) * (ze - beta) ** (m - 1)))
    if u:
        u1 = u[0]
        integrate_term(
            sp.together(sp.Integer(1) * c /
                        ((m - 1) * (ze - u1) * (ze - beta) ** (m - 1))),
            u[1:], out)


def _int_poly_times_G(poly, u, out):
    if poly.is_zero:
        return
    if not u:
        prim = poly.integrate(ze)
        rep_add(out, ((), ONE), prim.as_expr())
        return
    # integrate by parts, highest degree first:
    # int a*ze^k G(u) = a*ze^{k+1}/(k+1) G(u)
    #                   - a/(k+1) int ze^{k+1}/(ze-u1) G(u[1:])
    k = poly.degree()
    a = poly.nth(k)
    rep_add(out, (u, ONE), a * ze ** (k + 1) / (k + 1))
    integrate_term(-a * ze ** (k + 1) / ((k + 1) * (ze - u[0])), u[1:], out)
    rest = poly - sp.Poly(a * ze ** k, ze)
    _int_poly_times_G(rest, u, out)


# ---------------------------------------------------------------------------
# The engine
# ---------------------------------------------------------------------------
class Engine:
    """Sound fibration engine.  Constructing an Engine SETS the ambient
    mpmath precision to `dps` (single-purpose process model, as in the
    pilot) and RUNS the PSLQ positive control — construction fails loudly
    if the control battery is broken.

    val_points: >=2 distinct rational points in (0,1), all != the fit
    point 1/2.  Validation is mandatory; there is no off switch.
    """

    FIT_POINT = Fr(1, 2)

    def __init__(self, dps=130, val_points=(Fr(2, 5), Fr(3, 8))):
        if dps < 60:
            raise FibrateError(f"dps={dps} < 60: too low for the PSLQ "
                               "control battery")
        self.dps = dps
        self.dps2 = max(50, dps - 40)          # dual-precision partner
        assert self.dps2 != self.dps
        pts = [Fr(p) for p in val_points]
        if len(set(pts)) < 2:
            raise FibrateError("validation requires >= 2 DISTINCT points "
                               f"(got {val_points}) — validation is "
                               "mandatory, not tunable-away")
        for p in pts:
            if not (Fr(0) < p < Fr(1)) or p == self.FIT_POINT:
                raise FibrateError(f"validation point {p} invalid: need "
                                   "0 < p < 1 and p != fit point 1/2")
        self.val_points = pts
        self.val_tol_exp = dps - 40            # rel tol 10^-(dps-40)
        self.zero_tol_exp = dps - 25           # |C| below => constant is 0
        self.resid_tol_exp = dps - 30          # PSLQ residual gate

        mp.mp.dps = dps
        self.consts = build_consts(dps)
        self.cval = {nm: v for nm, v, _wt in self.consts}
        self.pslq_report = self._positive_control()

        self.zip12 = ZipEval()
        self.g12 = GEvalX(self.FIT_POINT)
        self.val_ctx = [(p, ZipEval(), GEvalX(p)) for p in pts]
        self.memo = {}
        self.meta = {}     # word -> {'used_pslq': bool, 'digits': {...}}

    # -- constants ----------------------------------------------------------
    def cval_tuple(self, k):
        v = mp.mpf(1)
        for nm in k:
            v *= self.cval[nm]
        return v

    # -- PSLQ with full control battery --------------------------------------
    def identify_const(self, x, maxweight=5):
        """Identify x as a rational combination of basis constants of
        weight <= maxweight.  Controls: dual-precision stability
        (dps, dps2), height cap (HEIGHT_CAP), residual gate
        10^-(dps-30).  Returns {} for numeric zero, dict name->Fraction
        on success, None on failure (caller must treat None as FATAL)."""
        cols = [c for c in self.consts if c[2] <= maxweight]
        if abs(x) < mp.mpf(10) ** (-self.zero_tol_exp):
            return {}
        rels = []
        for dps_try in (self.dps, self.dps2):
            with mp.workdps(dps_try):
                vec = [mp.mpf(str(x))] + [-mp.mpf(str(c[1])) for c in cols]
                rel = mp.pslq(vec, maxcoeff=HEIGHT_CAP,
                              maxsteps=PSLQ_MAXSTEPS)
            if rel is None:
                return None
            rels.append(rel)
        if rels[0] != rels[1]:
            return None                     # two-precision stability FAIL
        rel0 = rels[0]
        if rel0[0] == 0:
            return None                     # relation not involving x
        out = {}
        for c, (nm, _v, _wt) in zip(rel0[1:], cols):
            if c:
                out[nm] = Fr(c, rel0[0])
        resid = x - sum(mp.mpf(q.numerator) / q.denominator * self.cval[nm]
                        for nm, q in out.items())
        if abs(resid) > mp.mpf(10) ** (-self.resid_tol_exp):
            return None                     # residual gate FAIL
        return out

    def _positive_control(self):
        """Known nontrivial combination MUST be recovered exactly by the
        same identify_const path used in production (integrity rule:
        PSLQ needs a positive control).  Raises FibrateError on failure."""
        truth = {'z4': Fr(-7, 4), 'li4': Fr(3), 'z3l2': Fr(1, 2),
                 'l2p4': Fr(-5, 24)}
        x = mp.mpf(0)
        for nm, q in truth.items():
            x += mp.mpf(q.numerator) / q.denominator * self.cval[nm]
        got = self.identify_const(x, maxweight=4)
        if got != truth:
            raise FibrateError(
                f"PSLQ POSITIVE CONTROL FAILED: expected {truth}, got {got} "
                f"(dps={self.dps}/{self.dps2}, cap={HEIGHT_CAP}) — refusing "
                "to run")
        return {
            "positive_control": "pass",
            "dual_precision": [self.dps, self.dps2],
            "height_cap": HEIGHT_CAP,
            "residual_gate": f"1e-{self.resid_tol_exp}",
            "zero_gate": f"1e-{self.zero_tol_exp}",
        }

    def negative_control(self):
        """e = 2.718... is not in the basis span: identify_const must
        return None (selfcheck-only; not run per-word)."""
        with mp.workdps(self.dps):
            bad = self.identify_const(mp.e, maxweight=5)
        return bad is None

    # -- numeric evaluation of a rep -----------------------------------------
    def letters_fr(self, w, zeq):
        zr = sp.Rational(zeq.numerator, zeq.denominator)
        out = []
        for a in w:
            v = a.subs(ze, zr)
            out.append(Fr(int(sp.numer(v)), int(sp.denom(v))))
        return out

    def rep_eval(self, r, zeq, gx):
        tot = mp.mpf(0)
        zr = sp.Rational(zeq.numerator, zeq.denominator)
        for (u, K), c in r.items():
            cv = sp.N(c.subs(ze, zr), mp.mp.dps)
            tot += mp.mpf(str(cv)) * gx.G(u) * self.cval_tuple(K)
        return tot

    # -- the fibration --------------------------------------------------------
    def fib(self, w):
        """Fibrate ZIP_reg word w (tuple of exact sympy letters in ze).
        Returns the rep; memoized.  Steps: Goncharov derivative in ze
        (paper:1423-1426; peel-letter relation eq 3.41) -> recursive exact
        integration -> integration constant at ze=1/2 (PSLQ battery) ->
        MANDATORY numeric validation at self.val_points against the
        trusted step-1 ZIP semantics (zip_num)."""
        w = tuple(w)
        if w in self.memo:
            return self.memo[w]
        n = len(w)
        if n == 0:
            self.memo[w] = {((), ONE): sp.Integer(1)}
            self.meta[w] = {'used_pslq': False, 'digits': {}}
            return self.memo[w]
        # I(0; b1..bn; oo) with b = reversed word; total differential in ze:
        # d I = sum_i I(..b_{i-1}, b_{i+1}..) * dlog((b_{i+1}-b_i)/(b_{i-1}-b_i))
        b = [sp.Integer(0)] + [sp.nsimplify(a) for a in reversed(w)] + [None]
        acc = {}
        dep_pslq = False
        for i in range(1, n + 1):
            rho = sp.Integer(0)
            up = b[i + 1]
            if up is not None:
                d = sp.cancel(up - b[i])
                if d != 0:
                    rho += sp.cancel(sp.diff(d, ze) / d)
            dn = sp.cancel(b[i - 1] - b[i])
            if dn != 0:
                rho -= sp.cancel(sp.diff(dn, ze) / dn)
            if rho == 0:
                continue
            # delete b_i  <-> original position n-i
            sub = w[:n - i] + w[n - i + 1:]
            sr = self.fib(sub)
            dep_pslq = dep_pslq or self.meta[sub]['used_pslq']
            for (u, K), c in rep_scale(sr, rho).items():
                rep_add(acc, (u, K), c)
        # integrate exactly in ze
        res = {}
        for (u, K), c in acc.items():
            tmp = {}
            integrate_term(c, u, tmp)
            for (uu, K2), cc in tmp.items():
                rep_add(res, (uu, cmul(K, K2)), cc)
        # integration constant at ze = 1/2 (PSLQ battery — see class doc)
        zw = self.zip12.zip_value(self.letters_fr(w, self.FIT_POINT))
        rv = self.rep_eval(res, self.FIT_POINT, self.g12)
        C = zw - rv
        combo = self.identify_const(C, maxweight=n)
        if combo is None:
            raise FibrateError(
                f"PSLQ failed for word {w}, C={mp.nstr(C, 30)} "
                "(controls: dual-precision/height-cap/residual — no "
                "fabricated constants; REFUSING)")
        # ANY nonzero identified constant (rational 'one' included) is a
        # PSLQ-fitted quantity — stamp it (honest provenance downstream).
        used_pslq = bool(combo)
        for nm, q in combo.items():
            rep_add(res, ((), (nm,) if nm != 'one' else ONE),
                    sp.Rational(q.numerator, q.denominator))
        # MANDATORY per-word validation vs trusted step-1 ZIP semantics
        digits = {}
        for zeq, zipev, gx in self.val_ctx:
            t = zipev.zip_value(self.letters_fr(w, zeq))
            v = self.rep_eval(res, zeq, gx)
            err = abs(t - v)
            sc = max(abs(t), mp.mpf('1e-15'))
            rel = err / sc
            if rel > mp.mpf(10) ** (-self.val_tol_exp):
                raise FibrateError(
                    f"VALIDATION FAIL word {w} at ze={zeq}: "
                    f"zip={mp.nstr(t, 25)} rep={mp.nstr(v, 25)} "
                    f"(rel {mp.nstr(rel, 5)} > 1e-{self.val_tol_exp})")
            d = self.dps if rel == 0 else int(mp.floor(-mp.log10(rel)))
            digits[f"{zeq.numerator}/{zeq.denominator}"] = min(d, self.dps)
        self.memo[w] = res
        self.meta[w] = {'used_pslq': used_pslq or dep_pslq, 'digits': digits}
        return res


# ---------------------------------------------------------------------------
# Canonical serialization (shared by CLI output and fixture generation —
# scripts/fibrate/fixtures/make_regression.py uses these exact functions).
# ---------------------------------------------------------------------------
def serialize_rep(rep):
    """rep -> canonically ordered list of
    {"gword": [letter strs], "consts": [names], "coef": str}."""
    items = []
    for (u, K), c in rep.items():
        items.append({
            "gword": [sp.sstr(a) for a in u],
            "consts": list(K),
            "coef": sp.sstr(sp.cancel(c)),
        })
    items.sort(key=lambda t: (len(t["gword"]), t["gword"], t["consts"],
                              t["coef"]))
    return items


def run_wordlist(words, dps=130, val_points=(Fr(2, 5), Fr(3, 8))):
    """words: list of words (list of letter strings).  Returns the output
    document (dict) — every word fibrated, PSLQ-batteried, 2-point
    validated; raises FibrateError on ANY failure (no partial output)."""
    eng = Engine(dps=dps, val_points=val_points)
    parsed = []
    for wstr in words:
        if not isinstance(wstr, (list, tuple)):
            raise FibrateError(f"wordlist entry not a list: {wstr!r}")
        parsed.append(tuple(parse_letter(x) for x in wstr))
    for i in sorted(range(len(parsed)), key=lambda j: len(parsed[j])):
        try:
            eng.fib(parsed[i])
        except AssertionError as e:
            raise FibrateError(
                f"input refusal on word {words[i]}: {e}") from e
    out_words = []
    any_pslq = False
    for wstr, w in zip(words, parsed):
        meta = eng.meta[w]
        prov = ["goncharov_differential", "validated_2pt"]
        if meta['used_pslq']:
            prov.append("pslq_constants")
            any_pslq = True
        out_words.append({
            "word": list(wstr),
            "weight": len(w),
            "terms": serialize_rep(eng.memo[w]),
            "used_pslq": meta['used_pslq'],
            "provenance": prov,
            "validated_digits": meta['digits'],
        })
    top_prov = ["goncharov_differential", "validated_2pt"]
    if any_pslq:
        top_prov.append("pslq_constants")
    return {
        "tool": TOOL_NAME,
        "version": TOOL_VERSION,
        "source_lineage": SOURCE_LINEAGE,
        "dps": dps,
        "fit_point": "1/2",
        "validation": {
            "points": [f"{p.numerator}/{p.denominator}"
                       for p in eng.val_points],
            "rel_tol": f"1e-{eng.val_tol_exp}",
            "mandatory": True,
            "oracle": "zip_num.ZipEval (trusted step-1 semantics)",
        },
        "pslq_controls": eng.pslq_report,
        "provenance": top_prov,
        "words": out_words,
    }


# ---------------------------------------------------------------------------
# Selfcheck: control battery + invariant regression checks
# ---------------------------------------------------------------------------
def selfcheck(dps=130):
    """Runs (a) PSLQ positive control (Engine init), (b) PSLQ negative
    control, (c) adaptive-ZIP-series cross-dps check (fixed-length
    dps-ceiling guard), (d) GEval/GEvalX cross-check.  Returns a
    report dict; raises FibrateError on any failure."""
    eng = Engine(dps=dps)            # (a) runs positive control
    report = {"pslq_positive_control": "pass"}

    if not eng.negative_control():   # (b)
        raise FibrateError("PSLQ NEGATIVE CONTROL FAILED: e was 'identified' "
                           "in the constant basis")
    report["pslq_negative_control"] = "pass"

    # (c) adaptive series length: letter -9/26 maps to a Hoelder letter of
    # magnitude 9/17 (close to 1/2), where a fixed series length silently
    # caps accuracy independent of dps (the footgun this cap closes).  Independent
    # ground truth: ZIP_reg([-9/26]) = log(26/9) exactly (single-letter
    # word; the reg-stripped 1-letter piece vanishes).
    word = [Fr(-9, 26)]
    with mp.workdps(130):
        hi = ZipEval().zip_value(word)
        truth = mp.log(mp.mpf(26) / 9)
        rel = abs((hi - truth) / truth)
        d = int(mp.floor(-mp.log10(rel))) if rel > 0 else 130
    if d < 120:
        raise FibrateError(f"ADAPTIVE-SERIES CHECK FAILED: ZIP_reg([-9/26]) "
                           f"vs log(26/9) agrees to only {d} digits at "
                           "dps=130 (need >=120) — fixed-length ceiling "
                           "footgun has regressed")
    report["adaptive_series_vs_log_truth_digits"] = d

    # (d) GEval ({0,1,2} strings) vs GEvalX (exact sympy letters) on the
    # same word at x=2/5
    with mp.workdps(60):
        ga = GEval(Fr(2, 5)).G(('0', '1'))
        gb = GEvalX(Fr(2, 5)).G((sp.Integer(0), sp.Integer(1)))
        relg = abs((ga - gb) / gb)
        dg = int(mp.floor(-mp.log10(relg))) if relg > 0 else 60
    if dg < 50:
        raise FibrateError(f"GEval/GEvalX CROSS-CHECK FAILED: {dg} digits")
    report["geval_crosscheck_digits"] = dg

    report["selfcheck"] = "PASS"
    return report


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="fibrate.py",
        description="Sound fibration engine: ZIP_reg words with "
                    "ze-dependent letters -> G(...; ze) x constants.  "
                    "Per-word 2-point validation against step-1 ZIP "
                    "semantics is MANDATORY (no off switch).")
    ap.add_argument("wordlist", nargs="?",
                    help="JSON file: list of words (lists of letter "
                         "strings in ze)")
    ap.add_argument("--out", default=None,
                    help="output JSON path (default: stdout)")
    ap.add_argument("--dps", type=int, default=130,
                    help="working precision (decimal digits, >=60; "
                         "default 130)")
    ap.add_argument("--validation-points", nargs="+",
                    default=["2/5", "3/8"], metavar="P",
                    help=">=2 distinct rationals in (0,1), != 1/2 "
                         "(default: 2/5 3/8)")
    ap.add_argument("--selfcheck", action="store_true",
                    help="run the control battery + invariant checks "
                         "and exit")
    args = ap.parse_args(argv)

    try:
        if args.selfcheck:
            report = selfcheck(dps=args.dps)
            print(json.dumps(report, indent=1))
            return 0
        if not args.wordlist:
            ap.error("wordlist required (or --selfcheck)")
        with open(args.wordlist) as f:
            words = json.load(f)
        if not isinstance(words, list):
            raise FibrateError("wordlist JSON must be a list of words")
        pts = [Fr(p) for p in args.validation_points]
        doc = run_wordlist(words, dps=args.dps, val_points=pts)
        text = json.dumps(doc, indent=1)
        if args.out:
            with open(args.out, "w") as f:
                f.write(text + "\n")
            print(f"fibrate: {len(doc['words'])} words -> {args.out}")
        else:
            print(text)
        return 0
    except (FibrateError, ValueError, ZeroDivisionError, OSError,
            json.JSONDecodeError) as e:
        print(f"fibrate: FATAL: {e}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())

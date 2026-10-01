#!/usr/bin/env python3
"""
famhar.py — generic FAMILY EVALUATOR HARNESS (thin layer over wayfinder).

A family is a CONFIG STAMP: {connection A, anchor, alphabet/branch loci,
masters<->targets map}. The harness supplies, for any family config:

    eval(point_or_path, dps) -> value(s) + partial derivatives (A.J) +
                                sheet record (per-letter winding / continued logs)
    loop(closed_path, dps)   -> monodromy action on the transported vector

Paths are certified wayfinder marches: each straight chart segment is
pulled back to an inline DESystem in t on [0,1] (exact A_series from the
partial-fraction form of the dlog pullback; declared singular_points; the
complex-pole clearance law is wayfinder's own step control, PLUS an
explicit preflight clearance refusal here). eps is threaded through
(fixed-eps marches; the ladder calibration family is eps-independent with
eps0 = 0). Frobenius landings, if ever requested, inherit wayfinder's
REAL-eps0-only law (complex eps raises there by design — respected, never
worked around).

SHEET TRACKING (the homotopy data): the path spec itself is the homotopy
class; the harness derives and returns the abelianized sheet record — for
every letter, the CONTINUOUS winding (accumulated d arg) along the path and
the continued log(letter) at the endpoint. The tracker
unwraps principal args to the accumulated angle, never re-principal-izing
through a crossing. Per straight segment the
letter image is itself a straight segment, so |delta arg| < pi and the
principal-ratio increment is exact; the tracker still bisects to |.|<pi/2
per sub-step (numerical safety near the +-pi edge) and REFUSES any segment
passing within `clearance_min` of a branch locus (the seg-gating bug class:
nothing is ever evaluated across an uncertified region).

CONNECTION FORMAT: dlog form. letters are polynomials in the chart vars
with exact Fraction coefficients — AFFINE letters (letter = c0 + sum_i
c_i * x_i, the {"const", <var>} schema) evaluate their dlog pullback in
closed form per segment; POLYNOMIAL letters (the {"terms"} canonical-terms
schema, total degree >= 2) are composed EXACTLY along each straight chart
segment (Gaussian-rational arithmetic) and their t-polynomial is rooted
with mp.polyroots — dlog w(t) = sum_roots dt/(t - r) — so the same pole
terms, clearance certification (t-root distance to [0,1]) and per-factor
continuous winding apply. Root-finding failure REFUSES the segment
(clearance cannot be certified without the roots), and polynomial-letter
marches REFUSE floating waypoints (exact composition needs exact
rationals). The connection is
    dJ = [ sum_letters M_letter dlog(letter) ] J,
M_letter sparse with exact Fraction entries. Stored single-var eps-dependent
systems route through wayfinder's own loaders instead (load_monomial_json /
load_amatrix_json / load_kira_targets) — this layer does not duplicate them.
Elliptic families: the connection slot is the same; the anchor/reference
plugin carries the period machinery (Eichler-integral pattern) — see
README extension points.

No global mp.dps mutation anywhere (workdps only). All config numbers parse
as exact Fractions (the ambient-dps input-conversion trap is structurally
excluded).
"""
import json
import sys
import os
from fractions import Fraction

from mpmath import mp, mpf, mpc, pi

# repo layout: wayfinder ships beside this package in tools/
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from wayfinder import transport_fixed_eps  # noqa: E402

_HALF_PI = None


def _frac(s):
    """Parse 'p/q' or 'p' or 'p/q+r/qi' style rational-complex strings into
    (Fraction_re, Fraction_im)."""
    s = str(s).strip().replace(" ", "")
    if "j" in s or "i" in s:
        # form: re+imi / re-imi / imi
        s2 = s.replace("i", "j")
        # split on last +/- not at position 0 and not inside a fraction
        for k in range(len(s2) - 1, 0, -1):
            if s2[k] in "+-" and s2[k - 1] not in "+-/*e":
                re_part, im_part = s2[:k], s2[k:-1]
                if im_part in ("+", "-"):
                    im_part += "1"
                return Fraction(re_part), Fraction(im_part)
        return Fraction(0), Fraction(s2[:-1] if s2[:-1] not in ("", "+", "-") else s2[:-1] + "1")
    return Fraction(s), Fraction(0)


def frac_to_mpc(fr, dps):
    re_, im_ = fr
    with mp.workdps(dps):
        return mpc(mpf(re_.numerator) / re_.denominator,
                   mpf(im_.numerator) / im_.denominator)


class Letter:
    """Letter in the chart vars, exact Fraction coefficients. Two schemas:

    - affine: {"name", "const", "<var>": coeff, ...} -> c0 + sum_i ci * x_i
    - polynomial: {"name", "terms": [[[e_1..e_k], "p/q"], ...]} over the
      chart vars, in the same canonical-terms format as the rational
      connection block (term = [exponent tuple, coefficient string]).

    A "terms" letter of total degree <= 1 folds into the affine fast path;
    genuinely polynomial letters (degree >= 2) switch the segment machinery
    to per-segment exact composition + root-finding (Family._poly_seg_analysis)."""

    def __init__(self, name, spec, chart_vars):
        self.name = name
        self.vars = list(chart_vars)
        self.terms = None
        if "terms" in spec:
            extra = set(spec) - {"name", "terms"}
            if extra:
                raise ValueError(f"letter '{name}': 'terms' schema does not mix "
                                 f"with affine keys {sorted(extra)}")
            terms = []
            for t in spec["terms"]:
                exps = tuple(int(e) for e in t[0])
                if len(exps) != len(chart_vars):
                    raise ValueError(f"letter '{name}': exponent tuple length "
                                     f"{len(exps)} != {len(chart_vars)} chart vars")
                if any(e < 0 for e in exps):
                    raise ValueError(f"letter '{name}': negative exponent — refused")
                c = Fraction(t[1])
                if c:
                    terms.append((exps, c))
            if not terms:
                raise ValueError(f"letter '{name}': identically zero — refused")
            if max(sum(e) for e, _ in terms) >= 2:
                self.terms = tuple(terms)
            else:
                # total degree <= 1: fold into the affine fast path
                self.c0 = Fraction(0)
                self.coef = {v: Fraction(0) for v in chart_vars}
                for exps, c in terms:
                    if sum(exps) == 0:
                        self.c0 += c
                    else:
                        self.coef[chart_vars[exps.index(1)]] += c
        else:
            self.c0 = Fraction(spec.get("const", "0"))
            self.coef = {v: Fraction(spec.get(v, "0")) for v in chart_vars}

    @property
    def affine(self):
        return self.terms is None

    def value(self, point, dps):
        """point: dict var -> mpc."""
        with mp.workdps(dps):
            if self.terms is None:
                v = mpc(mpf(self.c0.numerator) / self.c0.denominator)
                for var, c in self.coef.items():
                    if c:
                        v += (mpf(c.numerator) / c.denominator) * point[var]
                return v
            v = mpc(0)
            for exps, c in self.terms:
                t = mpc(mpf(c.numerator) / c.denominator)
                for var, e in zip(self.vars, exps):
                    if e:
                        t *= point[var] ** e
                v += t
            return v

    def dvalue(self, point, dps, var):
        """d(letter)/d var at point (mpc); exact-coefficient derivative."""
        with mp.workdps(dps):
            if self.terms is None:
                c = self.coef[var]
                return mpc(mpf(c.numerator) / c.denominator)
            v = mpc(0)
            for exps, c in self.terms:
                k = self.vars.index(var)
                if not exps[k]:
                    continue
                t = exps[k] * mpc(mpf(c.numerator) / c.denominator)
                for var2, e2 in zip(self.vars, exps):
                    ee = e2 - 1 if var2 == var else e2
                    if ee:
                        t *= point[var2] ** ee
                v += t
            return v


class Family:
    def __init__(self, cfg, cfg_path=None):
        self.cfg = cfg
        self.cfg_path = cfg_path
        self.name = cfg["name"]
        self.chart_vars = cfg["chart_vars"]
        self.n = cfg["n_masters"]
        self.masters = cfg["masters"]
        self.letters = [Letter(l["name"], l, self.chart_vars) for l in cfg["letters"]]
        e = _frac(cfg.get("eps0", "0"))
        self.eps0_frac = e
        # sparse connection: letter name -> {(i,j): Fraction}
        self.conn = {}
        for lname, entries in cfg["connection"].items():
            m = {}
            for key, val in entries.items():
                i, j = (int(t) for t in key.split(","))
                m[(i, j)] = Fraction(val)
            self.conn[lname] = m
        for lname in self.conn:
            if lname not in {l.name for l in self.letters}:
                raise ValueError(f"connection letter {lname} not declared")
        self.targets = cfg.get("targets", {})
        self.clearance_min = Fraction(cfg.get("clearance_min", "1/50"))

    @classmethod
    def from_json(cls, path):
        with open(path) as f:
            return cls(json.load(f), cfg_path=path)

    # ---------------- points ------------------------------------------------
    def parse_point(self, pt, dps):
        """pt: dict var -> rational string / Fraction / (Fraction,Fraction) /
        mpc/number. Returns dict var->mpc (exact-rational parse when given
        exact forms)."""
        out = {}
        with mp.workdps(dps):
            for v in self.chart_vars:
                x = pt[v]
                if isinstance(x, str):
                    out[v] = frac_to_mpc(_frac(x), dps)
                elif isinstance(x, Fraction):
                    out[v] = frac_to_mpc((x, Fraction(0)), dps)
                elif isinstance(x, tuple):
                    out[v] = frac_to_mpc(x, dps)
                else:
                    out[v] = mpc(x)
        return out

    # ---------------- pullback DESystem ------------------------------------
    def segment_desys(self, p0, p1, dps, poly=None):
        """Inline DESystem for the straight chart segment p0 -> p1, t in [0,1].
        poly: _poly_seg_analysis output for this segment (required iff the
        connection uses polynomial letters): each root r of a polynomial
        letter's t-polynomial contributes a unit pole term dt/(t - r), since
        dlog w(t) = sum_roots dt/(t - r) (the leading coefficient is constant
        in t and drops out)."""
        fam = self

        class SegDE:
            n = fam.n
            var = "t"

            def __init__(self):
                wp = dps + 15
                with mp.workdps(wp):
                    self._terms = []   # (conn_matrix, pole t_p) ; dlog letter = 1/(t - t_p)
                    self.singular_points = []
                    for L in fam.letters:
                        if L.name not in fam.conn:
                            continue
                        if not L.affine:
                            if poly is None or L.name not in poly:
                                raise ValueError(
                                    f"segment_desys: polynomial letter '{L.name}' "
                                    "needs the per-segment root analysis "
                                    "(march supplies it; direct callers use "
                                    "Family._poly_seg_analysis)")
                            for r in poly[L.name]["roots"]:
                                t_p = mpc(r)
                                self._terms.append((fam.conn[L.name], t_p))
                                self.singular_points.append(t_p)
                            continue
                        w0 = L.value(p0, wp)
                        w1 = L.value(p1, wp)
                        slope = w1 - w0
                        if slope == 0:
                            continue   # letter constant along segment: no dlog
                        t_p = -w0 / slope
                        self._terms.append((fam.conn[L.name], t_p))
                        self.singular_points.append(t_p)
                self.meta = {"source": f"famhar segment {fam.name}", "format": "inline-dlog",
                             "singular_points": self.singular_points}

            def A(self, x, eps, dps_):
                with mp.workdps(dps_):
                    xv = mpc(x)
                    Amat = [[mpc(0)] * fam.n for _ in range(fam.n)]
                    for conn, t_p in self._terms:
                        g = 1 / (xv - t_p)
                        for (i, j), c in conn.items():
                            Amat[i][j] += (mpf(c.numerator) / c.denominator) * g
                    return Amat

            def A_series(self, z0, eps, dps_, M):
                with mp.workdps(dps_ + 10):
                    z0v = mpc(z0)
                    zero = mpc(0)
                    out = [[[zero] * fam.n for _ in range(fam.n)] for _ in range(M + 1)]
                    for conn, t_p in self._terms:
                        d = t_p - z0v
                        if d == 0:
                            raise ValueError("A_series requested AT a singular point")
                        inv = 1 / d
                        # 1/(t-t_p) = -sum_k (t-z0)^k / d^{k+1}
                        pk = inv
                        for k in range(M + 1):
                            g = -pk
                            for (i, j), c in conn.items():
                                cur = out[k][i][j]
                                out[k][i][j] = cur + (mpf(c.numerator) / c.denominator) * g
                            pk *= inv
                    return out

        return SegDE()

    # ---------------- clearance + winding ----------------------------------
    def _poly_seg_analysis(self, p0x, p1x, dps):
        """Per-segment data for POLYNOMIAL connection letters (degree >= 2):
        the letter is composed EXACTLY along the straight chart segment
        (Gaussian-rational arithmetic) and its t-polynomial rooted with
        mp.polyroots — the same root-finding + clearance-certification
        machinery as the rational-block segments. Returns
        {letter_name: {"poly_t": gq coeffs, "roots": [mpc]}} or None when
        the connection has no polynomial letters.

        REFUSALS (fail-closed): floating waypoints (exact composition
        impossible — pass rational strings/ints/Fractions); a letter
        identically zero along the segment (path ON a branch locus);
        root-finding failure (clearance cannot be certified without the
        roots)."""
        polyL = [L for L in self.letters if L.name in self.conn and not L.affine]
        if not polyL:
            return None
        if p0x is None or p1x is None:
            raise ValueError(
                "polynomial-letter march REFUSED: waypoints must be EXACT "
                "rationals (string 'p/q(+r/si)', int, Fraction or "
                "(Fraction,Fraction) pairs) so the letter t-polynomial can "
                "be composed exactly before root-finding")
        out = {}
        wp = dps + 25
        for L in polyL:
            pt = _letter_poly_t(L, p0x, p1x, self.chart_vars)
            if len(pt) == 1:
                if _gq_is0(pt[0]):
                    raise ValueError(
                        f"letter '{L.name}' identically ZERO along segment — "
                        "path lies ON a branch locus, REFUSED (seg-gating law)")
                out[L.name] = {"poly_t": pt, "roots": []}
                continue
            with mp.workdps(wp + 10):
                coeffs = [_gq_to_mpc(c, wp + 10) for c in reversed(pt)]
                try:
                    rts = mp.polyroots(coeffs, maxsteps=200, extraprec=60)
                except Exception as e:
                    raise ValueError(
                        f"letter '{L.name}': segment root-finding failed: "
                        f"{type(e).__name__}: {e} — REFUSED (cannot certify "
                        "clearance without the roots)")
                out[L.name] = {"poly_t": pt, "roots": list(rts)}
        return out

    def _normalize_waypoints(self, waypoints, dps):
        """Accept parsed waypoints (var -> mpc/number) or raw-exact waypoints
        (var -> rational string / int / Fraction / (Fraction,Fraction)).
        Returns (mpc_points, exact_points); exact_points is None unless EVERY
        waypoint parsed exactly (polynomial-letter marches require that)."""
        mpts, xpts = [], []
        for p in waypoints:
            if all(isinstance(p[v], (str, int, Fraction, tuple))
                   and not isinstance(p[v], bool) for v in self.chart_vars):
                px = self.parse_point_exact(p)
                xpts.append(px)
                mpts.append({v: frac_to_mpc(px[v], dps) for v in self.chart_vars})
            else:
                xpts.append(None)
                mpts.append(dict(p))
        return mpts, (xpts if all(x is not None for x in xpts) else None)

    def _seg_clearance(self, p0, p1, dps, poly=None):
        """min over letters (that appear in the connection) of the segment
        clearance. AFFINE letters: distance from the letter-image straight
        segment [w0,w1] to 0, normalized by the letter scale max(|w0|,|w1|).
        POLYNOMIAL letters (poly = _poly_seg_analysis output): distance of
        the nearest t-root to the real segment [0,1], in t-units (the
        rational-block clearance convention — same clearance_min bar).
        Returns (min_clear, letter_name)."""
        with mp.workdps(dps):
            worst = None
            for L in self.letters:
                if L.name not in self.conn or not L.affine:
                    continue
                w0 = L.value(p0, dps)
                w1 = L.value(p1, dps)
                scale = max(abs(w0), abs(w1))
                if scale == 0:
                    return mpf(0), L.name
                d = w1 - w0
                if abs(d) == 0:
                    c = abs(w0) / scale
                else:
                    # distance from segment w0 + s*d, s in [0,1], to origin
                    s = -(w0.real * d.real + w0.imag * d.imag) / abs(d) ** 2
                    s = min(max(s, mpf(0)), mpf(1))
                    c = abs(w0 + s * d) / scale
                if worst is None or c < worst[0]:
                    worst = (c, L.name)
            for name, dat in (poly or {}).items():
                for r in dat["roots"]:
                    re_, im_ = mp.re(r), mp.im(r)
                    if 0 <= re_ <= 1:
                        c = abs(im_)
                    else:
                        c = min(abs(r), abs(r - 1))
                    if worst is None or c < worst[0]:
                        worst = (mpf(c), name)
            return worst if worst is not None else (mp.inf, "(no gating letters)")

    @staticmethod
    def _delta_arg(w0, w1, depth=0):
        """Continuous d(arg) along the straight segment w0->w1 (avoiding 0).
        Bisect until |increment| <= pi/2; principal Im log of the ratio is
        then exact (never re-principal-ize)."""
        if depth > 60:
            raise ValueError("winding bisection depth cap: segment too close to a branch locus")
        if w0 == 0 or w1 == 0:
            raise ValueError("path endpoint ON branch locus — REFUSED (seg-gating law)")
        r = mp.arg(w1 / w0)
        if abs(r) <= mp.pi / 2:
            return r
        wm = (w0 + w1) / 2
        if wm == 0:
            raise ValueError("path through branch locus")
        return (Family._delta_arg(w0, wm, depth + 1) +
                Family._delta_arg(wm, w1, depth + 1))

    def segment_winding(self, p0, p1, dps, poly=None):
        """dict letter -> continuous delta-arg along the straight chart
        segment. AFFINE letters: the letter image is itself a straight
        segment, tracked directly. POLYNOMIAL letters (poly =
        _poly_seg_analysis output): w(t) = lc * prod_k (t - r_k) with lc
        constant in t, so the winding is the exact per-factor sum of the
        same tracker applied to each linear factor's straight image."""
        out = {}
        with mp.workdps(dps + 10):
            for L in self.letters:
                if L.name not in self.conn or not L.affine:
                    continue
                w0 = L.value(p0, dps + 10)
                w1 = L.value(p1, dps + 10)
                out[L.name] = self._delta_arg(w0, w1)
            for name, dat in (poly or {}).items():
                w = mpf(0)
                for r in dat["roots"]:
                    w += self._delta_arg(mpc(0) - r, mpc(1) - r)
                out[name] = w
        return out

    # ---------------- path march -------------------------------------------
    def march(self, waypoints, y0, dps, eps0=None, mtay=None):
        """March J along piecewise-straight chart path. waypoints: list of
        point dicts INCLUDING start — parsed (var -> mpc) or raw-exact
        (var -> rational string / int / Fraction / (Fraction,Fraction)).
        Families whose connection uses POLYNOMIAL letters require the
        raw-exact form (the letter t-polynomial is composed exactly per
        segment before root-finding; floats REFUSE). y0 = J at waypoints[0].

        Returns (J_end, record) with record = {
          'windings': letter -> total continuous delta-arg,
          'log_letters': letter -> continued log(letter) at endpoint,
          'diags': per-segment wayfinder diag dicts,
          'clearances': per-segment (min_clear, letter),
          'walls_s': per-segment wall seconds }."""
        import time
        wps, wpsx = self._normalize_waypoints(waypoints, dps + 15)
        eps_v = eps0 if eps0 is not None else frac_to_mpc(self.eps0_frac, dps)
        wind = {L.name: mpf(0) for L in self.letters if L.name in self.conn}
        diags, clears, walls = [], [], []
        y = list(y0)
        cmin_bar = mpf(self.clearance_min.numerator) / self.clearance_min.denominator
        with mp.workdps(dps + 10):
            for k in range(len(wps) - 1):
                p0, p1 = wps[k], wps[k + 1]
                pa = self._poly_seg_analysis(
                    wpsx[k] if wpsx else None,
                    wpsx[k + 1] if wpsx else None, dps)
                cl = self._seg_clearance(p0, p1, dps + 10, poly=pa)
                clears.append((float(cl[0]), cl[1]))
                if cl[0] < cmin_bar:
                    raise ValueError(
                        f"segment {k}: letter '{cl[1]}' clearance {float(cl[0]):.3g} < "
                        f"clearance_min {float(cmin_bar):.3g} — REFUSED (seg-gating law). "
                        f"Fix: insert a complex detour waypoint.")
                sw = self.segment_winding(p0, p1, dps + 10, poly=pa)
                for nm, d in sw.items():
                    wind[nm] += d
                de = self.segment_desys(p0, p1, dps, poly=pa)
                t0 = time.time()
                y, dg = transport_fixed_eps(de, eps_v, 0, 1, y, dps,
                                            mtay=mtay, return_diag=True)
                walls.append(round(time.time() - t0, 3))
                diags.append(dg)
            # continued letter logs at endpoint
            pe = wps[-1]
            p_start = wps[0]
            log_letters = {}
            for L in self.letters:
                if L.name not in self.conn:
                    continue
                w_end = L.value(pe, dps + 10)
                w_st = L.value(p_start, dps + 10)
                log_letters[L.name] = mp.log(abs(w_end)) + mpc(0, 1) * (mp.arg(w_st) + wind[L.name])
        return y, {"windings": wind, "log_letters": log_letters, "diags": diags,
                   "clearances": clears, "walls_s": walls}

    # ---------------- targets + derivatives ---------------------------------
    def _pref(self, tgt, point, dps, deriv_var=None):
        expr = tgt["prefactor"] if deriv_var is None else tgt["dprefactor"][deriv_var]
        ns = dict(point)
        with mp.workdps(dps):
            return eval(expr, {"__builtins__": {}}, ns)

    def target_value(self, tname, J, point, dps):
        tgt = self.targets[tname]
        with mp.workdps(dps):
            s = mpc(0)
            for key, val in tgt["coeffs"].items():
                c = Fraction(val)
                s += (mpf(c.numerator) / c.denominator) * J[int(key)]
            return self._pref(tgt, point, dps) * s

    def A_at(self, point, dps):
        """dict var -> dense A_var(point) (n x n), from the dlog connection:
        A_var = sum_letters M_letter * (d letter/d var) / letter(point).
        Rational => single-valued => sheet-safe by construction."""
        out = {}
        with mp.workdps(dps):
            for var in self.chart_vars:
                Amat = [[mpc(0)] * self.n for _ in range(self.n)]
                for L in self.letters:
                    if L.name not in self.conn:
                        continue
                    if L.affine:
                        cv = L.coef[var]
                        if not cv:
                            continue
                        g = (mpf(cv.numerator) / cv.denominator) / L.value(point, dps)
                    else:
                        dv = L.dvalue(point, dps, var)
                        if dv == 0:
                            continue
                        g = dv / L.value(point, dps)
                    for (i, j), c in self.conn[L.name].items():
                        Amat[i][j] += (mpf(c.numerator) / c.denominator) * g
                out[var] = Amat
        return out

    def target_derivs_chart(self, tname, J, point, dps):
        """dict var -> d(target)/d var via A.J + prefactor product rule."""
        tgt = self.targets[tname]
        Av = self.A_at(point, dps)
        with mp.workdps(dps):
            combo = mpc(0)
            cvec = {int(k): Fraction(v) for k, v in tgt["coeffs"].items()}
            for i, c in cvec.items():
                combo += (mpf(c.numerator) / c.denominator) * J[i]
            pref = self._pref(tgt, point, dps)
            out = {}
            for var in self.chart_vars:
                dJ = [mpc(0)] * self.n
                for i in range(self.n):
                    row = Av[var][i]
                    acc = mpc(0)
                    for j in range(self.n):
                        if row[j] != 0 and J[j] != 0:
                            acc += row[j] * J[j]
                    dJ[i] = acc
                dcombo = mpc(0)
                for i, c in cvec.items():
                    dcombo += (mpf(c.numerator) / c.denominator) * dJ[i]
                dpref = self._pref(tgt, point, dps, deriv_var=var)
                out[var] = dpref * combo + pref * dcombo
            return out

    def target_derivs_kin(self, tname, J, point, dps):
        """dict kin_var -> d(target)/d kin_var, via the config's kinematic
        Jacobian (2 chart vars x 2 kin vars only, v1)."""
        kin = self.cfg["kinematics"]
        kvars = kin["vars"]
        cvars = self.chart_vars
        if len(kvars) != 2 or len(cvars) != 2:
            raise NotImplementedError("v1: 2x2 kinematic maps only")
        dch = self.target_derivs_chart(tname, J, point, dps)
        ns = dict(point)
        with mp.workdps(dps):
            Jk = [[eval(kin["jacobian"][f"d{kv}_d{cv}"], {"__builtins__": {}}, ns)
                   for kv in kvars] for cv in cvars]   # rows: chart var, cols: kin var
            det = Jk[0][0] * Jk[1][1] - Jk[0][1] * Jk[1][0]
            r0, r1 = dch[cvars[0]], dch[cvars[1]]
            d_k0 = (r0 * Jk[1][1] - r1 * Jk[0][1]) / det
            d_k1 = (r1 * Jk[0][0] - r0 * Jk[1][0]) / det
            return {kvars[0]: d_k0, kvars[1]: d_k1}

    # ---------------- high-level eval ---------------------------------------
    def evaluate(self, path_spec, y0, dps, targets=None, derivs=True):
        """path_spec: list of raw point dicts (strings ok); y0 = J at start.
        Returns dict with values, derivatives (chart + kinematic), sheet
        record and transport diagnostics. Exactness of raw rational
        waypoints is preserved into the march (polynomial-letter families
        need it for the per-segment exact composition)."""
        J, rec = self.march(path_spec, y0, dps)
        pe = self.parse_point(path_spec[-1], dps + 15)
        out = {"J": J, "sheet": rec, "point": pe, "targets": {}}
        for tname in (targets or self.targets.keys()):
            entry = {"value": self.target_value(tname, J, pe, dps)}
            if derivs:
                entry["d_chart"] = self.target_derivs_chart(tname, J, pe, dps)
                if "kinematics" in self.cfg:
                    entry["d_kin"] = self.target_derivs_kin(tname, J, pe, dps)
            out["targets"][tname] = entry
        return out

    # ==================== RATIONAL-BLOCK MARCH =============================
    # Additive extension to the dlog engine above. It marches the EXACT
    # rational connection block (config slot
    # `connection_rational`: sha-pinned per-op A files, entries num/den
    # canonical_terms over the file's vars) through the UNCHANGED wayfinder
    # DESystem duck-typed contract (.n/.var/.A/.A_series/.singular_points).
    # Laws imported unchanged: fail-closed pin verification, the seg-gating
    # refusal law (R1 class), the case-3 continuous-arg tracker (R2 class,
    # via the SAME _delta_arg applied per exact linear factor t - root), and
    # the SQUARE-SYSTEM refusal: a block whose entry columns are not covered
    # by its DE rows has NO march defined — march_rational REFUSES with
    # SYSTEM_NOT_CLOSED (never integrates a partial block silently).

    def load_connection_rational(self, spec=None, verify_sha=True):
        """Load (and cache) the config's `connection_rational` block.

        spec: override dict with the same shape as the config slot
        ({"files": {"A_<op>": {"path":..., "sha256":...}, ...}}). Paths are
        resolved relative to the config file's directory. Every file's
        sha256 is verified against its pin (fail-closed) unless
        verify_sha=False (testing only; the loader records the mode).
        Returns a RationalConnection."""
        if spec is None:
            spec = self.cfg.get("connection_rational")
            if getattr(self, "_ratconn", None) is not None:
                return self._ratconn
        if not spec or "files" not in spec:
            raise ValueError("load_connection_rational REFUSED: no "
                             "connection_rational.files block in config")
        # path convention (i214 config + fullstamp-battery precedent): file
        # paths are relative to the HARNESS dir = parent of the config's
        # families/ dir; for configs not living in a families/ dir, the
        # config's own dir.
        if self.cfg_path:
            cdir = os.path.dirname(os.path.abspath(self.cfg_path))
            base = os.path.dirname(cdir) if os.path.basename(cdir) == "families" else cdir
        else:
            base = os.getcwd()
        rc = RationalConnection.from_files(spec["files"], base, verify_sha=verify_sha)
        if spec.get("basis_n") is not None and rc.n != spec["basis_n"]:
            raise ValueError(f"connection_rational basis_n pin {spec['basis_n']} "
                             f"!= file basis length {rc.n}")
        if spec is self.cfg.get("connection_rational"):
            self._ratconn = rc
        return rc

    def parse_point_exact(self, pt):
        """dict var -> (Fraction re, Fraction im). REFUSES non-exact input
        (floats would silently break the exact segment composition)."""
        out = {}
        for v in self.chart_vars:
            x = pt[v]
            if isinstance(x, str):
                out[v] = _frac(x)
            elif isinstance(x, Fraction):
                out[v] = (x, Fraction(0))
            elif isinstance(x, int):
                out[v] = (Fraction(x), Fraction(0))
            elif isinstance(x, tuple) and len(x) == 2 and all(isinstance(t, Fraction) for t in x):
                out[v] = x
            else:
                raise ValueError(
                    f"parse_point_exact REFUSED: waypoint var '{v}' has type "
                    f"{type(x).__name__}; rational-block marches require EXACT "
                    "rational(-complex) waypoints (string 'p/q(+r/si)', int, "
                    "Fraction, or (Fraction,Fraction))")
        return out

    def march_rational(self, path_spec, y0, dps, eps0=None, subs=None,
                       mtay=None, ratconn=None):
        """March J along a piecewise-straight chart path using the EXACT
        rational connection block: dJ = [sum_op A_op(vars) d op] J.

        path_spec: list of raw point dicts (exact rationals; floats refuse),
        INCLUDING the start. y0 = J at path_spec[0], length = basis n.
        subs: {var: exact rational} for every file var that is NOT a chart
        var (e.g. {"d": "28/5"} for a fixed-eps march with d = 4-2*eps).
        The caller owns the eps->var law; famhar only checks coverage.

        Fail-closed preconditions (REFUSED, never worked around):
          - closure: every entry column must be covered by the block's own
            DE rows (square system). A 474x706-class open block refuses with
            SYSTEM_NOT_CLOSED — marching it would fabricate a solution.
          - vars coverage: file vars == chart_vars + subs keys, exactly.
          - t-plane clearance: any denominator root within clearance_min of
            the segment [0,1] refuses (seg-gating law; fix = detour).

        Returns (J_end, record); record mirrors march(): 'windings' and
        'log_dens' are PER DISTINCT DENOMINATOR (canonical term-key), the
        winding computed exactly per linear factor (t - root) with the same
        case-3 _delta_arg tracker, plus 'clearances', 'diags', 'walls_s',
        'closure' (the measured census) and 'provenance' (file sha16s)."""
        import time
        rc = ratconn if ratconn is not None else self.load_connection_rational()
        cen = rc.closure_census()
        if not cen["closed"]:
            raise ValueError(
                "march_rational REFUSED: SYSTEM_NOT_CLOSED("
                f"{cen['n_rows_covered']}x{rc.n}): {cen['n_cols_outside_rows']} "
                "entry columns land outside the row-covered set (escape/deferred "
                "DE rows pending). No march is defined on a non-square block; "
                "integrating it would silently fabricate a solution.")
        if len(y0) != rc.n:
            raise ValueError(f"march_rational REFUSED: len(y0)={len(y0)} != basis n={rc.n}")
        subs = dict(subs or {})
        subs_f = {}
        for k, v in subs.items():
            subs_f[k] = _frac(v) if isinstance(v, str) else \
                ((v, Fraction(0)) if isinstance(v, (int, Fraction)) else v)
            if not (isinstance(subs_f[k], tuple) and len(subs_f[k]) == 2):
                raise ValueError(f"march_rational REFUSED: subs['{k}'] not exact rational")
        need = set(rc.vars) - set(self.chart_vars)
        if set(subs_f) != need:
            raise ValueError(
                f"march_rational REFUSED: file vars {rc.vars} need subs for "
                f"{sorted(need)}, got {sorted(subs_f)} — every non-chart var "
                "must be pinned exactly (no silent defaults)")
        wps = [self.parse_point_exact(p) for p in path_spec]
        eps_v = eps0 if eps0 is not None else frac_to_mpc(self.eps0_frac, dps)
        cmin = Fraction(self.cfg.get("clearance_min", "1/50")) if not hasattr(self, "clearance_min") \
            else self.clearance_min
        wind, diags, clears, walls = {}, [], [], []
        y = list(y0)
        seg = None
        with mp.workdps(dps + 10):
            cmin_f = mpf(cmin.numerator) / cmin.denominator
            start_vals = None
            for k in range(len(wps) - 1):
                seg = rc.segment_desys(wps[k], wps[k + 1], self.chart_vars,
                                       subs_f, dps)
                if start_vals is None:
                    start_vals = seg.den_values_at(mpf(0), dps + 10)
                cl_val, cl_key = seg.clearance()
                clears.append((float(cl_val), cl_key))
                if cl_val < cmin_f:
                    raise ValueError(
                        f"segment {k}: denominator '{cl_key}' t-root clearance "
                        f"{float(cl_val):.3g} < clearance_min {float(cmin_f):.3g} "
                        "— REFUSED (seg-gating law). Fix: insert a complex "
                        "detour waypoint.")
                for key, dw in seg.windings(dps + 10).items():
                    wind[key] = wind.get(key, mpf(0)) + dw
                t0 = time.time()
                y, dg = transport_fixed_eps(seg, eps_v, 0, 1, y, dps,
                                            mtay=mtay, return_diag=True)
                walls.append(round(time.time() - t0, 3))
                diags.append(dg)
            log_dens = {}
            if seg is not None:
                end_vals = seg.den_values_at(mpf(1), dps + 10)
                for key in wind:
                    w_end = end_vals.get(key)
                    w_st = start_vals.get(key)
                    if w_end is None or w_st is None:
                        continue
                    log_dens[key] = (mp.log(abs(w_end))
                                     + mpc(0, 1) * (mp.arg(w_st) + wind[key]))
        rec = {"windings": wind, "log_dens": log_dens, "diags": diags,
               "clearances": clears, "walls_s": walls, "closure": cen,
               "provenance": rc.provenance}
        return y, rec

    def A_rational_at(self, point, dps, subs=None, ratconn=None):
        """Pointwise sparse A_op matrices from the rational block:
        {op: {(i,j): mpc}}. Pointwise evaluation is legal on a non-closed
        block (it IS the stored A); only marching requires closure."""
        rc = ratconn if ratconn is not None else self.load_connection_rational()
        subs = dict(subs or {})
        pt = {}
        for v, x in list(point.items()) + list(subs.items()):
            if isinstance(x, str):
                pt[v] = frac_to_mpc(_frac(x), dps)
            elif isinstance(x, tuple):
                pt[v] = frac_to_mpc(x, dps)
            else:
                with mp.workdps(dps):
                    pt[v] = mpc(x)
        missing = set(rc.vars) - set(pt)
        if missing:
            raise ValueError(f"A_rational_at REFUSED: vars {sorted(missing)} unpinned")
        return rc.eval_at(pt, dps)


# ---------------------------------------------------------------------------
# v2 rational-block support classes (F3 leg). Exact Gaussian-rational
# polynomial composition; the per-segment DESystem satisfies the unchanged
# wayfinder contract with an EXACT-coefficient A_series fast path
# (numerator/denominator Taylor shift + inverse-series division).
# ---------------------------------------------------------------------------
def _gq_add(x, y):
    return (x[0] + y[0], x[1] + y[1])


def _gq_mul(x, y):
    return (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0])


def _gq_is0(x):
    return x[0] == 0 and x[1] == 0


_GQ0 = (Fraction(0), Fraction(0))
_GQ1 = (Fraction(1), Fraction(0))


def _poly_mul(p, q):
    out = [_GQ0] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        if _gq_is0(a):
            continue
        for j, b in enumerate(q):
            if _gq_is0(b):
                continue
            out[i + j] = _gq_add(out[i + j], _gq_mul(a, b))
    return out


def _poly_strip(p):
    while len(p) > 1 and _gq_is0(p[-1]):
        p = p[:-1]
    return p


def _gq_to_mpc(x, dps):
    with mp.workdps(dps):
        return mpc(mpf(x[0].numerator) / x[0].denominator,
                   mpf(x[1].numerator) / x[1].denominator)


def _sha256_file(path):
    import hashlib
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def _den_key(den_terms):
    """Canonical key for a denominator: sorted terms, json-compact."""
    return json.dumps(sorted([[list(e), str(c)] for e, c in
                              ((tuple(t[0]), t[1]) for t in den_terms)]),
                      separators=(",", ":"))


def _letter_poly_t(letter, p0x, p1x, chart_vars):
    """Exact Gaussian-rational t-polynomial of a polynomial letter along the
    straight chart segment p0x -> p1x (t in [0,1]); p0x/p1x: {var: (Fr,Fr)}.
    Returns ascending gq coefficients, trailing exact zeros stripped (degree
    drops from special segment directions are detected exactly — no
    numerical thresholds)."""
    aff = {}
    for v in chart_vars:
        a0 = p0x[v]
        aff[v] = (a0, _gq_add(p1x[v], (-a0[0], -a0[1])))
    pow_cache = {}

    def var_pow(v, e):
        key = (v, e)
        if key in pow_cache:
            return pow_cache[key]
        base = _poly_strip([aff[v][0], aff[v][1]])
        r = [_GQ1]
        for _ in range(e):
            r = _poly_mul(r, base)
        pow_cache[key] = r
        return r

    acc = [_GQ0]
    for exps, c in letter.terms:
        p = [(c, Fraction(0))]
        for v, e in zip(chart_vars, exps):
            if e:
                p = _poly_mul(p, var_pow(v, e))
        if len(p) > len(acc):
            acc = acc + [_GQ0] * (len(p) - len(acc))
        for k, a in enumerate(p):
            acc[k] = _gq_add(acc[k], a)
    return _poly_strip(acc)


class RationalConnection:
    """Exact rational connection block: per op var, sparse entries
    (i,j) -> num/den with canonical_terms over `vars`.

    File format (the i214 A-block format): {"basis": [...], "vars": [...],
    "entries": {"i,j": [[num_terms],[den_terms]]}} with term =
    [[e_1..e_nvars], coeff_string]."""

    def __init__(self, vars_, n_basis, ops, den_polys, provenance, zero_rows=()):
        self.vars = list(vars_)
        self.n = n_basis
        self.ops = ops            # op -> {(i,j): (num_terms, den_key)}
        self.den_polys = den_polys  # den_key -> term list
        self.provenance = provenance
        self.zero_rows = frozenset(zero_rows)

    @classmethod
    def from_files(cls, files_spec, base_dir, verify_sha=True):
        ops, den_polys, prov = {}, {}, {}
        vars_ = None
        n_basis = None
        zero_rows = None
        for fkey in sorted(files_spec):
            spec = files_spec[fkey]
            op = fkey[2:] if fkey.startswith("A_") else fkey
            path = spec["path"]
            if not os.path.isabs(path):
                path = os.path.normpath(os.path.join(base_dir, path))
            if verify_sha:
                got = _sha256_file(path)
                want = spec.get("sha256")
                if want is None or got != want:
                    raise ValueError(
                        f"RationalConnection REFUSED: sha256 mismatch/missing pin "
                        f"for {fkey}: file {got[:16]} vs pin "
                        f"{str(want)[:16]} ({path})")
                prov[op] = {"path": path, "sha16": got[:16], "sha_verified": True}
            else:
                prov[op] = {"path": path, "sha_verified": False}
            d = json.load(open(path))
            fv = list(d["vars"])
            if vars_ is None:
                vars_ = fv
            elif fv != vars_:
                raise ValueError(f"vars mismatch across files: {fv} vs {vars_}")
            nb = len(d["basis"])
            if n_basis is None:
                n_basis = nb
            elif nb != n_basis:
                raise ValueError(f"basis length mismatch: {nb} vs {n_basis}")
            # zero-row DECLARATION (fail-closed: a zero dJ row is
            # indistinguishable from a MISSING row in sparse format, so
            # coverage credit is given only to rows the emitting receipt
            # explicitly certifies as zero — identical across files)
            zr = frozenset(d.get("zero_rows", []))
            if zero_rows is None:
                zero_rows = zr
            elif zr != zero_rows:
                raise ValueError("zero_rows declarations differ across files")
            if any(not (0 <= r < nb) for r in zr):
                raise ValueError("zero_rows index outside basis range")
            ent = {}
            for key, (num, den) in d["entries"].items():
                i, j = (int(t) for t in key.split(","))
                if not (0 <= i < nb and 0 <= j < nb):
                    raise ValueError(f"entry {key} outside basis range [0,{nb})")
                if not den:
                    raise ValueError(f"entry {key}: empty denominator — refused")
                nt = tuple((tuple(t[0]), Fraction(t[1])) for t in num)
                dk = _den_key(den)
                if dk not in den_polys:
                    den_polys[dk] = tuple((tuple(t[0]), Fraction(t[1])) for t in den)
                ent[(i, j)] = (nt, dk)
            ops[op] = ent
            prov[op]["n_entries"] = len(ent)
            if spec.get("n_entries") is not None and spec["n_entries"] != len(ent):
                raise ValueError(f"{fkey}: n_entries pin {spec['n_entries']} != {len(ent)}")
        if vars_ is None:
            raise ValueError("RationalConnection REFUSED: no files given")
        for t in (t for dk in den_polys for t in den_polys[dk]):
            if any(e < 0 for e in t[0]):
                raise ValueError("negative exponent in denominator term — refused")
        rc = cls(vars_, n_basis, ops, den_polys, prov, zero_rows or ())
        for op, ent in ops.items():
            bad = {i for (i, _j) in ent} & rc.zero_rows
            if bad:
                raise ValueError(f"rows {sorted(bad)} declared zero_rows but "
                                 f"carry entries in A_{op} — contradiction, refused")
        return rc

    # -------------------------------------------------------------- census
    def closure_census(self):
        rows, cols = set(), set()
        for ent in self.ops.values():
            for (i, j) in ent:
                rows.add(i)
                cols.add(j)
        covered = rows | self.zero_rows
        outside = cols - covered
        return {"n_basis": self.n, "n_rows_covered": len(covered),
                "n_rows_with_entries": len(rows),
                "n_zero_rows_declared": len(self.zero_rows),
                "n_cols_touched": len(cols),
                "n_cols_outside_rows": len(outside),
                "closed": (not outside) and bool(rows)}

    # ---------------------------------------------------------- evaluation
    def _term_val(self, terms, pt, dps):
        with mp.workdps(dps):
            s = mpc(0)
            for exps, c in terms:
                v = mpc(mpf(c.numerator) / c.denominator)
                for var, e in zip(self.vars, exps):
                    if e:
                        v *= pt[var] ** e
                s += v
            return s

    def eval_at(self, pt, dps):
        out = {}
        with mp.workdps(dps + 5):
            dvals = {dk: self._term_val(self.den_polys[dk], pt, dps + 5)
                     for dk in self.den_polys}
            for op, ent in self.ops.items():
                m = {}
                for (i, j), (num, dk) in ent.items():
                    dv = dvals[dk]
                    if dv == 0:
                        raise ValueError(f"eval_at ON a denominator zero ({op} {i},{j})")
                    m[(i, j)] = self._term_val(num, pt, dps + 5) / dv
                out[op] = m
        return out

    # ----------------------------------------------------- segment pullback
    def segment_desys(self, p0, p1, chart_vars, subs_f, dps):
        """Exact pullback DESystem on the straight chart segment p0 -> p1
        (t in [0,1]); p0/p1: {var: (Fr,Fr)}; subs_f pins non-chart vars."""
        affine = {}
        for v in self.vars:
            if v in subs_f:
                affine[v] = (subs_f[v], _GQ0)
            else:
                a0 = p0[v]
                affine[v] = (a0, (_gq_add(p1[v], (-a0[0], -a0[1]))))
        pow_cache = {}

        def var_pow(v, e):
            key = (v, e)
            if key in pow_cache:
                return pow_cache[key]
            if e == 0:
                r = [_GQ1]
            else:
                base = _poly_strip([affine[v][0], affine[v][1]])
                r = base
                for _ in range(e - 1):
                    r = _poly_mul(r, base)
            pow_cache[key] = r
            return r

        def compose(terms, scale=_GQ1):
            acc = [_GQ0]
            for exps, c in terms:
                p = [(_gq_mul(scale, (c, Fraction(0))))]
                for v, e in zip(self.vars, exps):
                    if e:
                        p = _poly_mul(p, var_pow(v, e))
                if len(p) > len(acc):
                    acc = acc + [_GQ0] * (len(p) - len(acc))
                for k2, a in enumerate(p):
                    acc[k2] = _gq_add(acc[k2], a)
            return _poly_strip(acc)

        # denominators used by any entry on this segment
        den_t = {}
        for dk in self.den_polys:
            den_t[dk] = compose(self.den_polys[dk])
            if len(den_t[dk]) == 1 and _gq_is0(den_t[dk][0]):
                raise ValueError(f"segment denominator identically zero: {dk}")
        # contributions: (i, j, num_poly_t, den_key); dvar/dt folded into num
        contribs = []
        for op, ent in self.ops.items():
            dv = affine[op][1] if op in chart_vars else _GQ0
            if _gq_is0(dv):
                continue    # op var constant along segment: no contribution
            for (i, j), (num, dk) in ent.items():
                np_t = compose(num, scale=dv)
                if len(np_t) == 1 and _gq_is0(np_t[0]):
                    continue
                contribs.append((i, j, np_t, dk))
        used_keys = {c[3] for c in contribs}
        return _RatSegDE(self, den_t, contribs, used_keys, dps)

    def eval_poly_gq(self, p, x, dps):
        with mp.workdps(dps):
            acc = mpc(0)
            for c in reversed(p):
                acc = acc * x + _gq_to_mpc(c, dps)
            return acc


class _RatSegDE:
    """wayfinder DESystem for one rational-block segment (duck-typed
    contract: .n, .var, .A, .A_series, .singular_points, .meta)."""
    var = "t"

    def __init__(self, rc, den_t, contribs, used_keys, dps):
        self.rc = rc
        self.n = rc.n
        self.den_t = den_t
        self.contribs = contribs
        self.used_keys = used_keys
        wp = dps + 15
        self.roots = {}
        with mp.workdps(wp + 10):
            for dk in used_keys:
                p = den_t[dk]
                deg = len(p) - 1
                if deg == 0:
                    self.roots[dk] = []
                    continue
                coeffs = [_gq_to_mpc(c, wp + 10) for c in reversed(p)]
                try:
                    rts = mp.polyroots(coeffs, maxsteps=200, extraprec=60)
                except Exception as e:
                    raise ValueError(
                        f"segment denominator root-finding failed for {dk}: "
                        f"{type(e).__name__}: {e} — REFUSED (cannot certify "
                        "clearance without the roots)")
                self.roots[dk] = list(rts)
        self.singular_points = sorted({r for rl in self.roots.values() for r in rl},
                                      key=lambda z: (mp.re(z), mp.im(z)))
        self.meta = {"source": "famhar rational segment", "format": "rational-block",
                     "singular_points": self.singular_points}
        self._mpc_cache = {}

    def clearance(self):
        """(min distance of any denominator root to the real segment [0,1],
        den_key). Distance in t-units."""
        worst = (mp.inf, None)
        for dk, rl in self.roots.items():
            for r in rl:
                re_, im_ = mp.re(r), mp.im(r)
                if 0 <= re_ <= 1:
                    dd = abs(im_)
                else:
                    dd = min(abs(r), abs(r - 1))
                if dd < worst[0]:
                    worst = (dd, dk)
        return worst if worst[1] is not None else (mp.inf, "(no denominator roots)")

    def windings(self, dps):
        """Continuous delta-arg of each used denominator along t in [0,1]:
        exact per linear factor (t - r) via the case-3 tracker (each factor
        image IS a straight segment)."""
        out = {}
        with mp.workdps(dps):
            for dk in self.used_keys:
                w = mpf(0)
                for r in self.roots[dk]:
                    w += Family._delta_arg(mpc(0) - r, mpc(1) - r)
                out[dk] = w
        return out

    def den_values_at(self, t, dps):
        return {dk: self.rc.eval_poly_gq(self.den_t[dk], mpc(t), dps)
                for dk in self.used_keys}

    def _conv(self, dps):
        if dps in self._mpc_cache:
            return self._mpc_cache[dps]
        with mp.workdps(dps):
            dens = {dk: [_gq_to_mpc(c, dps) for c in self.den_t[dk]]
                    for dk in self.used_keys}
            cons = [(i, j, [_gq_to_mpc(c, dps) for c in np_t], dk)
                    for (i, j, np_t, dk) in self.contribs]
        self._mpc_cache[dps] = (dens, cons)
        return self._mpc_cache[dps]

    def A(self, x, eps, dps_):
        dens, cons = self._conv(dps_)
        with mp.workdps(dps_):
            xv = mpc(x)
            dvals = {}
            for dk, p in dens.items():
                acc = mpc(0)
                for c in reversed(p):
                    acc = acc * xv + c
                if acc == 0:
                    raise ValueError("A evaluated ON a denominator zero")
                dvals[dk] = 1 / acc
            Amat = [[mpc(0)] * self.n for _ in range(self.n)]
            for (i, j, np_c, dk) in cons:
                acc = mpc(0)
                for c in reversed(np_c):
                    acc = acc * xv + c
                Amat[i][j] += acc * dvals[dk]
            return Amat

    def A_series(self, z0, eps, dps_, M):
        """Taylor coefficients [A_0..A_M] about z0: exact-coefficient shift
        of num/den + inverse-series division per distinct denominator."""
        dens, cons = self._conv(dps_ + 10)
        with mp.workdps(dps_ + 10):
            z0v = mpc(z0)

            def shift(p):
                # Taylor coeffs of poly about z0: repeated synthetic
                # division by (t - z0) on ascending coefficients
                c = list(p)
                out = []
                for _ in range(len(p)):
                    rem = c[-1]
                    q = []
                    for k2 in range(len(c) - 2, -1, -1):
                        q.append(rem)
                        rem = rem * z0v + c[k2]
                    out.append(rem)
                    c = list(reversed(q))
                    if not c:
                        break
                return out

            inv = {}
            for dk, p in dens.items():
                a = shift(p)[:M + 1]
                if a[0] == 0:
                    raise ValueError("A_series requested AT a denominator zero")
                b = [1 / a[0]]
                for k2 in range(1, M + 1):
                    s = mpc(0)
                    for m2 in range(1, min(k2, len(a) - 1) + 1):
                        s += a[m2] * b[k2 - m2]
                    b.append(-b[0] * s)
                inv[dk] = b
            zero = mpc(0)
            out = [[[zero] * self.n for _ in range(self.n)] for _ in range(M + 1)]
            for (i, j, np_c, dk) in cons:
                ns = shift(np_c)[:M + 1]
                b = inv[dk]
                for k2 in range(M + 1):
                    s = mpc(0)
                    for m2 in range(0, min(k2, len(ns) - 1) + 1):
                        s += ns[m2] * b[k2 - m2]
                    out[k2][i][j] = out[k2][i][j] + s
            return out


# ---------------------------------------------------------------------------
# SELFTEST battery (public, self-contained). Synthetic families whose closed
# forms make every check free — no reference data needed, toy scale (n <= 4,
# dps 25, seconds). Every reference below is derived from the FUNCTION
# representation itself (the monodromy of log a is log a + 2*pi*i), never
# from exp(2*pi*i*Res A) — the path-ordering law.
#
# Legs:
#   B1 straight march: values + A.J chart/kinematic derivatives vs closed
#      forms (masters {1, log a, log^2 a / 2, log b}, letters {a, b})
#   B2 identity loop: zero windings, J returns componentwise
#   B3 branch-locus loop around a=0: winding exactly 2*pi; J_end vs the
#      exact continuation of the closed-form basis
#   B4 polynomial letter (w = 1 - a*b, canonical-terms schema): per-segment
#      exact composition + root-finding — value, A.J derivative, locus-loop
#      monodromy, degree<=1 fold
#   B5 refusal laws: near-locus / through-locus / endpoint-on-locus refuse
#      with the named message; polynomial marches refuse floating waypoints
#   B6 ladder sheet fixture: the vendored reference plugin's ell on the real
#      lambda>0 sheet against the fixture record (skipped by name if the
#      plugin or the fixture is absent)
#
# Verdict contract: each leg registers itself as it opens (leg) or is skipped
# (skip_leg); the verdict line's counts are read from that registry. PASS_ALL
# (rc 0) only at 0 skipped and 0 failed; FAIL (rc 1) on any failed check;
# PASS_EXECUTED_WITH_SKIPS (rc 2) when no check failed but a leg was skipped,
# the skipped legs named on the line. There are no expected-fail legs.
# ---------------------------------------------------------------------------
def _selftest():
    """Run the battery; return the process exit code.

    0 = PASS_ALL: every leg executed and every check passed (0 skipped, 0 failed).
    1 = FAIL: at least one check failed (each failing check is printed by name).
    2 = PASS_EXECUTED_WITH_SKIPS: no check failed but at least one leg was skipped
        (the skip printed as `SKIP: <reason>` under the leg's [Bn] banner and the
        skipped legs named on the verdict line) -- the battery did not run in full.
    The executed/skipped counts on the verdict line are read from the leg registry
    that every leg writes as it opens (leg) or is skipped (skip_leg); no count is
    typed. The battery has no expected-fail legs.
    """
    import time
    dps = 25
    bar = 1e-15
    ok = True
    t_all = time.time()
    legs = []   # the leg registry: [leg id, 'executed' | 'skipped'], one entry per leg()

    def leg(leg_id, title):
        """Open a leg: print its banner and register it (executed unless skip_leg follows)."""
        legs.append([leg_id, "executed"])
        print(f"[{leg_id}] {title}", flush=True)

    def skip_leg(reason):
        """Mark the open leg skipped, printing the reason by name (its checks did not run)."""
        legs[-1][1] = "skipped"
        print(f"  SKIP: {reason}", flush=True)

    def check(label, err, b=None):
        nonlocal ok
        good = float(abs(err)) < (b if b is not None else bar)
        print(f"  {label}: err={float(abs(err)):.2e}  {'PASS' if good else 'FAIL'}",
              flush=True)
        if not good:
            ok = False

    def expect_refusal(label, fn, needle):
        nonlocal ok
        try:
            fn()
        except ValueError as e:
            hit = needle in str(e)
            print(f"  {label}: refused, message {'named' if hit else 'WRONG'}  "
                  f"{'PASS' if hit else 'FAIL'}", flush=True)
            if not hit:
                ok = False
            return
        print(f"  {label}: NO refusal  FAIL", flush=True)
        ok = False

    # -- family A: masters J = [1, log a, log^2(a)/2, log b], letters {a, b}
    cfgA = {
        "name": "selftest_log", "chart_vars": ["a", "b"], "eps0": "0",
        "n_masters": 4, "masters": [{}, {}, {}, {}],
        "letters": [{"name": "a", "const": "0", "a": "1", "b": "0"},
                    {"name": "b", "const": "0", "a": "0", "b": "1"}],
        "connection": {"a": {"1,0": "1", "2,1": "1"}, "b": {"3,0": "1"}},
        "anchor": {"point": {"a": "1/4", "b": "1/9"}, "values_plugin": ""},
        "kinematics": {"vars": ["u", "v"],
                       "map": {"u": "a*b", "v": "a/b"},
                       "jacobian": {"du_da": "b", "du_db": "a",
                                    "dv_da": "1/b", "dv_db": "-a/b**2"}},
        "targets": {"T": {"coeffs": {"2": "1"}, "prefactor": "a*b",
                          "dprefactor": {"a": "b", "b": "a"}}},
        "clearance_min": "1/50"}
    famA = Family(cfgA)
    with mp.workdps(dps + 15):
        la0 = mp.log(mpf(1) / 4)
        lb0 = mp.log(mpf(1) / 9)
        y0A = [mpc(1), la0, la0 ** 2 / 2, lb0]

    leg("B1", "straight march: values + A.J derivatives vs closed forms")
    t0 = time.time()
    res = famA.evaluate([{"a": "1/4", "b": "1/9"}, {"a": "3/5", "b": "2/7"}],
                        y0A, dps)
    with mp.workdps(dps + 10):
        a1 = mpf(3) / 5
        b1 = mpf(2) / 7
        la1 = mp.log(a1)
        lb1 = mp.log(b1)
        check("J[log a]", res["J"][1] - la1)
        check("J[log^2 a/2]", res["J"][2] - la1 ** 2 / 2)
        check("J[log b]", res["J"][3] - lb1)
        T = res["targets"]["T"]
        check("target value", T["value"] - a1 * b1 * la1 ** 2 / 2)
        check("d_chart a", T["d_chart"]["a"] - (b1 * la1 ** 2 / 2 + b1 * la1))
        check("d_chart b", T["d_chart"]["b"] - a1 * la1 ** 2 / 2)
        # closed form in the invariants: T = u*(log u + log v)^2/8
        u1 = a1 * b1
        v1 = a1 / b1
        lsum = mp.log(u1) + mp.log(v1)
        check("d_kin u", T["d_kin"]["u"] - (lsum ** 2 / 8 + lsum / 4))
        check("d_kin v", T["d_kin"]["v"] - u1 * lsum / (4 * v1))
    print(f"  wall {time.time() - t0:.1f}s", flush=True)

    leg("B2", "identity loop: zero windings, J returns")
    t0 = time.time()
    loop0 = [{"a": "1/4", "b": "1/9"}, {"a": "3/8+1/8i", "b": "1/9"},
             {"a": "1/2", "b": "1/9"}, {"a": "1/4", "b": "1/9"}]
    J, rec = famA.march(loop0, y0A, dps)
    with mp.workdps(dps + 10):
        check("winding a", rec["windings"]["a"])
        check("winding b", rec["windings"]["b"])
        check("J return", max(abs(J[i] - y0A[i]) for i in range(4)))
    print(f"  wall {time.time() - t0:.1f}s", flush=True)

    leg("B3", "loop around a=0: winding 2*pi; monodromy vs closed-form continuation")
    t0 = time.time()
    loop1 = [{"a": "1/4", "b": "1/9"}, {"a": "1/4i", "b": "1/9"},
             {"a": "-1/4", "b": "1/9"}, {"a": "-1/4i", "b": "1/9"},
             {"a": "1/4", "b": "1/9"}]
    J, rec = famA.march(loop1, y0A, dps)
    with mp.workdps(dps + 10):
        tpi = 2 * mp.pi * mpc(0, 1)
        check("winding a - 2*pi", rec["windings"]["a"] - 2 * mp.pi)
        check("winding b", rec["windings"]["b"])
        # log a -> log a + 2*pi*i; log^2(a)/2 continues accordingly; log b fixed
        check("mono J[log a]", J[1] - (y0A[1] + tpi))
        check("mono J[log^2 a/2]", J[2] - (y0A[2] + tpi * y0A[1] + tpi ** 2 / 2))
        check("mono J[log b]", J[3] - y0A[3])
    print(f"  wall {time.time() - t0:.1f}s", flush=True)

    # -- family B: polynomial letter w = 1 - a*b; J = [1, log a, log(1-a*b)]
    cfgB = {
        "name": "selftest_poly", "chart_vars": ["a", "b"], "eps0": "0",
        "n_masters": 3, "masters": [{}, {}, {}],
        "letters": [{"name": "a", "const": "0", "a": "1", "b": "0"},
                    {"name": "w", "terms": [[[0, 0], "1"], [[1, 1], "-1"]]}],
        "connection": {"a": {"1,0": "1"}, "w": {"2,0": "1"}},
        "anchor": {"point": {"a": "1/8", "b": "1/8"}, "values_plugin": ""},
        "targets": {"P": {"coeffs": {"2": "1"}, "prefactor": "1+0*a",
                          "dprefactor": {"a": "0*a", "b": "0*a"}}},
        "clearance_min": "1/50"}
    famB = Family(cfgB)
    leg("B4", "polynomial letter: exact composition + per-segment root-finding")
    fold = Letter("lin", {"terms": [[[1, 0], "1"], [[0, 0], "-1/2"]]}, ["a", "b"])
    print(f"  degree<=1 'terms' letter folds to affine: "
          f"{'PASS' if fold.affine and fold.coef['a'] == 1 else 'FAIL'}", flush=True)
    if not (fold.affine and fold.coef["a"] == 1):
        ok = False
    with mp.workdps(dps + 15):
        y0B = [mpc(1), mp.log(mpf(1) / 8), mp.log(1 - mpf(1) / 64)]
    t0 = time.time()
    resB = famB.evaluate([{"a": "1/8", "b": "1/8"}, {"a": "1/2", "b": "1/2"}],
                         y0B, dps)
    with mp.workdps(dps + 10):
        check("J[log(1-a*b)]", resB["J"][2] - mp.log(mpf(3) / 4))
        check("J[log a]", resB["J"][1] - mp.log(mpf(1) / 2))
        # d/da log(1-a*b) = -b/(1-a*b)
        check("d_chart a", resB["targets"]["P"]["d_chart"]["a"]
              - (-(mpf(1) / 2) / (1 - mpf(1) / 4)))
        check("winding w (no crossing)", resB["sheet"]["windings"]["w"])
    # loop around the w-locus (b = 1/2 fixed, a circling a=2)
    with mp.workdps(dps + 15):
        y0L = [mpc(1), mp.log(mpf(3) / 2), mp.log(1 - mpf(3) / 4)]
    loopW = [{"a": "3/2", "b": "1/2"}, {"a": "2-1/2i", "b": "1/2"},
             {"a": "5/2", "b": "1/2"}, {"a": "2+1/2i", "b": "1/2"},
             {"a": "3/2", "b": "1/2"}]
    JW, recW = famB.march(loopW, y0L, dps)
    with mp.workdps(dps + 10):
        tpi = 2 * mp.pi * mpc(0, 1)
        check("winding w - 2*pi", recW["windings"]["w"] - 2 * mp.pi)
        check("winding a", recW["windings"]["a"])
        check("mono J[log(1-a*b)]", JW[2] - (y0L[2] + tpi))
        check("mono J[log a]", JW[1] - y0L[1])
    print(f"  wall {time.time() - t0:.1f}s", flush=True)

    leg("B5", "refusal laws")
    expect_refusal("near-locus segment",
                   lambda: famA.march([{"a": "1/4", "b": "1/9"},
                                       {"a": "1/1000", "b": "1/9"}], y0A, dps),
                   "REFUSED")
    expect_refusal("endpoint ON locus",
                   lambda: famA.march([{"a": "1/4", "b": "1/9"},
                                       {"a": "0", "b": "1/9"}], y0A, dps),
                   "REFUSED")
    expect_refusal("through polynomial locus",
                   lambda: famB.march([{"a": "3/2", "b": "1/2"},
                                       {"a": "5/2", "b": "1/2"}], y0L, dps),
                   "REFUSED")
    expect_refusal("float waypoints on polynomial family",
                   lambda: famB.march([{"a": 0.25, "b": 0.5},
                                       {"a": 0.5, "b": 0.5}], y0B, dps),
                   "EXACT")

    # -- B6: the ladder sheet fixture (the vendored closed-form plugin, if present)
    leg("B6", "ladder sheet fixture: the reference plugin's ell on the real lambda>0 sheet")
    t0 = time.time()
    here = os.path.dirname(os.path.abspath(__file__))
    fixture_path = os.path.join(here, "tests", "fixtures", "ladder_sheet_receipt.json")
    plugin_path = os.path.join(here, "ladder_reference.py")
    if not (os.path.exists(plugin_path) and os.path.exists(fixture_path)):
        missing = [p for p in (plugin_path, fixture_path) if not os.path.exists(p)]
        skip_leg(f"absent {', '.join(os.path.basename(p) for p in missing)}")
    else:
        import importlib.util
        spec = importlib.util.spec_from_file_location("ladder_reference", plugin_path)
        LR = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(LR)
        with open(fixture_path) as fh:
            fix = json.load(fh)

        def rel_digits(x, y):
            with mp.workdps(200):
                x = mpf(x) if isinstance(x, str) else x
                y = mpf(y) if isinstance(y, str) else y
                return 999.0 if x == y else float(-mp.log10(abs(x - y) / abs(y)))

        def check_digits(label, d, need):
            nonlocal ok
            good = d >= need
            print(f"  {label}: {d:.2f} digits (need >= {need})  "
                  f"{'PASS' if good else 'FAIL'}", flush=True)
            if not good:
                ok = False

        def check_eq(label, got, want):
            nonlocal ok
            good = got == want
            print(f"  {label}: {got!r} (want {want!r})  {'PASS' if good else 'FAIL'}", flush=True)
            if not good:
                ok = False

        pt = fix["points"][0]
        L6, dps6, bar6, selfbar6 = fix["L"], 50, 38, 48   # bar 38: the record's dps-40 cap (measured 40.3)
        X, Y = Fraction(pt["X"]), Fraction(pt["Y"])
        lam = (1 + X - Y) ** 2 - 4 * X
        check_eq("lambda_exact", str(lam), pt["lambda_exact"])
        with mp.workdps(dps6):
            Xm = mpf(X.numerator) / X.denominator
            Ym = mpf(Y.numerator) / Y.denominator
            delta = 2 * Ym / ((1 - Xm + Ym) + mp.sqrt(mpf(lam.numerator) / lam.denominator))
            a6 = -Xm / ((1 - Xm) - delta)      # w = z/(z-1), negative real on this sheet
            b6 = -(1 - delta) / delta          # wbar = zbar/(zbar-1), negative real
            check_eq("chart a = w (30 digits vs record w_str)", mp.nstr(a6, 30), pt["w_str"])
            check_eq("chart b = wbar (30 digits vs record wbar_str)", mp.nstr(b6, 30), pt["wbar_str"])
            v_k = LR.phi_value(L6, a6, b6, dps6, {"ell_k": -1})   # the record's exact call
            v_d = LR.phi_value(L6, a6, b6, dps6, None)            # the sheet-selected default
            v_0 = LR.phi_value(L6, a6, b6, dps6, {"ell_k": 0})    # the old cont=None branch
            check_digits("record call ell_k=-1 vs own_value_str", rel_digits(mp.re(v_k), pt["own_value_str"]), bar6)
            check("record call |Im|", mp.im(v_k), 1e-40)
            check_digits("default cont=None vs own_value_str", rel_digits(mp.re(v_d), pt["own_value_str"]), bar6)
            check("default |Im|", mp.im(v_d), 1e-40)
            check_digits("default vs record call", rel_digits(mp.re(v_d), mp.re(v_k)), selfbar6)
            check_eq("selected ell_k on the real lambda>0 sheet", LR.default_ell_k(a6, b6, dps6), -1)
            sheet = LR.ell_sheet(a6, b6, dps6)
            check_eq("ell_sheet selected_by", sheet["selected_by"], "default")
            check_digits("wrong-sheet control ell_k=0 vs record cont_none", rel_digits(mp.re(v_0), pt["ladder_reference_cont_none_re_str"]), bar6)
            check("wrong-sheet control |Im| vs record", abs(mp.im(v_0)) - pt["ladder_reference_cont_none_im_abs"], 1e-9 * pt["ladder_reference_cont_none_im_abs"])
            d_wrong = rel_digits(mp.re(v_d), pt["ladder_reference_cont_none_re_str"])
            good = d_wrong < 1 and abs(mp.im(v_d)) < 1e-40
            print(f"  default is NOT the wrong-sheet branch: {d_wrong:.2f} digits vs the old branch  "
                  f"{'PASS' if good else 'FAIL (DEFAULT RETURNED THE WRONG-SHEET BRANCH)'}", flush=True)
            if not good:
                ok = False
            # controls: Euclidean slice (b = conj(a)) and the anchor keep k = 0
            ae = mpc(mpf(1) / 6, -mpf(5) / 6)
            check_eq("Euclidean slice k", LR.default_ell_k(ae, mp.conj(ae), dps6), 0)
            check_eq("Euclidean slice default == log(a)+log(b)", LR.ell_sheet(ae, mp.conj(ae), dps6)["ell"] == mp.log(ae) + mp.log(mp.conj(ae)), True)
            check_eq("anchor k", LR.default_ell_k(mpf(1) / 8, mpf(1) / 9, dps6), 0)
    print(f"  wall {time.time() - t0:.1f}s", flush=True)

    executed = [i for i, s in legs if s == "executed"]
    skipped = [i for i, s in legs if s == "skipped"]
    if not ok:
        verdict, rc = "FAIL", 1
    elif skipped:
        verdict, rc = "PASS_EXECUTED_WITH_SKIPS", 2
    else:
        verdict, rc = "PASS_ALL", 0
    named = f" ({', '.join(skipped)})" if skipped else ""
    print(f"[famhar selftest] {verdict}  (legs {len(executed)} executed / "
          f"{len(skipped)} skipped{named}; dps {dps}, total wall "
          f"{time.time() - t_all:.1f}s)", flush=True)
    return rc


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(
        description="famhar — multi-sheet family evaluator harness "
                    "(import as a library; see README for the config schema)")
    ap.add_argument("--selftest", action="store_true",
                    help="run the self-contained synthetic-family battery "
                         "(exit 0 only if every leg passes)")
    args = ap.parse_args()
    if args.selftest:
        sys.exit(_selftest())
    ap.print_help()

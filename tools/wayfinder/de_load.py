#!/usr/bin/env python3
"""
de_load.py — load saved DE systems into one DESystem interface (wayfinder).

THE POINT
=========
Three on-disk DE formats keep getting re-parsed ad hoc per run. This
module gives them ONE evaluator contract:

    class DESystem:
        n: int                    # system size
        var: str                  # transport variable ('s', 'y', 'eta', ...)
        A(x, eps, dps) -> n x n list[list[mpc]]
        meta: dict                # {'source', 'sha256', 'format', ...}

    load_monomial_json(path, dps_check=50)   # exact-monomial json format
    load_amatrix_json(path)                  # A_of_y.json thiele-fit format
    load_kira_targets(path)                  # kira_target.m (eta,d) rules -> eta-DE

FORMATS (each pinned against the REAL producing/consuming code, not guessed)
============================================================================
monomial-json  written by wayfinder.flintexport (flintexport.py; exact
    fmpq_mpoly arithmetic); pinned against a production connection file
    (a reference fixture not shipped here). Matrix key 'A_<var>' is a
    dict 'i,j' -> {'num': [[k_eps, k_var, 'p/q'], ...], 'den': [...]} of exact
    rational monomials; entry = num(eps,x)/den(eps,x), ALREADY in eps (the
    exporter did the d = 4-2*eps substitution). Aux: companion closure-index
    maps can be attached via closure_json=.

amatrix-json   pinned against a 79x79 production A_of_y.json fit file and
    the code that produced and consumed it
    (schema is theirs, low->high coefficient
    order confirmed against the producer's own evaluator). entries['i,j'] with status 'OK':
        entry(y, d) = N(y; d) / Q(y; d),
        N(y; d) = sum_k n_k(d) y^k,   n_k(d) = P_num(d)/P_den(d)
    where e['N'][k]['num'/'den'] are Fraction-string coefficient lists in d,
    low->high, and same for 'Q' (monic-normalized in the fit). This file is in
    d; A(x, eps, dps) evaluates at d = 4 - 2*eps.

kira-target-m  AMFlow eta-DE reduction output:
    .../target_reduce/results/<family>/kira_target.m from an AMFlow
    working directory. Grammar observed across three production systems
    (up to 7 MB):
        { fam[a1,..,an] ->
           + fam[b1,..,bn]*((poly(d,eta))/(poly(d,eta)))
           + ...
        , ... }
    integer-coefficient polynomials, '^' powers, every term '+'-led. Targets
    are the masters with ONE index raised on an eta-carrying propagator; the
    eta-DE is assembled as
        d/deta I[a] = sum_p (-a_p * c_p) I[a + e_p],   D_p = (...) + c_p*eta,
    (AMFlow writes D_p = momenta^2 - eta, c_p = -1, so the factor is +a_p),
    then each I[a+e_p] is substituted by its kira_target.m rule (or is itself
    a master). Convention CHECK baked into tests/test_de_load.py: for a
    single-mass vacuum bubble subsystem this gives
    A[[1,0],[1,0]] = (d-2)/(2*eta), which matches d/deta log(eta^(d/2-1)) —
    the exact eta-scaling of the massive tadpole. Masters order = the
    results/<family>/masters file order (what AMFlow's boundary vectors use).
    Unhandled constructs raise NotImplementedError naming what is not
    supported — partial-but-correct by design, never guessy.

NUMERICS DISCIPLINE
===================
* Entries are evaluated EXACTLY (fractions.Fraction; python-flint fmpq fast
  path behind try/import for real-rational points) and rounded to dps only at
  the very end, inside mp.workdps — global mp.dps is NEVER touched
  (import-dps footgun).
* mpf/mpc inputs are dyadic rationals, hence converted EXACTLY to Fractions
  (mantissa * 2^exp) — waypoint calls from transport.py lose nothing.
* Complex rational points use exact Gaussian-rational arithmetic (_QQi).
* No sympy anywhere. No fabricated constants: this module only re-expresses
  what is in the source files.
"""
import json
import os
import re
import sys
from fractions import Fraction

import mpmath as mp

try:  # OPTIONAL accel only — everything works on pure Fraction
    import flint as _flint
    _HAVE_FLINT = True
except ImportError:  # pragma: no cover
    _flint = None
    _HAVE_FLINT = False

try:
    from .manifest import sha256_file
except ImportError:  # run outside package (tests, scripts)
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from manifest import sha256_file

try:
    from .ratfun import RF as _RF, padd as _rf_padd, pscale as _rf_pscale
except ImportError:  # flat mode
    from ratfun import RF as _RF, padd as _rf_padd, pscale as _rf_pscale

__all__ = ["DESystem", "load_monomial_json", "load_amatrix_json",
           "load_kira_targets"]

_F0 = Fraction(0)
_F1 = Fraction(1)
_POW_CAP = 10 ** 6  # sanity cap on exponents read from files


# ---------------------------------------------------------------------------
# exact scalar layer: Fraction / fmpq / Gaussian-rational (_QQi)
# ---------------------------------------------------------------------------

class _QQi:
    """Gaussian rational a + b*i with Fraction parts. Exact +,-,*,/,**."""
    __slots__ = ("re", "im")

    def __init__(self, re, im=_F0):
        self.re = Fraction(re)
        self.im = Fraction(im)

    def __add__(self, o):
        return _QQi(self.re + o.re, self.im + o.im)

    def __sub__(self, o):
        return _QQi(self.re - o.re, self.im - o.im)

    def __neg__(self):
        return _QQi(-self.re, -self.im)

    def __mul__(self, o):
        return _QQi(self.re * o.re - self.im * o.im,
                    self.re * o.im + self.im * o.re)

    def __truediv__(self, o):
        d = o.re * o.re + o.im * o.im
        if d == 0:
            raise ZeroDivisionError("QQi division by zero")
        return _QQi((self.re * o.re + self.im * o.im) / d,
                    (self.im * o.re - self.re * o.im) / d)

    def __pow__(self, k):
        if not isinstance(k, int):
            raise TypeError("QQi power must be int")
        if k < 0:
            return _QQi(_F1) / self.__pow__(-k)
        out = _QQi(_F1)
        base = self
        while k:
            if k & 1:
                out = out * base
            base = base * base
            k >>= 1
        return out

    def __eq__(self, o):
        return isinstance(o, _QQi) and self.re == o.re and self.im == o.im

    def __repr__(self):
        return f"_QQi({self.re!r}, {self.im!r})"

    def is_zero(self):
        return self.re == 0 and self.im == 0


def _mpf_to_fraction(v):
    """EXACT mpf -> Fraction (mpf is a dyadic rational: sign*man*2^exp)."""
    sign, man, exp, _bc = v._mpf_
    man = int(man)
    if man == 0:
        if exp == 0:
            return _F0
        raise ValueError(f"non-finite mpf {v!r} cannot be a DE evaluation point")
    fr = Fraction(man << exp) if exp >= 0 else Fraction(man, 1 << (-exp))
    return -fr if sign else fr


def _as_fraction(v, what):
    """Exact conversion of a REAL scalar to Fraction. Raises on lossy input."""
    if isinstance(v, bool):
        raise TypeError(f"{what}: bool is not a number")
    if isinstance(v, (int, Fraction)):
        return Fraction(v)
    if isinstance(v, str):
        try:
            return Fraction(v.strip())
        except (ValueError, ZeroDivisionError) as e:
            raise ValueError(f"{what}: cannot parse {v!r} as exact rational: {e}")
    if isinstance(v, float):
        return Fraction(v)  # exact: float is dyadic
    if isinstance(v, mp.mpf):
        return _mpf_to_fraction(v)
    raise TypeError(f"{what}: unsupported scalar type {type(v).__name__}")


def _as_qqi(v, what):
    """Exact conversion of a real/complex scalar to _QQi."""
    if isinstance(v, _QQi):
        return v
    if isinstance(v, complex):
        return _QQi(Fraction(v.real), Fraction(v.imag))
    if isinstance(v, mp.mpc):
        return _QQi(_mpf_to_fraction(v.real), _mpf_to_fraction(v.imag))
    return _QQi(_as_fraction(v, what))


def _fraction_to_mpf(fr, dps):
    """Fraction -> mpf correctly rounded at dps (computed at dps+10, then
    re-rounded inside workdps(dps)). Never touches global mp.dps."""
    with mp.workdps(dps + 10):
        v = mp.mpf(fr.numerator) / mp.mpf(fr.denominator)
    with mp.workdps(dps):
        return +v


def _qqi_to_mpc(q, dps):
    with mp.workdps(dps):
        return mp.mpc(_fraction_to_mpf(q.re, dps), _fraction_to_mpf(q.im, dps))


def _val_to_qqi(v):
    """Backend value (Fraction | fmpq | _QQi) -> _QQi."""
    if isinstance(v, _QQi):
        return v
    if isinstance(v, Fraction):
        return _QQi(v)
    return _QQi(Fraction(int(v.p), int(v.q)))  # flint fmpq


def _gpow(v, k):
    """Generic integer power for any backend value (fmpq/Fraction/_QQi)."""
    if not isinstance(k, int):
        raise TypeError("exponent must be int")
    if abs(k) > _POW_CAP:
        raise ValueError(f"exponent {k} exceeds sanity cap {_POW_CAP}")
    return v ** k


def _make_backend(xq, eq):
    """Choose the exact backend for one A() call.

    Returns (conv, xv, ev) with conv: Fraction -> backend value, and the
    point/eps already converted. Real-rational points get fmpq when
    python-flint is importable (accel only, bit-identical results);
    complex points fall back to Gaussian rationals.
    """
    real = (xq.im == 0 and eq.im == 0)
    if real and _HAVE_FLINT:
        def conv(fr):
            return _flint.fmpq(fr.numerator, fr.denominator)
        return conv, conv(xq.re), conv(eq.re)
    if real:
        return Fraction, xq.re, eq.re
    def conv(fr):  # noqa: E306
        return _QQi(fr)
    return conv, xq, eq


# ---------------------------------------------------------------------------
# entry evaluators (one class per source format)
# ---------------------------------------------------------------------------

def _poly_mono_eval(monos, ev, xv, conv, pcache_e, pcache_x):
    """sum coeff * eps^ke * x^kx over monomial list, exact, with power memo."""
    tot = None
    for ke, kx, coeff in monos:
        pe = pcache_e.get(ke)
        if pe is None:
            pe = pcache_e[ke] = _gpow(ev, ke)
        px = pcache_x.get(kx)
        if px is None:
            px = pcache_x[kx] = _gpow(xv, kx)
        term = conv(coeff) * pe * px
        tot = term if tot is None else tot + term
    return tot if tot is not None else conv(_F0)


class _MonoEntry:
    """monomial-json entry: num(eps,x)/den(eps,x), exact monomial lists."""
    __slots__ = ("num", "den")

    def __init__(self, num, den):
        self.num = num  # list[(k_eps int, k_x int, Fraction)]
        self.den = den

    def eval(self, xv, ev, dv, conv, pe, px):
        den = _poly_mono_eval(self.den, ev, xv, conv, pe, px)
        num = _poly_mono_eval(self.num, ev, xv, conv, pe, px)
        return num / den  # ZeroDivisionError propagates with context added by caller


def _horner(coeffs, x, conv):
    """Evaluate sum coeffs[k] * x^k (low->high), exact Horner."""
    r = conv(_F0)
    for c in reversed(coeffs):
        r = r * x + conv(c)
    return r


class _AmatrixEntry:
    """A_of_y entry: (sum_k n_k(d) y^k)/(sum_k q_k(d) y^k), n_k,q_k = num/den in d."""
    __slots__ = ("N", "Q")

    def __init__(self, N, Q):
        self.N = N  # list[(num_coeffs list[Fraction], den_coeffs list[Fraction])]
        self.Q = Q

    def _side(self, coeff_pairs, xv, dv, conv):
        tot = None
        xpow = conv(_F1)
        for num_c, den_c in coeff_pairs:
            ck = _horner(num_c, dv, conv) / _horner(den_c, dv, conv)
            term = ck * xpow
            tot = term if tot is None else tot + term
            xpow = xpow * xv
        return tot if tot is not None else conv(_F0)

    def eval(self, xv, ev, dv, conv, pe, px):
        return self._side(self.N, xv, dv, conv) / self._side(self.Q, xv, dv, conv)


class _KiraEntry:
    """eta-DE entry: sum of integer factors * coefficient-AST(d, eta)."""
    __slots__ = ("terms",)

    def __init__(self):
        self.terms = []  # list[(int factor, ast)]

    def eval(self, xv, ev, dv, conv, pe, px):
        env = {"d": dv, "eta": xv}
        tot = None
        for fac, ast in self.terms:
            term = conv(Fraction(fac)) * _ast_eval(ast, env, conv)
            tot = term if tot is None else tot + term
        return tot


# ---------------------------------------------------------------------------
# DESystem
# ---------------------------------------------------------------------------

class DESystem:
    """One loaded DE system: y'(x) = A(x, eps) y(x).

    Attributes (shared wayfinder contract):
        n     system size
        var   transport variable name ('s', 'y', 'eta', ...)
        meta  {'source', 'sha256', 'format', ...} + per-format extras
        A(x, eps, dps) -> n x n list[list[mpc]]

    A() evaluates every stored entry EXACTLY at (x, eps) — Fraction/fmpq for
    real-rational points, Gaussian rationals for complex ones — and rounds to
    dps only at the end, inside mp.workdps. Where the source file is written
    in d, the loader-installed entries use d = 4 - 2*eps. Unstored entries
    are exact zeros. mpf/mpc arguments are dyadic, so their conversion to
    rationals is exact; global mp.dps is never modified.
    """

    def __init__(self, n, var, entries, meta, uses_d):
        self.n = int(n)
        self.var = var
        self.meta = meta
        self._entries = entries  # dict[(i,j)] -> entry object
        self._uses_d = bool(uses_d)
        for (i, j) in entries:
            if not (0 <= i < self.n and 0 <= j < self.n):
                raise ValueError(f"entry index ({i},{j}) out of range for n={self.n}")

    def A(self, x, eps, dps):
        if not (isinstance(dps, int) and dps >= 2):
            raise ValueError(f"dps must be an int >= 2, got {dps!r}")
        xq = _as_qqi(x, f"{self.var} (transport point)")
        eq = _as_qqi(eps, "eps")
        conv, xv, ev = _make_backend(xq, eq)
        if self._uses_d:
            dv = conv(Fraction(4)) - conv(Fraction(2)) * ev
        else:
            dv = None
        pcache_e, pcache_x = {}, {}
        with mp.workdps(dps):
            zero = mp.mpc(0)
        rows = [[zero] * self.n for _ in range(self.n)]
        for (i, j), ent in self._entries.items():
            try:
                val = ent.eval(xv, ev, dv, conv, pcache_e, pcache_x)
            except ZeroDivisionError:
                raise ZeroDivisionError(
                    f"A[{i}][{j}] denominator vanishes at "
                    f"{self.var}={x!r}, eps={eps!r} "
                    f"(source {self.meta.get('source')})")
            rows[i][j] = _qqi_to_mpc(_val_to_qqi(val), dps)
        return rows

    def entry_exact(self, i, j, x, eps):
        """Exact _QQi value of A[i][j] (0 if unstored). Debug/test hook —
        not part of the shared contract, used by tests to check exactness."""
        ent = self._entries.get((i, j))
        if ent is None:
            return _QQi(_F0)
        xq = _as_qqi(x, "x")
        eq = _as_qqi(eps, "eps")
        conv, xv, ev = _make_backend(xq, eq)
        dv = conv(Fraction(4)) - conv(Fraction(2)) * ev if self._uses_d else None
        return _val_to_qqi(ent.eval(xv, ev, dv, conv, {}, {}))

    def enable_fast_path(self, eps, wp=250, backend="mpmath"):
        """Attach the OPTIONAL A_series/singular_points fast-path protocol.

        backend="mpmath" (default): EXACTLY the legacy behavior below —
        nothing changes for existing callers. backend="acb" additionally
        builds a python-flint acb entry table (acbfast.AcbFastTable, same
        RF entries converted mpc->acb exactly) and enables
        transport_fixed_eps(..., backend="acb") — the flint step kernel.
        The mpmath A_series /
        singular_points attached here are IDENTICAL in both modes.

        The fast path (ClosureSystem-style .A_series, entry_to_rf,
        denominator-root
        singular points) is attached here so
        transport.transport_fixed_eps takes its alias-safe
        fast path; this method is the package port, so any caller gets the
        same without monkey-patching. Both production gates passed at full
        stored ~70 digits through that layer.

        What it does: converts every stored entry to a fixed-eps rational
        function num(x)/den(x) with mpc coefficients (ratfun.RF, built at
        working precision `wp`), then sets ON THIS INSTANCE:

            self.A_series(x, eps, dps, M) -> [A_0..A_M] Taylor matrices
                (exact-rational route: shifted-Horner + series division —
                no Cauchy-circle aliasing, no two-radius cross-check cost)
            self.singular_points -> denominator roots (mp.polyroots per
                distinct denominator, deduplicated at 1e-20)

        transport_fixed_eps duck-types both (fast path + provably-capped
        circle radius). The plain .A() path is UNCHANGED — systems that
        never call enable_fast_path behave exactly as before, and
        hasattr(desys, 'A_series') stays False for them (instance
        attribute, not a class method — that is load-bearing for the
        transport duck-typing).

        Contracts / honesty:
        * FIXED eps: A_series refuses (ValueError) an eps that differs from
          the enabled one by more than 10^-dps relative (exact-rational
          comparison) — no silent wrong-eps coefficients.
        * dps <= wp required (coefficients only carry wp digits); raise
          rather than silently under-deliver. Re-enable with a bigger wp
          for deeper targets.
        * RF arithmetic does NOT reduce num/den to lowest terms, so
          singular_points may include APPARENT (cancelled) roots — exactly
          as in the driver. That is conservative for step control; if a
          landing point is spuriously flagged, prune self.singular_points
          after checking A() is finite there.
        * eps enters exactly: entries in d use d = 4 - 2*eps (same
          convention as .A()).

        Returns self (chainable). Calling again rebuilds (new eps/wp).
        """
        if not (isinstance(wp, int) and wp >= 30):
            raise ValueError(f"wp must be an int >= 30, got {wp!r}")
        eq = _as_qqi(eps, "eps (enable_fast_path)")
        with mp.workdps(wp):
            eps_v = _qqi_to_mpc(eq, wp)
            d_v = mp.mpc(4) - 2 * eps_v if self._uses_d else None
            rf_entries = {}
            for (i, j), ent in self._entries.items():
                rf = _entry_to_rf(ent, self.var, eps_v, d_v, wp)
                if not rf.is_zero():
                    rf_entries[(i, j)] = rf
            sings = []
            seen_dens = set()
            seen_roots = set()
            for rf in rf_entries.values():
                den = rf.den
                if len(den) <= 1:
                    continue
                key = tuple(str(c) for c in den)
                if key in seen_dens:
                    continue
                seen_dens.add(key)
                # Root-finding failure is NOT swallowed (unlike the driver):
                # declaring sings DISABLES transport's two-radius alias
                # cross-check, so an under-declared list is a silent-alias
                # risk. Fall back to numpy roots (the recorded precedent —
                # transport_yline.py poles_union uses np.roots for exactly
                # this: RF products carry REPEATED factors, which break
                # Durand-Kerner), polished by findroot where it converges;
                # refuse only if that fails too.
                roots = _den_roots(den)
                # O(N) bucket dedup at ~1e-10 (the numpy fallback route is
                # float64; nearby duplicates across entries collapse —
                # poles_union rounded to 12 decimals for the same reason)
                for r in roots:
                    r = mp.mpc(r)
                    if not (mp.isfinite(r.real) and mp.isfinite(r.imag)):
                        continue
                    key = (round(float(r.real), 10), round(float(r.imag), 10))
                    if key not in seen_roots:
                        seen_roots.add(key)
                        sings.append(r)
            sings = sorted(sings, key=lambda z: (abs(z), z.real, z.imag))

        def _a_series(x, eps_call, dps, M):
            if dps > wp:
                raise ValueError(
                    f"A_series: requested dps={dps} exceeds the fast-path "
                    f"build precision wp={wp} — re-run enable_fast_path "
                    f"with wp >= {dps} (silent under-precision refused)")
            eq_call = _as_qqi(eps_call, "eps (A_series)")
            dre = eq_call.re - eq.re
            dim = eq_call.im - eq.im
            if not (dre == 0 and dim == 0):
                scale2 = max(Fraction(1), eq.re * eq.re + eq.im * eq.im)
                tol2 = Fraction(1, 10 ** (2 * min(dps, wp))) * scale2
                if dre * dre + dim * dim > tol2:
                    raise ValueError(
                        f"A_series: eps={eps_call!r} differs from the "
                        f"enabled fast-path eps={eps!r} beyond 10^-{dps} — "
                        f"rebuild with enable_fast_path(eps) for this eps "
                        f"(fixed-eps layer; silent wrong-eps refused)")
            with mp.workdps(max(dps, wp)):
                xv = mp.mpc(_qqi_to_mpc(_as_qqi(x, "x (A_series)"),
                                     max(dps, wp)))
                out = [[[mp.mpc(0)] * self.n for _ in range(self.n)]
                       for _ in range(M + 1)]
                for (i, j), rf in rf_entries.items():
                    tc = rf.taylor(xv, M)
                    for m in range(M + 1):
                        out[m][i][j] = tc[m]
            return out

        self.A_series = _a_series
        self.singular_points = sings
        self._fastpath = {"eps": eps, "wp": wp, "backend": backend,
                          "n_rf_entries": len(rf_entries),
                          "n_singular_points": len(sings)}
        if backend == "acb":
            try:
                from . import acbfast
            except ImportError:
                import acbfast  # flat mode
            self._acb_fast = acbfast.AcbFastTable(rf_entries, wp)

            def _check_eps(eps_call, dps):
                """Fixed-eps refusal for the acb kernel — same tolerance
                rule as A_series above (silent wrong-eps refused)."""
                eq_call = _as_qqi(eps_call, "eps (acb backend)")
                dre = eq_call.re - eq.re
                dim = eq_call.im - eq.im
                if not (dre == 0 and dim == 0):
                    scale2 = max(Fraction(1),
                                 eq.re * eq.re + eq.im * eq.im)
                    tol2 = Fraction(1, 10 ** (2 * min(dps, wp))) * scale2
                    if dre * dre + dim * dim > tol2:
                        raise ValueError(
                            f"acb backend: eps={eps_call!r} differs from "
                            f"the enabled fast-path eps={eps!r} beyond "
                            f"10^-{dps} — rebuild with enable_fast_path"
                            f"(eps, backend='acb') for this eps")

            self._acb_check_eps = _check_eps
            self._fastpath["n_distinct_dens"] = \
                self._acb_fast.n_distinct_dens
        elif backend != "mpmath":
            raise ValueError(
                f"enable_fast_path: backend must be 'mpmath' or 'acb', "
                f"got {backend!r}")
        return self


# ---------------------------------------------------------------------------
# fast-path entry -> RF conversion (fixed eps)
# (_ast_to_rf + entry_to_rf, generalized: the transport
# variable name comes from the DESystem instead of a hard-coded 'eta',
# and the _MonoEntry/_AmatrixEntry formats — which that driver never loads —
# get their own direct constructors here).
# ---------------------------------------------------------------------------

def _frac_mpc(fr):
    """Fraction -> mpc at the CURRENT working precision (call inside
    mp.workdps; same rounding route as the driver's mp.mpf(num)/mp.mpf(den))."""
    return mp.mpc(mp.mpf(fr.numerator) / mp.mpf(fr.denominator))


def _ast_to_rf(node, d_rf, x_rf, var, wp):
    op = node[0]
    if op == "num":
        return _RF([_frac_mpc(node[1])])
    if op == "sym":
        if node[1] == "d":
            if d_rf is None:
                raise NotImplementedError(
                    "coefficient references 'd' but the system is not in d")
            return d_rf
        if node[1] == var:
            return x_rf
        raise NotImplementedError(
            f"fast path: symbol {node[1]!r} in coefficient (var={var!r})")
    if op == "neg":
        return -_ast_to_rf(node[1], d_rf, x_rf, var, wp)
    if op == "pow":
        return _ast_to_rf(node[1], d_rf, x_rf, var, wp) ** node[2]
    a = _ast_to_rf(node[1], d_rf, x_rf, var, wp)
    b = _ast_to_rf(node[2], d_rf, x_rf, var, wp)
    if op == "add":
        return a + b
    if op == "sub":
        return a - b
    if op == "mul":
        return a * b
    if op == "div":
        return a / b
    raise RuntimeError(f"corrupt AST {op!r}")


def _mono_to_poly(monos, eps_v):
    """Monomial list [(k_eps, k_x, Fraction)] -> coefficient list in x at
    fixed eps (call inside mp.workdps)."""
    deg = max((kx for _, kx, _ in monos), default=0)
    out = [mp.mpc(0)] * (deg + 1)
    ep_cache = {}
    for ke, kx, fr in monos:
        pe = ep_cache.get(ke)
        if pe is None:
            pe = ep_cache[ke] = eps_v ** ke
        out[kx] += _frac_mpc(fr) * pe
    return out


def _entry_to_rf(ent, var, eps_v, d_v, wp):
    """One stored DE entry -> ratfun.RF (num/den polys in the transport
    variable, mpc coefficients at fixed eps). Call inside mp.workdps(wp)."""
    if isinstance(ent, _KiraEntry):
        # driver's entry_to_rf: group terms by identical denominator to
        # curb degree blowup
        d_rf = _RF([d_v])
        x_rf = _RF([mp.mpc(0), mp.mpc(1)])
        groups = {}
        for fac, ast in ent.terms:
            rf = _ast_to_rf(ast, d_rf, x_rf, var, wp)
            num = _rf_pscale(rf.num, mp.mpc(fac))
            key = tuple(str(c) for c in rf.den)
            if key in groups:
                gn, gd = groups[key]
                groups[key] = (_rf_padd(gn, num), gd)
            else:
                groups[key] = (num, rf.den)
        tot = None
        for gn, gd in groups.values():
            t = _RF(gn, gd)
            tot = t if tot is None else tot + t
        return tot if tot is not None else _RF([mp.mpc(0)])
    if isinstance(ent, _MonoEntry):
        return _RF(_mono_to_poly(ent.num, eps_v),
                   _mono_to_poly(ent.den, eps_v))
    if isinstance(ent, _AmatrixEntry):
        def side(pairs):
            out = []
            for num_c, den_c in pairs:
                dv = _horner(den_c, d_v, _frac_mpc_conv)
                if dv == 0:
                    raise ZeroDivisionError(
                        "amatrix coefficient d-denominator vanishes at the "
                        "enabled eps — fast path unavailable at this eps")
                out.append(_horner(num_c, d_v, _frac_mpc_conv) / dv)
            return out if out else [mp.mpc(0)]
        return _RF(side(ent.N), side(ent.Q))
    raise NotImplementedError(
        f"enable_fast_path: entry type {type(ent).__name__} has no RF "
        f"conversion yet — extend _entry_to_rf (plain .A() path unaffected)")


def _frac_mpc_conv(fr):
    """Fraction -> mpc conversion hook for _horner (fast-path route)."""
    return _frac_mpc(Fraction(fr))


def _den_roots(den):
    """All roots of a denominator poly (mpc coeffs, low->high) for the
    singular-point declaration. Primary: mp.polyroots at 50 dps. Fallback on
    NoConvergence (repeated factors from unreduced RF products): np.roots —
    the recorded precedent, transport_yline.py poles_union — polished per-root
    by mp.findroot where it converges. Raises ValueError if every route
    fails (under-declaring sings would silently disable the alias
    cross-check). Root accuracy only steers step clearance / circle radii;
    the transport tail bound still gates every accepted step."""
    coeffs_rev = list(reversed([mp.mpc(c) for c in den]))
    if len(coeffs_rev) <= 6:  # low degree: polyroots, fast + 50-digit roots
        try:
            with mp.workdps(50):
                return mp.polyroots(coeffs_rev, maxsteps=200, extraprec=80)
        except (mp.libmp.NoConvergence, ValueError, ZeroDivisionError):
            pass
    # higher degree, or polyroots refused: numpy route (recorded precedent)
    try:
        import numpy as np
        cf = [complex(c) for c in coeffs_rev]
        raw = np.roots(cf)
    except Exception as exc:  # numpy missing or non-finite coefficients
        raise ValueError(
            f"enable_fast_path: root-finding failed on a degree-"
            f"{len(den) - 1} denominator (polyroots NoConvergence and numpy "
            f"fallback unavailable: {exc}) — refusing to under-declare "
            f"singular_points (that would disable the alias cross-check "
            f"with holes in the clearance map)")
    # float64 accuracy, unpolished — exactly the recorded poles_union
    # discipline (transport_yline.py): root LOCATIONS only steer step
    # clearance / circle radii; the transport tail bound gates every
    # accepted step regardless.
    return [mp.mpc(complex(r)) for r in raw
            if not (np.isnan(r.real) or np.isnan(r.imag)
                    or np.isinf(r.real) or np.isinf(r.imag))]


# ---------------------------------------------------------------------------
# loader 1: monomial-json (exact-monomial connection files)
# ---------------------------------------------------------------------------

def _detect_matrix_key(j, var):
    if var is not None:
        key = f"A_{var}"
        if key not in j:
            raise ValueError(f"no matrix key {key!r} in file (keys: {sorted(j)})")
        return key, var
    cands = [k for k, v in j.items()
             if k.startswith("A_") and isinstance(v, dict)
             and v and all("," in kk for kk in list(v)[:3])]
    if len(cands) != 1:
        raise ValueError(
            f"cannot auto-detect matrix key (candidates {cands}); pass var=")
    return cands[0], cands[0][2:]


def load_monomial_json(path, dps_check=50, var=None, closure_json=None):
    """Load a monomial-json exact-monomial connection.

    Format (exporter: wayfinder.flintexport, flintexport.py): matrix key
    'A_<var>' maps 'i,j' -> {'num': [[k_eps, k_var, 'p/q'], ...], 'den': [...]};
    entry = num(eps, x)/den(eps, x), already written in eps (NOT d). Missing
    keys are exact zeros. n comes from 'NC' (falls back to 1+max index).

    dps_check: if truthy, a load-time sanity pass runs — every denominator is
    checked non-identically-zero (exact), and A() is evaluated once at a
    benign exact-rational probe point at dps_check digits (finite result
    proves plumbing before a long run hangs on it). Set dps_check=0 to skip.

    closure_json: optional path to the closure index-map file
    stored verbatim under meta['closure'].

    The high-precision IC block ('ic'), base_point etc. are carried RAW
    (strings) in meta — parse them with mpf inside workdps at use time
    (import-dps footgun).
    """
    path = os.path.abspath(path)
    with open(path) as f:
        j = json.load(f)
    key, varname = _detect_matrix_key(j, var)
    raw = j[key]
    entries = {}
    maxidx = -1
    for k, rec in raw.items():
        i, jj = (int(t) for t in k.split(","))
        maxidx = max(maxidx, i, jj)
        num = [(int(m[0]), int(m[1]), Fraction(str(m[2]))) for m in rec["num"]]
        den = [(int(m[0]), int(m[1]), Fraction(str(m[2]))) for m in rec["den"]]
        if not den or all(c == 0 for _, _, c in den):
            raise ValueError(f"entry {k}: denominator identically zero")
        entries[(i, jj)] = _MonoEntry(num, den)
    n = int(j.get("NC", maxidx + 1))
    meta = {
        "source": path,
        "sha256": sha256_file(path),
        "format": "monomial-json",
        "matrix_key": key,
        "n_stored_entries": len(entries),
        "eps_convention": "entries already in eps (exporter did d=4-2*eps)",
    }
    for mk in ("description", "provenance", "t", "M", "mt2", "RA", "RB",
               "eps_orders", "A_range", "base_point", "ic", "checks_summary",
               "cone_nu", "cone_sec"):
        if mk in j:
            meta[mk] = j[mk]
    if closure_json is not None:
        cpath = os.path.abspath(closure_json)
        with open(cpath) as f:
            meta["closure"] = json.load(f)
        meta["closure_source"] = cpath
        meta["closure_sha256"] = sha256_file(cpath)
    desys = DESystem(n, varname, entries, meta, uses_d=False)
    if dps_check:
        _monomial_load_check(desys, int(dps_check))
    return desys


def _monomial_load_check(desys, dps_check):
    """One exact probe evaluation; probes shift if a denominator happens to
    vanish at the candidate point (legit pole, not an error)."""
    eps0 = Fraction(1, 97)
    for x0 in (Fraction(1, 3), Fraction(2, 7), Fraction(5, 11), Fraction(7, 13)):
        try:
            rows = desys.A(x0, eps0, dps_check)
        except ZeroDivisionError:
            continue
        with mp.workdps(dps_check):
            for r in rows:
                for v in r:
                    if not (mp.isfinite(v.real) and mp.isfinite(v.imag)):
                        raise ValueError("load check: non-finite entry")
        desys.meta["load_check"] = {
            "probe": {desys.var: str(x0), "eps": str(eps0)}, "dps": dps_check,
            "result": "finite",
        }
        return
    raise ValueError("load check: all probe points hit a denominator zero — "
                     "pass dps_check=0 and probe manually")


# ---------------------------------------------------------------------------
# loader 2: amatrix-json (A_of_y.json family)
# ---------------------------------------------------------------------------

def load_amatrix_json(path, n=None, var="y"):
    """Load an A_of_y.json thiele-fit DE matrix.

    Schema pinned against the production fit code that writes these files
    and the transport code that consumes them. entries['i,j'] (status 'OK') hold
        'N': [ {num: [c0,c1,..], den: [..]} , ...]   # y-coeffs, low->high
        'Q': [ ... ]                                 # same for denominator
    with each y-coefficient a rational function of d given by Fraction-string
    d-polynomial coefficient lists (low->high; fit_A_of_y.ev() is the
    reference Horner). Entry value: N(y;d)/Q(y;d). File is in d; DESystem.A
    evaluates at d = 4 - 2*eps (so the saved d_fit prime points d = 4 - 2/p
    correspond to eps = 1/p). 'zero_entries' and unlisted keys are exact
    zeros; a non-empty 'failures' list means the matrix is incomplete and the
    loader refuses it.

    n: override system size (default: 1 + max index seen; the source file
    does not store n — transport_yline reads it from its boundary file).
    """
    path = os.path.abspath(path)
    with open(path) as f:
        j = json.load(f)
    failures = j.get("failures") or []
    if failures:
        raise ValueError(
            f"{path}: {len(failures)} fit failures present "
            f"(first: {failures[0]}) — DE matrix incomplete, refusing to load")
    entries = {}
    maxidx = -1
    bad_status = []
    for k, rec in j.get("entries", {}).items():
        if rec.get("status") != "OK":
            bad_status.append((k, rec.get("status")))
            continue
        i, jj = (int(t) for t in k.split(","))
        maxidx = max(maxidx, i, jj)
        N = [([Fraction(c) for c in nd["num"]], [Fraction(c) for c in nd["den"]])
             for nd in rec["N"]]
        Q = [([Fraction(c) for c in qd["num"]], [Fraction(c) for c in qd["den"]])
             for qd in rec["Q"]]
        for side, name in ((N, "N"), (Q, "Q")):
            for num_c, den_c in side:
                if not den_c or all(c == 0 for c in den_c):
                    raise ValueError(f"entry {k}: {name}-coefficient with "
                                     f"identically-zero d-denominator")
        if not Q:
            raise ValueError(f"entry {k}: empty Q")
        entries[(i, jj)] = _AmatrixEntry(N, Q)
    if bad_status:
        raise ValueError(f"{path}: non-OK entries in 'entries' dict: "
                         f"{bad_status[:5]} — refusing to load")
    for k in j.get("zero_entries", []):
        i, jj = (int(t) for t in k.split(","))
        maxidx = max(maxidx, i, jj)
    if n is None:
        n = maxidx + 1
    meta = {
        "source": path,
        "sha256": sha256_file(path),
        "format": "amatrix-json",
        "d_convention": "file in d; A(x, eps) evaluated at d = 4 - 2*eps",
        "masters_hash": j.get("masters_hash"),
        "n_points": j.get("n_points"),
        "d_fit": j.get("d_fit"),
        "d_verify": j.get("d_verify"),
        "zero_entries": j.get("zero_entries", []),
        "rescue": j.get("rescue"),
        "n_stored_entries": len(entries),
        "n_inferred": (n == maxidx + 1),
    }
    return DESystem(n, var, entries, meta, uses_d=True)


# ---------------------------------------------------------------------------
# loader 3: kira_target.m -> eta-DE
# ---------------------------------------------------------------------------
# Scalar expression grammar (recursive descent; observed across three
# production kira_target.m files):
#   expr   := term { ('+'|'-') term }
#   term   := factor { ('*'|'/') factor }
#   factor := ('+'|'-') factor | primary [ '^' [sign] INT ]
#   primary:= INT | SYMBOL | '(' expr ')'
# Anything else -> NotImplementedError (partial-but-correct policy).

_TOKEN_RE = re.compile(r"\s*(?:(\d+)|([A-Za-z_]\w*)|(.))")


def _tokenize(s):
    toks = []
    pos = 0
    while pos < len(s):
        m = _TOKEN_RE.match(s, pos)
        if not m or m.end() == pos:
            break
        pos = m.end()
        if m.group(1) is not None:
            toks.append(("int", int(m.group(1))))
        elif m.group(2) is not None:
            toks.append(("sym", m.group(2)))
        else:
            ch = m.group(3)
            if ch in "+-*/^()":
                toks.append((ch, ch))
            else:
                raise NotImplementedError(
                    f"kira coefficient grammar: unexpected character {ch!r} "
                    f"in {s[:80]!r}")
    return toks


class _ExprParser:
    def __init__(self, toks, src):
        self.toks = toks
        self.i = 0
        self.src = src

    def peek(self):
        return self.toks[self.i][0] if self.i < len(self.toks) else None

    def take(self):
        t = self.toks[self.i]
        self.i += 1
        return t

    def expect(self, kind):
        if self.peek() != kind:
            raise NotImplementedError(
                f"kira coefficient grammar: expected {kind!r} at token "
                f"{self.i} of {self.src[:80]!r}")
        return self.take()

    def parse(self):
        node = self.expr()
        if self.i != len(self.toks):
            raise NotImplementedError(
                f"kira coefficient grammar: trailing tokens in {self.src[:80]!r}")
        return node

    def expr(self):
        node = self.term()
        while self.peek() in ("+", "-"):
            op = self.take()[0]
            rhs = self.term()
            node = ("add" if op == "+" else "sub", node, rhs)
        return node

    def term(self):
        node = self.factor()
        while self.peek() in ("*", "/"):
            op = self.take()[0]
            rhs = self.factor()
            node = ("mul" if op == "*" else "div", node, rhs)
        return node

    def factor(self):
        if self.peek() in ("+", "-"):
            op = self.take()[0]
            node = self.factor()
            return node if op == "+" else ("neg", node)
        node = self.primary()
        if self.peek() == "^":
            self.take()
            sign = 1
            if self.peek() in ("+", "-"):
                if self.take()[0] == "-":
                    sign = -1
            kind, val = self.expect("int")
            node = ("pow", node, sign * val)
        return node

    def primary(self):
        kind = self.peek()
        if kind == "int":
            return ("num", Fraction(self.take()[1]))
        if kind == "sym":
            return ("sym", self.take()[1])
        if kind == "(":
            self.take()
            node = self.expr()
            self.expect(")")
            return node
        raise NotImplementedError(
            f"kira coefficient grammar: unexpected token {kind!r} "
            f"in {self.src[:80]!r}")


def _ast_eval(node, env, conv):
    op = node[0]
    if op == "num":
        return conv(node[1])
    if op == "sym":
        try:
            return env[node[1]]
        except KeyError:
            raise NotImplementedError(
                f"kira coefficient references symbol {node[1]!r}; only "
                f"{sorted(env)} are supported. Not supported: symbolic "
                f"kinematic invariants (load_kira_targets would need a "
                f"kinematic-substitution dict for this family).")
    if op == "neg":
        return -_ast_eval(node[1], env, conv)
    if op == "pow":
        return _gpow(_ast_eval(node[1], env, conv), node[2])
    a = _ast_eval(node[1], env, conv)
    b = _ast_eval(node[2], env, conv)
    if op == "add":
        return a + b
    if op == "sub":
        return a - b
    if op == "mul":
        return a * b
    if op == "div":
        return a / b
    raise RuntimeError(f"corrupt AST node {op!r}")


def _split_toplevel(s, seps):
    """Split on chars in seps at zero ()/[] depth. Returns list of chunks."""
    parts, cur = [], []
    dp = db = 0
    for ch in s:
        if ch == "(":
            dp += 1
        elif ch == ")":
            dp -= 1
        elif ch == "[":
            db += 1
        elif ch == "]":
            db -= 1
        if ch in seps and dp == 0 and db == 0:
            parts.append("".join(cur))
            cur = []
        else:
            cur.append(ch)
    parts.append("".join(cur))
    return parts


_INTEGRAL_RE = re.compile(r"^\s*([A-Za-z_]\w*)\s*\[([^\]]*)\]\s*$", re.S)


def _parse_integral(s, where):
    m = _INTEGRAL_RE.match(s)
    if not m:
        raise NotImplementedError(
            f"kira grammar ({where}): cannot parse integral from {s[:120]!r}")
    name = m.group(1)
    try:
        idx = tuple(int(t) for t in m.group(2).split(","))
    except ValueError:
        raise NotImplementedError(
            f"kira grammar ({where}): non-integer index in {s[:120]!r} — "
            f"symbolic indices are not supported")
    return name, idx


def _split_terms_signed(rhs):
    """RHS -> list of (sign, chunk) split at top-level +/-."""
    out = []
    cur = []
    sign = 1
    dp = db = 0
    started = False
    for ch in rhs:
        if ch == "(":
            dp += 1
        elif ch == ")":
            dp -= 1
        elif ch == "[":
            db += 1
        elif ch == "]":
            db -= 1
        if ch in "+-" and dp == 0 and db == 0 and not _in_token(cur):
            if started and "".join(cur).strip():
                out.append((sign, "".join(cur)))
                cur = []
            sign = 1 if ch == "+" else -1
            started = True
            continue
        cur.append(ch)
        if ch.strip():
            started = True
    if "".join(cur).strip():
        out.append((sign, "".join(cur)))
    return out


def _in_token(cur):
    """True if the char stream is mid-exponent like '^-' (not seen in the
    wild for kira_target.m, but cheap to guard)."""
    tail = "".join(cur).rstrip()
    return tail.endswith("^")


def _parse_kira_rules(text, path):
    """kira_target.m text -> dict[idx_tuple] -> (family, [(idx_tuple, ast)])."""
    body = text.strip()
    if not (body.startswith("{") and body.endswith("}")):
        raise NotImplementedError(
            f"{path}: expected '{{ rules }}' wrapper — got "
            f"{body[:40]!r}...{body[-20:]!r}")
    body = body[1:-1]
    rules = {}
    family = None
    for chunk in _split_toplevel(body, ","):
        if not chunk.strip():
            continue
        if "->" not in chunk:
            raise NotImplementedError(
                f"{path}: rule without '->': {chunk[:120]!r}")
        lhs_s, rhs_s = chunk.split("->", 1)
        name, lidx = _parse_integral(lhs_s, "rule LHS")
        if family is None:
            family = name
        elif name != family:
            raise NotImplementedError(
                f"{path}: multiple families in one file ({family!r} vs "
                f"{name!r}) — not supported")
        terms = []
        rhs_clean = rhs_s.strip()
        if rhs_clean == "0":
            rules[lidx] = terms
            continue
        for sign, tchunk in _split_terms_signed(rhs_s):
            tchunk = tchunk.strip()
            if not tchunk:
                continue
            if tchunk == "0":
                continue
            # expected shape: fam[idx]*(coefficient-expression)
            parts = _split_toplevel(tchunk, "*")
            name2, tidx = _parse_integral(parts[0], "rule RHS term")
            if name2 != family:
                raise NotImplementedError(
                    f"{path}: RHS integral family {name2!r} != {family!r}")
            if len(parts) == 1:
                ast = ("num", _F1)
            else:
                coef_src = "*".join(parts[1:])
                ast = _ExprParser(_tokenize(coef_src), coef_src).parse()
            if sign == -1:
                ast = ("neg", ast)
            terms.append((tidx, ast))
        rules[lidx] = terms
    return family, rules


def _parse_masters_file(path):
    masters = []
    with open(path) as f:
        for line in f:
            line = line.split("#", 1)[0].strip()
            if not line:
                continue
            name, idx = _parse_integral(line, f"masters file {path}")
            masters.append((name, idx))
    if not masters:
        raise ValueError(f"{path}: no masters parsed")
    return masters


def _parse_family_yaml(path):
    """Minimal targeted parse of Kira integralfamilies.yaml: name +
    propagator strings. Raises NotImplementedError on constructs outside the
    AMFlow-generated shape (quoted propagators, empty cut_propagators)."""
    with open(path) as f:
        txt = f.read()
    m = re.search(r"-\s*name:\s*\"([^\"]+)\"", txt)
    if not m:
        raise NotImplementedError(f"{path}: no quoted family name found")
    name = m.group(1)
    props = []
    in_props = False
    for line in txt.splitlines():
        ls = line.strip()
        if ls.startswith("propagators:"):
            in_props = True
            continue
        if in_props:
            pm = re.match(r"-\s*\[\s*\"([^\"]+)\"\s*,\s*[^\]]*\]", ls)
            if pm:
                props.append(pm.group(1).strip())
            elif ls.startswith("-"):
                raise NotImplementedError(
                    f"{path}: unparsed propagator line {ls!r}")
            elif ls:
                in_props = False
                if ls.startswith("cut_propagators:") and not ls.endswith("[]"):
                    raise NotImplementedError(
                        f"{path}: non-empty cut_propagators — cut families "
                        f"are not supported (eta-derivative bookkeeping differs on cuts)")
    if not props:
        raise NotImplementedError(f"{path}: no propagators parsed")
    return name, props


_ETA_WORD_RE = re.compile(r"(?<![A-Za-z0-9_])eta(?![A-Za-z0-9_])")
# unit-coefficient additive "+/- eta" term ANYWHERE in the propagator
# (AMFlow SingleMass insertion on an already-
# massive line emits "(...)^2 - eta-1" — eta is NOT trailing). Reject
# multiplicative context after the word (would be a non-unit coefficient).
_ETA_TERM_RE = re.compile(
    r"([+-])\s*eta(?![A-Za-z0-9_])(?!\s*[*^/(])")


def _eta_coefficient(prop):
    """c_p with D_p = (eta-free part) + c_p*eta, for AMFlow-shaped
    propagators. 0 if eta absent. NotImplementedError on any other shape.

    Accepted shapes (a superset of a trailing-only rule):
    bare "eta", and a single signed additive unit term "+ eta"/"- eta"
    anywhere (e.g. "(l)^2 - eta", "(l)^2 - eta-1", "(l)^2 - eta - 1").
    Anything multiplicative ("2*eta", "eta*s", "eta^2") still raises."""
    hits = _ETA_WORD_RE.findall(prop)
    if not hits:
        return 0
    if len(hits) > 1:
        raise NotImplementedError(
            f"propagator {prop!r}: eta appears {len(hits)} times — not supported")
    s = prop.strip()
    if s == "eta":
        return 1
    m = _ETA_TERM_RE.search(s)
    if not m:
        raise NotImplementedError(
            f"propagator {prop!r}: eta not a signed unit additive "
            f"'+/- eta' term — general linear-in-eta propagators are not supported")
    # the signed term must sit at paren depth 0 (top-level additive term);
    # "(a - eta)*b" would otherwise be misread as an additive unit eta.
    depth = 0
    for ch in s[:m.start()]:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
    if depth != 0:
        raise NotImplementedError(
            f"propagator {prop!r}: eta inside parentheses — general "
            f"linear-in-eta propagators are not supported")
    return 1 if m.group(1) == "+" else -1


def load_kira_targets(path, masters_path=None, family_yaml=None, var="eta"):
    """kira_target.m exact-rational (eta, d) reduction rules -> eta-DE DESystem.

    path: .../target_reduce/results/<family>/kira_target.m from an AMFlow
    system_*_diffeq run (pinned on a production two-loop five-point
    run). Companion files located automatically (override via kwargs):
        masters file   <same dir>/masters.final  (fallback: masters)
        family yaml    <target_reduce>/config/integralfamilies.yaml

    Assembly (see module docstring for the derivation + tadpole-scaling
    convention check): with D_p = (eta-free) + c_p*eta from the family yaml,
        d/deta I[a] = sum_p (-a_p * c_p) * I[a + e_p]
    and each raised integral is either itself a master or looked up in the
    kira_target.m rules (which express it over the masters). Basis order =
    masters-file order. File is in d; A(x, eps) evaluates at d = 4 - 2*eps.

    Raises NotImplementedError (naming the unsupported construct) rather than guessing on:
    non-AMFlow eta placement, cut propagators, coefficient symbols other than
    d/eta, or a needed dotted target missing from the rules.
    """
    path = os.path.abspath(path)
    fam_dir = os.path.dirname(path)
    if masters_path is None:
        for cand in (os.path.join(fam_dir, "masters.final"),
                     os.path.join(fam_dir, "masters")):
            if os.path.isfile(cand):
                masters_path = cand
                break
        if masters_path is None:
            raise FileNotFoundError(
                f"no masters/masters.final next to {path}; pass masters_path=")
    if family_yaml is None:
        family_yaml = os.path.join(os.path.dirname(os.path.dirname(fam_dir)),
                                   "config", "integralfamilies.yaml")
        if not os.path.isfile(family_yaml):
            raise FileNotFoundError(
                f"no {family_yaml} (expected AMFlow target_reduce layout); "
                f"pass family_yaml=")
    with open(path) as f:
        family, rules = _parse_kira_rules(f.read(), path)
    masters = _parse_masters_file(masters_path)
    yname, props = _parse_family_yaml(family_yaml)
    if family is None:
        family = masters[0][0]
    for mname, _ in masters:
        if mname != family:
            raise NotImplementedError(
                f"masters file family {mname!r} != rules family {family!r}")
    if yname != family:
        raise NotImplementedError(
            f"integralfamilies.yaml family {yname!r} != rules family "
            f"{family!r} — wrong config dir?")
    eta_coeffs = [_eta_coefficient(p) for p in props]
    if not any(eta_coeffs):
        raise NotImplementedError(
            f"{family_yaml}: no propagator carries eta — not an eta-DE family")
    nprops = len(props)
    midx = {}
    for col, (_, idx) in enumerate(masters):
        if len(idx) != nprops:
            raise NotImplementedError(
                f"master {idx} has {len(idx)} indices but family has "
                f"{nprops} propagators")
        midx[idx] = col
    n = len(masters)
    entries = {}
    missing = []
    for row, (_, a) in enumerate(masters):
        for p in range(nprops):
            cp = eta_coeffs[p]
            if cp == 0 or a[p] == 0:
                continue
            factor = -a[p] * cp  # d/deta D_p^-a_p = (-a_p c_p) D_p^-(a_p+1)
            raised = tuple(a[q] + (1 if q == p else 0) for q in range(nprops))
            if raised in midx:
                contribs = [(midx[raised], ("num", _F1))]
            elif raised in rules:
                contribs = []
                for tidx, ast in rules[raised]:
                    if tidx not in midx:
                        raise NotImplementedError(
                            f"rule for {family}{list(raised)} references "
                            f"non-master {family}{list(tidx)} — incomplete "
                            f"reduction; chaining through intermediate rules is not supported")
                    contribs.append((midx[tidx], ast))
            else:
                missing.append((row, a, p, raised))
                continue
            for col, ast in contribs:
                ent = entries.get((row, col))
                if ent is None:
                    ent = entries[(row, col)] = _KiraEntry()
                ent.terms.append((factor, ast))
    if missing:
        lines = "; ".join(
            f"d/d{var} {family}{list(a)} needs {family}{list(r)} "
            f"(prop {p})" for _, a, p, r in missing[:6])
        raise NotImplementedError(
            f"{path}: {len(missing)} dotted target(s) have no reduction rule "
            f"and are not masters: {lines} — re-run the kira "
            f"target_reduce with these targets, or pass the right "
            f"kira_target.m")
    meta = {
        "source": path,
        "sha256": sha256_file(path),
        "format": "kira-target-m",
        "d_convention": "file in d; A(x, eps) evaluated at d = 4 - 2*eps",
        "family": family,
        "masters": [f"{family}[{','.join(map(str, idx))}]" for _, idx in masters],
        "masters_order": "masters-file order (AMFlow boundary-vector order)",
        "masters_path": os.path.abspath(masters_path),
        "masters_sha256": sha256_file(masters_path),
        "family_yaml": os.path.abspath(family_yaml),
        "family_yaml_sha256": sha256_file(family_yaml),
        "propagators": props,
        "eta_coefficients": eta_coeffs,
        "n_rules": len(rules),
        "n_stored_entries": len(entries),
    }
    return DESystem(n, var, entries, meta, uses_d=True)

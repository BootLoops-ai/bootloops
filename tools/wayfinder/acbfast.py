#!/usr/bin/env python3
"""
acbfast.py — OPT-IN python-flint (Arb) acb backend for the wayfinder
fixed-eps A_series/step kernel.

WHY (measured)
==============
The pure-mpmath Taylor-step transport on a reference 378-dim m2-line DE
(12,162 stored rational entries, mtay=60 at dps 38) measures >= 60-120
s/step (>21 min for a ~10-20-step class leg)
=> a 10-20 h full-wire projection. The whole cost is per-step mpmath series
arithmetic: (a) per-entry RF.taylor (shift + series division), (b) the
O(mtay^2 x nnz) Taylor convolution y_{m+1} = (1/(m+1)) sum_k A_k y_{m-k}.
Both are re-expressed here in flint acb (ball) arithmetic with the Python
call count batched down (measured primitive costs on one workstation, 245-bit:
acb scalar mul ~0.35 us from Python, acb_poly mul 61x61 ~99 us => ~27
ns/coefficient-pair in C; the design keeps per-step Python-level calls at
~nnz x mtay/B with B~4-8 instead of ~nnz x mtay).

DESIGN RULES
============
* Ball discipline — "M=50 +
  ftrim=acb.mid() per step" (arb ball blowup is contained by trimming to
  midpoints per accepted step): Taylor order stays the chassis dps-scaled
  mtay (max(60, 0.75*dps+25) — NEVER lowered here), and the state vector
  is trimmed to ball MIDPOINTS once per ACCEPTED step. Inside a step the
  arithmetic is genuine ball arithmetic; across steps the radii are
  dropped, so results are NOT ball-certified — accuracy is certified the
  same way as the mpmath backend: the chassis geometric tail bound +
  guard (trunc_worst), two-precision digit gates on top. Same honesty
  class as the default backend, ~x0-x00 faster.
* Exactness — the acb entry table is built from the SAME fixed-eps
  ratfun.RF entries as the mpmath fast path (de_load.enable_fast_path),
  converted coefficient-by-coefficient mpc -> acb EXACTLY (mpf mantissas
  are dyadic; no re-rounding). No partial fractions / no float roots: a
  wrong-but-finite A_series is the catalogued silent-wrong class, so the
  per-step Taylor build stays shift + series-division on the exact
  num/den data, in balls.
* Global-context footgun — flint's ctx.prec/ctx.cap are process-global
  (same trap class as mpmath import-dps). Every public entry point here
  saves and restores both in try/finally; nothing leaks.

WHAT LIVES HERE
===============
  mpc_to_acb / acb_to_mpc / acb_mag_mpf   exact dyadic converters
  AcbTrouble                              step-kernel "halve h" signal
  AcbFastTable                            per-DESystem entry table:
      .taylor_entries(z_acb, M) -> dict (i,j) -> (low_coeffs, acb_poly)
      Taylor coefficients of every stored entry about z to order M
      (grouped denominator inversion: distinct shifted dens inverted
      once via acb_series.inv, then one short-num x inv product per
      entry — measured 5,540 distinct dens / 12,162 entries at n=378).

Wired by de_load.DESystem.enable_fast_path(eps, wp, backend="acb") and
consumed by transport._march_leg_acb (transport_fixed_eps(...,
backend="acb")). The default mpmath paths are byte-identical-untouched.
"""

from fractions import Fraction

from mpmath import mp, mpf, mpc

try:
    import flint as _fl
    from flint import acb as _acb, acb_poly as _acb_poly, \
        acb_series as _acb_series, arb as _arb, ctx as _fctx
    HAVE_FLINT = True
except Exception:  # pragma: no cover — module import must stay optional
    _fl = None
    HAVE_FLINT = False

__all__ = ["HAVE_FLINT", "AcbTrouble", "AcbFastTable",
           "mpc_to_acb", "acb_to_mpc", "acb_mag_mpf", "prec_bits"]

def prec_bits(wp):
    """Working decimal digits -> acb bit precision (guarded).

    Exact integer sizing (BootLoops exactness standard): smallest b
    with 10**wp <= 2**b, via pure-integer bit_length — no float log2(10)
    literal in the runtime path. Never undershoots the old float sizing."""
    return (10 ** wp - 1).bit_length() + 20


class AcbTrouble(Exception):
    """acb Taylor build failed at this expansion point (denominator ball
    contains 0 / non-finite entry) — transport halves the step, exactly
    like _CircleTrouble on the generic path."""


# ---------------------------------------------------------------- converters
def _mpf_to_arb(x):
    """EXACT mpf -> arb (mpf is dyadic: sign*man*2^exp). No rounding."""
    sign, man, exp, _bc = x._mpf_
    man = int(man)
    if man == 0:
        if exp == 0:
            return _arb(0)
        raise ValueError(f"non-finite mpf {x!r} cannot enter the acb kernel")
    a = _arb(-man if sign else man)
    return a * (_arb(2) ** int(exp))


def mpc_to_acb(v):
    """EXACT mpc/mpf/int/Fraction-free scalar -> acb (dyadic, no rounding).
    Call with ctx.prec already set high enough that the 2^exp scaling is
    exact (power-of-two scaling is exact in arf at any prec)."""
    if isinstance(v, mpc):
        return _acb(_mpf_to_arb(v.real), _mpf_to_arb(v.imag))
    if isinstance(v, mpf):
        return _acb(_mpf_to_arb(v))
    if isinstance(v, (int,)):
        return _acb(v)
    raise TypeError(f"mpc_to_acb: unsupported type {type(v).__name__}")


def acb_to_mpc(v):
    """acb MIDPOINT -> mpc, exact (man_exp on the exact mid arf; rounding
    happens only when the caller's mp context re-rounds). Call inside
    mp.workdps(...)."""
    zm = v.mid()
    out = []
    for part in (zm.real, zm.imag):
        man, exp = part.man_exp()
        man = int(man)
        out.append(mp.ldexp(mp.mpf(man), int(exp)) if man else mp.mpf(0))
    return mpc(out[0], out[1])


def acb_mag_mpf(v):
    """|acb midpoint| as mpf (exponent-safe — no float over/underflow; deep
    ladder magnitudes reach 10^±1800). Call inside mp.workdps(...)."""
    zm = v.mid()
    re_m, re_e = zm.real.man_exp()
    im_m, im_e = zm.imag.man_exp()
    re = mp.ldexp(mp.mpf(int(re_m)), int(re_e)) if int(re_m) else mp.mpf(0)
    im = mp.ldexp(mp.mpf(int(im_m)), int(im_e)) if int(im_m) else mp.mpf(0)
    if im == 0:
        return abs(re)
    if re == 0:
        return abs(im)
    return mp.hypot(re, im)


# ---------------------------------------------------------------- table
class AcbFastTable:
    """Fixed-eps acb entry table for one DESystem.

    Built ONCE (enable_fast_path(backend="acb")) from the ratfun.RF
    entries: per stored entry the num/den mpc coefficient lists are
    converted exactly to acb_poly. Distinct denominators are grouped
    (identity by exact coefficient string key — same discipline as the
    singular-point dedup) so each shifted denominator is inverted once
    per step, not once per entry.

    taylor_entries(z, M) returns, per stored entry, the Taylor
    coefficients of num/den about z to order M as
        (low, poly): low = list[acb] (the first coefficients, for the
        step kernel's within-block scalar recursion), poly = acb_poly
        (full truncated series, for the cross-block products).
    Raises AcbTrouble when any shifted denominator's constant ball
    contains 0 (expansion point on/too near a pole at this precision) —
    the caller halves the step.
    """

    def __init__(self, rf_entries, wp):
        if not HAVE_FLINT:
            raise ImportError(
                "acb backend requested but python-flint is not importable — "
                "install python-flint or use backend='mpmath'")
        self.wp = int(wp)
        self.prec = prec_bits(self.wp)
        old = _fctx.prec
        _fctx.prec = self.prec
        try:
            dens = {}          # den key -> den index
            den_polys = []     # den index -> acb_poly
            entries = {}       # (i,j) -> (num_poly, den_index)
            for (i, j), rf in rf_entries.items():
                key = tuple(str(c) for c in rf.den)
                di = dens.get(key)
                if di is None:
                    di = dens[key] = len(den_polys)
                    den_polys.append(_acb_poly([mpc_to_acb(c)
                                                for c in rf.den]))
                entries[(i, j)] = (_acb_poly([mpc_to_acb(c)
                                              for c in rf.num]), di)
            self._den_polys = den_polys
            self._entries = entries
            self.n_distinct_dens = len(den_polys)
        finally:
            _fctx.prec = old

    # -- internals ---------------------------------------------------------
    @staticmethod
    def _shift(poly, lin):
        """poly(z + u) as acb_poly in u — Horner with the linear factor
        lin = (u + z). Exact-structure port of ratfun.pshift."""
        cs = poly.coeffs()
        if not cs:
            return _acb_poly([])
        acc = _acb_poly([cs[-1]])
        for c in reversed(cs[:-1]):
            acc = acc * lin + _acb_poly([c])
        return acc

    def taylor_entries(self, z, M, nlow=8, prec=None):
        """Taylor coefficients of every stored entry about z (acb) to order
        M. Returns dict (i,j) -> (low, poly); see class docstring. Work
        precision: `prec` bits if given (capped at the build precision),
        else the build precision — the transport kernel passes its
        wp-scaled prec (mpmath-path parity: that path also rounds every
        recursion op at the transport wp). ctx saved/restored."""
        oldp, oldc = _fctx.prec, _fctx.cap
        _fctx.prec = min(prec, self.prec) if prec else self.prec
        _fctx.cap = M + 1
        try:
            lin = _acb_poly([z, _acb(1)])
            invs = []
            for dp in self._den_polys:
                dsh = self._shift(dp, lin)
                d0 = dsh[0]
                if (not d0.is_finite()) or d0.contains(_acb(0)):
                    raise AcbTrouble(
                        "shifted denominator constant term contains 0 / "
                        "non-finite at this expansion point")
                inv = _acb_series(dsh.coeffs()).inv()
                invs.append(_acb_poly(inv.coeffs()))
            out = {}
            for (i, j), (np_, di) in self._entries.items():
                nsh = self._shift(np_, lin)
                q = nsh * invs[di]           # len <= M + deg_num; extra
                # orders beyond M are dead weight the step kernel ignores
                q = q.truncate(M + 1)
                cs = q.coeffs()
                for c in cs:
                    if not c.is_finite():
                        raise AcbTrouble(
                            f"non-finite Taylor coefficient in entry "
                            f"({i},{j}) at this expansion point")
                low = cs[:nlow]
                out[(i, j)] = (low, q)
            return out
        finally:
            _fctx.prec, _fctx.cap = oldp, oldc

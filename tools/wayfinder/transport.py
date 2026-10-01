#!/usr/bin/env python3
r"""
transport.py — fixed-eps Taylor endpoint transport for wayfinder DESystems.

DESIGN RULES (this module's own laws)
=====================================
* Numerics chassis: scalar-mpc Taylor stepping at FIXED
  eps0 — none of the windowed tool's Laurent aliasing; adaptive step fraction
  starting at ratio, halved on non-convergence and re-grown (x2, capped at
  ratio) after success; per-step term-decay convergence guard; final step
  lands EXACTLY on the target; regularity assert at the landing point;
  hard step cap (400/leg).
* Why fixed-eps at all: the windowed
  Laurent endpoint transport is UNSOUND on a pointwise
  basis (the connection's (d-4)-poles make the Laurent depth unbounded; the
  window aliases). Fixed-eps scalar transport + a Cauchy eps-fan (epsfan.py)
  is the calibrated replacement (R=1e-3 x 32-node circle; R=1e-17 loses ~50d
  to 1/eps pole leverage).
* Complex-waypoint detours (the MUM-crossing detour rule): a
  straight path through/near a singular point collapses the adaptive step
  (a Zeno hang). Fix = explicit detour path, e.g. [x0/4, 0.03i, x1]. This
  module raises (with that fix named) instead of hanging.
* Step control:
    (a) COMPLEX-pole clearance — ALL denominator
        roots (complex included) count, and steps go by the Euclidean
        clearance:
        declared singular points are parsed as mpc and the clearance is
        min |z - s| over them (complex-safe abs). Real-pole-only clearance
        is a catalogued silent-wrong footgun; the regression test
        (tests/test_transport.py::test_complex_pole_pair_off_axis) pins it.
    (b) GEOMETRIC TAIL-BOUND acceptance — the Taylor
        tail is bounded by (t1 + t0)/(1 - r) with the MEASURED ratio r = t1/t0 of the
        last two term magnitudes, clamped at 3/4 (r = 1/2 when t0 == 0),
        halving h until the bound clears tolerance. A single-small-term
        accept is silently wrong when a complex
        pole pair symmetric about the path makes alternate Taylor terms
        vanish (the dip-then-climb trap):
        one tiny term does NOT bound the
        tail. The 3/4 clamp is sound here because
        |h| <= 0.5 x dist-to-singularity when sings are declared, and when
        they are unknown the two-radius cross-check rejects any step with
        a singularity inside the r = 2|h| sampling circle — either way the
        true asymptotic term ratio is <= 1/2 < 3/4.
    (c) ZERO-RUN-AWARE acceptance: feeding the raw LAST TWO terms
        to the sliding early break is unsound on a local Taylor support
        lattice with
        gaps >= 2 (y' = x^k y stepped from x = 0, k >= 2 — exactly
        de_load's monomial-cone shape): two consecutive zero/roundoff terms
        make the bound 0 and ACCEPT a catastrophically truncated step
        while banking a false trunc_worst certificate. So
        (_tail_bound_nz) the bound
        is built from the last two ABOVE-FLOOR term magnitudes t_a, t_b
        (indices a < b) with the gap-normalized ratio
        r = (t_b/t_a)^(1/(b-a)); the early break additionally requires the
        pair to be FRESH (trailing zero-run since b shorter than the
        largest observed inter-term gap), and an all-zero/roundoff
        trailing window NEVER accepts early — it falls through to the
        fixed top order mtay (the fixed-top-order step rule). The
        t_prev == 0 -> r = 1/2 fallback survives ONLY at that final order.
        Regression: tests/test_transport.py::test_sparse_lattice_monomials
        (y'=x^k y, k in {1,2,3,5,10,12}, dps 60/80, certificate honesty).
* Import-dps footgun: every
  high-precision string is parsed inside `with mp.workdps(...)`; global
  mp.dps is NEVER touched.

WHAT IT DOES
============
transport_fixed_eps(desys, eps0, x0, x1, y0, dps, path=None, mtay=None,
                    ratio=0.5, guard_extra=8) -> list[mpc]

March y'(x) = A(x, eps0) y(x) from x0 to x1 (optionally via complex
waypoints) and return y(x1) rounded to dps. `desys` is any object with the
shared wayfinder contract: .n, .var, .A(x, eps, dps) -> n x n list[list].

LOCAL A-SERIES: the contract exposes A only pointwise, so per step the local
Taylor coefficients of A about the expansion point are recovered by a Cauchy
circle (radix-2 FFT over mpc samples). Discipline:
  - if the DESystem declares its singular points (meta['singular_points'] or
    attribute .singular_points), the circle radius is capped at 0.5x the
    distance to the nearest one, so the alias error is <= 2^-N (N >= 3.5*wp
    nodes, wp = dps+guard_extra+25 working digits): provably below the guard.
  - if singular points are UNKNOWN, every step's low-order coefficients are
    recomputed on a second circle at half the radius and cross-checked (the
    (r,delta) two-parameter discipline) — disagreement rejects the
    step. No silent trust.
  - a duck-typed fast path is honored: if desys has .A_series(x, eps, dps, M)
    returning Taylor coefficients [A_0..A_M] of A about x, it is used instead
    (hook for exact-rational loaders / python-flint acceleration).

python-flint is an OPTIONAL accel (behind try/import); pure mpmath otherwise.
No sympy anywhere.
"""

import math
from fractions import Fraction

from mpmath import mp, mpf, mpc

try:  # optional acceleration hook — NEVER required
    import flint as _flint  # noqa: F401
    HAVE_FLINT = True
except Exception:
    HAVE_FLINT = False

try:  # acb step-kernel backend (opt-in) —
    # OPT-IN via transport_fixed_eps(..., backend="acb"); import-guarded so
    # the default mpmath path never needs flint.
    from . import acbfast as _acbfast
except ImportError:  # flat mode (tests import transport directly)
    try:
        import acbfast as _acbfast
    except Exception:
        _acbfast = None
except Exception:
    _acbfast = None

__all__ = ["transport_fixed_eps", "HAVE_FLINT"]

# hard caps
_MAX_STEPS_PER_LEG = 400
# Assert when halving drives the step FRACTION below 1/64
# (`assert hfrac > 1/64`). NOTE this deliberately reads the contract's
# "assert if h/|x1-x0| < 1/64" as a fraction-of-distance-proxy bound, not an
# absolute-|h| bound: near a singular approach (Frobenius landing feed) |h|
# scales with dist-to-singularity and legitimately drops below |x1-x0|/64
# while hfrac stays healthy. The convergence-driven hfrac asymptote is
# 10^-((dps+g)/MTAY) >= ~1/32 for the chassis MTAY default, so hfrac < 1/64
# genuinely signals a near-path singularity, same as in ode_fixed_eps.py.
_MIN_HFRAC = mpf(1) / 64
# Zeno abort when the sing-limited |h| drops below 2^-16 x the REMAINING
# distance to the leg end (a leg-length-relative
# threshold spuriously aborted legitimate approaches on long legs, e.g. a
# frobenius.land() launched far from x_sing with a second declared
# singularity nearby capping the match radius).
_ZENO_FRAC = mpf(2) ** -16


def _to_mpc(v):
    """Parse v (str/int/float/complex/Fraction/mpf/mpc/(re,im)) to mpc at the
    CURRENT working precision. Call only inside `with mp.workdps(...)` —
    import-dps footgun discipline."""
    if isinstance(v, (list, tuple)) and len(v) == 2:
        return mpc(_to_mpc(v[0]).real, _to_mpc(v[1]).real)
    if isinstance(v, Fraction):
        return mpc(mpf(v.numerator) / mpf(v.denominator))
    if isinstance(v, str):
        return mpc(mp.mpmathify(v))
    return mpc(v)


def _round_vec(vec, dps):
    """Round a list of mpc to dps (unary + re-rounds at the context prec)."""
    with mp.workdps(dps):
        return [+mpc(v) for v in vec]


def _fft(vec, wroots):
    """In-place-ish radix-2 DIT FFT with supplied twiddle table.
    wroots[j] = w^j for j in 0..N-1 where w = exp(sign*2*pi*i/N); the caller
    picks the sign via the table. len(vec) must be a power of two."""
    n = len(vec)
    a = list(vec)
    j = 0
    for i in range(1, n):
        bit = n >> 1
        while j & bit:
            j ^= bit
            bit >>= 1
        j |= bit
        if i < j:
            a[i], a[j] = a[j], a[i]
    length = 2
    while length <= n:
        step = n // length
        half = length // 2
        for i in range(0, n, length):
            for k in range(half):
                w = wroots[k * step]
                u = a[i + k]
                v = a[i + k + half] * w
                a[i + k] = u + v
                a[i + k + half] = u - v
        length <<= 1
    return a


def _next_pow2(m):
    n = 1
    while n < m:
        n <<= 1
    return n


class _CircleTrouble(Exception):
    """Sampling circle hit a non-finite / wildly inconsistent A — halve h."""


def _tail_bound(t_prev, t_last):
    """FINAL-ORDER-ONLY geometric fallback bound from the last two RAW term
    magnitudes: (t_last + t_prev) / (1 - r), r = min(3/4, t_last/t_prev)
    (r = 1/2 when t_prev == 0). This is the fixed-top-order step rule
    (the tail is bounded ONLY at the fixed top order M with exactly
    this formula). It is UNSOUND as an early-break criterion on sparse
    Taylor lattices — _tail_bound(0, 0) = 0 —
    so the sliding early break uses _tail_bound_nz over the
    last two ABOVE-FLOOR terms, and this fallback is reached only at the
    fixed top order mtay when fewer than two above-floor terms exist or the
    trailing zero-run exceeds the observed lattice gap (series genuinely
    terminated — A == 0 fixed point, polynomial solutions — or decayed
    below the roundoff floor)."""
    if t_prev == 0 and t_last == 0:
        return mpf(0)
    if t_prev > 0:
        r_meas = min(mpf(3) / 4, t_last / t_prev)
    else:
        r_meas = mpf(1) / 2
    return (t_last + t_prev) / (1 - r_meas)


def _tail_bound_nz(t_a, t_b, gap):
    """Zero-run-aware geometric tail bound from the last two ABOVE-FLOOR
    Taylor term magnitudes t_a (index a) and t_b (index b), gap = b - a >= 1
    The per-order decay is estimated
    gap-NORMALIZED, r = (t_b/t_a)^(1/gap), clamped at 3/4 (same soundness
    argument as the dense rule: |h| <= 0.5 x dist-to-singularity, or the
    two-radius cross-check when sings are unknown, keeps the true asymptotic
    per-order ratio <= 1/2 < 3/4), and the tail from order b is bounded by
    (t_a + t_b) / (1 - r) — for gap == 1 this reduces EXACTLY to the dense
    two-term rule with t_prev > 0. Zero/roundoff lattice terms never enter:
    on local support gaps >= 2 (y' = x^k y from x = 0, k >= 2) the raw
    last-two-terms rule would see (0, 0) and bank a false zero bound."""
    r_meas = (t_b / t_a) ** (mpf(1) / gap)
    if r_meas > mpf(3) / 4:
        r_meas = mpf(3) / 4
    return (t_a + t_b) / (1 - r_meas)


_TWIDDLE_CACHE = {}


def _circle_tables(N, wp):
    """(nodes, wroots_inv) with nodes[k] = e^{+2 pi i k/N},
    wroots_inv[j] = e^{-2 pi i j/N}, cached per (N, wp)."""
    key = (N, wp)
    if key not in _TWIDDLE_CACHE:
        with mp.workdps(wp):
            inv = [mp.expjpi(mpf(-2) * j / N) for j in range(N)]
            nodes = [mp.conj(w) for w in inv]
        if len(_TWIDDLE_CACHE) > 8:
            _TWIDDLE_CACHE.clear()
        _TWIDDLE_CACHE[key] = (nodes, inv)
    return _TWIDDLE_CACHE[key]


def _sample_circle(desys, center, eps_v, r, N, wp):
    """A(center + r*w^k, eps) for k=0..N-1 at working dps wp.
    Returns samples[k] = n x n list. Raises _CircleTrouble on non-finite."""
    n = desys.n
    nodes, _ = _circle_tables(N, wp)
    samples = []
    for k in range(N):
        z = center + r * nodes[k]
        try:
            Ak = desys.A(z, eps_v, wp)
        except (ZeroDivisionError, ValueError, ArithmeticError) as exc:
            raise _CircleTrouble(f"A eval failed on circle at {mp.nstr(z, 8)}: {exc}")
        for i in range(n):
            for j in range(n):
                v = mpc(Ak[i][j])
                if not (mp.isfinite(v.real) and mp.isfinite(v.imag)):
                    raise _CircleTrouble(
                        f"non-finite A[{i}][{j}] on circle at {mp.nstr(z, 8)}")
        samples.append(Ak)
    return samples


def _taylor_from_circle(samples, r, M, wroots_inv):
    """Taylor coefficients c_m (m=0..M) of each entry from circle samples:
    c_m * r^m = DFT(samples)[m]/N. Returns dict (i,j) -> [c_0..c_M], skipping
    identically-zero entries. wroots_inv[j] = exp(-2*pi*i*j/N)."""
    N = len(samples)
    n = len(samples[0])
    rinv = 1 / r
    out = {}
    for i in range(n):
        for j in range(n):
            col = [mpc(samples[k][i][j]) for k in range(N)]
            if all(v == 0 for v in col):
                continue
            hat = _fft(col, wroots_inv)
            coeffs = []
            rp = mpf(1)
            for m in range(M + 1):
                coeffs.append(hat[m] / N * rp)
                rp *= rinv
            out[(i, j)] = coeffs
    return out


def _local_a_series(desys, z0, eps_v, r, M, wp, crosscheck):
    """Taylor coefficients of A about z0 to order M.

    Fast path: desys.A_series(z0, eps, dps, M) if the loader provides it.
    Generic path: Cauchy circle of radius r with N >= max(2M+16, 3.5*wp)
    nodes; if crosscheck, a second circle at r/2 must agree on the low-order
    coefficients (the (r,delta) discipline) — else _CircleTrouble.
    """
    if hasattr(desys, "A_series"):
        ser = desys.A_series(z0, eps_v, wp, M)
        if len(ser) < M + 1:
            raise ValueError(
                f"A_series returned {len(ser)} order(s) < requested M+1 = "
                f"{M + 1} — refusing to zero-pad: a truncated fast path would "
                f"corrupt every step silently. A_series(x, eps, dps, M) must "
                f"return the Taylor coefficients [A_0..A_M] of A about x.")
        out = {}
        n = desys.n
        for i in range(n):
            for j in range(n):
                col = [mpc(ser[m][i][j]) for m in range(M + 1)]
                if any(v != 0 for v in col):
                    out[(i, j)] = col
        return out
    N = _next_pow2(max(2 * M + 16, int(3.5 * wp) + 16))
    _, wroots_inv = _circle_tables(N, wp)
    samples = _sample_circle(desys, z0, eps_v, r, N, wp)
    ser = _taylor_from_circle(samples, r, M, wroots_inv)
    if crosscheck:
        samples2 = _sample_circle(desys, z0, eps_v, r / 2, N, wp)
        ser2 = _taylor_from_circle(samples2, r / 2, min(8, M), wroots_inv)
        scale = max([abs(c[0]) for c in ser.values()], default=mpf(0))
        scale = max(scale, mpf(1))
        tol = mpf(10) ** (-(wp - 12))
        for key, c1 in ser.items():
            c2 = ser2.get(key, [mpc(0)] * 9)
            for m in range(min(8, M) + 1):
                # compare contributions at the sampling scale r/2
                if abs(c1[m] - c2[m]) * (r / 2) ** m > tol * scale:
                    raise _CircleTrouble(
                        f"two-radius Taylor cross-check failed on entry {key} "
                        f"order {m} (unknown singularity inside circle?)")
    return ser


def _dist_to_sings(z, sings):
    """min |z - s| over declared singular points; None if none declared."""
    if not sings:
        return None
    return min(abs(z - s) for s in sings)


def _declared_sings(desys, wp):
    """Singular points from desys.meta['singular_points'] or
    desys.singular_points (optional, duck-typed), parsed at wp."""
    raw = None
    meta = getattr(desys, "meta", None)
    if isinstance(meta, dict):
        raw = meta.get("singular_points")
    if raw is None:
        raw = getattr(desys, "singular_points", None)
    if raw is None:
        return []
    return [_to_mpc(s) for s in raw]


def _assert_landing_regular(desys, x1, eps_v, sings, dps):
    """Chassis rule (ode_fixed_eps.py: 'entry singular at y=0 — system not
    eligible'): refuse to march toward a singular landing point. For a
    singular endpoint use frobenius.land() instead (RESULT.md item 12)."""
    d = _dist_to_sings(x1, sings)
    if d is not None and d < mpf(10) ** (-max(10, dps // 2)):
        raise AssertionError(
            f"landing point {mp.nstr(x1, 10)} is a declared singular point of "
            f"A — fixed-eps Taylor transport cannot land there; use "
            f"frobenius.land() (Frobenius landing)")
    try:
        Aend = desys.A(x1, eps_v, min(dps, 30))
    except (ZeroDivisionError, ValueError, ArithmeticError) as exc:
        raise AssertionError(
            f"A(x1) evaluation failed at landing point {mp.nstr(x1, 10)}: {exc} "
            f"— landing point looks singular; use frobenius.land()")
    for row in Aend:
        for v in row:
            v = mpc(v)
            if not (mp.isfinite(v.real) and mp.isfinite(v.imag)):
                raise AssertionError(
                    f"A(x1) non-finite at landing point {mp.nstr(x1, 10)} — "
                    f"singular landing; use frobenius.land()")


def _assert_landing_regular_acb(desys, x1, sings, dps):
    """acb-backend landing check: same distance-to-declared-sings refusal
    as _assert_landing_regular, with the A(x1) finiteness probe evaluated
    through the acb entry table (taylor_entries at order 0 — the same
    fixed-eps rational functions the march uses; the generic exact .A()
    probe costs seconds per call at n~378 and is skipped on this backend)."""
    d = _dist_to_sings(x1, sings)
    if d is not None and d < mpf(10) ** (-max(10, dps // 2)):
        raise AssertionError(
            f"landing point {mp.nstr(x1, 10)} is a declared singular point of "
            f"A — fixed-eps Taylor transport cannot land there; use "
            f"frobenius.land() (Frobenius landing)")
    try:
        desys._acb_fast.taylor_entries(_acbfast.mpc_to_acb(mpc(x1)), 0)
    except _acbfast.AcbTrouble as exc:
        raise AssertionError(
            f"A(x1) evaluation failed at landing point {mp.nstr(x1, 10)}: "
            f"{exc} — landing point looks singular; use frobenius.land()")


def _march_leg(desys, eps_v, za, zb, y, sings, dps, mtay, ratio_frac,
               guard_extra, wp):
    """Taylor-march y from za to zb along the straight segment.
    Returns (y(zb), leg_diag) with leg_diag = {'steps', 'trunc_worst',
    'min_hfrac'} (honesty diagnostics — surfaced through transport_fixed_eps).
    Chassis: the fixed-eps main loop, vector specialization (the
    original marched the [U|V] matrix; a vector is the k=1 column case)."""
    z = mpc(za)
    leg_diag = {"steps": 0, "trunc_worst": mpf(0), "min_hfrac": mpf(ratio_frac)}
    leg_len = abs(zb - za)
    if leg_len == 0:
        return y, leg_diag
    if all(v == 0 for v in y):
        # exactly-zero state stays zero under a linear DE — nothing to march
        return list(y), leg_diag
    hfrac = mpf(ratio_frac)
    guard = mpf(10) ** (-(dps + guard_extra))
    # tiny floor: ONLY to keep the relative threshold nonzero if the state
    # ever underflows to exactly 0 mid-leg; far below any meaningful scale.
    tiny = mpf(10) ** (-(4 * wp))
    step = 0
    while True:
        rem = abs(zb - z)
        if rem == 0:
            break
        unit = (zb - z) / rem
        dsing = _dist_to_sings(z, sings)
        landing_step = False
        # distance proxy: declared dist-to-singularity when known, else the
        # leg length
        dprox = dsing if dsing is not None else leg_len
        hmag = hfrac * dprox
        if hmag >= rem:
            hmag = rem
            landing_step = True
        if (dsing is not None) and (not landing_step) \
                and hmag < rem * _ZENO_FRAC:
            raise RuntimeError(
                f"Zeno abort: sing-limited step |h|={mp.nstr(hmag, 4)} < "
                f"2^-16 x leg length — the path passes (essentially) through "
                f"a singular point. Fix: pass a complex detour `path` "
                f"(MUM-crossing pattern [x0/4, 0.03i, x1])")
        h = hmag * unit
        # sampling circle radius: <= 0.5 x dist-to-sing when known (alias
        # provably < 2^-N), else 2|h| + two-radius cross-check
        if dsing is not None:
            r_c = min(2 * hmag, dsing / 2)
        else:
            r_c = 2 * hmag
        try:
            Aser = _local_a_series(desys, z, eps_v, r_c, mtay, wp,
                                   crosscheck=(dsing is None))
        except _CircleTrouble:
            hfrac = hfrac / 2
            if hfrac < _MIN_HFRAC:
                raise AssertionError(
                    f"step {step}: A-series circle keeps failing even at "
                    f"step fraction < 1/64 — undeclared singularity on/near "
                    f"the path; use a complex detour `path` (MUM-crossing "
                    f"detour fix)")
            continue
        # Taylor recursion y_{m+1} = (1/(m+1)) sum_k A_k y_{m-k}
        n = desys.n
        ys = [list(y)]
        # BOTH convergence tests are RELATIVE to the
        # current state scale ||y||_inf. A chassis that marches
        # W0 = [I|0] with norm >= 1 has absolute == relative; that
        # invariant does NOT survive arbitrary y0 (eps-fan boundary seeds
        # ~ eps^k, eta^lambda-premultiplied Tradition seeds, Frobenius-landing
        # feeds are all legitimately small-norm). An absolute early-break
        # (mx < guard) cuts the series at order ~9 and silently loses
        # ~log10(1/||y0||) digits (measured: y0=1e-60 -> 12.8 of 65 digits);
        # an acceptance floor max(||y||, 1) has the same asymmetry.
        ynorm = max(abs(v) for v in y)
        thresh = guard * max(ynorm, tiny)
        # Acceptance is by the
        # ZERO-RUN-AWARE geometric tail bound over the last two ABOVE-FLOOR
        # term magnitudes (gap-normalized ratio, _tail_bound_nz), never by
        # raw consecutive terms that may sit on the zero/roundoff lattice.
        # Floor: term magnitudes below flr_rel x (running max term) are
        # roundoff-lattice zeros — Cauchy-circle noise sits ~10^-wp relative,
        # the 10-order headroom keeps genuine terms at the acceptance
        # threshold (10^-(dps+guard_extra) relative, 15+ orders above the
        # floor at the chassis wp = dps+guard_extra+25) unmistakable.
        flr_rel = mpf(10) ** (-(wp - 10))
        tmag = [ynorm]        # tmag[m] = max|y_m| * hmag^m
        tmax_run = ynorm
        nz_prev = None        # (index, mag) of the 2nd-last above-floor term
        nz_last = (0, ynorm)  # (index, mag) of the last above-floor term
        gap_max = 0           # largest observed gap between above-floor terms
        bound = None
        for m in range(mtay):
            nxt = [mpc(0)] * n
            for (i, j), coeffs in Aser.items():
                top = min(m, len(coeffs) - 1)
                for k in range(top + 1):
                    c = coeffs[k]
                    if c == 0:
                        continue
                    v = ys[m - k][j]
                    if v != 0:
                        nxt[i] += c * v
            inv = mpf(1) / (m + 1)
            nxt = [v * inv for v in nxt]
            ys.append(nxt)
            idx = m + 1
            tm = max(abs(v) for v in nxt) * hmag ** idx
            tmag.append(tm)
            if tm > tmax_run:
                tmax_run = tm
            if tm > flr_rel * tmax_run:
                g = idx - nz_last[0]
                if g > gap_max:
                    gap_max = g
                nz_prev = nz_last
                nz_last = (idx, tm)
            # Early break ONLY on the zero-run-aware bound, and ONLY while
            # the above-floor pair is FRESH (trailing zero-run since the last
            # above-floor term strictly shorter than the largest observed
            # inter-term gap; dense series: gap_max = 1, so the current term
            # itself must be above floor — the pre-F1 dense behavior). An
            # all-zero/roundoff trailing window NEVER accepts early — fall
            # through to the fixed top order mtay
            # (soundness > speed).
            if m > 8 and nz_prev is not None \
                    and idx - nz_last[0] < gap_max:
                b_early = _tail_bound_nz(nz_prev[1], nz_last[1],
                                         nz_last[0] - nz_prev[0])
                if b_early < thresh:
                    bound = b_early
                    break
        if bound is None:
            # ran to the fixed top order mtay
            if nz_prev is not None and mtay - nz_last[0] <= gap_max:
                # the above-floor pair is current through the top (every
                # lattice slot up to mtay was inspected): honest
                # zero-run-aware bound — also closes the top-order
                # lattice-misalignment hole the raw-top-two rule has
                bound = _tail_bound_nz(nz_prev[1], nz_last[1],
                                       nz_last[0] - nz_prev[0])
            else:
                # fewer than two above-floor terms (A ~ 0 / polynomial
                # solution: series genuinely terminated) or the trailing run
                # decayed below the roundoff floor beyond the observed
                # lattice gap: the fixed-top-order rule. The
                # t_prev == 0 -> r = 1/2 fallback lives ONLY here.
                bound = _tail_bound(tmag[-2], tmag[-1])
        if bound >= thresh:
            # insufficient convergence: halve the step fraction, retry
            hfrac = hfrac / 2
            if hfrac < _MIN_HFRAC:
                raise AssertionError(
                    f"step {step} cannot converge even at step fraction "
                    f"< 1/64 (geometric tail bound {mp.nstr(bound, 3)}) — "
                    f"singularity too close to the path; use a complex "
                    f"detour `path` (MUM-crossing detour fix)")
            continue
        # honesty diagnostics (per accepted step): bank the geometric BOUND,
        # not the raw last term (the bound is what the acceptance certified)
        trunc_rel = bound / max(ynorm, tiny)
        if trunc_rel > leg_diag["trunc_worst"]:
            leg_diag["trunc_worst"] = trunc_rel
        if hfrac < leg_diag["min_hfrac"]:
            leg_diag["min_hfrac"] = hfrac
        # evaluate at t = h (Horner over stored orders)
        ynew = [mpc(0)] * n
        hp = mpc(1)
        for coeffs in ys:
            for i in range(n):
                if coeffs[i] != 0:
                    ynew[i] += coeffs[i] * hp
            hp *= h
        y = ynew
        z = z + h
        hfrac = min(mpf(ratio_frac), hfrac * 2)
        step += 1
        leg_diag["steps"] = step
        if step > _MAX_STEPS_PER_LEG:
            raise RuntimeError(
                f"too many steps ({step}) on leg {mp.nstr(za, 8)} -> "
                f"{mp.nstr(zb, 8)} — likely a near-path singularity; "
                f"use a complex detour `path`")
    return y, leg_diag


_ACB_BLOCK = 8  # cross-block batch width for the acb convolution kernel
# (measured on one workstation, n=378/dps38: B=8 11.2 s/step vs B=4 11.8 s/step,
# results bit-identical; per-(entry,order) flint calls are ~us-class
# Python-side, so blocking cuts calls to ~nnz*mtay/B with the within-block
# remainder done as ~B/2 scalar acb ops per (entry,order))


def _march_leg_acb(desys, eps_v, za, zb, y, sings, dps, mtay, ratio_frac,
                   guard_extra, wp):
    """acb-backend mirror of _march_leg. SAME chassis semantics step for step — step control,
    Zeno abort, zero-run-aware geometric tail acceptance, recorded bound
    certificates — with the two hot kernels in flint acb balls:

      (1) local A-series: AcbFastTable.taylor_entries (grouped shifted-
          denominator series inversion + one short product per entry)
          instead of per-entry mpmath RF.taylor;
      (2) Taylor convolution y_{m+1} = (1/(m+1)) sum_k A_k y_{m-k}:
          contributions of completed B-order blocks of y are folded with
          one acb_poly product per (entry, block); the within-block
          triangular remainder runs as short scalar acb ops.

    Ball discipline ("M=50 + ftrim=acb.mid()
    per step"): mtay stays the chassis dps-scaled order, and the state is
    trimmed to ball midpoints once per ACCEPTED step; radii never
    accumulate across steps. Results are therefore NOT ball-certified —
    accuracy is certified exactly as on the mpmath backend (geometric
    tail bound + guard, recorded as trunc_worst). Magnitude bookkeeping for
    the acceptance logic is done in mpf via exact man_exp conversion
    (float would over/underflow at deep-ladder scales ~10^±1800).

    AcbTrouble from the series build (shifted-den ball contains 0) is
    treated exactly like _CircleTrouble: halve the step fraction."""
    table = desys._acb_fast
    if wp > table.wp:
        raise ValueError(
            f"acb backend: working digits wp={wp} exceed the table build "
            f"precision {table.wp} — re-run enable_fast_path(backend='acb') "
            f"with wp >= {wp} (silent under-precision refused)")
    z = mpc(za)
    leg_diag = {"steps": 0, "trunc_worst": mpf(0), "min_hfrac": mpf(ratio_frac)}
    leg_len = abs(zb - za)
    if leg_len == 0:
        return y, leg_diag
    if all(v == 0 for v in y):
        return list(y), leg_diag
    hfrac = mpf(ratio_frac)
    guard = mpf(10) ** (-(dps + guard_extra))
    tiny = mpf(10) ** (-(4 * wp))
    n = desys.n
    B = _ACB_BLOCK
    old_prec = _acbfast._fctx.prec
    # recursion precision is wp-scaled (mpmath-path parity: its recursion
    # rounds at the transport wp); the entry table itself computes at its
    # own build precision inside taylor_entries
    _acbfast._fctx.prec = _acbfast.prec_bits(wp)
    try:
        yv = [_acbfast.mpc_to_acb(v) for v in y]   # exact
        _zero = _acbfast._acb(0)
        kprec = _acbfast.prec_bits(wp)
        flr_rel = mpf(10) ** (-(wp - 10))
        step = 0
        # Per-expansion-point cache (h-INDEPENDENT work): the local Taylor
        # coefficients of A about z and the full series recursion depend
        # only on z and the state — NOT on the trial step h. On a halving
        # retry only the acceptance sweep (magnitude re-weighting by the
        # new hmag) needs recomputing, at mpf cost ~mtay. The recursion is
        # therefore computed ONCE per expansion point to the FULL chassis
        # order mtay (no early compute break; the acceptance logic runs on
        # the cached per-order norms with EXACTLY the chassis bookkeeping,
        # so accept/reject decisions and recorded bounds match a fresh run
        # at that h; the step evaluation then uses all computed orders —
        # extra orders sit below the accepted tail bound by construction).
        cache_z = None
        cache_ys = None
        cache_rawnorm = None
        cache_ynorm = None
        while True:
            rem = abs(zb - z)
            if rem == 0:
                break
            unit = (zb - z) / rem
            dsing = _dist_to_sings(z, sings)
            landing_step = False
            dprox = dsing if dsing is not None else leg_len
            hmag = hfrac * dprox
            if hmag >= rem:
                hmag = rem
                landing_step = True
            if (dsing is not None) and (not landing_step) \
                    and hmag < rem * _ZENO_FRAC:
                raise RuntimeError(
                    f"Zeno abort: sing-limited step |h|={mp.nstr(hmag, 4)} < "
                    f"2^-16 x leg length — the path passes (essentially) "
                    f"through a singular point. Fix: pass a complex detour "
                    f"`path` (MUM-crossing pattern [x0/4, 0.03i, x1])")
            h = hmag * unit
            if cache_z is None or cache_z != z:
                try:
                    Aser = table.taylor_entries(_acbfast.mpc_to_acb(z), mtay,
                                                nlow=max(8, B), prec=kprec)
                except _acbfast.AcbTrouble:
                    # z-level failure (expansion point on/too near a pole
                    # ball): halving cannot fix it, but the bounded-halving
                    # exit gives the same terminal detour message as the
                    # generic path — no silent hang.
                    hfrac = hfrac / 2
                    if hfrac < _MIN_HFRAC:
                        raise AssertionError(
                            f"step {step}: acb A-series build keeps failing "
                            f"at step fraction < 1/64 — undeclared "
                            f"singularity on/near the path; use a complex "
                            f"detour `path`")
                    continue
                items = list(Aser.items())
                # ---- Taylor recursion in acb (block-batched convolution,
                # computed once per z to full mtay) ----
                ys = [list(yv)]
                ynorm = max(_acbfast.acb_mag_mpf(v) for v in yv)
                rawnorm = [ynorm]
                acc = [None] * n       # cross-block accumulator polys
                for m in range(mtay):
                    bstart = (m // B) * B
                    if m == bstart and m > 0:
                        # fold the completed y block [bstart-B, bstart);
                        # pre-truncate to the live orders <= mtay
                        bs = bstart - B
                        ypolys = {}
                        for jj in range(n):
                            blk = [ys[bs + t][jj] for t in range(B)]
                            if any(not v.is_zero() for v in blk):
                                ypolys[jj] = _acbfast._acb_poly(blk)
                        for (i, j), (_low, poly) in items:
                            yp = ypolys.get(j)
                            if yp is None:
                                continue
                            pt = poly if bs == 0 \
                                else poly.truncate(mtay + 1 - bs)
                            contrib = (pt * yp).left_shift(bs)
                            acc[i] = contrib if acc[i] is None \
                                else acc[i] + contrib
                    s = [_zero] * n
                    for (i, j), (low, _poly) in items:
                        t = None
                        for k in range(m - bstart + 1):
                            v = ys[m - k][j]
                            if v.is_zero():
                                continue
                            if k >= len(low):
                                break
                            t = low[k] * v if t is None else t + low[k] * v
                        if t is not None:
                            s[i] = s[i] + t
                    inv = _acbfast._acb(1) / _acbfast._acb(m + 1)
                    nxt = [None] * n
                    for i in range(n):
                        a = acc[i]
                        si = s[i] if a is None else s[i] + a[m]
                        nxt[i] = si * inv
                    ys.append(nxt)
                    # per-order state norm (mpf, exponent-safe)
                    mx_i, mx_a = 0, None
                    for i in range(n):
                        ai = abs(nxt[i])
                        if mx_a is None or ai > mx_a:
                            mx_a, mx_i = ai, i
                    rawnorm.append(_acbfast.acb_mag_mpf(nxt[mx_i]))
                cache_z = mpc(z)
                cache_ys = ys
                cache_rawnorm = rawnorm
                cache_ynorm = ynorm
            else:
                ys = cache_ys
                rawnorm = cache_rawnorm
                ynorm = cache_ynorm
            # ---- acceptance sweep at THIS h (chassis bookkeeping verbatim
            # on the cached per-order norms; tm = rawnorm[m]*hmag^m is
            # exactly what a fresh run computes) ----
            thresh = guard * max(ynorm, tiny)
            tmag = [ynorm]
            tmax_run = ynorm
            nz_prev = None
            nz_last = (0, ynorm)
            gap_max = 0
            bound = None
            for m in range(mtay):
                idx = m + 1
                tm = rawnorm[idx] * hmag ** idx
                tmag.append(tm)
                if tm > tmax_run:
                    tmax_run = tm
                if tm > flr_rel * tmax_run:
                    g = idx - nz_last[0]
                    if g > gap_max:
                        gap_max = g
                    nz_prev = nz_last
                    nz_last = (idx, tm)
                if m > 8 and nz_prev is not None \
                        and idx - nz_last[0] < gap_max:
                    b_early = _tail_bound_nz(nz_prev[1], nz_last[1],
                                             nz_last[0] - nz_prev[0])
                    if b_early < thresh:
                        bound = b_early
                        break
            if bound is None:
                if nz_prev is not None and mtay - nz_last[0] <= gap_max:
                    bound = _tail_bound_nz(nz_prev[1], nz_last[1],
                                           nz_last[0] - nz_prev[0])
                else:
                    bound = _tail_bound(tmag[-2], tmag[-1])
            if bound >= thresh:
                hfrac = hfrac / 2
                if hfrac < _MIN_HFRAC:
                    raise AssertionError(
                        f"step {step} cannot converge even at step fraction "
                        f"< 1/64 (geometric tail bound {mp.nstr(bound, 3)}) "
                        f"— singularity too close to the path; use a complex "
                        f"detour `path` (MUM-crossing detour fix)")
                continue
            trunc_rel = bound / max(ynorm, tiny)
            if trunc_rel > leg_diag["trunc_worst"]:
                leg_diag["trunc_worst"] = trunc_rel
            if hfrac < leg_diag["min_hfrac"]:
                leg_diag["min_hfrac"] = hfrac
            # evaluate at t = h; trim to ball midpoints (the ftrim rule)
            h_acb = _acbfast.mpc_to_acb(h)
            ynew = list(ys[0])
            hp = h_acb
            for coeffs in ys[1:]:
                for i in range(n):
                    ci = coeffs[i]
                    if not ci.is_zero():
                        ynew[i] = ynew[i] + ci * hp
                hp = hp * h_acb
            yv = [v.mid() for v in ynew]
            z = z + h
            cache_z = None      # state advanced: cached recursion is stale
            cache_ys = cache_rawnorm = cache_ynorm = None
            hfrac = min(mpf(ratio_frac), hfrac * 2)
            step += 1
            leg_diag["steps"] = step
            if step > _MAX_STEPS_PER_LEG:
                raise RuntimeError(
                    f"too many steps ({step}) on leg {mp.nstr(za, 8)} -> "
                    f"{mp.nstr(zb, 8)} — likely a near-path singularity; "
                    f"use a complex detour `path`")
        y_out = [+_acbfast.acb_to_mpc(v) for v in yv]
    finally:
        _acbfast._fctx.prec = old_prec
    return y_out, leg_diag


def transport_fixed_eps(desys, eps0, x0, x1, y0, dps, path=None, mtay=None,
                        ratio=0.5, guard_extra=8, return_diag=False,
                        backend=None):
    """Transport y' = A(x, eps0) y from x0 to x1; return y(x1) as list[mpc].

    Step control: h = ratio x dist-to-singularity,
    halve on non-convergence, exact final landing; MTAY default
    max(60, 0.75*dps+25); landing regularity assert.
    The convergence ACCEPTANCE rule is the ZERO-RUN-AWARE
    geometric tail bound
    (t_a + t_b)/(1 - min(3/4, (t_b/t_a)^(1/(b-a)))) over the last two
    ABOVE-FLOOR Taylor term magnitudes (gap-normalized ratio; early break
    only while that pair is fresh — an all-zero/roundoff trailing window
    falls through to the fixed top order), with the
    complex-safe pole clearance min |z - s| over the declared singular
    points (real-pole-only clearance / single-term acceptance / raw
    last-two-terms acceptance on sparse lattices are catalogued
    silent-wrong footguns).

    Args:
        desys: wayfinder DESystem (needs .n, .A(x, eps, dps); optional
            .meta['singular_points'] / .singular_points list, optional
            .A_series fast path).
        eps0: fixed epsilon (str/mpf/mpc/Fraction/...). Complex eps0 is
            allowed HERE (the transport itself is analytic in eps); the
            real-eps0-only policy applies to Frobenius matching (frobenius.py).
        x0, x1: endpoints (parsed at working precision).
        y0: initial vector at x0, length desys.n.
        dps: target decimal precision of the result.
        path: optional list of intermediate complex waypoints between x0 and
            x1 (detours around singularities — the MUM-crossing detour pattern
            [x0/4, 0.03i, x1]). x0/x1 need not be repeated in it (a repeated
            endpoint produces a zero-length leg, which is skipped).
        mtay: Taylor order cap (default max(60, int(0.75*dps)+25)).
        ratio: initial/maximal step fraction of the distance-to-singularity
            estimate (default 0.5, chassis value).
        guard_extra: convergence guard is tail < 10^-(dps+guard_extra)
            RELATIVE to the state norm ||y||_inf (an absolute
            guard is correct only for a
            norm >= 1 [I|0] state; small-norm y0 would silently lose digits).
        return_diag: if True, return (y, diag) with the chassis honesty
            diagnostics diag = {'steps': int (accepted steps, all legs),
            'trunc_worst': str (worst per-step relative geometric tail
            BOUND — what the acceptance certified, not the raw last term),
            'trunc_worst_log10': float or None,
            'min_hfrac': float (smallest accepted step fraction)} —
            callers should record these alongside results.

    Raises:
        AssertionError: singular landing point, or step halved below
            |x1-x0|/64 without convergence (use a detour path).
        RuntimeError: Zeno abort / step-count cap (use a detour path).

    backend: None/"mpmath" (default — byte-identical legacy behavior) or
        "acb" (OPT-IN flint acb step kernel; requires desys.enable_fast_path(eps, wp,
        backend="acb") and python-flint; same chassis semantics and
        certificates, state ball-mid-trimmed per accepted step — see
        _march_leg_acb / acbfast.py).
    """
    if backend not in (None, "mpmath", "acb"):
        raise ValueError(f"backend must be None|'mpmath'|'acb', "
                         f"got {backend!r}")
    use_acb = (backend == "acb")
    if use_acb:
        if _acbfast is None or not _acbfast.HAVE_FLINT:
            raise ImportError("backend='acb' requested but python-flint is "
                              "not importable")
        if getattr(desys, "_acb_fast", None) is None:
            raise ValueError(
                "backend='acb' requires the acb fast path — call "
                "desys.enable_fast_path(eps, wp, backend='acb') first")
        chk = getattr(desys, "_acb_check_eps", None)
        if chk is not None:
            chk(eps0, dps)   # fixed-eps layer: silent wrong-eps refused
    if mtay is None:
        mtay = max(60, int(0.75 * dps) + 25)
    wp = dps + guard_extra + 25
    with mp.workdps(wp):
        eps_v = _to_mpc(eps0)
        za = _to_mpc(x0)
        zb = _to_mpc(x1)
        y = [_to_mpc(v) for v in y0]
        if len(y) != desys.n:
            raise ValueError(f"y0 length {len(y)} != system size {desys.n}")
        sings = _declared_sings(desys, wp)
        if use_acb:
            _assert_landing_regular_acb(desys, zb, sings, dps)
        else:
            _assert_landing_regular(desys, zb, eps_v, sings, dps)
        waypoints = [za]
        for w in (path or []):
            waypoints.append(_to_mpc(w))
        waypoints.append(zb)
        diag_acc = {"steps": 0, "trunc_worst": mpf(0),
                    "min_hfrac": mpf(ratio)}
        marcher = _march_leg_acb if use_acb else _march_leg
        for k in range(len(waypoints) - 1):
            y, legd = marcher(desys, eps_v, waypoints[k], waypoints[k + 1],
                              y, sings, dps, mtay, ratio, guard_extra, wp)
            diag_acc["steps"] += legd["steps"]
            if legd["trunc_worst"] > diag_acc["trunc_worst"]:
                diag_acc["trunc_worst"] = legd["trunc_worst"]
            if legd["min_hfrac"] < diag_acc["min_hfrac"]:
                diag_acc["min_hfrac"] = legd["min_hfrac"]
        result = list(y)
        if return_diag:
            tw = diag_acc["trunc_worst"]
            diag_out = {
                "steps": int(diag_acc["steps"]),
                "trunc_worst": mp.nstr(tw, 6),
                "trunc_worst_log10": (float(mp.log10(tw)) if tw > 0 else None),
                "min_hfrac": float(diag_acc["min_hfrac"]),
            }
    if return_diag:
        return _round_vec(result, dps), diag_out
    return _round_vec(result, dps)

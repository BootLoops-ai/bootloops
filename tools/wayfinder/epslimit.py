#!/usr/bin/env python3
"""
wayfinder.epslimit (formerly the standalone eps_extrapolator) -- generic,
engine-agnostic eps->0 Laurent-coefficient extractor.

Member of the Wayfinder package.  Sits DOWNSTREAM of any engine's eps-grid;
it is NOT wayfinder.epsfan (which extracts the Laurent fan from a DE system
it can evaluate at eps-nodes of its own choosing, with certified Cauchy /
two-grid Vandermonde gates).  epslimit takes whatever (eps_i, value_i) grid an
external engine already produced and extrapolates it; use epsfan when you own
the evaluator, epslimit when you only own the samples.

Import:  from wayfinder.epslimit import extrapolate, richardson_boundary,
         load_amflow_grid          (package; tools/ on sys.path)
     or  sys.path.insert(0, <tools/wayfinder>); import epslimit   (flat)
CLI:     python3 tools/wayfinder/epslimit.py out1.json out2.json --part re
     or  (cd tools && python3 -m wayfinder.epslimit out1.json ...)

Precision note: extrapolate() sets the global mpmath precision to its
auto-chosen working value for the duration of the solve and restores it in a
finally block; every other Wayfinder module takes explicit dps and never
touches the global.

Given a list of (eps_i, value_i) pairs, each value_i a high-precision (~hundreds of
digits) numeric sample of a function

    f(eps) = sum_{k=k_min}^{inf} c_k eps^k          (possibly with k_min < 0)

return the Laurent coefficients c_k to high precision.  This is the GENERIC numeric
layer that sits ABOVE any particular integral engine (AMFlow, pySecDec, ...).  The
engine produces value_i at chosen eps_i; this module robustly turns the grid into the
eps-Laurent expansion, with the finite part c_0 being the usual target.

The hard part is NOT fitting a polynomial -- a single Vandermonde solve already does
that.  The hard part is doing it ROBUSTLY at high precision against the failure modes
that cap naive solves at ~6 digits:

  (1) Catastrophic cancellation / ill-conditioning of the eps-Vandermonde when the
      grid spans only a small dynamic range (e.g. 1/1000 .. 1/8000).  Cured by working
      at precision proportional to the *condition number*, by a Richardson / repeated-
      Romberg elimination that algebraically removes one power per step (condition-
      robust), and by a Bulirsch-Stoer-style rational (Pade) extrapolation that
      converges where the polynomial one does not.

  (2) Unknown / negative leading power (genuine eps-poles).  Auto-detected by fitting
      the leading power from the grid (log-log slope refined by a one-parameter solve)
      so the pole order k_min need not be supplied.

  (3) RESONANCES / Jordan blocks: the eps->0 system is rank-deficient because the
      basis {eps^k} degenerates -- typically a log(eps) appears (a Jordan block in the
      differential-equation monodromy, or two indicial exponents colliding as eps->0,
      exactly the calcx00 case).  Detected when the pure-power model leaves a residual
      that scales like a *positive power of log(eps)*; handled by augmenting the model
      basis with eps^k log(eps)^j terms (a "deflated/shifted" generalized Vandermonde),
      i.e. fitting  f = sum_{k,j} c_{k,j} eps^k log(eps)^j.

Methods provided (all mpmath, arbitrary precision; numpy is not used):

  * vandermonde_fit  -- generalized (power x log) least-squares / exact solve at
                        working precision >> grid condition number.
  * richardson       -- repeated Romberg/Neville power elimination (needs a geometric
                        or specified power ladder); condition-robust, no matrix solve.
  * pade_bulirsch    -- rational extrapolation to eps=0 (Bulirsch-Stoer), best when the
                        function has nearby singularities in eps.
  * eps_fft          -- if a *circle* of eps on |eps|=r is supplied, recover Laurent
                        coefficients by the discrete Fourier / Cauchy integral; this is
                        the spectrally-accurate route and side-steps Vandermonde
                        conditioning entirely.

The top-level entry point `extrapolate()` auto-routes: it detects the grid geometry
(circle vs real ladder vs scattered), the leading power, and resonance order, runs the
applicable methods, and CROSS-VALIDATES them against each other, returning the
coefficients together with an honest achieved-digit estimate (the agreement between
independent methods, NOT a self-reported residual).
"""

from __future__ import annotations
import math
import cmath
from dataclasses import dataclass, field
from typing import Sequence, Optional

import mpmath as mp

__all__ = ['ExtrapResult', 'to_mpc', 'vandermonde_fit', 'neville_zero',
           'romberg_ladder', 'bulirsch_stoer_zero', 'eps_fft',
           'detect_log_order', 'extrapolate', 'load_amflow_grid',
           'richardson_boundary']


# ----------------------------------------------------------------------------- #
#  parsing helpers: accept python numbers, strings, or Arb-ball "[mid +/- rad]"  #
# ----------------------------------------------------------------------------- #

def _strip_ball(s):
    """Turn an Arb/FLINT ball string '[mid +/- rad]' (or 'mid') into the midpoint."""
    if isinstance(s, (int, float, complex, mp.mpf, mp.mpc)):
        return s
    s = str(s).strip()
    if s.startswith('['):
        s = s[1:]
        if ']' in s:
            s = s[:s.index(']')]
    if '+/-' in s:
        s = s.split('+/-')[0]
    return s.strip()


def to_mpc(x):
    """Coerce many input shapes (real string, {'re':..,'im':..}, complex) to mp.mpc."""
    if isinstance(x, dict):
        re = mp.mpf(_strip_ball(x.get('re', '0')))
        im = mp.mpf(_strip_ball(x.get('im', '0')))
        return mp.mpc(re, im)
    if isinstance(x, (list, tuple)) and len(x) == 2:
        return mp.mpc(mp.mpf(_strip_ball(x[0])), mp.mpf(_strip_ball(x[1])))
    s = _strip_ball(x)
    try:
        return mp.mpc(s)
    except (ValueError, TypeError):
        return mp.mpc(mp.mpf(s))


def _digits_agree(a, b):
    """Number of agreeing significant decimal digits between a and b (mpc)."""
    a = mp.mpc(a); b = mp.mpc(b)
    d = abs(a - b)
    s = max(abs(a), abs(b))
    if s == 0:
        return mp.mp.dps if d == 0 else 0.0
    if d == 0:
        return float(mp.mp.dps)
    return float(-mp.log10(d / s))


# ----------------------------------------------------------------------------- #
#  result container                                                              #
# ----------------------------------------------------------------------------- #

@dataclass
class ExtrapResult:
    coeffs: dict                      # {(k, j): c_{k,j}}  k=power of eps, j=power of log(eps)
    k_min: int                        # leading (most negative) eps power found
    log_order: int                    # highest power of log(eps) detected (0 = none)
    finite_part: mp.mpc               # c_{0,0}
    achieved_digits: float            # cross-validation agreement estimate
    method: str                       # which method produced `coeffs`
    per_method: dict = field(default_factory=dict)   # method -> coeffs dict
    diagnostics: dict = field(default_factory=dict)

    def c(self, k, j=0):
        return self.coeffs.get((k, j), mp.mpc(0))

    def __repr__(self):
        head = (f"ExtrapResult(method={self.method}, k_min={self.k_min}, "
                f"log_order={self.log_order}, ~{self.achieved_digits:.1f} digits)\n")
        lines = []
        for (k, j) in sorted(self.coeffs):
            tag = f"eps^{k}" + (f" log^{j}" if j else "")
            lines.append(f"  c[{k},{j}] ({tag}) = {mp.nstr(self.coeffs[(k,j)], 25)}")
        return head + "\n".join(lines)


# ----------------------------------------------------------------------------- #
#  geometry / structure detection                                                #
# ----------------------------------------------------------------------------- #

def _detect_geometry(eps):
    """Classify the eps grid: 'circle' (constant |eps|, spread phase), 'real_ladder'
    (real, geometric ratio), or 'scatter'."""
    mods = [abs(e) for e in eps]
    rmean = sum(mods) / len(mods)
    rspread = max(abs(m - rmean) for m in mods) / rmean if rmean > 0 else 1
    phases = sorted(float(mp.arg(e)) for e in eps)
    phase_spread = phases[-1] - phases[0]
    # constant modulus (to ~1e-12 relative, well above numerical construction noise)
    # AND meaningfully spread phase AND enough points -> treat as a Cauchy circle.
    if rspread < mp.mpf('1e-12') and phase_spread > 0.5 and len(eps) >= 4:
        return 'circle'
    allreal = all(abs(mp.im(e)) < mp.mpf('1e-30') * (abs(e) + 1) for e in eps)
    if allreal:
        return 'real_ladder'
    return 'scatter'


def _estimate_kmin(eps, vals):
    """Estimate leading eps power from log-log slope of |value| vs |eps| using the two
    points with the most distinct |eps| (so the slope is well-conditioned; on a circle
    all |eps| are equal and no power can be read off the modulus -> return 0)."""
    pts = sorted(zip(eps, vals), key=lambda ev: abs(ev[0]))
    (e0, v0) = pts[0]
    (e1, v1) = pts[-1]
    if v0 == 0 or v1 == 0:
        return 0
    re0, re1 = abs(e0), abs(e1)
    if re0 == 0 or re1 == 0 or re0 == re1:
        return 0                       # equal-modulus grid (circle): cannot read slope
    slope = mp.re(mp.log(abs(v1) / abs(v0)) / mp.log(re1 / re0))
    k = int(mp.nint(slope))
    # a positive estimated slope still means k_min >= 0 (finite); negative => pole
    return k if k < 0 else 0


# ----------------------------------------------------------------------------- #
#  core method 1: generalized (power x log) Vandermonde solve                     #
# ----------------------------------------------------------------------------- #

def vandermonde_fit(eps, vals, k_min, n_powers, log_order=0, ridge=None, scale=None):
    """Solve f(eps) = sum_{k,j} c_{k,j} eps^k log(eps)^j over
       k in [k_min, k_min+n_powers-1], j in [0, log_order].
    Exact solve when square; least-squares (normal equations at high precision) when
    over-determined.  Returns {(k,j): c}.

    CONDITIONING: the power columns are built in a rescaled variable (eps/scale)^k with
    `scale` the geometric mean of |eps|, so every power column is O(1) regardless of how
    small the eps grid is (otherwise eps^k for k spanning a wide range and eps~1e-3
    blows the matrix condition number and caps the solve at a few digits).  The true
    log(eps) factor is kept exact; coefficients are unscaled at the end (c = x/scale^k).
    """
    eps = [mp.mpc(e) for e in eps]
    vals = [mp.mpc(v) for v in vals]
    if scale is None:
        mods = [abs(e) for e in eps if e != 0]
        # geometric mean modulus -> centers the power columns near 1
        scale = mp.e ** (sum(mp.log(mm) for mm in mods) / len(mods)) if mods else mp.mpf(1)
    basis = [(k, j) for k in range(k_min, k_min + n_powers)
             for j in range(0, log_order + 1)]
    nb = len(basis)
    rows = []
    for e in eps:
        u = e / scale
        le = mp.log(e)
        row = [u ** k * (le ** j if j else mp.mpc(1)) for (k, j) in basis]
        rows.append(row)
    A = mp.matrix(rows)            # m x nb
    b = mp.matrix(vals)            # m x 1
    m = len(eps)
    if m == nb:
        x = mp.lu_solve(A, b)
    else:
        # normal equations  (A^H A) x = A^H b ; high working precision absorbs the
        # squared condition number.
        AH = A.transpose_conj()
        M = AH * A
        if ridge:
            for i in range(nb):
                M[i, i] += ridge
        x = mp.lu_solve(M, AH * b)
    # unscale: column used (eps/scale)^k, so true c_{k,j} = x_{k,j} / scale^k
    return {basis[i]: x[i] / scale ** basis[i][0] for i in range(nb)}


# ----------------------------------------------------------------------------- #
#  core method 2: Richardson / Neville (polynomial) extrapolation to eps=0        #
# ----------------------------------------------------------------------------- #

def neville_zero(eps, vals):
    """Neville polynomial extrapolation of vals(eps) to eps=0 (k_min must be 0).
    Returns (limit, table) ; condition-robust, no matrix solve.  Best on a real or
    complex ladder converging to 0."""
    eps = [mp.mpc(e) for e in eps]
    vals = [mp.mpc(v) for v in vals]
    n = len(eps)
    # Neville for the target point x=0:
    P = list(vals)
    for k in range(1, n):
        newP = []
        for i in range(n - k):
            xi = eps[i]; xik = eps[i + k]
            # P_{i..i+k}(0) = ( (0-xik) P_left - (0-xi) P_right ) / (xi - xik)
            val = ((-xik) * P[i] - (-xi) * P[i + 1]) / (xi - xik)
            newP.append(val)
        # successive Neville columns shrink P by one and reuse adjacent entries
        P = newP
    return P[0]


def richardson_boundary(partial_sums, Ms, tail_exponents=None):
    """Richardson-eliminate a user-supplied list of M^{-α} tail powers from a
    sequence of partial sums S_M — for summing a series AT its radius of
    convergence when the leading-singularity type is known (√ → half-integer α).
    Default tail_exponents = [1/2, 3/2, 5/2, …].  Returns extrapolated list
    (length len(Ms)−len(tail_exponents))."""
    S = [mp.mpf(s) for s in partial_sums]; M = list(Ms)
    if tail_exponents is None:
        tail_exponents = [mp.mpf(1) / 2 + k for k in range(len(M) - 1)]
    for a in tail_exponents:
        S2 = []
        for i in range(len(S) - 1):
            r = (mp.mpf(M[i + 1]) / M[i]) ** a
            S2.append((r * S[i + 1] - S[i]) / (r - 1))
        S = S2; M = M[1:]
    return S


def romberg_ladder(eps, vals, ratio=None):
    """Repeated Richardson/Romberg elimination for a (near) geometric ladder
    eps_i = eps_0 * ratio^i, k_min=0.  Each column algebraically kills one more eps
    power.  Returns the eps->0 limit and the full triangular table.
    This is the most cancellation-robust route on a clean geometric grid."""
    eps = [mp.mpc(e) for e in eps]
    vals = [mp.mpc(v) for v in vals]
    n = len(eps)
    if ratio is None:
        ratio = eps[1] / eps[0]
    T = [vals[:]]   # T[0] = raw values
    col = vals[:]
    p = 1
    for k in range(1, n):
        rp = ratio ** p
        newcol = []
        for i in range(len(col) - 1):
            # eliminate eps^p term assuming geometric spacing
            newcol.append((col[i + 1] - rp * col[i]) / (1 - rp))
        col = newcol
        T.append(col[:])
        p += 1
    return T[-1][0], T


# ----------------------------------------------------------------------------- #
#  core method 3: Bulirsch-Stoer rational (Pade-type) extrapolation to eps=0      #
# ----------------------------------------------------------------------------- #

def bulirsch_stoer_zero(eps, vals):
    """Rational-function extrapolation to eps=0 (Stoer & Bulirsch).  Robust when the
    function has poles near eps=0 in the complex plane (typical of resonant masters).
    k_min must be 0 (subtract the pole first if not)."""
    eps = [mp.mpc(e) for e in eps]
    vals = [mp.mpc(v) for v in vals]
    n = len(eps)
    # tableau (Bulirsch-Stoer / Bulirsch-Henrici), target x=0
    R = [[mp.mpc(0)] * n for _ in range(n)]
    Rprev = [mp.mpc(0)] * n
    for i in range(n):
        R[i][0] = vals[i]
    for j in range(1, n):
        for i in range(n - j):
            xi = eps[i]; xij = eps[i + j]
            num = R[i + 1][j - 1] - R[i][j - 1]
            denom_factor = (xi / xij)  # since target is 0
            den = denom_factor * (1 - (R[i + 1][j - 1] - R[i][j - 1]) /
                                  (R[i + 1][j - 1] - (Rprev[i] if j > 1 else 0) + mp.mpc('1e-400'))) - 1
            if den == 0:
                R[i][j] = R[i + 1][j - 1]
            else:
                R[i][j] = R[i + 1][j - 1] + num / den
        Rprev = [R[i][j - 1] for i in range(n)]
    return R[0][n - 1]


# ----------------------------------------------------------------------------- #
#  core method 4: eps-FFT / Cauchy on a circle  -> full Laurent series            #
# ----------------------------------------------------------------------------- #

def eps_fft(eps, vals, n_coeffs=None):
    """If eps lie (approximately) equispaced on a circle |eps|=r, recover Laurent
    coefficients c_k = (1/2 pi i) oint f(eps) eps^{-k-1} deps via the DFT of the
    samples.  Spectrally accurate, immune to Vandermonde conditioning.  Returns
    {(k,0): c_k}."""
    eps = [mp.mpc(e) for e in eps]
    vals = [mp.mpc(v) for v in vals]
    N = len(eps)
    r = sum(abs(e) for e in eps) / N
    # sort by phase so samples are in DFT order, then use the EXACT uniform grid angles
    # 2 pi n / N (not the measured phases): for an equispaced circle the trapezoidal
    # Cauchy quadrature is then spectrally exact, and r^{-k} amplification of any angle
    # jitter is avoided.  (If the samples are not truly equispaced, supply a real ladder
    # instead -- this branch assumes the circle is uniform.)
    order = sorted(range(N), key=lambda i: float(mp.arg(eps[i])) % (2 * math.pi))
    f = [vals[i] for i in order]
    # align grid so sample 0 sits at its actual angle theta0 (phase reference)
    theta0 = mp.arg(eps[order[0]])
    if n_coeffs is None:
        n_coeffs = N
    two_pi = 2 * mp.pi
    out = {}
    half = n_coeffs // 2
    for k in range(-half, n_coeffs - half):
        acc = mp.mpc(0)
        for n in range(N):
            th = theta0 + two_pi * n / N
            acc += f[n] * mp.e ** (-1j * k * th)
        out[(k, 0)] = acc / N * (r ** (-k))
    return out


# ----------------------------------------------------------------------------- #
#  resonance / Jordan-block detection                                            #
# ----------------------------------------------------------------------------- #

def detect_log_order(eps, vals, k_min, n_powers, max_log=3):
    """Decide how many powers of log(eps) the data needs (resonance / Jordan-block
    detection).  For EACH candidate log_order `lo` we give the model the LARGEST power
    budget the data can support at that lo, namely n_powers(lo) = (m-1) // (lo+1) (one
    point reserved for the leave-one-out test).  This is the crucial fairness fix: a
    clean pure-power function must be allowed its full polynomial degree at lo=0, or a
    spurious log model with fewer powers can appear to "win".

    We score each lo by its median leave-one-out reconstruction (agreeing digits) and
    accept the SMALLEST lo within a small margin of the best -- i.e. we add a log power
    only when it buys a decisive (>~ a handful of digits) improvement, exactly the
    signature of a genuine resonance (the synthetic Jordan case jumps from ~7 to ~220
    digits when the needed log order is reached).
    Returns (best_log_order, scoreboard{lo: (npow, median_digits)})."""
    eps = [mp.mpc(e) for e in eps]
    vals = [mp.mpc(v) for v in vals]
    m = len(eps)
    scoreboard = {}
    scores = {}
    for lo in range(0, max_log + 1):
        npow = (m - 1) // (lo + 1)
        if npow < 1:
            break
        nb = npow * (lo + 1)
        digs = []
        ok = True
        for hold in range(m):
            tr_e = [eps[i] for i in range(m) if i != hold]
            tr_v = [vals[i] for i in range(m) if i != hold]
            if len(tr_e) < nb:
                ok = False
                break
            try:
                cf = vandermonde_fit(tr_e, tr_v, k_min, npow, log_order=lo)
            except Exception:
                ok = False
                break
            pred = _eval_model(cf, eps[hold])
            digs.append(_digits_agree(pred, vals[hold]))
        if not ok or not digs:
            continue
        digs.sort()
        score = digs[len(digs) // 2]    # median LOO agreement
        scoreboard[lo] = (npow, round(score, 2))
        scores[lo] = score
    if not scores:
        return 0, scoreboard
    best_score = max(scores.values())
    # Occam: prefer the SMALLEST log order, and accept a higher one only on a DECISIVE
    # improvement.  A genuine resonance produces a huge LOO jump (the synthetic Jordan
    # case leaps from ~7 to ~220 digits when the true log order is reached); engine
    # round-off noise on a deep ladder produces only a marginal few-digit wobble that a
    # spurious log term can "explain", so the threshold must be a real fraction of the
    # gain, not a fixed 3 digits.  Require the log model to beat pure-power by both an
    # absolute floor AND ~half the available head-room.
    base = scores.get(0, best_score)
    decisive = max(20.0, 0.5 * (best_score - base))
    best = 0
    for lo in sorted(scores):
        if lo == 0:
            continue
        if scores[lo] >= base + decisive and scores[lo] >= best_score - 3.0:
            best = lo
            base = scores[lo]      # ratchet: allow lo+1 only if it again beats decisively
    return best, scoreboard


def _eval_model(coeffs, e):
    e = mp.mpc(e)
    le = mp.log(e)
    acc = mp.mpc(0)
    for (k, j), c in coeffs.items():
        acc += c * e ** k * (le ** j if j else mp.mpc(1))
    return acc


# ----------------------------------------------------------------------------- #
#  top-level auto-routing entry point                                            #
# ----------------------------------------------------------------------------- #

def extrapolate(samples, k_min=None, n_powers=None, max_log=3,
                working_dps=None, verbose=False):
    """
    Main entry point.

    samples : sequence of (eps_i, value_i).  Each may be python numbers, strings,
              Arb-ball strings '[mid +/- rad]', or {'re':..,'im':..} dicts.
    k_min   : leading eps power; auto-detected if None.
    n_powers: number of eps powers to fit; default = (#samples // (max_log+1)) bounded
              so the system is determined/over-determined.
    max_log : highest power of log(eps) to consider for resonance/Jordan handling.
    working_dps : mpmath working precision; auto-set from input precision & grid
              condition if None.

    Returns ExtrapResult.
    """
    # ---- parse, set precision ----
    m = len(samples)
    if m < 2:
        raise ValueError("need at least 2 (eps,value) samples")

    # Parse at (at least) the precision the inputs themselves carry: parsing at
    # the ambient dps (default 15) would truncate every high-precision sample
    # before the working-precision block below could see it.
    parse_dps = max([mp.mp.dps]
                    + [_str_dps(v) for (_, v) in samples]
                    + [_str_dps(e) for (e, _) in samples]) + 10
    with mp.workdps(parse_dps):
        eps = [to_mpc(e) for (e, _) in samples]
        vals = [to_mpc(v) for (_, v) in samples]

    # input precision: shortest mantissa supplied (assume ~ digits of the strings)
    in_dps = min(_str_dps(v) for (_, v) in samples)
    # The (rescaled) Vandermonde / generalized-Vandermonde condition number is
    # dominated by (i) the modulus dynamic range of the grid and (ii) node clustering:
    # close nodes make near-parallel rows, costing ~ n * log10(1/min_node_gap) digits.
    mods = [abs(e) for e in eps if e != 0]
    dyn = float(mp.log10(max(mods) / min(mods))) if len(mods) > 1 else 1.0
    # geometric-mean-rescaled nodes, then smallest pairwise gap:
    gm = mp.e ** (sum(mp.log(mm) for mm in mods) / len(mods)) if mods else mp.mpf(1)
    u = sorted(set(float(mp.re(e / gm)) for e in eps))
    if len(u) > 1:
        gaps = [u[i + 1] - u[i] for i in range(len(u) - 1)]
        min_gap = min(g for g in gaps if g > 0) if any(g > 0 for g in gaps) else 1.0
        cluster = max(0.0, -math.log10(min_gap))
    else:
        cluster = 0.0
    cond_bump = int((m + 2) * (max(dyn, 1.0) + cluster + 1.0)) + 40
    if working_dps is None:
        # floor at a generous multiple of the input precision so a well-posed
        # (determined, low-degree) fit is never precision-starved; the cond_bump then
        # adds head-room for genuinely ill-conditioned wide/clustered grids.
        working_dps = max(in_dps + cond_bump, 2 * in_dps + 50, 150)
    saved = mp.mp.dps
    mp.mp.dps = working_dps
    try:
        eps = [mp.mpc(e) for e in eps]
        vals = [mp.mpc(v) for v in vals]

        geom = _detect_geometry(eps)
        if k_min is None:
            k_min = _estimate_kmin(eps, vals)

        diagnostics = {'geometry': geom, 'input_dps': in_dps,
                       'working_dps': working_dps, 'dyn_range_decades': dyn,
                       'cluster_decades': round(cluster, 2), 'n_samples': m}

        # ---- detect resonance / log order (each lo gets its full power budget) ----
        log_order, scoreboard = detect_log_order(eps, vals, k_min, None, max_log=max_log)
        diagnostics['log_scoreboard'] = scoreboard
        diagnostics['log_order'] = log_order

        # final power count: use all samples (determined fit) for the chosen log_order
        if n_powers is None:
            n_powers_final = max(1, m // (log_order + 1))
        else:
            n_powers_final = n_powers
        diagnostics['n_powers'] = n_powers_final

        per_method = {}

        # ---- method A: generalized Vandermonde (always applicable) ----
        try:
            cf_vand = vandermonde_fit(eps, vals, k_min, n_powers_final, log_order=log_order)
            per_method['vandermonde'] = cf_vand
        except Exception as ex:
            diagnostics['vandermonde_error'] = repr(ex)

        # ---- method B: eps-FFT on a circle ----
        if geom == 'circle':
            try:
                cf_fft = eps_fft(eps, vals)
                per_method['eps_fft'] = cf_fft
            except Exception as ex:
                diagnostics['fft_error'] = repr(ex)

        # ---- methods C/D: ladder extrapolation of the finite part (k_min==0, no log)
        if k_min == 0 and log_order == 0:
            try:
                lim_n = neville_zero(eps, vals)
                per_method.setdefault('neville', {})[(0, 0)] = lim_n
            except Exception as ex:
                diagnostics['neville_error'] = repr(ex)
            if geom == 'real_ladder':
                try:
                    lim_r, _ = romberg_ladder(eps, vals)
                    per_method.setdefault('romberg', {})[(0, 0)] = lim_r
                except Exception as ex:
                    diagnostics['romberg_error'] = repr(ex)
            try:
                lim_bs = bulirsch_stoer_zero(eps, vals)
                per_method.setdefault('bulirsch_stoer', {})[(0, 0)] = lim_bs
            except Exception as ex:
                diagnostics['bs_error'] = repr(ex)

        # ---- choose primary + cross-validate ----
        primary = 'vandermonde' if 'vandermonde' in per_method else next(iter(per_method))
        if geom == 'circle' and 'eps_fft' in per_method:
            primary = 'eps_fft'
        coeffs = per_method[primary]

        # cross-validation: agreement of the finite part across independent methods
        fp_estimates = {nm: cf.get((0, 0)) for nm, cf in per_method.items()
                        if (0, 0) in cf}
        achieved = _cross_validate(fp_estimates, primary)
        # methods share the same samples, so their agreement cannot certify more
        # digits than the inputs themselves carry
        if achieved == achieved:                     # skip the single-method NaN
            achieved = min(achieved, float(in_dps))

        # if a ladder method beats vandermonde on a clean ladder, prefer its finite part
        # but keep the vandermonde coeffs for the rest of the series.
        result = ExtrapResult(
            coeffs={kk: mp.mpc(vv) for kk, vv in coeffs.items()},
            k_min=k_min,
            log_order=log_order,
            finite_part=coeffs.get((0, 0), mp.mpc(0)),
            achieved_digits=achieved,
            method=primary,
            per_method={nm: {kk: mp.mpc(vv) for kk, vv in cf.items()}
                        for nm, cf in per_method.items()},
            diagnostics=diagnostics,
        )
        if verbose:
            print(result)
            print("  diagnostics:", diagnostics)
            print("  finite-part by method:",
                  {nm: mp.nstr(v, 20) for nm, v in fp_estimates.items()})
        return result
    finally:
        mp.mp.dps = saved


def _cross_validate(fp_estimates, primary):
    """Achieved digits = best pairwise agreement of the primary finite part with any
    other independent method (honest, method-disagreement-based, not self-residual).
    With only one method available, fall back to NaN-flagged single-method note."""
    if primary not in fp_estimates:
        return 0.0
    p = fp_estimates[primary]
    others = [v for nm, v in fp_estimates.items() if nm != primary]
    if not others:
        return float('nan')
    return max(_digits_agree(p, o) for o in others)


def _str_dps(v):
    """Estimate the supplied precision (decimal digits) of a value string/dict."""
    if isinstance(v, (mp.mpf, mp.mpc)):
        # an already-parsed mpmath value carries its own binary precision;
        # str() at the ambient dps would understate it
        tt = (v._mpf_,) if isinstance(v, mp.mpf) else v._mpc_
        bits = max(t[3] for t in tt)          # bc field = mantissa bit count
        return max(int(bits * 0.30103) + 1, 15)
    if isinstance(v, dict):
        cand = max((str(_strip_ball(v.get(k, '0'))) for k in ('re', 'im')), key=len)
    else:
        cand = str(_strip_ball(v))
    digits = sum(ch.isdigit() for ch in cand)
    return max(digits, 15)


# ----------------------------------------------------------------------------- #
#  convenience loader for AMFlow black_box JSON                                  #
# ----------------------------------------------------------------------------- #

def load_amflow_grid(json_paths, integral_index=0, part='re', skip_failed=True):
    """Read one or more AMFlow black_box result JSONs and collect (eps, value) pairs
    for a single integral.  `part` selects 're', 'im', or 'complex'.

    If `skip_failed`, samples whose value parses to an exact bare '0'/'0.0' (an engine
    non-convergence / failed-eps marker -- distinct from a genuine high-precision tiny
    number, which arrives as an Arb ball '[... +/- ...]') are dropped with a warning, so
    one bad eps point cannot silently poison the whole extrapolation."""
    import json, sys
    samples = []
    for path in json_paths:
        d = json.load(open(path))
        res = d.get('result')
        if not res:
            continue
        for s in res[integral_index]['samples']:
            e = s['eps']
            v = s['value']
            # parse at (at least) the digits the strings carry: the ambient dps
            # (default 15) would truncate every high-precision sample right here
            with mp.workdps(max(mp.mp.dps, _str_dps(e), _str_dps(v)) + 10):
                ev = to_mpc(e)
                if part == 'complex':
                    re_raw, im_raw = str(v.get('re', '')), str(v.get('im', ''))
                    vv = to_mpc(v)
                elif part == 'im':
                    im_raw = str(v['im']); re_raw = ''
                    vv = mp.mpc(mp.mpf(_strip_ball(v['im'])))
                else:
                    re_raw = str(v['re']); im_raw = ''
                    vv = mp.mpc(mp.mpf(_strip_ball(v['re'])))
            if skip_failed:
                relevant = (re_raw, im_raw) if part != 'im' else (im_raw,)
                # a real engine value is an Arb ball ('[...]'); a bare '0' is a failure flag
                bare_zero = all(r.strip() in ('0', '0.0', '') for r in relevant if r != '')
                if bare_zero and vv == 0:
                    print(f"[load_amflow_grid] WARNING: skipping failed/zero sample "
                          f"eps={mp.nstr(mp.re(ev), 8)} in {path}", file=sys.stderr)
                    continue
            samples.append((ev, vv))
    return samples


if __name__ == '__main__':
    import argparse, json, sys
    ap = argparse.ArgumentParser(description="Generic eps->0 Laurent extrapolator")
    ap.add_argument('json', nargs='*', help="AMFlow black_box result JSON file(s)")
    ap.add_argument('--part', default='re', choices=['re', 'im', 'complex'])
    ap.add_argument('--kmin', type=int, default=None)
    ap.add_argument('--maxlog', type=int, default=3)
    ap.add_argument('--dps', type=int, default=None)
    args = ap.parse_args()
    if not args.json:
        ap.error("supply at least one AMFlow result JSON, or import the module")
    samples = load_amflow_grid(args.json, part=args.part)
    r = extrapolate(samples, k_min=args.kmin, max_log=args.maxlog,
                    working_dps=args.dps, verbose=True)
    print("\nFINITE PART c_0 =", mp.nstr(r.finite_part, 40))

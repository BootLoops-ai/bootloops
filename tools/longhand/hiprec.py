# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision.
r"""
hiprec.py -- the engines of Longhand's route A: the arbitrary-precision
Feynman-parametric integral evaluator.

A SECOND, independent numerical ground-truth route for finite (eps^0) Feynman
integrals, built to break the ~8-9 digit double-precision ceiling of pySecDec
disteval / FIESTA double mode.  It evaluates the parametric representation

    I = Gamma(N - L d/2) * \int_{simplex} dx  delta(1-sum x)  U^{N-(L+1)d/2} / F^{N-Ld/2}

at the requested regulator order directly in arbitrary precision, using two
deterministic engines on a smooth, positive Euclidean integrand:

  * ENGINE "ts"   : nested arbitrary-precision tanh-sinh (mpmath.quad).  Spectral
                    for analytic integrands; the method of choice once the
                    effective dimension is <= 3 (e.g. the one-loop box after the
                    closed-form loop-parameter integration -> 40+ digits).

  * ENGINE "qmc"  : a rank-1 lattice rule whose generating vector is built by
                    fast component-by-component (CBC) construction, with a
                    Korobov periodizing transform and M independent random
                    shifts giving a statistical (median / std-error) uncertainty estimate.
                    This is the workhorse for the irreducible higher-dimensional
                    integrals (e.g. the 3-loop self-energy benchmark, which is
                    6-dimensional after one analytic loop-parameter
                    integration).  Multiprocess.

DIMENSION REDUCTION.  Before numerics, every Feynman parameter that enters the
*denominator polynomial linearly* is integrated out in closed form
(\int_0^inf dx/(a x + b)^2 = 1/(a b) in the Cheng-Wu chart; the analogous
1/(quadratic)^2 closed form when a variable is quadratic with smooth, sign-fixed
discriminant -- the "box" reduction).  Each such step removes one numerical
dimension exactly and is what lets the box reach 40 digits and the 3-loop
self-energy reach the QMC regime at all.

REGION.  Deep-Euclidean (s,t < 0, internal m^2 > 0), where F > 0 on the whole
open simplex and the integrand is real-analytic -- no contour deformation, no
threshold crossing.  This is the regime in which an independent >=30 d cross
check is actually needed (AMFlow's eps->0 extrapolation is precision-limited,
disteval caps at ~8 d).

See GUIDE.md (the package) and PARAMETRIC.md (this route's manual, with the
benchmark numbers).
"""
import multiprocessing as mp_proc
from math import comb, factorial, gcd

import gmpy2
from gmpy2 import mpfr


# --------------------------------------------------------------------------- #
#  CBC rank-1 lattice generating vector (product weights, fast Nuyens-Cools)
# --------------------------------------------------------------------------- #
def cbc_generating_vector(N, dim, alpha=2, gamma=0.8):
    """Component-by-component construction of a rank-1 lattice generating vector
    for prime N, product weights gamma, smoothness order alpha (kernel = scaled
    Bernoulli B_{2 alpha}).  Returns a python list of ints (z[0]=1).

    Pure-python/float.  NOT cheap at large N: cost is O(dim * N * phi(N))
    single-core (measured ~15 min at N=64007, dim=12 under load), so the
    result is cached per (N, dim, alpha, gamma) — it is deterministic.  The
    cache files read, in order: the vendored cbc_cache.json beside this
    module (read-only; ships the vectors the benchmarks and selftests use),
    then the user cache ($LONGHAND_CBC_CACHE, else the older name
    $HIPREC_CBC_CACHE, else $XDG_CACHE_HOME or ~/.cache under
    longhand/cbc_cache.json), then the cache location of this module's
    earlier packaging (hiprec-sectordecomp/cbc_cache.json under the same
    cache root, read-only, so vectors computed before the rename are still
    found); new vectors are written to the user cache only, so a run never
    modifies the package tree.
    """
    import numpy as np
    import json as _json
    import os as _os
    _shipped = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)),
                             'cbc_cache.json')
    _cache_root = (_os.environ.get('XDG_CACHE_HOME')
                   or _os.path.join(_os.path.expanduser('~'), '.cache'))
    _cache = (_os.environ.get('LONGHAND_CBC_CACHE')
              or _os.environ.get('HIPREC_CBC_CACHE')
              or _os.path.join(_cache_root, 'longhand', 'cbc_cache.json'))
    _legacy = _os.path.join(_cache_root, 'hiprec-sectordecomp', 'cbc_cache.json')
    _key = f"{N}_{dim}_{alpha}_{gamma}"
    _d = {}
    for _path in (_shipped, _cache, _legacy):
        try:
            with open(_path) as _fh:
                _got = _json.load(_fh)
            if _path == _cache:
                _d = _got
            if _key in _got:
                return list(_got[_key])
        except (OSError, ValueError):
            pass
    j = np.arange(N)
    x = j / N
    if alpha == 1:
        B = x**2 - x + 1.0/6.0
        c = (2*np.pi)**2 / 2.0
    elif alpha == 2:
        B = x**4 - 2*x**3 + x**2 - 1.0/30.0
        c = (2*np.pi)**4 / 24.0
    elif alpha == 3:
        B = x**6 - 3*x**5 + 2.5*x**4 - 0.5*x**2 + 1.0/42.0
        c = (2*np.pi)**6 / 720.0
    else:
        raise ValueError("alpha must be 1,2,3")
    omega = -c * B
    z = [1] * dim
    Pprev = np.ones(N)
    cand = [g for g in range(1, N) if gcd(g, N) == 1]
    for sdim in range(dim):
        if sdim == 0:
            z[0] = 1
            Pprev = Pprev * (1 + gamma * omega[(1 * j) % N])
            continue
        best, bestval = None, None
        for g in cand:
            val = float(np.mean(Pprev * (1 + gamma * omega[(g * j) % N])))
            if bestval is None or val < bestval:
                bestval, best = val, g
        z[sdim] = best
        Pprev = Pprev * (1 + gamma * omega[(best * j) % N])
    try:
        _d[_key] = [int(v) for v in z]
        _os.makedirs(_os.path.dirname(_cache), exist_ok=True)
        _tmp = _cache + '.tmp'
        with open(_tmp, 'w') as _fh:
            _json.dump(_d, _fh)
        _os.replace(_tmp, _cache)
    except OSError:
        pass
    return z


# --------------------------------------------------------------------------- #
#  Korobov periodizing transform  phi_p : (0,1) -> (0,1), weight w_p
# --------------------------------------------------------------------------- #
def _korobov(p):
    Bnorm = mpfr(factorial(p))**2 / factorial(2*p+1)
    coeffs = [mpfr(((-1)**k) * comb(p, k)) / (p+1+k) for k in range(p+1)]

    def phi(u):
        acc = mpfr(0)
        up = u**(p+1)
        for k in range(p+1):
            acc += coeffs[k]*up
            up *= u
        return acc / Bnorm

    def w(u):
        return (u*(1-u))**p / Bnorm
    return phi, w


# --------------------------------------------------------------------------- #
#  QMC engine
# --------------------------------------------------------------------------- #
def _qmc_one_shift(args):
    """One randomly-shifted lattice estimate over [0,inf)^dim.
    integrand_factory() must return a callable f(*x) -> gmpy2 mpfr (the FULL
    integrand on [0,inf)^dim, INCLUDING any Gamma prefactor).  We supply the
    Korobov tail/periodization Jacobian.
    """
    (N, gvec, shift, dim, p, prec, integrand_factory, fac_arg) = args
    gmpy2.get_context().precision = prec
    f = integrand_factory(fac_arg)
    phi, w = _korobov(p)
    one = mpfr(1)
    tot = mpfr(0)
    Nm = mpfr(N)
    sh = [mpfr(v) for v in shift]
    for jj in range(N):
        jac = one
        xv = []
        for k in range(dim):
            u = mpfr((gvec[k]*jj) % N) / Nm + sh[k]
            if u >= 1:
                u -= 1
            up = phi(u)
            ww = w(u)
            den = one - up
            xv.append(up/den)
            jac *= ww/den**2
        tot += jac * f(*xv)
    return tot / Nm


def integrate_qmc(integrand_factory, fac_arg, dim, N=64007, n_shifts=16,
                  korobov_p=3, dps=40, alpha=2, gamma=0.8, nproc=None, seed=1):
    """Median-of-random-shifts CBC lattice QMC over [0,inf)^dim in arbitrary precision.

    integrand_factory(fac_arg) -> f(*x): the full integrand on [0,inf)^dim, returning
    a gmpy2 mpfr.  (A factory is used so each worker process builds its own
    closed-over compiled integrand -- gmpy2/sympy lambdas do not always pickle.)

    Returns dict: value (mpfr), error (mpfr, std of the shift mean), digits (int),
                  shifts (list of per-shift estimates).
    """
    import random
    prec = int(dps * 3.3219) + 60
    gmpy2.get_context().precision = prec
    z = cbc_generating_vector(N, dim, alpha=alpha, gamma=gamma)
    rng = random.Random(seed)
    shifts = [[rng.randint(1, 10**9) / 10**9 for _ in range(dim)]
              for _ in range(n_shifts)]
    tasks = [(N, z, sh, dim, korobov_p, prec, integrand_factory, fac_arg)
             for sh in shifts]
    if nproc is None:
        nproc = min(n_shifts, mp_proc.cpu_count())
    if nproc > 1:
        with mp_proc.Pool(nproc) as pool:
            ests = pool.map(_qmc_one_shift, tasks)
    else:
        ests = [_qmc_one_shift(tk) for tk in tasks]
    gmpy2.get_context().precision = prec
    m = len(ests)
    mean = sum(ests) / m
    if m > 1:
        var = sum((e-mean)**2 for e in ests) / (m*(m-1))
        err = gmpy2.sqrt(var)
    else:
        err = mpfr('nan')
    digits = (-int(gmpy2.log10(err)) if err > 0 and gmpy2.is_finite(err) else 0)
    return {'value': mean, 'error': err, 'digits': digits, 'shifts': ests,
            'N': N, 'n_shifts': m, 'gvec': z}


# --------------------------------------------------------------------------- #
#  tanh-sinh engine (mpmath) -- for effective dimension <= ~3
# --------------------------------------------------------------------------- #
_GL_CACHE = {}


def _gl(n, work):
    """Cached Gauss-Legendre nodes/weights on [-1,1] at working precision `work`."""
    import mpmath as mp
    key = (n, work)
    xw = _GL_CACHE.get(key)
    if xw is None:
        old = mp.mp.dps
        mp.mp.dps = work
        xw = mp.gauss_quadrature(n, 'legendre')
        mp.mp.dps = old
        _GL_CACHE[key] = xw
    return xw


def integrate_gauss_product(f, dim, dps=40, nodes=None, limits=None):
    """FIXED (non-adaptive) tensored Gauss-Legendre product rule over [0,inf)^dim.

    Deterministic -- no adaptive error noise -- so it converges SPECTRALLY (and
    reproducibly) for integrands that are real-analytic and smooth on the whole
    [0,inf)^dim path (the deep-Euclidean regime).  This is the method that takes
    the one-loop box to 40+ digits; prefer it over the adaptive nested tanh-sinh
    for dim <= 3, whose nested error floor is ~17 d.

    Each axis k is mapped onto [lo_k, inf) <- [-1,1] by x = lo_k + L_k*(1+u)/(1-u).
    `limits(fixed) -> (lo_k, L_k)` may depend on the already-fixed OUTER variables;
    matching the lower limit and tail scale to where the integrand actually has
    support is what restores spectral convergence (e.g. for the box the inner
    variable lives on [x1+1, inf) with scale 1+x1, not [0,inf) with a fixed
    scale -- with the matched limits the rule reaches 40+ d, with [0,inf)/L=2 it
    plateaus at ~5 d).  Default limits = ([0,inf), L=2) on every axis.

    `nodes` per axis defaults to a spectral count linear in `dps`.  f takes `dim`
    mpmath.mpf args and must be smooth & decaying.  Returns an mpf good to ~dps d.
    """
    import mpmath as mp
    work = dps + 25
    mp.mp.dps = work
    if nodes is None:
        nodes = int(1.4*work) + 14
    one = mp.mpf(1)
    gl = _gl(nodes, work)
    if limits is None:
        def limits(fixed):
            return (mp.mpf(0), mp.mpf(2))

    def rec(level, fixed):
        if level == dim:
            return f(*fixed)
        lo, L = limits(fixed)
        lo = mp.mpf(lo); L = mp.mpf(L)
        tot = mp.mpf(0)
        for (tnode, w) in zip(*gl):
            x = lo + L*(one+tnode)/(one-tnode)
            jac = L*2/(one-tnode)**2
            tot += w*jac*rec(level+1, fixed+[x])
        return tot
    val = rec(0, [])
    mp.mp.dps = dps
    return +val


def integrate_tanhsinh(f, dim, dps=40, intervals=None):
    """Nested arbitrary-precision tanh-sinh (mpmath.quad) of f over the box
    `intervals` (default [0,inf)^dim).  Spectral for analytic integrands; use
    only for dim <= 3.  f takes `dim` mpmath.mpf arguments.
    """
    import mpmath as mp
    mp.mp.dps = dps + 15
    if intervals is None:
        intervals = [[0, mp.inf]] * dim

    def nest(level, fixed):
        if level == dim:
            return f(*fixed)
        return mp.quad(lambda x: nest(level+1, fixed+[x]), intervals[level])
    val = nest(0, [])
    mp.mp.dps = dps
    return +val

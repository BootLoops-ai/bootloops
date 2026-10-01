"""blade.ratrec — standalone multivariate rational reconstruction over
fflow primes, driven through Blade's recmod + dynamicrr binaries.

THE CATCH (measured in Phase 0): recmod needs evaluations AT fflow's own
structured sample points (Newton/Thiele grids derived from the degree info
with the SAME opt 7-tuple {0,0,0,pid,1,maxdeg,0} that recmod passes to
reconstruct_mod).  Free random points DO NOT work.  So reconstruct() first
dumps the required points per prime (dumppoints shim), hands each point to
the caller's evaluator, fabricates the eval files, then drives
recmod (per prime) -> dynamicrr (CRT + rational reconstruction).

SECOND CATCH (measured while building this): the degrees file recmod/
dumppoints consume must carry the EXACT total num/den degrees -- an
over-estimate makes recmod fail outright.  reconstruct() therefore first
LEARNS the degrees from the evaluator (scan_degrees: validated univariate
rational fits along random lines) unless the caller supplies exact ones.

API:
    reconstruct(evaluator, nvars, maxdeg, nfuns=1, primes=(0,1,2), ...)
        evaluator(coords: tuple[int], prime: int) -> int | sequence[int]
            values mod prime of the nfuns target functions at the point.
        primes: fflow prime IDS, must be consecutive 0..k-1 (dynamicrr
            indexes BIG_UINT_PRIMES by consecutive id, dynamicrr.cc:101-103),
            k >= 2 (prime-0 fit + fresh-prime verify).
        -> list of sympy expressions (or a single expression if nfuns == 1).

Every stage is artifact-gated (see pipeline.py); a corrupted evaluation
surfaces as BladeGateError (recmod pattern mismatch across primes, or
dynamicrr requesting a prime beyond those staged), never as a silently
wrong answer certified by rc/flags.
"""

from __future__ import annotations

import os
import tempfile
from typing import Callable, List, Optional, Sequence

from . import formats, pipeline
from .formats import (BladeFormatError, DegreeInfo, DegreesFile, EvalFile,
                      EvalList, PointsFile, RecCoeff, RecMono, RecPart, RRRes,
                      assemble_rational_functions, big_uint_primes, flags_size)
from .pipeline import (DEFAULT_BIN_DIR, PHASE0_DUMPPOINTS, BladeGateError,
                       run_dynamicrr_generic, run_recmod_generic, _run)


def _dumppoints_bin(bin_dir: str, dumppoints_bin: Optional[str]) -> str:
    if dumppoints_bin:
        return dumppoints_bin
    cand = os.path.join(bin_dir, "dumppoints")
    return cand if os.path.exists(cand) else PHASE0_DUMPPOINTS


def dump_points(nvars: int, nfuns: int, pid: int, maxdeg: int,
                degrees_path: str, points_path: str,
                bin_dir: str = DEFAULT_BIN_DIR,
                dumppoints_bin: Optional[str] = None) -> PointsFile:
    """Dump the fflow-internal structured sample points that reconstruction
    at prime id `pid` will request.  Gated like pipeline.run_dumppoints."""
    binp = _dumppoints_bin(bin_dir, dumppoints_bin)
    res = _run([binp, str(nvars), str(nfuns), str(pid), str(maxdeg),
                degrees_path, points_path],
               cwd=os.path.dirname(points_path) or ".", timeout=600.0,
               stage=f"dumppoints[p{pid}]")
    if not os.path.exists(points_path):
        raise BladeGateError(f"dumppoints: {points_path} missing. "
                             f"log tail: {res.log[-300:]}")
    pts = PointsFile.read(points_path)
    fs = flags_size(nfuns)
    if pts.n != nvars + fs or not pts.rows:
        raise BladeGateError(f"dumppoints: {points_path}: header n={pts.n} "
                             f"(want {nvars}+{fs}) with {len(pts.rows)} samples")
    p = big_uint_primes()[pid]
    if any(r[nvars] != p for r in pts.rows):
        raise BladeGateError(f"dumppoints: {points_path}: row prime != "
                             f"BIG_UINT_PRIMES[{pid}]={p}")
    return pts


def _evals_from_callback(pts: PointsFile, nvars: int, nfuns: int, prime: int,
                         evaluator: Callable) -> EvalFile:
    """Fabricate the ssolve-format eval file (iofflow.cpp:205) from evaluator
    values at the dumped points: row = point words + needed outputs compacted
    per the flags, zero-padded to nfuns."""
    fs = flags_size(nfuns)
    ev = EvalFile(nvars, nfuns + fs, [])
    for i in range(len(pts.rows)):
        coords, p, flags = pts.split_row(i, nvars)
        vals = evaluator(tuple(coords), prime)
        if isinstance(vals, int):
            vals = [vals]
        vals = list(vals)
        if len(vals) != nfuns:
            raise BladeGateError(f"evaluator returned {len(vals)} values, "
                                 f"expected nfuns={nfuns}")
        for j, v in enumerate(vals):
            if not (0 <= v < prime):
                raise BladeGateError(f"evaluator value {v} for output {j} not "
                                     f"reduced mod {prime}")
        row = list(coords) + [p] + list(flags)
        needed = [vals[j] for j in range(nfuns)
                  if (flags[j // 64] >> (j % 64)) & 1]
        row += needed + [0] * (nfuns - len(needed))
        ev.rows.append(row)
    return ev


# --------------------------------------------------------------------------
# degree scan (pure python, from the evaluator itself)
#
# MEASURED: recmod FAILS if the degrees file
# over-estimates the TOTAL numerator/denominator degrees (probe 'totplus'),
# while over-estimated PER-VAR max degrees are tolerated (probe 'varbump').
# DegreesFile.uniform(maxdeg) is therefore NOT a valid recmod input; the
# degrees must be learned.  fflowml_alg_all_degrees needs an evaluatable
# graph (not a dummy-graph + eval files), so we scan in python: univariate
# rational fits along a generic line (totals) and along each axis line
# (per-var max/min), every fit validated at held-out points.
# --------------------------------------------------------------------------

def _nullspace_vec_modp(A, ncols, p):
    """One nonzero nullspace vector of A (rows of length ncols) mod p, or
    None if the nullspace is trivial."""
    A = [row[:] for row in A]
    pivots = {}   # col -> row index
    r = 0
    for c in range(ncols):
        piv = None
        for i in range(r, len(A)):
            if A[i][c] % p:
                piv = i
                break
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        inv = pow(A[r][c], -1, p)
        A[r] = [v * inv % p for v in A[r]]
        for i in range(len(A)):
            if i != r and A[i][c]:
                f = A[i][c]
                A[i] = [(v - f * w) % p for v, w in zip(A[i], A[r])]
        pivots[c] = r
        r += 1
        if r == len(A):
            break
    free = [c for c in range(ncols) if c not in pivots]
    if not free:
        return None
    c0 = free[0]
    v = [0] * ncols
    v[c0] = 1
    for c, row in pivots.items():
        v[c] = (-A[row][c0]) % p
    return v


def _poly_deg(c):
    for i in range(len(c) - 1, -1, -1):
        if c[i]:
            return i
    return -1


def _poly_mod_gcd(a, b, p):
    a, b = a[:_poly_deg(a) + 1], b[:_poly_deg(b) + 1]
    while b and _poly_deg(b) >= 0:
        # a mod b
        da, db = _poly_deg(a), _poly_deg(b)
        if da < 0:
            break
        inv = pow(b[db], -1, p)
        r = a[:]
        for k in range(da - db, -1, -1):
            if _poly_deg(r) < db + k:
                continue
            if len(r) > db + k and r[db + k]:
                f = r[db + k] * inv % p
                for j in range(db + 1):
                    r[j + k] = (r[j + k] - f * b[j]) % p
        r = r[:_poly_deg(r) + 1]
        a, b = b, r
    return a if a else [0]


def _poly_divexact(a, g, p):
    da, dg = _poly_deg(a), _poly_deg(g)
    if dg < 0:
        raise BladeGateError("degree scan: division by zero polynomial")
    q = [0] * (da - dg + 1)
    r = a[:]
    inv = pow(g[dg], -1, p)
    for k in range(da - dg, -1, -1):
        coef = (r[dg + k] * inv) % p if len(r) > dg + k else 0
        q[k] = coef
        if coef:
            for j in range(dg + 1):
                r[j + k] = (r[j + k] - coef * g[j]) % p
    if any(x % p for x in r):
        raise BladeGateError("degree scan: inexact polynomial division")
    return q


def _poly_eval(c, t, p):
    acc = 0
    for v in reversed(c):
        acc = (acc * t + v) % p
    return acc


def _univ_ratfit(ts, fs, dmax, p):
    """Fit fs[i] = P(ts[i])/Q(ts[i]) with deg P, deg Q <= dmax; return the
    REDUCED (P, Q) coefficient lists.  Uses len(ts) >= 2*dmax+2+4 points; the
    last 4 are held out for validation.  Raises BladeGateError on failure."""
    need = 2 * dmax + 2
    fit_n = len(ts) - 4
    if fit_n < need:
        raise BladeGateError("degree scan: not enough points for the fit")
    if all(f == 0 for f in fs):
        return [0], [1]
    rows = []
    for t, f in zip(ts[:fit_n], fs[:fit_n]):
        tp = [1]
        for _ in range(dmax):
            tp.append(tp[-1] * t % p)
        rows.append([x % p for x in tp] + [(-f * x) % p for x in tp])
    v = _nullspace_vec_modp(rows, 2 * (dmax + 1), p)
    if v is None:
        raise BladeGateError(f"degree scan: no rational fit with deg <= {dmax} "
                             f"(raise maxdeg)")
    P, Q = v[:dmax + 1], v[dmax + 1:]
    if _poly_deg(Q) < 0:
        raise BladeGateError("degree scan: degenerate fit (Q == 0)")
    g = _poly_mod_gcd(P if _poly_deg(P) >= 0 else Q, Q, p)
    if _poly_deg(P) >= 0 and _poly_deg(g) > 0:
        P, Q = _poly_divexact(P, g, p), _poly_divexact(Q, g, p)
    P, Q = P[:_poly_deg(P) + 1] or [0], Q[:_poly_deg(Q) + 1]
    for t, f in zip(ts[fit_n:], fs[fit_n:]):
        qv = _poly_eval(Q, t, p)
        if qv == 0 or _poly_eval(P, t, p) * pow(qv, -1, p) % p != f:
            raise BladeGateError("degree scan: fit fails held-out validation "
                                 "(corrupt evaluations or maxdeg too small)")
    return P, Q


def scan_degrees(evaluator: Callable, nvars: int, nfuns: int, maxdeg: int,
                 pid: int = 0, seed: int = 20260708) -> DegreesFile:
    """Learn the EXACT degree info recmod needs, from the evaluator itself.
    Totals from a generic parametric line x = b + t*d; per-var max/min from
    axis lines (other coords random).  All fits validated at 4 held-out
    points each.  Evaluator cost: (1 + nvars) * (2*maxdeg + 6) calls.
    Degrees are exact w.h.p. (random 63-bit lines); a bad draw fails the
    held-out validation loudly rather than mis-reporting.
    MEASURED: when maxdeg exceeds the true degrees, ONE corrupt evaluation is
    ABSORBED as a spurious (t - t_bad) common factor and gcd-stripped -- the
    scan then reports the correct degrees (held-out-validated), it does not
    go silently wrong."""
    import random
    rng = random.Random(seed)
    p = big_uint_primes()[pid]
    npts = 2 * maxdeg + 2 + 4

    def sample_line(coords_of_t):
        ts, vals = [], []
        seen = set()
        while len(ts) < npts:
            t = rng.randrange(1, p)
            if t in seen:
                continue
            seen.add(t)
            out = evaluator(coords_of_t(t), p)
            out = [out] if isinstance(out, int) else list(out)
            if len(out) != nfuns:
                raise BladeGateError(f"evaluator returned {len(out)} values, "
                                     f"expected {nfuns}")
            ts.append(t)
            vals.append(out)
        return ts, vals

    # totals: generic line
    b = [rng.randrange(1, p) for _ in range(nvars)]
    d = [rng.randrange(1, p) for _ in range(nvars)]
    ts, vals = sample_line(lambda t: tuple((bi + t * di) % p
                                           for bi, di in zip(b, d)))
    totals = []
    for j in range(nfuns):
        P, Q = _univ_ratfit(ts, [v[j] for v in vals], maxdeg, p)
        totals.append((max(_poly_deg(P), 0), max(_poly_deg(Q), 0)))

    # per-var: axis lines
    pervar = [[None] * nvars for _ in range(nfuns)]
    for v_i in range(nvars):
        base = [rng.randrange(1, p) for _ in range(nvars)]

        def coords(t, v_i=v_i, base=base):
            c = list(base)
            c[v_i] = t
            return tuple(c)

        ts, vals = sample_line(coords)
        for j in range(nfuns):
            P, Q = _univ_ratfit(ts, [v[j] for v in vals], maxdeg, p)
            num_max, den_max = max(_poly_deg(P), 0), max(_poly_deg(Q), 0)
            num_min = next((i for i, c in enumerate(P) if c), 0)
            den_min = next((i for i, c in enumerate(Q) if c), 0)
            pervar[j][v_i] = (num_max, num_min, den_max, den_min)

    info = [DegreeInfo(totals[j][0], totals[j][1], pervar[j])
            for j in range(nfuns)]
    return DegreesFile(nvars, nfuns, info)


def reconstruct(evaluator: Callable, nvars: int, maxdeg: int, nfuns: int = 1,
                primes: Sequence[int] = (0, 1, 2), workdir: Optional[str] = None,
                nthreads: int = 1, degrees: Optional[DegreesFile] = None,
                bin_dir: str = DEFAULT_BIN_DIR,
                dumppoints_bin: Optional[str] = None,
                symbols=None, maxdeg_cli: Optional[int] = None,
                scan_seed: int = 20260708):
    """Reconstruct nfuns multivariate rational functions from a black-box
    modular evaluator.  See module docstring.  Returns sympy expression(s).

    degrees: optional EXACT DegreesFile (skips the scan; fewer evaluator
    calls).  Default: scan_degrees() learns them from the evaluator --
    NECESSARY, because recmod fails when the declared TOTAL degrees
    over-estimate the true ones (measured; see scan_degrees docstring).
    maxdeg must be >= the true num/den total degrees or the scan aborts.
    """
    primes = list(primes)
    if primes != list(range(len(primes))):
        raise BladeGateError(f"primes must be consecutive fflow prime ids "
                             f"0..k-1 (dynamicrr indexes by id); got {primes}")
    if len(primes) < 2:
        raise BladeGateError("need >= 2 primes (prime-0 fit + fresh-prime verify)")
    if nvars > 9:
        raise BladeGateError("ptsord digit encoding limits nvars <= 9")

    own_dir = workdir is None
    if own_dir:
        workdir = tempfile.mkdtemp(prefix="blade_ratrec_")
    os.makedirs(workdir, exist_ok=True)

    deg = degrees or scan_degrees(evaluator, nvars, nfuns, maxdeg,
                                  pid=primes[0], seed=scan_seed)
    if deg.nparsin != nvars or deg.nparsout != nfuns:
        raise BladeGateError("degrees file dims disagree with nvars/nfuns")
    degrees_path = os.path.join(workdir, "degrees.fflow")
    deg.write(degrees_path)
    cli_maxdeg = maxdeg_cli if maxdeg_cli is not None else maxdeg

    prefix = "rec"
    for pid in primes:
        prime = big_uint_primes()[pid]
        pts = dump_points(nvars, nfuns, pid, cli_maxdeg, degrees_path,
                          os.path.join(workdir, f"points_p{pid}.fflow"),
                          bin_dir, dumppoints_bin)
        ev = _evals_from_callback(pts, nvars, nfuns, prime, evaluator)
        eval_path = os.path.join(workdir, f"eval_p{pid}_1.fflow")
        ev.write(eval_path)
        want = EvalFile.expected_bytes(nvars, nfuns, len(pts.rows))
        got = os.path.getsize(eval_path)
        if got != want:
            raise BladeGateError(f"ratrec: {eval_path}: {got}B != expected {want}B")
        EvalList([os.path.abspath(eval_path)]).write(
            os.path.join(workdir, f"evallist_p{pid}.txt"))
        run_recmod_generic(os.path.join(bin_dir, "recmod"), nvars, nfuns, pid,
                           nthreads, degrees_path,
                           os.path.join(workdir, f"evallist_p{pid}.txt"),
                           workdir, prefix, cli_maxdeg)

    coeff0 = RecCoeff.read(os.path.join(workdir, f"{prefix}_0_coeff"))
    run_dynamicrr_generic(os.path.join(bin_dir, "dynamicrr"),
                          len(coeff0.tokens), nthreads, workdir, prefix,
                          "rrres", len(primes))

    part = RecPart.read(os.path.join(workdir, f"{prefix}_0_part"))
    mono = RecMono.read(os.path.join(workdir, f"{prefix}_0_mono"), nvars)
    rr = RRRes.read(os.path.join(workdir, "rrres"))
    funcs = assemble_rational_functions(part, mono, rr.values, symbols)
    return funcs[0] if nfuns == 1 else funcs


# --------------------------------------------------------------------------
# helpers for evaluators over exact rational functions (used by self-tests)
# --------------------------------------------------------------------------

def poly_eval_mod(terms, coords, prime: int) -> int:
    """terms: iterable of (Fraction_or_int_coeff, exponent_tuple)."""
    from fractions import Fraction
    acc = 0
    for c, exps in terms:
        c = Fraction(c)
        num = c.numerator % prime
        den = c.denominator % prime
        v = num * pow(den, -1, prime) % prime
        for x, e in zip(coords, exps):
            if e:
                v = v * pow(x % prime, e, prime) % prime
        acc = (acc + v) % prime
    return acc


def ratfun_evaluator(num_terms, den_terms):
    """-> evaluator(coords, prime) for num/den given as term lists.
    Raises BladeGateError if the denominator vanishes at a sample point
    (resample by perturbing degrees, or supply exact degrees)."""
    def ev(coords, prime):
        den = poly_eval_mod(den_terms, coords, prime)
        if den == 0:
            raise BladeGateError(f"denominator vanished at sample point {coords} "
                                 f"mod {prime}")
        return poly_eval_mod(num_terms, coords, prime) * pow(den, -1, prime) % prime
    return ev

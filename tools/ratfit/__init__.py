#!/usr/bin/env python3
r"""ratfit.py — exact rational-function fitting toolkit (python-flint + fractions).

One canonical implementation of the Thiele/Newton/screen family. No sympy
in the hot path.

Members:
  thiele, thiele_loo, newton_interp, ratrecon_with_denom, screen_modp
      the exact-rational fitting core.
  extract_stairs, row_lcm
      stair-ratio and row-LCM denominator-candidate builders
      (row-LCM family proven 31/31 on sub rows).
  thiele_loo_screened, degree_budget
      the screened fit + budget layer.
  parametric_rec_interp, mpoly_eval
      multivariate iterated-Newton lift of a parametric recurrence over a
      tensor-product param grid (off-grid cross-checked) + its evaluator.

Submodules:
  ratfit.thiele_gate   family-level byte-exact node gate (fit/validate/gate + CLI).
  ratfit.li2close      weight-2 rational x log closed-form closer.
  ratfit.degree_alias  failure-axis triage for reconstruction failures that
                       +prime passes never lift: grid alias (extend the grid)
                       vs CRT height (buy primes); pure functions + --selftest,
                       law receipt DEGREE_ALIAS_DISCRIMINATOR.json (package data).

Conventions:
  xs, ys      exact rationals (fractions.Fraction, int, or anything exposing
              .numerator/.denominator)
  polynomials flint.fmpq_poly, coefficients ascending
  qex         {(i,j): {m: [coeff_strs]}} — validated exact denominators per
              DE entry (i,j) and epsilon-order m, as cached by the production fits
"""
from fractions import Fraction
import flint

__all__ = ["coincidence_loci", "design_grid",
           "thiele", "thiele_loo", "thiele_loo_screened", "newton_interp",
           "newton_interp_fmpq", "ratrecon_with_denom", "screen_modp",
           "extract_stairs", "row_lcm", "degree_budget", "parametric_rec_interp",
           "mpoly_eval", "fq", "qpoly", "poly_lcm", "exact_div", "SCREEN_PRIMES"]

# mod-p screen primes: 2^61-1 (Mersenne) + a 2nd 61-bit prime for 0-div fallback.
SCREEN_PRIMES = [(1 << 61) - 1, 2305843009213693967]


def fq(x):
    """Exact rational -> flint.fmpq (canonical_255.py)."""
    if isinstance(x, flint.fmpq):
        return x
    return flint.fmpq(x.numerator, x.denominator)


def qpoly(cs):
    """Ascending coeff list (str/int/Fraction) -> fmpq_poly (probe_knownQ4.Qpoly)."""
    return flint.fmpq_poly([fq(Fraction(c)) for c in cs])


def poly_lcm(a, b):
    """Polynomial lcm via a*b//gcd (flint gcd is monic). (probe_knownQ4.poly_lcm)"""
    g = a.gcd(b)
    return (a * b) // g if g.degree() >= 0 else a * b


def exact_div(a, b):
    """a/b when the division is exact, else None. (probe_knownQ4.exact_div)"""
    q = a // b
    return q if q * b == a else None


def _loo_split(n, n_loo):
    """Deterministic held-out node selection shared by all LOO fits."""
    step = max(1, n // n_loo)
    test_idx = list(range(step // 2, n, step))[:n_loo]
    fit_idx = [i for i in range(n) if i not in test_idx]
    return test_idx, fit_idx


def thiele(xs, ys):
    """Exact Thiele continued-fraction rational interpolation -> (P, Q) fmpq_polys.
    Verbatim port of the origin implementation."""
    if all(y == 0 for y in ys):
        return flint.fmpq_poly([0]), flint.fmpq_poly([1])
    n = len(xs)
    xq = [fq(x) for x in xs]
    yq = [fq(y) for y in ys]
    phi = [list(yq)]
    for k in range(1, n):
        prev = phi[-1]; new = []; died = False
        for ii in range(len(prev) - 1):
            den = prev[ii + 1] - prev[ii]
            if den == 0:
                died = True; break
            new.append((xq[ii + k] - xq[ii]) / den
                       + (phi[-2][ii + 1] if k >= 2 else flint.fmpq(0)))
        if died:
            break
        phi.append(new)
        if len(set(str(v) for v in new)) == 1 and len(new) >= 2:
            break
    K = len(phi) - 1
    a = [phi[0][0]]
    if K >= 1:
        a.append(phi[1][0])
    for k in range(2, K + 1):
        a.append(phi[k][0] - phi[k - 2][0])
    Pm1 = flint.fmpq_poly([1]); P0 = flint.fmpq_poly([a[0]])
    Qm1 = flint.fmpq_poly([0]); Q0 = flint.fmpq_poly([1])
    for k in range(1, K + 1):
        tk = flint.fmpq_poly([-xq[k - 1], 1])
        Pk = flint.fmpq_poly([a[k]]) * P0 + tk * Pm1
        Qk = flint.fmpq_poly([a[k]]) * Q0 + tk * Qm1
        Pm1, P0 = P0, Pk; Qm1, Q0 = Q0, Qk
    return P0, Q0


def newton_interp(xv, yv):
    """Exact Newton-form polynomial interpolation over fmpq, O(n^2).
    Port: canonical_255.py::newton_interp_fmpq / probe_knownQ4.py::newton_interp."""
    m = len(xv)
    dd = list(yv)
    for k in range(1, m):
        for r in range(m - 1, k - 1, -1):
            dd[r] = (dd[r] - dd[r - 1]) / (xv[r] - xv[r - k])
    P = flint.fmpq_poly([dd[-1]])
    for k in range(m - 2, -1, -1):
        P = P * flint.fmpq_poly([-xv[k], 1]) + flint.fmpq_poly([dd[k]])
    return P


newton_interp_fmpq = newton_interp  # canonical_255.py name


def ratrecon_with_denom(xs, ys, Qknown, n_loo=6):
    """Fit Num(t)/Qknown(t) by exact Newton through the fit nodes, LOO-certify
    on held-out nodes, gcd-reduce on success. Returns (Pnum, Qden, ok, deg).
    ok=False also when Qknown vanishes at ANY node (fit or held-out): that
    node's sample would be absorbed (y*Qknown = 0 there) -- the fit-point D check.
    Verbatim port: canonical_255.py ([DENOM] route, mm0_route_223_v2 lineage)."""
    n = len(xs)
    xq = [fq(x) for x in xs]
    yQ = [fq(y) * Qknown(xq[k]) for k, y in enumerate(ys)]
    test_idx, fit_idx = _loo_split(n, n_loo)
    Pnum = newton_interp([xq[i] for i in fit_idx], [yQ[i] for i in fit_idx])
    # the fit-point D check: a node where Qknown vanishes has y*Q = 0 there
    # whatever its sample was -- that sample would be absorbed, so it is refused (ok=False).
    ok = all(Qknown(x) != 0 for x in xq) and all(Pnum(xq[r]) == yQ[r] for r in test_idx)
    if not ok:
        return Pnum, Qknown, False, Pnum.degree()
    g = Pnum.gcd(Qknown)
    if g.degree() > 0:
        Pnum, Qden = Pnum // g, Qknown // g
    else:
        Qden = Qknown
    return Pnum, Qden, True, max(Pnum.degree(), Qden.degree())
def _nmfr(fr, p):
    """Fraction-like -> nmod (canonical_255.py::_nmfr)."""
    return flint.nmod(fr.numerator, p) / flint.nmod(fr.denominator, p)


def _nmfq(c, p):
    """fmpq -> nmod (canonical_255.py::_nmfq)."""
    return flint.nmod(int(c.p), p) / flint.nmod(int(c.q), p)


def screen_modp(xs, ys, Q, n_loo=6):
    """GF(p) screen for a known-denominator fit: exact mod-p Newton through the
    fit nodes + held-out test-node check. ~1000x faster than the fmpq certify.
    No false rejects: arithmetic is exact mod p and any 0-division falls
    through to the 2nd prime (SCREEN_PRIMES); false-accept prob ~ n_loo/p ~
    1e-18 — ALWAYS follow a pass with the exact ratrecon_with_denom certify.
    Verbatim port of the origin implementation (validated standalone,
    sub-rows 18/18)."""
    n = len(xs)
    test_idx, fit_idx = _loo_split(n, n_loo)
    Qc = Q.coeffs()
    for p in SCREEN_PRIMES:
        try:
            xp = [_nmfr(x, p) for x in xs]
            qc = [_nmfq(c, p) for c in Qc]

            def Qv(x):
                a = qc[-1]
                for cc in reversed(qc[:-1]):
                    a = a * x + cc
                return a

            yQ = [_nmfr(y, p) * Qv(xp[k]) for k, y in enumerate(ys)]
            xv = [xp[i] for i in fit_idx]; dd = [yQ[i] for i in fit_idx]
            m = len(xv)
            for k in range(1, m):
                for r in range(m - 1, k - 1, -1):
                    dd[r] = (dd[r] - dd[r - 1]) / (xv[r] - xv[r - k])
            for r in test_idx:
                acc = dd[m - 1]
                for k in range(m - 2, -1, -1):
                    acc = acc * (xp[r] - xv[k]) + dd[k]
                if acc != yQ[r]:
                    return False
            return True
        except ZeroDivisionError:
            continue
    return False


def thiele_loo(xs, ys, n_loo=6):
    """Blind exact Thiele fit, LOO-validated over 3 node orderings; on success
    returns (P, Q, True, deg); on all-fail returns the single full-node exact
    Thiele with ok=False (deg report). Verbatim port: canonical_255.py."""
    n = len(xs)
    test_idx, fit_idx = _loo_split(n, n_loo)
    for od in [fit_idx, list(reversed(fit_idx)), fit_idx[0::2] + fit_idx[1::2]]:
        P, Q = thiele([xs[i] for i in od], [ys[i] for i in od])
        ok = True
        for r in test_idx + fit_idx:
            xr = fq(xs[r]); Qv = Q(xr)
            if Qv == 0 or P(xr) / Qv != fq(ys[r]):
                ok = False; break
        if ok:
            return P, Q, True, max(P.degree(), Q.degree())
    P, Q = thiele(xs, ys)
    return P, Q, False, max(P.degree(), Q.degree())


def _thiele_modp_ok(xs, ys, od, check_idx, p):
    """Mod-p replica of one thiele_loo ordering: build the Thiele reciprocal-
    difference CF over GF(p) on nodes od, then verify P/Q==y at check_idx via
    the 3-term recurrence (no poly construction). ZeroDivisionError propagates
    to the caller (ambiguous mod p -> certify exactly)."""
    if all(ys[i] == 0 for i in od):
        return all(ys[r] == 0 for r in check_idx)
    xq = [_nmfr(xs[i], p) for i in od]
    yq = [_nmfr(ys[i], p) for i in od]
    n = len(xq); phi = [list(yq)]
    for k in range(1, n):
        prev = phi[-1]; new = []; died = False
        for ii in range(len(prev) - 1):
            den = prev[ii + 1] - prev[ii]
            if den == 0:
                died = True; break
            new.append((xq[ii + k] - xq[ii]) / den
                       + (phi[-2][ii + 1] if k >= 2 else flint.nmod(0, p)))
        if died:
            break
        phi.append(new)
        if len(set(str(v) for v in new)) == 1 and len(new) >= 2:
            break
    K = len(phi) - 1
    a = [phi[0][0]]
    if K >= 1:
        a.append(phi[1][0])
    for k in range(2, K + 1):
        a.append(phi[k][0] - phi[k - 2][0])
    for r in check_idx:
        xr = _nmfr(xs[r], p)
        Pm1, P0 = flint.nmod(1, p), a[0]
        Qm1, Q0 = flint.nmod(0, p), flint.nmod(1, p)
        for k in range(1, K + 1):
            Pk = a[k] * P0 + (xr - xq[k - 1]) * Pm1
            Qk = a[k] * Q0 + (xr - xq[k - 1]) * Qm1
            Pm1, P0 = P0, Pk; Qm1, Q0 = Q0, Qk
        if Q0 == 0 or P0 / Q0 != _nmfr(ys[r], p):
            return False
    return True
def thiele_loo_screened(xs, ys, n_loo=6):
    """thiele_loo with a per-ordering mod-p prescreen:
    each candidate node ordering is first validated over GF(p); only
    screened-in orderings get the exact fmpq Thiele + exact certify, making
    this cheaper than thiele_loo when all orderings fail (measured 2.4x at
    deg(12,14)/n=62, production scale; the shared full-node fallback Thiele is the
    remaining floor). No practical
    false rejects: BOTH SCREEN_PRIMES must cleanly reject to skip an ordering
    (joint spurious-rejection prob ~1e-36) and any mod-p 0-division falls back
    to the exact fit. Semantics identical to thiele_loo: (P, Q, True, deg) on
    success; on all-fail the single full-node exact Thiele with ok=False
    (deg report)."""
    n = len(xs)
    test_idx, fit_idx = _loo_split(n, n_loo)
    check = test_idx + fit_idx
    for od in [fit_idx, list(reversed(fit_idx)), fit_idx[0::2] + fit_idx[1::2]]:
        screened = False
        for p in SCREEN_PRIMES:
            try:
                if _thiele_modp_ok(xs, ys, od, check, p):
                    screened = True; break
            except ZeroDivisionError:
                screened = True; break   # ambiguous mod p -> certify exactly
        if not screened:
            continue
        P, Q = thiele([xs[i] for i in od], [ys[i] for i in od])
        ok = True
        for r in check:
            xr = fq(xs[r]); Qv = Q(xr)
            if Qv == 0 or P(xr) / Qv != fq(ys[r]):
                ok = False; break
        if ok:
            return P, Q, True, max(P.degree(), Q.degree())
    P, Q = thiele(xs, ys)
    return P, Q, False, max(P.degree(), Q.degree())


def extract_stairs(qex, cap=12):
    """Per-row stair-ratio library R = Q_{m+1}/Q_m (origin Rlib
    implementation). qex: {(i,j): {m: [coeff_strs]}}. For every
    consecutive-m pair of validated denominators, keep the ratio iff the
    division is exact (r*Q_m == Q_{m+1} verified by exact_div) and deg r >= 1;
    dedupe by str(coeffs); cap `cap`=12 distinct ratios per row.
    Returns {row_i: [fmpq_poly, ...]}."""
    Rlib = {}
    for (i, j), md in qex.items():
        qp = {int(m): qpoly(cs) for m, cs in md.items()}
        ms = sorted(qp)
        for a, b in zip(ms, ms[1:]):
            if b != a + 1:
                continue
            r = exact_div(qp[b], qp[a])
            if r is None or r.degree() < 1:
                continue
            lst = Rlib.setdefault(i, [])
            if len(lst) < cap and str(r.coeffs()) not in {str(x.coeffs()) for x in lst}:
                lst.append(r)
    return Rlib


def row_lcm(qex):
    """Per-row LCM of all validated denominator polys, lcm via a*b//gcd
    (the 'C_rowLCM' candidate family, proven 31/31 on
    production sub rows). Returns {row_i: fmpq_poly}."""
    out = {}
    for (i, j), md in qex.items():
        for m, cs in md.items():
            q = qpoly(cs)
            out[i] = poly_lcm(out[i], q) if i in out else q
    return out


def degree_budget(qex, n, n_loo=6, project_extra=4):
    """Sample-budget analytics: 'how many samples do we
    need' as a one-liner. With n path nodes and n_loo held out, a blind Thiele
    resolves deg P + deg Q + 1 <= n - n_loo, and a known-Q Newton certify
    needs deg(y*Q) <= n - n_loo - 2 (probe_knownQ4 cap = n_good - 8 at
    n_loo=6). Returns per-row dict:
      blind_ceiling    (n - n_loo)//2 — max per-side degree a blind fit resolves
      certify_cap      n - n_loo - 2 — max candidate-Q degree worth screening
      stair_step_degs  sorted distinct stair-ratio degrees for the row
      projected_degQ   {m: deg} — observed max deg(Q_m) at validated m, then
                       projected project_extra more steps at +max(stair deg)
    A future order m is certifiable iff projected_degQ[m] <= certify_cap, and
    blind-findable iff projected_degQ[m] <= blind_ceiling."""
    Rlib = extract_stairs(qex)
    rows = {}
    for (i, j), md in qex.items():
        dm = rows.setdefault(i, {})
        for m, cs in md.items():
            m = int(m)
            dm[m] = max(dm.get(m, 0), qpoly(cs).degree())
    out = {}
    for i, dm in rows.items():
        steps = sorted({r.degree() for r in Rlib.get(i, [])})
        step = max(steps) if steps else 0
        proj = dict(sorted(dm.items()))
        mmax = max(dm)
        for m in range(mmax + 1, mmax + 1 + project_extra):
            proj[m] = proj[m - 1] + step
        out[i] = {"blind_ceiling": (n - n_loo) // 2,
                  "certify_cap": n - n_loo - 2,
                  "stair_step_degs": steps,
                  "projected_degQ": proj}
    return out


def _pfrac(x):
    """Exact rational -> Fraction (int, Fraction, or anything exposing
    .numerator/.denominator, fmpq included)."""
    if isinstance(x, Fraction):
        return x
    if isinstance(x, int):
        return Fraction(x)
    return Fraction(int(x.numerator), int(x.denominator))


def mpoly_eval(poly_coeffs, params):
    """Evaluate a {exponent_tuple: Fraction} multivariate polynomial at a
    param-tuple, Fraction-exact. Companion evaluator for the
    parametric_rec_interp return format."""
    ps = [_pfrac(p) for p in (params if isinstance(params, (tuple, list))
                              else (params,))]
    total = Fraction(0)
    for exps, c in poly_coeffs.items():
        term = _pfrac(c)
        for p, e in zip(ps, exps):
            if e:
                term *= p ** e
        total += term
    return total


def _tensor_lift(axes, values):
    """Iterated-Newton lift of exact leaf values on a tensor-product grid to a
    multivariate polynomial {exponent_tuple: fmpq}. axes: per-axis fmpq node
    lists; values: nested lists matching the axes shape, fmpq leaves.
    Interpolation is linear in the data, so lifting axis 0 coefficient-wise
    over the already-lifted tail axes is exact."""
    if len(axes) == 1:
        P = newton_interp(axes[0], values)
        return {(e,): c for e, c in enumerate(P.coeffs()) if c != 0}
    sub = [_tensor_lift(axes[1:], v) for v in values]
    keys = set()
    for s in sub:
        keys |= set(s)
    out = {}
    zero = flint.fmpq(0)
    for key in sorted(keys):
        P = newton_interp(axes[0], [s.get(key, zero) for s in sub])
        for e, c in enumerate(P.coeffs()):
            if c != 0:
                out[(e,) + key] = c
    return out


def _normalized_rec(rec_at_point, params, normalize_key):
    """Call rec_at_point(params) and normalize so coeffs[normalize_key] == 1.
    Fail-loud on a missing or vanishing normalizer. Returns {key: Fraction}."""
    rec = rec_at_point(params)
    if normalize_key not in rec:
        raise ValueError(
            f"parametric_rec_interp: normalize_key {normalize_key!r} absent "
            f"from rec_at_point({params!r})")
    v0 = _pfrac(rec[normalize_key])
    if v0 == 0:
        raise ValueError(
            f"parametric_rec_interp: coeffs[{normalize_key!r}] == 0 at "
            f"{params!r} — cannot normalize (redesign the grid off this "
            f"degeneracy)")
    return {k: _pfrac(v) / v0 for k, v in rec.items()}


def _fresh_offgrid_tuple(axes, seed_w):
    """One deterministic off-grid param-tuple: per axis with >=2 nodes take
    a0 + w*(a1-a0) for the first weight w (walking seed_w, 1-seed_w, then a
    fixed ladder) that lands off the axis node set; single-node axes keep
    their node."""
    ws = [seed_w, 1 - seed_w, Fraction(1, 5), Fraction(2, 5), Fraction(3, 7),
          Fraction(4, 7), Fraction(5, 11), Fraction(6, 11), Fraction(7, 13)]
    pt = []
    for ax in axes:
        if len(ax) == 1:
            pt.append(ax[0])
            continue
        nodes = set(ax)
        for w in ws:
            cand = ax[0] + w * (ax[1] - ax[0])
            if cand not in nodes:
                pt.append(cand)
                break
        else:  # 9 distinct candidates on >=2 nodes cannot all collide
            raise ValueError(
                "parametric_rec_interp: could not place an off-grid "
                "cross-check node")
    return tuple(pt)


def parametric_rec_interp(rec_at_point, param_grid, normalize_key, deg_bound=12):
    """
    Given rec_at_point(params)->{(j,k):Fraction} (a numeric recurrence at a param-tuple),
    evaluate on param_grid, normalize each so coeffs[normalize_key]==1, then for each
    (j,k) fit a multivariate polynomial in the params via iterated newton_interp.
    Returns {(j,k): poly_coeffs}. Cross-checks on 2 fresh held-out param-tuples.

    param_grid: a COMPLETE tensor-product grid given as a flat list of
    param-tuples (bare scalars accepted for one parameter); the per-axis node
    sets are inferred from the distinct coordinate values and every
    combination must be present — an incomplete grid raises ValueError.
    Coordinates and coefficient values follow the module convention (Fraction,
    int, or anything exposing .numerator/.denominator).

    deg_bound caps the per-axis interpolation degree (len(axis)-1 <= deg_bound,
    else ValueError) — the sample-budget guard, in the degree_budget spirit.

    poly_coeffs is {exponent_tuple: Fraction} with exponents aligned to the
    param-tuple order; evaluate with mpoly_eval. Keys missing at some grid
    points are treated as exact 0 there (structural-zero drop tolerance), but
    normalize_key must be present and nonzero everywhere.

    The cross-check calls rec_at_point at 2 deterministic OFF-GRID
    param-tuples, normalizes, and demands exact equality with the fitted
    polynomials on the union key set (a key appearing fresh at a check node is
    equally fatal). Any mismatch raises ValueError — a degenerate grid slice
    yields a deterministic wrong-but-valid fit (see FOOTGUNS), so the check is
    not optional.
    """
    if not param_grid:
        raise ValueError("parametric_rec_interp: empty param_grid")
    tuples = [tuple(t) if isinstance(t, (tuple, list)) else (t,)
              for t in param_grid]
    ndim = len(tuples[0])
    if any(len(t) != ndim for t in tuples):
        raise ValueError("parametric_rec_interp: mixed param-tuple lengths")
    key_tuples = [tuple(_pfrac(c) for c in t) for t in tuples]
    axes = [sorted({t[i] for t in key_tuples}) for i in range(ndim)]
    for i, ax in enumerate(axes):
        if len(ax) - 1 > deg_bound:
            raise ValueError(
                f"parametric_rec_interp: axis {i} needs degree {len(ax) - 1} "
                f"> deg_bound={deg_bound}")
    n_expect = 1
    for ax in axes:
        n_expect *= len(ax)
    by_node = dict(zip(key_tuples, tuples))
    if len(by_node) != len(tuples) or len(by_node) != n_expect:
        raise ValueError(
            f"parametric_rec_interp: param_grid is not a complete duplicate-"
            f"free tensor-product grid ({len(by_node)} distinct tuples, "
            f"axes want {n_expect})")

    normed = {kt: _normalized_rec(rec_at_point, by_node[kt], normalize_key)
              for kt in key_tuples}
    all_keys = set()
    for rec in normed.values():
        all_keys |= set(rec)

    def nest(prefix, rest, key):
        if not rest:
            return fq(normed[prefix].get(key, Fraction(0)))
        return [nest(prefix + (a,), rest[1:], key) for a in rest[0]]

    axes_fq = [[fq(a) for a in ax] for ax in axes]
    fit = {}
    for key in all_keys:
        lifted = _tensor_lift(axes_fq, nest((), axes, key))
        fit[key] = {exps: _pfrac(c) for exps, c in lifted.items()}

    for seed_w in (Fraction(1, 2), Fraction(1, 3)):
        chk = _fresh_offgrid_tuple(axes, seed_w)
        rec = _normalized_rec(rec_at_point, chk, normalize_key)
        for key in sorted(all_keys | set(rec)):
            want = rec.get(key, Fraction(0))
            got = mpoly_eval(fit.get(key, {}), chk)
            if got != want:
                raise ValueError(
                    f"parametric_rec_interp: held-out cross-check FAILED at "
                    f"{chk!r} for key {key!r}: fit gives {got}, "
                    f"rec_at_point gives {want} — grid under-resolves the "
                    f"parameter dependence (raise the per-axis node counts) "
                    f"or a slice is degenerate")
    return fit


# ===== degeneracy-safe grid designer (pathgrid; tools/pathgrid.py forwards here) =====



def _frac(x):
    """str 'p/q' / Fraction / int -> exact Fraction."""
    return x if isinstance(x, Fraction) else Fraction(x)


def coincidence_loci(PA, PB, vars=None):
    """All t in the OPEN interval (0,1) where two path variables collide:
    solve s_i(t)==s_j(t), i.e. (a_i-a_j) + t*((b_i-a_i)-(b_j-a_j)) = 0, for
    every unordered pair — Fraction-exact, no floats. PA, PB: {var: rational}.
    vars defaults to sorted(PA). Returns {Fraction(t): "si=sj"} (labels
    comma-joined when several pairs collide at the same t). Raises ValueError
    on a WHOLE-PATH coincidence (s_i==s_j at both endpoints): no 1-D grid on
    that path can be non-degenerate, so it must error loudly, not silently."""
    vs = list(vars) if vars is not None else sorted(PA)
    loci = {}
    for a in range(len(vs)):
        for b in range(a + 1, len(vs)):
            vi, vj = vs[a], vs[b]
            ai, aj = _frac(PA[vi]), _frac(PA[vj])
            di, dj = _frac(PB[vi]) - ai, _frac(PB[vj]) - aj
            if di == dj:
                if ai == aj:
                    raise ValueError(
                        f"whole-path coincidence {vi}={vj}: equal at both "
                        f"endpoints — every t on this path is degenerate")
                continue
            t = (aj - ai) / (di - dj)
            if 0 < t < 1:
                lab = f"{vi}={vj}"
                loci[t] = loci[t] + "," + lab if t in loci else lab
    return loci


def design_grid(PA, PB, n_points, denominators=(60, 120, 240), exclude=(),
                poles=(), vars=None):
    """Emit n_points FRESH exact-rational t in (0,1) for the path PA->PB,
    avoiding (1) the path's coincidence loci (coincidence_loci above), (2) a user-supplied pole list `poles` (known
    DE-denominator zeros on the path), (3) any t in `exclude` (existing grids
    — pass every previously farmed t so no node is ever resampled).

    Candidates are k/q for q in denominators, k=1..q-1, deduped as exact
    rationals (30/60 == 1/2 == 120/240 counts once). exclude/poles entries may
    be 'p/q' strings, Fractions, or ints. Raises ValueError LOUDLY if it
    cannot supply n_points — never silently returns a short grid.

    Each point carries provenance:
      {"t": "p/q", "checked": ["no-coincidence", "no-pole", "fresh"]}"""
    loci = coincidence_loci(PA, PB, vars)
    bad_pole = {_frac(x) for x in poles}
    bad_old = {_frac(x) for x in exclude}
    out, used = [], set()
    for q in denominators:
        for k in range(1, q):
            t = Fraction(k, q)
            if t in used or t in loci or t in bad_pole or t in bad_old:
                continue
            used.add(t)
            out.append({"t": f"{t.numerator}/{t.denominator}",
                        "checked": ["no-coincidence", "no-pole", "fresh"]})
            if len(out) == n_points:
                return out
    raise ValueError(
        f"design_grid: only {len(out)}/{n_points} admissible t values from "
        f"denominators={tuple(denominators)} after excluding {len(loci)} "
        f"coincidence loci, {len(bad_pole)} poles, {len(bad_old)} prior nodes "
        f"— widen denominators or lower n_points")

#!/usr/bin/env python3
r"""thiele_gate.py -- exact Thiele continued-fraction family reconstruction +
byte-exact gating of NEW sample points (Fraction-exact throughout).

A(u)-integrity
gate (ACHECK_GATE.json = 88 entries / 10032 byte-exact new-node
checks, 0 fails).

WHAT IT DOES
  Given verified sample points (x_k, y_k) of a quantity believed to be one
  fixed rational function y = R(x), with x_k, y_k exact `fractions.Fraction`s:
  build the Thiele continued fraction with SATURATION detection (the CF stops
  as soon as it already reproduces every remaining fit point -- its depth is
  then an honest degree certificate, not the point count), require byte-exact
  reproduction of all held-out points, then GATE new, unverified points:
  a new point passes iff it lies byte-exact on the reconstructed R.

THE PATTERN THIS ENCODES (config-forcing + external gate; the stage-1 skip)
  Use when each sample point is expensive because a per-point derivation /
  verification stage dominates the farm cost.  In the originating production
  farm, every node re-derived and stage-1-verified its config choice per node.  Cure: FORCE the
  expensive stage's config to the known-good one and skip its per-node
  verification, then restore integrity EXTERNALLY with this gate -- fit the
  family on the stage-1-verified nodes only, demand CF saturation + byte-exact
  held-out reproduction, and byte-exact gate every config-forced node.  A
  forced node with any defect (wrong config branch, corrupted export, wrong
  kinematic slice) has essentially zero chance of landing byte-exact on the
  correct rational function.  Measured on the originating farm: config-forced nodes
  ~9x cheaper than stage-1-verified ones; gate result 88/88 entries, 114 new
  nodes, 10032/10032 byte-exact checks (ACHECK_GATE.json).

  Point ordering: fit points are interleaved from both ends of the sorted x
  list (lo, hi, lo+1, hi-1, ...) -- keeps the CF well-conditioned (depth stays
  near-minimal, inverse-difference degeneracies rare).  On a degeneracy
  (CFDead: zero inverse difference, e.g. an even function on a symmetric grid)
  fit() evicts the colliding point to the held-out set and retries;
  correctness is guaranteed by the byte-exact validation, not by the ordering.

API
  cf = fit(xs, ys, max_fit=None, probe=4, interleave=True, retries=8)
      -> ThieleCF or None.  Fits on the first min(max_fit, n) ordered points;
      the remaining in-sample points are validated byte-exact (None on any
      mismatch or unresolvable degeneracy).
  validate(cf, holdout_pairs) -> bool          (byte-exact, all-or-nothing)
  gate(cf, new_pairs) -> {'ok', 'fail', 'fails': [{x, want, got}]}
  gate_family(banked, new, ...) -> acheck-style report over vector-valued
      samples (banked/new = [(x, [y_0..y_m]), ...]); fits each component,
      reports entries_ok/entries_fail/new_ok/new_fail/cf_depths/fails.

CLI (the common JSON case)
  python3 thiele_gate.py IN.json [-o REPORT.json]
  IN.json: {"banked": [[x, y-or-ylist], ...], "new": [[x, y-or-ylist], ...],
            "max_fit": 140, "labels": [...]}   (x, y as Fraction strings)
  Exit 1 on any entry-fail or new-point fail.
"""
import json
import sys
from fractions import Fraction

DEFAULT_PROBE = 4      # early-stop lookahead before full saturation check


class CFDead(Exception):
    """Zero inverse difference / zero CF tail -- ordering degeneracy.
    idx = offending point's position in the fit ordering (construction-time
    divided-difference collision only), else None."""

    def __init__(self, idx=None):
        super().__init__(idx)
        self.idx = idx


def _poly_add(p, q):
    n = max(len(p), len(q))
    p = p + [0] * (n - len(p))
    q = q + [0] * (n - len(q))
    return [a + b for a, b in zip(p, q)]


def _poly_eval(p, x):
    v = 0
    for c in reversed(p):
        v = v * x + c
    return v


class ThieleCF:
    """Thiele continued fraction a_0 + (x-u_0)/(a_1 + (x-u_1)/(...)).

    Evaluation: fast backward recursion; if that hits a zero tail (a removable
    artifact of the CF representation, not a true pole), falls back to the
    exact convergent polynomials P(x)/Q(x).  CFDead is only raised at a true
    pole (Q(x) == 0)."""
    __slots__ = ('a', 'used', '_pq')

    def __init__(self, a, used):
        self.a = list(a)
        self.used = list(used)
        self._pq = None

    @property
    def depth(self):
        return len(self.a)

    def _poly_pair(self):
        """Convergent polynomials (P, Q) with CF == P(x)/Q(x), exact."""
        if self._pq is None:
            Pm, Qm = [1], [0]
            P, Q = [self.a[0]], [1]
            for k in range(1, len(self.a)):
                u, ak = self.used[k - 1], self.a[k]
                # new = ak*cur + (x-u)*prev
                P, Pm = _poly_add([ak * c for c in P],
                                  _poly_add([0] + Pm, [-u * c for c in Pm])), P
                Q, Qm = _poly_add([ak * c for c in Q],
                                  _poly_add([0] + Qm, [-u * c for c in Qm])), Q
            self._pq = (P, Q)
        return self._pq

    def __call__(self, x):
        try:
            val = self.a[-1]
            for j in range(len(self.a) - 2, -1, -1):
                if val == 0:
                    raise CFDead()
                val = self.a[j] + (x - self.used[j]) / val
            return val
        except CFDead:
            P, Q = self._poly_pair()
            qv = _poly_eval(Q, x)
            if qv == 0:
                raise CFDead()
            return _poly_eval(P, x) / qv


def _interleave_order(n):
    order = []
    lo, hi = 0, n - 1
    while lo <= hi:
        order.append(lo)
        lo += 1
        if lo <= hi:
            order.append(hi)
            hi -= 1
    return order


def _thiele_core(xs, ys, probe):
    """Exact Thiele CF with saturation detection over (xs, ys) in this order.
    Returns (cf, verified_all): verified_all=True means cf was checked
    byte-exact on EVERY input point (saturated early); False means the CF
    consumed all points and the caller must validate it (incl. on xs itself:
    a degenerate CF is not guaranteed to interpolate its own support).
    Raises CFDead (with .idx) on a zero inverse difference."""
    n = len(xs)
    if all(y == ys[0] for y in ys):
        return ThieleCF([ys[0]], [xs[0]]), True
    g = list(ys)
    a, used = [], []
    for k in range(n):
        a.append(g[k])
        used.append(xs[k])
        if k + 1 < n:
            cf = ThieleCF(a, used)
            try:
                if all(cf(xs[t]) == ys[t]
                       for t in range(k + 1, min(k + 1 + probe, n))) \
                   and all(cf(xs[t]) == ys[t] for t in range(k + 1, n)) \
                   and all(cf(xs[t]) == ys[t] for t in range(k + 1)):
                    return cf, True
            except CFDead:
                pass
            new = []
            for i in range(k + 1, n):
                d = g[i] - a[k]
                if d == 0:
                    raise CFDead(idx=i)
                new.append((xs[i] - xs[k]) / d)
            g = g[:k + 1] + new
    return ThieleCF(a, used), False   # caller must validate


def fit(xs, ys, max_fit=None, probe=DEFAULT_PROBE, interleave=True, retries=8):
    """Fit + in-sample held-out validation.  Returns ThieleCF or None.

    xs, ys: equal-length sequences of Fractions (any exact type with ==, -, /).
    max_fit: fit on the first min(max_fit, n) ordered points; the rest are
             byte-exact validation points (None if any mismatches).
    interleave: sort by x and interleave from both ends (recommended).
    retries: on a CFDead degeneracy (zero inverse difference, e.g. an even
             function on a symmetric grid) the colliding point is EVICTED to
             the held-out set and the fit retried, up to `retries` times --
             safe because correctness comes from the byte-exact validation.
    """
    if len(xs) != len(ys):
        raise ValueError('xs/ys length mismatch')
    if len(xs) == 0:
        raise ValueError('no points')
    if len(set(xs)) != len(xs):
        raise ValueError('duplicate x values')
    pairs = list(zip(xs, ys))
    if interleave:
        pairs.sort(key=lambda p: p[0])
        pairs = [pairs[i] for i in _interleave_order(len(pairs))]
    nfit = len(pairs) if max_fit is None else min(max_fit, len(pairs))
    fit_pairs, rest = pairs[:nfit], pairs[nfit:]
    for _ in range(retries + 1):
        try:
            cf, verified = _thiele_core([p[0] for p in fit_pairs],
                                        [p[1] for p in fit_pairs], probe)
        except CFDead as e:
            if e.idx is not None and len(fit_pairs) > 2:
                rest = rest + [fit_pairs.pop(e.idx)]
                continue
            return None
        try:
            held = all(cf(x) == y for x, y in rest)
            if held and not verified:      # full-depth CF: self-check support
                held = all(cf(x) == y for x, y in fit_pairs)
        except CFDead:
            return None
        return cf if held else None
    return None


def validate(cf, holdout_pairs):
    """True iff cf reproduces every (x, y) byte-exact."""
    try:
        return all(cf(x) == y for x, y in holdout_pairs)
    except CFDead:
        return False


def gate(cf, new_pairs):
    """Byte-exact gate of NEW points against a validated cf."""
    rep = {'ok': 0, 'fail': 0, 'fails': []}
    for x, y in new_pairs:
        try:
            got = cf(x)
            ok = (got == y)
            got_s = str(got)
        except CFDead:
            ok, got_s = False, 'CFDead'
        if ok:
            rep['ok'] += 1
        else:
            rep['fail'] += 1
            rep['fails'].append({'x': str(x), 'want': str(y), 'got': got_s})
    return rep


def _fit_gate_component(payload):
    """Worker for gate_family (module-level so it pickles)."""
    xs, ys, newpairs, label, max_fit, probe, retries = payload
    cf = fit(xs, ys, max_fit=max_fit, probe=probe, retries=retries)
    if cf is None:
        return (label, None, None)
    return (label, cf.depth, gate(cf, newpairs))


def gate_family(banked, new, max_fit=140, probe=DEFAULT_PROBE, labels=None,
                retries=4, progress=None, workers=1):
    """acheck-style report over vector-valued samples.

    banked, new: lists of (x, [y_0..y_{m-1}]) with exact Fractions; every
    banked row must be trusted (stage-1-verified).  Fits each component,
    validates in-sample held-out points beyond max_fit, gates all new rows.
    Report keys mirror ACHECK_GATE.json: entries_ok / entries_fail / new_ok /
    new_fail / cf_depths / fails.  workers>1 parallelizes over components.
    """
    if not banked:
        raise ValueError('no banked rows')
    ncomp = len(banked[0][1])
    if labels is None:
        labels = [str(i) for i in range(ncomp)]
    rep = {'entries_ok': 0, 'entries_fail': 0, 'new_ok': 0, 'new_fail': 0,
           'cf_depths': {}, 'fails': []}
    xs = [row[0] for row in banked]
    payloads = [(xs, [row[1][j] for row in banked],
                 [(x, yv[j]) for x, yv in new],
                 labels[j], max_fit, probe, retries) for j in range(ncomp)]
    if workers > 1:
        from multiprocessing import Pool
        with Pool(workers) as pool:
            results = pool.map(_fit_gate_component, payloads)
    else:
        results = map(_fit_gate_component, payloads)
    for j, (label, depth, g) in enumerate(results):
        if depth is None:
            rep['entries_fail'] += 1
            rep['fails'].append({'kind': 'entry_no_family', 'entry': label})
            continue
        rep['entries_ok'] += 1
        rep['cf_depths'][label] = depth
        rep['new_ok'] += g['ok']
        rep['new_fail'] += g['fail']
        for f in g['fails']:
            rep['fails'].append({'kind': 'new_mismatch', 'entry': label,
                                 'x': f['x']})
        if progress:
            progress(j, ncomp, label, depth)
    return rep


# ----------------------------------------------------------------------- CLI
def _parse_rows(rows):
    out = []
    for x, y in rows:
        yv = [Fraction(v) for v in y] if isinstance(y, list) \
            else [Fraction(y)]
        out.append((Fraction(x), yv))
    return out


def main(argv):
    import argparse
    import time
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('input', help='JSON with banked/new rows (see docstring)')
    ap.add_argument('-o', '--output', default=None, help='report JSON path')
    args = ap.parse_args(argv)
    t0 = time.time()
    spec = json.load(open(args.input))
    banked = _parse_rows(spec['banked'])
    new = _parse_rows(spec.get('new', []))
    rep = gate_family(
        banked, new, max_fit=spec.get('max_fit', 140),
        probe=spec.get('probe', DEFAULT_PROBE), labels=spec.get('labels'),
        workers=spec.get('workers', 1),
        progress=lambda j, n, lab, d: print(
            f'entry {lab} ({j + 1}/{n}) depth {d} [{time.time() - t0:.0f}s]',
            flush=True))
    if args.output:
        json.dump(rep, open(args.output, 'w'), indent=1)
    print(f"entries: {rep['entries_ok']} ok / {rep['entries_fail']} no-family; "
          f"new-point checks: {rep['new_ok']} ok / {rep['new_fail']} FAIL "
          f"[{time.time() - t0:.0f}s]")
    return 1 if (rep['entries_fail'] or rep['new_fail']) else 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))

#!/usr/bin/env python3
"""verify_eras_adversarial.py — consumer-pluggable adversarial verification
battery for eras adaptations.

LAW: rerun THIS verifier against ANY adaptation of eras before that
adaptation is used, and CHECK the returned result. The battery is generic
over the protocol below and vector-p capable (length-1 vectors ARE the
univariate case; nothing is special-cased — the reference univariate arms
run as the built-in `reference` battery at full spec). The battery was
hardened against independent adversarial attack suites — adaptive
observed-hull enclosures, boundary-band-only corruptions, vacuous
enclosures laundered through coarse point balls, scale-dependent plant
self-destruction — and the laws below carry the measured outcomes.

THREAT MODEL (stated, binding on claims)
----------------------------------------
The battery defends against HONEST-BUT-BUGGY adaptations: wrong remainder
terms, dropped dimensions, narrowed/shifted enclosures, stale caches,
precision loss, boundary-band mistakes, and ADAPTIVE bugs (enclosures
accidentally derived from observed point queries — the battery commits
every enclosure before it evaluates any point, and adds a second fresh
post-enclose sample stream as belt and braces). It does NOT defend against
a consumer who maliciously edits the verifier itself or fakes its report;
no self-run harness can. It also cannot certify what finite sampling
cannot see: corruptions confined to sets the samples never hit are caught
only via the certified-envelope channel (B) or not at all — stated
honestly below, per arm.

CONSUMER PROTOCOL
-----------------
An *adaptation* is any object with:

  label : str                    [attribute] battery-report name.

  enclose(p0, radii, K) -> enclosure
      p0, radii: equal-length sequences (strings like "1/3"/"1e-4", or
      floats) — the box is p0 +/- radii per dimension. K: enclosure order;
      0 = naive interval evaluation, 1 = mean-value form, K >= 2 = order-K
      Taylor-Lagrange. An adaptation may cap K (returning its best rigorous
      form for the requested box); every returned enclosure MUST be
      rigorous for the CLOSED box (boundary included — the battery samples
      exact-face and exact-corner points AT p0 +/- radii). Returns one
      enclosure per output: an arb/acb ball, or a (lo, hi) TUPLE (outward
      hull is taken), or a LIST of these for multi-output adaptations
      (tuple = pair, list = outputs).

  point_eval(p, prec=None) -> value(s)
      High-precision evaluation at point p (sequence, len(p0)).
      QUANTITATIVE TIGHTNESS REQUIREMENT — 'tighter than any enclosure
      width under test' is NOT sufficient (measured: a point_eval with
      ball radius ~width/2000 inflates the sampled spread enough to
      convert a 330,000x-vacuous enclosure into informative PASS
      coverage): the max point-ball radius must be <= enclosure_width /
      (vacuous_factor * point_margin) for every cell you want judged
      (defaults: width / 1e4). Coarser point balls make the cell
      POINT-PRECISION INDET ('point_eval too coarse to judge vacuousness
      at this width') — never informative, never a pass, and plants
      cannot be verified from such cells. Honor the optional prec hint
      (the battery passes prec=fd_prec for derivative checks). Same
      output arity as enclose.

  deriv(p, dim) -> value(s)      [OPTIONAL]
      Partial derivative along dimension dim at point p, evaluated at
      >= fd_prec. If exposed, the battery checks it per dim against
      central finite differences.

  degraded() -> adaptation       [OPTIONAL]
      The same adaptation with its internal constants rebuilt at 53-bit
      precision. If exposed, the battery runs the PRECISION PROBE and the
      degradation MUST be caught (the documented pathology class: 53-bit
      s-constants degrade the reference unc K=16 r=1e-4 enclosure from
      relwidth 4.7e-3 to ~2e9 — measured, and regression-gated by the
      reference battery every run).

THE BATTERY (generic over the protocol)
---------------------------------------
(o)   ORDERING LAW — ALL enclosures for ALL shells and orders are
      collected BEFORE the battery evaluates a single point. An
      adaptation therefore commits its claims before it can observe any
      battery query. A second, fresh-stream random batch is sampled
      after collection as well (belt + braces). This closes the
      adaptive-hull hole: an adaptation that returns the padded hull of
      previously observed point values is caught at every seed.
(i)   CONTAINMENT, CLOSED-BOX — per radius shell: exact corners AT
      p0 +/- radii (all 2^dim for dim <= 6; for dim >= 7 a documented
      random-sign-vector corner sample of 2*dim + 64 draws — exact +/-r
      coordinates, sign patterns randomized, NOT all 2^dim), exact face
      points (dim >= 2: per dimension one point on each face, coordinate
      AT +/-r), stratified NEAR-FACE BAND points (per dimension, per
      sign, 4 points whose that-dim coordinate is drawn uniformly in
      +/-[0.9995, 1)*r, other dims uniform over the closed box — exact
      boundary points close only the band's endpoints, so the open band
      needs its own samples), plus >= --points (default 100) fresh
      random points drawn over the CLOSED box (no interior clamp), plus
      a second fresh closed-box batch of the same size on a ":post" seed
      stream. Every point value must lie inside every FINITE informative
      enclosure. ZERO violations tolerated. THIN-SET RESIDUAL (honest,
      measured): corruptions confined to thin sets are detected with
      probability scaling in points x seeds, NOT deterministically.
      Measured catch rates on the attack suite: an enclosure rigorous
      everywhere except the outer near-face band (escape window ~2e-4*r)
      is caught at 6/6 seeds via the stratified band points; a bump on a
      1e-4*r-thin face sliver over ~5% of one face (2D) at 9/100 fresh
      seeds — face-interior slivers remain a points-x-seeds detection,
      stated rather than pretended away; an interior thin-dip
      stream-prediction overfit at 2/12 fresh seeds (an overfit tuned to
      a specific sampler stream remains the documented single-seed
      residual). A single-seed pass is NECESSARY, NOT SUFFICIENT: rerun
      with several --seed values before relying on an adaptation
      (corruptions confined to unsampled interior sets are otherwise
      only visible to channel B).
(ii)  NON-FINITE = NO-INFORMATION — a nan/inf enclosure makes no claim;
      recorded as an INDET cell and NEVER counted as coverage.
(iii) VACUOUSNESS GATE — per box/K/output cell, the enclosure width is
      compared to the spread of the sampled value MIDPOINTS (midpoints
      ONLY — a spread computed from point-ball hulls including radii
      lets a coarse point_eval inflate the spread and launder a
      330,000x-vacuous enclosure into informative PASS coverage;
      measured): width > --vacuous-factor (default 1e3) * midpoint
      spread => the cell is VACUOUS = INDET (no-information; containment
      inside it proves nothing), reported separately and NEVER counted
      as coverage.
      POINT-PRECISION GATE (runs first): if the max point-ball radius
      > enclosure_width / (vacuous_factor * --point-margin, default 10)
      the cell is POINT-PRECISION INDET ('point_eval too coarse to judge
      vacuousness at this width') — no containment claim, no vacuousness
      claim, never informative. Exception: a point ball PROVABLY
      DISJOINT from the enclosure still proves a containment violation
      (FATAL) at any coarseness — an honest point ball contains the true
      value, so disjointness certifies escape; this cannot false-FAIL an
      honest adaptation. A fat-ball vacuousness launder therefore lands
      INDET (caught — no informative coverage, result.ok False), and its
      HONEST twin (honest enclosure + coarse point_eval) lands the same
      POINT-PRECISION INDET rather than a false FAIL. Callers may
      DECLARE expected-no-information cells (expect_indet={(shell_idx,
      K), ...}) — e.g. the reference battery's dependency-blowup
      regression rows, whose widths are separately regression-gated
      against the recorded reference values. Declared cells are still
      containment-checked when finite (a violation is FATAL regardless)
      and still reported INDET; they are excluded ONLY from the majority
      arithmetic below. A declaration can never create coverage.
(iv)  PLANTED ERRORS MUST FIRE — the battery wraps the adaptation in
      planted variants: enclosure narrowed (width x (1-1e-6) and x 0.5),
      center shifted (by width x 1e-3 and by a SCALE-ADAPTIVE absolute
      shift, max(1e-9, rad*2^-28, |mid|*2^-52), plant arithmetic pinned
      at 512 bits — a fixed 1e-9 shift is swallowed by 30-bit mag radius
      round-up whenever the enclosure radius exceeds ~0.5, and a harness
      that then hard-FAILs the HONEST adaptation as "HARNESS BLIND" is
      itself broken; regression: the ScaleRegression arm below,
      enclosure mid ~192 / rad ~21.5, must PASS), and one dimension's
      radius silently ignored (vector case; each dim in turn). Each
      plant must be flagged through at least one channel:
        (A) point containment on the fresh sample (exact corners/faces
            make this bite for gross plants on tight enclosures);
        (B) certified-envelope cross-check — the honest enclosure, once
            certified by (i) on an informative cell, must be contained in
            any claimed enclosure of the SAME box by the SAME adaptation;
            a planted enclosure that has lost certified content fails.
      Channel (B) exists because finite sampling cannot resolve sub-width
      corruptions (the 1e-6 narrowing); that is stated honestly rather
      than pretending 100 points can see a 1e-6 sliver. A non-finite
      planted enclosure counts as flagged (INDET can never pass). If a
      plant does not fire, the battery SPLITS the verdict: planted
      enclosure IDENTICAL to the honest one => "plant UNRESOLVABLE at
      this scale" (INDET for that plant, scale printed — the
      perturbation collapsed into the ball representation); planted
      enclosure DIFFERENT yet unflagged => THE HARNESS IS BLIND = hard
      error (FAIL). Plants anchor on an INFORMATIVE certified cell — an
      adaptation with only vacuous/non-finite cells at the plant shell
      cannot have its plants verified and FAILs (INDET coverage is not a
      pass, and plants unverifiable = no stamp). Exception: when the
      missing anchor is attributable to the POINT-PRECISION gate (coarse
      point_eval), plants-unverifiable is an INDET note, not a FATAL —
      the honest twin of a fat-ball launder must not become a
      false-rejection class; the run still cannot PASS (no informative
      coverage). Multi-output plants: the fire rule is ANY-output-fires,
      and the report prints per-output channel flags ([per-output flags:
      out0=A+B,out1=-,...]) so masking is visible; the ShiftPlant prints
      the applied shift PER OUTPUT.
(v)   FD DERIVATIVE AGREEMENT — per dim, central differences
      (h = 2^-fd_h_bits, prec fd_prec) vs deriv(); tolerance --fd-tol,
      default 1.1e-12 (the univariate spec class). The gate adds the FD
      and deriv ball radii OUTWARD ((|mid_fd - mid_deriv| + rad_fd +
      rad_deriv) / (|mid_deriv| - rad_deriv)); a derivative ball
      straddling zero at this scale is an INDET note, not a silent
      midpoint-only pass.
(vi)  PRECISION PROBE — where degraded() is exposed: the degraded
      enclosure must be caught — non-finite, or relwidth inflated by
      >= --degrade-factor (default 1e3) over the honest one.

VERDICT / EXIT CLASSES (class 2 is genuinely reachable)
-------------------------------------------------------
  0 = PASS   — no fatals, >= 1 informative coverage cell, and undeclared
               no-information cells do NOT outnumber informative cells.
  1 = FAIL   — any containment violation, harness-blind plant, FD
               mismatch, probe miss, or unverifiable plants.
  2 = INDET  — no fatal, but the run produced no informative coverage at
               all, OR its coverage is majority no-information (undeclared
               INDET cells > informative PASS cells). A majority-vacuous
               run can never be stamped PASS.
Importing consumers: run_battery(...) RETURNS a BatteryResult; check
result.ok (bool(result) == result.ok) — exiting 0 without checking it is
a consumer bug. result.ok = (no new fatals) and (>= 1 informative cell)
and (undeclared INDET cells <= informative cells).

BUILT-IN BATTERIES
------------------
  reference — the reference univariate arms (unc-collapsed@192,
      corr-4dim@192, corr-4dim@48) through the protocol, at full spec:
      radii {1e-4,1e-3,1e-2,5e-2,1e-1}; orders K in {0(naive),1(mv),
      2,4,8,16} for the @192 arms ({0,1} for corr@48, mv only); fresh
      points per shell as in (i), zero violations; FD gate 1.1e-12
      class; mid-exactness and K=1==mv probes; the 0001:+1
      planted-count error at anchor / fiber-equality / enclosure level;
      the generic plants of (iv); precision probe on the unc arm ONLY —
      the pathology lives in the collapsed-polynomial path; the corr
      arm's spec-radius enclosures are dependency-dominated (measured: a
      53-bit rebuild moves corr relwidth by < 1e-15 relative at K in
      {1,4,8}, invisible under its 5.8e-2 true range, and the corr@48
      arm's working precision is itself below 53 bits), so a probe there
      is no-information by construction. COVERAGE HONESTY: most
      full-spec cells are the recorded dependency-blowup rows — the
      documented pathology this tool exists to tame. Those cells are
      DECLARED expected-no-information (expect_indet) and their widths
      are regression-gated against the recorded numbers (naive @1e-4
      relwidth must reproduce the 8.84e22 class; unc K=16 @1e-4 must
      reproduce the 4.7e-3 class). The informative coverage of the
      reference battery is the working-configuration cells (unc K=16
      @1e-4; corr small-radius cells) — that is what "PASS" covers, and
      the report prints the split.
  selftest — MULTIVARIATE SELF-TEST ARM (TEST SCAFFOLDING, NOT SUBSTRATE):
      a small honest 3-parameter product integrand with rigorous naive
      interval and multivariate mean-value enclosures (flint/arb). Its
      only job is proving the harness end-to-end in vector mode: the
      correct adaptation passes, all plants fire (incl. dropped-dim, each
      dim), per-dim FD agrees. It is NOT a shipped enclosure form.
      PLUS the ScaleRegression arm: an HONEST univariate adaptation at
      enclosure scale mid ~192 / rad ~21.5 (the scale where a fixed
      absolute shift plant collapses into the ball representation) —
      must PASS with every plant firing or scale-INDET, never FAIL.

USAGE
  python3 verify_eras_adversarial.py [all|reference|selftest]
      [--seed N] [--points N] [--fd-tol X] [--degrade-factor X]
      [--vacuous-factor X] [--fd-h-bits N]
Consumers: import this module, implement the protocol, call
run_battery(adaptation, shells, Ks, cfg, ...) and CHECK THE RESULT;
MultivariateSelfTest at the bottom is a worked vector-mode example.
"""
import argparse, os, random, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)   # eras.py + QUARTET_COLLAPSE.json ship beside this file
os.environ.setdefault("ERAS_RESULTS", os.path.join(
    os.environ.get("TMPDIR", "/tmp"), "eras_verify_results"))
import eras as tp
from flint import arb, acb, ctx

FATALS, NOTES = [], []
N_PASS = [0]
CELLS_PASS = [0]      # informative coverage cells (box x K x output)
CELLS_INDET = [0]     # undeclared no-information cells (non-finite/vacuous)
CELLS_DECLARED = [0]  # declared-expected no-information cells

def fatal(msg):
    FATALS.append(msg)
    print("FATAL " + msg, flush=True)

def note(msg):
    NOTES.append(msg)

def ok(msg):
    N_PASS[0] += 1
    print("OK    " + msg, flush=True)

class Cfg:
    def __init__(self, seed=987654321202607, points=100, fd_tol="1.1e-12",
                 degrade_factor=1e3, point_prec=256, fd_prec=512,
                 fd_h_bits=40, vacuous_factor=1e3, plant_prec=512,
                 point_margin=10):
        self.seed, self.points, self.fd_tol = seed, points, fd_tol
        self.degrade_factor, self.point_prec, self.fd_prec = \
            degrade_factor, point_prec, fd_prec
        self.fd_h_bits = fd_h_bits
        self.vacuous_factor = vacuous_factor
        self.plant_prec = plant_prec
        self.point_margin = point_margin  # point-precision gate
        # margin: a cell is judgeable only if the max point-ball radius
        # <= enclosure_width / (vacuous_factor * point_margin)

class BatteryResult:
    """Returned by run_battery. bool(result) == result.ok.
    ok = no new fatals AND >= 1 informative coverage cell AND undeclared
    no-information cells do not outnumber informative cells. Importing
    consumers MUST check this (exiting 0 on a False result is a consumer
    bug — the old battery returned None and made that failure silent)."""
    def __init__(self, label, fatals, notes, cells_pass, cells_indet,
                 cells_declared):
        self.label = label
        self.fatals, self.notes = fatals, notes
        self.cells_pass, self.cells_indet = cells_pass, cells_indet
        self.cells_declared = cells_declared
        self.ok = (not fatals) and cells_pass > 0 \
            and not (cells_indet > cells_pass)

    def __bool__(self):
        return self.ok

    def __repr__(self):
        return (f"BatteryResult({self.label!r}, ok={self.ok}, "
                f"fatals={len(self.fatals)}, cells_pass={self.cells_pass}, "
                f"cells_indet={self.cells_indet}, "
                f"cells_declared_indet={self.cells_declared})")

# ------------------------------------------------- enclosure normalization
def _norm_one(o):
    """(lo, hi) tuple -> outward arb hull; balls pass through."""
    if isinstance(o, tuple) and len(o) == 2:
        lo = o[0] if isinstance(o[0], (arb, acb)) else arb(str(o[0]))
        hi = o[1] if isinstance(o[1], (arb, acb)) else arb(str(o[1]))
        return lo.union(hi)
    return o

def norm_out(x):
    """Normalize enclose()/point_eval() output to a list (tuple=pair,
    list=multiple outputs — see protocol docstring)."""
    if isinstance(x, list):
        return [_norm_one(o) for o in x]
    return [_norm_one(x)]

def finite(b):
    return bool(b.is_finite())

def _as_arb(b):
    return b.real if isinstance(b, acb) else b

def relwidth(b):
    b = _as_arb(b)
    if b.mid() == 0:
        return float("inf")
    return float(arb(b.rad()) / abs(arb(b.mid())))

def mid_exact(e):
    if isinstance(e, acb):
        return acb(arb(e.real.mid()), arb(e.imag.mid()))
    return arb(e.mid())

def rad_arb(e):
    return arb(e.real.rad()) if isinstance(e, acb) else arb(e.rad())

def _ball_same(a, b):
    """Exact ball identity (mid and rad both equal) — used to split
    'plant unresolvable at this scale' from 'harness blind'."""
    if isinstance(a, acb) or isinstance(b, acb):
        a, b = acb(a), acb(b)
        return _ball_same(a.real, b.real) and _ball_same(a.imag, b.imag)
    return bool(arb(a.mid()) == arb(b.mid())) \
        and bool(arb(a.rad()) == arb(b.rad()))

def fmt(e):
    try:
        return f"relw {relwidth(e):.2e}"
    except Exception:
        return "relw n/a"

# ----------------------------------------------------------- box sampling
BAND_PTS = 4   # near-face band points per dim per sign

def sample_box(p0, radii, n, rng, prec=256, boundary=True):
    """CLOSED-box sampling. Interior points are drawn over the CLOSED box
    (no clamp — a sampler that clamps its interior points leaves the open
    near-face band deterministically unsampled at every seed, and exact
    boundary points close only that band's endpoints), plus a stratified
    near-face band sample.
    boundary=True adds deterministic exact-boundary points:
      - EXACT CORNERS at p0 +/- radii (offset = the exact double
        float(radii[d])): all 2^dim for dim <= 6; for dim >= 7 a
        random-sign-vector corner sample of 2*dim + 64 draws (exact +/-r
        coordinates, sign patterns from rng — NOT all 2^dim; a
        documented truncation, never a silent one).
      - EXACT FACE POINTS (dim >= 2): per dimension one point on each of
        the two faces — coordinate d exactly at +/-r, other coordinates
        uniform over the closed box. (dim = 1: the corners ARE the faces.)
      - STRATIFIED NEAR-FACE BAND POINTS: per dimension, per sign,
        BAND_PTS (= 4) points whose dim-d coordinate is drawn uniformly
        in +/-[0.9995, 1)*r (the band a clamped sampler never hits; the
        exact corners/faces pin the u = 1 endpoint itself), other
        coordinates uniform over the closed box.
    Plus n fresh random points (each coordinate uniform in +/-r — CLOSED
    box, no clamp). All points are asserted inside the per-dim balls; the
    enclosure under test must be rigorous for the CLOSED box."""
    old = ctx.prec
    ctx.prec = prec
    try:
        dim = len(p0)
        mids = [arb(tp.frac(str(c)).mid()) for c in p0]
        exact_r = [arb(float(radii[d])) for d in range(dim)]  # exact doubles
        balls = [arb(mids[d].mid(), float(radii[d])) for d in range(dim)]
        pts = []
        if boundary:
            if dim <= 6:
                for mask in range(2 ** dim):
                    pts.append([mids[d] + (1 if (mask >> d) & 1 else -1)
                                * exact_r[d] for d in range(dim)])
            else:
                for _ in range(2 * dim + 64):
                    pts.append([mids[d] + (1 if rng.random() < 0.5 else -1)
                                * exact_r[d] for d in range(dim)])
            if dim >= 2:
                for d in range(dim):
                    for sgn in (-1, 1):
                        pt = []
                        for j in range(dim):
                            if j == d:
                                pt.append(mids[j] + sgn * exact_r[j])
                            else:
                                u = arb(str(rng.uniform(-1.0, 1.0)))
                                pt.append(mids[j] + u * arb(str(radii[j])))
                        pts.append(pt)
            for d in range(dim):           # near-face band points
                for sgn in (-1, 1):
                    for _ in range(BAND_PTS):
                        pt = []
                        for j in range(dim):
                            if j == d:
                                u = arb(str(rng.uniform(0.9995, 1.0)))
                                pt.append(mids[j] + sgn * u
                                          * arb(str(radii[j])))
                            else:
                                u = arb(str(rng.uniform(-1.0, 1.0)))
                                pt.append(mids[j] + u * arb(str(radii[j])))
                        pts.append(pt)
        for _ in range(n):
            pt = []
            for d in range(dim):
                u = arb(str(rng.uniform(-1.0, 1.0)))
                pt.append(mids[d] + u * arb(str(radii[d])))
            pts.append(pt)
        for pt in pts:
            for d in range(dim):
                assert balls[d].contains(pt[d]), "sample point left the box"
        return pts
    finally:
        ctx.prec = old

# ------------------------------------------------------------- the plants
PLANT_PREC = 512   # plant arithmetic pinned high, never at ambient
                   # (an ambient-53 path can lose the shift)

class Plant:
    """Wrapper corrupting inner.enclose(); point_eval passes through.
    Corruption arithmetic runs at pinned PLANT_PREC."""
    def __init__(self, inner, tag):
        self.inner, self.tag = inner, tag
        self.label = f"{inner.label}+PLANT[{tag}]"

    def point_eval(self, p, prec=None):
        return self.inner.point_eval(p, prec=prec)

    def enclose(self, p0, radii, K):
        inner_out = norm_out(self.inner.enclose(p0, radii, K))
        old = ctx.prec
        ctx.prec = PLANT_PREC
        try:
            self._begin()
            return [self._corrupt(e) for e in inner_out]
        finally:
            ctx.prec = old

    def _begin(self):
        pass

    def _corrupt(self, e):
        return e

class NarrowPlant(Plant):
    def __init__(self, inner, factor_str):
        super().__init__(inner, f"narrow x{factor_str}")
        self.s = arb(factor_str)

    def _corrupt(self, e):
        m = mid_exact(e)
        return m + (e - m) * self.s

class ShiftPlant(Plant):
    """rel: shift by width*rel (width = 2*rad). abs_: SCALE-ADAPTIVE
    absolute shift: applied shift = max(abs_, rad * 2^-28, |mid| *
    2^-52). A fixed abs shift is swallowed by the 30-bit mag radius
    round-up (relative quantum ~2^-29) whenever rad > ~0.5, and by
    53-bit midpoint rounding when |mid| > ~9e6 — a harness applying it
    would hard-FAIL HONEST adaptations. rad*2^-28 sits 2x above the mag
    quantum; |mid|*2^-52 clears double-ulp midpoint loss; the 1e-9 floor
    preserves the sub-slack probe at small scales. The applied shift is
    recorded PER OUTPUT (self.applied, a list — a scalar would show only
    the LAST output's scale for multi-output adaptations) and printed."""
    def __init__(self, inner, rel=None, abs_=None):
        tag = f"shift width*{rel}" if rel else f"shift abs {abs_} scaled"
        super().__init__(inner, tag)
        self.rel = arb(rel) * 2 if rel else None
        self.abs_ = arb(abs_) if abs_ else None
        self.applied = None

    def _begin(self):
        self.applied = []

    def _corrupt(self, e):
        if self.rel is not None:
            s = rad_arb(e) * self.rel
            self.applied.append(s)
            return e + s
        s = self.abs_
        for c in (rad_arb(e) * arb(2) ** -28,
                  abs(arb(_as_arb(e).mid())) * arb(2) ** -52):
            if bool(c > s):
                s = c
        self.applied.append(s)
        return e + s

class DropDimPlant(Plant):
    """One dimension's radius silently ignored (set to 0)."""
    def __init__(self, inner, dim):
        super().__init__(inner, f"drop-dim {dim}")
        self.dim = dim

    def enclose(self, p0, radii, K):
        r2 = list(radii)
        r2[self.dim] = 0.0
        return norm_out(self.inner.enclose(p0, r2, K))

# ---------------------------------------------------- generic battery core
def _vacuous(e, midhull, cfg):
    """(iii): enclosure width > vacuous_factor * sampled-value spread =>
    no-information cell. midhull is the hull of the point-ball MIDPOINTS
    ONLY — a hull including the point-ball radii lets a coarse point_eval
    inflate the spread and defeat the gate (the fat-ball launder). Exact
    arb comparison (widths can exceed float range). Returns (is_vacuous,
    ratio_str)."""
    old = ctx.prec
    ctx.prec = 256
    try:
        ew = 2 * rad_arb(e)
        sw = 2 * rad_arb(midhull)
        if bool(sw == 0):
            return (bool(ew > 0), "inf(spread=0)")
        vac = bool(ew > arb(cfg.vacuous_factor) * sw)
        try:
            return (vac, (ew / sw).str(3, radius=False) + "x")
        except Exception:
            return (vac, "n/a")
    finally:
        ctx.prec = old

def _too_coarse(e, vals, j, cfg):
    """Point-precision gate: the cell is judgeable ONLY if the max
    point-ball radius <= enclosure_width / (vacuous_factor * point_margin)
    (defaults: width / 1e4). Coarser point balls could both hide
    containment violations and inflate the sampled spread, so the cell is
    POINT-PRECISION INDET — no containment claim, no vacuousness claim,
    never informative. Non-finite point radii do NOT trip this gate (they
    fall through to containment, which fatals on non-contained values).
    Returns (too_coarse, maxr, ew)."""
    old = ctx.prec
    ctx.prec = 256
    try:
        ew = 2 * rad_arb(e)
        maxr = arb(0)
        for v in vals:
            r = rad_arb(v[j])
            if bool(r > maxr):
                maxr = r
        if not bool(maxr.is_finite()):
            return (False, maxr, ew)
        trip = bool(maxr * arb(cfg.vacuous_factor) * arb(cfg.point_margin)
                    > ew)
        return (trip, maxr, ew)
    finally:
        ctx.prec = old

def containment_shell(adapt, p0, radii, Ks, cfg, shell_idx, encmap,
                      expect_indet=frozenset()):
    """(i)+(ii)+(iii) for one shell. encmap comes PRE-COLLECTED by
    run_battery (ordering law (o): enclosures committed before any
    point_eval). Samples batch 1 (boundary + band + interior, seeded
    stream) and batch 2 (fresh closed-box, ':post' stream), evaluates
    truth, then tests every requested-K enclosure: non-finite -> INDET
    cell; point balls too coarse for this width (point-precision
    gate) -> POINT-PRECISION INDET cell; any point outside -> FATAL;
    vacuous vs sampled MIDPOINT spread -> INDET cell; else informative
    PASS cell. Returns (pts, vals, encmap, clean_Ks, pp_Ks) where clean_Ks
    are the INFORMATIVE certified orders (plant anchors) and pp_Ks are the
    orders with at least one POINT-PRECISION INDET cell."""
    t0 = time.time()
    rng = random.Random(f"{cfg.seed}:{adapt.label}:shell{shell_idx}")
    pts = sample_box(p0, radii, cfg.points, rng, prec=cfg.point_prec)
    rng2 = random.Random(f"{cfg.seed}:{adapt.label}:shell{shell_idx}:post")
    pts += sample_box(p0, radii, cfg.points, rng2, prec=cfg.point_prec,
                      boundary=False)
    vals = [norm_out(adapt.point_eval(p)) for p in pts]
    clean = []
    ppks = set()
    stat = []
    for K in Ks:
        enc = encmap[K]
        declared = (shell_idx, K) in expect_indet
        dsfx = "*" if declared else ""
        kclean = True
        for j, e in enumerate(enc):
            osfx = f"/out{j}" if len(enc) > 1 else ""
            if not finite(e):
                kclean = False
                (CELLS_DECLARED if declared else CELLS_INDET)[0] += 1
                stat.append(f"K={K}{osfx}:INDET{dsfx}")
                note(f"{adapt.label} shell{shell_idx} r={radii} K={K}{osfx}: "
                     f"non-finite enclosure (NO-INFORMATION, not a pass"
                     f"{'; declared-expected' if declared else ''})")
                continue
            coarse, maxr, ew = _too_coarse(e, vals, j, cfg)
            if coarse:
                # Even a coarse point ball PROVES a violation when it is
                # entirely DISJOINT from the enclosure: an honest point
                # ball contains the true value, so disjoint => the true
                # value is certainly outside. This cannot false-FAIL an
                # honest adaptation (its point balls contain values inside
                # the honest enclosure and therefore overlap it).
                missD = sum(0 if e.overlaps(v[j]) else 1 for v in vals)
                if missD:
                    kclean = False
                    fatal(f"CONTAINMENT VIOLATION {adapt.label} "
                          f"shell{shell_idx} r={radii} K={K}{osfx}: "
                          f"{missD}/{len(pts)} point balls PROVABLY DISJOINT "
                          f"from finite enclosure ({fmt(e)}) — proven escape "
                          f"despite coarse point_eval (max point-ball radius "
                          f"{maxr.str(3, radius=False)})")
                    continue
                kclean = False
                ppks.add(K)
                (CELLS_DECLARED if declared else CELLS_INDET)[0] += 1
                stat.append(f"K={K}{osfx}:PPREC{dsfx}")
                note(f"{adapt.label} shell{shell_idx} r={radii} K={K}{osfx}: "
                     f"POINT-PRECISION INDET — point_eval too coarse to "
                     f"judge vacuousness at this width (max point-ball "
                     f"radius {maxr.str(3, radius=False)} > enclosure width "
                     f"{ew.str(3, radius=False)} / (vacuous_factor "
                     f"{cfg.vacuous_factor:.0e} x point_margin "
                     f"{cfg.point_margin:g}); no containment claim either — "
                     f"NO-INFORMATION, not a pass"
                     f"{'; declared-expected' if declared else ''})")
                continue
            miss = sum(0 if e.contains(v[j]) else 1 for v in vals)
            if miss:
                kclean = False
                fatal(f"CONTAINMENT VIOLATION {adapt.label} shell{shell_idx} "
                      f"r={radii} K={K}{osfx}: {miss}/{len(pts)} points "
                      f"outside finite enclosure ({fmt(e)})")
                continue
            midhull = None
            for v in vals:
                m = mid_exact(v[j])
                midhull = m if midhull is None else midhull.union(m)
            vac, ratio = _vacuous(e, midhull, cfg)
            if vac:
                kclean = False
                (CELLS_DECLARED if declared else CELLS_INDET)[0] += 1
                stat.append(f"K={K}{osfx}:VAC{dsfx}({ratio})")
                note(f"{adapt.label} shell{shell_idx} r={radii} K={K}{osfx}: "
                     f"VACUOUS enclosure — width {ratio} the sampled spread "
                     f"(> gate {cfg.vacuous_factor:.0e}; NO-INFORMATION, "
                     f"not a pass{'; declared-expected' if declared else ''})")
            else:
                N_PASS[0] += 1
                CELLS_PASS[0] += 1
                stat.append(f"K={K}{osfx}:{len(pts)}/{len(pts)}({fmt(e)})")
        if kclean:
            clean.append(K)
    print(f"      {adapt.label} shell{shell_idx} r={list(radii)}: "
          + "  ".join(stat) + f"  [{time.time()-t0:.1f}s]", flush=True)
    return pts, vals, encmap, clean, ppks

def plants_fire(adapt, p0, radii, K, honest_enc, pts, vals, cfg):
    """(iv): every plant must be flagged via channel A (containment),
    channel B (certified-envelope cross-check), or non-finiteness. The
    honest enclosure at (p0, radii, K) must already be certified on an
    INFORMATIVE cell. Unfired plants split: identical planted ball =
    unresolvable-at-scale INDET; different = HARNESS BLIND hard FAIL."""
    dim = len(p0)
    plants = [NarrowPlant(adapt, "0.999999"),          # width x (1-1e-6)
              NarrowPlant(adapt, "0.5"),               # width x 0.5
              ShiftPlant(adapt, rel="1e-3"),           # center += width*1e-3
              ShiftPlant(adapt, abs_="1e-9")]          # scale-adaptive abs
    if dim > 1:
        plants += [DropDimPlant(adapt, d) for d in range(dim)]
    else:
        note(f"{adapt.label}: drop-dim plant N/A at dim=1 (vector case only)")
    for pl in plants:
        enc = pl.enclose(p0, radii, K)
        chans = []
        flags = []   # per-output fire flags (masking visible)
        multi = len(enc) > 1
        for j, e in enumerate(enc):
            pfx = f"out{j}:" if multi else ""
            fired = []
            if not finite(e):
                chans.append(pfx + "INDET(non-finite planted enclosure: no-"
                             "information can never pass)")
                flags.append("I")
                continue
            missA = sum(0 if e.contains(v[j]) else 1 for v in vals)
            if missA:
                chans.append(f"{pfx}A(containment: {missA}/{len(pts)} points "
                             f"escape)")
                fired.append("A")
            if not e.contains(honest_enc[j]):
                chans.append(pfx + "B(envelope: certified honest enclosure "
                             "lost)")
                fired.append("B")
            flags.append("+".join(fired) if fired else "-")
        sc = getattr(pl, "applied", None)
        if sc is None:
            sctxt = ""
        elif len(sc) == 1:
            sctxt = f" [applied shift {sc[0].str(3, radius=False)}]"
        else:
            sctxt = " [applied shifts " + ", ".join(
                f"out{j}={s.str(3, radius=False)}" for j, s in enumerate(sc)) \
                + "]"
        ftxt = (" [per-output flags: " + ",".join(
            f"out{j}={f}" for j, f in enumerate(flags)) + "]") if multi else ""
        if chans:
            ok(f"plant FIRES  {pl.label} @K={K} r={list(radii)} via "
               + "; ".join(chans) + sctxt + ftxt)
        elif all(_ball_same(e, honest_enc[j]) for j, e in enumerate(enc)):
            h0 = honest_enc[0]
            note(f"plant UNRESOLVABLE at this scale: {pl.label} @K={K} "
                 f"r={list(radii)} — perturbation collapsed into the ball "
                 f"representation (honest rad {rad_arb(h0).str(3, radius=False)}, "
                 f"|mid| {abs(arb(_as_arb(h0).mid())).str(3, radius=False)}"
                 f"{sctxt}); INDET for this plant, not counted as verified "
                 f"and not a harness-blind FAIL")
        else:
            fatal(f"HARNESS BLIND: plant SURVIVED {pl.label} @K={K} "
                  f"r={list(radii)} — planted enclosure differs from honest "
                  f"yet neither containment nor envelope channel flagged "
                  f"it{sctxt}")

def fd_check(adapt, fd_p0, fd_radii, cfg, npts=5):
    """(v): per-dim central FD (h=2^-fd_h_bits, prec fd_prec) vs deriv().
    Ball radii are added OUTWARD — gate on
    (|mid_fd - mid_dv| + rad_fd + rad_dv) / (|mid_dv| - rad_dv)."""
    if not hasattr(adapt, "deriv"):
        note(f"{adapt.label}: deriv() not exposed — FD arm skipped")
        return
    old = ctx.prec
    ctx.prec = cfg.fd_prec
    try:
        h = arb(2) ** -cfg.fd_h_bits
        rng = random.Random(f"{cfg.seed}:{adapt.label}:FD")
        dim = len(fd_p0)
        mids = [arb(tp.frac(str(c)).mid()) for c in fd_p0]
        for t in range(npts):
            p = [mids[d] + arb(str(rng.uniform(-0.9995, 0.9995)))
                 * arb(str(fd_radii[d])) for d in range(dim)]
            for d in range(dim):
                ctx.prec = cfg.fd_prec
                dv = norm_out(adapt.deriv(p, d))
                pp = list(p); pp[d] = p[d] + h
                pm = list(p); pm[d] = p[d] - h
                vp = norm_out(adapt.point_eval(pp, prec=cfg.fd_prec))
                vm = norm_out(adapt.point_eval(pm, prec=cfg.fd_prec))
                for j in range(len(dv)):
                    fd = (vp[j] - vm[j]) / (2 * h)
                    num = abs(arb(fd.mid()) - arb(dv[j].mid())) \
                        + arb(fd.rad()) + arb(dv[j].rad())
                    den = abs(arb(dv[j].mid())) - arb(dv[j].rad())
                    if not bool(den > 0):
                        note(f"FD {adapt.label} pt{t} dim{d}: derivative "
                             f"ball straddles zero at this scale — FD gate "
                             f"no-information (INDET, not a pass)")
                        continue
                    rd = num / den
                    if bool(rd < arb(cfg.fd_tol)):
                        ok(f"FD {adapt.label} pt{t} dim{d}: rel FD-deriv "
                           f"(radii outward) {rd.str(3)} < {cfg.fd_tol}")
                    else:
                        fatal(f"FD MISMATCH {adapt.label} pt{t} dim{d}: "
                              f"rel (radii outward) {rd.str(3)} >= "
                              f"{cfg.fd_tol}")
    finally:
        ctx.prec = old

def precision_probe(adapt, p0, radii, K, honest_enc, cfg):
    """(vi): degraded() (constants rebuilt at 53-bit) must be CAUGHT —
    non-finite, or relwidth inflated >= degrade_factor over honest."""
    if not hasattr(adapt, "degraded"):
        note(f"{adapt.label}: degraded() not exposed — precision probe "
             f"skipped (opt-in arm; see docstring (vi) for where it is "
             f"meaningful)")
        return
    deg = adapt.degraded()
    enc = norm_out(deg.enclose(p0, radii, K))
    for j, e in enumerate(enc):
        rh, rd_ = relwidth(honest_enc[j]), relwidth(e)
        if not finite(e):
            ok(f"precision probe {adapt.label} @K={K}: 53-bit rebuild -> "
               f"non-finite (caught; honest relw {rh:.2e})")
        elif rd_ >= cfg.degrade_factor * rh:
            ok(f"precision probe {adapt.label} @K={K}: 53-bit rebuild "
               f"caught — relwidth {rh:.2e} -> {rd_:.2e} "
               f"(x{rd_/rh:.1e} >= gate x{cfg.degrade_factor:.0e})")
        else:
            fatal(f"PRECISION PROBE BLIND {adapt.label} @K={K}: 53-bit "
                  f"rebuild NOT caught — relwidth {rh:.2e} -> {rd_:.2e} "
                  f"(< gate x{cfg.degrade_factor:.0e})")

def run_battery(adapt, shells, Ks, cfg, plant_at=None, fd_box=None,
                Ks_by_shell=None, expect_indet=frozenset()):
    """Full battery for one adaptation. RETURNS a BatteryResult — check
    result.ok (see docstring; exiting 0 without checking it is a consumer
    bug).
    shells: list of (p0, radii) — vectors (length-1 = univariate).
    Ks: enclosure orders to test per shell (Ks_by_shell overrides).
    plant_at: (shell_index, K) for plants + precision probe; default =
    (0, last INFORMATIVE clean K of shell 0). fd_box: (p0, radii) box for
    FD sampling; default = shell with the largest radii.
    expect_indet: set of (shell_idx, K) cells DECLARED expected-no-
    information (see docstring (iii)); they are still containment-checked
    when finite, still reported, and can never create coverage."""
    print(f"\n== battery: {adapt.label} — {len(shells)} shells, seed base "
          f"{cfg.seed}, {cfg.points} fresh closed-box points/shell x2 "
          f"batches (+exact corners/faces + near-face band), vacuous gate "
          f"x{cfg.vacuous_factor:.0e}, point-precision margin "
          f"x{cfg.point_margin:g}"
          + (f", {len(expect_indet)} cells declared expected-INDET"
             if expect_indet else "") + " ==", flush=True)
    f0, n0 = len(FATALS), len(NOTES)
    cp0, ci0, cd0 = CELLS_PASS[0], CELLS_INDET[0], CELLS_DECLARED[0]
    # ---- phase (o): EVERY enclosure committed before ANY point_eval
    encs = {}
    for si, (p0, radii) in enumerate(shells):
        ks = Ks_by_shell[si] if Ks_by_shell else Ks
        encs[si] = {K: norm_out(adapt.enclose(p0, radii, K)) for K in ks}
    # ---- phases (i)-(iii)
    cache = {}
    for si, (p0, radii) in enumerate(shells):
        ks = Ks_by_shell[si] if Ks_by_shell else Ks
        cache[si] = containment_shell(adapt, p0, radii, ks, cfg, si,
                                      encs[si], expect_indet)
    psi, pK = plant_at if plant_at else (0, None)
    pts, vals, encmap, clean, ppks = cache[psi]
    p0, radii = shells[psi]
    if pK is None:
        pK = clean[-1] if clean else None
    if pK is None or pK not in clean:
        # If the missing anchor is attributable to the POINT-PRECISION
        # gate (coarse point_eval), this is an INDET note, not a harness
        # FAIL — coarse points are the consumer's instrument being too
        # blunt to judge ANY claim at this width (the honest twin of a
        # fat-ball launder must not become a false-rejection class).
        # The run still cannot PASS: it has no informative coverage.
        pp_blocked = (pK in ppks) if pK is not None else bool(ppks)
        if pp_blocked:
            note(f"{adapt.label}: plants cannot be verified at shell {psi} "
                 f"— point_eval too coarse to judge vacuousness at this "
                 f"width (POINT-PRECISION INDET; no informative anchor). "
                 f"INDET, not a FAIL — and INDET can never pass")
        else:
            fatal(f"{adapt.label}: no certified INFORMATIVE finite enclosure "
                  f"at plant shell {psi} — plants cannot be verified (non-"
                  f"finite/vacuous coverage is no-information; INDET is not "
                  f"a pass)")
    else:
        plants_fire(adapt, p0, radii, pK, encmap[pK], pts, vals, cfg)
        precision_probe(adapt, p0, radii, pK, encmap[pK], cfg)
    if fd_box is None:
        fd_box = max(shells, key=lambda s: max(float(r) for r in s[1]))
    fd_check(adapt, fd_box[0], fd_box[1], cfg)
    res = BatteryResult(adapt.label, FATALS[f0:], NOTES[n0:],
                        CELLS_PASS[0] - cp0, CELLS_INDET[0] - ci0,
                        CELLS_DECLARED[0] - cd0)
    print(f"      -> {adapt.label} coverage: {res.cells_pass} informative "
          f"PASS cells, {res.cells_indet} no-information cells, "
          f"{res.cells_declared} declared-expected no-information cells; "
          f"result.ok={res.ok}", flush=True)
    return res

# --------------------------------------------------- reference adaptations
class ErasArm:
    """The reference univariate arms through the protocol.
    kind: 'unc' (y3-collapsed uncorrected, s3-closed) or 'corr' (4-dim
    corrected product). prec: working precision of the enclosure forms.
    sprec: precision at which the s-constants are built (53 => the
    degraded pathology variant). counts: mutated pattern counts
    (domain plant; None = pristine)."""
    def __init__(self, kind, prec, sprec=None, counts=None):
        self.kind, self.prec = kind, prec
        sp = sprec or max(prec, 192)
        old = ctx.prec
        ctx.prec = sp
        try:
            self.s1, self.s2, self.s3 = arb("0.1"), arb("0.2"), arb("0.3")
        finally:
            ctx.prec = old
        self.counts = counts
        self.groups = tp.grouped(counts) if counts is not None else tp.GROUPS
        self.label = f"{kind}@{prec}" + ("+s53" if sprec == 53 else "") \
            + ("+mut" if counts is not None else "")
        self._cache = {}
        if kind == "unc":
            # precision probe on the unc arm ONLY: the 53-bit pathology
            # lives in the collapsed-polynomial path (documented caveat,
            # measured 4.7e-3 -> ~2e9); the corr arm's spec-radius
            # enclosures are dependency-dominated (53-bit rebuild moves
            # relwidth < 1e-15 relative — measured, no-information) and
            # corr@48's working precision is itself below 53 bits.
            self.degraded = self._degraded

    # -- internal closures over the arm's own constants
    def _fnD(self, pD):
        if self.kind == "unc":
            return tp.G_unc_D(pD, self.s1, self.s2, groups=self.groups,
                              prec=self.prec)
        return tp.L_prod_D(pD, self.s1, self.s2, self.s3,
                           self.counts, corrected=True)

    def _fnS(self, pS):
        if self.kind == "unc":
            return tp.G_unc_series(pS, self.s1, self.s2, groups=self.groups)
        return tp.L_prod_S(pS, self.s1, self.s2, self.s3,
                           self.counts, corrected=True)

    def enclose(self, p0, radii, K):
        key = (str(p0[0]), str(radii[0]), K)
        if key in self._cache:
            return self._cache[key]
        old = ctx.prec
        ctx.prec = self.prec
        try:
            pball = arb(tp.frac(str(p0[0])).mid(), float(radii[0]))
            if K == 0:
                e = self._fnD(tp.D(pball, pball * 0 + 1)).v   # naive ball
            elif K == 1:
                e, _, _, _ = tp.mv_form(self._fnD, pball)     # mean-value
            else:
                e, _, _ = tp.taylor_encl(self._fnS, pball, K) # Taylor-Lagr.
            self._cache[key] = e
            return e
        finally:
            ctx.prec = old

    def point_eval(self, p, prec=None):
        old = ctx.prec
        pr = prec or 256
        ctx.prec = pr
        try:
            pD = tp.D(p[0], p[0] * 0 + 1)
            if self.kind == "unc":
                return tp.G_unc_D(pD, self.s1, self.s2, groups=self.groups,
                                  prec=pr).v
            return tp.L_prod_D(pD, self.s1, self.s2, self.s3,
                               self.counts, corrected=True).v
        finally:
            ctx.prec = old

    def deriv(self, p, dim):
        assert dim == 0
        old = ctx.prec
        ctx.prec = 512
        try:
            pD = tp.D(p[0], p[0] * 0 + 1)
            if self.kind == "unc":
                return tp.G_unc_D(pD, self.s1, self.s2, groups=self.groups,
                                  prec=512).d
            return tp.L_prod_D(pD, self.s1, self.s2, self.s3,
                               self.counts, corrected=True).d
        finally:
            ctx.prec = old

    def _degraded(self):
        return ErasArm(self.kind, self.prec, sprec=53, counts=self.counts)

def reference_battery(cfg):
    """UNIVARIATE REGRESSION: the reference arms at full spec through the
    consumer protocol. COVERAGE HONESTY: the recorded dependency-blowup
    cells are declared expected-no-information and their widths
    regression-gated against the recorded numbers; informative coverage =
    the working-configuration cells. See docstring."""
    print("\n#### REFERENCE BATTERY (univariate regression, full spec) "
          "####", flush=True)
    # soundness probe: mid() exactness (mv-form center must not leak)
    ctx.prec = 192
    pball = arb((arb(1) / 3).mid(), 1e-2)
    m = arb(pball.mid())
    assert float(arb(m.rad())) == 0.0 and pball.contains(m)
    ok("probe: arb(pball.mid()) exact (zero radius) and inside ball")

    shells = [(["1/3"], [r]) for r in ("1e-4", "1e-3", "1e-2", "5e-2", "1e-1")]
    fdbox = (["0.48"], ["0.40"])   # original FD domain p in (0.08, 0.88)
    arm_u192 = ErasArm("unc", 192)
    arm_c192 = ErasArm("corr", 192)
    arm_c48 = ErasArm("corr", 48)
    Ks_full = [0, 1, 2, 4, 8, 16]
    # Declared expected-no-information cells (the recorded blowup
    # pathology; regression-gated below — declarations are not
    # free-floating):
    #   unc@192: everything except the working configuration (0, K=16).
    #   corr@192: shell1 K<8 (recorded relw 1.9e5-2.6e6 vs true range
    #             ~1.15) and all cells at r >= 1e-2 (2.7e61 and beyond/inf).
    #   corr@48: all cells at r >= 1e-3 (same recorded rows as corr@192).
    exp_u = {(si, K) for si in range(5) for K in Ks_full} - {(0, 16)}
    exp_c192 = {(1, K) for K in (0, 1, 2, 4)} \
        | {(si, K) for si in (2, 3, 4) for K in Ks_full}
    exp_c48 = {(si, K) for si in (1, 2, 3, 4) for K in (0, 1)}
    run_battery(arm_u192, shells, Ks_full, cfg, plant_at=(0, 16),
                fd_box=fdbox, expect_indet=exp_u)
    run_battery(arm_c192, shells, Ks_full, cfg, plant_at=(0, 1),
                fd_box=fdbox, expect_indet=exp_c192)
    run_battery(arm_c48, shells, [0, 1], cfg, plant_at=(0, 1),
                fd_box=fdbox, expect_indet=exp_c48)

    # blowup regression gates: the declared-INDET cells must REPRODUCE
    # the recorded pathology numbers (fiber s=(0.1,0.2,0.3), p0=1/3 —
    # the TAYLOR_P_BLOWUP / TAYLOR_P_ORDERK rows eras.py's blowup and
    # taylork subcommands regenerate; the blowup at fixed radius is
    # steeply fiber-dependent, so the gate pins THIS fiber's numbers)
    rw_naive = relwidth(norm_out(arm_u192.enclose(["1/3"], ["1e-4"], 0))[0])
    if 8.84e21 < rw_naive < 8.84e23:
        ok(f"blowup regression: unc naive @1e-4 relwidth {rw_naive:.3e} "
           f"reproduces banked 8.84e22 (TAYLOR_P_BLOWUP, this fiber)")
    else:
        fatal(f"blowup regression BROKEN: unc naive @1e-4 relwidth "
              f"{rw_naive:.3e} vs banked 8.84e22")
    rw_k16 = relwidth(norm_out(arm_u192.enclose(["1/3"], ["1e-4"], 16))[0])
    if 2e-3 < rw_k16 < 1e-2:
        ok(f"working-config regression: unc K=16 @1e-4 relwidth "
           f"{rw_k16:.3e} reproduces banked 4.7e-3 (TAYLOR_P_ORDERK)")
    else:
        fatal(f"working-config regression BROKEN: unc K=16 @1e-4 relwidth "
              f"{rw_k16:.3e} vs banked 4.7e-3")

    # probe: K=1 Taylor form == mean-value form (same enclosure class)
    ctx.prec = 192
    s1, s2, s3 = arm_c192.s1, arm_c192.s2, arm_c192.s3
    pball = arb((arb(1) / 3).mid(), 1e-3)
    e1, _, SB = tp.taylor_encl(
        lambda pS: tp.L_prod_S(pS, s1, s2, s3), pball, 1)
    mv, _, _, FB = tp.mv_form(
        lambda pD: tp.L_prod_D(pD, s1, s2, s3), pball)
    if bool(e1.overlaps(mv)) and bool(SB.c[1].overlaps(FB.d)):
        ok("probe: K=1 series form == mv dual form at r=1e-3 (enclosure "
           "and f1(P) vs F'(P) overlap)")
    else:
        fatal("K=1 series form and mv dual form disagree at r=1e-3")

    # original domain plant: 0001:+1 planted count error
    print("\n== original domain plant 0001:+1 (anchor / fiber-equality / "
          "enclosure level) ==", flush=True)
    ctx.prec = 384
    mut = tp.mutate_counts("0001:+1")
    gmut = tp.grouped(mut)
    pfr = tp.frac("1/3")
    pD = tp.D(pfr, pfr * 0 + 1)
    (V0, _), beta = tp.upoly_pair(pD, s1, s2)
    (Vm, _), _ = tp.upoly_pair(pD, s1, s2, gmut)
    u = ((-(beta.v * s3)).exp()) ** 2
    sh_u = Vm(u).log() - arb(tp.ANCHORS["T1"][2])
    sh_c = (tp.L_prod_D(pD, s1, s2, s3, mut, corrected=True).v.log()
            - arb(tp.ANCHORS["T1"][1]))
    if bool(abs(sh_u) > arb("1e-3")) and bool(abs(sh_c) > arb("1e-3")):
        ok(f"domain plant anchor level: unc moves {sh_u.str(6)}, corr "
           f"moves {sh_c.str(6)}")
    else:
        fatal("domain plant 0001:+1 silent at anchor level")
    shc = Vm(u).log() - V0(u).log()
    shd = (tp.L_prod_D(pD, s1, s2, s3, mut, corrected=False).v.log()
           - tp.L_prod_D(pD, s1, s2, s3, corrected=False).v.log())
    if bool(abs(shc - shd) < arb("1e-40")) and bool(abs(shc) > arb("1e-2")):
        ok(f"domain plant fiber equality: {shc.str(15)} == {shd.str(15)}")
    else:
        fatal("domain plant fiber-shift equality broken")
    # enclosure level through the protocol: corr mv@1e-4, unc K=16@1e-4
    arm_c192m = ErasArm("corr", 192, counts=mut)
    arm_u192m = ErasArm("unc", 192, counts=mut)
    for arm0, armm, K in [(arm_c192, arm_c192m, 1), (arm_u192, arm_u192m, 16)]:
        e0 = norm_out(arm0.enclose(["1/3"], ["1e-4"], K))[0]
        em = norm_out(armm.enclose(["1/3"], ["1e-4"], K))[0]
        if e0.contains(arb(0)) or em.contains(arb(0)):
            fatal(f"domain plant enclosure level {arm0.label} K={K}: "
                  f"enclosure straddles 0, no log-shift claim")
            continue
        shift = em.log() - e0.log()
        if bool(not shift.contains(arb(0))) \
                and bool(abs(arb(shift.mid())) > arb("0.5")):
            ok(f"domain plant enclosure level {arm0.label} K={K} r=1e-4: "
               f"log-shift {shift.str(10)} excludes 0, |mid|>0.5")
        else:
            fatal(f"domain plant SILENT through {arm0.label} K={K} "
                  f"enclosure @1e-4")

# ------------------------- multivariate self-test arm (TEST SCAFFOLDING)
class MultivariateSelfTest:
    """TEST SCAFFOLDING — NOT SUBSTRATE. Small honest 3-parameter product
    integrand F(q1,q2,q3) = prod_{i=1..6} ( exp(-q1*a_i)*(q2+b_i) +
    q3*c_i ), strictly positive on the test boxes. Enclosures (flint/arb,
    correct math, modest object): K=0 naive interval evaluation; K>=1
    multivariate mean-value form F(m) + sum_d dF/dq_d(Box)*[-r_d, r_d]
    (mean value theorem along the segment, gradient enclosed over the
    convex box via per-dim forward-mode duals; K>1 returns the mv form —
    documented cap). Its ONLY job: prove the harness end-to-end in vector
    mode — correct adaptation passes, all plants fire (incl. dropped-dim,
    each dim), per-dim FD agrees. No degraded(): the precision probe's
    pathology class lives in the reference battery's 318-factor object,
    not in this 6-factor scaffold."""
    label = "mv3-selftest"
    dim = 3
    _A = ("0.3", "0.7", "1.1", "0.2", "0.9", "0.5")
    _B = ("0.15", "0.4", "0.25", "0.6", "0.35", "0.5")
    _C = ("0.2", "0.1", "0.3", "0.25", "0.15", "0.05")

    def __init__(self, prec=128):
        self.prec = prec
        old = ctx.prec
        ctx.prec = max(prec, 192)
        try:
            self.abc = [(arb(a), arb(b), arb(c))
                        for a, b, c in zip(self._A, self._B, self._C)]
        finally:
            ctx.prec = old

    def _F(self, x1, x2, x3):
        t = None
        for a, b, c in self.abc:
            f = (-(x1 * a)).exp() * (x2 + b) + x3 * c
            t = f if t is None else t * f
        return t

    def enclose(self, p0, radii, K):
        old = ctx.prec
        ctx.prec = self.prec
        try:
            mids = [arb(tp.frac(str(c)).mid()) for c in p0]
            balls = [arb(mids[d].mid(), float(radii[d])) for d in range(3)]
            if K == 0:
                return self._F(*balls)             # naive interval, rigorous
            E = self._F(*mids)                     # exact-center point term
            for d in range(3):
                args = [tp.D(balls[j], balls[j] * 0 + (1 if j == d else 0))
                        for j in range(3)]
                Fd = self._F(*args).d              # dF/dq_d enclosed on Box
                E = E + Fd * (balls[d] - mids[d])  # * [-r_d, +r_d]
            return E
        finally:
            ctx.prec = old

    def point_eval(self, p, prec=None):
        old = ctx.prec
        ctx.prec = prec or 256
        try:
            return self._F(*p)
        finally:
            ctx.prec = old

    def deriv(self, p, dim):
        old = ctx.prec
        ctx.prec = 512
        try:
            args = [tp.D(p[j], p[j] * 0 + (1 if j == dim else 0))
                    for j in range(3)]
            return self._F(*args).d
        finally:
            ctx.prec = old

class ScaleRegression:
    """REGRESSION against false-FAILs at large enclosure scales: an
    HONEST univariate adaptation whose enclosure at shell r=5e-2 sits at
    mid ~192 / rad ~21.5 — exactly where a fixed abs-1e-9 shift plant is
    swallowed by 30-bit mag radius round-up, so a harness without the
    scale-adaptive shift would hard-FAIL the honest consumer as 'HARNESS
    BLIND'. g(x) = 100*(exp(2x)+x^2), naive interval enclosure (rigorous,
    monotone-exact class). This arm must PASS with all plants firing or
    scale-INDET — NEVER FAIL."""
    label = "scale-regression"

    def _g(self, x):
        return ((x * 2).exp() + x * x) * 100

    def point_eval(self, p, prec=None):
        old = ctx.prec
        ctx.prec = prec or 256
        try:
            return self._g(p[0])
        finally:
            ctx.prec = old

    def enclose(self, p0, radii, K):
        old = ctx.prec
        ctx.prec = 128
        try:
            mid = arb(tp.frac(str(p0[0])).mid())
            return self._g(arb(mid.mid(), float(radii[0])))
        finally:
            ctx.prec = old

def selftest_battery(cfg):
    """MULTIVARIATE SELF-TEST: vector-mode proof of the harness, plus the
    ScaleRegression arm (see class docstrings)."""
    print("\n#### MULTIVARIATE SELF-TEST (vector mode; TEST SCAFFOLDING, "
          "not substrate) ####", flush=True)
    st = MultivariateSelfTest()
    center = ("1/3", "1/4", "1/5")
    shells = [(center, ("1e-3", "2e-3", "5e-4")),
              (center, ("5e-3", "1e-2", "2.5e-3")),
              (center, ("2e-2", "4e-2", "1e-2"))]
    run_battery(st, shells, [0, 1], cfg, plant_at=(0, 1))

    print("\n#### SCALE REGRESSION (honest adaptation at enclosure "
          "mid ~192 / rad ~21.5 — the large-scale false-FAIL point; must "
          "PASS) ####", flush=True)
    f0 = len(FATALS)
    res = run_battery(ScaleRegression(),
                      [(("3/10",), ("1e-3",)), (("3/10",), ("5e-2",))],
                      [0], cfg, plant_at=(1, 0))
    if len(FATALS) > f0 or not res.ok:
        fatal("SCALE REGRESSION FAILED: honest large-scale adaptation was "
              "rejected — the shift-plant scale adaptation has regressed")
    else:
        ok("scale regression: honest mid~192/rad~21.5 adaptation accepted "
           "(no fatals; plants fired or scale-INDET)")

# --------------------------------------------------------------- verdict
def main():
    ap = argparse.ArgumentParser(
        description="eras adversarial verification battery "
                    "(consumer-pluggable, vector-p)")
    ap.add_argument("battery", nargs="?", default="all",
                    choices=["all", "reference", "selftest"])
    ap.add_argument("--seed", type=int, default=987654321202607,
                    help="containment base seed (single-seed pass is "
                         "NECESSARY, not sufficient — rerun with several)")
    ap.add_argument("--points", type=int, default=100,
                    help="fresh random points per shell per batch "
                         "(>=100 for spec; two batches are drawn)")
    ap.add_argument("--fd-tol", default="1.1e-12",
                    help="FD-vs-deriv relative gate (univariate spec class)")
    ap.add_argument("--degrade-factor", type=float, default=1e3,
                    help="precision probe: min relwidth inflation to count "
                         "the 53-bit rebuild as caught")
    ap.add_argument("--vacuous-factor", type=float, default=1e3,
                    help="vacuousness gate: enclosure width > factor x "
                         "sampled midpoint spread => cell is INDET "
                         "no-information")
    ap.add_argument("--point-margin", type=float, default=10,
                    help="point-precision gate: max point-ball "
                         "radius must be <= enclosure width / "
                         "(vacuous-factor x this margin) or the cell is "
                         "POINT-PRECISION INDET")
    ap.add_argument("--fd-h-bits", type=int, default=40,
                    help="FD step h = 2^-bits (default 40: h^2 truncation "
                         "far below the fd-tol gate)")
    a = ap.parse_args()
    cfg = Cfg(seed=a.seed, points=a.points, fd_tol=a.fd_tol,
              degrade_factor=a.degrade_factor, fd_h_bits=a.fd_h_bits,
              vacuous_factor=a.vacuous_factor, point_margin=a.point_margin)
    t00 = time.time()
    print(f"eras adversarial battery: {a.battery}, seed {cfg.seed}, "
          f"{cfg.points} points/shell/batch, fd-tol {cfg.fd_tol}, "
          f"degrade-factor {cfg.degrade_factor:.0e}, vacuous-factor "
          f"{cfg.vacuous_factor:.0e}, point-margin {cfg.point_margin:g}",
          flush=True)
    if a.battery in ("all", "reference"):
        reference_battery(cfg)
    if a.battery in ("all", "selftest"):
        selftest_battery(cfg)
    print(f"\nNOTES ({len(NOTES)}):")
    for n in NOTES:
        print("  " + n)
    print(f"\nFATALS ({len(FATALS)}):")
    for f in FATALS:
        print("  " + f)
    import resource
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    if sys.platform == "darwin":
        rss /= 1024  # ru_maxrss is bytes on macOS, KiB on Linux
    cov = (f"coverage cells: {CELLS_PASS[0]} informative PASS, "
           f"{CELLS_INDET[0]} no-information, {CELLS_DECLARED[0]} "
           f"declared-expected no-information")
    if FATALS:
        verdict, code = "FAIL", 1
    elif N_PASS[0] == 0 or CELLS_PASS[0] == 0:
        verdict, code = "INDET (no informative-enclosure coverage)", 2
    elif CELLS_INDET[0] > CELLS_PASS[0]:
        verdict, code = ("INDET (majority no-information coverage: "
                         f"{CELLS_INDET[0]} INDET > {CELLS_PASS[0]} PASS "
                         "cells)"), 2
    else:
        verdict, code = "PASS", 0
    print(f"\nVERIFY {verdict} — {N_PASS[0]} checks passed, "
          f"{len(NOTES)} no-information notes, {len(FATALS)} fatals; {cov} "
          f"({time.time()-t00:.1f}s, maxrss {rss:.1f} MB)")
    return code

if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""
epsfan.py — Laurent-coefficient extraction around eps=0: Cauchy eps circle
(mode 1) and real geometric-grid Vandermonde solve (mode 2), plus an anchor
round-trip probe.

Mode 1 (cauchy_laurent):
    c_k = (1/2*pi*i) * oint f(eps)/eps^(k+1) deps
        = R^(-k)/N * sum_j f(R*w^j) * w^(-j*k),   w = exp(2*pi*i/N)

Mode 2 (vandermonde_laurent): real-eps geometric decade grid + square
Vandermonde solve with interleaved held-out verify nodes (a production
Laurent-gate node design — see the function docstring for the measured
design rule and provenance).

Design rules and their provenance
---------------------------------
* Default R='1e-3', nodes=32: measured calibration — a tiny radius destroys
  the extraction through 1/eps pole leverage (R=1e-17 was measured to lose
  ~50 digits); R~1e-3 keeps the leverage penalty at ~3 digits per Laurent
  order while keeping the alias suppression R^nodes astronomically strong.
  The same extractor design carried a production reference
  ("Cauchy circle r=1/32, npts=256, dps=220").
* Working precision: dividing by R^k amplifies absolute sample error by
  R^(kmin-k); the module works at wp = dps + ceil((kmax-kmin)*log10(1/R))
  + 15 so the requested dps survives the pole leverage.
* diag with a MEASURED alias/noise estimate: the eps-grid
  alias/noise tradeoff must be measured at two settings (r, delta) — never
  trust a single extraction.  cauchy_laurent always recomputes at (2R, nodes) (or
  (R, nodes/2) with check='half_nodes') and reports per-k agreement digits
  between the two extractions.  Aliasing (a coefficient c_m with
  m ≡ k mod N leaking into c_k with weight R^(m-k)) is R- and N-dependent,
  so it shows up as measured disagreement — no silent trust.
  tests/test_epsfan.py includes a deliberately-aliased case proving the
  diag catches it.
* Full circle vs conjugate mirror (a production full-chain finding): the
  17-of-32 conjugate-mirror shortcut is valid ONLY for Schwarz-real rows
  (ALL Laurent coefficients real).  Complex-Laurent rows — most
  physical-region rows — need the FULL node circle; a 17-node
  least-squares extraction was MEASURED truncation-limited to ~13 digits
  (catastrophically worse, not degraded).  Hence `schwarz_real=False` is
  the default (full circle); the mirror is opt-in and carries a loud
  docstring warning.
* Precision discipline: explicit dps argument, all parsing inside
  mp.workdps (the mpf import-dps footgun).  Global mp.dps is
  never touched.

f may return a scalar (mpc/mpf) or a list/tuple of them (vector Laurent
fan); coeffs entries are then lists.  Total cost: 2*nodes evaluations of f
(nodes + the check pass), halved-ish by schwarz_real=True.

Certified layer and extractor mode
----------------------------------
* vandermonde_laurent_extract runs E4 as an EXTRACTOR, not just a
  certifier.  The engine owns the fan: node counts n1/n2 come from the
  certified-layer f(dps) rule (n1 = max(nodes_seed, P+4, P+ceil(dps/25)),
  n2 = ceil(1.5*n1)), node VALUES are requested from a caller-supplied
  node generator (node_gen(n, exclude) -> exact 'p/q' strings), and
  samples are generated on demand through the caller's transport callback
  f at those nodes.  It is a pure additive wrapper that delegates to
  vandermonde_laurent_certified with the generated explicit grids, so the
  cond monitor / wp escalation / held-out verify / two-grid / k_certify
  machinery is shared by construction (a test asserts coefficient equality
  against a direct certified call on the same grids).
* The vandermonde_laurent_certified condition monitor ('vand-cond') is
  live on every path: if the wp escalation exhausts its rounds, the passes
  are re-run and cond is re-measured at the FINAL wp, and the monitor
  gates that post-escalation state (it fires under cond-estimator
  saturation, where measured cond tracks ~10^wp while the true cond
  outruns every escalation round; genuine rational grids do not saturate,
  they collapse to the named singular raise below, so the monitor is the
  backstop against estimator pathologies).  On the break path the
  escalation exit condition itself subsumes the monitor (margin >= 1e-20
  by construction).  A numerically SINGULAR fit matrix (node collapse at
  working precision) is reported as the charter-named
  EpsFanCertifyError('vand-cond', cond=inf, singular=True) with
  node-separation diagnostics rather than an anonymous ZeroDivisionError.
* vandermonde_laurent_certified takes two optional parameters
  `eps_nodes1` / `eps_nodes2`: explicit REAL node lists (exact 'p/q'
  strings preferred; the fixed-eps transported-values pattern: dict-backed
  f + 1/q prime grids).  Both-or-neither; per-list dedupe, cross-list
  disjointness (parsed values), positivity, and the f(dps)-seeded MINIMUM
  node counts (n1 >= max(nodes_seed, P+4, P+ceil(dps/25)),
  n2 >= ceil(1.5*n1)) are ENFORCED (ValueError, fail-closed); the max
  node must sit inside R0/10 when R0 is finite (with explicit nodes the
  grid cannot be shrunk, so this REFUSES instead of clamping).
  `k_certify=(ka, kb)` is an optional certification subwindow: the
  two-grid per-coefficient RAISE gate is enforced on ka..kb only (e.g. a
  quoted window whose higher fitted orders are tail-absorbers and are not
  certifiable at any claim); ALL orders still report per-k agreement in
  diag, and grid_min_agree_digits stays the full-window min.  Default
  k_certify=None = full-window gate.  Cond monitor / wp escalation /
  held-out verify gates are always full-fit.
* The E4 certified layer (cauchy_laurent_certified /
  vandermonde_laurent_certified / eps_analyticity_radius /
  desystem_eps_polys / EpsFanCertifyError) removes the alias floor of a
  frozen (r, M) choice: exact runtime analyticity radius R0 from the DE's
  eps-denominator roots (REFUSES to run without singularity data or a
  caller R0), control-circle Cauchy alias bound with geometric M
  escalation / r shrinking to caps (fail-closed raise naming
  bound/tol/r/M/R0), cross-radius r-vs-r/2 belt gate, two-grid (n vs 1.5n
  shifted) + condition-monitored Vandermonde variant, and working
  precision computed from each extraction's OWN conditioning (pole
  leverage + node count + sampled magnitude / measured cond(V)), not a
  constant.  The legacy entry points cauchy_laurent / vandermonde_laurent
  / anchor_roundtrip keep their contracts.
"""

import math
from fractions import Fraction

import mpmath as mp

__all__ = ["cauchy_laurent", "vandermonde_laurent", "anchor_roundtrip",
           "cauchy_laurent_certified", "vandermonde_laurent_certified",
           "vandermonde_laurent_extract", "extract_node_counts",
           "eps_analyticity_radius", "desystem_eps_polys",
           "EpsFanCertifyError"]


def _f_row(f, z, n_comp, who):
    """Evaluate f at z, normalize to a list of mpc, check vector length."""
    val = f(z)
    if isinstance(val, (list, tuple)):
        row = [mp.mpc(v) for v in val]
    else:
        row = [mp.mpc(val)]
    if n_comp is not None and len(row) != n_comp:
        raise ValueError("%s: f returned inconsistent vector length "
                         "(%d vs %d)" % (who, len(row), n_comp))
    return row


def _extract(f, Rv, N, kmin, kmax, schwarz_real=False, evals=None):
    """One Cauchy-fan pass at radius Rv with N nodes.

    Must be called inside `with mp.workdps(wp)`.  Returns (coeffs, n_comp)
    where coeffs[k] is a list of n_comp mpc values (n_comp=1 for scalar f).

    schwarz_real=True fills the lower half circle by conjugate mirroring
    (f(conj z) = conj f(z)) — evaluates only j = 0..N//2.  VALID ONLY for
    Schwarz-real f; see cauchy_laurent docstring.  `evals` (optional list)
    accumulates the number of f evaluations (measured-cost handle).
    """
    samples = []
    n_comp = None
    n_eval_j = (N // 2 + 1) if schwarz_real else N
    for j in range(N):
        if schwarz_real and j >= n_eval_j:
            # conjugate mirror: z_j = conj(z_{N-j}) on a real-R circle
            row = [mp.conj(v) for v in samples[N - j]]
        else:
            # exact rational angle: e^{2*pi*i*j/N} = expjpi(2j/N)
            z = Rv * mp.expjpi(mp.mpf(2 * j) / N)
            row = _f_row(f, z, n_comp, "cauchy_laurent")
            if evals is not None:
                evals.append(1)
        if n_comp is None:
            n_comp = len(row)
        samples.append(row)
    coeffs = {}
    for k in range(kmin, kmax + 1):
        acc = [mp.mpc(0)] * n_comp
        for j in range(N):
            ph = mp.expjpi(mp.mpf(-2 * j * k) / N)
            for c in range(n_comp):
                acc[c] += samples[j][c] * ph
        scale = mp.power(Rv, -k) / N
        coeffs[k] = [a * scale for a in acc]
    return coeffs, n_comp


def cauchy_laurent(f, kmin, kmax, dps, R='1e-3', nodes=32, check='R2',
                   schwarz_real=False):
    """Extract Laurent coefficients c_k, k = kmin..kmax, of f around eps=0.

    Args:
      f:     callable eps -> mpc (or list/tuple of mpc for a vector fan).
             Called at the module's internal working precision.
      kmin, kmax: integer Laurent-order window (kmin may be negative).
      dps:   target decimal precision of the returned coefficients.
      R:     circle radius (string preferred — parsed inside workdps).
             Default '1e-3' (see module docstring calibration note).
      nodes: number of equally spaced sample nodes N.  Must satisfy
             N >= kmax - kmin + 1 (window must fit in one DFT period).
      check: 'R2'          -> re-extract at (2R, nodes)      [default]
             'half_nodes'  -> re-extract at (R, nodes//2)
             The second pass feeds the measured (r,delta) diagnostics.
      schwarz_real: OPT-IN conjugate-mirror shortcut.  ***WARNING — valid
             ONLY for Schwarz-real rows***, i.e. f with ALL Laurent
             coefficients REAL (f(conj eps) = conj f(eps)).  Evaluates only
             the upper half circle (N//2 + 1 nodes, e.g. 17 of 32) and
             mirrors the rest by conjugation.  On a complex-Laurent row
             (most physical-region rows) the mirror forcibly symmetrizes
             the samples: the DFT of a conjugate-symmetric sequence is
             REAL, so the imaginary content of every coefficient is
             silently destroyed — and the corruption at the affected order
             is R-INDEPENDENT, so the (r,delta) two-pass diag does NOT
             reliably catch it there (tests/test_epsfan.py E7 demonstrates
             a ~50%-wrong coefficient with a clean-looking per-k diag at
             that order).  A production full-chain run measured the
             half-circle rescue attempt (17-node LS) truncation-limited to
             ~13 digits — catastrophically worse, not degraded
             (measured).  Default False = full circle, always safe.

    Returns:
      {'coeffs': {k: mpc}            (scalar f)
                 {k: [mpc, ...]}     (vector f),
       'diag':   {'R','R2','nodes','nodes2','check','wp',
                  'schwarz_real': bool,
                  'n_evals': int,               # measured f-evaluation count
                  'agree_digits': {k: float},   # per-k measured agreement
                  'min_agree_digits': float,    # min over k (and components)
                  'gscale': str,                # max |c_k| (magnitude scale)
                  'calibration': str}}

    agree_digits[k] = -log10(|c_k^(pass1) - c_k^(pass2)| / scale_k), capped
    at dps, with scale_k = max(|c1|,|c2|); a coefficient consistent with 0
    at the dps level (both passes below gscale*10^-dps) reports dps.  Treat
    min_agree_digits as the honest digit count of the extraction — if it is
    below dps the window/radius/nodes are aliased or noise-limited.
    """
    kmin = int(kmin)
    kmax = int(kmax)
    if kmax < kmin:
        raise ValueError("cauchy_laurent: kmax < kmin")
    nodes = int(nodes)
    width = kmax - kmin + 1
    if nodes < width:
        raise ValueError(
            "cauchy_laurent: nodes=%d cannot resolve window width %d "
            "(need nodes >= kmax-kmin+1)" % (nodes, width))
    if check not in ("R2", "half_nodes"):
        raise ValueError("cauchy_laurent: check must be 'R2' or 'half_nodes'")
    if check == "half_nodes" and nodes // 2 < width:
        raise ValueError(
            "cauchy_laurent: check='half_nodes' needs nodes//2 >= window "
            "width %d" % width)

    # pole-leverage working precision (see module docstring)
    with mp.workdps(30):
        R_probe = mp.mpmathify(R) if isinstance(R, str) else mp.mpf(R)
        if not (R_probe > 0):
            raise ValueError("cauchy_laurent: R must be > 0")
        leverage = (kmax - kmin) * max(-mp.log10(R_probe), mp.mpf(0))
        wp = int(dps + math.ceil(float(leverage)) + 15)

    evals = []
    with mp.workdps(wp):
        Rv = mp.mpmathify(R) if isinstance(R, str) else +mp.mpf(R)
        c1, n_comp = _extract(f, Rv, nodes, kmin, kmax,
                              schwarz_real=schwarz_real, evals=evals)
        if check == "R2":
            R2v, nodes2 = 2 * Rv, nodes
        else:
            R2v, nodes2 = Rv, nodes // 2
        c2, _ = _extract(f, R2v, nodes2, kmin, kmax,
                         schwarz_real=schwarz_real, evals=evals)

        # global magnitude scale for zero-consistency
        gscale = mp.mpf(0)
        for k in c1:
            for v in c1[k]:
                gscale = max(gscale, abs(v))
        zero_floor = gscale * mp.mpf(10) ** (-dps)

        agree = {}
        min_agree = float(dps)
        for k in range(kmin, kmax + 1):
            worst = float(dps)
            for c in range(n_comp):
                a, b = c1[k][c], c2[k][c]
                diff = abs(a - b)
                scale = max(abs(a), abs(b))
                if scale <= zero_floor:
                    ad = float(dps)      # consistent with 0 at the dps level
                elif diff == 0:
                    ad = float(dps)
                else:
                    ad = float(-mp.log10(diff / scale))
                    ad = max(0.0, min(ad, float(dps)))
                worst = min(worst, ad)
            agree[k] = worst
            min_agree = min(min_agree, worst)
        gscale_str = mp.nstr(gscale, 8)

    with mp.workdps(dps):
        if n_comp == 1:
            coeffs = {k: +c1[k][0] for k in c1}
        else:
            coeffs = {k: [+v for v in c1[k]] for k in c1}

    diag = {
        "R": str(R),
        "R2": ("2*R" if check == "R2" else str(R)),
        "nodes": nodes,
        "nodes2": nodes2,
        "check": check,
        "wp": wp,
        "schwarz_real": bool(schwarz_real),
        "n_evals": len(evals),
        "agree_digits": agree,
        "min_agree_digits": min_agree,
        "gscale": gscale_str,
        "calibration": ("default R=1e-3/nodes=32 (measured: R=1e-17 loses "
                        "~50d to 1/eps pole leverage); agreement digits are "
                        "MEASURED between two passes — trust "
                        "min_agree_digits, not dps"
                        + ("; schwarz_real MIRROR active — valid ONLY if "
                           "ALL Laurent coefficients are real, see "
                           "docstring warning" if schwarz_real else "")),
    }
    return {"coeffs": coeffs, "diag": diag}


def vandermonde_laurent(f, kmin, kmax, dps, eps_max='1e-3', span='10',
                        n_verify=7, eps_nodes=None):
    """Laurent window c_k, k = kmin..kmax, from REAL geometric eps nodes
    (second validated extraction mode, alongside the Cauchy circle).

    Node-design rule (MEASURED, copied with provenance from a production
    Laurent gate):
      * real nodes eps = 1/q, q = 1013..10037 (24 primes): a GEOMETRIC grid
        with total q-ratio ~ 10 (one decade; per-step ratio ~ 1.105) and
        eps_max = 1/1013 ~ 1e-3 (same leverage calibration as the Cauchy
        default R — going deeper multiplies the 1/eps^|kmin| pole leverage);
      * square Vandermonde solve on P = kmax-kmin+1 fit nodes for the whole
        window;
      * the remaining nodes are INTERLEAVED verify spares, spread across
        the band ("conditioning + honest check"): the held-out residual is
        the honest extraction error — reported digits = min(oracle match,
        self-consistency), never the fit's own claim.

    Args:
      f:        callable eps -> mpc (or list/tuple for a vector fan).  For
                precomputed fixed-eps transported values use a dict-backed
                callable + explicit `eps_nodes`.
      kmin, kmax: integer Laurent-order window.
      dps:      target decimal precision of the returned coefficients.
      eps_max:  largest node (string preferred), default '1e-3'.
      span:     total geometric span eps_max/eps_min, default '10' (one
                decade — the measured production design point).
      n_verify: held-out verify nodes (default 7, the measured design's
                24-17). Must be
                >= 1: with no verify nodes the extraction error is
                UNMEASURED and this module refuses to fabricate a digit
                claim.
      eps_nodes: optional explicit REAL node list (str/mpf, > 0) overriding
                the generated grid, e.g. exact rationals 1/q for
                precomputed transports.  len must be >= P + 1; all but P
                interleaved nodes are used as verify spares.

    Working precision: wp = dps + 4*P + 20.  Calibration measured on this
    box (known irrational fan, decade grid): worst-case relative
    coefficient loss of the square solve is ~ 3.6*P - 4 digits and is
    P-dominated (window-offset independent: P=5/7/17 lost 13.5/21/60
    digits at both dps 80 and 120/200).

    Returns:
      {'coeffs': {k: mpc} (scalar) / {k: [mpc,...]} (vector),
       'diag': {'mode': 'vandermonde_realgrid', 'eps_min','eps_max','span',
                'n_nodes','n_fit','n_verify','wp',
                'self_digits': float,   # MEASURED held-out digits — the
                                        # honest count, trust it, not dps
                'design': provenance str}}
    """
    kmin = int(kmin)
    kmax = int(kmax)
    if kmax < kmin:
        raise ValueError("vandermonde_laurent: kmax < kmin")
    P = kmax - kmin + 1
    n_verify = int(n_verify)
    if n_verify < 1:
        raise ValueError(
            "vandermonde_laurent: n_verify must be >= 1 — with no held-out "
            "verify nodes the extraction error is unmeasured (no fabricated "
            "digit claims; the reference design used 7 spares)")

    wp = int(dps) + 4 * P + 20
    with mp.workdps(wp):
        if eps_nodes is not None:
            eps_list = []
            for e in eps_nodes:
                try:
                    ev = (mp.mpmathify(e) if isinstance(e, str)
                          else mp.mpc(e))
                except Exception as exc:
                    raise ValueError("vandermonde_laurent: unparseable "
                                     "eps node %r (%s)" % (e, exc))
                eps_list.append(ev)
            if len(eps_list) < P + 1:
                raise ValueError(
                    "vandermonde_laurent: need >= P+1 = %d explicit nodes "
                    "(P fit + >=1 verify), got %d" % (P + 1, len(eps_list)))
        else:
            eps_maxv = (mp.mpmathify(eps_max) if isinstance(eps_max, str)
                        else +mp.mpf(eps_max))
            spanv = (mp.mpmathify(span) if isinstance(span, str)
                     else +mp.mpf(span))
            if not (eps_maxv > 0):
                raise ValueError("vandermonde_laurent: eps_max must be > 0")
            if not (spanv > 1):
                raise ValueError("vandermonde_laurent: span must be > 1")
            n = P + n_verify
            eps_list = [eps_maxv * spanv ** (-mp.mpf(i) / (n - 1))
                        for i in range(n)]
        for e in eps_list:
            if not (mp.im(e) == 0 and mp.re(e) > 0):
                raise ValueError("vandermonde_laurent: nodes must be real "
                                 "and > 0 (real-grid mode)")
        eps_list = sorted([mp.mpf(mp.re(e)) for e in eps_list])
        n = len(eps_list)
        k_spare = n - P
        # interleaved verify spares, spread across the band
        # ("conditioning + honest check")
        ver_idx = set(round((i + 1) * n / (k_spare + 1))
                      for i in range(k_spare))
        fit_eps = [e for i, e in enumerate(eps_list) if i not in ver_idx]
        ver_eps = [e for i, e in enumerate(eps_list) if i in ver_idx]
        if len(fit_eps) != P:
            raise ValueError(
                "vandermonde_laurent: interleave gave %d fit nodes, need "
                "square P=%d (increase n_verify or pass distinct eps_nodes)"
                % (len(fit_eps), P))

        # sample f on all nodes (fit + verify)
        n_comp = None
        fit_rows, ver_rows = [], []
        for e in fit_eps:
            row = _f_row(f, e, n_comp, "vandermonde_laurent")
            n_comp = len(row) if n_comp is None else n_comp
            fit_rows.append(row)
        for e in ver_eps:
            row = _f_row(f, e, n_comp, "vandermonde_laurent")
            ver_rows.append(row)

        V = mp.matrix(P, P)
        for r in range(P):
            for c in range(P):
                V[r, c] = mp.mpc(fit_eps[r]) ** (kmin + c)
        sols = []
        for comp in range(n_comp):
            b = mp.matrix([fit_rows[r][comp] for r in range(P)])
            sols.append(mp.lu_solve(V, b))

        # held-out self-consistency: the honest extraction error
        selferr = mp.mpf(0)
        for vi, e in enumerate(ver_eps):
            for comp in range(n_comp):
                pred = sum(sols[comp][c] * mp.mpc(e) ** (kmin + c)
                           for c in range(P))
                direct = ver_rows[vi][comp]
                scale = max(abs(direct), abs(pred))
                if scale > 0:
                    selferr = max(selferr, abs(pred - direct) / scale)
        if selferr > 0:
            self_digits = max(0.0, min(float(-mp.log10(selferr)),
                                       float(dps)))
        else:
            self_digits = float(dps)
        eps_min_str = mp.nstr(eps_list[0], 8)
        eps_max_str = mp.nstr(eps_list[-1], 8)

    with mp.workdps(dps):
        if n_comp == 1:
            coeffs = {kmin + c: +sols[0][c] for c in range(P)}
        else:
            coeffs = {kmin + c: [+sols[comp][c] for comp in range(n_comp)]
                      for c in range(P)}

    diag = {
        "mode": "vandermonde_realgrid",
        "eps_min": eps_min_str,
        "eps_max": eps_max_str,
        "span": str(span) if eps_nodes is None else "explicit nodes",
        "n_nodes": n,
        "n_fit": P,
        "n_verify": len(ver_eps),
        "wp": wp,
        "self_digits": self_digits,
        "design": ("production Laurent-gate node design: "
                   "geometric real grid, total span ~10 (one decade), "
                   "eps_max ~1e-3, square Vandermonde on P fit nodes + "
                   "interleaved verify spares; self_digits is MEASURED on "
                   "the held-out nodes — trust it, not dps"),
    }
    return {"coeffs": coeffs, "diag": diag}


def anchor_roundtrip(f, coeffs, dps, eps_anchor='6.1e-4'):
    """Anchor round-trip probe: resum extracted Laurent coefficients at a
    fresh anchor eps and compare against a direct evaluation there.

    Both sides of the comparison carry the SAME upstream boundary/anchor
    content (the transported value f(eps_anchor) and the coefficients
    extracted from transported values share every boundary constant), so
    the boundary CANCELS — the residual isolates extraction error, never
    the boundary.  On a window-limited Laurent (all orders inside
    [kmin, kmax], e.g. the known-fan unit test) the truncation term is zero
    and the probe measures PURE extraction error; on a generic f the
    reported digits are a lower bound, truncated at roughly
    |c_{kmax+1} * eps_anchor^{kmax+1}| relative to |f(eps_anchor)|.

    Args:
      f:          callable eps -> value(s), or a PRECOMPUTED value at
                  eps_anchor (mpf/mpc/str, or list/tuple for a vector fan)
                  for the fixed-eps transported-values use case.
      coeffs:     {k: mpc} or {k: [mpc, ...]} as returned by
                  cauchy_laurent / vandermonde_laurent.
      dps:        cap for the reported agreement digits (and rounding
                  precision of the reported strings).
      eps_anchor: probe point (string preferred; parsed inside workdps).
                  Default '6.1e-4' — inside the calibrated [1e-4, 1e-3]
                  band but on NEITHER default grid (not a Cauchy node of
                  R=1e-3/32 and not a decade-grid node).

    Returns:
      {'agree_digits': float,           # min over components
       'per_component': [float, ...],
       'anchor': str, 'resummed': [str, ...], 'direct': [str, ...],
       'note': str}
    """
    if not coeffs:
        raise ValueError("anchor_roundtrip: empty coeffs")
    ks = sorted(coeffs.keys())
    kmin = ks[0]
    with mp.workdps(30):
        ea = (mp.mpmathify(eps_anchor) if isinstance(eps_anchor, str)
              else mp.mpc(eps_anchor))
        leverage = max(-kmin, 0) * max(float(-mp.log10(abs(ea))), 0.0)
    wp = int(dps + math.ceil(leverage) + 15)

    with mp.workdps(wp):
        ea = (mp.mpmathify(eps_anchor) if isinstance(eps_anchor, str)
              else mp.mpc(eps_anchor))
        # normalize coeffs to lists
        rows = {}
        n_comp = None
        for k in ks:
            v = coeffs[k]
            row = ([mp.mpc(x) for x in v] if isinstance(v, (list, tuple))
                   else [mp.mpc(v)])
            if n_comp is None:
                n_comp = len(row)
            elif len(row) != n_comp:
                raise ValueError("anchor_roundtrip: inconsistent component "
                                 "count across coefficient orders")
            rows[k] = row
        resummed = [sum(rows[k][c] * ea ** k for k in ks)
                    for c in range(n_comp)]
        if callable(f):
            direct = _f_row(f, ea, n_comp, "anchor_roundtrip")
        else:
            direct = ([mp.mpc(mp.mpmathify(str(x)))
                       for x in f] if isinstance(f, (list, tuple))
                      else [mp.mpc(mp.mpmathify(str(f)))])
        if len(direct) != n_comp:
            raise ValueError("anchor_roundtrip: direct value has %d "
                             "components, coeffs have %d"
                             % (len(direct), n_comp))
        per = []
        for c in range(n_comp):
            diff = abs(resummed[c] - direct[c])
            scale = max(abs(resummed[c]), abs(direct[c]))
            if diff == 0 or scale == 0:
                per.append(float(dps))
            else:
                per.append(max(0.0, min(float(-mp.log10(diff / scale)),
                                        float(dps))))
        res_str = [mp.nstr(v, min(int(dps), 30)) for v in resummed]
        dir_str = [mp.nstr(v, min(int(dps), 30)) for v in direct]

    return {
        "agree_digits": min(per),
        "per_component": per,
        "anchor": str(eps_anchor),
        "resummed": res_str,
        "direct": dir_str,
        "note": ("boundary cancels — residual isolates extraction error "
                 "(plus window truncation unless f is window-limited)"),
    }


# ===========================================================================
# E4 certified layer
# ===========================================================================
#
# THE DEFECT this layer removes: DFT extraction on a frozen circle radius r
# with a frozen node count M aliases c_k <- c_k + sum_{j>=1} c_{k+jM} r^{jM}.
# With (r, M) frozen the alias term is a dps-INDEPENDENT error floor
# (measured in production: r=1/32, M frozen -> ~99d floor).  Two-dps reruns agree
# perfectly with each other and are both wrong past the floor — cranking dps
# is a lie there.  The certified layer makes (r, M) runtime quantities chosen
# against an explicit analyticity radius and a certified alias bound, and
# fails CLOSED (raise, never a silent number) whenever it cannot certify.


class EpsFanCertifyError(RuntimeError):
    """Fail-closed refusal from the certified eps-fan layer.

    kind:
      'alias-cap'   (r, M) escalation hit its caps with the certified alias
                    bound still above tolerance;
      'belt'        cross-radius (r vs r/2) per-coefficient gate failed —
                    a conditioning/contamination error the Cauchy bound
                    cannot see (e.g. non-analytic contamination);
      'vand-grid'   Vandermonde two-grid (n1 vs shifted n2) per-coefficient
                    gate failed;
      'vand-verify' Vandermonde held-out verify residual above ~10^-(dps-1);
      'vand-cond'   condition monitor (live on every path):
                    cond(V)*10^-wp > 10^-(dps+guard) with cond RE-MEASURED
                    at the FINAL post-escalation wp (on the break path the
                    wp escalation subsumes the monitor: margin >= 1e-20 by
                    its exit condition); also raised with cond=inf when the
                    Vandermonde fit matrix is numerically SINGULAR at wp
                    (node collapse at working precision — previously an
                    unnamed ZeroDivisionError);
      'vand-wp-cap' conditioning-driven working precision exceeded max_wp.

    .info carries the named numbers (bound, tol, r, M, R0, ...) as strings —
    the spec requires the raise to NAME them, so tests/pipelines can assert
    on them instead of re-parsing the message.
    """

    def __init__(self, kind, msg, **info):
        super().__init__(msg)
        self.kind = kind
        self.info = {k: str(v) for k, v in info.items()}


_R0_REFUSAL = (
    "%s: no eps-singularity data and no caller R0 — REFUSING to assume an "
    "analyticity radius (certified-layer rule: never silently assume). "
    "Pass R0=<radius> (str/Fraction/mpf) or eps_polys=<exact eps-polynomial "
    "coefficient lists> (e.g. desystem_eps_polys(desys, x) on a de_load "
    "DESystem, whose parsers keep the rationals exact).")


def _exact_fraction(v, what):
    """Exact-rational parse: int / Fraction / 'p/q' string.  Floats are
    REFUSED (inexact input to an exactness-critical path)."""
    if isinstance(v, Fraction):
        return v
    if isinstance(v, int):
        return Fraction(v)
    if isinstance(v, str):
        try:
            return Fraction(v.strip())
        except Exception as exc:
            raise ValueError("%s: unparseable exact rational %r (%s)"
                             % (what, v, exc))
    raise ValueError("%s: need int/Fraction/'p/q' string, got %r "
                     "(floats refused — exact rationals only)" % (what, v))


def _parse_radius(v, what):
    """Parse a positive radius (str preferred; Fraction/mpf/int accepted;
    mp.inf allowed).  Must be called inside mp.workdps."""
    if v is mp.inf:
        return mp.inf
    if isinstance(v, str):
        rv = mp.mpmathify(v)
    elif isinstance(v, Fraction):
        rv = mp.mpf(v.numerator) / mp.mpf(v.denominator)
    else:
        rv = mp.mpf(v)
    if not (rv > 0):
        raise ValueError("%s: radius must be > 0, got %r" % (what, v))
    return rv


def eps_analyticity_radius(eps_polys, dps=60):
    """Exact-data analyticity radius: R0 = min |nonzero eps-root| over the
    given eps-polynomials (certified-layer rule 1).

    eps_polys: iterable of polynomials in eps, each a coefficient list
    (low->high) of exact rationals (int / Fraction / 'p/q' string) — the
    DE's leading-coefficient/denominator polynomials in eps, e.g. from
    desystem_eps_polys().  Floats are refused.

    Roots via mp.polyroots at dps+20 with the returned error bound
    SUBTRACTED (certified lower bound on the singularity distance).  Roots
    at eps=0 are the expansion point and are excluded EXACTLY (low-order
    zero coefficients stripped before root-finding).  A numerically
    uncertifiable near-zero root (|root| <= 2*err with nonzero constant
    term) raises rather than guessing.

    Returns mp.mpf R0, or mp.inf if the data has no nonzero roots at all
    (entire in eps).  Raises ValueError on empty/absent data — the caller
    must then supply R0 explicitly (never silently assume).
    """
    if eps_polys is None:
        raise ValueError(_R0_REFUSAL % "eps_analyticity_radius")
    polys = list(eps_polys)
    if not polys:
        raise ValueError(
            "eps_analyticity_radius: EMPTY singularity data — refusing to "
            "assume an analyticity radius; pass R0 explicitly if the system "
            "is known to be eps-regular")
    best = None
    with mp.workdps(int(dps) + 20):
        for pi, p in enumerate(polys):
            cs = [_exact_fraction(c, "eps_polys[%d]" % pi) for c in p]
            while cs and cs[-1] == 0:
                cs.pop()
            lo = 0
            while lo < len(cs) and cs[lo] == 0:
                lo += 1          # exact eps^lo factor: roots AT 0 excluded
            cs = cs[lo:]
            if len(cs) <= 1:
                continue         # zero or constant poly: no roots off 0
            coeffs = [mp.mpf(c.numerator) / mp.mpf(c.denominator)
                      for c in reversed(cs)]
            roots, err = mp.polyroots(coeffs, error=True, maxsteps=200,
                                      extraprec=120)
            for rt in roots:
                mag = abs(rt)
                if mag <= 2 * err:
                    raise ValueError(
                        "eps_analyticity_radius: eps_polys[%d] has a root "
                        "at |eps|~%s with certified error %s — too close to "
                        "0 to certify a positive radius at dps=%d; raise "
                        "dps or inspect the polynomial" %
                        (pi, mp.nstr(mag, 8), mp.nstr(err, 8), dps))
                cand = mag - err     # certified-side deflation
                if best is None or cand < best:
                    best = cand
    if best is None:
        return mp.inf
    with mp.workdps(int(dps)):
        return +best


# ---- exact Fraction-polynomial helpers (for desystem_eps_polys) -----------

def _padd(a, b):
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else Fraction(0)) +
            (b[i] if i < len(b) else Fraction(0)) for i in range(n)]


def _pmul(a, b):
    out = [Fraction(0)] * (len(a) + len(b) - 1) if a and b else []
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j, bj in enumerate(b):
            out[i + j] += ai * bj
    return out


def _compose_d_to_eps(p):
    """p(d) -> p(4 - 2*eps) as an exact eps-coefficient list (low->high)."""
    res = [Fraction(0)]
    lin = [Fraction(4), Fraction(-2)]
    for c in reversed(p):
        res = _padd(_pmul(res, lin), [Fraction(c)])
    return res


def desystem_eps_polys(desys, x):
    """Collect the exact eps-polynomials whose roots bound the
    eps-analyticity of A(x, eps) at FIXED transport point x (exact rational:
    int / 'p/q' string / Fraction; floats refused).

    Duck-types the de_load entry classes (same-package engine hook; the
    loaders keep every coefficient an exact Fraction):
      * monomial entries (.num/.den = [(k_eps, k_x, Fraction), ...]):
        denominator collapsed to a poly in eps at x;
      * A_of_y entries (.N/.Q = [(num_coeffs, den_coeffs), ...] in d):
        Q's common-denominator numerator at y=x plus EVERY coefficient
        denominator (both sides), composed through d = 4 - 2*eps;
      * anything else raises NotImplementedError with the class name — pass
        R0 explicitly for such systems (partial-but-correct by design).

    SUPERSET-safe: the list may include apparent (cancelling) singularities,
    so the resulting R0 is conservative — never too large.  Returns a list
    of Fraction coefficient lists (low->high) for eps_analyticity_radius().
    """
    entries = getattr(desys, "_entries", None)
    if not isinstance(entries, dict) or not entries:
        raise ValueError("desystem_eps_polys: no stored entries — not a "
                         "loaded DESystem-like object")
    xq = _exact_fraction(x, "desystem_eps_polys: x")
    polys, seen = [], set()

    def add(poly):
        while poly and poly[-1] == 0:
            poly = poly[:-1]
        if len(poly) <= 1:
            return               # constant/zero: no eps roots
        key = tuple(poly)
        if key not in seen:
            seen.add(key)
            polys.append(poly)

    for (_i, _j), ent in sorted(entries.items()):
        if hasattr(ent, "den") and hasattr(ent, "num"):        # monomial
            acc = {}
            for (ke, kx, c) in ent.den:
                acc[ke] = acc.get(ke, Fraction(0)) + Fraction(c) * xq ** kx
            deg = max(acc)
            add([acc.get(k, Fraction(0)) for k in range(deg + 1)])
        elif hasattr(ent, "N") and hasattr(ent, "Q"):          # A_of_y (d)
            for side in (ent.N, ent.Q):
                for (_num_c, den_c) in side:
                    add(_compose_d_to_eps([Fraction(c) for c in den_c]))
            # Q common-denominator numerator at y = x
            pairs = [([Fraction(c) for c in num_c],
                      [Fraction(c) for c in den_c])
                     for (num_c, den_c) in ent.Q]
            common = [Fraction(0)]
            ypow = Fraction(1)
            for k, (num_c, _den_c) in enumerate(pairs):
                term = [c * ypow for c in num_c]
                for j, (_n2, den2) in enumerate(pairs):
                    if j != k:
                        term = _pmul(term, den2)
                common = _padd(common, term)
                ypow *= xq
            add(_compose_d_to_eps(common))
        else:
            raise NotImplementedError(
                "desystem_eps_polys: entry class %r not supported — pass R0 "
                "explicitly (never silently assume)" % type(ent).__name__)
    return polys


def _parse_node_list(nodes, who):
    """Parse an explicit REAL eps-node list (str preferred — exact 'p/q'
    parsed by mpmathify — or Fraction/mpf/int).  Must be called inside
    mp.workdps.  Positive, real, duplicate-free; returns sorted mpf list."""
    out = []
    for e in nodes:
        if isinstance(e, str):
            try:
                ev = mp.mpmathify(e)
            except Exception as exc:
                raise ValueError("%s: unparseable eps node %r (%s)"
                                 % (who, e, exc))
        elif isinstance(e, Fraction):
            ev = mp.mpf(e.numerator) / mp.mpf(e.denominator)
        else:
            ev = mp.mpc(e)
        if not (mp.im(ev) == 0 and mp.re(ev) > 0):
            raise ValueError("%s: explicit nodes must be real and > 0, "
                             "got %r" % (who, e))
        out.append(mp.mpf(mp.re(ev)))
    out = sorted(out)
    for a, b in zip(out, out[1:]):
        if a == b:
            raise ValueError("%s: duplicate explicit node %s"
                             % (who, mp.nstr(a, 12)))
    return out


def _circle_max(f, radius, n):
    """max |f| sampled at n nodes on |eps|=radius (half-step angle offset so
    the probe nodes never coincide with the fan nodes).  Must be called
    inside mp.workdps."""
    m = mp.mpf(0)
    for j in range(n):
        z = radius * mp.expjpi(mp.mpf(2 * j + 1) / n)
        for v in _f_row(f, z, None, "circle probe"):
            m = max(m, abs(v))
    return m


def _resolve_R0(R0, eps_polys, wp_probe, who):
    """Shared R0 resolution: caller R0 wins; else exact eps_polys roots;
    else REFUSE (certified-layer rule 1)."""
    if R0 is None and eps_polys is None:
        raise ValueError(_R0_REFUSAL % who)
    if R0 is not None:
        return _parse_radius(R0, who + ": R0"), "caller"
    return eps_analyticity_radius(eps_polys, dps=wp_probe), "eps_polys"


def cauchy_laurent_certified(f, kmin, kmax, dps, R0=None, eps_polys=None,
                             r_seed='1e-3', nodes_seed=32, guard=10,
                             max_nodes=8192, max_r_halvings=24,
                             ctrl_nodes=16, ctrl_safety=4,
                             schwarz_real=False):
    """Certified Cauchy Laurent fan: c_k, k = kmin..kmax, of f around eps=0
    with a RUNTIME-chosen (r, M) certified against an explicit analyticity
    radius — no frozen radius/grid, no dps-independent alias floor
    (certified-layer rules 1-3 + 5).

    Requires the caller to know where f stops being analytic:
      R0:        analyticity radius (nearest eps-singularity distance);
                 str/Fraction/mpf/mp.inf.  Caller-supplied R0 wins.
      eps_polys: exact eps-polynomial singularity data (see
                 desystem_eps_polys / eps_analyticity_radius) used to
                 compute R0 at runtime when R0 is None.
    With NEITHER, this function raises ValueError — never silently assume.

    Certification pipeline (all fail-closed):
      1. Control circle r' = sqrt(r*R0) (2r when R0=inf), r < r' < R0.
         M_{r'} = ctrl_safety * (max |f| sampled at ctrl_nodes offset nodes)
         gives the Cauchy coefficient bound |c_n| <= M_{r'}/r'^n and hence
         the certified alias bound on the M-node DFT at radius r:
             alias(c_k) <= M_{r'} (r/r')^M / (r'^k (1 - (r/r')^M)).
         Requires M > kmax-kmin so only upward orders alias (caller
         contract: kmin at or below f's true pole order).
      2. (r, M) chosen at RUNTIME from the seeds (r_seed, nodes_seed —
         "current values as seeds"): M escalates geometrically (x2) to
         max_nodes; on failure r shrinks (x1/2, which also shrinks
         q = r/r' = sqrt(r/R0)) up to max_r_halvings, control circle
         re-probed per r.  Accept when the worst-k alias bound <
         10^-(dps+guard) * max(1, sampled max|f| on the r circle).
         At the caps: EpsFanCertifyError('alias-cap') NAMING
         (bound, tol, r, M, R0).
      3. Working precision = dps + margin(M, r) with margin computed from
         THIS extraction's own conditioning — window pole leverage at the
         belt radius r/2, log10(M) DFT accumulation, sampled magnitude —
         not a constant.
      4. Cross-radius belt: extract at (r, M) AND (r/2, M) — the alias
         structures differ; per-coefficient |c(r) - c(r/2)| must be <
         10^-(dps+guard//2) * max(1, gscale) or
         EpsFanCertifyError('belt').  Cheap, and catches
         conditioning/contamination errors the Cauchy bound cannot see.

    Returns {'coeffs': {k: mpc} (scalar f) / {k: [mpc,...]} (vector f),
             'diag': {..., 'certified': True, 'R0','r','r_belt','nodes',
                      'alias_bound','alias_tol','wp','margin',
                      'belt_agree_digits','belt_min_agree_digits', ...}}.
    Legacy cauchy_laurent is unchanged; this is the additive certified
    entry point (consumers wire when their gates allow).
    """
    kmin = int(kmin)
    kmax = int(kmax)
    if kmax < kmin:
        raise ValueError("cauchy_laurent_certified: kmax < kmin")
    width = kmax - kmin + 1
    dps = int(dps)
    guard = int(guard)
    if guard < 2:
        raise ValueError("cauchy_laurent_certified: guard must be >= 2")
    nodes_seed = int(nodes_seed)
    max_nodes = int(max_nodes)
    ctrl_nodes = int(ctrl_nodes)
    if max_nodes <= width:
        raise ValueError("cauchy_laurent_certified: max_nodes must exceed "
                         "the window width %d" % width)

    wp_probe = max(60, dps // 2)
    with mp.workdps(wp_probe):
        R0v, r0_src = _resolve_R0(R0, eps_polys, wp_probe,
                                  "cauchy_laurent_certified")
        rv = _parse_radius(r_seed, "cauchy_laurent_certified: r_seed")
        if mp.isfinite(R0v) and rv >= R0v / 2:
            rv = R0v / 8          # seed must sit well inside the disk
        tol_rel = mp.mpf(10) ** (-(dps + guard))
        safety = mp.mpf(ctrl_safety)
        accepted = None
        last = {"bound": mp.inf, "tol": tol_rel, "M": nodes_seed}
        n_probe = 0
        for _h in range(max_r_halvings + 1):
            rp = mp.sqrt(rv * R0v) if mp.isfinite(R0v) else 2 * rv
            Mc = _circle_max(f, rp, ctrl_nodes) * safety
            Mr = _circle_max(f, rv, ctrl_nodes)
            n_probe += 2 * ctrl_nodes
            scale_hint = max(mp.mpf(1), Mr)
            tol_abs = tol_rel * scale_hint
            q = rv / rp
            M = max(nodes_seed, width + 1)
            while M <= max_nodes:
                qM = q ** M
                if qM >= 1:
                    bound = mp.inf
                else:
                    bound = (Mc * qM /
                             ((1 - qM) * min(rp ** kmin, rp ** kmax)))
                last = {"bound": bound, "tol": tol_abs, "r": rv, "M": M,
                        "ctrl_r": rp, "ctrl_max": Mc}
                if bound <= tol_abs:
                    accepted = (rv, rp, M, Mc, Mr, bound, tol_abs)
                    break
                M *= 2
            if accepted is not None:
                break
            rv = rv / 2
        if accepted is None:
            raise EpsFanCertifyError(
                'alias-cap',
                "cauchy_laurent_certified: certified alias bound did NOT "
                "reach tolerance at the caps: bound=%s tol=%s r=%s M=%d "
                "R0=%s (max_nodes=%d, max_r_halvings=%d, ctrl_r=%s, "
                "ctrl_max=%s) — raise the caps, shrink the window, or "
                "lower dps" % (
                    mp.nstr(last["bound"], 8), mp.nstr(last["tol"], 8),
                    mp.nstr(last.get("r", rv), 12), last["M"],
                    mp.nstr(R0v, 12), max_nodes, max_r_halvings,
                    mp.nstr(last.get("ctrl_r", mp.mpf(0)), 8),
                    mp.nstr(last.get("ctrl_max", mp.mpf(0)), 8)),
                bound=mp.nstr(last["bound"], 8), tol=mp.nstr(last["tol"], 8),
                r=mp.nstr(last.get("r", rv), 12), M=last["M"],
                R0=mp.nstr(R0v, 12))
        rv, rp, M, Mc, Mr, bound, tol_abs = accepted
        # item 5: working precision from THIS extraction's own conditioning:
        # window pole leverage at the belt radius r/2, log10(M) DFT
        # accumulation, and the sampled magnitude (belt-circle samples grow
        # like Mr * 2^{-kmin} for kmin < 0 — covered explicitly).
        lev = (kmax - kmin) * max(float(-mp.log10(rv / 2)), 0.0)
        mag = max(Mc, Mr * mp.mpf(2) ** max(0, -kmin), mp.mpf(1))
        amp = math.log10(M) + max(float(mp.log10(mag)), 0.0)
        margin = int(math.ceil(lev + amp)) + guard + 10
        r_str = mp.nstr(rv, 15)
        rb_str = mp.nstr(rv / 2, 15)
        rp_str = mp.nstr(rp, 15)
        R0_str = mp.nstr(R0v, 15) if mp.isfinite(R0v) else "inf"
        bound_str = mp.nstr(bound, 8)
        tol_str = mp.nstr(tol_abs, 8)
        Mc_str = mp.nstr(Mc, 8)

    wp = dps + margin
    evals = []
    with mp.workdps(wp):
        c1, n_comp = _extract(f, rv, M, kmin, kmax,
                              schwarz_real=schwarz_real, evals=evals)
        c2, _ = _extract(f, rv / 2, M, kmin, kmax,
                         schwarz_real=schwarz_real, evals=evals)
        gscale = mp.mpf(0)
        for k in c1:
            for v in c1[k]:
                gscale = max(gscale, abs(v))
        belt_scale = max(gscale, mp.mpf(1))
        tol_belt = mp.mpf(10) ** (-(dps + guard // 2)) * belt_scale
        agree = {}
        min_agree = None
        cap_d = float(dps + guard)
        for k in range(kmin, kmax + 1):
            worst = cap_d
            for c in range(n_comp):
                diff = abs(c1[k][c] - c2[k][c])
                if diff > tol_belt:
                    raise EpsFanCertifyError(
                        'belt',
                        "cauchy_laurent_certified: cross-radius belt gate "
                        "FAILED at k=%d (component %d): |c(r)-c(r/2)|=%s > "
                        "tol=%s (r=%s, r/2=%s, M=%d, R0=%s, alias bound "
                        "passed at %s vs %s) — conditioning or non-analytic "
                        "contamination the Cauchy bound cannot see; do NOT "
                        "trust either extraction" % (
                            k, c, mp.nstr(diff, 8), mp.nstr(tol_belt, 8),
                            r_str, rb_str, M, R0_str, bound_str, tol_str),
                        k=k, component=c, diff=mp.nstr(diff, 8),
                        tol=mp.nstr(tol_belt, 8), r=r_str, r_belt=rb_str,
                        M=M, R0=R0_str)
                ad = (cap_d if diff == 0 else
                      max(0.0, min(float(-mp.log10(diff / belt_scale)),
                                   cap_d)))
                worst = min(worst, ad)
            agree[k] = worst
            min_agree = worst if min_agree is None else min(min_agree, worst)
        gscale_str = mp.nstr(gscale, 8)

    with mp.workdps(dps):
        if n_comp == 1:
            coeffs = {k: +c1[k][0] for k in c1}
        else:
            coeffs = {k: [+v for v in c1[k]] for k in c1}

    diag = {
        "mode": "cauchy_certified",
        "certified": True,
        "R0": R0_str,
        "R0_source": r0_src,
        "r": r_str,
        "r_seed": str(r_seed),
        "r_belt": rb_str,
        "nodes": M,
        "nodes_seed": nodes_seed,
        "ctrl_radius": rp_str,
        "ctrl_max": Mc_str,
        "alias_bound": bound_str,
        "alias_tol": tol_str,
        "guard": guard,
        "wp": wp,
        "margin": margin,
        "schwarz_real": bool(schwarz_real),
        "n_evals": len(evals),
        "n_probe_evals": n_probe,
        "belt_agree_digits": agree,
        "belt_min_agree_digits": min_agree,
        "gscale": gscale_str,
        "calibration": ("E4 certified layer: (r,M) "
                        "chosen at runtime against R0=%s via the "
                        "control-circle Cauchy alias bound (%s < %s); "
                        "cross-radius belt gate passed at >= %.1f digits; "
                        "margin=%d computed from this extraction's own "
                        "conditioning (pole leverage + log10(M) + sampled "
                        "magnitude)" % (R0_str, bound_str, tol_str,
                                        min_agree, margin)),
    }
    return {"coeffs": coeffs, "diag": diag}


def _vand_pass(f, kmin, kmax, em, spanv, n, wp, eps_nodes=None):
    """One certified-variant Vandermonde pass: geometric grid of n nodes
    ending at em (span spanv) — or an EXPLICIT node list (eps_nodes,
    re-parsed inside this pass's workdps) — interleaved P fit + spares (the
    legacy production node design), square solve, MEASURED cond_1(V) and
    held-out residual.
    Returns (coeffs dict k->[mpc], n_comp, rel_selferr mpf, cond mpf)."""
    P = kmax - kmin + 1
    with mp.workdps(wp):
        if eps_nodes is not None:
            eps_list = _parse_node_list(eps_nodes,
                                        "vandermonde_laurent_certified")
            n = len(eps_list)
        else:
            eps_list = sorted(em * spanv ** (-mp.mpf(i) / (n - 1))
                              for i in range(n))
        k_spare = n - P
        ver_idx = set(round((i + 1) * n / (k_spare + 1))
                      for i in range(k_spare))
        fit_eps = [e for i, e in enumerate(eps_list) if i not in ver_idx]
        ver_eps = [e for i, e in enumerate(eps_list) if i in ver_idx]
        if len(fit_eps) != P:
            raise ValueError(
                "vandermonde_laurent_certified: interleave gave %d fit "
                "nodes, need square P=%d (raise nodes_seed)"
                % (len(fit_eps), P))
        n_comp = None
        fit_rows, ver_rows = [], []
        for e in fit_eps:
            row = _f_row(f, e, n_comp, "vandermonde_laurent_certified")
            n_comp = len(row) if n_comp is None else n_comp
            fit_rows.append(row)
        for e in ver_eps:
            ver_rows.append(_f_row(f, e, n_comp,
                                   "vandermonde_laurent_certified"))
        V = mp.matrix(P, P)
        for r in range(P):
            for c in range(P):
                V[r, c] = mp.mpc(fit_eps[r]) ** (kmin + c)
        try:
            cond = mp.mnorm(V, 1) * mp.mnorm(mp.inverse(V), 1)
            sols = []
            for comp in range(n_comp):
                b = mp.matrix([fit_rows[r][comp] for r in range(P)])
                sols.append(mp.lu_solve(V, b))
        except ZeroDivisionError as exc:
            # Name the ultra-ill-conditioned
            # path.  mpmath LU raises a bare ZeroDivisionError ("matrix is
            # numerically singular") when the fit nodes COLLAPSE at this
            # working precision (e.g. relative node separation ~1e-18 at
            # the default wp).  This is the
            # cond = inf limit of the condition monitor — raise it under
            # the charter name with the named numbers instead of leaking
            # an anonymous crash.
            sep = (min(abs(fit_eps[i + 1] - fit_eps[i])
                       for i in range(P - 1)) if P > 1 else None)
            raise EpsFanCertifyError(
                'vand-cond',
                "vandermonde_laurent_certified: Vandermonde fit matrix is "
                "NUMERICALLY SINGULAR at wp=%d (%s) — measured cond_1(V) "
                "is effectively infinite, so the condition monitor "
                "cond(V)*10^-wp > 10^-(dps+guard) fails maximally; the "
                "fit nodes collapse at this working precision (min "
                "adjacent fit-node separation %s, n=%d, P=%d, window "
                "%d..%d) — separate the nodes or raise the working "
                "precision" % (
                    wp, exc,
                    "n/a" if sep is None else mp.nstr(sep, 6),
                    n, P, kmin, kmax),
                cond='inf', wp=wp, singular=True,
                min_node_sep=("n/a" if sep is None else mp.nstr(sep, 6)),
                n=n, P=P) from exc
        selferr = mp.mpf(0)
        for vi, e in enumerate(ver_eps):
            for comp in range(n_comp):
                pred = sum(sols[comp][c] * mp.mpc(e) ** (kmin + c)
                           for c in range(P))
                direct = ver_rows[vi][comp]
                scale = max(abs(direct), abs(pred))
                if scale > 0:
                    selferr = max(selferr, abs(pred - direct) / scale)
        coeffs = {kmin + c: [sols[comp][c] for comp in range(n_comp)]
                  for c in range(P)}
        return coeffs, n_comp, selferr, cond


def vandermonde_laurent_certified(f, kmin, kmax, dps, R0=None,
                                  eps_polys=None, nodes_seed=24, guard=10,
                                  eps_max=None, span='10',
                                  grid_shift='0.93', max_wp=None,
                                  eps_nodes1=None, eps_nodes2=None,
                                  k_certify=None):
    """Certified real-grid Vandermonde Laurent window (certified-layer
    rules 4 + 5) — the frozen fixed 24-q fan disease removed.

    R0/eps_polys contract identical to cauchy_laurent_certified (REFUSES
    without singularity data or a caller R0); the grid is clamped inside
    the analyticity disk (eps_max <= R0/10 when R0 is finite).

    Node count is a function of dps SEEDED at 24:
        n1 = max(nodes_seed, P+4, P + ceil(dps/25)),   n2 = ceil(1.5*n1)
    (24 vs 36 at the seed), the second grid SHIFTED (eps_max * grid_shift)
    so the two node sets share no nodes — a genuine two-grid agreement gate,
    per-coefficient |c1 - c2| < 10^-(dps+guard//2) * max(1, gscale) or
    EpsFanCertifyError('vand-grid').

    Working precision from the extraction's OWN conditioning: seed
    wp = dps + 4P + 20 + guard (legacy measured calibration + guard), then
    escalate to wp = dps + guard + ceil(log10 cond_1(V)) + 20 from the
    MEASURED condition number; EpsFanCertifyError('vand-wp-cap') if that
    exceeds max_wp (default 8*dps + 400).  Condition monitor (live on
    every path): if the escalation loop exhausts its
    rounds, the passes are re-run and cond RE-MEASURED at the FINAL wp;
    raise ('vand-cond') if cond(V)*10^-wp > 10^-(dps+guard) on that
    post-escalation state.  On the break path the escalation exit
    condition itself guarantees monitor margin >= 1e-20 (wp >= dps +
    guard + ceil(log10 cond) + 20 with cond measured AT wp) — i.e. the
    wp-escalation subsumes the monitor there; the monitor is the
    fail-closed backstop for the exhaustion path, where cond-estimator
    saturation (measured cond tracking ~10^wp while the true cond outruns
    every round) would otherwise certify silently.  MEASURED note:
    on genuine rational grids
    mpmath's exact rounding makes deep node clusters collapse to EXACT
    singularity rather than saturate the estimator (measured cond is
    order-accurate whenever LU completes), so on genuine grids the
    wp-escalation + the named singular raise subsume the monitor; the
    live monitor is the backstop against estimator-saturation measurement
    pathologies (a forced-saturation witness raises here).  A
    numerically SINGULAR
    fit matrix (node collapse at wp) raises the same charter name with
    cond=inf instead of leaking an unnamed ZeroDivisionError.  Held-out
    verify
    residual (interleaved spares, legacy design) must clear ~10^-(dps-1)
    or EpsFanCertifyError('vand-verify') — window truncation and grid
    pathologies fail CLOSED instead of shipping a wrong digit claim.

    Explicit-node variant (additive): `eps_nodes1`/`eps_nodes2`
    (both or neither) replace the generated grids with caller node lists
    (exact 'p/q' strings preferred — the fixed-eps transported-values
    pattern: dict-backed f + prime 1/q grids).  Enforced fail-closed:
    dedupe + positivity per list, cross-list value disjointness, the
    f(dps)-seeded minimum counts (len1 >= n1 = max(nodes_seed, P+4,
    P+ceil(dps/25)); len2 >= ceil(1.5*len1)), and max node <= R0/10 when
    R0 finite (REFUSES — explicit grids cannot be shrunk).

    `k_certify=(ka, kb)` (optional): enforce the two-grid RAISE gate on the
    subwindow ka..kb only (quoted-window certification; higher fitted
    orders may be tail-absorbers with no certifiable content).  Per-k
    agreement is still reported for ALL orders in diag; the coefficients
    outside ka..kb are returned UNCERTIFIED (diag names the subwindow).
    Default None = full-window gate (legacy).

    Returns {'coeffs', 'diag': {..., 'certified': True, 'cond', 'wp',
             'n_grid1','n_grid2','grid_agree_digits', ...}}.  Legacy
    vandermonde_laurent is unchanged.
    """
    kmin = int(kmin)
    kmax = int(kmax)
    if kmax < kmin:
        raise ValueError("vandermonde_laurent_certified: kmax < kmin")
    P = kmax - kmin + 1
    dps = int(dps)
    guard = int(guard)
    if guard < 2:
        raise ValueError("vandermonde_laurent_certified: guard must be >= 2")
    n1 = max(int(nodes_seed), P + 4, P + int(math.ceil(dps / 25.0)))
    n2 = int(math.ceil(1.5 * n1))
    if max_wp is None:
        max_wp = 8 * dps + 400
    max_wp = int(max_wp)
    explicit = eps_nodes1 is not None or eps_nodes2 is not None
    if explicit and (eps_nodes1 is None or eps_nodes2 is None):
        raise ValueError("vandermonde_laurent_certified: eps_nodes1 and "
                         "eps_nodes2 must be given together (two-grid "
                         "gate needs both explicit grids)")
    if k_certify is not None:
        kc_lo, kc_hi = (int(k_certify[0]), int(k_certify[1]))
        if not (kmin <= kc_lo <= kc_hi <= kmax):
            raise ValueError("vandermonde_laurent_certified: k_certify "
                             "%r outside window [%d, %d]"
                             % (k_certify, kmin, kmax))
    else:
        kc_lo, kc_hi = kmin, kmax

    wp_probe = max(60, dps // 2)
    with mp.workdps(wp_probe):
        R0v, r0_src = _resolve_R0(R0, eps_polys, wp_probe,
                                  "vandermonde_laurent_certified")
        if explicit:
            g1 = _parse_node_list(eps_nodes1,
                                  "vandermonde_laurent_certified: "
                                  "eps_nodes1")
            g2 = _parse_node_list(eps_nodes2,
                                  "vandermonde_laurent_certified: "
                                  "eps_nodes2")
            if set(map(str, g1)) & set(map(str, g2)):
                raise ValueError("vandermonde_laurent_certified: explicit "
                                 "grids SHARE nodes — the two-grid gate "
                                 "needs disjoint node sets")
            if len(g1) < n1:
                raise ValueError(
                    "vandermonde_laurent_certified: eps_nodes1 has %d "
                    "nodes < required n1 = %d (f(dps) node-count rule: "
                    "max(nodes_seed=%d, P+4, P+ceil(dps/25)) at dps=%d) — "
                    "supply more nodes (refine) instead of lowering the "
                    "bar" % (len(g1), n1, int(nodes_seed), dps))
            if len(g2) < int(math.ceil(1.5 * len(g1))):
                raise ValueError(
                    "vandermonde_laurent_certified: eps_nodes2 has %d "
                    "nodes < ceil(1.5*len(eps_nodes1)) = %d (shifted-grid "
                    "rule)" % (len(g2), int(math.ceil(1.5 * len(g1)))))
            n1, n2 = len(g1), len(g2)
            em, em2 = g1[-1], g2[-1]
            if mp.isfinite(R0v) and max(em, em2) > R0v / 10:
                raise ValueError(
                    "vandermonde_laurent_certified: explicit grid reaches "
                    "eps=%s > R0/10=%s (R0=%s) — REFUSING (an explicit "
                    "grid cannot be shrunk; supply nodes inside the "
                    "analyticity disk)" % (mp.nstr(max(em, em2), 10),
                                           mp.nstr(R0v / 10, 10),
                                           mp.nstr(R0v, 10)))
            spanv = None
        else:
            em = (_parse_radius(eps_max, "vandermonde_laurent_certified: "
                                "eps_max") if eps_max is not None
                  else mp.mpf('1e-3'))
            if mp.isfinite(R0v):
                em = min(em, R0v / 10)   # grid well inside the disk
            spanv = _parse_radius(span,
                                  "vandermonde_laurent_certified: span")
            if not (spanv > 1):
                raise ValueError("vandermonde_laurent_certified: span "
                                 "must be > 1")
            shiftv = _parse_radius(grid_shift,
                                   "vandermonde_laurent_certified: "
                                   "grid_shift")
            if not (0 < shiftv < 1):
                raise ValueError("vandermonde_laurent_certified: "
                                 "grid_shift must be in (0,1)")
            em2 = em * shiftv
        R0_str = mp.nstr(R0v, 15) if mp.isfinite(R0v) else "inf"
        em_str = mp.nstr(em, 12)
        em2_str = mp.nstr(em2, 12)

    wp = dps + 4 * P + 20 + guard      # legacy measured calibration + guard
    wp = min(wp, max_wp)
    tol_cond = None
    wp_meas = None            # wp at which the passes were actually run
    for _esc in range(6):
        c1, n_comp, err1, cond1 = _vand_pass(f, kmin, kmax, em, spanv,
                                             n1, wp, eps_nodes=eps_nodes1)
        c2, _, err2, cond2 = _vand_pass(f, kmin, kmax, em2, spanv, n2, wp,
                                        eps_nodes=eps_nodes2)
        wp_meas = wp
        with mp.workdps(wp):
            cond = max(cond1, cond2)
            lcond = max(float(mp.log10(cond)), 0.0) if cond > 0 else 0.0
        wp_req = dps + guard + int(math.ceil(lcond)) + 20
        if wp >= wp_req:
            break
        if wp_req > max_wp:
            raise EpsFanCertifyError(
                'vand-wp-cap',
                "vandermonde_laurent_certified: conditioning-driven working "
                "precision wp_req=%d (measured cond_1(V)=%.3e -> margin "
                "ceil(log10 cond)+%d) exceeds max_wp=%d at dps=%d "
                "(R0=%s, n1=%d, n2=%d) — raise max_wp, shrink the window, "
                "or rescale the grid" % (wp_req, float(mp.mpf(mp.nstr(
                    cond, 8))), guard + 20, max_wp, dps, R0_str, n1, n2),
                cond=mp.nstr(cond, 8), wp=wp, wp_req=wp_req, max_wp=max_wp,
                R0=R0_str, tol=mp.mpf(10) ** (-(dps + guard)))
        wp = wp_req
    if wp != wp_meas:
        # The escalation loop EXHAUSTED its rounds with a final wp
        # assignment it never ran a pass at.  Comparing the STALE cond
        # (measured at the previous, lower wp) against the new wp would
        # make the monitor raise unreachable: every wp assignment satisfies
        # wp = dps+guard+ceil(log10 cond_meas)+20, so the stale margin is
        # >= 1e-20 BY CONSTRUCTION on every exit path (a forced saturating
        # cond of 1.0e394 would certify silently).  Re-extract and
        # RE-MEASURE at the FINAL wp so (a) the
        # monitor below gates post-escalation state — LIVE — and (b) the
        # shipped coefficients come from the wp quoted in diag, not the
        # last pre-escalation pass.
        c1, n_comp, err1, cond1 = _vand_pass(f, kmin, kmax, em, spanv,
                                             n1, wp, eps_nodes=eps_nodes1)
        c2, _, err2, cond2 = _vand_pass(f, kmin, kmax, em2, spanv, n2, wp,
                                        eps_nodes=eps_nodes2)
        with mp.workdps(wp):
            cond = max(cond1, cond2)
    with mp.workdps(wp):
        tol_cond = mp.mpf(10) ** (-(dps + guard))
        monitor = cond * mp.mpf(10) ** (-wp)
        if monitor > tol_cond:
            raise EpsFanCertifyError(
                'vand-cond',
                "vandermonde_laurent_certified: condition monitor FAILED "
                "(LIVE — cond re-measured at the final post-escalation "
                "wp): cond_1(V)*10^-wp = %s > 10^-(dps+guard) = %s "
                "(cond=%s, wp=%d, dps=%d, guard=%d, R0=%s) — the wp "
                "escalation exhausted its rounds with the measured "
                "conditioning still outrunning the working precision" % (
                    mp.nstr(monitor, 8), mp.nstr(tol_cond, 8),
                    mp.nstr(cond, 8), wp, dps, guard, R0_str),
                bound=mp.nstr(monitor, 8), tol=mp.nstr(tol_cond, 8),
                cond=mp.nstr(cond, 8), wp=wp, R0=R0_str)

        # held-out verify gate (fail-closed digit honesty)
        verify_tol = mp.mpf(10) ** (-(dps - 1))
        for tag, err, ngrid in (("grid1", err1, n1), ("grid2", err2, n2)):
            if err > verify_tol:
                raise EpsFanCertifyError(
                    'vand-verify',
                    "vandermonde_laurent_certified: held-out verify "
                    "residual on %s (n=%d) is %s > 10^-(dps-1)=%s — window "
                    "truncation or grid pathology; refusing the digit "
                    "claim (dps=%d, wp=%d, cond=%s, R0=%s)" % (
                        tag, ngrid, mp.nstr(err, 8), mp.nstr(verify_tol, 8),
                        dps, wp, mp.nstr(cond, 8), R0_str),
                    residual=mp.nstr(err, 8), tol=mp.nstr(verify_tol, 8),
                    grid=tag, wp=wp, R0=R0_str)

        # two-grid agreement gate
        gscale = mp.mpf(0)
        for k in c1:
            for v in c1[k]:
                gscale = max(gscale, abs(v))
        gsc = max(gscale, mp.mpf(1))
        tol_grid = mp.mpf(10) ** (-(dps + guard // 2)) * gsc
        agree = {}
        min_agree = None
        cap_d = float(dps + guard)
        for k in range(kmin, kmax + 1):
            worst = cap_d
            for comp in range(n_comp):
                diff = abs(c1[k][comp] - c2[k][comp])
                if diff > tol_grid and kc_lo <= k <= kc_hi:
                    raise EpsFanCertifyError(
                        'vand-grid',
                        "vandermonde_laurent_certified: two-grid agreement "
                        "gate FAILED at k=%d (component %d): |c_n1 - c_n2| "
                        "= %s > tol = %s (n1=%d @ eps_max=%s, n2=%d @ "
                        "eps_max=%s, wp=%d, cond=%s, R0=%s, certify "
                        "window %d..%d) — do NOT trust either grid" % (
                            k, comp, mp.nstr(diff, 8), mp.nstr(tol_grid, 8),
                            n1, em_str, n2, em2_str, wp,
                            mp.nstr(cond, 8), R0_str, kc_lo, kc_hi),
                        k=k, component=comp, diff=mp.nstr(diff, 8),
                        tol=mp.nstr(tol_grid, 8), n1=n1, n2=n2, wp=wp,
                        R0=R0_str)
                ad = (cap_d if diff == 0 else
                      max(0.0, min(float(-mp.log10(diff / gsc)), cap_d)))
                worst = min(worst, ad)
            agree[k] = worst
            min_agree = worst if min_agree is None else min(min_agree, worst)
        self_digits = float(min(
            -mp.log10(err1) if err1 > 0 else mp.mpf(wp),
            -mp.log10(err2) if err2 > 0 else mp.mpf(wp)))
        cond_str = mp.nstr(cond, 8)
        gscale_str = mp.nstr(gscale, 8)

    with mp.workdps(dps):
        if n_comp == 1:
            coeffs = {k: +c1[k][0] for k in c1}
        else:
            coeffs = {k: [+v for v in c1[k]] for k in c1}

    diag = {
        "mode": "vandermonde_certified",
        "certified": True,
        "R0": R0_str,
        "R0_source": r0_src,
        "eps_max1": em_str,
        "eps_max2": em2_str,
        "span": ("explicit-nodes" if explicit else str(span)),
        "n_grid1": n1,
        "n_grid2": n2,
        "nodes_seed": int(nodes_seed),
        "explicit_nodes": bool(explicit),
        "k_certify": [kc_lo, kc_hi],
        "certify_min_agree_digits": (min(agree[k]
                                         for k in range(kc_lo, kc_hi + 1))),
        "guard": guard,
        "wp": wp,
        "max_wp": max_wp,
        "cond": cond_str,
        "self_digits": self_digits,
        "grid_agree_digits": agree,
        "grid_min_agree_digits": min_agree,
        "gscale": gscale_str,
        "design": ("E4 certified variant: "
                   "node count f(dps) seeded at 24 (n1=%d vs n2=%d "
                   "shifted), two-grid per-coefficient gate + cond_1(V) "
                   "monitor; wp=%d escalated from the MEASURED condition "
                   "number, grid clamped inside R0=%s; held-out "
                   "verify residual gated fail-closed" % (n1, n2, wp,
                                                          R0_str)),
    }
    return {"coeffs": coeffs, "diag": diag}


def extract_node_counts(kmin, kmax, dps, nodes_seed=24):
    """The certified-layer node-count rule, exposed so consumers
    can COST a run before sampling (the extractor's on-demand legs may be
    expensive transports):  n1 = max(nodes_seed, P+4, P+ceil(dps/25)),
    n2 = ceil(1.5*n1).  Returns (n1, n2)."""
    P = int(kmax) - int(kmin) + 1
    if P < 1:
        raise ValueError("extract_node_counts: kmax < kmin")
    n1 = max(int(nodes_seed), P + 4, P + int(math.ceil(int(dps) / 25.0)))
    return n1, int(math.ceil(1.5 * n1))


def vandermonde_laurent_extract(f, node_gen, kmin, kmax, dps, R0=None,
                                eps_polys=None, nodes_seed=24, guard=10,
                                max_wp=None, k_certify=None,
                                n1_request=None):
    """E4 EXTRACTOR MODE: the
    engine generates its OWN fan samples on demand — it is the coefficient
    PRODUCER, not just the certifier of somebody else's fit.

    Division of labor:
      * ENGINE decides HOW MANY nodes: (n1, n2) from the certified-layer
        f(dps) rule (extract_node_counts; nodes_seed=24 is the production
        seed).
        Raising dps raises the demanded node count — eps-precision becomes
        genuinely dps-dialable, with the raising gates as the honesty
        backstop.
      * CALLER decides WHICH nodes: node_gen(n, exclude) must return >= n
        exact REAL eps nodes (exact 'p/q' strings preferred; Fraction/mpf
        accepted) disjoint from `exclude` (a tuple of already-issued
        nodes).  This keeps domain node design (e.g. prime 1/q grids that
        a rational-eps transport chain can hit exactly) with the caller.
      * CALLER supplies the samples ON DEMAND through the transport
        callback f(eps) -> value or vector, invoked at the issued nodes
        inside the certified passes (expensive transports should memoize /
        checkpoint per node — the engine may re-call f on wp escalation).

    n1_request: optional dial-up of the grid-1 node count ABOVE the f(dps)
    minimum (n1 = max(rule, n1_request)); requests BELOW the rule are
    refused — the bar cannot be lowered.

    Everything downstream DELEGATES to vandermonde_laurent_certified with
    the generated explicit grids: two-grid agreement gate, condition
    monitor, conditioning-driven wp escalation, held-out verify, R0
    refusal contract, k_certify subwindow — cond/wp machinery UNCHANGED by
    construction (tests assert coefficient equality vs a direct certified
    call on the same grids).  All failure modes remain the certified
    layer's fail-closed raises.

    Returns the vandermonde_laurent_certified result with diag extras:
    'extractor': True, 'n_demanded': [n1, n2] (the f(dps) demand),
    'n_issued': [len(g1), len(g2)] (what node_gen supplied).
    """
    n1_rule, _n2_rule = extract_node_counts(kmin, kmax, dps,
                                            nodes_seed=nodes_seed)
    if n1_request is not None:
        if int(n1_request) < n1_rule:
            raise ValueError(
                "vandermonde_laurent_extract: n1_request=%d is BELOW the "
                "f(dps) node-count rule n1=%d (max(nodes_seed=%d, P+4, "
                "P+ceil(dps/25)) at dps=%d) — the bar cannot be lowered"
                % (int(n1_request), n1_rule, int(nodes_seed), int(dps)))
        n1 = int(n1_request)
    else:
        n1 = n1_rule
    n2 = int(math.ceil(1.5 * n1))
    if not callable(node_gen):
        raise ValueError("vandermonde_laurent_extract: node_gen must be "
                         "callable (n, exclude) -> node list")
    g1 = list(node_gen(n1, ()))
    if len(g1) < n1:
        raise ValueError(
            "vandermonde_laurent_extract: node_gen(%d, ()) returned only "
            "%d nodes — the engine DEMANDS n1=%d (f(dps) rule at dps=%d); "
            "supply more nodes instead of lowering the bar"
            % (n1, len(g1), n1, int(dps)))
    # if node_gen over-supplied grid 1, the shifted-grid rule scales with
    # the ISSUED count (the certified layer enforces len2 >= 1.5*len1)
    n2 = max(n2, int(math.ceil(1.5 * len(g1))))
    g2 = list(node_gen(n2, tuple(g1)))
    if len(g2) < n2:
        raise ValueError(
            "vandermonde_laurent_extract: node_gen(%d, exclude) returned "
            "only %d nodes — the engine DEMANDS n2=ceil(1.5*n1)=%d"
            % (n2, len(g2), n2))
    res = vandermonde_laurent_certified(
        f, kmin, kmax, dps, R0=R0, eps_polys=eps_polys,
        nodes_seed=nodes_seed, guard=guard, max_wp=max_wp,
        eps_nodes1=g1, eps_nodes2=g2, k_certify=k_certify)
    res["diag"]["extractor"] = True
    res["diag"]["n_demanded"] = [n1, n2]
    res["diag"]["n_issued"] = [len(g1), len(g2)]
    return res

"""holonomic_transport.py — house layer over the ore_algebra certified-continuation engine.

EXTERNAL-ENGINE PORT. The engine — all of the certification mathematics — is
ore_algebra (Marc Mezzarobba's `ore_algebra.analytic`): certified numerical
analytic continuation of holonomic ODEs with Arb ball enclosures
(`numerical_transition_matrix`). This module is only the debugged recipe
around it, distilled from production use where transports run this way were
checked against an independent Arb evaluation to 48+ certified digits.

RUN CONTEXT (mandatory): import this module under Sage's Python —

    sage -python your_script.py

with ore_algebra installed into that Sage (`sage -pip install ore_algebra`).
Never `sage -python -c "import ore_algebra"` bare — `sage.all` must be
imported first (bare import fails with a Category ImportError in the
reference build); this module does that.

REFERENCE-BUILD FACTS (SAGEENV_PIN.md in this directory):
  * Reference build: SageMath 10.7 (2025-08-09) + ore_algebra 0.5.
  * algorithm='naive' is hard-wired. In the reference build the engine's
    compiled Cython extension is binary-incompatible with the Sage
    libraries: 'naive' warns "Cython extensions not found" then WORKS
    (pure-Python ball summation, certified, correct), while the default
    (binsplit) path RAISES `TypeError: C variable
    sage.rings.integer._small_primes_table has wrong signature (expected
    int [500], got int [0x1F4])` — hard failure, no fallback. Every
    recorded number in this package assumes the naive path; a build where
    the compiled path works is reported by the smoke as ENV-FACT-CHANGED.

INDEPENDENCE NOTE: this tool is the EXTERNAL certified oracle that the
house transport evaluators (wayfinder, famhar, baller) are spot-checked
against. It must never be merged into them, and none of their code may be
folded in here — see GUIDE.md.
"""
from sage.all import (QQ, ZZ, RR, PolynomialRing, RealBallField,  # noqa: E402 (sage.all FIRST — Category footgun)
                      ComplexBallField, QQbar, I)
from ore_algebra import OreAlgebra  # noqa: E402


def engine_versions():
    """(sage_version, ore_algebra_version) of the live process — record in
    every receipt; all recorded numbers assume 10.7 / 0.5."""
    from sage.version import version as sage_version
    try:
        from importlib.metadata import version as _v
        ore_v = _v("ore_algebra")
    except Exception:
        ore_v = "unknown"
    return str(sage_version), ore_v


def dyadic_point(x, bits=30):
    """Pin x to an exact dyadic k/2^bits (path-point convention: every leg
    then evaluates at IDENTICAL exact points; 30-bit dyadics are
    Float64-exact, so float-coercing consumers see the SAME point)."""
    k = ZZ((RR(2 ** bits) * RR(x)).round())
    return QQ(k) / QQ(2 ** bits)


def mat_vec(M, v):
    n = M.nrows()
    return [sum(M[i, j] * v[j] for j in range(n)) for i in range(n)]


def _to_ball(x, RBF):
    """Re-wrap a transported entry at RBF's precision. The engine returns
    complex balls even on all-real paths; a detoured continuation across a
    real singularity carries a genuinely COMPLEX value, which must survive.
    Exactly-zero imaginary part -> real ball (RBF); anything else -> complex
    ball at the same precision. Never discard an imaginary enclosure."""
    im = x.imag() if hasattr(x, "imag") else None
    if im is not None and not im.is_zero():
        return ComplexBallField(RBF.precision())(x)
    return RBF(x.real()) if hasattr(x, "real") else RBF(x)


def determine_ic_convention(eps=1e-30, prec=430):
    """DETERMINE (never assume) the transition-matrix IC convention —
    derivatives [f, f', f'', ...] vs Taylor coefficients [f, f', f''/2!, ...].

    The probe runs on Du^3 (solutions 1, u, u^2), which separates the two
    conventions in a single matrix-vector product. Re-determine per session
    — the convention is a property of the engine build, not of your script.
    Returns "taylor" or "derivatives".
    """
    Ru = PolynomialRing(QQ, "u")
    u = Ru.gen()
    A = OreAlgebra(Ru, "Du")
    Du = A.gen()
    RBF = RealBallField(prec)
    u0, u1 = QQ(1) / 3, QQ(1) / 2
    Mc = (Du ** 3).numerical_transition_matrix([u0, u1], eps,
                                               algorithm='naive')
    # f = u^2: derivatives at u0: (u0^2, 2u0, 2); Taylor coeffs: (u0^2, 2u0, 1)
    vder = [RBF(u0 ** 2), RBF(2 * u0), RBF(2)]
    vtay = [RBF(u0 ** 2), RBF(2 * u0), RBF(1)]
    tgt = RBF(u1 ** 2)
    od = _to_ball(mat_vec(Mc, vder)[0], RBF)
    ot = _to_ball(mat_vec(Mc, vtay)[0], RBF)
    if abs((ot - tgt).mid()) < 1e-20:
        return "taylor"
    if abs((od - tgt).mid()) < 1e-20:
        return "derivatives"
    raise RuntimeError("IC convention probe failed: %s / %s vs %s"
                       % (od, ot, tgt))


def ic_vector(vals_der, conv):
    """vals_der = [f, f', f'', ...] at the base point; convert to the
    DETERMINED convention."""
    if conv == "derivatives":
        return list(vals_der)
    if conv == "taylor":
        return [vals_der[i] / ZZ(i).factorial() for i in range(len(vals_der))]
    raise ValueError("conv must be 'derivatives' or 'taylor', got %r" % conv)


def plan_path(L, u0, u1, clearance=QQ(1) / 100, detour_imag=QQ(1) / 8):
    """Path law: check the REAL roots (not just rational) of the leading
    coefficient against the segment [u0, u1] (with `clearance` margin); if
    any root sits on the segment, detour through the complex midpoint
    (u0+u1)/2 + I*detour_imag.  Returns (path, roots_on_path)."""
    lead = L.leading_coefficient()
    lead_roots = [r for r, _m in lead.roots(RR)]
    lo = min(float(u0), float(u1)) - float(clearance)
    hi = max(float(u0), float(u1)) + float(clearance)
    onpath = [r for r in lead_roots if lo <= r <= hi]
    if onpath:
        path = [u0, (u0 + u1) / 2 + QQbar(I) * QQ(detour_imag), u1]
    else:
        path = [u0, u1]
    return path, onpath


def transition_matrix(L, path, eps=1e-48):
    """The engine call — algorithm='naive' HARD-WIRED (reference-build law;
    the default path hard-fails in that build, see module docstring)."""
    return L.numerical_transition_matrix(path, eps, algorithm='naive')


def enclosure_too_wide(ball, rel=1e-6):
    """True when the certified radius swamps the value (the measured
    high-rank near-singular failure mode: rank-48 near-singular eps=1e-40
    gave max radius 1.4e18 — an enclosure that certifies nothing).
    Accepts real or complex balls (|mid| is the magnitude either way)."""
    r = float(ball.rad())
    m = float(abs(ball.mid()))
    return r > rel * max(m, 1.0)


def certified_transport(L, u0, u1, ics_der, eps=1e-48, prec=430, conv=None):
    """The debugged production recipe, one call.

    L        : OreAlgebra operator in Du over QQ[u], annihilating f.
    u0, u1   : exact (rational/dyadic) base and target points — pin dyadics
               via dyadic_point() when a float-coercing leg must agree.
    ics_der  : [f(u0), f'(u0), ..., f^(r-1)(u0)] as balls (RBF), r = order.
    Returns dict: value (entry 0 re-wrapped as a ball — REAL ball when its
    imaginary part is exactly zero, COMPLEX ball otherwise; a detoured
    continuation across a real singularity is genuinely complex and the
    imaginary part is never discarded), vector (all transported entries,
    raw engine balls), path, detoured, conv, enclosure_ok.

    Every downstream use must still GATE the value against an independent
    oracle — this function certifies the continuation, not the inputs.
    """
    RBF = RealBallField(prec)
    if conv is None:
        conv = determine_ic_convention(prec=prec)
    path, onpath = plan_path(L, u0, u1)
    M = transition_matrix(L, path, eps)
    v0 = ic_vector(list(ics_der), conv)
    vec = mat_vec(M, v0)
    val = _to_ball(vec[0], RBF)
    return {"value": val,
            "vector": vec,
            "path": path,
            "detoured": bool(onpath),
            "roots_on_path": [str(r) for r in onpath],
            "conv": conv,
            "eps": eps,
            "enclosure_ok": not enclosure_too_wide(val)}

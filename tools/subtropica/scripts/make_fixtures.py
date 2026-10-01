#!/usr/bin/env python3
# make_fixtures.py -- Phase-A fixture generator for subtropica (C13).
#
# Builds the Phase-A fixture set as RAW EULER QUADRUPLES (DESIGN.md RT-x) in the
# JSON schema "subtropica-euler-quad-v1" documented in fixtures/README.md, matching
# src/types.jl EulerIntegrand field-for-field:
#
#   integrand = prefactor * prod_i x_i^(nu_i) * prod_j P_j^(a_j + b_j*eps)
#   over x in [0,inf)^n, FLAT measure dx; all integration variables and all
#   kinematic coefficients POSITIVE (declared-Euclidean, types.jl contract).
#
# Fixture set (task spec, DESIGN.md C13 + Gate A):
#   (1) moller_box_eq41      -- paper Sec 4.1.1 one-loop massive Moller box
#                               (PaperChecks.wl:254, Lstlisting 17), Schwinger/
#                               Feynman parametric rep via pySecDec
#                               LoopIntegralFromGraph, Cheng-Wu gauge x3=1.
#   (2) smirnov_tst2_5var    -- Smirnov-class convergent 5-var period integral
#                               (DESIGN.md Gate A: "eq 4.12 Smirnov period, HF
#                               test/Smirnov tst2-CLASS integrand"). Integrand
#                               DEFINITION extracted from the upstream HF test
#                               files tst0/tst1/tst2 (NOT any stored value) and
#                               the eps-twist exponents VERIFIED symbolically.
#   (3) synth_log2_1var, synth_zeta2_2var, synth_zeta3_beta_1var
#                            -- synthetic-truth POSITIVES with known closed
#                               forms (log2 / zeta(2) / zeta(2)+zeta(3) combos)
#                               derived independently IN THIS SCRIPT with
#                               sympy + mpmath (memory: synthetic-truth
#                               controls mandatory).
#   (4) divergent_xinv_1var  -- one DIVERGENT integrand, flagged
#                               "divergent": true. Phase-A gate (iv): the
#                               pipeline must refuse it LOUDLY
#                               ({"divergent":true} fatal; the in-repo
#                               divergent_xinv_1var fixture is the control).
#
# VALUES POLICY (DESIGN.md RT-iii, enforced structurally by src/verify.jl):
#   Fixture files contain NO result values of any kind. Known-truth values for
#   the synthetic fixtures are written as SEPARATE oracle-run artifacts
#   (schema "subtropica-oracle-run-v1", kind "oracle-run", full provenance stamp)
#   under fixtures/oracle/. verify.jl accepts comparison targets ONLY from such
#   artifacts -- fixture-embedded values are FORBIDDEN as comparison targets.
#
# Conventions transliterated (cite: DESIGN.md C2; SubTropica.wl
# STGetIntegrandData:2975 / STSymanzik:2892):
#   eU = nu_tot - (L+1)*D/2,  eF = -(nu_tot - L*D/2),  D = 4 - 2*eps,
#   prefactor Gamma(-eF)/prod_i Gamma(nu_i);  monomial exponents stored as
#   nu_i(flat measure) = nu_i(propagator power) - 1. The F polynomial is stored
#   in pySecDec's Euclidean-positive arrangement (pySecDec "F" == -F_paper,
#   positive for sb,tb,msq > 0) -- the (-F)^eF sign convention is thereby
#   already resolved to a positive base; recorded in provenance.
#
# Usage (ALWAYS under the house memory cap):
#   bash -c 'ulimit -v 32505856; python3 scripts/make_fixtures.py'                # full set
#   bash -c 'ulimit -v 32505856; python3 scripts/make_fixtures.py --synthetic-only'
#   bash -c 'ulimit -v 32505856; python3 scripts/make_fixtures.py \
#       --oracle fixtures/synth_log2_1var.json --dps 60 --out <dir>'              # mpmath oracle run
#
# The --oracle mode is the live mpmath tanh-sinh oracle adapter used by
# src/verify.jl run_oracle(:mpmath, ...): low-dim (<=2 var), convergent,
# eps-free fixtures only; writes an oracle-run artifact and prints its path.

import argparse
import datetime
import hashlib
import json
import os
import platform
import sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
FIXDIR_DEFAULT = os.path.join(ROOT, "fixtures")
ORACLE_SUBDIR = "oracle"

FIXTURE_SCHEMA = "subtropica-euler-quad-v1"
ORACLE_SCHEMA = "subtropica-oracle-run-v1"

# Upstream sources (READ-ONLY; reference/PROVENANCE.md pins).
# Set SUBTROPICA_UPSTREAM to your SubTropica checkout (required for fixture
# GENERATION; --oracle evaluation of existing fixtures does not read it —
# the refusal fires lazily at first upstream read).
UPSTREAM = os.environ.get("SUBTROPICA_UPSTREAM")


def _require_upstream():
    if not UPSTREAM:
        raise SystemExit("REFUSE: env var SUBTROPICA_UPSTREAM is unset -- point "
                         "it at your SubTropica checkout (read-only)")
    return UPSTREAM


SMIRNOV_DIR = os.path.join(UPSTREAM, "HyperFLINT", "test", "Smirnov") if UPSTREAM else ""
PAPERCHECKS = os.path.join(UPSTREAM, "PaperChecks.wl") if UPSTREAM else ""
WL_PIN = "adac2f722be64337aa095b2b2e7266628b03289b"


# --------------------------------------------------------------------------
# small helpers
# --------------------------------------------------------------------------
def utcnow():
    return datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def qstr(fr):
    """Canonical rational string: gcd=1, q>0, sign on numerator, 'p' if q==1
    (CONTRACTS.md section d content-canonical rule)."""
    fr = Fraction(fr)
    return str(fr.numerator) if fr.denominator == 1 else f"{fr.numerator}/{fr.denominator}"


def epsexp(a, b=0):
    """EpsExp a + b*eps as the JSON pair [qstr(a), qstr(b)] (types.jl EpsExp)."""
    return [qstr(a), qstr(b)]


def base_provenance(extra=None):
    p = {
        "generated_by": "subtropica scripts/make_fixtures.py",
        "argv": sys.argv,
        "timestamp_utc": utcnow(),
        "host": "redacted",
        "user": "redacted",
        "python": platform.python_version(),
        "memcap": "ulimit -v 32505856",
        "upstream_wl_pin": WL_PIN,
    }
    if extra:
        p.update(extra)
    return p


def write_json(path, obj):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        json.dump(obj, f, indent=1)
        f.write("\n")
    print(f"[make_fixtures] wrote {path}")
    return path


# --------------------------------------------------------------------------
# canonical Mma-grammar polynomial emitter (CONTRACTS.md section d)
# --------------------------------------------------------------------------
def mma_poly(expr, allvars, require_positive=True):
    """Serialize a sympy polynomial over `allvars` (integration vars first,
    then kinvars -- types.jl EulerIntegrand ring order) to a canonical
    Mma-grammar string: explicit '*', '^' powers, content-canonical rational
    coefficients, deterministic term order (exponent tuple, descending).
    require_positive enforces the declared-Euclidean all-positive-coefficients
    contract of types.jl EulerIntegrand."""
    import sympy as sp

    p = sp.Poly(sp.expand(expr), *allvars)
    terms = sorted(p.terms(), key=lambda t: t[0], reverse=True)
    if not terms:
        return "0"
    pieces = []
    for mono, coeff in terms:
        c = Fraction(sp.Rational(coeff).p, sp.Rational(coeff).q)
        if require_positive and c <= 0:
            raise ValueError(
                f"mma_poly: non-positive coefficient {c} in {expr} -- violates "
                "the declared-Euclidean positivity contract (types.jl EulerIntegrand)"
            )
        parts = []
        for v, k in zip(allvars, mono):
            if k == 1:
                parts.append(str(v))
            elif k >= 2:
                parts.append(f"{v}^{k}")
        if not parts or c != 1:
            parts.insert(0, qstr(c))
        pieces.append("*".join(parts))
    return " + ".join(pieces)


# --------------------------------------------------------------------------
# fixture writers
# --------------------------------------------------------------------------
def fixture_obj(fid, title, prefactor, nu, polys, variables, kinvars, divergent,
                notes, provenance):
    return {
        "schema": FIXTURE_SCHEMA,
        "id": fid,
        "title": title,
        "divergent": bool(divergent),
        "prefactor": prefactor,
        "nu": nu,
        "polys": polys,
        "vars": variables,
        "kinvars": kinvars,
        "notes": notes,
        "provenance": provenance,
    }


def trivial_prefactor():
    # Prefactor c*eps^k*prod Gamma(g)^p * e^(n*gammaE*eps) (types.jl Prefactor)
    return {"c": "1", "eps_power": 0, "gammas": [], "gammaE_eps": "0"}


# ---------------------------------------------------------------- fixture 1
def build_moller_box(fixdir):
    """Paper Sec 4.1.1 one-loop massive Moller box (PaperChecks.wl:254,
    Lstlisting 17): diag = {{{{1,2},m},{{2,3},0},{{3,4},m},{{1,4},0}},
    {{1,m},{2,m},{3,m},{4,m}}} -- alternating massive/massless edges, all
    external legs on-shell p_i^2 = m^2. Feynman/Schwinger parametric rep via
    pySecDec LoopIntegralFromGraph (the live-tested glue, DESIGN.md reuse
    notes), Cheng-Wu gauge x3 = 1 (last Feynman parameter).

    Mandelstams (all-incoming): (p1+p2)^2 = s, (p2+p3)^2 = t,
    (p1+p3)^2 = u = 4*msq - s - t. Euclidean chart: s = -sb, t = -tb with
    sb, tb > 0, so every polynomial coefficient is positive (checked)."""
    import sympy as sp
    from pySecDec.loop_integral import LoopIntegralFromGraph

    li = LoopIntegralFromGraph(
        internal_lines=[["m", [1, 2]], [0, [2, 3]], ["m", [3, 4]], [0, [4, 1]]],
        external_lines=[["p1", 1], ["p2", 2], ["p3", 3], ["p4", 4]],
        replacement_rules=[
            ("p1*p1", "msq"), ("p2*p2", "msq"), ("p3*p3", "msq"), ("p4*p4", "msq"),
            # (p1+p2)^2 = s ; (p3+p4)^2 = s   (momentum conservation checked below)
            ("p1*p2", "s/2-msq"), ("p3*p4", "s/2-msq"),
            # (p2+p3)^2 = t ; (p1+p4)^2 = t   -- the massive-edge cut carries t
            ("p2*p3", "t/2-msq"), ("p1*p4", "t/2-msq"),
            # (p1+p3)^2 = (p2+p4)^2 = u = 4 msq - s - t
            ("p1*p3", "msq-s/2-t/2"), ("p2*p4", "msq-s/2-t/2"),
            ("m*m", "msq"),
        ],
    )

    x0, x1, x2, x3 = sp.symbols("x0 x1 x2 x3", positive=True)
    s, t, msq, sb, tb, eps = sp.symbols("s t msq sb tb eps", positive=True)

    U = sp.sympify(str(li.U), locals={"x0": x0, "x1": x1, "x2": x2, "x3": x3})
    F = sp.sympify(str(li.F), locals={"x0": x0, "x1": x1, "x2": x2, "x3": x3,
                                      "s": s, "t": t, "msq": msq})
    # exponent conventions (DESIGN.md C2; SubTropica.wl STGetIntegrandData:2975):
    # eU = nu - (L+1)D/2 = 2 eps, eF = -(nu - L D/2) = -2 - eps, D = 4 - 2eps
    eU = sp.sympify(str(li.exponent_U), locals={"eps": eps})
    eF = sp.sympify(str(li.exponent_F), locals={"eps": eps})
    assert sp.simplify(eU - 2 * eps) == 0, f"unexpected exponent_U: {eU}"
    assert sp.simplify(eF - (-2 - eps)) == 0, f"unexpected exponent_F: {eF}"
    assert str(li.Gamma_factor) == "gamma(eps + 2)", str(li.Gamma_factor)
    assert list(li.powerlist) == [1, 1, 1, 1], li.powerlist

    # momentum-conservation consistency: p4 = -(p1+p2+p3) so p4^2 must equal msq
    dot = {"p1p2": s / 2 - msq, "p1p3": msq - s / 2 - t / 2, "p2p3": t / 2 - msq}
    p4sq = 3 * msq + 2 * (dot["p1p2"] + dot["p1p3"] + dot["p2p3"])
    assert sp.simplify(p4sq - msq) == 0, "replacement rules violate momentum conservation"

    # Euclidean chart + Cheng-Wu x3 -> 1 (SubTropica Gauge option semantics;
    # PaperChecks.wl:120 'bare symbols default to [0,Infinity)')
    subs = {s: -sb, t: -tb, x3: 1}
    U1 = sp.expand(U.subs(subs))
    F1 = sp.expand(F.subs(subs))

    intvars = [x0, x1, x2]
    kins = [sb, tb, msq]
    allv = intvars + kins
    U_str = mma_poly(U1, allv)          # raises if any coefficient non-positive
    F_str = mma_poly(F1, allv)

    prefactor = {
        # Gamma(-eF)/prod Gamma(nu_i) = Gamma(2+eps)/1 (DESIGN.md C2)
        "c": "1",
        "eps_power": 0,
        "gammas": [{"arg": epsexp(2, 1), "power": 1}],
        "gammaE_eps": "0",
    }
    nu = [epsexp(0), epsexp(0), epsexp(0)]  # flat-measure exponents nu_i - 1 = 0
    polys = [
        {"poly": U_str, "exp": epsexp(0, 2)},    # U^(2 eps)
        {"poly": F_str, "exp": epsexp(-2, -1)},  # F^(-2-eps), Euclidean-positive base
    ]
    notes = (
        "Paper Sec 4.1.1 / eq 4.1-class one-loop massive Moller box "
        "(PaperChecks.wl:254 Lstlisting 17; upstream pin " + WL_PIN[:8] + "). "
        "pySecDec disteval normalization (int d^D l/(i pi^(D/2)) per loop) -- same "
        "convention as tools/longhand route A (its PARAMETRIC.md benchmarks this very "
        "topology). Euclidean chart sb=-s>0, tb=-t>0, msq>0; u=4*msq+sb+tb is then "
        "ABOVE the u-channel threshold, which is invisible in this parametric rep "
        "(F stays coefficient-positive) but matters for momentum-space oracles. "
        "RAY WARNING: the scaling x0,x2 ~ lam, "
        "x1 ~ lam^2 (soft massless edge {1,4} gauged to 1) gives F ~ lam^2 and "
        "measure lam^4 dlam/lam, i.e. integral ~ int lam^(-1-2*eps) dlam: a "
        "MARGINAL (log) ray at eps=0, regulated only for eps<0. Expect the "
        "Phase-A standalone run (check_divergences hard-ON) to flag this "
        "integrand divergent at eps^0; Gate A(i) may need the Phase-B1 "
        "subtraction route or an approved off-shell deformation. "
        "Flagged in the fixture notes; divergent=false records only that the Euler "
        "integral converges for eps in a strip (dim-reg sense), per the DESIGN "
        "Gate-A classification of this benchmark."
    )
    prov = base_provenance({
        "pysecdec": __import__("pySecDec").__version__,
        "sympy": sp.__version__,
        "graph": "internal [m,[1,2]],[0,[2,3]],[m,[3,4]],[0,[4,1]]; external p1..p4; p_i^2=msq",
        "mandelstams": "(p1+p2)^2=s=-sb, (p2+p3)^2=t=-tb, u=4*msq-s-t",
        "gauge": "Cheng-Wu x3=1 (last Feynman parameter, massless edge {1,4})",
        "exponent_convention": "eU=nu-(L+1)D/2=2eps, eF=-(nu-LD/2)=-2-eps, "
                               "prefactor Gamma(-eF)/prod Gamma(nu_i) "
                               "(DESIGN.md C2; SubTropica.wl STGetIntegrandData:2975)",
        "euclidean_positivity_check": "all monomial coefficients of U|x3=1, F|x3=1 "
                                      "positive over (sb,tb,msq) -- enforced by emitter",
        "paperchecks_citation": "PaperChecks.wl:254 (Lstlisting 17, Sec 4.1.1)",
    })
    return write_json(
        os.path.join(fixdir, "moller_box_eq41.json"),
        fixture_obj("moller_box_eq41",
                    "one-loop massive Moller box (paper Sec 4.1.1), Cheng-Wu x3=1",
                    prefactor, nu, polys,
                    ["x0", "x1", "x2"], ["sb", "tb", "msq"],
                    False, notes, prov))


# ---------------------------------------------------------------- fixture 2
def build_smirnov(fixdir):
    """Smirnov-class convergent 5-var period integral (DESIGN.md Gate A:
    'eq 4.12 Smirnov period, HF test/Smirnov tst2-CLASS integrand').

    Extraction (integrand DEFINITION only -- no stored values touched):
    tst0.txt is the eps^0 base rational integrand; tst1.txt / tst2.txt are its
    eps^1 / eps^2 twist expansions: tstk = base * L^k / k! with
    L = sum_i c_i*Log[f_i]. Reading the Log coefficients off tst1 gives the
    eps-twist exponent c_i of each factor f_i; we VERIFY the reading
    symbolically (sympy): cancel(tst1/tst0) == L and cancel(tst2/tst0) == L^2/2
    with Log[f_i] treated as independent atoms. The Euler quadruple is then
        prod t_i^(m_i + c_i*eps) * prod P_j^(e_j + c_j*eps)
    with the base monomial/poly exponents read off tst0."""
    import sympy as sp

    t1, t2, t3, t4, t5 = sp.symbols("t1 t2 t3 t4 t5", positive=True)
    loc = {"t1": t1, "t2": t2, "t3": t3, "t4": t4, "t5": t5, "log": sp.log}

    def read_tst(name):
        _require_upstream()
        path = os.path.join(SMIRNOV_DIR, name)
        with open(path) as f:
            txt = f.read()
        # Mma -> sympy: only Log[...] brackets occur in these files
        txt = txt.replace("Log[", "log(").replace("]", ")").replace("^", "**")
        return sp.sympify(txt, locals=loc), path

    tst0, p0 = read_tst("tst0.txt")
    tst1, p1 = read_tst("tst1.txt")
    tst2, p2 = read_tst("tst2.txt")

    # factor list and eps-twist coefficients read off tst1's Log coefficients
    P1 = t3 + t1 * t3 + t2 * t3 + t2 * t4 * t5
    P2 = t3 + t1 * t3 + t2 * t3 + t2 * t4 + t2 * t4 * t5
    P3 = t1 * t3 + t2 * t3 + t3 ** 2 + t2 * t4 * t5 + t3 * t4 * t5
    Q1 = t1 + t2 + t3
    Q2 = 1 + t5
    factors = [t1, t2, t3, t4, t5, Q1, Q2, P1, P2, P3]
    twist = [3, 5, -87, -12, 19, -17, 9, 24, 50, -3]
    # base exponents read off tst0: t2^2 t3^4 t4^2 t5 / (Q1^2 Q2 P1^2 P2^2 P3^2)
    base_mono = [0, 2, 4, 2, 1]
    base_poly = {Q1: -2, Q2: -1, P1: -2, P2: -2, P3: -2}
    base = (t2 ** 2 * t3 ** 4 * t4 ** 2 * t5) / (Q1 ** 2 * Q2 * P1 ** 2 * P2 ** 2 * P3 ** 2)
    assert sp.cancel(base - tst0) == 0, "tst0 base integrand mismatch"

    # symbolic verification of the twist reading: replace log(f_i) -> L_i atoms
    Ls = sp.symbols("L0:10")
    logmap = {}
    for f, L in zip(factors, Ls):
        logmap[sp.log(f)] = L
        logmap[sp.log(sp.expand(f))] = L
        logmap[sp.log(sp.factor(f))] = L

    def sub_logs(e):
        out = e
        for k, v in logmap.items():
            out = out.subs(k, v)
        assert not out.atoms(sp.log), f"unmapped Log atoms: {out.atoms(sp.log)}"
        return out

    Lsum = sum(c * L for c, L in zip(twist, Ls))
    r1 = sp.cancel(sp.together(sub_logs(tst1)) / base)
    assert sp.expand(r1 - Lsum) == 0, "tst1 != base * sum(c_i Log f_i) -- twist reading wrong"
    r2 = sp.cancel(sp.together(sub_logs(tst2)) / base)
    assert sp.expand(r2 - sp.expand(Lsum ** 2 / 2)) == 0, \
        "tst2 != base * L^2/2 -- twist reading wrong"
    print("[make_fixtures] smirnov twist exponents VERIFIED against tst1/tst2 (sympy)")

    allv = [t1, t2, t3, t4, t5]
    nu = [epsexp(m, c) for m, c in zip(base_mono, twist[:5])]
    polys = []
    for f, c in zip([Q1, Q2, P1, P2, P3], twist[5:]):
        polys.append({"poly": mma_poly(f, allv), "exp": epsexp(base_poly[f], c)})

    notes = (
        "Smirnov-class convergent 5-variable period integral: base rational "
        "integrand of the upstream HyperFLINT regression family test/Smirnov "
        "(tst0), eps-twisted per the tst1/tst2 expansion structure; the twist "
        "exponents were inferred from the Log coefficients and "
        "VERIFIED symbolically (tst1 == base*L, tst2 == base*L^2/2). DESIGN.md "
        "Gate A names this the 'eq 4.12 Smirnov period' convergent benchmark. "
        "No kinematic symbols; all polynomial coefficients positive. No stored "
        "upstream value was read or copied (RT-iii)."
    )
    prov = base_provenance({
        "sympy": sp.__version__,
        "sources": [
            {"path": p0, "sha256": sha256_file(p0)},
            {"path": p1, "sha256": sha256_file(p1)},
            {"path": p2, "sha256": sha256_file(p2)},
        ],
        "extraction": "integrand definition only; twist exponents "
                      f"{twist} on factors [t1..t5, t1+t2+t3, 1+t5, P1, P2, P3]; "
                      "verified: cancel(tst1/tst0)==L, cancel(tst2/tst0)==L^2/2 "
                      "with Log atoms symbolized",
    })
    return write_json(
        os.path.join(fixdir, "smirnov_tst2_5var.json"),
        fixture_obj("smirnov_tst2_5var",
                    "Smirnov-class 5-var convergent period (HF test/Smirnov tst2-class)",
                    trivial_prefactor(), nu, polys,
                    ["t1", "t2", "t3", "t4", "t5"], [],
                    False, notes, prov))


# ---------------------------------------------------------------- fixture 3
def _mp_str(x, dps):
    import mpmath as mp
    return mp.nstr(x, dps, strip_zeros=False)


def _floored_digits(a, b):
    import mpmath as mp
    if a == b:
        return mp.mp.dps
    scale = max(abs(a), abs(b))
    if scale == 0:
        return mp.mp.dps
    rel = abs(a - b) / scale
    if rel == 0:
        return mp.mp.dps
    return max(int(mp.floor(-mp.log10(rel))), 0)


def oracle_artifact(fixture_id, method, laurent, cross_check, dps, point=None,
                    oracle="mpmath", sanity_only=False):
    return {
        "schema": ORACLE_SCHEMA,
        "kind": "oracle-run",
        "fixture_id": fixture_id,
        "oracle": oracle,
        "sanity_only": bool(sanity_only),
        "dps": dps,
        "point": point or {},
        "laurent": laurent,
        "cross_check": cross_check,
        "provenance": base_provenance({
            "mpmath": __import__("mpmath").__version__,
            "sympy": __import__("sympy").__version__,
            "method": method,
        }),
    }


def build_synthetic(fixdir, dps=60):
    """The three synthetic-truth POSITIVE fixtures with closed forms derived
    independently here (sympy symbolic + mpmath numeric cross-check), plus
    their oracle-run artifacts under fixtures/oracle/. The negative control is
    NOT a fixture: it is produced at verify time by
    SubTropica.mutate_coefficient on a SCRATCH COPY (standing rule: mutation tests
    never run in place)."""
    import sympy as sp
    import mpmath as mp

    x, y, e = sp.symbols("x y e", positive=True)
    outdir = os.path.join(fixdir, ORACLE_SUBDIR)
    paths = []

    def defint(f, v):
        """Definite [0,oo) integral via antiderivative + limits (sympy's
        direct meijerg path returns unevaluated meijerg forms here)."""
        F = sp.integrate(f, v)
        return sp.simplify(sp.limit(F, v, sp.oo) - sp.limit(F, v, 0, dir="+"))

    # ---- (3a) int_0^inf dx / ((1+x)(2+x)) = log 2 --------------------------
    closed = defint(1 / ((1 + x) * (2 + x)), x)
    assert sp.simplify(closed - sp.log(2)) == 0, f"sympy closed form: {closed}"
    mp.mp.dps = dps + 20
    val = mp.log(2)
    q = mp.quad(lambda u: 1 / ((1 + u) * (2 + u)), [0, mp.inf], method="tanh-sinh")
    d_q = _floored_digits(val, q)
    assert d_q >= dps - 5, f"log2 quadrature only matched {d_q} digits"
    fx = fixture_obj(
        "synth_log2_1var", "synthetic truth: int dx/((1+x)(2+x)) = log 2",
        trivial_prefactor(), [epsexp(0)],
        [{"poly": "1 + x", "exp": epsexp(-1)}, {"poly": "2 + x", "exp": epsexp(-1)}],
        ["x"], [], False,
        "Synthetic-truth positive control (Gate A(ii)). Closed form log(2) "
        "derived in make_fixtures.py via sympy.integrate and cross-checked by "
        "mpmath tanh-sinh quadrature; value lives ONLY in the oracle-run "
        "artifact. Verified vs the "
        "live HF binary.",
        base_provenance({"sympy": sp.__version__, "mpmath": mp.__version__,
                         "derivation": "sympy integrate == log(2)"}))
    paths.append(write_json(os.path.join(fixdir, "synth_log2_1var.json"), fx))
    art = oracle_artifact(
        "synth_log2_1var",
        "closed form log(2): sympy.integrate(1/((1+x)*(2+x)),(x,0,oo)) == log(2), "
        "evaluated with mpmath at dps; independent mpmath tanh-sinh quadrature "
        "cross-check recorded",
        {"minorder": 0, "coeffs": [{
            "order": 0, "value_re": _mp_str(val, dps), "value_im": "0",
            "closed_form": "log(2)", "claimed_digits": dps}]},
        {"tanh_sinh_dps": mp.mp.dps, "matched_digits_vs_quadrature": d_q},
        dps)
    paths.append(write_json(os.path.join(outdir, "synth_log2_1var.mpmath.json"), art))

    # ---- (3b) int_0^inf dxdy / ((1+x)(1+y)(1+x+y)) = zeta(2) ---------------
    # inner x-integral -> log(1+y)/y (antiderivative + limits)
    inner = defint(1 / ((1 + x) * (1 + x + y)), x)
    assert sp.simplify(inner - sp.log(1 + y) / y) == 0, f"inner integral: {inner}"
    # outer: substitute u = y/(1+y) (0,1); integrand -> -log(1-u)/u exactly.
    # sympy leaves log(-1/(u-1)) unexpanded (branch caution), so verify the
    # transform by (i) exp-equality (algebraic) and (ii) exact evaluation at
    # u = 1/2 -- on 0<u<1 the difference is a continuous real function with
    # exp(D)=1, so D is a constant in 2*pi*I*Z and D(1/2)=0 pins D == 0.
    uu = sp.Symbol("u", positive=True)
    yu = uu / (1 - uu)
    transformed = sp.simplify((sp.log(1 + y) / (y * (1 + y))).subs(y, yu)
                              * sp.diff(yu, uu))
    D = sp.simplify(transformed * uu + sp.log(1 - uu))   # log(...) + log(1-u)
    assert sp.simplify(sp.exp(D)) == 1, f"substitution check (exp): {D}"
    assert sp.simplify(D.subs(uu, sp.Rational(1, 2))) == 0, f"substitution check at 1/2: {D}"
    outer = sp.integrate(-sp.log(1 - uu) / uu, (uu, 0, 1))
    assert sp.simplify(outer - sp.pi ** 2 / 6) == 0, f"outer integral: {outer}"
    mp.mp.dps = dps + 20
    val2 = mp.zeta(2)
    # nested adaptive tanh-sinh has a well-known ~15-17 digit error floor
    # (tools/longhand PARAMETRIC.md); it is a CROSS-CHECK here, the value
    # itself comes from the symbolic derivation above.
    mp.mp.dps = 30
    q2 = mp.quad(lambda v: mp.quad(lambda u: 1 / ((1 + u) * (1 + v) * (1 + u + v)),
                                   [0, mp.inf], method="tanh-sinh"),
                 [0, mp.inf], method="tanh-sinh")
    mp.mp.dps = dps + 20
    d_q2 = _floored_digits(val2, mp.mpf(q2))
    assert d_q2 >= 12, f"zeta2 nested quadrature only matched {d_q2} digits"
    fx = fixture_obj(
        "synth_zeta2_2var", "synthetic truth: int dxdy/((1+x)(1+y)(1+x+y)) = zeta(2)",
        trivial_prefactor(), [epsexp(0), epsexp(0)],
        [{"poly": "1 + x", "exp": epsexp(-1)}, {"poly": "1 + y", "exp": epsexp(-1)},
         {"poly": "1 + x + y", "exp": epsexp(-1)}],
        ["x", "y"], [], False,
        "Synthetic-truth positive control, genuinely coupled 2-var integrand. "
        "Closed form zeta(2) derived in make_fixtures.py by exact iterated "
        "sympy integration (inner -> log(1+y)/y, outer -> pi^2/6) and "
        "cross-checked by nested mpmath tanh-sinh quadrature (~15-digit nested "
        "floor, recorded honestly).",
        base_provenance({"sympy": sp.__version__, "mpmath": mp.__version__,
                         "derivation": "iterated sympy integrate == pi^2/6"}))
    paths.append(write_json(os.path.join(fixdir, "synth_zeta2_2var.json"), fx))
    art = oracle_artifact(
        "synth_zeta2_2var",
        "closed form zeta(2): exact iterated sympy integration "
        "(inner log(1+y)/y, outer pi^2/6), evaluated with mpmath at dps; nested "
        "tanh-sinh quadrature cross-check at dps=30 (nested adaptive floor)",
        {"minorder": 0, "coeffs": [{
            "order": 0, "value_re": _mp_str(val2, dps), "value_im": "0",
            "closed_form": "zeta(2)", "claimed_digits": dps}]},
        {"tanh_sinh_dps": 30, "matched_digits_vs_quadrature": d_q2},
        dps)
    paths.append(write_json(os.path.join(outdir, "synth_zeta2_2var.mpmath.json"), art))

    # ---- (3c) int_0^inf x^(2eps) (1+x)^(-2-3eps) dx = B(1+2eps, 1+eps) -----
    # Laurent: B = (1+3e)^(-1) * exp(-2 zeta2 e^2 + 6 zeta3 e^3 + O(e^4))
    #        = 1 - 3e + (9 - 2 zeta2) e^2 + (-27 + 6 zeta2 + 6 zeta3) e^3 + ...
    z2s, z3s = sp.zeta(2), sp.zeta(3)
    mine = [sp.Integer(1), sp.Integer(-3), 9 - 2 * z2s, -27 + 6 * z2s + 6 * z3s]
    ser = sp.series(sp.beta(1 + 2 * e, 1 + e), e, 0, 4).removeO()
    for k in range(4):
        diff = sp.expand_func(ser.coeff(e, k) - mine[k])
        diff = sp.simplify(sp.expand(diff.rewrite(sp.zeta)))
        assert diff == 0, f"beta series order {k}: sympy gives {ser.coeff(e, k)}"
    print("[make_fixtures] zeta3_beta Laurent coefficients VERIFIED vs sympy beta series")
    mp.mp.dps = dps + 20
    vals3 = [mp.mpf(1), mp.mpf(-3), 9 - 2 * mp.zeta(2), -27 + 6 * mp.zeta(2) + 6 * mp.zeta(3)]
    # numeric cross-checks: mp.beta and direct quadrature at small rational eps
    checks = []
    for ev in [mp.mpf(1) / 128, -mp.mpf(1) / 128]:
        bval = mp.beta(1 + 2 * ev, 1 + ev)
        ssum = sum(v * ev ** k for k, v in enumerate(vals3))
        mp.mp.dps = 30
        qv = mp.quad(lambda u: u ** (2 * ev) * (1 + u) ** (-2 - 3 * ev),
                     [0, mp.inf], method="tanh-sinh")
        mp.mp.dps = dps + 20
        checks.append({
            "eps": qstr(Fraction(1, 128) if ev > 0 else Fraction(-1, 128)),
            "series_vs_beta_digits": _floored_digits(bval, ssum),      # ~ eps^4 truncation
            "quad_vs_beta_digits": _floored_digits(bval, mp.mpf(qv)),  # quadrature control
        })
        assert checks[-1]["series_vs_beta_digits"] >= 7, checks[-1]
        assert checks[-1]["quad_vs_beta_digits"] >= 20, checks[-1]
    fx = fixture_obj(
        "synth_zeta3_beta_1var",
        "synthetic truth: int x^(2eps)(1+x)^(-2-3eps) dx = B(1+2eps,1+eps); "
        "zeta(3) at eps^3",
        trivial_prefactor(), [epsexp(0, 2)],
        [{"poly": "1 + x", "exp": epsexp(-2, -3)}],
        ["x"], [], False,
        "Synthetic-truth positive control with eps-DEPENDENT exponents "
        "(exercises the EpsExp schema and multi-order Laurent comparison). "
        "Closed form B(1+2eps,1+eps); Laurent coefficients "
        "[1, -3, 9-2zeta(2), -27+6zeta(2)+6zeta(3)] derived from the standard "
        "logGamma psi-expansion and verified against sympy's beta series; "
        "numeric controls vs mp.beta and direct quadrature at eps=+-1/128.",
        base_provenance({"sympy": sp.__version__, "mpmath": mp.__version__,
                         "derivation": "logGamma series; sympy beta series match"}))
    paths.append(write_json(os.path.join(fixdir, "synth_zeta3_beta_1var.json"), fx))
    art = oracle_artifact(
        "synth_zeta3_beta_1var",
        "closed form B(1+2eps,1+eps) Laurent coefficients "
        "[1, -3, 9-2zeta(2), -27+6zeta(2)+6zeta(3)]: derived via the "
        "psi-expansion of logGamma, verified against sympy.series(beta(...)); "
        "per-order values via mpmath.zeta at dps; controls: 4-term series vs "
        "mp.beta at eps=+-1/128 (>=7 digits, eps^4 truncation) and tanh-sinh "
        "quadrature vs mp.beta (>=20 digits at dps=30)",
        {"minorder": 0, "coeffs": [
            {"order": k, "value_re": _mp_str(v, dps), "value_im": "0",
             "closed_form": cf, "claimed_digits": dps}
            for k, (v, cf) in enumerate(zip(vals3, [
                "1", "-3", "9 - 2*zeta(2)", "-27 + 6*zeta(2) + 6*zeta(3)"]))]},
        {"eps_probes": checks},
        dps)
    paths.append(write_json(os.path.join(outdir, "synth_zeta3_beta_1var.mpmath.json"), art))
    return paths


# ---------------------------------------------------------------- fixture 3d
def build_synth7(fixdir):
    """Synthetic 7-var positive control with a factorizing closed form:

        I(eps) = int_{[0,inf)^7} P^(-2-eps) dx,  P = prod_{i=0}^{6} (1+x_i)
               = prod_i int_0^inf (1+x_i)^(-2-eps) dx_i = (1+eps)^(-7).

    P expands to 128 monomials, every coefficient 1 -- a deterministic
    7-variable Euler quadruple that any consumer can regenerate and verify
    independently (the Julia suite reconstructs P and requires exact
    equality). The closed form is symbolic documentation only; no numeric
    value is stored (RT-iii)."""
    import sympy as sp

    xs = sp.symbols("x0 x1 x2 x3 x4 x5 x6", positive=True)
    e = sp.symbols("e", positive=True)
    # per-factor closed form, derived here: int_0^inf (1+x)^(-2-e) dx = 1/(1+e)
    x = xs[0]
    F = sp.integrate((1 + x) ** (-2 - e), x)
    val = sp.simplify(sp.limit(F, x, sp.oo) - F.subs(x, 0))
    assert sp.simplify(val - 1 / (1 + e)) == 0, f"factor integral: {val}"
    P = sp.expand(sp.prod([1 + xi for xi in xs]))
    poly = mma_poly(P, list(xs))
    assert poly.count("+") == 127          # 128 monomials, unit coefficients
    fx = fixture_obj(
        "synth_prod7_7var",
        "synthetic 7-var product control: int P^(-2-eps), P = prod(1+x_i), "
        "= (1+eps)^(-7) by per-variable factorization",
        trivial_prefactor(), [epsexp(0)] * 7,
        [{"poly": poly, "exp": epsexp(-2, -1)}],
        [str(v) for v in xs], [], False,
        "Synthetic 7-variable positive control. The integrand factorizes per "
        "variable, so I(eps) = (1+eps)^(-7) exactly (each factor integrates "
        "to 1/(1+eps); derived in make_fixtures.py via sympy and asserted at "
        "generation time). The closed form is documentation only -- no "
        "numeric value is stored in this file (RT-iii). The Julia suite "
        "reconstructs P = prod(1+x_i) independently and requires exact "
        "equality.",
        base_provenance({"sympy": sp.__version__,
                         "derivation": "int (1+x)^(-2-e) dx = 1/(1+e); "
                                       "product over 7 variables"}))
    return [write_json(os.path.join(fixdir, "synth_prod7_7var.json"), fx)]


# ---------------------------------------------------------------- fixture 4
def build_divergent(fixdir):
    """The mandatory DIVERGENT fixture (Phase-A gate iv):
    int_0^inf dx x^(-1) (1+x)^(-1), log-divergent at x -> 0. The pipeline must
    refuse it LOUDLY: HF standalone with check_divergences hard-ON returns
    {"divergent":true,...} (silent-0-on-divergent control, live-verified;
    policy: SubTropica.wl:12374-12384)."""
    fx = fixture_obj(
        "divergent_xinv_1var",
        "DIVERGENT control: int dx/(x(1+x)) -- pipeline must refuse loudly",
        trivial_prefactor(), [epsexp(-1)],
        [{"poly": "1 + x", "exp": epsexp(-1)}],
        ["x"], [], True,
        "Deliberately divergent (log at x->0; NOT regulated by eps -- the "
        "integrand carries no eps). Phase-A gate (iv): standalone run must "
        "surface {'divergent':true} as FATAL; any numeric or symbolic 'result' "
        "for this fixture is a pipeline bug. verify.jl refuses to compare "
        "anything against it.",
        base_provenance({
            "control_citation": "in-repo control (this fixture + its oracle run): "
                                "1/(x(1+x)) returns {'result':[]} silently without "
                                "check_divergences, loud {'divergent':true} with it",
        }))
    return write_json(os.path.join(fixdir, "divergent_xinv_1var.json"), fx)


# --------------------------------------------------------------------------
# --oracle mode: live mpmath tanh-sinh oracle run on a fixture (verify.jl's
# run_oracle(:mpmath, ...) adapter shells to this)
# --------------------------------------------------------------------------
def run_mpmath_oracle(fixture_path, dps, outdir, point_str=None):
    import sympy as sp
    import mpmath as mp

    with open(fixture_path) as f:
        fx = json.load(f)
    if fx.get("schema") != FIXTURE_SCHEMA:
        raise SystemExit(f"--oracle: {fixture_path} is not a {FIXTURE_SCHEMA} fixture")
    if fx.get("divergent", False):
        raise SystemExit("--oracle: REFUSING divergent fixture "
                         f"{fx['id']} (Phase-A gate iv: divergent inputs are fatal, "
                         "never evaluated)")
    variables = fx["vars"]
    if len(variables) > 2:
        raise SystemExit(f"--oracle: {len(variables)} vars > 2 -- direct tanh-sinh "
                         "oracle is low-dim only (DESIGN.md RT-iii); use "
                         "tools/longhand (hiprec_sectordecomp oracle) / AMFlow-port instead")

    def frac(s):
        return Fraction(s)

    # eps-free only: this mode evaluates the single order eps^0
    for nu in fx["nu"]:
        if frac(nu[1]) != 0:
            raise SystemExit("--oracle: eps-dependent monomial exponent -- generic "
                             "quadrature mode is eps-free/order-0 only; synthetic "
                             "Laurent artifacts come from the dedicated derivations")
    for p in fx["polys"]:
        if frac(p["exp"][1]) != 0:
            raise SystemExit("--oracle: eps-dependent poly exponent -- see above")
    pf = fx["prefactor"]
    if pf["eps_power"] != 0 or frac(pf["gammaE_eps"]) != 0:
        raise SystemExit("--oracle: eps-dependent prefactor unsupported in this mode")

    point = {}
    if fx["kinvars"]:
        if not point_str:
            raise SystemExit(f"--oracle: fixture has kinvars {fx['kinvars']}; "
                             "pass --point 'k1=p/q,k2=p/q,...' (exact rationals only "
                             "-- dyadic/float basepoints poison PSLQ)")
        for item in point_str.split(","):
            k, v = item.split("=")
            point[k.strip()] = Fraction(v.strip())
        missing = [k for k in fx["kinvars"] if k not in point]
        if missing:
            raise SystemExit(f"--oracle: missing point values for {missing}")

    syms = {v: sp.Symbol(v, positive=True) for v in variables}
    ksyms = {k: sp.Rational(point[k]) for k in fx["kinvars"]}
    loc = dict(syms)
    loc.update({k: sp.Symbol(k) for k in fx["kinvars"]})

    def parse_poly(s):
        return sp.sympify(s.replace("^", "**"), locals=loc).subs(
            {loc[k]: ksyms[k] for k in fx["kinvars"]})

    integrand = sp.Integer(1)
    for v, nu in zip(variables, fx["nu"]):
        integrand *= syms[v] ** sp.Rational(str(frac(nu[0])))
    for p in fx["polys"]:
        integrand *= parse_poly(p["poly"]) ** sp.Rational(str(frac(p["exp"][0])))

    mp.mp.dps = dps + 15
    fn = sp.lambdify([syms[v] for v in variables], integrand, modules="mpmath")
    prefval = mp.mpf(Fraction(pf["c"]).numerator) / Fraction(pf["c"]).denominator
    for g in pf["gammas"]:
        if frac(g["arg"][1]) != 0:
            raise SystemExit("--oracle: eps-dependent Gamma argument unsupported here")
        a = frac(g["arg"][0])
        prefval *= mp.gamma(mp.mpf(a.numerator) / a.denominator) ** g["power"]

    t0 = datetime.datetime.now()
    if len(variables) == 1:
        val, err = mp.quad(fn, [0, mp.inf], method="tanh-sinh", error=True)
    else:
        # nested adaptive tanh-sinh: known ~15-17 digit floor (hiprec README);
        # claimed_digits below is set from the measured error estimate, so a
        # 30-digit gate can NEVER pass off this 2-var path -- honest by design.
        val, err = mp.quad(lambda u, v: fn(u, v), [0, mp.inf], [0, mp.inf],
                           method="tanh-sinh", error=True)
    val = prefval * val
    wall = (datetime.datetime.now() - t0).total_seconds()
    claimed = min(dps, _floored_digits(val, val + err * abs(prefval))
                  if err != 0 else dps)
    art = oracle_artifact(
        fx["id"],
        f"direct mpmath tanh-sinh quadrature of the eps^0 integrand at "
        f"working dps={mp.mp.dps}; claimed_digits from mp.quad error estimate "
        f"(wall {wall:.1f}s)",
        {"minorder": 0, "coeffs": [{
            "order": 0, "value_re": _mp_str(val, dps), "value_im": "0",
            "claimed_digits": int(claimed)}]},
        {"quad_error_estimate": _mp_str(mp.mpf(err), 3), "wall_s": wall},
        dps,
        point={k: qstr(v) for k, v in point.items()})
    os.makedirs(outdir, exist_ok=True)
    out = os.path.join(outdir, f"{fx['id']}.mpmath.json")
    write_json(out, art)
    print(out)  # LAST LINE = artifact path (verify.jl run_oracle contract)
    return out


# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fixdir", default=FIXDIR_DEFAULT)
    ap.add_argument("--synthetic-only", action="store_true",
                    help="only fixtures (3)+(4): no pySecDec, no upstream reads")
    ap.add_argument("--dps", type=int, default=60)
    ap.add_argument("--oracle", metavar="FIXTURE_JSON",
                    help="oracle-run mode: mpmath tanh-sinh on a low-dim "
                         "convergent eps-free fixture; writes + prints artifact path")
    ap.add_argument("--out", default=None, help="--oracle output dir")
    ap.add_argument("--point", default=None,
                    help="--oracle kinematic point 'k=p/q,...' (exact rationals)")
    args = ap.parse_args()

    if args.oracle:
        outdir = args.out or os.path.join(args.fixdir, ORACLE_SUBDIR)
        run_mpmath_oracle(args.oracle, args.dps, outdir, args.point)
        return

    build_synthetic(args.fixdir, args.dps)
    build_synth7(args.fixdir)
    build_divergent(args.fixdir)
    if not args.synthetic_only:
        build_smirnov(args.fixdir)
        build_moller_box(args.fixdir)
    print("[make_fixtures] DONE")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
r"""
geotriage.py — per-graph maximal-cut GEOMETRY + VALUE-FITTABILITY auto-triage.

One entry point for the per-graph triage.  The pieces it combines:

  - Baikov → genus via square-free degree (hyperelliptic Riemann–Hurwitz)
  - sector sweep / degree-2 fibration reduction of multi-variable cuts
  - j-invariant, Kodaira fiber configuration, congruence level
  - CM j-table, Sym²-of-sunrise structure, Beauville fiber-config matching
  - Livné CM test via a_p / χ_D(p) point counts
  - loop-by-loop localization verdicts fed in through spec["maxcut"] (proxy mode)
  - the A/B/C value-fittability criterion (period-matrix Jordan structure)

THE DECISION TREE (encoded in closure_route()):

  graph spec → maximal-cut Baikov polynomial → genus/dim of residual variety
   ├─ genus 0 (rational)           → POLYLOG; route = A (value-fit)
   ├─ genus 1 (elliptic curve)     → j, Kodaira config, congruence level N
   │     ├─ j ∈ {0,1728,CM-13}    → CM elliptic; route = A or B (value-fit, ring known)
   │     └─ generic j              → check congruence-modular (Γ₁(N))
   │           ├─ on X₁(N)        → route = B (Eichler integral on Γ₁(N), bootstrappable)
   │           └─ no              → route = C (DE-transport mandatory; e.g. the
   │                                            ice-cream-cone graph)
   ├─ K3 (PF order 3)             → Sym²-of-elliptic check; a_p point-count → CM via Livné
   │     ├─ rank-2 motive         → CM (Livné theorem); route = B
   │     └─ rank≥3                → non-CM; route = C/D (frontier)
   ├─ CY3+ (PF order ≥4)          → route = D (CY transport; new module)
   └─ genus≥2 hyperelliptic       → route = E (Siegel; hard limit)

INPUT graph_spec JSON  (super-set of tools/landau-alphabet/graph_specs/*.json):
  {
    "name": "...",
    "edges": [[v1,v2,"m2"], ...],
    "nodes": {"<vtx>":[<leg>,...], ...},
    "kinematics": {"invariants":[...], "external_masses":{...}, "channels":{...}},
    # OPTIONAL maximal-cut shortcut (used when automatic Baikov is overkill / known):
    "maxcut": {
       "poly": "<sympy string in vars + kinematics>",   # y² = P(var) curve, or F=0 hypersurface
       "vars": ["z"],                                   # residual integration vars
       "modulus": "t",                                  # the kinematic base variable (for j(t))
       "source": "literature ref / derivation script"
    },
    # OPTIONAL kinematic specialization for arithmetic tests (CM / a_p):
    "fiber_point": {"s": -1, "t": 7, ...}
  }

CLI:
    python geotriage.py graph_specs/icc.json --out icc_classify.json

HONESTY (every report carries):
    genus_method        : "exact" | "degree-bound" | "numeric-probe"
    cm_evidence         : list of (p, a_p, χ_D(p)) tuples or j-value match
    pf_order_certified  : bool — did we BUILD the PF op, or just count masters/dim
    caveats             : ["X assumed because Y", ...]
"""
from __future__ import annotations
import sympy as sp
import json, sys, os, argparse
from itertools import combinations
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))

# Reuse the graph-agnostic GraphSpec / Symanzik builder from the Landau Alphabet tool.
sys.path.insert(0, os.path.join(HERE, "..", "landau-alphabet"))
try:
    from landau_alphabet import GraphSpec, symanzik_subgraph    # noqa: E402
    _HAVE_LB = True
except Exception:
    _HAVE_LB = False


# ===========================================================================
#  CONSTANT TABLES  (CM j-invariants, Beauville/Sebbar fiber configs, X₁(N))
# ===========================================================================

# Class-number-1 CM j-invariants (the 13 imaginary-quadratic orders with h(D)=1).
# Verified against Hecke.
CM_J_TABLE = {
    0:                   -3,
    1728:                -4,
    -3375:               -7,
    8000:                -8,
    -32768:              -11,
    54000:               -12,
    287496:              -16,
    -884736:             -19,
    -12288000:           -27,
    16581375:            -28,
    -884736000:          -43,
    -147197952000:       -67,
    -262537412640768000: -163,
}

# Beauville / Sebbar list: the six rational elliptic surfaces with 4 cusps,
# Σnᵢ = 12 (extremal, semistable).  fiber-config (sorted) → congruence subgroup.
BEAUVILLE_4CUSP = {
    (1, 1, 1, 9): ("Gamma_0(9)∩Gamma_1(3)", 9),
    (1, 1, 2, 8): ("Gamma_0(8)∩Gamma_1(4)", 8),
    (1, 1, 5, 5): ("Gamma_1(5)",            5),
    (1, 2, 3, 6): ("Gamma_1(6)",            6),
    (2, 2, 4, 4): ("Gamma_1(4)∩Gamma(2)",   4),
    (3, 3, 3, 3): ("Gamma(3)",              3),
}
# 3-cusp extremal (Σ=12, with one I*-type / additive — incomplete list; the
# configurations covered here):
BEAUVILLE_3CUSP = {
    (1, 4, 1): ("Gamma_0(4)", 4),    # I1* + I4 + I1 type config (gg→Zγ deformed)
    (2, 2, 2): ("Gamma(2)",   2),
}

# X₁(N) is genus-0 for N ∈ {1..10, 12}.  For these N a one-parameter family
# whose j(t) is a rational function of degree = [PSL₂(ℤ):Γ₁(N)] is congruence
# on the nose (Hauptmodul exists).
X1N_GENUS0 = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12}
PSL2_INDEX_GAMMA1 = {1: 1, 2: 3, 3: 4, 4: 6, 5: 12, 6: 12, 7: 24, 8: 24,
                     9: 36, 10: 36, 12: 48}


# ===========================================================================
#  1.  baikov_maxcut(graph_spec)  →  Symanzik / maximal-cut polynomial
# ===========================================================================

def _load_spec(path_or_dict):
    if isinstance(path_or_dict, dict):
        return dict(path_or_dict)
    with open(path_or_dict) as f:
        return json.load(f)


def baikov_maxcut(spec: dict):
    """
    Return the maximal-cut polynomial data of the TOP sector.

    Two routes:
      (i)  spec["maxcut"] explicitly given  →  parse & return verbatim
           (used for graphs whose loop-by-loop localization was already
            computed elsewhere; the proxy mode).
      (ii) generic: build the second-Symanzik F via GraphSpec/symanzik_subgraph
           on ALL edges; the maximal-cut variety is {F=0} ⊂ ℙ^{n-1} in
           Feynman parameters.  (For banana / sunrise / vertex this IS the
           Baikov maximal-cut variety up to birational equivalence.)

    Returns dict
       { 'P': sympy expr,  'vars': [Symbol,...],  'modulus': Symbol|None,
         'kin': [Symbol,...], 'method': str, 'caveats': [...] }
    """
    caveats = []
    kin_names = spec.get("kinematics", {}).get("invariants", [])
    kin = [sp.Symbol(s) for s in kin_names]

    # ---- (i) explicit maxcut polynomial ---------------------------------
    if "maxcut" in spec:
        mc = spec["maxcut"]
        loc = {s: sp.Symbol(s) for s in (mc.get("vars", []) + kin_names
                                         + ([mc["modulus"]] if mc.get("modulus") else []))}
        P = sp.sympify(mc["poly"], locals=loc)
        vars_ = [loc[v] for v in mc["vars"]]
        mod = loc.get(mc.get("modulus")) if mc.get("modulus") else None
        caveats.append(f"maxcut polynomial taken from spec (source: "
                       f"{mc.get('source','user-supplied')}); not re-derived here")
        return {'P': P, 'vars': vars_, 'modulus': mod, 'kin': kin,
                'method': "spec-supplied", 'caveats': caveats}

    # ---- (ii) generic Symanzik route ------------------------------------
    if not _HAVE_LB:
        raise RuntimeError("landau_alphabet.GraphSpec unavailable and no spec['maxcut'] given")
    gs = GraphSpec(spec)
    top = symanzik_subgraph(gs, gs.ALLE)
    if top is None:
        raise RuntimeError("symanzik_subgraph returned None on top sector")
    F = top['F']; U = top['U']; L = top['L']
    fp = [gs.ev[e] for e in gs.ALLE]
    # affine chart: set last Feynman param = 1
    chart = {fp[-1]: sp.Integer(1)}
    Fa = sp.expand(F.subs(chart))
    vars_ = fp[:-1]
    # pick a kinematic modulus (the first non-mass invariant) for j(t) etc.
    # IMPORTANT: use gs.invariants (which carry real=True) so the symbol
    # actually appears in F.free_symbols — sp.Symbol('t') != sp.Symbol('t',real=True).
    kin = list(gs.invariants)
    mod = next((k for k in kin if str(k) not in ("m2", "m", "M")), kin[0] if kin else None)
    caveats.append("maximal-cut variety modelled as {F=0} in Feynman parameters "
                   "(2nd-Symanzik); birational to Baikov maximal cut for "
                   "banana/vertex topologies.")
    return {'P': Fa, 'vars': vars_, 'modulus': mod, 'kin': kin,
            'L': L, 'n_edges': len(fp), 'n_verts': len(gs.verts), 'U': U,
            'method': "symanzik-F", 'caveats': caveats}


# ===========================================================================
#  2.  classify_variety(P, vars)  →  {dim, genus, pf_order, type}
# ===========================================================================

def _squarefree_degree(P, var):
    """Square-free degree in one variable: product of factors with odd mult."""
    num = sp.together(P).as_numer_denom()[0]
    fl = sp.factor_list(num, var)
    odd = [f for f, m in fl[1] if m % 2 == 1
           and var in sp.sympify(f).free_symbols]
    sqf = sp.Mul(*odd) if odd else sp.Integer(1)
    dsq = sp.Poly(sqf, var).degree() if odd else 0
    return dsq, sqf


def _genus_from_sqfree(dsq):
    """y² = P(x) genus from square-free degree (hyperelliptic Riemann–Hurwitz)."""
    if dsq <= 2:
        return 0
    return (dsq - 1) // 2


def classify_variety(mc: dict):
    """
    Determine {dim, genus, pf_order, type} of the maximal-cut variety.

    Strategy:
      *  1 residual var: y² = P(z) → genus from square-free degree (EXACT for
         hyperelliptic).
      *  2 residual vars: try fibration in var₁; if
         square-free deg ≤ 2 in var₁, curve localizes on disc_{var₁} → recurse.
         Else: 2-dim variety = K3 candidate (PF order 3).
      *  ≥3 residual vars: CY (n−1)-fold; PF order = #vars + 1 (banana ladder).

    pf_order_certified = False unless we actually built the operator (we don't
    here — the master count is the proxy, flagged honestly).
    """
    P, vars_ = mc['P'], list(mc['vars'])
    caveats = list(mc.get('caveats', []))
    out = {'pf_order_certified': False, 'caveats': caveats}

    # ---- Symanzik-F hypersurface: {F=0} ⊂ ℙ^{n-1} has dim = n−2. -----------
    # For BANANA topologies (2 vertices, n=L+1 edges) this is the CY (L−1)-fold
    # ladder; PF order = L.  For non-banana
    # multi-var Symanzik we attempt a deg-≤2 fibration reduction below; if that
    # fails, the user should supply spec['maxcut'] (loop-by-loop localization).
    if mc.get('method') == 'symanzik-F' and mc.get('n_verts') == 2:
        n = mc['n_edges']; L = mc.get('L', n - 1)
        cy_dim = n - 2
        if cy_dim <= 1:
            pass  # sunrise: fall through to the curve / fibration logic below
        else:
            tp = "K3" if cy_dim == 2 else ("CY3" if cy_dim == 3 else f"CY{cy_dim}")
            out.update(dim=cy_dim, type=tp, genus=None, pf_order=L,
                       genus_method="degree-bound",
                       hypersurface_degree=sp.Poly(sp.expand(P), *vars_).total_degree())
            caveats.append(f"{tp} classed by banana-ladder rule (2-vertex graph, "
                           f"{n} edges ⇒ CY{cy_dim}, PF order {L}); PF op NOT built.")
            return out

    dim = len(vars_)
    out['dim'] = dim

    # ---------- 0-dim (everything localized) -----------------------------
    if dim == 0:
        out.update(type="polylog", genus=0, pf_order=1,
                   genus_method="exact", squarefree_part=str(sp.factor(P)))
        return out

    # ---------- 1 residual var: hyperelliptic curve y²=P -----------------
    if dim == 1:
        z = vars_[0]
        dsq, sqf = _squarefree_degree(P, z)
        g = _genus_from_sqfree(dsq)
        if g == 0:
            tp = "polylog"
        elif g == 1:
            tp = "elliptic"
        elif g >= 2:
            tp = "hyperelliptic"
        out.update(type=tp, genus=g, pf_order=max(1, g + 1),
                   genus_method="exact",
                   squarefree_degree=dsq,
                   squarefree_part=str(sp.factor(sqf)))
        return out

    # ---------- 2 residual vars ------------------------------------------
    if dim == 2:
        u, v = vars_
        # First try: degree-2 fibration → reduce to curve.
        for (fib, base) in ((u, v), (v, u)):
            try:
                pf = sp.Poly(sp.expand(P), fib)
            except sp.PolynomialError:
                continue
            if pf.degree() <= 2:
                # disc_fib(P) is the residual curve in `base`
                disc = sp.factor(sp.discriminant(pf, fib)) if pf.degree() == 2 \
                       else sp.factor(pf.coeff_monomial(1)**2)  # linear → trivial
                dsq, sqf = _squarefree_degree(disc, base)
                g = _genus_from_sqfree(dsq)
                tp = ("polylog" if g == 0 else "elliptic" if g == 1 else "hyperelliptic")
                out.update(type=tp, genus=g, pf_order=max(1, g + 1),
                           genus_method="exact",
                           reduction=f"disc_{fib}",
                           squarefree_degree=dsq,
                           squarefree_part=str(sp.factor(sqf)),
                           reduced_vars=[str(base)])
                # stash the reduced 1-var curve for elliptic_invariants()
                out['_reduced_P'] = sqf
                out['_reduced_var'] = base
                return out
        # Neither fibration is ≤2: genuine surface.  For a degree-d hypersurface
        # in ℙ² (3 Feynman params, one set to 1): d=3 → cubic = elliptic curve
        # (sunrise!); d=4 → quartic K3.
        deg = sp.Poly(sp.expand(P), *vars_).total_degree()
        n_fp = mc.get('n_edges', dim + 1)
        if n_fp == 3 or deg == 3:
            # ternary cubic ⇒ elliptic curve (genus 1).  Reduce via disc-in-one-var
            # to get a binary quartic.
            try:
                pf = sp.Poly(sp.expand(P), v)
                if pf.degree() == 2:
                    A2, A1, A0 = pf.all_coeffs()
                    quartic = sp.expand(A1**2 - 4 * A2 * A0)
                else:
                    quartic = sp.discriminant(pf, v)
                dsq, sqf = _squarefree_degree(quartic, u)
                g = _genus_from_sqfree(dsq)
            except Exception:
                g, dsq, sqf = 1, 3, P
            out.update(type="elliptic", genus=1, pf_order=2,
                       genus_method="exact" if dsq in (3, 4) else "degree-bound",
                       squarefree_degree=dsq,
                       squarefree_part=str(sp.factor(sqf)),
                       reduction="ternary-cubic→binary-quartic")
            out['_reduced_P'] = sqf
            out['_reduced_var'] = u
            return out
        # degree-4 in ℙ³ chart → K3; otherwise generic surface
        out.update(type="K3", genus=None, pf_order=3,
                   genus_method="degree-bound",
                   hypersurface_degree=deg)
        caveats.append("2-dim residual variety classed as K3 by hypersurface "
                       "degree; PF order NOT certified (master-count proxy).")
        return out

    # ---------- ≥3 residual vars: CY ladder ------------------------------
    # Banana ladder fact: l-loop banana → CY (l−1)-fold, PF order = l.
    # We use n_vars = l−1 (after chart) so PF order = n_vars + 1.
    cy_dim = dim
    pf = dim + 1
    if cy_dim == 2:
        tp = "K3"
    elif cy_dim == 3:
        tp = "CY3"
    else:
        tp = f"CY{cy_dim}"
    out.update(type=tp, genus=None, pf_order=pf,
               genus_method="degree-bound", hypersurface_dim=cy_dim)
    caveats.append(f"{tp} classed by residual-var count (={dim}); PF order = "
                   f"#vars+1 (banana-ladder rule), NOT certified.")
    return out


# ===========================================================================
#  3.  elliptic_invariants(P)  →  {j, disc, kodaira, level_N, on_X1N}
# ===========================================================================

def _binary_quartic_invariants(P, x):
    """Return (I, J, Δ_bq, j) of y² = P(x) via the binary-quartic invariants
       (validated: equal-mass sunrise → Γ₁(6))."""
    Qp = sp.Poly(sp.expand(P), x)
    deg = Qp.degree()
    cc = [Qp.nth(k) for k in (4, 3, 2, 1, 0)]
    a, b, c, d, e = cc
    I = sp.factor(sp.expand(12*a*e - 3*b*d + c**2))
    J = sp.factor(sp.expand(72*a*c*e + 9*b*c*d - 27*a*d**2 - 27*e*b**2 - 2*c**3))
    Disc = sp.factor(sp.expand(4*I**3 - J**2))   # ∝ discriminant of the quartic
    j = sp.cancel(1728 * sp.factor(4*I**3) / Disc)
    return I, J, Disc, j


def elliptic_invariants(P, var, modulus=None):
    """
    Weierstrass-model invariants of the genus-1 maximal-cut curve y² = P(var).

    Returns {j, disc, kodaira:{loc:order}, level_N, on_X1N, congruence_group,
             genus_method, caveats}.

    Kodaira fiber config = orders of poles of j(modulus) (semistable Iₙ fibers);
    level via Beauville/Sebbar match on the sorted fiber signature.
    """
    caveats = []
    I, J, Disc, j = _binary_quartic_invariants(P, var)
    out = {'j': str(sp.factor(j)), 'disc': str(Disc), 'I': str(I), 'J': str(J)}

    if modulus is None or modulus not in sp.sympify(j).free_symbols:
        # j is a constant → single curve, not a family
        out.update(kodaira={}, level_N=None, on_X1N=None,
                   congruence_group=None,
                   caveats=["j is modulus-independent; no fiber config computed"])
        return out

    t = modulus
    jn, jd = sp.fraction(sp.together(j))
    config = {}
    for fac, mult in sp.factor_list(sp.Poly(jd, t))[1]:
        fs = sp.factor(fac)
        if t in fs.free_symbols:
            # only an Iₙ cusp if numerator does NOT also vanish there
            try:
                rts = sp.roots(sp.Poly(fs, t))
                cancels = any(sp.simplify(jn.subs(t, r)) == 0 for r in rts)
            except Exception:
                cancels = False
            if not cancels:
                config[str(fs)] = int(mult)
    # behaviour at ∞
    try:
        degn = sp.degree(sp.Poly(jn, t)); degd = sp.degree(sp.Poly(jd, t))
        if degn > degd:
            config["oo"] = int(degn - degd)
        degj = max(degn, degd)
    except Exception:
        degj = None

    sig = tuple(sorted(config.values()))
    sumn = sum(sig)
    out['kodaira'] = config
    out['fiber_signature'] = list(sig)
    out['sum_fiber_orders'] = sumn
    out['deg_j'] = degj

    # Congruence test (Riemann–Hurwitz):
    level_N = None; group = None; on_X1N = False
    if sumn == 12 and sig in BEAUVILLE_4CUSP:
        group, level_N = BEAUVILLE_4CUSP[sig]
        on_X1N = (level_N in X1N_GENUS0)
    elif sumn == 12 and len(sig) == 3 and sig in BEAUVILLE_3CUSP:
        group, level_N = BEAUVILLE_3CUSP[sig]
        on_X1N = (level_N in X1N_GENUS0)
    elif degj in PSL2_INDEX_GAMMA1.values():
        # heuristic: deg(j) matches some [PSL₂:Γ₁(N)] for genus-0 N
        cands = [N for N, idx in PSL2_INDEX_GAMMA1.items() if idx == degj]
        if cands:
            level_N = min(cands); on_X1N = True
            group = f"index-{degj} genus-0 (candidate N∈{cands})"
            caveats.append(f"level_N={level_N} from deg(j)={degj} index match only "
                           f"(NOT fiber-config certified)")
    else:
        caveats.append(f"fiber config {sig} (Σ={sumn}) is NOT in the Beauville/"
                       f"Sebbar 4-cusp list and deg(j)={degj} matches no genus-0 "
                       f"Γ₁(N) index ⇒ likely NON-congruence / non-modular family.")

    out.update(level_N=level_N, on_X1N=on_X1N, congruence_group=group,
               caveats=caveats)
    return out


# ===========================================================================
#  4.  cm_check(curve_or_K3)  →  {is_cm, discriminant_D, evidence}
# ===========================================================================

def _kronecker(a, p):
    """(a|p) for odd prime p."""
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def _primes_upto(n):
    s = bytearray([1]) * (n + 1); s[0] = s[1] = 0
    for i in range(2, int(n**0.5) + 1):
        if s[i]:
            s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(2, n + 1) if s[i]]


def _ap_elliptic_over_Fp(a4, a6, p):
    """a_p = p + 1 − #E(𝔽_p) for y² = x³ + a4 x + a6 (naive O(p) point-count)."""
    a4 %= p; a6 %= p; cnt = 1  # point at infinity
    for x in range(p):
        rhs = (x*x % p * x + a4*x + a6) % p
        if rhs == 0:
            cnt += 1
        elif pow(rhs, (p - 1)//2, p) == 1:
            cnt += 2
    return p + 1 - cnt


def _elliptic_model_from_j(j):
    """
    Short Weierstrass model (a4, a6) over ℚ with the given j ∈ ℚ — SOME twist
    of the curve (from j alone the curve is only defined up to twist).  That is
    enough for the CM test: twisting multiplies a_p by a root of unity, so the
    VANISHING pattern a_p = 0 ⟺ χ_D(p) = −1 is twist-invariant, and that
    pattern is all the callers consume.  Self-verifies j(a4, a6) == j exactly.
    """
    j = Fraction(j)
    if j == 0:
        a4, a6 = Fraction(0), Fraction(1)      # y² = x³ + 1
    elif j == 1728:
        a4, a6 = Fraction(-1), Fraction(0)     # y² = x³ − x
    else:
        k = 1728 - j
        a4, a6 = 3 * j * k, 2 * j * k * k      # standard j-line model
    disc = 4 * a4**3 + 27 * a6**2
    if disc == 0:
        raise ValueError(f"singular j-line model at j={j}")
    if a4 != 0 and Fraction(1728) * 4 * a4**3 / disc != j:
        raise RuntimeError(f"j-line model self-check failed at j={j}")
    return a4, a6


def _ap_from_j(j, p):
    """
    a_p of (a twist of) the elliptic curve with j-invariant j ∈ ℚ over 𝔽_p,
    via the wired naive point counter.  Returns None when the model has bad
    reduction at p (p | denominator of j, or singular reduction).  Only the
    twist-invariant content (a_p = 0 or not) is meaningful to callers; the
    Hasse bound a_p² ≤ 4p is asserted as an internal sanity check.
    """
    j = Fraction(j)
    if p < 5 or j.denominator % p == 0:
        return None
    jp = (j.numerator * pow(j.denominator, p - 2, p)) % p
    if jp == 0:
        a4, a6 = 0, 1
    elif jp == 1728 % p:
        a4, a6 = p - 1, 0
    else:
        k = (1728 - jp) % p
        a4 = (3 * jp * k) % p
        a6 = (2 * jp * k * k) % p
    if (4 * a4 * a4 * a4 + 27 * a6 * a6) % p == 0:
        return None                      # singular reduction of the model
    ap = _ap_elliptic_over_Fp(a4, a6, p)
    if ap * ap > 4 * p:
        raise RuntimeError(f"Hasse bound violated at p={p} (a_p={ap}) — counter bug")
    return ap


def _supersingular_pattern(aps, D, min_primes=8):
    """
    (p, a_p, χ_D(p)) triples over the χ_D-unramified primes in `aps`, plus
    whether the CM supersingular pattern  a_p = 0 ⟺ χ_D(p) = −1  holds on
    ALL of them (Deuring; requires ≥ min_primes usable primes to count).
    """
    trip = [(p, aps[p], _kronecker(D, p)) for p in sorted(aps)
            if _kronecker(D, p) != 0]
    ok = len(trip) >= min_primes and \
        all((a == 0) == (c == -1) for _, a, c in trip)
    return trip, ok


def cm_check_elliptic(j_expr, modulus=None, fiber_point=None, pmax=100):
    """
    For an elliptic curve / family:
      (a) if j is a CONSTANT rational, decide CM by CM_J_TABLE — complete for
          j ∈ ℚ (a rational j is CM iff it is one of the 13 class-number-1
          values; higher-class-number CM j's are irrational) — and CORROBORATE
          by point count with the wired a_p counter:
          a_p = 0 ⟺ χ_D(p) = −1 for the primes 5 ≤ p ≤ pmax (twist-invariant).
          `evidence` carries the (p, a_p, χ_D(p)) triples the module header
          promises; `certified` records whether the pattern check succeeded.
      (b) else specialize the modulus at fiber_point (if given) and run (a);
      (c) else: modulus-dependent j and no fiber_point → no verdict.
    """
    j = sp.nsimplify(sp.sympify(j_expr), rational=True) \
        if not hasattr(j_expr, 'free_symbols') else sp.sympify(j_expr)
    ev = []
    # try specialization
    if j.free_symbols and fiber_point:
        # match by NAME (GraphSpec symbols carry real=True assumptions)
        fs_by_name = {str(s): s for s in j.free_symbols}
        subs = {fs_by_name[k]: sp.nsimplify(v) for k, v in fiber_point.items()
                if k in fs_by_name}
        if subs:
            j = sp.simplify(j.subs(subs))
            ev.append(f"j specialized at {fiber_point} → {j}")
    if j.free_symbols:
        return {'is_cm': None, 'discriminant_D': None, 'certified': False,
                'evidence': ev + ["j is modulus-dependent and no fiber_point given; "
                                  "CM test skipped (a generic fiber is non-CM)."]}
    if not j.is_rational:
        return {'is_cm': None, 'discriminant_D': None, 'certified': False,
                'evidence': ev + [f"j = {j} is not a rational number at this fiber "
                                  "(cusp/singular specialization?); CM-13 table "
                                  "does not apply."]}
    jr = sp.Rational(j)
    jq = Fraction(int(jr.p), int(jr.q))
    aps = {}
    for p in _primes_upto(pmax):
        a_ = _ap_from_j(jq, p)
        if a_ is not None:
            aps[p] = a_
    if jq.denominator == 1 and int(jq) in CM_J_TABLE:
        D = CM_J_TABLE[int(jq)]
        trip, ok = _supersingular_pattern(aps, D)
        out = {'is_cm': True, 'discriminant_D': D, 'certified': bool(ok),
               'method': "CM-13 table + a_p ⟺ χ_D supersingular pattern "
                         "(wired point counter)",
               'evidence': ev + [f"j = {int(jq)} ∈ CM-13 table (h(D)=1, D={D})"]
                              + trip}
        if not ok:
            out['evidence'].append(
                f"WARNING: supersingular pattern check did not certify over "
                f"{len(trip)} usable primes ≤ {pmax}; verdict stands on the "
                f"j-table, but inspect the model reduction.")
        return out
    # Rational j NOT in the table ⇒ non-CM (the table is complete for j ∈ ℚ).
    # Corroborate: no candidate D may satisfy the full supersingular pattern.
    matches = [D for D in sorted(set(CM_J_TABLE.values()))
               if _supersingular_pattern(aps, D)[1]]
    scan_note = (f"a_p scan over {len(aps)} good primes ≤ {pmax}: no CM "
                 f"discriminant matches the supersingular pattern" if not matches
                 else f"a_p scan matched D={matches} at pmax={pmax} — finite-"
                      f"sample coincidence; the table verdict (non-CM) stands")
    return {'is_cm': False, 'discriminant_D': None, 'certified': True,
            'method': "CM-13 table completeness for rational j + a_p scan "
                      "(wired point counter)",
            'evidence': ev + [f"j = {j} NOT in class-number-1 CM table",
                              scan_note, ('ap', sorted(aps.items()))]}


def _toric_banana_count(a, tval, p):
    """
    #{x ∈ (𝔽_p^×)^L : (a₀ + Σᵢ aᵢxᵢ)(a₀ + Σᵢ aᵢ/xᵢ) = t},  L = len(a) − 1.
    Walks the small torus in the first L−1 variables and solves the resulting
    quadratic in the last variable exactly.
    """
    from itertools import product
    inv = [0]*p
    for k in range(1, p):
        inv[k] = pow(k, p-2, p)
    a0 = a[0] % p; rest = [ai % p for ai in a[1:]]
    tt = tval % p
    L = len(a) - 1
    cnt = 0
    rng = range(1, p)
    for xs in product(rng, repeat=L-1):
        S1 = a0; S2 = a0
        for ai, xi in zip(rest[:-1], xs):
            S1 = (S1 + ai*xi) % p
            S2 = (S2 + ai*inv[xi]) % p
        aL = rest[-1]
        # (S1 + aL z)(S2 + aL/z) = t  ⇒  aL S2 z² + (S1 S2 + aL² − t) z + aL S1 = 0
        A = (aL * S2) % p
        B = (S1*S2 + aL*aL - tt) % p
        C = (aL * S1) % p
        if A == 0:
            if B == 0:
                cnt += (p-1) if C == 0 else 0
            else:
                u = (-C * inv[B]) % p
                if u != 0:
                    cnt += 1
        else:
            disc = (B*B - 4*A*C) % p
            if disc == 0:
                u = (-B * inv[(2*A) % p]) % p
                cnt += (1 if u != 0 else 0)
            elif pow(disc, (p-1)//2, p) == 1:
                cnt += 2   # both roots; z=0 impossible since C = aL S1, drop if 0
                if C == 0:
                    cnt -= 1
    return cnt


def _fit_quadratic_exact(pts):
    """Exact rational (c2, c1, c0) with c2 p² + c1 p + c0 through the three
       (p, count) points, or None when degenerate."""
    (p1, r1), (p2, r2), (p3, r3) = pts
    den = (p1 - p2) * (p1 - p3) * (p2 - p3)
    if den == 0:
        return None
    c2 = Fraction(r1 * (p2 - p3) - r2 * (p1 - p3) + r3 * (p1 - p2), den)
    c1 = (Fraction(r1 - r2) - c2 * (p1 * p1 - p2 * p2)) / (p1 - p2)
    c0 = Fraction(r1) - c2 * p1 * p1 - c1 * p1
    return (c2, c1, c0)


def k3_ap_probe(raw, D_candidates=(-3, -4, -7, -8, -11, -15, -20, -24)):
    """
    Exact background centring + Livné χ_D scan on raw K3 point counts
    (the background-subtraction step the docstring of cm_check_K3 promises).

    The algebraic (Néron–Severi) part of the count is an integer polynomial
    c₂p² + c₁p + c₀ for the models this probe accepts; a CM (rank-2)
    transcendental motive makes the residual a_p VANISH on the χ_D-inert
    primes (density ½), which LOCKS the background: some exact quadratic
    through 3 of the counts must reproduce the count at ≥ max(4, n/3) primes,
    with every residual inside the weight-3 Deligne bound |a_p| ≤ 2p.  The
    centring is exact arithmetic throughout — never a least-squares fit.

    Returns {'locked', 'background', 'support', 'n_primes', 'resid',
             'matches', 'zero_density', 'all_zero'} where matches lists every
    candidate D with  a_p = 0 ⟺ χ_D(p) = −1  on ALL unramified primes.
    No lock ⇒ either non-CM (a_p zeros too sparse to centre) or a
    character-twisted background — the caller must stay inconclusive.
    """
    primes = sorted(raw)
    n = len(primes)
    out = {'locked': False, 'background': None, 'support': 0, 'n_primes': n,
           'resid': None, 'matches': [], 'zero_density': None, 'all_zero': False}
    if n < 8:
        return out
    # enumerate exact quadratics through triples; keep integer-coefficient ones
    cands = set()
    for trip in combinations(primes, 3):
        c = _fit_quadratic_exact([(p, raw[p]) for p in trip])
        if c is not None and all(x.denominator == 1 for x in c):
            cands.add(tuple(int(x) for x in c))
    best = None
    for c in sorted(cands):
        c2, c1, c0 = c
        resid = {p: raw[p] - (c2*p*p + c1*p + c0) for p in primes}
        if any(abs(r) > 2*p for p, r in resid.items()):
            continue                      # violates the weight-3 Deligne bound
        support = sum(1 for r in resid.values() if r == 0)
        if best is None or support > best[1]:
            best = (c, support, resid)
    if best is None or best[1] < max(4, (n + 2) // 3):
        if best is not None:
            out['support'] = best[1]
        return out
    c, support, resid = best
    out.update(locked=True, background=c, support=support, resid=resid,
               zero_density=support / n, all_zero=(support == n))
    for D in D_candidates:
        trip = [(p, resid[p], _kronecker(D, p)) for p in primes
                if _kronecker(D, p) != 0]
        if len(trip) >= 8 and all((a == 0) == (chi == -1) for _, a, chi in trip):
            out['matches'].append(D)
    return out


def cm_check_K3(spec, mc, pmax=60, D_candidates=(-3, -4, -7, -8, -11, -15, -20, -24)):
    """
    Livné CM test for a maximal-cut K3:

      Count #X_t(𝔽_p) on the toric banana model for the primes 5 ≤ p ≤ pmax,
      CENTRE the counts by their exact integer-polynomial background
      (k3_ap_probe — locked on the a_p = 0 primes, Deligne-bounded, exact
      arithmetic), and test
            a_p = 0  ⟺  χ_D(p) = −1      over the candidate D < 0.
      Zero density ≈ ½ with a single D matching ⇒ CM (rank-2) PROBE verdict:
      is_cm=True with certified=False — a finite prime scan on a proxy model,
      to be gated before compute moves.  No exact background lock ⇒ honest
      is_cm=None: either non-CM (a_p zeros too sparse to centre) or a
      character-twisted background; the probe refuses rather than guesses.

    The point-count model is the banana toric model
        (Σ aᵢ xᵢ)(Σ aᵢ/xᵢ) = t   over (𝔽_p^×)^{L}.
    For non-banana K3s this model is a HEURISTIC PROXY and is flagged as such.

    CERTIFIED special case: the equal-mass 4-edge banana goes through the
    eta-quotient a_p of the level-15 weight-3 form instead (Livné) — that path
    returns certified=True.
    """
    caveats = []
    primes = [p for p in _primes_upto(pmax) if p >= 5]
    # Masses² from edges; t from fiber_point or 0 (Broadhurst on-shell point).
    edges = spec.get("edges", [])
    a = []
    for e in edges:
        m2 = e[2] if len(e) >= 3 else 0
        try:
            a.append(int(sp.nsimplify(m2)))
        except Exception:
            a.append(1)
    if len(set(a)) > 1 and 0 in a:
        caveats.append("mixed massive/massless edges in K3 toric model; "
                       "background subtraction may be off.")
    if not a:
        a = [1, 1, 1, 1]
    fp = spec.get("fiber_point", {})
    tval = int(fp.get("t", fp.get("s", 0)))

    # Certified path for the EQUAL-MASS 3-loop banana (the case that matters):
    # build f3 = 15.3.b.a a_p from the eta product and run the χ_{−15} test.
    if sorted(a) == [1]*len(a) and len(a) == 4:
        ap = _f3_level15_ap(pmax)
        triples = [(p, ap[p], _kronecker(-15, p)) for p in primes if p not in (3, 5)]
        zeros = [p for p, a_, _ in triples if a_ == 0]
        match = all((a_ == 0) == (chi == -1) for _, a_, chi in triples)
        return {'is_cm': bool(match), 'discriminant_D': -15 if match else None,
                'evidence': triples, 'certified': bool(match),
                'zero_density': len(zeros)/max(1, len(triples)),
                'method': "eta-quotient a_p vs χ_{-15} (Livné, certified)",
                'caveats': caveats}

    # Generic K3: count on the toric model, centre the exact polynomial
    # background, and run the Livné χ_D scan (k3_ap_probe).
    raw = {p: _toric_banana_count(a, tval, p) for p in primes}
    probe = k3_ap_probe(raw, D_candidates)
    base = {'method': "toric point-count + exact background centring + "
                      "Livné χ_D scan (probe)",
            'certified': False,
            'background': list(probe['background']) if probe['background'] else None,
            'support': f"{probe['support']}/{probe['n_primes']}",
            'raw_counts': sorted(raw.items())}
    if mc.get('method') != 'symanzik-F' or mc.get('n_verts') != 2:
        caveats.append("non-banana K3 probed on the banana toric model — "
                       "HEURISTIC PROXY; treat any verdict as a hypothesis.")
    if not probe['locked']:
        caveats.append(
            f"no exact integer-quadratic background locks the counts (best "
            f"support {probe['support']}/{probe['n_primes']}): either non-CM "
            f"(rank ≥ 3 — a_p zeros too sparse to centre) or a character-"
            f"twisted background; probe inconclusive at pmax={pmax}.")
        return {'is_cm': None, 'discriminant_D': None,
                'evidence': sorted(raw.items()), **base, 'caveats': caveats}
    resid = probe['resid']
    if probe['all_zero']:
        caveats.append(
            "centred residuals vanish at EVERY sampled prime — the count is "
            "purely polynomial in p (degenerate/singular fiber, or a fully "
            "algebraic count model); transcendental a_p invisible, probe "
            "inconclusive.")
        return {'is_cm': None, 'discriminant_D': None,
                'evidence': sorted(resid.items()),
                'zero_density': probe['zero_density'], **base, 'caveats': caveats}
    zeros = probe['support']
    nonzeros = probe['n_primes'] - zeros
    if len(probe['matches']) == 1 and zeros >= 3 and nonzeros >= 3:
        D = probe['matches'][0]
        triples = [(p, resid[p], _kronecker(D, p)) for p in sorted(resid)
                   if _kronecker(D, p) != 0]
        caveats.append(
            f"K3 CM verdict is a point-count PROBE (background "
            f"{base['background']} exact on {base['support']} primes, "
            f"Deligne-bounded residuals) — NOT theorem-certified; "
            f"cross-certify the motive rank before committing the route.")
        return {'is_cm': True, 'discriminant_D': D, 'evidence': triples,
                'zero_density': probe['zero_density'], **base, 'caveats': caveats}
    if len(probe['matches']) >= 2:
        caveats.append(
            f"background locks but MULTIPLE candidate discriminants match "
            f"({probe['matches']}) at pmax={pmax} — undecidable at this depth; "
            f"raise pmax.")
    else:
        caveats.append(
            f"background locks (zero density {probe['zero_density']:.2f}) but "
            f"the zero pattern matches no candidate discriminant "
            f"{tuple(D_candidates)} — D outside the candidate list, or a "
            f"non-CM lock; probe inconclusive.")
    return {'is_cm': None, 'discriminant_D': None,
            'evidence': sorted(resid.items()),
            'zero_density': probe['zero_density'], **base, 'caveats': caveats}


def _f3_level15_ap(N):
    """a_p of f3 = (η₃η₅)³ + (η₁η₁₅)³ ∈ S₃(Γ₀(15), χ), by direct q-expansion
       of the eta products."""
    def prodpow(mult, power, M):
        base = [0]*(M+1); base[0] = 1; n = 1
        while mult*n <= M:
            m = mult*n
            for k in range(M, m-1, -1):
                base[k] -= base[k-m]
            n += 1
        res = [0]*(M+1); res[0] = 1
        for _ in range(power):
            new = [0]*(M+1)
            for i in range(M+1):
                if res[i] == 0:
                    continue
                ri = res[i]
                for j in range(M+1-i):
                    new[i+j] += ri*base[j]
            res = new
        return res
    def mul(a, b, M):
        c = [0]*(M+1)
        for i in range(M+1):
            if a[i] == 0:
                continue
            ai = a[i]
            for j in range(M+1-i):
                c[i+j] += ai*b[j]
        return c
    M = N + 4
    A = mul(prodpow(3, 3, M), prodpow(5, 3, M), M)
    B = mul(prodpow(1, 3, M), prodpow(15, 3, M), M)
    a = [0]*(N+1)
    for k in range(M+1):
        if k+1 <= N: a[k+1] += A[k]
        if k+2 <= N: a[k+2] += B[k]
    return a


# ===========================================================================
#  5.  closure_route(classification)  →  {route, reason, ring_conductors, honesty}
# ===========================================================================

def closure_route(cls: dict):
    """
    Map the geometry classification to a closure route ∈ {A,B,C,D,E}.
    Encodes the decision tree at the top of this file and the value-fittability
    criterion (memory/value-fittability-criterion.md).
    """
    tp = cls.get('type')
    honesty = {
        'genus_method': cls.get('genus_method', 'unknown'),
        'pf_order_certified': cls.get('pf_order_certified', False),
        'cm_evidence': cls.get('cm', {}).get('evidence', []),
        'caveats': list(cls.get('caveats', [])),
    }

    # ---- POLYLOG --------------------------------------------------------
    if tp == "polylog":
        return {'route': 'A', 'reason': "genus 0 ⇒ dlog leading singularity ⇒ "
                "pure polylog; value-fit in MPL/MZV ring closes.",
                'ring_conductors': ['MZV'], 'honesty': honesty}

    # ---- ELLIPTIC -------------------------------------------------------
    if tp == "elliptic":
        ell = cls.get('elliptic', {})
        cm = cls.get('cm', {})
        N = ell.get('level_N')
        on_X1N = ell.get('on_X1N')
        grp = ell.get('congruence_group')
        if cm.get('is_cm'):
            D = cm.get('discriminant_D')
            return {'route': 'A',
                    'reason': f"genus-1 CM curve (D={D}, j∈CM-13) ⇒ periods in "
                              f"Chowla–Selberg Γ-ring; value-fit closes.",
                    'ring_conductors': [f"Q(sqrt({D}))", f"Gamma_1({N})" if N else "CM"],
                    'honesty': honesty}
        if on_X1N and N:
            # ring conductor: the Dirichlet L-value generator at this level.
            # Γ₁(6)→χ_{-3}; Γ₀(4)→χ_{-4} (Catalan); Γ₁(5)→χ_5; default → level.
            chi_map = {6: -3, 3: -3, 4: -4, 8: -4, 2: -4, 5: 5, 7: -7, 12: -3}
            D_ring = chi_map.get(N, -N)
            return {'route': 'B',
                    'reason': f"genus-1 congruence-modular family on {grp} "
                              f"(level {N}, X₁({N}) genus-0); Eichler integral on "
                              f"Γ₁({N}) bootstrappable (eMPL/iterated-Eisenstein).",
                    'ring_conductors': [f"Gamma_1({N})", f"L(chi_{{{D_ring}}},k)",
                                        f"D={D_ring}"],
                    'honesty': honesty}
        # non-congruence / generic j
        honesty['caveats'].append(
            "non-congruence verdict from fiber-config / deg(j) only; a "
            "Hauptmodul search on higher-N X₁(N) was NOT run.")
        return {'route': 'C',
                'reason': "genus-1, j NOT on any genus-0 X₁(N) ⇒ non-congruence-"
                          "modular; "
                          "inhomogeneous kernel is OUTSIDE the modular/eMPL span "
                          "⇒ DE-transport mandatory (value-fittability class C).",
                'ring_conductors': [], 'honesty': honesty}

    # ---- K3 -------------------------------------------------------------
    if tp == "K3":
        cm = cls.get('cm', {})
        if cm.get('is_cm') is not None and cm.get('certified') is False:
            honesty['caveats'].append(
                "K3 CM verdict from the point-count probe (exact background "
                "centring + Livné χ_D scan) — NOT theorem-certified; gate the "
                "route before committing compute.")
        if cm.get('is_cm'):
            D = cm.get('discriminant_D')
            tag = "" if cm.get('certified', True) else \
                  " [point-count PROBE verdict — gate before compute]"
            return {'route': 'B',
                    'reason': f"K3, transcendental rank 2 ⇒ CM by Livné "
                              f"(D={D}); boundary in Chowla–Selberg Γ-ring "
                              f"(row-19 K3 banana → A on the boundary, B on the bulk)."
                              f"{tag}",
                    'ring_conductors': [f"Q(sqrt({D}))", "L(f3,s)"],
                    'honesty': honesty}
        if cm.get('is_cm') is False:
            return {'route': 'D',
                    'reason': "K3 with non-CM transcendental motive (rank ≥ 3) "
                              "⇒ Hilbert/Siegel frontier (NONCM_FEYNMAN_REALIZATION); "
                              "no named ring; CY-transport.",
                    'ring_conductors': [], 'honesty': honesty}
        honesty['caveats'].append("K3 CM status undetermined (background "
                                  "centring probe inconclusive).")
        return {'route': 'C',
                'reason': "K3, CM status undetermined; default to DE-transport.",
                'ring_conductors': [], 'honesty': honesty}

    # ---- CY3+ -----------------------------------------------------------
    if tp and tp.startswith("CY"):
        return {'route': 'D',
                'reason': f"{tp} (PF order ≥ 4); Calabi–Yau transport module "
                          f"(banana ladder / MUM-cusp); no value-fit.",
                'ring_conductors': [], 'honesty': honesty}

    # ---- HYPERELLIPTIC --------------------------------------------------
    if tp == "hyperelliptic":
        g = cls.get('genus')
        return {'route': 'E',
                'reason': f"genus {g} hyperelliptic; Siegel modular / "
                          f"Abel–Jacobi on Jac (hard pipeline limit, cf. H+jet g=2).",
                'ring_conductors': [], 'honesty': honesty}

    honesty['caveats'].append(f"unclassified variety type {tp!r}")
    return {'route': '?', 'reason': "unclassified", 'ring_conductors': [],
            'honesty': honesty}


# ===========================================================================
#  6.  classify(graph_spec)  →  full report
# ===========================================================================

def classify(spec_or_path, verbose=False):
    spec = _load_spec(spec_or_path)
    name = spec.get("name", "graph")
    report = {'name': name, 'spec_source': spec.get("description", "")}

    # 1. maximal-cut polynomial
    mc = baikov_maxcut(spec)
    report['maxcut'] = {'method': mc['method'],
                        'poly': str(sp.factor(mc['P'])),
                        'vars': [str(v) for v in mc['vars']],
                        'modulus': str(mc['modulus']) if mc['modulus'] else None}

    # 2. variety classification
    cls = classify_variety(mc)
    report['variety'] = {k: v for k, v in cls.items() if not k.startswith('_')}

    # 3. type-specific invariants
    if cls['type'] == "elliptic":
        Pred = cls.get('_reduced_P', mc['P'])
        vred = cls.get('_reduced_var', mc['vars'][0])
        ell = elliptic_invariants(Pred, vred, mc.get('modulus'))
        cls['elliptic'] = ell
        report['elliptic'] = ell
        # CM check (specialize at fiber_point if given)
        cm = cm_check_elliptic(ell['j'], mc.get('modulus'),
                               spec.get('fiber_point'))
        cls['cm'] = cm
        report['cm'] = cm

    elif cls['type'] == "K3":
        cm = cm_check_K3(spec, mc)
        cls['cm'] = cm
        report['cm'] = cm
        # Sym²-of-elliptic probe: for the equal-mass banana this is exact;
        # recorded as structural evidence.
        if sorted([e[2] if len(e) >= 3 else 0 for e in spec.get("edges", [])]) \
           and len(spec.get("edges", [])) == 4:
            report['sym2_root'] = {'group': 'Gamma_1(6)',
                                   'note': "3-loop banana K3 PF op = Sym² of "
                                           "sunrise Γ₁(6) op (Verrill/Joyce)"}

    # 4. route
    route = closure_route(cls)
    report['route'] = route['route']
    report['reason'] = route['reason']
    report['ring_conductors'] = route['ring_conductors']
    report['honesty'] = route['honesty']
    # merge caveats
    report['honesty']['caveats'] = sorted(set(
        route['honesty'].get('caveats', []) + mc.get('caveats', [])
        + cls.get('caveats', [])))

    if verbose:
        print(json.dumps(report, indent=2, default=str))
    return report


# ===========================================================================
#  CLI
# ===========================================================================

def main():
    ap = argparse.ArgumentParser(description=__doc__.split("THE DECISION")[0])
    ap.add_argument("spec", help="graph_spec JSON path")
    ap.add_argument("--out", help="write full report JSON here", default=None)
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()
    rep = classify(args.spec, verbose=args.verbose)
    if args.out:
        with open(args.out, "w") as f:
            json.dump(rep, f, indent=2, default=str)
        print(f"[geotriage] wrote {args.out}")
    else:
        print(json.dumps(rep, indent=2, default=str))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
r"""
linear_landau.py — Landau / symbol-alphabet engine for MIXED quadratic + LINEAR
(worldline / eikonal / HQET) propagator families.

Extends landau_alphabet.py to the post-Minkowskian / WQFT integral families,
where some propagators are LINEAR in the loop momenta:

    D_quad = ℓ_j(k)²  - m_j²                  (graviton / standard)
    D_lin  = 2 u_a · ℓ_j(k)                   (worldline / eikonal)

with ℓ_j(k) = Σ_l c_{jl} k_l + r_j a fixed loop-momentum routing (r_j a linear
combination of external momenta, typically 0 or ±q for PM).

THEORY (see linear_landau.md for the derivation):

  Schwinger / Symanzik:  Σ_j α_j D_j  =  k^T M(α) k  +  2 Q(α)·k  +  J(α)
    M_{ll'}(α) = Σ_{j∈quad} α_j c_{jl} c_{jl'}                 (linear-prop α's ABSENT)
    Q_l(α)^μ   = Σ_{j∈quad} α_j c_{jl} r_j^μ  +  Σ_{j∈lin} α_j c_{jl} u_{a(j)}^μ
    J(α)       = Σ_{j∈quad} α_j (r_j² - m_j²)

  After Gaussian k-integration:

      U = det M(α)                                            (degree L in α_quad ONLY)
      F = U·J  -  U· Q^T M^{-1} Q                             (degree L+1; ≤ QUADRATIC in each α_lin)

  ⇒ Linear-prop Feynman parameters live in F at degree ≤ 2, NOT in U.  This is
    the parametric counterpart of the eikonal Landau loop equation:

      Σ_{j∈quad on loop l} α_j ℓ_j  +  Σ_{j∈lin on loop l} α_j u_{a(j)}  =  0      (∂/∂k_l)
      α_j D_j = 0                                                                  (on-shell)

    A linear prop contributes a CONSTANT vector α_j u_a to the Landau loop sum,
    exactly like a heavy on-shell external line (Coleman–Norton: an eikonal prop
    is a static worldline of an infinitely heavy particle along u_a).

  PM kinematic ring:  u₁²=u₂²=1, u₁·u₂=y, u_i·q=0, q²=q2 (typically -1).
  Single physical scale y → x via y=(1+x²)/(2x).

API
---
  symanzik_mixed(quad, lin, loops, ext, repl)       → (U, F, αsyms)
  landau_locus_mixed(quad, lin, loops, ext, repl)   → {face: [letters]}, full alphabet
  pm_alphabet_x(letters_y)                          → letters in x via y=(1+x²)/(2x)

INPUT FORMAT
------------
  quad = [(coeffs:dict[loop->int], r:"ext-mom expr", m2), ...]
  lin  = [(coeffs:dict[loop->int], r:"ext-mom expr", u_label), ...]
  loops = ["k1","k2",...]
  ext   = ["q"]                                       (external momenta names)
  repl  = {"u1.u1":1, "u2.u2":1, "u1.u2":"y",
           "u1.q":0, "u2.q":0, "q.q":"q2"}            (scalar-product table)

This module is self-contained sympy; PLD.jl is NOT a dependency.  See
landau_alphabet.py for the pure-quadratic graph-driven engine that this
parallels.

For graphs with linear/eikonal propagators (PM-gravity bootstrap class).
"""
from __future__ import annotations
import sympy as sp
from itertools import combinations, product
import json, sys, argparse, signal, time


class _FaceTimeout(Exception):
    pass


def _alarm(signum, frame):  # pragma: no cover
    raise _FaceTimeout()


# ===========================================================================
#  Scalar-product expander.  We never carry explicit Lorentz indices; instead
#  every "vector" is a formal sympy linear combination over the basis
#  {ext momenta} ∪ {u_a}, and dot(A,B) distributes and looks up the
#  replacement table.
# ===========================================================================
class DotRing:
    def __init__(self, basis_names, repl: dict):
        self.basis = {n: sp.Symbol(f"_v_{n}") for n in basis_names}
        # build a symmetric lookup of basis·basis → expr
        self.tab = {}
        for k, v in repl.items():
            a, b = k.split(".")
            self.tab[frozenset((a, b))] = sp.sympify(v)
        # any unspecified self-product / cross-product defaults to a fresh symbol
        for a in basis_names:
            for b in basis_names:
                key = frozenset((a, b))
                if key not in self.tab:
                    nm = f"{a}D{b}" if a == b else f"{min(a,b)}D{max(a,b)}"
                    self.tab[key] = sp.Symbol(nm)

    def vec(self, expr_str):
        """Parse 'q' or 'k1-q' or '0' into a dict basis_name -> sympy coeff."""
        if expr_str in (0, "0", None, ""):
            return {}
        e = sp.sympify(expr_str, locals={n: s for n, s in self.basis.items()})
        e = sp.expand(e)
        out = {}
        for n, s in self.basis.items():
            c = e.coeff(s, 1)
            if c != 0:
                out[n] = c
        # sanity: residual must be 0
        return out

    def dot(self, A: dict, B: dict):
        s = sp.Integer(0)
        for a, ca in A.items():
            for b, cb in B.items():
                s += ca * cb * self.tab[frozenset((a, b))]
        return sp.expand(s)


# ===========================================================================
#  symanzik_mixed — build (U, F) for an arbitrary mixed family by Gaussian
#  completion (the matrix-tree route, not the spanning-tree enumeration; the
#  latter has no obvious graph for linear props).
# ===========================================================================
def symanzik_mixed(quad, lin, loops, ext, repl, *, alpha_prefix="a"):
    """
    quad : list of (coeffs:dict loop->int, r:str ext-routing, m2)
    lin  : list of (coeffs:dict loop->int, r:str ext-routing, u_label:str)
    loops: list of loop-momentum names ['k1','k2',...]
    ext  : list of external-vector names ['q', ...]  (u-labels added automatically)
    repl : dict 'a.b' -> expr  (scalar-product replacement table)

    Returns (U, F, alphas, info) with alphas a list of sympy Symbols in the
    SAME order as quad++lin.
    """
    L = len(loops)
    nq, nl = len(quad), len(lin)
    alphas = [sp.Symbol(f"{alpha_prefix}{i+1}", positive=True) for i in range(nq + nl)]
    a_q, a_l = alphas[:nq], alphas[nq:]

    u_labels = sorted({u for *_, u in lin})
    dr = DotRing(list(ext) + u_labels, repl)

    # routing vectors r_j for both prop types, as basis-dicts
    r_q = [dr.vec(r) for (_, r, _) in quad]
    r_l = [dr.vec(r) for (_, r, _) in lin]
    u_l = [dr.vec(u) for (*_, u) in lin]

    # --- M matrix (L×L) -----------------------------------------------------
    M = sp.zeros(L, L)
    for j, (c, _, _) in enumerate(quad):
        for il, kl in enumerate(loops):
            for ilp, klp in enumerate(loops):
                M[il, ilp] += a_q[j] * c.get(kl, 0) * c.get(klp, 0)
    M = sp.Matrix(L, L, lambda i, j: sp.expand(M[i, j]))
    U = sp.factor(sp.expand(M.det()))

    # --- Q_l (vector-valued, length-L list of basis-dicts with α-coeffs) ----
    Q = [dict() for _ in range(L)]
    def _acc(D, vec, coeff):
        for n, c in vec.items():
            D[n] = sp.expand(D.get(n, 0) + coeff * c)
    for j, (c, _, _) in enumerate(quad):
        for il, kl in enumerate(loops):
            cj = c.get(kl, 0)
            if cj:
                _acc(Q[il], r_q[j], a_q[j] * cj)
    for j, (c, _, _) in enumerate(lin):
        for il, kl in enumerate(loops):
            cj = c.get(kl, 0)
            if cj:
                _acc(Q[il], u_l[j], a_l[j] * cj)
                if r_l[j]:
                    _acc(Q[il], r_l[j], 0)   # ext routing on a linear prop: 2u·(k+r) → 2u·k + const; the 2u·r piece goes to J, NOT Q
    # NOTE: a linear prop D = 2u·(Σc k + r) contributes 2α u·(Σc k) to the
    # k-linear part (→ Q) and 2α u·r to the constant part (→ J).  In PM with
    # u·q=0 the latter is zero anyway.

    # --- J scalar -----------------------------------------------------------
    J = sp.Integer(0)
    for j, (_, _, m2) in enumerate(quad):
        J += a_q[j] * (dr.dot(r_q[j], r_q[j]) - sp.sympify(m2))
    for j, (_, _, u) in enumerate(lin):
        if r_l[j]:
            J += 2 * a_l[j] * dr.dot(u_l[j], r_l[j])

    # --- F = U·J - U · Q^T M^{-1} Q  ---------------------------------------
    # adj(M) = U · M^{-1} is polynomial in α; compute it once.
    adj = M.adjugate()
    QTMQ = sp.Integer(0)
    for il in range(L):
        for ilp in range(L):
            QTMQ += adj[il, ilp] * dr.dot(Q[il], Q[ilp])
    F = sp.expand(U * J - QTMQ)
    F = sp.factor(sp.expand(F))

    return U, F, alphas, {"M": M, "adjM": adj, "Q": Q, "J": J, "dr": dr}


# ===========================================================================
#  Landau locus by face: enumerate faces (subsets of propagators with α_j→0
#  for the complement), solve  ∂F_face/∂α_i = 0  ∧  F_face = 0  by resultant
#  elimination, return the kinematic-only factors.
# ===========================================================================
def _kin_factors(expr, alphas, kinsyms):
    """Irreducible factors depending on kinematics only (no α's)."""
    facs = set()
    aset = set(alphas)
    for f in sp.Mul.make_args(sp.factor(expr)):
        b = f.as_base_exp()[0]
        fs = b.free_symbols
        if (fs & kinsyms) and not (fs & aset):
            facs.add(sp.factor(b))
    return facs


def _strip_alpha_content(g, aset):
    """Drop pure-alpha monomial factors and numeric coefficients from a
    generator.  Shared pure-alpha monomial content across generators makes
    subsequent pairwise resultants identically 0 (chain self-annihilation),
    silently losing every letter harvested only at the end.  Returns 0 if
    nothing but alpha-monomial/numeric content remains."""
    out = sp.Integer(1)
    for f in sp.Mul.make_args(sp.factor(g)):
        b, e = f.as_base_exp()
        if b.is_number or (b in aset):
            continue
        out *= f
    return out if out != 1 else sp.Integer(0)


def _discriminant_letters(P, fvars, kinsyms, alphas, face_timeout=30):
    """Bounded resultant elimination of {P, ∂P/∂α_i} over α's.  Returns
    (letters:set, status).  Mirrors landau_alphabet.discriminant_letters
    with the same SIGALRM hard cap (the 2-loop H leading face already
    overflows a naive resultant chain).

    LETTER-LOSS FIX (regression reproducer shipped in tests): on deep
    faces the surviving resultants after one elimination round can share
    pure-alpha monomial content, so subsequent pairwise resultants are
    identically 0 -> cur=[] -> letters harvested only at the end were DROPPED
    while the face still reported 'resolved'.  Reproducer face:
    P = -a2*a5*(a7^2 + 2*y*a7*a8 + a8^2)  ->  pre-fix: resolved, NO letters,
    though (y-1)(y+1) is in every generator after eliminating a7.  Two
    conservative fixes: (i) alpha-monomial content stripped from every
    generator each round; (ii) kinematic factors harvested from EVERY
    intermediate generator set, not only the final one.  Over-generation is
    policed downstream by verify_letter (saturated Groebner)."""
    aset = set(alphas)
    fvars = [v for v in fvars if v in P.free_symbols]
    if not fvars:
        return _kin_factors(P, alphas, kinsyms), "no-free-params"
    chart = fvars[-1]
    Pc = sp.expand(P.subs(chart, 1))
    elim = fvars[:-1]
    gens = [Pc] + [sp.diff(Pc, v) for v in elim]
    cur = [g for g in gens if g != 0]
    remaining = list(elim)
    letters = set()
    signal.signal(signal.SIGALRM, _alarm)
    signal.alarm(face_timeout)
    try:
      for g in cur:
          letters |= _kin_factors(g, alphas, kinsyms)
      while remaining:
        v = remaining.pop()
        withv = [g for g in cur if v in g.free_symbols]
        if not withv:
            continue
        withv.sort(key=lambda g: sp.Poly(g, v).degree())
        piv = withv[0]
        new = []
        for g in cur:
            if g is piv:
                continue
            if v in g.free_symbols:
                r = sp.resultant(piv, g, v)
                if r != 0:
                    new.append(sp.expand(r))
            else:
                new.append(g)
        stripped = set()
        for g in new:
            if g == 0:
                continue
            gs = _strip_alpha_content(g, aset)
            if gs != 0:
                stripped.add(gs)
        cur = list(stripped)
        for g in cur:
            letters |= _kin_factors(g, alphas, kinsyms)
        if not cur:
            break
      status = "resolved"
      for g in cur:
        if g.free_symbols & aset:
            status = "bounded-partial"
      signal.alarm(0)
      return letters, status
    except _FaceTimeout:
      signal.alarm(0)
      letters |= _kin_factors(sp.factor(Pc), alphas, kinsyms)
      return letters, "bounded-timeout"
    finally:
      signal.alarm(0)


def landau_locus_mixed(quad, lin, loops, ext, repl, *,
                       kin_invariants=("y", "q2"),
                       min_face=1, include_second_type=True,
                       face_timeout=30, verbose=False):
    """
    Full face-by-face Landau alphabet for a mixed family.

    Returns dict:
       { 'U': U, 'F': F, 'alphabet': [...], 'faces': [...], 'kin': [...] }
    """
    U, F, alphas, info = symanzik_mixed(quad, lin, loops, ext, repl)
    kinsyms = {sp.Symbol(s) for s in kin_invariants}
    nP = len(alphas)
    faces = []
    raw = set()
    for k in range(nP, min_face - 1, -1):
        for keep in combinations(range(nP), k):
            sub = {alphas[j]: 0 for j in range(nP) if j not in keep}
            Ff = sp.expand(F.subs(sub))
            Uf = sp.expand(U.subs(sub))
            if Ff == 0:
                continue
            fvars = [alphas[j] for j in keep]
            lf, stf = _discriminant_letters(Ff, fvars, kinsyms, alphas, face_timeout)
            lu, stu = (set(), "skip")
            if include_second_type and Uf not in (0, 1) and (set(fvars) & Uf.free_symbols):
                lu, stu = _discriminant_letters(Uf, fvars, kinsyms, alphas, face_timeout)
            raw |= (lf | lu)
            rec = {"keep": keep, "n": k, "letters_F": sorted(map(str, lf)),
                   "letters_U": sorted(map(str, lu)),
                   "status_F": stf, "status_U": stu}
            faces.append(rec)
            if verbose:
                sys.stderr.write(f"[face {keep}] F:{stf} {sorted(map(str,lf))}\n")
    alphabet = _canonical(raw)
    return {"U": U, "F": F, "alphas": alphas, "alphabet": alphabet,
            "faces": faces, "kin": list(kin_invariants)}


def _canonical(raw):
    seen, out = [], []
    for f in sorted(raw, key=lambda z: (sp.count_ops(z), str(z))):
        ff = sp.factor(sp.expand(f))
        if ff.is_number:
            continue
        c, _ = ff.as_coeff_Mul()
        if c.is_number and c not in (1, 0):
            ff = sp.factor(sp.expand(ff / c))
        if any(sp.cancel(ff / g).is_number for g in seen):
            continue
        seen.append(ff); out.append(ff)
    return out


# ===========================================================================
#  y → x re-expression:  y = (1+x²)/(2x).
# ===========================================================================
def pm_alphabet_x(letters_y, *, q2_value=-1):
    x = sp.Symbol("x")
    y = sp.Symbol("y"); q2 = sp.Symbol("q2")
    sub = {y: (1 + x**2) / (2 * x), q2: q2_value}
    out = set()
    for L in letters_y:
        Lx = sp.factor(sp.together(sp.sympify(L).subs(sub)))
        # split numerator and denominator into irreducibles over Q[x]
        num, den = sp.fraction(Lx)
        for poly in (num, den):
            for f in sp.Mul.make_args(sp.factor(poly)):
                b = f.as_base_exp()[0]
                if b.free_symbols:
                    out.add(sp.factor(b))
    return _canonical(out)


# ===========================================================================
#  PRESET PM FAMILIES (used by the test driver and as documentation).
#  Kinematics: u1²=u2²=1, u1·u2=y, u1·q=u2·q=0, q²=q2.
# ===========================================================================
PM_REPL = {"u1.u1": 1, "u2.u2": 1, "u1.u2": "y",
           "u1.q": 0, "u2.q": 0, "q.q": "q2"}


def family_2PM_tri():
    """1-loop {k², (k-q)², 2u1·k} — single-worldline triangle."""
    quad = [({"k": 1}, "0", 0),
            ({"k": 1}, "-q", 0)]
    lin  = [({"k": 1}, "0", "u1")]
    return quad, lin, ["k"], ["q"], PM_REPL


def family_2PM_box():
    """1-loop {k², (k-q)², 2u1·k, 2u2·k} — two-worldline 'box'."""
    quad = [({"k": 1}, "0", 0),
            ({"k": 1}, "-q", 0)]
    lin  = [({"k": 1}, "0", "u1"),
            ({"k": 1}, "0", "u2")]
    return quad, lin, ["k"], ["q"], PM_REPL


def family_3PM_H():
    """2-loop H/double-box: {k1²,k2²,(k1-q)²,(k2-q)²,(k1-k2)², 2u1·k1, 2u2·k2}."""
    quad = [({"k1": 1}, "0", 0),
            ({"k2": 1}, "0", 0),
            ({"k1": 1}, "-q", 0),
            ({"k2": 1}, "-q", 0),
            ({"k1": 1, "k2": -1}, "0", 0)]
    lin  = [({"k1": 1}, "0", "u1"),
            ({"k2": 1}, "0", "u2")]
    return quad, lin, ["k1", "k2"], ["q"], PM_REPL


def family_3PM_IY():
    """2-loop 'IY'/nested triangle (3PM 0SF): {k1²,(k1-q)²,(k1-k2)²,k2², 2u1·k1, 2u1·k2, 2u2·k2}.
    Three worldline props (2 on line u1, 1 on u2) — the kind that produces
    cross-terms u1·u2 in F at degree 2."""
    quad = [({"k1": 1}, "0", 0),
            ({"k1": 1}, "-q", 0),
            ({"k1": 1, "k2": -1}, "0", 0),
            ({"k2": 1}, "0", 0)]
    lin  = [({"k1": 1}, "0", "u1"),
            ({"k2": 1}, "0", "u1"),
            ({"k2": 1}, "0", "u2")]
    return quad, lin, ["k1", "k2"], ["q"], PM_REPL


def family_4PM_K3():
    """3-loop 4PM K3.0 sector (Parra-Martinez et al.).  Propagators per
    2401.07899 fig. K3.0 / 2211.16357: 3 graviton rungs + 2 horizontal +
    3 worldline props on u1, 1 on u2 (1SF)."""
    # k1,k2,k3 routed so the rungs carry k_i, horizontals (k_i - k_{i+1}),
    # last rung (k3 - q); worldline u1 carries k1,k2,k3; u2 the cut.
    quad = [({"k1": 1}, "0", 0),
            ({"k1": 1, "k2": -1}, "0", 0),
            ({"k2": 1, "k3": -1}, "0", 0),
            ({"k3": 1}, "-q", 0),
            ({"k2": 1}, "-q", 0)]
    lin  = [({"k1": 1}, "0", "u1"),
            ({"k2": 1}, "0", "u1"),
            ({"k3": 1}, "0", "u1"),
            ({"k2": 1}, "0", "u2")]
    return quad, lin, ["k1", "k2", "k3"], ["q"], PM_REPL


# ===========================================================================
#  γ=3 (x = 3-2√2) singularity — discriminant of the topology-40 K3′
#  leading-singularity surface  Q₆(t1,t2;x)=0  from 2505.10274 eq.(3.26):
#      Q₆ = (x(1+t1)(t1+t2²) - (1+x²) t1 t2)² - 4 x² t1² t2²
# ===========================================================================
def apery_K3_landau_locus():
    """Return the discriminant in x of the Apéry-K3 LS surface Q₆(t₁,t₂;x)=0.

    Q₆ is the maximal-cut leading-singularity surface of topology #40 from
    2505.10274 eq.(3.26).  Its t-discriminant is the projection of the
    leading Landau variety onto kinematic x — i.e. the locus where the K3
    period integral pinches.
    """
    x, t1, t2 = sp.symbols("x t1 t2")
    Q6 = sp.expand((x*(1+t1)*(t1+t2**2) - (1+x**2)*t1*t2)**2 - 4*x**2*t1**2*t2**2)
    # Step 1: Res_{t1}(Q6, ∂_{t1}Q6) → product of factors in (t2, x)
    r1 = sp.factor(sp.resultant(Q6, sp.diff(Q6, t1), t1))
    facs = set()
    for f in sp.Mul.make_args(r1):
        b = f.as_base_exp()[0]
        fs = b.free_symbols
        if fs == {x}:
            facs.add(sp.factor(b))
        elif t2 in fs and x in fs:
            # Step 2: disc_{t2} of each (t2,x)-factor
            d = sp.factor(sp.discriminant(b, t2))
            for g in sp.Mul.make_args(d):
                gb = g.as_base_exp()[0]
                if gb.free_symbols == {x}:
                    facs.add(sp.factor(gb))
    return _canonical(facs), Q6


# ===========================================================================
#  verify_letter — re-check a candidate letter against the SATURATED Groebner
#  elimination ideal of {F, ∂F/∂α_i} on the torus (∏α_i ≠ 0).  The bounded
#  resultant chain over-generates (e.g. spurious y on the 3PM_H (0,3,4,5,6)
#  face); this is the honesty filter.
# ===========================================================================
def verify_letter(F, alphas_face, letter, kinsyms, *, timeout=30):
    """Returns True iff `letter` divides the saturated elimination ideal of
    the leading Landau system on this face (i.e. is a GENUINE Landau
    component, not a resultant artifact)."""
    fvars = [a for a in alphas_face if a in F.free_symbols]
    if not fvars:
        return False
    Fc = sp.expand(F.subs(fvars[-1], 1))
    elim = fvars[:-1]
    t = sp.Symbol("_tRab")
    eqs = [Fc] + [sp.diff(Fc, v) for v in elim]
    sat = 1 - t * sp.prod(elim) if elim else sp.Integer(0)
    signal.signal(signal.SIGALRM, _alarm)
    signal.alarm(timeout)
    try:
        gb = sp.groebner(eqs + ([sat] if sat != 0 else []),
                         t, *elim, *kinsyms, order="lex")
        signal.alarm(0)
    except _FaceTimeout:
        signal.alarm(0)
        return None  # unknown
    except Exception:
        signal.alarm(0)
        return None
    # the elimination ideal in kinsyms is generated by the kin-only polys in gb
    for g in gb:
        if g.free_symbols <= set(kinsyms) and g != 0:
            if sp.gcd(sp.Poly(g, *kinsyms), sp.Poly(letter, *kinsyms)).total_degree() > 0:
                return True
    return False


# ===========================================================================
#  CLI / quick test driver
# ===========================================================================
def _run_preset(name):
    fams = {"2PM_tri": family_2PM_tri, "2PM_box": family_2PM_box,
            "3PM_H": family_3PM_H, "3PM_IY": family_3PM_IY,
            "4PM_K3": family_4PM_K3}
    fam = fams[name]()
    res = landau_locus_mixed(*fam, kin_invariants=("y", "q2"))
    return res


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("preset", nargs="?", default="2PM_box",
                    choices=["2PM_tri","2PM_box","3PM_H","3PM_IY","4PM_K3","apery"])
    ap.add_argument("--json-out", default=None)
    args = ap.parse_args(argv)
    if args.preset == "apery":
        facs, Q6 = apery_K3_landau_locus()
        print("Apery-K3 LS-surface discriminant factors in x:")
        for f in facs:
            print("  ", f)
        return 0
    res = _run_preset(args.preset)
    print(f"[{args.preset}]  U = {res['U']}")
    print(f"[{args.preset}]  F = {sp.factor(res['F'])}")
    print(f"[{args.preset}]  alphabet (in y, q2): {res['alphabet']}")
    ax = pm_alphabet_x(res["alphabet"])
    print(f"[{args.preset}]  alphabet (in x, q2=-1): {ax}")
    if args.json_out:
        with open(args.json_out, "w") as fh:
            json.dump({"preset": args.preset,
                       "U": str(res["U"]), "F": str(res["F"]),
                       "alphabet_y": [str(a) for a in res["alphabet"]],
                       "alphabet_x": [str(a) for a in ax],
                       "n_faces": len(res["faces"])}, fh, indent=2)
    return 0


if __name__ == "__main__":
    sys.exit(main())

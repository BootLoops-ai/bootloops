#!/usr/bin/env python3
r"""
canonical_form.py — canonical-ε-form / UT-basis builder.

Sits BETWEEN the IBP cache (Kira master basis + DE matrix) and the value-fit
(gatekeeper.heldout_certify).  Implements the principle: the fit must sample
the LS-defined canonical UT basis, NOT the raw Kira basis.  In practice that
means

    J_i = ε^{p_i} · (1/LS_i(x)) · I_i^raw      (diagonal leading-order rotation)

where LS_i is the maximal-cut leading singularity of the i-th master, computed
loop-by-loop in the Baikov representation.  This module produces that diagonal
rotation directly from the family definition (the same `family` block
amflow-cpp already consumes), and exports it in the format

    {ut_tag: {raw_tag: 'p/q' | mpf}}      or      (matrix, row_names, col_names)

that `gatekeeper.heldout_cv.heldout_certify(ut_rotation=...)` accepts.

The full Henn algorithm (fuchsify → normalize → ε-factorize → dlog) is
DELEGATED to the in-house Julia engine `Counterweight` (Lee 1411.0911 over
Nemo/Flint), which this module wraps via a subprocess bridge.  No Mathematica
(CANONICA / Libra) is required, but if those are present the bridge can be
swapped.  Where no engine is available the heavy steps raise a documented
NotImplemented with a recipe pointer.

PUBLIC API
==========
    Family(spec)                                  — parse the amflow `family` block
    leading_singularities(family, masters)        → {tag: LSResult}
    ut_rotation(family, masters, point=None)      → heldout_certify-format rotation
    fuchsify(A, x, eps)                           → wrap_epsfactor (Julia)
    factor_epsilon(A, x, eps)                     → wrap_epsfactor (Julia)
    to_dlog(Atilde, x, alphabet)                  → {letter: constant ℚ-matrix}
    elliptic_period(quartic, var, point)          → ϖ₀(point) numerically

ENGINES
=======
    ε-factorization engine            : ../ — the Counterweight engine in this package (Fuchsia.jl)
    heldout_certify ut_rotation slot  : tools/gatekeeper/heldout_cv.py:_apply_ut_rotation
"""
from __future__ import annotations
import json
import os
import re
import subprocess
import sys
import tempfile
from dataclasses import dataclass, field, asdict
from fractions import Fraction
from itertools import combinations

import mpmath as mp
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
COUNTERWEIGHT_ROOT = os.environ.get(
    "COUNTERWEIGHT_ROOT", os.path.normpath(os.path.join(HERE, "..")))
COUNTERWEIGHT_BRIDGE = os.path.join(HERE, "counterweight_bridge.jl")

eps = sp.Symbol('eps')


# ===========================================================================
#  Family — parse the amflow-cpp `family` block into a Baikov rep.
# ===========================================================================
class Family:
    """Parse a `family` block (the same JSON sub-object amflow-cpp consumes)
    into the data the Baikov LS engine needs:

        loops, ext (independent external momenta), scalar-product symbols,
        the kinematic replacement table, and the propagator polynomials in
        scalar products.

    `spec` keys:
        loops        : ["k1","k2",...]       loop momenta names
        legs         : ["p1",...,"pN"]       external leg names
        conservation : {"pN":"-p1-...-p_{N-1}"}   eliminated leg(s)
        replacement  : {"p1^2":"0", "(p1+p2)^2":"s", ...}  kinematics
        propagators  : ["k1^2", "(k1-p1)^2 - m2", ...]
    """

    _MOM_TOK = re.compile(r'[A-Za-z_]\w*')

    def __init__(self, spec: dict):
        self.spec  = spec
        self.loops = list(spec["loops"])
        self.legs  = list(spec.get("legs", []))
        # eliminated legs via momentum conservation
        cons = spec.get("conservation", {}) or {}
        self.elim = {k: self._parse_vec(v) for k, v in cons.items()}
        self.ext  = [p for p in self.legs if p not in self.elim]
        # independent external momenta actually entering the Gram (drop the last
        # one if all p_i^2 are fixed AND conservation already removed one — we
        # keep ALL independent ext; the Gram is over loops + ext[:-?])
        # Standard Baikov uses E = #ext - 1 independent external momenta after
        # conservation; `ext` here already has conservation applied.
        self.E = len(self.ext)
        self.L = len(self.loops)

        # ---- scalar-product symbols ------------------------------------
        # loop-loop and loop-ext are the Baikov SPs (the integration variables
        # before changing to z_i = D_i); ext-ext are fixed by `replacement`.
        self._sp_sym: dict[frozenset, sp.Symbol] = {}
        self._sp_loop: list[sp.Symbol] = []
        for i, a in enumerate(self.loops):
            for b in self.loops[i:]:
                s = sp.Symbol(f'sp_{a}_{b}')
                self._sp_sym[frozenset((a, b))] = s
                self._sp_loop.append(s)
        for a in self.loops:
            for b in self.ext:
                s = sp.Symbol(f'sp_{a}_{b}')
                self._sp_sym[frozenset((a, b))] = s
                self._sp_loop.append(s)
        self.nSP = len(self._sp_loop)      # = L(L+1)/2 + L·E

        # ---- kinematic replacement table  (ext·ext) --------------------
        # Build a complete  p_i·p_j  table from the `replacement` rules by
        # expanding each "(Σ p)^2 = X" linearly.
        self.kin_syms = sorted(
            {sp.Symbol(str(s)) for v in (spec.get("replacement") or {}).values()
             for s in sp.sympify(v).free_symbols},
            key=str)
        self._ext_dot = self._build_ext_dot_table(spec.get("replacement") or {})

        # ---- propagators as polys in scalar products -------------------
        self.props_raw = list(spec["propagators"])
        self.props = [self._expand_propagator(p) for p in self.props_raw]
        self.nP = len(self.props)

        # ---- Baikov variables z_i and the SP→z linear map --------------
        self.z = sp.symbols(f'z1:{self.nP + 1}')
        # solve  D_i(SP) = z_i  for the loop SPs.  This is LINEAR in SP.
        eqs = [sp.Eq(self.props[i], self.z[i]) for i in range(self.nP)]
        sol = sp.solve(eqs, self._sp_loop, dict=True)
        if not sol:
            # under-determined (fewer props than SPs) — solve for the first
            # len(props) SPs that actually appear; remaining SPs stay free
            # (they are the ISP / residual Baikov directions).
            present = [s for s in self._sp_loop
                       if any(p.has(s) for p in self.props)]
            sol = sp.solve(eqs, present[:self.nP], dict=True)
        if not sol:
            raise ValueError("Family: could not invert propagator → Baikov "
                             "linear map (propagators not independent?)")
        self._sp_to_z = sol[0]
        self._free_sp = [s for s in self._sp_loop if s not in self._sp_to_z]

        # ---- Baikov polynomial F = det Gram(loops, ext) in z -----------
        self._gram = self._build_gram()
        self._F_sp = sp.expand(self._gram.det())
        self._Fz   = sp.expand(self._F_sp.subs(self._sp_to_z))

    # -------------------------------------------------------------------
    # vector / dot-product helpers (a "vector" is a dict {name: coeff})
    # -------------------------------------------------------------------
    def _parse_vec(self, expr_str: str) -> dict[str, sp.Expr]:
        """'-p1 - p2 - p3'  →  {'p1': -1, 'p2': -1, 'p3': -1}."""
        e = sp.sympify(expr_str, evaluate=False)
        # collect over the leg/loop names appearing
        names = set(self._MOM_TOK.findall(expr_str)) & set(self.legs + self.loops)
        out: dict[str, sp.Expr] = {}
        # represent each name as a fresh non-commutative? No — momenta are
        # vectors so the *expression* is purely a ℚ-linear combo of names.
        syms = {n: sp.Symbol(n) for n in (self.legs + self.loops)}
        e2 = sp.sympify(expr_str, locals=syms)
        for n, s in syms.items():
            c = sp.expand(e2).coeff(s)
            if c != 0:
                out[n] = c
        return out

    def _resolve(self, name: str) -> dict[str, sp.Expr]:
        """Expand a momentum name through `conservation`."""
        if name in self.elim:
            return dict(self.elim[name])
        return {name: sp.Integer(1)}

    def _dot(self, A: dict, B: dict) -> sp.Expr:
        """Minkowski dot of two ℚ-linear momentum combinations."""
        out = sp.Integer(0)
        for a, ca in A.items():
            ra = self._resolve(a)
            for b, cb in B.items():
                rb = self._resolve(b)
                for u, cu in ra.items():
                    for v, cv in rb.items():
                        out += ca * cb * cu * cv * self._sp(u, v)
        return sp.expand(out)

    def _sp(self, u: str, v: str) -> sp.Expr:
        """Scalar product of two BASIC (already-resolved) momentum names."""
        key = frozenset((u, v))
        if key in self._sp_sym:
            return self._sp_sym[key]
        # both external — look up the kinematic table
        return self._ext_dot.get(key, sp.Integer(0))

    def _build_ext_dot_table(self, repl: dict) -> dict[frozenset, sp.Expr]:
        """Turn  {"p1^2":0, "(p1+p2)^2":"s", ...}  into a complete  p_i·p_j
        table over the INDEPENDENT external momenta `self.ext`.  The system
        is linear in the unknown dot products."""
        unk = {}
        for i, a in enumerate(self.ext):
            for b in self.ext[i:]:
                unk[frozenset((a, b))] = sp.Symbol(f'_dd_{a}_{b}')
        eqs = []
        syms = {n: sp.Symbol(n) for n in self.legs}
        for lhs, rhs in repl.items():
            # parse "(p1+p2)^2" or "p1*p3" etc. into a polynomial in the legs,
            # then expand into dot-products (using `unk` for ext·ext) and set
            # equal to rhs.
            # We support "(...)^2" and "Pi^2" forms (the only forms the family
            # block uses).
            m = re.fullmatch(r'\s*\(([^)]*)\)\s*\^\s*2\s*', lhs)
            if m:
                vec = self._parse_vec(m.group(1))
            else:
                m2 = re.fullmatch(r'\s*([A-Za-z_]\w*)\s*\^\s*2\s*', lhs)
                if m2:
                    vec = {m2.group(1): sp.Integer(1)}
                else:
                    # generic "p1*p2" etc.
                    e = sp.sympify(lhs, locals=syms)
                    # only support a*b monomial
                    fac = sp.Mul.make_args(e)
                    if len(fac) == 2 and all(str(f) in self.legs for f in fac):
                        a, b = str(fac[0]), str(fac[1])
                        eqs.append(sp.Eq(self._dot_ext({a: 1}, {b: 1}, unk),
                                         sp.sympify(rhs)))
                        continue
                    raise ValueError(f"Family.replacement: cannot parse '{lhs}'")
            eqs.append(sp.Eq(self._dot_ext(vec, vec, unk), sp.sympify(rhs)))
        sol = sp.solve(eqs, list(unk.values()), dict=True)
        tab = {k: sp.Integer(0) for k in unk}
        if sol:
            for k, sym in unk.items():
                tab[k] = sp.simplify(sol[0].get(sym, sp.Integer(0)))
        return tab

    def _dot_ext(self, A, B, unk):
        out = sp.Integer(0)
        for a, ca in A.items():
            ra = self._resolve(a)
            for b, cb in B.items():
                rb = self._resolve(b)
                for u, cu in ra.items():
                    for v, cv in rb.items():
                        out += ca * cb * cu * cv * unk[frozenset((u, v))]
        return sp.expand(out)

    # -------------------------------------------------------------------
    def _expand_propagator(self, pstr: str) -> sp.Expr:
        """'(k1-p1)^2 - m2'  →  polynomial in loop/ext scalar products and
        kinematic symbols.  The mass term (anything not a momentum square)
        is sympified directly."""
        # split additive top level: "(...)² [- m2]" — the mass piece is
        # everything that does NOT contain a loop/ext symbol after the square.
        # Simplest robust route: replace each momentum name by a sympy
        # MatrixSymbol-like placeholder?  We instead exploit that propagators
        # are always  (Σ ± q_a)^2  -  m²  in the family blocks we consume.
        m = re.match(r'\s*\(([^)]*)\)\s*\^\s*2\s*(.*)$', pstr)
        if m:
            vec = self._parse_vec(m.group(1))
            rest = m.group(2).strip()
        else:
            m2 = re.match(r'\s*([A-Za-z_]\w*)\s*\^\s*2\s*(.*)$', pstr)
            if not m2:
                raise ValueError(f"Family.propagator: cannot parse '{pstr}'")
            vec = {m2.group(1): sp.Integer(1)}
            rest = m2.group(2).strip()
        sq = self._dot(vec, vec)
        if rest:
            sq = sq + sp.sympify(rest)
        return sp.expand(sq)

    def _build_gram(self):
        mom = self.loops + self.ext
        n = len(mom)
        return sp.Matrix(n, n, lambda i, j:
                         self._dot({mom[i]: 1}, {mom[j]: 1}))

    # -------------------------------------------------------------------
    def baikov_poly(self) -> sp.Expr:
        """Full Baikov polynomial F(z₁..z_nP, free_SP, kin)."""
        return self._Fz

    def baikov_exponent(self) -> sp.Expr:
        """The Baikov measure exponent  γ = (d - L - E - 1)/2  with d=4-2ε."""
        return sp.Rational(4 - self.L - self.E - 1, 2) - eps

    def sector_of(self, indices) -> list[int]:
        """1-based propagator indices with positive power."""
        return [i + 1 for i, a in enumerate(indices) if a > 0]


# ===========================================================================
#  Leading singularities  (loop-by-loop Baikov, iterated residues)
# ===========================================================================
@dataclass
class LSResult:
    tag: str
    sector: list[int]
    kind: str                      # 'rational' | 'sqrt' | 'elliptic' | 'unknown'
    expr: sp.Expr | None = None    # the LS itself, e.g. 1/(s*t) or 1/sqrt(R)
    radicand: sp.Expr | None = None
    quartic: sp.Expr | None = None
    residual_vars: list = field(default_factory=list)
    note: str = ""

    def to_json(self):
        d = {"tag": self.tag, "sector": self.sector, "kind": self.kind,
             "expr": str(self.expr) if self.expr is not None else None,
             "radicand": str(self.radicand) if self.radicand is not None else None,
             "quartic": str(self.quartic) if self.quartic is not None else None,
             "residual_vars": [str(v) for v in self.residual_vars],
             "note": self.note}
        return d


def _maxcut_ls(Fcut: sp.Expr, residual: list[sp.Symbol],
               gamma0: sp.Rational | None = None) -> LSResult:
    """Iterated-residue maximal-cut leading singularity.

    Input is the Baikov polynomial restricted to the cut, F_cut(v_k; kin),
    and the list of residual integration variables v_k.  At d=4 the measure
    is  ∫ ∏ dv_k · F_cut^{γ₀}  with γ₀ = (4-L-E-1)/2 ∈ ½ℤ.  The LS is the
    multivariate residue / period of that integrand.

    What this routine returns is the KINEMATIC CONTENT of the LS — the
    rational/algebraic prefactor that the UT rotation divides out.  Overall
    numerical constants and signs are dropped; they are fixed downstream by
    the ε⁻ᴸ unit-normalization test (see verify_unit_pole()).

    Localization rules (one variable at a time, greedy on lowest degree):

        deg 0       → free ISP direction; drop.
        deg 1       → F = a·v + b ; residue at v=-b/a gives Jacobian 1/a and
                      the next-stage polynomial is `a` (the only surviving
                      kinematic carrier).  Track under sqrt as a².
        deg 2       → F = a v² + b v + c ; ∫ dv/√F ∝ 1/√a, and the dlog
                      residue (= leading singularity) carries 1/√(b²-4ac).
                      Track √a, continue with F → disc = b²-4ac.
                      (This sqrt-class extraction has been verified against
                      numeric purity checks.)
        deg 3,4 (last var) → ELLIPTIC; return the quartic (LS = ϖ₀).
        otherwise   → UNKNOWN; route decision deferred to geotriage.

    For 0 residual variables (full localization, e.g. 1-loop box): LS is just
    F_cut^{γ₀}, returned as 'rational' if γ₀∈ℤ or with the appropriate sqrt.
    """
    F = sp.expand(Fcut)
    res = [v for v in residual if sp.expand(F).has(v)]
    LS_rad = sp.Integer(1)       # accumulates under the sqrt
    chain = []
    while res:
        # pick lowest-degree variable
        degs = {w: sp.Poly(F, w).degree() for w in res}
        v = min(res, key=lambda w: (degs[w], str(w)))
        d = degs[v]
        res.remove(v)
        P = sp.Poly(F, v)
        chain.append((str(v), d))
        if d == 0:
            continue
        if d == 1:
            a = P.nth(1)
            LS_rad *= a ** 2       # 1/a = 1/sqrt(a²) → keep everything under sqrt
            F = sp.expand(a)
            res = [w for w in res if F.has(w)]
            continue
        if d == 2:
            a, b, c = P.nth(2), P.nth(1), P.nth(0)
            disc = sp.expand(b * b - 4 * a * c)
            LS_rad *= a
            F = sp.factor(disc)
            res = [w for w in res if F.has(w)]
            continue
        if d in (3, 4) and not res:
            # elliptic: leading singularity is the holomorphic period of
            # y² = F(v).  Return the quartic and let elliptic_period() handle
            # numerics.
            return LSResult(tag="", sector=[], kind="elliptic",
                            expr=None,
                            radicand=sp.factor(LS_rad),
                            quartic=sp.factor(F),
                            residual_vars=[v],
                            note=f"elliptic curve y^2 = {sp.factor(F)} in {v}")
        # higher / multivariate residual
        return LSResult(tag="", sector=[], kind="unknown",
                        expr=None, radicand=sp.factor(LS_rad),
                        quartic=None,
                        residual_vars=[v] + res,
                        note=f"residual deg {d} in {v} with {len(res)} vars "
                             f"left; CY/hyperelliptic — use geotriage")
    # fully localized → polylog.  Any leftover kinematic-constant factor in F
    # (the last discriminant / the 0-residual F_cut itself) is the final
    # radicand, weighted by the Baikov exponent.
    rad = sp.factor(LS_rad)
    Ffin = sp.factor(F)
    if gamma0 is not None and not chain:
        # 0-residual case: LS = F_cut^{γ₀} on the nose.
        e2 = sp.Rational(-2) * gamma0       # so LS = 1/sqrt(F^{e2})
        rad = rad * Ffin ** int(e2) if e2.is_integer else rad * Ffin
    else:
        rad = rad * Ffin
    sq, rem = _split_square(rad)
    if rem == 1 or rem == -1:
        return LSResult(tag="", sector=[], kind="rational",
                        expr=sp.together(1 / sq), radicand=None,
                        residual_vars=[],
                        note=f"polylog (genus 0); LS ∝ 1/({sq}); chain={chain}")
    return LSResult(tag="", sector=[], kind="sqrt",
                    expr=1 / (sq * sp.sqrt(rem)),
                    radicand=rem, residual_vars=[],
                    note=f"polylog (genus 0, algebraic letter); "
                         f"LS ∝ 1/({sq}·sqrt({rem})); chain={chain}")


def _split_square(expr: sp.Expr):
    """Factor `expr` and split off the perfect-square part:
    expr = sq² · rem  with rem square-free."""
    f = sp.factor(expr)
    sq = sp.Integer(1)
    rem = sp.Integer(1)
    if f.is_Mul:
        args = f.args
    else:
        args = (f,)
    for a in args:
        b, e = (a.base, a.exp) if a.is_Pow else (a, sp.Integer(1))
        if e.is_Integer:
            sq  *= b ** (int(e) // 2)
            rem *= b ** (int(e) % 2)
        else:
            rem *= a
    # constant sign/magnitude → push into sq² where possible
    if rem.is_number:
        r = sp.nsimplify(rem)
        rt = sp.sqrt(r)
        if rt.is_rational:
            sq *= rt; rem = sp.Integer(1)
    return sp.factor(sq), sp.factor(rem)


def leading_singularities(family, masters, *, maxcut_override=None):
    """Compute the maximal-cut leading singularity of each master.

    Parameters
    ----------
    family  : Family  OR  the dict spec (will be wrapped in Family())
    masters : list of {"tag": str, "indices": [int,...]}  (amflow-cpp format)
    maxcut_override : optional dict {tag: {"poly": str, "vars": [str,...]}}
              Per-master shortcut in the SAME schema as
              tools/geotriage/graph_specs/*.json["maxcut"].
              Use when the automatic Baikov is degenerate (e.g. equal-mass
              sunrise, where γ₀=0 in the standard rep and one must use the
              loop-by-loop / Symanzik rep instead) or when the curve is
              already known from the literature.

    Returns
    -------
    dict {tag: LSResult}
    """
    fam = family if isinstance(family, Family) else Family(family)
    Fz = fam.baikov_poly()
    z = fam.z
    g0 = fam.baikov_exponent().subs(eps, 0)
    mco = maxcut_override or {}
    out = {}
    for m in masters:
        idx = m["indices"] if isinstance(m, dict) else list(m)
        tag = m.get("tag") if isinstance(m, dict) else None
        tag = tag or "I[" + ",".join(map(str, idx)) + "]"
        sec = fam.sector_of(idx)
        # per-master override (geotriage "maxcut" schema)
        mo = mco.get(tag) or (m.get("maxcut") if isinstance(m, dict) else None)
        if mo:
            vrs = [sp.Symbol(v) for v in mo.get("vars", [])]
            r = _maxcut_ls(sp.sympify(mo["poly"]), vrs, gamma0=None)
            r.note += f" [override: {mo.get('source','user')}]"
        else:
            sub = {z[i - 1]: 0 for i in sec}
            Fcut = sp.expand(Fz.subs(sub))
            residual = [z[i] for i in range(fam.nP) if (i + 1) not in sec] \
                       + list(fam._free_sp)
            r = _maxcut_ls(Fcut, residual, gamma0=g0)
        r.tag, r.sector = tag, sec
        out[tag] = r
    return out


# ===========================================================================
#  Elliptic period  ϖ₀  (numerical, via complete elliptic integral)
# ===========================================================================
def _mp_num(x, dps: int) -> mp.mpc:
    """Convert a sympy number to mp.mpc at `dps` digits — never through a
    Python float/complex (which would silently cap accuracy at 53 bits)."""
    xe = sp.N(sp.sympify(x), dps)
    re_, im_ = xe.as_real_imag()
    return mp.mpc(mp.mpf(sp.Float(re_, dps)._mpf_),
                  mp.mpf(sp.Float(im_, dps)._mpf_))


def elliptic_period(quartic: sp.Expr, var: sp.Symbol, point: dict,
                    dps: int = 60) -> mp.mpc:
    """Holomorphic period ϖ₀ of the curve y² = P(z) (deg 3 or 4) at the
    given kinematic `point`.

    Standard reduction to Legendre form (Byrd & Friedman 1954), generic in
    the polynomial.  Quartic with roots e₁..e₄ (sorted by real part):
    ϖ₀ = 2 K(k²) / √(lc·(e₄-e₂)(e₃-e₁)) with k² the cross ratio.  Cubic with
    roots e₁..e₃ (B&F 235.00 — an exact reduction, no root-at-infinity
    trick): ϖ₀ = 4 K(k²) / √(lc·(e₃-e₁)) with k² = (e₂-e₁)/(e₃-e₁), which is
    2 ∫_{e₁}^{e₂} dz/√P(z) for real sorted roots.

    Root-finding and all arithmetic run at the requested `dps` (mp.polyroots
    on full-precision coefficients — no float64 hop anywhere).
    """
    mp.mp.dps = dps
    P = sp.Poly(quartic.subs({sp.Symbol(k): v for k, v in point.items()}), var)
    deg = P.degree()
    if deg not in (3, 4):
        raise ValueError(f"elliptic_period: need cubic/quartic, got deg {deg}")
    coeffs = [_mp_num(c, dps) for c in P.all_coeffs()]
    lc = coeffs[0]
    roots = mp.polyroots([c / lc for c in coeffs], maxsteps=200, extraprec=dps)
    # sort by real part for a canonical pairing
    roots = sorted(roots, key=lambda r: (mp.re(r), mp.im(r)))
    if deg == 3:
        e1, e2, e3 = roots
        k2 = (e2 - e1) / (e3 - e1)
        return 4 * mp.ellipk(k2) / mp.sqrt(lc * (e3 - e1))
    e1, e2, e3, e4 = roots
    k2 = ((e3 - e2) * (e4 - e1)) / ((e4 - e2) * (e3 - e1))
    pre = 2 / mp.sqrt((e4 - e2) * (e3 - e1))
    return pre * mp.ellipk(k2) / mp.sqrt(lc)


# ===========================================================================
#  ut_rotation — the deliverable that plugs into heldout_certify
# ===========================================================================
def ut_rotation(family, masters, point=None, *, eps_power=None, dps=60,
                include_eps_shift=False):
    """Build the diagonal LS rotation T₀ = T(x, ε=0) that maps raw Kira
    masters → UT basis:

        J_tag  =  (1 / LS_tag(point)) · I_tag

    Returned in the format `gatekeeper.heldout_cv.heldout_certify`
    accepts:

        {ut_tag: {raw_tag: coeff}}

    If `point` is given, coeffs are mpmath numbers (so the dict can be passed
    straight to heldout_certify).  If `point` is None, coeffs are sympy
    strings in the kinematic variables (for serialization / per-point
    evaluation by the caller).

    `include_eps_shift`: if True, also return a sibling dict
    {tag: p_tag} of recommended ε-power shifts (so the UT master is
    ε^{p}·(1/LS)·I_raw).  These are NOT applied to the rotation matrix
    itself (heldout_certify rotates per ε-order; the ε-shift is selected by
    the `weight` argument there).
    """
    fam = family if isinstance(family, Family) else Family(family)
    LS = leading_singularities(fam, masters)
    rot: dict[str, dict[str, object]] = {}
    eps_shift: dict[str, int] = {}
    for tag, r in LS.items():
        ut = f"ut_{tag}"
        if r.kind in ("rational", "sqrt"):
            inv_ls = sp.simplify(1 / r.expr)       # this is LS⁻¹⁻¹? no:
            # r.expr IS the LS (e.g. 1/(s t)); the UT normalization multiplies
            # the raw master by 1/LS? — convention: canonical g = ε^p · LS · I_raw
            # (the LS is the Jacobian, so multiplying BY it makes the maximal
            # cut = 1).  But "LS" in the literature is usually the residue
            # ITSELF (1/(st) for the box), and the UT master is  s t · I_box.
            # Our r.expr = 1/(s t) = LS, hence the rotation coeff = 1/LS = s t.
            coeff = sp.simplify(1 / r.expr)
            if point is None:
                rot[ut] = {tag: str(coeff)}
            else:
                val = coeff.subs({sp.Symbol(k): sp.nsimplify(v)
                                  for k, v in point.items()})
                rot[ut] = {tag: _to_coeff(val, dps)}
        elif r.kind == "elliptic":
            if point is None:
                rot[ut] = {tag: f"1/varpi0[{r.quartic}]"}
            else:
                w0 = elliptic_period(r.quartic, r.residual_vars[0], point, dps)
                # account for any rational pre-sqrt accumulated before the
                # elliptic step:
                pre = sp.sqrt(r.radicand or 1).subs(
                    {sp.Symbol(k): sp.nsimplify(v) for k, v in point.items()})
                rot[ut] = {tag: mp.mpc(1) / (_mp_num(pre, dps) * w0)}
        else:
            rot[ut] = {tag: "NotImplemented:" + r.note}
        eps_shift[tag] = (eps_power or {}).get(tag, fam.L)  # default ε^L
    if include_eps_shift:
        return rot, eps_shift
    return rot


def _to_coeff(val: sp.Expr, dps):
    """Render a sympy value into the format `_apply_ut_rotation` accepts:
    'p/q' for rationals, mpf/mpc otherwise."""
    v = sp.nsimplify(val, rational=True)
    if v.is_Rational:
        return f"{v.p}/{v.q}" if v.q != 1 else str(v.p)
    mp.mp.dps = dps
    return _mp_num(val, dps)


def verify_unit_pole(values: dict, ut_rot: dict, *, L: int, tol_digits=8):
    """UT-normalization sanity check on amflow oracle data.

    `values` is the gatekeeper `load_oracle_values` dict
    {raw_tag: {ptkey: {order: mpc}}}.  After applying `ut_rot` (per the
    heldout_cv convention), the leading ε⁻ᴸ coefficient of each UT master
    must be a kinematic CONSTANT — the same rational number at every point.

    Returns {ut_tag: {"const": mpf, "spread_digits": int, "ok": bool}}.
    This is the purity check in operational form: the rotated pole
    coefficient must reproduce the same constant at every point.
    """
    # local import to avoid a hard dep at module load
    sys.path.insert(0, os.path.join(HERE, "..", "gatekeeper"))
    # repo layout: the gatekeeper package sits beside this package
    sys.path.insert(0, os.path.join(HERE, "..", "..", "gatekeeper"))
    from heldout_cv import _apply_ut_rotation     # type: ignore
    rotated = _apply_ut_rotation(values, ut_rot)
    out = {}
    for ut, pts in rotated.items():
        if not ut.startswith("ut_"):
            continue
        lead = []
        for pk, orders in pts.items():
            if -L in orders:
                lead.append(orders[-L])
        if not lead:
            out[ut] = {"const": None, "spread_digits": 0, "ok": False,
                       "note": f"no eps^{-L} data"}
            continue
        ref = lead[0]
        spread = max((mp.fabs(v - ref) for v in lead[1:]), default=mp.mpf(0))
        sd = int(-mp.log10(spread / max(mp.fabs(ref), mp.mpf(1)))) \
             if spread > 0 else 999
        out[ut] = {"const": ref, "spread_digits": sd,
                   "ok": sd >= tol_digits, "n_points": len(lead)}
    return out


def ut_rotation_matrix(family, masters, point):
    """Same as ut_rotation(..., point=point) but returned as
    (matrix, row_names, col_names) — the alternative format heldout_certify
    accepts."""
    rot = ut_rotation(family, masters, point)
    rows = sorted(rot)
    cols = sorted({c for r in rot.values() for c in r})
    M = mp.matrix(len(rows), len(cols))
    for i, r in enumerate(rows):
        for j, c in enumerate(cols):
            v = rot[r].get(c, 0)
            M[i, j] = v if isinstance(v, (mp.mpf, mp.mpc)) else _to_mp(v)
    return M, rows, cols


def _to_mp(x):
    if isinstance(x, str) and "/" in x:
        a, b = x.split("/"); return mp.mpf(a) / mp.mpf(b)
    return mp.mpf(str(x))


# ===========================================================================
#  fuchsify / factor_epsilon / to_dlog
# ===========================================================================
def _have_counterweight():
    return os.path.isdir(os.path.join(COUNTERWEIGHT_ROOT, "src"))


def wrap_epsfactor(A: sp.Matrix, x: sp.Symbol, *, action="epsform",
                   timeout=600) -> dict:
    """Call the in-house Julia engine `Counterweight` (Lee 1411.0911 over Nemo).

    `A` is a sympy Matrix with entries rational in (x, eps).  The bridge
    serializes A as strings, runs `counterweight_bridge.jl`, and returns the
    parsed JSON result:

        {"ok": bool, "fuchsian": bool, "epsform": bool,
         "T": [[str,...],...], "Atilde": [[str,...],...],
         "singular_points": [...], "report": {...}}

    Actions:
        "certify"  — only run certify_epsform (cheap)
        "epsform"  — full try_epsfactor (fuchsify+normalize+decouple)
    """
    if not _have_counterweight():
        raise NotImplementedError(_NO_ENGINE_MSG)
    n = A.shape[0]
    payload = {
        "action": action,
        "n": n, "x": str(x), "eps": "eps",
        "A": [[str(sp.together(A[i, j])) for j in range(n)] for i in range(n)],
    }
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as fh:
        json.dump(payload, fh)
        inp = fh.name
    cmd = ["julia", "--project=" + COUNTERWEIGHT_ROOT, COUNTERWEIGHT_BRIDGE, inp]
    try:
        out = subprocess.run(cmd, capture_output=True, text=True,
                             timeout=timeout, check=False)
    finally:
        try: os.unlink(inp)
        except OSError: pass
    if out.returncode != 0:
        raise RuntimeError(f"Counterweight bridge failed:\n{out.stderr}")
    return json.loads(out.stdout)


def fuchsify(A: sp.Matrix, x: sp.Symbol, eps_: sp.Symbol = eps):
    """Reduce A to Fuchsian form (all poles simple).  Delegated to the Julia
    `Counterweight` engine (Fuchsia.jl).  Returns (A_fuchsian, T, report)."""
    r = wrap_epsfactor(A, x, action="epsform")
    if not r.get("fuchsian"):
        raise NotImplementedError(
            "Counterweight reports non-Fuchsian input.  The Julia engine's "
            "Barkatou-Moser reduction (Fuchsia.moser_reduce via "
            "Normalize.try_epsfactor) runs on such input, but this Python "
            "wrapper keys on the pre-reduction certificate and does not "
            "consume the reduced matrix -- non-Fuchsian input is a STOP here. "
            "Recipe: (i) strip the LS prefactor from the raw masters FIRST "
            "(leading_singularities()); the LS-normalized DE is Fuchsian for "
            "dlog sectors by construction. (ii) For genuinely irregular "
            "(Moser-irreducible) singular points no automatic route exists "
            "(see e.g. github.com/magv/fuchsia.git, `fuchsify`).")
    T = _strmat(r["T"]) if r.get("T") else None
    Anew = _strmat(r["Anew"]) if r.get("Anew") else _strmat(r.get("Atilde"))
    return Anew, T, r.get("report", {})


def factor_epsilon(A: sp.Matrix, x: sp.Symbol, eps_: sp.Symbol = eps):
    """Find T(x,ε) such that  T A T⁻¹ + (∂ₓT) T⁻¹ = ε · Ã(x).

    Delegated to the Julia `Counterweight` engine (Normalize.jl, Lee's algorithm).
    Returns (Atilde, T, report).  Raises NotImplementedError with a recipe
    pointer if no engine is available or the engine STOPs."""
    r = wrap_epsfactor(A, x, action="epsform")
    if not r.get("ok"):
        raise NotImplementedError(
            "Counterweight STOP: " + str(r.get("report", {}).get("stop", "?")) +
            "\n" + _NO_ENGINE_MSG)
    return _strmat(r["Atilde"]), _strmat(r["T"]), r.get("report", {})


def _strmat(M):
    return sp.Matrix([[sp.sympify(e) for e in row] for row in M])


_NO_ENGINE_MSG = (
    "Full ε-factorization requires an external engine.  Available wrapped:\n"
    f"  • Counterweight (Julia, Lee 1411.0911) at {COUNTERWEIGHT_ROOT}\n"
    "Not found / not requested.  Alternatives (NOT auto-wrapped):\n"
    "  • CANONICA (Meyer, 1611.01087)  — Mathematica, Get[\"CANONICA.m\"]\n"
    "  • Libra    (Lee,   2012.00279)  — Mathematica\n"
    "  • Fuchsia  (Gituliar–Magerya, 1701.04269) — Python2/Maple,\n"
    "             github.com/gituliar/fuchsia\n"
    "  • epsilon  (Prausa, 1701.00725) — C++, github.com/mprausa/epsilon\n"
    "Recipe: export A as a Mathematica/Maple matrix and run the tool there;\n"
    "import T back and verify with to_dlog().  The MINIMUM-VIABLE rotation\n"
    "for UT-basis sampling is leading_singularities() → ut_rotation(), which\n"
    "needs NO engine.")


def to_dlog(Atilde: sp.Matrix, x: sp.Symbol, alphabet: list):
    """Express  Ã(x) = Σ_a  A_a · ∂ₓ log a(x)  and return {a: A_a}.

    `alphabet` is a list of sympy expressions (or strings) in x.  We
    partial-fraction Ã in x, collect the residue at each simple pole, and
    match against the pole structure of the dlog forms.  Any residue not
    spanned by the alphabet's poles is reported in the '_unmatched' key.
    """
    n = Atilde.shape[0]
    letters = [sp.sympify(a) for a in alphabet]
    # pole → which letters carry it (with their residue in ∂log a)
    pole_map = {}
    for k, a in enumerate(letters):
        dl = sp.cancel(sp.diff(sp.log(a), x))
        num, den = sp.fraction(sp.together(dl))
        for r, mult in sp.roots(sp.Poly(den, x), x).items():
            res = sp.limit(dl * (x - r), x, r)
            pole_map.setdefault(sp.nsimplify(r), {})[k] = res
    # residues of Ã at every pole appearing in any entry
    poles = set()
    for i in range(n):
        for j in range(n):
            _, d = sp.fraction(sp.together(Atilde[i, j]))
            for r in sp.roots(sp.Poly(d, x), x):
                poles.add(sp.nsimplify(r))
    # solve  Σ_k A_a[k] · res_k(p) = ResÃ(p)  for each pole p, jointly over
    # all letters touching p.
    Aa = {str(a): sp.zeros(n) for a in letters}
    unmatched = {}
    for p in poles:
        R = sp.Matrix(n, n, lambda i, j:
                      sp.limit(sp.together(Atilde[i, j]) * (x - p), x, p))
        if p not in pole_map:
            unmatched[str(p)] = R
            continue
        contrib = pole_map[p]             # {letter_idx: residue_of_dlog}
        if len(contrib) == 1:
            (k, c), = contrib.items()
            Aa[str(letters[k])] += R / c
        else:
            # multiple letters share this pole — under-determined at a single
            # pole; record and let the caller resolve via a second pole.
            unmatched[str(p)] = (R, contrib)
    out = {a: sp.simplify(M) for a, M in Aa.items()}
    if unmatched:
        out["_unmatched"] = unmatched
    # sanity: reconstruct and compare
    recon = sp.zeros(n)
    for a, M in Aa.items():
        recon += M * sp.cancel(sp.diff(sp.log(sp.sympify(a)), x))
    out["_residual"] = sp.simplify(Atilde - recon)
    return out


# ===========================================================================
#  CLI
# ===========================================================================
def _main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("target", help="target.json with a `family` block (or bare family.json)")
    ap.add_argument("--point", help='JSON kinematic point, e.g. \'{"s":-3,"t":-7}\'')
    ap.add_argument("--out", help="write LS + ut_rotation JSON here")
    a = ap.parse_args(argv)
    spec = json.load(open(a.target))
    fam = Family(spec.get("family", spec))
    masters = spec.get("masters") or []
    if not masters:
        print("[warn] no 'masters' in target — nothing to rotate", file=sys.stderr)
    LS = leading_singularities(fam, masters)
    pt = json.loads(a.point) if a.point else None
    rot = ut_rotation(fam, masters, pt)
    rec = {"family": spec.get("name", fam.spec.get("name", "?")),
           "point": pt,
           "leading_singularities": {t: r.to_json() for t, r in LS.items()},
           "ut_rotation": {u: {k: (str(v) if not isinstance(v, str) else v)
                               for k, v in d.items()}
                           for u, d in rot.items()}}
    js = json.dumps(rec, indent=2)
    if a.out:
        open(a.out, "w").write(js)
    print(js)


if __name__ == "__main__":
    _main()

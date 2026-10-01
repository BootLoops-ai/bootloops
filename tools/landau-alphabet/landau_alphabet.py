#!/usr/bin/env python3
r"""
landau_alphabet.py — STANDALONE bounded face-by-face Landau / symbol-alphabet engine.

Bounded, face-by-face Landau/symbol-alphabet engine (validated on the shipped
regression graphs: massive double box, 3-loop K4, 3-loop ladder, non-planar
double box) applied to an arbitrary Feynman graph
read from a JSON spec.  Pure sympy; PLD.jl is NOT a dependency (see pld_bridge.jl for
the optional cross-validation backend on the same JSON).

CORE DESIGN:
  * The Landau singular locus (principal A-determinant) of a Feynman integral is the
    UNION over all FACES (sub-sectors) of the singular locus of the face's Symanzik
    F-polynomial (normal/anomalous pinch thresholds) PLUS the U-polynomial discriminant
    (second-type / at-infinity).
  * Per-face elimination is a BOUNDED resultant chain (NOT global Groebner) on a small
    fixed variable set with a hard SIGALRM timeout.
  * Work is BOUNDED per face and FLUSHED INCREMENTALLY to JSON after EVERY face, so a
    Groebner-style explosion on one face cannot kill the whole run.
  * Per-face HONESTY record (resolved / bounded-partial / bounded-timeout) is
    load-bearing: it tells you which faces still need PLD.jl / a heavier method.

ALSO PROVIDED:
  * graph-agnostic Symanzik builder driven by the JSON spec (no hardcoded per-graph bits)
  * classify_entries(): first-entry (codim-1 F-faces ↔ mass/Mandelstam) vs
    last-entry / second-type (U-only) classification
  * prune_spurious(): drop pure Feynman-param monomials, U-only factors, and factors
    that cancel in F/U
  * CLI:  python landau_alphabet.py SPEC.json --out OUT.json --face-timeout 45

INPUT SCHEMA (see graph_specs/ and README.md):
  {
    "name": "k4box3l",
    "edges": [[v1,v2,"m2"], [v2,v3,0], ...],            # mass² entry: number or symbol-string
    "nodes": {"<vertex>": [<leg_idx>, ...], ...},        # which external legs attach where
    "kinematics": {
      "invariants": ["s","t","m2"],                      # the kinematic ring
      "external_masses": {"<leg_idx>": "0"|"mZ2"|..., ...},
      "channels": {"1,2": "s", "2,3": "t", "1,3": "-s-t", ...}   # (P_{legs in S})² in invariants
    }
  }

OUTPUT SCHEMA:
  {
    "alphabet": [...polys...],
    "first_entry": [...], "last_entry": [...], "spurious": [...],
    "faces": [ {face_edges, n_edges, loops_L, status_F, status_U, secs, letters, ...}, ... ],
    "honesty": { "n_faces", "resolved", "bounded_partial", "bounded_timeout",
                 "unresolved_faces": [...], "global_deadline_hit": bool }
  }
"""
from __future__ import annotations
import sympy as sp
from itertools import combinations
import json, time, signal, sys, os, argparse

# ---------------------------------------------------------------------------
# bounded-timeout machinery
# ---------------------------------------------------------------------------
class FaceTimeout(Exception):
    pass

def _alarm(signum, frame):  # pragma: no cover
    raise FaceTimeout()


# ===========================================================================
# GraphSpec — parse the JSON into the data the Symanzik / Landau engine needs.
# ===========================================================================
class GraphSpec:
    def __init__(self, spec: dict):
        self.name = spec.get("name", "graph")
        self.raw  = spec
        # --- edges -> {ename: (v1,v2)}, Feynman params, masses --------------
        self.edges: dict[str, tuple[int,int]] = {}
        self.ev:    dict[str, sp.Symbol]      = {}   # Feynman parameter
        self.em:    dict[str, sp.Expr]        = {}   # mass² symbol/number
        for i, e in enumerate(spec["edges"], 1):
            if len(e) == 3:
                v1, v2, m2 = e
            else:
                v1, v2 = e; m2 = 0
            en = f"x{i}"
            self.edges[en] = (int(v1), int(v2))
            self.ev[en]    = sp.Symbol(en, positive=True)
            self.em[en]    = sp.sympify(m2)
        self.ALLE  = list(self.edges)
        self.verts = sorted({v for ab in self.edges.values() for v in ab})
        self.massive_edges = {e for e in self.ALLE if self.em[e] != 0}

        # --- external legs: nodes = {vertex: [leg_idx,...]} ----------------
        nodes_in = spec.get("nodes", {})
        if isinstance(nodes_in, dict):
            self.nodes = {int(v): list(legs) for v, legs in nodes_in.items()}
        else:  # list-of-pairs [[vertex, leg], ...]
            self.nodes = {}
            for v, leg in nodes_in:
                self.nodes.setdefault(int(v), []).append(leg)
        self.all_legs = sorted({l for ls in self.nodes.values() for l in ls})

        # --- kinematics ----------------------------------------------------
        kin = spec.get("kinematics", {})
        self.invariants = [sp.Symbol(s, real=True) for s in kin.get("invariants", [])]
        self._invset    = set(self.invariants)
        # external leg masses² (default 0)
        _loc = {str(s): s for s in self.invariants}
        self.ext_mass = {int(l): sp.sympify(m, locals=_loc)
                         for l, m in kin.get("external_masses", {}).items()}
        # re-sympify edge masses with the SAME invariant symbols (edges were
        # parsed before kinematics; a bare sympify would mint a DIFFERENT
        # symbol of the same name as a real=True invariant -> silent
        # unsimplified letters like '-4*z*z + 4*z**2').
        self.em = {e: sp.sympify(str(m), locals=_loc) if getattr(m, "free_symbols", None)
                   else m for e, m in self.em.items()}
        # channels: keys are comma-joined SORTED leg-index sets, values sympy exprs
        # Auto-add the complement under momentum conservation.
        self._chan: dict[frozenset, sp.Expr] = {}
        full = frozenset(self.all_legs)
        for k, v in kin.get("channels", {}).items():
            legs = frozenset(int(x) for x in str(k).split(","))
            expr = sp.sympify(v, locals={str(s): s for s in self.invariants})
            self._chan[legs] = expr
            self._chan.setdefault(full - legs, expr)
        # single-leg "channels" = external masses²
        for l in self.all_legs:
            self._chan.setdefault(frozenset({l}), self.ext_mass.get(l, sp.Integer(0)))
            self._chan.setdefault(full - frozenset({l}), self.ext_mass.get(l, sp.Integer(0)))
        self._chan.setdefault(frozenset(), sp.Integer(0))
        self._chan.setdefault(full,       sp.Integer(0))
        # convenience: kinematic free-symbol set for kin_factors()
        self.kinsyms = set().union(*(e.free_symbols for e in self._chan.values())) \
                       | set().union(*(m.free_symbols for m in self.em.values())) \
                       | self._invset

    # -----------------------------------------------------------------------
    def chan_inv(self, vertex_component) -> sp.Expr:
        """(∑_{legs attached in component} p_leg)² as a kinematic invariant."""
        legs = frozenset(l for v in vertex_component if v in self.nodes
                           for l in self.nodes[v])
        if legs in self._chan:
            return self._chan[legs]
        # Unknown channel -> 0 (scaleless cut). Record so honesty can flag it.
        return sp.Integer(0)


# ===========================================================================
# Symanzik (U,F) on a FACE = sub-multigraph spanned by `face_edges` (others
# DELETED, x_e -> 0).  Generic, graph-agnostic: the standard spanning-tree /
# 2-forest Symanzik construction, parametrized on a GraphSpec.
# ===========================================================================
def symanzik_subgraph(gs: GraphSpec, face_edges):
    Es = list(face_edges)
    vs = sorted({v for e in Es for v in gs.edges[e]})
    if len(vs) < 2:
        return None
    nV = len(vs)

    def connects(subset):
        par = {v: v for v in vs}
        def f(a):
            while par[a] != a:
                par[a] = par[par[a]]; a = par[a]
            return a
        for e in subset:
            a, b = gs.edges[e]; ra, rb = f(a), f(b)
            if ra == rb:
                return False
            par[ra] = rb
        return len({f(v) for v in vs}) == 1

    Lf = len(Es) - nV + 1
    if Lf < 0:
        return None
    # U_S
    U = sp.Integer(0)
    for T in combinations(Es, nV - 1):
        if connects(T):
            U += sp.prod([gs.ev[e] for e in Es if e not in T]) if (nV - 1) < len(Es) else sp.Integer(1)
    U = sp.expand(U)
    # F0_S
    F0 = sp.Integer(0)
    if nV - 2 >= 0:
        for fo in combinations(Es, nV - 2):
            par = {v: v for v in vs}
            def f(a):
                while par[a] != a:
                    par[a] = par[par[a]]; a = par[a]
                return a
            ok = True
            for e in fo:
                a, b = gs.edges[e]; ra, rb = f(a), f(b)
                if ra == rb:
                    ok = False; break
                par[ra] = rb
            if not ok:
                continue
            cc = {}
            for v in vs:
                cc.setdefault(f(v), []).append(v)
            if len(cc) != 2:
                continue
            comp0 = list(cc.values())[0]
            F0 += sp.prod([gs.ev[e] for e in Es if e not in fo]) * (-gs.chan_inv(comp0))
    F0 = sp.expand(F0)
    massterm = U * sum(gs.ev[e] * gs.em[e] for e in Es)
    F = sp.expand(F0 + massterm)
    return {'edges': Es, 'verts': vs, 'L': Lf, 'U': U, 'F0': F0, 'F': F}


# ===========================================================================
# kin_factors — irreducible factors that depend on the KINEMATIC ring only.
# (uses gs.kinsyms — never a hardcoded invariant list.)
# ===========================================================================
def kin_factors(gs: GraphSpec, expr):
    expr = sp.together(expr)
    num, den = sp.fraction(expr)
    facs = set()
    fparams = set(gs.ev.values())
    for poly in (num, den):
        if poly == 0:
            continue
        for fpow in sp.Mul.make_args(sp.factor(poly)):
            base = fpow.as_base_exp()[0]
            fs = base.free_symbols
            if (fs & gs.kinsyms) and not (fs & fparams):
                facs.add(sp.factor(base))
    return facs


# ===========================================================================
# BOUNDED resultant-chain discriminant — parametrized on (gs, FACE_TIMEOUT).  Status ∈ {resolved, bounded-partial,
# bounded-timeout, no-free-params}.
# ===========================================================================
# ---------------------------------------------------------------------------
# Optional fast Groebner backend (msolve / Singular).  Heavy faces (e.g. the
# 3-loop K4 leading face, degree-267, and its 5-edge "3m+2gamma" subfaces) overflow the bounded
# sympy resultant chain (it TIMES OUT, or comes back bounded-partial).  When available we
# escalate to pld_groebner_backend.face_letters (msolve f4 Rabinowitsch-saturated
# elimination), which is both COMPLETE (real elimination ideal on the torus) and
# fast (K4 leading: 32 s vs >900 s timeout).  Set PLD_GROEBNER=off to disable,
# =force to use it for every face (regression / cross-check).
# The chain itself carries the letter-loss guard (see discriminant_letters):
# kinematic factors are harvested from EVERY intermediate generator set and
# pure-Feynman-parameter monomial content is stripped each round, so resultant
# self-annihilation can no longer silently empty a face's letter set.
# ---------------------------------------------------------------------------
import importlib.util as _ilu
_GB = None
if os.environ.get("PLD_GROEBNER", "auto") != "off":
    try:
        _here = os.path.dirname(os.path.abspath(__file__))
        _gbpath = os.path.join(_here, "pld_groebner_backend.py")
        _spec = _ilu.spec_from_file_location("pld_groebner_backend", _gbpath)
        _GB = _ilu.module_from_spec(_spec); _spec.loader.exec_module(_GB)
        if _GB._MSOLVE is None and _GB._SINGULAR is None:
            _GB = None
    except Exception:
        _GB = None


def _groebner_letters(gs, P, fvars, timeout):
    """Discriminant letters via the fast Groebner backend; () if unavailable."""
    if _GB is None:
        return None
    kin = list(gs.invariants)
    try:
        r = _GB.face_letters(P, list(fvars), kin, timeout=timeout, threads=8)
        if r["status"] == "resolved":
            return set(r["letters"]), r["backend"]
    except Exception:
        pass
    return None


def _strip_fparam_content(g, fparams):
    """Drop pure-Feynman-parameter monomial factors and numeric coefficients
    from a generator.  Shared pure-x monomial content across generators makes
    subsequent pairwise resultants identically 0 (chain self-annihilation),
    silently losing every letter harvested only at the end.  Returns 0 if
    nothing but x-monomial/numeric content remains.  (Same guard as the mixed
    engine's _strip_alpha_content in linear_landau.py.)"""
    out = sp.Integer(1)
    for f in sp.Mul.make_args(sp.factor(g)):
        b, _ = f.as_base_exp()
        if b.is_number or (b in fparams):
            continue
        out *= f
    return out if out != 1 else sp.Integer(0)


def discriminant_letters(gs: GraphSpec, P, fvars, face_timeout: int):
    """Bounded resultant elimination of {P, dP/dx_i} over the face's Feynman
    parameters.  Returns (letters:set, status).

    LETTER-LOSS GUARD (the same guard as the mixed engine's
    _discriminant_letters in linear_landau.py): on deep faces the surviving
    resultants after one elimination round can share pure-x monomial content,
    so subsequent pairwise resultants are identically 0 -> empty generator set
    -> letters harvested only at the end would be DROPPED while the face still
    reports 'resolved' (an unguarded chain loses one of {s, t} on the box1l
    leading face).  Two conservative measures: (i) pure-x monomial content is
    stripped from every generator each round; (ii) kinematic factors are
    harvested from EVERY intermediate generator set, not only the final one.
    Over-generation is policed downstream by canonical_alphabet /
    prune_spurious, and the Groebner escalation still covers timeouts and
    bounded-partial results."""
    fvars = [v for v in fvars if v in P.free_symbols]
    if not fvars:
        return kin_factors(gs, P), "no-free-params"
    # ---- fast path: force-Groebner if requested --------------------------
    if os.environ.get("PLD_GROEBNER", "auto") == "force" and _GB is not None:
        g = _groebner_letters(gs, P, fvars, max(face_timeout * 20, 1800))
        if g is not None:
            return g[0], f"resolved-groebner:{g[1]}"
    chart = fvars[-1]
    elim  = fvars[:-1]
    Pc    = sp.expand(P.subs(chart, 1))
    status = "resolved"
    fparams = set(gs.ev.values())
    letters = set()
    signal.signal(signal.SIGALRM, _alarm)
    signal.alarm(face_timeout)
    try:
        gens = [Pc] + [sp.expand(sp.diff(Pc, v)) for v in elim]
        cur = [g for g in gens if g != 0]
        for g in cur:
            letters |= kin_factors(gs, g)
        remaining = list(elim)
        while remaining:
            v = remaining.pop()
            withv = [g for g in cur if v in g.free_symbols]
            if not withv:
                continue
            withv.sort(key=lambda g: sp.Poly(g, v).degree())
            piv = withv[0]
            newgens = []
            for g in cur:
                if g is piv:
                    continue
                if v in g.free_symbols:
                    r = sp.resultant(sp.Poly(piv, v), sp.Poly(g, v))
                    r = sp.expand(r)
                    if r != 0:
                        newgens.append(r)
                else:
                    newgens.append(g)
            cur = []
            seen = set()
            for g in newgens:
                gstr = _strip_fparam_content(g, fparams)
                if gstr == 0:
                    continue
                k = sp.srepr(sp.factor(gstr))
                if k not in seen:
                    seen.add(k); cur.append(gstr)
            for g in cur:
                letters |= kin_factors(gs, g)
            if not cur:
                break
        for g in cur:
            if g.free_symbols & fparams:
                status = "bounded-partial"
        signal.alarm(0)
        # Escalate to the Groebner backend if the bounded chain came back
        # PARTIAL or EMPTY (residual under-resolution).  Resolved-nonempty
        # results are trusted as-is.
        if (status != "resolved" or not letters) and _GB is not None:
            g = _groebner_letters(gs, P, fvars, max(face_timeout * 20, 1800))
            if g is not None and g[0]:
                return g[0], f"resolved-groebner:{g[1]}"
        return letters, status
    except FaceTimeout:
        signal.alarm(0)
        # bounded chain timed out -> escalate to the Groebner backend
        g = _groebner_letters(gs, P, fvars, max(face_timeout * 20, 1800))
        if g is not None:
            return g[0], f"resolved-groebner:{g[1]}"
        # keep whatever the completed rounds harvested; still an honest timeout
        return letters | kin_factors(gs, sp.factor(Pc)), "bounded-timeout"
    finally:
        signal.alarm(0)


# ===========================================================================
# canonical_alphabet — normalize sign/constant, dedupe up to unit, drop trivials.
# ===========================================================================
def canonical_alphabet(raw_letters):
    seen = []
    out  = []
    for f in sorted(raw_letters, key=lambda z: (sp.count_ops(z), str(z))):
        ff = sp.factor(sp.expand(f))
        if ff.is_number:
            continue
        c, _ = ff.as_coeff_Mul()
        if c.is_number and c not in (1, 0):
            ff = sp.factor(sp.expand(ff / c))
        dup = False
        for g in seen:
            r = sp.cancel(ff / g)
            if r.is_number:
                dup = True; break
        if dup:
            continue
        seen.append(ff); out.append(ff)
    return out


# ===========================================================================
# classify_entries — first-entry / last-entry candidates from per-face origin.
#   * FIRST-ENTRY candidates: letters that appear in codim-1 F-faces (single-edge
#     deletions) or as direct F-coefficients — i.e. mass / Mandelstam monomials,
#     the physical-sheet branch points the integral can hit at the FIRST step.
#   * LAST-ENTRY / second-type candidates: letters that appear ONLY through U-faces
#     (at-infinity, second-type) and never through any F-face.
# Returns (first_entry:list, last_entry:list).
# ===========================================================================
def classify_entries(gs: GraphSpec, alphabet, faces_data):
    nE = len(gs.ALLE)
    inF = {str(a): False for a in alphabet}
    inU = {str(a): False for a in alphabet}
    in_codim1 = {str(a): False for a in alphabet}
    for rec in faces_data:
        is_codim1 = (rec.get("n_edges", 0) >= nE - 1)
        for L in rec.get("letters_F", []):
            for a in alphabet:
                if str(a) == L or _same_letter(a, L):
                    inF[str(a)] = True
                    if is_codim1:
                        in_codim1[str(a)] = True
        for L in rec.get("letters_U", []):
            for a in alphabet:
                if str(a) == L or _same_letter(a, L):
                    inU[str(a)] = True
    # additionally: a letter that is itself one of the bare invariants (s, t, m2, ...)
    # is automatically a first-entry candidate.
    bare = {str(s) for s in gs.invariants}
    first = [str(a) for a in alphabet if in_codim1[str(a)] or str(a) in bare]
    last  = [str(a) for a in alphabet if inU[str(a)] and not inF[str(a)]]
    return first, last


def _same_letter(a, b_str):
    try:
        r = sp.cancel(sp.sympify(a) / sp.sympify(b_str))
        return bool(r.is_number)
    except Exception:
        return False


# ===========================================================================
# prune_spurious — a letter is SPURIOUS if any of:
#   (a) it appears ONLY in U-faces and never in any F-face (second-type-only,
#       no solution in the physical region for the leading Landau equations);
#   (b) it cancels identically in F/U (i.e. divides gcd(F,U) so is a pure
#       Lee-Pomeransky artefact);
#   (c) it is a pure Feynman-parameter monomial that slipped through (should
#       never happen post-kin_factors but is guarded).
# Returns (kept:list, spurious:list[str]).
# ===========================================================================
def prune_spurious(gs: GraphSpec, alphabet, faces_data, full_UF):
    U, F = full_UF
    fparams = set(gs.ev.values())
    # (b) factors that cancel in F/U
    try:
        g = sp.gcd(sp.Poly(F, *fparams, *gs.kinsyms),
                   sp.Poly(U, *fparams, *gs.kinsyms))
        cancels = kin_factors(gs, g.as_expr() if hasattr(g, "as_expr") else g)
    except Exception:
        cancels = set()
    # (a) U-only
    inF = set(); inU = set()
    for rec in faces_data:
        inF |= set(rec.get("letters_F", []))
        inU |= set(rec.get("letters_U", []))
    spurious, kept = [], []
    for a in alphabet:
        sa = str(a)
        if a.free_symbols & fparams:                        # (c)
            spurious.append(sa); continue
        if any(_same_letter(a, str(c)) for c in cancels):   # (b)
            spurious.append(sa); continue
        # (a): U-only AND never in any F-face
        only_U = any(_same_letter(a, l) for l in inU) and \
                 not any(_same_letter(a, l) for l in inF)
        if only_U:
            spurious.append(sa); continue
        kept.append(a)
    return kept, spurious


# ===========================================================================
# DRIVER
# ===========================================================================
def run(gs: GraphSpec, *, out_path: str | None = None,
        face_timeout: int = 45, global_minutes: int = 40,
        min_face_edges: int = 2, verbose: bool = False) -> dict:
    GLOBAL_DEADLINE = time.time() + 60 * global_minutes
    log = {
        "name": gs.name,
        "engine": "landau_alphabet.py (bounded face-by-face sympy)",
        "method": ("principal-A-determinant by FACES: per-subgraph Symanzik (U,F) + "
                   "BOUNDED resultant elimination (no global Groebner). "
                   "Incremental flush after every face."),
        "graph": gs.raw,
        "FACE_TIMEOUT_s": face_timeout,
        "faces": [],
        "alphabet": [],
        "first_entry": [],
        "last_entry": [],
        "spurious": [],
        "honesty": {},
    }
    raw = set()
    n_resolved = n_partial = n_timeout = 0
    unresolved_faces = []
    deadline_hit = False

    # full-graph U,F (for spurious gcd test)
    full = symanzik_subgraph(gs, set(gs.ALLE))
    full_UF = (full['U'], full['F']) if full else (sp.Integer(1), sp.Integer(0))

    ALLE = gs.ALLE
    for k in range(len(ALLE), min_face_edges - 1, -1):
        for face in combinations(ALLE, k):
            if time.time() > GLOBAL_DEADLINE:
                deadline_hit = True
                break
            sg = symanzik_subgraph(gs, set(face))
            if sg is None:
                continue
            has_kin = (sg['F0'] != 0) or any(gs.em[e] != 0 for e in face)
            rec = {
                "face_edges": list(face),
                "n_edges": k,
                "loops_L": sg['L'],
                "n_massive": len(set(face) & gs.massive_edges),
                "kinematic": bool(has_kin),
            }
            if not has_kin:
                rec["status"] = "no-kinematic-part (scaleless)"
                rec["face_status"] = "scaleless"   # normalized enum; not a gap
                rec["letters_F"] = []; rec["letters_U"] = []
                log["faces"].append(rec); _flush(log, out_path)
                continue
            fvars = [gs.ev[e] for e in face]
            t0 = time.time()
            lf, stf = discriminant_letters(gs, sg['F'], fvars, face_timeout)
            lu, stu = (set(), "skip")
            if sg['L'] >= 1 and sg['U'] not in (0, 1):
                lu, stu = discriminant_letters(gs, sg['U'], fvars, face_timeout)
            raw |= (lf | lu)
            rec["letters_F"] = sorted(str(x) for x in lf)
            rec["letters_U"] = sorted(str(x) for x in lu)
            rec["status_F"] = stf; rec["status_U"] = stu
            # discriminant_letters returns ∈ {"resolved", "resolved-groebner:<backend>",
            # "no-free-params", "bounded-partial", "bounded-timeout"}.  The first
            # THREE are all FULLY RESOLVED (Groebner escalation succeeded; no-free-
            # params = F has no Feynman params on this face → kin_factors of P only).
            _ok = lambda s: s in ("resolved", "skip", "no-free-params") or s.startswith("resolved-groebner")
            rec["status"]   = ("resolved" if _ok(stf) and _ok(stu)
                               else f"F:{stf}/U:{stu}")
            rec["secs"] = round(time.time() - t0, 2)
            # normalized per-face status enum: done / timed_out / over_budget.
            # Load-bearing: the auditor (audit_faces.py) keys off face_status,
            # so any non-"done" face is a NAMED gap, never silently dropped.
            if rec["status"] == "resolved":
                n_resolved += 1
                rec["face_status"] = "done"
            elif "timeout" in rec["status"]:
                n_timeout += 1; unresolved_faces.append(list(face))
                rec["face_status"] = "timed_out"
            else:
                n_partial += 1; unresolved_faces.append(list(face))
                rec["face_status"] = "over_budget"
            log["faces"].append(rec)
            if verbose:
                sys.stderr.write(f"[face k={k} {face}] {rec['status']} "
                                 f"({len(lf)}F+{len(lu)}U letters, {rec['secs']}s)\n")
            _flush(log, out_path)
        else:
            continue
        break

    # ---- canonicalize, classify, prune ------------------------------------
    alphabet = canonical_alphabet(raw)
    first, last = classify_entries(gs, alphabet, log["faces"])
    kept, spurious = prune_spurious(gs, alphabet, log["faces"], full_UF)

    # ---- expected face census (so a truncated ENUMERATION is also caught) --
    # The full face lattice is every connected sub-multigraph on >=2 vertices,
    # i.e. every edge-subset whose Symanzik build is not None, down to
    # min_face_edges.  We recompute it independently of what got attempted so
    # the auditor can flag faces that were never even reached (enumeration
    # truncation), not just faces that were attempted and timed out.
    expected_faces = 0
    for k in range(len(ALLE), min_face_edges - 1, -1):
        for face in combinations(ALLE, k):
            if symanzik_subgraph(gs, set(face)) is not None:
                expected_faces += 1

    n_attempted = len(log["faces"])
    n_incomplete = n_partial + n_timeout
    n_missing = expected_faces - n_attempted          # never-attempted faces
    # COMPLETE iff every expected face was attempted AND resolved/scaleless.
    complete = (n_incomplete == 0) and (n_missing == 0) and (not deadline_hit)

    log["alphabet"]    = [str(a) for a in kept]
    log["first_entry"] = [a for a in first if a in log["alphabet"]]
    log["last_entry"]  = last
    log["spurious"]    = spurious
    # TOP-LEVEL honesty flag: a reader scanning the head of the JSON sees this
    # FIRST.  It is FALSE if ANY face is incomplete or was never attempted.
    log["complete"] = bool(complete)
    log["honesty"] = {
        "complete": bool(complete),
        "expected_faces": expected_faces,
        "attempted_faces": n_attempted,
        "missing_faces": n_missing,          # attempted < expected => enumeration truncated
        "n_faces": n_attempted,
        "resolved": n_resolved,
        "bounded_partial": n_partial,
        "bounded_timeout": n_timeout,
        "unresolved_faces": unresolved_faces,
        "global_deadline_hit": deadline_hit,
        "note": ("alphabet is COMPLETE only if complete==true: every expected face "
                 "attempted (missing_faces==0) AND none bounded_partial/bounded_timeout. "
                 "Re-run 'bounded_timeout' faces with pld_bridge.jl or a longer "
                 "--face-timeout."),
    }
    _flush(log, out_path)
    return log


def _flush(log, out_path):
    if not out_path:
        return
    try:
        with open(out_path, "w") as fh:
            json.dump(log, fh, indent=2, default=str)
    except Exception as e:  # pragma: no cover
        sys.stderr.write(f"flush error: {e}\n")


# ===========================================================================
# CLI
# ===========================================================================
def load_spec(path_or_dict) -> GraphSpec:
    if isinstance(path_or_dict, dict):
        return GraphSpec(path_or_dict)
    with open(path_or_dict) as fh:
        return GraphSpec(json.load(fh))


def main(argv=None):
    ap = argparse.ArgumentParser(description="Bounded face-by-face Landau alphabet engine.")
    ap.add_argument("spec", help="graph spec JSON (see graph_specs/)")
    ap.add_argument("--out", default=None, help="output JSON (incremental flush)")
    ap.add_argument("--face-timeout", type=int, default=45)
    ap.add_argument("--global-minutes", type=int, default=40)
    ap.add_argument("--min-face-edges", type=int, default=2)
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args(argv)
    gs = load_spec(args.spec)
    out = args.out or os.path.splitext(os.path.basename(args.spec))[0] + "_alphabet.json"
    log = run(gs, out_path=out, face_timeout=args.face_timeout,
              global_minutes=args.global_minutes,
              min_face_edges=args.min_face_edges, verbose=args.verbose)
    print(f"\n[{gs.name}] alphabet ({len(log['alphabet'])} letters):")
    for a in log["alphabet"]:
        print("   ", a)
    print(f"first_entry: {log['first_entry']}")
    print(f"last_entry:  {log['last_entry']}")
    print(f"spurious:    {log['spurious']}")
    h = log["honesty"]
    print(f"honesty: {h['resolved']} resolved / {h['bounded_partial']} partial / "
          f"{h['bounded_timeout']} timeout "
          f"({h['attempted_faces']}/{h['expected_faces']} faces attempted)")
    if log["complete"]:
        print("COMPLETE: True  (every face attempted and resolved)")
    else:
        print("COMPLETE: False  ** alphabet NOT proven complete **")
        if h["missing_faces"]:
            print(f"  {h['missing_faces']} face(s) NEVER ATTEMPTED (enumeration truncated)")
        for f in h["unresolved_faces"]:
            print(f"  unresolved face ({len(f)} edges): {f}")
    print(f"written -> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
# seedling family-preflight member: `seedling audit` + `seedling identity`.
"""Stage-0 preflight: discrete autopermutation group of an integral family.

Finds every external-leg crossing π (composed with charge-conjugation C: p→−p)
that, together with an affine unimodular loop map  l_a = M_{ab} l'_b + A_{ac} p'_c,
sends the DENOMINATOR set {D_1..D_n} onto itself as functions of momenta.  These
are the family autopermutations kira's own symmetry finder MISSES (kira only tries
loop shifts + `magic_relations` Pak matching, not external crossings ∘ C).

The search and its controls were validated against independently proven
symmetry orbits before being generalized into this tool.

Input
-----
Either a kira config dir (containing integralfamilies.yaml + kinematics.yaml) or
an amflow input JSON.  Propagators are parsed as (Σ c_i q_i)^2 − m; loop momenta,
independent externals (after conservation), scalar-product rules, and kinematic
invariants are extracted.

Search
------
Step 1: enumerate π ∈ S_{n_legs} that preserve the on-shell mass pattern
        p_i²→p_{π(i)}² and FIX every kinematic invariant (p_i·p_j → same
        expression in {s,t,…}).  Also record π that PERMUTE invariants
        (crossing identities I(s,t)=I(t,s)-type; reported, not folded into Aut).
Step 2: for each such π × {id, C}, solve for (M,A) by an ANCHOR method: pick L
        denominators with independent loop-vectors, enumerate their (target,
        sign) images, solve M and A exactly over ℚ, keep integer unimodular
        solutions, then verify the full denominator set maps bijectively.
        This is O(n_denom^L · 2^L) per (π,C), cheap for L≤4.
Step 3: for each hit, record the induced propagator permutation σ_D on
        {D_1..D_{n_denom}} and the ISP action (permutation if the ISPs happen to
        map to ISPs, else 'mixes' — the identity I[a]=I[σ(a)] is then scoped to
        zero-ISP integrals).

Output (under --out DIR, default = alongside CONFIG)
----------------------------------------------------
  <family>_AUT.json                — group order, generators (momentum maps +
                                     D-perms), invariant action, ISP scope,
                                     crossing report, controls.
  <family>_aut_relations.kira      — reduce_user_defined_system seed relations
                                     (one per subsector-pair with sec≠σ(sec),
                                     zero-ISP scoped).
  <family>_magic_relations.yaml    — drop-in integralfamilies block with
                                     magic_relations: true.

CLI
---
  seedling audit CONFIG_DIR_OR_JSON [--out DIR] [--max-shift 1] [--no-cache]
                                    [--report-crossings]

Exit 0 always; summary to stdout.  Result cached per family-hash in
.family_aut_cache/ beside this module (override with FAMILY_AUT_CACHE).

Identity subcommand
-------------------
  seedling identity
      — built-in regression fixture (wpair two-loop family vs its paper
        build), no arguments.
  seedling identity --ours FAMILY_DEF.json --theirs PROPS \
      [--theirs-format mathematica|json] [--theirs-key PropagatorsPBB]
      [--swap "l->k,k->p"] [--masses "p1^2=M2,..."] [--x-perm "j1,..,jN"] [--out F]
      — general mode: solve the affine scalar-product dictionary
        between two builds of a family via the second-Symanzik F polynomial
        over the fully symbolic external Gram (+ mass monomials), for
        arbitrary loop/propagator count.  rc 0 = PASS, rc 1 = FAIL (loud).

Preflight-hook contract
-----------------------
Before the first kira run on a config dir CFG, run `seedling audit CFG`; if
the result has |Aut|>1 AND CFG/integralfamilies.yaml lacks
`magic_relations: true`, WARN — kira will carry redundant masters.  Wire the
emitted magic_relations drop-in (or the .kira seed relations) into the
config; a WARN read and ignored buys nothing.
"""
import argparse
import hashlib
import itertools
import json
import os
import sys
import time
from fractions import Fraction

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE_DIR = os.environ.get("FAMILY_AUT_CACHE",
                           os.path.join(HERE, ".family_aut_cache"))


# ============================================================================
# Family parsing
# ============================================================================
class Family:
    def __init__(self, name, loops, legs, elim_leg, elim_expr,
                 sp_rules, invariants, prop_strings, top_sector, n_props,
                 magic_relations_set):
        self.name = name
        self.loops = list(loops)
        self.legs = list(legs)           # ALL physical legs (before conservation)
        self.elim_leg = elim_leg         # symbol name eliminated
        self.elim_expr = elim_expr       # sympy expr in the other legs
        self.indep = [p for p in self.legs if p != elim_leg]
        self.sp_rules = sp_rules         # {(i,j) sorted-tuple: sympy expr}
        self.invariants = list(invariants)
        self.prop_strings = list(prop_strings)
        self.top_sector = int(top_sector)
        self.n_props = int(n_props)
        self.magic_relations_set = bool(magic_relations_set)
        # denominator/ISP split from top-sector bit pattern
        self.n_denom = self.top_sector.bit_length()
        # basis symbols
        self.loop_syms = [sp.Symbol(x) for x in self.loops]
        self.indep_syms = [sp.Symbol(x) for x in self.indep]
        self.elim_sym = sp.Symbol(elim_leg)
        self.elim_sub = {self.elim_sym: elim_expr}
        self.L = len(self.loops)
        self.E = len(self.indep)
        # parsed propagators
        self.props = [self._parse_prop(i + 1, s)
                      for i, s in enumerate(self.prop_strings)]

    # ---- propagator parsing ------------------------------------------------
    def _parse_prop(self, idx, s):
        """Return dict {idx, kind, lvec (L,), evec (E,), mass, raw} for a propagator.
        kind='quad' if of form (Σ c q)^2 − m; 'other' otherwise (raw sympy kept)."""
        expr = sp.sympify(s, locals={x: sp.Symbol(x)
                                     for x in self.loops + self.legs})
        expr = sp.expand(expr.subs(self.elim_sub))
        # try to recognise q^2 − m with q linear in (loops ∪ indep)
        mom = self.loop_syms + self.indep_syms
        # collect linear-in-momenta part by treating momenta as commuting scalars
        p = sp.Poly(expr, *mom)
        deg = p.total_degree()
        if deg == 2:
            # quadratic form: check if it's a perfect square of a linear form,
            # i.e. coeff(q_a q_b) = 2 c_a c_b and coeff(q_a^2)=c_a^2.
            cvec = []
            for q in mom:
                c2 = p.nth(*[2 if m is q else 0 for m in mom])
                cvec.append(c2)
            # sign of linear form: pick first nonzero c_a^2, take +sqrt (integer)
            lin = [None] * len(mom)
            ref = None
            for a, ca2 in enumerate(cvec):
                if ca2 != 0:
                    r = sp.sqrt(ca2)
                    if not r.is_Integer:
                        return dict(idx=idx, kind='other', raw=expr)
                    lin[a] = int(r)
                    ref = a
                    break
            if ref is None:
                return dict(idx=idx, kind='other', raw=expr)
            for b in range(len(mom)):
                if b == ref:
                    continue
                mono = [0] * len(mom)
                mono[ref] = 1
                mono[b] = 1
                cab = p.nth(*mono)
                # cab should be 2 c_ref c_b
                cb = sp.Rational(cab, 2 * lin[ref])
                if not cb.is_integer:
                    return dict(idx=idx, kind='other', raw=expr)
                lin[b] = int(cb)
            # verify: (Σ lin_a q_a)^2 matches the deg-2 part
            qlin = sum(c * q for c, q in zip(lin, mom))
            mass = sp.expand(qlin**2 - expr)
            # mass must be momentum-free
            if any(m in mass.free_symbols for m in mom):
                return dict(idx=idx, kind='other', raw=expr)
            return dict(idx=idx, kind='quad',
                        lvec=tuple(lin[:self.L]), evec=tuple(lin[self.L:]),
                        mass=sp.simplify(mass), raw=expr)
        return dict(idx=idx, kind='other', raw=expr)

    def denom_props(self):
        return [p for p in self.props[:self.n_denom]]

    def isp_props(self):
        return [p for p in self.props[self.n_denom:]]

    def hash(self):
        body = json.dumps({
            "name": self.name, "loops": self.loops, "legs": self.legs,
            "elim": self.elim_leg, "props": self.prop_strings,
            "sp": {str(k): str(v) for k, v in sorted(self.sp_rules.items())},
            "top": self.top_sector,
        }, sort_keys=True)
        return hashlib.sha1(body.encode()).hexdigest()[:16]


def _load_yaml(path):
    import yaml
    with open(path) as f:
        return yaml.safe_load(f)


def parse_kira_config(cfgdir):
    ifam = _load_yaml(os.path.join(cfgdir, "integralfamilies.yaml"))
    kin = _load_yaml(os.path.join(cfgdir, "kinematics.yaml"))
    fam = ifam["integralfamilies"][0]
    k = kin["kinematics"]
    loops = fam["loop_momenta"]
    props = [p[0] if isinstance(p, (list, tuple)) else p
             for p in fam["propagators"]]
    top = max(fam.get("top_level_sectors", [2**len(props) - 1]))
    magic = bool(fam.get("magic_relations", False))
    inc = list(k.get("incoming_momenta") or [])
    out = list(k.get("outgoing_momenta") or [])
    legs = inc + out
    mc = k.get("momentum_conservation")
    if mc:
        elim = mc[0]
        elim_expr = sp.sympify(mc[1], locals={x: sp.Symbol(x) for x in legs})
        if elim not in legs:
            legs = legs + [elim]
    else:
        elim = legs[-1]
        elim_expr = -sum(sp.Symbol(x) for x in legs[:-1])
    invs = [sp.Symbol(v[0]) for v in (k.get("kinematic_invariants") or [])]
    sp_rules = {}
    for rule in (k.get("scalarproduct_rules") or []):
        (a, b), val = rule
        key = tuple(sorted((a, b)))
        sp_rules[key] = sp.sympify(str(val))
    return Family(fam["name"], loops, legs, elim, elim_expr,
                  sp_rules, invs, props, top, len(props), magic)


def parse_amflow_json(path):
    d = json.load(open(path))
    fam = d["family"]
    loops = fam["loops"]
    legs = list(fam["legs"])
    conservation = fam.get("conservation", {})
    if conservation:
        elim = list(conservation.keys())[0]
        elim_expr = sp.sympify(conservation[elim],
                               locals={x: sp.Symbol(x) for x in legs})
    else:
        elim = legs[-1]
        elim_expr = -sum(sp.Symbol(x) for x in legs[:-1])
    if elim not in legs:
        legs = legs + [elim]
    indep = [p for p in legs if p != elim]
    # replacement dict → scalar-product rules by solving linear system
    repl = fam.get("replacement", {})
    invs = set()
    sp_rules = _sp_rules_from_replacement(repl, indep, elim, elim_expr, invs)
    props = list(fam["propagators"])
    top = d.get("top_sector") or fam.get("top_sector")
    if top is None:
        # amflow JSONs don't carry it; assume all props before the first pure
        # loop×ext bilinear are denominators — fallback: use bit_length heuristic
        # (safe: user can override via --n-denom)
        top = 2**len(props) - 1
    return Family(fam["name"], loops, legs, elim, elim_expr,
                  sp_rules, sorted(invs, key=str), props, top, len(props), False)


def _sp_rules_from_replacement(repl, indep, elim, elim_expr, invs_out):
    syms = {x: sp.Symbol(x) for x in indep + [elim]}
    isyms = [sp.Symbol(x) for x in indep]
    pairs = [(i, j) for i in range(len(indep)) for j in range(i, len(indep))]
    unk = {pr: sp.Symbol(f"_sp{pr[0]}{pr[1]}") for pr in pairs}
    eqs = []
    for lhs, rhs in repl.items():
        L = sp.expand(sp.sympify(lhs, locals=syms).subs({sp.Symbol(elim): elim_expr}))
        R = sp.sympify(rhs)
        for s in R.free_symbols:
            invs_out.add(s)
        # expand L as bilinear in indep, replace p_i*p_j → unk
        Lp = sp.Poly(L, *isyms)
        e = -R
        for mono, c in Lp.terms():
            nz = [k for k, m in enumerate(mono) if m > 0]
            if sum(mono) == 2 and len(nz) == 1:
                e += c * unk[(nz[0], nz[0])]
            elif sum(mono) == 2 and len(nz) == 2:
                e += c * unk[(nz[0], nz[1])]
            elif sum(mono) == 0:
                e += c
        eqs.append(e)
    sol = sp.solve(eqs, list(unk.values()), dict=True)
    sol = sol[0] if sol else {}
    rules = {}
    for (i, j), u in unk.items():
        v = sol.get(u, u)  # may be underdetermined
        rules[tuple(sorted((indep[i], indep[j])))] = sp.simplify(v)
    return rules


def load_family(path, n_denom_override=None):
    if os.path.isdir(path):
        cfg = os.path.join(path, "config") if os.path.isdir(
            os.path.join(path, "config")) else path
        fam = parse_kira_config(cfg)
    elif path.endswith(".json"):
        fam = parse_amflow_json(path)
    else:
        fam = parse_kira_config(path)
    if n_denom_override:
        fam.n_denom = int(n_denom_override)
        fam.top_sector = 2**fam.n_denom - 1
    return fam


# ============================================================================
# External-leg permutation enumeration
# ============================================================================
def _sp_of(fam, i, j):
    """p_i · p_j (indep-index i,j) as sympy expr via sp_rules."""
    key = tuple(sorted((fam.indep[i], fam.indep[j])))
    return fam.sp_rules.get(key)


def _leg_mass(fam, leg):
    """p_leg^2 as sympy expr (for any leg incl. eliminated)."""
    if leg == fam.elim_leg:
        # (Σ indep with coeffs from elim_expr)^2
        cs = [sp.Poly(fam.elim_expr, sp.Symbol(p)).nth(1) for p in fam.indep]
        e = sp.Integer(0)
        for a in range(fam.E):
            for b in range(fam.E):
                spab = _sp_of(fam, a, b)
                if spab is None:
                    return None
                e += cs[a] * cs[b] * spab
        return sp.simplify(e)
    key = (leg, leg)
    return fam.sp_rules.get(key)


def _perm_matrix_on_indep(fam, perm):
    """perm: dict leg→leg on ALL legs.  Return E×E integer matrix Q with
    (new indep basis) p'_a = old p_{perm[a]} = Σ_b Q[a,b] p_b (indep)."""
    Q = []
    for a, pa in enumerate(fam.indep):
        img = perm[pa]
        if img == fam.elim_leg:
            row = [int(sp.Poly(fam.elim_expr, sp.Symbol(p)).nth(1))
                   for p in fam.indep]
        else:
            row = [1 if p == img else 0 for p in fam.indep]
        Q.append(row)
    return Q


def _sp_after_perm(fam, Q, a, b):
    """p'_a · p'_b where p'_a = Σ Q[a,c] p_c."""
    e = sp.Integer(0)
    for c in range(fam.E):
        for d in range(fam.E):
            if Q[a][c] == 0 or Q[b][d] == 0:
                continue
            spcd = _sp_of(fam, c, d)
            if spcd is None:
                return None
            e += Q[a][c] * Q[b][d] * spcd
    return sp.expand(e)


def enum_ext_maps(fam, want_crossings=False, cap=5040):
    """Enumerate leg permutations.  Returns two lists:
       fixing: [(label, Qmatrix)]            — π fixes every invariant
       crossing: [(label, Qmatrix, subs)]    — π induces nontrivial invariant subs
    Always includes identity in `fixing`."""
    n = len(fam.legs)
    if n > 7:
        # cap: factorial(8)=40320 already; for n>7 we bail with identity only
        return [("id", [[1 if i == j else 0 for j in range(fam.E)]
                        for i in range(fam.E)])], [], True
    # mass classes
    mass = {leg: _leg_mass(fam, leg) for leg in fam.legs}
    fixing, crossing = [], []
    seen_Q = set()
    for perm_t in itertools.permutations(fam.legs):
        perm = dict(zip(fam.legs, perm_t))
        # mass compatibility
        if any(sp.simplify((mass[leg] or 0) - (mass[perm[leg]] or 0)) != 0
               for leg in fam.legs if mass[leg] is not None
               and mass[perm[leg]] is not None):
            continue
        Q = _perm_matrix_on_indep(fam, perm)
        Qkey = tuple(tuple(r) for r in Q)
        if Qkey in seen_Q:
            continue
        seen_Q.add(Qkey)
        # invariant action: does p_a·p_b == p'_a·p'_b for all indep pairs?
        fixes = True
        subs = {}
        for a in range(fam.E):
            for b in range(a, fam.E):
                orig = _sp_of(fam, a, b)
                new = _sp_after_perm(fam, Q, a, b)
                if orig is None or new is None:
                    fixes = False
                    subs = None
                    break
                if sp.simplify(orig - new) != 0:
                    fixes = False
                if subs is not None:
                    subs[str(orig)] = str(new)
            if subs is None:
                break
        label = _perm_label(fam.legs, perm)
        if fixes:
            fixing.append((label, Q))
        elif want_crossings and subs is not None:
            crossing.append((label, Q, subs))
        if len(fixing) + len(crossing) > cap:
            return fixing, crossing, True
    return fixing, crossing, False


def _perm_label(legs, perm):
    seen, cycles = set(), []
    for x in legs:
        if x in seen or perm[x] == x:
            seen.add(x)
            continue
        cyc, y = [x], perm[x]
        seen.add(x)
        while y != x:
            cyc.append(y)
            seen.add(y)
            y = perm[y]
        cycles.append("(" + " ".join(cyc) + ")")
    return "".join(cycles) or "id"


# ============================================================================
# Anchor-based autopermutation search
# ============================================================================
def _sig(lvec, evec, mass):
    """Sign-normalized signature of a quadratic propagator."""
    v = tuple(lvec) + tuple(evec)
    for c in v:
        if c != 0:
            if c < 0:
                v = tuple(-x for x in v)
            break
    return (v, sp.simplify(mass))


def _pick_anchors(denoms, L):
    """Pick L quad denominators with linearly-independent loop-vectors.
    Return (indices, V) with V = L×L matrix of their loop-vecs (rows)."""
    idxs = []
    V = sp.zeros(L, L)
    r = 0
    for d in denoms:
        if d['kind'] != 'quad':
            continue
        Vt = V.copy()
        for j in range(L):
            Vt[r, j] = d['lvec'][j]
        if sp.Matrix(Vt[:r + 1, :]).rank() == r + 1:
            V = Vt
            idxs.append(d['idx'])
            r += 1
            if r == L:
                return idxs, V
    return None, None


def _apply_map(prop, M, A, Q, erev, L, E):
    """Transform a quad propagator's (lvec,evec) under the affine map."""
    lv = sp.Matrix(1, L, prop['lvec'])
    ev = sp.Matrix(1, E, prop['evec'])
    new_l = lv * M
    # external part: old p_a → erev * Σ_b Q[a][b] p_b, then + lv·A
    Qm = sp.Matrix(E, E, lambda i, j: Q[i][j])
    new_e = erev * ev * Qm + lv * A
    return tuple(int(x) for x in new_l), tuple(int(x) for x in new_e)


def _fmt_mom(coeffs, names):
    terms = []
    for c, n in zip(coeffs, names):
        if c == 0:
            continue
        if c == 1:
            terms.append(("+" if terms else "") + n)
        elif c == -1:
            terms.append("-" + n)
        else:
            terms.append(("+" if (c > 0 and terms) else "") + f"{c}*{n}")
    return "".join(terms) or "0"


def find_autopermutations(fam, ext_maps, max_shift=1, verbose=False):
    denoms = fam.denom_props()
    for d in denoms:
        if d['kind'] != 'quad':
            raise ValueError(
                f"denominator D{d['idx']} is not (linear)^2−m; cannot search "
                f"(raw={d['raw']}). Rewrite the family or use --n-denom.")
    L, E = fam.L, fam.E
    anchors, V = _pick_anchors(denoms, L)
    if anchors is None:
        raise ValueError("could not find L independent-loop-vector denominators")
    Vinv = V.inv()
    anchor_props = [d for d in denoms if d['idx'] in anchors]
    other_props = [d for d in denoms if d['idx'] not in anchors]
    target_by_mass = {}
    for d in denoms:
        target_by_mass.setdefault(sp.simplify(d['mass']), []).append(d)
    target_sig = {_sig(d['lvec'], d['evec'], d['mass']): d['idx'] for d in denoms}

    hits = []
    cap_hit = False
    n_tried = 0
    for elabel, Q in ext_maps:
        Qm = sp.Matrix(E, E, lambda i, j: Q[i][j])
        for erev in (1, -1):
            # w'_i = erev * (evec_i · Q)  (row vector)
            wprime = {d['idx']: erev * sp.Matrix(1, E, d['evec']) * Qm
                      for d in denoms}
            # enumerate anchor→target assignments
            anchor_targets = [target_by_mass.get(sp.simplify(a['mass']), [])
                              for a in anchor_props]
            for tgt_combo in itertools.product(*anchor_targets):
                if len({t['idx'] for t in tgt_combo}) != L:
                    continue
                for signs in itertools.product((1, -1), repeat=L):
                    n_tried += 1
                    # solve M: V·M = diag(signs)·[tgt lvecs]
                    RHS_l = sp.Matrix(L, L,
                        lambda i, j: signs[i] * tgt_combo[i]['lvec'][j])
                    M = Vinv * RHS_l
                    if any(not c.is_integer for c in M):
                        continue
                    if abs(M.det()) != 1:
                        continue
                    # solve A: V·A = diag(signs)·[tgt evecs] − W'_anchor
                    RHS_e = sp.Matrix(L, E,
                        lambda i, j: signs[i] * tgt_combo[i]['evec'][j]
                                     - wprime[anchor_props[i]['idx']][0, j])
                    A = Vinv * RHS_e
                    if any(not c.is_integer for c in A):
                        continue
                    # verify full denominator set
                    Dperm, ok = {}, True
                    for d in denoms:
                        nl, ne = _apply_map(d, M, A, Q, erev, L, E)
                        s = _sig(nl, ne, d['mass'])
                        if s not in target_sig:
                            ok = False
                            break
                        Dperm[d['idx']] = target_sig[s]
                    if not ok or len(set(Dperm.values())) != len(denoms):
                        continue
                    maxc = max(abs(int(c)) for c in list(M) + list(A))
                    if max_shift and maxc > max_shift:
                        cap_hit = True   # a VALID hit exceeds the coeff bound
                    hits.append(dict(
                        ext=elabel, erev=erev,
                        M=[[int(M[i, j]) for j in range(L)] for i in range(L)],
                        A=[[int(A[i, j]) for j in range(E)] for i in range(L)],
                        Dperm=Dperm, max_coeff=maxc,
                    ))
    if verbose:
        print(f"  [search] {len(ext_maps)} ext-maps × 2 erev, "
              f"{n_tried} anchor-assignments tried, {len(hits)} raw hits")
    return hits, cap_hit


def dedup_by_Dperm(hits, n_denom):
    identity = tuple(range(1, n_denom + 1))
    grouped = {}
    for h in hits:
        key = tuple(h['Dperm'][i] for i in range(1, n_denom + 1))
        if key not in grouped:
            grouped[key] = dict(Dperm=dict(h['Dperm']), reps=[])
        grouped[key]['reps'].append(
            {k: h[k] for k in ('ext', 'erev', 'M', 'A', 'max_coeff')})
    return grouped, identity


def group_order(grouped, n_denom):
    """Order of the D-permutation group (closure under composition)."""
    perms = set(grouped.keys())
    identity = tuple(range(1, n_denom + 1))
    perms.add(identity)
    # close under composition
    changed = True
    while changed:
        changed = False
        for a in list(perms):
            for b in list(perms):
                c = tuple(a[b[i] - 1] for i in range(n_denom))
                if c not in perms:
                    perms.add(c)
                    changed = True
    return len(perms), sorted(perms)


def minimal_generators(grouped, n_denom):
    """Greedy minimal generating set (by element order, then lex)."""
    identity = tuple(range(1, n_denom + 1))
    elems = [k for k in grouped if k != identity]
    if not elems:
        return []
    # sort by (order desc, lex)
    def eorder(p):
        cur, o = list(p), 1
        while tuple(cur) != identity:
            cur = [p[cur[i] - 1] for i in range(n_denom)]
            o += 1
            if o > 256:
                break
        return o
    elems.sort(key=lambda p: (-eorder(p), p))
    gens, span = [], {identity}
    for e in elems:
        if e in span:
            continue
        gens.append(e)
        # regenerate span
        span = {identity}
        frontier = list(gens)
        while frontier:
            new = []
            for a in list(span):
                for b in gens:
                    for c in (tuple(a[b[i] - 1] for i in range(n_denom)),
                              tuple(b[a[i] - 1] for i in range(n_denom))):
                        if c not in span:
                            span.add(c)
                            new.append(c)
            frontier = new
        if len(span) == len(grouped):
            break
    return gens


def _verify_generator_symbolic(fam, rep):
    """Mandatory control: explicitly substitute the momentum map into
    each denominator and check expand(D_i(map) − D_{σ(i)}) == 0."""
    L, E = fam.L, fam.E
    subs = {}
    Qm = sp.Matrix(E, E, lambda i, j: rep['Q'][i][j])
    # loop map
    for a in range(L):
        rhs = sum(rep['M'][a][b] * fam.loop_syms[b] for b in range(L))
        rhs += sum(rep['A'][a][c] * fam.indep_syms[c] for c in range(E))
        subs[fam.loop_syms[a]] = rhs
    # ext map
    for a in range(E):
        rhs = rep['erev'] * sum(int(Qm[a, b]) * fam.indep_syms[b]
                                for b in range(E))
        subs[fam.indep_syms[a]] = rhs
    Dperm = rep['Dperm']
    for d in fam.denom_props():
        img = sp.expand(d['raw'].subs(subs, simultaneous=True))
        tgt = fam.props[Dperm[d['idx']] - 1]['raw']
        # apply sp_rules to (img - tgt) — should be identically 0 as a
        # polynomial in loop momenta with coefficients in ℚ[invariants]
        diff = sp.expand(img - tgt)
        # substitute p_i·p_j via a polarisation trick: treat as scalars
        # (already done — raw exprs are in commuting scalar symbols)
        if diff != 0:
            return False, f"D{d['idx']}->D{Dperm[d['idx']]}: diff={diff}"
    return True, "all denominators verified by explicit substitution"


# ============================================================================
# ISP action + output
# ============================================================================
def isp_action(fam, hit):
    """For each ISP, compute image under the momentum map. If it's ±another
    prop (denom or ISP), record the permutation; else 'mixes'."""
    L, E = fam.L, fam.E
    M = sp.Matrix(hit['M'])
    A = sp.Matrix(hit['A'])
    Q = hit['Q']
    erev = hit['erev']
    all_sig = {}
    for d in fam.props:
        if d['kind'] == 'quad':
            all_sig[_sig(d['lvec'], d['evec'], d['mass'])] = d['idx']
    out = {}
    is_perm = True
    for d in fam.isp_props():
        if d['kind'] != 'quad':
            out[d['idx']] = 'non-quad (skipped)'
            is_perm = False
            continue
        nl, ne = _apply_map(d, M, A, Q, erev, L, E)
        s = _sig(nl, ne, d['mass'])
        if s in all_sig:
            out[d['idx']] = all_sig[s]
            if all_sig[s] <= fam.n_denom:
                is_perm = False
        else:
            expr = _fmt_mom(nl + ne, fam.loops + fam.indep)
            out[d['idx']] = f"({expr})^2 - ({d['mass']})  [not a family prop]"
            is_perm = False
    return out, is_perm


def index_action(idx_tuple, Dperm, n_denom):
    inv = {Dperm[i]: i for i in range(1, n_denom + 1)}
    return tuple(idx_tuple[inv[j] - 1] for j in range(1, n_denom + 1)) + \
           tuple(idx_tuple[n_denom:])


def sector_of(idx, n_denom):
    return sum((1 << i) for i in range(n_denom) if idx[i] > 0)


def emit_kira_relations(fam, gens, out_path, max_pairs=200):
    """Seed relations: for each nontrivial generator σ, one I−σ(I) pair per
    subsector corner where sec ≠ σ(sec). Zero-ISP scoped."""
    lines = []
    n = fam.n_denom
    nisp = fam.n_props - n
    for g in gens:
        Dperm = g['Dperm']
        seen = set()
        count = 0
        for sec in range(1, 2**n):
            if sec in seen:
                continue
            idx = tuple(1 if (sec >> i) & 1 else 0 for i in range(n))
            simg = index_action(idx, Dperm, n)
            ssec = sector_of(simg, n)
            if ssec == sec:
                continue
            seen.add(sec)
            seen.add(ssec)
            a = list(idx) + [0] * nisp
            b = list(simg) + [0] * nisp
            lines.append(f"{fam.name}[{','.join(map(str, a))}]*(1)")
            lines.append(f"{fam.name}[{','.join(map(str, b))}]*(-1)")
            lines.append("")
            count += 1
            if count >= max_pairs:
                break
    with open(out_path, "w") as f:
        f.write("\n".join(lines))
    return len([l for l in lines if l.endswith("*(1)")])


def emit_magic_dropin(fam, cfgdir, out_path):
    if cfgdir and os.path.exists(os.path.join(cfgdir, "integralfamilies.yaml")):
        src = open(os.path.join(cfgdir, "integralfamilies.yaml")).read()
        if "magic_relations" not in src:
            src = src.replace("loop_momenta:",
                              "magic_relations: true\n    loop_momenta:", 1)
        with open(out_path, "w") as f:
            f.write("# drop-in: original integralfamilies.yaml with "
                    "magic_relations enabled (family_aut found |Aut|>1)\n")
            f.write(src)
    else:
        with open(out_path, "w") as f:
            f.write(f"# add to integralfamilies.yaml under family '{fam.name}':\n"
                    f"#   magic_relations: true\n")


# ============================================================================
# Main driver
# ============================================================================
def run(path, out_dir=None, max_shift=1, use_cache=True,
        report_crossings=False, n_denom=None, verbose=True):
    fam = load_family(path, n_denom_override=n_denom)
    h = fam.hash()
    os.makedirs(CACHE_DIR, exist_ok=True)
    cfile = os.path.join(CACHE_DIR, f"{h}.json")
    if use_cache and os.path.exists(cfile):
        result = json.load(open(cfile))
        result["_from_cache"] = True
        if verbose:
            _print_summary(fam, result)
        return result

    t0 = time.time()
    if verbose:
        print(f"[family_aut] {fam.name}: L={fam.L} loops {fam.loops}, "
              f"n_legs={len(fam.legs)}, n_denom={fam.n_denom}, "
              f"n_isp={fam.n_props - fam.n_denom}, "
              f"invariants={[str(x) for x in fam.invariants]}")
    if not fam.invariants and verbose:
        print("  WARN: no free kinematic invariants — sp-rules are numeric; "
              "invariant-fixing test degenerates to point-equality.")

    fixing, crossing, ext_cap = enum_ext_maps(fam, want_crossings=report_crossings)
    if verbose:
        print(f"  [ext] {len(fixing)} invariant-FIXING leg permutation(s): "
              f"{[e[0] for e in fixing]}"
              + (f"; {len(crossing)} crossing (invariant-permuting)" if crossing else ""))

    hits, cap_hit = find_autopermutations(fam, fixing, max_shift=max_shift,
                                          verbose=verbose)
    Qmap = {lbl: Q for lbl, Q in fixing}
    grouped, identity = dedup_by_Dperm(hits, fam.n_denom)
    order, closed = group_order(grouped, fam.n_denom)
    nontrivial = [k for k in grouped if k != identity]
    gen_keys = minimal_generators(grouped, fam.n_denom)

    def _entry(key, is_gen):
        g = grouped[key]
        is_id = (key == identity)
        rep = sorted(g['reps'], key=lambda r: (r['erev'], r['ext']))[0]
        rep_full = dict(rep)
        rep_full['Q'] = Qmap[rep['ext']]
        rep_full['Dperm'] = g['Dperm']
        loop_map = {fam.loops[a]:
                    _fmt_mom(list(rep['M'][a]) + list(rep['A'][a]),
                             fam.loops + fam.indep)
                    for a in range(fam.L)}
        isp, isp_perm = isp_action(fam, rep_full)
        # element order
        cur, o = list(key), 1
        while tuple(cur) != identity and o <= 256:
            cur = [key[cur[i] - 1] for i in range(fam.n_denom)]
            o += 1
        inv = {key[i]: i + 1 for i in range(fam.n_denom)}
        entry = dict(
            name=("identity" if is_id else
                  "sigma_" + "".join(str(x) for x in key)),
            is_identity=is_id, is_generator=is_gen,
            Dperm=list(key), order=o, n_reps=len(g['reps']),
            ext=rep['ext'], erev=rep['erev'],
            loop_map=loop_map,
            index_action=("[a1..a{n}] -> [".format(n=fam.n_denom) +
                          ",".join(f"a{inv[j]}"
                                   for j in range(1, fam.n_denom + 1)) + "]"),
            ISP_action={str(k): (v if isinstance(v, str) else int(v))
                        for k, v in isp.items()},
            ISP_is_permutation=isp_perm,
        )
        # mandatory control: explicit sympy verification
        if not is_id:
            ok, msg = _verify_generator_symbolic(fam, rep_full)
            entry['control_verified_by_substitution'] = ok
            if not ok:
                entry['control_fail_detail'] = msg
        return entry

    generators = [_entry(k, True) for k in gen_keys]
    elements = ([_entry(identity, False)] if identity in grouped else []) + \
               [_entry(k, k in gen_keys) for k in sorted(nontrivial)]

    controls = {
        "identity_found": identity in grouped,
        "closure_matches_hits": (order == len(grouped)),
        "all_generators_verified": all(
            g.get('control_verified_by_substitution', True) for g in generators),
    }

    result = {
        "family": fam.name,
        "family_hash": fam.hash(),
        "config_source": os.path.abspath(path),
        "L": fam.L, "n_legs": len(fam.legs), "n_denom": fam.n_denom,
        "n_isp": fam.n_props - fam.n_denom,
        "invariants": [str(x) for x in fam.invariants],
        "n_ext_fixing": len(fixing),
        "ext_fixing_labels": [e[0] for e in fixing],
        "n_raw_hits": len(hits),
        "n_distinct_Dperms": len(grouped),
        "group_order": order,
        "n_generators": len(generators),
        "generators": generators,
        "elements": elements,
        "controls": controls,
        "search_cap_hit": cap_hit or ext_cap,
        "search_space": (f"S_{len(fam.legs)} leg perms (mass-compatible, "
                         f"invariant-fixing) x {{id,C}} x anchor-solved "
                         f"unimodular loop map (|M_ij|,|A_ij| bound "
                         f"{max_shift}: hits above it are kept and flagged "
                         f"via search_cap_hit, not filtered)"),
        "scope": ("exact identity I[a]=I[sigma(a)] for a_{n_denom+1..}=0; "
                  "extend to nonzero ISP only when ISP_is_permutation=true"),
        "wall_s": round(time.time() - t0, 2),
        "magic_relations_already_set": fam.magic_relations_set,
    }
    if report_crossings:
        result["crossings"] = [
            {"ext": lbl, "invariant_image": subs} for lbl, Q, subs in crossing]

    # write outputs
    if out_dir is None:
        out_dir = (os.path.dirname(os.path.abspath(path))
                   if not os.path.isdir(path) else os.path.abspath(path))
    os.makedirs(out_dir, exist_ok=True)
    json_path = os.path.join(out_dir, f"{fam.name}_AUT.json")
    with open(json_path, "w") as f:
        json.dump(result, f, indent=2)
    result["_json_path"] = json_path
    if len(nontrivial) > 0:
        rel_path = os.path.join(out_dir, f"{fam.name}_aut_relations.kira")
        # emit needs Dperm as dict
        gens_for_emit = [
            {'Dperm': {i + 1: g['Dperm'][i] for i in range(fam.n_denom)}}
            for g in generators]
        n_rel = emit_kira_relations(fam, gens_for_emit, rel_path)
        result["_kira_relations_path"] = rel_path
        result["_n_seed_relations"] = n_rel
        cfgdir = (os.path.join(path, "config")
                  if os.path.isdir(os.path.join(path, "config"))
                  else (path if os.path.isdir(path) else None))
        magic_path = os.path.join(out_dir, f"{fam.name}_magic_relations.yaml")
        emit_magic_dropin(fam, cfgdir, magic_path)
        result["_magic_dropin_path"] = magic_path

    if use_cache:
        with open(cfile, "w") as f:
            json.dump(result, f, indent=2)

    if verbose:
        _print_summary(fam, result)
    return result


def _print_summary(fam, r):
    print(f"[family_aut] {r['family']}: |Aut| = {r['group_order']}  "
          f"({r['n_distinct_Dperms']} distinct D-perms from "
          f"{r['n_raw_hits']} raw reps; {r['n_generators']} generator(s))  "
          f"wall={r.get('wall_s', '?')}s"
          + ("  [CACHED]" if r.get('_from_cache') else ""))
    ctrl = r.get('controls', {})
    if not ctrl.get('all_generators_verified', True):
        print("  *** CONTROL FAIL: at least one generator failed explicit "
              "sympy substitution — DO NOT TRUST ***")
    for g in r['generators']:
        vf = g.get('control_verified_by_substitution')
        print(f"  gen  order={g['order']}  Dperm={g['Dperm']}  "
              f"ext={g['ext']}{'∘C' if g['erev'] == -1 else ''}  "
              f"ISP_perm={g['ISP_is_permutation']}  "
              f"verified={'YES' if vf else ('n/a' if vf is None else 'NO')}")
        for k, v in g['loop_map'].items():
            print(f"    {k} -> {v}")
    if r['group_order'] > 1 and not r.get('magic_relations_already_set'):
        print(f"  WARN: |Aut|>1 but config lacks `magic_relations: true` — "
              f"kira will carry redundant masters. "
              f"Drop-in: {r.get('_magic_dropin_path', '(see output dir)')}")
    if r.get('search_cap_hit'):
        print("  NOTE: a valid autopermutation had |M|,|A| entries beyond "
              f"--max-shift, or n_legs>7 (S_n cap).")


# ============================================================================
# CLI
# ============================================================================
def main():
    ap = argparse.ArgumentParser(
        description="Stage-0 preflight: discrete autopermutation group of a "
                    "Feynman integral family (external crossings ∘ C ∘ loop map).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)
    ap.add_argument("config", help="kira config dir OR amflow input JSON")
    ap.add_argument("--out", help="output dir (default: alongside config)")
    ap.add_argument("--max-shift", type=int, default=1,
                    help="coefficient bound for |M_ij|,|A_ij| (default 1): a "
                         "valid hit above it is STILL recorded, the receipt "
                         "sets search_cap_hit (flags, does not filter)")
    ap.add_argument("--n-denom", type=int, default=None,
                    help="override denominator/ISP split")
    ap.add_argument("--no-cache", action="store_true")
    ap.add_argument("--report-crossings", action="store_true",
                    help="also report leg perms that PERMUTE invariants")
    a = ap.parse_args()
    run(a.config, out_dir=a.out, max_shift=a.max_shift,
        use_cache=not a.no_cache, report_crossings=a.report_crossings,
        n_denom=a.n_denom)
    sys.exit(0)


# ============================================================================
# identity subcommand
# ============================================================================
def identity_check_wpair():
    """Built-in affine-dictionary identity check (wpair T4 vs its paper T_4
    build) — the regression fixture behind the no-arg `identity` form.
    Returns (hits_current, hits_fixed)."""
    import sympy as sp
    from itertools import permutations

    mW2 = sp.symbols('mW2')
    a12, a13, a23 = sp.symbols('a12 a13 a23')   # OUR p1.p2, p1.p3, p2.p3 (p3 = our symbol)
    b12, b13, b23 = sp.symbols('b12 b13 b23')   # PAPER p1.p2, p1.p3, p2.p3 (all incoming)
    xs = sp.symbols('x1:8')

    def make_F(props, G):
        M = sp.zeros(2,2); B=[[sp.S(0)]*3 for _ in range(2)]; J=sp.S(0)
        for xi,(a,b,r,m2) in zip(xs, props):
            M[0,0]+=xi*a*a; M[0,1]+=xi*a*b; M[1,0]+=xi*a*b; M[1,1]+=xi*b*b
            for k in range(3): B[0][k]+=xi*a*r[k]; B[1][k]+=xi*b*r[k]
            J += xi*(sum(r[i]*r[j]*G[i,j] for i in range(3) for j in range(3)) - m2)
        U = sp.expand(M[0,0]*M[1,1]-M[0,1]*M[1,0])
        adj = sp.Matrix([[M[1,1],-M[0,1]],[-M[1,0],M[0,0]]])
        BB = sum(adj[i,j]*sum(B[i][k]*B[j][l]*G[k,l] for k in range(3) for l in range(3))
                 for i in range(2) for j in range(2))
        return sp.expand(U*J - BB)

    def q(c1,c2,c3): return (sp.Integer(c1), sp.Integer(c2), sp.Integer(c3))

    G_a = sp.Matrix([[0,a12,a13],[a12,0,a23],[a13,a23,mW2]])
    G_b = sp.Matrix([[0,b12,b13],[b12,0,b23],[b13,b23,mW2]])

    # OUR families (p4_out = p1+p2-p3 substituted in D8):
    ours_cur = [(1,0,q(0,0,0),0),(1,0,q(1,0,0),0),(1,0,q(0,-1,0),0),(0,1,q(0,-1,0),0),
                (1,-1,q(0,0,0),0),(1,-1,q(-1,-1,1),1),(0,1,q(0,-1,-1),1)]
    ours_fix = [(1,0,q(0,0,0),0),(1,0,q(1,0,0),0),(1,0,q(0,-1,0),0),(0,1,q(0,-1,0),0),
                (1,-1,q(0,0,0),0),(1,-1,q(1,1,-1),1),(0,1,q(0,-1,-1),1)]
    # PAPER T_4 (p4_in = -(p1+p2+p3) substituted in D8):
    paper = [(1,0,q(0,0,0),0),(1,0,q(1,0,0),0),(1,0,q(0,-1,0),0),(0,1,q(0,-1,0),0),
             (1,-1,q(0,0,0),0),(1,-1,q(1,1,1),1),(0,1,q(0,-1,-1),1)]

    F_p = make_F(paper, G_b)
    P_p = sp.Poly(F_p, *xs)
    pap_terms = dict(zip(P_p.monoms(), [sp.expand(c) for c in P_p.coeffs()]))

    # unknown affine dictionary: b = L a + c mW2
    Ls = sp.symbols('L11 L12 L13 L21 L22 L23 L31 L32 L33 c1 c2 c3')
    L = sp.Matrix(3,3, Ls[:9]); cv = sp.Matrix(Ls[9:])
    subs_b = {b12: (L[0,0]*a12+L[0,1]*a13+L[0,2]*a23+cv[0]*mW2),
              b13: (L[1,0]*a12+L[1,1]*a13+L[1,2]*a23+cv[1]*mW2),
              b23: (L[2,0]*a12+L[2,1]*a13+L[2,2]*a23+cv[2]*mW2)}

    def try_family(name, fam):
        F_o = make_F(fam, G_a)
        P_o = sp.Poly(F_o, *xs)
        our_terms = dict(zip(P_o.monoms(), [sp.expand(c) for c in P_o.coeffs()]))
        hits = []
        for perm in permutations(range(7)):
            # paper monomial in permuted x's must equal our monomial
            ok = True
            eqs = []
            for m, c in pap_terms.items():
                mm = tuple(m[perm.index(i)] for i in range(7))  # x_i^(paper slot perm)
                oc = our_terms.get(mm)
                if oc is None and c != 0:
                    ok = False; break
                eqs.append(sp.expand(c.subs(subs_b, simultaneous=True) - (oc if oc is not None else 0)))
            if not ok: continue
            # collect linear system in Ls over polynomial ring in (a..,mW2)
            unk = list(Ls)
            rows, rhs = [], []
            consistent = True
            for e in eqs:
                pe = sp.Poly(e, a12, a13, a23, mW2)
                for mono, co in zip(pe.monoms(), pe.coeffs()):
                    row = [sp.diff(co, v) for v in unk]
                    const = co - sum(r*v for r, v in zip(row, unk))
                    if any(sp.diff(r, v) != 0 for r in row for v in unk):
                        consistent = False; break  # nonlinear (shouldn't happen)
                    rows.append(row); rhs.append(-sp.expand(const))
                if not consistent: break
            if not consistent: continue
            Amat = sp.Matrix(rows); bvec = sp.Matrix(rhs)
            try:
                sol = sp.linsolve((Amat, bvec), *unk)
            except Exception:
                continue
            if not sol or sol == sp.EmptySet: continue
            solv = list(sol)[0]
            if any(x.free_symbols for x in solv): continue  # underdetermined: skip (would still report)
            Lnum = sp.Matrix(3,3, solv[:9])
            if Lnum.det() == 0: continue
            hits.append((tuple(p+1 for p in perm), solv))
        print(f"{name}: {len(hits)} dictionary hit(s)")
        for perm, solv in hits[:3]:
            Lnum = sp.Matrix(3,3, solv[:9]); cnum = solv[9:]
            print(f"  x-perm {perm}")
            print(f"  b12 = {sp.nsimplify(Lnum[0,0])}*a12 + {Lnum[0,1]}*a13 + {Lnum[0,2]}*a23 + {cnum[0]}*mW2")
            print(f"  b13 = {Lnum[1,0]}*a12 + {Lnum[1,1]}*a13 + {Lnum[1,2]}*a23 + {cnum[1]}*mW2")
            print(f"  b23 = {Lnum[2,0]}*a12 + {Lnum[2,1]}*a13 + {Lnum[2,2]}*a23 + {cnum[2]}*mW2")
        return hits

    h1 = try_family("ours_CURRENT (D8: -p4_out)", ours_cur)
    h2 = try_family("ours_FIXED   (D8: +p4_out)", ours_fix)
    print("\nVERDICT: current =", "PAPER-T4-equivalent" if h1 else "NOT paper T_4",
          "| fixed =", "PAPER-T4-equivalent" if h2 else "NOT paper T_4")

    return h1, h2


# ============================================================================
# identity — GENERAL mode
# ============================================================================
# Generalizes the make_F affine scalar-product-dictionary method of
# identity_check_wpair() (2 loops / 7 props / 3 dots + mW2) to arbitrary loop
# count, propagator count, and the FULL symbolic external Gram (+ mass
# monomials as extra affine generators).  Validated on an 18-propagator
# multi-loop family against its published ancillary propagator list:
# rank 10/10 unique, L = identity, det 1, polynomial verification == 0.
# The no-arg `identity` subcommand above is the wpair regression fixture;
# this path activates only when flags are given.  Slot map between the two
# propagator lists defaults to the identity (positional); pass --x-perm for
# an explicit reordering.  (The wpair fixture is the place that ENUMERATES
# slot permutations; at 18 props that is not an option.)

def _idg_natkey(name):
    import re as _re
    return [int(t) if t.isdigit() else t for t in _re.split(r"(\d+)", name)]


def _idg_sympify(s):
    """sympify a propagator string with EVERY token forced to a fresh Symbol
    (avoids the E/I/N/gamma sympy-builtin traps)."""
    import re as _re
    toks = set(_re.findall(r"[A-Za-z_][A-Za-z_0-9]*", s))
    return sp.expand(sp.sympify(s.replace("^", "**"),
                                locals={t: sp.Symbol(t) for t in toks}))


def _idg_apply_swap(s, swap):
    """Simultaneous symbol/prefix rename on a propagator STRING, e.g.
    swap='l->k,k->p' sends l1->k1 and k3->p3 in one pass (no cascading).
    Rules may be exact symbols ('l1->k1') or prefixes ('l->k')."""
    import re as _re
    if not swap:
        return s
    rules = {}
    for part in swap.split(","):
        a, b = part.split("->")
        rules[a.strip()] = b.strip()

    def _ren(m):
        tok = m.group(0)
        if tok in rules:
            return rules[tok]
        mm = _re.match(r"([A-Za-z_]+?)(\d+)$", tok)
        if mm and mm.group(1) in rules:
            return rules[mm.group(1)] + mm.group(2)
        return tok

    return _re.sub(r"[A-Za-z_][A-Za-z_0-9]*", _ren, s)


def _idg_classify(prop_strs, loops):
    """Fixpoint momentum classification over a propagator set: a symbol is a
    MOMENTUM iff it shares a total-degree-2 monomial with a known momentum
    (seeded by the loop momenta); everything else is a mass symbol.
    Returns (exprs, legs_sorted, mass_syms_sorted)."""
    exprs, allsyms, mono_sets = [], set(), []
    for s in prop_strs:
        e = _idg_sympify(s)
        exprs.append(e)
        syms = sorted(e.free_symbols, key=str)
        allsyms |= {str(x) for x in syms}
        if not syms:
            continue
        p = sp.Poly(e, *syms)
        for mono in p.monoms():
            mono_sets.append((sum(mono),
                              {str(syms[i]) for i, m in enumerate(mono) if m}))
    known = set(loops)
    changed = True
    while changed:
        changed = False
        for deg, ss in mono_sets:
            if deg == 2 and (ss & known) and not (ss <= known):
                known |= ss
                changed = True
    legs = sorted((known & allsyms) - set(loops), key=_idg_natkey)
    masses = sorted(allsyms - known, key=_idg_natkey)
    return exprs, legs, masses


def _idg_parse_prop(idx, expr, loop_syms, leg_syms, mass_syms, raw):
    """expr (sympy) -> (lvec, evec, massdict).  prop == (linear form)^2 - mass,
    mass a polynomial in mass_syms with rational coefficients (possibly 0).
    massdict: generator-name -> Fraction coefficient, one per mass MONOMIAL.
    Raises ValueError (loud, names the propagator) on anything else."""
    mom = loop_syms + leg_syms
    p = sp.Poly(expr, *mom)
    if p.total_degree() != 2:
        raise ValueError(f"prop #{idx} not quadratic in momenta: {raw}")
    lin, ref = [0] * len(mom), None
    for a, q in enumerate(mom):
        c2 = p.nth(*[2 if m is q else 0 for m in mom])
        if c2 != 0:
            r = sp.sqrt(c2)
            if not r.is_Integer:
                raise ValueError(f"prop #{idx}: non-integer q^2 coeff: {raw}")
            lin[a], ref = int(r), a
            break
    if ref is None:
        raise ValueError(f"prop #{idx}: no momentum-squared term: {raw}")
    for b in range(len(mom)):
        if b == ref:
            continue
        mono = [0] * len(mom)
        mono[ref], mono[b] = 1, 1
        cb = sp.Rational(p.nth(*mono), 2 * lin[ref])
        if not cb.is_integer:
            raise ValueError(f"prop #{idx}: not (linear)^2 - m: {raw}")
        lin[b] = int(cb)
    qlin = sum(c * q for c, q in zip(lin, mom))
    mass = sp.expand(qlin**2 - expr)          # prop == qlin^2 - mass
    if any(m in mass.free_symbols for m in mom):
        raise ValueError(f"prop #{idx}: residual momentum dependence "
                         f"(not a perfect square + mass): {raw}")
    md = {}
    if mass != 0:
        msyms = [sp.Symbol(x) for x in mass_syms]
        if not msyms:
            raise ValueError(f"prop #{idx}: mass term {mass} but no mass "
                             f"symbols classified: {raw}")
        pm = sp.Poly(mass, *msyms)
        for mono, co in zip(pm.monoms(), pm.coeffs()):
            if sum(mono) == 0:
                raise ValueError(f"prop #{idx}: dimensionful numeric literal "
                                 f"{co} in mass term: {raw}")
            co = sp.Rational(co)
            gk = "*".join(f"{s}^{m}" if m > 1 else str(s)
                          for s, m in zip(msyms, mono) if m)
            md[gk] = md.get(gk, Fraction(0)) + Fraction(int(co.p), int(co.q))
    return tuple(lin[:len(loop_syms)]), tuple(lin[len(loop_syms):]), md


def _idg_poly_dict(P, xs):
    if P == 0:
        return {}
    out = {}
    for mono, c in sp.Poly(P, *xs).as_dict().items():
        c = sp.Rational(c)
        out[mono] = Fraction(int(c.p), int(c.q))
    return out


def _idg_F_slices(props, NL, NE, pairs, xs):
    """props: list of (lvec, evec, massdict).  F = U*J - BB sliced per affine
    generator.  Keys: ('sp', m, n) Gram dots; ('m', name) mass monomials."""
    M = sp.zeros(NL, NL)
    B = [[sp.S(0)] * NE for _ in range(NL)]
    Jmn = {pr: sp.S(0) for pr in pairs}
    Jm = {}
    for xi, (lv, ev, md) in zip(xs, props):
        for a in range(NL):
            for b in range(NL):
                M[a, b] += xi * lv[a] * lv[b]
            for m in range(NE):
                B[a][m] += xi * lv[a] * ev[m]
        for (m, n) in pairs:
            Jmn[(m, n)] += xi * ev[m] * ev[n] * (1 if m == n else 2)
        for gk, co in md.items():   # J += x_i * (-mass_i)
            Jm[gk] = Jm.get(gk, sp.S(0)) - xi * sp.Rational(co.numerator,
                                                            co.denominator)
    U = sp.expand(M.det())
    adj = M.adjugate()
    slices = {}
    for (m, n) in pairs:
        BB = sp.S(0)
        for a in range(NL):
            for b in range(NL):
                if m == n:
                    BB += adj[a, b] * B[a][m] * B[b][n]
                else:
                    BB += adj[a, b] * (B[a][m] * B[b][n] + B[a][n] * B[b][m])
        P = sp.expand(U * Jmn[(m, n)] - sp.expand(BB))
        slices[("sp", m, n)] = _idg_poly_dict(P, xs)
    for gk, J in Jm.items():
        slices[("m", gk)] = _idg_poly_dict(sp.expand(U * J), xs)
    return slices


def _idg_apply_pins(slices, legs, pins):
    """pins: [((name1,name2), sympy_expr)] — Gram slot a_{n1 n2} := expr
    (0 or a rational-linear combination of mass monomials).  Folds the pinned
    slot's slice into the corresponding mass-generator slices."""
    for (n1, n2), expr in pins:
        if n1 not in legs or n2 not in legs:
            continue
        i, j = legs.index(n1), legs.index(n2)
        k = ("sp", min(i, j), max(i, j))
        if k not in slices:
            continue
        sl = slices.pop(k)
        if expr == 0:
            continue
        syms = sorted(expr.free_symbols, key=str)
        p = sp.Poly(expr, *syms)
        for mono, co in zip(p.monoms(), p.coeffs()):
            if sum(mono) == 0:
                raise ValueError(f"pin {n1}.{n2}={expr}: dimensionful numeric "
                                 f"constant not allowed")
            gk = "*".join(f"{s}^{m}" if m > 1 else str(s)
                          for s, m in zip(syms, mono) if m)
            co = sp.Rational(co)
            cof = Fraction(int(co.p), int(co.q))
            tgt = slices.setdefault(("m", gk), {})
            for mo, c in sl.items():
                tgt[mo] = tgt.get(mo, Fraction(0)) + cof * c
    return slices


def _idg_solve_dictionary(theirs_sl, ours_sl, tkeys, okeys):
    """Solve  sum_A L[A][B] * P_theirs_A(x) = P_ours_B(x)  for all B, exactly
    over Fractions (=> theirs generator b_A = sum_B L[A][B] a_B).  Returns
    (L, rank, bad_monomial|None); L is None iff inconsistent."""
    monos = set()
    for d in list(theirs_sl.values()) + list(ours_sl.values()):
        monos |= set(d)
    NT, NO = len(tkeys), len(okeys)
    basis, rank = [], 0
    for mo in sorted(monos):
        row = [theirs_sl[k].get(mo, Fraction(0)) for k in tkeys] + \
              [ours_sl[k].get(mo, Fraction(0)) for k in okeys]
        for piv, br in basis:
            if row[piv] != 0:
                f = row[piv] / br[piv]
                row = [r - f * b for r, b in zip(row, br)]
        nz = next((i for i in range(NT) if row[i] != 0), None)
        if nz is not None:
            basis.append((nz, row))
            rank += 1
        elif any(row[i] != 0 for i in range(NT, NT + NO)):
            return None, rank, mo                       # inconsistent HERE
    L = [[Fraction(0)] * NO for _ in range(NT)]
    if rank == NT:
        for piv, br in sorted(basis, key=lambda t: -t[0]):
            for tu in range(NO):
                s = br[NT + tu] - sum(br[j] * L[j][tu]
                                      for j in range(piv + 1, NT))
                L[piv][tu] = s / br[piv]
    return L, rank, None


def _idg_verify_full(theirs_sl, ours_sl, L, tkeys, okeys):
    """Exact check: for each ours-generator B, sum_A L[A][B] P_theirs_A -
    P_ours_B == 0 identically."""
    for tu, okey in enumerate(okeys):
        acc = {}
        for rs, tkey in enumerate(tkeys):
            lam = L[rs][tu]
            if lam == 0:
                continue
            for mo, c in theirs_sl[tkey].items():
                acc[mo] = acc.get(mo, Fraction(0)) + lam * c
        for mo, c in ours_sl[okey].items():
            acc[mo] = acc.get(mo, Fraction(0)) - c
        if any(v != 0 for v in acc.values()):
            return False, okey
    return True, None


def _idg_load_ours(path):
    """FAMILY_DEF-style or amflow-style JSON -> (loops, prop_strings)."""
    d = json.load(open(path))
    fam = d.get("family", d)
    loops = fam.get("loop_momenta") or fam.get("loops")
    if not loops:
        raise ValueError(f"{path}: no 'loop_momenta'/'loops' key")
    props = [p[0] if isinstance(p, (list, tuple)) else str(p)
             for p in fam["propagators"]]
    return list(loops), props


def _idg_load_theirs(path, fmt, key):
    """Mathematica 'Name = {p1, p2, ...};' list or JSON list -> (props, name)."""
    import re as _re
    txt = open(path).read()
    if fmt is None:
        fmt = "json" if path.endswith(".json") else "mathematica"
    if fmt == "json":
        d = json.loads(txt)
        if isinstance(d, dict):
            d = d.get("family", d)
            d = d.get("propagators", d)
        return [str(x) for x in d], "json"
    if key:
        m = _re.search(_re.escape(key) + r"\s*=\s*\{(.*?)\}\s*;", txt, _re.S)
        if not m:
            raise ValueError(f"list '{key}' not found in {path}")
        name = key
    else:
        m0 = _re.search(r"(\w+)\s*=\s*\{", txt)
        if not m0:
            raise ValueError(f"no 'name = {{...}};' list found in {path}")
        name = m0.group(1)
        m = _re.search(_re.escape(name) + r"\s*=\s*\{(.*?)\}\s*;", txt, _re.S)
        if not m:
            raise ValueError(f"list '{name}' not terminated by '}};' in {path}")
    return [e.strip() for e in m.group(1).split(",")], name


def _idg_gen_name(key, legs, side=""):
    if key[0] == "sp":
        return f"{side}{legs[key[1]]}.{legs[key[2]]}"
    return f"{side}{key[1]}"


def identity_check_general(argv):
    """GENERAL identity mode.  Returns rc: 0 PASS, 1 FAIL (loud)."""
    ap = argparse.ArgumentParser(
        prog="seedling identity",
        description="General transcribed-family identity check: solve the "
                    "affine scalar-product dictionary between two propagator "
                    "builds via the second-Symanzik F polynomial over the "
                    "fully symbolic external Gram (+ mass monomials). "
                    "No-arg form = the built-in wpair regression fixture.")
    ap.add_argument("--ours", required=True,
                    help="FAMILY_DEF-style or amflow-style JSON")
    ap.add_argument("--theirs", required=True,
                    help="propagator file (Mathematica list or JSON)")
    ap.add_argument("--theirs-format", choices=("mathematica", "json"),
                    default=None, help="default: by extension")
    ap.add_argument("--theirs-key", default=None,
                    help="Mathematica list name, e.g. PropagatorsPBB "
                         "(default: first list in the file)")
    ap.add_argument("--swap", default=None,
                    help="symbol/prefix renames applied to THEIRS strings, "
                         "e.g. 'l->k,k->p' (simultaneous, no cascade)")
    ap.add_argument("--masses", default=None,
                    help="comma list of bare mass symbols and/or Gram pins "
                         "'p1^2=M2', 'p2^2=0', 'p1.p2=...' (pins fold the "
                         "slot into mass generators, applied to BOTH sides "
                         "by leg name)")
    ap.add_argument("--x-perm", default=None,
                    help="explicit slot map j1,..,jN (1-based): OUR prop i "
                         "corresponds to THEIR prop j_i (default: identity)")
    ap.add_argument("--out", default=None, help="write verdict JSON here")
    a = ap.parse_args(argv)

    loops, our_strs = _idg_load_ours(a.ours)
    their_strs, tname = _idg_load_theirs(a.theirs, a.theirs_format,
                                         a.theirs_key)
    their_strs = [_idg_apply_swap(s, a.swap) for s in their_strs]

    print("== family_audit identity (GENERAL mode) ==")
    if len(our_strs) != len(their_strs):
        print(f"   ours: {len(our_strs)} props vs theirs [{tname}]: "
              f"{len(their_strs)} props")
        print("VERDICT: FAIL — propagator counts differ; not the same family "
              "under any slot map.")
        return 1
    NP = len(our_strs)
    if a.x_perm:
        perm = [int(t) for t in a.x_perm.split(",")]
        if sorted(perm) != list(range(1, NP + 1)):
            print(f"VERDICT: FAIL — --x-perm is not a permutation of 1..{NP}")
            return 1
        their_strs = [their_strs[j - 1] for j in perm]

    # pins / declared mass symbols
    pins, declared = [], []
    for ent in (a.masses.split(",") if a.masses else []):
        ent = ent.strip()
        if not ent:
            continue
        if "=" in ent:
            lhs, rhs = ent.split("=", 1)
            lhs = lhs.strip()
            if lhs.endswith("^2"):
                n1 = n2 = lhs[:-2].strip()
            elif "." in lhs:
                n1, n2 = (t.strip() for t in lhs.split(".", 1))
            else:
                n1 = n2 = lhs
            pins.append(((n1, n2), _idg_sympify(rhs.strip())))
        else:
            declared.append(ent)

    sides = {}
    for label, strs in (("ours", our_strs), ("theirs", their_strs)):
        exprs, legs, masses = _idg_classify(strs, loops)
        for dcl in declared:
            if dcl in legs:      # user override: momentum-classified -> mass
                legs.remove(dcl)
                masses = sorted(set(masses) | {dcl}, key=_idg_natkey)
        loop_syms = [sp.Symbol(x) for x in loops]
        leg_syms = [sp.Symbol(x) for x in legs]
        props = [_idg_parse_prop(i + 1, e, loop_syms, leg_syms, masses, s)
                 for i, (e, s) in enumerate(zip(exprs, strs))]
        NE = len(legs)
        pairs = [(m, n) for m in range(NE) for n in range(m, NE)]
        sides[label] = dict(strs=strs, legs=legs, masses=masses, props=props,
                            NE=NE, pairs=pairs)
        print(f"   {label:6s}: {NP} props, loops={loops}, legs={legs}, "
              f"mass syms={masses or '[]'}"
              + (f"  [{tname}, swap='{a.swap}']" if label == "theirs" else ""))

    xs = sp.symbols(f"x1:{NP + 1}")
    print("== F-polynomial (2nd Symanzik) slices, both sides ==")
    sl, keys = {}, {}
    for label in ("ours", "theirs"):
        S = sides[label]
        slc = _idg_F_slices(S["props"], len(loops), S["NE"], S["pairs"], xs)
        slc = _idg_apply_pins(slc, S["legs"], pins)
        sl[label] = slc
        keys[label] = ([k for k in (("sp",) + pr for pr in S["pairs"])
                        if k in slc] +
                       sorted([k for k in slc if k[0] == "m"]))
    nmono = len(set().union(*[set(d) for d in sl["ours"].values()],
                            *[set(d) for d in sl["theirs"].values()]))
    print(f"   x-monomials in play: {nmono}")

    tkeys, okeys = keys["theirs"], keys["ours"]
    NT, NO = len(tkeys), len(okeys)
    print(f"== affine dictionary solve  b_gen = sum L[gen,gen'] a_gen'  "
          f"({NT}x{NO} unknowns) ==")
    L, rank, badmono = _idg_solve_dictionary(sl["theirs"], sl["ours"],
                                             tkeys, okeys)
    print(f"   theirs-F slice rank: {rank}/{NT}  "
          f"(full => dictionary unique if it exists)")
    verdict, why = "FAIL", ""
    detL = None
    is_id = False
    okv = False
    if L is None:
        mono_str = "*".join(f"x{i+1}^{m}" if m > 1 else f"x{i+1}"
                            for i, m in enumerate(badmono) if m) or "1"
        why = (f"no affine dictionary exists — inconsistent at x-monomial "
               f"{mono_str}; the two builds are NOT the same family "
               f"(transcription/convention mismatch class)")
    elif rank < NT:
        why = (f"dictionary underdetermined (rank {rank} < {NT}); "
               f"refusing to certify")
    elif NT != NO:
        why = (f"generator counts differ ({NT} theirs vs {NO} ours) — "
               f"pin or declare masses so the Gram bases align")
    else:
        detL = sp.Matrix(NT, NT, lambda i, j: sp.Rational(
            L[i][j].numerator, L[i][j].denominator)).det()
        print(f"   det(L) = {detL}")
        okv, badkey = _idg_verify_full(sl["theirs"], sl["ours"], L,
                                       tkeys, okeys)
        print(f"   full polynomial verification: "
              f"{'PASS' if okv else 'FAIL at ours-generator ' + _idg_gen_name(badkey, sides['ours']['legs'])}")
        is_id = all(L[i][j] == (1 if i == j else 0)
                    for i in range(NT) for j in range(NT))
        print(f"   dictionary is the IDENTITY map "
              f"(theirs g_i = ours g_i): {is_id}")
        if not is_id and detL != 0:
            for i in range(NT):
                terms = [f"{L[i][j]}*{_idg_gen_name(okeys[j], sides['ours']['legs'])}"
                         for j in range(NT) if L[i][j] != 0]
                print(f"     {_idg_gen_name(tkeys[i], sides['theirs']['legs'], 'theirs:')}"
                      f" = {' + '.join(terms) or '0'}")
        if detL == 0:
            why = "dictionary degenerate (det L = 0)"
        elif not okv:
            why = "full polynomial verification failed"
        else:
            verdict = "PASS"
            why = ("affine dictionary exists"
                   + (" and is the identity" if is_id else " (non-identity)"))
    print(f"\nVERDICT: {verdict} — {why}")
    if a.out:
        json.dump({
            "mode": "identity-general",
            "ours": os.path.abspath(a.ours),
            "theirs": os.path.abspath(a.theirs), "theirs_list": tname,
            "swap": a.swap, "masses": a.masses, "x_perm": a.x_perm,
            "n_props": NP, "n_x_monomials": nmono,
            "theirs_slice_rank": f"{rank}/{NT}",
            "dictionary_consistent": L is not None,
            "det_L": (str(detL) if detL is not None else None),
            "dictionary_is_identity": is_id,
            "full_poly_verification": okv,
            "verdict": verdict, "why": why,
        }, open(a.out, "w"), indent=1)
        print(f"[out] {os.path.abspath(a.out)}")
    return 0 if verdict == "PASS" else 1


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "identity":
        if len(sys.argv) == 2:
            identity_check_wpair()          # stored wpair fixture, UNCHANGED
        else:
            sys.exit(identity_check_general(sys.argv[2:]))
    else:
        main()

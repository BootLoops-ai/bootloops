#!/usr/bin/env python3
r"""
deriv_m2.py — OFF-SHELL-LEG mass-derivative operator (d/dm2 -> shifted-G rows).

CANONICAL PATH: tools/maxcut/deriv_m2.py (the flat tools/deriv_m2.py is an
attribute-complete shim, so `import deriv_m2` keeps working; the --out
default writes to CWD).
The engine is family-generic: feed any Family(props, sp_rules, invariants,
mass_invariant, mass_leg); PBB1m below is the bundled reference instantiation
+ executable audit, not a hard-coded assumption.

WHAT IT DOES
  Given a family definition (props as momentum strings, kinematic replacement
  table), builds the first-order momentum operator
        D = sum_{i,j} c_{ij} p_i . d/dp_j        (Gauss-relation chain rule)
  such that on the on-shell surface D == d/d(m2) at FIXED Mandelstams and fixed
  on-shell conditions, i.e.
        D[p1^2] = 1,  D[p_i^2] = 0 (i=2..5 INCLUDING the conserved leg),
        D[s12] = D[s23] = D[s34] = D[s45] = D[s15] = 0.
  Then expresses the integrand-level derivative of any G(nu) as a FINITE SUM of
  shifted G's:  d/dm2 G(nu) = sum (-nu_a) G(nu+e_a) * D[D_a],  with D[D_a]
  re-expanded over the family's own 18 propagators/ISPs (EVERY p1-carrying
  slot is walked, including ISP numerators — D2,D3,D4,D6,D8,D9,N1..N6 for PBB1m).

MANDATORY SELF-CHECKS (executed, not comments):
  * assert D[p_i^2] = delta_{i1} for EVERY external leg, INCLUDING the conserved
    one (operators that preserve only p1^2..p3^2 can still violate the
    conserved-leg condition — a hidden contamination class).
  * all-slots walk audit: symbolic p-dependence of ALL 18 slots is computed
    and compared against the operator's chain-rule coverage; for pbb1m the
    p1-carrying set must be exactly {D2,D3,D4,D6,D8,D9,N1..N6}.
  Smoking-gun class: prop-shift-only recipes give A_diag=(d-3)/u where the
  truth is -eps/u — the chain-rule terms are mandatory.

Usage:
  python3 -m maxcut.deriv_m2 [--audit] # build operator for PBB1m at X0 + audits
  (the audit ALWAYS runs; --audit is accepted for compatibility and not consulted)
  (the flat shim tools/deriv_m2.py forwards script mode via runpy)
  (as module)  fam = pbb1m_family(); op = MassDerivOperator(fam);
               terms = op.deriv_terms(nu)   # dict shift-tuple -> sympy coeff

Suggested limits: `ulimit -v 32505856`, nice -5.
"""
import argparse
import json
import os
import re
import sys
from itertools import combinations_with_replacement

import sympy as sp

# famdef default: package-local FAMILY_DEF_pbb1m.json or DERIV_M2_FAMDEF
# (the reference family-definition JSON is not shipped in this repo — point
# DERIV_M2_FAMDEF at your own).


# ---------------------------------------------------------------------------
# momentum parsing:  "(k1-p1-p2)^2" -> {'k1':1,'p1':-1,'p2':-1}
# ---------------------------------------------------------------------------
_MOM_TOK = re.compile(r'([+-]?)\s*([A-Za-z]\w*)')


def parse_mom(s):
    s = s.strip()
    if s.endswith('^2'):
        s = s[:-2]
    s = s.strip()
    if s.startswith('(') and s.endswith(')'):
        s = s[1:-1]
    out = {}
    pos = 0
    for m in _MOM_TOK.finditer(s):
        sign = -1 if m.group(1) == '-' else 1
        name = m.group(2)
        out[name] = out.get(name, 0) + sign
        pos = m.end()
    # sanity: nothing but tokens/signs/space
    rest = _MOM_TOK.sub('', s).replace(' ', '')
    assert rest == '', f"unparsed momentum junk in {s!r}: {rest!r}"
    return {k: v for k, v in out.items() if v != 0}


class Family:
    """Loop-integral family: props (squared momenta, all massless internal),
    independent external momenta + one conserved leg, invariant replacement
    table sp_rules[(pi,pj)] -> sympy expr in the invariants."""

    def __init__(self, name, loops, ext_indep, conserved, props, sp_rules,
                 invariants, mass_invariant, mass_leg):
        self.name = name
        self.loops = loops                    # e.g. ['k1','k2','k3']
        self.ext = ext_indep                  # e.g. ['p1','p2','p3','p4']
        self.conserved = conserved            # ('p5', {'p1':-1,...}) or None
        self.props_str = props
        self.props = [parse_mom(p) for p in props]
        self.n = len(props)
        self.sp_rules = sp_rules              # {(a,b): expr} a<=b in ext
        self.invariants = invariants          # {name: combo dict over ext}
        self.mass_invariant = mass_invariant  # 'm2'
        self.mass_leg = mass_leg              # 'p1'
        # eliminate conserved leg from props if present
        if conserved is not None:
            cname, crepl = conserved
            for pr in self.props:
                if cname in pr:
                    c = pr.pop(cname)
                    for b, v in crepl.items():
                        pr[b] = pr.get(b, 0) + c * v
        # bilinear symbol table
        self.T = {}
        for a, b in combinations_with_replacement(range(len(loops)), 2):
            self.T[(a, b)] = sp.Symbol(f'T_{loops[a]}_{loops[b]}')
        self.W = {}
        for l in range(len(loops)):
            for i in range(len(self.ext)):
                self.W[(l, i)] = sp.Symbol(f'W_{loops[l]}_{self.ext[i]}')
        self._build_basis()

    def v(self, i, j):
        """p_i . p_j (0-based indices into self.ext), from sp_rules."""
        a, b = self.ext[i], self.ext[j]
        key = (a, b) if (a, b) in self.sp_rules else (b, a)
        return sp.sympify(self.sp_rules[key])

    def expand_sq(self, mom):
        """(sum c_x x)^2 -> expr in T,W symbols + invariant constants."""
        loops, ext = self.loops, self.ext
        lc = [mom.get(k, 0) for k in loops]
        pc = [mom.get(p, 0) for p in ext]
        expr = sp.Integer(0)
        nl, ne = len(loops), len(ext)
        for a in range(nl):
            for b in range(a, nl):
                co = lc[a] * lc[b] * (1 if a == b else 2)
                if co:
                    expr += co * self.T[(a, b)]
        for l in range(nl):
            for i in range(ne):
                co = 2 * lc[l] * pc[i]
                if co:
                    expr += co * self.W[(l, i)]
        for i in range(ne):
            for j in range(i, ne):
                co = pc[i] * pc[j] * (1 if i == j else 2)
                if co:
                    expr += co * self.v(i, j)
        return expr

    def _build_basis(self):
        """Invert props -> bilinears: each T,W symbol as affine expr in the
        prop 'x_a' symbols. Requires the family to span the full bilinear
        space (rank n_bilinear)."""
        ysyms = list(self.T.values()) + list(self.W.values())
        xsyms = [sp.Symbol(f'x{a+1}') for a in range(self.n)]
        self.xsyms = xsyms
        eqs = []
        for a, pr in enumerate(self.props):
            eqs.append(self.expand_sq(pr) - xsyms[a])
        # linear in ysyms; solve
        sol = sp.solve(eqs, ysyms, dict=True)
        assert sol, (f"family {self.name}: prop set does not span the "
                     f"bilinear space (need {len(ysyms)} indep, have {self.n})")
        assert len(sol) == 1
        self.bilinear_of_x = sol[0]   # {T/W symbol: affine expr in xsyms + invs}
        # rank check: all y solved
        assert set(self.bilinear_of_x) == set(ysyms), \
            f"underdetermined bilinear basis for {self.name}"


class MassDerivOperator:
    r"""D = sum_j w_j . d/dp_j,  w_j = sum_i c_ij p_i, with D == d/dm2 on shell.

    js: which d/dp_j are allowed (default: only the massive leg — sufficient
    whenever the reduced condition system is nonsingular; falls back to all)."""

    def __init__(self, fam: Family, js=None, verbose=True):
        self.fam = fam
        self.verbose = verbose
        ne = len(fam.ext)
        mass_j = fam.ext.index(fam.mass_leg)
        tried = []
        for js_try in ([js] if js is not None else [[mass_j], list(range(ne))]):
            c = self._solve_c(js_try)
            tried.append(js_try)
            if c is not None:
                self.js = js_try
                self.c = c   # dict (i,j)->sympy expr
                break
        else:
            raise RuntimeError(f"no operator solution; tried js={tried}")
        self._self_check_legs()
        self._precompute_prop_action()

    # -- linear solve for c_ij ------------------------------------------------
    def _apply_D_to_combo(self, c, q):
        """D[(q.p)^2] with unknown/known c: 2 sum_j q_j sum_i c_ij sum_a q_a v_ia"""
        fam = self.fam
        ne = len(fam.ext)
        expr = sp.Integer(0)
        for (i, j), cij in c.items():
            qj = q.get(fam.ext[j], 0)
            if qj == 0:
                continue
            for a in range(ne):
                qa = q.get(fam.ext[a], 0)
                if qa == 0:
                    continue
                expr += 2 * qj * qa * cij * fam.v(min(i, a), max(i, a))
        return expr

    def _solve_c(self, js):
        fam = self.fam
        ne = len(fam.ext)
        unk = {}
        syms = []
        for j in js:
            for i in range(ne):
                s = sp.Symbol(f'c_{i}_{j}')
                unk[(i, j)] = s
                syms.append(s)
        eqs = []
        for nameQ, q in fam.invariants.items():
            target = sp.Integer(1) if nameQ == fam.mass_invariant else sp.Integer(0)
            eqs.append(sp.together(self._apply_D_to_combo(unk, q) - target))
        sol = sp.solve(eqs, syms, dict=True)
        if not sol:
            return None
        sol = sol[0]
        # any leftover free symbols among unknowns -> set to 0 (valid: they span
        # LI directions annihilating every invariant)
        c = {}
        for key, s in unk.items():
            val = sol.get(s, s)
            val = val.subs({t: 0 for t in syms if t not in sol})
            val = sp.cancel(sp.together(val))
            if val != 0:
                c[key] = val
        # verify exactly (the solve may return a parametric family)
        for nameQ, q in fam.invariants.items():
            target = sp.Integer(1) if nameQ == fam.mass_invariant else sp.Integer(0)
            r = sp.cancel(sp.together(self._apply_D_to_combo(c, q) - target))
            if sp.simplify(r) != 0:
                return None
        return c

    # -- MANDATORY self-check: every leg incl. conserved ---------------------
    def _self_check_legs(self):
        fam = self.fam
        legs = {fam.ext[i]: {fam.ext[i]: 1} for i in range(len(fam.ext))}
        if fam.conserved is not None:
            cname, crepl = fam.conserved
            legs[cname] = dict(crepl)   # p5 as combo of independents
        results = {}
        for lname, q in legs.items():
            val = sp.simplify(sp.cancel(self._apply_D_to_combo(self.c, q)))
            want = 1 if lname == fam.mass_leg else 0
            results[lname] = (val, want)
            assert val == want, (
                f"SELF-CHECK FAILED: D[{lname}^2] = {val}, "
                f"expected {want}. Operator MUST preserve every on-shell "
                f"condition including the conserved leg.")
        if self.verbose:
            print(f"[deriv_m2] self-check PASS: D[p_i^2]=delta_i,mass for "
                  f"{list(results)} (incl. conserved leg)")
        self.leg_check = {k: str(v[0]) for k, v in results.items()}

    # -- v4: walk ALL slots; expand D[D_a] over the prop basis ----------------
    def _precompute_prop_action(self):
        fam = self.fam
        ne = len(fam.ext)
        self.p_dependence = []   # per slot: sorted list of ext names with lam!=0
        self.action = []         # per slot: None or dict {None: const, b: coeff}
        for a, pr in enumerate(fam.props):
            lam = {j: pr.get(fam.ext[j], 0) for j in range(ne)}
            dep = sorted(fam.ext[j] for j, l in lam.items() if l != 0)
            self.p_dependence.append(dep)
            # D[D_a] = 2 sum_j lam_aj sum_i c_ij (p_i . M_a)
            expr = sp.Integer(0)
            for (i, j), cij in self.c.items():
                if lam.get(j, 0) == 0:
                    continue
                # p_i . M_a
                pim = sp.Integer(0)
                for l, kname in enumerate(fam.loops):
                    cl = pr.get(kname, 0)
                    if cl:
                        pim += cl * fam.W[(l, i)]
                for b in range(ne):
                    cb = pr.get(fam.ext[b], 0)
                    if cb:
                        pim += cb * fam.v(min(i, b), max(i, b))
                expr += 2 * lam[j] * cij * pim
            if expr == 0:
                self.action.append(None)
                continue
            expr = sp.expand(expr.subs(fam.bilinear_of_x))
            poly = {}
            const = expr
            for b, xs in enumerate(fam.xsyms):
                co = sp.cancel(sp.together(sp.expand(expr.coeff(xs, 1))))
                if co != 0:
                    poly[b] = co
                const = const.coeff(xs, 0)
            const = sp.cancel(sp.together(const))
            if const != 0:
                poly[None] = const
            self.action.append(poly if poly else None)

    def audit(self, expect_mass_dep=None):
        """v4 executable audit: report per-slot p-dependence + chain-rule
        coverage. expect_mass_dep: set of 1-based slot indices expected to
        carry the massive leg (assert exact match if given)."""
        fam = self.fam
        rows = []
        mass = fam.mass_leg
        walked = set()
        for a in range(fam.n):
            dep = self.p_dependence[a]
            act = self.action[a]
            covered = act is not None
            if covered:
                walked.add(a + 1)
            rows.append({'slot': a + 1, 'prop': fam.props_str[a],
                         'p_dependence': dep, 'chain_rule_active': covered,
                         'n_shift_terms': (len(act) if act else 0)})
        mass_dep = {a + 1 for a in range(fam.n) if mass in self.p_dependence[a]}
        js_ext = {fam.ext[j] for j in self.js}
        # every slot depending on any operator direction must be walked
        for a in range(fam.n):
            if set(self.p_dependence[a]) & js_ext:
                assert self.action[a] is not None, \
                    f"v4 AUDIT FAIL: slot {a+1} depends on {js_ext} but has no action"
        if expect_mass_dep is not None:
            assert mass_dep == set(expect_mass_dep), (
                f"v4 AUDIT FAIL: {mass}-carrying slots {sorted(mass_dep)} != "
                f"expected {sorted(expect_mass_dep)}")
        if self.verbose:
            print(f"[deriv_m2] v4 audit: {len(mass_dep)} slots carry {mass}: "
                  f"{sorted(mass_dep)}; operator directions {sorted(js_ext)}; "
                  f"walked slots with action: {sorted(walked)}")
        return {'per_slot': rows, 'mass_carrying_slots': sorted(mass_dep),
                'operator_directions': sorted(js_ext),
                'walked_slots': sorted(walked),
                'coefficients_c_ij': {f"{fam.ext[i]}.d/d{fam.ext[j]}": str(v)
                                      for (i, j), v in self.c.items()},
                'leg_self_check': self.leg_check}

    # -- the deliverable ------------------------------------------------------
    def deriv_terms(self, nu):
        """d/dm2 G(nu) = sum_shift coeff * G(shift). Returns {tuple: sympy}."""
        fam = self.fam
        nu = tuple(nu)
        assert len(nu) == fam.n
        terms = {}
        for a, na in enumerate(nu):
            if na == 0:
                continue
            poly = self.action[a]
            if poly is None:
                continue
            for key, co in poly.items():
                c = -na * co
                sh = list(nu)
                sh[a] += 1
                if key is not None:
                    sh[key] -= 1
                sh = tuple(sh)
                terms[sh] = sp.cancel(sp.together(terms.get(sh, 0) + c))
        return {k: v for k, v in terms.items() if v != 0}


# ---------------------------------------------------------------------------
# PBB1m instantiation
# ---------------------------------------------------------------------------
X0_SIJ = {'s12': sp.Rational(-22, 5), 's23': sp.Rational(241, 25),
          's34': sp.Rational(-377, 100), 's45': sp.Rational(13, 50),
          's15': sp.Rational(249, 50)}
X0_M2 = sp.Rational(137, 50)   # anchor value; production keeps m2 SYMBOLIC


def pbb1m_family(numeric_x0=True, m2_symbolic=True, famdef_path=None):
    famdef_path = famdef_path or os.environ.get('DERIV_M2_FAMDEF') or os.path.join(
        os.path.dirname(os.path.abspath(__file__)), 'FAMILY_DEF_pbb1m.json')
    fd = json.load(open(famdef_path))
    props = fd['propagators']
    m2, s12, s23, s34, s45, s15 = sp.symbols('m2 s12 s23 s34 s45 s15')
    subs = {}
    if numeric_x0:
        subs = {sp.Symbol(k): v for k, v in X0_SIJ.items()}
        if not m2_symbolic:
            subs[m2] = X0_M2
    # replacement table == the family definition's 'replacements' block
    sp_rules = {
        ('p1', 'p1'): m2, ('p2', 'p2'): 0, ('p3', 'p3'): 0, ('p4', 'p4'): 0,
        ('p1', 'p2'): (s12 - m2) / 2, ('p2', 'p3'): s23 / 2,
        ('p3', 'p4'): s34 / 2, ('p1', 'p3'): (s45 - s12 - s23) / 2,
        ('p1', 'p4'): (s23 - s45 - s15) / 2, ('p2', 'p4'): (s15 - s23 - s34) / 2,
    }
    sp_rules = {k: sp.sympify(v).subs(subs) for k, v in sp_rules.items()}
    invariants = {
        'm2':  {'p1': 1},
        'p2sq': {'p2': 1}, 'p3sq': {'p3': 1}, 'p4sq': {'p4': 1},
        'p5sq': {'p1': 1, 'p2': 1, 'p3': 1, 'p4': 1},
        's12': {'p1': 1, 'p2': 1}, 's23': {'p2': 1, 'p3': 1},
        's34': {'p3': 1, 'p4': 1}, 's45': {'p1': 1, 'p2': 1, 'p3': 1},
        's15': {'p2': 1, 'p3': 1, 'p4': 1},
    }
    fam = Family('pbb1m', ['k1', 'k2', 'k3'], ['p1', 'p2', 'p3', 'p4'],
                 ('p5', {'p1': -1, 'p2': -1, 'p3': -1, 'p4': -1}),
                 props, sp_rules, invariants, 'm2', 'p1')
    # closure sanity: p5^2 == 0 exactly under the table
    p5sq = sum(fam.v(i, j) * (1 if i == j else 2)
               for i in range(4) for j in range(i, 4))
    assert sp.simplify(p5sq) == 0, f"p5^2 != 0 under replacement table: {p5sq}"
    return fam


# expected p1-carrying slots (1-based) for the pbb1m family:
PBB1M_P1_SLOTS = {2, 3, 4, 6, 8, 9, 12, 13, 14, 15, 16, 17}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--audit', action='store_true',
                    help='accepted for compatibility; the audit always runs '
                         '(flag not consulted)')
    ap.add_argument('--symbolic', action='store_true',
                    help='full symbolic s_ij (slower); default X0-numeric')
    ap.add_argument('--out',
                    default=os.path.join(os.getcwd(), 'DERIV_M2_AUDIT.json'))
    a = ap.parse_args()
    fam = pbb1m_family(numeric_x0=not a.symbolic)
    op = MassDerivOperator(fam)
    audit = op.audit(expect_mass_dep=PBB1M_P1_SLOTS)
    # demo: term list for the target corner T
    T = tuple([1] * 11 + [0] * 7)
    terms = op.deriv_terms(T)
    print(f"[deriv_m2] d/dm2 T: {len(terms)} shifted integrals")
    audit['corner_T'] = {'nu': list(T), 'n_terms': len(terms),
                         'terms': {str(list(k)): str(v) for k, v in terms.items()}}
    json.dump(audit, open(a.out, 'w'), indent=1)
    print(f"[deriv_m2] audit written -> {a.out}")


if __name__ == '__main__':
    main()

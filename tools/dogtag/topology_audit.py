#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision.
# =============================================================================
#  topology_audit.py  --  Dogtag: integral-family identity checks; the one-line
#  provenance/topology auditor for Feynman integral families (Kira
#  integralfamilies.yaml  or  AmflowFamily family.jl).  Package tools/dogtag
#  (formerly topology-audit; the script keeps its historical file name).
#
#  WHY THIS EXISTS
#  ---------------
#  A real case: a family file labeled "nonplanar 2L double box" turned out to be the
#  PLANAR Smirnov double box mislabeled.  Under k2 -> -k2 its propagator SET
#  maps exactly onto Smirnov's planar set; amflow gave Im=0 at u>0 (the planar
#  signature; a true crossed box has an open u-cut -> nonzero Im).  Catching
#  that cost a manual investigation.  This tool makes it a one-liner.
#
#  WHAT IT DOES, given a family file:
#    (a) LOOP-RELABELING ISOMORPHISM.  Tries every signed permutation of the
#        loop momenta k_i -> +-k_{sigma(i)} (and, with --leg-perms, external-leg
#        permutations).  A propagator q^2 is invariant under q -> -q, so two
#        families that are the SAME GRAPH under such a relabeling have identical
#        canonical propagator sets.  We compare against a CANONICAL catalog and
#        cross-check pairwise.  The set is LABELED: every line carries its
#        mass and its multiplicity (propagator power), every leg its
#        virtuality class, and a relabeling must carry them along -- a match
#        on the momenta alone is reported as MASS-BLIND, never as identity.
#        Two routings of ONE graph differ by a loop-momentum basis change
#        (k1 -> k2 - k1, k -> k + p_i) that no signed permutation reaches:
#        the AFFINE class k -> U k + c.p (U unimodular, c solved) is searched
#        second and describes the map; the identity test that is routing-
#        independent by construction is the isomorphism of the REALIZED
#        labeled graphs (VF2), reported beside it with a canonical hash of
#        that graph (section 2b).
#    (b) PLANARITY.  Reconstructs the Feynman graph from momentum-flow incidence
#        and runs a graph-theoretic planarity test (networkx.check_planarity).
#    (c) CUT SIGNATURE.  For each Mandelstam channel (s,t,u, as 2-particle
#        external partitions) counts how many propagators a cut must sever
#        (the channel's minimal edge cut separating that external pair from the
#        rest).  The (s,t,u) multiset is a label-independent fingerprint.
#
#  OUTPUT: a fingerprint dict + human-readable WARNINGS, e.g.
#     "graph-isomorphic to <known family> under k2 -> -k2"
#     "labeled 'nonplanar' but the reconstructed graph IS planar"
#
#  Pure sympy + networkx.  No amflow, no Kira, no heavy compute.
#
#  USAGE
#    python3 topology_audit.py FAMILY_FILE [--name NAME] [--leg-perms]
#    python3 topology_audit.py FAMILY_FILE --drawn DRAWN_GRAPH [--name NAME] [--verdict]
#      (--name: the family of that name in a multi-family file is the side
#       compared; an unknown name is refused, rc 2, the file's names listed)
#    python3 topology_audit.py --self-test
#    python3 -m loomcheck --selftest | --basis FILE   (the loomcheck member,
#      loomcheck/: Yangian/loom/fishnet applicability screen for position-
#      space conformal integrals; its battery is the self-test's [loomcheck] leg)
#  FAMILY_FILE may be
#    * a Kira integralfamilies.yaml                       (load_kira_yaml)
#    * an AmflowFamily *.jl                                (load_amflow_jl)
#    * an AMFlow-port JSON config, or its bare family block written by hand
#      as a propagator list in the same keys               (load_amflow_json)
#    * a pySecDec script calling LoopIntegralFromGraph     (load_pysecdec_graph)
#    * a drawn-graph edge list (JSON, a "drawn_graph:" yaml block, or the
#      "# edge list:" / "# legs:" header of a routed yaml) (load_edge_list)
# =============================================================================
from __future__ import annotations
import sys, os, re, json, argparse, itertools
import sympy as sp

try:
    import networkx as nx
    _HAVE_NX = True
except Exception:                                   # pragma: no cover
    _HAVE_NX = False


# ---------------------------------------------------------------------------
#  0.  PROPAGATOR PARSING
# ---------------------------------------------------------------------------
#  Every propagator is  q^2 - mass^2  with q a linear combination of loop
#  momenta k_i and external momenta p_j.  We represent q by its coefficient
#  vector over the basis (k_1..k_L, p_1..p_E) after applying external
#  substitutions (e.g. p4 = -p1-p2-p3).  ISPs that are pure dot products
#  (k_i * p_j) are NOT squared propagators and are dropped from the graph.

class Family:
    def __init__(self, name, loops, exts, ext_subs, propagators, source,
                 physical=None, nu=None, leg_virt=None, mass_values=None,
                 kinematics_source=None, cyclic_leg_order=None,
                 source_kind=None, record_notes=None):
        self.name = name                # str
        self.loops = list(loops)        # ['k1','k2',...]
        self.exts = list(exts)          # ['p1','p2','p3','p4']
        self.ext_subs = dict(ext_subs)  # {'p4':'-p1-p2-p3'}
        self.propagators = propagators  # list of (raw_expr_str, mass_str|0)
        self.source = source            # file path
        # physical[i] = True if propagator i is in the top sector (a real edge);
        # False = ISP / numerator auxiliary.  None entries default to True.
        if physical is None:
            physical = [True] * len(propagators)
        self.physical = list(physical)
        # nu[i] = propagator power of the integral that defines the top sector
        # (the record's index vector: >0 a line with that multiplicity, 0 an
        # ISP, <0 a numerator power).  None = not declared by the record.
        self.nu = list(nu) if nu is not None else None
        # leg_virt = {leg: virtuality string p^2} for the legs the record
        # declares it for ("0" on shell); dependent legs are filled by the
        # kinematics reader.
        self.leg_virt = dict(leg_virt or {})
        # mass_values = {mass symbol: value string} from the record's numeric
        # values (an AMFlow config's numeric_values); propagator strings stay
        # verbatim, resolve_mass() applies them.
        self.mass_values = dict(mass_values or {})
        # kinematics_source: where the external-leg set came from ("json:legs+
        # conservation", "jl:exts+ext_subs", "graph:legs"); None when the leg
        # set was inferred from the propagator symbols (the Kira reader).
        self.kinematics_source = kinematics_source
        # cyclic_leg_order: the legs in the drawn cyclic order (drawn graphs only).
        self.cyclic_leg_order = (list(cyclic_leg_order)
                                 if cyclic_leg_order is not None else None)
        self.source_kind = source_kind  # "kira_yaml" | "amflow_jl" | "amflow_json"
                                        # | "pysecdec_graph" | "edge_list"
        # record_notes: reader annotations (chosen integral index, replacement
        # rules, routing tree/chords, vertex-conservation witness ...).
        self.record_notes = dict(record_notes or {})
        self.kinematics_path = None     # --kinematics file recorded by load_family

    def resolve_mass(self, m):
        """Apply the record's numeric mass values to a mass^2 expression."""
        expr = sp.sympify(m)
        if not self.mass_values:
            return expr
        subs = {sp.Symbol(k): sp.nsimplify(v) for k, v in self.mass_values.items()}
        return sp.nsimplify(expr.subs(subs))

    def __repr__(self):
        return (f"Family({self.name!r}, loops={self.loops}, "
                f"exts={self.exts}, nprop={len(self.propagators)}, "
                f"nphys={sum(self.physical)})")


def _symbols(loops, exts):
    syms = {}
    for s in list(loops) + list(exts):
        syms[s] = sp.Symbol(s)
    return syms


def _resolve_ext_subs(exts, ext_subs, syms):
    """Return {ext_symbol: sympy expr in the *independent* externals}."""
    out = {}
    for e in exts:
        if e in ext_subs:
            out[syms[e]] = sp.sympify(ext_subs[e], locals=syms)
        else:
            out[syms[e]] = syms[e]
    return out


def _parse_momentum(expr_str, loops, exts, syms, ext_resolved, mass_str):
    """
    Parse one propagator into (momentum_vector, is_squared_propagator, mass_sym).

    Handles the two raw conventions seen in this repo:
      * full squared form:   "k1^2", "(k1+p1)^2 - m2"
      * bare-momentum form (kira [expr, mass]):  "k1", "k1 + k2 - p1 - p2"
        meaning  expr^2 - mass.
      * pure dot products (ISP):  "k1*p2", "k1.p2"  -> NOT a propagator.

    momentum_vector: tuple of sympy coefficients over (loops..., exts...),
                     with externals already substituted to the independent set.
    """
    s = expr_str.strip()
    s = s.replace('.', '*')             # k1.p2 -> k1*p2 dot-product spelling
    # Detect a pure dot product (ISP): a single product of two distinct momenta,
    # i.e. it is linear in NO single momentum but bilinear.  Easiest robust test:
    # build the sympy expr; if it is NOT of the form (linear)^2 - const, treat as ISP.
    full = sp.sympify(s, locals=syms)

    # Did the raw string already contain a square?  ("^2" / "**2")
    has_explicit_square = bool(re.search(r'\*\*\s*2|\^\s*2', expr_str))

    if has_explicit_square:
        # form: <linear>^2 [ - mass ].  The squared part is the propagator momentum.
        # Strip any additive non-squared mass term (already separated for kira;
        # for jl strings the "- m2" rides along -> remove it).
        q2 = full
        mass = sp.Integer(0)
        # peel off additive terms that are NOT a perfect square of a linear form
        # (these are the explicit "- m2").  We do it by collecting the Pow(...,2).
        terms = sp.Add.make_args(sp.expand(full)) if full.is_Add else [full]
        # Simpler & robust: re-extract the linear momentum from the ORIGINAL string
        # by grabbing the parenthesised / bare token raised to 2.
        msq_match = re.search(r'(.*?)(?:\)\s*\^?\s*\*{0,2}2|\^\s*2|\*\*\s*2)\s*([-+].*)?$', expr_str.strip())
        # fall back to symbolic route below
        lin = _extract_linear_from_square(expr_str, syms)
        if lin is None:
            # could not isolate; treat whole thing minus mass as ISP-ish
            return None
        mass_extra = _trailing_mass(expr_str, syms)
        mass = mass_extra + _mass_from_field(mass_str, syms)
        vec = _linear_to_vec(lin, loops, exts, syms, ext_resolved)
        if vec is None:
            return None
        return (vec, True, mass)
    else:
        # bare-momentum kira form: expr is the momentum, mass field is mass^2.
        lin = full
        # Is it actually linear in the momenta?  (ISP dot products are not.)
        vec = _linear_to_vec(lin, loops, exts, syms, ext_resolved)
        if vec is None:
            return None                 # nonlinear -> ISP dot product, skip
        mass = _mass_from_field(mass_str, syms)
        return (vec, True, mass)


def _extract_linear_from_square(expr_str, syms):
    """From 'X^2 - m2' or '(X)^2' return sympy linear form X (or None)."""
    s = expr_str.strip()
    # normalise ^2 / **2
    s = s.replace('**2', '^2')
    # match the first  (....)^2  or  token^2
    m = re.match(r'\s*\(([^()]*(?:\([^()]*\)[^()]*)*)\)\s*\^2', s)
    if m:
        inner = m.group(1)
    else:
        m = re.match(r'\s*([A-Za-z]\w*)\s*\^2', s)
        if not m:
            return None
        inner = m.group(1)
    try:
        return sp.sympify(inner, locals=syms)
    except Exception:
        return None


def _trailing_mass(expr_str, syms):
    """Return the additive mass term after the '^2', as sympy (mass^2 quantity)."""
    s = expr_str.strip().replace('**2', '^2')
    m = re.search(r'\^2\s*(.+)$', s)
    if not m:
        return sp.Integer(0)
    tail = m.group(1).strip()
    if not tail:
        return sp.Integer(0)
    try:
        # tail like "- m2" -> the propagator is q^2 - m2, mass^2 = m2 = -tail
        return -sp.sympify(tail, locals=syms)
    except Exception:
        return sp.Integer(0)


def _mass_from_field(mass_str, syms):
    if mass_str in (0, "0", None, ""):
        return sp.Integer(0)
    try:
        return sp.sympify(str(mass_str), locals=syms)
    except Exception:
        return sp.Symbol(str(mass_str))


def _linear_to_vec(lin, loops, exts, syms, ext_resolved):
    """
    Express linear momentum form 'lin' as a coefficient vector over
    (loops..., independent-exts...).  Returns None if 'lin' is not linear
    (i.e. it is an ISP dot product like k1*p2).
    """
    expr = sp.expand(sp.sympify(lin).subs(ext_resolved))
    basis = [syms[k] for k in loops] + [syms[e] for e in exts if e not in
                                        # only independent externals survive subst
                                        []]
    # Determine which external symbols remain after substitution (independent set)
    indep_ext = [e for e in exts if syms[e] in expr.free_symbols
                 or True]   # keep full ext order for a stable basis
    basis = [syms[k] for k in loops] + [syms[e] for e in exts]
    coeffs = []
    poly_ok = True
    for b in basis:
        c = expr.coeff(b, 1)
        # linearity: c must be free of every momentum symbol
        if any(m in c.free_symbols for m in basis):
            poly_ok = False
            break
        coeffs.append(sp.nsimplify(c))
    if not poly_ok:
        return None
    # check there is no constant or quadratic remainder beyond the linear part
    reconstructed = sum(c * b for c, b in zip(coeffs, basis))
    remainder = sp.expand(expr - reconstructed)
    if any(m in remainder.free_symbols for m in basis):
        return None
    return tuple(coeffs)


# ---------------------------------------------------------------------------
#  1.  FILE LOADERS
# ---------------------------------------------------------------------------
def load_kira_yaml(path, kinematics=None):
    """
    Parse a Kira integralfamilies.yaml (minimal hand-parser; no pyyaml dep).

    The external-leg set comes from DATA when the record carries it, in this
    order of precedence:
      1. `kinematics` -- a Kira kinematics.yaml path (read_kinematics_yaml:
         incoming_momenta / outgoing_momenta, the momentum_conservation rule,
         the scalar-product rules -> leg virtualities);
      2. the family block's own external_momenta / momentum_conservation /
         leg_virtualities keys (the routed-drawing form);
    and only when neither exists from the FALLBACK _infer_externals /
    _infer_ext_subs (the propagator symbols), which marks the family
    kinematics_source = "inferred" and records the LEG-SET-INFERRED note that
    audit() prints and acts on (planarity withheld).  A family whose top
    sector never spells its dependent leg is read by the fallback with one leg
    too few and the wrong conservation rule -- that is the trap the data
    routes exist for.  Dependent-leg completion (single-current second leg,
    implicit p_N, composite legs) is complete_legs(), applied by load_family.
    """
    txt = open(path).read()
    fams = []
    kin = read_kinematics_yaml(kinematics) if kinematics else None
    # split into family blocks at '- name:'
    blocks = re.split(r'(?m)^\s*-\s*name\s*:', txt)
    header = blocks[0]
    for blk in blocks[1:]:
        nm = re.match(r'\s*"?([^"\n]+)"?', blk)
        name = nm.group(1).strip().strip('"') if nm else "unknown"
        lm = re.search(r'loop_momenta\s*:\s*\[([^\]]*)\]', blk)
        loops = [x.strip() for x in lm.group(1).split(',')] if lm else []
        # top_level_sectors: [N] -> bitmask of which propagators are physical
        tls = re.search(r'top_level_sectors\s*:\s*\[\s*([0-9]+)', blk)
        top_sector = int(tls.group(1)) if tls else None
        # propagators: list of [ "expr", mass ]
        props = []
        pblock = re.search(r'propagators\s*:(.*?)(?:cut_propagators|top_level_sectors|\Z|\n\s*-\s*name)',
                           blk, re.S)
        seg = pblock.group(1) if pblock else blk
        for m in re.finditer(r'\[\s*"([^"]*)"\s*,\s*([^\]\s]+)\s*\]', seg):
            props.append((m.group(1), m.group(2)))
        # the family block's own leg keys (routed-drawing yamls carry them)
        yaml_legs = None
        me = re.search(r'(?m)^\s*external_momenta\s*:\s*\[([^\]]*)\]', blk)
        if me:
            yaml_legs = [x.strip().strip('"\'') for x in me.group(1).split(",") if x.strip()]
        yaml_subs = {}
        mc = re.search(r'(?m)^\s*momentum_conservation\s*:\s*\{([^}]*)\}', blk)
        if mc:
            for a, b in re.findall(r'(\w+)\s*:\s*"([^"]*)"', mc.group(1)):
                yaml_subs[a] = b.strip()
        yaml_virt = {}
        mv = re.search(r'(?m)^\s*leg_virtualities\s*:\s*\{([^}]*)\}', blk)
        if mv:
            yaml_virt = {a: b for a, b in re.findall(r'(\w+)\s*:\s*"([^"]*)"', mv.group(1))}
        leg_virt = {}
        notes = {}
        if kin is not None:
            exts = list(kin["legs"])
            ext_subs = dict(kin["conservation"])
            leg_virt = dict(kin["virtualities"])
            kinematics_source = "kinematics:" + os.path.basename(str(kinematics))
            notes["leg_set"] = {"source": kinematics_source, "declared_legs": list(kin["incoming"]),
                                "declared_outgoing": list(kin["outgoing"]),
                                "conservation_declared": dict(kin["conservation"]),
                                "kinematics_file": os.path.basename(str(kinematics))}
            mass_values = ({kin["symbol_to_replace_by_one"]: "1"}
                           if kin.get("symbol_to_replace_by_one") else {})
        elif yaml_legs is not None:
            exts = yaml_legs
            ext_subs = yaml_subs
            leg_virt = yaml_virt
            kinematics_source = "yaml:external_momenta+momentum_conservation"
            notes["leg_set"] = {"source": kinematics_source, "declared_legs": list(exts),
                                "conservation_declared": dict(ext_subs)}
            mass_values = {}
        else:
            # FALLBACK: externals are not declared in this Kira yaml; infer them
            # from the propagator symbols and say so (LEG-SET-INFERRED).
            exts = _infer_externals(props, loops)
            ext_subs = _infer_ext_subs(props, exts)
            kinematics_source = "inferred"
            notes["leg_set"] = {
                "source": "inferred",
                "declared_legs": list(exts),
                "conservation_declared": dict(ext_subs),
                "warning": _leg_set_inferred_note(exts, ext_subs)}
            mass_values = {}
        # physical mask from top sector bitmask (bit i set -> propagator i physical)
        if top_sector is not None:
            physical = [bool((top_sector >> i) & 1) for i in range(len(props))]
        else:
            physical = None
        fams.append(Family(name, loops, exts, ext_subs, props, path, physical,
                           leg_virt=leg_virt, mass_values=mass_values,
                           kinematics_source=kinematics_source,
                           source_kind="kira_yaml", record_notes=notes))
    return fams


#  symbols that are masses (never momenta) -- excluded from the external set.
_MASS_RE = re.compile(r'(?:m|M)(?:sq|2|t|b|w|z|h|H|W|Z)?\d*$|^msq$|^mm\d*$')


def _leg_set_inferred_note(exts, ext_subs):
    """The LEG-SET-INFERRED sentence for an inferred external-leg set."""
    rule = ", ".join(f"{k} = {v}" for k, v in ext_subs.items()) or "no conservation rule"
    return ("LEG-SET-INFERRED: the external-leg set %s (%s) was inferred from the "
            "propagator symbols, not read from the record (no --kinematics file, no "
            "kinematics.yaml beside the family file, no external_momenta / "
            "momentum_conservation keys in the yaml); a family whose top sector never "
            "spells its dependent leg is read this way with one leg too few and the "
            "wrong conservation rule, so no planarity verdict is given and the cut "
            "signature / canonical set below are those of the inferred set"
            % (exts, rule))


def _infer_externals(props, loops):
    """
    FALLBACK (LEG-SET-INFERRED), used by load_kira_yaml only when the record
    carries no kinematics data: External momenta = symbols appearing in
    propagators that are not loop momenta and not masses.  Handles p1,p2,...;
    bare p / q / Q (single-current); and non-standard loop names (u1,u2) with
    external Q.  A leg that the top sector never spells (the eliminated p_N of
    a four-point family written in p1..p3) is invisible to this rule; the
    caller marks the family kinematics_source = "inferred" and audit() prints
    the LEG-SET-INFERRED line and withholds the planarity verdict.
    """
    loopset = set(loops)
    syms = set()
    for expr, _ in props:
        for tok in re.findall(r'[A-Za-z]\w*', expr):
            syms.add(tok)
    exts = []
    for s in sorted(syms):
        if s in loopset:
            continue
        if _MASS_RE.match(s):
            continue
        exts.append(s)
    # stable ordering: p1,p2,... first (numeric), then others
    def keyf(e):
        m = re.fullmatch(r'p(\d+)', e)
        return (0, int(m.group(1))) if m else (1, e)
    return sorted(exts, key=keyf)


def _infer_ext_subs(props, exts):
    """FALLBACK (LEG-SET-INFERRED) heuristic, paired with _infer_externals: if
    p1..pn present and momentum conservation p1+..+pn=0 is the convention,
    substitute the highest-index leg.  We only do this when exactly the
    contiguous set p1..pN is present (the all-incoming 4-pt convention).  This
    is the rule that turns a four-point family written in p1..p3 into a
    three-point one; a record's own conservation rule (kinematics.yaml, the
    JSON "conservation" block, the .jl ext_subs) always takes precedence."""
    nums = sorted(int(e[1:]) for e in exts if re.fullmatch(r'p\d+', e))
    if nums and nums == list(range(1, len(nums) + 1)) and len(nums) >= 3:
        top = f"p{nums[-1]}"
        rest = " - ".join(f"p{n}" for n in nums[:-1])
        return {top: "-" + rest}
    return {}


def load_amflow_jl(path):
    """Parse an AmflowFamily(...) definition from a *.jl file."""
    txt = open(path).read()
    # locate AmflowFamily( ... ) with balanced parens
    i = txt.find("AmflowFamily(")
    if i < 0:
        raise ValueError("no AmflowFamily( in " + path)
    j = i + len("AmflowFamily(")
    depth = 1
    k = j
    while k < len(txt) and depth:
        if txt[k] == '(':
            depth += 1
        elif txt[k] == ')':
            depth -= 1
        k += 1
    body = txt[j:k - 1]
    # split top-level commas
    args = _split_top_commas(body)
    # arg0 name, arg1 loops list, arg2 exts list, arg3 ext-subs Dict,
    # arg4 kinematics Dict, arg5 propagators list
    name = args[0].strip().strip('"')
    loops = _jl_string_list(args[1])
    exts = _jl_string_list(args[2])
    ext_subs = _jl_dict(args[3])
    props_raw = _jl_string_list(args[5])
    props = [(p, 0) for p in props_raw]    # mass rides inside the string ("- msq")
    # physical index vector:  const <PREFIX>_INDICES = [1,1,...,0,0]
    physical = None
    mind = re.search(r'_INDICES\s*=\s*\[([^\]]*)\]', txt)
    if mind:
        idx = [int(x) for x in re.findall(r'-?\d+', mind.group(1))]
        if len(idx) == len(props):
            physical = [v != 0 for v in idx]
    return [Family(name, loops, exts, ext_subs, props, path, physical)]


def _split_top_commas(s):
    out, depth, cur = [], 0, ""
    for ch in s:
        if ch in "([{":
            depth += 1
        elif ch in ")]}":
            depth -= 1
        if ch == ',' and depth == 0:
            out.append(cur)
            cur = ""
        else:
            cur += ch
    if cur.strip():
        out.append(cur)
    return out


def _jl_string_list(s):
    return re.findall(r'"([^"]*)"', s)


def _jl_dict(s):
    out = {}
    for m in re.finditer(r'"([^"]*)"\s*=>\s*"([^"]*)"', s):
        out[m.group(1)] = m.group(2)
    return out


def load_amflow_json(path, integral_index=None):
    """
    Parse an AMFlow-port JSON config:
        {"family": {"name", "loops", "legs", "conservation", "replacement",
                    "propagators"}, "integrals": [{"indices": [...]}, ...],
         "amf_options": {"blackbox": {"numeric_values": {...}}}}
    or the bare family block alone (a propagator list written by hand in the
    same keys, with optional "indices" / "nu" and "numeric_values").

    Legs and momentum conservation are taken from the file, never inferred:
    a family block with two or more legs and NO "conservation" key is refused
    by name ("conservation": {} declares the listed legs independent).  The
    top sector is the integral with the most positive indices (first among
    ties) unless integral_index picks one; its index vector is Family.nu.
    Mass symbols = numeric_values keys that are not kinematic symbols (those
    named in the replacement rules); they are carried in Family.mass_values
    and the propagator strings stay verbatim.
    """
    d = json.load(open(path))
    if isinstance(d, dict) and "family" in d:
        fb = d["family"]
    elif isinstance(d, dict) and "propagators" in d:
        fb = d
    else:
        raise ValueError("REFUSED: no 'family' block and no 'propagators' "
                         "list in " + path)
    for key in ("name", "loops", "propagators"):
        if key not in fb:
            raise ValueError(f"REFUSED: family block lacks '{key}' in {path}")
    loops = [str(x) for x in fb["loops"]]
    legs = [str(x) for x in fb.get("legs", [])]
    if not legs:
        raise ValueError("REFUSED: family block declares no 'legs' in " + path)
    if "conservation" not in fb:
        if len(legs) >= 2:
            raise ValueError(
                "REFUSED: 'conservation' is not declared for the %d legs %s of "
                "family %r in %s; momentum conservation is never inferred "
                "(write \"conservation\": {} to declare the listed legs "
                "independent)" % (len(legs), legs, fb["name"], path))
        cons = {}
    else:
        cons = {str(k): str(v) for k, v in dict(fb["conservation"]).items()}
    for dep in cons:
        if dep not in legs:
            raise ValueError(f"REFUSED: conservation names {dep!r}, not one of "
                             f"the legs {legs} in {path}")
    props = []
    for p in fb["propagators"]:
        if isinstance(p, str):
            props.append((p, 0))
        elif isinstance(p, (list, tuple)) and len(p) == 2:
            props.append((str(p[0]), p[1]))
        else:
            raise ValueError(f"REFUSED: propagator entry {p!r} is neither a "
                             f"string nor [expr, mass] in {path}")
    integrals = []
    if fb is d:
        for key in ("indices", "nu"):
            if fb.get(key) is not None:
                integrals.append([int(a) for a in fb[key]])
    else:
        for it in (d.get("integrals") or []):
            if isinstance(it, dict) and it.get("indices") is not None:
                integrals.append([int(a) for a in it["indices"]])
    for v in integrals:
        if len(v) != len(props):
            raise ValueError(f"REFUSED: an index vector has {len(v)} entries "
                             f"for {len(props)} propagators in {path}")
    nu = None
    physical = None
    integral_index_given = integral_index is not None
    if integrals:
        if integral_index is None:
            npos = [sum(1 for a in v if a > 0) for v in integrals]
            integral_index = npos.index(max(npos))
        if not 0 <= integral_index < len(integrals):
            raise ValueError(f"REFUSED: integral index {integral_index} out of "
                             f"range ({len(integrals)} integrals) in {path}")
        nu = list(integrals[integral_index])
        physical = [a > 0 for a in nu]
    rep = {str(k): str(v) for k, v in dict(fb.get("replacement", {})).items()}
    leg_virt = {}
    for p in legs:
        for key, val in rep.items():
            if key.replace(" ", "") == p + "^2":
                leg_virt[p] = val
    if fb is d:
        nv = fb.get("numeric_values") or {}
    else:
        nv = ((d.get("amf_options") or {}).get("blackbox") or {}).get(
            "numeric_values") or {}
    kin_syms = set()
    for val in rep.values():
        kin_syms |= set(re.findall(r'[A-Za-z_]\w*', str(val)))
    mass_values = {str(k): str(v) for k, v in nv.items() if k not in kin_syms}
    notes = {"integral_index": integral_index, "n_integrals": len(integrals),
             # every index vector the record lists, so the compare's evidence can
             # name the one it read beside the others when there are several
             "integrals": [list(v) for v in integrals],
             "integral_index_given": integral_index_given,
             "replacement": rep, "numeric_values": {str(k): str(v) for k, v in nv.items()},
             "kinematic_symbols": sorted(kin_syms)}
    fam = Family(str(fb["name"]), loops, legs, cons, props, path, physical,
                 nu=nu, leg_virt=leg_virt, mass_values=mass_values,
                 kinematics_source="json:legs+conservation",
                 source_kind="amflow_json", record_notes=notes)
    return [fam]


def route_edges(vertices, edges, legs, loop_names=None, dependent_leg=None):
    """
    Spanning-tree routing of a drawn graph: propagator momenta from an edge
    list.  edges = [(u, v, mass[, nu])] (momentum flows u -> v), legs =
    [(name, vertex[, virtuality])] all incoming; the dependent leg (default:
    the last one) is -(sum of the others).  The spanning tree is grown in
    edge order (union-find); the chords carry the loop momenta in order; the
    tree edges follow from conservation at every vertex.

    Returns a dict: props [(momentum string, mass, nu)], loops, exts,
    ext_subs, tree, chords, vertex_conservation (every vertex must read "0"),
    degrees, E, V, N, L.  Refuses a disconnected graph, a leg at an unknown
    vertex, or a loop count that disagrees with loop_names.
    """
    verts = [str(v) for v in vertices]
    E = [(str(e[0]), str(e[1]), (e[2] if len(e) > 2 else 0),
          (int(e[3]) if len(e) > 3 else 1)) for e in edges]
    for u, v, _m, _n in E:
        for x in (u, v):
            if x not in verts:
                verts.append(x)
    L = [(str(lg[0]), str(lg[1]), (str(lg[2]) if len(lg) > 2 else None))
         for lg in legs]
    for name, at, _v in L:
        if at not in verts:
            raise ValueError(f"REFUSED: leg {name} sits at unknown vertex {at!r}")
    parent = {v: v for v in verts}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    tree, chords = [], []
    for i, (u, v, _m, _n) in enumerate(E):
        ru, rv = find(u), find(v)
        if ru == rv:
            chords.append(i)
        else:
            parent[ru] = rv
            tree.append(i)
    if len({find(v) for v in verts}) != 1:
        raise ValueError("REFUSED: the drawn graph is disconnected")
    n_loops = len(E) - len(verts) + 1
    if loop_names is None:
        loop_names = [f"k{i + 1}" for i in range(n_loops)]
    loop_names = [str(x) for x in loop_names]
    if len(loop_names) != n_loops:
        raise ValueError(f"REFUSED: E - V + 1 = {n_loops} loops for the edge "
                         f"list but loop names {loop_names} were given")
    exts = [name for name, _at, _v in L]
    if dependent_leg is None:
        dependent_leg = exts[-1] if exts else None
    if dependent_leg is not None and dependent_leg not in exts:
        raise ValueError(f"REFUSED: dependent leg {dependent_leg!r} is not one "
                         f"of the legs {exts}")
    syms = {n: sp.Symbol(n) for n in loop_names + exts}
    indep = [e for e in exts if e != dependent_leg]
    ext_subs = {}
    dep_sub = {}
    if dependent_leg is not None and len(exts) > 1:
        ext_subs = {dependent_leg: "-" + " - ".join(indep)}
        dep_sub = {syms[dependent_leg]: -sum(syms[e] for e in indep)}
    q = {}
    for li, i in enumerate(chords):
        q[i] = syms[loop_names[li]]
    tsym = {i: sp.Symbol("Qtree%d" % i) for i in tree}
    for i in tree:
        q[i] = tsym[i]
    eqs = []
    for v in verts:
        acc = sp.Integer(0)
        for i, (a, b, _m, _n) in enumerate(E):
            if b == v:
                acc += q[i]
            if a == v:
                acc -= q[i]
        for name, at, _v in L:
            if at == v:
                acc += syms[name]
        eqs.append(sp.expand(acc.subs(dep_sub)))
    sol = sp.solve(eqs, list(tsym.values()), dict=True) if tree else [{}]
    if not sol:
        raise ValueError("REFUSED: momentum conservation has no solution on "
                         "the drawn incidence")
    sol = sol[0]
    props = []
    for i, (u, v, m, n) in enumerate(E):
        expr = sp.expand(q[i].subs(sol).subs(dep_sub))
        if any(str(s).startswith("Qtree") for s in expr.free_symbols):
            raise ValueError(f"REFUSED: tree edge {u}-{v} left unsolved")
        props.append((sp.sstr(expr), m, n))
    checks, degrees = {}, {}
    for v in verts:
        acc = sp.Integer(0)
        deg = 0
        for i, (a, b, _m, _n) in enumerate(E):
            if b == v:
                acc += q[i].subs(sol)
                deg += 1
            if a == v:
                acc -= q[i].subs(sol)
                deg += 1
        for name, at, _v in L:
            if at == v:
                acc += syms[name]
                deg += 1
        checks[v] = sp.sstr(sp.expand(acc.subs(dep_sub)))
        degrees[v] = deg
    if any(c != "0" for c in checks.values()):
        raise ValueError("REFUSED: a routed vertex does not conserve momentum: "
                         + str(checks))
    return {"props": props, "loops": loop_names, "exts": exts,
            "ext_subs": ext_subs, "dependent_leg": dependent_leg,
            "leg_virt": {name: virt for name, _at, virt in L if virt is not None},
            "legs_at": {name: at for name, at, _v in L},
            "tree": [E[i][:2] for i in tree], "chords": [E[i][:2] for i in chords],
            "vertex_conservation": checks, "degrees": degrees,
            "E": len(E), "V": len(verts), "N": len(exts), "L": n_loops}


def load_pysecdec_graph(path):
    """
    Parse a pySecDec script's LoopIntegralFromGraph(internal_lines=[[mass,
    [u, v]], ...], external_lines=[[name, vertex], ...], replacement_rules=
    [(expr, value), ...]) call (literal arguments only) and route it with
    route_edges.  pySecDec masses are MASSES (the propagator is q^2 - m^2):
    the mass^2 field is written as (m)**2.  The optional powerlist gives nu.
    The family name is loop_package(name=...) when present, else the file
    stem.
    """
    import ast
    tree = ast.parse(open(path).read(), filename=path)
    call = None
    name = os.path.splitext(os.path.basename(path))[0]
    for node in ast.walk(tree):
        if isinstance(node, ast.Call):
            fn = node.func
            fname = fn.attr if isinstance(fn, ast.Attribute) else getattr(fn, "id", "")
            if fname == "LoopIntegralFromGraph" and call is None:
                call = node
            if fname == "loop_package":
                for kw in node.keywords:
                    if kw.arg == "name":
                        try:
                            name = str(ast.literal_eval(kw.value))
                        except Exception:
                            pass
    if call is None:
        raise ValueError("REFUSED: no LoopIntegralFromGraph(...) call in " + path)
    kw = {}
    for k in call.keywords:
        try:
            kw[k.arg] = ast.literal_eval(k.value)
        except Exception:
            raise ValueError(f"REFUSED: LoopIntegralFromGraph argument "
                             f"{k.arg!r} is not a literal in {path}")
    if len(call.args) >= 1 and "internal_lines" not in kw:
        kw["internal_lines"] = ast.literal_eval(call.args[0])
    if len(call.args) >= 2 and "external_lines" not in kw:
        kw["external_lines"] = ast.literal_eval(call.args[1])
    for key in ("internal_lines", "external_lines"):
        if key not in kw:
            raise ValueError(f"REFUSED: LoopIntegralFromGraph lacks {key} in {path}")
    powers = kw.get("powerlist")
    edges = []
    for i, line in enumerate(kw["internal_lines"]):
        mtok, (u, v) = str(line[0]), line[1]
        msq = 0 if mtok.strip() in ("0", "") else f"({mtok})**2"
        nu = int(powers[i]) if powers is not None else 1
        edges.append((str(u), str(v), msq, nu))
    verts = []
    for u, v, _m, _n in edges:
        for x in (u, v):
            if x not in verts:
                verts.append(x)
    virt = {}
    for rule in kw.get("replacement_rules", []) or []:
        key = str(rule[0]).replace(" ", "")
        m = re.fullmatch(r'([A-Za-z_]\w*)\*\1|([A-Za-z_]\w*)\*\*2|([A-Za-z_]\w*)\^2', key)
        if m:
            virt[m.group(1) or m.group(2) or m.group(3)] = str(rule[1])
    legs = [(str(lg[0]), str(lg[1]), virt.get(str(lg[0]))) for lg in kw["external_lines"]]
    if legs and all(re.fullmatch(r'p\d+', n) for n, _at, _v in legs):
        legs.sort(key=lambda t: int(t[0][1:]))     # basis convention p1, p2, ...
    r = route_edges(verts, edges, legs)
    props = [(expr, m) for expr, m, _n in r["props"]]
    nu = [n for _e, _m, n in r["props"]]
    notes = {"routing": {"tree": r["tree"], "chords": r["chords"]},
             "vertex_conservation": r["vertex_conservation"],
             "degrees": r["degrees"], "E": r["E"], "V": r["V"], "N": r["N"],
             "L": r["L"], "legs_at": r["legs_at"],
             "replacement_rules": [[str(a), str(b)] for a, b in
                                   (kw.get("replacement_rules") or [])]}
    fam = Family(name, r["loops"], r["exts"], r["ext_subs"], props, path,
                 [n > 0 for n in nu], nu=nu, leg_virt=r["leg_virt"],
                 kinematics_source="graph:external_lines",
                 cyclic_leg_order=list(r["exts"]), source_kind="pysecdec_graph",
                 record_notes=notes)
    return [fam]


def load_edge_list(path):
    """
    Read a drawn graph as an edge list and route it (route_edges).  Forms:
      * JSON {"name", "edges": [[u, v, mass(, nu)] | {"u","v","mass","nu"}],
        "legs": [[name, vertex(, virtuality)] | {"name","vertex","virtuality"}]
        or "nodes": {vertex: [leg numbers]} with "kinematics":
        {"external_masses": {leg number: p^2}}, optional "vertices",
        "loop_momenta", "cyclic_leg_order", "dependent_leg"}.
      * a yaml "drawn_graph:" block with the same keys, one flow-style item
        per line (- [a1, a2, msq]).
      * the routed-yaml header form: "# edge list: a1-a2[msq]; ..." and
        "# legs: p2@a1, ..." lines (loop_momenta, leg_virtualities,
        momentum_conservation and a "cyclic leg order ... [..]" comment are
        read from the same file when present).
    """
    txt = open(path).read()
    spec = None
    if path.lower().endswith(".json"):
        spec = json.load(open(path))
    else:
        blk = re.search(r'(?ms)^drawn_graph\s*:\s*\n(.*?)(?=^\S|\Z)', txt)
        if blk:
            spec = {}
            cur = None
            for line in blk.group(1).splitlines():
                s = line.strip()
                if not s or s.startswith("#"):
                    continue
                m = re.match(r'-\s*\[(.*)\]\s*$', s)
                if m and cur is not None:
                    spec[cur].append([x.strip().strip('"\'') for x in m.group(1).split(",")])
                    continue
                m = re.match(r'([A-Za-z_]\w*)\s*:\s*(.*)$', s)
                if not m:
                    raise ValueError(f"REFUSED: unparsed drawn_graph line {s!r} in {path}")
                key, val = m.group(1), m.group(2).strip()
                if val == "":
                    spec[key] = []
                    cur = key
                elif val.startswith("["):
                    spec[key] = [x.strip().strip('"\'') for x in val.strip("[]").split(",")
                                 if x.strip()]
                    cur = None
                else:
                    spec[key] = val.strip('"\'')
                    cur = None
        else:
            me = re.search(r'(?m)^#\s*edge list\s*:\s*(.+)$', txt)
            ml = re.search(r'(?m)^#\s*legs\s*:\s*(.+)$', txt)
            if not (me and ml):
                raise ValueError("REFUSED: no drawn_graph block, no JSON edge list "
                                 "and no '# edge list:' / '# legs:' header in " + path)
            edges = []
            for item in me.group(1).split(";"):
                item = item.strip()
                if not item:
                    continue
                m = re.fullmatch(r'([A-Za-z_]\w*)\s*-\s*([A-Za-z_]\w*)\s*\[([^\]]*)\]', item)
                if not m:
                    raise ValueError(f"REFUSED: unparsed edge {item!r} in {path}")
                edges.append([m.group(1), m.group(2), m.group(3).strip()])
            legs = []
            for item in ml.group(1).split(","):
                item = item.strip()
                if not item:
                    continue
                m = re.fullmatch(r'([A-Za-z_]\w*)\s*@\s*([A-Za-z_]\w*)', item)
                if not m:
                    raise ValueError(f"REFUSED: unparsed leg {item!r} in {path}")
                legs.append([m.group(1), m.group(2)])
            spec = {"edges": edges, "legs": legs}
            mn = re.search(r'(?m)^\s*-?\s*name\s*:\s*"?([^"\n]+)"?\s*$', txt)
            if mn:
                spec["name"] = mn.group(1).strip()
            mlp = re.search(r'loop_momenta\s*:\s*\[([^\]]*)\]', txt)
            if mlp:
                spec["loop_momenta"] = [x.strip() for x in mlp.group(1).split(",") if x.strip()]
            mv = re.search(r'leg_virtualities\s*:\s*\{([^}]*)\}', txt)
            if mv:
                spec["leg_virtualities"] = {a: b for a, b in
                                            re.findall(r'(\w+)\s*:\s*"([^"]*)"', mv.group(1))}
            mc = re.search(r'momentum_conservation\s*:\s*\{\s*(\w+)\s*:', txt)
            if mc:
                spec["dependent_leg"] = mc.group(1)
            mo = re.search(r'cyclic leg order[^\[\n]*\[([^\]]*)\]', txt)
            if mo:
                spec["cyclic_leg_order"] = [x.strip().strip('"\'') for x in
                                            mo.group(1).split(",") if x.strip()]
    if not isinstance(spec, dict) or "edges" not in spec:
        raise ValueError("REFUSED: no 'edges' in the drawn-graph spec of " + path)
    edges = []
    for e in spec["edges"]:
        if isinstance(e, dict):
            edges.append((e["u"], e["v"], e.get("mass", 0), e.get("nu", 1)))
        else:
            e = list(e)
            edges.append((e[0], e[1], e[2] if len(e) > 2 else 0,
                          e[3] if len(e) > 3 else 1))
    legs = []
    if "legs" in spec:
        for lg in spec["legs"]:
            if isinstance(lg, dict):
                legs.append((lg["name"], lg["vertex"], lg.get("virtuality")))
            else:
                lg = list(lg)
                legs.append((lg[0], lg[1], lg[2] if len(lg) > 2 else None))
    elif "nodes" in spec:
        ext_m = ((spec.get("kinematics") or {}).get("external_masses") or {})
        for at, nums in spec["nodes"].items():
            for n in nums:
                legs.append((f"p{n}", str(at), (str(ext_m[str(n)])
                                                if str(n) in ext_m else None)))
        legs.sort(key=lambda t: int(t[0][1:]))
    else:
        raise ValueError("REFUSED: no 'legs' (or 'nodes') in the drawn-graph "
                         "spec of " + path)
    virt = {str(k): str(v) for k, v in (spec.get("leg_virtualities") or {}).items()}
    legs = [(n, at, (v if v is not None else virt.get(str(n)))) for n, at, v in legs]
    if legs and all(re.fullmatch(r'p\d+', str(n)) for n, _at, _v in legs):
        # the tool's basis convention: p1, p2, ... in index order (the last is
        # the dependent leg unless the spec names one)
        legs.sort(key=lambda t: int(str(t[0])[1:]))
    verts = [str(v) for v in (spec.get("vertices") or [])]
    loops = spec.get("loop_momenta") or spec.get("loop_names")
    r = route_edges(verts, edges, legs, loop_names=loops,
                    dependent_leg=spec.get("dependent_leg"))
    props = [(expr, m) for expr, m, _n in r["props"]]
    nu = [n for _e, _m, n in r["props"]]
    order = spec.get("cyclic_leg_order")
    notes = {"routing": {"tree": r["tree"], "chords": r["chords"]},
             "vertex_conservation": r["vertex_conservation"],
             "degrees": r["degrees"], "E": r["E"], "V": r["V"], "N": r["N"],
             "L": r["L"], "legs_at": r["legs_at"]}
    name = str(spec.get("name") or os.path.splitext(os.path.basename(path))[0])
    fam = Family(name, r["loops"], r["exts"], r["ext_subs"], props, path,
                 [n > 0 for n in nu], nu=nu, leg_virt=r["leg_virt"],
                 kinematics_source="graph:legs",
                 cyclic_leg_order=(list(order) if order else None),
                 source_kind="edge_list", record_notes=notes)
    return [fam]


def load_family(path, name=None, kinematics=None, integral_index=None):
    """Dispatch on the record form (see the header): .yaml/.yml -> Kira yaml,
    or an edge list when the file carries a drawn_graph: block or only an
    '# edge list:' header; .jl -> AmflowFamily; .json -> AMFlow-port config
    (an 'edges' list -> edge list); .py -> pySecDec graph.

    Kinematics: a Kira yaml reads its legs as DATA from `kinematics` (a
    kinematics.yaml path, the --kinematics flag) or, when none is given, from
    a file named kinematics.yaml beside the family file (Kira's own config
    layout); the path used is recorded on every family
    (Family.kinematics_path).  The other readers carry their record's own
    legs (JSON legs + conservation, .jl exts + ext_subs, a drawn graph's
    legs); a kinematics file given with one of those is recorded, not
    applied.  Every family then passes through complete_legs(): the
    single-current second leg, the implicit p_N, composite legs."""
    low = path.lower()
    kin_used = None
    if low.endswith(".yaml") or low.endswith(".yml"):
        txt = open(path).read()
        if re.search(r'(?m)^drawn_graph\s*:', txt) or (
                not re.search(r'(?m)^integralfamilies\s*:', txt)
                and re.search(r'(?m)^#\s*edge list\s*:', txt)):
            fams = load_edge_list(path)
        else:
            kin_used = kinematics
            if kin_used is None:
                beside = os.path.join(os.path.dirname(os.path.abspath(path)), "kinematics.yaml")
                if os.path.exists(beside):
                    kin_used = beside
            fams = load_kira_yaml(path, kinematics=kin_used)
            for f in fams:
                f.source_kind = "kira_yaml"
                if kin_used is not None and kinematics is None:
                    f.kinematics_source = "kinematics.yaml beside the family file"
                    f.record_notes.setdefault("leg_set", {})["source"] = f.kinematics_source
    elif low.endswith(".jl"):
        fams = load_amflow_jl(path)
        for f in fams:
            f.source_kind = "amflow_jl"
            f.kinematics_source = "jl:exts+ext_subs"
    elif low.endswith(".json"):
        d = json.load(open(path))
        if isinstance(d, dict) and "edges" in d and "family" not in d:
            fams = load_edge_list(path)
        else:
            fams = load_amflow_json(path, integral_index=integral_index)
    elif low.endswith(".py"):
        fams = load_pysecdec_graph(path)
    else:
        raise ValueError("unknown family file type: " + path)
    if name:
        fams = [f for f in fams if f.name == name]
    kin = read_kinematics_yaml(kin_used) if kin_used else None
    for f in fams:
        if kinematics is not None:
            f.kinematics_path = kinematics
        elif kin_used is not None:
            f.kinematics_path = kin_used
        if kinematics is not None and f.source_kind != "kira_yaml":
            f.record_notes["kinematics_file_not_applied"] = (
                "%s: the %s record declares its own legs; the kinematics file is "
                "recorded, not applied" % (os.path.basename(str(kinematics)), f.source_kind))
        complete_legs(f, kin)
    return fams


# ---------------------------------------------------------------------------
#  1b. KINEMATICS: the external-leg set as DATA, and its completion
# ---------------------------------------------------------------------------
def read_kinematics_yaml(path):
    """
    Read a Kira kinematics.yaml (hand parser, no pyyaml):

      kinematics:
        incoming_momenta: [p1, p2, p3]        outgoing_momenta: [p4]
        momentum_conservation: [p4, -p1-p2-p3]   (or [] : all legs independent)
        scalarproduct_rules:
          - [[p1, p1], 0]     - [[p1, p2], "s/2"]     - [ "p1*p1", 0 ]
        symbol_to_replace_by_one: msq

    Returns {"incoming", "outgoing", "legs", "conservation": {dep: expr},
    "rules": {(a, b): expr}, "invariants", "symbol_to_replace_by_one",
    "virtualities": {leg: p^2 expr}, "path"}.  "legs" lists every external
    leg: the declared incoming and outgoing ones in declared order, then a
    dependent leg that the conservation rule names but the lists do not (the
    AuxLeg form: `momentum_conservation: [fooAuxLeg, -(p)]`), which keeps
    its declared name.  Virtualities come from the rules (leg_virtualities);
    a leg with no rule for one of its terms is left out of the dict.  The
    file's in/out lists are NOT taken as the momentum-flow orientation of the
    legs (Kira uses them only to name the symbols); the conservation rule is.
    """
    txt = open(path).read()
    body = re.sub(r'(?m)#.*$', '', txt)

    def lst(key):
        m = re.search(r'(?m)^\s*%s\s*:\s*\[([^\]]*)\]' % key, body)
        if not m:
            return []
        return [x.strip().strip('"\'') for x in m.group(1).split(",") if x.strip()]
    incoming = lst("incoming_momenta")
    outgoing = lst("outgoing_momenta")
    cons = {}
    m = re.search(r'(?m)^\s*momentum_conservation\s*:\s*\[\s*([A-Za-z_]\w*)\s*,\s*(.*?)\s*\]\s*$', body)
    if m:
        cons[m.group(1)] = m.group(2).strip().strip('"\'')
    rules = {}
    for mm in re.finditer(r'-\s*\[\s*\[\s*(\w+)\s*,\s*(\w+)\s*\]\s*,\s*(?:"([^"]*)"|\'([^\']*)\'|([^\]\s]+))\s*\]', body):
        a, b, v1, v2, v3 = mm.groups()
        rules[(a, b)] = (v1 if v1 is not None else (v2 if v2 is not None else v3)).strip()
    for mm in re.finditer(r'-\s*\[\s*"(\w+)\s*\*\s*(\w+)"\s*,\s*(?:"([^"]*)"|\'([^\']*)\'|([^\]\s]+))\s*\]', body):
        a, b, v1, v2, v3 = mm.groups()
        rules[(a, b)] = (v1 if v1 is not None else (v2 if v2 is not None else v3)).strip()
    inv = []
    mi = re.search(r'(?ms)^\s*kinematic_invariants\s*:\s*\n(.*?)(?=^\s*\w+\s*:|\Z)', body)
    if mi:
        inv = re.findall(r'-\s*\[\s*(\w+)\s*,', mi.group(1))
    one = None
    mo = re.search(r'(?m)^\s*symbol_to_replace_by_one\s*:\s*(\w+)', body)
    if mo:
        one = mo.group(1)
    legs = list(incoming) + [x for x in outgoing if x not in incoming]
    for dep in cons:
        if dep not in legs:
            legs.append(dep)
    kin = {"incoming": incoming, "outgoing": outgoing, "legs": legs, "conservation": cons,
           "rules": rules, "invariants": inv, "symbol_to_replace_by_one": one,
           "path": path}
    kin["virtualities"] = leg_virtualities(kin, legs, cons)
    return kin


def leg_virtualities(kin, legs=None, ext_subs=None, leg_exprs=None):
    """
    p^2 of every external leg from the kinematics' scalar-product rules.
    legs / ext_subs default to the kinematics' own; a dependent leg's p^2 is
    the expansion of its rule ((sum c_i p_i)^2 through the rules), and
    leg_exprs = {leg: sympy expression in the declared symbols} overrides the
    momentum of any leg (used for a completed implicit leg).  A leg is left
    out when a rule its expansion needs is absent (never guessed).  Values are
    sympy strings ("0" on shell).
    """
    rules = {}
    for (a, b), v in kin.get("rules", {}).items():
        try:
            rules[frozenset([a, b])] = sp.sympify(v)
        except Exception:
            continue
    legs = list(legs if legs is not None else kin.get("legs", []))
    ext_subs = dict(ext_subs if ext_subs is not None else kin.get("conservation", {}))
    syms = {}
    for lg in legs:
        syms[lg] = sp.Symbol(lg)
    for expr in list(ext_subs.values()) + list((leg_exprs or {}).values()):
        for tok in re.findall(r'[A-Za-z_]\w*', str(expr)):
            syms.setdefault(tok, sp.Symbol(tok))
    out = {}
    for lg in legs:
        if leg_exprs and lg in leg_exprs:
            e = sp.expand(sp.sympify(leg_exprs[lg], locals=syms))
        elif lg in ext_subs:
            e = sp.expand(sp.sympify(ext_subs[lg], locals=syms))
        else:
            e = syms[lg]
        terms = [(c, s) for s, c in e.as_coefficients_dict().items()]
        tot = sp.Integer(0)
        ok = True
        for c1, s1 in terms:
            for c2, s2 in terms:
                key = frozenset([str(s1), str(s2)])
                if key not in rules:
                    ok = False
                    break
                tot += c1 * c2 * rules[key]
            if not ok:
                break
        if ok:
            out[lg] = str(sp.nsimplify(sp.simplify(tot)))
    return out


def _kin_from_replacement(rep, legs):
    """
    A kinematics-like dict ({"rules": {(a, b): expr}}) from an AMFlow-port
    JSON "replacement" block: "p^2" -> (p, p); "(pi+pj)^2" / "(pi-pj)^2" ->
    the dot product (value - pi^2 - pj^2)/2 (resp. -(...)/2) when both
    squares are declared; "pi*pj" -> (pi, pj).  Only what the block states
    is derived; leg_virtualities leaves out any leg whose expansion needs a
    rule that is absent.
    """
    rules = {}
    squares = {}
    pairs = []
    for k, v in (rep or {}).items():
        ks = str(k).replace(" ", "")
        m = re.fullmatch(r'([A-Za-z_]\w*)\^2', ks) or re.fullmatch(r'\(([A-Za-z_]\w*)\)\^2', ks)
        if m:
            squares[m.group(1)] = sp.sympify(str(v))
            rules[(m.group(1), m.group(1))] = str(v)
            continue
        m = re.fullmatch(r'([A-Za-z_]\w*)\*([A-Za-z_]\w*)', ks)
        if m:
            rules[(m.group(1), m.group(2))] = str(v)
            continue
        m = re.fullmatch(r'\(([A-Za-z_]\w*)([+-])([A-Za-z_]\w*)\)\^2', ks)
        if m:
            pairs.append((m.group(1), m.group(2), m.group(3), sp.sympify(str(v))))
    for a, sign, b, val in pairs:
        if a in squares and b in squares:
            dot = (val - squares[a] - squares[b]) / 2
            if sign == "-":
                dot = -dot
            rules[(a, b)] = str(sp.expand(dot))
    return {"rules": rules, "legs": list(legs), "conservation": {}}


def _leg_implicit_name(exts):
    """Name of the implicit (completing) leg: pB for a single current p (xB
    for a single current x), the first gap of a p-numbered list else p<K+1>
    (p1..p4 -> p5; p1, p2, p4 -> p3), the missing digit pair for composite
    P<ij> legs (P12, P56 -> P34), else Pdep."""
    if len(exts) == 1:
        return exts[0] + "B"
    if all(re.fullmatch(r'p\d+', e) for e in exts):
        nums = sorted(int(e[1:]) for e in exts)
        gaps = [n for n in range(1, nums[-1] + 1) if n not in nums]
        return "p%d" % (gaps[0] if gaps else nums[-1] + 1)   # p3 for [p1, p2, p4]
    if all(re.fullmatch(r'P\d+', e) for e in exts):
        digs = [int(ch) for e in exts for ch in e[1:]]
        if len(set(digs)) == len(digs) and digs:
            missing = [d for d in range(1, max(digs) + 1) if d not in digs]
            if missing:
                return "P" + "".join(str(d) for d in missing)
    return "Pdep"


def _realizes(fam):
    """True when the top sector of fam admits a connected vertex realization
    (build_graph ok) -- used to resolve leg orientations, never printed."""
    try:
        G, ok = build_graph(fam)
    except Exception:
        return False
    return bool(ok)


def _rewrite_outgoing(props, outgoing):
    """Rewrite the propagator strings for legs that flow OUT: the token X of
    each outgoing leg becomes (-X), so that the leg symbol X denotes the
    INCOMING momentum of that leg (the realizer's all-incoming convention)
    while every propagator momentum stays what the record wrote."""
    if not outgoing:
        return list(props), {}
    out = []
    rew = {}
    for expr, mass in props:
        s = expr
        for x in outgoing:
            s2 = re.sub(r'(?<![\w.])%s(?![\w.])' % re.escape(x), "(-%s)" % x, s)
            if s2 != s:
                rew[x] = "%s -> (-%s): the record's %s flows out; the leg symbol now names its incoming momentum" % (x, x, x)
            s = s2
        out.append((s, mass))
    return out, rew


def complete_legs(fam, kin=None):
    """
    Complete the external-leg set of a family IN PLACE so that every leg the
    graph has is present with its momentum-conservation relation in the
    tool's all-incoming convention, and record what was done in
    fam.record_notes["leg_set"]:

      * a declared conservation rule (kinematics.yaml, JSON "conservation",
        .jl ext_subs, a drawn graph) is DATA: the dependent leg it names is
        added to the leg list when the lists omit it (the AuxLeg form; a
        dummy "...AuxLeg" name that no propagator spells is renamed to the
        completion name, pB for a single current p) and the rule is kept;
        nothing is second-guessed, a wrong rule fails the vertex realization
        by name downstream;
      * no rule and K declared legs (all independent): the graph has one more
        leg, the IMPLICIT leg, = -(sum of the declared legs) in the
        all-incoming convention -- pB = -p for a single current p (a two-point
        family), p_{K+1} for p1..pK (a five-point family written in p1..p4),
        the missing pair for composite legs (P12, P56 -> P34).  Which declared
        legs flow IN or OUT is not a datum Kira stores (its in/out lists only
        name symbols), so when the all-incoming completion does not realize
        the top sector the orientation is RESOLVED from the propagator
        momenta: the declared legs are flipped one subset at a time (fewest
        flips first, earliest declared leg first; the last declared leg is
        never flipped -- a global flip is the same graph) until the top
        sector realizes, and the flip is reported ("P12 outgoing"); if no
        pattern realizes, the all-incoming completion stands and the
        realizer reports the failure by name.  An outgoing leg X has its
        token rewritten to (-X) in the propagator strings (_rewrite_outgoing)
        so that X names the leg's incoming momentum; the rewrite is recorded.
      * composite (non-p<digits>) legs are put in the tool's basis order
        [incoming legs in declared order, the implicit leg, the outgoing
        legs], the last leg eliminated, so one graph has one canonical set
        whichever declared form the record uses (_normalize_composite_basis);
        p<digits> legs keep their declared order and elimination.
      * the canonical hash of the DECLARED leg set (the record's own basis,
        before completion) is kept as leg_set["declared_canonical_hash"];
      * leg virtualities of the completed / dependent legs are expanded from
        the kinematics' scalar-product rules when `kin` is given.
    """
    notes = dict(fam.record_notes.get("leg_set") or {})
    notes.setdefault("source", fam.kinematics_source)
    notes["declared_legs"] = list(notes.get("declared_legs") or fam.exts)
    exts = list(fam.exts)
    subs = dict(fam.ext_subs)
    if not exts:
        notes["completion"] = "none (no legs declared)"
        fam.record_notes["leg_set"] = notes
        return fam
    try:
        vd, _md, _i = momentum_vectors(fam)
        flat = sorted(tuple(str(c) for c in v) for v in canonical_set(vd))
        import hashlib as _hl
        notes["declared_canonical_hash"] = _hl.sha1(json.dumps(flat).encode()).hexdigest()[:12]
    except Exception:
        pass
    is_p = all(re.fullmatch(r'p\d+', e) for e in exts)
    spelled = set()
    for expr, _m in fam.propagators:
        spelled |= set(re.findall(r'[A-Za-z_]\w*', expr))
    unspelled = [e for e in exts if e not in subs and e not in spelled]
    if unspelled:
        raise ValueError(
            "REFUSED: LEG-NOT-SPELLED: the declared independent leg(s) %s of family %r "
            "appear in no propagator of the record (declared legs %s, source %s): the "
            "kinematics does not belong to this family, or the leg list names a momentum "
            "the family never carries; no leg set is completed from it"
            % (unspelled, fam.name, exts, notes.get("source")))
    if kin is None and fam.record_notes.get("replacement"):
        kin = _kin_from_replacement(fam.record_notes["replacement"], exts)
    # (a) a declared rule: keep it; add the dependent leg to the list
    if subs:
        added = []
        for dep in list(subs):
            if dep not in exts:
                exts.append(dep)
                added.append(dep)
        renamed = {}
        for dep in list(subs):
            if dep.endswith("AuxLeg") and dep not in spelled:
                others = [e for e in exts if e != dep]
                new = _leg_implicit_name(others) if others else dep
                while new in exts:
                    new += "B"
                exts = [new if e == dep else e for e in exts]
                subs = {(new if k == dep else k): v for k, v in subs.items()}
                if dep in fam.leg_virt:
                    fam.leg_virt[new] = fam.leg_virt.pop(dep)
                renamed[dep] = new
        notes["completion"] = ("declared rule %s" % subs
                               + ("; dependent leg %s added to the leg list" % added if added else "")
                               + ("; dummy dependent leg renamed %s" % renamed if renamed else ""))
        is_p = all(re.fullmatch(r'p\d+', e) for e in exts)
        if not is_p and len(subs) == 1:
            order, new_subs, outgoing, norm = _normalize_composite_basis(exts, subs)
            if norm:
                props, rew = _rewrite_outgoing(fam.propagators, outgoing)
                if rew:
                    notes["propagator_rewrite"] = rew
                    fam.propagators = props
                exts, subs = order, new_subs
                notes["basis_normalization"] = norm
                notes["orientation"] = {"resolved_from_propagators": False, "from_rule": True,
                                        "outgoing": outgoing,
                                        "note": "orientation read from the declared conservation rule"}
        fam.exts, fam.ext_subs = exts, subs
        if kin is not None:
            fam.leg_virt.update({k: v for k, v in
                                 leg_virtualities(kin, exts, subs).items()
                                 if k not in fam.leg_virt})
        fam.record_notes["leg_set"] = notes
        return fam
    # (b) no rule: the implicit leg, orientation resolved from the propagators
    K = len(exts)
    imp = _leg_implicit_name(exts)
    while imp in exts or imp in spelled:
        imp += "B"
    all_in = "-" + " - ".join(exts)
    tried = []
    chosen = None
    patterns = [()]
    if K >= 2:
        for r in range(1, K):
            for combo in itertools.combinations(range(K - 1), r):
                patterns.append(combo)
    for flips in patterns:
        outgoing = [exts[i] for i in flips]
        props_c, rew = _rewrite_outgoing(fam.propagators, outgoing)
        cand = Family(fam.name, fam.loops, exts + [imp], {imp: all_in}, props_c,
                      fam.source, fam.physical, nu=fam.nu, mass_values=fam.mass_values)
        ok = _realizes(cand)
        tried.append({"outgoing": outgoing, "realizes": ok})
        if ok:
            chosen = (outgoing, props_c, rew)
            break
    resolved = chosen is not None
    if not resolved:
        chosen = ([], list(fam.propagators), {})
    outgoing, props_c, rew = chosen
    notes["implicit_leg"] = {imp: all_in}
    notes["orientation"] = {
        "resolved_from_propagators": resolved,
        "outgoing": outgoing,
        "patterns_tried": tried,
        "note": ("all-incoming completion realizes the top sector" if resolved and not outgoing
                 else ("declared leg(s) %s flow OUT: the top sector conserves momentum at "
                       "every vertex only with that orientation (LEG-ORIENTATION-RESOLVED)"
                       % outgoing if resolved
                       else "no orientation of the declared legs realizes the top sector; "
                            "the all-incoming completion is kept and the realizer reports "
                            "the failure by name"))}
    notes["completion"] = ("implicit leg %s = %s (%s)"
                           % (imp, all_in, "single-current second leg" if K == 1 else
                              ("p_N of p1..p%d" % K if is_p else "composite legs")))
    if rew:
        notes["propagator_rewrite"] = rew
        fam.propagators = props_c
    new_exts = exts + [imp]
    new_subs = {imp: all_in}
    if not is_p and K >= 2:
        order, nsubs, out2, norm = _normalize_composite_basis(new_exts, new_subs, outgoing)
        if norm:
            new_exts, new_subs = order, nsubs
            notes["basis_normalization"] = norm
    fam.exts, fam.ext_subs = new_exts, new_subs
    if kin is not None:
        # the implicit leg's momentum in the record's symbols: -(sum of the
        # incoming declared legs) + (sum of the outgoing ones)
        syms = {e: sp.Symbol(e) for e in exts}
        q = -sum((-syms[e] if e in outgoing else syms[e]) for e in exts)
        leg_exprs = {imp: sp.sstr(sp.expand(q))}
        fam.leg_virt.update({k: v for k, v in
                             leg_virtualities(kin, new_exts, new_subs, leg_exprs).items()
                             if k not in fam.leg_virt})
    elif K == 1 and exts[0] in fam.leg_virt:
        fam.leg_virt[imp] = fam.leg_virt[exts[0]]
    fam.record_notes["leg_set"] = notes
    return fam


def _normalize_composite_basis(exts, subs, outgoing_known=None):
    """
    For composite (non-p<digits>) legs with ONE conservation rule dep = expr:
    write the rule as sum_i c_i P_i = 0 (c_dep = +1); when every c_i is +-1 it
    is a momentum-conservation relation, the all-incoming momenta are
    q_i = c_i P_i and the legs with c = -1 flow OUT.  Returns (order, subs,
    outgoing, note) in the tool's basis order: legs with c = +1 in their
    given order, then those with c = -1 (their symbols reinterpreted as the
    incoming momenta, see _rewrite_outgoing), the LAST leg eliminated through
    q_last = -(sum of the others).  With outgoing_known (the orientation
    already resolved and the strings already rewritten) the rule is read in
    the incoming symbols directly.  Any other rule (a coefficient that is not
    +-1, a symbol outside the leg list) is returned as given with note None.
    """
    (dep, expr), = subs.items()
    syms = {e: sp.Symbol(e) for e in exts}
    try:
        rel = sp.expand(syms[dep] - sp.sympify(expr, locals=syms))
    except Exception:
        return exts, subs, [], None
    coeffs = {}
    for e in exts:
        c = rel.coeff(syms[e])
        if c not in (1, -1, 0):
            return exts, subs, [], None
        coeffs[e] = int(c)
    rest = sp.expand(rel - sum(coeffs[e] * syms[e] for e in exts))
    if rest != 0 or coeffs[dep] != 1:
        return exts, subs, [], None
    if outgoing_known is not None:
        # the strings are already in incoming symbols: every c must be +1
        if any(coeffs[e] != 1 for e in exts):
            return exts, subs, [], None
        inc = [e for e in exts if e not in outgoing_known]
        out = [e for e in exts if e in outgoing_known]
    else:
        inc = [e for e in exts if coeffs[e] == 1]
        out = [e for e in exts if coeffs[e] == -1]
    if len(inc) + len(out) < 2:
        return exts, subs, [], None
    order = inc + out
    last = order[-1]
    # all incoming symbols: q_last = -(sum of the other q)
    others = [e for e in order[:-1]]
    new_subs = {last: "-" + " - ".join(others)}
    note = {"basis_order": order, "eliminated": last, "rule": new_subs,
            "outgoing": out,
            "why": "composite legs: basis = incoming legs in declared order, the implicit "
                   "leg, then outgoing legs; the last leg eliminated, so both declared forms "
                   "of one record hash alike"}
    return order, new_subs, out, note


def two_point_cut(fam):
    """
    For a family with exactly two external legs: the minimal number of
    internal lines a cut must sever to separate the one leg from the other
    (the two-point cut; the channel notion of cut_signature is empty here
    because {leg, leg} vs 'the rest' has no rest).  None when the graph does
    not realize.  Sunrise 3, kite 2, an L-loop banana L+1.
    """
    if len(fam.exts) != 2 or not _HAVE_NX:
        return None
    G, ok = build_graph(fam)
    if not ok or G is None:
        return None
    a, b = ("ext_%s" % fam.exts[0]), ("ext_%s" % fam.exts[1])
    if a not in G or b not in G:
        return None
    H = nx.Graph()
    INF = 10 ** 6
    for u, v, k in G.edges(keys=True):
        w = INF if str(k).startswith("ext_") else 1
        if H.has_edge(u, v):
            H[u][v]["capacity"] += w
        else:
            H.add_edge(u, v, capacity=w)
    try:
        val, _ = nx.minimum_cut(H, a, b)
    except Exception:
        return None
    return None if val >= INF else int(val)


# ---------------------------------------------------------------------------
#  2.  CANONICAL PROPAGATOR-VECTOR SET  +  ISOMORPHISM
# ---------------------------------------------------------------------------
class MomentumVectors(tuple):
    """The value of momentum_vectors(): a 3-tuple (vecs, masses, n_isp) that
    unpacks as before, carrying beside it the per-line LABELS the record
    declares -- .nu (the propagator power of each physical line: the index
    vector of the top-sector integral for a JSON / edge-list record, 1 for a
    line the record only flags physical, as a Kira top_level_sectors bitmask
    does), .index (the propagator index of each line in the record's list)
    and .nu_source ("record" | "assumed 1 (the record declares no index
    vector)")."""

    def __new__(cls, vecs, masses, isps, nu=None, index=None, nu_source=None):
        self = super().__new__(cls, (vecs, masses, isps))
        self.nu = list(nu) if nu is not None else [1] * len(vecs)
        self.index = list(index) if index is not None else list(range(len(vecs)))
        self.nu_source = nu_source or "assumed 1 (the record declares no index vector)"
        return self


def momentum_vectors(fam):
    """
    Return (vecs, masses, n_isp) for the PHYSICAL (top-sector) propagators only
    -- a MomentumVectors tuple whose .nu / .index / .nu_source carry the
    multiplicity (propagator power) and record index of every kept line.

    A propagator is excluded from the graph when it is (i) flagged non-physical
    by the top-sector mask / index vector, (ii) a pure dot-product ISP, or
    (iii) carries no loop momentum.  These are the numerator auxiliaries; they
    are not graph edges and must not enter the isomorphism / planarity / cut.
    """
    syms = _symbols(fam.loops, fam.exts)
    ext_resolved = _resolve_ext_subs(fam.exts, fam.ext_subs, syms)
    nloop = len(fam.loops)
    vecs, masses, isps, nus, idx = [], [], 0, [], []
    declared = getattr(fam, "nu", None)
    for i, (expr, mass) in enumerate(fam.propagators):
        if i < len(fam.physical) and not fam.physical[i]:
            isps += 1
            continue
        parsed = _parse_momentum(expr, fam.loops, fam.exts, syms, ext_resolved, mass)
        if parsed is None:
            isps += 1
            continue
        vec, _is_prop, m = parsed
        if all(c == 0 for c in vec[:nloop]):
            isps += 1            # pure external "propagator" -> not a graph edge
            continue
        vecs.append(vec)
        masses.append(m)
        nus.append(int(declared[i]) if declared is not None and i < len(declared) else 1)
        idx.append(i)
    return MomentumVectors(vecs, masses, isps, nus, idx,
                           "record" if declared is not None
                           else "assumed 1 (the record declares no index vector)")


def _canon_vec(vec):
    """Canonicalise a momentum up to overall sign (q^2 == (-q)^2): fix the sign
    of the first nonzero coefficient to be positive."""
    v = [sp.nsimplify(c) for c in vec]
    for c in v:
        if c != 0:
            if c < 0:
                v = [-x for x in v]
            break
    return tuple(v)


def canonical_set(vecs, masses=None, nu=None):
    """The canonical propagator set as a Counter (a multiset).  With `masses`
    and/or `nu` (lists aligned with `vecs`: the mass label string of each line
    and its multiplicity) the keys are LABELED, (canonical momentum, mass, nu);
    without labels the keys are the bare canonical momenta, so a caller that
    iterates the set or takes its length sees exactly what the former
    frozenset gave.  Two labeled sets are equal only when the same momenta
    carry the same masses with the same multiplicities."""
    import collections
    vecs = list(vecs)
    if masses is None and nu is None:
        return collections.Counter(_canon_vec(v) for v in vecs)
    ms = list(masses) if masses is not None else [""] * len(vecs)
    ns = list(nu) if nu is not None else [1] * len(vecs)
    if len(ms) != len(vecs) or len(ns) != len(vecs):
        raise ValueError("canonical_set: labels are not aligned with the vectors "
                         f"({len(vecs)} vectors, {len(ms)} masses, {len(ns)} multiplicities)")
    return collections.Counter((_canon_vec(v), str(m), int(n))
                               for v, m, n in zip(vecs, ms, ns))


def _mass_label(m, fam=None, mass_map=None):
    """The label string of a line's mass^2: the record's numeric values
    applied (Family.resolve_mass), then an optional {symbol: value} map (the
    other side's dictionary for a compare), nsimplified and printed -- "0",
    "1", "2", "7/25", "msq", "mt2"."""
    expr = sp.sympify(m)
    if fam is not None:
        expr = fam.resolve_mass(expr)
    if mass_map:
        expr = sp.sympify(expr).subs({sp.Symbol(str(k)): sp.nsimplify(v)
                                      for k, v in mass_map.items()})
    return str(sp.nsimplify(expr))


def _leg_class(virt):
    """The virtuality class of a leg label: "0" on shell, "offshell" for any
    non-zero virtuality (a symbol, a number, the literal "offshell"), None
    when the record declares nothing for the leg."""
    if virt is None:
        return None
    s = str(virt).strip()
    if s == "":
        return None
    return "0" if s in ("0", "0.0") else "offshell"


def labeled_set(fam, mass_map=None):
    """
    The LABELED canonical set of a family: the census's (momentum, mass,
    multiplicity) triples with the leg labels beside them.  Returns a dict:
      "lines"        Counter over (canonical momentum, mass label, nu)
      "masses"       mass label per physical line (aligned with momentum_vectors)
      "nu"           multiplicity per physical line, "nu_source" where it came from
      "mass_multiset", "nu_multiset"   the sorted multisets
      "n_massive_lines", "n_dotted_lines"   derived counts
      "leg_labels"   {leg: class} ("0" / "offshell" / None) from Family.leg_virt
      "leg_classes"  the class per leg in basis order, or None when any leg
                     has no declared virtuality (then no leg restriction applies)
    mass_map = {symbol: value} applied after the record's own numeric values
    (the dictionary that identifies a drawing's `msq` with the record's 1).
    """
    mv = momentum_vectors(fam)
    vecs, masses, _isps = mv
    labels = [_mass_label(m, fam, mass_map) for m in masses]
    lines = canonical_set(vecs, labels, mv.nu)
    leg_labels = {e: _leg_class(fam.leg_virt.get(e)) for e in fam.exts}
    classes = [leg_labels[e] for e in fam.exts]
    return {
        "lines": lines,
        "vecs": list(vecs),
        "masses": labels,
        "nu": list(mv.nu),
        "nu_source": mv.nu_source,
        "index": list(mv.index),
        "mass_multiset": sorted(labels),
        "nu_multiset": sorted(mv.nu),
        "n_massive_lines": sum(1 for m in labels if m != "0"),
        "n_dotted_lines": sum(1 for n in mv.nu if n > 1),
        "leg_labels": leg_labels,
        "leg_classes": classes if fam.exts and all(c is not None for c in classes) else None,
    }


def mass_multiset(fam, mass_map=None):
    """The sorted multiset of mass labels of the family's physical lines (the
    census lint's "mass_multiset"), e.g. ['0', '0', '1', '1', '2']."""
    return labeled_set(fam, mass_map)["mass_multiset"]


def labeled_isomorphism(famA, famB, try_leg_perms=True, mass_mapA=None, mass_mapB=None,
                        leg_labelsA=None, leg_labelsB=None, find_all=False):
    """
    The label-preserving isomorphism of two families under the tool's signed
    loop relabelings (and, with try_leg_perms, external-leg permutations with
    B's momentum conservation re-imposed), reported at three levels so that a
    mass-blind or leg-label-blind match is never mistaken for identity:
      "momentum_match"  the bare canonical momentum sets (the former notion)
      "mass_match"      (momentum, mass, nu) carried along, any leg permutation
      "labeled_match"   (momentum, mass, nu) AND only the leg permutations that
                        preserve the leg virtuality classes
    each a relabeling string or None (a list with find_all), with the mass
    and multiplicity multisets of both sides, whether they agree, the leg
    classes compared and how many leg permutations the labels admit.
    leg_labelsA / leg_labelsB: {leg: virtuality} overriding Family.leg_virt
    (the census's leg dictionary for a drawing whose yaml carries none).
    This is the census's labeled_isomorphism (census_lib_r21 L89) on the
    tool's own primitives; the class of maps is the signed-permutation one,
    so two routings of one graph that differ by a loop-momentum basis change
    are NOT identified here (the realized-graph test is the routing-
    independent one).
    """
    LA = labeled_set(famA, mass_mapA)
    LB = labeled_set(famB, mass_mapB)
    nloop, n_ext = len(famA.loops), len(famA.exts)
    out = {
        "n_lines": (len(LA["vecs"]), len(LB["vecs"])),
        "n_loops": (nloop, len(famB.loops)),
        "n_legs": (n_ext, len(famB.exts)),
        "mass_multiset_A": LA["mass_multiset"], "mass_multiset_B": LB["mass_multiset"],
        "mass_multisets_equal": LA["mass_multiset"] == LB["mass_multiset"],
        "nu_multiset_A": LA["nu_multiset"], "nu_multiset_B": LB["nu_multiset"],
        "nu_multisets_equal": LA["nu_multiset"] == LB["nu_multiset"],
        "nu_source": (LA["nu_source"], LB["nu_source"]),
        "momentum_match": None, "mass_match": None, "labeled_match": None,
        "leg_classes_A": None, "leg_classes_B": None, "leg_labels_compared": False,
        "n_leg_perms_admissible": None, "n_leg_perms_any": None,
        "n_relabelings_tried": {}, "reason": None,
    }
    if leg_labelsA is not None:
        out["leg_classes_A"] = [_leg_class(leg_labelsA.get(e)) for e in famA.exts]
    else:
        out["leg_classes_A"] = LA["leg_classes"]
    if leg_labelsB is not None:
        out["leg_classes_B"] = [_leg_class(leg_labelsB.get(e)) for e in famB.exts]
    else:
        out["leg_classes_B"] = LB["leg_classes"]
    if len(famB.loops) != nloop or len(famB.exts) != n_ext or len(LA["vecs"]) != len(LB["vecs"]):
        out["reason"] = ("line, loop or leg counts differ (%d/%d lines, %d/%d loops, %d/%d legs)"
                         % (len(LA["vecs"]), len(LB["vecs"]), nloop, len(famB.loops),
                            n_ext, len(famB.exts)))
        return out
    dep = _dependent_leg_rows(famB)
    cA, cB = out["leg_classes_A"], out["leg_classes_B"]
    labels_ok = (cA is not None and cB is not None and
                 all(c is not None for c in cA) and all(c is not None for c in cB))
    out["leg_labels_compared"] = labels_ok
    n_perms_any = len(list(itertools.permutations(range(n_ext)))) if try_leg_perms else 1
    out["n_leg_perms_any"] = n_perms_any
    if labels_ok and try_leg_perms:
        out["n_leg_perms_admissible"] = sum(
            1 for lp in itertools.permutations(range(n_ext))
            if all(cA[j] == cB[lp[j]] for j in range(n_ext)))
    else:
        out["n_leg_perms_admissible"] = n_perms_any if try_leg_perms else 1
    det = {}
    out["momentum_match"] = isomorphism_to(LA["vecs"], LB["vecs"], nloop, n_ext,
                                          try_leg_perms=try_leg_perms, dep_rows=dep,
                                          find_all=find_all, details=det)
    out["n_relabelings_tried"]["momentum"] = det.get("n_relabelings_tried")
    det = {}
    out["mass_match"] = isomorphism_to(LA["vecs"], LB["vecs"], nloop, n_ext,
                                      try_leg_perms=try_leg_perms, dep_rows=dep,
                                      find_all=find_all, details=det,
                                      labels=(list(zip(LA["masses"], LA["nu"])),
                                              list(zip(LB["masses"], LB["nu"]))))
    out["n_relabelings_tried"]["mass"] = det.get("n_relabelings_tried")
    if det.get("reason"):
        out["reason"] = det["reason"]
    det = {}
    out["labeled_match"] = isomorphism_to(LA["vecs"], LB["vecs"], nloop, n_ext,
                                         try_leg_perms=try_leg_perms, dep_rows=dep,
                                         find_all=find_all, details=det,
                                         leg_labels=(cA, cB) if labels_ok else None,
                                         labels=(list(zip(LA["masses"], LA["nu"])),
                                                 list(zip(LB["masses"], LB["nu"]))))
    out["n_relabelings_tried"]["labeled"] = det.get("n_relabelings_tried")
    out["labeled_matches"] = list(det.get("matches", []))
    found = out["labeled_match"] if not find_all else bool(out["labeled_match"])
    if not found and out["reason"] is None:
        if not labels_ok:
            out["reason"] = ("no signed loop relabeling + leg permutation maps the (momentum, "
                             "mass, multiplicity) sets onto each other (leg classes not "
                             "declared on both sides: no leg-label restriction applied)")
        elif (out["mass_match"] if not find_all else bool(out["mass_match"])):
            out["reason"] = ("same lines and masses under some relabeling, but no such "
                             "relabeling carries the off-shell-leg labels onto each other")
        elif out["n_leg_perms_admissible"] == 0:
            out["reason"] = ("the leg virtuality classes admit no leg permutation (%s vs %s)"
                             % (cA, cB))
        else:
            out["reason"] = ("no signed loop relabeling + label-preserving leg permutation "
                             "maps the labeled sets")
    return out


def _dependent_leg_rows(fam):
    """For every external leg the family eliminates by momentum conservation
    (Family.ext_subs), the coefficient row of its expression over the
    family's external block: {ext index: [coefficient of each ext]}.

    A momentum vector is written over (loops..., exts...) with the dependent
    legs already substituted away, so its dependent slots are zero.  A leg
    permutation can move a coefficient INTO a dependent slot; the row is what
    _reimpose_conservation needs to substitute it back out (the eliminated
    leg re-expressed through the independent ones).  Family.exts / ext_subs
    are read as data; nothing is inferred here."""
    syms = _symbols(fam.loops, fam.exts)
    rows = {}
    for j, e in enumerate(fam.exts):
        if e in fam.ext_subs:
            expr = sp.expand(sp.sympify(fam.ext_subs[e], locals=syms))
            rows[j] = [sp.nsimplify(expr.coeff(syms[e2], 1)) for e2 in fam.exts]
    return rows


def _reimpose_conservation(vec, nloop, dep):
    """Re-impose momentum conservation on a permuted momentum vector: every
    coefficient sitting in a dependent-leg slot (dep = _dependent_leg_rows of
    the TARGET family, whose basis the vector now lives in) is moved onto the
    independent legs through that leg's conservation row.  Without this step
    a permutation that lands a leg in the eliminated slot (p4 = -p1-p2-p3) is
    compared unreduced and a true relabeling is missed."""
    v = list(vec)
    for j, row in dep.items():
        c = v[nloop + j]
        if c != 0:
            v[nloop + j] = sp.Integer(0)
            for j2, r in enumerate(row):
                v[nloop + j2] += c * r
    return tuple(v)


def _apply_signed_perm_to_vec(vec, loop_perm, loop_signs, leg_perm, nloop, next_ext,
                              dep_rows=None):
    """Relabel a momentum vector: k_i -> sign_i * k_{perm(i)}, p_j -> p_{legperm(j)},
    then re-impose the target family's momentum conservation (dep_rows =
    _dependent_leg_rows of the family whose basis the result lives in) so a
    coefficient moved into an eliminated-leg slot is substituted back onto
    the independent legs.  dep_rows None = no re-imposition (a caller that
    has no target family).  A leg permutation that indexes past the vector's
    external block is refused by name (ValueError), never an IndexError."""
    new = [sp.Integer(0)] * len(vec)
    n_ext_vec = len(vec) - nloop
    # loop block
    for i in range(nloop):
        c = vec[i]
        if c == 0:
            continue
        tgt = loop_perm[i]
        new[tgt] += loop_signs[i] * c
    # external block
    for j in range(nloop, len(vec)):
        c = vec[j]
        if c == 0:
            continue
        if j - nloop >= len(leg_perm) or leg_perm[j - nloop] >= n_ext_vec:
            raise ValueError(
                f"leg permutation {tuple(leg_perm)} indexes past the vector's "
                f"external block ({n_ext_vec} legs): the two families do not "
                "have the same number of external legs")
        tgt = nloop + leg_perm[j - nloop]
        new[tgt] += c
    if dep_rows:
        return _reimpose_conservation(tuple(new), nloop, dep_rows)
    return tuple(new)


def isomorphism_to(vecsA, vecsB, nloop, n_ext, try_leg_perms=False,
                   dep_rows=None, leg_labels=None, find_all=False, details=None,
                   labels=None, affine=False, names=None):
    """
    Search for a signed loop-permutation (and optional leg permutation) that
    maps the canonical propagator SET of A onto that of B.  Returns a
    description string of the first match found, else None; with find_all
    the list of every matching relabeling (possibly empty).

    affine     = True: when no signed permutation maps the sets, search the
                 census's AFFINE class k_i -> sum_m U_im k_m + sum_j c_ij p_j
                 (U in unimodular_maps(nloop), the shift c solved from the
                 labeled set; _affine_search_vectors) over the same leg
                 permutations and labels, and return its describer
                 ('l1 -> k1 -k2 +p3, l2 -> -k2 +p1 +p2 +p3; p2 -> p4, p4 ->
                 p2') -- details["relabeling_class"] then reads "affine
                 unimodular loop redefinition" and details["affine"] carries
                 U, c, the leg map and the search bound; a signed match keeps
                 details["relabeling_class"] = "signed loop permutation".
                 With find_all the affine matches are appended after the
                 signed ones.  names = (loops_A, loops_B, exts_B) for the
                 describer's symbols (default k1.., p1..).

    dep_rows   = _dependent_leg_rows(family B): after every leg permutation
                 B's momentum conservation is re-imposed on the permuted
                 vectors, so permutations that move a leg into B's eliminated
                 slot are compared reduced (the row-21 p2 <-> p4 class).
                 None = compare the permuted coordinates as they are.
    leg_labels = (labels_A, labels_B), one class string per external index
                 ("0" on shell / "offshell"); only label-preserving leg
                 permutations are searched.  None = every permutation (no
                 labels present on either side).
    labels     = (labels_A, labels_B), one (mass label, multiplicity) pair per
                 VECTOR of each side (labeled_set(fam)["masses"] / ["nu"]):
                 the sets compared are then the LABELED canonical sets
                 (momentum, mass, nu), so a relabeling must carry every
                 line's mass and multiplicity along; when the mass or
                 multiplicity multisets of the two sides differ no relabeling
                 can exist and None / [] is returned with details["reason"]
                 naming the multiset (nothing is searched).  None = the bare
                 momentum sets (mass- and multiplicity-blind).
    details    = a dict the caller supplies to receive n_relabelings_tried,
                 n_leg_perms, leg_perm (the tuple of the match) and the
                 leg-count guard's verdict.

    Leg-count guard: the external blocks of A and B (len(vec) - nloop) must
    both equal n_ext; when A and B carry different leg counts no relabeling
    exists and None / [] is returned with details["leg_count_skip"] set --
    the permutation never indexes past a vector.
    """
    det = details if isinstance(details, dict) else {}
    det.update({"n_relabelings_tried": 0, "n_leg_perms": 0, "leg_perm": None,
                "leg_count_skip": None, "matches": []})
    empty = [] if find_all else None
    if len(vecsA) != len(vecsB):
        det["reason"] = f"propagator counts differ ({len(vecsA)} vs {len(vecsB)})"
        return empty
    n_extA = len(vecsA[0]) - nloop if vecsA else n_ext
    n_extB = len(vecsB[0]) - nloop if vecsB else n_ext
    if n_extA != n_extB:
        det["leg_count_skip"] = (f"external-leg counts differ ({n_extA} vs {n_extB}): "
                                 "no leg relabeling exists, comparison skipped")
        det["reason"] = det["leg_count_skip"]
        return empty
    if n_extA != n_ext:
        raise ValueError(f"n_ext={n_ext} does not match the vectors' external "
                         f"block ({n_extA} legs)")
    massA = nuA = massB = nuB = None
    if labels is not None:
        labA, labB = labels
        if len(labA) != len(vecsA) or len(labB) != len(vecsB):
            raise ValueError("labels are not aligned with the vectors "
                             f"({len(labA)}/{len(vecsA)} on A, {len(labB)}/{len(vecsB)} on B)")
        massA = [str(m) for m, _n in labA]
        nuA = [int(n) for _m, n in labA]
        massB = [str(m) for m, _n in labB]
        nuB = [int(n) for _m, n in labB]
        if sorted(massA) != sorted(massB):
            det["reason"] = ("mass multiset differs (%s vs %s): no label-preserving "
                             "relabeling exists" % (sorted(massA), sorted(massB)))
            return empty
        if sorted(nuA) != sorted(nuB):
            det["reason"] = ("multiplicity (propagator power) multiset differs (%s vs %s): "
                             "no label-preserving relabeling exists"
                             % (sorted(nuA), sorted(nuB)))
            return empty
    targetB = canonical_set(vecsB, massB, nuB)
    loop_perms = list(itertools.permutations(range(nloop)))
    sign_choices = list(itertools.product([1, -1], repeat=nloop))
    if try_leg_perms:
        leg_perms = list(itertools.permutations(range(n_ext)))
        if leg_labels is not None:
            labA, labB = leg_labels
            leg_perms = [lp for lp in leg_perms
                         if all(labA[j] == labB[lp[j]] for j in range(n_ext))]
    else:
        leg_perms = [tuple(range(n_ext))]
    det["n_leg_perms"] = len(leg_perms)
    found = []
    tried = 0
    for lp in loop_perms:
        for signs in sign_choices:
            for legp in leg_perms:
                tried += 1
                mapped = canonical_set(
                    [_apply_signed_perm_to_vec(v, lp, signs, legp, nloop, n_ext,
                                               dep_rows=dep_rows)
                     for v in vecsA], massA, nuA)
                if mapped == targetB:
                    desc = _describe_relabel(lp, signs, legp, nloop, n_ext)
                    found.append(desc)
                    det["matches"].append({"relabeling": desc, "loop_perm": tuple(lp),
                                           "signs": tuple(signs), "leg_perm": tuple(legp),
                                           "relabeling_class": "signed loop permutation"})
                    if det["leg_perm"] is None:
                        det["leg_perm"] = tuple(legp)
                        det["n_relabelings_tried"] = tried
                        det["relabeling_class"] = "signed loop permutation"
                    if not find_all:
                        return desc
    det["n_relabelings_tried"] = tried
    if affine:
        # the affine stage: the same leg permutations and labels, the shift
        # solved (never enumerated), the bound reported by name
        aff = _affine_search_vectors(
            vecsA, vecsB, nloop, n_ext, leg_perms, dep_rows=dep_rows,
            labelsA=(list(zip(massA, nuA)) if labels is not None else None),
            labelsB=(list(zip(massB, nuB)) if labels is not None else None),
            names=names, find_all=find_all)
        det["affine"] = aff
        for m in aff["matches"]:
            found.append(m["relabeling"])
            det["matches"].append({"relabeling": m["relabeling"],
                                   "loop_map_matrix": m["loop_map_matrix"],
                                   "shift_matrix": m["shift_matrix"],
                                   "leg_perm": m["leg_perm"], "leg_map": m["leg_map"],
                                   "relabeling_class": "affine unimodular loop redefinition"})
        if aff["found"]:
            if det["leg_perm"] is None:
                det["leg_perm"] = tuple(aff["matches"][0]["leg_perm"])
                det["relabeling_class"] = "affine unimodular loop redefinition"
            if not find_all:
                return aff["matches"][0]["relabeling"]
    if find_all:
        return found
    return None


def _describe_relabel(loop_perm, signs, leg_perm, nloop, n_ext):
    parts = []
    for i in range(nloop):
        sgn = "" if signs[i] == 1 else "-"
        if loop_perm[i] != i or signs[i] != 1:
            parts.append(f"k{i+1} -> {sgn}k{loop_perm[i]+1}")
    legparts = [f"p{j+1} -> p{leg_perm[j]+1}" for j in range(n_ext)
                if leg_perm[j] != j]
    desc = ", ".join(parts) if parts else "identity"
    if legparts:
        desc += "; " + ", ".join(legparts)
    return desc


# ---------------------------------------------------------------------------
#  2b. THE AFFINE CLASS  k -> U k + c.p  (loop-momentum basis changes) and the
#      ROUTING-INDEPENDENT identity test on the realized graph
# ---------------------------------------------------------------------------
#  Two routings of one graph differ by a change of loop-momentum basis
#  (k1 -> k2 - k1, a shift k -> k + p_i) that a signed permutation never
#  reaches: the propagator SET comparison above then says "different" for
#  the same graph.  The identity census met this on its double boxes (rows
#  1, 2, 7: a spanning-tree routing of the drawing vs the record's textbook
#  routing), on the sunrises, the kite and the bananas, and on the zgamma
#  box (row 20).  The two routes that resolve it:
#    * the AFFINE search (famlib.iso_search): loops k_i(A) -> sum_m U_im k_m(B)
#      + sum_j c_ij p_j(B) with U in unimodular_maps(L) -- every integer LxL
#      matrix with entries in {-1,0,1} and det +-1 for L <= 3, the signed
#      permutations composed with the two-vertex line reflections for L >= 4
#      (the census's class) -- legs permuted (label-preserving), masses and
#      multiplicities carried as labels.  The shift c is SOLVED, never
#      enumerated: for a given U the loop block of every mapped line is fixed
#      (U alone), so U survives only when the multiset of (loop block, mass,
#      nu) equals the target's; then L lines with independent loop blocks are
#      matched to their candidate target lines and the linear system a.c =
#      s b_ext - e determines c, and the full labeled set is verified.  The
#      bound is printed with every result (matrices, leg permutations,
#      matrices surviving the loop-block filter, shift systems solved).
#    * the REALIZED-GRAPH isomorphism (rows_01_02_04_05_07 graph_iso,
#      census_lib_r21 graph_isomorphism, famlib.nx_iso): VF2 on the realized
#      labeled multigraphs, three ways -- bare / with (mass, multiplicity) on
#      the lines / with the leg virtuality classes on the external nodes as
#      well -- plus a strict by-leg-name level; this is routing-independent
#      by construction and is the identity test to trust; the affine map is
#      its DESCRIBER (a map may not exist in the searched class where VF2
#      says isomorphic: then the isomorphism stands and the describer says
#      "none in the class").
#    * a CANONICAL FORM of the realized labeled graph (an individualization-
#      refinement canonical labeling, exact, no hash collisions by
#      construction): two routings of one graph give one hash.  The fingerprint
#      prints it beside the routed canonical_hash under its own name.

def unimodular_maps(L, entries=(-1, 0, 1)):
    """
    The census's class of loop-momentum basis changes (famlib.unimodular_maps):
    for L <= 3 every integer LxL matrix with entries in `entries` and det +-1
    (L = 1: 2, L = 2: 40, L = 3: 6960 matrices); for L >= 4 the products of the
    signed permutation matrices with the L+1 'line reflections' R_i (row i =
    (-1, ..., -1): the map that exchanges the role of loop line i with the
    (L+1)-th line of a two-vertex graph), in both orders, de-duplicated -- the
    routing changes that occur for the banana graphs.  Returns (matrices as
    lists of rows, class description).  The matrices are generated in a fixed
    order (itertools.product over the flat entries; signed permutations in
    itertools order times the reflections), so "the first map found" is
    reproducible.
    """
    mats = []
    if L <= 3:
        for flat in itertools.product(entries, repeat=L * L):
            M = [list(flat[i * L:(i + 1) * L]) for i in range(L)]
            if _int_det(M) in (1, -1):
                mats.append(M)
        kind = "all integer matrices with entries in {-1,0,1} and det = +-1"
    else:
        sps = []
        for sigma in itertools.permutations(range(L)):
            for signs in itertools.product([1, -1], repeat=L):
                M = [[0] * L for _ in range(L)]
                for i in range(L):
                    M[i][sigma[i]] = signs[i]
                sps.append(M)
        refl = [[[1 if i == j else 0 for j in range(L)] for i in range(L)]]
        for i in range(L):
            Rm = [[1 if r == cc else 0 for cc in range(L)] for r in range(L)]
            for m in range(L):
                Rm[i][m] = -1
            refl.append(Rm)
        seen = set()
        for S in sps:
            for Rm in refl:
                for M in (_int_matmul(S, Rm), _int_matmul(Rm, S)):
                    key = tuple(x for row in M for x in row)
                    if key not in seen:
                        seen.add(key)
                        mats.append(M)
        kind = "signed permutations composed with the two-vertex line reflections R_i (L >= 4)"
    return mats, kind


def _int_det(M):
    n = len(M)
    if n == 0:
        return 1
    if n == 1:
        return M[0][0]
    if n == 2:
        return M[0][0] * M[1][1] - M[0][1] * M[1][0]
    s = 0
    for j in range(n):
        if M[0][j] == 0:
            continue
        minor = [row[:j] + row[j + 1:] for row in M[1:]]
        s += (-1) ** j * M[0][j] * _int_det(minor)
    return s


def _int_matmul(A, B):
    n, m, p = len(A), len(B), len(B[0])
    return [[sum(A[i][k] * B[k][j] for k in range(m)) for j in range(p)] for i in range(n)]


def _solve_rational(M, R):
    """Solve M c = R exactly (Fractions): M is n x n (rows = the loop blocks of
    the basis lines), R is n x w; returns c (n x w) or None when M is
    singular."""
    from fractions import Fraction
    n = len(M)
    w = len(R[0]) if R else 0
    A = [[Fraction(x) for x in M[i]] + [Fraction(x) for x in R[i]] for i in range(n)]
    for col in range(n):
        piv = None
        for r in range(col, n):
            if A[r][col] != 0:
                piv = r
                break
        if piv is None:
            return None
        A[col], A[piv] = A[piv], A[col]
        pv = A[col][col]
        A[col] = [x / pv for x in A[col]]
        for r in range(n):
            if r != col and A[r][col] != 0:
                f = A[r][col]
                A[r] = [x - f * y for x, y in zip(A[r], A[col])]
    return [row[n:n + w] for row in A]


def _fr(x):
    from fractions import Fraction
    q = sp.Rational(x)
    return Fraction(int(q.p), int(q.q))


def _canon_fr(v):
    for c in v:
        if c != 0:
            if c < 0:
                return tuple(-x for x in v)
            break
    return tuple(v)


def _describe_affine(U, c, loops_a, loops_b, exts_b):
    """The census's describer strings, one per loop of A: 'l1 -> k1 -k2 +p3'
    (famlib.iso_search): the U row's +-loop names of B, then the shift's
    external momenta with their coefficients (+p3, -p1, +2*p2)."""
    out = []
    for i, la in enumerate(loops_a):
        terms = []
        for m, lb in enumerate(loops_b):
            if U[i][m]:
                terms.append(("-" if U[i][m] < 0 else "+") + lb)
        for j, e in enumerate(exts_b):
            x = c[i][j]
            if x:
                if x in (1, -1):
                    terms.append(("-" if x < 0 else "+") + e)
                else:
                    terms.append(("-" if x < 0 else "+") + str(abs(x)) + "*" + e)
        out.append(f"{la} -> " + (" ".join(terms).lstrip("+") if terms else "0"))
    return out


def _affine_search_vectors(vecsA, vecsB, nloop, n_ext, leg_perms, dep_rows=None,
                           labelsA=None, labelsB=None, names=None, find_all=False,
                           umaps=None):
    """
    The affine search on momentum vectors (the primitive isomorphism_to calls
    with affine=True and affine_iso_search wraps for two families):
    find (leg permutation, U, c) with  k_i(A) -> sum_m U_im k_m(B) + sum_j
    c_ij p_j(B)  and  p_j(A) -> p_legperm(j)(B)  (B's momentum conservation
    re-imposed through dep_rows) mapping the LABELED canonical set of A
    ((momentum, mass, nu) Counter; labels None = bare momenta) onto B's.

    leg_perms: the leg permutations to search (already label-restricted by
    the caller).  names: (loops_A, loops_B, exts_B) for the describer
    strings (default k1.., k1.., p1..).  umaps: (matrices, class) to search;
    default unimodular_maps(nloop).

    Returns a dict: found, matches (list of {loop_relabeling, loop_map_matrix,
    shift_matrix, leg_perm, leg_map, relabeling}), the first match's fields
    at top level (loop_relabeling, loop_map_matrix, shift_matrix, leg_map,
    relabeling), loop_map_class, n_unimodular, n_leg_perms,
    n_U_surviving_loop_block_filter, n_shift_systems_solved, tried
    (= verifications of the full set), bound (a string), shift_in_census_range
    (every |c_ij| <= 1, the census's enumerated range) and reason when nothing
    was found.
    """
    import collections
    from fractions import Fraction
    if names:
        loops_a, loops_b, exts_b = names[0], names[1], names[2]
        exts_a = names[3] if len(names) > 3 else exts_b
    else:
        loops_a = loops_b = [f"k{i+1}" for i in range(nloop)]
        exts_a = exts_b = [f"p{j+1}" for j in range(n_ext)]
    indep = [j for j in range(n_ext) if j not in (dep_rows or {})]
    mats, kind = umaps if umaps is not None else unimodular_maps(nloop)
    out = {"found": False, "matches": [], "loop_map_class": kind,
           "n_unimodular": len(mats), "n_leg_perms": len(leg_perms),
           "n_U_surviving_loop_block_filter": 0, "n_shift_systems_solved": 0,
           "tried": 0, "reason": None,
           "loop_relabeling": None, "loop_map_matrix": None, "shift_matrix": None,
           "shift_legs": None, "leg_map": None, "relabeling": None,
           "shift_in_census_range": None}
    if len(vecsA) != len(vecsB):
        out["reason"] = f"propagator counts differ ({len(vecsA)} vs {len(vecsB)})"
        out["bound"] = "0 (line counts differ)"
        return out
    if not vecsA:
        out["reason"] = "no lines"
        out["bound"] = "0"
        return out
    labA = [(str(m), int(n)) for m, n in labelsA] if labelsA is not None else [("", 1)] * len(vecsA)
    labB = [(str(m), int(n)) for m, n in labelsB] if labelsB is not None else [("", 1)] * len(vecsB)
    if sorted(labA) != sorted(labB):
        out["reason"] = ("label multisets differ (%s vs %s): no label-preserving map exists"
                         % (sorted(labA), sorted(labB)))
        out["bound"] = "0 (label multisets differ)"
        return out
    A = [tuple(_fr(c) for c in v) for v in vecsA]
    B = [tuple(_fr(c) for c in v) for v in vecsB]
    targetB = collections.Counter((_canon_fr(b), lab) for b, lab in zip(B, labB))
    B_loop_counter = collections.Counter((_canon_fr(b[:nloop]), lab) for b, lab in zip(B, labB))
    B_by_loop = collections.defaultdict(list)
    for b, lab in zip(B, labB):
        B_by_loop[(_canon_fr(b[:nloop]), lab)].append(_canon_fr(b))
    # (i) the loop-block filter: U alone fixes every mapped loop block (integer
    #     arithmetic when the loop blocks are integral, which they are for
    #     every record read here; a Fraction and its integer value hash alike)
    aL = [tuple((int(x) if x.denominator == 1 else x) for x in a[:nloop]) for a in A]
    rng = range(nloop)
    surviving = []
    for U in mats:
        cols = [[U[i][m] for i in rng] for m in rng]
        keys = []
        for a, lab in zip(aL, labA):
            w = tuple(sum(a[i] * col[i] for i in rng) for col in cols)
            keys.append((_canon_fr(w), lab))
        if collections.Counter(keys) == B_loop_counter:
            surviving.append((U, keys))
    out["n_U_surviving_loop_block_filter"] = len(surviving)
    ident = tuple(range(nloop))
    ones = tuple([1] * nloop)
    solved = 0
    tried = 0
    found = []

    def rank_rows(rows):
        return int(sp.Matrix([[sp.Rational(x.numerator, x.denominator) for x in r]
                              for r in rows]).rank()) if rows else 0

    for legp in leg_perms:
        # A's external coordinates in B's basis (B's conservation re-imposed)
        Ap = [tuple(_fr(c) for c in _apply_signed_perm_to_vec(v, ident, ones, legp, nloop, n_ext,
                                                              dep_rows=dep_rows))
              for v in vecsA]
        for U, keys in surviving:
            cands = [B_by_loop[k] for k in keys]
            # basis: lines with independent (original) loop blocks, fewest candidates first
            order = sorted(range(len(A)), key=lambda r: (len(cands[r]), r))
            basis = []
            for r in order:
                if rank_rows([A[q][:nloop] for q in basis] + [A[r][:nloop]]) > len(basis):
                    basis.append(r)
                if len(basis) == nloop:
                    break
            if len(basis) < nloop:
                continue
            M = [list(A[r][:nloop]) for r in basis]
            signs = []
            for r in basis:
                w = keys[r][0]
                w_raw = tuple(sum(A[r][i] * U[i][m] for i in range(nloop)) for m in range(nloop))
                signs.append(1 if w_raw == w else -1)
            for combo in itertools.product(*[cands[r] for r in basis]):
                R = []
                for s_r, b, r in zip(signs, combo, basis):
                    e = Ap[r][nloop:]
                    R.append([s_r * b[nloop + j] - e[j] for j in range(n_ext)])
                c = _solve_rational(M, R)
                solved += 1
                if c is None:
                    continue
                tried += 1
                mapped = collections.Counter()
                for a, ap, lab in zip(A, Ap, labA):
                    w = [sum(a[i] * U[i][m] for i in range(nloop)) for m in range(nloop)]
                    e2 = [ap[nloop + j] + sum(a[i] * c[i][j] for i in range(nloop))
                          for j in range(n_ext)]
                    mapped[(_canon_fr(tuple(w + e2)), lab)] += 1
                if mapped == targetB:
                    cc = [[(int(x) if x.denominator == 1 else x) for x in row] for row in c]
                    desc = _describe_affine(U, cc, loops_a, loops_b, exts_b)
                    leg_map = {(f"p{j+1}" if len(exts_a) <= j else exts_a[j]):
                               (f"p{legp[j]+1}" if len(exts_b) <= legp[j] else exts_b[legp[j]])
                               for j in range(n_ext)}
                    legparts = [f"{a} -> {b}" for a, b in leg_map.items() if a != b]
                    rel = ", ".join(desc) + (("; " + ", ".join(legparts)) if legparts else "")
                    # the shift over B's INDEPENDENT legs (the dependent slot is
                    # always 0 after conservation), the census's S shape
                    m = {"loop_relabeling": desc, "loop_map_matrix": [list(r) for r in U],
                         "shift_matrix": [[(cc[i][j] if isinstance(cc[i][j], int) else str(cc[i][j]))
                                           for j in indep] for i in range(nloop)],
                         "shift_legs": [(f"p{j+1}" if len(exts_b) <= j else exts_b[j]) for j in indep],
                         "leg_perm": tuple(legp), "leg_map": leg_map, "relabeling": rel,
                         "shift_in_census_range": all(isinstance(x, int) and abs(x) <= 1
                                                      for row in cc for x in row)}
                    found.append(m)
                    if not find_all:
                        break
            if found and not find_all:
                break
        if found and not find_all:
            break
    out.update({"tried": tried, "n_shift_systems_solved": solved, "matches": found,
                "found": bool(found)})
    out["bound"] = ("%d unimodular matrices (%s) x %d leg permutation(s); %d matrices survive the "
                    "loop-block filter; %d shift systems solved, %d full-set verifications"
                    % (len(mats), kind, len(leg_perms), len(surviving), solved, tried))
    if found:
        for k in ("loop_relabeling", "loop_map_matrix", "shift_matrix", "shift_legs", "leg_map",
                  "relabeling", "shift_in_census_range"):
            out[k] = found[0][k]
    else:
        out["reason"] = ("no map in the class (%s) with any leg permutation of the %d searched "
                         "carries the labeled set: %d matrices survive the loop-block filter, "
                         "%d shift systems solved, %d verifications"
                         % (kind, len(leg_perms), len(surviving), solved, tried))
    return out


def affine_iso_search(famA, famB, mass_mapA=None, mass_mapB=None, leg_labelsA=None,
                      leg_labelsB=None, try_leg_perms=True, respect_leg_classes=True,
                      with_mass=True, find_all=False):
    """
    The census's iso_search on two families (famlib L473): a relabeling
    p_j(A) -> p_pi(j)(B) (a bijection; the leg virtuality classes respected
    when respect_leg_classes and both sides declare them) and k_i(A) ->
    sum_m U_im k_m(B) + sum_j c_ij p_j(B) (U in unimodular_maps(L), c solved)
    with masses and multiplicities carried as labels (with_mass; False = the
    bare momentum sets).  mass_mapA / mass_mapB and leg_labelsA / leg_labelsB
    as in labeled_isomorphism.  Returns the census's dict shape: found,
    loop_relabeling (['l1 -> k1 -k2 +p3', ...]), loop_map_matrix, shift_matrix,
    leg_map, tried, with_mass, respect_virtuality, n_leg_perms_admissible,
    loop_map_class, plus the bound fields of _affine_search_vectors and
    "signed_permutation_match" (the tool's former notion on the same labeled
    set, for the record).
    """
    LA = labeled_set(famA, mass_mapA)
    LB = labeled_set(famB, mass_mapB)
    nloop, n_ext = len(famA.loops), len(famA.exts)
    out = {"found": False, "with_mass": bool(with_mass),
           "respect_virtuality": bool(respect_leg_classes),
           "n_lines": (len(LA["vecs"]), len(LB["vecs"])),
           "n_loops": (nloop, len(famB.loops)), "n_legs": (n_ext, len(famB.exts)),
           "mass_multiset_A": LA["mass_multiset"], "mass_multiset_B": LB["mass_multiset"],
           "nu_multiset_A": LA["nu_multiset"], "nu_multiset_B": LB["nu_multiset"]}
    if len(famB.loops) != nloop or len(famB.exts) != n_ext or len(LA["vecs"]) != len(LB["vecs"]):
        out["reason"] = ("line/loop/leg counts differ: A %d/%d/%d vs B %d/%d/%d"
                         % (len(LA["vecs"]), nloop, n_ext, len(LB["vecs"]), len(famB.loops),
                            len(famB.exts)))
        out["n_leg_perms_admissible"] = 0
        out["tried"] = 0
        return out
    cA = ([_leg_class(leg_labelsA.get(e)) for e in famA.exts] if leg_labelsA is not None
          else LA["leg_classes"])
    cB = ([_leg_class(leg_labelsB.get(e)) for e in famB.exts] if leg_labelsB is not None
          else LB["leg_classes"])
    labels_ok = (cA is not None and cB is not None and all(c is not None for c in cA)
                 and all(c is not None for c in cB))
    out["leg_classes_A"], out["leg_classes_B"] = cA, cB
    out["leg_labels_compared"] = bool(labels_ok and respect_leg_classes)
    if try_leg_perms:
        leg_perms = list(itertools.permutations(range(n_ext)))
        if respect_leg_classes and labels_ok:
            leg_perms = [lp for lp in leg_perms if all(cA[j] == cB[lp[j]] for j in range(n_ext))]
    else:
        leg_perms = [tuple(range(n_ext))]
    out["n_leg_perms_admissible"] = len(leg_perms)
    labels = ((list(zip(LA["masses"], LA["nu"])), list(zip(LB["masses"], LB["nu"])))
              if with_mass else (None, None))
    dep = _dependent_leg_rows(famB)
    det = {}
    signed = isomorphism_to(LA["vecs"], LB["vecs"], nloop, n_ext, try_leg_perms=try_leg_perms,
                            dep_rows=dep, leg_labels=((cA, cB) if (labels_ok and respect_leg_classes
                                                                    and try_leg_perms) else None),
                            labels=(labels if with_mass else None), details=det)
    out["signed_permutation_match"] = signed
    out["n_signed_relabelings_tried"] = det.get("n_relabelings_tried")
    res = _affine_search_vectors(LA["vecs"], LB["vecs"], nloop, n_ext, leg_perms, dep_rows=dep,
                                 labelsA=labels[0], labelsB=labels[1],
                                 names=(list(famA.loops), list(famB.loops), list(famB.exts),
                                        list(famA.exts)),
                                 find_all=find_all)
    out.update(res)
    return out


# ---- the realized labeled graph, its isomorphism (VF2) and its canonical form

def _realized_labeled_graphs(fam, mass_map=None, leg_labels=None):
    """Every distinct leg-labelled realization of the family (build_graph's
    first plus its alternatives) as an annotated networkx MultiGraph: nodes
    kind "int" / "ext" (ext: leg = the leg name, cls = its virtuality class
    or None), internal edges kind "line" with mass (label string, the record's
    numeric values and mass_map applied) and nu (multiplicity), leg edges
    kind "leg".  Returns (graphs, report) -- graphs empty when the family
    does not realize (report["cause"] names why)."""
    if not _HAVE_NX:
        return [], {"verdict": "UNAVAILABLE", "cause": "networkx unavailable"}
    G, ok = build_graph(fam)
    rep = G.graph.get("realization", {}) if G is not None else {"cause": "no graph"}
    if not ok or G is None:
        return [], rep
    mv = momentum_vectors(fam)
    vecs, masses, _isps = mv
    labels = [_mass_label(m, fam, mass_map) for m in masses]
    nus = list(mv.nu)
    virt = dict(fam.leg_virt)
    if leg_labels is not None:
        virt.update({str(k): v for k, v in leg_labels.items()})
    graphs = []
    for src in [G] + list(G.graph.get("alternative_realizations", []) or []):
        H = nx.MultiGraph()
        for n in src.nodes():
            ext = str(n).startswith("ext_")
            if ext:
                leg = str(n)[4:]
                H.add_node(n, kind="ext", leg=leg, cls=_leg_class(virt.get(leg)))
            else:
                H.add_node(n, kind="int", leg="", cls="")
        for u, v, k in src.edges(keys=True):
            k = str(k)
            if k.startswith("ext_"):
                H.add_edge(u, v, key=k, kind="leg", mass="leg", nu=0)
            else:
                i = int(k[1:])
                H.add_edge(u, v, key=k, kind="line", mass=labels[i], nu=nus[i])
        graphs.append(H)
    return graphs, rep


def realized_graph_iso(famA, famB, mass_mapA=None, mass_mapB=None, leg_labelsA=None,
                       leg_labelsB=None):
    """
    Graph isomorphism of the two families' REALIZED Feynman graphs (networkx
    VF2 on the annotated multigraphs), the census's routing-independent
    identity test (rows_01_02_04_05_07 graph_iso L502; census_lib_r21
    graph_isomorphism L335; famlib.nx_iso L564), four ways:
      "bare"                     internal / external nodes, lines / legs
      "with_masses"              + (mass label, multiplicity) on every line
      "with_masses_and_legs"     + the leg virtuality class on every external
                                 node (None with a reason when either side
                                 declares no class for some leg -- nothing is
                                 inferred)
      "with_masses_and_leg_names"  + the external nodes matched by leg NAME
                                 (census_lib_r21 strict_legs)
    each True / False / None with "<level>_leg_map" ({A leg: B leg} of the
    VF2 mapping) when True.  Every distinct leg-labelled realization of each
    side is tried (a symmetric insertion gives several); "realization_pair"
    names the (A index, B index) that matched, "n_realizations" both counts.
    V_A / E_A / V_B / E_B are the annotated graphs' node and edge counts
    (external nodes and leg edges included, as the census counted: the row-4
    'V,E 10,12 vs 12,14'); V_internal_* / E_lines_* the internal ones.
    "error" names a side that does not realize (every level None).
    mass_mapA / mass_mapB: {symbol: value} applied to the line masses after
    the record's own numeric values; leg_labelsA / leg_labelsB: {leg:
    virtuality} overriding Family.leg_virt.
    """
    out = {"bare": None, "with_masses": None, "with_masses_and_legs": None,
           "with_masses_and_leg_names": None, "error": None, "realization_pair": None,
           "n_realizations": None, "V_A": None, "E_A": None, "V_B": None, "E_B": None,
           "V_internal_A": None, "E_lines_A": None, "V_internal_B": None, "E_lines_B": None}
    if not _HAVE_NX:
        out["error"] = "networkx unavailable: no realized-graph isomorphism"
        return out
    GA, repA = _realized_labeled_graphs(famA, mass_mapA, leg_labelsA)
    GB, repB = _realized_labeled_graphs(famB, mass_mapB, leg_labelsB)
    if not GA or not GB:
        side = "A" if not GA else "B"
        out["error"] = ("vertex realization failed for %s: %s"
                        % (side, (repA if not GA else repB).get("cause")))
        return out
    from networkx.algorithms.isomorphism import MultiGraphMatcher

    def counts(G):
        Vi = sum(1 for n in G.nodes if G.nodes[n]["kind"] == "int")
        El = sum(1 for _u, _v, d in G.edges(data=True) if d["kind"] == "line")
        return G.number_of_nodes(), G.number_of_edges(), Vi, El
    out["V_A"], out["E_A"], out["V_internal_A"], out["E_lines_A"] = counts(GA[0])
    out["V_B"], out["E_B"], out["V_internal_B"], out["E_lines_B"] = counts(GB[0])
    out["n_realizations"] = (len(GA), len(GB))
    import collections

    def em_kind(ea, eb):
        return (collections.Counter(d["kind"] for d in ea.values())
                == collections.Counter(d["kind"] for d in eb.values()))

    def em_mass(ea, eb):
        return (collections.Counter((d["kind"], d["mass"], d["nu"]) for d in ea.values())
                == collections.Counter((d["kind"], d["mass"], d["nu"]) for d in eb.values()))

    def nm_kind(a, b):
        return a["kind"] == b["kind"]

    def nm_cls(a, b):
        return a["kind"] == b["kind"] and a["cls"] == b["cls"]

    def nm_name(a, b):
        return a["kind"] == b["kind"] and a["leg"] == b["leg"]
    legs_declared = (all(GA[0].nodes[n]["cls"] is not None for n in GA[0].nodes
                         if GA[0].nodes[n]["kind"] == "ext")
                     and all(GB[0].nodes[n]["cls"] is not None for n in GB[0].nodes
                             if GB[0].nodes[n]["kind"] == "ext"))
    levels = [("bare", nm_kind, em_kind), ("with_masses", nm_kind, em_mass),
              ("with_masses_and_legs", nm_cls, em_mass),
              ("with_masses_and_leg_names", nm_name, em_mass)]
    for name, nm, em in levels:
        if name == "with_masses_and_legs" and not legs_declared:
            out[name] = None
            out[name + "_reason"] = ("leg virtuality classes not declared for every leg on both "
                                     "sides; nothing inferred")
            continue
        hit = None
        for ia, HA in enumerate(GA):
            for ib, HB in enumerate(GB):
                gm = MultiGraphMatcher(HA, HB, node_match=nm, edge_match=em)
                if gm.is_isomorphic():
                    hit = (ia, ib, gm.mapping)
                    break
            if hit:
                break
        out[name] = hit is not None
        if hit:
            ia, ib, mapping = hit
            out[name + "_leg_map"] = {str(a)[4:]: str(b)[4:] for a, b in mapping.items()
                                      if str(a).startswith("ext_")}
            out[name + "_realization_pair"] = (ia, ib)
            if out["realization_pair"] is None:
                out["realization_pair"] = (ia, ib)
    return out


def _canonical_form(node_labels, edges, max_leaves=20000):
    """
    Exact canonical form of a labeled multigraph by individualization-
    refinement (the search tree is built by isomorphism-invariant rules and
    every leaf is visited, so the minimum leaf encoding is a canonical form;
    no hash-collision class as a Weisfeiler-Lehman hash would have).
    node_labels: {node: label string}; edges: [(u, v, label string)].
    Returns (encoding, n_leaves, capped): encoding = (node labels in canonical
    order, sorted (i, j, label) edge list) or None when the leaf count exceeds
    max_leaves (capped True).
    """
    nodes = sorted(node_labels, key=lambda n: (str(node_labels[n]), str(n)))
    idx = {n: i for i, n in enumerate(nodes)}
    n = len(nodes)
    inc = [[] for _ in range(n)]
    elist = []
    for u, v, lab in edges:
        iu, iv = idx[u], idx[v]
        elist.append((iu, iv, str(lab)))
        inc[iu].append((str(lab), iv))
        inc[iv].append((str(lab), iu))
    labs = [str(node_labels[x]) for x in nodes]

    def canon_colors(pre):
        distinct = sorted(set(pre))
        rank = {p: i for i, p in enumerate(distinct)}
        return [rank[p] for p in pre]

    def refine(col):
        col = list(col)
        while True:
            sig = [(col[i], tuple(sorted((lab, col[j]) for lab, j in inc[i]))) for i in range(n)]
            new = canon_colors(sig)
            if len(set(new)) == len(set(col)):
                return new
            col = new
    state = {"leaves": 0, "best": None, "capped": False}

    def leaf(col):
        order = sorted(range(n), key=lambda i: col[i])
        pos = {v: k for k, v in enumerate(order)}
        enc = (tuple(labs[v] for v in order),
               tuple(sorted((min(pos[u], pos[v]), max(pos[u], pos[v]), lab) for u, v, lab in elist)))
        state["leaves"] += 1
        if state["best"] is None or enc < state["best"]:
            state["best"] = enc

    def search(col):
        if state["capped"]:
            return
        if len(set(col)) == n:
            leaf(col)
            if state["leaves"] > max_leaves:
                state["capped"] = True
            return
        cells = {}
        for i, c in enumerate(col):
            cells.setdefault(c, []).append(i)
        target = min(c for c, members in cells.items() if len(members) > 1)
        for v in cells[target]:
            pre = [(col[i], 1 if i == v else 0) for i in range(n)]
            search(refine(canon_colors(pre)))
            if state["capped"]:
                return
    search(refine(canon_colors(labs)))
    if state["capped"]:
        return None, state["leaves"], True
    return state["best"], state["leaves"], False


def realized_canonical_hash(fam, with_leg_classes=False, mass_map=None, leg_labels=None):
    """
    The routing-independent canonical hash of the family's realized labeled
    graph: sha1[:12] of _canonical_form on the annotated realization (nodes
    "int" / "ext"; lines "line:m=<mass label>;nu=<n>"; legs "leg").  With
    with_leg_classes the external nodes carry their virtuality class
    ("ext:0" / "ext:offshell"); when any leg has no declared class the hash is
    None with the reason (nothing inferred).  Two routings of one graph (a
    spanning-tree routing of a drawing, a record's textbook routing) give one
    hash; the leg NAMES do not enter (a drawing's p3 may be the record's p1).
    Returns a dict: hash, method, reason, n_leaves, V, E (annotated graph
    counts), alternatives (the hashes of the other distinct leg-labelled
    realizations when they differ from the first: the family's realization is
    then ambiguous and the fingerprint says so).
    """
    method = ("individualization-refinement canonical form of the realized labeled multigraph "
              "(exact; nodes int/ext%s, lines (mass, multiplicity), legs unlabeled), sha1[:12]"
              % (" with the leg virtuality class" if with_leg_classes else ""))
    out = {"hash": None, "method": method, "reason": None, "n_leaves": None,
           "V": None, "E": None, "alternatives": []}
    graphs, rep = _realized_labeled_graphs(fam, mass_map, leg_labels)
    if not graphs:
        out["reason"] = "realization failed: %s" % rep.get("cause")
        return out
    hashes = []
    import hashlib
    for H in graphs:
        nl = {}
        for x in H.nodes:
            d = H.nodes[x]
            if d["kind"] == "ext":
                if with_leg_classes:
                    if d["cls"] is None:
                        out["reason"] = ("leg virtuality class not declared for leg %s; "
                                         "nothing inferred" % d["leg"])
                        return out
                    nl[x] = "ext:" + d["cls"]
                else:
                    nl[x] = "ext"
            else:
                nl[x] = "int"
        el = []
        for u, v, d in H.edges(data=True):
            el.append((u, v, "leg" if d["kind"] == "leg" else "line:m=%s;nu=%d" % (d["mass"], d["nu"])))
        enc, n_leaves, capped = _canonical_form(nl, el)
        if capped:
            out["reason"] = "canonical form search capped at %d leaves" % n_leaves
            return out
        h = hashlib.sha1(json.dumps(enc).encode()).hexdigest()[:12]
        hashes.append((h, n_leaves, H.number_of_nodes(), H.number_of_edges()))
    out["hash"], out["n_leaves"], out["V"], out["E"] = hashes[0]
    out["alternatives"] = sorted({h for h, _n, _v, _e in hashes[1:] if h != hashes[0][0]})
    return out


# ---------------------------------------------------------------------------
#  3.  GRAPH RECONSTRUCTION  +  PLANARITY  +  CUT SIGNATURE
# ---------------------------------------------------------------------------
#  We build the Feynman graph from momentum-flow incidence.  Each internal edge
#  carries a momentum q_e (the propagator vector).  External legs carry p_j.
#  We reconstruct vertices by momentum conservation: this is the standard dual
#  of "the loop momenta are independent cycles".  Concretely we use the
#  incidence-from-routing construction validated for these families.

def build_graph(fam):
    """
    Reconstruct the Feynman graph from the physical-propagator momenta by
    momentum-conservation vertex realization, and attach the external legs.

    Returns (G, ok) where G is a networkx.MultiGraph with internal edges keyed
    "e<i>" (attr id="e<i>") and external legs keyed "ext_<leg>" joining a
    degree-1 external node "ext_<leg>" to its attachment vertex.  ok is False
    when no realization exists; G is then an EMPTY MultiGraph.  In both cases
    G.graph["realization"] carries the realization report (non_graph_test +
    connected_covers): the count test 3V <= 2E + N with V = E - L + 1, the
    loop-momentum rank, repeated line momenta, the number of exact covers of
    the half-edges with V = E - L + 1, how many of them are connected, how many
    distinct connected graphs survive (ambiguity), legs per vertex, and a named
    "cause" when ok is False.

    METHOD.  Every internal edge contributes two half-edges (+q_e at its head,
    -q_e at its tail) and every external leg one inward half-edge (+p_j, with
    the dependent leg resolved through ext_subs).  A vertex is a MINIMAL
    zero-sum set of half-edges of size 3..6 from distinct edges, not all
    external.  A realization is an exact cover of the half-edges by such
    vertices with exactly V = E - L + 1 vertices whose graph is CONNECTED (one
    independent cycle per loop momentum).  All covers are enumerated (bounded
    by a step budget and a cover cap, both reported); the first zero-sum
    partition alone is NOT accepted: for a two-particle-reducible insertion it
    closes the insertion on itself as a detached vacuum piece, and for a
    non-textbook routing it merges vertices (three lines meeting at a vertex
    have far ends that also sum to zero) and under-counts V.
    """
    return _realize_graph(fam)


def _realize_graph(fam, step_budget=None, max_covers=None):
    """Momentum-conservation vertex realization (see _exact_vertex_realizer)."""
    if not _HAVE_NX:
        return None, False
    vecs, masses, isps = momentum_vectors(fam)
    nloop = len(fam.loops)
    n_ext = len(fam.exts)
    syms = _symbols(fam.loops, fam.exts)
    ext_resolved = _resolve_ext_subs(fam.exts, fam.ext_subs, syms)
    kw = {"masses": masses}
    if step_budget is not None:
        kw["step_budget"] = step_budget
    if max_covers is not None:
        kw["max_covers"] = max_covers
    return _exact_vertex_realizer(vecs, fam, nloop, n_ext, ext_resolved, syms, **kw)


#  Step cap of the exact-cover search (Algorithm X nodes).  The cap the first
#  realizer carried (400000) is kept: on the largest family this tool has been
#  timed on (10 lines, 4 legs, 3 loops: 24 half-edges) the full enumeration
#  finishes far below it; a family that exhausts it is reported by name
#  ("budget_exhausted": True, cause "step budget exhausted") and never given a
#  planarity or cut verdict.  REALIZER_MAX_COVERS bounds the number of exact
#  covers collected for the ambiguity count ("enumeration_capped": True).
REALIZER_STEP_BUDGET = 400000
REALIZER_MAX_COVERS = 200
REALIZER_MAX_VERTEX_DEGREE = 6


def non_graph_test(vecs, fam):
    """
    Arithmetic non-graph tests on the top-sector lines (no search):

      * count test: with V = E - L + 1 (a connected graph with one independent
        cycle per loop momentum) and every vertex at least trivalent,
        3V <= 2E + N must hold (each internal line has two ends, each leg one);
      * loop rank: the loop-momentum block of the line momenta must have rank
        L, else the top sector does not depend on every loop momentum and
        V = E - L + 1 is miscounted;
      * distinct lines: two top-sector lines with the same momentum (up to sign)
        are a repeated line -- a dot belongs in the propagator power, a second
        line in series would need a 2-valent vertex, which no realization here
        admits;
      * leg sum: the external legs (each resolved through ext_subs) must sum
        to zero -- every half-edge is covered by a zero-sum vertex and the
        internal half-edges cancel pairwise, so no cover exists otherwise.  A
        failure here is a LEG-SET defect of the input (a leg unspelled by the
        propagators, no conservation rule), not a property of the graph: it is
        reported as NOT-CHECKABLE, never as NON-GRAPH.

    Returns a dict with E, L, N, V, three_V, two_E_plus_N, the four booleans,
    "passes" (all four) and "cause" (None, or the first failing test by name).
    """
    E = len(vecs)
    L = len(fam.loops)
    N = len(fam.exts)
    V = E - L + 1
    three_V = 3 * V
    two_E_plus_N = 2 * E + N
    count_ok = three_V <= two_E_plus_N
    syms = _symbols(fam.loops, fam.exts)
    ext_resolved = _resolve_ext_subs(fam.exts, fam.ext_subs, syms)
    leg_sum = sp.expand(sum((ext_resolved[syms[e]] for e in fam.exts), sp.Integer(0)))
    legs_ok = (leg_sum == 0)
    if vecs and L > 0:
        loop_rank = int(sp.Matrix([[c for c in v[:L]] for v in vecs]).rank())
    else:
        loop_rank = 0
    rank_ok = (loop_rank == L)
    seen = {}
    repeated = []
    for i, v in enumerate(vecs):
        key = _canon_vec(v)
        if key in seen:
            repeated.append([seen[key], i])
        else:
            seen[key] = i
    distinct_ok = not repeated
    cause = None
    if not legs_ok:
        cause = ("NOT-CHECKABLE: the external legs %s do not sum to zero (sum = %s; "
                 "ext_subs = %s): a leg is missing or no momentum-conservation rule "
                 "was given, so no vertex cover can close -- fix the leg set before "
                 "reading any graph verdict" % (fam.exts, leg_sum, dict(fam.ext_subs)))
    elif not count_ok:
        cause = ("NON-GRAPH: 3V = %d > 2E + N = %d with V = E - L + 1 = %d "
                 "(E = %d lines, L = %d loops, N = %d legs)"
                 % (three_V, two_E_plus_N, V, E, L, N))
    elif not rank_ok:
        cause = ("NON-GRAPH: the top-sector lines span only %d of the %d "
                 "declared loop momenta (loop rank %d < L = %d)"
                 % (loop_rank, L, loop_rank, L))
    elif not distinct_ok:
        cause = ("NON-GRAPH: repeated line momentum in the top sector (line "
                 "pairs %s carry the same momentum up to sign); a dot is a "
                 "propagator power, not a second line" % repeated)
    return {"E": E, "L": L, "N": N, "V": V,
            "three_V": three_V, "two_E_plus_N": two_E_plus_N,
            "count_test_3V_le_2E_plus_N": bool(count_ok),
            "loop_momentum_rank": loop_rank, "loop_rank_test": bool(rank_ok),
            "repeated_line_pairs": repeated, "distinct_lines_test": bool(distinct_ok),
            "leg_sum": str(leg_sum), "legs_sum_to_zero_test": bool(legs_ok),
            "passes": bool(legs_ok and count_ok and rank_ok and distinct_ok),
            "cause": cause}


def connected_covers(he, n_loops, target_V=None, step_budget=REALIZER_STEP_BUDGET,
                     max_covers=REALIZER_MAX_COVERS, maxdeg=REALIZER_MAX_VERTEX_DEGREE,
                     masses=None):
    """
    Enumerate the exact covers of the half-edges by minimal zero-sum vertices
    and keep the CONNECTED ones with E - V + 1 = n_loops.

    he: list of (momentum_vector, (edge_id, end)) -- edge_id "e<i>" for an
        internal line (ends "head"/"tail") or "ext_<leg>" (end "in").
    target_V: when given, only covers with exactly that many vertices are
        collected (the search is pruned to them); None collects every cover
        (used to report what exists when the count test already failed).
    masses: per-line mass^2 (indexed like the lines), stored on the internal
        edges as attribute "mass" ("0" for massless) and used as a label in
        the distinctness test.

    Three notions of "distinct" are reported, coarse to fine.  "graphs":
    distinct as mass-labelled graphs with the legs unlabelled (legs told apart
    only as legs) -- the ambiguity of the graph as an object; "realization_unique"
    is derived from it.  "graphs_leg_labelled": masses on the lines AND leg
    names on the legs -- two covers with the same leg-labelled graph give the
    same leg-order planarity, cut signature and legs per vertex, so this is
    the ambiguity that can change a verdict; the graphs after the first are
    handed on as the alternatives to check a verdict against.
    "graphs_edge_labelled": edge identities (which propagator attaches where)
    -- finest; a symmetric insertion (a kite on a rung) has several of these
    for one labelled graph.

    Returns a dict: "graphs", "graphs_leg_labelled", "graphs_edge_labelled"
    (MultiGraph, in the order found), "covers" (every cover collected), "connected_covers",
    "legs_per_vertex_by_realization", "n_candidate_vertices",
    "n_half_edges", "steps", "budget_exhausted", "enumeration_capped".
    Vertex candidates are enumerated with integer arithmetic (the rational
    coefficients are scaled by their common denominator).
    """
    from fractions import Fraction
    from math import lcm
    n = len(he)
    dim = len(he[0][0]) if he else 0
    fr = [[Fraction(int(sp.Rational(c).p), int(sp.Rational(c).q)) for c in v]
          for (v, _) in he]
    den = 1
    for row in fr:
        for c in row:
            den = lcm(den, c.denominator)
    vec = [tuple(int(c * den) for c in row) for row in fr]
    edge_of = [d[0] for (_, d) in he]
    is_ext = [d[0].startswith("ext_") for (_, d) in he]
    zero = tuple([0] * dim)

    # -- vertex candidates: zero-sum sets of size 3..maxdeg, distinct edges,
    #    not all external; minimal (no zero-sum subset of size >= 3).
    cands = []
    chosen = []
    edges_used = set()

    def dfs(start, acc):
        for k in range(start, n):
            if edge_of[k] in edges_used:
                continue
            acc2 = tuple(a + b for a, b in zip(acc, vec[k]))
            chosen.append(k)
            edges_used.add(edge_of[k])
            if len(chosen) >= 3 and acc2 == zero and not all(is_ext[j] for j in chosen):
                cands.append(tuple(chosen))
            if len(chosen) < maxdeg:
                dfs(k + 1, acc2)
            chosen.pop()
            edges_used.discard(edge_of[k])

    dfs(0, zero)
    cands.sort(key=len)
    cand_set = set()
    minimal = []
    for c in cands:
        ok = True
        for s2 in range(3, len(c)):
            for sub in itertools.combinations(c, s2):
                if sub in cand_set:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            minimal.append(c)
            cand_set.add(c)
    cands = minimal
    masks = [sum(1 << k for k in c) for c in cands]
    by_he = [[ci for ci, c in enumerate(cands) if k in c] for k in range(n)]
    full = (1 << n) - 1

    covers = []
    state = {"steps": 0, "exhausted": False, "capped": False}

    def solve(covered, picked):
        state["steps"] += 1
        if state["steps"] > step_budget:
            state["exhausted"] = True
            return
        if len(covers) >= max_covers:
            state["capped"] = True
            return
        if covered == full:
            if target_V is None or len(picked) == target_V:
                covers.append(list(picked))
            return
        if target_V is not None:
            rem_groups = target_V - len(picked)
            rem_he = n - bin(covered).count("1")
            if rem_groups <= 0 or rem_he < 3 * rem_groups or rem_he > maxdeg * rem_groups:
                return
        best, bestc = None, None
        for k in range(n):
            if (covered >> k) & 1:
                continue
            opts = [ci for ci in by_he[k] if not (masks[ci] & covered)]
            if best is None or len(opts) < len(bestc):
                best, bestc = k, opts
            if not opts:
                return
        for ci in bestc:
            picked.append(cands[ci])
            solve(covered | masks[ci], picked)
            picked.pop()
            if state["exhausted"] or state["capped"]:
                return

    if n:
        solve(0, [])

    by_edge = {}
    for k, (mom, (eid, end)) in enumerate(he):
        by_edge.setdefault(eid, {})[end] = k
    mass_of = {}
    for i, m in enumerate(masses or []):
        mass_of["e%d" % i] = "0" if m == 0 else str(m)

    def build(cover):
        G = nx.MultiGraph()
        hv = {}
        for vi, group in enumerate(cover):
            for k in group:
                hv[k] = "v%d" % vi
        for eid, ends in by_edge.items():
            if eid.startswith("ext_"):
                G.add_node(eid, external=True, leg=eid[4:])
                G.add_edge(eid, hv[ends["in"]], key=eid, id=eid, mass="ext")
            elif "head" in ends and "tail" in ends:
                G.add_edge(hv[ends["head"]], hv[ends["tail"]], key=eid, id=eid,
                           mass=mass_of.get(eid, "0"))
        for v in G.nodes:
            if "external" not in G.nodes[v]:
                G.nodes[v]["external"] = False
        return G

    def legs_as_legs(a, b):
        return a.get("external", False) == b.get("external", False)

    def legs_by_name(a, b):
        return legs_as_legs(a, b) and a.get("leg") == b.get("leg")

    def by_ids(a, b):
        return sorted(d["id"] for d in a.values()) == sorted(d["id"] for d in b.values())

    def by_masses(a, b):
        return sorted(d["mass"] for d in a.values()) == sorted(d["mass"] for d in b.values())

    def same(A, B, nm, em):
        return nx.is_isomorphic(A, B, node_match=nm, edge_match=em)

    graphs, graphs_ll, graphs_el, connected, legs_per_vertex = [], [], [], [], []
    for cover in covers:
        G = build(cover)
        V = sum(1 for v in G.nodes if not G.nodes[v].get("external"))
        E = sum(1 for _, _, k in G.edges(keys=True) if not str(k).startswith("ext_"))
        if not nx.is_connected(G) or (E - V + 1) != n_loops:
            continue
        connected.append(cover)
        legs_per_vertex.append(sorted(
            [sum(1 for k in group if is_ext[k]) for group in cover], reverse=True))
        if not any(same(G, H, legs_by_name, by_ids) for H in graphs_el):
            graphs_el.append(G)
            if not any(same(G, H, legs_by_name, by_masses) for H in graphs_ll):
                graphs_ll.append(G)
                if not any(same(G, H, legs_as_legs, by_masses) for H in graphs):
                    graphs.append(G)
    return {"graphs": graphs, "graphs_leg_labelled": graphs_ll,
            "graphs_edge_labelled": graphs_el, "covers": covers,
            "connected_covers": connected,
            "legs_per_vertex_by_realization": legs_per_vertex,
            "n_candidate_vertices": len(cands), "n_half_edges": n,
            "steps": state["steps"], "step_budget": step_budget,
            "budget_exhausted": state["exhausted"],
            "max_covers": max_covers, "enumeration_capped": state["capped"]}


def _exact_vertex_realizer(vecs, fam, nloop, n_ext, ext_resolved, syms, masses=None,
                           step_budget=REALIZER_STEP_BUDGET,
                           max_covers=REALIZER_MAX_COVERS):
    """
    Exact vertex realization by momentum conservation (see build_graph).

      * half-edges: for edge e, head carries +q_e, tail carries -q_e; for each
        external leg one half-edge carrying +p_j (inward);
      * non_graph_test on the lines (count, loop rank, distinct lines);
      * connected_covers: every exact cover of the half-edges by minimal
        zero-sum vertices with V = E - L + 1 (all covers when the count test
        already failed, to report what exists), kept only when connected;
        distinct realizations counted (ambiguity).

    Returns (G, ok).  G is the first distinct connected realization (or an
    empty MultiGraph when there is none); G.graph["realization"] is the report
    with "verdict" ("GRAPH" / "NON-GRAPH" / "NOT-CHECKABLE" for a leg set that
    cannot close / "UNDECIDED" when the search budget ran out), a named
    "cause" when not ok, V, E, L, N, the count numbers, the cover counts,
    "n_distinct_connected_realizations" (mass-labelled graphs, legs
    unlabelled), "n_distinct_leg_labelled_realizations" (masses + leg names),
    "n_distinct_edge_labelled_realizations", "realization_unique",
    "legs_per_vertex", "vertex_degrees", "conservation_at_every_vertex",
    "connected", and the search-budget fields.  G.graph["alternative_realizations"]
    lists the other distinct leg-labelled realizations, for consumers to check
    a verdict against (cut_signature does).
    """
    he = []   # list of (momentum_vector, descriptor)
    for e, v in enumerate(vecs):
        full = tuple(sp.nsimplify(c) for c in v)
        he.append((full, ("e%d" % e, "head")))
        he.append((tuple(-c for c in full), ("e%d" % e, "tail")))
    # external legs (independent + dependent, all original legs) inward
    for j, name in enumerate(fam.exts):
        ev = sp.expand(ext_resolved[syms[name]])
        vec = [sp.Integer(0)] * nloop + [sp.nsimplify(ev.coeff(syms[e2], 1))
                                         for e2 in fam.exts]
        he.append((tuple(vec), ("ext_%s" % name, "in")))

    ngt = non_graph_test(vecs, fam)
    target_V = ngt["V"] if ngt["count_test_3V_le_2E_plus_N"] else None
    if ngt["legs_sum_to_zero_test"]:
        res = connected_covers(he, nloop, target_V=target_V, step_budget=step_budget,
                               max_covers=max_covers, masses=masses)
    else:   # no cover can close; do not search
        res = {"graphs": [], "graphs_leg_labelled": [], "graphs_edge_labelled": [],
               "covers": [], "connected_covers": [],
               "legs_per_vertex_by_realization": [], "n_candidate_vertices": 0,
               "n_half_edges": len(he), "steps": 0, "step_budget": step_budget,
               "budget_exhausted": False, "max_covers": max_covers,
               "enumeration_capped": False}
    n_covers = len(res["covers"])
    n_conn = len(res["connected_covers"])
    n_distinct = len(res["graphs"])
    ok = ngt["passes"] and n_distinct > 0
    report = dict(ngt)
    report.update({
        "method": ("exact cover of the half-edges by minimal zero-sum vertices "
                   "(size 3..%d, distinct edges, not all external); V = E - L + 1; "
                   "connected; distinct realizations as labelled graphs (masses, "
                   "leg names) and as edge-identity-labelled graphs"
                   % REALIZER_MAX_VERTEX_DEGREE),
        "n_half_edges": res["n_half_edges"],
        "n_candidate_vertices": res["n_candidate_vertices"],
        "n_exact_covers_found": n_covers,
        "covers_restricted_to_V_eq_E_minus_L_plus_1": target_V is not None,
        "n_covers_with_V_eq_E_minus_L_plus_1": sum(
            1 for c in res["covers"] if len(c) == ngt["V"]),
        "n_connected_covers": n_conn,
        "n_distinct_connected_realizations": n_distinct,
        "n_distinct_leg_labelled_realizations": len(res["graphs_leg_labelled"]),
        "n_distinct_edge_labelled_realizations": len(res["graphs_edge_labelled"]),
        "realization_unique": (n_distinct == 1) if n_distinct else None,
        "legs_per_vertex_by_realization": res["legs_per_vertex_by_realization"],
        "steps": res["steps"], "step_budget": res["step_budget"],
        "budget_exhausted": res["budget_exhausted"],
        "max_covers": res["max_covers"], "enumeration_capped": res["enumeration_capped"],
    })
    if ok:
        G = res["graphs_leg_labelled"][0]
        G.graph["alternative_realizations"] = list(res["graphs_leg_labelled"][1:])
        cover = res["connected_covers"][0]
        dim = nloop + n_ext
        cons = []
        for group in cover:
            acc = [sp.Integer(0)] * dim
            for k in group:
                for c in range(dim):
                    acc[c] += he[k][0][c]
            cons.append(all(x == 0 for x in acc))
        report.update({
            "verdict": "GRAPH", "cause": None, "connected": True,
            "vertex_degrees": sorted([len(g) for g in cover], reverse=True),
            "legs_per_vertex": res["legs_per_vertex_by_realization"][0],
            "conservation_at_every_vertex": all(cons),
        })
        n_ll = len(res["graphs_leg_labelled"])
        if n_ll > 1:
            report["ambiguity"] = ("%d distinct connected realizations as leg-labelled "
                                   "graphs (%d with the legs unlabelled) with V = %d; "
                                   "the first is reported, the others are carried as "
                                   "alternatives" % (n_ll, n_distinct, ngt["V"]))
    else:
        G = nx.MultiGraph()
        cause = ngt["cause"]
        if cause is None:
            if res["budget_exhausted"]:
                cause = ("realization not decided: step budget %d exhausted in the "
                         "exact-cover search (%d covers collected)"
                         % (res["step_budget"], n_covers))
            elif n_covers == 0:
                cause = ("NON-GRAPH: no exact cover of the %d half-edges by zero-sum "
                         "vertices with V = E - L + 1 = %d exists (E = %d, L = %d, N = %d)"
                         % (res["n_half_edges"], ngt["V"], ngt["E"], ngt["L"], ngt["N"]))
            elif res["enumeration_capped"]:
                cause = ("realization not decided: cover cap %d reached with no "
                         "connected cover among the first %d" % (res["max_covers"], n_covers))
            else:
                cause = ("NON-GRAPH: %d exact cover(s) with V = E - L + 1 = %d found, "
                         "none connected with E - V + 1 = L" % (n_covers, ngt["V"]))
        elif n_covers:
            cause += ("; %d exact cover(s) of the half-edges found, %d with V = %d"
                      % (n_covers, report["n_covers_with_V_eq_E_minus_L_plus_1"], ngt["V"]))
        report.update({
            "verdict": ("NON-GRAPH" if cause.startswith("NON-GRAPH") else
                        "NOT-CHECKABLE" if cause.startswith("NOT-CHECKABLE") else
                        "UNDECIDED"),
            "cause": cause, "connected": None, "vertex_degrees": None,
            "legs_per_vertex": None, "conservation_at_every_vertex": None,
        })
    G.graph["realization"] = report
    return G, ok


#  Planarity of the realized graph with the external legs on the outer face.
#
#  The closure-cycle test (planarity) only DISCRIMINATES planar from crossed
#  when at least four legs sit on the boundary: a two- or three-leg closure
#  is a label check on the realized graph, never a planar/crossed verdict,
#  and the field says so by name ("not discriminating (N < 4 legs)").  The
#  cyclic order the legs are closed in is a DATUM when the record carries one
#  (a drawn graph's cyclic_leg_order, a caller's leg_order); the canonical
#  index order p1, p2, ..., pN is the fallback convention and is named as
#  such.  Beside the closure test the legs-joined-at-infinity criterion of
#  the identity census (planarity_at_infinity) is always reported: it is the
#  Feynman-graph criterion (an embedding with EVERY leg on the outer face
#  exists) and discriminates for any number of legs.  The record's Mandelstam
#  declaration (s = (p1+p2)^2, t = (p1+p3)^2, ...) is read and REPORTED as a
#  candidate cyclic order with its own closure verdict, never enforced: a
#  record's channel names need not follow the drawn boundary order (the
#  lbl3se record names its crossed channel t while its box is drawn
#  p1,p2,p3,p4; Smirnov's double box names its adjacent channel t and sits
#  p1,p2,p4,p3), so which order decides is the reader's call (the compare
#  mode's, given the drawing).
_PLANARITY_CANONICAL_ORDER_SOURCE = ("canonical leg index order (no cyclic order "
                                     "declared by the record: a convention, not a datum)")
_PLANARITY_NOT_DISCRIMINATING = "not discriminating (N < 4 legs)"
_PLANARITY_ENUMERATION_MAX_LEGS = 7


def _cyclic_leg_orders(legs):
    """Every cyclic order of `legs` up to rotation and reflection: the first
    leg fixed, the rest permuted, one of each mirror pair kept -- (N-1)!/2
    orders for N >= 4, the single order for N <= 3."""
    legs = list(legs)
    if len(legs) <= 3:
        return [legs]
    first, rest = legs[0], legs[1:]
    out = []
    for perm in itertools.permutations(rest):
        if str(perm[0]) < str(perm[-1]):      # drops the reversed (mirror) order
            out.append([first] + list(perm))
    return out


def _same_cyclic_order(a, b):
    """True when the two leg lists are the same cyclic order up to rotation
    and reflection."""
    a, b = [str(x) for x in a], [str(x) for x in b]
    if len(a) != len(b) or sorted(a) != sorted(b):
        return False
    n = len(a)
    if n <= 3:
        return True
    for r in range(n):
        rot = b[r:] + b[:r]
        if rot == a or rot[::-1] == a:
            return True
    return False


def _closure_leg_order(fam, attached, leg_order=None):
    """The cyclic order the closure test closes the legs in, and where it came
    from: the caller's `leg_order` (refused by name unless it is a
    permutation of the attached legs), else the record's Family.cyclic_leg_order
    (a drawn graph's boundary order, a pySecDec graph's leg list) when it
    names the attached legs, else the canonical index order.  Returns
    (order, source, note)."""
    attached = [str(x) for x in attached]
    if leg_order is not None:
        lo = [str(x) for x in leg_order]
        if sorted(lo) != sorted(attached):
            raise ValueError(
                "REFUSED: PLANARITY-LEG-ORDER: the leg order %s is not a permutation "
                "of the attached external legs %s of family %r" % (lo, attached, fam.name))
        return lo, "leg_order argument (the caller's declared cyclic order)", None
    cyc = getattr(fam, "cyclic_leg_order", None)
    if cyc:
        cyc = [str(x) for x in cyc]
        if sorted(cyc) == sorted(attached):
            return cyc, "family.cyclic_leg_order (the record's drawn / declared cyclic order)", None
        return attached, _PLANARITY_CANONICAL_ORDER_SOURCE, (
            "the record's cyclic_leg_order %s does not name the attached legs %s; "
            "the canonical index order is used" % (cyc, attached))
    return attached, _PLANARITY_CANONICAL_ORDER_SOURCE, None


def planarity_at_infinity(fam_or_graph, legs_in_order=None):
    """
    Planarity three ways on the realized graph, in the identity census's
    planarity_closure form (census_lib_rows06_24_27_30.planarity_closure):

      abstract_planar                  the graph alone (legs ignored);
      planar_canonical_leg_order       the external legs closed by a cycle in
                                       the given cyclic order (the tool's
                                       closure method; fewer than three
                                       attached legs: the abstract value);
      planar_legs_joined_at_infinity   every leg joined to ONE vertex at
                                       infinity -- the Feynman-graph
                                       criterion: an embedding with every
                                       external leg on the outer face exists;
                                       it discriminates for any number of legs;
      closure_order                    the order used;
      discriminating                   whether the closure test can separate
                                       planar from crossed here (>= 4 legs).

    fam_or_graph is a Family (realized through build_graph) or a networkx
    graph carrying the 'ext_<leg>' nodes; legs_in_order defaults to the
    family's leg list as the census used it.  Returns None when the family
    has no realization.
    """
    if not _HAVE_NX:
        return None
    if isinstance(fam_or_graph, Family):
        G, ok = build_graph(fam_or_graph)
        if not ok or G is None or G.number_of_edges() == 0:
            return None
        if legs_in_order is None:
            legs_in_order = list(fam_or_graph.exts)
    else:
        G = fam_or_graph
    base = nx.Graph()
    if G.is_multigraph():
        base.add_edges_from((u, v) for u, v, k in G.edges(keys=True))
    else:
        base.add_edges_from(G.edges())
    if legs_in_order is None:
        legs_in_order = [str(n)[4:] for n in base.nodes() if str(n).startswith("ext_")]
    legs_in_order = [str(x) for x in legs_in_order]
    abstract = nx.check_planarity(base)[0]
    ext = ["ext_" + l for l in legs_in_order if ("ext_" + l) in base]
    if len(ext) >= 3:
        H = base.copy()
        for a in range(len(ext)):
            H.add_edge(ext[a], ext[(a + 1) % len(ext)])
        canon = nx.check_planarity(H)[0]
    else:
        canon = abstract
    Hinf = base.copy()
    for e in ext:
        Hinf.add_edge("INFINITY", e)
    inf_planar = nx.check_planarity(Hinf)[0]
    return {"planar_canonical_leg_order": bool(canon), "abstract_planar": bool(abstract),
            "planar_legs_joined_at_infinity": bool(inf_planar),
            "closure_order": list(legs_in_order), "discriminating": len(ext) >= 4}


def _mandelstam_declaration(fam, legs):
    """
    The record's declared two-particle channels and the cyclic leg order they
    imply, REPORTED (never enforced by planarity).  Sources, in order: the
    kinematics.yaml recorded on the family (scalarproduct_rules), a .jl
    AmflowFamily's kinematics Dict (its fourth argument), an AMFlow-port
    JSON's replacement block.  A rule pi.pj = expr declares the channel X for
    the pair {pi, pj} when 2*expr + pi^2 + pj^2 is the single symbol X (so
    p1.p2 = s/2 with massless legs declares s = (p1+p2)^2; p2.p3 =
    (7/25 - s - t)/2 declares nothing).  A cyclic order is implied when
    exactly one cyclic order (up to rotation and reflection) has every
    declared pair adjacent; for four legs a pair and its complement are the
    same channel and are adjacent together, so s = (p1+p2)^2, t = (p1+p3)^2
    implies p1,p2,p4,p3.  Returns None when the record carries no
    declaration; otherwise {"source", "channels": {X: [[a, b], ...]},
    "n_candidate_orders", "cyclic_order" (or None), "note"}.
    """
    legs = [str(x) for x in legs]
    rules, source = None, None
    try:
        kp = getattr(fam, "kinematics_path", None)
        if kp and os.path.exists(str(kp)):
            rules = read_kinematics_yaml(str(kp)).get("rules") or {}
            source = "kinematics.yaml scalarproduct_rules (%s)" % os.path.basename(str(kp))
        elif (getattr(fam, "source_kind", None) == "amflow_jl" or
              str(fam.source).lower().endswith(".jl")) and os.path.exists(str(fam.source)):
            txt = open(fam.source).read()
            i = txt.find("AmflowFamily(")
            if i >= 0:
                j = i + len("AmflowFamily(")
                depth, k = 1, j
                while k < len(txt) and depth:
                    if txt[k] == '(':
                        depth += 1
                    elif txt[k] == ')':
                        depth -= 1
                    k += 1
                args = _split_top_commas(txt[j:k - 1])
                if len(args) > 4:
                    rules = _kin_from_replacement(_jl_dict(args[4]), legs)["rules"]
                    source = ("AmflowFamily kinematics Dict (%s)"
                              % os.path.basename(str(fam.source)))
        elif fam.record_notes.get("replacement"):
            rules = _kin_from_replacement(fam.record_notes["replacement"], legs)["rules"]
            source = "AMFlow-port JSON replacement block"
    except Exception as exc:                       # a declaration that does not parse is named
        return {"source": source, "channels": {}, "n_candidate_orders": None,
                "cyclic_order": None, "note": "declaration not readable: %s: %s"
                % (type(exc).__name__, exc)}
    if not rules:
        return None
    squares = {}
    for (a, b), v in rules.items():
        if a == b:
            try:
                squares[a] = sp.sympify(str(v))
            except Exception:
                continue
    channels = {}
    for (a, b), v in rules.items():
        if a == b or a not in legs or b not in legs or a not in squares or b not in squares:
            continue
        try:
            tot = sp.expand(2 * sp.sympify(str(v)) + squares[a] + squares[b])
        except Exception:
            continue
        if isinstance(tot, sp.Symbol):
            channels.setdefault(str(tot), []).append(sorted([a, b]))
    out = {"source": source, "channels": channels, "n_candidate_orders": None,
           "cyclic_order": None, "note": None}
    if not channels:
        out["note"] = ("no two-particle channel is declared as a single invariant "
                       "(no pi.pj rule of the form (X - pi^2 - pj^2)/2)")
        return out
    if len(legs) < 4:
        out["note"] = "fewer than four legs: one cyclic order only"
        return out
    pairs = {frozenset(p) for ps in channels.values() for p in ps}

    def adjacent(order, pair):
        a, b = tuple(pair)
        i, j = order.index(a), order.index(b)
        return (j - i) % len(order) in (1, len(order) - 1)
    if len(legs) > _PLANARITY_ENUMERATION_MAX_LEGS:
        out["note"] = "more than %d legs: cyclic orders not enumerated" % _PLANARITY_ENUMERATION_MAX_LEGS
        return out
    cands = [o for o in _cyclic_leg_orders(legs) if all(adjacent(o, p) for p in pairs)]
    out["n_candidate_orders"] = len(cands)
    if len(cands) == 1:
        out["cyclic_order"] = list(cands[0])
        out["note"] = ("the declared channels %s are adjacent in exactly one cyclic order"
                       % {k: v for k, v in channels.items()})
    elif not cands:
        out["note"] = ("the declared channels %s are not all adjacent in any cyclic order "
                       "(a crossed channel is declared as a single invariant)" % channels)
    else:
        out["note"] = ("the declared channels %s leave %d cyclic orders: no order implied"
                       % (channels, len(cands)))
    return out


def planarity(fam, leg_order=None):
    """
    PHYSICAL planarity of the Feynman diagram.

    A double box and a *crossed* double box have the SAME abstract graph up to
    the cyclic order in which the external legs sit on the outer boundary, so
    abstract-graph planarity does NOT distinguish them.  The physical question
    is: can the diagram be drawn in the plane with the external legs in a
    given cyclic order?  Equivalently, is the graph planar after we ADD a
    closure cycle linking the external legs in that order?

    planar == True   <=>  drawable with legs in the closure order   (planar box)
    planar == False  <=>  forces a crossing                          (crossed box)

    The closure order is `leg_order` (the caller's; refused by name unless it
    is a permutation of the attached legs), else the record's
    Family.cyclic_leg_order, else the canonical index order p1..pN -- named
    in "closure_order" / "closure_order_source" and in the method string.
    The closure test DISCRIMINATES only with >= 4 attached legs
    ("discriminating"); with two or three the method string reads
    "not discriminating (N < 4 legs): ..." and `planar` is the closure /
    abstract value as a label check.  Reported beside it:
    `abstract_planar`; `planar_legs_joined_at_infinity` (every leg joined to
    one vertex at infinity -- the census's Feynman-graph criterion, valid
    for any N); `planar_in_some_leg_order` with the cyclic orders ENUMERATED
    (up to 7 legs; `n_leg_orders`, `n_planar_leg_orders`,
    `planar_leg_orders`); and `mandelstam_declaration`, the record's channel
    naming as a candidate order with its own closure verdict
    ("planar_in_declared_order", "agrees_with_closure_order"), never
    enforced.  The report is cached on the family for the same data.
    """
    import copy
    key = (tuple(fam.loops), tuple(fam.exts), tuple(sorted(fam.ext_subs.items())),
           tuple((str(e), str(m)) for e, m in fam.propagators), tuple(fam.physical),
           tuple(fam.cyclic_leg_order or ()), tuple(leg_order or ()) if leg_order is not None else None,
           str(getattr(fam, "kinematics_path", None)), str(fam.source),
           tuple(sorted(fam.mass_values.items())))
    cache = getattr(fam, "_planarity_cache", None)
    if cache is not None and cache[0] == key:
        return copy.deepcopy(cache[1])
    G, ok = build_graph(fam)
    if not _HAVE_NX or not ok or G is None or G.number_of_edges() == 0:
        return {"planar": None, "method": "realization-failed-or-unavailable",
                "n_external_legs": len(fam.exts), "discriminating": None,
                "closure_order": None, "closure_order_source": None,
                "planar_legs_joined_at_infinity": None}
    base = nx.Graph()
    base.add_edges_from((u, v) for u, v, k in G.edges(keys=True))
    abstract_planar = nx.check_planarity(base)[0]

    extnodes = [n for n in base.nodes() if str(n).startswith("ext_")]
    # canonical cyclic order by leg index (a leg without an index keeps the
    # record's declared position: P56, P34, P12 as the family lists them)
    def legidx(n):
        m = re.search(r'p(\d+)', str(n))
        return int(m.group(1)) if m else 0

    def canon_key(n):
        name = str(n)[4:]
        if re.fullmatch(r'p\d+', name):
            return (0, int(name[1:]), "")
        return (1, fam.exts.index(name) if name in fam.exts else len(fam.exts), name)
    canonical = [str(n)[4:] for n in sorted(extnodes, key=canon_key)]
    order_names, source, order_note = _closure_leg_order(fam, canonical, leg_order)
    order = ["ext_" + l for l in order_names]
    all_p = all(re.fullmatch(r'p\d+', l) for l in order_names)
    order_txt = (",".join(str(legidx(n)) for n in order) if all_p
                 else ",".join(order_names))
    declared = source != _PLANARITY_CANONICAL_ORDER_SOURCE
    where = (("canonical order " + order_txt) if not declared
             else ("the declared cyclic order " + ",".join(order_names) + " (" + source + ")"))

    def closure_planar(o):
        H = base.copy()
        m = len(o)
        for a in range(m):
            H.add_edge(o[a], o[(a + 1) % m])
        return nx.check_planarity(H)[0]

    n = len(order)
    planar_orders, n_orders, enumerated = [], None, False
    if n >= 4:
        canon_planar = closure_planar(order)
        discriminating = True
        # planar for SOME cyclic order? (a crossed box IS planar in a crossed
        # order; a truly non-planar graph is False for every order) -- the
        # cyclic orders are enumerated, the abstract shortcut is wrong on the
        # (2,1,1) theta graph (abstract planar, no order planar)
        if n <= _PLANARITY_ENUMERATION_MAX_LEGS:
            orders = _cyclic_leg_orders(order)
            planar_orders = [[str(x)[4:] for x in o] for o in orders if closure_planar(o)]
            n_orders = len(orders)
            some_order_planar = bool(planar_orders)
            enumerated = True
        else:
            some_order_planar = abstract_planar   # abstract planar => some order works
        method = ("networkx.check_planarity with external-leg closure cycle in " + where)
    elif n == 3:
        canon_planar = closure_planar(order)
        discriminating = False
        some_order_planar = canon_planar       # one cyclic order of three legs
        planar_orders, n_orders, enumerated = ([list(order_names)] if canon_planar else []), 1, True
        method = (_PLANARITY_NOT_DISCRIMINATING + ": networkx.check_planarity with "
                  "external-leg closure cycle in " + where)
    else:
        canon_planar = abstract_planar
        discriminating = False
        some_order_planar = abstract_planar
        planar_orders, n_orders, enumerated = ([list(order_names)] if abstract_planar else []), (1 if n else 0), True
        method = (_PLANARITY_NOT_DISCRIMINATING + ": networkx.check_planarity "
                  "(abstract; <3 external legs)")
    at_inf = planarity_at_infinity(G, order_names)
    md = _mandelstam_declaration(fam, canonical)
    if md and md.get("cyclic_order"):
        mo = ["ext_" + l for l in md["cyclic_order"]]
        md["planar_in_declared_order"] = bool(closure_planar(mo) if n >= 3 else abstract_planar)
        md["agrees_with_closure_order"] = _same_cyclic_order(md["cyclic_order"], order_names)
    out = {"planar": bool(canon_planar),
           "abstract_planar": bool(abstract_planar),
           "planar_in_some_leg_order": bool(some_order_planar),
           "method": method,
           "n_vertices": base.number_of_nodes(),
           "n_edges": G.number_of_edges(),
           "n_external_legs": n,
           "discriminating": discriminating,
           "not_discriminating_note": (None if discriminating else
                                       _PLANARITY_NOT_DISCRIMINATING + ": the leg-order closure "
                                       "test cannot separate planar from crossed with fewer than "
                                       "four legs on the boundary; planar is a label check on the "
                                       "realized graph"),
           "closure_order": list(order_names),
           "closure_order_source": source,
           "closure_order_note": order_note,
           "planar_legs_joined_at_infinity": at_inf["planar_legs_joined_at_infinity"],
           "leg_orders_enumerated": enumerated,
           "n_leg_orders": n_orders,
           "n_planar_leg_orders": len(planar_orders) if enumerated else None,
           "planar_leg_orders": planar_orders if enumerated else None,
           "mandelstam_declaration": md}
    fam._planarity_cache = (key, copy.deepcopy(out))
    return out


class _CutSignature(dict):
    """
    Channel -> minimal internal cut (int), or None.  A plain mapping for JSON
    and for every consumer that reads the channel values; two attributes name
    what a None entry means:

      cause        None when every channel has a value; otherwise the named
                   reason (the realization's NON-GRAPH cause, or, per channel,
                   the two legs that attach to one vertex so no internal cut
                   separates them);
      realization  the realization report of build_graph (G.graph["realization"]);
      alternatives the cut signatures of the other distinct connected
                   realizations (labelled-graph notion), when there are any;
      ambiguity    None, or a sentence naming the alternatives whose cut
                   signature differs from the reported one.
    """
    cause = None
    realization = None
    alternatives = None
    ambiguity = None


def cut_signature(fam):
    """
    For each Mandelstam channel, the minimal number of INTERNAL propagators a
    cut must sever to separate that 2-particle external cluster from the rest.

    Computed as a graph min-cut (networkx) between the two external clusters on
    the realized graph (build_graph: connected, V = E - L + 1), counting only
    internal edges.  The (s,t,u) multiset is the label-independent fingerprint:
    planar dbox -> t != u (asymmetric); genuine crossed box -> t == u (t<->u
    symmetric).  Returns a _CutSignature: a dict whose None entries carry a
    named cause (.cause) -- no realization (the NON-GRAPH cause by name), or
    two legs of opposite clusters attached to the same vertex.
    """
    G, ok = build_graph(fam)
    nums = sorted(int(e[1:]) for e in fam.exts if re.fullmatch(r'p\d+', e))
    if nums == [1, 2, 3, 4]:
        channels = {('p1', 'p2'): 's', ('p2', 'p3'): 't', ('p1', 'p3'): 'u'}
    else:
        channels = {}
        for a, b in itertools.combinations(fam.exts, 2):
            channels[(a, b)] = f"{a}{b}"
    out = _CutSignature()
    report = G.graph.get("realization") if G is not None else None
    out.realization = report
    if not _HAVE_NX or not ok or G is None:
        for v in channels.values():
            out[v] = None
        out.cause = ((report or {}).get("cause")
                     or "realization unavailable (networkx not importable)")
        return out
    attach = {}
    for e in fam.exts:
        node = "ext_%s" % e
        if node in G:
            attach[e] = sorted(str(v) for v in G.neighbors(node))
    causes = []
    for (a, b), cname in channels.items():
        out[cname] = _min_internal_cut(G, fam, a, b)
        if out[cname] is None:
            shared = [(x, y) for x in (a, b) for y in fam.exts if y not in (a, b)
                      and attach.get(x) and attach.get(x) == attach.get(y)]
            if shared:
                causes.append("%s: legs %s attach to one vertex (legs per vertex %s), "
                              "no internal cut separates {%s,%s} from the rest"
                              % (cname, ", ".join("%s/%s" % s for s in shared),
                                 (report or {}).get("legs_per_vertex"), a, b))
            else:
                causes.append("%s: min-cut undefined on the realized graph" % cname)
    out.cause = "; ".join(causes) if causes else None
    alts = G.graph.get("alternative_realizations") or []
    if alts:
        out.alternatives = [{cname: _min_internal_cut(H, fam, a, b)
                             for (a, b), cname in channels.items()} for H in alts]
        differing = [i + 2 for i, alt in enumerate(out.alternatives) if alt != dict(out)]
        if differing:
            out.ambiguity = ("cut signature differs between the %d distinct connected "
                             "realizations: realization(s) %s give %s"
                             % (len(alts) + 1, differing,
                                [out.alternatives[i - 2] for i in differing]))
    return out


def _min_internal_cut(G, fam, a, b):
    """Min number of internal edges separating ext cluster {a,b} from the rest."""
    na, nb = "ext_%s" % a, "ext_%s" % b
    others = ["ext_%s" % e for e in fam.exts if e not in (a, b)]
    if na not in G or nb not in G or any(o not in G for o in others):
        return None
    # contract {a,b} into a super-source, others into a super-sink, then
    # min edge cut counting only internal edges (give ext legs infinite weight).
    H = nx.Graph()
    INF = 10 ** 6
    for u, v, k in G.edges(keys=True):
        w = INF if str(k).startswith("ext_") else 1
        if H.has_edge(u, v):
            H[u][v]["capacity"] += w
        else:
            H.add_edge(u, v, capacity=w)
    # merge sources / sinks
    src, snk = "SRC", "SNK"
    H.add_node(src); H.add_node(snk)
    for node, supr in [(na, src), (nb, src)] + [(o, snk) for o in others]:
        if node in H:
            for nbr in list(H.neighbors(node)):
                cap = H[node][nbr]["capacity"]
                tgt = supr
                if H.has_edge(tgt, nbr):
                    H[tgt][nbr]["capacity"] += cap
                else:
                    H.add_edge(tgt, nbr, capacity=cap)
            H.remove_node(node)
    try:
        val, _ = nx.minimum_cut(H, src, snk)
    except Exception:
        return None
    if val >= INF:
        return None
    return int(val)


# ---------------------------------------------------------------------------
#  4.  CANONICAL CATALOG OF KNOWN FAMILIES
# ---------------------------------------------------------------------------
#  A small set of textbook representatives.  Each entry stores the canonical
#  propagator-vector set so a new family can be matched against it under signed
#  loop relabeling.  Loaded lazily; extend by editing CATALOG_SOURCES below
#  (loops, exts, ext_subs, propagators, note, legs).  The "planar" key
#  declares the expected leg-order planarity verdict; the self-test checks
#  every entry against it (leg [6]).  "planar_in_some_leg_order" declares
#  whether ANY cyclic order of the legs on the outer face is planar (a crossed
#  box is False for every order; a planar box written in a crossed leg
#  convention is True).  "legs" declares the virtuality class of every
#  external leg, the eliminated one included: "0" = massless on shell,
#  "offshell" = generic virtuality; a mass symbol would name a leg at that
#  mass (none in this catalog).  Off-shell legs are part of a family's
#  identity (a three-point off-shell ladder is not the on-shell four-point
#  ladder with a leg pair merged), so every entry carries the classes.
#
#  Retired entries (CATALOG_RETIRED, below) are propagator lists that are easy
#  to mistake for a catalog family.  They are kept by name as negative
#  controls, NOT searched, and the catalog self-test leg proves what each
#  actually is.

CATALOG_SOURCES = {
    "planar_smirnov_dbox": {
        "loops": ["k1", "k2"], "exts": ["p1", "p2", "p3", "p4"],
        "ext_subs": {"p4": "-p1 - p2 - p3"},
        # massless planar double box, Smirnov hep-ph/9905323 routing
        "propagators": [
            ("k1^2", 0), ("(k1 + p1)^2", 0), ("(k1 + p1 + p2)^2", 0),
            ("k2^2", 0), ("(k2 + p1 + p2)^2", 0), ("(k2 + p1 + p2 + p3)^2", 0),
            ("(k1 - k2)^2", 0),
        ],
        "note": "massless planar double box (Smirnov hep-ph/9905323)",
        "planar": True,
        "planar_in_some_leg_order": True,
        "legs": {"p1": "0", "p2": "0", "p3": "0", "p4": "0"},
    },
    "crossed_tausk_dbox": {
        "loops": ["k1", "k2"], "exts": ["p1", "p2", "p3", "p4"],
        "ext_subs": {"p4": "-p1 - p2 - p3"},
        # The genuine nonplanar crossed double box (Tausk hep-ph/9909506) is
        # the (2,1,1) theta graph: two 3-valent hub vertices mt, mb joined by
        # three paths, mt-a1-b1-mb (legs p2 at a1, p1 at b1), mt-ct-mb (leg p3
        # at ct) and mt-cb-mb (leg p4 at cb); six vertices, seven massless
        # lines, four on-shell legs.  Non-planar with the legs on the outer
        # face in EVERY cyclic leg order (not only the canonical one), cut
        # signature {s: 2, t: 3, u: 3}, and not isomorphic as a bare graph to
        # the planar double box.  Three lines carry both loop momenta in this
        # routing (routing-dependent, informational).  Its easy impostor, the
        # planar box with legs 3, 4 interchanged, is the negative control
        # CATALOG_RETIRED["planar_smirnov_dbox_legs34_interchanged"].
        "propagators": [
            ("(-k1 - k2 + p3 + p4)^2", 0), ("(-k1 - k2 - p1)^2", 0),
            ("(-k1 - k2)^2", 0), ("(k1 - p3)^2", 0), ("k1^2", 0),
            ("(k2 + p1 + p2 + p3)^2", 0), ("k2^2", 0),
        ],
        "note": "genuine nonplanar crossed double box, the (2,1,1) theta graph "
                "(Tausk hep-ph/9909506)",
        "planar": False,
        "planar_in_some_leg_order": False,
        "legs": {"p1": "0", "p2": "0", "p3": "0", "p4": "0"},
    },
    "massless_box_1l": {
        "loops": ["k1"], "exts": ["p1", "p2", "p3", "p4"],
        "ext_subs": {"p4": "-p1 - p2 - p3"},
        "propagators": [
            ("k1^2", 0), ("(k1 + p1)^2", 0), ("(k1 + p1 + p2)^2", 0),
            ("(k1 + p1 + p2 + p3)^2", 0),
        ],
        "note": "massless on-shell one-loop box",
        "planar": True,
        "planar_in_some_leg_order": True,
        "legs": {"p1": "0", "p2": "0", "p3": "0", "p4": "0"},
    },
    "massless_triangle_1l": {
        "loops": ["k1"], "exts": ["p1", "p2", "p3"],
        "ext_subs": {"p3": "-p1 - p2"},
        "propagators": [
            ("k1^2", 0), ("(k1 + p1)^2", 0), ("(k1 + p1 + p2)^2", 0),
        ],
        # 3 legs: the closure cycle is trivially planar, matching is the value.
        # The massless triangle with every leg on shell is scaleless; the
        # textbook object has the merged leg p3 = -p1-p2 off shell.
        "note": "massless one-loop triangle (three-point)",
        "planar": True,
        "planar_in_some_leg_order": True,
        "legs": {"p1": "0", "p2": "0", "p3": "offshell"},
    },
    "massless_bubble_1l": {
        "loops": ["k1"], "exts": ["p1", "p2"],
        "ext_subs": {"p2": "-p1"},
        "propagators": [("k1^2", 0), ("(k1 + p1)^2", 0)],
        # two-point massless objects are scaleless at p1^2 = 0: the leg is off shell
        "note": "massless one-loop bubble (two-point)",
        "planar": True,
        "planar_in_some_leg_order": True,
        "legs": {"p1": "offshell", "p2": "offshell"},
    },
    "sunrise_2l": {
        "loops": ["k1", "k2"], "exts": ["p1", "p2"],
        "ext_subs": {"p2": "-p1"},
        "propagators": [("k1^2", 0), ("k2^2", 0), ("(k1 + k2 - p1)^2", 0)],
        "note": "massless two-loop sunrise (two-point)",
        "planar": True,
        "planar_in_some_leg_order": True,
        "legs": {"p1": "offshell", "p2": "offshell"},
    },
    "ladder_vertex_2l": {
        "loops": ["k1", "k2"], "exts": ["p1", "p2", "p3"],
        "ext_subs": {"p3": "-p1 - p2"},
        # the planar double box with the p3/p4 corner legs merged into the
        # single off-shell leg p3 = -p1-p2: the two-loop form-factor ladder.
        "propagators": [
            ("k1^2", 0), ("(k1 + p1)^2", 0), ("(k1 + p1 + p2)^2", 0),
            ("k2^2", 0), ("(k2 + p1 + p2)^2", 0), ("(k1 - k2)^2", 0),
        ],
        "note": "massless planar two-loop ladder vertex (three-point)",
        "planar": True,
        "planar_in_some_leg_order": True,
        "legs": {"p1": "0", "p2": "0", "p3": "offshell"},
    },
    "planar_triple_box_3l": {
        "loops": ["k1", "k2", "k3"], "exts": ["p1", "p2", "p3", "p4"],
        "ext_subs": {"p4": "-p1 - p2 - p3"},
        # three-rung extension of the planar Smirnov double box routing
        "propagators": [
            ("k1^2", 0), ("(k1 + p1)^2", 0), ("(k1 + p1 + p2)^2", 0),
            ("k2^2", 0), ("(k2 + p1 + p2)^2", 0),
            ("k3^2", 0), ("(k3 + p1 + p2)^2", 0), ("(k3 + p1 + p2 + p3)^2", 0),
            ("(k1 - k2)^2", 0), ("(k2 - k3)^2", 0),
        ],
        "note": "massless planar triple box (three-rung ladder)",
        "planar": True,
        "planar_in_some_leg_order": True,
        "legs": {"p1": "0", "p2": "0", "p3": "0", "p4": "0"},
    },
    "ud_ladder_3pt_offshell_3l": {
        "loops": ["k1", "k2", "k3"], "exts": ["p1", "p2", "p3"],
        "ext_subs": {"p3": "-p1 - p2"},
        # the L=3 Usyukina-Davydychev three-point ladder C^(3)(p1^2, p2^2, p3^2)
        # (Broadhurst–Davydychev arXiv:1007.0237 Fig. 1a): nine massless lines, three
        # loops, three legs; apex carries p3, the far rung carries p1 and p2.
        # Legs p1 and p3 off shell, p2 on shell (the slice y = p2^2/p3^2 -> 0
        # whose holomorphic block is the ladder function Phi^(3)).  It is NOT
        # the on-shell four-point triple box (V,E 10,12 here vs 12,14 there).
        "propagators": [
            ("(k1 + k2 + k3 + p2 + p3)^2", 0), ("(-k1 - k2 - k3 - p2)^2", 0),
            ("(k2 + k3 - p1)^2", 0), ("(k3 - p1)^2", 0),
            ("(-k2 - k3 - p2)^2", 0), ("(-k3 - p2)^2", 0),
            ("k1^2", 0), ("k2^2", 0), ("k3^2", 0),
        ],
        "note": "massless three-loop three-point off-shell ladder, Usyukina-"
                "Davydychev Phi^(3) (legs p1 and p3 off shell, p2 on shell)",
        "planar": True,
        "planar_in_some_leg_order": True,
        "legs": {"p1": "offshell", "p2": "0", "p3": "offshell"},
    },
}

#  Retired entries: negative controls kept by name, never searched.  Each
#  records what the object actually is; the catalog self-test leg asserts it.
CATALOG_RETIRED = {
    "planar_smirnov_dbox_legs34_interchanged": {
        "former_name": "crossed_tausk_dbox",
        "loops": ["k1", "k2"], "exts": ["p1", "p2", "p3", "p4"],
        "ext_subs": {"p4": "-p1 - p2 - p3"},
        # An easy impostor for the crossed Tausk box.  It is the planar
        # Smirnov double box with legs p3 and p4 interchanged: substituting
        # p3 <-> p4 in these strings (with p4 = -p1-p2-p3 re-imposed) gives
        # planar_smirnov_dbox line for line, and as a bare graph it is that
        # box (leg map p4->p1, p3->p2, p1->p3, p2->p4 onto the drawn planar
        # box).  Non-planar only in the canonical leg order (planar in 8 of
        # the 24 orders), cut signature {s: 2, t: 4, u: 3} -- the planar box's
        # multiset -- which is why a "planar is False" test alone accepts it,
        # and why any family in the t = (p1+p3)^2 leg convention would match it
        # and be reported crossed if it were searched.
        "propagators": [
            ("k1^2", 0), ("(k1 + p1)^2", 0), ("(k1 + p1 + p2)^2", 0),
            ("k2^2", 0), ("(k2 + p1 + p2)^2", 0), ("(k2 - p3)^2", 0),
            ("(k1 - k2)^2", 0),
        ],
        "note": "NEGATIVE CONTROL: the planar double box with legs 3 and 4 interchanged "
                "(an easy impostor for crossed_tausk_dbox)",
        "is": "planar_smirnov_dbox",
        "is_under_leg_map": {"p3": "p4", "p4": "p3"},
        "planar": False,
        "planar_in_some_leg_order": True,
        "legs": {"p1": "0", "p2": "0", "p3": "0", "p4": "0"},
    },
}


def _catalog_families():
    out = []
    for name, spec in CATALOG_SOURCES.items():
        fam = Family(name, spec["loops"], spec["exts"], spec["ext_subs"],
                     spec["propagators"], "catalog:" + name)
        out.append((fam, spec.get("note", "")))
    return out


# ---------------------------------------------------------------------------
#  5.  TOP-LEVEL AUDIT
# ---------------------------------------------------------------------------
def fingerprint(fam):
    vecs, masses, isps = momentum_vectors(fam)
    cset = canonical_set(vecs)
    # the LABELS of the lines and legs: the mass multiset, the multiplicity
    # (propagator power) multiset and the leg virtuality classes -- the
    # census carried these beside the momentum set; n_massive_lines is the
    # count derived from the mass multiset
    lab = labeled_set(fam)
    nmass = lab["n_massive_lines"]
    csig = cut_signature(fam)
    plan = planarity(fam)
    # a label-free hash of the canonical set (routing-dependent), and beside
    # it the hash of the labeled set (momentum, mass, multiplicity)
    flat = sorted(tuple(str(c) for c in v) for v in cset)
    import hashlib
    h = hashlib.sha1(json.dumps(flat).encode()).hexdigest()[:12]
    flat_l = sorted((tuple(str(c) for c in v), str(m), int(n)) for v, m, n in lab["lines"])
    h_l = hashlib.sha1(json.dumps(flat_l).encode()).hexdigest()[:12]
    # the ROUTING-INDEPENDENT hash beside the routed ones: the canonical form
    # of the realized labeled graph (masses and multiplicities on the lines,
    # the legs as legs) and, when every leg declares a virtuality class, the
    # same with the classes on the external nodes; None with the reason when
    # the family does not realize (the realizer's verdict) or a class is undeclared
    rh = realized_canonical_hash(fam)
    rhl = realized_canonical_hash(fam, with_leg_classes=True)
    return {
        "name": fam.name,
        "source": fam.source,
        "n_loops": len(fam.loops),
        "n_genuine_propagators": len(vecs),
        "n_isp_or_dotproduct": isps,
        "n_massive_lines": nmass,
        "mass_multiset": lab["mass_multiset"],
        "nu_multiset": lab["nu_multiset"],
        "nu_source": lab["nu_source"],
        "n_dotted_lines": lab["n_dotted_lines"],
        "leg_labels": lab["leg_labels"],
        "cut_signature": csig,
        "cut_signature_multiset": sorted(v for v in csig.values()
                                         if v is not None),
        "planar": plan["planar"],
        "planarity_method": plan["method"],
        "canonical_hash": h,
        "canonical_hash_labeled": h_l,
        "canonical_set_size": len(cset),
        "canonical_hash_realized": rh["hash"],
        "canonical_hash_realized_legs": rhl["hash"],
        "canonical_hash_realized_method": rh["method"],
        "canonical_hash_realized_reason": (rh["reason"] if rh["hash"] is None else
                                          (rhl["reason"] if rhl["hash"] is None else None)),
        "realized_hash_alternatives": rh["alternatives"],
        "realized_V": rh["V"],
        "realized_E": rh["E"],
    }


def audit(fam, try_leg_perms=False):
    fp = fingerprint(fam)
    warnings = []
    vecsA, _, _ = momentum_vectors(fam)
    nloop = len(fam.loops)
    n_ext = len(fam.exts)

    # (a) isomorphism vs catalog.  A catalog entry with a different number of
    #     external legs than the family is SKIPPED by name (no leg relabeling
    #     between different leg counts exists; the former code indexed past
    #     the vector under --leg-perms -- rc 1 on the row-20 class).  Every
    #     searched entry has the target entry's momentum conservation
    #     re-imposed after each leg permutation (Family.ext_subs read as data);
    #     leg permutations are restricted to the label-preserving ones only
    #     when BOTH sides carry a virtuality class for every leg.  A match on
    #     the momentum set alone is a MASS-BLIND match: every match also
    #     reports whether the family's line masses, multiplicities and leg
    #     classes are carried by some relabeling onto the entry's (the
    #     catalog entries are massless, undotted, with declared leg classes),
    #     and the mass pattern is stated with the match, never a bare identity.
    #     The relabeling searched is the signed loop permutation first and,
    #     when none maps the sets, the census's AFFINE class k -> U k + c.p
    #     (isomorphism_to(..., affine=True): a second routing of the same
    #     graph, the row-4 drawn triple box vs planar_triple_box_3l); the
    #     match names its class ("relabeling_class") and, for an affine
    #     match, U, c and the leg map ("affine").  Beside the set search the
    #     ROUTING-INDEPENDENT identity test is run on every entry of the
    #     family's loop and leg count: realized_graph_iso (VF2 on the
    #     realized labeled graphs) -- an entry isomorphic there but reached by
    #     no map in the searched class is still reported (relabeling None,
    #     class "realized-graph isomorphism only"), the census's rule that
    #     identity is VF2's verdict and the map string its describer.
    matches = []
    skipped_leg_count = []
    lab = labeled_set(fam)
    fam_labels = lab["leg_classes"]
    fam_line_labels = list(zip(lab["masses"], lab["nu"]))
    for cfam, note in _catalog_families():
        if len(cfam.loops) != nloop:
            continue
        vecsB, _, _ = momentum_vectors(cfam)
        if len(cfam.exts) != n_ext:
            skipped_leg_count.append(
                {"family": cfam.name, "legs": len(cfam.exts),
                 "n_genuine_propagators": len(vecsB)})
            continue
        spec = CATALOG_SOURCES.get(cfam.name, {})
        entry_classes = None
        if isinstance(spec.get("legs"), dict) and all(e in spec["legs"] for e in cfam.exts):
            entry_classes = [_leg_class(spec["legs"][e]) for e in cfam.exts]
            cfam.leg_virt = dict(spec["legs"])
        labels = None
        if fam_labels is not None and entry_classes is not None:
            labels = (fam_labels, entry_classes)
        names = (list(fam.loops), list(cfam.loops), list(cfam.exts), list(fam.exts))
        det = {}
        rel = isomorphism_to(vecsA, vecsB, nloop, n_ext,
                             try_leg_perms=try_leg_perms,
                             dep_rows=_dependent_leg_rows(cfam),
                             leg_labels=labels, details=det, affine=True, names=names)
        rgi = None
        vf2_only = False
        if len(vecsB) == len(vecsA):
            rgi = realized_graph_iso(fam, cfam)
            # the realized-graph verdict under the flag's own leg semantics:
            # legs fixed -> matched by NAME; --leg-perms -> by virtuality
            # class when the family declares one for every leg, else as legs
            if not try_leg_perms:
                vf2_only = rgi.get("with_masses_and_leg_names") is True
            elif rgi.get("with_masses_and_legs") is not None:
                vf2_only = rgi.get("with_masses_and_legs") is True
            else:
                vf2_only = rgi.get("with_masses") is True
        if rel is None and not vf2_only:
            continue
        # the labeled search: masses and multiplicities carried, leg classes preserved
        clab = labeled_set(cfam)
        det_l = {}
        rel_l = isomorphism_to(vecsA, vecsB, nloop, n_ext,
                               try_leg_perms=try_leg_perms,
                               dep_rows=_dependent_leg_rows(cfam),
                               leg_labels=labels, details=det_l,
                               labels=(fam_line_labels, list(zip(clab["masses"], clab["nu"]))),
                               affine=True, names=names)
        why = []
        if lab["mass_multiset"] != clab["mass_multiset"]:
            if clab["n_massive_lines"] == 0:
                why.append("massless catalog entry; family carries %d massive lines "
                           "(mass multiset %s)" % (lab["n_massive_lines"], lab["mass_multiset"]))
            else:
                why.append("mass multiset %s vs the entry's %s"
                           % (lab["mass_multiset"], clab["mass_multiset"]))
        if lab["nu_multiset"] != clab["nu_multiset"]:
            why.append("family carries %d dotted lines (multiplicities %s); the entry has none"
                       % (lab["n_dotted_lines"], lab["nu_multiset"]))
        if fam_labels is None:
            why.append("leg classes not declared by the record (the entry declares %s)"
                       % (dict(zip(cfam.exts, entry_classes)) if entry_classes else "none"))
        elif rel_l is None and not why:
            why.append("no relabeling carries the masses onto the entry's lines with the "
                       "leg classes %s preserved" % dict(zip(fam.exts, fam_labels)))
        if rel_l is not None and not why:
            pattern = ("labels carried: masses %s, multiplicities %s and leg classes %s "
                       "match the entry under [%s]"
                       % (lab["mass_multiset"], lab["nu_multiset"],
                          dict(zip(fam.exts, fam_labels)), rel_l))
        elif rel_l is not None:
            pattern = ("masses and multiplicities carried under [%s]; %s"
                       % (rel_l, "; ".join(why)))
        else:
            pattern = "MASS-BLIND match (momentum set only): " + "; ".join(why)
        rel_class = det.get("relabeling_class") if rel is not None else None
        if rel is None:
            rel_class = "realized-graph isomorphism only (no map in the searched class)"
            pattern = ("realized-graph isomorphism (VF2) only: no signed permutation and no "
                       "map in the class %s carries the set; %s"
                       % ((det.get("affine") or {}).get("loop_map_class"), pattern))
        aff = None
        if rel is not None and rel_class == "affine unimodular loop redefinition":
            a0 = (det.get("affine") or {})
            aff = {"loop_relabeling": a0.get("loop_relabeling"),
                   "loop_map_matrix": a0.get("loop_map_matrix"),
                   "shift_matrix": a0.get("shift_matrix"), "leg_map": a0.get("leg_map"),
                   "loop_map_class": a0.get("loop_map_class"), "bound": a0.get("bound")}
        matches.append({"family": cfam.name, "relabeling": rel, "note": note,
                        "labels_match": rel_l is not None and not why,
                        "labeled_relabeling": rel_l,
                        "mass_pattern": pattern,
                        "relabeling_class": rel_class,
                        "affine": aff,
                        "realized_graph_iso": ({k: rgi.get(k) for k in
                                                ("bare", "with_masses", "with_masses_and_legs",
                                                 "with_masses_and_legs_reason",
                                                 "with_masses_leg_map", "V_A", "E_A", "V_B", "E_B",
                                                 "error")} if rgi is not None else None)})
    fp["catalog_matches"] = matches
    fp["catalog_skipped_leg_count"] = skipped_leg_count
    for m in matches:
        if m["relabeling"] is None:
            warnings.append(
                f"graph-isomorphic to catalog family '{m['family']}' as REALIZED GRAPHS (VF2 with "
                f"masses and multiplicities; no relabeling in the searched class describes it)  "
                f"({m['note']}); {m['mass_pattern']}")
        elif m["relabeling_class"] == "affine unimodular loop redefinition":
            warnings.append(
                f"graph-isomorphic to catalog family '{m['family']}' under the affine loop "
                f"redefinition [{m['relabeling']}] (a second routing of the same graph: "
                f"{m['affine']['loop_map_class']})  ({m['note']}); {m['mass_pattern']}")
        else:
            warnings.append(
                f"graph-isomorphic to catalog family '{m['family']}' under [{m['relabeling']}]  "
                f"({m['note']}); {m['mass_pattern']}")
    if try_leg_perms:
        same_size = [s for s in skipped_leg_count
                     if s["n_genuine_propagators"] == len(vecsA)]
        if same_size:
            warnings.append(
                "LEG-COUNT SKIP: catalog entries with the same loop and line count "
                f"but a different number of external legs were not searched (this "
                f"family reads {n_ext} legs): "
                + ", ".join(f"{s['family']} ({s['legs']} legs)" for s in same_size)
                + " -- no leg relabeling exists between different leg counts; "
                "if the family's leg set was inferred from its propagator strings, "
                "supply the record's legs")

    # (b) label vs planarity cross-check
    label = fam.name.lower() + " " + os.path.basename(fam.source).lower()
    src_txt = ""
    try:
        src_txt = open(fam.source).read().lower() if os.path.exists(fam.source) else ""
    except Exception:
        pass
    looks_nonplanar = any(w in (label + " " + src_txt[:2000])
                          for w in ("nonplanar", "non-planar", "crossed", "np2l",
                                    "np_", "tausk"))
    # planarity-with-leg-order only DISCRIMINATES when there are >= 4 external
    # legs on the boundary; a 3-cycle closure is trivially planar for any graph,
    # so we do not raise a mismatch in that degenerate case.  Count ALL external
    # legs (including the momentum-conservation-eliminated one, p4), since the
    # closure cycle is built over every attached external node.
    n_boundary_legs = len(fam.exts)
    if fp["planar"] is True and looks_nonplanar and n_boundary_legs >= 4:
        warnings.append(
            "LABEL/PROVENANCE MISMATCH: file is labeled nonplanar/crossed but the "
            "reconstructed graph is PLANAR in the canonical external-leg order "
            "(networkx). This is the NP-dbox-style mislabel. Cross-check the cut "
            f"signature {fp['cut_signature']}.")
    elif fp["planar"] is True and looks_nonplanar and n_boundary_legs < 4:
        warnings.append(
            f"NOTE: labeled nonplanar/crossed and reads planar, but only "
            f"{n_boundary_legs} boundary external legs -> the leg-order planarity "
            "test cannot discriminate planar vs crossed here (3-pt/2-pt closure is "
            "always planar). Verdict INCONCLUSIVE from topology alone.")
    # (b2) the planarity report by name (planarity's own keys): the closure
    #      order fp["planar"] was computed in and where it came from, the
    #      non-discriminating flag when fewer than four legs are attached (a
    #      two- or three-leg closure is a label check, not a planar/crossed
    #      verdict), the legs-joined-at-infinity criterion beside it, the
    #      enumerated leg orders, and the record's Mandelstam declaration as a
    #      REPORTED candidate order (not enforced: a record's channel naming
    #      need not follow the drawn boundary order).  On an inferred leg set
    #      (block (d) withholds the verdict) the verdict-class fields are
    #      withheld here too; the computed values stay under
    #      fp["planarity_report"].
    prep = planarity(fam)
    inferred_set = (fam.kinematics_source
                    or (fam.record_notes.get("leg_set") or {}).get("source")) == "inferred"
    fp["planarity_discriminating"] = prep.get("discriminating")
    fp["planarity_closure_order"] = prep.get("closure_order")
    fp["planarity_closure_order_source"] = prep.get("closure_order_source")
    fp["planar_legs_joined_at_infinity"] = (None if inferred_set
                                            else prep.get("planar_legs_joined_at_infinity"))
    report = dict(prep)
    if inferred_set:
        report["withheld"] = ("LEG-SET-INFERRED: the leg set these values were computed on "
                              "was inferred from the propagator symbols; no planarity "
                              "verdict is given on it (see planarity_on_inferred_leg_set)")
    fp["planarity_report"] = report
    if not inferred_set and prep.get("planar") is not None:
        if prep.get("discriminating") is False:
            warnings.append(
                "PLANARITY-NOT-DISCRIMINATING: %d external leg(s) attached (N < 4): the "
                "leg-order closure test cannot separate a planar leg order from a non-planar "
                "one with fewer than four legs on the boundary, so planar=%s is a label check "
                "on the realized graph, not a leg-order verdict (abstract planar %s; legs joined "
                "at one vertex at infinity planar %s -- the Feynman-graph criterion, reported "
                "beside it)"
                % (prep.get("n_external_legs"), prep["planar"], prep.get("abstract_planar"),
                   prep.get("planar_legs_joined_at_infinity")))
        md = prep.get("mandelstam_declaration") or {}
        if md.get("cyclic_order") and md.get("agrees_with_closure_order") is False \
                and md.get("planar_in_declared_order") is not prep["planar"]:
            warnings.append(
                "PLANARITY-LEG-ORDER: the record's Mandelstam declaration (%s; %s) puts the "
                "legs in the cyclic order %s, where the graph reads planar=%s; the verdict "
                "above (planar=%s) closes them in %s %s. A record's channel names do not fix "
                "the drawn boundary order (t = (p1+p3)^2 names the diagonal channel of a box "
                "drawn p1,p2,p3,p4 as well as the adjacent channel of one drawn p1,p2,p4,p3), "
                "so neither order is enforced: give the drawn order (planarity(fam, "
                "leg_order=...)) to decide"
                % (", ".join("%s = (%s)^2" % (k, "+".join(v[0])) for k, v in
                             sorted(md.get("channels", {}).items())),
                   md.get("source"), md["cyclic_order"], md["planar_in_declared_order"],
                   prep["planar"],
                   "the canonical index order" if not prep.get("closure_order_source", "").startswith(
                       ("leg_order", "family.")) else "the declared order",
                   prep.get("closure_order")))

    # (c) cut-signature symmetry note (informational; planarity is the rigorous
    #     discriminator -- the t<->u cut symmetry is routing-dependent and only a
    #     hint, so it is reported, not used to raise an independent warning).
    t = fp["cut_signature"].get("t")
    u = fp["cut_signature"].get("u")
    if t is not None and u is not None:
        fp["cut_symmetry"] = ("t<->u symmetric" if t == u
                              else "t/u asymmetric")

    # (d) the external-leg set: where it came from, how it was completed, and
    #     the LEG-SET-INFERRED consequence -- an inferred set (the Kira-yaml
    #     fallback) gets no planarity verdict; the numbers computed on it are
    #     kept under a named key.
    leg_set = dict(fam.record_notes.get("leg_set") or {})
    # a Family built in memory with exts + ext_subs (the catalog, a control)
    # declares its legs by construction; only the Kira-yaml fallback infers
    leg_set["source"] = fam.kinematics_source or leg_set.get("source") or "constructor:exts+ext_subs"
    leg_set["legs"] = list(fam.exts)
    leg_set["ext_subs"] = dict(fam.ext_subs)
    leg_set["virtualities"] = dict(fam.leg_virt)
    declared = list(leg_set.get("declared_legs", fam.exts))
    if leg_set.get("declared_canonical_hash") and (
            declared != list(fam.exts) or leg_set.get("propagator_rewrite")):
        # the record's own basis (before completion / normalization): the hash
        # the CLI printed on this file before the leg set was completed
        fp["canonical_hash_declared_legs"] = leg_set["declared_canonical_hash"]
        fp["declared_legs"] = declared
    if len(fam.exts) == 2:
        fp["two_point_cut_lines_between_the_legs"] = two_point_cut(fam)
    orient = (leg_set.get("orientation") or {})
    if orient.get("resolved_from_propagators") and orient.get("outgoing"):
        warnings.append("LEG-ORIENTATION-RESOLVED: declared leg(s) %s flow OUT (the top sector "
                        "conserves momentum at every vertex only with that orientation); the "
                        "record's kinematics declares no conservation rule, the implicit leg "
                        "%s was completed with it" % (orient["outgoing"], leg_set.get("implicit_leg")))
    elif orient.get("patterns_tried") and not orient.get("resolved_from_propagators"):
        warnings.append("LEG-COMPLETION-UNREALIZED: no orientation of the declared legs %s "
                        "realizes the top sector with the implicit leg %s; the all-incoming "
                        "completion is kept and the realizer's verdict names the failure"
                        % (leg_set.get("declared_legs"), leg_set.get("implicit_leg")))
    if leg_set.get("implicit_leg"):
        (imp, expr), = leg_set["implicit_leg"].items()
        warnings.append("LEG-SET-COMPLETED: the record declares the independent leg(s) %s and "
                        "no conservation rule; the graph's remaining leg %s = %s (all legs "
                        "incoming) was completed from momentum conservation%s"
                        % (leg_set.get("declared_legs"), imp, expr,
                           " with %s read as outgoing (see LEG-ORIENTATION-RESOLVED)"
                           % orient.get("outgoing") if orient.get("outgoing") else ""))
    if leg_set["source"] == "inferred":
        note = leg_set.get("warning") or _leg_set_inferred_note(
            leg_set.get("declared_legs", fam.exts), leg_set.get("conservation_declared", {}))
        warnings.append(note)
        fp["planarity_on_inferred_leg_set"] = {"planar": fp["planar"],
                                               "method": fp["planarity_method"]}
        fp["planar"] = None
        fp["planarity_method"] = ("withheld: LEG-SET-INFERRED (the external-leg set was "
                                  "inferred from the propagator symbols; give the record's "
                                  "kinematics for a verdict)")
    fp["leg_set"] = leg_set

    fp["warnings"] = warnings
    return fp


# ---------------------------------------------------------------------------
#  5b.  COMPARE MODE  --  a family of record against a drawn graph: the
#       identity census's verdict (IDENTITY-PASS / FAMILY-MISMATCH /
#       NON-GRAPH / NOT-CHECKABLE) with its evidence, from the units above.
# ---------------------------------------------------------------------------
#  The verdict emitter compare_drawn(famA, famB, kinA, kinB):
#    NOT-CHECKABLE   a side cannot be read, its legs cannot be completed (a
#                    leg set that does not sum to zero, a leg the propagators
#                    never spell), a realization that exhausts its budget, or
#                    -- with three or more legs -- a side that declares no
#                    virtuality class for some leg (the identity test needs
#                    the labels; nothing is inferred).  The side and the
#                    reason are named.
#    NON-GRAPH       either side fails non_graph_test (3V <= 2E + N, the loop
#                    rank, distinct lines) or admits no connected cover with
#                    V = E - L + 1 (build_graph): the realizer's units.
#    IDENTITY-PASS   the realized labeled graphs are isomorphic WITH masses,
#                    multiplicities and leg virtuality classes
#                    (realized_graph_iso "with_masses_and_legs"), the legs
#                    mapped by class -- the leg map is printed with, when the
#                    two propagator sets differ by a loop redefinition, the
#                    affine describer (affine_iso_search).  Two legs (N = 2)
#                    carry one p^2 by conservation, so an undeclared class
#                    there is vacuous and the mass level decides, said by name.
#    FAMILY-MISMATCH otherwise -- the FIRST failing level named (bare /
#                    masses / legs) with the legs-per-vertex lists of both
#                    sides (legs_per_vertex) and V, E of both.
#  Symbolic mass labels are IDENTIFIED across the two sides before the
#  labeled compare (a drawing's m2 against a record's resolved 1; the
#  census's mass_map): every count-preserving bijection of the symbolic
#  labels is tried and the one that carries the labels is printed; numeric
#  labels are values and never identified with one another.
#  CYCLIC ORDER: when BOTH sides declare a cyclic leg order (a drawn graph's
#  boundary order; Family.cyclic_leg_order) the leg map is checked to preserve
#  it up to rotation and reflection and a failure is the named sub-finding
#  "CYCLIC-ORDER MISMATCH".  The sub-finding does NOT change the verdict:
#  that is the census's rule (census_lib_r21.labeled_isomorphism reports
#  "leg_perm_cyclic_order_preserved" and census_specs_r21.make_verdict never
#  reads it; rows_01_02_04_05_07_census.std_verdict: "the checker's
#  canonical-order planarity flag is reported but not used here: it assumes
#  the family's legs are labelled in the drawn cyclic order, which a family
#  in Smirnov's convention t=(p1+p3)^2 is not" -- row 2 is IDENTITY-PASS
#  under the leg map {p2: p1, p1: p2, p3: p3, p4: p4}).  A record's leg index
#  order is a convention, not a declared cyclic order, so a family that
#  declares none gets the drawn order reported in its own leg names instead.

VERDICT_STRINGS = ("IDENTITY-PASS", "FAMILY-MISMATCH", "NON-GRAPH", "NOT-CHECKABLE")
CYCLIC_ORDER_MISMATCH = "CYCLIC-ORDER MISMATCH"


def _parse_dict_literal(s):
    """A JSON or Python dict literal of strings -> {str: str}; None when it
    is neither."""
    s = str(s).strip()
    try:
        d = json.loads(s)
        if isinstance(d, dict):
            return {str(k): str(v) for k, v in d.items()}
    except Exception:
        pass
    pairs = re.findall(r'''['"]?([A-Za-z_]\w*)['"]?\s*:\s*['"]([^'"]*)['"]''', s)
    if pairs:
        return {k: v for k, v in pairs}
    return None


def _parse_list_literal(s):
    return [x.strip().strip('"\'') for x in str(s).strip().strip("[]").split(",")
            if x.strip()]


def read_drawn_header(path):
    """
    The leg data a routed drawn-graph yaml carries in its HEADER COMMENTS
    (the identity census's four forms), for the compare mode's drawn side --
    a routed yaml whose family block declares no external_momenta key is
    otherwise read with an inferred leg set.  Returns None when the file
    carries none of the forms, else a dict:
      "form"             which form was read (named below)
      "legs"             the external legs in the header's order (all incoming)
      "conservation"     {dependent leg: expression}
      "virtualities"     {leg: p^2 string} for the legs the header labels
      "identity_labels"  {leg: "on" | "off"} when the header carries them
      "cyclic_leg_order" the drawn boundary order when stated
      "nu"               the per-line multiplicity list when stated
    Forms (a line of each, verbatim shapes):
      "external legs"   # external legs (all incoming): p1, p2, p3; p3 = -(p1+p2)
                        # leg virtualities: {"p1": "0", "p2": "0", "p3": "s"}
      "legs@vertex"     # legs: p1@TL (p^2=0), p2@TR (p^2=0), ...
                        # loop momenta [...]; {'p4': '-p1-p2-p3'}; cyclic leg order as drawn: [...]
                        # ...; edge multiplicity (nu) = [1, 1, 1, 2, ...]
      "leg dictionary"  # leg labels (identity): {"p1": "on", ...}
                        # leg virtualities (kinematics): {"p1": "p56^2 = -1", ...}
                        # convention: all legs incoming, p_N = -(p_1 + ... + p_{N-1})
      "External legs"   # External legs (all incoming; last = -(sum of the others)): p at L p^2=pp; pout at R p^2=pp
    A two-point routing note ("p in, p out") declares no legs: None.
    """
    try:
        txt = open(path).read()
    except Exception:
        return None
    out = {"form": None, "legs": None, "conservation": {}, "virtualities": {},
           "identity_labels": {}, "cyclic_leg_order": None, "nu": None}
    mnu = re.search(r'edge multiplicity \(nu\)\s*=\s*\[([^\]]*)\]', txt)
    if mnu:
        try:
            out["nu"] = [int(x) for x in _parse_list_literal(mnu.group(1))]
        except ValueError:
            out["nu"] = None
    mo = re.search(r'cyclic leg order[^\[\n]*\[([^\]]*)\]', txt)
    if mo:
        out["cyclic_leg_order"] = _parse_list_literal(mo.group(1))
    # form "external legs"
    m = re.search(r'(?m)^#\s*external legs \(all incoming\)\s*:\s*(.+)$', txt)
    if m:
        parts = [p.strip() for p in m.group(1).split(";")]
        legs = [x.strip() for x in parts[0].split(",") if x.strip()]
        cons = {}
        for p in parts[1:]:
            mc = re.fullmatch(r'(\w+)\s*=\s*(.+)', p.strip())
            if mc:
                cons[mc.group(1)] = mc.group(2).strip()
        mv = re.search(r'(?m)^#\s*leg virtualities\s*:\s*(\{.*\})\s*$', txt)
        virt = _parse_dict_literal(mv.group(1)) if mv else None
        out.update({"form": "external legs", "legs": legs, "conservation": cons,
                    "virtualities": {k: v for k, v in (virt or {}).items() if k in legs}})
        return out
    # form "legs@vertex"
    m = re.search(r'(?m)^#\s*legs\s*:\s*(.+)$', txt)
    if m and "@" in m.group(1):
        legs, virt = [], {}
        for item in m.group(1).split(","):
            mi = re.match(r'\s*(\w+)\s*@\s*(\w+)(?:\s*\(\s*p\^2\s*=\s*([^)]*)\))?', item)
            if not mi:
                continue
            legs.append(mi.group(1))
            if mi.group(3) is not None:
                virt[mi.group(1)] = mi.group(3).strip()
        cons = {}
        md = re.search(r'(?m)^#\s*loop momenta[^\n]*?(\{[^}]*\})', txt)
        if md:
            cons = _parse_dict_literal(md.group(1)) or {}
        if virt:
            out.update({"form": "legs@vertex", "legs": legs, "conservation": cons,
                        "virtualities": virt})
            return out
    # form "leg dictionary"
    mv = re.search(r'(?m)^#\s*leg virtualities \(kinematics\)\s*:\s*(\{.*\})\s*$', txt)
    if mv:
        raw = _parse_dict_literal(mv.group(1)) or {}
        legs = list(raw)
        virt = {}
        for k, v in raw.items():
            virt[k] = (v.split("=", 1)[1] if "=" in v else v).strip()
        mi = re.search(r'(?m)^#\s*leg labels \(identity\)\s*:\s*(\{.*\})\s*$', txt)
        ident = _parse_dict_literal(mi.group(1)) if mi else {}
        cons = {}
        if len(legs) >= 2:
            cons[legs[-1]] = "-(" + "+".join(legs[:-1]) + ")"
        out.update({"form": "leg dictionary", "legs": legs, "conservation": cons,
                    "virtualities": virt, "identity_labels": ident or {}})
        return out
    # form "External legs"
    m = re.search(r'(?m)^#\s*External legs \(all incoming; last = -\(sum of the others\)\)\s*:\s*(.+)$',
                  txt)
    if m:
        legs, virt = [], {}
        for item in m.group(1).split(";"):
            mi = re.match(r'\s*(\w+)\s+at\s+\w+\s+p\^2\s*=\s*(\S+)', item)
            if mi:
                legs.append(mi.group(1))
                virt[mi.group(1)] = mi.group(2).strip()
        if legs:
            cons = {legs[-1]: "-(" + "+".join(legs[:-1]) + ")"} if len(legs) >= 2 else {}
            out.update({"form": "External legs", "legs": legs, "conservation": cons,
                        "virtualities": virt})
            return out
    if out["nu"] is not None or out["cyclic_leg_order"] is not None:
        out["form"] = "nu / cyclic order only"
        return out
    return None


def _jl_replacement_virtualities(path, exts, ext_subs):
    """The kinematics Dict of an AmflowFamily(...) record (its fourth
    argument: "p1^2" => "0", "(p1 + p2)^2" => "s", ...) -> {leg: p^2} for
    every leg the rules determine (through _kin_from_replacement and
    leg_virtualities); {} when the file carries no such Dict."""
    try:
        txt = open(path).read()
    except Exception:
        return {}
    i = txt.find("AmflowFamily(")
    if i < 0:
        return {}
    j = i + len("AmflowFamily(")
    depth, k = 1, j
    while k < len(txt) and depth:
        if txt[k] == '(':
            depth += 1
        elif txt[k] == ')':
            depth -= 1
        k += 1
    args = _split_top_commas(txt[j:k - 1])
    if len(args) < 5:
        return {}
    rep = _jl_dict(args[4])
    if not rep:
        return {}
    kin = _kin_from_replacement(rep, exts)
    out = leg_virtualities(kin, legs=list(exts), ext_subs=dict(ext_subs))
    # a leg's square the Dict states DIRECTLY ("p4^2" => "0") is the record's
    # own label for it, before any expansion through the other rules
    for k, v in rep.items():
        m = re.fullmatch(r'\(?([A-Za-z_]\w*)\)?\^2', str(k).replace(" ", ""))
        if m and m.group(1) in exts:
            out[m.group(1)] = str(sp.nsimplify(sp.sympify(str(v))))
    return out


class UnknownFamilyName(ValueError):
    """A family name the file does not carry (load_compare_side with `name`,
    the CLI's --name): refused by name, the message listing the families the
    file carries.  Raised through compare_drawn, never folded into a
    NOT-CHECKABLE verdict -- the record is readable, the caller's name is not
    in it."""


def load_compare_side(path, kinematics=None, name=None, integral_index=None):
    """
    One side of a compare from a file: load_family(...) and then, ONLY where
    the record itself left the leg set to inference, the data the file
    carries elsewhere -- a routed drawn yaml's header comments
    (read_drawn_header: legs, conservation, virtualities, cyclic order,
    multiplicities), the AmflowFamily kinematics Dict of a .jl
    (_jl_replacement_virtualities) -- so that the compare runs on declared
    legs and labels or says it cannot.  A kinematics file, when given, wins
    (load_family's rule).  With `name` the family of that name is the side
    (a multi-family file); a name the file does not carry raises
    UnknownFamilyName listing the names it does; without a name the file's
    first family is the side.  Returns (Family, notes) with notes naming what
    was applied: "leg_set_source", "header", "header_applied",
    "leg_virtualities_source", "nu_source", "cyclic_leg_order_source", and
    "family_selected_by" ("name" / "first family in the file") beside
    "family_name" / "n_families_in_file".
    """
    fams_all = load_family(path, None, kinematics=kinematics, integral_index=integral_index)
    if not fams_all:
        raise ValueError("no family found in " + path)
    if name:
        fams = [f for f in fams_all if f.name == name]
        if not fams:
            raise UnknownFamilyName("no family found named %r in %s; the file carries: %s"
                                    % (name, path, ", ".join(f.name for f in fams_all)))
    else:
        fams = fams_all
    fam = fams[0]
    notes = {"file": os.path.basename(path), "n_families_in_file": len(fams_all),
             "family_name": fam.name,
             "family_selected_by": ("name" if name else "first family in the file"),
             "kinematics": os.path.basename(str(kinematics)) if kinematics else None,
             "leg_set_source": fam.kinematics_source, "header": None,
             "header_applied": [], "leg_virtualities_source": None,
             "nu_source": "record" if fam.nu is not None else None,
             "cyclic_leg_order_source": ("record" if fam.cyclic_leg_order else None)}
    if fam.leg_virt:
        notes["leg_virtualities_source"] = fam.kinematics_source
    hdr = read_drawn_header(path) if path.lower().endswith((".yaml", ".yml")) else None
    if hdr:
        notes["header"] = hdr["form"]
        if hdr.get("legs") and fam.kinematics_source == "inferred" and kinematics is None:
            src = "drawn-header:" + str(hdr["form"])
            legs, cons = list(hdr["legs"]), dict(hdr["conservation"])
            # a routing that spells the header's dependent leg and leaves an
            # independent one unspelled (row 21's k - p4 - b with p4 = -p1-p2-p3)
            # declares the same conservation law in another basis: re-solve the
            # rule for the unspelled leg, never guess a leg set
            spelled = set()
            for expr, _m in fam.propagators:
                spelled |= set(re.findall(r'[A-Za-z_]\w*', expr))
            unspelled = [e for e in legs if e not in cons and e not in spelled]
            if len(unspelled) == 1 and len(cons) == 1:
                (dep, rule), = cons.items()
                if dep in spelled:
                    syms = {e: sp.Symbol(e) for e in legs}
                    try:
                        sol = sp.solve(sp.Eq(syms[dep], sp.sympify(rule, locals=syms)),
                                       syms[unspelled[0]])
                    except Exception:
                        sol = []
                    if len(sol) == 1:
                        cons = {unspelled[0]: str(sp.expand(sol[0]))}
                        notes["conservation_resolved"] = (
                            "the header eliminates %s (%s = %s) but the routing spells %s and not "
                            "%s; the same rule re-solved: %s = %s"
                            % (dep, dep, rule, dep, unspelled[0], unspelled[0], cons[unspelled[0]]))
            new = Family(fam.name, fam.loops, legs, cons, fam.propagators,
                         fam.source, fam.physical, nu=fam.nu, leg_virt=hdr["virtualities"],
                         mass_values=fam.mass_values, kinematics_source=src,
                         cyclic_leg_order=(hdr.get("cyclic_leg_order") or fam.cyclic_leg_order),
                         source_kind=fam.source_kind,
                         record_notes={"leg_set": {"source": src,
                                                   "declared_legs": list(hdr["legs"]),
                                                   "conservation_declared": dict(hdr["conservation"]),
                                                   "header_form": hdr["form"]}})
            new.kinematics_path = fam.kinematics_path
            complete_legs(new, None)
            fam = new
            notes["leg_set_source"] = src
            notes["header_applied"].append("legs")
            if hdr["virtualities"]:
                notes["leg_virtualities_source"] = src
                notes["header_applied"].append("virtualities")
        elif hdr.get("virtualities") and not fam.leg_virt:
            fam.leg_virt = {k: v for k, v in hdr["virtualities"].items() if k in fam.exts}
            notes["leg_virtualities_source"] = "drawn-header:" + str(hdr["form"])
            notes["header_applied"].append("virtualities")
        if hdr.get("nu") is not None and fam.nu is None:
            if len(hdr["nu"]) == len(fam.propagators):
                fam.nu = list(hdr["nu"])
                notes["nu_source"] = "drawn-header: edge multiplicity (nu)"
                notes["header_applied"].append("nu")
            else:
                notes["nu_header_ignored"] = ("header nu list has %d entries, the record %d "
                                              "propagators" % (len(hdr["nu"]), len(fam.propagators)))
        if hdr.get("cyclic_leg_order") and not fam.cyclic_leg_order:
            fam.cyclic_leg_order = list(hdr["cyclic_leg_order"])
            notes["cyclic_leg_order_source"] = "drawn-header: cyclic leg order as drawn"
            notes["header_applied"].append("cyclic_leg_order")
    if fam.source_kind == "amflow_jl" and not fam.leg_virt:
        v = _jl_replacement_virtualities(path, fam.exts, fam.ext_subs)
        if v:
            fam.leg_virt = v
            notes["leg_virtualities_source"] = "jl:kinematics Dict (replacement rules)"
    # a kinematics file that lists its dependent leg under outgoing_momenta
    # (momentum_conservation [p3, p1+p2] with p3 the outgoing Higgs) states
    # that leg's OUTGOING momentum; the kinematics reader keeps the rule as
    # written (the in/out lists are not an orientation datum for the symbols)
    # and the all-incoming leg set then does not close.  The file's own
    # outgoing declaration is applied here, the tool's way (_rewrite_outgoing:
    # the leg symbol names the incoming momentum, = minus the rule), ONLY when
    # the declared reading does not close and the declared flip closes it.
    ls = fam.record_notes.get("leg_set") or {}
    dep_out = [e for e in (ls.get("declared_outgoing") or []) if e in fam.ext_subs]
    if kinematics is not None and dep_out and len(dep_out) == 1 and fam.exts:
        syms = _symbols(fam.loops, fam.exts)
        ext_resolved = _resolve_ext_subs(fam.exts, fam.ext_subs, syms)
        s_in = sp.expand(sum((ext_resolved[syms[e]] for e in fam.exts), sp.Integer(0)))
        if s_in != 0:
            x = dep_out[0]
            flipped = dict(fam.ext_subs)
            flipped[x] = "-(%s)" % fam.ext_subs[x]
            ext2 = _resolve_ext_subs(fam.exts, flipped, syms)
            s_flip = sp.expand(sum((ext2[syms[e]] for e in fam.exts), sp.Integer(0)))
            if s_flip == 0:
                props, rew = _rewrite_outgoing(fam.propagators, [x])
                fam.propagators = props
                fam.ext_subs = flipped
                kin = read_kinematics_yaml(kinematics)
                fam.leg_virt = leg_virtualities(kin, legs=list(fam.exts), ext_subs=flipped)
                ls["orientation"] = {"resolved_from_propagators": False, "from_kinematics_outgoing": True,
                                     "outgoing": [x],
                                     "note": ("the kinematics file lists %s under outgoing_momenta and "
                                              "names its outgoing momentum in momentum_conservation; "
                                              "the leg symbol now names the incoming momentum %s"
                                              % (x, flipped[x]))}
                ls["propagator_rewrite"] = rew
                fam.record_notes["leg_set"] = ls
                notes["outgoing_leg_applied"] = ls["orientation"]["note"]
    return fam, notes


def legs_per_vertex(fam):
    """
    The number of external legs attached to each realized vertex, sorted
    descending (the census's "photons per vertex": the retired row-27
    three-point family reads [2, 1, 1, 0] against the drawn box's
    [1, 1, 1, 1]), with V (internal vertices), E (lines), N (legs) and the
    realizer's cause when the family does not realize.  Counted on
    build_graph's first realization; the realizer's own count
    ("realizer_legs_per_vertex", the realizer leg) is carried beside it.
    """
    out = {"ok": False, "legs_per_vertex": None, "V": None, "E": None, "N": len(fam.exts),
           "cause": None, "n_distinct_connected_realizations": None,
           "realizer_legs_per_vertex": None}
    if not _HAVE_NX:
        out["cause"] = "networkx unavailable"
        return out
    G, ok = build_graph(fam)
    rep = (G.graph.get("realization", {}) if G is not None else {}) or {}
    out["cause"] = rep.get("cause")
    out["n_distinct_connected_realizations"] = rep.get("n_distinct_connected_realizations")
    out["realizer_legs_per_vertex"] = rep.get("legs_per_vertex")
    out["V"], out["E"] = rep.get("V"), rep.get("E")
    if not ok or G is None:
        return out
    counts = []
    n_lines = 0
    for n in G.nodes():
        if str(n).startswith("ext_"):
            continue
        counts.append(sum(1 for _u, _v, k in G.edges(n, keys=True) if str(k).startswith("ext_")))
    for _u, _v, k in G.edges(keys=True):
        if not str(k).startswith("ext_"):
            n_lines += 1
    out.update({"ok": True, "legs_per_vertex": sorted(counts, reverse=True),
                "V": len(counts), "E": n_lines})
    return out


def _is_numeric_label(lbl):
    try:
        return bool(sp.nsimplify(str(lbl)).is_number)
    except Exception:
        return False


def _mass_identifications(LA, LB):
    """
    Candidate identifications of the two sides' non-zero mass labels
    (labeled_set outputs): [(mass_mapA, mass_mapB, description)], the
    identity first when the label multisets already agree.  A SYMBOLIC label
    on one side may be identified with a label of the same line count on the
    other (the map is applied to the symbolic side); two different NUMERIC
    labels are values and are never identified.  Empty when no
    count-preserving identification exists.
    """
    import collections
    cA = collections.Counter(m for m in LA["masses"] if m != "0")
    cB = collections.Counter(m for m in LB["masses"] if m != "0")
    out = []
    if cA == cB:
        out.append(({}, {}, "identical mass labels"))
    if sorted(cA.values()) != sorted(cB.values()):
        return out
    labsA = sorted(cA)
    for perm in itertools.permutations(sorted(cB)):
        pairs = list(zip(labsA, perm))
        if any(cA[a] != cB[b] for a, b in pairs):
            continue
        mA, mB, desc, ok = {}, {}, [], True
        for a, b in pairs:
            if a == b:
                continue
            if not _is_numeric_label(b):
                mB[b] = a
                desc.append("%s (drawn) = %s (family)" % (b, a))
            elif not _is_numeric_label(a):
                mA[a] = b
                desc.append("%s (family) = %s (drawn)" % (a, b))
            else:
                ok = False
                break
        if ok and (mA or mB):
            out.append((mA, mB, "; ".join(desc)))
    return out


def _cyclic_order_finding(famA, famB, leg_map):
    """The cyclic-order sub-finding of a compare (see the section header):
    {"A_declared", "B_declared", "preserved", "finding", "drawn_order_in_A_names",
    "note"}; preserved is None when a side declares no order or no leg map
    exists."""
    oa = list(famA.cyclic_leg_order) if famA.cyclic_leg_order else None
    ob = list(famB.cyclic_leg_order) if famB.cyclic_leg_order else None
    out = {"A_declared": oa, "B_declared": ob, "leg_map": dict(leg_map) if leg_map else None,
           "preserved": None, "finding": None, "drawn_order_in_A_names": None, "note": None,
           "rule": ("a named sub-finding, reported beside the verdict and never changing it "
                    "(the census's rule: leg_perm_cyclic_order_preserved is reported, the "
                    "verdict functions never read it)")}
    if not leg_map:
        out["note"] = "no leg map (the sides are not isomorphic at the level the map needs)"
        return out
    inv = {b: a for a, b in leg_map.items()}
    if ob and all(b in inv for b in ob):
        out["drawn_order_in_A_names"] = [inv[b] for b in ob]
    if oa and ob:
        if sorted(oa) != sorted(famA.exts) or sorted(ob) != sorted(famB.exts):
            out["note"] = "a declared cyclic order does not name every leg of its side"
            return out
        image = [leg_map[a] for a in oa]
        out["preserved"] = bool(_same_cyclic_order(image, ob))
        if not out["preserved"]:
            out["finding"] = ("%s: the leg map carries the family's declared cyclic order %s to "
                              "%s, not a rotation or reflection of the drawn order %s"
                              % (CYCLIC_ORDER_MISMATCH, oa, image, ob))
    elif ob:
        out["note"] = ("the family declares no cyclic leg order (a record's leg index order is a "
                       "convention, not a datum); the drawn order %s reads %s in the family's "
                       "leg names" % (ob, out["drawn_order_in_A_names"]))
    elif oa:
        out["note"] = "the drawn side declares no cyclic leg order"
    else:
        out["note"] = "neither side declares a cyclic leg order"
    return out


def _side_report(fam, notes, tag):
    L = labeled_set(fam)
    vecs, _m, _i = momentum_vectors(fam)
    ng = non_graph_test(vecs, fam)
    lpv = legs_per_vertex(fam)
    # which integral of the record defined this side: an AMFlow-port JSON lists
    # index vectors (integrals[].indices) and the reader takes one of them for
    # the top sector and the multiplicities -- integral_index (the record's
    # call through --integral, else the reader's default: the integral with
    # the most positive indices, first among ties), n_integrals, the index
    # vector read (integral_indices = Family.nu) and the others by index, so a
    # multiplicity-level mismatch is readable as "the dotted integral" without
    # the source.  A record listing one integral or none carries n_integrals
    # 1 / 0 and no others.
    rn = fam.record_notes or {}
    ints = [list(v) for v in (rn.get("integrals") or [])]
    idx = rn.get("integral_index")
    n_int = int(rn.get("n_integrals") or 0) if rn.get("n_integrals") is not None else len(ints)
    others = [[i, v] for i, v in enumerate(ints) if i != idx]
    idx_source = None
    if n_int:
        given = rn.get("integral_index_given")
        if given is None:
            given = (notes or {}).get("integral_index") is not None
        idx_source = ("--integral N (the record's call)" if given
                      else "reader default: the integral with the most positive indices "
                           "(first among ties)")
    return {"tag": tag, "name": fam.name, "file": (notes or {}).get("file"),
            "loops": list(fam.loops), "legs": list(fam.exts), "ext_subs": dict(fam.ext_subs),
            "leg_virtualities": dict(fam.leg_virt),
            "leg_classes": dict(L["leg_labels"]), "classes_declared": L["leg_classes"] is not None,
            "n_lines": len(vecs), "mass_multiset": L["mass_multiset"], "nu_multiset": L["nu_multiset"],
            "nu_source": L["nu_source"], "cyclic_leg_order": (list(fam.cyclic_leg_order)
                                                              if fam.cyclic_leg_order else None),
            "integral_index": idx, "n_integrals": n_int,
            "integral_indices": (list(fam.nu) if fam.nu is not None else None),
            "other_integrals": others, "integral_index_source": idx_source,
            "non_graph_test": ng, "realization": lpv, "notes": notes or {}}


def compare_drawn(famA, famB, kinA=None, kinB=None, integral_indexA=None, nameA=None):
    """
    THE VERDICT EMITTER (see the section header for the contract).  famA is
    the family of record, famB the drawn graph: Family objects, or file paths
    read by load_compare_side with kinA / kinB as their kinematics files
    (integral_indexA: an AMFlow-port JSON's integrals[N] defining the family's
    top sector and multiplicities, the CLI's --integral; nameA: the family of
    that name in a multi-family famA path, the CLI's --name -- a name the
    file does not carry raises UnknownFamilyName, not a verdict).
    Returns a dict: "verdict" (one of VERDICT_STRINGS), "reason", "line" (the
    one-line verdict), "first_failing_level" (bare / masses / legs / None),
    "levels" (the VF2 booleans), "leg_map", "mass_identification",
    "relabeling" (signed / affine describer or the reason none was found),
    "cyclic_order" (the sub-finding), "planarity_in_drawn_order",
    "legs_per_vertex" / "V" / "E" (both sides), "sides" (the two side
    reports) and "not_checkable" (side, reason) when that is the verdict.
    """
    res = {"verdict": None, "reason": None, "line": None, "first_failing_level": None,
           "levels": None, "leg_map": None, "mass_identification": None, "relabeling": None,
           "cyclic_order": None, "planarity_in_drawn_order": None,
           "legs_per_vertex": {"A": None, "B": None}, "V": {"A": None, "B": None},
           "E": {"A": None, "B": None}, "sides": {}, "not_checkable": None, "leg_class_note": None}

    def not_checkable(side, why):
        res.update({"verdict": "NOT-CHECKABLE", "reason": "%s side: %s" % (side, why),
                    "not_checkable": {"side": side, "reason": why}})
        res["line"] = "VERDICT: NOT-CHECKABLE (%s side: %s)" % (side, why)
        return res

    sides = {}
    for tag, obj, kin in (("family", famA, kinA), ("drawn", famB, kinB)):
        try:
            if isinstance(obj, Family):
                fam, notes = obj, {"file": os.path.basename(str(obj.source)),
                                   "kinematics": os.path.basename(str(kin)) if kin else None,
                                   "leg_set_source": obj.kinematics_source}
            else:
                fam, notes = load_compare_side(str(obj), kinematics=kin,
                                               name=(nameA if tag == "family" else None),
                                               integral_index=(integral_indexA if tag == "family"
                                                               else None))
                if tag == "family" and integral_indexA is not None:
                    notes["integral_index"] = integral_indexA
        except UnknownFamilyName:
            raise
        except Exception as exc:
            return not_checkable(tag, "cannot be read: %s: %s" % (type(exc).__name__, exc))
        sides[tag] = (fam, notes)
    A, nA = sides["family"]
    B, nB = sides["drawn"]
    if not _HAVE_NX:
        return not_checkable("either", NX_SKIP_REASON)
    rep = {"family": _side_report(A, nA, "family"), "drawn": _side_report(B, nB, "drawn")}
    res["sides"] = rep
    for key, tag in (("A", "family"), ("B", "drawn")):
        res["legs_per_vertex"][key] = rep[tag]["realization"]["legs_per_vertex"]
        res["V"][key] = rep[tag]["realization"]["V"]
        res["E"][key] = rep[tag]["realization"]["E"]
    # (1) legs completed / (2) the non-graph tests, either side
    for tag in ("family", "drawn"):
        ng = rep[tag]["non_graph_test"]
        cause = ng.get("cause") or ""
        if cause.startswith("NOT-CHECKABLE"):
            return not_checkable(tag, cause)
    for tag in ("family", "drawn"):
        ng = rep[tag]["non_graph_test"]
        if not ng["passes"]:
            res.update({"verdict": "NON-GRAPH", "reason": "%s side: %s" % (tag, ng["cause"])})
            res["line"] = "VERDICT: NON-GRAPH (%s side: %s)" % (tag, ng["cause"])
            return res
    for tag in ("family", "drawn"):
        rl = rep[tag]["realization"]
        if not rl["ok"]:
            cause = rl.get("cause") or "no realization"
            if str(cause).startswith("NON-GRAPH"):
                res.update({"verdict": "NON-GRAPH", "reason": "%s side: %s" % (tag, cause)})
                res["line"] = "VERDICT: NON-GRAPH (%s side: %s)" % (tag, cause)
                return res
            return not_checkable(tag, "vertex realization did not complete: %s" % cause)
    # (3) the realized-graph isomorphism with the masses identified
    LA, LB = labeled_set(A), labeled_set(B)
    cands = _mass_identifications(LA, LB) or [({}, {}, "no count-preserving identification of the "
                                                       "mass labels exists (mass multisets %s vs %s)"
                                                       % (LA["mass_multiset"], LB["mass_multiset"]))]
    rgi, chosen = None, None
    for mA, mB, desc in cands:
        r = realized_graph_iso(A, B, mass_mapA=mA, mass_mapB=mB)
        if rgi is None:
            rgi, chosen = r, (mA, mB, desc)
        if r.get("with_masses") is True:
            rgi, chosen = r, (mA, mB, desc)
            break
    mA, mB, mdesc = chosen
    res["mass_identification"] = {"mass_mapA": mA, "mass_mapB": mB, "description": mdesc,
                                  "n_candidates": len(cands)}
    res["levels"] = {k: rgi.get(k) for k in ("bare", "with_masses", "with_masses_and_legs",
                                             "with_masses_and_legs_reason",
                                             "with_masses_and_leg_names", "V_A", "E_A", "V_B", "E_B",
                                             "n_realizations", "error")}
    if rgi.get("error"):
        return not_checkable("either", rgi["error"])
    N = len(A.exts)
    declA, declB = rep["family"]["classes_declared"], rep["drawn"]["classes_declared"]
    lpv = "legs per vertex family %s vs drawn %s; V,E family %s,%s vs drawn %s,%s" % (
        res["legs_per_vertex"]["A"], res["legs_per_vertex"]["B"], res["V"]["A"], res["E"]["A"],
        res["V"]["B"], res["E"]["B"])

    def mismatch(level, why):
        res.update({"verdict": "FAMILY-MISMATCH", "first_failing_level": level,
                    "reason": "first failing level: %s -- %s; %s" % (level, why, lpv)})
        res["line"] = "VERDICT: FAMILY-MISMATCH (first failing level: %s -- %s; %s)" % (level, why, lpv)
        return res

    if rgi.get("bare") is not True:
        return mismatch("bare", "the realized graphs are not isomorphic as bare graphs "
                                "(%d vs %d lines, %d vs %d legs)" % (len(LA["vecs"]), len(LB["vecs"]),
                                                                     N, len(B.exts)))
    if rgi.get("with_masses") is not True:
        return mismatch("masses", "isomorphic as bare graphs but no isomorphism carries the line "
                                  "masses and multiplicities (mass multisets family %s vs drawn %s; "
                                  "multiplicities %s vs %s; %s)"
                        % (LA["mass_multiset"], LB["mass_multiset"], LA["nu_multiset"],
                           LB["nu_multiset"], mdesc))
    legs_ok = rgi.get("with_masses_and_legs")
    leg_map = None
    if legs_ok is True:
        leg_map = rgi.get("with_masses_and_legs_leg_map")
    elif legs_ok is False:
        return mismatch("legs", "isomorphic with masses but no isomorphism preserves the leg "
                                "virtuality classes (family %s vs drawn %s)"
                        % (rep["family"]["leg_classes"], rep["drawn"]["leg_classes"]))
    else:
        undeclared = [t for t, d in (("family", declA), ("drawn", declB)) if not d]
        if N == 2:
            leg_map = rgi.get("with_masses_leg_map")
            res["leg_class_note"] = ("leg virtuality classes not declared on the %s side; N = 2: "
                                     "both legs carry one p^2 by momentum conservation, so the "
                                     "class level has nothing to distinguish and the mass level "
                                     "decides (said by name, nothing inferred)" % " and ".join(undeclared))
        else:
            missing = {t: [e for e, c in rep[t]["leg_classes"].items() if c is None]
                       for t in undeclared}
            return not_checkable(" and ".join(undeclared),
                                 "leg virtuality classes undeclared for legs %s with N = %d: the "
                                 "identity test needs the labels (give the kinematics or the "
                                 "drawn legs' p^2); the graphs are isomorphic with masses under "
                                 "the leg map %s"
                                 % (missing, N, rgi.get("with_masses_leg_map")))
    res["leg_map"] = dict(leg_map) if leg_map else None
    # (4) the describer: signed permutation, else the affine class
    try:
        aff = affine_iso_search(A, B, mass_mapA=mA, mass_mapB=mB,
                                respect_leg_classes=bool(declA and declB), with_mass=True)
    except Exception as exc:                                              # pragma: no cover
        aff = {"found": False, "reason": "%s: %s" % (type(exc).__name__, exc)}
    rel = {"signed_permutation": aff.get("signed_permutation_match"),
           "affine_found": bool(aff.get("found")),
           "affine": (aff.get("loop_relabeling") if aff.get("found") else None),
           "affine_leg_map": (aff.get("leg_map") if aff.get("found") else None),
           "loop_map_class": aff.get("loop_map_class"), "bound": aff.get("bound"),
           "reason": aff.get("reason"), "text": None}
    if rel["signed_permutation"]:
        rel["text"] = "signed loop relabeling [%s]" % rel["signed_permutation"]
    elif rel["affine_found"]:
        rel["text"] = ("the propagator sets differ by a loop redefinition: affine map %s with leg map %s"
                       % (rel["affine"], rel["affine_leg_map"]))
    else:
        rel["text"] = ("no relabeling in the searched class (%s) maps the propagator sets; the "
                       "identity is the realized-graph isomorphism's" % aff.get("loop_map_class"))
    res["relabeling"] = rel
    # (5) the cyclic-order sub-finding and planarity in the drawn order
    cyc = _cyclic_order_finding(A, B, leg_map)
    res["cyclic_order"] = cyc
    pl = {}
    try:
        pB = planarity(B)
        pl["drawn"] = {k: pB.get(k) for k in ("planar", "closure_order", "closure_order_source",
                                               "planar_legs_joined_at_infinity", "discriminating")}
        order_A = cyc.get("drawn_order_in_A_names")
        pA = planarity(A, leg_order=order_A) if order_A else planarity(A)
        pl["family"] = {k: pA.get(k) for k in ("planar", "closure_order", "closure_order_source",
                                                "planar_legs_joined_at_infinity", "discriminating")}
    except Exception as exc:                                              # pragma: no cover
        pl["error"] = "%s: %s" % (type(exc).__name__, exc)
    res["planarity_in_drawn_order"] = pl
    parts = ["leg map %s" % res["leg_map"]]
    if mA or mB:
        parts.append("masses identified: %s" % mdesc)
    parts.append(rel["text"])
    if cyc.get("finding"):
        parts.append(cyc["finding"])
    elif cyc.get("preserved") is True:
        parts.append("cyclic order preserved (drawn %s)" % cyc["B_declared"])
    elif cyc.get("drawn_order_in_A_names"):
        parts.append("drawn cyclic order %s = %s in the family's leg names (the family declares none)"
                     % (cyc["B_declared"], cyc["drawn_order_in_A_names"]))
    if res["leg_class_note"]:
        parts.append(res["leg_class_note"])
    res.update({"verdict": "IDENTITY-PASS",
                "reason": "realized graphs isomorphic with masses, multiplicities and leg classes"})
    res["line"] = "VERDICT: IDENTITY-PASS (%s)" % "; ".join(parts)
    return res


def verdict_lines(res):
    """The compare's printable evidence: the verdict line first, then the
    sides, the realized graphs, the isomorphism levels, the describer, the
    cyclic-order sub-finding and the planarity in the drawn order."""
    out = [res.get("line") or "VERDICT: %s (%s)" % (res.get("verdict"), res.get("reason"))]
    sides = res.get("sides") or {}
    for tag in ("family", "drawn"):
        s = sides.get(tag)
        if not s:
            continue
        rl = s["realization"]
        out.append("   %-6s %s (%s): loops %s, legs %s, %s; leg set from %s; classes %s; masses %s; "
                   "multiplicities %s (%s)"
                   % (tag, s["name"], s["file"], s["loops"], s["legs"], s["ext_subs"],
                      (s["notes"] or {}).get("leg_set_source"), s["leg_classes"], s["mass_multiset"],
                      s["nu_multiset"], s["nu_source"]))
        if (s.get("n_integrals") or 0) > 1:
            # the record lists several integrals: name the one read, its index
            # vector, and the others (a one-integral file prints nothing here)
            out.append("          integral %s of %s: %s; others: %s  (%s)"
                       % (s["integral_index"], s["n_integrals"], s["integral_indices"],
                          ", ".join("%s %s" % (i, v) for i, v in (s.get("other_integrals") or [])),
                          s.get("integral_index_source")))
        out.append("          realized: V=%s E=%s N=%s legs per vertex %s; distinct connected "
                   "realizations %s%s"
                   % (rl["V"], rl["E"], rl["N"], rl["legs_per_vertex"],
                      rl["n_distinct_connected_realizations"],
                      ("; " + str(rl["cause"])) if rl.get("cause") else ""))
        ng = s["non_graph_test"]
        out.append("          non-graph tests: 3V=%s <= 2E+N=%s %s; loop rank %s/%s; distinct lines %s; "
                   "legs sum to zero %s"
                   % (ng["three_V"], ng["two_E_plus_N"], ng["count_test_3V_le_2E_plus_N"],
                      ng["loop_momentum_rank"], ng["L"], ng["distinct_lines_test"],
                      ng["legs_sum_to_zero_test"]))
    lv = res.get("levels")
    if lv:
        out.append("   isomorphism of the realized labeled graphs (VF2): bare %s, with masses and "
                   "multiplicities %s, with leg classes %s%s; by leg name %s"
                   % (lv["bare"], lv["with_masses"], lv["with_masses_and_legs"],
                      (" (%s)" % lv["with_masses_and_legs_reason"]) if lv.get("with_masses_and_legs_reason")
                      else "", lv["with_masses_and_leg_names"]))
    mi = res.get("mass_identification")
    if mi:
        out.append("   mass labels: %s (candidates tried %s)" % (mi["description"], mi["n_candidates"]))
    rel = res.get("relabeling")
    if rel:
        out.append("   relabeling: %s" % rel["text"])
        if rel.get("bound"):
            out.append("   affine search bound: %s" % rel["bound"])
    cyc = res.get("cyclic_order")
    if cyc:
        out.append("   cyclic order: family declares %s, drawn declares %s; preserved %s%s"
                   % (cyc["A_declared"], cyc["B_declared"], cyc["preserved"],
                      ("; " + (cyc["finding"] or cyc["note"] or "")) if (cyc.get("finding") or cyc.get("note")) else ""))
    pl = res.get("planarity_in_drawn_order")
    if pl:
        for tag in ("family", "drawn"):
            p = pl.get(tag)
            if p:
                out.append("   planarity %-6s legs closed in %s (%s): planar %s; legs joined at infinity %s; "
                           "discriminating %s"
                           % (tag, p["closure_order"], p["closure_order_source"], p["planar"],
                              p["planar_legs_joined_at_infinity"], p["discriminating"]))
        if pl.get("error"):
            out.append("   planarity: %s" % pl["error"])
    if res.get("not_checkable"):
        out.append("   not checkable: %s" % res["not_checkable"])
    return out


# ---------------------------------------------------------------------------
#  6.  SELF-TEST  (NP-dbox must flag planar; planar/crossed/sunrise must not
#      false-flag)
# ---------------------------------------------------------------------------
#  self_test() is a dispatcher over SELF_TEST_LEGS.  A leg is a function of no
#  arguments returning {"status": "PASS" | "FAIL" | "SKIP", "detail": str,
#  "results": {...}}; a leg whose fixture is absent returns SKIP with the
#  fixture named.  Legs register with register_self_test_leg(id, title, fn);
#  a new battery appends its own leg function and one registration line and
#  never edits self_test.  The summary counts "N executed / M skipped" and
#  prints ALL PASS only when nothing failed, nothing was skipped AND no
#  expected-fail sub-leg is open; the exit code is 0 for ALL PASS (or with
#  declared expected-fails only), 1 when any leg FAILED, 2 when any leg was
#  SKIPPED (a skip is an unknown, never a pass).

SELF_TEST_LEGS = []

_FIXTURE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            "tests", "fixtures")


def register_self_test_leg(leg_id, title, fn):
    """Append a self-test leg (id, title, function) to SELF_TEST_LEGS."""
    SELF_TEST_LEGS.append((str(leg_id), str(title), fn))
    return fn


LEGACY_ENV_PREFIX = "TOPOLOGY_AUDIT_"   # the package's former name; still read
ENV_PREFIX = "DOGTAG_"                  # the primary override names


def _env_names(env_var):
    """The override names read for one fixture, primary first: DOGTAG_<X>,
    then the legacy TOPOLOGY_AUDIT_<X> (the package was renamed from
    topology-audit to Dogtag; both spellings keep working).  A name passed
    with either prefix, or with none, resolves to the same pair."""
    stem = env_var
    for pre in (ENV_PREFIX, LEGACY_ENV_PREFIX):
        if stem.startswith(pre):
            stem = stem[len(pre):]
            break
    return ENV_PREFIX + stem, LEGACY_ENV_PREFIX + stem


def _fixture_env(env_var):
    """(name, value) of the first override that is set among _env_names(),
    or (primary_name, "") when neither is."""
    names = _env_names(env_var)
    for name in names:
        val = os.environ.get(name, "")
        if val:
            return name, val
    return names[0], ""


def _fixture_or_skip(env_var, vendored_rel, what):
    """Resolve a fixture: the env override if set (DOGTAG_* first, then the
    legacy TOPOLOGY_AUDIT_* name), else the vendored copy under
    tests/fixtures; returns (path, None) or (None, skip_reason)."""
    name, override = _fixture_env(env_var)
    if override:
        if os.path.exists(override):
            return override, None
        return None, (f"{what}: {name}={override} does not exist")
    vendored = os.path.join(_FIXTURE_DIR, vendored_rel) if vendored_rel else ""
    if vendored and os.path.exists(vendored):
        return vendored, None
    named = (f"vendored copy tests/fixtures/{vendored_rel} absent"
             if vendored_rel else "no vendored copy")
    primary, legacy = _env_names(env_var)
    return None, f"{what}: fixture absent ({named}; set {primary} (or {legacy}) to run)"


NX_SKIP_REASON = "networkx unavailable (declared dependency; pip install networkx)"


def _require_nx(what=""):
    """Dependency gate for a self-test leg or sub-leg that needs networkx (the
    realizer, planarity, the cut signature, VF2, the realized hash): returns
    None when networkx imported, else the skip reason naming the dependency.
    The leg returns SKIP with it -- the missing-fixture pattern: counted,
    named on the summary line, rc 2 -- never PASS, never FAIL.  The tool
    itself still runs without networkx (every value that needs the realized
    graph is None with the reason); only the battery's verdicts are withheld."""
    if _HAVE_NX:
        return None
    return (f"{what}: " if what else "") + NX_SKIP_REASON


def _selftest_npdbox():
    # 1. THE NP-DBOX (the mislabel).  Must flag = planar Smirnov.
    #    Vendored fixture tests/fixtures/npdbox/family.jl (the historically
    #    mislabeled AmflowFamily file, one path-valued constant line
    #    rewritten; PROVENANCE.json beside it pins the record's and the
    #    vendored sha256); DOGTAG_NPDBOX_FAMILY (or the legacy
    #    TOPOLOGY_AUDIT_NPDBOX_FAMILY) overrides.  The
    #    fingerprint this leg printed on the record file is pinned: 7 lines,
    #    2 ISPs, 0 massive, planar True in the canonical leg order, cut
    #    {s: 2, t: 3, u: 4}, canonical hash b6577ea5d3bf, isomorphic to
    #    planar_smirnov_dbox under [k2 -> -k2].
    why = _require_nx()
    if why:
        return {"status": "SKIP", "detail": why, "results": {}}
    path, why = _fixture_or_skip("DOGTAG_NPDBOX_FAMILY",
                                 os.path.join("npdbox", "family.jl"),
                                 "NP-dbox family file")
    if path is None:
        return {"status": "SKIP", "detail": why, "results": {}}
    print(f"   file: {path}")
    fp = audit(load_family(path)[0])
    _print_audit(fp)
    iso = [m for m in fp["catalog_matches"] if "planar_smirnov" in m["family"]]
    ok = (fp["planar"] is True) and bool(iso)
    pinned = {"n_genuine_propagators": 7, "n_isp_or_dotproduct": 2, "n_massive_lines": 0,
              "cut_signature": {"s": 2, "t": 3, "u": 4}, "canonical_hash": "b6577ea5d3bf"}
    got = {k: (dict(fp[k]) if k == "cut_signature" else fp[k]) for k in pinned}
    pin_ok = got == pinned and any(m["relabeling"] == "k2 -> -k2" for m in iso)
    print(f"   EXPECT: planar + iso to planar_smirnov  ->  {'PASS' if ok else 'FAIL'}")
    print(f"   EXPECT: fingerprint == the leg's record print {pinned} under [k2 -> -k2]"
          f"  ->  {'PASS' if pin_ok else 'FAIL (got %s, %s)' % (got, [m['relabeling'] for m in iso])}")
    return {"status": "PASS" if (ok and pin_ok) else "FAIL", "detail": "",
            "results": {"npdbox": fp}}


def _selftest_planar_control():
    # 2. Genuine planar Smirnov dbox catalog rep (must NOT false-flag nonplanar)
    why = _require_nx()
    if why:
        return {"status": "SKIP", "detail": why, "results": {}}
    cf = dict(CATALOG_SOURCES["planar_smirnov_dbox"])
    fam = Family("planar_dbox_control", cf["loops"], cf["exts"],
                 cf["ext_subs"], cf["propagators"], "memory:control")
    fp = audit(fam)
    _print_audit(fp)
    ok = fp["planar"] is True and not any("MISMATCH" in w for w in fp["warnings"])
    print(f"   EXPECT: planar, no mismatch warning  ->  {'PASS' if ok else 'FAIL'}")
    return {"status": "PASS" if ok else "FAIL", "detail": "",
            "results": {"planar_control": fp}}


def _selftest_crossed_control():
    # 3. Every catalog entry declared "planar": False (the crossed boxes) must
    #    read nonplanar in the canonical leg order; no false 'planar' flag.
    #    The entries are found by their declared planarity, not by name, so
    #    the catalog can change without touching this leg.
    why = _require_nx()
    if why:
        return {"status": "SKIP", "detail": why, "results": {}}
    names = [n for n, spec in CATALOG_SOURCES.items() if spec.get("planar") is False]
    if not names:
        print("   no catalog entry declares planar False")
        return {"status": "FAIL", "detail": "no crossed entry in the catalog",
                "results": {}}
    ok = True
    results = {}
    for n in names:
        cf = dict(CATALOG_SOURCES[n])
        fam = Family(n + "_control", cf["loops"], cf["exts"], cf["ext_subs"],
                     cf["propagators"], "memory:control")
        fp = audit(fam)
        results["crossed_control_" + n] = fp
        _print_audit(fp)
        leg_ok = fp["planar"] is False
        ok = ok and leg_ok
        print(f"   EXPECT: {n} nonplanar  ->  {'PASS' if leg_ok else 'FAIL (see notes)'}")
    return {"status": "PASS" if ok else "FAIL", "detail": "", "results": results}


def _selftest_c3_dbox():
    # 4. C3 planar dbox (one internal mass) -- a real family of record,
    #    planar, must NOT false-flag.  Vendored under tests/fixtures/row07/
    #    (the record's bytes with one path-valued constant line rewritten;
    #    PROVENANCE.json beside it pins both the record's and the vendored
    #    sha256 and names the line).
    why = _require_nx()
    if why:
        return {"status": "SKIP", "detail": why, "results": {}}
    path, why = _fixture_or_skip("DOGTAG_C3_FAMILY",
                                 os.path.join("row07", "family_of_record.jl"),
                                 "C3 dbox family.jl")
    if path is None:
        return {"status": "SKIP", "detail": why, "results": {}}
    print(f"   file: {path}")
    fp = audit(load_family(path)[0])
    _print_audit(fp)
    ok = fp["planar"] is True and not any("MISMATCH" in w for w in fp["warnings"])
    print(f"   EXPECT: planar, no mismatch  ->  {'PASS' if ok else 'FAIL'}")
    return {"status": "PASS" if ok else "FAIL", "detail": "", "results": {"c3_dbox": fp}}


def _selftest_banana():
    # 5. A 3-loop banana (sunrise-type) -- different topology, must not match
    #    the dbox catalog.  Vendored fixture: the row-19 three-loop banana
    #    family of record (tests/fixtures/row19/); the file this leg
    #    originally named is no longer present, DOGTAG_BANANA_FAMILY (or the
    #    legacy TOPOLOGY_AUDIT_BANANA_FAMILY) overrides.
    path, why = _fixture_or_skip("DOGTAG_BANANA_FAMILY",
                                 os.path.join("row19", "family_of_record.yaml"),
                                 "3-loop banana family")
    if path is None:
        return {"status": "SKIP", "detail": why, "results": {}}
    print(f"   file: {path}")
    fp = audit(load_family(path)[0])
    _print_audit(fp)
    ok = fp["n_loops"] == 3 and not fp["catalog_matches"]
    print(f"   EXPECT: 3 loops, no dbox catalog match  ->  {'PASS' if ok else 'FAIL'}")
    return {"status": "PASS" if ok else "FAIL", "detail": "", "results": {"banana": fp}}


def _selftest_catalog_integrity():
    # 6. Catalog integrity sweep: every catalog entry must parse to its declared
    #    propagator count, match ITSELF under the identity relabeling, read its
    #    declared "planar" verdict, and match NO other entry (in particular the
    #    planar/crossed dbox pair must not cross-match: their propagator sets
    #    differ under every signed loop relabeling with legs fixed).
    why = _require_nx()
    if why:
        return {"status": "SKIP", "detail": why, "results": {}}
    ok = True
    results = {}
    for cname, spec in CATALOG_SOURCES.items():
        fam = Family(cname, spec["loops"], spec["exts"], spec["ext_subs"],
                     spec["propagators"], "memory:catalog")
        fp = audit(fam)
        results["catalog_" + cname] = fp
        nprops_ok = fp["n_genuine_propagators"] == len(spec["propagators"])
        selfmatch = any(m["family"] == cname and m["relabeling"] == "identity"
                        for m in fp["catalog_matches"])
        expect_planar = spec.get("planar")
        planar_ok = (expect_planar is None) or (fp["planar"] is expect_planar)
        crossmatch = [m["family"] for m in fp["catalog_matches"]
                      if m["family"] != cname]
        leg_ok = nprops_ok and selfmatch and planar_ok and not crossmatch
        ok = ok and leg_ok
        print(f"   {cname}: props={fp['n_genuine_propagators']}/"
              f"{len(spec['propagators'])}  planar={fp['planar']} "
              f"(expect {expect_planar})  self-match={'yes' if selfmatch else 'NO'}  "
              f"cross-match={crossmatch if crossmatch else 'none'}  ->  "
              f"{'PASS' if leg_ok else 'FAIL'}", flush=True)
    print(f"   EXPECT: every entry self-matches, planarity as declared, no "
          f"cross-matches  ->  {'PASS' if ok else 'FAIL'}", flush=True)
    return {"status": "PASS" if ok else "FAIL", "detail": "", "results": results}


def _selftest_signed_relabel():
    # 7. Signed-relabeling positive control (in-memory): the one-loop box
    #    written with k1 -> -k1 must be caught as graph-isomorphic to
    #    'massless_box_1l' under the sign flip -- the q^2 == (-q)^2 relabeling
    #    catch this tool exists for, exercised without any fixture file.
    fam = Family("relabeled_box_control", ["k1"], ["p1", "p2", "p3", "p4"],
                 {"p4": "-p1 - p2 - p3"},
                 [("k1^2", 0), ("(k1 - p1)^2", 0), ("(k1 - p1 - p2)^2", 0),
                  ("(k1 - p1 - p2 - p3)^2", 0)], "memory:control")
    fp = audit(fam)
    _print_audit(fp)
    ok = any(m["family"] == "massless_box_1l" and "k1 -> -k1" in m["relabeling"]
             for m in fp["catalog_matches"])
    print(f"   EXPECT: ISO to massless_box_1l under [k1 -> -k1]  ->  "
          f"{'PASS' if ok else 'FAIL'}", flush=True)
    return {"status": "PASS" if ok else "FAIL", "detail": "",
            "results": {"relabeled_box_control": fp}}


def _selftest_readers():
    # 8. The reader battery (tests/test_readers.py): every vendored record
    #    fixture under tests/fixtures/ must load through load_family and
    #    reproduce its hand-checked object (fingerprint hash, counts, leg set,
    #    (momentum, mass) list against the transcription of record), and each
    #    planted negative control must flip.  Fixtures are sha-pinned and
    #    refused on drift.  SKIP by name when the battery or its fixtures are
    #    absent.
    here = os.path.dirname(os.path.abspath(__file__))
    mod_path = os.path.join(here, "tests", "test_readers.py")
    if not os.path.exists(mod_path):
        return {"status": "SKIP", "detail": "tests/test_readers.py absent", "results": {}}
    if not os.path.isdir(_FIXTURE_DIR):
        return {"status": "SKIP", "detail": "tests/fixtures/ absent", "results": {}}
    import importlib.util
    spec = importlib.util.spec_from_file_location("topology_audit_test_readers", mod_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    rows = mod.run_battery(verbose=True)
    n_pass = sum(1 for r in rows if r["status"] == "PASS")
    n_fail = [r["name"] for r in rows if r["status"] == "FAIL"]
    n_skip = [r["name"] + ": " + r["detail"] for r in rows if r["status"] == "SKIP"]
    print(f"   battery: {len(rows)} cases, {n_pass} PASS, {len(n_fail)} FAIL, "
          f"{len(n_skip)} SKIP", flush=True)
    for name in n_fail:
        print(f"   FAIL {name}")
    if n_skip:
        return {"status": "SKIP", "detail": "; ".join(n_skip),
                "results": {"reader_battery": rows}}
    return {"status": "PASS" if not n_fail else "FAIL",
            "detail": ", ".join(n_fail), "results": {"reader_battery": rows}}


register_self_test_leg("1", "NP-dbox (the mislabeled family, must flag planar + iso)",
                       _selftest_npdbox)
register_self_test_leg("2", "catalog planar Smirnov dbox (control, must read planar)",
                       _selftest_planar_control)
register_self_test_leg("3", "catalog crossed box entries (declared planar False, "
                            "must read NONplanar)", _selftest_crossed_control)
register_self_test_leg("4", "C3 dbox (rung mass), family of record, planar, no mismatch",
                       _selftest_c3_dbox)
register_self_test_leg("5", "3-loop banana (control, different topology, no dbox match)",
                       _selftest_banana)
register_self_test_leg("6", "catalog integrity (parse, self-match, declared planarity, "
                            "no cross-match)", _selftest_catalog_integrity)
register_self_test_leg("7", "signed-relabeling positive control (one-loop box, k1 -> -k1)",
                       _selftest_signed_relabel)
register_self_test_leg("8", "reader battery (tests/test_readers.py on the vendored "
                            "record fixtures + planted controls)", _selftest_readers)


def _run_extra_leg(fn):
    """Adapt an EXTRA_SELF_TEST_LEGS callable to the dispatcher protocol.

    The callable returns its own counts: {"ok": bool, "executed": int,
    "skipped": int | [names], "expected_fail": int (optional),
    "skipped_names": [...] (optional), "results": {...} (optional)}.  The
    leg is FAIL when ok is False, SKIP (named) when any sub-leg was skipped,
    PASS otherwise; expected-fail sub-legs are counted and named, never
    folded into PASS silently.  Returns {"status", "detail", "results",
    "counts": {"executed", "skipped", "expected_fail", "skipped_names"}}."""
    out = fn() or {}
    skipped = out.get("skipped", 0)
    names = list(out.get("skipped_names", []) or [])
    if isinstance(skipped, (list, tuple)):
        names = list(skipped) + [n for n in names if n not in skipped]
        n_skip = len(skipped)
    else:
        n_skip = int(skipped or 0)
    n_exec = int(out.get("executed", 0) or 0)
    n_xfail = int(out.get("expected_fail", 0) or 0)
    ok = bool(out.get("ok", False))
    counts = {"executed": n_exec, "skipped": n_skip, "expected_fail": n_xfail,
              "skipped_names": names}
    detail = f"{n_exec} sub-legs executed, {n_skip} skipped"
    if names:
        detail += f": {names}"
    if n_xfail:
        detail += f"; {n_xfail} EXPECTED-FAIL by name (see the leg's lines)"
    status = "FAIL" if not ok else ("SKIP" if n_skip else "PASS")
    return {"status": status, "detail": detail,
            "results": out.get("results", {}) or {}, "counts": counts}


def self_test():
    """Run every registered leg -- SELF_TEST_LEGS (the dispatcher protocol)
    followed by EXTRA_SELF_TEST_LEGS (count-returning callables, adapted by
    _run_extra_leg) -- print each leg PASS / SKIP (named) / FAIL, count
    executed vs skipped, and print ALL PASS only when no leg failed, none
    was skipped and no expected-fail sub-leg is open.  Returns (results, rc)
    with rc = 1 when any leg FAILED, 2 when any leg was SKIPPED (a skip is an
    unknown), else 0 (declared expected-fail sub-legs keep rc 0 but are named
    on the summary line and hold back ALL PASS)."""
    print("=" * 74)
    print("TOPOLOGY AUDITOR SELF-TEST")
    print("=" * 74)
    results = {}
    statuses = []
    sub = {"executed": 0, "skipped": 0, "expected_fail": 0, "skipped_names": []}
    legs = list(SELF_TEST_LEGS)
    legs += [(leg_id, title, (lambda fn=fn: _run_extra_leg(fn)))
             for leg_id, title, fn in globals().get("EXTRA_SELF_TEST_LEGS", [])]
    for leg_id, title, fn in legs:
        print(f"\n[{leg_id}] {title}:", flush=True)
        try:
            out = fn() or {}
        except Exception as exc:      # a crashing leg is a FAIL, never a skip
            out = {"status": "FAIL", "detail": f"{type(exc).__name__}: {exc}",
                   "results": {}}
        status = str(out.get("status", "FAIL")).upper()
        if status not in ("PASS", "FAIL", "SKIP"):
            status = "FAIL"
        detail = out.get("detail", "")
        results.update(out.get("results", {}) or {})
        for k in ("executed", "skipped", "expected_fail"):
            sub[k] += int((out.get("counts") or {}).get(k, 0) or 0)
        sub["skipped_names"] += list((out.get("counts") or {}).get("skipped_names", []))
        statuses.append((leg_id, title, status, detail))
        tail = f"  ({detail})" if detail else ""
        print(f"   [{leg_id}] {status}{tail}", flush=True)
    executed = [s for s in statuses if s[2] in ("PASS", "FAIL")]
    skipped = [s for s in statuses if s[2] == "SKIP"]
    failed = [s for s in statuses if s[2] == "FAIL"]
    print("\n" + "=" * 74)
    print(f"SELF-TEST: {len(executed)} executed / {len(skipped)} skipped "
          f"({len(executed) - len(failed)} PASS, {len(failed)} FAIL) of "
          f"{len(statuses)} legs")
    if sub["executed"] or sub["skipped"] or sub["expected_fail"]:
        print(f"SELF-TEST: registered sub-legs {sub['executed']} executed / "
              f"{sub['skipped']} skipped / {sub['expected_fail']} expected-fail"
              + (f" (skipped: {sub['skipped_names']})" if sub["skipped_names"] else ""))
    for leg_id, title, _s, detail in skipped:
        print(f"SKIPPED [{leg_id}] {title} -- {detail}")
    for leg_id, title, _s, detail in failed:
        print(f"FAILED  [{leg_id}] {title}" + (f" -- {detail}" if detail else ""))
    rc = 1 if failed else (2 if skipped else 0)
    if failed:
        print("SELF-TEST: SOME FAILED -- see above (rc 1)")
    elif skipped or sub["expected_fail"]:
        print("SELF-TEST: PASS on the executed legs; NOT ALL PASS "
              f"({len(skipped)} skipped{', named above' if skipped else ''}; "
              f"{sub['expected_fail']} open EXPECTED-FAIL sub-legs by name) (rc {rc})")
    else:
        print("SELF-TEST: ALL PASS (rc 0)")
    print("=" * 74)
    return results, rc


def _selftest_realizer():
    """
    Self-test leg of the vertex realizer (build_graph / connected_covers /
    non_graph_test / cut_signature).  Registered in EXTRA_SELF_TEST_LEGS for the
    self-test dispatcher; also runnable on its own.

    Sub-legs: in-memory controls (no fixture needed) and vendored fixtures
    under tests/fixtures/ (each pinned by sha256 in its PROVENANCE.json and
    refused on drift; an absent fixture is SKIPPED by name and counted).
    Returns {"ok": bool, "executed": int, "skipped": [names]}.
    """
    import hashlib
    why = _require_nx("[realizer]")
    if why:
        print(f"   [realizer] SKIP -- {why}", flush=True)
        return {"ok": True, "executed": 0, "skipped": [why]}
    here = os.path.dirname(os.path.abspath(__file__))
    fixdir = os.path.join(here, "tests", "fixtures")
    results = []
    skipped = []

    def leg(name, ok, detail=""):
        results.append(bool(ok))
        print(f"   [realizer] {name}: {'PASS' if ok else 'FAIL'}  {detail}", flush=True)

    def pinned(row, rel):
        """Path of a vendored fixture, or None (SKIP) when absent; a sha256
        mismatch against PROVENANCE.json is a FAIL by name (refused)."""
        path = os.path.join(fixdir, row, rel)
        prov = os.path.join(fixdir, row, "PROVENANCE.json")
        if not (os.path.exists(path) and os.path.exists(prov)):
            skipped.append("%s/%s" % (row, rel))
            print(f"   [realizer] fixture {row}/{rel} SKIPPED (absent)", flush=True)
            return None
        want = None
        for ent in json.load(open(prov)).get("files", []):
            if ent.get("path") == rel:
                want = ent.get("sha256")
        have = hashlib.sha256(open(path, "rb").read()).hexdigest()
        if want != have:
            leg("fixture pin %s/%s" % (row, rel), False,
                "sha256 drift: PROVENANCE %s != file %s -- refused" % (want, have))
            return None
        return path

    def realized(fam):
        G, ok = build_graph(fam)
        return ok, G.graph.get("realization", {})

    print("\n[realizer] vertex realizer: connected covers, V = E - L + 1, non-graph verdict:",
          flush=True)

    # realizer.1  textbook planar double box (in memory): six vertices, unique, cuts.
    cf = CATALOG_SOURCES["planar_smirnov_dbox"]
    fam = Family("g6_planar_dbox", cf["loops"], cf["exts"], cf["ext_subs"],
                 cf["propagators"], "memory:control")
    ok, rep = realized(fam)
    csig = cut_signature(fam)
    leg("planar double box realizes: V = 6, connected, unique, 3V = 18 <= 2E + N = 18, "
        "cuts {s: 2, t: 3, u: 4}",
        ok and rep.get("V") == 6 and rep.get("connected") is True
        and rep.get("realization_unique") is True and rep.get("three_V") == 18
        and rep.get("two_E_plus_N") == 18 and rep.get("conservation_at_every_vertex") is True
        and dict(csig) == {"s": 2, "t": 3, "u": 4} and csig.cause is None,
        "V=%s unique=%s cuts=%s" % (rep.get("V"), rep.get("realization_unique"), dict(csig)))

    # realizer.2  seven massless lines of a two-loop three-point family all in the
    #       denominator (the ladder with its ISP promoted to a line): 3V = 18 >
    #       2E + N = 17, no realization -> NON-GRAPH by name, no planar verdict.
    fam7 = Family("g6_sevenline_3pt", ["k1", "k2"], ["p1", "p2", "p3"],
                  {"p3": "-p1 - p2"},
                  [("k1^2", 0), ("k2^2", 0), ("(k1 + k2)^2", 0), ("(k1 + p1)^2", 0),
                   ("(k2 + p1)^2", 0), ("(k1 + p1 + p2)^2", 0), ("(k2 + p1 + p2)^2", 0)],
                  "memory:control")
    ok7, rep7 = realized(fam7)
    csig7 = cut_signature(fam7)
    plan7 = planarity(fam7)
    leg("seven-line three-point object is NON-GRAPH (3V = 18 > 2E + N = 17), cut "
        "entries None with the cause named, planar None",
        (not ok7) and rep7.get("verdict") == "NON-GRAPH"
        and rep7.get("three_V") == 18 and rep7.get("two_E_plus_N") == 17
        and "3V = 18 > 2E + N = 17" in (rep7.get("cause") or "")
        and all(v is None for v in csig7.values())
        and "3V = 18 > 2E + N = 17" in (csig7.cause or "")
        and plan7.get("planar") is None,
        "verdict=%s cause=%s" % (rep7.get("verdict"), rep7.get("cause")))

    # realizer.3  planted control: the planar double box with one ring line listed
    #       twice in the top sector -> refused by name (repeated line / count).
    fam_dup = Family("g6_dup_line", cf["loops"], cf["exts"], cf["ext_subs"],
                     list(cf["propagators"]) + [cf["propagators"][1]], "memory:control")
    okd, repd = realized(fam_dup)
    leg("duplicated ring line refused by name",
        (not okd) and repd.get("verdict") == "NON-GRAPH"
        and ("repeated line" in (repd.get("cause") or "")
             or "3V =" in (repd.get("cause") or "")),
        "cause=%s" % repd.get("cause"))

    # realizer.4  planted control: the step budget is reported by name when exhausted.
    Gb, okb = _realize_graph(fam, step_budget=1)
    repb = Gb.graph.get("realization", {})
    leg("step budget exhaustion reported by name (budget 1)",
        (not okb) and repb.get("budget_exhausted") is True
        and "step budget" in (repb.get("cause") or "") and repb.get("verdict") == "UNDECIDED",
        "cause=%s" % repb.get("cause"))

    # realizer.5  vendored fixture: box with a kite self-energy on one rung (row 28
    #       of the identity census; the record declares four legs, p4 = -p1-p2-p3).
    p = pinned("row28", "family_of_record.yaml")
    if p:
        f28 = load_family(p)[0]
        f28.exts = ["p1", "p2", "p3", "p4"]
        f28.ext_subs = {"p4": "-p1 - p2 - p3"}
        ok28, rep28 = realized(f28)
        c28 = cut_signature(f28)
        leg("row28 kite-on-a-rung box realizes CONNECTED with V = 6, cuts {s: 3, t: 2, u: 5}",
            ok28 and rep28.get("connected") is True and rep28.get("V") == 6
            and rep28.get("E") == 8 and rep28.get("three_V") == 18
            and rep28.get("two_E_plus_N") == 20 and dict(c28) == {"s": 3, "t": 2, "u": 5},
            "V=%s cuts=%s" % (rep28.get("V"), dict(c28)))

    # realizer.6  vendored fixtures: row 6 retired seven-line object (NON-GRAPH) and
    #       the cured ladder (a graph with V = 5).
    p = pinned("row06", "retired_sevenline_vertex2L/family_of_record_transcribed.yaml")
    if p:
        f6r = load_family(p)[0]
        f6r.exts = ["p1", "p2", "p3"]
        f6r.ext_subs = {"p3": "-p1 - p2"}
        ok6r, rep6r = realized(f6r)
        leg("row06 retired seven-line vertex is NON-GRAPH (3V = 18 > 2E + N = 17)",
            (not ok6r) and rep6r.get("verdict") == "NON-GRAPH"
            and rep6r.get("three_V") == 18 and rep6r.get("two_E_plus_N") == 17,
            "cause=%s" % rep6r.get("cause"))
    p = pinned("row06", "cured_ladder_sudakovPR0/family_of_record.yaml")
    if p:
        f6c = load_family(p)[0]
        f6c.exts = ["p1", "p2", "p3"]
        f6c.ext_subs = {"p3": "-p1 - p2"}
        ok6c, rep6c = realized(f6c)
        c6c = cut_signature(f6c)
        leg("row06 cured ladder is a graph: V = 5, unique, cuts {p1p2: 2, p1p3: 2, p2p3: 2}",
            ok6c and rep6c.get("V") == 5 and rep6c.get("realization_unique") is True
            and dict(c6c) == {"p1p2": 2, "p1p3": 2, "p2p3": 2},
            "V=%s cuts=%s" % (rep6c.get("V"), dict(c6c)))

    # realizer.7  vendored fixture: a drawn double box routed by a spanning tree (row
    #       1 drawing): six vertices, not five, unique realization.
    p = pinned("row01", "drawn_graph.yaml")
    if p:
        fd = load_family(p)[0]
        okd1, repd1 = realized(fd)
        cd1 = cut_signature(fd)
        leg("row01 drawn double box realizes with V = 6 (not 5), unique, cuts {s: 2, t: 3, u: 4}",
            okd1 and repd1.get("V") == 6 and repd1.get("realization_unique") is True
            and dict(cd1) == {"s": 2, "t": 3, "u": 4},
            "V=%s unique=%s cuts=%s" % (repd1.get("V"), repd1.get("realization_unique"), dict(cd1)))

    allok = all(results)
    print(f"   [realizer] {len(results)} executed / {len(skipped)} skipped  ->  "
          f"{'PASS' if allok and not skipped else ('FAIL' if not allok else 'PASS WITH SKIPS')}",
          flush=True)
    return {"ok": allok, "executed": len(results), "skipped": skipped}


# ---------------------------------------------------------------------------
#  6b. CATALOG SELF-TEST LEG  (registered in EXTRA_SELF_TEST_LEGS)
#      [3]  the crossed-box entry is the (2,1,1) theta graph: cut {s:2,t:3,u:3},
#           non-planar with legs on the outer face in EVERY leg order, no
#           planar catalog entry reachable by any signed loop relabeling with
#           any leg permutation (conservation re-imposed), self-match identity;
#           the vendored fixture reproduces the same object by hash.
#      [3b] the RETIRED former entry is the named negative: it FAILS the [3]
#           criteria and IS planar_smirnov_dbox under p3 <-> p4 (proved by the
#           textual leg substitution and by the tool's own --leg-perms search,
#           which re-imposes the target's momentum conservation after the
#           permutation).
#      [3c] the fixture family in the t = (p1+p3)^2 convention no longer
#           matches a crossed entry.
#      [8]  the three-loop three-point off-shell ladder entry; every entry
#           declares "legs" and "planar_in_some_leg_order", the latter checked
#           by enumerating every leg order on the realized graph.
#      Fixtures live under tests/fixtures/<row>/ beside this file and are
#      pinned by sha256: an absent fixture is a named SKIP (counted), a
#      drifted fixture is a FAIL.
# ---------------------------------------------------------------------------
def _selftest_catalog(verbose=True):
    """Return {"ok", "executed", "skipped", "expected_fail", "skipped_names",
    "results"}; prints every sub-leg by name when verbose."""
    import hashlib
    say = print if verbose else (lambda *a, **k: None)
    why = _require_nx("[catalog]")
    if why:
        say(f"   [catalog] SKIP -- {why}")
        return {"ok": True, "executed": 0, "skipped": 1, "expected_fail": 0,
                "skipped_names": [why], "results": {}}
    here = os.path.dirname(os.path.abspath(__file__))
    fxroot = os.path.join(here, "tests", "fixtures")
    # fixture relpath -> sha256 of the vendored copy (PROVENANCE.json "sha256";
    # the record's own sha256 is pinned there as record_sha256)
    FIX = {
        "row02/genuine_crossed_box_theta211.yaml":
            "05d1300140c4b9cad4929783b25df5225b1e35623491a0b96903097d5c18690e",
        "row02/family_of_record.jl":
            "ddc75ae24954b08103d01e996f614cae6a6b74a8bbc485c80711cc7e17434917",
        "row04/family_of_record.yaml":
            "c61b80c50011b900a1144d396d50d9c1b41e6c77dc40566da6201abefa0d633c",
        "row04/drawn_graph.yaml":
            "87bcbde125b06c812b48408a7b0cc2b5a54ca4700aa1730b8341f337d096915b",
    }
    THETA_CUT = {"s": 2, "t": 3, "u": 3}
    THETA_HASH = "7790fb92864c"
    RETIRED_CUT = {"s": 2, "t": 4, "u": 3}
    RETIRED_HASH = "a542a386080e"
    RETIRED_LEG_MAP = {"p3": "p4", "p4": "p3"}
    UD_HASH = "fd256fcd39ad"
    UD_CUT = {"p1p2": 2, "p1p3": 2, "p2p3": 2}
    UD_LEGS = {"p1": "offshell", "p2": "0", "p3": "offshell"}
    UD_VE = (10, 12)
    LEG_CLASSES = ("0", "offshell")

    out = {"ok": True, "executed": 0, "skipped": 0, "expected_fail": 0,
           "skipped_names": [], "results": {}}

    def leg(name, ok, detail):
        out["executed"] += 1
        out["ok"] = out["ok"] and bool(ok)
        out["results"][name] = {"ok": bool(ok), "detail": detail}
        say(f"   {name}: {detail}  ->  {'PASS' if ok else 'FAIL'}")

    def xfail(name, ok, detail, reason):
        # an expectation the tool cannot meet yet: PASS if it does, else
        # EXPECTED-FAIL by name (counted separately, does not clear ALL PASS
        # on its own but is reported)
        out["executed"] += 1
        if ok:
            out["results"][name] = {"ok": True, "detail": detail}
            say(f"   {name}: {detail}  ->  PASS (expected-fail cleared)")
        else:
            out["expected_fail"] += 1
            out["results"][name] = {"ok": None, "detail": detail,
                                    "expected_fail": reason}
            say(f"   {name}: {detail}  ->  EXPECTED-FAIL ({reason})")

    def fixture(rel):
        path = os.path.join(fxroot, rel)
        if not os.path.exists(path):
            out["skipped"] += 1
            out["skipped_names"].append(rel)
            say(f"   fixture tests/fixtures/{rel} SKIPPED (absent)")
            return None
        got = hashlib.sha256(open(path, "rb").read()).hexdigest()
        if got != FIX[rel]:
            leg(f"fixture pin {rel}", False,
                f"sha256 {got[:16]} != pinned {FIX[rel][:16]} (drift: refused)")
            return None
        return path

    def fam_from(name, spec):
        return Family(name, spec["loops"], spec["exts"], spec["ext_subs"],
                      spec["propagators"], "memory:" + name)

    def leg_orders(fam):
        """(abstract_planar, n_planar_orders, n_orders) on the realized graph,
        closing the legs in every cyclic order."""
        G, ok = build_graph(fam)
        if not ok or G is None:
            return None, 0, 0
        base = nx.Graph()
        base.add_edges_from((u, v) for u, v, k in G.edges(keys=True))
        ext = [n for n in base.nodes() if str(n).startswith("ext_")]
        n_ok, n_all = 0, 0
        for perm in itertools.permutations(ext):
            H = base.copy()
            for a in range(len(perm)):
                H.add_edge(perm[a], perm[(a + 1) % len(perm)])
            n_all += 1
            n_ok += 1 if nx.check_planarity(H)[0] else 0
        return bool(nx.check_planarity(base)[0]), n_ok, n_all

    def bare_graph(fam):
        G, ok = build_graph(fam)
        if not ok or G is None:
            return None
        base = nx.Graph()
        base.add_edges_from((u, v) for u, v, k in G.edges(keys=True))
        return base

    def leg_subst(props, legmap, loops, exts):
        """Rewrite every propagator string under the leg map p_a -> p_b
        (simultaneous); momentum conservation is re-imposed when the family
        built from the result resolves its ext_subs."""
        syms = {s: sp.Symbol(s) for s in list(loops) + list(exts)}
        res = []
        for e, m in props:
            lin = (_extract_linear_from_square(e, syms)
                   if re.search(r'\^\s*2|\*\*\s*2', e)
                   else sp.sympify(e, locals=syms))
            lin2 = lin.subs({syms[a]: syms[b] for a, b in legmap.items()},
                            simultaneous=True)
            res.append((f"({sp.sstr(lin2)})^2", m))
        return res

    def complete_search(specA, specB):
        """Every leg permutation applied to A's strings (conservation
        re-imposed), then the signed loop relabeling search with legs fixed.
        Returns [(legmap, relabeling)] -- the complete search the catalog
        claims need."""
        fb = fam_from("b", specB)
        vb, _, _ = momentum_vectors(fb)
        found = []
        exts = list(specA["exts"])
        for perm in itertools.permutations(exts):
            lm = dict(zip(exts, perm))
            fa = Family("a", specA["loops"], exts, specA["ext_subs"],
                        leg_subst(specA["propagators"], lm, specA["loops"], exts),
                        "memory:a")
            va, _, _ = momentum_vectors(fa)
            rel = isomorphism_to(va, vb, len(specA["loops"]), len(exts))
            if rel is not None:
                found.append((lm, rel))
        return found

    def crossed_box_criteria(spec):
        """The [3] criteria on one 2-loop 4-leg entry; returns (all_ok, dict)."""
        fam = fam_from("crossed_box_candidate", spec)
        fp = audit(fam)
        ab, n_ok, n_all = leg_orders(fam)
        planar_entries = [n for n, s in CATALOG_SOURCES.items()
                          if s.get("planar") is True
                          and len(s["loops"]) == 2 and len(s["exts"]) == 4]
        reach = {n: complete_search(spec, CATALOG_SOURCES[n]) for n in planar_entries}
        d = {"cut_signature": fp["cut_signature"], "planar": fp["planar"],
             "canonical_hash": fp["canonical_hash"], "abstract_planar": ab,
             "n_planar_leg_orders": n_ok, "n_leg_orders": n_all,
             "planar_entries_reachable": {n: r for n, r in reach.items() if r}}
        ok = (fp["cut_signature"] == THETA_CUT and fp["planar"] is False
              and ab is True and n_all > 0 and n_ok == 0
              and not d["planar_entries_reachable"])
        return ok, d

    say("\n[3] catalog crossed box = the (2,1,1) theta graph (object of record):")
    theta = CATALOG_SOURCES["crossed_tausk_dbox"]
    ok3, d3 = crossed_box_criteria(theta)
    leg("[3] theta entry criteria", ok3,
        f"cut={d3['cut_signature']} (expect {THETA_CUT}) planar={d3['planar']} "
        f"abstract_planar={d3['abstract_planar']} planar_leg_orders="
        f"{d3['n_planar_leg_orders']}/{d3['n_leg_orders']} (expect 0) "
        f"planar entries reachable under any leg map={d3['planar_entries_reachable']}")
    fpt = audit(fam_from("theta_check", theta))
    leg("[3] theta hash + self-match", fpt["canonical_hash"] == THETA_HASH and
        any(m["family"] == "crossed_tausk_dbox" and m["relabeling"] == "identity"
            for m in fpt["catalog_matches"]),
        f"canonical_hash={fpt['canonical_hash']} (expect {THETA_HASH}); "
        f"matches={[(m['family'], m['relabeling']) for m in fpt['catalog_matches']]}")
    pth = fixture("row02/genuine_crossed_box_theta211.yaml")
    if pth:
        fpy = audit(load_family(pth)[0], try_leg_perms=True)
        leg("[3] fixture theta yaml", fpy["canonical_hash"] == THETA_HASH and
            fpy["cut_signature"] == THETA_CUT and fpy["planar"] is False and
            [(m["family"], m["relabeling"]) for m in fpy["catalog_matches"]]
            == [("crossed_tausk_dbox", "identity")],
            f"hash={fpy['canonical_hash']} cut={fpy['cut_signature']} "
            f"planar={fpy['planar']} matches="
            f"{[(m['family'], m['relabeling']) for m in fpy['catalog_matches']]}")
    # the planarity() field for the theta entry is recorded, not asserted:
    # the closure test's "planar_in_some_leg_order" shortcut reads abstract
    # planarity and disagrees with the enumeration on this graph
    pl = planarity(fam_from("theta_planarity", theta))
    out["results"]["[3] planarity() field vs enumeration"] = {
        "planarity_planar_in_some_leg_order": pl.get("planar_in_some_leg_order"),
        "enumeration_planar_in_some_leg_order": d3["n_planar_leg_orders"] > 0}
    say(f"   (planarity()['planar_in_some_leg_order']={pl.get('planar_in_some_leg_order')} "
        f"vs enumeration {d3['n_planar_leg_orders'] > 0}: recorded)")

    say("\n[3b] retired former entry (named negative: NOT the crossed box):")
    ret = CATALOG_RETIRED["planar_smirnov_dbox_legs34_interchanged"]
    okr, dr = crossed_box_criteria(ret)
    leg("[3b] retired entry FAILS the crossed-box criteria", (not okr) and
        dr["cut_signature"] == RETIRED_CUT and dr["canonical_hash"] == RETIRED_HASH
        and dr["n_planar_leg_orders"] > 0,
        f"cut={dr['cut_signature']} (expect {RETIRED_CUT}) hash={dr['canonical_hash']} "
        f"(expect {RETIRED_HASH}) planar_leg_orders={dr['n_planar_leg_orders']}/"
        f"{dr['n_leg_orders']} reachable={list(dr['planar_entries_reachable'])}")
    smi = CATALOG_SOURCES[ret["is"]]
    fsw = Family("retired_p3p4", ret["loops"], ret["exts"], ret["ext_subs"],
                 leg_subst(ret["propagators"], ret["is_under_leg_map"],
                           ret["loops"], ret["exts"]), "memory:retired_p3p4")
    rel_sw = isomorphism_to(momentum_vectors(fsw)[0], momentum_vectors(fam_from("smi", smi))[0],
                            2, 4)
    found = complete_search(ret, smi)
    full_map = {e: RETIRED_LEG_MAP.get(e, e) for e in ret["exts"]}
    leg("[3b] retired == planar_smirnov_dbox under p3<->p4", rel_sw is not None and
        any(lm == full_map for lm, _ in found),
        f"textual p3<->p4 then signed-relabeling search: {rel_sw}; complete leg-map "
        f"search finds {found}")
    gt, gr, gs = bare_graph(fam_from("t", theta)), bare_graph(ret and fam_from("r", ret)), bare_graph(fam_from("s", smi))
    leg("[3b] bare-graph isomorphism (VF2)", gt is not None and gr is not None and gs is not None and
        nx.is_isomorphic(gr, gs) and not nx.is_isomorphic(gt, gr) and not nx.is_isomorphic(gt, gs),
        f"retired~planar_box={gr is not None and gs is not None and nx.is_isomorphic(gr, gs)} "
        f"theta~retired={gt is not None and gr is not None and nx.is_isomorphic(gt, gr)} "
        f"theta~planar_box={gt is not None and gs is not None and nx.is_isomorphic(gt, gs)}")
    smi_fam = fam_from("s", smi)
    rel_tool = isomorphism_to(momentum_vectors(fam_from("r", ret))[0],
                              momentum_vectors(smi_fam)[0], 2, 4,
                              try_leg_perms=True, dep_rows=_dependent_leg_rows(smi_fam))
    leg("[3b] tool leg-permutation search finds p3<->p4",
        rel_tool is not None and "p3 -> p4" in rel_tool and "p4 -> p3" in rel_tool,
        f"isomorphism_to(..., try_leg_perms=True, dep_rows=<planar box>) = {rel_tool!r}")
    leg("[3b] retired entry is not searched",
        ret["former_name"] not in CATALOG_SOURCES or
        CATALOG_SOURCES[ret["former_name"]]["propagators"] != ret["propagators"],
        f"former name {ret['former_name']!r} now carries a different object; "
        f"CATALOG_RETIRED keys are not CATALOG_SOURCES keys: "
        f"{not set(CATALOG_RETIRED) & set(CATALOG_SOURCES)}")

    say("\n[3c] fixture family in the t=(p1+p3)^2 convention (must not read crossed):")
    p2 = fixture("row02/family_of_record.jl")
    if p2:
        fp2 = audit(load_family(p2)[0], try_leg_perms=True)
        mts = [(m["family"], m["relabeling"]) for m in fp2["catalog_matches"]]
        leg("[3c] no crossed match; any match is the planar box",
            fp2["canonical_hash"] == RETIRED_HASH and fp2["cut_signature"] == RETIRED_CUT
            and all("crossed" not in f for f, _ in mts)
            and all(f == "planar_smirnov_dbox" for f, _ in mts),
            f"hash={fp2['canonical_hash']} cut={fp2['cut_signature']} planar={fp2['planar']} "
            f"matches={mts}")

    say("\n[8] off-shell three-point ladder entry + leg virtuality classes on every entry:")
    ud = CATALOG_SOURCES["ud_ladder_3pt_offshell_3l"]
    fpu = audit(fam_from("ud_check", ud))
    plu = planarity(fam_from("ud_planarity", ud))
    leg("[8] ud_ladder_3pt_offshell_3l fingerprint",
        fpu["canonical_hash"] == UD_HASH and fpu["cut_signature"] == UD_CUT and
        fpu["planar"] is True and ud["legs"] == UD_LEGS and
        (plu.get("n_vertices"), plu.get("n_edges")) == UD_VE and
        [(m["family"], m["relabeling"]) for m in fpu["catalog_matches"]]
        == [("ud_ladder_3pt_offshell_3l", "identity")],
        f"hash={fpu['canonical_hash']} (expect {UD_HASH}) cut={fpu['cut_signature']} "
        f"planar={fpu['planar']} V,E={(plu.get('n_vertices'), plu.get('n_edges'))} "
        f"(expect {UD_VE}) legs={ud['legs']} matches="
        f"{[(m['family'], m['relabeling']) for m in fpu['catalog_matches']]}")
    bad = {}
    for cname, spec in list(CATALOG_SOURCES.items()) + list(CATALOG_RETIRED.items()):
        legs = spec.get("legs")
        why = []
        if not isinstance(legs, dict) or list(legs) != list(spec["exts"]):
            why.append("legs keys != exts")
        elif any(v not in LEG_CLASSES for v in legs.values()):
            why.append("leg class outside " + str(LEG_CLASSES))
        if spec.get("planar") not in (True, False):
            why.append("planar undeclared")
        if spec.get("planar_in_some_leg_order") not in (True, False):
            why.append("planar_in_some_leg_order undeclared")
        else:
            ab, n_ok, n_all = leg_orders(fam_from(cname, spec))
            if n_all == 0 or (n_ok > 0) != spec["planar_in_some_leg_order"]:
                why.append(f"planar_in_some_leg_order declared "
                           f"{spec['planar_in_some_leg_order']} but enumeration "
                           f"gives {n_ok}/{n_all} planar orders")
        if why:
            bad[cname] = why
    leg("[8] every entry declares legs / planar / planar_in_some_leg_order (enumerated)",
        not bad, f"{len(CATALOG_SOURCES)} entries + {len(CATALOG_RETIRED)} retired checked; "
        f"defects={bad if bad else 'none'}")
    p4 = fixture("row04/family_of_record.yaml")
    if p4:
        fp4 = audit(load_family(p4)[0], try_leg_perms=True)
        mts = [(m["family"], m["relabeling"]) for m in fp4["catalog_matches"]]
        src = open(p4).read()
        m = re.search(r'leg_virtualities:\s*\{([^}]*)\}', src)
        lv = dict(re.findall(r'(p\d+):\s*"([^"]*)"', m.group(1))) if m else None
        leg("[8] fixture ud ladder matches the ud entry, not the triple box",
            mts == [("ud_ladder_3pt_offshell_3l", "identity")] and
            fp4["canonical_hash"] == UD_HASH and fp4["cut_signature"] == UD_CUT and
            lv == ud["legs"],
            f"matches={mts} hash={fp4['canonical_hash']} cut={fp4['cut_signature']} "
            f"fixture leg_virtualities={lv} == entry legs {ud['legs']}: {lv == ud['legs']}")
    p4d = fixture("row04/drawn_graph.yaml")
    if p4d:
        f4d = load_family(p4d)[0]
        fp4d = audit(f4d, try_leg_perms=True)
        mts = [(m["family"], m["relabeling"]) for m in fp4d["catalog_matches"]]
        src = open(p4d).read()
        m = re.search(r'leg_virtualities:\s*\{([^}]*)\}', src)
        lv = dict(re.findall(r'(p\d+):\s*"([^"]*)"', m.group(1))) if m else None
        leg("[8] drawn four-point triple box: not the ud entry; legs all on shell",
            all(f != "ud_ladder_3pt_offshell_3l" for f, _ in mts) and
            lv == CATALOG_SOURCES["planar_triple_box_3l"]["legs"] and lv != ud["legs"] and
            isomorphism_to(momentum_vectors(f4d)[0],
                           momentum_vectors(fam_from("ud", ud))[0], 3, 4) is None,
            f"matches={mts} drawn leg_virtualities={lv}")
        # the two routings differ by a loop-momentum basis change: the affine
        # search (isomorphism unit) describes it and the realized-graph
        # isomorphism certifies it -- a plain sub-leg since that unit landed
        leg("[8] drawn four-point triple box matches planar_triple_box_3l",
            any(f == "planar_triple_box_3l" for f, _ in mts),
            f"matches={mts}")

    say(f"\n   catalog leg: {out['executed']} executed / {out['skipped']} skipped / "
        f"{out['expected_fail']} expected-fail"
        + (f" (skipped: {out['skipped_names']})" if out["skipped_names"] else ""))
    return out


# ---------------------------------------------------------------------------
#  leg-permutation self-test leg: --leg-perms with momentum conservation re-imposed, the
#  leg-count guard, and the label-preserving restriction.  Every value is the
#  identity census's hand-checked object for that row (quoted in the sub-leg
#  from rows/<row>/iso.json, ROW.md or catalog_check.json); every planted
#  control is the one the plan names.
#      [legperms.1] the retired former crossed entry == planar_smirnov_dbox under
#             p3 <-> p4 by the tool's OWN --leg-perms search (CAT leg [3b]).
#      [legperms.2] row 21 (lbl3x, drawn K4 vs the AMFlow-port record): the drawn
#             graph maps onto the record under "k1 -> -k1, k3 -> -k3; p2 ->
#             p4, p4 -> p2" (leg_perm {p1: p1, p2: p4, p3: p3, p4: p2}, 126
#             relabelings tried, 24 leg permutations, cyclic order preserved);
#             edited-line control "k - p4 - b" -> "k - p3 - b" is a K4 with
#             the legs re-seated: every relabeling found breaks the drawn
#             cyclic leg order (the census field flips true -> false); the
#             graph-breaking control "k - p4 - b" -> "k - p4 + b" matches
#             nothing; positive control p1 <-> p2 exchanged in every drawn
#             string matches under the composed map.
#      [legperms.3] row 5 (pentagon, p5 dependent): exactly two relabelings map the
#             family onto the drawing, identity and the p5-moving reflection
#             "k1 -> -k1; p1 -> p5, p2 -> p4, p4 -> p2, p5 -> p1"; without
#             re-imposition only identity is found (the former search);
#             negative control "k1 + p1 + p2" -> "k1 + p1 + p3" matches nothing.
#      [legperms.4] the four rc-1 specimens (row 20 family + drawn, row 6 retired,
#             Smirnov 2box as printed) run rc 0 under --leg-perms with
#             catalog_matches []; where the family reads fewer legs than the
#             4-leg entries (row 20: p1..p3 spelled; row 6 retired: p1, p2
#             spelled and the implicit p3 completed) the LEG-COUNT SKIP
#             warning names the 4-leg entries; the Smirnov yaml's own label
#             keys give four legs as data, so the 4-leg entries are searched
#             (no match) and nothing of the same line count is skipped; the
#             fingerprint equals the run without --leg-perms.
#      [legperms.5] rows 15 / 16 (ice-cream cones, composite legs P12 / P56 with
#             P34 eliminated): the reader alone gives [P12, P56] and no rule;
#             normalized to the drawn dictionary (p1 = P56, p2 = P34, p3 =
#             -P12) "identity" / "k1 -> -k1, k2 -> -k2", 48 relabelings in
#             the full search -- and the same from the family as load_family
#             completes it (P34 added, P12 resolved outgoing), no rewriting.
#      [legperms.6] unit controls: _dependent_leg_rows / _reimpose_conservation by
#             value; a leg permutation past the vector's block is refused by
#             name (ValueError), never IndexError; the leg-count guard returns
#             None with the skip named; label-preserving restriction counts.
#      Fixtures live under tests/fixtures/<row>/ (PROVENANCE.json beside
#      them); an absent fixture is a named SKIP, a drifted one a FAIL.
# ---------------------------------------------------------------------------
def _selftest_legperms(verbose=True):
    """Return {"ok", "executed", "skipped", "expected_fail", "skipped_names",
    "results"}; prints every sub-leg by name when verbose."""
    import hashlib
    import subprocess
    say = print if verbose else (lambda *a, **k: None)
    here = os.path.dirname(os.path.abspath(__file__))
    fxroot = os.path.join(here, "tests", "fixtures")
    FIX = {
        "row21/family_of_record.json":
            "209d2d1741147557a6b7529415606084e3a57660beb6b57b89b1c8bec1dfe42b",
        "row21/drawn_graph.yaml":
            "c307f8008add78e64099b7fb9e7948f91a81c5646bb75d393d9a47a7efea6acf",
        "row05/family_of_record.yaml":
            "0ec539a5fec61bdf5f3c0e5c59f9cda7f8044469200559ec21aadd32d8a8bccc",
        "row05/drawn_graph.yaml":
            "33984b814f95f701aef84f7965acbdf029247cda5e06b73a6ec2343f92e6599d",
        "row20/family_of_record.yaml":
            "d4aef7e8609a1f653af4fa6ba04f1d870c50be86d60ce87a20f1690d4bb8a993",
        "row20/drawn_graph.yaml":
            "c756f4b316e8953f8736a4c3b81cb88b4d12a4d6603fed6da2ecb65a2cebfae9",
        "row06/retired_sevenline_vertex2L/family_of_record_transcribed.yaml":
            "e759606f59374bd60bf8593867ac650a44990107466b46802f550fcb5c1649d5",
        "row02/smirnov_2box_as_printed.yaml":
            "93717fc0180fd92bd72e8956fcbfe745d527afec60073d171da8d347e200d522",
        "row15/family_of_record_kira_icc_eqmass.yaml":
            "9d0e526a983173045d5571fd4dd05148c306ca176b414d79353ead17054b82e2",
        "row15/drawn_graph.yaml":
            "2673aab9b205d47eaf2378bb50385b8b6b17cce690f4f4a1925f7a8c580aac34",
        "row16/family_of_record_kira_iccg_generic.yaml":
            "5be2fe7689e8dee6aab40bd7ab9d21f1efa2e7655d80333f0d87c86dfca1e581",
        "row16/drawn_graph.yaml":
            "0d09c246ad19868783db9bb3ff8eeae721aea2f48abcf79da3fb500fed806139",
    }
    # census objects (tests/fixtures/row21/iso.json "labeled_isomorphism_drawn_vs_record";
    # row05 "family_of_record_vs_drawn"; row15, row16 iso.json "pairs"[0]; the row-6
    # replaced family's leg reading; tests/fixtures/row02/catalog_check.json)
    R21_RELABEL = "k1 -> -k1, k3 -> -k3; p2 -> p4, p4 -> p2"
    R21_LEG_PERM = {"p1": "p1", "p2": "p4", "p3": "p3", "p4": "p2"}
    R21_TRIED, R21_LEG_PERMS = 126, 24
    R21_LEGS = (["p1", "p2", "p3", "p4"], {"p4": "-p1-p2-p3"})
    R05_N_FOUND, R05_FIRST = 2, "identity"
    R05_P5_MAP = "k1 -> -k1; p1 -> p5, p2 -> p4, p4 -> p2, p5 -> p1"
    R15_RELABEL, R16_RELABEL, R1516_TESTED = "identity", "k1 -> -k1, k2 -> -k2", 48
    RETIRED_MAP = ("p3 -> p4", "p4 -> p3")
    FOUR_LEG_ENTRIES = ["planar_smirnov_dbox", "crossed_tausk_dbox"]

    out = {"ok": True, "executed": 0, "skipped": 0, "expected_fail": 0,
           "skipped_names": [], "results": {}}

    def leg(name, ok, detail):
        out["executed"] += 1
        out["ok"] = out["ok"] and bool(ok)
        out["results"][name] = {"ok": bool(ok), "detail": detail}
        say(f"   {name}: {detail}  ->  {'PASS' if ok else 'FAIL'}")

    def fixture(rel):
        path = os.path.join(fxroot, rel)
        if not os.path.exists(path):
            out["skipped"] += 1
            out["skipped_names"].append(rel)
            say(f"   fixture tests/fixtures/{rel} SKIPPED (absent)")
            return None
        got = hashlib.sha256(open(path, "rb").read()).hexdigest()
        if got != FIX[rel]:
            leg(f"fixture pin {rel}", False,
                f"sha256 {got[:16]} != pinned {FIX[rel][:16]} (drift refused)")
            return None
        return path

    def fam_from(name, spec):
        return Family(name, spec["loops"], spec["exts"], spec["ext_subs"],
                      spec["propagators"], "memory:" + name)

    def with_legs(fam, legs):
        f2 = Family(fam.name, fam.loops, legs[0], legs[1], fam.propagators, fam.source,
                    fam.physical, nu=fam.nu, leg_virt=fam.leg_virt,
                    mass_values=fam.mass_values)
        return f2

    def edited(fam, old, new):
        props = [(e.replace(old, new), m) for e, m in fam.propagators]
        assert props != list(fam.propagators), f"edit {old!r} -> {new!r} changed nothing"
        return Family(fam.name + "_edited", fam.loops, fam.exts, fam.ext_subs, props,
                      "memory:edited", fam.physical)

    def search(A, B, find_all=False, dep=True):
        det = {}
        r = isomorphism_to(momentum_vectors(A)[0], momentum_vectors(B)[0], len(A.loops),
                           len(A.exts), try_leg_perms=True,
                           dep_rows=_dependent_leg_rows(B) if dep else None,
                           find_all=find_all, details=det)
        return r, det

    def cyclic_ok(legp):
        n = len(legp)
        seq = list(legp)
        for r in range(n):
            rot = seq[r:] + seq[:r]
            if rot == list(range(n)) or rot == list(range(n))[::-1]:
                return True
        return False

    say("\n[legperms.1] retired former crossed entry by the tool's own --leg-perms search:")
    ret = CATALOG_RETIRED["planar_smirnov_dbox_legs34_interchanged"]
    smi = fam_from("planar_smirnov_dbox", CATALOG_SOURCES["planar_smirnov_dbox"])
    rel, det = search(fam_from("retired", ret), smi)
    leg("[legperms.1] retired == planar_smirnov_dbox under p3<->p4 (conservation re-imposed)",
        rel is not None and all(m in rel for m in RETIRED_MAP),
        f"relabeling={rel!r} tried={det['n_relabelings_tried']} leg_perms={det['n_leg_perms']}")
    rel_old, _ = search(fam_from("retired", ret), smi, dep=False)
    leg("[legperms.1] control: without re-imposition the p3<->p4 map is NOT found (the former search)",
        rel_old is None, f"relabeling without dep_rows = {rel_old!r}")

    say("\n[legperms.2] row 21 drawn K4 vs the AMFlow-port record (legs by object):")
    p_rec, p_dr = fixture("row21/family_of_record.json"), fixture("row21/drawn_graph.yaml")
    if p_rec and p_dr:
        rec = load_family(p_rec)[0]
        drawn = with_legs(load_family(p_dr)[0], R21_LEGS)
        rel, det = search(drawn, rec)
        legmap = ({drawn.exts[j]: rec.exts[det["leg_perm"][j]] for j in range(4)}
                  if det["leg_perm"] else None)
        leg("[legperms.2] relabeling drawn -> record == census string",
            rel == R21_RELABEL and legmap == R21_LEG_PERM and
            det["n_relabelings_tried"] == R21_TRIED and det["n_leg_perms"] == R21_LEG_PERMS
            and cyclic_ok(det["leg_perm"]),
            f"relabeling={rel!r} (expect {R21_RELABEL!r}) leg_perm={legmap} "
            f"(expect {R21_LEG_PERM}) tried={det['n_relabelings_tried']} (expect {R21_TRIED}) "
            f"leg_perms={det['n_leg_perms']} (expect {R21_LEG_PERMS}) "
            f"cyclic_order_preserved={cyclic_ok(det['leg_perm']) if det['leg_perm'] else None}")
        rel_old, _ = search(drawn, rec, dep=False)
        leg("[legperms.2] control: without re-imposition no relabeling is found (the former search)",
            rel_old is None, f"relabeling without dep_rows = {rel_old!r}")
        # the plan's edited-line control 'k - p4 - b' -> 'k - p3 - b' is itself a
        # K4 with the legs re-seated (it realizes, masses match): the field that
        # flips is the census's leg_perm_cyclic_order_preserved (true -> false
        # on every match); the graph-breaking edit 'k - p4 - b' -> 'k - p4 + b'
        # flips the match itself.
        neg = edited(drawn, "k - p4 - b", "k - p3 - b")
        allm, detn = search(neg, rec, find_all=True)
        cyc = [cyclic_ok(m["leg_perm"]) for m in detn["matches"]]
        leg("[legperms.2] edited-line control 'k - p4 - b' -> 'k - p3 - b': every relabeling found "
            "breaks the drawn cyclic leg order (census field true -> false)",
            len(allm) >= 1 and not any(cyc),
            f"matches={allm} cyclic_order_preserved={cyc} tried={detn['n_relabelings_tried']}")
        neg2 = edited(drawn, "k - p4 - b", "k - p4 + b")
        allm2, detn2 = search(neg2, rec, find_all=True)
        leg("[legperms.2] negative control 'k - p4 - b' -> 'k - p4 + b' (not a graph): matches under "
            "no signed relabeling + leg permutation",
            allm2 == [] and build_graph(neg2)[1] is False,
            f"matches={allm2} tried={detn2['n_relabelings_tried']} realizes={build_graph(neg2)[1]}")
        swapped = Family(drawn.name + "_p1p2", drawn.loops, drawn.exts, drawn.ext_subs,
                         [(e.replace("p1", "@").replace("p2", "p1").replace("@", "p2"), m)
                          for e, m in drawn.propagators], "memory:p1p2", drawn.physical)
        allm, detp = search(swapped, rec, find_all=True)
        composed = {"p1": "p4", "p2": "p1", "p3": "p3", "p4": "p2"}   # (p1<->p2) then the map above
        maps = [{swapped.exts[j]: rec.exts[m["leg_perm"][j]] for j in range(4)}
                for m in detp["matches"]]
        leg("[legperms.2] positive control: p1<->p2 exchanged in every drawn string matches under "
            "the composed map",
            composed in maps and len(allm) >= 1,
            f"matches={allm} leg maps={maps} composed map {composed} present: {composed in maps}")

    say("\n[legperms.3] row 5 pentagon (p5 dependent): every relabeling family -> drawing:")
    p_f5, p_d5 = fixture("row05/family_of_record.yaml"), fixture("row05/drawn_graph.yaml")
    if p_f5 and p_d5:
        f5, d5 = load_family(p_f5)[0], load_family(p_d5)[0]
        allm, det = search(f5, d5, find_all=True)
        leg("[legperms.3] exactly two relabelings: identity and the p5-moving reflection",
            len(allm) == R05_N_FOUND and allm[0] == R05_FIRST and R05_P5_MAP in allm and
            f5.exts == ["p1", "p2", "p3", "p4", "p5"] and "p5" in f5.ext_subs,
            f"matches={allm} (expect {R05_N_FOUND}: {R05_FIRST!r}, {R05_P5_MAP!r}) "
            f"tried={det['n_relabelings_tried']} legs={f5.exts} ext_subs={f5.ext_subs}")
        allm_old, _ = search(f5, d5, find_all=True, dep=False)
        leg("[legperms.3] control: without re-imposition only identity is found (the former search)",
            allm_old == ["identity"], f"matches without dep_rows = {allm_old}")
        neg = edited(d5, "k1 + p1 + p2", "k1 + p1 + p3")
        allm, detn = search(f5, neg, find_all=True)
        leg("[legperms.3] negative control: 'k1 + p1 + p2' -> 'k1 + p1 + p3' matches under no leg permutation",
            allm == [], f"matches={allm} tried={detn['n_relabelings_tried']}")

    say("\n[legperms.4] the four rc-1 specimens under --leg-perms (CLI):")
    tool = os.path.abspath(__file__)
    # (fixture, legs read after leg completion, LEG-COUNT SKIP expected): row 6
    # retired reads p1, p2 and the completed implicit p3; the Smirnov yaml's
    # label keys give four legs as data (the 4-leg entries searched, not skipped)
    specimens = [("row20/family_of_record.yaml", 3, True), ("row20/drawn_graph.yaml", 3, True),
                 ("row06/retired_sevenline_vertex2L/family_of_record_transcribed.yaml", 3, True),
                 ("row02/smirnov_2box_as_printed.yaml", 4, False)]
    for rel, n_legs, skip_expected in specimens:
        p = fixture(rel)
        if not p:
            continue
        r1 = subprocess.run([sys.executable, tool, p, "--leg-perms", "--json"],
                            capture_output=True, text=True)
        r0 = subprocess.run([sys.executable, tool, p, "--json"],
                            capture_output=True, text=True)
        try:
            j1 = json.loads(r1.stdout)[0]
            j0 = json.loads(r0.stdout)[0]
        except Exception as exc:
            j1 = j0 = None
            err = f"{type(exc).__name__}: {exc}; stderr tail {r1.stderr.strip().splitlines()[-1:]}"
        if j1 is None or j0 is None:
            leg(f"[legperms.4] {rel} runs rc 0 with JSON", False, err)
            continue
        warn = [w for w in j1["warnings"] if w.startswith("LEG-COUNT SKIP")]
        skipped = [s["family"] for s in j1.get("catalog_skipped_leg_count", [])]
        same = all(j1.get(k) == j0.get(k) for k in
                   ("planar", "cut_signature", "canonical_hash", "canonical_set_size",
                    "n_genuine_propagators", "n_massive_lines"))
        legs_read = len((j1.get("leg_set") or {}).get("legs") or [])
        if skip_expected:
            skip_ok = (len(warn) == 1 and all(e in warn[0] for e in FOUR_LEG_ENTRIES) and
                       f"reads {n_legs} legs" in warn[0] and
                       all(e in skipped for e in FOUR_LEG_ENTRIES))
            what = "LEG-COUNT SKIP names the 4-leg entries"
        else:
            skip_ok = warn == [] and not any(e in skipped for e in FOUR_LEG_ENTRIES)
            what = "four legs read as data, the 4-leg entries searched (no skip warning)"
        leg(f"[legperms.4] {rel}: rc 0, no catalog match, {what}, "
            "fingerprint == the run without --leg-perms",
            r1.returncode == 0 and r0.returncode == 0 and j1["catalog_matches"] == [] and
            legs_read == n_legs and skip_ok and same and "IndexError" not in r1.stderr,
            f"rc={r1.returncode}/{r0.returncode} matches={j1['catalog_matches']} "
            f"legs_read={legs_read} (expect {n_legs}) skipped={skipped} "
            f"warning={'present' if warn else 'ABSENT'} fingerprint_equal={same} "
            f"hash={j1['canonical_hash']}")

    say("\n[legperms.5] rows 15 / 16 ice-cream cones (composite legs, normalized to the drawn dictionary):")
    for rel_f, rel_d, expect in (("row15/family_of_record_kira_icc_eqmass.yaml", "row15/drawn_graph.yaml", R15_RELABEL),
                                 ("row16/family_of_record_kira_iccg_generic.yaml", "row16/drawn_graph.yaml", R16_RELABEL)):
        pf, pd = fixture(rel_f), fixture(rel_d)
        if not (pf and pd):
            continue
        Fr, D = load_kira_yaml(pf)[0], load_family(pd)[0]     # Fr: the reader alone
        # drawn header dictionary: p1 = P56, p2 = P34, p3 = -P12, all incoming,
        # p3 = -p1 - p2  =>  P56 = p1, P12 = p1 + p2 (P34 = P12 - P56 eliminated)
        Fn = Family(Fr.name + "_norm", Fr.loops, ["p1", "p2", "p3"], {"p3": "-p1 - p2"},
                    [(e.replace("P56", "(p1)").replace("P12", "(p1 + p2)"), m)
                     for e, m in Fr.propagators], "memory:norm", Fr.physical)
        Dn = with_legs(D, (["p1", "p2", "p3"], {"p3": "-p1 - p2"}))
        rel, det = search(Dn, Fn)
        allm, deta = search(Dn, Fn, find_all=True)
        leg(f"[legperms.5] {rel_d} -> {rel_f}: relabeling == census, {R1516_TESTED} relabelings in the full search",
            rel == expect and allm == [expect] and deta["n_relabelings_tried"] == R1516_TESTED
            and Fr.exts == ["P12", "P56"] and Fr.ext_subs == {},
            f"relabeling={rel!r} (expect {expect!r}) all={allm} tried_full={deta['n_relabelings_tried']} "
            f"family legs as the reader gives them={Fr.exts} ext_subs={Fr.ext_subs}")
        # the completed family (load_family: P34 added, P12 resolved outgoing,
        # basis P56 = p1, P34 = p2, -P12 = p3) against the drawing as read;
        # the completion resolves the outgoing leg through the realizer
        why = _require_nx(f"[legperms.5] {rel_d} -> {rel_f} as load_family completes it")
        if why:
            out["skipped"] += 1
            out["skipped_names"].append(why)
            say(f"   {why}  ->  SKIPPED")
            continue
        F = load_family(pf)[0]
        rel2, _det2 = search(D, F)
        allm2, deta2 = search(D, F, find_all=True)
        leg(f"[legperms.5] {rel_d} -> {rel_f} as load_family completes it: the same relabeling, no rewriting",
            rel2 == expect and allm2 == [expect] and deta2["n_relabelings_tried"] == R1516_TESTED
            and F.exts == ["P56", "P34", "P12"] and F.ext_subs == {"P12": "-P56 - P34"},
            f"relabeling={rel2!r} (expect {expect!r}) all={allm2} tried_full={deta2['n_relabelings_tried']} "
            f"completed legs={F.exts} ext_subs={F.ext_subs}")

    say("\n[legperms.6] unit controls:")
    four = Family("four", ["k1", "k2"], ["p1", "p2", "p3", "p4"], {"p4": "-p1 - p2 - p3"}, [], "m")
    five = Family("five", ["k1"], ["p1", "p2", "p3", "p4", "p5"], {"p5": "-p1 - p2 - p3 - p4"}, [], "m")
    comp = Family("comp", ["l1", "l2"], ["P12", "P56", "P34"], {"P34": "P12 - P56"}, [], "m")
    rows4, rows5, rowsc = _dependent_leg_rows(four), _dependent_leg_rows(five), _dependent_leg_rows(comp)
    v = (sp.Integer(1), sp.Integer(0), sp.Integer(0), sp.Integer(0), sp.Integer(0), sp.Integer(1))
    re_v = _reimpose_conservation(v, 2, rows4)
    leg("[legperms.6] _dependent_leg_rows / _reimpose_conservation by value",
        rows4 == {3: [-1, -1, -1, 0]} and rows5 == {4: [-1, -1, -1, -1, 0]} and
        rowsc == {2: [1, -1, 0]} and re_v == (1, 0, -1, -1, -1, 0),
        f"rows4={rows4} rows5={rows5} rowsc={rowsc} k1+p4 -> {re_v}")
    try:
        # a 2-leg vector (k1 + p1) under a 4-leg permutation sending p1 -> p4
        _apply_signed_perm_to_vec((sp.Integer(1), sp.Integer(0), sp.Integer(1), sp.Integer(0)),
                                  (0, 1), (1, 1), (3, 2, 1, 0), 2, 4)
        named = "no error"
    except ValueError as exc:
        named = "ValueError: " + str(exc)
    except IndexError as exc:
        named = "IndexError: " + str(exc)
    leg("[legperms.6] a leg permutation past the vector's block is refused by name (ValueError)",
        named.startswith("ValueError") and "external legs" in named, named)
    det = {}
    guard = isomorphism_to([(sp.Integer(1), sp.Integer(1), sp.Integer(0))],
                           [(sp.Integer(1), sp.Integer(1), sp.Integer(0), sp.Integer(0))],
                           1, 2, try_leg_perms=True, details=det)
    leg("[legperms.6] leg-count guard: different leg counts -> None with the skip named, no IndexError",
        guard is None and det["leg_count_skip"] is not None and "differ" in det["leg_count_skip"],
        f"result={guard!r} skip={det['leg_count_skip']!r}")
    det = {}
    isomorphism_to(momentum_vectors(smi)[0], momentum_vectors(smi)[0], 2, 4, try_leg_perms=True,
                   dep_rows=_dependent_leg_rows(smi),
                   leg_labels=(["0", "0", "offshell", "offshell"], ["0", "0", "offshell", "offshell"]),
                   find_all=True, details=det)
    leg("[legperms.6] label-preserving restriction: 4 of 24 leg permutations searched for labels (0,0,off,off)",
        det["n_leg_perms"] == 4, f"n_leg_perms={det['n_leg_perms']} (expect 4 = 2! * 2!)")

    say(f"\n   legperms leg: {out['executed']} executed / {out['skipped']} skipped / "
        f"{out['expected_fail']} expected-fail"
        + (f" (skipped: {out['skipped_names']})" if out["skipped_names"] else ""))
    return out


#  Self-test legs contributed outside self_test(): (leg id, description,
#  callable of no arguments).  Each callable returns its own counts,
#  {"ok": bool, "executed": int, "skipped": int | [names],
#   "expected_fail": int (optional), "skipped_names": [...] (optional),
#   "results": {...} (optional)}; the dispatcher folds executed / skipped /
#  expected-fail sub-legs into the summary and prints ALL PASS only when
#  nothing failed, nothing was skipped and no expected-fail sub-leg is open
#  (rc 2 on any skip, rc 0 with declared expected-fails).  ONE list for the
#  whole package:
#  the realizer leg (_selftest_realizer), the catalog leg (_selftest_catalog),
#  the kinematics leg (_selftest_kinematics), the leg-perms leg
#  (_selftest_legperms: momentum conservation re-imposed after the
#  permutation) and the later legs register here.
def _selftest_kinematics():
    """
    Self-test leg of the kinematics reader and the leg completion
    (read_kinematics_yaml / leg_virtualities / complete_legs / the
    LEG-SET-INFERRED line): runs tests/test_kinematics.py's battery -- every
    case reproduces a census object BY VALUE on a vendored, sha-pinned
    fixture (row 20 planar False / {s: 2, t: 3, u: 3} / size 7 with its
    kinematics and the named trap without; rows 19/31/32/34 {ppB: 0} and
    173439ca7eee / 778036c88aae; rows 8-14 two-point cuts 3 / 2; rows 15-17
    8a73188ccb5a with P12 resolved outgoing; row 23 bb2ae2aad8f4; row 24 five
    legs with the census's virtualities; rows 1/2/4/7 label keys read) and
    every planted control flips or is refused by name (wrong-sign composite
    rule -> NON-GRAPH; [p, q] for a single current -> REFUSED; p5 rule
    deleted -> LEG-SET-COMPLETED; off-shell leg moved -> virtualities change;
    label keys removed -> LEG-SET-INFERRED).  Registered in
    EXTRA_SELF_TEST_LEGS; returns {"ok", "executed", "skipped": [names]}.
    """
    here = os.path.dirname(os.path.abspath(__file__))
    mod_path = os.path.join(here, "tests", "test_kinematics.py")
    if not os.path.exists(mod_path):
        print("   [kinematics] tests/test_kinematics.py absent -> SKIP", flush=True)
        return {"ok": True, "executed": 0, "skipped": ["tests/test_kinematics.py absent"]}
    if not os.path.isdir(_FIXTURE_DIR):
        print("   [kinematics] tests/fixtures/ absent -> SKIP", flush=True)
        return {"ok": True, "executed": 0, "skipped": ["tests/fixtures/ absent"]}
    import importlib.util
    spec = importlib.util.spec_from_file_location("topology_audit_test_kinematics", mod_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    print("\n[kinematics] kinematics: legs as data, completion, LEG-SET-INFERRED:", flush=True)
    rows = mod.run_battery(verbose=True)
    n_pass = sum(1 for r in rows if r["status"] == "PASS")
    fails = [r["name"] for r in rows if r["status"] == "FAIL"]
    skips = [r["name"] + ": " + r["detail"] for r in rows if r["status"] == "SKIP"]
    print(f"   [kinematics] battery: {len(rows)} cases, {n_pass} PASS, {len(fails)} FAIL, "
          f"{len(skips)} SKIP  ->  "
          f"{'PASS' if not fails and not skips else ('FAIL' if fails else 'PASS WITH SKIPS')}",
          flush=True)
    for name in fails:
        print(f"   [kinematics] FAIL {name}", flush=True)
    return {"ok": not fails, "executed": n_pass + len(fails), "skipped": skips,
            "results": {"kinematics_battery": rows}}


def _selftest_labels(verbose=True):
    """
    Self-test leg of the LABELS (masses on lines, multiplicities = propagator
    powers, leg virtuality classes) in the canonical set and the isomorphism
    (momentum_vectors .nu / canonical_set labeled Counter / labeled_set /
    mass_multiset / labeled_isomorphism / the fingerprint's mass, nu and leg
    label fields / the catalog match's mass pattern): runs tests/test_labels.py's
    battery -- every case reproduces a census object BY VALUE on a vendored,
    sha-pinned fixture (the row-12 3x3 panel matrix: mass multisets diagonal,
    every off-diagonal cell fails the labeled compare, every cell passes the
    mass-blind one; row 14's record kite labeled identity and the equal-mass
    control failing; row 18 family 2 identity with the leg labels and family 5
    mass-blind only; row 4 reading (b) momentum-isomorphic true / labelled
    false; row 35 identity with the W legs off shell and the moved-leg control
    failing; rows 21 / 28 / 29 / 33 with the record's multiplicities carried
    and the nu-mutated copies failing; rows 19 / 34 labeled identity and the
    row-31 mass-permuted control failing; the row-24 cured drawn pair
    separating on leg virtualities alone; rows 1 / 7 catalog matches naming
    the mass pattern) and the identities the signed-permutation class does not
    reach (the row-12 diagonal, rows 20 / 31 / 32, the row-14 half-turn
    control), which the affine / realized-graph isomorphism (the [isomorphism] unit:
    affine_iso_search, realized_graph_iso) lifts -- each of those cases asserts
    that the signed class finds nothing and then the census's map from the
    affine search and the realized-graph test, so no EXPECTED-FAIL sub-leg is
    open (no BATTERY entry of tests/test_labels.py carries a strict mark; the
    mechanism stays for any later one).  Registered in EXTRA_SELF_TEST_LEGS;
    returns {"ok", "executed", "skipped": [names], "expected_fail"}.
    """
    here = os.path.dirname(os.path.abspath(__file__))
    mod_path = os.path.join(here, "tests", "test_labels.py")
    say = print if verbose else (lambda *a, **k: None)
    if not os.path.exists(mod_path):
        say("   [labels] tests/test_labels.py absent -> SKIP", flush=True)
        return {"ok": True, "executed": 0, "skipped": ["tests/test_labels.py absent"]}
    if not os.path.isdir(_FIXTURE_DIR):
        say("   [labels] tests/fixtures/ absent -> SKIP", flush=True)
        return {"ok": True, "executed": 0, "skipped": ["tests/fixtures/ absent"]}
    import importlib.util
    spec = importlib.util.spec_from_file_location("topology_audit_test_labels", mod_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    say("\n[labels] labels: masses, multiplicities and leg classes in the canonical set and "
        "the isomorphism:", flush=True)
    rows = mod.run_battery(verbose=verbose)
    n_pass = sum(1 for r in rows if r["status"] == "PASS")
    fails = [r["name"] for r in rows if r["status"] == "FAIL"]
    skips = [r["name"] + ": " + r["detail"] for r in rows if r["status"] == "SKIP"]
    xfails = [r["name"] for r in rows if r["status"] == "EXPECTED-FAIL"]
    say(f"   [labels] battery: {len(rows)} cases, {n_pass} PASS, {len(fails)} FAIL, "
        f"{len(skips)} SKIP, {len(xfails)} EXPECTED-FAIL  ->  "
        f"{'PASS' if not fails and not skips else ('FAIL' if fails else 'PASS WITH SKIPS')}",
        flush=True)
    for name in fails:
        say(f"   [labels] FAIL {name}", flush=True)
    for name in xfails:
        say(f"   [labels] open EXPECTED-FAIL {name}: a strict mark declared in "
            "tests/test_labels.py BATTERY that still raises ExpectedFail (not lifted "
            "by the labeled / affine isomorphism); counted, holds back ALL PASS",
            flush=True)
    return {"ok": not fails, "executed": n_pass + len(fails) + len(xfails), "skipped": skips,
            "expected_fail": len(xfails),
            "results": {"labels_battery": rows}}


def _selftest_planarity():
    """
    Self-test leg of the planarity surface (planarity / planarity_at_infinity
    / the closure order and its source / the Mandelstam declaration report):
    runs tests/test_planarity.py's battery -- every case reproduces a census
    object BY VALUE on a vendored, sha-pinned fixture (row 24 retired and
    cured: the census's planarity_on_realized_graph dict {closure False,
    abstract True, legs joined at infinity False, order p1..p5,
    discriminating True}; rows 15/16/17 with their kinematics planar True and
    row 23 planar True with the "not discriminating (N < 4 legs)" flag by
    name in the method string, the PLANARITY-NOT-DISCRIMINATING line and the
    legs-at-infinity field; row 18 fam2 (three legs, inferred) withheld with
    the report; row 20 with its kinematics planar False {s: 2, t: 3, u: 3}
    unchanged; row 2's family planar False in the canonical order and planar
    True in its own t = (p1+p3)^2 order p1,p2,p4,p3, reported under
    PLANARITY-LEG-ORDER; the theta crossed box non-planar in every enumerated
    leg order and at infinity; row 28 with its kinematics planar True in the
    drawn order p1,p2,p3,p4 while its t = (p1+p3)^2 declaration is reported,
    not enforced; row 5 planar with its pySecDec leg order as the closure
    order) and every planted control flips or is refused by name (a declared
    cyclic order p1,p2,p4,p3 on the one-loop box -> planar False; a leg order
    that is not a permutation -> REFUSED; K3,3 and prism three-point graphs
    where the closure test decides False under the flag; the catalog's 2- and
    3-leg entries flagged, its 4-leg entries not).  Registered in
    EXTRA_SELF_TEST_LEGS; returns {"ok", "executed", "skipped": [names]}.
    """
    here = os.path.dirname(os.path.abspath(__file__))
    mod_path = os.path.join(here, "tests", "test_planarity.py")
    if not os.path.exists(mod_path):
        print("   [planarity] tests/test_planarity.py absent -> SKIP", flush=True)
        return {"ok": True, "executed": 0, "skipped": ["tests/test_planarity.py absent"]}
    if not os.path.isdir(_FIXTURE_DIR):
        print("   [planarity] tests/fixtures/ absent -> SKIP", flush=True)
        return {"ok": True, "executed": 0, "skipped": ["tests/fixtures/ absent"]}
    import importlib.util
    spec = importlib.util.spec_from_file_location("topology_audit_test_planarity", mod_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    print("\n[planarity] planarity: closure order by name, the N < 4 flag, legs at infinity, "
          "enumerated orders, the Mandelstam declaration reported:", flush=True)
    rows = mod.run_battery(verbose=True)
    n_pass = sum(1 for r in rows if r["status"] == "PASS")
    fails = [r["name"] for r in rows if r["status"] == "FAIL"]
    skips = [r["name"] + ": " + r["detail"] for r in rows if r["status"] == "SKIP"]
    print(f"   [planarity] battery: {len(rows)} cases, {n_pass} PASS, {len(fails)} FAIL, "
          f"{len(skips)} SKIP  ->  "
          f"{'PASS' if not fails and not skips else ('FAIL' if fails else 'PASS WITH SKIPS')}",
          flush=True)
    for name in fails:
        print(f"   [planarity] FAIL {name}", flush=True)
    return {"ok": not fails, "executed": n_pass + len(fails), "skipped": skips,
            "results": {"planarity_battery": rows}}


def _selftest_isomorphism():
    """
    Self-test leg of the ISOMORPHISM beyond signed loop permutations
    (unimodular_maps / _affine_search_vectors / affine_iso_search /
    isomorphism_to(affine=True) / realized_graph_iso / _canonical_form /
    realized_canonical_hash / the fingerprint's canonical_hash_realized and
    the catalog match's relabeling class): runs tests/test_isomorphism.py's
    battery -- every case reproduces a census object BY VALUE on a vendored,
    sha-pinned fixture (rows 1 / 2 / 7 graph-isomorphic to their drawings with
    the realized hash equal on both sides and the routed hashes different; row
    20's affine map ['l1 -> k1 -k2 +p3', 'l2 -> -k2 +p1 +p2 +p3']; rows 19 /
    31 / 32 / 34 and 9-12 with their maps; the row-12 diagonal, the row-14
    half-turn control and rows 31 / 32 labeled maps that were EXPECTED-FAIL
    under the signed class; row 30 adj vs opp NOT isomorphic; row 4 NOT
    isomorphic (V,E 10,12 vs 12,14); the drawn four-point triple box matching
    planar_triple_box_3l) and every planted control flips or is refused by
    name (a moved rung / a different mass on one line -> not isomorphic under
    every U, the search exhausts and says so; a planted shift that does not
    preserve the labeled set -> refused; a relabeled copy of a graph -> the
    same canonical hash).  Registered in EXTRA_SELF_TEST_LEGS; returns
    {"ok", "executed", "skipped": [names]}.
    """
    here = os.path.dirname(os.path.abspath(__file__))
    mod_path = os.path.join(here, "tests", "test_isomorphism.py")
    if not os.path.exists(mod_path):
        print("   [isomorphism] tests/test_isomorphism.py absent -> SKIP", flush=True)
        return {"ok": True, "executed": 0, "skipped": ["tests/test_isomorphism.py absent"]}
    if not os.path.isdir(_FIXTURE_DIR):
        print("   [isomorphism] tests/fixtures/ absent -> SKIP", flush=True)
        return {"ok": True, "executed": 0, "skipped": ["tests/fixtures/ absent"]}
    import importlib.util
    spec = importlib.util.spec_from_file_location("topology_audit_test_isomorphism", mod_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    print("\n[isomorphism] isomorphism: the affine class k -> U k + c.p as the describer, the realized-"
          "graph isomorphism (VF2) as the identity test, the routing-independent canonical "
          "hash:", flush=True)
    rows = mod.run_battery(verbose=True)
    n_pass = sum(1 for r in rows if r["status"] == "PASS")
    fails = [r["name"] for r in rows if r["status"] == "FAIL"]
    skips = [r["name"] + ": " + r["detail"] for r in rows if r["status"] == "SKIP"]
    print(f"   [isomorphism] battery: {len(rows)} cases, {n_pass} PASS, {len(fails)} FAIL, "
          f"{len(skips)} SKIP  ->  "
          f"{'PASS' if not fails and not skips else ('FAIL' if fails else 'PASS WITH SKIPS')}",
          flush=True)
    for name in fails:
        print(f"   [isomorphism] FAIL {name}", flush=True)
    return {"ok": not fails, "executed": n_pass + len(fails), "skipped": skips,
            "results": {"isomorphism_battery": rows}}


def _selftest_compare():
    """
    Self-test leg of the COMPARE MODE (read_drawn_header / load_compare_side
    / legs_per_vertex / _mass_identifications / _cyclic_order_finding /
    compare_drawn / verdict_lines): runs tests/census_battery.py -- every
    family-level entry of the identity-census corpus (30 published families
    and their drawn graphs; manifest tests/fixtures/CENSUS_BATTERY.json) is
    compared with the tool and the verdict STRING asserted equal to the
    census's, never typed; each row-level rollup the manifest lists ("see
    families: ...") is reproduced by joining the family verdicts in the
    manifest's order; the dual cases are the manifest's `extras` and `also`
    pairs, each with its own expected string (a family against another
    reading of the drawing, or against its own edge list), and the planted
    controls (a mass moved, a leg class moved, a line dropped) flip or are
    refused by name.  A corpus row the census did not compare (no drawn
    figure) has no entry, no rollup and no verdict here: its fixture files
    serve the reader and isomorphism tests, not battery cases.  The case
    count is the manifest's, printed by the battery's own summary line,
    never written here.  A census verdict the tool cannot reproduce is a
    FAIL of the named entry.  Registered in EXTRA_SELF_TEST_LEGS; returns
    {"ok", "executed", "skipped": [names]}.
    """
    here = os.path.dirname(os.path.abspath(__file__))
    mod_path = os.path.join(here, "tests", "census_battery.py")
    if not os.path.exists(mod_path):
        print("   [compare] tests/census_battery.py absent -> SKIP", flush=True)
        return {"ok": True, "executed": 0, "skipped": ["tests/census_battery.py absent"]}
    if not os.path.isdir(_FIXTURE_DIR):
        print("   [compare] tests/fixtures/ absent -> SKIP", flush=True)
        return {"ok": True, "executed": 0, "skipped": ["tests/fixtures/ absent"]}
    import importlib.util
    spec = importlib.util.spec_from_file_location("topology_audit_census_battery", mod_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    print("\n[compare] the compare mode against the identity census: every family-level "
          "verdict string of the corpus, the rollups, the dual cases, "
          "the planted controls:", flush=True)
    rows = mod.run_battery(verbose=True)
    n_pass = sum(1 for r in rows if r["status"] == "PASS")
    fails = [r["name"] for r in rows if r["status"] == "FAIL"]
    skips = [r["name"] + ": " + r["detail"] for r in rows if r["status"] == "SKIP"]
    print(f"   [compare] battery: {len(rows)} cases, {n_pass} PASS, {len(fails)} FAIL, "
          f"{len(skips)} SKIP  ->  "
          f"{'PASS' if not fails and not skips else ('FAIL' if fails else 'PASS WITH SKIPS')}",
          flush=True)
    for name in fails:
        print(f"   [compare] FAIL {name}", flush=True)
    return {"ok": not fails, "executed": n_pass + len(fails), "skipped": skips,
            "results": {"census_battery": rows}}


def _load_loomcheck_member():
    """Import the loomcheck member (loomcheck/loomcheck.py beside this file)
    by path, so the leg runs from any cwd; returns (module, None) or
    (None, skip_reason)."""
    here = os.path.dirname(os.path.abspath(__file__))
    mod_path = os.path.join(here, "loomcheck", "loomcheck.py")
    if not os.path.exists(mod_path):
        return None, "loomcheck/loomcheck.py absent"
    why = _require_nx("loomcheck (planarity + face traversal)")
    if why:
        return None, why
    import importlib.util
    spec = importlib.util.spec_from_file_location("topology_audit_loomcheck", mod_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod, None


def _selftest_loomcheck():
    """
    Self-test leg of the LOOMCHECK member (loomcheck/: the Yangian / loom /
    fishnet applicability screen for position-space conformal integrals --
    (1) the conformal weight D at every internal vertex, numerator edges
    counted as power -1; (2) planarity of the FULL graph with the numerator
    edges as edges; (3) the level-one-momentum face condition sum_face a_e =
    (k-2) D/2 per k-gon face): runs the member's structural known-answer
    battery (loomcheck.selftest_cases) -- K5 with one numerator edge is
    non-planar and excluded, the four-point star's internal vertex weighs
    4 = D with unit powers and 5 with one propagator squared, a unit 4-cycle's
    two quadrilateral faces satisfy the face sum and a unit 5-cycle's two
    pentagon faces both violate it -- one hand-checkable case per screened
    condition, no data files; then the screen's summary invariant on the same
    graphs (yangian_class is True only when all three conditions hold).  The
    member absent or networkx absent -> SKIP by name.  Registered in
    EXTRA_SELF_TEST_LEGS; returns {"ok", "executed", "skipped": [names]}.
    """
    mod, why = _load_loomcheck_member()
    if mod is None:
        print(f"   [loomcheck] {why} -> SKIP", flush=True)
        return {"ok": True, "executed": 0, "skipped": [why]}
    print("\n[loomcheck] the applicability screen's structural battery (internal-vertex "
          "weight, full-graph planarity with numerator edges, the P-hat face sums), one "
          "hand-checkable graph per condition:", flush=True)
    rows = []
    for name, got, want in mod.selftest_cases(D=4):
        ok = (got == want)
        rows.append({"name": name, "status": "PASS" if ok else "FAIL",
                     "detail": "" if ok else f"got {got!r}, expected {want!r}"})
    # the summary invariant: yangian_class == (weight and planar and 0 violations)
    graphs = {
        "K5+numerator": {**{(i, j): 1 for i in range(1, 6) for j in range(i + 1, 6)}, (1, 2): -1},
        "star": {(1, 5): 1, (2, 5): 1, (3, 5): 1, (4, 5): 1},
        "4-cycle": {(1, 2): 1, (2, 3): 1, (3, 4): 1, (1, 4): 1},
        "5-cycle": {(1, 2): 1, (2, 3): 1, (3, 4): 1, (4, 5): 1, (1, 5): 1},
    }
    for gname, powers in graphs.items():
        r = mod.screen_graph(powers, D=4)
        want = bool(r["conformal_internal_weight"] and r["planar_full_graph"]
                    and r["phat_face_violations"] == 0)
        ok = (r["yangian_class"] == want) and set(r) >= {
            "conformal_internal_weight", "planar_full_graph", "numerator_edges",
            "faces", "phat_face_violations", "yangian_class"}
        rows.append({"name": f"{gname}: yangian_class == all three conditions, keys complete",
                     "status": "PASS" if ok else "FAIL",
                     "detail": "" if ok else f"screen_graph -> {r!r}"})
    for r in rows:
        print(f"   {r['status']:4s} {r['name']}" + (f"  ({r['detail']})" if r["detail"] else ""),
              flush=True)
    n_pass = sum(1 for r in rows if r["status"] == "PASS")
    fails = [r["name"] for r in rows if r["status"] == "FAIL"]
    print(f"   [loomcheck] battery: {len(rows)} checks, {n_pass} PASS, {len(fails)} FAIL  ->  "
          f"{'PASS' if not fails else 'FAIL'}", flush=True)
    return {"ok": not fails, "executed": len(rows), "skipped": [],
            "results": {"loomcheck_battery": rows}}


EXTRA_SELF_TEST_LEGS = [
    ("realizer", "vertex realizer: connected covers, V = E - L + 1, non-graph verdict",
     _selftest_realizer),
    ("catalog", "catalog: theta crossed box, retired entry as the named negative, "
                "UD three-point ladder entry, per-leg virtuality classes",
     _selftest_catalog),
    ("kinematics", "kinematics: the external-leg set as data (kinematics.yaml, yaml keys, JSON "
           "legs), completion (pB, p_N, composite legs), LEG-SET-INFERRED fallback",
     _selftest_kinematics),
    ("legperms", "leg permutations: momentum conservation re-imposed after the permutation, "
           "leg-count guard (the rc-1 specimens), label-preserving restriction",
     _selftest_legperms),
    ("labels", "labels: masses, multiplicities (propagator powers) and leg virtuality classes "
           "in the canonical set, the isomorphism and the catalog match's mass pattern",
     _selftest_labels),
    ("planarity", "planarity: closure order by name (declared cyclic order / canonical index), "
             "the 'not discriminating (N < 4 legs)' flag, legs joined at infinity, "
             "enumerated leg orders, the Mandelstam declaration reported",
     _selftest_planarity),
    ("isomorphism", "isomorphism beyond signed permutations: the affine class k -> U k + c.p (U "
           "unimodular, c solved) as the describer, the realized-graph isomorphism (VF2 with "
           "masses, multiplicities, leg classes) as the identity test, the routing-independent "
           "canonical hash beside the routed one",
     _selftest_isomorphism),
    ("compare", "compare mode: a family against its drawn graph -> the identity census's "
                "verdict string (IDENTITY-PASS / FAMILY-MISMATCH / NON-GRAPH / NOT-CHECKABLE), "
                "every family-level entry of the corpus reproduced exactly, "
                "the rollups, the dual cases, the cyclic-order sub-finding, legs per vertex",
     _selftest_compare),
    ("loomcheck", "loomcheck member: the Yangian/loom/fishnet applicability screen's structural "
                  "battery (internal-vertex conformal weight, full-graph planarity with numerator "
                  "edges, the P-hat face condition), one hand-checkable graph per condition",
     _selftest_loomcheck),
]


def _print_audit(fp, verdict=None):
    if verdict is not None:
        # the compare mode's verdict line first, then its evidence, then the audit
        for ln in verdict_lines(verdict):
            print(ln)
    print(f"   name={fp['name']}  loops={fp['n_loops']}  "
          f"genuine_props={fp['n_genuine_propagators']}  "
          f"isp={fp['n_isp_or_dotproduct']}  massive={fp['n_massive_lines']}")
    print(f"   planar={fp['planar']}  ({fp['planarity_method']})")
    print(f"   cut_signature={fp['cut_signature']}  "
          f"multiset={fp['cut_signature_multiset']}  "
          f"{fp.get('cut_symmetry','')}")
    print(f"   canonical_hash={fp['canonical_hash']}  (routed)   "
          f"canonical_hash_realized={fp.get('canonical_hash_realized')}  (realized graph, "
          f"routing-independent"
          + (f"; {fp['canonical_hash_realized_reason']}" if fp.get('canonical_hash_realized_reason')
             and fp.get('canonical_hash_realized') is None else "") + ")")
    if fp["catalog_matches"]:
        for m in fp["catalog_matches"]:
            print(f"   >> ISO to '{m['family']}' under [{m['relabeling']}]")
    if fp["warnings"]:
        for w in fp["warnings"]:
            print(f"   !! {w}")
    else:
        print("   (no warnings)")


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Dogtag: integral-family identity checks (catch mislabeled / duplicate Feynman families; compare a family with its drawn graph).")
    ap.add_argument("family_file", nargs="?",
                    help="integralfamilies.yaml, family.jl, AMFlow-port config.json, "
                         "pySecDec graph .py, or a drawn-graph edge list")
    ap.add_argument("--name", help="select a single family by name (multi-family yaml)")
    ap.add_argument("--leg-perms", action="store_true",
                    help="also try external-leg permutations in the isomorphism search")
    ap.add_argument("--json", action="store_true", help="emit JSON only")
    ap.add_argument("--kinematics", metavar="KINEMATICS_YAML",
                    help="Kira kinematics.yaml recorded on the family "
                         "(Family.kinematics_path) for the kinematics reader")
    ap.add_argument("--drawn", metavar="DRAWN_GRAPH",
                    help="drawn-graph file (edge list, routed yaml or graphspec JSON) "
                         "compared with the family: the compare mode (compare_drawn) "
                         "emits the verdict IDENTITY-PASS / FAMILY-MISMATCH / NON-GRAPH / "
                         "NOT-CHECKABLE first, then its evidence, then the family's audit")
    ap.add_argument("--kinematics-drawn", metavar="KINEMATICS_YAML",
                    help="Kira kinematics.yaml of the drawn graph (its legs as data), "
                         "read with --drawn")
    ap.add_argument("--verdict", action="store_true",
                    help="with --drawn: print the verdict line and its evidence only "
                         "(no audit); --json carries the compare under 'verdict'")
    ap.add_argument("--integral", type=int, metavar="N",
                    help="AMFlow-port JSON: take the top sector from integrals[N] "
                         "(default: the integral with the most positive indices)")
    ap.add_argument("--self-test", "--selftest", action="store_true",
                    help="run validation suite")
    args = ap.parse_args()

    if args.self_test:
        results, rc = self_test()
        sys.exit(rc)

    if not args.family_file:
        ap.error("provide a family file or --self-test")
    if (args.kinematics_drawn or args.verdict) and not args.drawn:
        ap.error("--kinematics-drawn and --verdict need --drawn DRAWN_GRAPH")

    if args.drawn:
        # the compare mode: both sides through load_compare_side (the drawn
        # yaml's header legs, the .jl replacement virtualities; --name picks
        # the family side of a multi-family file, the battery's route), the
        # verdict emitter, then the family's audit unless --verdict
        try:
            res = compare_drawn(args.family_file, args.drawn, kinA=args.kinematics,
                                kinB=args.kinematics_drawn, integral_indexA=args.integral,
                                nameA=args.name)
        except UnknownFamilyName as exc:
            print("topology_audit: --name refused: %s" % exc, file=sys.stderr)
            sys.exit(2)
        if args.json:
            out = {"verdict": res}
            if not args.verdict:
                fams = load_family(args.family_file, args.name, kinematics=args.kinematics,
                                   integral_index=args.integral)
                out["audit"] = [audit(fam, try_leg_perms=args.leg_perms) for fam in fams]
            print(json.dumps(out, indent=2, default=str))
            sys.exit(0)
        if args.verdict:
            for ln in verdict_lines(res):
                print(ln)
            sys.exit(0)
        fams = load_family(args.family_file, args.name, kinematics=args.kinematics,
                           integral_index=args.integral)
        if not fams:
            print("no families found in", args.family_file)
            sys.exit(2)
        for i, fam in enumerate(fams):
            fp = audit(fam, try_leg_perms=args.leg_perms)
            _print_audit(fp, verdict=(res if i == 0 else None))
            print()
        sys.exit(0)

    fams = load_family(args.family_file, args.name, kinematics=args.kinematics,
                       integral_index=args.integral)
    if not fams:
        if args.name:
            # the audit path's refusal by name lists the families the file carries
            carried = [f.name for f in load_family(args.family_file, None, kinematics=args.kinematics,
                                                   integral_index=args.integral)]
            print("topology_audit: --name refused: no family found named %r in %s; the file carries: %s"
                  % (args.name, args.family_file, ", ".join(carried)), file=sys.stderr)
        else:
            print("no families found in", args.family_file)
        sys.exit(2)
    if args.kinematics:
        print(f"NOTE: --kinematics {args.kinematics} recorded on the family; its "
              "leg completion is applied by the kinematics reader, not here",
              file=sys.stderr)
    out = []
    for fam in fams:
        fp = audit(fam, try_leg_perms=args.leg_perms)
        out.append(fp)
        if not args.json:
            _print_audit(fp)
            print()
    if args.json:
        print(json.dumps(out, indent=2, default=str))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# galois member — sv-HPL vector-space layer.
"""vspace.py — allowed-space (V) machinery for the coaction cut.

Parses SOLVED conformal-integral results (published ancillary files) into
per-LS-block PURE single-valued functions over the {0,1}-letter sv spaces,
saturates the collection under the Galois derivations D3/D5/D7, and exposes
weight+parity fibers with exact annihilator functionals.

Solved feed (published sources only — quarantined I173 entries are NOT read
here; the standalone published fourloopI173.txt enters only in production
mode, never in the control):
  - 4Loop_integral_results_by_hyperlog.m : clean entries (no FAIL, no
    epsilon/gamma/Pi, no extended zz-letters)
  - fourloopI120.txt (clean standalone; FAIL in the hyperlog file)
  - fourloopI173.txt (production only; the G5 control excludes it)
  - 3Loop_integral_insvmpl.m : clean entries
Entries with extended letters (I[z, zz, ...]) live outside the {0,1} sv
space and are skipped (counted in the receipt).
"""
import json
import os
import re
import sys
from fractions import Fraction

import sympy as sp

TOOLS = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ANC = os.environ.get("GALOIS_ANC_DIR")  # solved-feed anc dir (arXiv:2607.11645 ancillary files)


def _anc():
    if not ANC:
        raise RuntimeError(
            "GALOIS_ANC_DIR is not set. vspace's solved feed = the ancillary "
            "files of arXiv:2607.11645 (4Loop_integral_results_by_hyperlog.m, "
            "fourloopI120.txt, fourloopI173.txt, 3Loop_integral_insvmpl.m) — "
            "download them and point GALOIS_ANC_DIR at that directory.")
    return ANC
if TOOLS not in sys.path:
    sys.path.insert(0, TOOLS)

from galois import core as sg
from fractions import Fraction

z, zz = sp.symbols("z zz")

FTOK = {  # f[...] -> zeta token of core.ZTOK
    (3,): "f3", (5,): "f5", (7,): "f7",
    (3, 3): "f33", (5, 3): "f53", (3, 5): "f35",
}

TARGETS = {144, 163, 167, 168, 169, 212, 214, 233, 260}


def split_top(body):
    out, depth, cur = [], 0, ""
    for ch in body:
        if ch in "{[(":
            depth += 1
        elif ch in "}])":
            depth -= 1
        if ch == "," and depth == 0:
            out.append(cur)
            cur = ""
        else:
            cur += ch
    if cur.strip():
        out.append(cur)
    return out


def load_entries(path):
    s = open(path).read()
    return split_top(s[s.index("{") + 1:s.rindex("}")])


# ---------------- entry parsing --------------------------------------------

_IPAT = re.compile(r"I\[z((?:,\s*[01])+),\s*0\]")
_FPAT = re.compile(r"f\[(\d(?:,\s*\d)*)\]")


class SkipEntry(Exception):
    pass


def parse_entry(txt):
    """Mathematica expr string -> list of blocks
    [{'prefactor': sympy-ratfunc (normalized), 'vec': {(word,mono): Fr},
      'weight': int}].  Raises SkipEntry for zz-letters / eps / FAIL."""
    t = txt.strip()
    if t == "FAIL" or "FAIL" in t and len(t) < 10:
        raise SkipEntry("FAIL")
    for bad in ("epsilon", "gamma", "Pi"):
        if bad in t:
            raise SkipEntry(bad)
    if "I[z, zz" in t or "I[zz" in t:
        raise SkipEntry("zz-letter")
    # replace I[...] and f[...] by symbols
    iw_syms = {}

    def _irep(m):
        w = tuple(int(x) for x in m.group(1).replace(" ", "").strip(",").split(","))
        name = "IW_" + "".join(map(str, w))
        iw_syms[sp.Symbol(name)] = w
        return name

    f_syms = {}

    def _frep(m):
        idx = tuple(int(x) for x in m.group(1).replace(" ", "").split(","))
        if idx not in FTOK:
            raise SkipEntry(f"unknown f{idx}")
        name = "FF_" + "_".join(map(str, idx))
        f_syms[sp.Symbol(name)] = FTOK[idx]
        return name

    t2 = _IPAT.sub(_irep, t)
    if "I[" in t2:
        raise SkipEntry("unparsed I[ form")
    t2 = _FPAT.sub(_frep, t2)
    if "f[" in t2:
        raise SkipEntry("unparsed f[ form")
    t2 = t2.replace("^", "**")
    expr = sp.sympify(t2, locals={s.name: s for s in
                                  list(iw_syms) + list(f_syms) + [z, zz]})
    expr = sp.expand(expr)
    syms = set(iw_syms) | set(f_syms)
    blocks = []          # [(numP, denP, scale-tracking via poly content)]

    def canon_pref(r):
        r = sp.cancel(sp.together(r))
        num, den = sp.fraction(r)
        pn = sp.Poly(num, z, zz)
        pd = sp.Poly(den, z, zz)
        cn = pn.LC()
        r_scale = sp.Rational(cn) / sp.Rational(pd.LC())
        key = (tuple(sorted(zip(map(tuple, pn.monoms()),
                                [sp.Rational(c / cn) for c in pn.coeffs()]))),
               tuple(sorted(zip(map(tuple, pd.monoms()),
                                [sp.Rational(c / pd.LC())
                                 for c in pd.coeffs()]))))
        return key, r_scale

    groups = {}
    for term in sp.Add.make_args(expr):
        if term == 0:
            continue
        rest, spart = term.as_independent(*syms)
        if rest == 0 or spart == 0:
            continue
        # spart: product of one optional IW and f-symbols with powers
        word = ()
        mono = {(0, 0, 0, 0, 0): Fraction(1)}
        for fac in sp.Mul.make_args(spart):
            base, exp = fac.as_base_exp()
            assert exp.is_Integer and exp > 0, (fac,)
            if base in iw_syms:
                assert word == (), "nonlinear in I-symbols"
                assert exp == 1, "I-symbol power > 1"
                word = iw_syms[base]
            elif base in f_syms:
                tokp = sg.ZTOK[f_syms[base]]
                for _ in range(int(exp)):
                    mono = sg.zp_mul(mono, tokp)
            else:
                assert base is sp.S.One, (fac,)
        key, scale = canon_pref(rest)
        g = groups.setdefault(key, {})
        sign = Fraction((-1) ** sg.nones(word))
        for k, q in mono.items():
            kk = (word, k)
            g[kk] = g.get(kk, Fraction(0)) + \
                Fraction(scale.p, scale.q) * sign * q
    for key, vec in groups.items():
        vec = {k: v for k, v in vec.items() if v != 0}
        if not vec:
            continue
        # split by weight (the weight grading is Galois-stable; mixed-weight
        # prefactor groups occur in reducible entries)
        byw = {}
        for (w, k), q in vec.items():
            byw.setdefault(len(w) + sg.zp_weight(k), {})[(w, k)] = q
        for wt, v2 in sorted(byw.items()):
            blocks.append({"vec": v2, "weight": wt})
    return blocks


# ---------------- feed collection ------------------------------------------

def collect_feed(include_i173):
    """returns (vectors, receipt): vectors = list of dicts with
    src, weight, vec; receipt = counts."""
    rec = {"clean": 0, "skipped": {}, "sources": []}
    vecs = []

    def eat(txt, src):
        try:
            blocks = parse_entry(txt)
        except SkipEntry as e:
            rec["skipped"][str(e)] = rec["skipped"].get(str(e), 0) + 1
            return
        rec["clean"] += 1
        for b in blocks:
            vecs.append({"src": src, "weight": b["weight"], "vec": b["vec"]})

    ents = load_entries(f"{_anc()}/4Loop_integral_results_by_hyperlog.m")
    assert len(ents) == 412
    for i, e in enumerate(ents, start=1):
        if i in TARGETS:
            rec["skipped"]["target"] = rec["skipped"].get("target", 0) + 1
            continue
        eat(e, f"hyperlog4L#{i}")
    eat(open(f"{_anc()}/fourloopI120.txt").read(), "fourloopI120")
    if include_i173:
        eat(open(f"{_anc()}/fourloopI173.txt").read(), "fourloopI173")
    for i, e in enumerate(load_entries(f"{_anc()}/3Loop_integral_insvmpl.m"),
                          start=1):
        eat(e, f"threeloop#{i}")
    rec["n_block_vectors"] = len(vecs)
    return vecs, rec


# ---------------- Galois closure + fibers ----------------------------------

def vkey(vec):
    return tuple(sorted((w, k, str(q)) for (w, k), q in vec.items()))


def saturate(vecs):
    """input [{'weight': w, 'vec': v}, ...]; returns dict
    (weight,parity) -> list of parity-pure vectors (Galois closure under
    D3/D5/D7, conjugation-split)."""
    seen = set()
    fibers = {}
    queue = []

    def push(w, v, src):
        if not v:
            return
        cv = sg.conj_vector(v)
        for par, comp in (("even", {k: (v.get(k, Fraction(0)) +
                                        cv.get(k, Fraction(0))) / 2
                                    for k in set(v) | set(cv)}),
                          ("odd", {k: (v.get(k, Fraction(0)) -
                                       cv.get(k, Fraction(0))) / 2
                                   for k in set(v) | set(cv)})):
            comp = {k: q for k, q in comp.items() if q != 0}
            if not comp:
                continue
            key = (w, par, vkey(comp))
            if key in seen:
                continue
            seen.add(key)
            fibers.setdefault((w, par), []).append(comp)
            queue.append((w, comp))

    for it in vecs:
        push(it["weight"], it["vec"], it["src"])
    while queue:
        w, v = queue.pop()
        for m in (3, 5, 7):
            if w - m < 0:
                continue
            img = sg.Dm_vector(v, w, m)
            push(w - m, img, "D")
    return fibers


# ---------------- fiber annihilators ---------------------------------------

def fiber_annihilator(fibers, wt, par):
    """exact annihilator functionals of span(V[(wt,par)]) INSIDE the parity
    subspace of the weight-wt sv coordinate space.
    Returns (functionals, dim_par_subspace, dim_fiber_span):
    functionals = list of dicts {(word,mono): Fraction} such that f.v = 0
    for all v in the fiber span AND f supported on the parity subspace
    (f applied to any vector's parity component)."""
    basis = sg.sv_basis(wt)
    # parity projector basis: for each sv basis element e, parity component
    # P e = (e + sign*conj e)/2; the parity subspace is spanned by {P e}.
    sign = Fraction(1) if par == "even" else Fraction(-1)
    # coordinates: represent subspace vectors over the full basis index
    idx = {bm: i for i, bm in enumerate(basis)}
    d = len(basis)

    def tovec(vd):
        out = [Fraction(0)] * d
        for k, q in vd.items():
            out[idx[k]] += q
        return out

    # parity-subspace basis via RREF of projected unit vectors
    proj = []
    for bm in basis:
        e = {bm: Fraction(1)}
        ce = sg.conj_vector(e)
        pv = {k: (e.get(k, Fraction(0)) + sign * ce.get(k, Fraction(0))) / 2
              for k in set(e) | set(ce)}
        pv = {k: q for k, q in pv.items() if q != 0}
        if pv:
            proj.append(tovec(pv))
    pbasis = _row_reduce(proj)          # independent rows spanning parity sub
    span = [tovec(v) for v in fibers.get((wt, par), [])]
    span_r = _row_reduce(span)
    # annihilator of span inside parity subspace: functionals f expressed in
    # the DUAL of the parity subspace; solve f.v = 0 for v in span_r, f in
    # row-space of pbasis (parity-supported functionals: the pairing uses the
    # standard coordinate dual; parity subspaces are conj-orthogonal here
    # because conj permutes the coordinate basis up to lower-triangular
    # corrections — we simply demand f in span(pbasis) and f.span = 0).
    # unknowns: coefficients a_i over pbasis rows
    nb = len(pbasis)
    rows = []
    for v in span_r:
        row = [sum(pb[j] * v[j] for j in range(d)) for pb in pbasis]
        rows.append(row)
    null = _null_space(rows, nb)
    funcs = []
    for coeffs in null:
        f = [Fraction(0)] * d
        for a, pb in zip(coeffs, pbasis):
            if a:
                for j in range(d):
                    if pb[j]:
                        f[j] += a * pb[j]
        funcs.append({basis[j]: f[j] for j in range(d) if f[j] != 0})
    return funcs, len(pbasis), len(span_r)


def _row_reduce(rows):
    """exact RREF; returns independent reduced rows (dense lists)."""
    out = []
    piv = []
    for r in rows:
        r = list(r)
        for p, pr in zip(piv, out):
            if r[p] != 0:
                f = r[p] / pr[p]
                for j in range(len(r)):
                    if pr[j]:
                        r[j] -= f * pr[j]
        nz = next((j for j, x in enumerate(r) if x != 0), None)
        if nz is None:
            continue
        piv.append(nz)
        out.append(r)
    return out


def _null_space(rows, n):
    """null space basis of the linear map Q^n -> Q^k given by rows
    (each row length n): returns list of length-n Fraction vectors."""
    red = _row_reduce(rows)
    pivs = []
    for r in red:
        pivs.append(next(j for j, x in enumerate(r) if x != 0))
    free = [j for j in range(n) if j not in pivs]
    basis = []
    for fcol in free:
        vec = [Fraction(0)] * n
        vec[fcol] = Fraction(1)
        for r, p in reversed(list(zip(red, pivs))):
            s = sum(r[j] * vec[j] for j in range(p + 1, n))
            vec[p] = -s / r[p]
        basis.append(vec)
    return basis

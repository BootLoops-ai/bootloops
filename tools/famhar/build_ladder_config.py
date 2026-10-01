#!/usr/bin/env python3
"""
build_ladder_config.py — emit families/ladder_L4.json (the 4L ladder family
config stamp) with EXACT Fraction entries throughout.

Basis/connection derivation (receipt): masters {ell^r Li_s(a), ell^r Li_s(b),
ell^r}, ell = log(a b). d ell = dlog a + dlog b; d Li_s(x) = Li_{s-1}(x) dlog x
(s>=2), d Li_1(x) = -dlog(1-x) * 1 (coefficient on the ell^r master, since
Li_0 would be rational — the s=1 row closes onto the ell^r tower through the
letter (1-x)). All entries integer. Validated by gate 0 (A.J vs finite
differences of the closed-form basis) before any transport is trusted.
"""
import json
import os
import sys
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ladder_reference import masters_spec, f_coeffs  # noqa: E402

L = 4
masters = masters_spec(L)
idx = {(m["type"], m["r"], m.get("s")): i for i, m in enumerate(masters)}
n = len(masters)

conn = {"a": {}, "one_minus_a": {}, "b": {}, "one_minus_b": {}}


def add(letter, i, j, val):
    key = f"{i},{j}"
    cur = Fraction(conn[letter].get(key, "0"))
    conn[letter][key] = str(cur + val)


for i, m in enumerate(masters):
    r = m["r"]
    if m["type"] == "lA":
        s = m["s"]
        if r >= 1:
            add("a", i, idx[("lA", r - 1, s)], r)
            add("b", i, idx[("lA", r - 1, s)], r)
        if s >= 2:
            add("a", i, idx[("lA", r, s - 1)], 1)
        else:  # s == 1: d Li_1(a) = -dlog(1-a)
            add("one_minus_a", i, idx[("l", r, None)], -1)
    elif m["type"] == "lB":
        s = m["s"]
        if r >= 1:
            add("a", i, idx[("lB", r - 1, s)], r)
            add("b", i, idx[("lB", r - 1, s)], r)
        if s >= 2:
            add("b", i, idx[("lB", r, s - 1)], 1)
        else:
            add("one_minus_b", i, idx[("l", r, None)], -1)
    else:  # ell^r
        if r >= 1:
            add("a", i, idx[("l", r - 1, None)], r)
            add("b", i, idx[("l", r - 1, None)], r)

targets = {}
for Lt in range(1, L + 1):
    cs = f_coeffs(Lt)
    coeffs = {}
    for r in range(Lt + 1):
        s = 2 * Lt - r
        coeffs[str(idx[("lA", r, s)])] = str(cs[r])
        coeffs[str(idx[("lB", r, s)])] = str(-cs[r])
    targets[f"Phi{Lt}"] = {
        "coeffs": coeffs,
        "prefactor": "((a-1)*(b-1))/(a-b)",
        "dprefactor": {"a": "-((b-1)**2)/((a-b)**2)",
                       "b": "((a-1)**2)/((a-b)**2)"},
        "doc": f"Phi^({Lt})(u,v) of 1303.6909 eq (2.2)/(2.3); "
               "= -f/(z-zbar) with chart a=z/(z-1), b=zbar/(zbar-1)."}

cfg = {
    "name": "ladder_L4",
    "doc": ("4L ladder calibration family (Usyukina-Davydychev Phi^(1..4)). "
            "Chart (a,b)=(w,wbar) independent; Euclidean slice b=conj(a). "
            "u=a*b/((a-1)*(b-1)), v=1/((a-1)*(b-1)). Branch loci in chart = "
            "letter zeros {a=0, a=1, b=0, b=1}; kinematic-space loci u=0 "
            "(a=0 or b=0 sheet-locus), v=0 (chart infinity), lambda=0 (a=b; "
            "NOT a DE singularity — prefactor pole only, cancels in Phi). "
            "ORIENTATION LAW (frame stamp, receipt gate_ladder.py matches): "
            "frame ladder integral = -(1/u) * Phi^(L)(1/u, v/u)."),
    "chart_vars": ["a", "b"],
    "eps0": "0",
    "eps_doc": "family is eps-independent; harness threads eps0 through detransport (fixed-eps marches; real-eps0 Frobenius law inherited).",
    "n_masters": n,
    "masters": masters,
    "letters": [
        {"name": "a", "const": "0", "a": "1", "b": "0"},
        {"name": "one_minus_a", "const": "1", "a": "-1", "b": "0"},
        {"name": "b", "const": "0", "a": "0", "b": "1"},
        {"name": "one_minus_b", "const": "1", "a": "0", "b": "-1"},
    ],
    "connection": conn,
    "anchor": {
        "point": {"a": "1/8", "b": "1/9"},
        "values_plugin": "ladder_reference.anchor_values",
        "doc": "anchor J = principal-branch basis values at a small real chart point (calibration family: closed basis functions; production families: AMFlow or stored anchor vector)."},
    "kinematics": {
        "vars": ["u", "v"],
        "map": {"u": "a*b/((a-1)*(b-1))", "v": "1/((a-1)*(b-1))"},
        "jacobian": {"du_da": "-b/(((a-1)**2)*(b-1))",
                     "du_db": "-a/((a-1)*((b-1)**2))",
                     "dv_da": "-1/(((a-1)**2)*(b-1))",
                     "dv_db": "-1/((a-1)*((b-1)**2))"}},
    "targets": targets,
    "clearance_min": "1/50",
    "calibration_points": {
        "PA": {"z": "2/5+3/5i", "a": "1/6-5/6i", "b": "1/6+5/6i", "u": "13/25", "v": "18/25"},
        "PB": {"z": "1/4+1/2i", "a": "1/13-8/13i", "b": "1/13+8/13i", "u": "5/16", "v": "13/16"},
        "refs": "phi_values.json: the 100-digit reference values at the calibration points (two-rep crosschecked), a record of the calibration run not shipped with this package"},
}

out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "families", "ladder_L4.json")
with open(out, "w") as f:
    json.dump(cfg, f, indent=1)
nnz = sum(len(v) for v in conn.values())
print(f"wrote {out}: n={n} masters, {nnz} nonzero connection entries, {len(targets)} targets")

#!/usr/bin/env python3
"""X(3,6) Grassmannian string integral pilot -- symbolic setup.

Positive parametrization: Gimenez Umbert--Sturmfels arXiv:2501.10805, eq.(eq:36):
    X = [[0, 0, 1, 1, 1,        1                    ],
         [0,-1, 0, 1, 1+x11,    1+x11*(1+x12)        ],
         [1, 0, 0, 1, 1+x11+x21, m33                 ]]
    m33 = 1 + x11*(1+x12) + x21*(1+x12+x22)
Integral (AHL 1912.08707, sec. "Grassmannian string integrals"):
    I_{3,6}(s;a') = (a')^4 int_{R+^4} prod dx/x  prod_{a<b<c} X_abc(x)^{a' s_abc}
Momentum conservation: for each a: sum_{b<c, b,c != a} s_abc = 0.

Outputs (JSON): minors, supports, conservation solve, split point, generic point.
"""
import sympy as sp
import json, itertools, sys

x11, x12, x21, x22 = sp.symbols('x11 x12 x21 x22', positive=True)
XV = [x11, x12, x21, x22]

m33 = 1 + x11*(1+x12) + x21*(1+x12+x22)
M = sp.Matrix([
    [0,  0, 1, 1, 1,          1],
    [0, -1, 0, 1, 1+x11,      1+x11*(1+x12)],
    [1,  0, 0, 1, 1+x11+x21,  m33]])

triples = list(itertools.combinations(range(1, 7), 3))  # 20 triples, 1-indexed
minors = {}
for t in triples:
    d = sp.expand(M[:, [t[0]-1, t[1]-1, t[2]-1]].det())
    minors[t] = d

# ---- checks: positivity of coefficients, classification -------------------
report = {}
bad = []
for t, d in minors.items():
    poly = sp.Poly(d, *XV)
    coeffs = poly.coeffs()
    if any(c <= 0 for c in coeffs):
        bad.append((t, str(d)))
report['all_minors_positive_coeffs'] = (len(bad) == 0)
report['bad'] = [str(b) for b in bad]

const_minors, mono_minors, poly_minors = [], [], []
for t, d in minors.items():
    if d.is_number:
        const_minors.append(t)
    elif sp.Poly(d, *XV).is_monomial:
        mono_minors.append(t)
    else:
        poly_minors.append(t)
report['const_minors'] = [str(t) for t in const_minors]
report['mono_minors'] = [(str(t), str(minors[t])) for t in mono_minors]
report['poly_minors'] = [(str(t), str(minors[t])) for t in poly_minors]
report['counts'] = dict(const=len(const_minors), mono=len(mono_minors),
                        poly=len(poly_minors))

# variable support of each minor
support = {t: sorted(str(v) for v in minors[t].free_symbols) for t in triples}
report['support'] = {str(t): support[t] for t in triples}

# ---- GS eq:recovery cross-ratio check --------------------------------------
def X(*t): return minors[tuple(t)]
rec = {
    'x11': sp.simplify(X(1,2,3)*X(1,4,5)/(X(1,2,5)*X(1,3,4)) - x11),
    'x12': sp.simplify(X(1,2,4)*X(1,5,6)/(X(1,2,6)*X(1,4,5)) - x12),
    'x21': sp.simplify(X(1,2,3)*X(1,2,4)*X(3,4,5)/(X(1,2,5)*X(1,3,4)*X(2,3,4)) - x21),
    'x22': sp.simplify(X(1,2,5)*X(1,3,4)*X(4,5,6)/(X(1,2,6)*X(1,4,5)*X(3,4,5)) - x22),
}
report['recovery_check_all_zero'] = all(v == 0 for v in rec.values())

# ---- GS split-kinematics admissible minors check ---------------------------
gs_split = {(1,2,3):1, (1,2,4):1, (1,2,5):1, (1,2,6):1, (1,3,4):1, (2,3,4):1,
            (1,4,5):x11, (1,5,6):x11*x12, (2,3,5):1+x11+x21, (3,4,5):x21,
            (3,4,6):x21*(1+x12+x22), (4,5,6):x11*x21*x22}
report['gs_split_minors_match'] = all(
    sp.simplify(minors[t]-e) == 0 for t, e in gs_split.items())

# ---- conservation: solve 6 constant-minor exponents from the other 14 ------
svars = {t: sp.Symbol('s_%d%d%d' % t) for t in triples}
cons = [sp.Add(*[svars[t] for t in triples if a in t]) for a in range(1, 7)]
sol = sp.solve(cons, [svars[t] for t in const_minors], dict=True)
report['conservation_solvable_for_const_exponents'] = (len(sol) == 1)
report['conservation_solution'] = {str(k): str(v) for k, v in sol[0].items()}

with open('setup_report.json', 'w') as f:
    json.dump(report, f, indent=1, default=str)

print(json.dumps({k: report[k] for k in
      ['all_minors_positive_coeffs', 'counts', 'recovery_check_all_zero',
       'gs_split_minors_match', 'conservation_solvable_for_const_exponents']},
      indent=1))
print("const:", report['const_minors'])
print("mono :", report['mono_minors'])
for t, d in report['poly_minors']:
    print("poly :", t, "=", d, "   support:", report['support'][t])

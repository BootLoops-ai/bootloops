#!/usr/bin/env python3
"""ibp_gen.py — independent IBP identity generator (fresh-identity
source; independent of kira's Generate AND of the Route-B solver, so the
coverage gate tests real IBP validity, not agreement with the staging tool).

Scope: L-loop families with E external momenta whose external scalar
products are FIXED by kinematics scalarproduct_rules (the whole vac3 corpus:
pure vacuum plus the AMFlow _sm0 chain families where a decoupled loop
became an external leg). Masses and rule values are linear in eta.

With unknown scalar products sp (>= one loop momentum) and props
P_j = q_j^2 - m_j = sum_s A_{js} sp_s + ee_j - m_j:

    sp = B (P + shift),  B = A^{-1},  shift_j = m_j - ee_j

IBP identity for seed a, operator (d/dk_i) . v  (i loop, v any momentum):

  0 = d*delta_{iv} I(a) - sum_j a_j u_j[i] * [
        sum_m C_{jvm} ( I(a+e_j-e_m) + shift_m I(a+e_j) )
        + K_{jv} I(a+e_j) ]

  C_{jvm} = sum_c 2 u_j[c] B[slot(c,v), m]   (over unknown slots (c,v))
  K_{jv}  = sum_{c: (c,v) both external} 2 u_j[c] rule(c,v)
"""
import itertools
import re
from fractions import Fraction


def _lin_eta(s):
    """eta-linear expression -> (c0, ceta) Fractions. Fast token path for
    the vac3-style forms; exact Fraction-eval fallback for general rational
    forms like '(-1)/(2)' (k2disp kinematics). Linearity asserted."""
    t = str(s).replace(" ", "")
    if not t:
        return Fraction(0), Fraction(0)
    toks = re.findall(r"([+-]?)(\d+\*?eta|eta|\d+)", t)
    if "".join(sg + v for sg, v in toks) == t:
        c0 = ceta = Fraction(0)
        for sg, tok in toks:
            f = Fraction(-1 if sg == "-" else 1)
            if "eta" in tok:
                ceta += f * Fraction(
                    tok.replace("*eta", "").replace("eta", "") or 1)
            else:
                c0 += f * Fraction(tok)
        return c0, ceta
    # general exact fallback: evaluate at eta = 0, 1, 2 with Fractions
    assert re.match(r"^[0-9eta+\-*/^()]+$", t), f"unparsed eta-linear: {s!r}"
    expr = re.sub(r"(\d+)", r"Fraction(\1)", t.replace("^", "**"))

    def ev(e0):
        return Fraction(eval(  # noqa: S307 (restricted charset, no builtins)
            expr, {"eta": Fraction(e0), "Fraction": Fraction,
                   "__builtins__": {}}))

    f0, f1, f2 = ev(0), ev(1), ev(2)
    assert f2 - f1 == f1 - f0, f"non-eta-linear kinematic value: {s!r}"
    return f0, f1 - f0


def parse_family_spec(fam_yaml, kin_yaml=""):
    """-> dict(loops, exts, mom, u[j], mass[j]=(c0,ceta), rules{(a,b):(c0,ceta)})."""
    loops = [x.strip() for x in re.search(
        r"loop_momenta:\s*\[([\w, ]+)\]", fam_yaml).group(1).split(",")]
    raw = []
    symbols = []
    for mom, _zero in re.findall(r'-\s*\[\s*"([^"]+)"\s*,\s*(\S+?)\s*\]',
                                 fam_yaml):
        m = re.match(r"^\(?([^()]*)\)?\s*\^\s*2\s*(.*)$", mom.strip())
        assert m, f"unparsed propagator: {mom!r}"
        core, tail = m.group(1).strip(), m.group(2).strip()
        terms = re.findall(r"([+-]?)\s*(\w+)", core)
        for _sg, var in terms:
            if var not in symbols:
                symbols.append(var)
        c0, ceta = _lin_eta(tail)
        raw.append((terms, (-c0, -ceta)))        # P = q^2 + tail => m = -tail
    exts = [v for v in symbols if v not in loops]
    mom = loops + exts
    u = []
    masses = []
    for terms, mm in raw:
        vec = [0] * len(mom)
        for sg, var in terms:
            vec[mom.index(var)] += -1 if sg == "-" else 1
        u.append(tuple(vec))
        masses.append(mm)
    rules = {}
    for a, b, val in re.findall(
            r"-\s*\[\[\s*(\w+)\s*,\s*(\w+)\s*\]\s*,\s*([^\]]+?)\s*\]",
            kin_yaml or ""):
        rules[(a, b)] = rules[(b, a)] = _lin_eta(val)
    return {"loops": loops, "exts": exts, "mom": mom, "u": u,
            "masses": masses, "rules": rules}


def build_maps(spec):
    """Invert the unknown-sp <-> prop map. Returns (slots, B, shifts, Kconst)
    with shifts_j = m_j - ee_j (eta-linear pairs) and Kconst[j][v] eta-linear."""
    mom, loops, exts = spec["mom"], spec["loops"], spec["exts"]
    nl, nm = len(loops), len(mom)

    def rule(c, v):
        key = (mom[c], mom[v])
        assert key in spec["rules"], f"missing scalarproduct_rule for {key}"
        return spec["rules"][key]

    slots = [(c, v) for c in range(nm) for v in range(c, nm) if c < nl]
    slot_of = {s: k for k, s in enumerate(slots)}
    n = len(slots)
    assert len(spec["u"]) == n, \
        f"family has ISPs or over-complete props: {len(spec['u'])} props vs {n} sps"
    A, shifts = [], []
    for uj, (m0, meta) in zip(spec["u"], spec["masses"]):
        row = [Fraction(0)] * n
        ee0, eeeta = Fraction(0), Fraction(0)
        for ci in range(nm):
            if not uj[ci]:
                continue
            for di in range(nm):
                if not uj[di]:
                    continue
                a_, b_ = min(ci, di), max(ci, di)
                if a_ < nl:
                    row[slot_of[(a_, b_)]] += Fraction(uj[ci] * uj[di])
                else:
                    r0, reta = rule(a_, b_)
                    ee0 += uj[ci] * uj[di] * r0
                    eeeta += uj[ci] * uj[di] * reta
        A.append(row)
        shifts.append((m0 - ee0, meta - eeeta))
    # invert A over Q
    M = [list(A[i]) + [Fraction(int(i == j)) for j in range(n)]
         for i in range(n)]
    for col in range(n):
        piv = next(r for r in range(col, n) if M[r][col] != 0)
        M[col], M[piv] = M[piv], M[col]
        pv = M[col][col]
        M[col] = [x / pv for x in M[col]]
        for r in range(n):
            if r != col and M[r][col]:
                f = M[r][col]
                M[r] = [x - f * y for x, y in zip(M[r], M[col])]
    B = [row[n:] for row in M]
    # constant (external-external) parts of 2 q_j . v
    K = [[(Fraction(0), Fraction(0)) for _ in range(nm)]
         for _ in range(len(spec["u"]))]
    for j, uj in enumerate(spec["u"]):
        for v in range(nl, nm):                 # v external
            k0 = keta = Fraction(0)
            for c in range(nl, nm):             # c external too
                if uj[c]:
                    r0, reta = rule(c, v)
                    k0 += 2 * uj[c] * r0
                    keta += 2 * uj[c] * reta
            K[j][v] = (k0, keta)
    return slots, slot_of, B, shifts, K


def identities(spec, seeds, p, d0, eta0):
    """Yield (seed, i, v, {index_tuple: coeff mod p}) for all L(L+E) ops."""
    mom, loops = spec["mom"], spec["loops"]
    nl, nm = len(loops), len(mom)
    nprop = len(spec["u"])
    slots, slot_of, B, shifts, K = build_maps(spec)

    def val(fr):
        return fr.numerator % p * pow(fr.denominator % p, p - 2, p) % p

    def lin(pair):
        return (val(pair[0]) + val(pair[1]) * eta0) % p

    shift_val = [lin(s) for s in shifts]
    K_val = [[lin(K[j][v]) for v in range(nm)] for j in range(nprop)]
    # C[j][v][m] over unknown slots
    C = [[[0] * nprop for _ in range(nm)] for _ in range(nprop)]
    for j, uj in enumerate(spec["u"]):
        for v in range(nm):
            for c in range(nm):
                if not uj[c]:
                    continue
                a_, b_ = min(c, v), max(c, v)
                if a_ >= nl:                    # both external -> in K
                    continue
                srow = B[slot_of[(a_, b_)]]
                for m in range(nprop):
                    if srow[m]:
                        C[j][v][m] = (C[j][v][m] + 2 * uj[c] * val(srow[m])) % p

    for a in seeds:
        for i in range(nl):
            for v in range(nm):
                row = {}

                def add(t, c):
                    if c % p:
                        row[t] = (row.get(t, 0) + c) % p
                        if not row[t]:
                            del row[t]

                if i == v:
                    add(a, d0)
                for j in range(nprop):
                    if not a[j] or not spec["u"][j][i]:
                        continue
                    pref = (-a[j] * spec["u"][j][i]) % p
                    for m in range(nprop):
                        cf = C[j][v][m]
                        if not cf:
                            continue
                        up = list(a)
                        up[j] += 1
                        dn = list(up)
                        dn[m] -= 1
                        add(tuple(dn), pref * cf)
                        if shift_val[m]:
                            add(tuple(up), pref * cf * shift_val[m])
                    if K_val[j][v]:
                        up = list(a)
                        up[j] += 1
                        add(tuple(up), pref * K_val[j][v])
                if row:
                    yield (a, i, v, row)


def seed_box(sectors_nontrivial, n_idx, rmax, smax, dmax):
    """Enumerate seed tuples with all IBP images inside the (rmax,smax,dmax)
    box: positive part dots<=dmax-1 and r<=rmax-1, negative part s<=smax-1."""
    seeds = []
    for sec in sectors_nontrivial:
        pos = [i for i in range(n_idx) if sec >> i & 1]
        neg = [i for i in range(n_idx) if not sec >> i & 1]
        t = len(pos)
        for dots in range(0, max(0, min(dmax - 1, rmax - 1 - t)) + 1):
            for dist in itertools.combinations_with_replacement(pos, dots):
                for stot in range(0, smax):     # s <= smax-1
                    for sdist in itertools.combinations_with_replacement(
                            neg, stot):
                        a = [0] * n_idx
                        for i in pos:
                            a[i] = 1
                        for i in dist:
                            a[i] += 1
                        for i in sdist:
                            a[i] -= 1
                        seeds.append(tuple(a))
    return sorted(set(seeds))

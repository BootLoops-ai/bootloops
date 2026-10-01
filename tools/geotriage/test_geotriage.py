#!/usr/bin/env python3
"""
Regression tests for geotriage.py on the 5 KNOWN cases.

Run:  python test_geotriage.py
"""
import os, sys, json
from fractions import Fraction
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from geotriage import classify, elliptic_invariants, cm_check_elliptic, \
    _binary_quartic_invariants, CM_J_TABLE, _ap_from_j, _supersingular_pattern, \
    k3_ap_probe, _toric_banana_count, _f3_level15_ap, _primes_upto  # noqa: E402

PASS = []; FAIL = []


def _check(label, cond, detail=""):
    if cond:
        PASS.append(label)
        print(f"  [PASS] {label}  {detail}")
    else:
        FAIL.append(label)
        print(f"  [FAIL] {label}  {detail}")


# ---------------------------------------------------------------------------
def test_box1l():
    print("\n=== TEST 1: 1-loop box ===")
    rep = classify(os.path.join(HERE, "graph_specs", "box1l.json"))
    _check("box1l type=polylog", rep['variety']['type'] == "polylog",
           f"got {rep['variety']['type']}")
    _check("box1l genus=0", rep['variety']['genus'] == 0)
    _check("box1l route=A", rep['route'] == 'A', f"got {rep['route']}")
    _check("box1l genus_method=exact",
           rep['honesty']['genus_method'] == 'exact')
    return rep


# ---------------------------------------------------------------------------
def test_sunrise_eq():
    print("\n=== TEST 2: equal-mass sunrise ===")
    rep = classify(os.path.join(HERE, "graph_specs", "sunrise_eq.json"))
    _check("sunrise type=elliptic", rep['variety']['type'] == "elliptic",
           f"got {rep['variety']['type']}")
    _check("sunrise genus=1", rep['variety']['genus'] == 1)
    ell = rep.get('elliptic', {})
    sig = tuple(ell.get('fiber_signature', []))
    _check("sunrise fiber config (1,2,3,6)", sig == (1, 2, 3, 6),
           f"got {sig}")
    _check("sunrise on X1(6)", ell.get('level_N') == 6 and ell.get('on_X1N'),
           f"level_N={ell.get('level_N')}, on_X1N={ell.get('on_X1N')}, "
           f"group={ell.get('congruence_group')}")
    # The equal-mass sunrise FAMILY is congruence-modular (route B); the
    # D=-3 in the spec is the ring conductor (L(χ_{-3},2)), not a generic-fiber
    # CM-j hit.  CM-j only fires at SPECIAL interior points (e.g. y=±3 of the
    # banana map), which are cusps in the t-variable — so cm.is_cm is None here.
    _check("sunrise ring conductor D=-3",
           any("D=-3" in str(c) or "chi_{-3}" in str(c)
               for c in rep.get('ring_conductors', [])),
           f"ring_conductors={rep.get('ring_conductors')}")
    _check("sunrise route=B", rep['route'] == 'B',
           f"got {rep['route']}: {rep['reason']}")
    return rep


# ---------------------------------------------------------------------------
def test_icc():
    print("\n=== TEST 3: ice-cream-cone ===")
    rep = classify(os.path.join(HERE, "graph_specs", "icc.json"))
    _check("icc type=elliptic", rep['variety']['type'] == "elliptic",
           f"got {rep['variety']['type']}")
    _check("icc genus=1", rep['variety']['genus'] == 1)
    ell = rep.get('elliptic', {})
    sig = tuple(ell.get('fiber_signature', []))
    on_small = ell.get('on_X1N')
    _check("icc NOT on small X1(N)", not on_small,
           f"level_N={ell.get('level_N')}, sig={sig}, sum={ell.get('sum_fiber_orders')}")
    cm = rep.get('cm', {})
    _check("icc non-CM at generic fiber", cm.get('is_cm') is not True,
           f"is_cm={cm.get('is_cm')}, ev={cm.get('evidence')}")
    _check("icc route=C", rep['route'] == 'C',
           f"got {rep['route']}: {rep['reason']}")
    return rep


# ---------------------------------------------------------------------------
def test_banana3l_eq():
    print("\n=== TEST 4: equal-mass 3-loop banana (K3) ===")
    rep = classify(os.path.join(HERE, "graph_specs", "banana3l_eq.json"))
    _check("banana3l type=K3", rep['variety']['type'] == "K3",
           f"got {rep['variety']['type']}")
    _check("banana3l pf_order=3", rep['variety']['pf_order'] == 3,
           f"got {rep['variety']['pf_order']}")
    cm = rep.get('cm', {})
    _check("banana3l CM by Q(√-15) (Livné)",
           cm.get('is_cm') is True and cm.get('discriminant_D') == -15,
           f"is_cm={cm.get('is_cm')}, D={cm.get('discriminant_D')}, "
           f"zero_density={cm.get('zero_density')}")
    _check("banana3l Sym² root recorded", 'sym2_root' in rep,
           f"{rep.get('sym2_root')}")
    _check("banana3l route=B", rep['route'] == 'B',
           f"got {rep['route']}: {rep['reason']}")
    # spot-check the a_p evidence tuples
    ev = cm.get('evidence', [])
    if ev and isinstance(ev[0], tuple) and len(ev[0]) == 3:
        bad = [(p, a, c) for (p, a, c) in ev if (a == 0) != (c == -1)]
        _check("banana3l a_p=0 ⇔ χ_{-15}(p)=-1 (all primes)", len(bad) == 0,
               f"mismatches: {bad[:5]}")
    return rep


# ---------------------------------------------------------------------------
def test_crossedbox3l():
    print("\n=== TEST 5: 3-loop crossed box (light-by-light) ===")
    rep = classify(os.path.join(HERE, "graph_specs", "crossedbox3l.json"))
    _check("crossedbox3l type=polylog", rep['variety']['type'] == "polylog",
           f"got {rep['variety']['type']}")
    _check("crossedbox3l genus=0", rep['variety']['genus'] == 0)
    _check("crossedbox3l route=A", rep['route'] == 'A', f"got {rep['route']}")
    # honesty: spec-supplied caveat MUST be present
    cav = rep['honesty']['caveats']
    _check("crossedbox3l caveat: maxcut spec-supplied",
           any("spec" in c for c in cav), f"caveats={cav}")
    return rep


# ---------------------------------------------------------------------------
def test_unit_ap_counter():
    """UNIT: wired a_p point counter — CM supersingular patterns (twist-invariant)."""
    print("\n=== UNIT: wired a_p counter (supersingular patterns) ===")
    for jv, D, label in [(1728, -4, "j=1728"), (0, -3, "j=0")]:
        aps = {}
        for p in _primes_upto(60):
            a = _ap_from_j(Fraction(jv), p)
            if a is not None:
                aps[p] = a
        trip, ok = _supersingular_pattern(aps, D)
        _check(f"{label}: a_p=0 ⟺ χ_{{{D}}}(p)=−1 (all primes ≤ 60)",
               ok and len(trip) >= 10, f"n={len(trip)}, first={trip[:4]}")
        _check(f"{label}: Hasse bound a_p² ≤ 4p",
               all(a*a <= 4*p for p, a in aps.items()))
        if jv == 1728:
            # wrong-discriminant control: the SAME a_p against χ_{-3} must fail
            _, okw = _supersingular_pattern(aps, -3)
            _check("control: j=1728 pattern vs wrong D=-3 FAILS", not okw)


# ---------------------------------------------------------------------------
def test_cm_elliptic_wired():
    """UNIT: cm_check_elliptic emits the promised (p, a_p, χ_D(p)) evidence."""
    print("\n=== UNIT: cm_check_elliptic wired evidence ===")
    r = cm_check_elliptic('1728')
    trips = [e for e in r['evidence'] if isinstance(e, tuple) and len(e) == 3]
    _check("j=1728 → CM, D=-4",
           r['is_cm'] is True and r['discriminant_D'] == -4,
           f"is_cm={r['is_cm']}, D={r['discriminant_D']}")
    _check("j=1728 certified by the a_p pattern", r.get('certified') is True)
    _check("j=1728 evidence carries (p,a_p,χ_D) triples",
           len(trips) >= 10, f"n={len(trips)}")
    _check("j=1728 evidence pattern consistent",
           all((a == 0) == (c == -1) for _, a, c in trips))
    r2 = cm_check_elliptic('1729')
    _check("j=1729 → non-CM (table complete for rational j)",
           r2['is_cm'] is False and r2.get('certified') is True,
           f"is_cm={r2['is_cm']}, certified={r2.get('certified')}")
    _check("j=1729 a_p scan matches no χ_D",
           any('no CM discriminant matches' in str(e) for e in r2['evidence']),
           f"ev={r2['evidence'][:2]}")


# ---------------------------------------------------------------------------
def test_unit_k3_probe():
    """UNIT: K3 exact background centring + Livné scan on synthetic truth."""
    print("\n=== UNIT: K3 background centring + Livné scan (synthetic truth) ===")
    primes = [p for p in _primes_upto(60) if p >= 5]
    ap = _f3_level15_ap(60)
    # planted truth: background p²+3p+1 + the level-15 weight-3 CM a_p
    bank = {p: p*p + 3*p + 1 + ap[p] for p in primes}
    pr = k3_ap_probe(bank)
    _check("synthetic CM bank locks", pr['locked'] is True,
           f"support={pr['support']}/{pr['n_primes']}")
    _check("planted background (1,3,1) recovered exactly",
           pr['background'] == (1, 3, 1), f"got {pr['background']}")
    _check("unique Livné match D=-15", pr['matches'] == [-15],
           f"got {pr['matches']}")
    _check("zero density ≈ 1/2", 0.3 <= pr['zero_density'] <= 0.7,
           f"got {pr['zero_density']:.2f}")
    # mutation control: corrupt ONE inert-prime count → the match must die
    mut = dict(bank); mut[29] += 1
    prm = k3_ap_probe(mut)
    _check("mutation control: one flipped count kills the D=-15 match",
           -15 not in prm['matches'], f"matches={prm['matches']}")
    # degenerate control: equal-mass t=0 count is purely polynomial in p
    raw0 = {p: _toric_banana_count([1, 1, 1, 1], 0, p) for p in primes}
    pr0 = k3_ap_probe(raw0)
    _check("degenerate fiber: all-zero residuals flagged, no D claimed",
           pr0['locked'] and pr0['all_zero'] and pr0['matches'] == [],
           f"bg={pr0['background']}")
    # non-CM-shaped control: scattered nonzero residuals → probe refuses
    vals = [1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 6, -6, 7, -7, 8]
    rawn = {p: p*p + 3*p + 1 + v for p, v in zip(primes, vals)}
    prn = k3_ap_probe(rawn)
    _check("scattered residuals: probe refuses to lock", prn['locked'] is False,
           f"support={prn['support']}")


# ---------------------------------------------------------------------------
def test_k3_generic_probe():
    """END-TO-END: unequal-mass 4-edge banana through the generic K3 probe."""
    print("\n=== TEST 6: unequal-mass 4-edge banana (generic K3 probe) ===")
    spec = {"name": "banana3l_1112",
            "description": "4-edge banana, model coefficients (1,1,1,2), t=1 — "
                           "exercises the GENERIC K3 probe path (not the "
                           "equal-mass eta-quotient shortcut).",
            "edges": [[1, 2, "1"], [1, 2, "1"], [1, 2, "1"], [1, 2, "2"]],
            "nodes": {"1": [1], "2": [2]},
            "kinematics": {"invariants": ["t"],
                           "external_masses": {"1": "t", "2": "t"},
                           "channels": {"1": "t"}},
            "fiber_point": {"t": 1}}
    rep = classify(spec)
    _check("banana1112 type=K3", rep['variety']['type'] == 'K3',
           f"got {rep['variety']['type']}")
    cm = rep['cm']
    _check("banana1112 probe verdict CM, D=-3",
           cm['is_cm'] is True and cm['discriminant_D'] == -3,
           f"is_cm={cm['is_cm']}, D={cm['discriminant_D']}, "
           f"bg={cm.get('background')}, support={cm.get('support')}")
    _check("banana1112 verdict flagged NOT certified (probe)",
           cm.get('certified') is False)
    _check("banana1112 route=B with probe caveat",
           rep['route'] == 'B' and
           any('NOT theorem-certified' in c for c in rep['honesty']['caveats']),
           f"route={rep['route']}")
    ev = cm['evidence']
    _check("banana1112 evidence: a_p=0 ⟺ χ_{-3}(p)=−1 on every prime",
           len(ev) >= 10 and all((a == 0) == (chi == -1) for (p, a, chi) in ev),
           f"n={len(ev)}")
    # independent cross-check (computed here, nothing hand-typed): on split
    # primes the residual a_p must equal the weight-3 Hecke trace 2(a²−3b²)
    # for the decomposition p = a² + 3b².
    split = [(p, a) for (p, a, chi) in ev if chi == 1]
    ok = len(split) >= 5
    for p, apv in split:
        dec = []
        b = 1
        while 3*b*b < p:
            x2 = p - 3*b*b
            x = int(x2**0.5)
            if x*x == x2:
                dec.append((x, b))
            b += 1
        ok = ok and bool(dec) and all(2*(x*x - 3*b*b) == apv for (x, b) in dec)
    _check("banana1112 split-prime a_p = 2(a²−3b²) with p=a²+3b² (Hecke trace)",
           ok, f"split={split}")
    return rep


# ---------------------------------------------------------------------------
def test_unit_jinvariant():
    """Unit: binary-quartic j-invariant pipeline reproduces j=1728 on y²=x⁴−1."""
    print("\n=== UNIT: j-invariant pipeline ===")
    x = sp.Symbol('x')
    I, J, D, j = _binary_quartic_invariants(x**4 - 1, x)
    _check("j(y²=x⁴−1) = 1728", sp.simplify(j - 1728) == 0, f"got j={j}")
    # y² = x³ - 1 (deg-3 → a=0): j = 0
    I, J, D, j = _binary_quartic_invariants(x**3 - 1, x)
    _check("j(y²=x³−1) = 0", sp.simplify(j) == 0, f"got j={j}")


# ---------------------------------------------------------------------------
def test_member_exactj():
    """MEMBER: exactj — exact j-minpoly for algebraic fibers, two-route COPAIR."""
    print("\n=== MEMBER: exactj (exact j-minpoly, two flint routes COPAIRed) ===")
    try:
        from exactj import exact_j_resultant
    except ImportError as e:
        print(f"  [SKIP] exactj legs skipped by name: python-flint unavailable ({e})")
        return
    # x²+1 (x0 = i, lam = −1): j = 256·3³/4 = 1728, CM table D = −4
    r = exact_j_resultant([1, 0, 1])
    _check("exactj x²+1 → rational j = 1728",
           r['j_minpoly_deg'] == 1 and r.get('j_exact') == '1728',
           f"got deg={r['j_minpoly_deg']}, j={r.get('j_exact')}")
    _check("exactj x²+1 in CM table with D=-4",
           r.get('j_in_CM_table') is True and r.get('D_from_table') == -4,
           f"in_table={r.get('j_in_CM_table')}, D={r.get('D_from_table')}")
    _check("exactj x²+1 all engine gates green",
           all(r['engine_gates'].values()), f"gates={r['engine_gates']}")
    # x⁴−x²+1 (x0 a primitive 12th root of unity, lam a primitive 6th root):
    # lam²−lam+1 = 0 → j = 0, CM table D = −3
    r0 = exact_j_resultant([1, 0, -1, 0, 1])
    _check("exactj x⁴−x²+1 → rational j = 0",
           r0['j_minpoly_deg'] == 1 and r0.get('j_exact') == '0',
           f"got deg={r0['j_minpoly_deg']}, j={r0.get('j_exact')}")
    _check("exactj x⁴−x²+1 in CM table with D=-3",
           r0.get('j_in_CM_table') is True and r0.get('D_from_table') == -3,
           f"in_table={r0.get('j_in_CM_table')}, D={r0.get('D_from_table')}")
    _check("exactj x⁴−x²+1 all engine gates green", all(r0['engine_gates'].values()))
    # deg-2 cyclotomic x²+x+1 (lam = ω a primitive cube root of unity):
    # independent oracle computed HERE in Q(ω) with pure Fraction arithmetic
    # (ω² = −1−ω), no flint anywhere — a genuine cross-route COPAIR.
    def _mul(u, v):          # (a+bω)(c+dω) with ω² = −1−ω
        a, b = u; c, d = v
        return (a*c - b*d, a*d + b*c - b*d)
    lam = (Fraction(0), Fraction(1))
    l2 = _mul(lam, lam)
    nbase = (l2[0] - lam[0] + 1, l2[1] - lam[1])            # lam²−lam+1
    num = _mul(_mul(nbase, nbase), nbase)                   # (lam²−lam+1)³
    lm1 = (lam[0] - 1, lam[1])
    den = _mul(l2, _mul(lm1, lm1))                          # lam²(lam−1)²
    _check("oracle sanity: N and D land in Q (ω-parts vanish)",
           num[1] == 0 and den[1] == 0, f"num={num}, den={den}")
    j_oracle = 256 * num[0] / den[0]
    r3 = exact_j_resultant([1, 1, 1])
    _check("exactj x²+x+1 COPAIR vs independent Q(ω) oracle",
           r3['j_minpoly_deg'] == 1 and r3.get('j_exact') == str(j_oracle),
           f"engine j={r3.get('j_exact')}, oracle j={j_oracle}")
    _check("exactj x²+x+1 route-A/route-B agreement gate green",
           r3['engine_gates']['copair_routes_agree'] and
           r3['engine_gates']['annihilating_factor_unique'],
           f"gates={r3['engine_gates']}")
    _check("exactj x²+x+1 non-CM certificate fires (lc≠1, j not an algebraic integer)",
           r3['j_minpoly_lc'] != 1 and r3['j_is_algebraic_integer'] is False and
           r3.get('j_in_CM_table') is False,
           f"lc={r3['j_minpoly_lc']}")


# ---------------------------------------------------------------------------
def test_member_fiberstack():
    """MEMBER: fiberstack — a_p COPAIR vs naive count + the a_{p²} identity."""
    print("\n=== MEMBER: fiberstack (a_p COPAIR + a_{p²} = a_p²−2p identity) ===")
    from fiberstack import chi_table, legendre_ap_fast, ap2_charsum, \
        rational_heights_menu
    # a_p COPAIR: character-sum engine vs brute-force point count
    # (a_p = p − #{(u,y) ∈ F_p²: y² = u(u−1)(u−lam)})
    pairs, bad = 0, []
    for p in (11, 13, 17, 19):
        chi = chi_table(p)
        for lam in (2, 3, 5, 7):
            if lam % p in (0, 1):
                continue
            naive = 0
            for u in range(p):
                f = u * (u - 1) % p * ((u - lam) % p) % p
                naive += sum(1 for y in range(p) if y * y % p == f)
            a_naive = p - naive
            a_fast = legendre_ap_fast(lam, p, chi)
            pairs += 1
            if a_fast != a_naive or a_fast * a_fast > 4 * p:
                bad.append((lam, p, a_fast, a_naive))
    _check(f"fiberstack a_p COPAIR engine==naive count on {pairs} (lam,p), "
           "Hasse-bounded", pairs >= 12 and not bad, f"bad={bad}")
    # a_{p²} identity: independent F_p² norm-character sum must satisfy
    # a_{p²} = a_p² − 2p for every good-reduction (lam, p)
    pairs2, bad2 = 0, []
    for p in (5, 7, 11, 13):
        chi = chi_table(p)
        for lam in (2, 3, 4):
            if lam % p in (0, 1):
                continue
            ap = legendre_ap_fast(lam, p, chi)
            ap2 = ap2_charsum(lam, p, chi)
            pairs2 += 1
            if ap2 != ap * ap - 2 * p:
                bad2.append((lam, p, ap, ap2))
    _check(f"fiberstack a_{{p²}} = a_p²−2p identity on {pairs2} (lam,p) COPAIRs",
           pairs2 >= 10 and not bad2, f"bad={bad2}")
    # mutation control: a wrong a_{p²} must violate the identity
    chi7 = chi_table(7)
    _check("control: corrupted a_{p²} breaks the identity",
           ap2_charsum(3, 7, chi7) + 1 != legendre_ap_fast(3, 7, chi7)**2 - 14)
    # menu sanity: deterministic, reduced, height-ordered, no (1,1), both signs
    from math import gcd
    menu = rational_heights_menu(3)
    hts = [max(abs(a), b) for a, b in menu]
    _check("fiberstack menu: reduced, no (1,1), signs paired, heights ascending",
           all(gcd(abs(a), b) == 1 for a, b in menu) and
           (1, 1) not in menu and
           all(menu[i + 1] == (-menu[i][0], menu[i][1])
               for i in range(0, len(menu), 2)) and
           hts == sorted(hts) and menu == rational_heights_menu(3),
           f"n={len(menu)}, first={menu[:4]}")


# ---------------------------------------------------------------------------
if __name__ == "__main__":
    reps = {}
    test_unit_jinvariant()
    test_unit_ap_counter()
    test_cm_elliptic_wired()
    test_unit_k3_probe()
    test_member_exactj()
    test_member_fiberstack()
    reps['box1l'] = test_box1l()
    reps['sunrise_eq'] = test_sunrise_eq()
    reps['icc'] = test_icc()
    reps['banana3l_eq'] = test_banana3l_eq()
    reps['crossedbox3l'] = test_crossedbox3l()
    reps['banana3l_1112'] = test_k3_generic_probe()

    print("\n" + "="*72)
    print(f"SUMMARY: {len(PASS)} passed, {len(FAIL)} failed")
    if FAIL:
        print("FAILED:", FAIL)
    print("="*72)
    print("\nROUTE TABLE:")
    print(f"{'graph':14s} {'type':12s} {'genus':6s} {'route':6s} reason")
    for name, rep in reps.items():
        v = rep['variety']
        print(f"  {name:12s} {v['type']:12s} {str(v.get('genus','-')):6s} "
              f"{rep['route']:6s} {rep['reason'][:60]}")

    # write combined output: GEOTRIAGE_TEST_OUT redirects the record off the
    # tracked tree (e.g. to a mktemp file); unset, it lands beside this script
    out_path = (os.environ.get("GEOTRIAGE_TEST_OUT")
                or os.path.join(HERE, "test_output.json"))
    with open(out_path, "w") as f:
        json.dump({'pass': PASS, 'fail': FAIL,
                   'reports': {k: v for k, v in reps.items()}},
                  f, indent=2, default=str)
    sys.exit(1 if FAIL else 0)

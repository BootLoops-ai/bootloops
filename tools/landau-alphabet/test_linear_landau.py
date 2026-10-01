#!/usr/bin/env python3
"""
test_linear_landau.py — regression tests for the mixed quadratic+linear
(worldline / eikonal) Landau engine.

Each test is a sympy-EXACT check (no numerics).  See P2_LANDAU_VERDICT.md
for the physics interpretation.
"""
import sympy as sp
import linear_landau as LL


y, q2, x = sp.symbols("y q2 x")


def _aset(letters):
    return {sp.factor(l) for l in letters}


# ---------------------------------------------------------------------------
# 1.  Symanzik (U,F) for the 1-loop PM box, by hand vs. engine.
# ---------------------------------------------------------------------------
def test_symanzik_1L_box():
    U, F, a, _ = LL.symanzik_mixed(*LL.family_2PM_box())
    a1, a2, a3, a4 = a
    assert sp.expand(U - (a1 + a2)) == 0
    F_expected = a1*a2*q2 - a3**2 - 2*y*a3*a4 - a4**2
    assert sp.expand(F - F_expected) == 0, f"F mismatch: {sp.expand(F-F_expected)}"


def test_symanzik_1L_tri():
    U, F, a, _ = LL.symanzik_mixed(*LL.family_2PM_tri())
    a1, a2, a3 = a
    assert sp.expand(U - (a1 + a2)) == 0
    assert sp.expand(F - (a1*a2*q2 - a3**2)) == 0


# ---------------------------------------------------------------------------
# 2.  1-loop PM box Landau alphabet  ==  {q2, y-1, y+1}  →  {x, 1-x, 1+x}
# ---------------------------------------------------------------------------
def test_alphabet_2PM_box():
    res = LL.landau_locus_mixed(*LL.family_2PM_box())
    assert _aset(res["alphabet"]) == {q2, y - 1, y + 1}
    ax = _aset(LL.pm_alphabet_x(res["alphabet"]))
    assert ax == {x, x - 1, x + 1}, ax


# ---------------------------------------------------------------------------
# 3.  2-loop 3PM IY:  alphabet  ==  {q2, y-1, y+1}  (same as 1-loop)
# ---------------------------------------------------------------------------
def test_alphabet_3PM_IY():
    res = LL.landau_locus_mixed(*LL.family_3PM_IY(), face_timeout=20)
    assert {q2, y - 1, y + 1} <= _aset(res["alphabet"])
    # no spurious high-degree letters
    assert all(sp.Poly(l, y, q2).total_degree() <= 1 for l in res["alphabet"])


# ---------------------------------------------------------------------------
# 4.  2-loop 3PM H:  bounded resultant produces a SPURIOUS `y` and a deg-8
#     polynomial on face (0,3,4,5,6).  Saturated Groebner refutes both.
#     Genuine alphabet = {q2, y-1, y+1}.
# ---------------------------------------------------------------------------
def test_3PM_H_spurious_y_refuted():
    U, F, a, _ = LL.symanzik_mixed(*LL.family_3PM_H())
    Ff = sp.expand(F.subs({a[1]: 0, a[2]: 0}))
    fvars = [a[0], a[3], a[4], a[5], a[6]]
    # y is NOT a genuine Landau component on this face:
    assert LL.verify_letter(Ff, fvars, y, [y, q2]) is False
    # y-1 IS genuine (on the rung+worldline sub-face {a5,a6,a7}):
    Fwl = sp.expand(F.subs({a[0]: 0, a[1]: 0, a[2]: 0, a[3]: 0}))
    assert LL.verify_letter(Fwl, [a[4], a[5], a[6]], y - 1, [y, q2]) is True
    assert LL.verify_letter(Fwl, [a[4], a[5], a[6]], y + 1, [y, q2]) is True


# ---------------------------------------------------------------------------
# 5.  γ=3 (Apéry-K3) singularity:  the leading-Landau discriminant of the
#     topology-#40 LS surface contains  (x²−6x+1)(x²+6x+1) = x⁴−34x²+1,
#     with root x = 3−2√2  ⇔  y = 3.
# ---------------------------------------------------------------------------
def test_apery_gamma3():
    facs, Q6 = LL.apery_K3_landau_locus()
    fset = _aset(facs)
    # the two quadratic factors must be present
    assert (x**2 - 6*x + 1) in fset, fset
    assert (x**2 + 6*x + 1) in fset, fset
    # and (x±1) (the standard thresholds)
    assert (x - 1) in fset and (x + 1) in fset
    # product identity
    assert sp.expand((x**2 - 6*x + 1)*(x**2 + 6*x + 1) - (x**4 - 34*x**2 + 1)) == 0
    # root → γ
    r = 3 - 2*sp.sqrt(2)
    assert sp.simplify((1 + r**2)/(2*r) - 3) == 0


# ---------------------------------------------------------------------------
# 6.  Letter-loss regression:
#     alpha-monomial resultant self-annihilation must NOT lose letters.
#     Face P = -a2*a5*(a7²+2y·a7·a8+a8²): pre-fix the chain returned
#     'resolved' with NO letters although (y−1)(y+1) is in every generator
#     after eliminating a7.
# ---------------------------------------------------------------------------
def test_letterloss_alpha_content_face():
    a2, a5, a7, a8 = sp.symbols("a2 a5 a7 a8", positive=True)
    P = -a2*a5*(a7**2 + 2*y*a7*a8 + a8**2)
    letters, status = LL._discriminant_letters(
        P, [a2, a5, a7, a8], {y, q2}, [a2, a5, a7, a8], face_timeout=30)
    got = _aset(letters)
    assert (y - 1) in got and (y + 1) in got, (status, got)
    # and the saturated-Groebner filter agrees these are genuine here
    assert LL.verify_letter(P, [a2, a5, a7, a8], y - 1, [y, q2]) is True
    assert LL.verify_letter(P, [a2, a5, a7, a8], y + 1, [y, q2]) is True


# ---------------------------------------------------------------------------
# 7.  U is independent of linear-prop α's  (the structural fact).
# ---------------------------------------------------------------------------
def test_U_independent_of_linear_alphas():
    for fam in (LL.family_2PM_box, LL.family_3PM_H, LL.family_3PM_IY):
        quad, lin, *rest = fam()
        U, F, a, _ = LL.symanzik_mixed(quad, lin, *rest)
        nq = len(quad)
        for al in a[nq:]:
            assert al not in U.free_symbols, f"{fam.__name__}: U depends on {al}"


if __name__ == "__main__":
    import sys
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    fail = 0
    for f in fns:
        try:
            f()
            print(f"PASS  {f.__name__}")
        except AssertionError as e:
            fail += 1
            print(f"FAIL  {f.__name__}: {e}")
        except Exception as e:
            fail += 1
            print(f"ERROR {f.__name__}: {type(e).__name__}: {e}")
    print(f"\n{len(fns)-fail}/{len(fns)} passed")
    sys.exit(1 if fail else 0)

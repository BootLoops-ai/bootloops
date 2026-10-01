#!/usr/bin/env python3
"""
test_frobenius.py — unit self-checks for wayfinder/frobenius.py.

Reference values are HAND-DERIVED closed forms, re-derivable from the ODE
(derivations inline below) — no fabricated digits. The resonant test system

    y' = ( M/u + B0 ) y,   M = [[0,1],[0,0]],  B0 = [[0,0],[1,0]]

has a regular-singular point at u=0 with indicial roots {0,0} (Jordan block)
and a genuine log layer: Y(u) = P(u) (I + M log u).

Closed forms used (all verified by plugging the series into the ODE):
  * Sylvester recursion q P_q - M P_q + P_q M = C_q, C_1 = B0 gives
        P_1 = [[1,-2],[1,-1]],   P_2 = [[1/4,-3/4],[1/2,-5/4]].
    ODE cross-check of P_1 (independent of the recursion): column 1 is the
    log-free solution y = [1,0] + u[p,r] + O(u^2); y2' = y1 forces r = 1,
    y1' = y2/u forces p = 1. Column 2 non-log part a(u) = [0,1] + u a_1
    satisfies a_1 - M a_1 = -P_1[:,0] = -[1,1], so a_1 = [-2,-1]. Both match
    P_1 = [[1,-2],[1,-1]].
  * k=1 log-layer coefficient (order u^1 of the log u part of Y):
        P_1 M = [[0,1],[0,1]]  =>  column-2 log coefficient = [1, 1].

SHEARING closed forms (all hand-derived by direct ODE integration; the
Moser/Turrittin shearing chain must reproduce each up to a
constant right factor C — checked by the two-point mixing test
Y(u2) = Yhat(u2) * [Yhat(u1)^{-1} Y(u1)]):
  * GapDiag      A = [[1/u,0],[0,0]]:   Y = [[u,0],[0,1]]            (1 shear)
  * GapResonant  A = [[1/u,1],[0,0]]:   y2=c2; y1'-y1/u=c2 has the log
    particular solution c2*u*log(u)  =>  Y = [[u, u log u],[0,1]]    (1 shear,
    reduced residue = the Jordan block [[0,1],[0,0]] up to scaling)
  * GapConjugated = R * GapResonant * R^{-1}, R=[[2,1],[1,1]] (det 1):
    A = [[2/u-2, -2/u+4],[1/u-1, -1/u+2]],  Y = R*[[u, u log u],[0,1]]
    (nontrivial spectral projector/T; residue [[2,-2],[1,-1]], eig {1,0})
  * Gap2         A = [[2/u,1],[0,0]]:   y1'-2y1/u=c2, particular -c2*u
    =>  Y = [[u^2, -u],[0,1]]                                       (2 shears)

Runnable directly (python3 tests/test_frobenius.py) or under pytest.
"""

import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from mpmath import mp, mpf, mpc  # noqa: E402
from frobenius import frobenius_basis, land  # noqa: E402


class ResonantLog:
    """A(u) = M/u + B0 with M=[[0,1],[0,0]], B0=[[0,0],[1,0]]."""
    n = 2
    var = "u"
    meta = {"singular_points": ["0"]}

    def A(self, u, eps, dps):
        with mp.workdps(dps):
            u = mpc(u)
            return [[mpc(0), 1 / u], [mpc(1), mpc(0)]]


class IntegerGap:
    """M = diag(1,0): nonzero-integer indicial gap -> sheared, Y=[[u,0],[0,1]]."""
    n = 2
    var = "u"
    meta = {"singular_points": ["0"]}

    def A(self, u, eps, dps):
        with mp.workdps(dps):
            u = mpc(u)
            return [[1 / u, mpc(0)], [mpc(0), mpc(0)]]


class GapResonant:
    """A = [[1/u,1],[0,0]]: gap 1 + resonant coupling -> Y = [[u, u log u],[0,1]]."""
    n = 2
    var = "u"
    meta = {"singular_points": ["0"]}

    def A(self, u, eps, dps):
        with mp.workdps(dps):
            u = mpc(u)
            return [[1 / u, mpc(1)], [mpc(0), mpc(0)]]


class GapConjugated:
    """R*(GapResonant)*R^{-1}, R=[[2,1],[1,1]]: nontrivial spectral split.
    Y = [[2u, 2u log u + 1],[u, u log u + 1]]."""
    n = 2
    var = "u"
    meta = {"singular_points": ["0"]}

    def A(self, u, eps, dps):
        with mp.workdps(dps):
            u = mpc(u)
            return [[2 / u - 2, -2 / u + 4], [1 / u - 1, -1 / u + 2]]


class Gap2:
    """A = [[2/u,1],[0,0]]: gap 2 -> TWO shears; Y = [[u^2, -u],[0,1]]."""
    n = 2
    var = "u"
    meta = {"singular_points": ["0"]}

    def A(self, u, eps, dps):
        with mp.workdps(dps):
            u = mpc(u)
            return [[2 / u, mpc(1)], [mpc(0), mpc(0)]]


class DoublePole:
    """A = 1/u^2: irregular singular point -> ValueError."""
    n = 1
    var = "u"
    meta = {"singular_points": ["0"]}

    def A(self, u, eps, dps):
        with mp.workdps(dps):
            u = mpc(u)
            return [[1 / u ** 2]]


def _absdiff(a, b, wp):
    with mp.workdps(wp):
        return abs(mpc(a) - mpc(b))


def test_indicial_and_log_layer():
    """Indicial roots {0,0}, log_max=1, P_1/P_2 and the k=1 log coefficient
    against the hand-derived closed forms above."""
    dps_saved = mp.dps
    dps, kmax = 60, 30
    b = frobenius_basis(ResonantLog(), "0", "0", dps, kmax)
    tol = mpf(10) ** (-(dps - 15))
    # indicial roots: double root at 0
    assert len(b["indicial"]) == 2
    for lam in b["indicial"]:
        assert abs(lam) < tol, f"indicial root {lam} != 0"
    assert b["log_max"] == 1, f"log_max {b['log_max']} != 1 (Jordan block)"
    # P_1, P_2 exact rational closed forms
    P1_ref = [[mpf(1), mpf(-2)], [mpf(1), mpf(-1)]]
    P2_ref = [[mpf(1) / 4, mpf(-3) / 4], [mpf(1) / 2, mpf(-5) / 4]]
    for i in range(2):
        for j in range(2):
            assert _absdiff(b["P"][1][i][j], P1_ref[i][j], dps) < tol, \
                f"P1[{i}][{j}] = {b['P'][1][i][j]} != {P1_ref[i][j]}"
            assert _absdiff(b["P"][2][i][j], P2_ref[i][j], dps) < tol, \
                f"P2[{i}][{j}] = {b['P'][2][i][j]} != {P2_ref[i][j]}"
    # k=1 log-layer coefficient: P_1 M, column 2 => [1, 1]
    with mp.workdps(dps):
        M = b["M"]
        P1 = b["P"][1]
        # (P_1 M) column 2 entries, rows i=0,1
        c2 = [sum(P1[i][k] * M[k][1] for k in range(2)) for i in range(2)]
    for i, ref in enumerate([mpf(1), mpf(1)]):
        assert _absdiff(c2[i], ref, dps) < tol, \
            f"k=1 log coeff col2[{i}] = {c2[i]} != {ref}"
    assert mp.dps == dps_saved, "global mp.dps was mutated (import-dps footgun)"
    print("  indicial {0,0}, log_max=1, P1/P2 and k=1 log coefficient OK")


def test_eval_log_structure():
    """Y(u) column 2 must contain [1,0]*log(u) at leading order: check
    Y(u)[0][1] - (P(u)[:,1] + log(u)*P(u)[:,0]) ~ 0 at a small u via the
    series definition, using only M, P from the dict (no eval internals)."""
    dps_saved = mp.dps
    dps, kmax = 60, 30
    b = frobenius_basis(ResonantLog(), "0", "0", dps, kmax)
    with mp.workdps(dps):
        u = mpf(1) / 32
        Y = b["eval"](u, dps)
        # reconstruct: Y = P(u) (I + M log u); column 2 = P(u)[:,1] + log u P(u)[:,0]
        lu = mp.log(u)
        Pu = [[mpc(0)] * 2 for _ in range(2)]
        up = mpf(1)
        for q in range(kmax + 1):
            for i in range(2):
                for j in range(2):
                    Pu[i][j] += b["P"][q][i][j] * up
            up *= u
        ref_col2 = [Pu[i][1] + lu * Pu[i][0] for i in range(2)]
        for i in range(2):
            d = abs(Y[i][1] - ref_col2[i])
            assert d < mpf(10) ** (-(dps - 12)), \
                f"eval col2[{i}] mismatch {d} vs series reconstruction"
    assert mp.dps == dps_saved
    print("  eval() log structure matches P(u)(I + M log u) reconstruction OK")


def test_land_residual_and_recovery():
    """land(): transport a known basis combination from x_from=1/4 and match
    back; residual < 10^-(dps-10) and kappa recovered."""
    dps_saved = mp.dps
    dps, kmax = 60, 40
    b = frobenius_basis(ResonantLog(), "0", "0", dps, kmax)
    with mp.workdps(dps + 30):
        Y0 = b["eval"](mpf(1) / 4, dps + 30)
        kt = [mpf(2) / 7, mpf(-3) / 5]
        yf = [Y0[i][0] * kt[0] + Y0[i][1] * kt[1] for i in range(2)]
    r = land(ResonantLog(), "0", "0.25", yf, "0", dps, kmax)
    with mp.workdps(dps + 10):
        assert mpf(r["residual_rel"]) < mpf(10) ** (-(dps - 10)), \
            f"land residual_rel {r['residual_rel']} >= 1e-{dps - 10}"
        for i in range(2):
            d = abs(r["kappa"][i] - kt[i])
            assert d < mpf(10) ** (-(dps - 15)), \
                f"kappa[{i}] off by {d}"
    assert mp.dps == dps_saved, "global mp.dps was mutated (import-dps footgun)"
    print(f"  land(): residual_rel={mp.nstr(mpf(r['residual_rel']), 3)} "
          f"< 1e-{dps - 10}, kappa recovered OK")


def test_policy_and_tripwires():
    """Real-eps0 policy (the flagged Im-linear defect) and
    double-pole ValueError."""
    dps_saved = mp.dps
    try:
        frobenius_basis(ResonantLog(), "0.1+0.2j", "0", 50, 10)
        raise SystemExit("complex eps0 accepted — policy violation")
    except ValueError as e:
        assert "Im-linear" in str(e)
    try:
        frobenius_basis(DoublePole(), "0", "0", 50, 10)
        raise SystemExit("double pole accepted — not regular-singular")
    except ValueError as e:
        assert "regular-singular" in str(e)
    assert mp.dps == dps_saved
    print("  policy tripwires (complex-eps / double-pole) OK")


def _check_sheared_closed_form(desys, yhat_fn, n_shears_ref, log_max_ref,
                               ind_ref, dps=60, kmax=30):
    """Shearing correctness up to a constant right factor C (any fundamental
    basis of the same system satisfies Y(u) = Yhat(u)*C): fit C at u1, check
    at u2. Also pins n_shears, log_max and both indicial multisets."""
    b = frobenius_basis(desys, "0", "0", dps, kmax)
    assert b["sheared"] is True
    assert b["n_shears"] == n_shears_ref, \
        f"n_shears {b['n_shears']} != {n_shears_ref}"
    assert b["log_max"] == log_max_ref, \
        f"log_max {b['log_max']} != {log_max_ref}"
    tol = mpf(10) ** (-(dps - 15))
    # original indicial multiset (sorted by real part)
    got = sorted(b["indicial"], key=lambda z: (z.real, z.imag))
    for g, rf in zip(got, sorted(ind_ref, key=lambda z: (z.real, z.imag))):
        assert abs(g - rf) < tol, f"indicial {g} != {rf}"
    # reduced residue must be gap-free: all integer differences are 0
    for ei in b["indicial_reduced"]:
        for ej in b["indicial_reduced"]:
            k = int(mp.nint((ei - ej).real))
            assert not (k != 0 and abs(ei - ej - k) < mpf(10) ** (-15)), \
                f"reduced residue still gapped: {ei} vs {ej}"
    with mp.workdps(dps + 20):
        u1, u2 = mpf(1) / 32, mpf(1) / 48
        Y1 = mp.matrix(b["eval"](u1, dps + 20))
        C = mp.inverse(mp.matrix(yhat_fn(u1))) * Y1
        R = mp.matrix(yhat_fn(u2)) * C - mp.matrix(b["eval"](u2, dps + 20))
        err = max(abs(R[i, j]) for i in range(2) for j in range(2))
        assert err < mpf(10) ** (-(dps - 15)), \
            f"closed-form mixing error {mp.nstr(err, 3)}"
    return b


def test_shearing_closed_forms():
    """The Moser/Turrittin shearing chain vs the four HAND-DERIVED 2x2
    closed forms in the module docstring (gap 1 diagonal / resonant-log /
    conjugated resonant-log, and the two-shear gap-2 case)."""
    dps_saved = mp.dps
    _check_sheared_closed_form(
        IntegerGap(),
        lambda u: [[u, mpc(0)], [mpc(0), mpc(1)]],
        n_shears_ref=1, log_max_ref=0, ind_ref=[mpc(0), mpc(1)])
    _check_sheared_closed_form(
        GapResonant(),
        lambda u: [[u, u * mp.log(u)], [mpc(0), mpc(1)]],
        n_shears_ref=1, log_max_ref=1, ind_ref=[mpc(0), mpc(1)])
    _check_sheared_closed_form(
        GapConjugated(),
        lambda u: [[2 * u, 2 * u * mp.log(u) + 1],
                   [u, u * mp.log(u) + 1]],
        n_shears_ref=1, log_max_ref=1, ind_ref=[mpc(0), mpc(1)])
    _check_sheared_closed_form(
        Gap2(),
        lambda u: [[u ** 2, -u], [mpc(0), mpc(1)]],
        n_shears_ref=2, log_max_ref=0, ind_ref=[mpc(0), mpc(2)])
    assert mp.dps == dps_saved, "global mp.dps was mutated (import-dps footgun)"
    print("  shearing vs 4 hand-derived closed forms (1- and 2-shear) OK")


def test_land_on_sheared_basis():
    """land() end-to-end THROUGH the shearing path: seed the exact solution
    [u + u log u, 1] of GapResonant at x=1/4, land at u=0, then reproduce
    the solution at an independent point u=1/40 from kappa."""
    dps_saved = mp.dps
    dps, kmax = 55, 40
    with mp.workdps(dps + 30):
        x0 = mpf(1) / 4
        yf = [x0 + x0 * mp.log(x0), mpc(1)]   # Yhat(1/4) . [1, 1]
    r = land(GapResonant(), "0", "0.25", yf, "0", dps, kmax)
    assert r["basis"]["sheared"] is True
    with mp.workdps(dps + 20):
        assert mpf(r["residual_rel"]) < mpf(10) ** (-(dps - 10)), \
            f"land residual_rel {r['residual_rel']}"
        uc = mpf(1) / 40
        Yc = mp.matrix(r["basis"]["eval"](uc, dps + 20))
        yp = Yc * mp.matrix([[r["kappa"][0]], [r["kappa"][1]]])
        yt = [uc + uc * mp.log(uc), mpc(1)]
        err = max(abs(yp[i] - yt[i]) for i in range(2))
        assert err < mpf(10) ** (-(dps - 12)), \
            f"sheared-basis solution reproduction error {mp.nstr(err, 3)}"
    assert mp.dps == dps_saved, "global mp.dps was mutated (import-dps footgun)"
    print(f"  land() through shearing: residual_rel="
          f"{mp.nstr(mpf(r['residual_rel']), 3)}, reproduction OK")


def main():
    t0 = time.time()
    print("test_frobenius.py self-checks:")
    test_indicial_and_log_layer()
    test_eval_log_structure()
    test_land_residual_and_recovery()
    test_policy_and_tripwires()
    test_shearing_closed_forms()
    test_land_on_sheared_basis()
    print(f"ALL PASS ({time.time() - t0:.1f}s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())

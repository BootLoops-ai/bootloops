#!/usr/bin/env python3
"""popcorn acceptance battery.

Default tier: fast real checks (~40 s on a laptop) — pins, the hardened CLI
selftest and two-route check, exact/ball SFS routes vs the per-i hyp1f1
oracle, the dominance oracle's h=1/2 collapse gate, the q-series evaluator vs
pinned reference values, a small certified DFE kernel with a finite-difference
gradient control, the LP/positivity certificate family on toy hulls with
planted negative controls, the region-certificate selftest (which carries
its own planted MUST-FAIL controls), the exact Lambda-coalescent identities and
n = 20 reference grid, the two-locus chain against exact Kingman moments and
closed forms, the transient PDE engine against the certified stationary
references and its self-convergence table, the EPO-join fixture (library
+ CLI), the exact linked 2-SFS engine with its Kingman kernel identities, the
shipped class certificates re-derived exactly, a live class verdict and the
Monte Carlo table, the folded two-site spectra (exact fold map,
convex/conic/diagonal certificates re-verified, live deciders at n = 4), and
the two-site projection operator against explicit enumeration plus the
two-site spectrum estimator against its reference (library + CLI), the
certified transient-SFS enclosures (live M-ladder cells, the certified
stationary objects against popcorn.sfs's exact route to >= 30 digits, the
six shipped cells and the independent mpmath route), the planted truths
against the closed form, the exact Lambda spectra, the two-locus engine and
a seeded reduced simulation, and the two-window simulation harness against
the exact linked Kingman spectrum and the two-locus moments.

--full adds the heavier gates: the held-out two-route gate (gate_phase1.py,
minutes), the dominance oracle's full selftest (minutes), the n=20
certified kernel selftest (minutes), the M = 200 certified-cell
reproduction (~25 s) and the CLI regeneration of two planted truths and the
gene-conversion arm (~10 s, msprime).

Every leg that needs an optional engine SKIPS BY NAME, stating what to
install, when that engine is absent. Exit 0 iff no leg FAILs (skips do not
fail the battery). Writes nothing into the package tree: scratch goes to a
temp directory.

Usage:  python3 tools/popcorn/selftest.py [--full]
"""
import argparse
import os
import shutil
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
PARENT = os.path.dirname(HERE)
for p in (PARENT, HERE):
    if p not in sys.path:
        sys.path.insert(0, p)

PY = sys.executable or "python3"

RESULTS = []


def leg(name):
    def deco(fn):
        fn._leg_name = name
        return fn
    return deco


def run_leg(fn, *args):
    name = fn._leg_name
    t0 = time.time()
    try:
        out = fn(*args)
        status, detail = ("SKIP", out[5:]) if isinstance(out, str) and out.startswith("SKIP:") \
            else ("PASS", out if isinstance(out, str) else "")
    except AssertionError as e:
        status, detail = "FAIL", str(e)
    except Exception as e:
        status, detail = "FAIL", f"{type(e).__name__}: {e}"
    wall = time.time() - t0
    print(f"LEG {name:<22s} {status}  {wall:6.1f}s  {detail}")
    RESULTS.append((name, status))
    return status


def _has(modname):
    try:
        __import__(modname)
        return True
    except Exception:
        return False


# --------------------------------------------------------------------- legs
@leg("pins")
def leg_pins():
    import popcorn
    bad = popcorn.verify()
    assert not bad, f"pin mismatches: {bad}"
    return "all shipped pins byte-identical"


@leg("cli_selftest")
def leg_cli_selftest():
    r = subprocess.run([PY, os.path.join(HERE, "certsfs.py"),
                        "selftest", "--dps", "40"],
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                       text=True, timeout=300)
    assert r.returncode == 0 and "SELFTEST PASS" in r.stdout, \
        f"rc={r.returncode}\n{r.stdout[-800:]}"
    return "sha pins + 2 reference values digit-for-digit"


@leg("cli_check")
def leg_cli_check():
    r = subprocess.run([PY, os.path.join(HERE, "certsfs.py"),
                        "check", "100", "50", "-1000", "--dps", "40"],
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                       text=True, timeout=300)
    assert r.returncode == 0 and "agree" in r.stdout, \
        f"rc={r.returncode}\n{r.stdout[-800:]}"
    return "two independent routes at n=100, S=-1000"


@leg("sfs_exact")
def leg_sfs_exact():
    from fractions import Fraction
    from mpmath import mp, mpf, hyp1f1, fabs, log10
    from popcorn.sfs import solve_M_exact, solve_dM_exact, assemble_M
    n, S = 50, Fraction(-10)          # the |S|~10 directional-leak regime
    alpha, beta = solve_M_exact(n, S)
    M = assemble_M(alpha, beta, S, 70)
    ders = solve_dM_exact(n, S, alpha, beta, 1)
    G = assemble_M(ders[0][0], ders[0][1], S, 70)
    with mp.workdps(95):
        worstM = worstG = 9999.0
        for i in (1, 10, 25, 40, 49):
            truth = hyp1f1(n - i, n, mpf(10))
            if M[i] != truth:
                worstM = min(worstM, float(-log10(fabs((M[i] - truth) / truth))))
        for i in (1, 25, 49):
            truthG = -(mpf(n - i) / n) * hyp1f1(n - i + 1, n + 1, mpf(10))
            if G[i] != truthG:
                worstG = min(worstG, float(-log10(fabs((G[i] - truthG) / truthG))))
    assert worstM >= 50, f"M vs hyp1f1 only {worstM:.1f}d (bar 50)"
    assert worstG >= 50, f"grad vs closed form only {worstG:.1f}d (bar 50)"
    return f"n=50 S=-10: M {worstM:.0f}d, grad {worstG:.0f}d vs hyp1f1"


@leg("sfs_arb")
def leg_sfs_arb():
    if not _has("flint"):
        return "SKIP:python-flint absent — pip install python-flint"
    from fractions import Fraction
    from mpmath import mp, mpf, hyp1f1, fabs, log10
    from popcorn.sfs import arb_route
    n, S = 300, Fraction(-50)
    sol = arb_route.solve_vector_arb(n, float(S), 40, want_grad=True, S_frac=S)
    with mp.workdps(70):
        worst = 9999.0
        for i in (1, 75, 150, 299):
            got = mp.mpmathify(sol["M"][i].mid().str(60, radius=False))
            truth = hyp1f1(n - i, n, mpf(50))
            if got != truth:
                worst = min(worst, float(-log10(fabs((got - truth) / truth))))
    assert sol["achieved_dps"] >= 40, f"achieved {sol['achieved_dps']:.1f} < 40"
    assert worst >= 35, f"ball mid vs hyp1f1 only {worst:.1f}d (bar 35)"
    return f"n=300 ball route: {worst:.0f}d vs hyp1f1"


@leg("dominance_oracle")
def leg_dominance_oracle():
    from mpmath import mp, mpf, hyp1f1, expm1, fabs, log10
    from popcorn.dominance import E_dominant
    with mp.workdps(80):
        got, sc = E_dominant(20, 3, -100.0, 0.5, dps=50)
        Sm = mpf('-100')
        truth = mpf(20) / (mpf(3) * 17) * (1 - hyp1f1(17, 20, -Sm)) / (-expm1(-Sm))
        d = float(-log10(fabs((got - truth) / truth))) if got != truth else 9999
    assert d >= 45, f"h=1/2 collapse only {d:.1f}d (bar 45)"
    got2, sc2 = E_dominant(20, 3, -5.0, 0.05, dps=50)
    assert sc2 >= 45, f"hostile-point selfcons only {sc2:.1f}d (bar 45)"
    return f"collapse {d:.0f}d, hostile selfcons {min(sc2, 9999):.0f}d"


@leg("dominance_qseries")
def leg_dominance_qseries():
    import json
    from mpmath import mp, mpf, fabs, log10
    from popcorn.dominance import E_dominant_qseries
    bank = json.load(open(os.path.join(HERE, "REFERENCE_VALUES.json")))
    worst = 9999.0
    for (n, k, S, h) in [(20, 3, -5.0, 0.05), (20, 10, 10.0, 0.3)]:
        rec = [p for p in bank['points']
               if p['n'] == n and p['k'] == k and p['S'] == S and p['h'] == h]
        assert len(rec) == 1, f"reference point {(n, k, S, h)} not unique"
        E = E_dominant_qseries(n, k, S, h, dps=30)
        with mp.workdps(45):
            oracle = mpf(rec[0]['oracle'])
            if E != oracle:
                worst = min(worst, float(-log10(fabs((E - oracle) / oracle))))
    assert worst >= 25, f"q-series vs pinned oracle only {worst:.1f}d (bar 25)"
    return f"2 points vs pinned oracle: {worst:.0f}d"


@leg("dfe_kernel")
def leg_dfe_kernel():
    if not _has("flint"):
        return "SKIP:python-flint absent — pip install python-flint"
    from mpmath import mp, mpf, fabs, log10
    from popcorn.dfe import CertifiedKernel
    K = CertifiedKernel(8, dps=25)
    K.build(-1)
    K.build(1)
    Ev, sc, tail = K.mix({'b': '0.4', 'Sd': '-100', 'pb': 0.02, 'Sb': 10.0})
    assert sc >= 20, f"degree-pair selfcons only {sc:.1f}d (bar 20)"
    assert tail < 1e-30, f"tail bound {tail:.1e} not < 1e-30"
    g = K.grad_mix({'b': 0.4, 'Sd': -100.0, 'pb': 0.02, 'Sb': 10.0})
    with mp.workdps(45):
        h = mpf('1e-10')
        E1, _, _ = K.mix({'b': mpf('0.4') + h, 'Sd': -100.0, 'pb': 0.02,
                          'Sb': 10.0}, dps_check=False)
        E0, _, _ = K.mix({'b': mpf('0.4') - h, 'Sd': -100.0, 'pb': 0.02,
                          'Sb': 10.0}, dps_check=False)
        fd = (E1[4] - E0[4]) / (2 * h)
        gd = float(-log10(fabs((g['b'][4] - fd) / fd))) if fd != 0 else 0
    assert gd >= 8, f"analytic grad vs FD only {gd:.1f}d (bar 8)"
    return f"n=8 kernel: selfcons {min(sc, 9999):.0f}d, grad vs FD {gd:.0f}d"


@leg("certificates_lp")
def leg_certificates_lp():
    if not _has("gmpy2"):
        return "SKIP:gmpy2 absent — pip install gmpy2"
    from gmpy2 import mpq
    from popcorn import certificates as C
    pts = [(mpq(1, 3), mpq(1, 3), mpq(1, 3)), (mpq(1, 2), mpq(1, 4), mpq(1, 4)),
           (mpq(1, 4), mpq(1, 2), mpq(1, 4)), (mpq(1, 4), mpq(1, 4), mpq(1, 2)),
           (mpq(2, 5), mpq(2, 5), mpq(1, 5))]
    q_in = (mpq(3, 8), mpq(3, 8), mpq(1, 4))          # inside the toy hull
    q_out = (mpq(9, 10), mpq(1, 20), mpq(1, 20))      # far outside
    done = []
    st, _ = C.exact_lp.membership(pts, q_in)
    assert st == 'FEASIBLE', f"exact_lp: {st} for interior point"
    st, _ = C.exact_lp.membership(pts, q_out)
    assert st == 'INFEASIBLE', f"exact_lp: {st} for planted exterior point"
    done.append("exact_lp")
    if _has("numpy"):
        st, _ = C.fast_lp.membership_fast(pts, q_in)
        assert st == 'FEASIBLE', f"fast_lp: {st}"
        st, _ = C.fast_lp.membership_fast(pts, q_out)
        assert st == 'INFEASIBLE', f"fast_lp: {st}"
        done.append("fast_lp")
    if _has("flint"):
        r = C.cone_lp_b.membership_exact(pts, q_in)
        assert r[0] == 'FEASIBLE', f"cone_lp_b: {r[0]}"
        r = C.cone_lp_b.membership_exact(pts, q_out)
        assert r[0] == 'INFEASIBLE', f"cone_lp_b: {r[0]}"
        done.append("cone_lp_b")
        if _has("sympy"):
            pos = {0: mpq(1, 4), 1: mpq(-1), 2: mpq(1)}   # (t-1/2)^2 >= 0
            v, info, _dips = C.positivity.nonneg_route_R(pos)
            assert v is True, f"positivity route R: {v} ({info})"
            vb, _ib = C.positivity.nonneg_route_B(pos, touch_points=(mpq(1, 2),))
            assert vb is True, f"positivity route B: {vb}"
            neg = {0: mpq(-1, 8), 1: mpq(1)}              # planted negative
            vn, _i, _d = C.positivity.nonneg_route_R(neg)
            vb2, _i2 = C.positivity.nonneg_route_B(neg)
            assert vn is False and vb2 is False, \
                "planted negative poly not refused"
            done.append("positivity")
    else:
        done.append("(cone_lp_b, positivity skipped: python-flint absent)")
    if _has("scipy") and _has("flint"):
        H = C.CertHull(pts)
        st, _pl, _b, route = H.decide(q_in)
        assert st == 'FEASIBLE', f"CertHull: {st} via {route}"
        st, _pl, _b, route = H.decide(q_out)
        assert st == 'INFEASIBLE', f"CertHull: {st} via {route}"
        done.append("cert_lp")
    else:
        done.append("(cert_lp skipped: scipy or python-flint absent)")
    return "; ".join(done)


@leg("region_certs")
def leg_region_certs():
    if not _has("sympy"):
        return "SKIP:sympy absent — pip install sympy"
    r = subprocess.run([PY, os.path.join(HERE, "region_certificates.py")],
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                       text=True, timeout=300)
    assert r.returncode == 0 and "REGION CERTS PASS" in r.stdout, \
        f"rc={r.returncode}\n{r.stdout[-800:]}"
    return "incl. planted MUST-FAIL controls"


@leg("lambda_kingman")
def leg_lambda_kingman():
    from fractions import Fraction
    from popcorn.lambda_coalescent import expected_lengths, xi_hat, kingman_rate
    for n in (20, 60):
        h = expected_lengths(n, kingman_rate())[n]
        bad = [i for i in range(1, n) if h[i] != Fraction(2, i)]
        assert not bad, f"n={n}: h(n,i) != 2/i at i={bad[:5]}"
    H = sum(Fraction(1, j) for j in range(1, 20))
    x = xi_hat(20, kingman_rate())
    assert all(x[i - 1] == Fraction(1, i) / H for i in range(1, 20)), \
        "normalized Kingman spectrum != (1/i)/H_19"
    assert sum(x) == 1, "xi_hat does not sum to 1 exactly"
    return "h(n,i) == 2/i exact at n=20,60; xi_hat == (1/i)/H_19"


@leg("lambda_identities")
def leg_lambda_identities():
    import json
    from fractions import Fraction
    from math import factorial
    from popcorn.lambda_coalescent import identity_checks, beta_rate
    ref = json.load(open(os.path.join(HERE, "reference", "lambda_coalescent",
                                       "spectra_n20.json")))
    c = identity_checks(20)
    names = ("check1_kingman_closed_form", "check2_beta_alpha2_is_kingman",
             "check3_dirac_psi0_is_kingman",
             "check4_beta_alpha1_is_bolthausen_sznitman",
             "check5_dirac_psi1_is_star")
    for k in names:
        assert c[k] is True, f"{k} live: {c[k]!r}"
        assert ref["checks"][k] is True, f"{k} reference: {ref['checks'][k]!r}"
    k6 = "check6_beta_alpha199_100_linf_to_kingman_float"
    assert c[k6] == ref["checks"][k6], \
        f"check 6 float {c[k6]!r} != reference {ref['checks'][k6]!r}"
    l1 = beta_rate(1)
    assert all(l1(20, k) == Fraction(factorial(k - 2) * factorial(20 - k),
                                    factorial(19)) for k in range(2, 21)), \
        "beta_rate(1)(20,k) != (k-2)!(20-k)!/19!"
    return f"checks 1-5 exact True; check 6 == {ref['checks'][k6]!r}; BS rates at b=20"


@leg("lambda_spectra_ref")
def leg_lambda_spectra_ref():
    import json
    from popcorn.lambda_coalescent import xi_hat, reference_grid
    ref = json.load(open(os.path.join(HERE, "reference", "lambda_coalescent",
                                       "spectra_n20.json")))
    rows = {(r["family"], r["value"]): r for r in ref["grid"]}
    grid = reference_grid()
    assert len(rows) == len(ref["grid"]) == len(grid) == 14, \
        f"reference grid size {len(ref['grid'])}, live {len(grid)}"
    assert ("beta", "3/2") in rows, "Beta(3/2) reference row missing"
    for fam, _pname, pval, lam in grid:
        key = (fam, str(pval) if pval is not None else None)
        assert key in rows, f"no reference row for {key}"
        x = xi_hat(20, lam)
        assert len(x) == 19 and sum(x) == 1, f"{key}: len {len(x)} / sum != 1"
        assert [str(v) for v in x] == rows[key]["xi_hat"], \
            f"{key}: exact spectrum differs from reference"
    return "14 rows (Kingman, 7 Beta incl. 3/2, 6 Dirac) string-equal, each sums to 1"


@leg("lambda_functional")
def leg_lambda_functional():
    import json
    from fractions import Fraction
    from popcorn.lambda_coalescent import (xi_hat, beta_rate, kingman_rate,
                                          load_witness, witness_margin)
    ref = json.load(open(os.path.join(HERE, "reference", "lambda_coalescent",
                                       "spectra_n20.json")))
    w, w0 = load_witness(os.path.join(HERE, "reference", "lambda_coalescent",
                                      "kingman_class_functional_n20.json"))
    assert len(w) == 19 and w0 == Fraction(-4898883, 1 << 26)
    rb = [r for r in ref["grid"] if r["family"] == "beta" and r["value"] == "3/2"]
    rk = [r for r in ref["grid"] if r["family"] == "kingman"]
    assert len(rb) == 1 and len(rk) == 1
    mb = witness_margin(w, w0, xi_hat(20, beta_rate(Fraction(3, 2))))
    mk = witness_margin(w, w0, xi_hat(20, kingman_rate()))
    assert str(mb) == rb[0]["witness_margin_exact"], "Beta(3/2) margin != reference"
    assert str(mk) == rk[0]["witness_margin_exact"], "Kingman margin != reference"
    assert mb < 0 and rb[0]["witness_sign"] == "<0", "Beta(3/2) margin not < 0"
    assert mk >= 0 and rk[0]["witness_sign"] == ">=0", "Kingman margin not >= 0"
    return f"exact margins == reference; Beta(3/2) F={float(mb):.3e} <0, Kingman F={float(mk):.3e} >=0"


@leg("twolocus_kingman")
def leg_twolocus_kingman():
    if not (_has("numpy") and _has("scipy")):
        return "SKIP:numpy/scipy absent — pip install numpy scipy"
    import json
    from fractions import Fraction as Fr
    import numpy as np
    from popcorn.twolocus import TwoLocus, pool_pairs
    ref = json.load(open(os.path.join(HERE, "reference", "twolocus",
                                      "kingman_n8_exact.json")))
    n = int(ref["n"])
    E = TwoLocus(n)
    assert E.N == 8405, f"n=8 reachable configurations {E.N} != 8405"
    M0, a0, b0 = E.moments(0.0, [(0.0, 1.0)])          # rho = 0, constant N
    two_over_i = 2.0 / np.arange(1, n)
    d1 = float(max(np.max(np.abs(a0 - two_over_i)), np.max(np.abs(b0 - two_over_i))))
    assert d1 < 1e-12, f"max|E[T_i]-2/i| = {d1:.2e} (bar 1e-12)"
    EL2 = float(Fr(ref["E_L2_analytic"]))
    dL2 = abs(float(M0.sum()) - EL2) / EL2              # sum over ordered (i,j) = E[L^2]
    assert dL2 < 1e-12, f"E[L^2] rel dev {dL2:.2e} (bar 1e-12)"
    keys = [tuple(k) for k in ref["keys"]]
    ex = {(i, j): float(Fr(ref["E_LiLj"][f"{i},{j}"])) for (i, j) in keys}
    draw = max(abs(M0[i - 1, j - 1] - ex[(i, j)]) / ex[(i, j)] for (i, j) in keys)
    assert draw < 1e-11, f"raw M vs exact E[L_i L_j]: max rel {draw:.2e} (bar 1e-11)"
    q = np.array([float(Fr(ref["pool_pairs"][f"{i},{j}"])) for (i, j) in keys])
    dq = float(np.max(np.abs(pool_pairs(M0) - q)))
    assert dq < 1e-12, f"pooled vector vs exact: max abs {dq:.2e} (bar 1e-12)"
    dsym = float(np.max(np.abs(M0 - M0.T)))
    assert dsym < 1e-12, f"max|M-M^T| = {dsym:.2e} (bar 1e-12)"
    return (f"n=8 rho=0 vs exact Kingman: raw rel {draw:.1e}, pooled {dq:.1e}, "
            f"|ET_i-2/i| {d1:.1e}, E[L^2] rel {dL2:.1e}")


@leg("twolocus_rho")
def leg_twolocus_rho():
    if not (_has("numpy") and _has("scipy")):
        return "SKIP:numpy/scipy absent — pip install numpy scipy"
    from fractions import Fraction as Fr
    from math import exp
    import numpy as np
    from popcorn.twolocus import TwoLocus
    # (a) n = 2, constant N, all rho: M = 4*E[T^A T^B], E[T^A T^B] = 1 + (rho+18)/(rho^2+13 rho+18)
    E2 = TwoLocus(2)
    w2 = 0.0
    for rho in (Fr(0), Fr(1, 2), Fr(1), Fr(3), Fr(10), Fr(100)):
        M, a, b = E2.moments(float(rho), [(0.0, 1.0)])
        exact = float(4 * (1 + (rho + 18) / (rho * rho + 13 * rho + 18)))
        w2 = max(w2, abs(M[0, 0] - exact) / exact)
        assert abs(a[0] - 2) < 1e-12 and abs(b[0] - 2) < 1e-12, f"n=2 E[T_1] != 2 at rho={rho}"
    assert w2 < 1e-12, f"n=2 closed form: max rel {w2:.2e} (bar 1e-12)"
    # (b) n = 2, two-epoch N(t) [(0, e1), (t1, e2)]: rho = 0 second moment and any-rho first moment in closed form
    e1, t1, e2 = 3.0, 0.4, 0.5
    q = exp(-e1 * t1)
    ET = (1 - q) / e1 + q / e2                       # int P(T>t) dt, hazard e1 on [0,t1), e2 after
    ET2 = 2 * (1 / e1 ** 2 - q * (t1 / e1 + 1 / e1 ** 2)) + q * (2 * t1 / e2 + 2 / e2 ** 2)   # int 2t P(T>t) dt
    Mh, ah, _ = E2.moments(0.0, [(0.0, e1), (t1, e2)])
    dh2 = abs(Mh[0, 0] - 4 * ET2) / (4 * ET2)
    _Mr, ar, _ = E2.moments(5.0, [(0.0, e1), (t1, e2)])
    dh1 = max(abs(ah[0] - 2 * ET), abs(ar[0] - 2 * ET)) / (2 * ET)
    assert dh2 < 1e-11 and dh1 < 1e-11, f"n=2 two-epoch: E[T^2] rel {dh2:.2e}, E[T] rel {dh1:.2e} (bar 1e-11)"
    # (c) n = 6 self-consistency at rho > 0
    n, rho, c = 6, 1.5, 2.5
    E = TwoLocus(n)
    assert E.N == 1042, f"n=6 reachable configurations {E.N} != 1042"
    two_over_i = 2.0 / np.arange(1, n)
    M1, a1, b1 = E.moments(rho, [(0.0, 1.0)])
    scale = float(np.max(np.abs(M1)))
    dm = float(max(np.max(np.abs(a1 - two_over_i)), np.max(np.abs(b1 - two_over_i))))
    assert dm < 1e-12, f"Kingman marginals at rho={rho}: {dm:.2e} (bar 1e-12)"
    dsym = float(np.max(np.abs(M1 - M1.T))) / scale
    assert dsym < 1e-12, f"A/B symmetry: {dsym:.2e} (bar 1e-12)"
    Ms, _, _ = E.moments(rho, [(0.0, 1.0), (0.7, 1.0)])     # constant epoch split in two: expm path == LU path
    dsplit = float(np.max(np.abs(Ms - M1))) / scale
    assert dsplit < 1e-11, f"epoch-split invariance: {dsplit:.2e} (bar 1e-11)"
    Mc, ac, _ = E.moments(rho * c, [(0.0, c)])             # time rescaling: M(c*rho; eta=c) = M(rho; 1)/c^2
    dscale = max(float(np.max(np.abs(Mc * c * c - M1))) / scale, float(np.max(np.abs(ac * c - a1))))
    assert dscale < 1e-11, f"time-rescaling covariance: {dscale:.2e} (bar 1e-11)"
    Mhi, ahi, bhi = E.moments(40.0, [(0.0, 1.0)])           # toward unlinked loci: M -> outer(m1A, m1B)
    gap1 = float(np.max(np.abs(M1 - np.outer(a1, b1))))
    gap2 = float(np.max(np.abs(Mhi - np.outer(ahi, bhi))))
    assert gap2 < 0.2 * gap1, f"linkage decay: |M-outer| {gap1:.3f} (rho={rho}) -> {gap2:.3f} (rho=40)"
    return (f"n=2 closed form {w2:.1e} (6 rho), n=2 two-epoch {max(dh1, dh2):.1e}; n=6 rho={rho}: "
            f"marginals {dm:.1e}, sym {dsym:.1e}, split {dsplit:.1e}, rescale {dscale:.1e}, "
            f"linkage gap {gap1:.2f}->{gap2:.3f}")


@leg("transient_stationary")
def leg_transient_stationary():
    if not (_has("numpy") and _has("scipy")):
        return "SKIP:numpy/scipy absent — pip install numpy scipy"
    import json
    import numpy as np
    from mpmath import mp, mpf, hyp1f1, expm1, fabs, log10
    from popcorn.transient import TransientSFSEngine
    ref_dir = os.path.join(HERE, "reference", "transient")
    tab = json.load(open(os.path.join(ref_dir, "eq_valid_n1000.json")))
    bands = {}
    for b in tab["bands"]:                      # "i=1", "i=2-5", ... -> slices
        lo, _, hi = b[2:].partition("-")
        bands[b] = (int(lo) - 1, int(hi or lo))
    eng = TransientSFSEngine(n=1000)
    assert eng.G == tab["G"], f"default grid G={eng.G} != {tab['G']}"
    # the shipped S=-5 reference vs the certified closed form (mpmath), spot i
    ref5 = json.load(open(os.path.join(ref_dir, "certref_n1000_S-5.json")))
    with mp.workdps(40):
        worst_d = 9999.0
        for i in (1, 2, 500, 998, 999):
            truth = mpf(1000) / (mpf(i) * (1000 - i)) \
                * (1 - hyp1f1(1000 - i, 1000, mpf(5))) / (-expm1(mpf(5)))
            got = mpf(ref5["E"][i - 1])
            if got != truth:
                worst_d = min(worst_d, float(-log10(fabs((got - truth) / truth))))
    assert worst_d >= 15, f"shipped S=-5 reference vs 1F1 only {worst_d:.1f}d"
    BOUND = 1e-4
    worstA = worstB = 0.0
    worst_repro = 0.0
    for cell in tab["cells"]:
        S = cell["S"]
        ref = np.array([float(v) for v in json.load(open(os.path.join(
            ref_dir, f"certref_n1000_S{S}.json")))["E"]])
        EA = eng.equilibrium_sfs(float(S))
        EB, meta = eng.expected_sfs(float(S), [(1.0, 0.5)])
        assert meta["negative_entries"] == 0, \
            f"S={S}: negative_entries={meta['negative_entries']}"
        for b, (lo, hi) in bands.items():
            a = float(np.max(np.abs((EA - ref)[lo:hi] / ref[lo:hi])))
            h = float(np.max(np.abs((EB - ref)[lo:hi] / ref[lo:hi])))
            assert a <= BOUND and h <= BOUND, \
                f"S={S} band {b}: A={a:.3e} B={h:.3e} > {BOUND:.0e}"
            worstA, worstB = max(worstA, a), max(worstB, h)
            for got, want in ((a, cell["A_eqproj"][b]), (h, cell["B_held"][b])):
                if want > 1e-6:                    # quadrature-dominated bands
                    worst_repro = max(worst_repro, abs(got / want - 1.0))
                else:                              # mid classes stay tiny
                    assert got <= 2e-6, f"S={S} band {b}: {got:.2e} (table {want:.2e})"
    assert worst_repro <= 1e-3, \
        f"band table not reproduced: worst rel diff {worst_repro:.1e}"
    return (f"n=1000, 5 S x 7 bands vs certified refs: proj {worstA:.2e}, "
            f"held {worstB:.2e} (<= 1e-4); table repro {worst_repro:.0e}; "
            f"ref vs 1F1 {worst_d:.0f}d")


@leg("transient_selfconv")
def leg_transient_selfconv():
    if not (_has("numpy") and _has("scipy")):
        return "SKIP:numpy/scipy absent — pip install numpy scipy"
    import json
    import numpy as np
    from popcorn.transient import TransientSFSEngine, project_down, fold_sfs
    ref_dir = os.path.join(HERE, "reference", "transient")
    sc = json.load(open(os.path.join(ref_dir, "selfconv_n1000.json")))
    eng = TransientSFSEngine(n=1000)
    i = np.arange(1, 1000)
    E0 = eng.equilibrium_sfs(0.0)                       # neutral: 1/i
    r1 = float(np.max(np.abs(E0 * i - 1.0)))
    E2, _ = eng.expected_sfs(0.0, [(2.0, 30.0)], dt0=2e-3)   # long hold: 2/i
    r2 = float(np.max(np.abs(E2 * i / 2.0 - 1.0)))
    E3a = eng.equilibrium_sfs(-5.0)                     # stationarity under CN
    E3b, _ = eng.expected_sfs(-5.0, [(1.0, 0.5)])
    r3 = float(np.max(np.abs(E3b / E3a - 1.0)))
    assert r1 < 2e-4, f"neutral eq vs 1/i: {r1:.2e}"
    assert r2 < 5e-3, f"nu=2 hold vs 2/i: {r2:.2e}"
    assert r3 < 5e-5, f"S=-5 hold-still: {r3:.2e}"
    # dt and grid refinement at S=-5 on the shipped two-epoch history
    hist = [tuple(e) for e in sc["hist"]]
    row = [r for r in sc["rows"] if r["S"] == -5.0][0]
    Ec, meta = eng.expected_sfs(-5.0, hist, dt0=4e-4)
    Ef, _ = eng.expected_sfs(-5.0, hist, dt0=1e-4)
    assert meta["negative_entries"] == 0
    m = Ec > sc["mass_filter_abs"]
    dt_rel = float(np.max(np.abs(Ef[m] / Ec[m] - 1.0)))
    assert dt_rel <= 5e-4, f"dt refinement 4e-4->1e-4: {dt_rel:.2e}"
    assert abs(dt_rel / row["dt_maxrel_massfilt"] - 1.0) <= 1e-2, \
        f"dt selfconv {dt_rel:.3e} vs table {row['dt_maxrel_massfilt']:.3e}"
    eng_f = TransientSFSEngine(n=1000, grid_kw=dict(per_decade=224, dx_mid=5.6e-4))
    assert eng_f.G == sc["G_fine"], f"fine grid G={eng_f.G} != {sc['G_fine']}"
    Eg, _ = eng_f.expected_sfs(-5.0, hist, dt0=4e-4)
    gr_rel = float(np.max(np.abs(Eg[m] / Ec[m] - 1.0)))
    assert gr_rel <= 1e-4, f"grid refinement: {gr_rel:.2e}"
    assert abs(gr_rel / row["grid_maxrel_massfilt"] - 1.0) <= 1e-2, \
        f"grid selfconv {gr_rel:.3e} vs table {row['grid_maxrel_massfilt']:.3e}"
    # exact utilities: hypergeometric down-sampling keeps 1/i; folding keeps mass
    k = np.arange(1, 100)
    pd = project_down(1.0 / i, 1000, 100)
    r4 = float(np.max(np.abs(pd * k - 1.0)))
    assert r4 < 1e-10, f"project_down(1/i) vs 1/k: {r4:.1e}"
    F = fold_sfs(Ec)
    assert len(F) == 500 and abs(F.sum() / Ec.sum() - 1.0) < 1e-13
    return (f"1/i {r1:.1e}, 2/i {r2:.1e}, hold {r3:.1e}; S=-5 two-epoch "
            f"dt-ref {dt_rel:.2e}, grid-ref {gr_rel:.2e} (= table); "
            f"down-proj {r4:.0e}")


@leg("ancestral_join")
def leg_ancestral_join():
    import gzip
    import io
    import json
    from popcorn.ancestral import (join_sites, polarization_rates, load_fasta,
                                   classify, orient, epo_join)
    refd = os.path.join(HERE, "reference", "ancestral")
    exp = json.load(open(os.path.join(refd, "epo_mini_expected.json")))
    fasta = os.path.join(refd, exp["fasta"])
    sites = os.path.join(refd, exp["sites"])
    hdrs, seq = load_fasta(fasta)
    assert seq == "ACgtN-.aTGCA" and len(hdrs) == 2, \
        f"load_fasta: seq={seq!r} headers={hdrs!r}"
    assert [classify(c) for c in "ACGTacgtN-.n"] == \
        ["high"] * 4 + ["low"] * 4 + ["none"] * 4, "classify convention"
    assert orient("A", "G", "A") == ("anc_eq_ref", "A", "G")
    assert orient("T", "A", "t") == ("anc_eq_ref", "T", "A")
    assert orient("T", "C", "C") == ("anc_eq_alt", "C", "T")
    assert orient("A", "C", "g") == ("anc_neither", None, None)
    assert orient("A", "G", "N") == ("unpolarizable", None, None)
    assert orient("AT", "A", "a") == (None, None, None)
    scratch = tempfile.mkdtemp(prefix="popcorn_anc_")
    try:
        out = os.path.join(scratch, "mini.anc.tsv.gz")
        counts = join_sites(fasta, sites, out, chrom=exp["chrom"])
        with gzip.open(out, "rt") as fh:
            assert fh.readline().split() == ["pos", "ref", "alt", "anc", "conf"]
            rows = [[int(f[0])] + f[1:] for f in
                    (l.rstrip("\n").split("\t") for l in fh)]
        assert rows == exp["expected_rows"], f"rows {rows}"
        for k, v in exp["expected_counts"].items():
            assert counts[k] == v, f"counts[{k}] {counts[k]!r} != {v!r}"
        assert counts["skipped_other_chrom"] == 1 and counts["layout"] == "chrom", \
            f"skipped_other_chrom {counts['skipped_other_chrom']} layout {counts['layout']}"
        rates = polarization_rates(counts)
        for k, v in exp["expected_rates"].items():
            assert rates[k] == v, f"rates[{k}] {rates[k]!r} != {v!r}"
        # 'pos ref alt' layout (header row, no chromosome column) and the
        # in-memory tuple form must give the same table and register.
        kept = [l.split("\t") for l in open(sites).read().splitlines()
                if not l.startswith("#") and l.split("\t")[0] in ("chr21", "21")]
        buf = io.StringIO("pos\tref\talt\textra\n" + "".join(
            f"{p}\t{r}\t{a}\tx\n" for _c, p, r, a in kept))
        sink = io.StringIO()
        c2 = join_sites(fasta, buf, sink)
        rows2 = [[int(f[0])] + f[1:] for f in
                 (l.split("\t") for l in sink.getvalue().splitlines()[1:])]
        c3 = join_sites(fasta, [(c, int(p), r, a) for c, p, r, a in kept]
                        + [("chr22", 9, "T", "G")], None, chrom="21")
        assert rows2 == exp["expected_rows"] and c2["layout"] == "nochrom", "nochrom layout"
        for k in ("total_sites", "out_of_fasta_range", "conf_counts_all_sites",
                  "biallelic_snv", "mismatch_class_biallelic_snv_anc_neither"):
            assert c2[k] == counts[k] == c3[k], f"layout disagreement on {k}"
        assert c3["skipped_other_chrom"] == 1, "tuple chrom filter"
    finally:
        shutil.rmtree(scratch, ignore_errors=True)
    b = counts["biallelic_snv"]
    return (f"12-bp fixture: {len(rows)} rows, conf "
            f"{tuple(counts['conf_counts_all_sites'].values())}, biallelic "
            f"{b['n_biallelic_snv']} neither "
            f"{counts['mismatch_class_biallelic_snv_anc_neither']} unpol "
            f"{b['unpolarizable']}; 3 input forms agree")


@leg("ancestral_cli")
def leg_ancestral_cli():
    import gzip
    import json
    refd = os.path.join(HERE, "reference", "ancestral")
    exp = json.load(open(os.path.join(refd, "epo_mini_expected.json")))
    scratch = tempfile.mkdtemp(prefix="popcorn_anc_cli_")
    try:
        out = os.path.join(scratch, "mini.anc.tsv.gz")
        summ = os.path.join(scratch, "mini.summary.json")
        r = subprocess.run([PY, os.path.join(HERE, "epo_join.py"),
                            "--fasta", os.path.join(refd, exp["fasta"]),
                            "--sites", os.path.join(refd, exp["sites"]),
                            "--out", out, "--chrom", exp["chrom"],
                            "--summary", summ],
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                           text=True, timeout=60)
        assert r.returncode == 0 and r.stdout.startswith("epo_join chr21: 9 sites"), \
            f"rc={r.returncode}\n{r.stdout[-800:]}"
        doc = json.load(open(summ))
        for k, v in exp["expected_counts"].items():
            assert doc[k] == v, f"summary[{k}] {doc[k]!r} != {v!r}"
        for k, v in exp["expected_rates"].items():
            assert doc["rates"][k] == v, f"summary rates[{k}]"
        with gzip.open(out, "rt") as fh:
            fh.readline()
            rows = [[int(f[0])] + f[1:] for f in
                    (l.rstrip("\n").split("\t") for l in fh)]
        assert rows == exp["expected_rows"], f"rows {rows}"
        # stdin form ('--sites -'), as in  cut -f1,2,4,5 body.vcf | epo_join.py
        out2 = os.path.join(scratch, "mini2.anc.tsv")
        r2 = subprocess.run([PY, os.path.join(HERE, "epo_join.py"),
                             "--fasta", os.path.join(refd, exp["fasta"]),
                             "--sites", "-", "--out", out2, "--chrom", "21"],
                            input=open(os.path.join(refd, exp["sites"])).read(),
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            text=True, timeout=60)
        assert r2.returncode == 0, f"rc={r2.returncode}\n{r2.stdout[-800:]}"
        rows2 = [[int(f[0])] + f[1:] for f in
                 (l.split("\t") for l in open(out2).read().splitlines()[1:])]
        assert rows2 == exp["expected_rows"], "stdin form rows"
    finally:
        shutil.rmtree(scratch, ignore_errors=True)
    return "CLI file + stdin forms: table and --summary JSON equal the reference"


@leg("twosfs_identities")
def leg_twosfs_identities():
    import json
    from fractions import Fraction
    from popcorn.twosfs import identity_checks, lambda_moments, kingman_rate
    from popcorn.lambda_coalescent import expected_lengths, beta_rate
    ref = json.load(open(os.path.join(HERE, "reference", "twosfs",
                                      "engine_checks.json")))["checks"]
    bad_ref = [k for k, v in ref.items() if v is not True]
    assert not bad_ref, f"reference record not all true: {bad_ref[:5]}"
    for n in (3, 4, 5, 6):
        c = identity_checks(n)
        bad = [k for k, v in c.items() if v is not True]
        assert not bad, f"n={n}: {bad}"
        for k in c:
            pref = k.split("_")[0]
            assert any(r.startswith(pref) for r in ref), f"{k}: no reference record"
    u, v = lambda_moments(3, kingman_rate())
    assert (u, v) == ({1: Fraction(2), 2: Fraction(1)},
                      {(1, 1): Fraction(6), (1, 2): Fraction(3), (2, 2): Fraction(2)})
    u5, _ = lambda_moments(5, beta_rate(Fraction(3, 2)))
    h5 = expected_lengths(5, beta_rate(Fraction(3, 2)))[5]
    assert all(u5[i] == h5[i] for i in range(1, 5)), "E[L_i] != lambda_coalescent at n=5"
    return ("checks 1-7 exact True at n=3..6 (7 families vs expected_lengths, "
            "set-partition route, Kingman 2nd route, limits, semigroup, leaf "
            f"conservation); record {len(ref)} true; n=3 hand values")


@leg("twosfs_hull_ref")
def leg_twosfs_hull_ref():
    import json
    from fractions import Fraction
    from popcorn.twosfs import (lambda_moments, component_moments, normalize_matrix,
                                project_hull, rate_function, pairs)
    doc = json.load(open(os.path.join(HERE, "reference", "twosfs", "support_hull.json")))
    rows = doc["rows"]
    assert len(rows) == 18 and sum(r["family"] == "kingman" for r in rows) == 1
    npos = 0
    for r in rows:
        n, key = r["n"], (r["n"], r["family"], r["param"])
        u, v = lambda_moments(n, rate_function(r["family"], r["param"]))
        totu = sum(u.values())
        xi = [u[i] / totu for i in range(1, n)]
        assert [str(x) for x in xi] == r["xi_hat"], f"{key}: xi_hat"
        q, _tot = normalize_matrix(v, n)
        assert all(str(q[ij]) == r["Q_lambda"][f"{ij[0]},{ij[1]}"] for ij in pairs(n)), \
            f"{key}: Q_lambda"
        w = [Fraction(x) for x in r["sfs_mixture_weights"]]
        assert sum(w) == 1 and all(x >= 0 for x in w) and len(w) == len(r["support_t"])
        comps, mix = [], [Fraction(0)] * (n - 1)
        for wk, t in zip(w, r["support_t"]):
            m1, m2 = component_moments(n, Fraction(t))
            s1 = sum(m1.values())
            mix = [mix[i - 1] + wk * m1[i] / s1 for i in range(1, n)]
            comps.append(normalize_matrix(m2, n)[0])
        assert mix == xi, f"{key}: 1-SFS is not the stated exact mixture of atoms"
        best = project_hull(n, q, comps)
        assert best is not None and str(best["d2"]) == r["dist2_sq_exact"], f"{key}: dist2"
        assert best["d2"] > 0
        npos += 1
        assert list(best["S"]) == r["proj_active_set"] and \
            [str(x) for x in best["lam"]] == r["proj_weights"], f"{key}: active set/weights"
        winf = max(abs(x) for x in best["w"].values())
        assert str(best["d2"] / winf) == r["delta_l1_lb_exact"], f"{key}: l1 bound"
    wl = doc["within_lambda_n3"]
    qa, _ = normalize_matrix(lambda_moments(3, rate_function("beta", "1/2"))[1], 3)
    qb, _ = normalize_matrix(lambda_moments(3, rate_function("dirac", "3/4"))[1], 3)
    assert qa == qb and wl["equal"] is True and \
        all(str(qa[ij]) == wl["Q_beta_1/2"][f"{ij[0]},{ij[1]}"] for ij in pairs(3))
    r4 = [r for r in rows if r["n"] == 4 and r["family"] == "beta" and r["param"] == "1/2"][0]
    assert r4["dist2_sq_exact"] == "347711/15859712"
    return (f"18 rows n=3..6: xi_hat, exact 1-SFS atom-mixture certificates, Q, "
            f"hull projections (d2>0 in {npos}, n4 Beta(1/2) d2=347711/15859712) "
            "string-equal; n=3 Beta(1/2)==Dirac(3/4) 2-SFS")


@leg("twosfs_kernel")
def leg_twosfs_kernel():
    if not _has("sympy"):
        return "SKIP:sympy absent — pip install sympy"
    import json
    import sympy as sp
    from fractions import Fraction
    from popcorn.twosfs import (kingman_kernel, kernel_constant_size_check,
                                twotime_kernel_poly, component_moments, pairs)
    ref = json.load(open(os.path.join(HERE, "reference", "twosfs",
                                      "engine_checks.json")))["checks"]
    for n in (3, 4, 5, 6):
        kern, m2p, s_, u_ = kingman_kernel(n)          # asserts diag identity
        assert len(kern) == len(pairs(n)) == n * (n - 1) // 2
        assert kernel_constant_size_check(n) is True, f"n={n}: dx/x reconstruction"
        assert ref[f"kernel_constant_size_reconstruction_n{n}"] is True
        assert ref[f"kernel_diagonal_is_m2_n{n}"] is True
    # kernel at (s, u=1) equals component_moments m2 at a rational level, n=5
    kern, m2p, s_, u_ = kingman_kernel(5)
    t = Fraction(373, 512)
    _m1, m2 = component_moments(5, t)
    for ij in pairs(5):
        val = kern[ij].subs({s_: sp.Rational(373, 512), u_: 1})
        assert Fraction(int(sp.numer(val)), int(sp.denom(val))) == m2[ij], f"n=5 {ij}"
    k2, (s2, u2) = twotime_kernel_poly(3)
    assert sp.expand(k2[(1, 1)].subs({s2: 1, u2: 1}) - 9) == 0   # n=3 star: N_1^2 = 9
    return ("n=3..6: Csym(s,1) == m2 polys (symbolic 0); int int Csym dx/x dy/y == "
            "constant-size Kingman E[L_i L_j] exact; Csym(373/512,1) == m2(373/512) at n=5")


@leg("twosfs_class_certs")
def leg_twosfs_class_certs():
    if not _has("sympy"):
        return "SKIP:sympy absent — pip install sympy"
    import json
    import sympy as sp
    from fractions import Fraction as Fr
    from popcorn.twosfs import (certify_witness, load_reference, rate_function,
                                lambda_moments, pairs, bern_1d_nonneg, bern_2d_nonneg,
                                poly_coeffs_2d, kingman_kernel, solve_exact)
    from popcorn import certificates
    R = certificates.region
    ref = load_reference(os.path.join(HERE, "reference", "twosfs", "class_certificates.json"))
    rows = ref["rows"]
    assert len(rows) == 19
    nfull = nfresh = 0
    for r in rows:
        n, key = r["n"], (r["n"], r["family"], r["param"])
        lam = rate_function(r["family"], r["param"])
        W = {tuple(int(x) for x in k.split(",")): Fr(v) for k, v in r["witness_W"].items()}
        assert set(W) == set(pairs(n)) and max(abs(x) for x in W.values()) == 1, f"{key}: W"
        _u, v = lambda_moments(n, lam)
        GM = sum(W[ij] * v[ij] for ij in pairs(n))
        tot = sum((1 if i == j else 2) * v[(i, j)] for i, j in pairs(n))
        assert GM < 0 and str(GM) == r["G_M_lambda_exact"], f"{key}: <W,M>"
        assert str(-GM / tot) == r["norm_margin_exact"] == r["delta_l1_lb_exact"], f"{key}: margin"
        if r["family"] == "kingman":
            assert r["verdict"] == "OUT_ATOMPOOLINGS", "control must never be OUT_FULLCLASS"
        elif n >= 4:
            assert r["verdict"] == "OUT_FULLCLASS", f"{key}: {r['verdict']}"
            nfull += 1
        if n <= 5:                                   # fresh exact certificates
            c = certify_witness(n, lam, r["witness_W"])
            assert c["verdict"] == r["verdict"], f"{key}: fresh {c['verdict']} != {r['verdict']}"
            assert str(c["G_M_exact"]) == r["G_M_lambda_exact"]
            if r["verdict"] == "OUT_FULLCLASS":
                assert c["full_class"] is True
            else:
                assert c["full_class"] is not True and c["atom_poolings"] is True
            nfresh += 1
    assert nfull == 13
    # planted MUST-FAIL / must-bisect controls (three-valued contract)
    kern, m2p, s_, u_ = kingman_kernel(4)
    assert bern_2d_nonneg(s_ * u_ - sp.Rational(1, 64), s_, u_) is False
    assert bern_1d_nonneg({0: Fr(-1, 64), 1: Fr(1)}) is False
    assert bern_1d_nonneg({0: Fr(1, 4) + Fr(1, 64), 1: Fr(-1), 2: Fr(1)}) is True
    assert bern_2d_nonneg((s_ - sp.Rational(1, 2)) ** 2 + sp.Rational(1, 64), s_, u_) is True
    # twin: the same K_W through popcorn.certificates.region certifies too
    r4 = [r for r in rows if r["n"] == 4 and r["family"] == "beta" and r["param"] == "1/2"][0]
    W4 = {tuple(int(x) for x in k.split(",")): Fr(v) for k, v in r4["witness_W"].items()}
    expr = sum(sp.Rational(W4[ij].numerator, W4[ij].denominator) * kern[ij] for ij in pairs(4))
    cf = poly_coeffs_2d(sp.expand(expr), s_, u_)
    ds, du = max(a for a, _ in cf), max(b for _, b in cf)
    assert R.bern_nonneg(R.to_bernstein(cf, ds, du), 7) is True, "region twin disagrees"
    neg = poly_coeffs_2d(s_ * u_ - sp.Rational(1, 64), s_, u_)
    assert R.bern_nonneg(R.to_bernstein(neg, 1, 1), 7) is False, "region twin accepted a negative"
    # square-solve contract vs the rectangular region.solve_exact
    A, b = [[Fr(1), Fr(1)], [Fr(2), Fr(2)]], [Fr(1), Fr(2)]
    assert solve_exact(A, b) is None and R.solve_exact(A, b) is not None
    assert solve_exact([[Fr(2), Fr(1)], [Fr(1), Fr(3)]], [Fr(3), Fr(5)]) == [Fr(4, 5), Fr(7, 5)]
    return (f"19 rows: exact <W,M><0 and margins string-equal; {nfresh} fresh certificates "
            f"(n<=5) reproduce the verdicts ({nfull} OUT_FULLCLASS at n>=4, controls "
            "atom-poolings only); planted negatives refuted; region twin agrees")


@leg("twosfs_verdict")
def leg_twosfs_verdict():
    if not (_has("sympy") and _has("numpy") and _has("scipy")):
        return "SKIP:sympy/numpy/scipy absent — pip install sympy numpy scipy"
    from fractions import Fraction as Fr
    from popcorn.twosfs import (class_verdict, load_reference, beta_rate, kingman_rate,
                                kingman_kernel, component_moments, hunt)
    ref = {(r["n"], r["family"], r["param"]): r for r in load_reference()["rows"]}
    rb = class_verdict(4, beta_rate(Fr(1, 2)))
    assert rb["verdict"] == "OUT_FULLCLASS", f"n=4 Beta(1/2): {rb['verdict']}"
    m = Fr(rb["norm_margin_exact"])
    assert m > 0 and Fr(rb["G_M_lambda_exact"]) < 0 and rb["bernstein"] == "certified"
    same = rb["witness_W"] == ref[(4, "beta", "1/2")]["witness_W"]
    rk = class_verdict(4, kingman_rate())
    assert rk["verdict"] != "OUT_FULLCLASS", "constant-size Kingman certified OUT of its own class"
    assert rk["verdict"] == "OUT_ATOMPOOLINGS", f"n=4 Kingman control: {rk['verdict']}"
    r3 = class_verdict(3, beta_rate(Fr(1, 2)))
    assert r3["verdict"] != "OUT_FULLCLASS" and \
        r3.get("full_class_search") == "no_certified_witness"
    kern, m2p, s_, u_ = kingman_kernel(4)               # a class member: one atom at t=1/2
    _m1, m2 = component_moments(4, Fr(1, 2))
    ra = hunt(4, m2, kern, m2p, s_, u_)
    assert not ra["verdict"].startswith("OUT"), f"atom certified out: {ra['verdict']}"
    assert ra["verdict"] == "INSIDE_2SFS" and "1/2" in ra["inside_support"], f"atom: {ra}"
    return (f"live: n=4 Beta(1/2) OUT_FULLCLASS margin {float(m):.4e} (witness "
            f"{'==' if same else '!='} reference), Kingman control OUT_ATOMPOOLINGS only, "
            f"n=3 no full-class witness, atom t=1/2 INSIDE_2SFS exactly")


@leg("twosfs_montecarlo")
def leg_twosfs_montecarlo():
    import json
    from popcorn.twosfs import lambda_moments, normalize_matrix, rate_function
    doc = json.load(open(os.path.join(HERE, "reference", "twosfs", "msprime_moments.json")))
    worst, nstat = 0.0, 0
    for b in doc["batteries"]:
        n = b["n"]
        u, v = lambda_moments(n, rate_function(b["family"], b["param"]))
        q, _ = normalize_matrix(v, n)
        z = {}
        if "raw" in b["mode"]:
            for i in range(1, n):
                z[f"EL_{i}"] = (b["first"][i - 1] - float(u[i])) / b["first_sd"][i - 1]
                for j in range(i, n):
                    z[f"ELL_{i}{j}"] = (b["second"][i - 1][j - 1] - float(v[(i, j)])) \
                        / b["second_sd"][i - 1][j - 1]
        for i in range(1, n):
            for j in range(i, n):
                z[f"q_{i}{j}"] = (b["q"][i - 1][j - 1] - float(q[(i, j)])) / b["q_sd"][i - 1][j - 1]
        assert set(z) == set(b["z"]), f"{b['family']} n={n}: statistic set"
        for k, val in z.items():
            assert abs(round(val, 3) - b["z"][k]) <= 0.0011, \
                f"{b['family']} n={n} {k}: {val:.4f} vs {b['z'][k]}"
            worst = max(worst, abs(val))
            nstat += 1
    assert nstat == 72 and worst < doc["acceptance_abs_z"] and \
        abs(worst - doc["worst_abs_z"]) < 2e-3
    return (f"{nstat} z-scores vs stored simulator moments (Kingman n=4,6 raw+norm; Beta(3/2) "
            f"n=4,5; Dirac-mixture n=4): equal the shipped table, worst |z| = {worst:.3f} < 4")


@leg("foldgate_map")
def leg_foldgate_map():
    import random
    from fractions import Fraction
    from popcorn.foldgate import (fold_coeffs, apply_fold, pullback, rank_report,
                                  PAIRS, FPAIRS, make_object, load_reference_rows)
    ref = load_reference_rows(os.path.join(HERE, "reference", "foldgate",
                                           "folded_reference_n456.json"))
    # the n = 4 map written out: classes C_1 = {1,3}, C_2 = {2}
    assert fold_coeffs(4) == {(1, 1): {(1, 1): 1, (1, 3): 2, (3, 3): 1},
                              (1, 2): {(1, 2): 1, (2, 3): 1},
                              (2, 2): {(2, 2): 1}}, f"fold_coeffs(4) = {fold_coeffs(4)}"
    # every unfolded cell of the full table lands in exactly one folded cell:
    # total integer weight over stored coords, counted with multiplicity, is (n-1)^2
    for n in range(3, 11):
        Phi = fold_coeffs(n)
        assert sorted(Phi) == FPAIRS(n), f"n={n}: folded coordinate set"
        tot = sum((1 if a == b else 2) * c for (a, b), co in Phi.items() for c in co.values())
        assert tot == (n - 1) ** 2, f"n={n}: total weight {tot} != {(n - 1) ** 2}"
    # exact ranks equal the reference table (Phi is onto; kernel 3, 7, 9)
    assert [r["n"] for r in ref["ranks"]] == [4, 5, 6]
    for r in ref["ranks"]:
        assert rank_report(r["n"]) == r, f"rank_report({r['n']}) = {rank_report(r['n'])} != {r}"
    ranks = ",".join(f"{r['rank_sym']}/{r['kernel_dim_sym']}" for r in ref["ranks"])
    # pullback identity <v, Phi x> == <Phi^T v, x> exactly on random rationals
    rng = random.Random(4096)
    cnt = 0
    for n in (4, 5, 6, 7):
        for _ in range(40):
            x = {ij: Fraction(rng.randrange(-999, 1000), rng.randrange(1, 97)) for ij in PAIRS(n)}
            v = {fp: Fraction(rng.randrange(-999, 1000), rng.randrange(1, 97)) for fp in FPAIRS(n)}
            lhs = sum(v[fp] * val for fp, val in apply_fold(n, x).items())
            W = pullback(n, v)
            cnt += (lhs == sum(W[ij] * x[ij] for ij in PAIRS(n)))
    assert cnt == 160, f"pullback identity held on {cnt}/160 vectors"
    # mass conservation and the shipped folded spectra, all 16 targets
    assert len(ref["rows"]) == 16
    for row in ref["rows"]:
        n = row["n"]
        obj = make_object(n, row["M"])
        assert sum((1 if i == j else 2) * q for (i, j), q in obj["q"].items()) == 1
        mass = sum((1 if a == b else 2) * obj["phi"][(a, b)] for a, b in FPAIRS(n))
        assert mass == 1, f"n{n} {row['family']} {row['param']}: folded mass {mass}"
        assert {f"{a},{b}": str(obj["phi"][(a, b)]) for a, b in FPAIRS(n)} == row["phi_folded"], \
            f"n{n} {row['family']} {row['param']}: folded spectrum != reference"
    return (f"n=4 map literal; onto with rank/kernel n4,5,6 = {ranks}; pullback identity "
            f"160/160 exact; 16 folded spectra string-equal, mass 1")


@leg("foldgate_annihilated")
def leg_foldgate_annihilated():
    from fractions import Fraction
    from popcorn.foldgate import (FPAIRS, make_object, load_reference_rows,
                                  load_ingredients, component_q, apply_fold, tgrid_of)
    refd = os.path.join(HERE, "reference", "foldgate")
    ref = load_reference_rows(os.path.join(refd, "folded_reference_n456.json"))
    ings = {n: load_ingredients(n, os.path.join(refd, "kingman_class_polys_n456.json"))
            for n in (4, 5, 6)}
    tg = set(tgrid_of())
    rows = [r for r in ref["rows"] if r["verdict"] == "ANNIHILATED"]
    assert len(rows) == 10 and all(r["inside_lp"] == "exact" for r in rows), len(rows)
    assert sorted((r["n"], r["family"]) for r in rows if r["family"] == "kingman") == \
        [(4, "kingman"), (5, "kingman")], "Kingman sanity rows at n=4,5 must be ANNIHILATED"
    first = None
    for row in rows:
        n = row["n"]
        phi = make_object(n, row["M"])["phi"]
        sw = [(Fraction(t), Fraction(w)) for t, w in row["inside_support_weights"]]
        assert len(sw) == 3 and all(w > 0 for _, w in sw) and sum(w for _, w in sw) == 1, \
            f"n{n} {row['family']} {row['param']}: weights {sw}"
        assert all(t in tg for t, _ in sw), "support nodes must lie on tgrid_of()"
        fcs = [apply_fold(n, component_q(n, ings[n], t)) for t, _ in sw]
        acc = {fp: sum(w * fc[fp] for (_, w), fc in zip(sw, fcs)) for fp in FPAIRS(n)}
        assert acc == phi, f"n{n} {row['family']} {row['param']}: convex identity fails"
        if first is None:
            first = (n, row, phi, sw, fcs)
    # planted MUST-FAIL: perturbing two weights (sum still 1) breaks the exact identity
    n, row, phi, sw, fcs = first
    eps = Fraction(1, 10 ** 12)
    bad = [(sw[0][0], sw[0][1] + eps), (sw[1][0], sw[1][1] - eps), sw[2]]
    acc = {fp: sum(w * fc[fp] for (_, w), fc in zip(bad, fcs)) for fp in FPAIRS(n)}
    assert sum(w for _, w in bad) == 1 and acc != phi, "planted perturbation not detected"
    # component_q at t = 0 is the directional limit: equals the t -> 0+ trend of t = 2^-15
    q0, q1 = component_q(6, ings[6], 0), component_q(6, ings[6], Fraction(1, 1 << 15))
    d = max(abs(q0[ij] - q1[ij]) for ij in q0)
    assert d < Fraction(1, 1000), f"t=0 limit vs t=2^-15: {float(d):.2e}"
    return ("10 rows: phi == sum_k lambda_k Phi(q(t_k)) exact from shipped polynomials, "
            "lambda > 0, sum 1, 3 nodes each (Kingman n=4,5 included); planted 1e-12 "
            "perturbation breaks it")


@leg("foldgate_obstructions")
def leg_foldgate_obstructions():
    from fractions import Fraction
    from popcorn.foldgate import (PAIRS, FPAIRS, fold_coeffs, apply_fold, pullback,
                                  make_object, load_reference_rows, load_ingredients,
                                  bern_1d_nonneg)
    refd = os.path.join(HERE, "reference", "foldgate")
    ref = load_reference_rows(os.path.join(refd, "folded_reference_n456.json"))
    ings = {n: load_ingredients(n, os.path.join(refd, "kingman_class_polys_n456.json"))
            for n in (4, 6)}
    rows = [r for r in ref["rows"] if r["verdict"] == "NO_LINEAR_WITNESS"]
    assert len(rows) == 6 and all(r["inside_lp"] == "lp_infeasible" for r in rows), len(rows)
    assert sorted(r["n"] for r in rows) == [4, 6, 6, 6, 6, 6]
    margins = []
    for row in rows:
        n = row["n"]
        tag = f"n{n} {row['family']} {row['param']}"
        obj = make_object(n, row["M"])
        phiM = apply_fold(n, obj["M"])
        # (a) exact conic certificate: Phi(M) == sum_k w_k Phi(Csym(x_k, y_k/x_k)), w_k >= 0
        om = row["outer_membership"]
        ws = [Fraction(w) for w in om["weights"]]
        assert len(ws) == len(om["generators"]) == len(FPAIRS(n)) and all(w >= 0 for w in ws), tag
        acc = {fp: Fraction(0) for fp in FPAIRS(n)}
        for g, w in zip(om["generators"], ws):
            xs, ys = g.strip("()").split(",")
            x, y = Fraction(xs), Fraction(ys)
            u = y / x if x != 0 else Fraction(0)
            col = {ij: sum(c * x ** es * u ** eu for (es, eu), c in p.items())
                   for ij, p in ings[n]["twotime_kernel"].items()}
            fcol = apply_fold(n, col)
            for fp in FPAIRS(n):
                acc[fp] += w * fcol[fp]
        assert acc == phiM, f"{tag}: conic identity fails"
        # (b) exact diagonal witness: margin string-equal and < 0, pullback cross-route,
        #     and a FRESH Bernstein certificate of p_v(t) = sum_ab v_ab (Phi m2)_ab(t) >= 0 on [0,1]
        dw = row["diagonal_witness"]
        v = {tuple(map(int, k.split(","))): Fraction(x) for k, x in dw["witness_v_folded"].items()}
        G = sum(v[fp] * phiM[fp] for fp in FPAIRS(n))
        assert str(G) == dw["G_phiM_exact"] and G < 0, f"{tag}: margin {G}"
        W = pullback(n, v)
        assert sum(W[ij] * obj["M"][ij] for ij in PAIRS(n)) == G, f"{tag}: pullback route"
        assert -G / obj["tot"] == Fraction(dw["norm_margin_exact"]), f"{tag}: normalized margin"
        assert Fraction(dw["standoff_l1_lb_exact"]) == \
            -G / obj["tot"] / max(abs(x) for x in v.values())
        pv = {}
        for fp, co in fold_coeffs(n).items():
            for ij, c in co.items():
                for e, ce in ings[n]["m2_diag"][ij].items():
                    pv[e] = pv.get(e, Fraction(0)) + v[fp] * c * ce
        assert bern_1d_nonneg(pv) is True, f"{tag}: diagonal polynomial not certified"
        margins.append(float(Fraction(dw["norm_margin_exact"])))
    # planted MUST-FAIL controls: a polynomial negative inside (0,1) is refused (False);
    # a functional with all-positive coordinates has a positive margin (no witness)
    assert bern_1d_nonneg({0: Fraction(-1, 64), 1: Fraction(1), 2: Fraction(-1)}) is False
    assert bern_1d_nonneg({1: Fraction(1), 2: Fraction(-1)}) is True      # t(1-t) >= 0
    obj4 = make_object(4, rows[0]["M"])
    assert sum(apply_fold(4, obj4["M"]).values()) > 0
    assert ref["realization_n4_beta_1_2"] == {"status": "no_2atom_realization_found_on_scan"}
    return ("6 rows: exact conic certificates re-verified from the shipped kernel; 6 diagonal "
            "witnesses: margin string-equal < 0, pullback route, fresh Bernstein PASS "
            f"(normalized margins {min(margins):.2e}..{max(margins):.2e}); "
            "planted negative poly refused")


@leg("foldgate_deciders")
def leg_foldgate_deciders():
    if not (_has("numpy") and _has("scipy") and _has("sympy")):
        return "SKIP:numpy/scipy/sympy absent — pip install numpy scipy sympy"
    from fractions import Fraction
    import sympy as sp
    from popcorn.foldgate import (PAIRS, FPAIRS, make_object, load_reference_rows,
                                  load_ingredients, folded_components, kernel_sympy,
                                  m2_sympy, fold_exprs, tgrid_of, pullback, apply_fold,
                                  inside_hunt, witness_hunt, outer_membership,
                                  diagonal_witness_hunt, realization_hunt_n4, decide,
                                  bern_2d_nonneg)
    refd = os.path.join(HERE, "reference", "foldgate")
    ref = load_reference_rows(os.path.join(refd, "folded_reference_n456.json"))
    ing = load_ingredients(4, os.path.join(refd, "kingman_class_polys_n456.json"))
    rows = {(r["family"], r["param"]): r for r in ref["rows"] if r["n"] == 4}
    tg = tgrid_of()
    assert len(tg) == 83
    fcomps = folded_components(4, ing, tg)
    kern, (s_, u_) = kernel_sympy(ing)
    kernF = fold_exprs(4, kern)
    s1 = sp.Symbol("s", nonnegative=True)
    m2pF = fold_exprs(4, m2_sympy(ing, s1))
    # (1) inside_hunt re-derives the reference 3-node weights of Beta(1) at n=4 from that hull
    rb = rows[("beta", "1")]
    objb = make_object(4, rb["M"])
    sw = [(Fraction(t), Fraction(w)) for t, w in rb["inside_support_weights"]]
    hull = [fcomps[tg.index(t)] for t, _ in sw]
    sol, status = inside_hunt(4, objb["phi"], hull, [t for t, _ in sw])
    assert status == "exact" and {Fraction(t): Fraction(w) for t, w in sol} == dict(sw), \
        f"inside_hunt on the 3-node hull: {status} {sol}"
    # (2) decide() on the full 83-node component set: Kingman -> ANNIHILATED (sanity law),
    #     Beta(1/2) -> NO_LINEAR_WITNESS with an exact conic certificate re-verified here
    reck = decide(4, make_object(4, rows[("kingman", None)]["M"]), fcomps, tg, kernF, kern, s_, u_)
    assert reck["verdict"] == "ANNIHILATED" and reck["inside_lp"] == "exact", reck["verdict"]
    obj = make_object(4, rows[("beta", "1/2")]["M"])
    rec = decide(4, obj, fcomps, tg, kernF, kern, s_, u_)
    assert rec["verdict"] == "NO_LINEAR_WITNESS" and rec["inside_lp"] == "lp_infeasible", \
        rec["verdict"]
    om = rec["outer_membership"]
    phiM = apply_fold(4, obj["M"])
    acc = {fp: Fraction(0) for fp in FPAIRS(4)}
    for g, w in zip(om["generators"], om["weights"]):
        xs, ys = g.strip("()").split(",")
        x, y = Fraction(xs), Fraction(ys)
        for fp in FPAIRS(4):
            val = kernF[fp].subs({s_: sp.Rational(x), u_: sp.Rational(y / x) if x != 0 else 0})
            acc[fp] += Fraction(w) * Fraction(int(sp.numer(val)), int(sp.denom(val)))
    assert all(Fraction(w) >= 0 for w in om["weights"]) and acc == phiM, "live conic certificate"
    same_gen = om["generators"] == rows[("beta", "1/2")]["outer_membership"]["generators"]
    # (3) witness_hunt alone finds no full-class witness for Beta(1/2) (consistent with (2))
    win, lp0 = witness_hunt(4, obj, kernF, kern, s_, u_)
    assert win is None and lp0 is not None and lp0 <= 1e-8, f"witness_hunt: {win} lp_delta0={lp0}"
    # (4) diagonal_witness_hunt certifies Beta(1/2) outside the single-atom-pooling hull;
    #     its exact margin re-derives through the pullback route
    w = diagonal_witness_hunt(4, obj, m2pF, s1, tg)
    assert isinstance(w, dict) and w["bernstein_1d"] == "certified", f"diagonal witness: {w}"
    v = {tuple(map(int, k.split(","))): Fraction(x) for k, x in w["witness_v_folded"].items()}
    G = sum(pullback(4, v)[ij] * obj["M"][ij] for ij in PAIRS(4))
    assert str(G) == w["G_phiM_exact"] and G < 0, f"margin {G}"
    same_dw = w["G_phiM_exact"] == rows[("beta", "1/2")]["diagonal_witness"]["G_phiM_exact"]
    # (5) the two-atom realization scan at n=4, Beta(1/2) comes back negative, as shipped
    r = realization_hunt_n4(obj, kern, s_, u_, budget_s=20)
    assert r == ref["realization_n4_beta_1_2"], f"realization scan: {r}"
    # (6) bern_2d_nonneg planted controls: exact refusal of a kernel negative at a corner,
    #     certification of a strictly positive one
    assert bern_2d_nonneg(s_ * u_ - sp.Rational(1, 64), s_, u_) is False
    assert bern_2d_nonneg(1 - s_ * u_ / 2, s_, u_) is True
    return (f"n=4 live: inside_hunt re-derives reference weights; decide(): Kingman ANNIHILATED, "
            f"Beta(1/2) NO_LINEAR_WITNESS (conic cert exact, generators "
            f"{'==' if same_gen else '!='} reference), "
            f"no full-class witness (lp_delta0={lp0:.0e}), diagonal witness margin "
            f"{float(-G / obj['tot']):.3e} ({'==' if same_dw else '!='} reference), "
            "realization scan negative")


@leg("twosite_projection")
def leg_twosite_projection():
    import itertools
    import json
    import math
    import random
    from fractions import Fraction
    from popcorn.twosite import (nine_cell, project_grid, brute_grid,
                                 flip_grid, seg_matrix, fold_matrix,
                                 add_into, add_scaled, zeros, cells_key,
                                 GridCache, brute_check, load_geno_tsv,
                                 twosite_projection)
    refd = os.path.join(HERE, "reference", "twosite")
    exp = json.load(open(os.path.join(refd, "brutecheck_expected.json")))
    sites, n_ind, _c = load_geno_tsv(os.path.join(refd, exp["geno"]))
    sites.sort(key=lambda s: s[0])
    dos = [s[3] for s in sites]
    assert (n_ind, len(sites)) == (exp["n_individuals"], exp["n_sites"]) == (12, 24)
    # 1. DP operator vs explicit enumeration of ALL subsets, n = 4, 6, 8.
    r = brute_check(dos, n_ind, (4, 6, 8))
    assert r["EXACT_EQUAL_ALL"] is True and exp["EXACT_EQUAL_ALL"] is True
    for k in ("n4", "n6", "n8"):
        assert r["per_n"][k] == exp["per_n"][k], f"{k}: {r['per_n'][k]} != {exp['per_n'][k]}"
    assert [r["per_n"][k]["cells_checked"] for k in ("n4", "n6", "n8")] == [6900, 13524, 22356]
    # 2. Hand case: N = 2 individuals, dosages (1,2) and (0,1), m = 1:
    #    subset {ind 0} -> (1,0), subset {ind 1} -> (2,1), each 1/2.
    G = project_grid(nine_cell((1, 2), (0, 1)), 2, 1)
    want = [[Fraction(0)] * 3 for _ in range(3)]
    want[1][0] = want[2][1] = Fraction(1, 2)
    assert G == want, f"hand case {G}"
    # 3. Independent oracle: row/column marginals of the pair grid equal the
    #    single-site projection over individuals (trinomial-hypergeometric
    #    on the dosage classes), exactly, for every fixture pair at m = 3.
    def single_site(d, N, m):
        c0, c1, c2 = (sum(1 for x in d if x == v) for v in (0, 1, 2))
        P = [Fraction(0)] * (2 * m + 1)
        for k1 in range(0, min(c1, m) + 1):
            for k2 in range(0, min(c2, m - k1) + 1):
                k0 = m - k1 - k2
                if k0 <= c0:
                    P[k1 + 2 * k2] += Fraction(math.comb(c1, k1) * math.comb(c2, k2)
                                               * math.comb(c0, k0), math.comb(N, m))
        return P
    m = 3
    cache = GridCache()
    acc_pairs, groups, n_marg = zeros(m), {}, 0
    scaled_bad = memo_bad = fold_bad = flip_bad = 0
    for i in range(len(dos)):
        for j in range(i + 1, len(dos)):
            cells = nine_cell(dos[i], dos[j])
            G = project_grid(cells, n_ind, m)
            rows = [sum(row) for row in G]
            cols = [sum(G[a][b] for a in range(2 * m + 1)) for b in range(2 * m + 1)]
            assert rows == single_site(dos[i], n_ind, m), f"row marginal, pair {i},{j}"
            assert cols == single_site(dos[j], n_ind, m), f"col marginal, pair {i},{j}"
            n_marg += 1
            # memo twin == direct; fold/seg mass conservation; flip is an involution
            ent = cache.entry(cells_key(cells), n_ind, m)
            memo_bad += (ent[0] != G) or (cache.fold(ent) != fold_matrix(G))
            memo_bad += cache.umat(ent, True, False) != seg_matrix(flip_grid(G, True, False))
            F, U = fold_matrix(G), seg_matrix(G)
            fold_bad += sum(map(sum, F)) != sum(map(sum, U))
            fold_bad += fold_matrix(flip_grid(G, True, True)) != F
            flip_bad += flip_grid(flip_grid(G, True, True), True, True) != G
            flip_bad += seg_matrix(flip_grid(G, True, False)) != [row for row in U[::-1]]
            add_into(acc_pairs, F)
            k = cells_key(cells)
            groups[k] = groups.get(k, 0) + 1
            if (i + j) % 37 == 0:
                for mult in (1, 3, 7):
                    a1, a2 = zeros(2 * m + 1), zeros(2 * m + 1)
                    add_scaled(a1, G, mult)
                    for _ in range(mult):
                        add_into(a2, G)
                    scaled_bad += a1 != a2
    acc_groups = zeros(m)
    for k, mult in groups.items():
        add_scaled(acc_groups, fold_matrix(project_grid(k, n_ind, m)), mult)
    assert acc_groups == acc_pairs, "sum over pairs != sum over distinct-table groups"
    assert memo_bad == 0 and fold_bad == 0 and flip_bad == 0 and scaled_bad == 0, \
        f"memo {memo_bad} fold {fold_bad} flip {flip_bad} scaled {scaled_bad}"
    assert len(cache) == len(groups) and cache.report()["grid_misses"] == len(groups)
    cache.clear()
    assert len(cache) == 0
    # 4. A larger single pair: N = 40, m = 4 (n = 8), all C(40,4) = 91390
    #    subsets tallied directly (not via brute_grid) vs the DP.
    rng = random.Random(12)
    N, m = 40, 4
    d1 = tuple(rng.choice((0, 0, 1, 1, 2)) for _ in range(N))
    d2 = tuple(rng.choice((0, 1, 0, 2, 0)) for _ in range(N))
    G = project_grid(nine_cell(d1, d2), N, m)
    cnt = {}
    for combo in itertools.combinations(range(N), m):
        key = (sum(d1[x] for x in combo), sum(d2[x] for x in combo))
        cnt[key] = cnt.get(key, 0) + 1
    ncomb = math.comb(N, m)
    assert sum(cnt.values()) == ncomb == 91390
    bad = sum(G[a][b] != Fraction(cnt.get((a, b), 0), ncomb)
              for a in range(2 * m + 1) for b in range(2 * m + 1))
    assert bad == 0 and sum(map(sum, G)) == 1, f"N=40 m=4: {bad} cell mismatches"
    ncells = sum(v['cells_checked'] for v in r['per_n'].values())
    return (f"DP == enumeration on 276 pairs x n=4,6,8 ({ncells} cells, exact); "
            f"marginals == single-site hypergeometric on {n_marg} pairs; "
            f"groups {len(groups)}/276; C(40,4)={ncomb} subsets, 0 mismatches")


@leg("twosite_spectrum")
def leg_twosite_spectrum():
    import json
    from fractions import Fraction
    from popcorn.twosite import (estimate, orient_site, load_geno_tsv,
                                 load_anc, load_map, annotate_cm,
                                 nine_cell, project_grid, fold_matrix,
                                 add_into, zeros, PAIR_CLASSES)
    refd = os.path.join(HERE, "reference", "twosite")
    exp = json.load(open(os.path.join(refd, "spectrum_expected.json")))
    prm = exp["params"]
    geno = os.path.join(refd, exp["geno"])
    anc = os.path.join(refd, exp["anc"])
    gmap = os.path.join(refd, exp["map"])
    kw = dict(gmap=gmap, chrom=prm["chrom"], n=prm["n"], bins=prm["bins"],
              max_bp=prm["max_bp"], mnv_bp=prm["mnv_bp"],
              hotspot_cm_per_mb=prm["hotspot_cm_per_mb"],
              block_mb=prm["block_mb"], fmt=prm["format"])

    def same(got, want, path="", tol=1e-12):
        """exact == on strings/ints/bools/structure; floats to rel tol."""
        if isinstance(want, dict):
            assert isinstance(got, dict) and set(got) == set(want), \
                (f"{path}: keys {sorted(got) if isinstance(got, dict) else type(got)}"
                 f" != {sorted(want)}")
            for k in want:
                same(got[k], want[k], f"{path}.{k}", tol)
        elif isinstance(want, list):
            assert isinstance(got, list) and len(got) == len(want), f"{path}: length"
            for i, (g, w) in enumerate(zip(got, want)):
                same(g, w, f"{path}[{i}]", tol)
        elif isinstance(want, float):
            assert abs(got - want) <= tol * max(1.0, abs(want)), f"{path}: {got!r} != {want!r}"
        else:
            assert got == want, f"{path}: {got!r} != {want!r}"

    doc = estimate(geno, anc=anc, **kw)
    ref = exp["default"]
    assert doc["dist_unit"] == "cM" and doc["n_sites"] == 24
    for k in ("n_individuals", "site_counts", "site_orientation_counts",
              "pair_class_counts", "mnv_control_pairs",
              "signal_pairs_unfolded_eligible", "results"):
        same(json.loads(json.dumps(doc[k])), ref[k], k)
    # dict order of the controls = first-occurrence order of the scan
    for n in ref["results"]:
        for reg in ref["results"][n]:
            assert list(doc["results"][n][reg]["controls"]) == \
                list(ref["results"][n][reg]["controls"])
    cc = doc["pair_class_counts"]
    assert sum(cc.values()) == 276 and all(cc[c] > 0 for c in PAIR_CLASSES), cc
    assert doc["aggregation"]["pairs_grouped"] == \
        276 - cc["na_gap_excluded"] - cc["out_of_bin_range"]
    # high-confidence-only orientation
    doc_hc = estimate(geno, anc=anc, high_conf_only=True, **kw)
    hc = exp["high_conf_only"]
    same(doc_hc["site_orientation_counts"], hc["site_orientation_counts"], "hc.orient")
    assert doc_hc["signal_pairs_unfolded_eligible"] == hc["signal_pairs_unfolded_eligible"]
    for n in hc["results"]:
        same(json.loads(json.dumps(doc_hc["results"][n]["unfolded"])),
             hc["results"][n]["unfolded"], f"hc.{n}.unfolded")
        same(json.loads(json.dumps(doc_hc["results"][n]["folded"])),
             ref["results"][n]["folded"], f"hc.{n}.folded")
    # pre-loaded inputs (sites tuple, anc dict, map dict) == path form
    sites, n_ind, counts = load_geno_tsv(geno)
    doc2 = estimate((sites, n_ind, counts), anc=load_anc(anc),
                    **dict(kw, gmap=load_map(gmap, "21")))
    same(json.loads(json.dumps(doc2["results"])), ref["results"], "preloaded")
    # bp mode, no polarization, no map: folded only; every retained pair is
    # mnv/signal/out_of_range; class totals == an independent per-pair sum.
    docb = estimate(geno, n=(4,), bins=[501, 5000, 20000, 60000], max_bp=70000,
                    mnv_bp=500, fmt="tsv")
    assert docb["dist_unit"] == "bp" and docb["site_orientation_counts"] == {"no_anc_row": 24}
    cb = docb["pair_class_counts"]
    assert cb["na_gap_excluded"] == 0 and cb["hotspot_control"] == 0 and cb["mnv_control"] == 2
    assert docb["results"]["n4"]["unfolded"] == {"bins": {}, "controls": {}}
    sites.sort(key=lambda s: s[0])
    tot_sig, tot_mnv, n_sig, n_oor = zeros(2), zeros(2), 0, 0
    for i in range(len(sites)):
        for j in range(i + 1, len(sites)):
            d = sites[j][0] - sites[i][0]
            if d > 70000:
                continue
            F = fold_matrix(project_grid(nine_cell(sites[i][3], sites[j][3]), n_ind, 2))
            if d <= 500:
                add_into(tot_mnv, F)
            elif 501 <= d < 60000:
                add_into(tot_sig, F); n_sig += 1
            else:
                n_oor += 1
    got_sig = zeros(2)
    for b in docb["results"]["n4"]["folded"]["bins"].values():
        add_into(got_sig, [[Fraction(x) for x in row] for row in b["phi_exact"]])
    assert got_sig == tot_sig and cb["signal"] == n_sig and \
        cb["out_of_bin_range"] == n_oor, "bp-mode signal total"
    got_mnv = [[Fraction(x) for x in row] for row in
               docb["results"]["n4"]["folded"]["controls"]["mnv_control"]["phi_exact"]]
    assert got_mnv == tot_mnv, "bp-mode mnv total"
    # orientation convention: case-folded ancestral base, conf governs
    am = {5: ("A", "G", "a", "low"), 6: ("A", "G", "G", "high"), 7: ("A", "G", "c", "low"),
          8: ("A", "G", "N", "none"), 9: ("A", "G", ".", "none"), 10: ("C", "T", "C", "high")}
    tt = [orient_site(p, "A", "G", am, False) for p in (4, 5, 6, 7, 8, 9, 10)]
    assert tt == [(None, "no_anc_row"), (1, "ok"), (-1, "ok"), (None, "anc_neither"),
                  (None, "unpolarizable"), (None, "unpolarizable"), (None, "refalt_mismatch")], tt
    assert orient_site(5, "A", "G", am, True) == (None, "low_conf_skipped")
    # genetic positions: gap sites have none; cumulative cM integrates covered intervals only
    mp = load_map(gmap, "chr21")
    assert mp["n_intervals"] == 6, mp["n_intervals"]
    ann = annotate_cm(sites, mp)
    assert sum(1 for c, _k in ann if c is None) == 2
    assert abs(ann[0][0] - (sites[0][0] - 14400000) * 1.0 / 1e6) < 1e-15
    return (f"cM mode == reference exactly (276 pairs: {cc}); high-conf-only, pre-loaded and "
            f"bp/folded-only forms consistent; bp totals == independent per-pair sums")


@leg("twosite_cli")
def leg_twosite_cli():
    import json
    refd = os.path.join(HERE, "reference", "twosite")
    exp = json.load(open(os.path.join(refd, "spectrum_expected.json")))
    bexp = json.load(open(os.path.join(refd, "brutecheck_expected.json")))
    prm = exp["params"]
    scratch = tempfile.mkdtemp(prefix="popcorn_twosite_")
    try:
        # projection engine CLI on the shipped fixture (default --geno)
        run_kw = dict(stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                      timeout=120)
        spec = os.path.join(HERE, "twosite_spectrum.py")
        r0 = subprocess.run([PY, os.path.join(HERE, "twosite_projection.py"), "--n", "4,6"],
                            **run_kw)
        assert r0.returncode == 0 and "EXACT_EQUAL_ALL = True" in r0.stdout, \
            f"rc={r0.returncode}\n{r0.stdout[-600:]}"
        # estimator CLI, --brute-check form with a JSON out
        bout = os.path.join(scratch, "brute.json")
        r1 = subprocess.run([PY, spec, "--brute-check", "--geno",
                             os.path.join(refd, bexp["geno"]), "--n", "4,6,8", "--out", bout],
                            **run_kw)
        assert r1.returncode == 0 and "EXACT_EQUAL_ALL = True" in r1.stdout, \
            f"rc={r1.returncode}\n{r1.stdout[-600:]}"
        assert json.load(open(bout))["per_n"] == bexp["per_n"], "CLI brute-check per_n != reference"
        # estimator CLI, cM mode, TSV input -> the reference spectra
        out = os.path.join(scratch, "spectrum.json")
        argv = ["--anc", os.path.join(refd, exp["anc"]), "--map", os.path.join(refd, exp["map"]),
                "--n", ",".join(str(x) for x in prm["n"]), "--bins", prm["bins"],
                "--max-bp", str(prm["max_bp"]), "--mnv-bp", str(prm["mnv_bp"]),
                "--hotspot-cm-per-mb", str(prm["hotspot_cm_per_mb"]),
                "--block-mb", str(prm["block_mb"])]
        r2 = subprocess.run([PY, spec, "--geno", os.path.join(refd, exp["geno"]),
                             "--chrom", prm["chrom"], "--format", "tsv", "--out", out] + argv,
                            **run_kw)
        assert r2.returncode == 0 and \
            r2.stdout.startswith("[twosite_spectrum] N=12 sites=24 unit=cM"), \
            f"rc={r2.returncode}\n{r2.stdout[-600:]}"
        doc = json.load(open(out))
        ref = exp["default"]
        for k in ("pair_class_counts", "site_orientation_counts", "mnv_control_pairs"):
            assert doc[k] == ref[k], k
        for n in ref["results"]:
            for reg in ref["results"][n]:
                for b, bd in ref["results"][n][reg]["bins"].items():
                    g = doc["results"][n][reg]["bins"][b]
                    assert g["phi_exact"] == bd["phi_exact"] and \
                        g["qhat_exact"] == bd["qhat_exact"] and \
                        g["n_pairs"] == bd["n_pairs"] and g["edges"] == bd["edges"], \
                        f"{n}.{reg}.bin{b}"
                for c, cd in ref["results"][n][reg]["controls"].items():
                    g = doc["results"][n][reg]["controls"][c]
                    assert g["phi_exact"] == cd["phi_exact"] and g["n_pairs"] == cd["n_pairs"] \
                        and g["bins_unbinned"] == cd["bins_unbinned"], f"{n}.{reg}.{c}"
        # the same panel as a VCF (mixed '/' and '|' GT separators, an extra
        # FORMAT field, plus one multiallelic, one indel, one missing-GT and
        # one other-chromosome record that must be skipped and counted),
        # read with --chrom 21 against 'chr21' records -> identical spectra.
        rows = [l.split("\t") for l in
                open(os.path.join(refd, exp["geno"])).read().splitlines()[1:]]
        vcf = os.path.join(scratch, "tiny.vcf")
        with open(vcf, "w") as fh:
            fh.write("##fileformat=VCFv4.2\n")
            fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                     + "\t".join(f"s{i}" for i in range(12)) + "\n")
            gt = {"0": ("0/0", "0|0"), "1": ("0/1", "1|0"), "2": ("1/1", "1|1")}
            for r_i, (pos, refb, alt, gts) in enumerate(rows):
                calls = "\t".join(gt[c][(r_i + k) % 2] + ":30" for k, c in enumerate(gts))
                fh.write(f"chr21\t{pos}\t.\t{refb}\t{alt}\t50\tPASS\t.\tGT:DP\t{calls}\n")
                if r_i == 3:
                    fh.write(f"chr21\t{int(pos)+1}\t.\tA\tC,T\t50\tPASS\t.\tGT:DP\t{calls}\n")
                    fh.write(f"chr21\t{int(pos)+2}\t.\tAT\tA\t50\tPASS\t.\tGT:DP\t{calls}\n")
                    miss = calls.replace("0/0:30", "./.:0", 1) if "0/0:30" in calls \
                        else "./.:0" + calls[6:]
                    fh.write(f"chr21\t{int(pos)+3}\t.\tA\tG\t50\tPASS\t.\tGT:DP\t{miss}\n")
                    fh.write(f"chr22\t{int(pos)+4}\t.\tA\tG\t50\tPASS\t.\tGT:DP\t{calls}\n")
        out2 = os.path.join(scratch, "spectrum_vcf.json")
        r3 = subprocess.run([PY, spec, "--geno", vcf, "--chrom", "21", "--out", out2] + argv,
                            **run_kw)
        assert r3.returncode == 0, f"rc={r3.returncode}\n{r3.stdout[-600:]}"
        doc2 = json.load(open(out2))
        assert doc2["inputs"]["format"] == "vcf" and doc2["n_sites"] == 24
        assert doc2["site_counts"] == {"nonbiallelic_or_multiallelic": 1, "non_snv": 1,
                                      "missing_gt": 1, "other_chrom": 1, "retained": 24}, \
            doc2["site_counts"]
        assert doc2["results"] == doc["results"] and \
            doc2["pair_class_counts"] == doc["pair_class_counts"], "VCF form != TSV form"
        # stdout form ('--out -'): JSON on stdout, summary on stderr; bp mode without --anc/--map
        pipe_kw = dict(run_kw, stderr=subprocess.PIPE)
        r4 = subprocess.run([PY, spec, "--geno", vcf, "--n", "4", "--bins",
                             "501,5000,20000,60000", "--mnv-bp", "500", "--max-bp", "70000"],
                            **pipe_kw)
        assert r4.returncode == 0 and "unit=bp" in r4.stderr, \
            f"rc={r4.returncode}\n{r4.stderr[-600:]}"
        d4 = json.loads(r4.stdout)
        assert d4["dist_unit"] == "bp" and d4["results"]["n4"]["unfolded"]["bins"] == {}
        # an odd n is a usage error (rc 2), not a traceback
        r5 = subprocess.run([PY, spec, "--geno", vcf, "--n", "5"], **pipe_kw)
        assert r5.returncode == 2 and "n must be even" in r5.stderr, f"rc={r5.returncode}"
    finally:
        shutil.rmtree(scratch, ignore_errors=True)
    return ("projection CLI + estimator CLI (--brute-check, TSV cM mode == reference, "
            "VCF form == TSV form with 4 skipped records counted, bp mode on stdout, "
            "usage error rc=2)")


@leg("enclosure_cell")
def leg_enclosure_cell():
    if not _has("flint"):
        return "SKIP:python-flint absent — pip install python-flint"
    import popcorn.enclosure as EN
    had_moments = "moments" in sys.modules
    cert = EN.certified_cell                    # first touch loads the engine
    eng = sys.modules.get("transient_enclosure_engine")
    assert eng is not None and EN.transient_enclosure_engine is eng, \
        "engine not aliased by identity"
    r40 = cert(40, -5, 20, 256, verbose=False)
    r60 = cert(60, -5, 20, 256, verbose=False)
    for tag, r in (("M=40", r40), ("M=60", r60)):
        g = r["gates"]
        assert g["stationarity_interior_contains_zero"] is True, f"{tag}: stationarity check"
        assert g["eq_projection_overlaps_closed_form"] is True, f"{tag}: eq-projection check"
        assert len(r["entry_balls"]) == 19 and len(r["enclosure"]) == 19, f"{tag}: 19 entries"
    w40 = max(e["rel_width"] for e in r40["enclosure"])
    w60 = max(e["rel_width"] for e in r60["enclosure"])
    assert w40 <= 1e-13, f"M=40 max rel width {w40:.2e} (bar 1e-13)"
    assert w60 <= 1e-34, f"M=60 max rel width {w60:.2e} (bar 1e-34)"
    assert w60 < 1e-15 * w40, "width does not shrink with M (truncation not self-reporting)"
    nov = sum(a.overlaps(b) for a, b in zip(r40["entry_balls"], r60["entry_balls"]))
    nin = sum(a.contains(b) for a, b in zip(r40["entry_balls"], r60["entry_balls"]))
    assert nov == 19, f"M=40 and M=60 enclosures overlap at {nov}/19 entries"
    mono = all(float(b.mid()) > 0 for b in r60["entry_balls"])
    assert mono, "non-positive entry"
    assert ("moments" in sys.modules) == had_moments, "engine imported moments eagerly"
    fdev = None
    if _has("numpy") and _has("scipy"):
        arb = eng.arb
        fl = EN.float_cell(60, -5, 20)["sfs"]
        fdev = 0.0
        for b, f in zip(r60["entry_balls"], fl):
            mid, rad = float(b.mid()), float(b.rad())
            fdev = max(fdev, abs(f - mid) / abs(mid))
            assert abs(f - mid) <= rad + 1e-8 * abs(mid), \
                f"float route {f!r} vs enclosure mid {mid!r}: rel {abs(f - mid) / abs(mid):.2e} (bar 1e-8)"
    ftxt = f"float route within {fdev:.1e} rel (bar 1e-8)" if fdev is not None \
        else "(float route skipped: numpy/scipy absent)"
    return (f"S=-5 n=20: checks true; rel width M=40 {w40:.1e} -> M=60 {w60:.1e}; "
            f"overlap 19/19 (nested {nin}/19); {ftxt}; moments not imported")


@leg("enclosure_stationary")
def leg_enclosure_stationary():
    if not _has("flint"):
        return "SKIP:python-flint absent — pip install python-flint"
    from fractions import Fraction
    from mpmath import mp, mpf, fabs, log10
    from popcorn.sfs import solve_M_exact, sfs_engine
    import popcorn.enclosure as EN
    eng = EN.transient_enclosure_engine
    eng.ctx.prec = 256
    worst_e = worst_p = 9999.0
    checked = 0
    for n, S, do_proj in ((20, -5, True), (20, -50, True), (60, -20, False)):
        exact = sfs_engine.sfs_from_M(n, S, solve_M_exact(n, Fraction(S)), 72)   # exact route, integer S
        weq = [eng.weq_ball(j, S) for j in range(1, n)] if do_proj else None
        with mp.workdps(90):
            for i in range(1, n):
                truth = exact[i - 1]
                e = eng.equilibrium_entry(n, i, S)
                got = mpf(e.mid().str(80, radius=False))
                d = float(-log10(fabs((got - truth) / truth))) if got != truth else 9999.0
                worst_e = min(worst_e, d)
                tball = eng.arb(f"[{mp.nstr(truth, 50)} +/- {float(fabs(truth)) * 1e-48:.2e}]")
                assert e.overlaps(tball), f"n={n} i={i} S={S}: ball misses the popcorn.sfs value (50 digits)"
                checked += 1
                if do_proj:
                    pj = eng.project_entry(n, i, weq)
                    gotp = mpf(pj.mid().str(80, radius=False))
                    dp = float(-log10(fabs((gotp - truth) / truth))) if gotp != truth else 9999.0
                    worst_p = min(worst_p, dp)
    assert worst_e >= 30, f"equilibrium_entry vs popcorn.sfs exact route only {worst_e:.1f}d (bar 30)"
    assert worst_p >= 30, f"projected weq moments vs popcorn.sfs only {worst_p:.1f}d (bar 30)"
    for i in range(1, 20):                                  # neutral: exactly 1/i
        assert eng.equilibrium_entry(20, i, 0).overlaps(eng.qarb(eng.fmpq(1, i)))
        pj = eng.project_entry(20, i, [eng.weq_ball(j, 0) for j in range(1, 20)])
        assert pj.overlaps(eng.qarb(eng.fmpq(1, i))), f"S=0 projection at i={i} misses 1/i"
    return (f"{checked} stationary entries vs popcorn.sfs exact route (theta=1, same S): "
            f"closed form {worst_e:.0f}d, projected moments {worst_p:.0f}d (bar 30); S=0 -> 1/i")


@leg("enclosure_reference")
def leg_enclosure_reference():
    if not _has("flint"):
        return "SKIP:python-flint absent — pip install python-flint"
    import popcorn.enclosure as EN
    eng = EN.transient_enclosure_engine
    cells = {(M, S): EN.load_reference_cell(M, S) for (M, S) in EN.REFERENCE_CELLS}
    assert len(cells) == 6, f"{len(cells)} reference cells"
    eng.ctx.prec = 256
    worst_w = 0.0
    for (M, S), c in cells.items():
        assert (c["M"], c["S"], c["n"], c["prec_bits"]) == (M, S, 20, 256), f"header {(M, S)}"
        assert c["epochs"] == [[str(r), str(t)] for r, t in eng.EPOCHS], f"{(M, S)}: history"
        g = c["gates"]
        assert g["stationarity_interior_contains_zero"] is True and \
            g["eq_projection_overlaps_closed_form"] is True, f"{(M, S)}: in-cell checks"
        assert g["stationarity_max_residual_radius"] < 1e-70, f"{(M, S)}: residual radius"
        assert len(c["enclosure"]) == 19 and [e["i"] for e in c["enclosure"]] == list(range(1, 20))
        worst_w = max(worst_w, max(e["rel_width"] for e in c["enclosure"]))
        for e, f in zip(c["enclosure"], c["float_scipy"]["sfs"]):
            b = EN.entry_ball(e)
            mid = float(b.mid())
            assert abs(f - mid) <= 1e-7 * abs(mid), f"{(M, S)} i={e['i']}: float route off by {abs(f - mid) / abs(mid):.1e}"
    assert worst_w <= 1e-60, f"max rel width {worst_w:.1e} (bar 1e-60)"
    same30 = nov = 0
    for S in (0, -5, -50):
        for e2, e5 in zip(cells[(200, S)]["enclosure"], cells[(500, S)]["enclosure"]):
            same30 += int(e2["mid30"] == e5["mid30"])
            nov += int(EN.entry_ball(e2).overlaps(EN.entry_ball(e5)))
    assert same30 == 57 and nov == 57, f"M=200 vs M=500: identical mid30 {same30}/57, overlap {nov}/57"
    assert cells[(200, 0)]["moment_diag"]["bracket_sep_mid_jM"] == "0", "S=0 closure not exact"
    # live engine vs shipped S=-5 cells: small-M enclosures must overlap them
    r60 = EN.certified_cell(60, -5, 20, 256, verbose=False)
    live = sum(b.overlaps(EN.entry_ball(e2)) and b.overlaps(EN.entry_ball(e5))
               for b, e2, e5 in zip(r60["entry_balls"], cells[(200, -5)]["enclosure"],
                                     cells[(500, -5)]["enclosure"]))
    assert live == 19, f"live M=60 enclosure overlaps shipped M=200/500 at {live}/19"
    # FLOAT-HP independent route, live (128 substeps) vs the shipped instance;
    # worst_rel_dev is measured against the hull midpoint, so at M = 40 it is
    # the half LOWER/UPPER bracket (1.5e-32) and substep-independent
    ref = EN.load_reference_xcheck()
    rx = EN.xcheck(M=ref["M"], S=ref["S"], n=20, prec=256, dps=ref["dps"], nsub=128, verbose=False)
    assert rx["worst_rel_dev"] < 1e-30, f"xcheck worst rel dev {rx['worst_rel_dev']:.2e} (bar 1e-30)"
    assert abs(rx["worst_rel_dev"] / ref["worst_rel_dev"] - 1.0) <= 1e-3, \
        f"xcheck {rx['worst_rel_dev']:.6e} vs shipped {ref['worst_rel_dev']:.6e}"
    return (f"6 cells: checks true, rel width <= {worst_w:.1e}, float route <= 1e-7; "
            f"M=200 vs M=500 mid30 identical 57/57 + overlap; live M=60 overlaps shipped 19/19; "
            f"mpmath route {rx['worst_rel_dev']:.4e} (= shipped, bit-identical "
            f"{rx['worst_rel_dev'] == ref['worst_rel_dev']})")


def _planted_mode():
    """('BIT-IDENTICAL' | 'VERSION-DRIFT', drift dict): whether the live
    msprime/tskit/numpy equal the versions the reference set was generated
    with (reference/planted/pins.json 'versions')."""
    import json
    from popcorn.planted import versions
    want = json.load(open(os.path.join(HERE, "reference", "planted",
                                       "pins.json")))["versions"]
    live = versions()
    drift = {k: (live.get(k), want[k]) for k in ("msprime", "tskit", "numpy")
             if live.get(k) != want[k]}
    return ("VERSION-DRIFT" if drift else "BIT-IDENTICAL"), drift


@leg("planted_analytic")
def leg_planted_analytic():
    if not _has("numpy"):
        return "SKIP:numpy absent — pip install numpy"
    import json
    import numpy as np
    from popcorn.planted import (expected_branch_sfs, tavare_coefs,
                                 analytic_selfchecks, reference_checks,
                                 normalized_sfs, load_pins, TK1_EPOCHS,
                                 TK2_TIMES, TRUTH_NS)
    assert analytic_selfchecks() is True
    # Tavare transition function at tau = 0: P(A_n(0) = k) = 1{k = n}
    for n in (2, 4, 6, 9):
        co = tavare_coefs(n)
        for k in range(2, n + 1):
            p0 = sum(c for c, _lam in co[k])
            assert abs(p0 - (1.0 if k == n else 0.0)) < 1e-9, (n, k, p0)
    # constant N: E[L_i] = 2N/i; epoch splitting is inert; (t, N) -> (cN, ct)
    # rescales E[L] by c
    N = 8000.0
    for n in (3, 5, 8):
        EL = expected_branch_sfs([(0.0, N)], n)
        assert np.allclose(EL, [2 * N / i for i in range(1, n)], rtol=1e-12, atol=0)
        split = expected_branch_sfs([(0.0, N), (300.0, N), (2500.0, N)], n)
        assert np.allclose(split, EL, rtol=1e-12, atol=0), "epoch split not inert"
    base = expected_branch_sfs(TK1_EPOCHS, 6)
    scaled = expected_branch_sfs([(3 * t, 3 * s) for t, s in TK1_EPOCHS], 6)
    assert np.allclose(scaled, 3 * base, rtol=1e-12, atol=0), "time rescaling"
    # shipped reference: TK2 fit and exact E[L_i](TK1, 6) reproduced to 1e-12
    dev = reference_checks()
    assert set(dev) == {"tk2_match_rel", "exact_EL_rel"}, dev
    pins = load_pins()
    ach = pins["tk2_match"]["achieved"]
    worst_match = max(ach[f"n{n}"]["max_abs_rel_diff"] for n in TRUTH_NS)
    assert worst_match <= 1e-4, f"TK1/TK2 SFS match only {worst_match:.2e}"
    # planted Monte Carlo xi (2e6 replicates) vs the closed form, |z| <= 4
    tk2 = [(TK2_TIMES[m], pins["tk2_fitted_sizes"][m]) for m in range(len(TK2_TIMES))]
    refd = os.path.join(HERE, "reference", "planted")
    zmax = 0.0
    for cls, ep in (("TK1", TK1_EPOCHS), ("TK2", tk2)):
        t = json.load(open(os.path.join(refd, f"truth_{cls}.json")))
        for n in TRUTH_NS:
            b = t["mc"][f"n{n}"]
            assert b["reps"] == 2000000, f"{cls} n{n} reps {b['reps']}"
            xi = normalized_sfs(ep, n)
            z = max(abs(m - e) / s for m, e, s in
                    zip(b["xi_1sfs"], xi, b["xi_1sfs_mc_se"]))
            assert z <= 4.0, f"{cls} n={n}: planted xi vs closed form |z|={z:.2f}"
            zmax = max(zmax, z)
    return (f"Tavare identities; reference fit/E[L] reproduced "
            f"{dev['tk2_match_rel']:.0e}/{dev['exact_EL_rel']:.0e}; TK1~TK2 SFS "
            f"match {worst_match:.1e}; planted xi(TK1,TK2) vs closed form max|z| {zmax:.2f}")


@leg("planted_vs_lambda")
def leg_planted_vs_lambda():
    import json
    from fractions import Fraction
    from popcorn.lambda_coalescent import (xi_hat, beta_rate, dirac_rate,
                                          msprime_dirac_rate)
    refd = os.path.join(HERE, "reference", "planted")
    rows = (("B13", beta_rate(Fraction(13, 10))), ("B17", beta_rate(Fraction(17, 10))),
            ("D005", msprime_dirac_rate(Fraction(1, 20), 10)),
            ("D020", msprime_dirac_rate(Fraction(1, 5), 10)))
    zmax = 0.0
    cells = 0
    for cls, lam in rows:
        t = json.load(open(os.path.join(refd, f"truth_{cls}.json")))
        for nk, b in t["mc"].items():
            exact = [float(v) for v in xi_hat(b["n"], lam)]
            assert len(exact) == len(b["xi_1sfs"]) == b["n"] - 1
            for m, e, s in zip(b["xi_1sfs"], exact, b["xi_1sfs_mc_se"]):
                z = abs(m - e) / s
                assert z <= 4.0, f"{cls} {nk}: planted xi vs exact Lambda spectrum |z|={z:.2f}"
                zmax = max(zmax, z)
                cells += 1
    # planted control: the one-atom Dirac(psi) law is NOT msprime's model and
    # must be rejected by the same comparison
    t = json.load(open(os.path.join(refd, "truth_D020.json")))
    b = t["mc"]["n6"]
    wrong = [float(v) for v in xi_hat(6, dirac_rate(Fraction(1, 5)))]
    zbad = max(abs(m - e) / s for m, e, s in zip(b["xi_1sfs"], wrong, b["xi_1sfs_mc_se"]))
    assert zbad > 10.0, f"control: one-atom Dirac not rejected (|z|={zbad:.1f})"
    return (f"planted xi of B13,B17,D005,D020 (n=4,6; {cells} cells) vs exact "
            f"Beta/Dirac spectra max|z| {zmax:.2f} <= 4; one-atom Dirac control "
            f"rejected |z| {zbad:.0f}")


@leg("planted_vs_twolocus")
def leg_planted_vs_twolocus():
    if not (_has("numpy") and _has("scipy")):
        return "SKIP:numpy/scipy absent — pip install numpy scipy"
    import json
    import numpy as np
    from popcorn.planted import load_pins, TK1_EPOCHS, TK2_TIMES, TRUTH_NS
    from popcorn.twolocus import moments
    pins = load_pins()
    tk2 = [(TK2_TIMES[m], pins["tk2_fitted_sizes"][m]) for m in range(len(TK2_TIMES))]
    refd = os.path.join(HERE, "reference", "planted")
    zq = zx = 0.0
    for cls, ep in (("TK1", TK1_EPOCHS), ("TK2", tk2)):
        N0 = ep[0][1]                     # (t_gen, N_chrom) -> (tau, eta)
        hist = tuple((t / N0, N0 / s) for t, s in ep)
        t = json.load(open(os.path.join(refd, f"truth_{cls}.json")))
        for n in TRUTH_NS:
            b = t["mc"][f"n{n}"]
            M, a, _b = moments(n, 0.0, hist)
            M = np.array(M)
            q = M / M.sum()
            z = np.abs((np.array(b["q_ij"]) - q) / np.array(b["q_ij_mc_se"]))
            assert z.max() <= 4.0, f"{cls} n={n}: planted q_ij vs two-locus rho=0 |z|={z.max():.2f}"
            zq = max(zq, float(z.max()))
            xi = np.array(a) / np.sum(a)
            w = np.abs((np.array(b["xi_1sfs"]) - xi) / np.array(b["xi_1sfs_mc_se"]))
            assert w.max() <= 4.0, f"{cls} n={n}: planted xi vs two-locus E[T_i] |z|={w.max():.2f}"
            zx = max(zx, float(w.max()))
    # the recorded TK2-vs-TK1 two-site distance is consistent with the two
    # shipped truths (recomputed here)
    t1 = json.load(open(os.path.join(refd, "truth_TK1.json")))
    t2 = json.load(open(os.path.join(refd, "truth_TK2.json")))
    sep = 0.0
    for n in TRUTH_NS:
        b1, b2 = t1["mc"][f"n{n}"], t2["mc"][f"n{n}"]
        d = np.abs(np.array(b2["q_ij"]) - np.array(b1["q_ij"]))
        se = np.sqrt(np.array(b2["q_ij_mc_se"]) ** 2 + np.array(b1["q_ij_mc_se"]) ** 2)
        mz = float((d / se).max())
        assert abs(mz - t2["vs_TK1_2sfs"][f"n{n}"]["max_z"]) <= 1e-9, f"vs_TK1_2sfs n{n}"
        sep = max(sep, mz)
    return (f"planted q_ij of TK1,TK2 (n=4,6) vs deterministic two-locus rho=0 "
            f"max|z| {zq:.2f}, xi {zx:.2f} <= 4; SFS-matched pair 2-SFS distance "
            f"{sep:.2f} SE (recorded block consistent)")


@leg("planted_sim")
def leg_planted_sim():
    for m in ("numpy", "msprime", "tskit"):
        if not _has(m):
            return f"SKIP:{m} absent — pip install msprime (simulation legs only)"
    import contextlib
    import io
    import json
    from popcorn.planted import selftest
    exp = json.load(open(os.path.join(HERE, "reference", "planted",
                                       "selftest_expected.json")))
    mode, drift = _planted_mode()
    with contextlib.redirect_stdout(io.StringIO()):   # generators print paths
        rec = selftest()                      # temp dir, removed by selftest
    assert rec["pass"] and rec["out_dir"] is None, "selftest verdict"
    z, ze = rec["zcheck_exact_vs_mc"], exp["zcheck_exact_vs_mc"]
    assert z["seed"] == ze["seed"] and z["max_abs_z"] < 4.0, \
        f"zcheck seed {z['seed']} / max|z| {z['max_abs_z']:.3f}"
    rel = max(abs(a - b) / abs(b) for a, b in zip(z["exact_EL_gen"], ze["exact_EL_gen"]))
    assert rel <= 1e-12, f"exact E[L_i] rel {rel:.2g}"
    st, se_ = rec["subscale_truth_TK1"], exp["subscale_truth_TK1"]
    assert st["seeds"] == se_["seeds"], "subscale truth seeds"
    arm, ae = rec["subscale_arm_TK1_rON_GC"], exp["subscale_arm_TK1_rON_GC"]
    assert arm["seeds"] == ae["seeds"], "subscale arm seeds"
    dz = abs(z["max_abs_z"] - ze["max_abs_z"])
    dq = max(abs(x - y) for nk in ("n4", "n6")
             for xr, yr in zip(st["q_ij"][nk], se_["q_ij"][nk]) for x, y in zip(xr, yr))
    same_counts = arm["counts"] == ae["counts"]
    same_payload = arm["payload_sha256"] == ae["payload_sha256"]
    if mode == "BIT-IDENTICAL":
        assert dz <= 1e-9, f"zcheck max|z| {z['max_abs_z']!r} != {ze['max_abs_z']!r}"
        assert dq <= 1e-12, f"subscale q_ij max|dq| {dq:.2g}"
        assert same_counts, f"arm counts {arm['counts']} != {ae['counts']}"
        assert same_payload, "arm payload hashes differ"
    else:
        assert arm["counts"]["kept_biallelic_snv"] > 0
    return (f"[{mode}{'' if not drift else ' ' + str(drift)}] reduced selftest: "
            f"zcheck max|z| {z['max_abs_z']:.6f} (ref {ze['max_abs_z']:.6f}); "
            f"TK1@300 max|dq| {dq:.1e}; GC arm 40 dip/300 kb kept_snv "
            f"{arm['counts']['kept_biallelic_snv']} trees {arm['counts']['num_trees']} "
            f"counts {'==' if same_counts else '!='} payload {'==' if same_payload else '!='} "
            f"({rec['wall_s_total']:.1f}s)")


@leg("twowindow_laws")
def leg_twowindow_laws():
    if not _has("numpy"):
        return "SKIP:numpy absent — pip install numpy"
    from fractions import Fraction as Fr
    import numpy as np
    import popcorn.twowindow as TW
    # geometry law (pure arithmetic): pa = margin, pb = pa + d, L = pb + 1 + margin
    assert TW.layout(10_000, TW.margin_for(500)) == (5000, 15000, 20001), TW.layout(10_000, TW.margin_for(500))
    assert TW.layout(10_000, TW.margin_for(55)) == (3000, 13000, 16001) == TW.layout(10_000), "default-margin layout"
    assert (TW.margin_for(), TW.margin_for(300), TW.margin_for(300.5), TW.margin_for(None, 100)) \
        == (3000, 3000, 3005, 100), "margin law max(margin_bp, 10 * tract)"
    e, ra, rb = TW.window_edges(*TW.layout(100, 300))
    assert e.tolist() == [0, 300, 301, 400, 401, 701] and (ra, rb) == (1, 3), f"edges {e.tolist()} rows {(ra, rb)}"
    assert TW.probe_weights((0, 100, 300)).tolist() == [100, 150, 200] \
        and TW.probe_weights((-5000, 0, 5000)).tolist() == [5000, 5000, 5000], "trapezoid probe weights"
    assert TW.probe_positions(20_000, (-5000, 0, 5000)).tolist() == [5000, 10000, 15000]
    # seed law
    assert TW.seed_from_name("popcorn.twowindow") == 1208489443 == TW.root_seed("popcorn.twowindow")
    assert TW.root_seed(290993229) == 290993229
    for bad in (True, 1.5):
        try:
            TW.root_seed(bad)
            raise AssertionError(f"root_seed accepted {bad!r}")
        except TypeError:
            pass
    # pooling convention (diagonal doubled) on hand values; refusals
    keys = TW.pair_keys(8)
    assert len(keys) == 28 and keys[0] == (1, 1) and keys[7] == (2, 2) and keys[-1] == (7, 7) \
        and keys == TW.keys28(), "pair_keys(8)"
    Q = np.array([[1.0, 2.0], [4.0, 3.0]])
    assert TW.pair_counts(Q).tolist() == [2, 6, 6] and np.allclose(TW.pool_pairs(Q), [1 / 7, 3 / 7, 3 / 7], atol=1e-15)
    assert TW.pool_pairs(Q[None, :, :]).tolist() == TW.v28_from_M(Q).tolist(), "1 x m x m form / alias"
    try:
        TW.pool_pairs(np.ones((2, 3)))
        raise AssertionError("non-square matrix accepted")
    except ValueError:
        pass
    # delete-one-block jackknife on hand values: blocks with pair counts (1,1,1) and (3,1,1)
    A = np.array([[0.5, 0.5], [0.5, 0.5]])
    B = np.array([[1.5, 0.5], [0.5, 0.5]])
    v, se = TW.block_jackknife([A, B])
    assert np.allclose(v, [0.5, 0.25, 0.25], atol=1e-15) and np.allclose(se, [2 / 15, 1 / 15, 1 / 15], atol=1e-15), (v, se)
    qh, Sig, nb = TW.jackknife_cov({"b0": [1, 1, 1], "b1": [3, 1, 1]})
    assert nb == 2 and np.allclose(qh, v, atol=1e-15) and np.allclose(np.sqrt(np.diag(Sig)), se, atol=1e-15)
    assert abs(Sig[0, 1] + 2 / 225) < 1e-15, "jackknife off-diagonal"
    v1, se1 = TW.block_jackknife_v28([A, B])
    assert v1.tolist() == v.tolist() and se1.tolist() == se.tolist(), "alias"
    try:
        TW.block_jackknife([A])
        raise AssertionError("single block accepted")
    except ValueError:
        pass
    # sim_kwargs_for: gene conversion needs a tract; batch/seed keys refused
    kw = TW.sim_kwargs_for(701, gc_rate=1e-8, gc_tract=300)
    assert kw["gene_conversion_tract_length"] == 300.0 and kw["recombination_rate"] == 0.0 and kw["samples"] == 4
    for bad in (dict(gc_rate=1e-8), dict(extra_sim_kwargs={"random_seed": 1})):
        try:
            TW.sim_kwargs_for(701, **bad)
            raise AssertionError(f"sim_kwargs_for accepted {bad}")
        except ValueError:
            pass
    # shipped exact comparand (n = 8 Kingman, completely linked): structure + live exact engine
    ref = TW.load_reference(8)
    assert ref["keys"] == keys and sum(ref["pooled"]) == 1 and ref["E_Li"] == [Fr(2, i) for i in range(1, 8)]
    tot = sum(ref["E_LiLj"].values())
    assert all(ref["pooled"][k] == ref["E_LiLj"][keys[k]] / tot for k in range(28)), "pooled != E_LiLj / sum"
    try:
        k2, ex = TW.exact_linked_pooled(8)
        assert k2 == keys and ex == ref["pooled"], "live popcorn.twosfs != shipped reference"
        T = TW._exact_engine()
        u2, v2 = T.kingman_moments_route2(8)
        assert all(v2[k] == ref["E_LiLj"][k] for k in keys) and [u2[i] for i in range(1, 8)] == ref["E_Li"], \
            "route-2 Kingman moments != shipped reference"
        live = f"live popcorn.twosfs (two exact routes) == shipped n=8 reference via {T.__name__}"
    except ImportError as e:
        live = f"(live popcorn.twosfs cross-check skipped: {e})"
    return f"geometry/margin/seed/pooling/jackknife laws on literals; {live}"


@leg("twowindow_linked")
def leg_twowindow_linked():
    if not _has("numpy"):
        return "SKIP:numpy absent — pip install numpy"
    if not _has("msprime"):
        return "SKIP:msprime absent — pip install msprime (brings tskit)"
    import numpy as np
    import popcorn.twowindow as TW
    ref = TW.load_reference(8)
    exact = np.array([float(x) for x in ref["pooled"]])
    NREP, B = 16_000, 40
    r = TW.run_two_locus(0.0, NREP, seed=290993229, n_blocks=B, batch=400)
    assert r["done"] == NREP and r["failed"] == 0, f"done {r['done']} failed {r['failed']}: {r['last_error']}"
    assert r["n_afs_distinct"] == 0, f"vacuity check: n_afs_distinct = {r['n_afs_distinct']} at rho = 0"
    M = r["M"]
    assert M.shape == (7, 7) and np.array_equal(M, M.T), "rho = 0 matrix not exactly symmetric"
    assert M.min() > 0 and np.isfinite(M).all(), "simulator returned zeros / non-finite"
    v, se = r["pooled"], r["pooled_se"]
    assert r["keys"] == ref["keys"] and v.shape == (28,) and abs(v.sum() - 1) < 1e-12
    assert (v > 0).all() and (se > 0).all(), "empty cell or zero SE"
    z = (v - exact) / se
    maxz, chi2 = float(np.max(np.abs(z))), float((z ** 2).mean())
    assert maxz <= 4.0, f"pooled 28-vector vs exact Kingman: max|z| {maxz:.2f} > 4 ({B} blocks; z = {np.round(z, 2).tolist()})"
    i = np.arange(1, 8)
    indep = TW.pool_pairs(np.outer(2.0 / i, 2.0 / i))          # unlinked limit outer(E[L_i], E[L_j])
    power = float(np.max(np.abs(indep - exact) / se))
    assert power > 20, f"no power: unlinked limit only {power:.1f} SE away"
    # contiguous two-window geometry: no process -> same tree; gene conversion -> decorrelation
    w0 = TW.run_two_window(100, 1_000, seed="selftest.twowindow.noproc", margin=300, batch=250)
    assert w0["done"] == 1_000 and w0["failed"] == 0 and w0["n_afs_distinct"] == 0 \
        and np.array_equal(w0["M"], w0["M"].T), f"two-window no-process arm: distinct {w0['n_afs_distinct']}"
    assert (w0["params"]["window_a_bp"], w0["params"]["window_b_bp"], w0["params"]["sequence_length"]) == (300, 400, 701)
    wg = TW.run_two_window(100, 1_000, seed="selftest.twowindow.gc", margin=300, batch=250,
                           gc_rate=5e-7, gc_tract=30.0)
    assert wg["done"] == 1_000 and wg["failed"] == 0 and wg["n_afs_distinct"] > 100, \
        f"gene conversion arm: distinct {wg['n_afs_distinct']}/1000"
    return (f"n=8 rho=0, {NREP} reps: max|z| {maxz:.2f}, chi2/28 {chi2:.2f} vs exact Kingman "
            f"({B}-block jackknife; unlinked limit {power:.0f} SE away); n_afs_distinct 0; "
            f"two-window d=100: distinct 0 -> {wg['n_afs_distinct']}/1000 with gene conversion")


@leg("twowindow_recomb")
def leg_twowindow_recomb():
    if not _has("numpy"):
        return "SKIP:numpy absent — pip install numpy"
    if not _has("msprime"):
        return "SKIP:msprime absent — pip install msprime (brings tskit)"
    import numpy as np
    import popcorn.twowindow as TW
    RHO, NREP, B = 1.0, 8_000, 40
    r = TW.run_two_locus(RHO, NREP, seed=290993229, samples=3, n_blocks=B, batch=200)
    assert r["done"] == NREP and r["failed"] == 0, f"done {r['done']} failed {r['failed']}: {r['last_error']}"
    assert r["params"]["n_haploid"] == 6 and r["M"].shape == (5, 5) and len(r["keys"]) == 15
    frac = r["afs_distinct_fraction"]
    assert frac > 0.5, f"rho = {RHO}: afs_distinct_fraction {frac:.3f} — sites not decorrelating (rate-map orientation?)"
    v, se = r["pooled"], r["pooled_se"]
    assert abs(v.sum() - 1) < 1e-12 and (v > 0).all() and (se > 0).all()
    parts = [f"n=6 rho={RHO:g}, {NREP} reps: distinct fraction {frac:.3f}"]
    try:
        k0, ex0 = TW.exact_linked_pooled(6)
        z0 = float(np.max(np.abs((v - np.array([float(x) for x in ex0])) / se)))
        assert k0 == r["keys"] and z0 > 10, f"arm indistinguishable from the linked spectrum (max|z| {z0:.1f})"
        parts.append(f"{z0:.0f} SE from the exact linked spectrum")
    except ImportError as e:
        parts.append(f"(linked-spectrum contrast skipped: {e})")
    if _has("scipy"):
        import popcorn.twolocus as TL
        q = TL.pool_pairs(TL.moments(6, RHO)[0])
        z = (v - q) / se
        maxz, chi2 = float(np.max(np.abs(z))), float((z ** 2).mean())
        assert maxz <= 4.5, (f"pooled vector vs popcorn.twolocus at rho={RHO:g}: max|z| {maxz:.2f} > 4.5 "
                             f"({B} blocks; z = {np.round(z, 2).tolist()})")
        zh = float(np.max(np.abs((v - TL.pool_pairs(TL.moments(6, RHO / 2)[0])) / se)))
        zd = float(np.max(np.abs((v - TL.pool_pairs(TL.moments(6, 2 * RHO)[0])) / se)))
        assert min(zh, zd) > 5, f"no power against a factor 2 in rho: {zh:.1f}, {zd:.1f}"
        parts.append(f"vs popcorn.twolocus E[T_i^A T_j^B]: max|z| {maxz:.2f}, chi2/15 {chi2:.2f} "
                     f"({B} blocks); rho/2 and 2 rho rejected at {zh:.0f}, {zd:.0f} SE")
    else:
        parts.append("(popcorn.twolocus cross-check skipped: scipy absent — pip install scipy)")
    return "; ".join(parts)


# --------------------------------------------------------------- full tier
@leg("gate_phase1")
def leg_gate_phase1(scratch):
    if not _has("flint"):
        return "SKIP:python-flint absent — pip install python-flint"
    cwd = os.path.join(scratch, "gate_phase1")
    os.makedirs(cwd, exist_ok=True)
    r = subprocess.run([PY, os.path.join(HERE, "gate_phase1.py")],
                       cwd=cwd, stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, text=True, timeout=3600)
    assert r.returncode == 0 and '"PASS": true' in r.stdout, \
        f"rc={r.returncode}\n{r.stdout[-1200:]}"
    return "held-out two-route gate >=40d"


@leg("dominance_full")
def leg_dominance_full(scratch):
    r = subprocess.run([PY, os.path.join(HERE, "dominance_oracle.py")],
                       cwd=scratch, stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, text=True, timeout=3600)
    assert r.returncode == 0, f"rc={r.returncode}\n{r.stdout[-1200:]}"
    return "full oracle selftest (collapse x3 + hostile x6)"


@leg("dfe_full")
def leg_dfe_full(scratch):
    if not _has("flint"):
        return "SKIP:python-flint absent — pip install python-flint"
    r = subprocess.run([PY, os.path.join(HERE, "dfe_layer.py")],
                       cwd=scratch, stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, text=True, timeout=3600)
    assert r.returncode == 0, f"rc={r.returncode}\n{r.stdout[-1200:]}"
    return "n=20 dps=40 kernel selftest"


@leg("enclosure_full")
def leg_enclosure_full(scratch):
    if not _has("flint"):
        return "SKIP:python-flint absent — pip install python-flint"
    import popcorn.enclosure as EN
    ref = EN.load_reference_cell(200, -5)
    rec = EN.certified_cell(200, -5, 20, 256, verbose=False)      # ~25 s
    g = rec["gates"]
    assert g["stationarity_interior_contains_zero"] is True and \
        g["eq_projection_overlaps_closed_form"] is True, "in-cell checks"
    same = sum(a["mid30"] == b["mid30"] for a, b in zip(rec["enclosure"], ref["enclosure"]))
    nov = sum(b.overlaps(EN.entry_ball(e)) for b, e in zip(rec["entry_balls"], ref["enclosure"]))
    assert same == 19, f"mid30 strings identical at {same}/19"
    assert nov == 19, f"balls overlap shipped at {nov}/19"
    for k in ("bracket_sep_mid_j1", "bracket_sep_mid_jn", "bracket_sep_mid_jM"):
        assert rec["moment_diag"][k] == ref["moment_diag"][k], f"{k} differs from shipped"
    xr = EN.load_reference_xcheck()
    rx = EN.xcheck(verbose=False)                                 # shipped instance, 512 substeps
    assert all(rx[k] == xr[k] for k in ("M", "S", "nsub", "dps")), "xcheck instance"
    assert rx["worst_rel_dev"] < 1e-30 and \
        abs(rx["worst_rel_dev"] / xr["worst_rel_dev"] - 1.0) <= 1e-3, \
        f"xcheck {rx['worst_rel_dev']:.6e} vs shipped {xr['worst_rel_dev']:.6e}"
    return (f"M=200 S=-5 cell reproduced: mid30 19/19 identical, balls overlap, separations "
            f"identical (jM {rec['moment_diag']['bracket_sep_mid_jM']}); xcheck "
            f"{rx['worst_rel_dev']:.6e} bit-identical {rx['worst_rel_dev'] == xr['worst_rel_dev']}")


@leg("planted_regen")
def leg_planted_regen(scratch):
    for m in ("numpy", "msprime", "tskit"):
        if not _has(m):
            return f"SKIP:{m} absent — pip install msprime (simulation legs only)"
    import json
    import math
    refd = os.path.join(HERE, "reference", "planted")
    mode, drift = _planted_mode()
    out = os.path.join(scratch, "planted")
    os.makedirs(out, exist_ok=True)
    eng = os.path.join(HERE, "synth_truths.py")
    parts = []
    for cls in ("B13", "D020"):
        ref = json.load(open(os.path.join(refd, f"truth_{cls}.json")))
        reps = ref["mc"]["n4"]["reps"]
        r = subprocess.run([PY, eng, "truth", "--cls", cls, "--out", out,
                            "--reps", str(reps)], stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT, text=True, timeout=600)
        assert r.returncode == 0, f"rc={r.returncode}\n{r.stdout[-800:]}"
        got = json.load(open(os.path.join(out, f"truth_{cls}.json")))
        dq = zq = 0.0
        for nk in ("n4", "n6"):
            b, g = ref["mc"][nk], got["mc"][nk]
            assert b["seed"] == g["seed"] and b["reps"] == g["reps"], f"{cls} {nk} seed/reps"
            for key in ("q_ij", "q_ij_mc_se"):
                for xr, yr, sr in zip(b[key], g[key], b["q_ij_mc_se"]):
                    for x, y, s in zip(xr, yr, sr):
                        dq = max(dq, abs(x - y))
                        zq = max(zq, abs(x - y) / (s * math.sqrt(2)))
            for x, y in zip(b["xi_1sfs"] + b["xi_1sfs_mc_se"], g["xi_1sfs"] + g["xi_1sfs_mc_se"]):
                dq = max(dq, abs(x - y))
        if mode == "BIT-IDENTICAL":
            assert dq <= 1e-12, f"{cls}: regenerated truth differs, max|d| {dq:.2g}"
        else:
            assert zq <= 4.0, f"{cls}: regenerated truth off by {zq:.1f} SE"
        parts.append(f"{cls}@{reps} max|d| {dq:.1e}")
    ref = json.load(open(os.path.join(refd, "arm_TK1_rON_GC.json")))
    r = subprocess.run([PY, eng, "arm", "--arm", "TK1_rON_GC", "--out", out],
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                       text=True, timeout=600)
    assert r.returncode == 0, f"rc={r.returncode}\n{r.stdout[-800:]}"
    got = json.load(open(os.path.join(out, "arm_TK1_rON_GC.json")))
    assert got["seeds"] == ref["seeds"], "arm seeds"
    assert os.path.isfile(got["files"]["gt"]["path"]) and got["files"]["gt"]["path"].startswith(out)
    same = {k: got[k] == ref[k] for k in ("params", "counts")}
    same["payload"] = all(got["files"][k]["payload_sha256"] == ref["files"][k]["payload_sha256"]
                          for k in ("gt", "anc"))
    if mode == "BIT-IDENTICAL":
        assert all(same.values()), f"arm TK1_rON_GC differs: {same}"
    else:
        assert got["counts"]["kept_biallelic_snv"] > 0
    return (f"[{mode}] CLI regeneration into scratch: {'; '.join(parts)}; arm TK1_rON_GC "
            f"1 Mb/500 dip kept_snv {got['counts']['kept_biallelic_snv']} "
            + " ".join(f"{k}{'==' if v else '!='}" for k, v in same.items()))


def main():
    ap = argparse.ArgumentParser(description="popcorn acceptance battery")
    ap.add_argument("--full", action="store_true",
                    help="also run the heavy gates (minutes)")
    a = ap.parse_args()

    if not _has("mpmath"):
        print("FAIL: mpmath is required for every leg — pip install mpmath")
        return 1

    t0 = time.time()
    for fn in (leg_pins, leg_cli_selftest, leg_cli_check, leg_sfs_exact,
               leg_sfs_arb, leg_dominance_oracle, leg_dominance_qseries,
               leg_dfe_kernel, leg_certificates_lp, leg_region_certs,
               leg_lambda_kingman, leg_lambda_identities, leg_lambda_spectra_ref,
               leg_lambda_functional, leg_twolocus_kingman, leg_twolocus_rho,
               leg_transient_stationary, leg_transient_selfconv,
               leg_ancestral_join, leg_ancestral_cli,
               leg_twosfs_identities, leg_twosfs_hull_ref, leg_twosfs_kernel,
               leg_twosfs_class_certs, leg_twosfs_verdict, leg_twosfs_montecarlo,
               leg_foldgate_map, leg_foldgate_annihilated,
               leg_foldgate_obstructions, leg_foldgate_deciders,
               leg_twosite_projection, leg_twosite_spectrum, leg_twosite_cli,
               leg_enclosure_cell, leg_enclosure_stationary, leg_enclosure_reference,
               leg_planted_analytic, leg_planted_vs_lambda, leg_planted_vs_twolocus,
               leg_planted_sim,
               leg_twowindow_laws, leg_twowindow_linked, leg_twowindow_recomb):
        run_leg(fn)
    if a.full:
        scratch = tempfile.mkdtemp(prefix="popcorn_battery_")
        try:
            for fn in (leg_gate_phase1, leg_dominance_full, leg_dfe_full,
                       leg_enclosure_full, leg_planted_regen):
                run_leg(fn, scratch)
        finally:
            shutil.rmtree(scratch, ignore_errors=True)

    n_pass = sum(1 for _, s in RESULTS if s == "PASS")
    n_skip = sum(1 for _, s in RESULTS if s == "SKIP")
    failed = [n for n, s in RESULTS if s == "FAIL"]
    print(f"{n_pass} pass, {n_skip} skip, {len(failed)} fail "
          f"({time.time() - t0:.1f}s)")
    if failed:
        print(f"OVERALL FAIL: {', '.join(failed)}")
        return 1
    print("OVERALL PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())

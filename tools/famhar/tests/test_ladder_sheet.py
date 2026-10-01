"""The ladder sheet battery for the vendored closed-form plugin `ladder_reference`.

Fixture: tests/fixtures/ladder_sheet_receipt.json, the sheet_checks record of an identity
receipt (its sha256 pinned in the fixture; every string copied from the record by object, never
retyped).  The record evaluated Phi^(3) on the real lambda > 0 sheet (z, zbar real, u + v < 1)
at (X, Y) = (1/4, 1/10) and (1/4, 1e-20) three independent ways and found the plugin's
cont=None call wrong there by one turn of ell: log(a) + log(b) on two negative reals carries an
extra 2*pi*i.  The cure is the sheet-selected default (cont=None means the principal log(a*b));
these tests pin it, keep the explicit ell_k semantics, and carry the planted controls.

Digit thresholds are set from measured agreements (the numbers in the docstrings are the
measurements at the time the thresholds were set; the asserts carry the margin).  Agreement
against the record's own value is capped by the record's own working precision (its
sheet_checks were computed at dps 40 and printed at 50 digits), so the bar against the record
is 38 digits; the 50-digit-class statements are self-consistency of the plugin at dps 50 vs 80.
"""
import copy
import hashlib
import importlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from fractions import Fraction

import pytest
from mpmath import mp, mpf, mpc

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
for p in (PKG, os.path.dirname(PKG)):
    if p not in sys.path:
        sys.path.insert(0, p)

import ladder_reference as LR  # noqa: E402
import famhar  # noqa: E402

FIXTURE_PATH = os.path.join(HERE, "fixtures", "ladder_sheet_receipt.json")
FIX = json.load(open(FIXTURE_PATH))
DPS = 50
RECORD_BAR = 38        # measured 40.30 (point 0) and 40.21-class (record's own dps-40 cap); margin 2
SELF_BAR = 48          # dps 50 vs dps 80 self-consistency, measured >= 50
IM_BOUND = mpf(10) ** -40   # measured |Im| 1.1e-68 at dps 50


# ----------------------------------------------------------------------------- helpers
def fmp(fr):
    """Exact Fraction -> mpf at the working precision (never through a float)."""
    fr = Fraction(fr)
    return mpf(fr.numerator) / fr.denominator


def rel_digits(x, y, dps_cmp=200):
    """Relative agreement digits of two numbers or decimal strings (999 = identical)."""
    with mp.workdps(dps_cmp):
        x = mpf(x) if isinstance(x, str) else x
        y = mpf(y) if isinstance(y, str) else y
        if x == y:
            return 999.0
        return float(-mp.log10(abs(x - y) / abs(y)))


def chart_point(X, Y, dps):
    """The record's chart point from exact (X, Y) on the real lambda > 0 sheet.

    lambda = (1+X-Y)^2 - 4X = (1-X+Y)^2 - 4Y exactly (Fraction).  z, zbar are the real roots of
    t^2 - (1+X-Y) t + X = 0 with z the smaller; a = w = z/(z-1), b = wbar = zbar/(zbar-1).
    Formed through delta = 1 - zbar = 2Y/((1-X+Y) + sqrt(lambda)) (the small root of
    delta^2 - (1-X+Y) delta + Y = 0), which keeps every digit when Y is tiny; the plain
    (s +- sqrt(lambda))/2 form is asserted equal at the first point below.
    """
    lam = (1 + X - Y) ** 2 - 4 * X
    assert lam == (X + Y - 1) ** 2 - 4 * X * Y
    assert lam > 0, "not on the real sheet: lambda <= 0"
    with mp.workdps(dps):
        Xm = mpf(X.numerator) / X.denominator
        Ym = mpf(Y.numerator) / Y.denominator
        r = mp.sqrt(mpf(lam.numerator) / lam.denominator)
        delta = 2 * Ym / ((1 - Xm + Ym) + r)
        zbar = 1 - delta
        z = Xm / (1 - delta)
        a = -Xm / ((1 - Xm) - delta)      # z/(z-1)
        b = -(1 - delta) / delta          # zbar/(zbar-1)
        return dict(lam=lam, z=z, zbar=zbar, a=a, b=b)


def point(i, dps=DPS):
    pt = FIX["points"][i]
    X, Y = Fraction(pt["X"]), Fraction(pt["Y"])
    assert str(chart_point(X, Y, dps)["lam"]) == pt["lambda_exact"], "lambda_exact"
    return pt, chart_point(X, Y, dps)


def check_record_point(pt, cp, L, dps, bar=RECORD_BAR):
    """The record comparisons at one fixture point; every failure names the fixture field."""
    with mp.workdps(dps):
        a, b = cp["a"], cp["b"]
        assert a < 0 and b < 0, "chart point must have both coordinates negative real"
        v_k = LR.phi_value(L, a, b, dps, {"ell_k": -1})       # the record's exact call
        v_d = LR.phi_value(L, a, b, dps, None)                # the sheet-selected default
        v_0 = LR.phi_value(L, a, b, dps, {"ell_k": 0})        # the old cont=None branch
    d_k = rel_digits(mp.re(v_k), pt["own_value_str"])
    assert d_k >= bar, f"own_value_str: ell_k=-1 call agrees to {d_k:.2f} digits < {bar}"
    assert abs(mp.im(v_k)) < IM_BOUND, f"own_value_str: |Im| of the ell_k=-1 call {abs(mp.im(v_k))}"
    d_d = rel_digits(mp.re(v_d), pt["own_value_str"])
    assert d_d >= bar, f"own_value_str: default call agrees to {d_d:.2f} digits < {bar}"
    assert abs(mp.im(v_d)) < IM_BOUND, f"own_value_str: |Im| of the default call {abs(mp.im(v_d))}"
    d_kk = rel_digits(mp.re(v_k), pt["ladder_reference_ellk_minus1_re_str"])
    assert d_kk >= bar, f"ladder_reference_ellk_minus1_re_str: {d_kk:.2f} digits < {bar}"
    d_00 = rel_digits(mp.re(v_0), pt["ladder_reference_cont_none_re_str"])
    assert d_00 >= bar, f"ladder_reference_cont_none_re_str: ell_k=0 reproduces the old branch to {d_00:.2f} digits < {bar}"
    im0 = float(abs(mp.im(v_0)))
    assert abs(im0 - pt["ladder_reference_cont_none_im_abs"]) <= 1e-9 * abs(im0), \
        f"ladder_reference_cont_none_im_abs: {im0} vs record {pt['ladder_reference_cont_none_im_abs']}"
    return dict(d_ellk_minus1_vs_own=d_k, d_default_vs_own=d_d, d_default_vs_ellk_minus1=rel_digits(mp.re(v_d), mp.re(v_k)),
                d_ellk0_vs_record_cont_none=d_00, v_k=v_k, v_d=v_d, v_0=v_0)


# ----------------------------------------------------------------------------- the fixture
def test_fixture_record_pinned():
    assert FIX["record_sha256_16"] == "3c50aeda13a5c9a2"
    assert FIX["record_sha256"].startswith(FIX["record_sha256_16"])
    assert FIX["L"] == 3 and FIX["record_dps"] == 40
    assert len(FIX["points"]) == 2
    for pt in FIX["points"]:
        assert pt["lambda_positive"] and pt["X_plus_Y_below_1"]
        assert len(pt["own_value_str"].replace(".", "").replace("-", "")) >= 50


def test_chart_point_matches_record_strings():
    """z, zbar, a = w, b = wbar reproduce the record's 30-digit strings at both points."""
    for i in range(2):
        pt, cp = point(i)
        with mp.workdps(DPS):
            for key, val in (("z_str", cp["z"]), ("zbar_str", cp["zbar"]), ("w_str", cp["a"]), ("wbar_str", cp["b"])):
                assert mp.nstr(val, 30) == pt[key], f"{key}: {mp.nstr(val, 30)} vs record {pt[key]}"
            # the plain root form agrees with the stable form at the first point
            if i == 0:
                s = fmp(Fraction(pt["X"]) + 1 - Fraction(pt["Y"]))
                r = mp.sqrt(fmp(Fraction(pt["lambda_exact"])))
                z_plain, zb_plain = (s - r) / 2, (s + r) / 2
                assert rel_digits(z_plain, cp["z"]) >= 45
                assert rel_digits(zb_plain / (zb_plain - 1), cp["b"]) >= 45


def test_point0_record_call_and_default_agree_with_record():
    """(i)+(ii): the record's exact call cont={'ell_k': -1} and the default cont=None both
    reproduce own_value_str (measured 40.30 digits, capped by the record's dps 40), Im ~ 0."""
    pt, cp = point(0)
    r = check_record_point(pt, cp, FIX["L"], DPS)
    assert r["d_default_vs_ellk_minus1"] >= SELF_BAR, r


def test_point0_selected_k_is_minus_one():
    """(iii): the sheet-selected k at the real lambda > 0 chart point is -1, reported by object."""
    pt, cp = point(0)
    assert LR.default_ell_k(cp["a"], cp["b"], DPS) == -1
    sh = LR.ell_sheet(cp["a"], cp["b"], DPS)
    assert sh["ell_k"] == -1 and sh["selected_by"] == "default"
    with mp.workdps(DPS):
        # ell is the real master ln(X/Y) = ln(u/v): the record's ell_str, and Im ell = 0
        assert mp.nstr(mp.re(sh["ell"]), 30) == pt["ell_str"]
        assert abs(mp.im(sh["ell"])) < IM_BOUND
        assert rel_digits(mp.re(sh["ell"]), mp.log(fmp(Fraction(pt["X"]) / Fraction(pt["Y"])))) >= SELF_BAR
    exp = LR.ell_sheet(cp["a"], cp["b"], DPS, {"ell_k": -1})
    assert exp["ell_k"] == -1 and exp["selected_by"] == "explicit"
    assert exp["ell"] == sh["ell"]


def test_point0_wrong_sheet_control_and_default_is_not_it():
    """(iv): cont={'ell_k': 0} is the old cont=None branch (real part -383.41..., |Im| 113.78,
    the record's ladder_reference_cont_none_* fields) and the DEFAULT does not return it."""
    pt, cp = point(0)
    r = check_record_point(pt, cp, FIX["L"], DPS)
    assert r["d_ellk0_vs_record_cont_none"] >= RECORD_BAR
    assert abs(mp.im(r["v_0"])) > 100
    d_wrong = rel_digits(mp.re(r["v_d"]), pt["ladder_reference_cont_none_re_str"])
    assert d_wrong < 1, ("DEFAULT RETURNED THE WRONG-SHEET BRANCH: cont=None reproduces "
                         f"ladder_reference_cont_none_re_str to {d_wrong:.2f} digits")
    assert abs(mp.im(r["v_d"])) < IM_BOUND, "DEFAULT RETURNED THE WRONG-SHEET BRANCH: Im != 0"


def test_point0_tampered_fixture_fails_by_name():
    """(v): one digit changed in the fixture's own_value_str fails naming the field."""
    pt, cp = point(0)
    bad = copy.deepcopy(pt)
    s = bad["own_value_str"]
    k = 20                                    # a digit well inside the 38-digit bar
    assert s[k].isdigit()
    bad["own_value_str"] = s[:k] + str((int(s[k]) + 1) % 10) + s[k + 1:]
    with pytest.raises(AssertionError, match="own_value_str"):
        check_record_point(bad, cp, FIX["L"], DPS)
    # and the untouched fixture passes the same check
    check_record_point(pt, cp, FIX["L"], DPS)


def test_point0_two_precision_self_consistency():
    """The dps-50 default value agrees with the dps-80 default value to >= 48 digits (measured >= 50)."""
    pt, cp50 = point(0, 50)
    _, cp80 = point(0, 80)
    v50 = LR.phi_value(FIX["L"], cp50["a"], cp50["b"], 50, None)
    v80 = LR.phi_value(FIX["L"], cp80["a"], cp80["b"], 80, None)
    d = rel_digits(mp.re(v50), mp.re(v80))
    assert d >= SELF_BAR, d
    assert LR.default_ell_k(cp80["a"], cp80["b"], 80) == -1


def test_point0_derivatives_vs_finite_differences():
    """(vi): phi_derivs_chart (the same default sheet) vs central differences of the default
    phi_value in a and b at dps 50; measured agreement recorded in the assertion margin."""
    pt, cp = point(0)
    L = FIX["L"]
    a, b = cp["a"], cp["b"]
    da, db = LR.phi_derivs_chart(L, a, b, DPS)
    with mp.workdps(DPS + 10):
        h = mpf(10) ** -20
        fd_a = (LR.phi_value(L, a + h, b, DPS) - LR.phi_value(L, a - h, b, DPS)) / (2 * h)
        fd_b = (LR.phi_value(L, a, b + h, DPS) - LR.phi_value(L, a, b - h, DPS)) / (2 * h)
    d_a = rel_digits(mp.re(da), mp.re(fd_a))
    d_b = rel_digits(mp.re(db), mp.re(fd_b))
    assert d_a >= 30 and d_b >= 30, (d_a, d_b)      # measured ~40 (h = 1e-20: truncation h^2, rounding 1e-60/h)
    assert abs(mp.im(da)) < IM_BOUND and abs(mp.im(db)) < IM_BOUND


def test_point1_tiny_Y_second_case():
    """(vii): the record's second point (Y = 1e-20) at dps 50 through the stable chart form
    (runs in well under a second).  Measured 36.82 digits against the record's own value, the
    same figure the record itself reports for its ell_k=-1 call there (its own dps-40 cap at a
    point where delta = 1 - zbar is 1e-20): bar 34."""
    pt, cp = point(1)
    r = check_record_point(pt, cp, FIX["L"], DPS, bar=34)
    assert LR.default_ell_k(cp["a"], cp["b"], DPS) == -1
    assert r["d_default_vs_ellk_minus1"] >= SELF_BAR


def test_euclidean_slice_control():
    """(viii): b = conj(a), Im a != 0 -> k = 0 and the default equals log(a) + log(b): no change
    of behaviour on the Euclidean slice (the calibration points of the gate battery)."""
    for re_fr, im_fr in (("1/6", "-5/6"), ("1/13", "-8/13"), ("-1/4", "-1/6")):
        with mp.workdps(DPS):
            a = mpc(fmp(re_fr), fmp(im_fr))
            b = mp.conj(a)
            assert mp.im(a) != 0
            assert LR.default_ell_k(a, b, DPS) == 0
            sh = LR.ell_sheet(a, b, DPS)
            assert sh["ell_k"] == 0 and sh["selected_by"] == "default"
            assert sh["ell"] == mp.log(a) + mp.log(b)
            v_d = LR.phi_value(4, a, b, DPS, None)
            v_0 = LR.phi_value(4, a, b, DPS, {"ell_k": 0})
            assert v_d == v_0


def test_anchor_control():
    """(ix): real 0 < a0, b0 < 1 -> k = 0; anchor_values unchanged vs the explicit ell_k = 0 basis."""
    ms = LR.masters_spec(4)
    with mp.workdps(DPS):
        for a0, b0 in ((mpf(1) / 8, mpf(1) / 9), (mpf(1) / 4, mpf(1) / 9), (mpf(1) / 2, mpf(3) / 4)):
            assert LR.default_ell_k(a0, b0, DPS) == 0
            J_anchor = LR.anchor_values({"a": a0, "b": b0}, DPS, ms)
            J_expl = LR.basis_values(a0, b0, DPS, ms, cont={"ell_k": 0})
            assert J_anchor == J_expl
            assert all(mp.im(x) == 0 for x in J_anchor)


def test_explicit_ell_k_keeps_its_meaning():
    """An explicit ell_k is applied relative to log(a) + log(b) on every sheet (k = 1 at the
    anchor shifts ell by exactly 2*pi*i, as the branch-locus gate of the battery expects)."""
    with mp.workdps(DPS):
        a0, b0 = mpf(1) / 8, mpf(1) / 9
        e0 = LR.ell_sheet(a0, b0, DPS)["ell"]
        e1 = LR.ell_sheet(a0, b0, DPS, {"ell_k": 1})["ell"]
        assert rel_digits(mp.im(e1 - e0), 2 * mp.pi) >= SELF_BAR
        pt, cp = point(0)
        e_expl = LR.ell_sheet(cp["a"], cp["b"], DPS, {"ell_k": 0})
        assert e_expl["selected_by"] == "explicit" and e_expl["ell_k"] == 0
        assert rel_digits(mp.im(e_expl["ell"]), 2 * mp.pi) >= SELF_BAR


# ----------------------------------------------------------------------------- membership
def readme_members(text):
    """Member names the README's Files section lists: backticked names with a suffix on the
    section's bullet lines (the unbulleted 'Not shipped' paragraph names records, not members)."""
    sec = text.split("## Files", 1)[1]
    names = set()
    for line in sec.splitlines():
        if line.startswith("- ") or line.startswith("  "):
            names.update(re.findall(r"`([A-Za-z0-9_./-]+\.(?:py|json|md))`", line))
    return sorted(names)


def test_readme_members_exist_in_package():
    text = open(os.path.join(PKG, "README.md")).read()
    names = readme_members(text)
    assert "ladder_reference.py" in names and "build_ladder_config.py" in names and "famhar.py" in names
    missing = [n for n in names if not os.path.exists(os.path.join(PKG, n))]
    assert not missing, f"README Files section names members absent from the package: {missing}"
    # planted control: a member the tree does not hold is caught by name
    fake = readme_members(text + "\n- `no_such_module.py` — planted\n")
    missing_fake = [n for n in fake if not os.path.exists(os.path.join(PKG, n))]
    assert missing_fake == ["no_such_module.py"]
    # gates_ladder.py is named as a record, not a member: it must not be listed as a file
    assert "gates_ladder.py" not in names


def test_guide_names_the_members_and_the_leg():
    text = open(os.path.join(PKG, "GUIDE.md")).read()
    for n in ("ladder_reference.py", "build_ladder_config.py", "families/ladder_L4.json", "B6"):
        assert n in text, n


# ----------------------------------------------------------------------------- the vendored builder
def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        h.update(fh.read())
    return h.hexdigest()


def test_provenance_pins_match_the_tree():
    prov = json.load(open(os.path.join(PKG, "PROVENANCE.json")))
    for entry in prov["vendored"]:
        path = os.path.join(PKG, entry["vendored_as"])
        assert os.path.exists(path), entry["vendored_as"]
        assert sha256_file(path) == entry["vendored_sha256"], entry["vendored_as"]
        if entry.get("changed"):
            assert entry["vendored_sha256"] != entry["record_sha256"]
        else:
            assert entry["vendored_sha256"] == entry["record_sha256"]
    assert sha256_file(FIXTURE_PATH) == prov["fixture"]["sha256"]


def test_builder_reproduces_the_committed_config(tmp_path):
    """build_ladder_config.py run from a copy of the package emits families/ladder_L4.json
    byte-identical to the committed one (the drift guard on the calibration family config)."""
    for n in ("build_ladder_config.py", "ladder_reference.py"):
        shutil.copy(os.path.join(PKG, n), tmp_path / n)
    (tmp_path / "families").mkdir()
    r = subprocess.run([sys.executable, "build_ladder_config.py"], cwd=tmp_path, capture_output=True, text=True, timeout=120)
    assert r.returncode == 0, r.stderr
    assert "n=65 masters, 156 nonzero connection entries, 4 targets" in r.stdout
    assert sha256_file(tmp_path / "families" / "ladder_L4.json") == sha256_file(os.path.join(PKG, "families", "ladder_L4.json"))


def aj_vs_fd_worst(fam, ptraw, dps=40):
    """Worst matched digit of A.J vs central finite differences of the closed-form basis at one
    chart point (the connection-validation gate of the calibration battery, one point)."""
    with mp.workdps(dps + 20):
        h = mpf(10) ** -18
        pt = fam.parse_point(ptraw, dps + 20)
        J = LR.basis_values(pt["a"], pt["b"], dps + 20, fam.masters)
        Av = fam.A_at(pt, dps + 20)
        scale = max(abs(x) for x in J)
        worst = 999.0
        for var in ("a", "b"):
            Jp = LR.basis_values(pt["a"] + (h if var == "a" else 0), pt["b"] + (h if var == "b" else 0), dps + 20, fam.masters)
            Jm = LR.basis_values(pt["a"] - (h if var == "a" else 0), pt["b"] - (h if var == "b" else 0), dps + 20, fam.masters)
            for i in range(65):
                fd = (Jp[i] - Jm[i]) / (2 * h)
                aj = sum(Av[var][i][j] * J[j] for j in range(65))
                err = abs(fd - aj) / scale
                d = 999.0 if err == 0 else float(-mp.log10(err))
                worst = min(worst, d)
    return worst


def test_config_loads_plugin_resolves_and_AJ_matches_finite_differences():
    """The committed calibration config loads in the engine, its anchor plugin resolves to the
    vendored module, and A.J at a fresh complex point matches central finite differences of the
    closed-form basis (measured worst digit recorded in the assertion margin)."""
    fam = famhar.Family.from_json(os.path.join(PKG, "families", "ladder_L4.json"))
    assert fam.n == 65 and len(fam.masters) == 65
    mod, fn = fam.cfg["anchor"]["values_plugin"].rsplit(".", 1)
    plugin = getattr(importlib.import_module(mod), fn)
    assert plugin is LR.anchor_values
    dps = 40
    anchor = fam.parse_point(fam.cfg["anchor"]["point"], dps + 15)
    y0 = plugin(anchor, dps + 15, fam.masters)
    assert len(y0) == 65 and y0[fam.masters.index({"type": "l", "r": 0})] == 1
    worst = aj_vs_fd_worst(fam, {"a": "1/3+1/7i", "b": "2/7-1/5i"}, dps)
    assert worst >= 25, worst     # measured 34.2 (h = 1e-18 at dps 60)


# ----------------------------------------------------------------------------- the registered self-test
VERDICT_RE = re.compile(r"^\[famhar selftest\] (\S+)  \(legs (\d+) executed / (\d+) skipped(?: \(([^)]*)\))?;")
BANNER_RE = re.compile(r"^\[(B\d+)\] ")


def run_selftest(pkg_dir):
    """Run `python3 famhar.py --selftest` from pkg_dir with an environment of PATH only (the
    registered line's form); return (rc, verdict token, executed, skipped, skipped names,
    banner ids in print order, stdout).  Exactly one verdict line is required."""
    env = {"PATH": os.environ.get("PATH", "")}
    r = subprocess.run([sys.executable, "famhar.py", "--selftest"], cwd=pkg_dir, env=env,
                       capture_output=True, text=True, timeout=600)
    lines = r.stdout.splitlines()
    banners = [m.group(1) for m in (BANNER_RE.match(l) for l in lines) if m]
    verdicts = [m for m in (VERDICT_RE.match(l) for l in lines) if m]
    assert len(verdicts) == 1, (r.stdout[-1500:], r.stderr[-1500:])
    m = verdicts[0]
    names = [n.strip() for n in m.group(4).split(",")] if m.group(4) else []
    return r.returncode, m.group(1), int(m.group(2)), int(m.group(3)), names, banners, r.stdout


def package_copy(tmp_path):
    """A scratch copy of the package beside a copy of the transport core it imports."""
    root = tmp_path / "tools"
    ignore = shutil.ignore_patterns("__pycache__", ".pytest_cache")
    shutil.copytree(PKG, root / "famhar", ignore=ignore)
    shutil.copytree(os.path.join(os.path.dirname(PKG), "wayfinder"), root / "wayfinder", ignore=ignore)
    return root / "famhar"


def test_selftest_verdict_counts_follow_the_legs_and_a_skip_is_never_pass_all(tmp_path):
    """The verdict line's executed + skipped count equals the number of [Bn] banners printed
    (the counts are read from the leg registry, not typed); PASS_ALL only at 0 skipped with
    rc 0.  Control: the plugin removed from a copy -> B6 is skipped by name, the token is the
    distinct PASS_EXECUTED_WITH_SKIPS and rc is 2 (PASS_ALL nowhere in the output)."""
    rc, tok, n_exec, n_skip, names, banners, out = run_selftest(PKG)
    assert (rc, tok) == (0, "PASS_ALL"), out[-1500:]
    assert (n_skip, names) == (0, [])
    assert n_exec == len(banners), (n_exec, banners)
    assert banners == [f"B{i}" for i in range(1, len(banners) + 1)] and "B6" in banners
    assert "SKIP:" not in out
    cp = package_copy(tmp_path)
    (cp / "ladder_reference.py").unlink()
    rc, tok, n_exec, n_skip, names, banners, out = run_selftest(cp)
    assert (rc, tok) == (2, "PASS_EXECUTED_WITH_SKIPS"), out[-1500:]
    assert (n_skip, names) == (1, ["B6"])
    assert n_exec + n_skip == len(banners) and "B6" in banners
    assert "  SKIP: absent ladder_reference.py" in out
    assert "PASS_ALL" not in out

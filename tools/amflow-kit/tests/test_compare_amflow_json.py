"""The battery of compare_amflow_json.py: the NUMERICAL ZERO rule and the MISSING-BY-WINDOW rule.

The reference fixture pair (tests/fixtures/pair_goal50.json, pair_goal70.json; PINS.json beside them): two solve_integrals
outputs of one finite eight-propagator integral at goals 50 and 70, whose pole-order coefficients are numerical zeros
(|re| ~ 1e-116 .. 1e-77 at goal 50, 1e-142 .. 1e-113 at goal 70) and whose eps window differs by goal (eps^-6 and eps^-5
present at goal 50 only).  pair_compare_previous.out is the reading of that pair by an earlier version of the tool, which
had neither rule: OVERALL FAIL (the two window orders as MISSING -> FAIL, the four pole-order re parts as 0.0 digits -> FAIL,
eps^0 / eps^1 / eps^2 re at their digit counts -> PASS).  The current tool reads the same pair OVERALL PASS with the four
numerical zeros and the two window orders named and counted; the eps^0 / eps^1 / eps^2 digit lines are read from that .out
and asserted equal.
Every count asserted against the reference pair is COUNTED from the fixture JSONs here (the orders present on one side only;
the components of the common orders whose two members are both below 10^-(50-2) and not both exactly zero), never typed.
The earlier tool is run as the positive control when COMPARE_AMFLOW_JSON_PREVIOUS names its file (skipped by name otherwise).
Synthetic pairs cover each rule in isolation: a numerical-zero pair below T (excluded, counted), the same pair with one member
above T (zero-vs-nonzero FAIL), the radii governing T when larger than the goal term, an im-part numerical zero, a window-only
order (not compared; --strict-window FAIL; the one-sided exact zero PASS), exact zeros (PASS, not counted as
numerical zeros), a planted digit in eps^0 (FAIL by digits), --goal G / GREF,GTEST parsing with the lower goal governing, the
goal read from a file key, a malformed --goal refused by name (rc 2), --n-integrals and --min-digits, and the exit
codes.  The tool runs as a subprocess with PATH, HOME, PYTHONDONTWRITEBYTECODE and the user site (PYTHONUSERBASE) passed through.
Synthetic coefficient maps with integer order keys are built as {**PAIR_OK, -1: ...} (a keyword mapping refuses integer keys).
A second fixture pair, tests/fixtures/pair_q1b_goal30.json / pair_q1b_goal60.json -- 98 integrals of one family
solved at goals 30 and 60 -- covers the imaginary parts: 46 im components of real coefficients are numerical zeros
that a radii-only comparison reads as `0.0 digits -> FAIL` (OVERALL FAIL) and that --goal 30,60 reads
`OVERALL: PASS (1358 compared, 46 numerical zeros excluded, 0 orders not compared)`; the counts are COUNTED from the
fixture JSONs by the numerical-zero rule and asserted equal to the reference figures below.
The goal from the CONFIG -- three sources in order: --goal; --goal-from-config REF_CFG[,TEST_CFG] reading the amflow_cli
input config's top-level goal_digits (or goal, first found); with neither flag the file's own key, else the config looked up
BY NAME beside the out (out/<name>.json <-> configs/<name>.json, every location tried named on the goal: line); no derivable
goal -> the run REFUSES by name rc 2 and compares nothing, --radii-only being the explicit opt-in to the radii-alone threshold.
Fixtures: the two configs of the reference pair vendored RE-CUT as pair_goal50.config.json / pair_goal70.config.json (goal_digits
kept, the machine paths and producer bookkeeping replaced by <re-cut> or removed, the re-cut fields listed in PINS.json) and
pair_radii_only_previous.out = the no-goal (radii-alone) reading of the pair, which --radii-only reproduces byte for
byte.  The legs that run the reference pair or a synthetic pair with no goal pass --radii-only (named at each site).
Further legs: PASS via --goal-from-config; PASS by name on a copy laid out as out/ + configs/ (one side by name + the other
absent named as tried); no goal -> rc 2 by name, nothing compared; --radii-only byte for byte and refused beside --goal;
--goal wins over a config; a config without goal_digits / non-numeric -> rc 2 by name; a missing config path -> rc 4 by
name; --help states the three sources and the refusal; the goal: line names the key read beside the path in one token
form for both sources (`--goal-from-config <path> <key> (SIDE)` / `config by name <path> <key> (SIDE)`).
"""
import hashlib
import json
import os
import re
import site
import subprocess
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
FIX = os.path.join(HERE, "fixtures")
PY = sys.executable
CLI = os.path.join(PKG, "compare_amflow_json.py")
sys.path.insert(0, PKG)
import compare_amflow_json as cmp_mod  # noqa: E402

from mpmath import mp, mpf  # noqa: E402

mp.dps = 250
G50, G70 = 50, 70                       # the goals of the fixture pair (the file names carry them)
OVERALL_RE = re.compile(r"^OVERALL: (PASS|FAIL) \((\d+) compared, (\d+) numerical zeros excluded, (\d+) orders not compared\)$")
DIGIT_LINE_RE = re.compile(r"^(int\d+ eps\^-?\d+ (?:re|im)): (-?\d+\.\d) digits -> PASS$")


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def fx(n):
    return os.path.join(FIX, n)


def base_env():
    """the subprocess env: PATH, HOME when set, PYTHONDONTWRITEBYTECODE=1 and the user site (PYTHONUSERBASE when set, else the
    interpreter's own user base read at run time) so a scratch HOME never hides the user-site mpmath."""
    env = {"PATH": os.environ.get("PATH", "/usr/bin:/bin"), "PYTHONDONTWRITEBYTECODE": "1"}
    if "HOME" in os.environ:
        env["HOME"] = os.environ["HOME"]
    ub = os.environ.get("PYTHONUSERBASE") or site.getuserbase()
    if ub:
        env["PYTHONUSERBASE"] = ub
    return env


def run(args, cli=CLI, cwd=None):
    p = subprocess.run([PY, cli] + list(args), cwd=cwd or HERE, capture_output=True, text=True, env=base_env())
    return p.returncode, p.stdout + p.stderr


def overall(out):
    """(verdict, compared, numerical zeros, orders not compared) from the OVERALL line; asserts there is exactly one."""
    lines = [l for l in out.splitlines() if l.startswith("OVERALL:")]
    assert len(lines) == 1, out
    m = OVERALL_RE.match(lines[0])
    assert m, lines[0]
    return m.group(1), int(m.group(2)), int(m.group(3)), int(m.group(4))


def digit_lines(text):
    return {m.group(1): m.group(2) for m in (DIGIT_LINE_RE.match(l) for l in text.splitlines()) if m}


def ball(mid, rad=None):
    return f"[{mid} +/- {rad}]" if rad is not None else str(mid)


def make_out(path, integrals):
    """a solve_integrals output in the reference form: integrals = [ {order: (re_str, im_str)}, ... ]."""
    result = []
    for k, coeffs in enumerate(integrals):
        result.append({"coefficients": [{"order": o, "value": {"re": re_s, "im": im_s}} for o, (re_s, im_s) in sorted(coeffs.items())],
                       "integral": {"family": "synthetic", "indices": [1] * 3 + [k]},
                       "leading_order": min(coeffs) if coeffs else 0})
    doc = {"mode": "solve_integrals", "options": {"working_pre": 100}, "result": result}
    json.dump(doc, open(path, "w"), indent=1)
    return path


PAIR_OK = {0: (ball("0.52495777678114463233299641528264604243617097224067", "1e-60"), "0"),
           1: (ball("-0.95204595894350870285622457294698144896291911395519", "1e-60"), "0")}


# ----------------------------------------------------------------------------- pins and the record pair

def test_fixtures_pinned():
    pins = json.load(open(fx("PINS.json")))
    files = {k: v for k, v in pins.items() if not k.startswith("_")}
    assert set(files) == {"pair_goal50.json", "pair_goal70.json", "pair_compare_previous.out",
                          "pair_q1b_goal30.json", "pair_q1b_goal60.json",
                          "pair_goal50.config.json", "pair_goal70.config.json", "pair_radii_only_previous.out"}
    for name, want in files.items():
        assert sha(fx(name)) == want, name
    vend = pins["_vendored_from"]
    for name in files:
        if name in CONFIGS:                                       # re-cut: the source sha differs, goal_digits kept, the fields listed
            assert vend[name]["recut"] is True and vend[name]["sha256"] != files[name], name
            assert vend[name]["recut_fields"] and vend[name]["goal_digits_kept"] == json.load(open(fx(name)))["goal_digits"], name
            assert "<re-cut>" in open(fx(name)).read()
        else:
            assert vend[name]["recut"] is False and vend[name]["sha256"] == files[name], name


CONFIGS = ("pair_goal50.config.json", "pair_goal70.config.json")     # the configs of the record pair (goal_digits 50 / 70 kept)
NO_GOAL = "no goal: pass --goal G|GREF,GTEST, --goal-from-config REF_CFG[,TEST_CFG], or --radii-only"


def config_goals():
    """the goal_digits of the two vendored configs, read from the files (the file names carry them)"""
    return tuple(json.load(open(fx(n)))["goal_digits"] for n in CONFIGS)


def after_goal_line(out):
    return "\n".join(out.splitlines()[1:])


def record_counts():
    """the expectations of the record pair COUNTED from the fixture JSONs: the orders on one side only, and the components of the
    common orders whose two members are both below T = 10^-(50-2) and not both exactly zero (the tool's rule, re-read here)."""
    a = cmp_mod.coeffs_of(json.load(open(fx("pair_goal50.json")))["result"][0])
    b = cmp_mod.coeffs_of(json.load(open(fx("pair_goal70.json")))["result"][0])
    one_side = sorted(set(a) ^ set(b))
    t = mpf(10) ** (-(min(G50, G70) - 2))
    zeros, exact, compared = [], [], []
    for o in sorted(set(a) & set(b)):
        for part, (x, y) in zip(("re", "im"), zip(a[o], b[o])):
            xm, xr = x
            ym, yr = y
            name = f"int0 eps^{o} {part}"
            if xm == 0 and ym == 0:
                exact.append(name)
                compared.append(name)
            elif abs(xm) < max(xr, yr, t) and abs(ym) < max(xr, yr, t):
                zeros.append(name)
            else:
                compared.append(name)
    return one_side, zeros, exact, compared


def test_record_pair_reads_pass_with_zeros_and_window_named():
    one_side, zeros, exact, compared = record_counts()
    assert one_side and zeros and exact and compared               # the pair exercises every rule
    prev = open(fx("pair_compare_previous.out")).read()
    assert prev.rstrip().endswith("OVERALL: FAIL")          # the earlier tool's reading
    for o in one_side:
        assert f"int0 eps^{o}: MISSING on one side, nonzero -> FAIL" in prev
    for name in zeros:
        assert f"{name}: 0.0 digits (ref=" in prev            # the previous tool's reading of each numerical zero
    prev_digits = digit_lines(prev)
    assert set(prev_digits) == set(compared) - set(exact)     # the previous PASS digit lines are the compared re parts

    rc, out = run([fx("pair_goal50.json"), fx("pair_goal70.json"), "--goal", f"{G50},{G70}"])
    assert rc == 0, out
    assert out.splitlines()[0] == (f"goal: REF {G50}, TEST {G70} (from --goal); zero threshold T = max(radii, "
                                   f"10^-({min(G50, G70)}-2) = {cmp_mod.fmt(mpf(10) ** (-(min(G50, G70) - 2)))})")
    verdict, n_cmp, n_zero, n_win = overall(out)
    assert (verdict, n_cmp, n_zero, n_win) == ("PASS", len(compared), len(zeros), len(one_side))
    assert f"numerical zeros: {len(zeros)} ({', '.join(zeros)})" in out
    assert f"orders not compared: {len(one_side)} ({', '.join('int0 eps^%d' % o for o in one_side)})" in out
    for o in one_side:
        assert f"int0 eps^{o}: not compared (present at goal {G50} (REF) only; outside the other member's window)" in out
    for name in zeros:
        line = [l for l in out.splitlines() if l.startswith(name + ":")]
        assert len(line) == 1 and "NUMERICAL ZERO (|ref| " in line[0] and "-> not compared" in line[0], line
        assert f"below T {cmp_mod.fmt(mpf(10) ** (-(min(G50, G70) - 2)))}" in line[0]
    for name in exact:
        assert f"{name}: both exactly 0 -> PASS" in out
    assert digit_lines(out) == prev_digits                    # eps^0 / eps^1 / eps^2 re at the digits the previous reading shows
    assert "MISSING" not in out and "FAIL" not in out

    rc30, out30 = run([fx("pair_goal50.json"), fx("pair_goal70.json"), "--goal", f"{G50},{G70}", "--min-digits", "30"])
    assert rc30 == 0 and overall(out30)[0] == "PASS"


def test_record_pair_radii_only_reproduces_the_previous_reading_byte_for_byte():
    """the pair with no goal under --radii-only = the explicit opt-in to the radii-alone threshold;
    the output is the reference no-goal reading, asserted byte for byte"""
    one_side, zeros, exact, compared = record_counts()
    rc, out = run([fx("pair_goal50.json"), fx("pair_goal70.json"), "--radii-only"])
    assert rc == 1
    assert out == open(fx("pair_radii_only_previous.out")).read()   # byte for byte the previous reading
    assert out.splitlines()[0] == "goal: none (zero threshold = the two radii alone)"
    assert overall(out) == ("FAIL", len(compared) + len(zeros), 0, len(one_side))
    for name in zeros:
        assert f"{name}: 0.0 digits (ref=" in out               # above the radii: compared, and 0.0 digits as before
    for o in one_side:
        assert f"int0 eps^{o}: not compared (present in REF only; outside the other member's window)" in out
    # --radii-only beside a goal source is refused by name (argparse, rc 2)
    for extra in (["--goal", "50"], ["--goal-from-config", fx(CONFIGS[0])]):
        rc_x, out_x = run([fx("pair_goal50.json"), fx("pair_goal70.json"), "--radii-only"] + extra)
        assert rc_x == 2 and "--radii-only cannot be combined with --goal or --goal-from-config" in out_x, out_x
        assert "OVERALL" not in out_x and "Traceback" not in out_x


def test_record_pair_reads_pass_via_goal_from_config():
    g50, g70 = config_goals()
    assert (g50, g70) == (G50, G70)
    rc_g, out_g = run([fx("pair_goal50.json"), fx("pair_goal70.json"), "--goal", f"{G50},{G70}"])
    rc, out = run([fx("pair_goal50.json"), fx("pair_goal70.json"), "--goal-from-config", f"{fx(CONFIGS[0])},{fx(CONFIGS[1])}"])
    assert rc == 0, out
    assert out.splitlines()[0] == (f"goal: REF {g50}, TEST {g70} (from --goal-from-config {fx(CONFIGS[0])} goal_digits (REF), "
                                   f"--goal-from-config {fx(CONFIGS[1])} goal_digits (TEST)); zero threshold T = max(radii, "
                                   f"10^-({min(g50, g70)}-2) = {cmp_mod.fmt(mpf(10) ** (-(min(g50, g70) - 2)))})")
    assert after_goal_line(out) == after_goal_line(out_g)      # the same reading as --goal 50,70 below the goal: line
    assert overall(out) == overall(out_g) and overall(out)[0] == "PASS"
    # one config path applies to both files
    rc1, out1 = run([fx("pair_goal50.json"), fx("pair_goal70.json"), "--goal-from-config", fx(CONFIGS[0])])
    assert rc1 == 0 and out1.splitlines()[0].startswith(
        f"goal: REF {g50}, TEST {g50} (from --goal-from-config {fx(CONFIGS[0])} goal_digits (REF), "
        f"--goal-from-config {fx(CONFIGS[0])} goal_digits (TEST))")


def test_config_carrying_goal_instead_of_goal_digits_is_read_and_the_key_named(tmp_path):
    """the job-config form that carries `goal` (the smoke fixtures' form) is read as the goal too; the key is named"""
    p = tmp_path / "job.json"; json.dump({"mode": "solve_integrals", "goal": 60, "eps_order": 12}, open(p, "w"))
    rc, out = run([fx("pair_goal50.json"), fx("pair_goal70.json"), "--goal-from-config", str(p)])
    assert rc == 0, out
    assert out.splitlines()[0] == (f"goal: REF 60, TEST 60 (from --goal-from-config {p} goal (REF), --goal-from-config {p} goal (TEST)); "
                                   "zero threshold T = max(radii, 10^-(60-2) = 1.0e-58)")
    # goal_digits wins when both keys are present
    p2 = tmp_path / "both.json"; json.dump({"goal_digits": 50, "goal": 60}, open(p2, "w"))
    rc2, out2 = run([fx("pair_goal50.json"), fx("pair_goal70.json"), "--goal-from-config", str(p2)])
    assert rc2 == 0 and out2.splitlines()[0].startswith(f"goal: REF 50, TEST 50 (from --goal-from-config {p2} goal_digits (REF)")


def test_record_pair_reads_pass_by_name_beside_the_out(tmp_path):
    """a copy laid out as out/<name>.json + configs/<name>.json: no flag; the config is found by name and named on the goal: line"""
    import shutil
    g50, g70 = config_goals()
    (tmp_path / "out").mkdir(); (tmp_path / "configs").mkdir()
    for src, cfg, name in (("pair_goal50.json", CONFIGS[0], "solve_g50.json"), ("pair_goal70.json", CONFIGS[1], "solve_g70.json")):
        shutil.copyfile(fx(src), tmp_path / "out" / name); shutil.copyfile(fx(cfg), tmp_path / "configs" / name)
    o50, o70 = str(tmp_path / "out" / "solve_g50.json"), str(tmp_path / "out" / "solve_g70.json")
    c50, c70 = str(tmp_path / "configs" / "solve_g50.json"), str(tmp_path / "configs" / "solve_g70.json")
    rc_g, out_g = run([o50, o70, "--goal", f"{G50},{G70}"])
    rc, out = run([o50, o70])
    assert rc == 0, out
    assert out.splitlines()[0] == (f"goal: REF {g50}, TEST {g70} (from config by name {c50} goal_digits (REF), config by name "
                                   f"{c70} goal_digits (TEST)); zero threshold T = max(radii, 10^-({min(g50, g70)}-2) = "
                                   f"{cmp_mod.fmt(mpf(10) ** (-(min(g50, g70) - 2)))})")
    assert after_goal_line(out) == after_goal_line(out_g) and overall(out)[0] == "PASS"
    # the config of one side absent: that location is named as tried, the other side's goal governs
    (tmp_path / "configs" / "solve_g70.json").rename(tmp_path / "configs" / "solve_g70.json.aside")
    rc2, out2 = run([o50, o70])
    assert rc2 == 0, out2
    assert out2.splitlines()[0] == (f"goal: REF {g50}, TEST none (from config by name {c50} goal_digits (REF); tried by name: "
                                    f"{c70} (TEST, absent)); zero threshold T = max(radii, 10^-({g50}-2) = "
                                    f"{cmp_mod.fmt(mpf(10) ** (-(g50 - 2)))})")
    assert "int0 eps^-6: not compared (present at goal 50 (REF) only; outside the other member's window)" in out2
    # --goal wins over the config by name
    rc3, out3 = run([o50, o70, "--goal", "20"])
    assert rc3 == 0 and out3.splitlines()[0].startswith("goal: REF 20, TEST 20 (from --goal); ")


def test_no_goal_refuses_by_name_rc2_and_compares_nothing():
    """the fixture pair with no flag: no key in the outs, no configs/ beside tests/fixtures -> rc 2 by name; nothing compared"""
    rc, out = run([fx("pair_goal50.json"), fx("pair_goal70.json")])
    assert rc == 2, out
    tried = [cmp_mod.config_by_name(fx(n)) for n in ("pair_goal50.json", "pair_goal70.json")]
    assert not any(os.path.exists(t) for t in tried)
    assert f"goal: none (tried by name: {tried[0]} (REF, absent), {tried[1]} (TEST, absent))" in out
    assert NO_GOAL in out and "Traceback" not in out
    assert "OVERALL" not in out and "int0 " not in out and "numerical zeros" not in out


def test_goal_wins_over_a_config():
    rc, out = run([fx("pair_goal50.json"), fx("pair_goal70.json"), "--goal", "20", "--goal-from-config", fx(CONFIGS[0])])
    assert rc == 0, out
    assert out.splitlines()[0] == ("goal: REF 20, TEST 20 (from --goal (--goal-from-config not read: --goal wins)); "
                                   "zero threshold T = max(radii, 10^-(20-2) = 1.0e-18)")
    # a config that would be refused is not even read when --goal is given
    rc2, out2 = run([fx("pair_goal50.json"), fx("pair_goal70.json"), "--goal", "20", "--goal-from-config", fx("pair_goal50.json")])
    assert rc2 == 0 and out2.splitlines()[0] == out.splitlines()[0]


def test_config_without_goal_digits_or_non_numeric_refused_by_name_rc2(tmp_path):
    base = json.load(open(fx(CONFIGS[0])))
    cases = []
    d = dict(base); del d["goal_digits"]; cases.append(("nokey", d, "carries none of goal_digits, goal at its top level"))
    for tag, val in (("str", "50"), ("bool", True), ("zero", 0), ("neg", -5), ("none", None), ("list", [50])):
        d = dict(base); d["goal_digits"] = val; cases.append((tag, d, f"carries goal_digits = {val!r}, not a positive number"))
    for tag, doc, reason in cases:
        p = tmp_path / f"cfg_{tag}.json"; json.dump(doc, open(p, "w"))
        rc, out = run([fx("pair_goal50.json"), fx("pair_goal70.json"), "--goal-from-config", str(p)])
        assert rc == 2, (tag, out)
        assert f"no goal: config {p} {reason} (REF, --goal-from-config)" in out, (tag, out)
        assert "OVERALL" not in out and "Traceback" not in out
    # the TEST side named when only it is unusable
    p = tmp_path / "cfg_nokey.json"
    rc, out = run([fx("pair_goal50.json"), fx("pair_goal70.json"), "--goal-from-config", f"{fx(CONFIGS[0])},{p}"])
    assert rc == 2 and f"no goal: config {p} carries none of goal_digits, goal at its top level (TEST, --goal-from-config)" in out
    # not JSON
    p2 = tmp_path / "cfg_bad.json"; open(p2, "w").write("{not json")
    rc, out = run([fx("pair_goal50.json"), fx("pair_goal70.json"), "--goal-from-config", str(p2)])
    assert rc == 2 and f"no goal: config {p2} is not JSON" in out and "Traceback" not in out
    # the same refusal for a config found BY NAME beside the out
    import shutil
    (tmp_path / "out").mkdir(); (tmp_path / "configs").mkdir()
    shutil.copyfile(fx("pair_goal50.json"), tmp_path / "out" / "s.json"); shutil.copyfile(p, tmp_path / "configs" / "s.json")
    rc, out = run([str(tmp_path / "out" / "s.json"), fx("pair_goal70.json")])
    assert rc == 2, out
    assert f"goal: none (REF: config by name {tmp_path / 'configs' / 's.json'} found; config {tmp_path / 'configs' / 's.json'} carries none of goal_digits, goal at its top level)" in out
    assert "(REF, config by name)" in out and "OVERALL" not in out


def test_missing_config_file_refused_by_name_rc4(tmp_path):
    missing = str(tmp_path / "absent.json")
    rc, out = run([fx("pair_goal50.json"), fx("pair_goal70.json"), "--goal-from-config", missing])
    assert rc == 4, out
    assert f"no config: {missing} (REF, --goal-from-config) does not exist" in out
    assert "OVERALL" not in out and "Traceback" not in out
    rc2, out2 = run([fx("pair_goal50.json"), fx("pair_goal70.json"), "--goal-from-config", f"{fx(CONFIGS[0])},{missing}"])
    assert rc2 == 4 and f"no config: {missing} (TEST, --goal-from-config) does not exist" in out2
    # an empty member of the pair is an argparse refusal (rc 2) by name
    rc3, out3 = run([fx("pair_goal50.json"), fx("pair_goal70.json"), "--goal-from-config", f"{fx(CONFIGS[0])},"])
    assert rc3 == 2 and "expected REF_CFG or REF_CFG,TEST_CFG as file paths" in out3


def test_help_states_the_three_sources_and_the_refusal():
    rc, out = run(["--help"])
    assert rc == 0
    for s in ("--goal G", "--goal-from-config REF_CFG[,TEST_CFG]", "--radii-only", "goal_digits", '"goal" when a config',
              "out/<name>.json <-> configs/<name>.json", NO_GOAL, "exits 4 by name", "exits 2 by name"):
        assert s in out, s


def test_record_pair_strict_window_fails_as_before():
    one_side, zeros, exact, compared = record_counts()
    rc, out = run([fx("pair_goal50.json"), fx("pair_goal70.json"), "--goal", f"{G50},{G70}", "--strict-window"])
    assert rc == 1
    for o in one_side:
        assert f"int0 eps^{o}: MISSING on one side, nonzero -> FAIL" in out
    assert overall(out) == ("FAIL", len(compared), len(zeros), 0)
    assert f"orders missing on one side (--strict-window): {len(one_side)} (" in out


def test_record_pair_on_the_previous_tool_is_the_positive_control():
    prev_cli = os.environ.get("COMPARE_AMFLOW_JSON_PREVIOUS")
    if not prev_cli:
        pytest.skip("COMPARE_AMFLOW_JSON_PREVIOUS unset: the previous tool's file is not named")
    assert os.path.isfile(prev_cli), prev_cli
    rc, out = run([fx("pair_goal50.json"), fx("pair_goal70.json")], cli=prev_cli)
    assert rc == 1
    assert out == open(fx("pair_compare_previous.out")).read()   # byte for byte the reference reading


# ----------------------------------------------------------------------------- parse_arb keeps the radius

def test_parse_arb_returns_mid_and_radius():
    assert cmp_mod.parse_arb("[1.5 +/- 2e-10]") == (mpf("1.5"), mpf("2e-10"))
    assert cmp_mod.parse_arb("[-3e-100 +/- 4.57e-407]") == (mpf("-3e-100"), mpf("4.57e-407"))
    assert cmp_mod.parse_arb("[+/- 3e-5]") == (mpf(0), mpf("3e-5"))
    assert cmp_mod.parse_arb("[]") == (mpf(0), mpf(0))
    assert cmp_mod.parse_arb("[1.25]") == (mpf("1.25"), mpf(0))
    assert cmp_mod.parse_arb("0") == (mpf(0), mpf(0))
    assert cmp_mod.parse_arb("1.5") == (mpf("1.5"), mpf(0))
    assert cmp_mod.parse_arb("[1.5 +/- -2e-10]")[1] == mpf("2e-10")    # a radius is a magnitude


# ----------------------------------------------------------------------------- synthetic pairs

def test_numerical_zero_pair_below_threshold_is_excluded_and_counted(tmp_path):
    a = make_out(tmp_path / "a.json", [{**PAIR_OK, -1: (ball("1e-60", "1e-80"), "0")}])
    b = make_out(tmp_path / "b.json", [{**PAIR_OK, -1: (ball("3e-70", "1e-90"), "0")}])
    rc, out = run([a, b, "--goal", "40"])
    assert rc == 0, out
    assert "int0 eps^-1 re: NUMERICAL ZERO (|ref| 1.0e-60, |test| 3.0e-70 below T 1.0e-38) -> not compared" in out
    assert "numerical zeros: 1 (int0 eps^-1 re)" in out
    assert overall(out) == ("PASS", 5, 1, 0)                      # eps^-1 im (exact zeros) + eps^0 re/im + eps^1 re/im


def test_one_member_above_threshold_is_zero_vs_nonzero_fail(tmp_path):
    a = make_out(tmp_path / "a.json", [{**PAIR_OK, -1: (ball("1e-60", "1e-80"), "0")}])
    b = make_out(tmp_path / "b.json", [{**PAIR_OK, -1: (ball("1e-30", "1e-80"), "0")}])
    rc, out = run([a, b, "--goal", "40"])
    assert rc == 1
    line = [l for l in out.splitlines() if l.startswith("int0 eps^-1 re:")][0]
    assert line.endswith("zero-vs-nonzero (|ref| below T 1.0e-38, |test| above) -> FAIL") and " digits (ref=" in line
    assert "numerical zeros: 0 ()" in out
    assert overall(out) == ("FAIL", 6, 0, 0)
    # the other way round names the sides the other way
    rc2, out2 = run([b, a, "--goal", "40"])
    assert rc2 == 1 and "zero-vs-nonzero (|test| below T 1.0e-38, |ref| above) -> FAIL" in out2


def test_radii_govern_the_threshold_when_larger_than_the_goal_term(tmp_path):
    a = make_out(tmp_path / "a.json", [{**PAIR_OK, -1: (ball("2e-25", "1e-20"), "0")}])
    b = make_out(tmp_path / "b.json", [{**PAIR_OK, -1: (ball("-7e-26", "3e-22"), "0")}])
    rc, out = run([a, b, "--goal", "100"])                          # 10^-(100-2) is far below the radii
    assert rc == 0, out
    assert "int0 eps^-1 re: NUMERICAL ZERO (|ref| 2.0e-25, |test| 7.0e-26 below T 1.0e-20) -> not compared" in out
    # the same mids as plain numbers (radius 0) with no goal: compared by digits, as before
    a2 = make_out(tmp_path / "a2.json", [{**PAIR_OK, -1: ("2e-25", "0")}])
    b2 = make_out(tmp_path / "b2.json", [{**PAIR_OK, -1: ("-7e-26", "0")}])
    rc2, out2 = run([a2, b2, "--radii-only"])                       # --radii-only: no goal is derivable here
    assert rc2 == 1 and "goal: none (zero threshold = the two radii alone)" in out2
    assert re.search(r"^int0 eps\^-1 re: -?\d+\.\d digits \(ref=.*\) -> FAIL$", out2, re.M), out2
    assert "zero-vs-nonzero" not in out2 and "NUMERICAL ZERO" not in out2


def test_im_parts_follow_the_same_rules(tmp_path):
    a = make_out(tmp_path / "a.json", [{0: (PAIR_OK[0][0], ball("4e-50", "1e-70"))}])
    b = make_out(tmp_path / "b.json", [{0: (PAIR_OK[0][0], ball("-1e-55", "1e-70"))}])
    rc, out = run([a, b, "--goal", "40"])
    assert rc == 0, out
    assert "int0 eps^0 im: NUMERICAL ZERO (|ref| 4.0e-50, |test| 1.0e-55 below T 1.0e-38) -> not compared" in out
    assert "numerical zeros: 1 (int0 eps^0 im)" in out
    assert overall(out) == ("PASS", 1, 1, 0)
    b2 = make_out(tmp_path / "b2.json", [{0: (PAIR_OK[0][0], ball("-1e-20", "1e-70"))}])
    rc2, out2 = run([a, b2, "--goal", "40"])
    assert rc2 == 1 and "int0 eps^0 im: " in out2 and "zero-vs-nonzero (|ref| below T 1.0e-38, |test| above) -> FAIL" in out2


def test_window_only_order_is_not_compared_unless_strict(tmp_path):
    a = make_out(tmp_path / "a.json", [{**PAIR_OK, -2: (ball("1.5e-3", "1e-80"), "0")}])
    b = make_out(tmp_path / "b.json", [PAIR_OK])
    rc, out = run([a, b, "--goal", "40"])
    assert rc == 0, out
    assert "int0 eps^-2: not compared (present at goal 40 (REF) only; outside the other member's window)" in out
    assert "orders not compared: 1 (int0 eps^-2)" in out
    assert overall(out) == ("PASS", 4, 0, 1)
    rc_t, out_t = run([b, a, "--radii-only"])                        # on the TEST side, no goal known (--radii-only)
    assert rc_t == 0 and "int0 eps^-2: not compared (present in TEST only; outside the other member's window)" in out_t
    rc_s, out_s = run([a, b, "--goal", "40", "--strict-window"])
    assert rc_s == 1
    assert "int0 eps^-2: MISSING on one side, nonzero -> FAIL" in out_s
    assert "orders missing on one side (--strict-window): 1 (int0 eps^-2)" in out_s
    assert overall(out_s) == ("FAIL", 4, 0, 0)
    # the one-sided order exactly zero: PASS as before under --strict-window, not compared without it
    a0 = make_out(tmp_path / "a0.json", [{**PAIR_OK, -2: ("0", "0")}])
    rc_z, out_z = run([a0, b, "--goal", "40", "--strict-window"])
    assert rc_z == 0 and "int0 eps^-2: only one side, but zero -> PASS" in out_z
    rc_z2, out_z2 = run([a0, b, "--goal", "40"])
    assert rc_z2 == 0 and "int0 eps^-2: not compared (present at goal 40 (REF) only" in out_z2


def test_exact_zero_pair_passes_and_is_not_a_numerical_zero(tmp_path):
    a = make_out(tmp_path / "a.json", [{**PAIR_OK, -1: ("0", "0")}])
    b = make_out(tmp_path / "b.json", [{**PAIR_OK, -1: ("[0 +/- 1e-50]", "0")}])
    rc, out = run([a, b, "--goal", "40"])
    assert rc == 0, out
    assert "int0 eps^-1 re: both exactly 0 -> PASS" in out and "int0 eps^-1 im: both exactly 0 -> PASS" in out
    assert "numerical zeros: 0 ()" in out
    assert overall(out) == ("PASS", 6, 0, 0)
    # exact zero against a value below T is a numerical zero; against a value above T a zero-vs-nonzero FAIL
    c = make_out(tmp_path / "c.json", [{**PAIR_OK, -1: (ball("1e-45", "1e-80"), "0")}])
    rc_c, out_c = run([a, c, "--goal", "40"])
    assert rc_c == 0 and "int0 eps^-1 re: NUMERICAL ZERO (|ref| 0.0, |test| 1.0e-45 below T 1.0e-38) -> not compared" in out_c
    d = make_out(tmp_path / "d.json", [{**PAIR_OK, -1: (ball("1e-5", "1e-80"), "0")}])
    rc_d, out_d = run([a, d, "--goal", "40"])
    assert rc_d == 1 and "int0 eps^-1 re: 0.0 digits (ref=0.0, test=" in out_d and "zero-vs-nonzero" in out_d


def test_planted_digit_in_eps0_fails_by_digits(tmp_path):
    a = make_out(tmp_path / "a.json", [PAIR_OK])
    planted = dict(PAIR_OK)
    planted[0] = (ball("0.52495777678114463233299641528264604243617097224068", "1e-60"), "0")   # last digit 7 -> 8
    b = make_out(tmp_path / "b.json", [planted])
    rc, out = run([a, b, "--goal", "40"])
    assert rc == 0, out                                              # ~50 digits agree: above the default bar 25
    m = re.search(r"^int0 eps\^0 re: (\d+\.\d) digits -> PASS$", out, re.M)
    assert m and 49 <= float(m.group(1)) <= 51, out
    rc2, out2 = run([a, b, "--goal", "40", "--min-digits", "60"])
    assert rc2 == 1
    assert re.search(r"^int0 eps\^0 re: \d+\.\d digits \(ref=.*, test=.*\) -> FAIL$", out2, re.M), out2
    assert overall(out2)[0] == "FAIL"


def test_goal_parsing_one_value_or_two_the_lower_governing(tmp_path):
    a = make_out(tmp_path / "a.json", [PAIR_OK])
    for flag, want in (("60", "goal: REF 60, TEST 60 (from --goal); zero threshold T = max(radii, 10^-(60-2) = 1.0e-58)"),
                       ("50,70", "goal: REF 50, TEST 70 (from --goal); zero threshold T = max(radii, 10^-(50-2) = 1.0e-48)"),
                       ("70,50", "goal: REF 70, TEST 50 (from --goal); zero threshold T = max(radii, 10^-(50-2) = 1.0e-48)"),
                       (" 40 , 45 ", "goal: REF 40, TEST 45 (from --goal); zero threshold T = max(radii, 10^-(40-2) = 1.0e-38)")):
        rc, out = run([a, a, "--goal", flag])
        assert rc == 0 and out.splitlines()[0] == want, (flag, out)


def test_goal_read_from_a_file_key_when_present(tmp_path):
    a = make_out(tmp_path / "a.json", [{**PAIR_OK, -1: (ball("1e-45", "1e-80"), "0")}])
    b = make_out(tmp_path / "b.json", [{**PAIR_OK, -1: (ball("1e-45", "1e-80"), "0")}])
    da, db = json.load(open(a)), json.load(open(b))
    da["goal_digits"] = 40
    db["options"]["goal"] = 60
    json.dump(da, open(a, "w")); json.dump(db, open(b, "w"))
    rc, out = run([a, b])
    assert rc == 0, out
    assert out.splitlines()[0] == ("goal: REF 40, TEST 60 (from file key goal_digits (REF), options.goal (TEST)); "
                                   "zero threshold T = max(radii, 10^-(40-2) = 1.0e-38)")
    assert "int0 eps^-1 re: NUMERICAL ZERO (|ref| 1.0e-45, |test| 1.0e-45 below T 1.0e-38) -> not compared" in out
    # the flag overrides the keys
    rc2, out2 = run([a, b, "--goal", "20"])
    assert rc2 == 0 and out2.splitlines()[0].startswith("goal: REF 20, TEST 20 (from --goal)")
    # one side only carries a key: the other reads none, the one governs
    dc = json.load(open(b)); del dc["options"]["goal"]; json.dump(dc, open(b, "w"))
    rc3, out3 = run([a, b])
    assert rc3 == 0 and out3.splitlines()[0].startswith(          # the side without a key is looked up by name; the location tried is named
        f"goal: REF 40, TEST none (from file key goal_digits (REF); tried by name: {cmp_mod.config_by_name(str(b))} (TEST, absent))")
    # the reference fixture pair carries none of the keys
    for name in ("pair_goal50.json", "pair_goal70.json"):
        assert cmp_mod.goal_of_file(json.load(open(fx(name)))) == (None, None)


@pytest.mark.parametrize("bad", ["abc", "50,", ",50", "0", "-5", "50,70,90", "5.5", "50;70", ""])
def test_malformed_goal_refused_by_name_rc2(tmp_path, bad):
    a = make_out(tmp_path / "a.json", [PAIR_OK])
    rc, out = run([a, a, "--goal", bad])
    assert rc == 2, (bad, out)
    assert f"argument --goal: expected G or GREF,GTEST as positive integers, got {bad!r}" in out
    assert "Traceback" not in out


def test_n_integrals_and_min_digits_unchanged(tmp_path):
    a = make_out(tmp_path / "a.json", [PAIR_OK, PAIR_OK])
    bad_second = dict(PAIR_OK); bad_second[0] = (ball("0.6", "1e-60"), "0")
    b = make_out(tmp_path / "b.json", [PAIR_OK, bad_second])
    rc, out = run([a, b, "--goal", "40"])
    assert rc == 1 and "int1 eps^0 re: " in out and overall(out)[0] == "FAIL"
    rc1, out1 = run([a, b, "--goal", "40", "--n-integrals", "1"])
    assert rc1 == 0 and "int1 " not in out1 and overall(out1) == ("PASS", 4, 0, 0)
    rc2, out2 = run([a, a, "--goal", "40", "--min-digits", "300"])   # identical values read 250.0 digits, below 300
    assert rc2 == 1 and "250.0 digits (ref=" in out2 and overall(out2)[0] == "FAIL"


def test_exit_codes():
    """0 iff every compared component >= --min-digits; 1 otherwise; 2 for a refused option, no derivable goal or a config without
    a usable goal_digits; 4 for a --goal-from-config path that does not exist."""
    rc0, _ = run([fx("pair_goal50.json"), fx("pair_goal70.json"), "--goal", f"{G50},{G70}"])
    rc1, _ = run([fx("pair_goal50.json"), fx("pair_goal70.json"), "--goal", f"{G50},{G70}", "--min-digits", "80"])
    rc2, _ = run([fx("pair_goal50.json"), fx("pair_goal70.json"), "--goal", "x"])
    rc2b, _ = run([fx("pair_goal50.json"), fx("pair_goal70.json")])
    rc4, _ = run([fx("pair_goal50.json"), fx("pair_goal70.json"), "--goal-from-config", fx("absent.config.json")])
    assert (rc0, rc1, rc2, rc2b, rc4) == (0, 1, 2, 2, 4)

# ----------------------------------------------------------------------------- the second fixture pair: Q1b (goals 30 / 60)

G30, G60 = 30, 60                       # the goals of the second pair (the file names carry them)
Q1B_N_INTEGRALS = 98
Q1B_IM_ZEROS = 46                     # im components of real coefficients that are numerical zeros (the previous tool's FAIL lines)
Q1B_COMPARED = 1358
Q1B_PREVIOUS_RE_FAIL = 0
Q1B_IM_ZERO_NAMES = (
    'int2 eps^-1 im',
    'int4 eps^-1 im',
    'int7 eps^-1 im',
    'int8 eps^-2 im',
    'int9 eps^-2 im',
    'int10 eps^-2 im',
    'int11 eps^-2 im',
    'int13 eps^-2 im',
    'int13 eps^-1 im',
    'int14 eps^-2 im',
    'int14 eps^-1 im',
    'int15 eps^-2 im',
    'int15 eps^-1 im',
    'int18 eps^-2 im',
    'int18 eps^-1 im',
    'int19 eps^-2 im',
    'int19 eps^-1 im',
    'int20 eps^-2 im',
    'int20 eps^-1 im',
    'int21 eps^-2 im',
    'int22 eps^-2 im',
    'int23 eps^-2 im',
    'int24 eps^-2 im',
    'int61 eps^-4 im',
    'int62 eps^-4 im',
    'int63 eps^-4 im',
    'int64 eps^-4 im',
    'int71 eps^-4 im',
    'int72 eps^-4 im',
    'int73 eps^-4 im',
    'int74 eps^-4 im',
    'int83 eps^-4 im',
    'int84 eps^-4 im',
    'int85 eps^-4 im',
    'int86 eps^-4 im',
    'int87 eps^-4 im',
    'int88 eps^-4 im',
    'int89 eps^-4 im',
    'int90 eps^-4 im',
    'int91 eps^-4 im',
    'int92 eps^-4 im',
    'int93 eps^-4 im',
    'int94 eps^-4 im',
    'int95 eps^-4 im',
    'int96 eps^-4 im',
    'int97 eps^-4 im',
)


def q1b_counts():
    """the expectations of the second pair COUNTED from the fixture JSONs by the numerical-zero rule: over every integral and every common
    order, the components whose two members are both below T = max(radii, 10^-(30-2)) and not both exactly zero (numerical zeros),
    the exact-zero pairs, the compared components, and the orders on one side only."""
    a = json.load(open(fx("pair_q1b_goal30.json")))["result"]
    b = json.load(open(fx("pair_q1b_goal60.json")))["result"]
    assert len(a) == len(b) == Q1B_N_INTEGRALS
    t = mpf(10) ** (-(min(G30, G60) - 2))
    zeros, exact, compared, one_side = [], [], [], []
    for i, (ra, rb) in enumerate(zip(a, b)):
        ca, cb = cmp_mod.coeffs_of(ra), cmp_mod.coeffs_of(rb)
        one_side += [f"int{i} eps^{o}" for o in sorted(set(ca) ^ set(cb))]
        for o in sorted(set(ca) & set(cb)):
            for part, (x, y) in zip(("re", "im"), zip(ca[o], cb[o])):
                (xm, xr), (ym, yr) = x, y
                name = f"int{i} eps^{o} {part}"
                if xm == 0 and ym == 0:
                    exact.append(name)
                    compared.append(name)
                elif abs(xm) < max(xr, yr, t) and abs(ym) < max(xr, yr, t):
                    zeros.append(name)
                else:
                    compared.append(name)
    return zeros, exact, compared, one_side


def test_q1b_pair_im_numerical_zeros_read_pass():
    rc, out = run([fx("pair_q1b_goal30.json"), fx("pair_q1b_goal60.json"), "--goal", f"{G30},{G60}"])
    assert rc == 0, out.splitlines()[-1]                             # the tool's verdict first (the previous tool reads OVERALL FAIL here)
    zeros, exact, compared, one_side = q1b_counts()
    assert (len(zeros), len(compared), len(one_side)) == (Q1B_IM_ZEROS, Q1B_COMPARED, 0)   # the fixture's own figures
    assert tuple(zeros) == Q1B_IM_ZERO_NAMES
    assert all(n.endswith(" im") for n in zeros) and exact           # every numerical zero is an im part; exact-zero im parts exist beside
    t_str = cmp_mod.fmt(mpf(10) ** (-(min(G30, G60) - 2)))
    assert out.splitlines()[0] == (f"goal: REF {G30}, TEST {G60} (from --goal); zero threshold T = max(radii, "
                                   f"10^-({min(G30, G60)}-2) = {t_str})")
    assert overall(out) == ("PASS", Q1B_COMPARED, Q1B_IM_ZEROS, 0)
    assert f"numerical zeros: {Q1B_IM_ZEROS} ({', '.join(Q1B_IM_ZERO_NAMES)})" in out
    assert "orders not compared: 0 ()" in out
    for name in zeros:
        line = [l for l in out.splitlines() if l.startswith(name + ":")]
        assert len(line) == 1 and "NUMERICAL ZERO (|ref| " in line[0] and f"below T {t_str}) -> not compared" in line[0], line
    for name in exact:
        assert f"{name}: both exactly 0 -> PASS" in out
    assert "FAIL" not in out and "MISSING" not in out
    rc30, out30 = run([fx("pair_q1b_goal30.json"), fx("pair_q1b_goal60.json"), "--goal", f"{G30},{G60}", "--min-digits", "30"])
    assert rc30 == 0 and overall(out30) == ("PASS", Q1B_COMPARED, Q1B_IM_ZEROS, 0)
    # under --radii-only the im zeros are compared on the radii alone and read 0.0 digits -> FAIL
    rc_ng, out_ng = run([fx("pair_q1b_goal30.json"), fx("pair_q1b_goal60.json"), "--radii-only"])
    assert rc_ng == 1 and overall(out_ng) == ("FAIL", Q1B_COMPARED + Q1B_IM_ZEROS, 0, 0)


def test_q1b_pair_on_the_previous_tool_fails_on_the_im_components():
    prev_cli = os.environ.get("COMPARE_AMFLOW_JSON_PREVIOUS")
    if not prev_cli:
        pytest.skip("COMPARE_AMFLOW_JSON_PREVIOUS unset: the previous tool's file is not named")
    assert os.path.isfile(prev_cli), prev_cli
    rc, out = run([fx("pair_q1b_goal30.json"), fx("pair_q1b_goal60.json"), "--min-digits", "30"], cli=prev_cli)
    assert rc == 1
    assert out.rstrip().endswith("OVERALL: FAIL")
    fails = re.findall(r"^(int\d+ eps\^-?\d+ (re|im)): .* -> FAIL$", out, re.M)
    assert len([n for n, p in fails if p == "re"]) == Q1B_PREVIOUS_RE_FAIL
    im_fails = [n for n, p in fails if p == "im"]
    assert len(im_fails) == Q1B_IM_ZEROS and tuple(im_fails) == Q1B_IM_ZERO_NAMES   # the previous FAIL lines are the numerical zeros
    for name in im_fails:
        assert f"{name}: 0.0 digits (ref=" in out


def test_boolean_goal_key_is_not_a_goal(tmp_path):
    """a JSON boolean under a goal key (goal: true) is NOT a goal: with no other
    source the run refuses by name rc 2 and compares nothing (it must never read the boolean as goal 1 and call every
    order a numerical zero); with --goal it compares as usual."""
    a = make_out(str(tmp_path / "a.json"), [PAIR_OK])
    b = make_out(str(tmp_path / "b.json"), [PAIR_OK])
    for p in (a, b):
        doc = json.load(open(p))
        doc["goal"] = True
        json.dump(doc, open(p, "w"), indent=1)
    rc, out = run([a, b])
    assert rc == 2 and "no goal" in out and "OVERALL" not in out and "REF 1" not in out and "Traceback" not in out, out
    rc2, out2 = run([a, b, "--goal", "60"])
    assert rc2 == 0 and overall(out2)[0] == "PASS", out2

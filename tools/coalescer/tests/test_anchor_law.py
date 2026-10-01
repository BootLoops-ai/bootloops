"""Anchor-law battery for coalescer.py: the s_b floor computed from the operator's own singularities.

Conventions (pf_engine header): z = 1/p^2 is the MUM coordinate, s = -p^2, so z = -1/s; the Route A seed is
the BFKNS multinomial series sum_n a_n z_b^n evaluated at z_b = -1/s_b (make_period_state).  The seed
converges for |z_b| < R, i.e. s_b > 1/R.  Two readings of the README's "s_b > 1/min|z_sing|":

  seed   R = 1/thr, thr = (sum_i sqrt(m_i^2))^2 -- the leading threshold z = -1/thr is the nearest true
         singularity of the seed series (exact: thr^n/(n+1)^(L-1) <= |a_n| <= thr^n by the multinomial
         identity), floor s_b > thr.  The CLI default.
  roots  R = min|z| over ALL roots of the leading polynomial P_order(z) (the same list Route B's auto-anchor
         ceil(1.6/min|z_sing|) uses); roots inside |z| < 1/thr are apparent for the seed.  Conservative.

Recorded values pinned below, each from an object (never retyped from memory):
  K3  (1,1,1,9):      thr 36;  leading roots 5;  min|z_sing| = 1/36 (the leading threshold itself; the
                      comment in routeB_monoproj.jl: "for K3 |z|=1/36 -> s_b > 36"); both floors 36.
  CY3 (1,1,1,1,16):   thr 64;  leading roots 9;  min|z_sing| = 0.005598521790772446, the real root of the
                      quintic factor 331776 z^5 - 161280 z^4 + 125840 z^3 + 5484 z^2 + 144 z - 1 (Route B's log
                      of the CY3 gate: "min |z_sing| ~ 0.005599 (=> s>178.6); anchor s_b=286"); roots floor
                      178.61857779105415; seed floor 64.  The reference Route A CY3 gate ran at s_b = 100 with the
                      seed agreeing to 80.0 of 80 digits -- the root at 0.0056 is apparent for the seed.
  CY4 (1,1,1,1,1,25): thr 100; leading roots 14; min|z_sing| = 0.0013248607864452342 (the complex pair
                      0.0002678142236327703 +/- 0.0012975097861211919 i of the degree-9 factor), roots floor
                      754.7962851879132 -- the figures of the CY4 probe's ZSING_CHECK record; Route B auto 1208;
                      the probe ran --sb 1210 --sb2 2420; the CLI defaults 2.2*thr = 220.00000000000003 /
                      4.2*thr = 420.0; seed floor 100.
The in-process legs build each operator once (K3 ~1 s, CY3 ~5 s, CY4 ~12 s) and exercise anchor_law()
directly; the CLI legs are K3 (seconds) plus the one CY4 refusal, and the generic --op door on a K3
operator + seed JSON freshly written by dump_pf_data.dump (check-only line, refusal rc 3, usage rc 2,
loader errors rc 2).  The K3 --test masked-identical leg (both doors + the Gauss 2F1 example) takes ~35 s
and is skipped by name unless COALESCER_SLOW=1.
"""
import importlib.util
import os
import re
import subprocess
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
CLI = os.path.join(PKG, "coalescer.py")
GOLDEN = os.path.join(HERE, "k3_selftest_masked.txt")

RECORD = {
    "K3": dict(msq=(1, 1, 1, 9), thr=36.0, n_roots=5, min_abs="0.027777777777777776",
               floor_roots="36.0", floor_seed="36.0", inside=0,
               sb_default="79.2", sb2_default="151.20000000000002", sb_auto_B=58),
    "CY3": dict(msq=(1, 1, 1, 1, 16), thr=64.0, n_roots=9, min_abs="0.005598521790772446",
                floor_roots="178.61857779105415", floor_seed="64.0", inside=1,
                sb_default="140.8", sb2_default="268.8", sb_auto_B=286),
    "CY4": dict(msq=(1, 1, 1, 1, 1, 25), thr=100.0, n_roots=14, min_abs="0.0013248607864452342",
                floor_roots="754.7962851879132", floor_seed="100.0", inside=2,
                sb_default="220.00000000000003", sb2_default="420.0", sb_auto_B=1208),
}


def _load_cli():
    if PKG not in sys.path:
        sys.path.insert(0, PKG)
    spec = importlib.util.spec_from_file_location("coalescer_cli", CLI)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def cli():
    return _load_cli()


@pytest.fixture(scope="module")
def ops(cli):
    """The three operators, built exactly as routeA_monoproj.run builds them (find_minimal_op)."""
    import mpmath as mp
    out = {}
    for name, rec in RECORD.items():
        order, degz, Pj, resid, _ = cli.find_minimal_op(list(rec["msq"]), ord_max=10, degz_max=18)
        thr = sum(mp.sqrt(m) for m in rec["msq"]) ** 2          # as main() computes it (mp.dps 15)
        out[name] = dict(order=order, Pj=Pj, thr=thr, resid=resid,
                         sb=float(2.2 * thr), sb2=float(4.2 * thr))
    return out


def _law(cli, ops, name, sb=None, sb2=None, law="seed"):
    o = ops[name]
    sb = o["sb"] if sb is None else float(sb)
    sb2 = o["sb2"] if sb2 is None else float(sb2)
    return cli.anchor_law(RECORD[name]["msq"], o["Pj"], o["order"], o["thr"], sb, sb2, law)


def test_floors_from_the_operators_match_the_record(cli, ops):
    for name, rec in RECORD.items():
        r = _law(cli, ops, name)
        assert ops[name]["resid"] == 0
        assert r["n_roots"] == rec["n_roots"]
        assert repr(r["min_abs"]) == rec["min_abs"]
        assert repr(r["floor_roots"]) == rec["floor_roots"]
        assert repr(r["floor_seed"]) == rec["floor_seed"]
        assert r["thr"] == rec["thr"] and r["thr_is_root"] is True
        assert r["roots_inside_seed_disk"] == rec["inside"]
        assert r["sb_auto_routeB"] == rec["sb_auto_B"]
        assert repr(ops[name]["sb"]) == rec["sb_default"] and repr(ops[name]["sb2"]) == rec["sb2_default"]


def test_written_defaults_under_each_law(cli, ops):
    # seed law: every written default is above thr (2.2 > 1 by construction)
    for name in RECORD:
        r = _law(cli, ops, name, law="seed")
        assert r["law_in_force"] == "seed" and r["sb_ok"] and r["sb2_ok"]
    # roots law: K3's floors coincide (the leading threshold is its nearest root); CY3's sb and both CY4
    # anchors sit below the all-roots floor
    r = _law(cli, ops, "K3", law="roots"); assert r["sb_ok"] and r["sb2_ok"]
    r = _law(cli, ops, "CY3", law="roots"); assert (r["sb_ok"], r["sb2_ok"]) == (False, True)
    r = _law(cli, ops, "CY4", law="roots"); assert (r["sb_ok"], r["sb2_ok"]) == (False, False)


def test_overrides_against_the_cy4_roots_floor(cli, ops):
    ok = lambda r: (r["sb_ok"], r["sb2_ok"])
    assert ok(_law(cli, ops, "CY4", 1210, 2420, "roots")) == (True, True)     # the recorded probe's anchors
    assert ok(_law(cli, ops, "CY4", 754, 1508, "roots")) == (False, True)     # just below the floor
    assert ok(_law(cli, ops, "CY4", 755, 1510, "roots")) == (True, True)      # just above
    assert ok(_law(cli, ops, "CY4", 1210, 420, "roots")) == (True, False)     # sb2 below its floor
    assert ok(_law(cli, ops, "CY4", 90, 180, "seed")) == (False, True)        # below thr: the seed diverges
    assert ok(_law(cli, ops, "CY4", 101, 202, "seed")) == (True, True)


def test_refusal_message_names_everything(cli, ops):
    r = _law(cli, ops, "CY4", law="roots")
    src = dict(sb="default 2.2*thr", sb2="default 4.2*thr", lift="default thr/4", lift_val=25.0)
    msg = cli.anchor_law_refusal(r, src)
    assert msg.startswith("[anchor-law] REFUSED (rc 3): masses=(1,1,1,1,1,25): law=roots floor s_b > 1/min|z_sing| = 754.7962851879132")
    assert "min|z_sing| = 0.0013248607864452342 over the 14 roots of the order-8 operator's leading polynomial" in msg
    assert "sb = 220.00000000000003 (default 2.2*thr) and sb2 = 420.0 (default 4.2*thr) are not above the floor" in msg
    assert "Nothing transported. Override: --sb S --sb2 S2 with S, S2 > 754.7962851879132, or --anchor-law seed (floor thr = 100.0)." in msg
    r = _law(cli, ops, "K3", 30, 60, "seed")
    msg = cli.anchor_law_refusal(r, dict(sb="override", sb2="override", lift="default thr/4", lift_val=9.0))
    assert "law=seed floor s_b > thr = 36.0" in msg and "sb = 30.0 (override) is not above the floor" in msg
    line = cli.anchor_law_line(_law(cli, ops, "CY4"), src, "A")
    assert "floor_roots=754.7962851879132 floor_seed=100.0" in line and "verdict=OK" in line
    assert cli.anchor_law_line(_law(cli, ops, "CY4"), src, "B").rstrip().find("sb_auto(routeB)=1208") > 0


def _run(*args):
    env = dict(PATH=os.environ.get("PATH", ""), PYTHONDONTWRITEBYTECODE="1")
    p = subprocess.run([sys.executable, CLI, *args], capture_output=True, text=True, env=env, cwd=HERE)
    return p.returncode, p.stdout, p.stderr


def test_cli_k3_check_only_is_rc0():
    rc, out, err = _run("--masses", "1,1,1,9", "--check-only")
    assert rc == 0, (out, err)
    assert "[anchor-law] masses=(1,1,1,9) thr=36.0 order=4 leading_roots=5 min|z_sing|=0.027777777777777776 floor_roots=36.0 floor_seed=36.0" in out
    assert "verdict=OK" in out and "check only: nothing transported" in out


def test_cli_k3_sb_below_thr_is_refused_rc3():
    rc, out, err = _run("--masses", "1,1,1,9", "--sb", "30", "--check-only")
    assert rc == 3, (out, err)
    assert "[anchor-law] REFUSED (rc 3): masses=(1,1,1,9): law=seed floor s_b > thr = 36.0" in out
    assert "sb = 30.0 (override) is not above the floor" in out and "sb2=60.0 (derived 2*sb)" in out


def test_cli_cy4_roots_law_refuses_the_written_defaults_rc3():
    rc, out, err = _run("--masses", "1,1,1,1,1,25", "--anchor-law", "roots")
    assert rc == 3, (out, err)
    assert "floor s_b > 1/min|z_sing| = 754.7962851879132 (min|z_sing| = 0.0013248607864452342" in out
    assert "sb = 220.00000000000003 (default 2.2*thr) and sb2 = 420.0 (default 4.2*thr) are not above the floor" in out
    assert "Route A monodromy-projection" not in out          # nothing transported


def test_cli_usage_errors_are_rc2():
    rc, out, err = _run("--dps", "30")
    assert rc == 2 and "--masses required" in err
    rc, out, err = _run("--masses", "1,1,1,9", "--route", "B", "--sb", "58.5", "--check-only")
    assert rc == 2 and "integer --sb" in err


def _mask(text):
    text = re.sub(r"\[\d+\.\d+s( total)?\]", lambda m: "[Ts%s]" % (m.group(1) or ""), text)
    text = re.sub(r"(built|loaded) in \d+\.\d+s\]", lambda m: "%s in Ts]" % m.group(1), text)
    return re.sub(r"wall = \d+\.\d+s", "wall = Ts", text)


# ---------------------------------------------------------------- the generic --op door (fast legs)
@pytest.fixture(scope="module")
def k3_op_json(tmp_path_factory):
    """K3 operator + seed re-exported through dump_pf_data.dump into the --op JSON layout."""
    if PKG not in sys.path:
        sys.path.insert(0, PKG)
    from dump_pf_data import dump
    d = tmp_path_factory.mktemp("k3op")
    _, js, summ = dump((1, 1, 1, 9), nthr=20, nser=120, out=str(d / "pf_k3_op.jl"), verbose=False)
    assert summ["order"] == 4 and summ["degz"] == 5 and summ["residual"] == "0"
    assert summ["frac_exps"] == ["3/2"] and summ["seed_radius"] == "1/36"
    assert len(summ["bfkns_a"]) == 121 and summ["bfkns_a"][:3] == [1, -12, 204]
    return js


def test_op_loader_round_trips_the_operator(cli, ops, k3_op_json):
    import sympy as sp
    OP = cli.load_operator_json(k3_op_json)
    assert OP["order"] == 4 and OP["degz"] == 5 and OP["alpha"] is None
    assert OP["seed_radius"] == pytest.approx(1 / 36) and len(OP["seed"]) == 121
    for j in range(5):                      # same polynomials find_minimal_op built in-process
        assert sp.expand(OP["Pj"][j] - ops["K3"]["Pj"][j]) == 0
    r = cli.anchor_law(None, OP["Pj"], 4, None, 79.2, 151.2, "roots", what="op=x", seed_radius=OP["seed_radius"])
    assert repr(r["floor_roots"]) == RECORD["K3"]["floor_roots"] and r["law_in_force"] == "roots"
    assert r["sb_ok"] and r["sb2_ok"] and r["thr_is_root"] and r["roots_inside_seed_disk"] == 0
    r = cli.anchor_law(None, OP["Pj"], 4, None, 79.2, 151.2, "seed", what="op=x", seed_radius=None)
    assert r["law_in_force"] == "roots" and "declares no seed_radius" in r["note"]


def test_cli_op_check_only_prints_the_line_and_transports_nothing(k3_op_json):
    rc, out, err = _run("--op", k3_op_json, "--check-only")
    assert rc == 0, (out, err)
    assert ("[anchor-law] op=pf_k3_op.json (banana msq=(1, 1, 1, 9)) seed_radius=0.027777777777777776 order=4 "
            "leading_roots=5 min|z_sing|=0.027777777777777776 floor_roots=36.0 floor_seed=36.0") in out
    assert "law=roots floor=36.0 sb=79.2 (default 2.2*floor) sb2=151.20000000000002 (default 4.2*floor)" in out
    assert "lift=9.0 (default floor/4) verdict=OK" in out and "seed terms=120" in out
    assert "check only: nothing transported" in out and "Route A monodromy-projection" not in out


def test_cli_op_sb_below_floor_is_refused_rc3(k3_op_json):
    rc, out, err = _run("--op", k3_op_json, "--sb", "30", "--check-only")
    assert rc == 3, (out, err)
    assert "[anchor-law] REFUSED (rc 3): op=pf_k3_op.json (banana msq=(1, 1, 1, 9)): law=roots floor s_b > 1/min|z_sing| = 36.0" in out
    assert "sb = 30.0 (override) is not above the floor" in out
    assert "or --anchor-law seed (floor 1/seed_radius = 36.0)" in out


def test_cli_op_usage_errors_are_rc2(k3_op_json, tmp_path):
    rc, out, err = _run("--op", k3_op_json, "--masses", "1,1,1,9", "--check-only")
    assert rc == 2 and "not allowed with argument" in err
    rc, out, err = _run("--op", k3_op_json, "--route", "B", "--check-only")
    assert rc == 2 and "--op runs Route A" in err
    bad = tmp_path / "noseed.json"
    bad.write_text('{"theta_coeffs": [[0, 1], [1]], "label": "toy"}', encoding="utf-8")
    rc, out, err = _run("--op", str(bad), "--check-only")
    assert rc == 2 and "no seed" in err
    bad2 = tmp_path / "floatseed.json"
    bad2.write_text('{"Pj": [[0, 1], [1]], "seed": [1, 0.5]}', encoding="utf-8")
    rc, out, err = _run("--op", str(bad2), "--check-only")
    assert rc == 2 and "inexact float seed coefficient" in err


@pytest.mark.skipif(not os.environ.get("COALESCER_SLOW"),
                    reason="K3 --test probe selftest takes ~55-70 s: set COALESCER_SLOW=1 to run it")
def test_k3_selftest_output_masked_identical_to_golden():
    rc, out, err = _run("--test")
    assert rc == 0, (out, err)
    assert _mask(out) == open(GOLDEN, encoding="utf-8").read()
    assert "[selftest] match = 24.73 d" in out

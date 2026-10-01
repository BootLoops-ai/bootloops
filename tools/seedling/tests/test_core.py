"""Unit tests for the core pipeline (synthetic fixtures only)."""
import os, sys, tempfile
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from seedling import support, pinner, bulkmodel

FIXTURE = """\
fam[1, 0, 2, -1, 0]

fam[1, 1, 1, 0, -2]
fam[0, 0, 1, 0, 0]
"""


def test_support():
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as fh:
        fh.write(FIXTURE); path = fh.name
    tg = support.parse_targets(path)
    assert len(tg) == 3
    a = tg[0]  # fam[1,0,2,-1,0]: t=2, r=3, s=1, d=1, sector=0b101=5
    assert (a.t, a.r, a.s, a.d, a.sector) == (2, 3, 1, 1, 5)
    g = support.global_support(tg)
    assert (g["r"], g["s"], g["d"], g["t"], g["n"]) == (3, 2, 1, 3, 3)
    tab = support.support_table(tg)
    assert len(tab) == 3 and tab[5]["n"] == 1
    os.unlink(path)


def test_symmetry_closure():
    from seedling import rankcheck
    # sectorSymmetries fixture (3 props): 5 -> 3 (prop0->pos0, prop2->pos1)
    # and the inverse 3 -> 5, so the closure walks a 2-cycle to a fixpoint
    lines = ("5 0 0 0 0 0 0 1 3 2 0 {r} {a} 0 -1 1 1 0\n"
             "3 0 0 0 0 0 0 1 5 2 0 {r} {a} 0 2 -1 1 0\n")
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as fh:
        fh.write(lines); mpath = fh.name
    maps = rankcheck.parse_sector_maps(mpath, 3)
    os.unlink(mpath)
    assert maps[0] == {"src": 5, "tgt": 3, "perm": [0, -1, 1]}
    # a sector-5 target's demand lands in sector 3 with the same support
    tg = [support.Target("fam", (1, 0, 2))]         # sector 5: r=3, s=0, d=1
    closed, rep = support.symmetry_closure(tg, maps)
    assert rep["n_added"] == 1 and not rep["undecidable"]
    assert not rep["sector_mismatch"]
    tab = support.support_table(closed)
    assert set(tab) == {5, 3}
    assert (tab[3]["r"], tab[3]["s"], tab[3]["d"]) == (3, 0, 1)
    # fixpoint: closing the closed list adds nothing (2-cycle terminates)
    closed2, rep2 = support.symmetry_closure(closed, maps)
    assert rep2["n_added"] == 0 and len(closed2) == len(closed)
    # ISP power at an unmapped position: image undecidable -> REPORTED,
    # never fabricated and never silently dropped
    isp = [support.Target("fam", (1, -1, 2))]
    closed3, rep3 = support.symmetry_closure(isp, maps)
    assert len(closed3) == 1
    assert rep3["undecidable"] == [{"indices": [1, -1, 2],
                                    "src": 5, "tgt": 3}]


def test_pinner_shape():
    txt = pinner.render_jobs_yaml("fam", [31], 7, 1, 2)
    assert "{topologies: [fam], sectors: [31], r: 7, s: 1, d: 2}" in txt
    assert "truncate_sp" not in txt
    txt2 = pinner.render_jobs_yaml("fam", [31], 7, 1, 2, truncate_sp_l=3)
    assert "truncate_sp:" in txt2 and "l: 3}" in txt2
    # truncate block must precede select_integrals (kira job order)
    assert txt2.index("truncate_sp") < txt2.index("select_integrals")


def test_bulkmodel():
    # single sector t=2 in a 5-sp family, box r=3,s=1,d=1:
    # dhat=min(1,1)=1 -> C(3,2)=3 dots; m=3 -> C(4,3)=4 ISP: 12 seeds x nops
    assert bulkmodel.seeds_per_sector(2, 5, 3, 1, 1) == 12
    assert bulkmodel.bulk([2], 5, 6, 3, 1, 1) == 72
    # r below t -> zero seeds
    assert bulkmodel.seeds_per_sector(4, 5, 3, 0, 0) == 0
    pf = bulkmodel.preflight([2], 5, 6, 3, 1, 1)
    assert pf["bulk_eqns"] == 72 and pf["generate_peak_GB"] > bulkmodel.A_GEN




def test_chain_s_floor():
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as fh:
        fh.write(FIXTURE); path = fh.name
    tg = support.parse_targets(path); os.unlink(path)
    # fixture has ISP-bearing targets -> floor >= their max s = 2
    assert pinner.chain_s_floor(tg, [7]) == 2
    # all-top-sector, s=0 targets -> floor 0 allowed
    top = [support.Target("fam", (1, 1, 1, 0, 0))]
    assert pinner.chain_s_floor(top, [7]) == 0
    # s=0 target BELOW top sector -> floor 1 (chain-floor lesson)
    sub = [support.Target("fam", (1, 1, 0, 0, 0))]
    assert pinner.chain_s_floor(sub, [7]) == 1
    # pinned_cell applies the floor: support s=0 but sub-sector target
    g = {"r": 2, "s": 0, "d": 0}
    assert pinner.pinned_cell(g, sub, [7]) == (2, 1, 0)
    # explicit unsafe override honored
    assert pinner.pinned_cell(g, sub, [7], s_floor=0) == (2, 0, 0)


def test_runner_gate():
    from seedling import runner
    b = {"r1": 5, "r2": 7, "r3": 11}
    ok = dict(b)
    bad = {"r1": 5, "r2": 99, "r3": 12}
    # clean pass
    assert runner.two_slice_gate(ok, b)["verdict"] == runner.GATE_PASS
    # fail without second slice -> demand it
    assert runner.two_slice_gate(bad, b)["verdict"] == "NEED_SECOND_SLICE"
    # both primes fail -> pin at fault (chain-loss pattern)
    v = runner.two_slice_gate(bad, b, bad, b)
    assert v["verdict"] == runner.GATE_FAIL_PIN and v["fails1"] == ["r2", "r3"]
    # second prime passes -> slice artifact
    assert runner.two_slice_gate(bad, b, ok, b)["verdict"] == runner.GATE_RERUN_SLICE


def test_runner_dryrun():
    from seedling import runner
    with tempfile.TemporaryDirectory() as d:
        ledger = os.path.join(d, "ledger.jsonl")
        rec, script = runner.run_phase(
            "stage", d, "kira jobs.yaml", 
            {"name": "testcap", "mem_max": "12G", "cpu_max": "400000 100000"},
            ledger, dry_run=True)
        assert rec["exit_status"] is None and rec["dry_run"]
        assert "memory.max" in script and "prlimit" not in script
        # cgroup.procs write-discipline guard tokens ride every script
        assert "${CG:?}" in script and 'case "$CG"' in script \
            and "grep -qx" in script
        assert os.path.exists(ledger)




def test_de_census_gate():
    from seedling import runner
    m = [(1, 1, 1, 0), (1, 1, 1, -1)]
    g = runner.de_census_gate(m, [0, 1, 2], gated_rows=[(2, 1, 1, 0), (1, 2, 1, 0)])
    assert g["census_n"] == 6 and g["gated_n"] == 2 and g["master_self_n"] == 0
    assert len(g["uncovered"]) == 4 and (2, 1, 1, -1) in g["uncovered"]
    spec = runner.minikira_spec("fam", [15], g["uncovered"])
    assert spec["cell"] == (5, 2, 2)          # support (4,1,1) + one unit everywhere
    assert "fam[2, 1, 1, -1]" in spec["target"]
    assert "r: 5, s: 2, d: 2" in spec["jobs_yaml"]


def test_escalate_ladders():
    from seedling import escalate as E
    st = E.EscalationState()
    # staging incompleteness: r -> d -> s order, then r again? no: each dim to cap
    a, m = E.next_margins(st, E.LEFTOVER_TARGETS)
    assert (a, m) == ("ESCALATE", (1, 0, 0))
    # chain-loss gate failure: s FIRST (chain-floor lesson)
    a, m = E.next_margins(st, E.GATE_FAIL_CHAIN)
    assert (a, m) == ("ESCALATE", (1, 1, 0))
    # rank problems are NOT margin problems
    a, m = E.next_margins(st, E.RANK_DEFICIT)
    assert a == E.REFER and m == (1, 1, 0)
    # de-census -> minikira, no margin change
    a, m = E.next_margins(st, E.DE_CENSUS_UNCOVERED)
    assert a == E.RUN_MINIKIRA and m == (1, 1, 0)
    # bounded: exhaust everything -> GIVE_UP
    st2 = E.EscalationState(margins=(3, 2, 3), max_margins=(3, 2, 3))
    a, m = E.next_margins(st2, E.REGENERATE_LOOP)
    assert a == E.GIVE_UP
    assert len(st.history) == 4


def test_classify_gate_failure():
    from seedling import escalate as E
    assert E.classify_gate_failure([(1, 1, 1, -1)]) == E.GATE_FAIL_CHAIN
    assert E.classify_gate_failure([(2, 1, 1, 0)]) == E.GATE_FAIL_OTHER


def test_certify_zero_normalize():
    from seedling import certify
    # the numeric-d artifact: reference carries a term whose coefficient
    # vanishes at this d; fresh slice dropped it. Must NOT flag as DIFF.
    ref = {"T[1]": {"m1": 3, "m2": 0}, "T[2]": {"m1": 1}}
    fresh = {"T[1]": {"m1": 3}, "T[2]": {"m1": 1}}
    rep = certify.verify_slice(fresh, ref)
    assert rep["fails"] == [] and rep["n_common"] == 2
    # a REAL diff still flags
    bad = {"T[1]": {"m1": 4}, "T[2]": {"m1": 1}}
    assert certify.verify_slice(bad, ref)["fails"] == ["T[1]"]


def test_certify_two_slice():
    from seedling import certify, runner
    ref = {"T[1]": {"m1": 3, "m2": 0}}
    ok = {"T[1]": {"m1": 3}}
    bad = {"T[1]": {"m1": 5}}
    assert certify.two_slice_certify(ok, ref)["verdict"] == runner.GATE_PASS
    v = certify.two_slice_certify(bad, ref, bad, ref)
    assert v["verdict"] == runner.GATE_FAIL_PIN
    assert certify.two_slice_certify(bad, ref, ok, ref)["verdict"] == runner.GATE_RERUN_SLICE


def test_certify_syzygy_degrade():
    from seedling import certify
    assert certify.syzygy_certify({}, oracle=None)["status"] == "SKIPPED"
    oracle = lambda tw: {"T[1]": {"m1": 3}}
    rep = certify.syzygy_certify({"T[1]": {"m1": 3, "m2": 0}}, oracle, towers=["tw0"])
    assert rep["status"] == "PASS"


def test_kira_target_m_splitter():
    from seedling import certify
    import tempfile, os
    text = "fam[1,1,-1] -> \n  fam[1,1,0]*(3/2) + fam[1,0,0]*(d-2)\n\nfam[2,1,0] -> fam[1,1,0]*(1),\n"
    with tempfile.NamedTemporaryFile("w", suffix=".m", delete=False) as fh:
        fh.write(text); p = fh.name
    rows = certify.split_kira_target_m(p); os.unlink(p)
    assert set(rows) == {"fam[1,1,-1]", "fam[2,1,0]"}
    assert rows["fam[2,1,0]"] == "fam[1,1,0]*(1)"


def test_schedule_v1():
    from seedling import pinner
    tab = {7: {"r": 4, "s": 1, "d": 1, "t": 3, "n": 2},
           3: {"r": 2, "s": 0, "d": 0, "t": 2, "n": 1}}
    census = {7: 3, 3: 2, 5: 2, 6: 2}
    g = pinner.build_schedule(tab, census, closure_cell=lambda t: (t + 2, 1, 2))
    assert g[(4, 1, 1)] == [7]        # target sector: own support
    assert g[(2, 1, 0)] == [3]        # target sector: s floored to 1 (chain floor)
    assert g[(4, 1, 2)] == [5, 6]     # closure sectors grouped by cell
    y = pinner.render_schedule_jobs_yaml("fam", g)
    assert y.count("- {topologies") == 3
    assert "sectors: [5, 6], r: 4, s: 1, d: 2" in y
    assert "select_mandatory_list" in y and "run_initiate: true" in y


def test_rankcheck():
    from seedling import rankcheck
    import tempfile, os
    log = ("Number of master integrals: 2\n"
           "fam[2,1,0]  #  3\n"
           "fam[1,2,0]  #  3\n"
           "Number of selected equations to reduce: 5 equations\n")
    with tempfile.NamedTemporaryFile("w", suffix=".log", delete=False) as fh:
        fh.write(log); p = fh.name
    masters = rankcheck.parse_kira_masters(p); os.unlink(p)
    assert masters == [((2, 1, 0), 3), ((1, 2, 0), 3)]
    # self-map swapping props 0<->1 in sector 3 relates the two dotted masters
    maps = [{"src": 3, "tgt": 3, "perm": [1, 0, -1]}]
    rep = rankcheck.rank_report(masters, maps=maps)
    assert rep["verdict"] == "LONG" and len(rep["redundant_pairs"]) == 1
    # image undecidable (ISP power at unmapped position) -> no fabricated pair
    assert rankcheck.apply_perm((1, 1, -1), [1, 0, -1]) is None
    # count directions
    assert rankcheck.rank_report(masters, expected_count=3)["verdict"] == "SHORT"
    assert rankcheck.rank_report(masters, expected_count=2)["verdict"] == "OK"
    # nothing independent supplied -> SUSPECT, never silent OK
    assert rankcheck.rank_report(masters)["verdict"] == "SUSPECT"


def test_run_const_dry_selfcontained():
    """`--margins const` runs the full dry pipeline with NO model/labels:
    the const arm must never touch the frozen model, and --sector-maps must
    fold symmetry images into the schedule (ledgered, target-grade cells)."""
    import json
    from seedling import cli, margins

    def _boom(*a, **k):
        raise AssertionError("const arm must not load the model/labels")

    orig_fit, orig_ml = margins.fit_frozen, margins._ml
    margins.fit_frozen = _boom
    margins._ml = _boom
    try:
        with tempfile.TemporaryDirectory() as td:
            cfg = os.path.join(td, "config")
            os.makedirs(cfg)
            with open(os.path.join(cfg, "integralfamilies.yaml"), "w") as fh:
                fh.write("integralfamilies:\n - name: \"toy\"\n"
                         "   loop_momenta: [l]\n   top_level_sectors: [3]\n"
                         "   propagators:\n    - [\"l^2\", 0]\n"
                         "    - [\"(l+p1)^2\", 0]\n")
            tpath = os.path.join(td, "targets")
            with open(tpath, "w") as fh:
                fh.write("toy[1,0]\n")
            npath = os.path.join(td, "nts")
            with open(npath, "w") as fh:
                fh.write("3 2\n2 1\n1 1\n")
            mpath = os.path.join(td, "maps")
            with open(mpath, "w") as fh:      # sector 1 -> 2 (prop0->pos1)
                fh.write("1 0 0 0 0 0 0 1 2 1 0 {r} {a} 1 -1 1 0\n")
            out = os.path.join(td, "run")
            rc = cli.main(["run", "--config", cfg, "--targets", tpath,
                           "--out", out, "--nts", npath, "--margins", "const",
                           "--sector-maps", mpath])
            assert rc == 0
            led = [json.loads(l)
                   for l in open(os.path.join(out, "ledger.jsonl"))]
            clo = [r for r in led if r["phase"] == "support_closure"]
            assert clo and clo[0]["n_added"] == 1 \
                and clo[0]["added"] == [["toy", [0, 1], 2]]
            mrow = [r for r in led if r["phase"] == "margins"][0]
            assert mrow["mode"] == "const"
            assert all(v["margin"] == 6
                       for v in mrow["per_class"].values())
            # image sector 2 rides the target-grade cell beside its source
            # (closure sector 3 keeps the chain cell t+6 = (8, 1, 6))
            sched = [r for r in led if r["phase"] == "schedule"][0]
            assert sched["groups"]["(1, 1, 0)"] == [1, 2]
            assert sched["groups"]["(8, 1, 6)"] == [3]
    finally:
        margins.fit_frozen = orig_fit
        margins._ml = orig_ml


def test_cli_wired():
    import json
    from seedling.cli import main
    with tempfile.TemporaryDirectory() as d:
        ref = os.path.join(d, "ref.json"); ok = os.path.join(d, "ok.json")
        json.dump({"T1": {"m1": 3, "m2": 0}}, open(ref, "w"))
        json.dump({"T1": {"m1": 3}}, open(ok, "w"))
        assert main(["certify", "--fresh", ok, "--ref", ref]) == 0
        bad = os.path.join(d, "bad.json")
        json.dump({"T1": {"m1": 9}}, open(bad, "w"))
        assert main(["certify", "--fresh", bad, "--ref", ref,
                     "--fresh2", bad, "--ref2", ref]) == 1
        # ledger on empty rundir
        assert main(["ledger", "--rundir", d]) == 1


if __name__ == "__main__":
    test_support(); test_symmetry_closure(); test_pinner_shape()
    test_bulkmodel()
    test_chain_s_floor(); test_runner_gate(); test_runner_dryrun()
    test_de_census_gate()
    test_escalate_ladders(); test_classify_gate_failure()
    test_certify_zero_normalize(); test_certify_two_slice()
    test_certify_syzygy_degrade(); test_kira_target_m_splitter()
    test_schedule_v1(); test_rankcheck()
    test_run_const_dry_selfcontained(); test_cli_wired()
    print("tests OK (18)")

#!/usr/bin/env python3
"""emitall test battery — planted truths, negative controls, incident replays.

Run:  python3 tests/test_all.py   (from tools/emitall/)  or  pytest -q

House rules: every behavior claim in README.md has
an executable receipt here; negative controls prove the checks can fail;
the shipped example deployment (specs/toy_stars.json) is exercised
end-to-end, clean and with planted corruptions. No bare asserts in package
code (tests use unittest assertions, which survive -O by design of this
suite being the receipt, not the guard).
"""
import json
import os
import subprocess
import sys
import tempfile
import textwrap
import time
import unittest
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
PKG_PARENT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, PKG_PARENT)

from emitall import (run_spec, compare_printed, compare_outward_band,  # noqa: E402
                     fmt_like, evaluate, resolve, resolve_pattern, Env,
                     ExprError, SpecError)
from emitall.lint import lint_source  # noqa: E402

# Shipped example deployment: specs/toy_stars.json + specs/toy_receipts/ —
# an INVENTED program whose numbers are synthetic, exercised end-to-end below.
SPECS = os.path.join(PKG_PARENT, "emitall", "specs")
RUN_PY = os.path.join(PKG_PARENT, "emitall", "run.py")
LINT_PY = os.path.join(PKG_PARENT, "emitall", "lint.py")
BATTERY_PY = os.path.join(PKG_PARENT, "emitall", "battery.py")


def _lint(src):
    return lint_source(textwrap.dedent(src), "<corpus>")


def _env(**receipts):
    return Env(receipts)


class TestPointers(unittest.TestCase):
    def test_nested_and_index(self):
        obj = {"a": {"b": [10, {"c": 7}]}, "x/y": 3, "t~": 4}
        self.assertEqual(resolve(obj, "/a/b/0"), 10)
        self.assertEqual(resolve(obj, "/a/b/1/c"), 7)
        self.assertEqual(resolve(obj, "a/b/1/c"), 7)      # leading / optional
        self.assertEqual(resolve(obj, "/x~1y"), 3)        # ~1 -> /
        self.assertEqual(resolve(obj, "/t~0"), 4)         # ~0 -> ~
        self.assertEqual(resolve(obj, ""), obj)
        with self.assertRaises(ExprError):
            resolve(obj, "/a/nope")
        with self.assertRaises(ExprError):
            resolve(obj, "/a/b/9")

    def test_pattern_star_and_alternation(self):
        obj = {"g": {"p1": {"v": 1}, "p2": {"v": 2}, "p3": {"v": 3}},
               "L": [{"s": 5}, {"s": 6}]}
        self.assertEqual(resolve_pattern(obj, "/g/*/v"), [1, 2, 3])
        self.assertEqual(resolve_pattern(obj, "/g/{p3|p1}/v"), [3, 1])  # order
        self.assertEqual(resolve_pattern(obj, "/L/*/s"), [5, 6])
        self.assertEqual(
            resolve_pattern(obj, "/g/{p1|p2}/{v|v}"), [1, 1, 2, 2])
        with self.assertRaises(ExprError):
            resolve_pattern(obj, "/g/p1/v/*")


class TestExpressions(unittest.TestCase):
    def setUp(self):
        self.env = _env(r={"n": {"a": 3, "b": 4}, "xs": [1, 2, 3, 4],
                           "s": "993/1067",
                           "rows": [{"k": "p", "v": 1}, {"k": "q", "v": 0},
                                    {"k": "p", "v": 5}],
                           "flag": True})

    def ev(self, e, **vars_):
        return evaluate(e, self.env, vars_ or None)

    def test_arithmetic_and_fields(self):
        self.assertEqual(self.ev("f('r','/n/a') + f('r','/n/b')"), 7)
        self.assertEqual(self.ev("100 * f('r','/n/a') / f('r','/n/b')"), 75.0)
        self.assertEqual(self.ev("len(f('r','/xs'))"), 4)
        self.assertEqual(self.ev("sum(vals('r','/xs/*'))"), 10)
        self.assertEqual(self.ev("max(vals('r','/xs/*') + [9])"), 9)
        self.assertEqual(self.ev("min(where(vals('r','/xs/*'), '>', 2))"), 3)

    def test_bool_comparison_arithmetic(self):
        # counting via summed comparisons — the "unique optima" pattern
        self.assertEqual(
            self.ev("(f('r','/n/a') == 3) + (f('r','/n/b') == 3)"), 1)
        self.assertEqual(self.ev("str(f('r','/flag'))"), "True")

    def test_split_join_jsondump_dict(self):
        self.assertEqual(self.ev("int(split(f('r','/s'), '/', 0))"), 993)
        self.assertEqual(self.ev("join(vals('r','/xs/*'), '/')"), "1/2/3/4")
        self.assertEqual(self.ev("jsondump({'b': 1, 'a': 2})"),
                         '{"a": 2, "b": 1}')
        self.assertEqual(
            self.ev("jsondump(f('r','/n'))"), '{"a": 3, "b": 4}')

    def test_joinrows_and_lets(self):
        self.assertEqual(
            self.ev("joinrows('r','/rows','{k}={v}','; ')"), "p=1; q=0; p=5")
        self.assertEqual(self.ev("L + 1", L=41), 42)

    def test_count_where(self):
        self.assertEqual(self.ev("cw('r','/rows','[[\"k\",\"==\",\"p\"]]')"), 2)
        self.assertEqual(
            self.ev('cw(\'r\',\'/rows\',\'[["k","==","p"],["v",">",1]]\')'), 1)
        # short-circuit guard before coercion (the CSV NA pattern)
        env = _env(c=[{"x": "NA"}, {"x": "3"}, {"x": "1"}])
        n = evaluate('cw(\'c\',\'\',\'[["x","not_in",["NA",""]],'
                     '["x",">=",2,"int"]]\')', env)
        self.assertEqual(n, 1)
        # unguarded coercion of "NA" is LOUD, not silent
        with self.assertRaises(ValueError):
            evaluate('cw(\'c\',\'\',\'[["x",">=",2,"int"]]\')', env)

    def test_cw_ref_membership(self):
        env = _env(rows={"detail": [{"id": "A", "m": 1}, {"id": "B", "m": 1},
                                    {"id": "C", "m": 2}]},
                   other={"erased": ["B"]})
        n = evaluate('cw(\'rows\',\'/detail\',\'[["m","==",1],'
                     '["id","not_in",{"$ref": "other:/erased"}]]\')', env)
        self.assertEqual(n, 1)

    def test_safety_rejections(self):
        for bad in ["__import__('os')", "open('/etc/passwd')",
                    "f('r','/n').__class__", "f('r','/xs')[0]",
                    "(lambda: 1)()", "f(alias='r', pointer='/n/a')",
                    "unknownname", "exec('x=1')"]:
            with self.assertRaises(ExprError, msg=bad):
                self.ev(bad)


class TestPrintedMode(unittest.TestCase):
    def test_fmt_like_bytecompat(self):
        # byte-compatible with the audited emitting scripts' fmt_like
        self.assertEqual(fmt_like("3182", 3182), "3182")
        self.assertEqual(fmt_like("1.0", 1.0), "1.0")
        self.assertEqual(fmt_like("6.7%", 6.652513268810491), "6.7")
        self.assertEqual(fmt_like("-0.69", -0.6845146732002659), "-0.68")
        self.assertEqual(fmt_like("2,978", 2978), "2978")
        self.assertEqual(fmt_like("93", 92.97533562285358), "93")
        self.assertEqual(fmt_like("54", 53.9, decimals=0), "54")
        self.assertEqual(fmt_like("x", 1.23456, decimals=2), "1.23")

    def test_compare_printed(self):
        self.assertEqual(compare_printed(6.652513268810491, "6.7%"),
                         (True, "6.7"))
        self.assertEqual(compare_printed(-0.6845146732002659, "-0.69"),
                         (False, "-0.68"))     # the printed-precision-drift class
        self.assertEqual(compare_printed("213/0", "213/0"), (True, "213/0"))
        self.assertEqual(compare_printed("95/241 = 39.4%", "95/241 = 39.4%"),
                         (True, "95/241 = 39.4"))
        self.assertEqual(compare_printed(165782, "165,782"), (True, "165782"))


class TestOutwardBand(unittest.TestCase):
    LO, HI = 66312800, 99469200          # certified_8 gross lo/hi, USD

    def band(self, qlo, qhi, g=1, scale="1e-6", lo=None, hi=None):
        return compare_outward_band(self.LO if lo is None else lo,
                                    self.HI if hi is None else hi,
                                    qlo, qhi, g, scale)

    def test_planted_outward_match(self):
        ok, emitted, viol = self.band(66, 100)
        self.assertTrue(ok)
        self.assertIsNone(viol)
        self.assertIn("[66-100]", emitted)

    def test_inward_edges_caught(self):
        ok, _, viol = self.band(67, 100)
        self.assertFalse(ok)
        self.assertIn("inward-edge", viol)
        ok, _, viol = self.band(66, 99)
        self.assertFalse(ok)
        self.assertIn("inward-edge", viol)

    def test_loose_band_caught(self):
        ok, _, viol = self.band(65, 101)
        self.assertFalse(ok)
        self.assertIn("outward-loose", viol)

    def test_planted_replay_3of6_inward(self):
        # 3-of-6 planted replay: near_tier_7
        # displayed [154.9, 232.3] against computed [154.8676, 232.3014]
        ok, _, viol = compare_outward_band(154867600, 232301400,
                                           "154.9", "232.3", "0.1", "1e-6")
        self.assertFalse(ok)
        self.assertIn("inward-edge", viol)
        self.assertIn("lo edge", viol)
        self.assertIn("hi edge", viol)     # both edges inward
        # corrected display passes
        ok, _, viol = compare_outward_band(154867600, 232301400,
                                           "154.8", "232.4", "0.1", "1e-6")
        self.assertTrue(ok)

    def test_on_grid_boundary_is_not_inward(self):
        # computed lo exactly on the grid: floor == value, quoted edge == lo
        ok, _, viol = compare_outward_band(66000000, 99469200, 66, 100,
                                           1, "1e-6")
        self.assertTrue(ok)

    def test_exactness_no_float_boundary_artifacts(self):
        # 0.1 is not a binary double; Fraction arithmetic must not misgrade
        ok, _, _ = compare_outward_band("3.3", "3.3", "3.3", "3.3", "0.1", 1)
        self.assertTrue(ok)

    def test_malformed_inputs_loud(self):
        with self.assertRaises(SpecError):
            compare_outward_band(5, 4, 5, 5, 1, 1)          # lo > hi
        with self.assertRaises(SpecError):
            compare_outward_band(1, 2, 1, 2, 0, 1)          # granularity 0
        with self.assertRaises(SpecError):
            compare_outward_band("abc", 2, 1, 2, 1, 1)      # non-number


class TestEngineEndToEnd(unittest.TestCase):
    """Planted-truth specs over synthetic receipts in a tempdir."""

    def setUp(self):
        self.td = tempfile.TemporaryDirectory()
        self.base = self.td.name
        self.receipt = {
            "total": 213, "rated": 3203,
            "peer_groups": {"p3": {"d": 55, "u": 1}, "p4": {"d": 68, "u": 1},
                            "p5": {"d": 90, "u": 2}},
            "band": {"lo": 66312800, "hi": 99469200},
            "flag": True,
        }
        self._write("main.json", self.receipt)
        with open(os.path.join(self.base, "table.csv"), "w") as fh:
            fh.write("star,q,cnt\nNA,1,NA\n4,1,3\n3,1,2\n5,0,7\n4,1,9\n")

    def tearDown(self):
        self.td.cleanup()

    def _write(self, name, obj):
        with open(os.path.join(self.base, name), "w") as fh:
            json.dump(obj, fh)

    def _spec(self, claims, receipts=None, path="spec.json"):
        spec = {"campaign": "planted", "base": self.base,
                "receipts": receipts or {"r": {"path": "main.json"},
                                         "c": {"path": "table.csv"}},
                "claims": claims}
        p = os.path.join(self.base, path)
        with open(p, "w") as fh:
            json.dump(spec, fh)
        return p

    def _run(self, claims, **kw):
        return run_spec(self._spec(claims, **kw), quiet=True)

    def test_all_match_exit_zero(self):
        rep, code = self._run([
            {"section": "S", "label": "total", "quoted": "213",
             "expr": "f('r','/total')"},
            {"section": "S", "label": "pct", "quoted": "6.7%",
             "expr": "100 * f('r','/total') / f('r','/rated')"},
            {"section": "S", "label": "per-peer", "quoted": "55/68/90",
             "format": "{a}/{b}/{c}",
             "exprs": {"a": "f('r','/peer_groups/p3/d')",
                       "b": "f('r','/peer_groups/p4/d')",
                       "c": "f('r','/peer_groups/p5/d')"}},
            {"section": "S", "label": "uniq", "quoted": "2",
             "lets": {"L": "vals('r','/peer_groups/{p3|p4|p5}/u')"},
             "expr": "len(where(L, '==', 1))"},
            {"section": "S", "label": "csv count", "quoted": "2",
             "expr": "cw('c','','[[\"cnt\",\"not_in\",[\"NA\",\"\"]],"
                     "[\"q\",\"==\",\"1\"],[\"cnt\",\">=\",3,\"int\"]]')"},
            {"section": "S", "label": "band", "mode": "outward-band",
             "quoted_lo": 66, "quoted_hi": 100, "granularity": 1,
             "scale": "1e-6",
             "lo_expr": "f('r','/band/lo')", "hi_expr": "f('r','/band/hi')"},
            {"section": "S", "label": "flag", "quoted": "True",
             "expr": "str(f('r','/flag'))"},
            {"section": "S", "label": "emitted only",
             "expr": "jsondump(f('r','/peer_groups/p3'))"},
        ])
        self.assertEqual(code, 0)
        self.assertEqual(rep["findings_count"], 0)
        self.assertEqual(rep["matches"], 7)
        self.assertEqual(rep["headlines_emitted"], 8)
        self.assertEqual(rep["table"][-1]["status"], "EMITTED")

    def test_planted_mismatch_is_finding_and_nonzero(self):
        rep, code = self._run([
            {"section": "S", "label": "total wrong", "quoted": "214",
             "expr": "f('r','/total')"},
            {"section": "S", "label": "total right", "quoted": "213",
             "expr": "f('r','/total')"},
        ])
        self.assertEqual(code, 1)
        self.assertEqual(rep["findings_count"], 1)
        self.assertEqual(rep["FINDINGS"][0]["kind"], "mismatch")
        self.assertEqual(rep["FINDINGS"][0]["label"], "total wrong")
        self.assertEqual(rep["matches"], 1)

    def test_planted_inward_band_is_finding(self):
        rep, code = self._run([
            {"section": "S", "label": "band inward", "mode": "outward-band",
             "quoted_lo": 67, "quoted_hi": 100, "granularity": 1,
             "scale": "1e-6",
             "lo_expr": "f('r','/band/lo')", "hi_expr": "f('r','/band/hi')"},
        ])
        self.assertEqual(code, 1)
        self.assertEqual(rep["FINDINGS"][0]["kind"], "outward-band")
        self.assertIn("inward-edge", rep["FINDINGS"][0]["note"])

    def test_exact_mode(self):
        rep, code = self._run([
            {"section": "S", "label": "exact ok", "mode": "exact",
             "quoted": "213", "expr": "f('r','/total')"},
            {"section": "S", "label": "exact catches register drift",
             "mode": "exact", "quoted": "213.0", "expr": "f('r','/total')"},
        ])
        self.assertEqual(code, 1)
        self.assertEqual([r["status"] for r in rep["table"]],
                         ["MATCH", "FINDING"])

    def test_missing_receipt(self):
        rep, code = self._run(
            [{"section": "S", "label": "gone", "quoted": "1",
              "expr": "f('gone','/x')"},
             {"section": "S", "label": "still evaluated", "quoted": "213",
              "expr": "f('r','/total')"}],
            receipts={"r": {"path": "main.json"},
                      "gone": {"path": "absent.json"}})
        self.assertEqual(code, 1)
        kinds = [f["kind"] for f in rep["FINDINGS"]]
        self.assertIn("missing-receipt", kinds)
        self.assertEqual(rep["table"][0]["status"], "MISSING-RECEIPT")
        self.assertEqual(rep["table"][1]["status"], "MATCH")

    def test_stale_receipt_mtime(self):
        self._write("derived.json", {"v": 1})
        self._write("input.json", {"v": 1})
        old = time.time() - 3600
        os.utime(os.path.join(self.base, "derived.json"), (old, old))
        rep, code = self._run(
            [{"section": "S", "label": "v", "quoted": "1",
              "expr": "f('d','/v')"}],
            receipts={"d": {"path": "derived.json",
                            "newer_than": ["input.json"]}})
        self.assertEqual(code, 1)
        self.assertEqual(rep["FINDINGS"][0]["kind"], "stale-receipt")
        # claim itself still evaluates and matches
        self.assertEqual(rep["table"][0]["status"], "MATCH")
        # negative control: fresh ordering -> clean
        now = time.time()
        os.utime(os.path.join(self.base, "derived.json"), (now, now))
        rep, code = self._run(
            [{"section": "S", "label": "v", "quoted": "1",
              "expr": "f('d','/v')"}],
            receipts={"d": {"path": "derived.json",
                            "newer_than": ["input.json"]}})
        self.assertEqual(code, 0)

    def test_stale_receipt_max_age(self):
        self._write("aged.json", {"v": 1})
        old = time.time() - 10 * 86400
        os.utime(os.path.join(self.base, "aged.json"), (old, old))
        rep, code = self._run(
            [{"section": "S", "label": "v", "quoted": "1",
              "expr": "f('a','/v')"}],
            receipts={"a": {"path": "aged.json", "max_age_days": 5}})
        self.assertEqual(code, 1)
        self.assertEqual(rep["FINDINGS"][0]["kind"], "stale-receipt")

    def test_eval_error_is_loud_finding(self):
        rep, code = self._run([
            {"section": "S", "label": "schema drifted", "quoted": "1",
             "expr": "f('r','/no_such_field')"},
            {"section": "S", "label": "unaffected", "quoted": "213",
             "expr": "f('r','/total')"},
        ])
        self.assertEqual(code, 1)
        self.assertEqual(rep["FINDINGS"][0]["kind"], "eval-error")
        self.assertEqual(rep["table"][1]["status"], "MATCH")

    def test_note_format(self):
        rep, _ = self._run([
            {"section": "S", "label": "n", "quoted": "213",
             "expr": "f('r','/total')",
             "note_format": "flag is {h}",
             "note_exprs": {"h": "str(f('r','/flag'))"}},
        ])
        self.assertEqual(rep["table"][0]["note"], "flag is True")

    def test_report_written(self):
        out = os.path.join(self.base, "report.json")
        run_spec(self._spec([{"section": "S", "label": "t", "quoted": "213",
                              "expr": "f('r','/total')"}]),
                 out_path=out, quiet=True)
        with open(out) as fh:
            rep = json.load(fh)
        self.assertEqual(rep["run"], "emitall")
        self.assertEqual(rep["matches"], 1)
        self.assertTrue(rep["date_utc"].endswith("Z"))


class TestCLI(unittest.TestCase):
    def _spec_file(self, td, quoted):
        with open(os.path.join(td, "r.json"), "w") as fh:
            json.dump({"x": 5}, fh)
        p = os.path.join(td, "s.json")
        with open(p, "w") as fh:
            json.dump({"base": td, "receipts": {"r": {"path": "r.json"}},
                       "claims": [{"section": "S", "label": "x",
                                   "quoted": quoted,
                                   "expr": "f('r','/x')"}]}, fh)
        return p

    def test_exit_codes(self):
        with tempfile.TemporaryDirectory() as td:
            ok = subprocess.run([sys.executable, RUN_PY,
                                 self._spec_file(td, "5"), "--quiet"])
            self.assertEqual(ok.returncode, 0)
        with tempfile.TemporaryDirectory() as td:
            bad = subprocess.run([sys.executable, RUN_PY,
                                  self._spec_file(td, "6"), "--quiet"])
            self.assertEqual(bad.returncode, 1)
        with tempfile.TemporaryDirectory() as td:
            p = os.path.join(td, "malformed.json")
            with open(p, "w") as fh:
                json.dump({"claims": []}, fh)   # no receipts key
            mal = subprocess.run([sys.executable, RUN_PY, p, "--quiet"],
                                 capture_output=True)
            self.assertEqual(mal.returncode, 2)


class TestToySpecReference(unittest.TestCase):
    """End-to-end on the SHIPPED example deployment (specs/toy_stars.json +
    specs/toy_receipts/): the spec must recount clean as committed, and
    planted corruptions of it must fire — proving the shipped example is a
    live check, not decoration."""

    def _patched(self, mutate):
        """load the toy spec, apply mutate(spec), run from a tempdir with
        base pinned back at specs/ so the receipts still resolve."""
        with open(os.path.join(SPECS, "toy_stars.json")) as fh:
            spec = json.load(fh)
        spec["base"] = SPECS
        mutate(spec)
        with tempfile.TemporaryDirectory() as td:
            p = os.path.join(td, "patched.json")
            with open(p, "w") as fh:
                json.dump(spec, fh)
            return run_spec(p, quiet=True)

    def test_toy_spec_recounts_clean(self):
        rep, code = run_spec(os.path.join(SPECS, "toy_stars.json"),
                             quiet=True)
        self.assertEqual(code, 0)
        self.assertEqual(rep["findings_count"], 0)
        self.assertEqual(rep["headlines_emitted"], 12)
        self.assertEqual(rep["matches"], 11)
        self.assertEqual(rep["table"][-1]["status"], "EMITTED")

    def test_planted_transcription_drift_fires(self):
        def mutate(spec):
            n = 0
            for c in spec["claims"]:
                if c["label"] == "star changes":
                    c["quoted"] = "6"          # receipt says 5
                    n += 1
            assert n == 1
        rep, code = self._patched(mutate)
        self.assertEqual(code, 1)
        self.assertEqual(rep["findings_count"], 1)
        self.assertEqual(rep["FINDINGS"][0]["kind"], "mismatch")
        self.assertEqual(rep["FINDINGS"][0]["label"], "star changes")
        self.assertEqual(rep["matches"], 10)

    def test_planted_inward_bands_fire(self):
        # nearest-rounding replay: print each band's edges nearest instead
        # of outward; the engine must fire inward-edge findings on both bands
        def mutate(spec):
            n = 0
            for c in spec["claims"]:
                if c["label"] == "payout band, integer $K":
                    c["quoted_lo"], c["quoted_hi"] = 1235, 2765
                    n += 1
                if c["label"] == "payout band, 0.1 $M":
                    c["quoted_lo"], c["quoted_hi"] = "1.3", "2.7"
                    n += 1
            assert n == 2
        rep, code = self._patched(mutate)
        self.assertEqual(code, 1)
        bands = [f for f in rep["FINDINGS"] if f["kind"] == "outward-band"]
        self.assertEqual(len(bands), 2)
        for f in bands:
            self.assertIn("inward-edge", f["note"])


class TestLintInterpolatedLiterals(unittest.TestCase):
    """Interpolated-literal planted truth: every other number in the headline
    is interpolated; the one typed '5' understates the receipt's $0.09."""

    PENNY = '''
        rec["emitted_headlines"] = [
            "DS2017_updated (2021 reissue) carries CENT-SCALE internal identity "
            "violations (after - before != adjustment) in "
            f"{L['DS2017_updated']['nonmerged']['identity_adj_eq_after_minus_before']['n_cent_scale']}"
            " nonmerged + "
            f"{L['DS2017_updated']['exiting_2017']['identity_adj_eq_after_minus_before']['n_cent_scale']}"
            " exiting-issuer rows (max 5 cents) -- the reissued file is "
            "internally inconsistent at its own money register; the original "
            "edition has zero such violations",
        ]
    '''

    def test_typed_digits_caught_placeholders_exempt(self):
        findings, allowed = _lint(self.PENNY)
        self.assertEqual(allowed, [])
        self.assertTrue(all(f["kind"] == "digit-literal" for f in findings))
        runs = sorted(f["detail"].split("'")[1] for f in findings)
        # exactly the three typed runs: the two year labels and THE '5'
        self.assertEqual(runs, ["2017", "2021", "5"])
        self.assertTrue(any("max 5 cents" in f["detail"] for f in findings))
        # interpolated receipt fields (subscript keys) are never flagged
        self.assertFalse(any("n_cent_scale" in f["detail"] for f in findings))

    def test_fully_interpolated_is_clean(self):
        findings, _ = _lint('''
            rec["emitted_headlines"] = [
                f"{ed} ({yr} reissue) carries CENT-SCALE internal identity "
                "violations (after - before != adjustment) in "
                f"{n_nm} nonmerged + {n_ex} exiting-issuer rows "
                f"(max {mx} cents) -- the reissued file is internally "
                "inconsistent at its own money register",
            ]
        ''')
        self.assertEqual(findings, [])


class TestLintUnseededSampler(unittest.TestCase):
    """Unseeded-sampler planted truth: PYTHONHASHSEED-dependent emission
    order — a set-intersection comprehension feeding a headline, plus an
    unseeded random.sample; the lint catches the pattern at the source."""

    SAMPLER = '''
        import random
        ch = {k: (o[k], u[k]) for k in set(o) & set(u) if o[k] != u[k]}
        rec["emitted_headlines"] = [
            "changed cells: " + ", ".join(f"{k}->{v}" for k, v in ch.items()),
            "witness pools: " + ", ".join(str(p) for p in set(o) & set(u)),
        ]
        examples = random.sample(list(pool_ids), k)
        rec["emitted_headlines"].append("examples: " + ", ".join(examples))
    '''

    def test_hash_order_and_unseeded_sampler_caught(self):
        findings, _ = _lint(self.SAMPLER)
        kinds = sorted(f["kind"] for f in findings)
        self.assertEqual(kinds, ["unordered-iteration", "unordered-iteration",
                                 "unseeded-sampling"])
        details = " | ".join(f["detail"] for f in findings)
        self.assertIn("dict-view", details)      # ch.items() built from sets
        self.assertIn("set iteration", details)  # set(o) & set(u) directly

    def test_sorted_and_seeded_is_clean(self):
        findings, _ = _lint('''
            import random
            ch = {k: (o[k], u[k]) for k in set(o) & set(u) if o[k] != u[k]}
            rec["emitted_headlines"] = [
                "changed cells: " + ", ".join(
                    f"{k}->{v}" for k, v in sorted(ch.items())),
                "witness pools: " + ", ".join(
                    str(p) for p in sorted(set(o) & set(u))),
            ]
            rng = random.Random(seed)
            examples = rng.sample(sorted(pool_ids), k)
            rec["emitted_headlines"].append("examples: " + ", ".join(examples))
        ''')
        self.assertEqual(findings, [])

    def test_seed_must_precede_the_draw(self):
        clean, _ = _lint('''
            import random
            random.seed(seed)
            xs = random.sample(pool, k)
        ''')
        self.assertEqual(clean, [])
        late, _ = _lint('''
            import random
            xs = random.sample(pool, k)
            random.seed(seed)
        ''')
        self.assertEqual([f["kind"] for f in late], ["unseeded-sampling"])
        bare_rng, _ = _lint('''
            import random
            rng = random.Random()
            xs = rng.choice(pool)
        ''')
        self.assertEqual([f["kind"] for f in bare_rng], ["unseeded-sampling"])


class TestLintSinksAndTemplates(unittest.TestCase):
    def test_non_emission_strings_not_flagged(self):
        findings, _ = _lint('''
            path = "runs/certified_at_132.json"
            cutoff_note = "strict > pinned in 6 of 6"
            table[3] = load(path)
        ''')
        self.assertEqual(findings, [])

    def test_dict_key_and_keyword_sinks(self):
        findings, _ = _lint('''
            summary = {"headline_dp": f"narrower in 379 of {n} raw"}
            bank(receipt=dict(emitted_headlines=["66 pools closed"]))
        ''')
        self.assertEqual([f["kind"] for f in findings],
                         ["digit-literal", "digit-literal"])

    def test_emit_function_args_are_sinks(self):
        # strict per the emission rule: EVERY digit run in a literal fed to a
        # registered emit function fires — including section labels like
        # "CP2" and transcribed quote strings ("213"); the fix is a pragma
        # with a reason or migration of the quotes into a claims spec (JSON
        # specs are data, not linted Python)
        findings, _ = _lint('''
            emit("CP2", "2026 census differs", "213", cc["total"])
        ''')
        runs = sorted(f["detail"].split("'")[1] for f in findings)
        self.assertEqual(runs, ["2", "2026", "213"])

    def test_format_and_percent_placeholders_exempt(self):
        findings, _ = _lint('''
            headlines.append("{} pools (max {} cents), share {v:.4f}".format(a, b, v=v))
            headlines.append("%05d rows in %.3f s" % (n, t))
            headlines.append(f"{x:.3f} share of {n}")
        ''')
        self.assertEqual(findings, [])
        findings, _ = _lint('''
            headlines.append("max 5 cents over %d rows" % n)
            headlines.append("max 5 cents over {} rows".format(n))
        ''')
        self.assertEqual([f["kind"] for f in findings],
                         ["digit-literal", "digit-literal"])
        self.assertTrue(all("'5'" in f["detail"] for f in findings))

    def test_for_loop_over_set_writing_sink(self):
        findings, _ = _lint('''
            for pool in set(pools):
                emitted.append(f"pool {pool} closed")
        ''')
        self.assertEqual([f["kind"] for f in findings],
                         ["unordered-iteration"])
        findings, _ = _lint('''
            for pool in sorted(set(pools)):
                emitted.append(f"pool {pool} closed")
        ''')
        self.assertEqual(findings, [])
        # same loop NOT writing a sink is out of scope for the emission rule
        findings, _ = _lint('''
            for pool in set(pools):
                totals[pool] = compute(pool)
        ''')
        self.assertEqual(findings, [])

    def test_order_insensitive_reducers_exempt(self):
        findings, _ = _lint('''
            summary["headline_n"] = f"{sum(1 for p in set(pools))} pools, "
            summary["headline_m"] = f"max {max(set(vals))}"
        ''')
        self.assertEqual(findings, [])

    def test_return_from_headline_named_function(self):
        findings, _ = _lint('''
            def build_headline(n):
                return f"{n} of 117 off-optimum"
        ''')
        self.assertEqual([f["kind"] for f in findings], ["digit-literal"])
        self.assertIn("'117'", findings[0]["detail"])

    def test_augassign_and_extend_sinks(self):
        findings, _ = _lint('''
            headlines += [f"tail note ({yr}): {v} rows kept 9dp"]
            emitted.extend(["all 3 reissuances characterized"])
        ''')
        runs = sorted(f["detail"].split("'")[1] for f in findings)
        self.assertEqual(runs, ["3", "9"])


class TestLintPragma(unittest.TestCase):
    def test_pragma_with_reason_honored_and_recorded(self):
        findings, allowed = _lint('''
            headlines.append(
                "guardrail cap 3182 hospitals"  # emitall-lint: allow-literal cap constant quoted from the method pin, not a result
            )
        ''')
        self.assertEqual(findings, [])
        self.assertEqual(len(allowed), 1)
        self.assertEqual(allowed[0]["runs"], ["3182"])
        self.assertIn("method pin", allowed[0]["reason"])

    def test_pragma_without_reason_is_a_finding_and_suppresses_nothing(self):
        findings, allowed = _lint('''
            headlines.append("guardrail cap 3182 hospitals")  # emitall-lint: allow-literal
        ''')
        self.assertEqual(allowed, [])
        kinds = sorted(f["kind"] for f in findings)
        self.assertEqual(kinds, ["digit-literal", "pragma-missing-reason"])


class TestLintCLI(unittest.TestCase):
    def _write(self, td, name, body):
        p = os.path.join(td, name)
        with open(p, "w") as fh:
            fh.write(textwrap.dedent(body))
        return p

    def test_exit_codes_and_json_report(self):
        with tempfile.TemporaryDirectory() as td:
            dirty = self._write(td, "dirty.py", '''
                headlines = ["max 5 cents residue"]
            ''')
            clean = self._write(td, "clean.py", '''
                headlines = [f"max {mx} cents residue"]
            ''')
            out = os.path.join(td, "lint.json")
            r = subprocess.run([sys.executable, LINT_PY, dirty,
                                "--json", out, "--quiet"],
                               capture_output=True, text=True)
            self.assertEqual(r.returncode, 1)
            with open(out) as fh:
                rep = json.load(fh)
            self.assertEqual(rep["run"], "emitall-lint")
            self.assertEqual(rep["counts"], {"digit-literal": 1})
            self.assertTrue(rep["date_utc"].endswith("Z"))
            r = subprocess.run([sys.executable, LINT_PY, clean, "--quiet"],
                               capture_output=True, text=True)
            self.assertEqual(r.returncode, 0)

    def test_directory_scan(self):
        with tempfile.TemporaryDirectory() as td:
            self._write(td, "a.py", 'headlines = ["9dp era"]\n')
            self._write(td, "b.py", 'x = "not an emission 123"\n')
            os.mkdir(os.path.join(td, "__pycache__"))
            self._write(os.path.join(td, "__pycache__"), "junk.py",
                        'headlines = ["55 skipped"]\n')
            r = subprocess.run([sys.executable, LINT_PY, td],
                               capture_output=True, text=True)
            self.assertEqual(r.returncode, 1)
            self.assertIn("1 digit-literal", r.stdout)
            self.assertIn("2 files scanned", r.stdout)

    def test_parse_error_is_loud(self):
        with tempfile.TemporaryDirectory() as td:
            bad = self._write(td, "bad.py", "def broken(:\n")
            r = subprocess.run([sys.executable, LINT_PY, bad, "--quiet"],
                               capture_output=True, text=True)
            self.assertEqual(r.returncode, 1)
            self.assertIn("parse-error", r.stdout)


class TestV11Expr(unittest.TestCase):
    """zipjoin, contains predicate, frac()."""

    def setUp(self):
        self.env = Env({"r": {"yrs": ["FY2019", "FY2020", "FY2021"],
                              "regs": ["9dp", "9dp", "raw"],
                              "rows": [{"tag": "restated-final"},
                                       {"tag": "original"},
                                       {"tag": "restated-draft"}]}})

    def ev(self, e):
        return evaluate(e, self.env)

    def test_zipjoin_rowwise_pairing(self):
        self.assertEqual(
            self.ev("zipjoin(vals('r','/yrs/*'), vals('r','/regs/*'),"
                    " '{a}={b}', '; ')"),
            "FY2019=9dp; FY2020=9dp; FY2021=raw")

    def test_zipjoin_length_mismatch_loud(self):
        with self.assertRaises(ExprError):
            self.ev("zipjoin(vals('r','/yrs/*'), ['x'], '{a}{b}', ',')")
        with self.assertRaises(ExprError):
            self.ev("zipjoin(f('r','/yrs/0'), vals('r','/regs/*'), '{a}', ',')")

    def test_contains_predicate_where_and_cw(self):
        self.assertEqual(
            self.ev("len(where(vals('r','/yrs/*'), 'contains', 'FY202'))"), 2)
        self.assertEqual(
            self.ev("where(vals('r','/regs/*'), 'not_contains', 'dp')"),
            ["raw"])
        self.assertEqual(
            self.ev('cw(\'r\',\'/rows\',\'[["tag","contains","restated"]]\')'), 2)
        with self.assertRaises(ExprError):
            self.ev("where([1, 2], 'contains', 'x')")   # non-container: LOUD

    def test_frac_exact_interval_widths(self):
        self.assertEqual(self.ev("str(frac('232.4') - frac('154.8'))"),
                         "388/5")
        self.assertIs(self.ev("frac('1/3') + frac('1/6') == frac('1/2')"),
                      True)
        # decimal strings parse exactly, never through a binary double
        self.assertIs(self.ev("frac('0.1') == frac(1) / frac(10)"), True)
        self.assertEqual(self.ev("float(frac('77.6'))"), 77.6)
        with self.assertRaises(ExprError):
            self.ev("frac('not a number')")

    def test_printed_mode_formats_fractions(self):
        self.assertEqual(compare_printed(Fraction(388, 5), "77.6"),
                         (True, "77.6"))
        self.assertEqual(compare_printed(Fraction(1, 3), "0.33"),
                         (True, "0.33"))
        # exact mode remains exact p/q via str()
        self.assertEqual(str(Fraction(388, 5)), "388/5")


class TestPaperClaimsBattery(unittest.TestCase):
    """battery.py --selftest rides the registered pytest battery: the gated
    fixture battery (every row class catches its planted failure and passes
    its clean twin, comment stripping, wrap-proof join, stale-render refusal)
    must go 8/8 green under both python3 and python3 -O."""

    def _run_selftest(self, extra_flags=()):
        return subprocess.run(
            [sys.executable, *extra_flags, BATTERY_PY, "--selftest"],
            capture_output=True, text=True)

    def test_selftest_green(self):
        p = self._run_selftest()
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        self.assertIn("8/8 gates green", p.stdout)

    def test_selftest_green_under_O(self):
        p = self._run_selftest(extra_flags=("-O",))
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        self.assertIn("8/8 gates green", p.stdout)


if __name__ == "__main__":
    unittest.main(verbosity=2)

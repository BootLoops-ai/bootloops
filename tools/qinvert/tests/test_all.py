"""qinvert executable receipts.

Every test is a planted-truth or negative-control receipt:
- evaluator receipts against EXTERNAL reference values (R's documented
  quantile(1:10, ...) outputs per Hyndman-Fan type; SAS QNTLDEF hand
  computations; the originating analysis's validated qntldef5 read directly
  from its source tree (env QINVERT_ORIGIN_ROOT, read-only; not included);
- planted-truth recovery batteries at every supported quantile definition
  (plant a known add/removal, publish the perturbed statistics, require the
  solver to recover the plant among its minimal families);
- negative controls (unreachable targets must yield empty results);
- parity-constraint units (the n mod 4 averaging argument, generalized);
- an end-to-end stacking scenario (3 datasets sharing 2 planted entities,
  plus an exact dataset guarded by the no-regression gate);
- a REGRESSION test reproducing the originating analysis's 2024 D02/MA-PD
  +1-addition family vs the published 1.41, read-only from that tree
  (receipt: its membership_families.json run record, not included).

Run:  python -m pytest tests/test_all.py -q   (from the qinvert directory)
  or  python tests/test_all.py
"""

import csv
import json
import os
import random
import sys
import unittest
from fractions import Fraction as F

_HERE = os.path.dirname(os.path.abspath(__file__))
_TOOLS = os.path.dirname(os.path.dirname(_HERE))
if _TOOLS not in sys.path:
    sys.path.insert(0, _TOOLS)

from qinvert import (ALL_DEFINITIONS, Dataset, Interval, MaxStat, Mean,
                     MinStat, OrderStat, Quantile, RangeStat, TukeyFence,
                     invert_fence_pair, solve, solve_many, stack,
                     targets_exact)

ORIGIN = __import__("os").environ.get("QINVERT_ORIGIN_ROOT", "")  # originating-analysis tree (env; unset/absent -> those tests SKIP)


def _import_origin_symbol(module, symbol):
    """Import from the origin tree STRICTLY read-only: bytecode writing is
    disabled for the duration so that tree cannot be touched."""
    src = os.path.join(ORIGIN, "src")
    prev = sys.dont_write_bytecode
    sys.dont_write_bytecode = True
    sys.path.insert(0, src)
    try:
        mod = __import__(module)
        return getattr(mod, symbol)
    finally:
        sys.path.remove(src)
        sys.dont_write_bytecode = prev


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def window_contains(win, v):
    lo, hi, lo_s, hi_s = win
    if lo is not None and (v < lo or (lo_s and v == lo)):
        return False
    if hi is not None and (v > hi or (hi_s and v == hi)):
        return False
    return True


def family_covers(fam, v):
    """Does this family's add set cover the value v (pinned or in-window)?"""
    for s in fam["adds"]:
        if s["window"] is None:
            if s["value"] == v:
                return True
        elif window_contains(s["window"], v):
            return True
    return False


def tie_free_data(seed, n, lo=500, hi=9500):
    """n distinct 2dp values in [lo/100, hi/100]."""
    rng = random.Random(seed)
    return sorted(F(x, 100) for x in rng.sample(range(lo, hi + 1), n))


def materialize(base_sorted, fam):
    vals = list(base_sorted)
    for r in fam["removed"]:
        vals.remove(r)
    return sorted(vals + [s["value"] for s in fam["adds"]])


# ---------------------------------------------------------------------------
# 1. evaluator receipts
# ---------------------------------------------------------------------------

class TestQuantileEvaluators(unittest.TestCase):

    def test_r_reference_values_1_to_10(self):
        """R's quantile(1:10, 0.25, type=1..9) documented outputs
        (Hyndman & Fan 1996 taxonomy; R ?quantile).  External planted truth
        for all nine R types."""
        xs = [F(i) for i in range(1, 11)]
        expected = {
            "r1": F(3), "r2": F(3), "r3": F(2), "r4": F(5, 2), "r5": F(3),
            "r6": F(11, 4), "r7": F(13, 4), "r8": F(35, 12), "r9": F(47, 16),
        }
        for d, want in expected.items():
            got = Quantile(F(1, 4), d).value(xs)
            self.assertEqual(got, want, f"{d}: {got} != {want}")
        # medians: r7 -> 5.5 ; r2 -> 5.5 ; r1 -> 5
        self.assertEqual(Quantile(F(1, 2), "r7").value(xs), F(11, 2))
        self.assertEqual(Quantile(F(1, 2), "r2").value(xs), F(11, 2))
        self.assertEqual(Quantile(F(1, 2), "r1").value(xs), F(5))

    def test_sas_hand_values(self):
        """SAS QNTLDEF hand computations (Base SAS 9.4 Proc Guide,
        'Calculating Percentiles') + the origin tree's tukey.py spot checks."""
        xs4 = [F(1), F(2), F(3), F(4)]
        # QNTLDEF=5, n=4: np=1 integer -> (x1+x2)/2 (tukey.py receipt)
        self.assertEqual(Quantile(F(1, 4), "sas5").value(xs4), F(3, 2))
        self.assertEqual(Quantile(F(3, 4), "sas5").value(xs4), F(7, 2))
        xs5 = [F(i) for i in range(1, 6)]
        self.assertEqual(Quantile(F(1, 4), "sas5").value(xs5), F(2))
        self.assertEqual(Quantile(F(3, 4), "sas5").value(xs5), F(4))
        # n=10, p=1/4: np=2.5, j=2, g=1/2
        xs = [F(i) for i in range(1, 11)]
        self.assertEqual(Quantile(F(1, 4), "sas1").value(xs), F(5, 2))
        self.assertEqual(Quantile(F(1, 4), "sas2").value(xs), F(2))  # tie->even j
        self.assertEqual(Quantile(F(1, 4), "sas3").value(xs), F(3))
        self.assertEqual(Quantile(F(1, 4), "sas4").value(xs), F(11, 4))
        self.assertEqual(Quantile(F(1, 4), "sas5").value(xs), F(3))

    def test_documented_equivalences(self):
        """SAS1=R4, SAS2=R3, SAS3=R1, SAS4=R6, SAS5=R2 on random exact data
        (H&F sec. 2 / SAS doc)."""
        pairs = [("sas1", "r4"), ("sas2", "r3"), ("sas3", "r1"),
                 ("sas4", "r6"), ("sas5", "r2")]
        rng = random.Random(7)
        ps = [F(1, 4), F(1, 2), F(3, 4), F(1, 3), F(2, 5), F(9, 10)]
        for n in (7, 8, 12, 13, 40):
            xs = sorted(F(rng.randint(0, 5000), 100) for _ in range(n))
            for p in ps:
                for a, b in pairs:
                    va = Quantile(p, a).value(xs)
                    vb = Quantile(p, b).value(xs)
                    self.assertEqual(va, vb, f"{a} vs {b} at n={n} p={p}")

    def test_definitions_not_all_aliased(self):
        """Negative control on the evaluator battery itself: distinct
        definitions must DISAGREE somewhere (guards against one
        implementation masquerading as fourteen)."""
        xs = [F(i) for i in range(1, 11)]
        vals = {d: Quantile(F(1, 4), d).value(xs) for d in ALL_DEFINITIONS}
        self.assertGreater(len(set(vals.values())), 4, vals)

    def test_origin_qntldef5_crosscheck(self):
        """Our sas5 vs the origin tree's validated qntldef5, read-only
        (origin receipt: tukey.py, validated against CMS K-5/K-6 tables)."""
        src = os.path.join(ORIGIN, "src")
        if not os.path.isdir(src):
            self.skipTest(f"origin tree not present at {src}")
        qntldef5 = _import_origin_symbol("tukey", "qntldef5")
        rng = random.Random(11)
        for trial in range(60):
            n = rng.randint(4, 90)
            xs = sorted(F(rng.randint(0, 9000), 100) for _ in range(n))
            for p in (F(1, 4), F(1, 2), F(3, 4), F(2, 5)):
                self.assertEqual(Quantile(p, "sas5").value(xs),
                                 qntldef5(xs, p), f"n={n} p={p}")


# ---------------------------------------------------------------------------
# 2. parity constraints (the n mod 4 argument, generalized per definition)
# ---------------------------------------------------------------------------

class TestParity(unittest.TestCase):

    def test_sas5_quartile_averaging_iff_n_mod_4(self):
        """QNTLDEF=5 at p=1/4: averaging (two-index support with weights
        1/2,1/2 at n/4, n/4+1) fires iff n % 4 == 0.
        Provenance: the origin inverse_membership.py ('averaging
        occurs iff n % 4 == 0'), validated there by its own battery."""
        q = Quantile(F(1, 4), "sas5")
        for n in range(8, 81):
            sup = q.support(n)
            if n % 4 == 0:
                self.assertEqual(len(sup), 2, n)
                self.assertEqual(sup[0], (n // 4, F(1, 2)), n)
                self.assertEqual(sup[1], (n // 4 + 1, F(1, 2)), n)
            else:
                self.assertEqual(len(sup), 1, n)
                self.assertEqual(sup[0][0], -(-n // 4), n)  # ceil(n/4)

    def test_parity_generalized_step_and_avg(self):
        """For p = a/b in lowest terms, the boundary event (np integer)
        happens iff b | n: sas5/r2 then average, sas3/r1 then step down to
        x_(np).  This is the per-definition generalization of the origin's
        n mod 4 argument."""
        cases = [F(1, 4), F(3, 4), F(1, 2), F(1, 3), F(2, 5)]
        for p in cases:
            b = p.denominator
            for n in range(b + 1, 90):
                np_int = (n * p).denominator == 1
                self.assertEqual(np_int, n % b == 0, (p, n))
                for d in ("sas5", "r2"):
                    sup = Quantile(p, d).support(n)
                    self.assertEqual(len(sup) == 2, np_int and
                                     1 <= n * p < n, (d, p, n, sup))
                for d in ("sas3", "r1"):
                    sup = Quantile(p, d).support(n)
                    j = sup[0][0]
                    if np_int:
                        self.assertEqual(j, n * p, (d, p, n))
                    else:
                        num = (n * p).numerator
                        den = (n * p).denominator
                        self.assertEqual(j, -(-num // den), (d, p, n))

    def test_r7_boundary(self):
        """Continuous type: r7 interpolates iff (n-1)p is not an integer."""
        p = F(1, 4)
        for n in range(5, 60):
            sup = Quantile(p, "r7").support(n)
            h = (n - 1) * p + 1
            self.assertEqual(len(sup) == 1, h.denominator == 1, (n, sup))

    def test_averaging_value(self):
        xs = [F(i) for i in range(1, 13)]  # n=12, 4 | 12
        self.assertEqual(Quantile(F(1, 4), "sas5").value(xs), F(7, 2))


# ---------------------------------------------------------------------------
# 3. fence algebra
# ---------------------------------------------------------------------------

class TestFenceAlgebra(unittest.TestCase):

    def test_inversion_round_trip(self):
        rng = random.Random(3)
        for mult in (F(3), F(3, 2)):
            for _ in range(50):
                q1 = F(rng.randint(0, 4000), 100)
                q3 = q1 + F(rng.randint(1, 3000), 100)
                lo = (1 + mult) * q1 - mult * q3
                hi = (1 + mult) * q3 - mult * q1
                r1, r3, mode = invert_fence_pair(lo, hi, mult, None, None)
                self.assertEqual(mode, "both")
                self.assertEqual((r1, r3), (q1, q3))

    def test_capped_modes(self):
        # published lo equals cap -> only hi equation survives
        self.assertEqual(invert_fence_pair("0", "1.41", 3, F(0), None)[2],
                         "hi_only")
        self.assertEqual(invert_fence_pair("57", "100", 3, F(0), F(100))[2],
                         "lo_only")
        self.assertEqual(invert_fence_pair("0", "100", 3, F(0), F(100))[2],
                         "none")
        self.assertEqual(invert_fence_pair(None, "88", 3, F(0), F(100))[2],
                         "hi_only")

    def test_fence_value_matches_manual(self):
        xs = tie_free_data(5, 41)
        q1 = Quantile(F(1, 4), "sas5").value(xs)
        q3 = Quantile(F(3, 4), "sas5").value(xs)
        lo = TukeyFence("lo", 3, "sas5", cap=0).value(xs)
        hi = TukeyFence("hi", 3, "sas5").value(xs)
        self.assertEqual(lo, max(4 * q1 - 3 * q3, F(0)))
        self.assertEqual(hi, 4 * q3 - 3 * q1)


# ---------------------------------------------------------------------------
# 4. planted-truth recovery battery: every supported definition
# ---------------------------------------------------------------------------

class TestPlantedAddAllDefinitions(unittest.TestCase):

    def _run_definition(self, definition):
        xs = tie_free_data(101, 57)
        rng = random.Random(hash(definition) & 0xFFFF)
        planted = None
        for _ in range(200):
            cand = F(rng.randint(500, 9500), 100)
            if cand in xs:
                continue
            full = sorted(xs + [cand])
            q1p = Quantile(F(1, 4), definition).value(full)
            q3p = Quantile(F(3, 4), definition).value(full)
            q1b = Quantile(F(1, 4), definition).value(xs)
            q3b = Quantile(F(3, 4), definition).value(xs)
            if (q1p, q3p) != (q1b, q3b):
                planted = cand
                break
        self.assertIsNotNone(planted, f"{definition}: no moving plant found")
        targets = [(Quantile(F(1, 4), definition), q1p),
                   (Quantile(F(3, 4), definition), q3p)]
        res = solve(xs, targets, register=F(1, 100), value_lo=0,
                    value_hi=100, k_cap=2, rem_cap=0, max_witnesses=64)
        self.assertEqual(res["k"], 1, f"{definition}: k={res['k']} "
                                      f"caps={res['caps_hit']}")
        # every family must re-verify independently
        for fam in res["families"]:
            after = materialize(xs, fam)
            self.assertTrue(targets_exact(after, targets), (definition, fam))
        # the planted delta must be covered by some minimal family
        covered = any(family_covers(f, planted) for f in res["families"])
        self.assertTrue(covered, f"{definition}: planted {planted} not in "
                                 f"{[(f['adds']) for f in res['families']]}")

    def test_all_14_definitions(self):
        for d in ALL_DEFINITIONS:
            with self.subTest(definition=d):
                self._run_definition(d)


class TestPlantedRemoval(unittest.TestCase):

    def _run(self, definition):
        xs = tie_free_data(202, 61)
        victim = xs[45]
        reduced = [v for v in xs if v != victim]
        q1 = Quantile(F(1, 4), definition).value(reduced)
        q3 = Quantile(F(3, 4), definition).value(reduced)
        # ensure removal moves at least one target
        self.assertNotEqual(
            (q1, q3),
            (Quantile(F(1, 4), definition).value(xs),
             Quantile(F(3, 4), definition).value(xs)),
            f"{definition}: victim did not move the targets; pick another")
        targets = [(Quantile(F(1, 4), definition), q1),
                   (Quantile(F(3, 4), definition), q3)]
        res = solve(xs, targets, register=F(1, 100), value_lo=0,
                    value_hi=100, k_cap=2, rem_cap=1, max_witnesses=32)
        self.assertEqual(res["k"], 1, definition)
        hit = [f for f in res["families"]
               if f["k_rem"] == 1 and f["removed"] == [victim]
               and f["k_add"] == 0]
        self.assertTrue(hit, f"{definition}: planted removal not recovered: "
                             f"{[(f['removed'], len(f['adds'])) for f in res['families']]}")

    def test_sas5(self):
        self._run("sas5")

    def test_r7(self):
        self._run("r7")


# ---------------------------------------------------------------------------
# 5. negative controls
# ---------------------------------------------------------------------------

class TestNegativeControls(unittest.TestCase):

    def test_unreachable_quantiles_empty(self):
        xs = tie_free_data(303, 45, lo=4000, hi=9900)
        targets = [(Quantile(F(1, 4), "sas5"), F(10)),
                   (Quantile(F(3, 4), "sas5"), F(200))]
        res = solve(xs, targets, register=F(1, 100), value_lo=0,
                    value_hi=100, k_cap=3, rem_cap=1)
        self.assertIsNone(res["k"])
        self.assertEqual(res["families"], [])
        # loudness contract: a witness-free search always records its depth
        self.assertIn("k_cap-exhausted", res["caps_hit"])

    def test_inverted_quartiles_empty(self):
        xs = tie_free_data(304, 45)
        targets = [(Quantile(F(1, 4), "sas5"), F(80)),
                   (Quantile(F(3, 4), "sas5"), F(20))]
        res = solve(xs, targets, register=F(1, 100), value_lo=0,
                    value_hi=100, k_cap=2, rem_cap=1)
        self.assertIsNone(res["k"])
        self.assertIn("k_cap-exhausted", res["caps_hit"])

    def test_unreachable_fence_empty(self):
        """The origin tree's negative control: fences (0, 99999) unreachable
        by any small delta (inverse_membership.py __main__)."""
        xs = tie_free_data(305, 121, lo=4000, hi=9900)
        targets = [(TukeyFence("lo", 3, "sas5", cap=0), F(0)),
                   (TukeyFence("hi", 3, "sas5"), F(99999))]
        res = solve(xs, targets, register=F(1, 100), value_lo=0,
                    value_hi=100, k_cap=2, rem_cap=0)
        self.assertIsNone(res["k"])
        self.assertEqual(res["families"], [])
        self.assertIn("k_cap-exhausted", res["caps_hit"])


# ---------------------------------------------------------------------------
# 5b. same-value blocks (equal pinned values): a scenario with two blocks
#     pinned to the same value must be found — a strict chaining inequality
#     in _assemble.choose would silently prune it. Same-value blocks merge
#     into one span constraint.  Planted truths + negative controls.
# ---------------------------------------------------------------------------

class TestSameValueBlocks(unittest.TestCase):

    def test_adjacent_orderstats_same_value(self):
        """Two adjacent order stats pinned to the SAME absent value: truth is
        two added copies (k=2), with caps_hit honest."""
        base = [F(i, 10) for i in range(1, 41) if i != 20]   # 2.0 absent
        v = F(2)
        full = sorted(base + [v, v])                         # n' = 41
        j = full.index(v) + 1                                # 1-based
        targets = [(OrderStat(j), v), (OrderStat(j + 1), v)]
        self.assertEqual(OrderStat(j).value(full), v)
        self.assertEqual(OrderStat(j + 1).value(full), v)
        res = solve(base, targets, register=F(1, 10), value_lo=0, value_hi=10,
                    k_cap=4, rem_cap=0, max_witnesses=16)
        self.assertEqual(res["k"], 2, res["caps_hit"])
        self.assertEqual(res["caps_hit"], [])
        self.assertTrue(any(
            sorted(s["value"] for s in f["adds"]) == [v, v]
            for f in res["families"]), res["families"])
        for fam in res["families"]:
            self.assertTrue(targets_exact(materialize(base, fam), targets))

    def test_equal_quantiles_adjacent_supports_sas3(self):
        """Two published sas3 quantiles with adjacent supports and EQUAL
        values on tie-heavy 2dp data (the CMS complaint-measure shape)."""
        lows = [F(i, 100) for i in range(10, 100, 10)]                # 9
        highs = [F(i, 100) for i in range(200, 200 + 29 * 5, 5)]     # 29
        base = sorted(lows + highs)                                  # n=38
        v = F(150, 100)                                              # absent
        full = sorted(base + [v, v])                                 # n'=40
        qa = Quantile(F(1, 4), "sas3")     # x_(10) at n'=40
        qb = Quantile(F(11, 40), "sas3")   # x_(11) at n'=40
        self.assertEqual(qa.value(full), v)
        self.assertEqual(qb.value(full), v)
        targets = [(qa, v), (qb, v)]
        res = solve(base, targets, register=F(1, 100), value_lo=0,
                    value_hi=100, k_cap=4, rem_cap=0, max_witnesses=16)
        self.assertEqual(res["k"], 2, res["caps_hit"])
        self.assertTrue(any(
            sorted(s["value"] for s in f["adds"]) == [v, v]
            for f in res["families"]), res["families"])

    def test_q1_equals_q3_k1(self):
        """The README's own Q1 == Q3 case on tie-heavy data, with a k=1
        witness: base [1,5,5,9] + one more 5 gives sas5 Q1 = Q3 = 5."""
        v = F(5)
        base = [F(1), v, v, F(9)]
        full = sorted(base + [v])          # n'=5: Q1 = x_2 = 5, Q3 = x_4 = 5
        q1, q3 = Quantile(F(1, 4), "sas5"), Quantile(F(3, 4), "sas5")
        self.assertEqual(q1.value(full), v)
        self.assertEqual(q3.value(full), v)
        targets = [(q1, v), (q3, v)]
        self.assertFalse(targets_exact(base, targets))
        res = solve(base, targets, register=F(1), value_lo=0, value_hi=10,
                    k_cap=4, rem_cap=0, max_witnesses=16)
        self.assertEqual(res["k"], 1, res["caps_hit"])
        self.assertTrue(any(
            [s["value"] for s in f["adds"]] == [v]
            for f in res["families"]), res["families"])

    def test_nonadjacent_equal_value_blocks_k1(self):
        """Non-adjacent point blocks pinned to the same value (gap indices
        forced equal by monotonicity), k=1 truth."""
        lows = [F(i, 100) for i in range(10, 100, 10)]                # 9
        v = F(150, 100)
        highs = [F(i, 100) for i in range(300, 300 + 28 * 7, 7)]     # 28
        base = sorted(lows + [v, v, v] + highs)                      # n=40
        full = sorted(base + [v])                                    # n'=41
        qa = Quantile(F(10, 41), "sas3")   # x_(10) = v at n'=41
        qb = Quantile(F(13, 41), "sas3")   # x_(13) = v at n'=41
        self.assertEqual(qa.value(full), v)
        self.assertEqual(qb.value(full), v)
        targets = [(qa, v), (qb, v)]
        self.assertFalse(targets_exact(base, targets))
        res = solve(base, targets, register=F(1, 100), value_lo=0,
                    value_hi=100, k_cap=3, rem_cap=0, max_witnesses=16)
        self.assertEqual(res["k"], 1, res["caps_hit"])
        self.assertTrue(any(
            [s["value"] for s in f["adds"]] == [v]
            for f in res["families"]), res["families"])

    def test_wavg_merged_with_point_same_value(self):
        """A two-index averaging block (sas5 Q1 at 4 | n') and a point block
        (sas3) pinned to the same value merge into one span: x_11 = x_12 =
        x_13 = v needs THREE added copies of the absent v (k=3)."""
        lows = [F(i, 100) for i in range(10, 110, 10)]                # 10
        highs = [F(i, 100) for i in range(300, 300 + 31 * 7, 7)]     # 31
        base = sorted(lows + highs)                                  # n=41
        v = F(2)                                                     # absent
        full = sorted(base + [v, v, v])                              # n'=44
        q1 = Quantile(F(1, 4), "sas5")     # avg(x_11, x_12) at n'=44
        qp = Quantile(F(13, 44), "sas3")   # x_(13) at n'=44
        self.assertEqual(q1.value(full), v)
        self.assertEqual(qp.value(full), v)
        targets = [(q1, v), (qp, v)]
        res = solve(base, targets, register=F(1, 100), value_lo=0,
                    value_hi=100, k_cap=4, rem_cap=0, max_witnesses=16)
        self.assertEqual(res["k"], 3, res["caps_hit"])
        self.assertTrue(any(
            sorted(s["value"] for s in f["adds"]) == [v, v, v]
            for f in res["families"]), res["families"])

    def test_negative_control_nonmonotone_pins(self):
        """Order stats pinned to DECREASING values are genuinely infeasible:
        must stay empty (the merge must not over-accept), with only the depth
        note in caps_hit (infeasible-order pruning stays silent)."""
        base = [F(1), F(5), F(5), F(9)]
        targets = [(OrderStat(2), F(5)), (OrderStat(3), F(4))]
        res = solve(base, targets, register=F(1), value_lo=0, value_hi=10,
                    k_cap=3, rem_cap=0)
        self.assertIsNone(res["k"])
        self.assertEqual(res["families"], [])
        self.assertEqual(res["caps_hit"], ["k_cap-exhausted"])

    def test_negative_control_span_unreachable(self):
        """Equal pinned values BELOW the allowed add range and absent from
        the base: the span constraint is unsatisfiable — empty result."""
        base = [F(6), F(7), F(8), F(9)]
        targets = [(OrderStat(2), F(5)), (OrderStat(3), F(5))]
        res = solve(base, targets, register=F(1), value_lo=6, value_hi=10,
                    k_cap=3, rem_cap=0)
        self.assertIsNone(res["k"])
        self.assertEqual(res["families"], [])


# ---------------------------------------------------------------------------
# 6. parity at solver level (averaging argument as constraints)
# ---------------------------------------------------------------------------

class TestSolverParity(unittest.TestCase):

    def test_averaging_straddle_from_existing(self):
        """n=55 base, +1 add -> n'=56 = 0 mod 4: Q1 must be the AVERAGE of
        order stats 14, 15.  The solver must derive the straddle constraint
        and recover the planted low add."""
        xs = tie_free_data(404, 55)
        planted = xs[0] - F(1, 100)   # below the minimum: pure index shift
        full = sorted(xs + [planted])
        self.assertEqual(len(full) % 4, 0)
        q1 = Quantile(F(1, 4), "sas5").value(full)
        q3 = Quantile(F(3, 4), "sas5").value(full)
        # averaging really is active and non-trivial:
        self.assertEqual(q1, (full[13] + full[14]) / 2)
        self.assertNotEqual(full[13], full[14])
        targets = [(Quantile(F(1, 4), "sas5"), q1),
                   (Quantile(F(3, 4), "sas5"), q3)]
        res = solve(xs, targets, register=F(1, 100), value_lo=0,
                    value_hi=100, k_cap=2, rem_cap=0, max_witnesses=64)
        self.assertEqual(res["k"], 1)
        self.assertTrue(any(family_covers(f, planted)
                            for f in res["families"]))

    def test_averaging_added_straddle_partner(self):
        """Plant the add INSIDE the averaged pair: the new Q1 averages an
        existing value with the added value; the solver must pin the added
        value exactly (crit)."""
        xs = tie_free_data(405, 55)
        # place the add just above order stat 14 of the base (wide gap check)
        lo_anchor = xs[13]
        planted = lo_anchor + F(1, 100)
        if planted in xs:
            planted += F(1, 100)
        full = sorted(xs + [planted])
        q1 = Quantile(F(1, 4), "sas5").value(full)
        q3 = Quantile(F(3, 4), "sas5").value(full)
        if q1 != (lo_anchor + planted) / 2:
            self.skipTest("layout shifted; plant not inside the pair")
        targets = [(Quantile(F(1, 4), "sas5"), q1),
                   (Quantile(F(3, 4), "sas5"), q3)]
        res = solve(xs, targets, register=F(1, 100), value_lo=0,
                    value_hi=100, k_cap=2, rem_cap=0, max_witnesses=64)
        self.assertEqual(res["k"], 1)
        pinned = [f for f in res["families"]
                  if any(s["role"] == "crit" and s["value"] == planted
                         for s in f["adds"])]
        self.assertTrue(pinned, res["families"])

    def test_both_equal_two_adds(self):
        """Required Q1 not present in the base at all: the averaging branch
        'both order stats equal v' needs TWO added copies -> minimal k=2."""
        xs = tie_free_data(406, 54)
        # v in the gap between x_13 and x_14 (1-based): with two copies
        # inserted there the merged order stats 14 and 15 are both v
        gap_lo, gap_hi = xs[12], xs[13]
        self.assertGreater(gap_hi - gap_lo, F(2, 100))
        v = gap_lo + F(1, 100)
        full = sorted(xs + [v, v])           # n'=56, Q1 = avg(y14, y15)
        q1 = Quantile(F(1, 4), "sas5").value(full)
        q3 = Quantile(F(3, 4), "sas5").value(full)
        self.assertEqual(q1, v)              # both averaged stats are v
        targets = [(Quantile(F(1, 4), "sas5"), q1),
                   (Quantile(F(3, 4), "sas5"), q3)]
        res = solve(xs, targets, register=F(1, 100), value_lo=0,
                    value_hi=100, k_cap=3, rem_cap=0, max_witnesses=64)
        self.assertIsNotNone(res["k"])
        if res["k"] == 2:
            self.assertTrue(any(
                sorted(s["value"] for s in f["adds"]) == [v, v] or
                (family_covers(f, v) and f["k_add"] == 2)
                for f in res["families"]), res["families"])
        else:
            # a smaller family exists; it must still verify exactly
            for fam in res["families"]:
                self.assertTrue(targets_exact(materialize(xs, fam), targets))


# ---------------------------------------------------------------------------
# 7. capped fences (health-like), min/max/range, mean, multi-target
# ---------------------------------------------------------------------------

class TestCappedFence(unittest.TestCase):

    def test_hi_only_recovery(self):
        """Synthetic complaints-like measure: lo fence capped at 0 (published
        '0'), hi fence published exactly; plant one high add."""
        rng = random.Random(77)
        xs = sorted(F(x, 100) for x in rng.sample(range(0, 300), 87))
        planted = F(285, 100)
        if planted in xs:
            planted = F(286, 100)
        full = sorted(xs + [planted])
        lo_stat = TukeyFence("lo", 3, "sas5", cap=0)
        hi_stat = TukeyFence("hi", 3, "sas5")
        self.assertLess(lo_stat.raw_value(full), 0)   # genuinely capped
        pub_lo = lo_stat.value(full)
        pub_hi = hi_stat.value(full)
        self.assertEqual(pub_lo, 0)
        targets = [(lo_stat, pub_lo), (hi_stat, pub_hi)]
        self.assertFalse(targets_exact(xs, targets))
        res = solve(xs, targets, register=F(1, 100), value_lo=0,
                    k_cap=2, rem_cap=0, max_witnesses=64)
        self.assertEqual(res["k"], 1, res["caps_hit"])
        self.assertTrue(any(family_covers(f, planted)
                            for f in res["families"]))


class TestMinMaxRangeMean(unittest.TestCase):

    def test_min_planted(self):
        xs = tie_free_data(505, 40, lo=1000, hi=9000)
        planted = xs[0] - F(7, 100)
        targets = [(MinStat(), planted), (MaxStat(), xs[-1])]
        res = solve(xs, targets, register=F(1, 100), value_lo=0,
                    value_hi=100, k_cap=2)
        self.assertEqual(res["k"], 1)
        self.assertTrue(any(family_covers(f, planted)
                            for f in res["families"]))

    def test_range_planted(self):
        xs = tie_free_data(506, 40, lo=1000, hi=9000)
        planted = xs[-1] + F(2, 100)
        pub_range = planted - xs[0]
        targets = [(RangeStat(), pub_range)]
        res = solve(xs, targets, register=F(1, 100), value_lo=0,
                    value_hi=100, k_cap=2, max_witnesses=32)
        self.assertEqual(res["k"], 1)
        self.assertTrue(any(family_covers(f, planted)
                            for f in res["families"]),
                        res["families"])

    def test_mean_only(self):
        xs = tie_free_data(507, 33)
        planted = F(4242, 100)
        full = sorted(xs + [planted])
        targets = [(Mean(), Mean().value(full))]
        res = solve(xs, targets, register=F(1, 100), value_lo=0,
                    value_hi=100, k_cap=2)
        self.assertEqual(res["k"], 1)
        vals = [s["value"] for s in res["families"][0]["adds"]]
        self.assertEqual(vals, [planted])   # sum constraint pins it exactly

    def test_mean_plus_median(self):
        """Mean pins the free filler inside its order-stat window."""
        xs = tie_free_data(508, 40)
        planted = xs[30] + F(1, 100)  # above the median, in a gap
        if planted in xs:
            planted += F(1, 100)
        full = sorted(xs + [planted])
        targets = [(Quantile(F(1, 2), "sas5"), Quantile(F(1, 2), "sas5").value(full)),
                   (Mean(), Mean().value(full))]
        res = solve(xs, targets, register=F(1, 100), value_lo=0,
                    value_hi=100, k_cap=2, max_witnesses=32)
        self.assertEqual(res["k"], 1)
        hit = [f for f in res["families"]
               if sorted(s["value"] for s in f["adds"]) == [planted]]
        self.assertTrue(hit, res["families"])


class TestMultiTargetTriple(unittest.TestCase):

    def test_q1_median_q3_two_adds(self):
        """Three simultaneous pinned targets (Q1 + median + Q3) with two
        planted adds: exercises 3-block assembly."""
        xs = tie_free_data(606, 41)
        p_lo = xs[5] + F(1, 100)
        p_hi = xs[35] + F(1, 100)
        for p in (p_lo, p_hi):
            self.assertNotIn(p, xs)
        full = sorted(xs + [p_lo, p_hi])
        targets = [
            (Quantile(F(1, 4), "sas5"), Quantile(F(1, 4), "sas5").value(full)),
            (Quantile(F(1, 2), "sas5"), Quantile(F(1, 2), "sas5").value(full)),
            (Quantile(F(3, 4), "sas5"), Quantile(F(3, 4), "sas5").value(full)),
        ]
        self.assertFalse(targets_exact(xs, targets))
        res = solve(xs, targets, register=F(1, 100), value_lo=0,
                    value_hi=100, k_cap=3, rem_cap=0, max_witnesses=64)
        self.assertIsNotNone(res["k"])
        self.assertLessEqual(res["k"], 2)
        if res["k"] == 2:
            self.assertTrue(any(family_covers(f, p_lo) and
                                family_covers(f, p_hi)
                                for f in res["families"]),
                            res["families"])


class TestIntervalTargets(unittest.TestCase):

    def test_interval_probe_path(self):
        """Interval-only targets exercise the probe fallback: bring Q3
        inside a published interval."""
        xs = tie_free_data(707, 41)
        q3 = Quantile(F(3, 4), "sas5")
        cur = q3.value(xs)
        iv = Interval(cur + F(1, 100), cur + F(200, 100))
        targets = [(q3, iv)]
        self.assertFalse(targets_exact(xs, targets))
        res = solve(xs, targets, register=F(1, 100), value_lo=0,
                    value_hi=100, k_cap=2, max_witnesses=8)
        self.assertIsNotNone(res["k"], res["caps_hit"])
        self.assertIn("probe-fallback-used", res["caps_hit"])
        for fam in res["families"]:
            self.assertTrue(targets_exact(materialize(xs, fam), targets))


class TestWindowReceipt(unittest.TestCase):

    def test_filler_window_is_real(self):
        """The window claim is executable: any register point inside a
        filler window keeps every order-stat target exact."""
        xs = tie_free_data(808, 57)
        planted = xs[-1] - F(3, 100)   # high region
        if planted in xs:
            planted -= F(1, 100)
        full = sorted(xs + [planted])
        targets = [(Quantile(F(1, 4), "sas5"), Quantile(F(1, 4), "sas5").value(full)),
                   (Quantile(F(3, 4), "sas5"), Quantile(F(3, 4), "sas5").value(full))]
        res = solve(xs, targets, register=F(1, 100), value_lo=0,
                    value_hi=100, k_cap=2, max_witnesses=64)
        self.assertEqual(res["k"], 1)
        windowed = [(f, s) for f in res["families"] for s in f["adds"]
                    if s["window"] is not None]
        self.assertTrue(windowed, "no filler-window family returned")
        fam, spec = windowed[0]
        # pick a different register point strictly inside the window
        alt = None
        for step in range(1, 200):
            for cand in (spec["value"] + F(step, 100),
                         spec["value"] - F(step, 100)):
                if cand != spec["value"] and window_contains(spec["window"],
                                                             cand):
                    alt = cand
                    break
            if alt is not None:
                break
        self.assertIsNotNone(alt, spec)
        # swap exactly one copy of the windowed filler for the alternate
        swapped = [dict(s) for s in fam["adds"]]
        for s in swapped:
            if s["value"] == spec["value"] and s["window"] == spec["window"]:
                s["value"] = alt
                break
        after = sorted(list(xs) + [s["value"] for s in swapped])
        self.assertTrue(targets_exact(after, targets), (spec, alt))


class TestLargeK(unittest.TestCase):

    def test_true_minimal_k_30(self):
        """30 planted adds all below the minimum force a genuinely large
        minimal family (index-shift arithmetic: fewer adds cannot reproduce
        the quartile pair).  The solver must reach k=30 within its caps."""
        rng = random.Random(43)
        xs = sorted(F(x, 100) for x in rng.sample(range(2000, 9500), 301))
        plant = [F(x, 100) for x in rng.sample(range(100, 1900), 30)]
        full = sorted(xs + plant)
        q1, q3 = Quantile(F(1, 4), "sas5"), Quantile(F(3, 4), "sas5")
        targets = [(q1, q1.value(full)), (q3, q3.value(full))]
        res = solve(xs, targets, register=F(1, 100), value_lo=0,
                    value_hi=100, k_cap=30, rem_cap=0, max_witnesses=3)
        self.assertEqual(res["k"], 30, res["caps_hit"])
        for fam in res["families"]:
            self.assertTrue(targets_exact(materialize(xs, fam), targets))


# ---------------------------------------------------------------------------
# 7b. OrderStat index beyond the current array size: evaluates not-exact,
#     never raises out of solve() — an index exceeding n' at some small k
#     must not kill a stack() run when larger k makes the instance solvable.
# ---------------------------------------------------------------------------

class TestOrderStatOutOfRange(unittest.TestCase):

    def test_solves_beyond_current_n(self):
        """OrderStat(43) on n=40 needs k_add >= 3; the k loop must reach it
        (k=1,2 are infeasible SIZES, not errors)."""
        base = [F(i, 10) for i in range(10, 50)]        # n=40
        targets = [(OrderStat(43), F(9))]
        res = solve(base, targets, register=F(1, 10), value_lo=0,
                    value_hi=10, k_cap=5, rem_cap=0, max_witnesses=8)
        self.assertEqual(res["k"], 3, res["caps_hit"])
        for fam in res["families"]:
            after = materialize(base, fam)
            self.assertEqual(len(after), 43)
            self.assertEqual(after[42], F(9))
            self.assertTrue(targets_exact(after, targets))

    def test_negative_index_beyond_current_n(self):
        """Negative (from-the-top) index out of range at small k: same
        infeasible-size semantics.  OrderStat(-43) is x'_(1) once n'=43."""
        base = [F(i, 10) for i in range(10, 50)]        # n=40, min 1.0
        v = F(1, 2)
        targets = [(OrderStat(-43), v)]
        res = solve(base, targets, register=F(1, 10), value_lo=0,
                    value_hi=10, k_cap=5, rem_cap=0, max_witnesses=8)
        self.assertEqual(res["k"], 3, res["caps_hit"])
        for fam in res["families"]:
            after = materialize(base, fam)
            self.assertEqual(after[0], v)

    def test_unreachable_within_k_cap_is_loud_not_crash(self):
        """k_cap below the needed size: empty result with the depth note —
        never a ValueError."""
        base = [F(i, 10) for i in range(10, 50)]
        targets = [(OrderStat(43), F(9))]
        res = solve(base, targets, register=F(1, 10), value_lo=0,
                    value_hi=10, k_cap=2, rem_cap=0)
        self.assertIsNone(res["k"])
        self.assertEqual(res["families"], [])
        self.assertIn("k_cap-exhausted", res["caps_hit"])

    def test_evaluate_out_of_range_is_not_exact(self):
        """targets_exact on an array too small for the OrderStat: False,
        not a crash (evaluation semantics: undefined statistic matches no
        published value)."""
        base = [F(i, 10) for i in range(10, 50)]
        self.assertFalse(targets_exact(base, [(OrderStat(43), F(9))]))
        self.assertFalse(targets_exact(base, [(OrderStat(-43), F(1))]))

    def test_out_of_range_with_removals_no_crash(self):
        """rem_cap > 0 exercises _removal_candidates on the ORIGINAL n; the
        out-of-range component must be skipped there too."""
        base = [F(i, 10) for i in range(10, 50)]
        targets = [(OrderStat(43), F(9))]
        res = solve(base, targets, register=F(1, 10), value_lo=0,
                    value_hi=10, k_cap=4, rem_cap=1, max_witnesses=4)
        self.assertEqual(res["k"], 3, res["caps_hit"])

    def test_stack_survives_out_of_range_orderstat(self):
        """In stack(), one dataset whose OrderStat target exceeds its row
        count must not kill the run: it is a miss, solved at the k that
        makes the index valid; the exact companion dataset stays exact."""
        rng = random.Random(1301)
        rows = [(f"S{i:02d}", F(x, 100)) for i, x in
                enumerate(rng.sample(range(100, 900), 20))]
        vals = sorted(v for _, v in rows)
        v_new = vals[-1] + F(50, 100)
        ds = Dataset("OOR", rows, [(OrderStat(22), v_new)],
                     register=F(1, 100), value_lo=0, value_hi=100)
        good_rows = [(f"G{i:02d}", F(x, 100)) for i, x in
                     enumerate(rng.sample(range(100, 900), 15))]
        good_vals = sorted(v for _, v in good_rows)
        good = Dataset("GOOD", good_rows,
                       [(Quantile(F(1, 2), "sas5"),
                         Quantile(F(1, 2), "sas5").value(good_vals))],
                       register=F(1, 100))
        res = stack([ds, good], k_cap=3, rem_cap=0, workers=0)
        self.assertEqual(res["unsolved"], [], res["notes"])
        self.assertEqual(res["exact_after"], 2, res["battery"])
        self.assertEqual(res["regressions"], [])


# ---------------------------------------------------------------------------
# 7c. loudness receipts: cap-exhausted searches emit a witness, never
#     silence — caps_hit is the certification currency: every truncation and
#     every structural skip must leave a note.
# ---------------------------------------------------------------------------

class TestLoudness(unittest.TestCase):

    def test_k_cap_exhausted_note(self):
        """Planted k=3 truth searched at k_cap=2: empty result MUST carry
        k_cap-exhausted, never silence; at k_cap=3 the witness is
        found and the note is absent."""
        base = [F(i, 10) for i in range(20, 60)]
        plant = [F(1, 10), F(2, 10), F(3, 10)]
        full = sorted(base + plant)
        q1, q3 = Quantile(F(1, 4), "sas5"), Quantile(F(3, 4), "sas5")
        targets = [(q1, q1.value(full)), (q3, q3.value(full))]
        res2 = solve(base, targets, register=F(1, 10), value_lo=0,
                     value_hi=10, k_cap=2, rem_cap=0)
        self.assertIsNone(res2["k"])
        self.assertIn("k_cap-exhausted", res2["caps_hit"])
        res3 = solve(base, targets, register=F(1, 10), value_lo=0,
                     value_hi=10, k_cap=3, rem_cap=0)
        self.assertEqual(res3["k"], 3)
        self.assertNotIn("k_cap-exhausted", res3["caps_hit"])

    def test_exact_already_no_depth_note(self):
        xs = tie_free_data(1401, 41)
        q = Quantile(F(1, 2), "sas5")
        res = solve(xs, [(q, q.value(xs))], register=F(1, 100))
        self.assertEqual(res["k"], 0)
        self.assertTrue(res["exact_already"])
        self.assertNotIn("k_cap-exhausted", res["caps_hit"])

    def test_register_none_fence_grid_is_loud(self):
        """register=None one-sided fence modes grid the free quantile on
        data values only — a structural limit that must leave a cap note,
        never pass silently."""
        base = [F(i, 7) for i in range(10, 51)]          # off-lattice, n=41
        lo = TukeyFence("lo", 3, "sas5", cap=0)
        hi = TukeyFence("hi", 3, "sas5")
        plant = F(40, 7) + F(1, 113)
        full = sorted(base + [plant])
        pub_lo, pub_hi = lo.value(full), hi.value(full)
        self.assertEqual(pub_lo, 0)                      # capped -> hi_only
        targets = [(lo, pub_lo), (hi, pub_hi)]
        self.assertFalse(targets_exact(base, targets))
        res = solve(base, targets, register=None, value_lo=0, value_hi=100,
                    k_cap=2, rem_cap=0, max_witnesses=8)
        self.assertIn("register-none-grid-data-values-only", res["caps_hit"])

    def test_touching_straddle_is_loud(self):
        """Residual structural limit, now loud: a straddle end touching a
        NEIGHBORING pinned value (different pinned values sharing a boundary
        value) is skipped with touching-realizations-skipped.  The planted
        k=1 witness [2] is missed (the solver returns the loud k=2), so the
        note is what keeps the miss honest."""
        base = [F(1), F(1), F(2), F(8), F(9), F(10), F(11)]     # n=7
        targets = [(OrderStat(3), F(2)),
                   (Quantile(F(1, 2), "sas5"), F(5))]
        # planted truth: one added 2 -> n'=8, x_3=2, median=(2+8)/2=5
        self.assertTrue(targets_exact(sorted(base + [F(2)]), targets))
        res = solve(base, targets, register=F(1), value_lo=0, value_hi=20,
                    k_cap=2, rem_cap=0, max_witnesses=8)
        self.assertIn("touching-realizations-skipped", res["caps_hit"])
        # any families returned must still verify exactly
        for fam in res["families"]:
            self.assertTrue(targets_exact(materialize(base, fam), targets))


# ---------------------------------------------------------------------------
# 8. spawn pool
# ---------------------------------------------------------------------------

class TestSpawnPool(unittest.TestCase):

    def test_solve_many_spawn_matches_serial(self):
        jobs = []
        expected = []
        for seed in (901, 902):
            xs = tie_free_data(seed, 45)
            planted = xs[40] + F(1, 100)
            full = sorted(xs + [planted])
            targets = [(Quantile(F(1, 4), "sas5"),
                        Quantile(F(1, 4), "sas5").value(full)),
                       (Quantile(F(3, 4), "sas5"),
                        Quantile(F(3, 4), "sas5").value(full))]
            jobs.append({"job_id": str(seed), "values": xs,
                         "targets": targets, "register": F(1, 100),
                         "value_lo": F(0), "value_hi": F(100),
                         "k_cap": 2, "rem_cap": 0, "max_witnesses": 4})
            expected.append(solve(xs, targets, register=F(1, 100),
                                  value_lo=0, value_hi=100, k_cap=2,
                                  max_witnesses=4))
        results = solve_many(jobs, workers=2)   # spawn context
        by_id = {r["job_id"]: r for r in results}
        for seed, exp in zip(("901", "902"), expected):
            self.assertEqual(by_id[seed]["k"], exp["k"])
            self.assertEqual(
                sorted(tuple(sorted(s["value"] for s in f["adds"]))
                       for f in by_id[seed]["families"]),
                sorted(tuple(sorted(s["value"] for s in f["adds"]))
                       for f in exp["families"]))


# ---------------------------------------------------------------------------
# 9. stacking end-to-end
# ---------------------------------------------------------------------------

class TestStackerE2E(unittest.TestCase):

    def _mkrows(self, seed, n, prefix):
        rng = random.Random(seed)
        vals = [F(x, 100) for x in rng.sample(range(500, 9500), n)]
        return [(f"{prefix}{i:03d}", v) for i, v in enumerate(vals)]

    def test_three_datasets_two_shared_entities(self):
        """Plant 2 shared entities across 3 datasets (D1 super / D2 sub
        value-coupled; D3 independent), publish targets from the full data,
        hand the stacker the reduced rows.  A 4th dataset is exact and must
        stay exact (no-regression gate)."""
        rows1 = self._mkrows(1001, 53, "A")
        rows2 = self._mkrows(1002, 45, "B")
        rows3 = self._mkrows(1003, 61, "C")
        rows4 = self._mkrows(1004, 41, "D")
        rng = random.Random(1005)

        def fresh(vals_rows, lo, hi):
            have = {v for _, v in vals_rows}
            while True:
                c = F(rng.randint(lo, hi), 100)
                if c not in have:
                    return c

        # entity E1: low values; E2: high values.  D1/D2 coupled => equal.
        e1_d12 = fresh(rows1 + rows2, 600, 2500)
        e2_d12 = fresh(rows1 + rows2, 7500, 9400)
        e1_d3 = fresh(rows3, 600, 2500)
        e2_d3 = fresh(rows3, 7500, 9400)

        def q_targets(rows, extra):
            vals = sorted([v for _, v in rows] + extra)
            return [(Quantile(F(1, 4), "sas5"),
                     Quantile(F(1, 4), "sas5").value(vals)),
                    (Quantile(F(3, 4), "sas5"),
                     Quantile(F(3, 4), "sas5").value(vals))]

        d1 = Dataset("D1", rows1, q_targets(rows1, [e1_d12, e2_d12]),
                     register=F(1, 100), value_lo=0, value_hi=100)
        d2 = Dataset("D2", rows2, q_targets(rows2, [e1_d12, e2_d12]),
                     register=F(1, 100), value_lo=0, value_hi=100)
        d3 = Dataset("D3", rows3, q_targets(rows3, [e1_d3, e2_d3]),
                     register=F(1, 100), value_lo=0, value_hi=100)
        d4 = Dataset("D4", rows4, q_targets(rows4, []),
                     register=F(1, 100), value_lo=0, value_hi=100)
        for d, miss in ((d1, True), (d2, True), (d3, True), (d4, False)):
            self.assertEqual(not d.exact(), miss, d.name)

        classes = ("SHARED", "SOLO")
        # SOLO entities may not participate in the sub dataset D2
        participation = lambda cls, name: not (cls == "SOLO" and name == "D2")
        res = stack([d1, d2, d3, d4], classes=classes,
                    participation=participation, couplings=[("D1", "D2")],
                    k_cap=3, rem_cap=0, max_witnesses=16, workers=0)
        self.assertEqual(res["unsolved"], [], res["notes"])
        self.assertEqual(res["exact_after"], 4, res["battery"])
        self.assertEqual(res["regressions"], [])
        self.assertLessEqual(len(res["entities"]), 2, res["entities"])
        # coupling: any entity with a D2 value carries the SAME D1 value
        for e in res["entities"]:
            if "D2" in e["values"]:
                self.assertIn("D1", e["values"])
                self.assertEqual(e["values"]["D1"], e["values"]["D2"])
        # planted delta itself is (independently) a valid solution
        self.assertTrue(d1.exact(adds=[e1_d12, e2_d12]))
        self.assertTrue(d2.exact(adds=[e1_d12, e2_d12]))
        self.assertTrue(d3.exact(adds=[e1_d3, e2_d3]))

    def test_unreachable_dataset_reported_unsolved(self):
        rows = self._mkrows(1101, 41, "U")
        vals = sorted(v for _, v in rows)
        bad = Dataset("BAD", rows,
                      [(Quantile(F(1, 4), "sas5"), F(1)),
                       (Quantile(F(3, 4), "sas5"), F(9999))],
                      register=F(1, 100), value_lo=0, value_hi=100)
        good_rows = self._mkrows(1102, 37, "G")
        good = Dataset("GOOD", good_rows,
                       [(Quantile(F(1, 2), "sas5"),
                         Quantile(F(1, 2), "sas5").value(
                             sorted(v for _, v in good_rows)))],
                       register=F(1, 100))
        res = stack([bad, good], k_cap=2, rem_cap=1, workers=0)
        self.assertIn("BAD", res["unsolved"])
        self.assertEqual(res["regressions"], [])
        self.assertEqual(
            [r["after_exact"] for r in res["battery"] if r["name"] == "GOOD"],
            [True])

    def test_removal_gate_blocks_shared_entity(self):
        """Dataset A is fixable ONLY by removing its minimum entity, but that
        entity also anchors exact dataset B's published minimum: the global
        no-regression gate must refuse, leaving A unsolved and B exact."""
        rng = random.Random(1201)
        shared_min = F(111, 100)
        rows_a = [("X000", shared_min)] + \
            [(f"A{i:03d}", F(x, 100))
             for i, x in enumerate(rng.sample(range(300, 8000), 40))]
        rows_b = [("X000", shared_min)] + \
            [(f"B{i:03d}", F(x, 100))
             for i, x in enumerate(rng.sample(range(300, 8000), 35))]
        vals_a = sorted(v for _, v in rows_a)
        second_min = vals_a[1]
        a = Dataset("A", rows_a, [(MinStat(), second_min)],
                    register=F(1, 100), value_lo=second_min, value_hi=100)
        b = Dataset("B", rows_b, [(MinStat(), shared_min)],
                    register=F(1, 100), value_lo=0, value_hi=100)
        self.assertFalse(a.exact())
        self.assertTrue(b.exact())
        res = stack([a, b], k_cap=2, rem_cap=1, workers=0)
        self.assertIn("A", res["unsolved"])
        self.assertEqual(res["regressions"], [])
        self.assertEqual(res["removed_entities"], [])
        self.assertEqual(
            [r["after_exact"] for r in res["battery"] if r["name"] == "B"],
            [True])


# ---------------------------------------------------------------------------
# 10. REGRESSION: the originating analysis's 2024 D02/MA-PD +1-addition family
# ---------------------------------------------------------------------------

class TestOriginRegression(unittest.TestCase):
    """Reproduce, read-only from the reference tree (QINVERT_ORIGIN_ROOT), the
    result recorded in its membership_families.json: star-year 2024, measure
    D02 (Complaints about the Drug Plan), MA-PD org split, hypothesis
    cost+d60r+emp, scores vintage 'original': n=534 contracts, published
    outer fences (0, 1.41) with the low side capped at 0, and a UNIQUE
    minimal profile: ONE addition in the high region with witness value
    21/50 = 0.42 (receipt fields asserted below)."""

    FAMS = os.path.join(ORIGIN, "runs", "membership_families.json")
    SCORES = os.path.join(ORIGIN, "data", "parsed",
                          "scores_2024_original.csv")
    SUMMARY = os.path.join(ORIGIN, "data", "parsed",
                           "summary_2024_original.csv")

    def _lane_receipt(self):
        with open(self.FAMS) as f:
            fams = json.load(f)
        for r in fams:
            if (r.get("year"), r.get("mid"), r.get("org")) == \
                    ("2024", "D02", "MA-PD"):
                return r
        return None

    def _rebuild_values(self):
        """Replays load_year semantics from the origin membership_stacker.py
        (read-only): MA-PD = org_type without 'PDP'; hypothesis
        cost+d60r+emp with the 2024 disaster rule active."""
        meta = {}
        with open(self.SUMMARY, newline="", encoding="utf-8-sig") as f:
            for r in csv.DictReader(f):
                vals = []
                for i in (1, 2):
                    v = (r.get(f"disaster_col{i}_raw") or "").strip()
                    try:
                        vals.append(int(v))
                    except ValueError:
                        vals.append(None)
                recent = vals[1] if vals[1] is not None else 0
                meta[r["contract_id"]] = (r.get("org_type", ""), recent)
        kept = []
        with open(self.SCORES, newline="", encoding="utf-8-sig") as f:
            for r in csv.DictReader(f):
                if r["measure_id"] != "D02" or not r["score"].strip():
                    continue
                org, recent = meta.get(r["contract_id"], ("", 0))
                if "PDP" in org or "1876" in org or "Employer" in org:
                    continue
                if recent >= 60:
                    continue
                kept.append(F(r["score"].strip()))
        return sorted(kept)

    def test_d02_2024_one_addition_family(self):
        for p in (self.FAMS, self.SCORES, self.SUMMARY):
            if not os.path.exists(p):
                self.skipTest(f"origin-tree file missing: {p}")
        # (a) the recorded origin receipt says: 1 profile, k_add=1, witness 21/50
        rec = self._lane_receipt()
        self.assertIsNotNone(rec, "2024 D02 MA-PD not in origin receipt")
        self.assertEqual(rec["n"], 534)
        self.assertEqual(rec["pub"], ["0", "1.41"])
        self.assertEqual(rec["n_profiles"], 1)
        prof = rec["profiles"][0]
        self.assertEqual((prof["k_add"], prof["k_rem"]), (1, 0))
        self.assertIn("21/50", prof["constraints"])
        # (b) rebuild the dataset read-only and reproduce with qinvert
        xs = self._rebuild_values()
        self.assertEqual(len(xs), 534)
        lo_stat = TukeyFence("lo", 3, "sas5", cap=0)
        hi_stat = TukeyFence("hi", 3, "sas5")
        targets = [(lo_stat, F(0)), (hi_stat, F("1.41"))]
        self.assertFalse(targets_exact(xs, targets))   # it is a miss
        res = solve(xs, targets, register=F(1, 100), value_lo=0,
                    k_cap=2, rem_cap=0, max_witnesses=12)
        self.assertEqual(res["k"], 1, res["caps_hit"])
        add_only = [f for f in res["families"] if f["k_rem"] == 0]
        self.assertTrue(add_only)
        self.assertTrue(any(family_covers(f, F(21, 50)) for f in add_only),
                        [(f["adds"]) for f in add_only])
        # (c) double gate through the origin tree's own fence code (read-only)
        tukey_fences = _import_origin_symbol("tukey", "tukey_fences")
        lo, hi = tukey_fences(sorted(xs + [F(21, 50)]), cap_lo="0")
        self.assertEqual((lo, hi), (F(0), F("1.41")))


if __name__ == "__main__":
    unittest.main(verbosity=2)

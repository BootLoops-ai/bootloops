"""Regression: the lockpick lattice fitter's synthetic house control (3f0-7f1)/5
carries the SPAN TEST.  On a linearly dependent member set the lattice returns
a SHORTER representative of the same function; the control's entry then reads
ok = True with reason = equivalent_representative when got - want lies in the
span of the member set's basis-internal relations (the exact rank test
rank(R) == rank(R + [got - want]) on the fit's own basis_relations — integer
elimination, no floating rank), ok = False with reason = not_in_span
otherwise, and reason = exact on a match.

The regression guarded: a 281-member set with 85 basis-internal relations (at
most 196 independent functions) once printed ok = False on both lattice routes
(flint and fplll proved mode) for a representative one unit shorter than the
planted vector, while got - want lay in the span of the relations (rank
85 -> 85 with the difference adjoined).  The two records' house_controls
(want / got) and main.basis_relations are vendored beside this file as slim
extracts pinned to the records' sha256 (mplll_span_fixtures/); the verbatim
records are cross-checked when MPLLL_SPAN_FIXTURES names the directory that
holds them (skipped by name otherwise).

Legs: an INDEPENDENT synthetic set (ok = True, reason exact) and a DEPENDENT
synthetic set F = [f0, f1, f2 = f0 - f1, f3] — the target (3f0-7f1)/5 equals
(3f2-4f1)/5, whose relation vector (5, 0, 4, -3, 0) has norm^2 50 against the
planted (5, -3, 7, 0, 0)'s 83, so the reduced lattice returns it: got != want
and ok = True with reason equivalent_representative — on every LLL backend
present (pure always; fplll-cli when the binary is on PATH, skipped by name
otherwise); a PLANTED non-span difference (got perturbed off the span ->
ok = False, not_in_span) and a planted in-span representative with a
different c0 scaling; the two item-75 fixtures through house_verdict
(ok = True, equivalent_representative, rank 85 -> 85, heights 6 vs 7); the
exact rank on hand-built matrices, one of which a floating rank gets wrong; a
source guard that controls() fills its entry from house_verdict.
"""
import hashlib
import inspect
import json
import os
import shutil

import mpmath as mp
import pytest

from lockpick import mplll as ML
from lockpick.pslq_gate import canonicalize

D = 40                               # fit dps of the synthetic legs (seconds-class)
FIT, HO = [0, 1, 2, 3, 4], [5, 6]    # 7 points: 5 fit / 2 withheld
W = 10 ** 8
HERE = os.path.dirname(os.path.abspath(__file__))
FIXDIR = os.path.join(HERE, "mplll_span_fixtures")
FIXTURES = {  # slim extract -> (verbatim source file, its sha256)
    "item75_flint.json": ("FIT_item75_dps100_flint_20260906T152846Z.json",
                          "08ca4e65bd6ec7e088d4e84ea0ea7b6e41433675aba89569fbdabbb4143af91e"),
    "item75_fplll.json": ("FIT_item75_dps100_fplll_20260906T153013Z.json",
                          "8529cb86510a72091fbf1d7847e3c40b9ffb5b1e929735e18c9a2e686b82c87b"),
}
BACKENDS = ["pure",
            pytest.param("fplll-cli", marks=pytest.mark.skipif(
                not shutil.which("fplll"), reason="fplll binary not on PATH"))]
PLANTED_RELATION = (1, -1, -1, 0)    # f0 - f1 - f2 = 0 in the dependent set


def _build(dependent):
    """(D+30)-digit value strings for the member set: f0 = pi^x, f1 = log(2+x),
    [f2 = f0 - f1,] f3 = exp(x) on seven rational points."""
    with mp.workdps(D + 60):
        xs = [mp.mpf(3) / 10 + mp.mpf(j) / 7 for j in range(7)]
        f0 = [mp.power(mp.pi, x) for x in xs]
        f1 = [mp.log(2 + x) for x in xs]
        f3 = [mp.exp(x) for x in xs]
        F = [f0, f1, [a - b for a, b in zip(f0, f1)], f3] if dependent else [f0, f1, f3]
        s = lambda v: mp.nstr(v, D + 30, strip_zeros=False)
        return [[s(v) for v in r] for r in F]


def _controls(F, backend):
    names = ["f%d" % i for i in range(len(F))]
    return ML.controls(["0"] * len(F[0]), F, names, FIT, HO, D, W, backend=backend, method="rawlll")


def _synthetic_fit(F, backend):
    with mp.workdps(D + 30):
        Isyn = [mp.nstr((3 * ML._mpc(F[0][j]) - 7 * ML._mpc(F[1][j])) / 5, D + 20)
                for j in range(len(F[0]))]
    names = ["f%d" % i for i in range(len(F))]
    return ML.mplll_fit(Isyn, F, names, FIT, HO, D, W, two_prec=False, backend=backend, method="rawlll")


@pytest.mark.parametrize("backend", BACKENDS)
def test_independent_set_is_exact(backend):
    c = _controls(_build(False), backend)
    syn = c["detail"][0]
    assert syn["control"] == "synthetic (3f0-7f1)/5"
    assert syn["got"] == syn["want"] == canonicalize([5, -3, 7, 0])
    assert syn["reason"] == "exact"
    assert syn["ok"] is True and c["ok"] is True
    assert syn["relations_n"] == 0


@pytest.mark.parametrize("backend", BACKENDS)
def test_dependent_set_shorter_representative_is_equivalent(backend):
    c = _controls(_build(True), backend)
    syn = c["detail"][0]
    assert syn["want"] == canonicalize([5, -3, 7, 0, 0])
    assert syn["got"] is not None and syn["got"] != syn["want"], \
        "the reduced lattice must return the shorter representative, not the planted vector"
    assert syn["reason"] == "equivalent_representative"
    assert syn["ok"] is True and c["ok"] is True
    assert (syn["relations_n"], syn["relations_rank"], syn["rank_with_difference"]) == (1, 1, 1)
    assert syn["height_got"] < syn["height_want"] == 7
    assert canonicalize(syn["difference"]) == PLANTED_RELATION
    neg = c["detail"][1]
    assert neg["control"] == "negative (random target)" and neg["ok"] is True


@pytest.mark.parametrize("backend", BACKENDS)
def test_planted_non_span_difference_is_not_in_span(backend):
    r = _synthetic_fit(_build(True), backend)
    rels = r["basis_relations"]
    assert [canonicalize(b) for b in rels] == [PLANTED_RELATION]
    want = canonicalize([5, -3, 7, 0, 0])
    got = canonicalize(r["relation"])
    assert got != want
    # perturb got by a unit vector on f3, which no relation touches -> off the span
    planted = list(got)
    planted[-1] += 1
    v = ML.house_verdict(canonicalize(planted), want, rels)
    assert v["ok"] is False and v["reason"] == "not_in_span"
    assert (v["relations_n"], v["relations_rank"], v["rank_with_difference"]) == (1, 1, 2)
    # the unperturbed representative stays equivalent under the same relations
    assert ML.house_verdict(got, want, rels)["reason"] == "equivalent_representative"


def test_house_verdict_forms():
    want = canonicalize([5, -3, 7, 0, 0])
    rels = [list(PLANTED_RELATION)]
    assert ML.house_verdict(None, want, rels) == {"ok": False, "reason": "null"}
    assert ML.house_verdict(want, want, rels) == {"ok": True, "reason": "exact"}
    assert ML.house_verdict(list(want), want, []) == {"ok": True, "reason": "exact"}
    # a representative with c0 = 10 (2*want + one relation, gcd 1): the difference
    # want0*got - got0*want = (0, 5, -5, -5, 0) lies in the span
    got10 = canonicalize([10, -5, 13, -1, 0])
    v = ML.house_verdict(got10, want, rels)
    assert got10[0] == 10 and v["reason"] == "equivalent_representative" and v["ok"] is True
    assert canonicalize(v["difference"]) == PLANTED_RELATION
    assert (v["height_got"], v["height_want"]) == (13, 7)
    # no relations at all: any difference is off the span
    v0 = ML.house_verdict(canonicalize([5, -3, 7, 1, 0]), want, [])
    assert v0["reason"] == "not_in_span" and v0["ok"] is False
    assert (v0["relations_n"], v0["relations_rank"], v0["rank_with_difference"]) == (0, 0, 1)


def test_exact_rank_on_hand_built_matrices():
    assert ML._rank_exact([]) == 0
    assert ML._rank_exact([[0, 0, 0]]) == 0
    assert ML._rank_exact([[1, 2, 3], [2, 4, 6]]) == 1
    assert ML._rank_exact([[1, 2, 3], [2, 4, 6], [0, 0, 1]]) == 2
    assert ML._rank_exact([[1, 0], [0, 1], [3, 5]]) == 2
    assert ML._rank_exact([[2, 4], [3, 6], [5, 10]]) == 1
    big = 10 ** 40                    # a floating rank sees two equal rows here
    assert ML._rank_exact([[big, 1], [big + 1, 1]]) == 2


def _load_fixture(name):
    src, sha = FIXTURES[name]
    d = json.load(open(os.path.join(FIXDIR, name)))
    assert (d["source_file"], d["source_sha256"]) == (src, sha), "fixture pin"
    return d


@pytest.mark.parametrize("name", sorted(FIXTURES))
def test_item75_fixture_is_an_equivalent_representative(name):
    d = _load_fixture(name)
    want, got, R = d["want"], d["got"], d["basis_relations"]
    assert d["K"] == 281 and len(want) == len(got) == 282
    assert len(R) == 85 and all(len(r) == 281 for r in R)
    assert d["recorded_ok"] is False and got != want          # the verdict recorded before the span test
    v = ML.house_verdict(got, want, R)
    assert v["ok"] is True and v["reason"] == "equivalent_representative"
    assert (v["relations_n"], v["relations_rank"], v["rank_with_difference"]) == (85, 85, 85)
    assert (v["height_got"], v["height_want"]) == (6, 7)
    assert sum(1 for x in v["difference"] if x) == 3
    # off the span by one unit on a member the relations do not reach: not_in_span
    planted = list(got)
    planted[-1] += 1
    v2 = ML.house_verdict(canonicalize(planted), want, R)
    assert v2["reason"] == "not_in_span" and v2["rank_with_difference"] == 86


@pytest.mark.parametrize("name", sorted(FIXTURES))
def test_item75_extract_matches_the_verbatim_object(name):
    root = os.environ.get("MPLLL_SPAN_FIXTURES")
    if not root:
        pytest.skip("MPLLL_SPAN_FIXTURES unset: the verbatim item-75 objects are not at hand "
                    "(the slim extracts beside this file carry their sha256 pins)")
    src, sha = FIXTURES[name]
    raw = open(os.path.join(root, src), "rb").read()
    assert hashlib.sha256(raw).hexdigest() == sha
    o, d = json.loads(raw), _load_fixture(name)
    hc = o["house_controls"]["detail"][0]
    assert hc["control"] == d["control"] == "synthetic (3f0-7f1)/5"
    assert [int(x) for x in hc["want"]] == d["want"]
    assert [int(x) for x in hc["got"]] == d["got"]
    assert hc["ok"] is d["recorded_ok"] is False
    assert [[int(x) for x in r] for r in o["main"]["basis_relations"]] == d["basis_relations"]
    assert o["main"]["K"] == d["K"] and o["main"]["backend"] == d["backend"]


def test_controls_fills_its_entry_from_house_verdict():
    src = inspect.getsource(ML.controls)
    assert "house_verdict(got, want, rels)" in src
    assert '"control": "negative (random target)"' in src

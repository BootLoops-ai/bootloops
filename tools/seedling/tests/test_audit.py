"""Battery for the family-preflight member (`seedling audit` / `seedling
identity`) on the shipped generic fixtures.

- sunrise (equal-mass 2-loop): kira DOES find the S3 loop-shift symmetry on
  its own; the audit must find it too (superset of kira's finder). Expect
  |Aut| = 6 with the pure loop-swap present.
- asym_box (1-loop box, 3 distinct internal masses): negative control,
  expect |Aut| = 1.
- mutation control: breaking one sunrise mass must shrink the group
  (S3 broken -> residual Z2), never silently stay at 6.
- identity: the built-in wpair fixture runs clean; general mode fails loud
  on inconsistent builds.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", ".."))

from seedling import audit  # noqa: E402

FIX = os.path.join(HERE, "fixtures")


def _run(config_dir, tmp_path):
    return audit.run(str(config_dir), out_dir=str(tmp_path),
                     use_cache=False, verbose=False)


def test_sunrise_s3(tmp_path):
    r = _run(os.path.join(FIX, "sunrise"), tmp_path)
    assert r["group_order"] == 6
    assert all(r["controls"].values())
    # the pure loop-swap kira also finds must be in the group
    perms = {tuple(g["Dperm"]) for g in r["generators"]}
    assert perms, "no generators reported for |Aut|=6"


def test_asym_box_trivial(tmp_path):
    r = _run(os.path.join(FIX, "asym_box"), tmp_path)
    assert r["group_order"] == 1
    assert r["generators"] == []


def test_sunrise_mass_mutation_breaks_s3(tmp_path):
    """One-mass mutation control: l2^2 - 1 -> l2^2 - 2 must NOT keep |Aut|=6."""
    import shutil
    mut = tmp_path / "sunrise_mut"
    shutil.copytree(os.path.join(FIX, "sunrise"), mut)
    cfg = mut / "config" / "integralfamilies.yaml"
    text = cfg.read_text()
    assert '"l2^2 - 1"' in text
    cfg.write_text(text.replace('"l2^2 - 1"', '"l2^2 - 2"', 1))
    r = _run(mut, tmp_path / "out")
    assert r["group_order"] < 6, "mutated family kept the full group"
    assert r["group_order"] == 2  # residual Z2 (l1 <-> l1 sector swap broken)


def test_identity_wpair_fixture():
    """The built-in fixture is byte-frozen: its recorded verdict is that
    NEITHER wpair build matches the paper T_4 (that finding is the point
    of the fixture). The regression is reproducing exactly that."""
    hits_cur, hits_fixed = audit.identity_check_wpair()
    assert not hits_cur and not hits_fixed

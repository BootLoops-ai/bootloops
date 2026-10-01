"""Tests for tools/ratfit/degree_alias.py (ratfit.degree_alias) — the
grid-alias vs CRT-height failure-axis classifier for rational reconstruction.
Replays the two-sided proof shapes carried by the shipped receipt
(DEGREE_ALIAS_DISCRIMINATOR.json), the two fail-closed negatives, the receipt
cross-check, and both CLI entry forms."""
import os
import subprocess
import sys

import pytest

import ratfit.degree_alias as da

TOOLS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def test_side_a_alias_pinned_at_ceiling_trivial_den():
    r = da.classify((50, 100), [[((49, 99), (0, 0))] * 7] * 5)
    assert r['verdict'] == 'ALIAS' and r['axis'] == 'd'
    assert r['counts']['n_fits'] == 5 and r['counts']['alias'] == 5


def test_side_a_axis_eta_when_only_eta_pinned():
    r = da.classify((50, 100), [[((10, 99), (0, 0))] * 3] * 4)
    assert r['verdict'] == 'ALIAS' and r['axis'] == 'eta'


def test_side_b_height_below_ceiling_nontrivial_den():
    r = da.classify((29, 101), [[((12, 46), (9, 47))] * 5,
                                [((13, 45), (10, 48))] * 5,
                                [((11, 47), (9, 46))] * 5])
    assert r['verdict'] == 'HEIGHT' and r['axis'] is None


def test_cross_prime_disagreement_refuses():
    r = da.classify((50, 100), [[((49, 99), (0, 0)), ((30, 99), (0, 0))]])
    assert r['verdict'] == 'INCONCLUSIVE'


def test_pinned_num_with_nontrivial_den_is_not_alias():
    assert da.classify_cell((50, 100), [((49, 20), (3, 2))] * 5) == 'INCONCLUSIVE'
    assert da.classify((50, 100), [[((49, 20), (3, 2))] * 5])['verdict'] == 'INCONCLUSIVE'


def test_one_unpinned_cell_blocks_alias_verdict():
    # reference rule: ALIAS needs EVERY fitted cell pinned; a HEIGHT cell among
    # ALIAS cells (no inconclusives) yields HEIGHT, never ALIAS.
    cells = [[((49, 99), (0, 0))] * 3] * 4 + [[((12, 46), (9, 47))] * 3]
    assert da.classify((50, 100), cells)['verdict'] == 'HEIGHT'


def test_empty_inputs_refuse():
    assert da.classify_cell((50, 100), []) == 'INCONCLUSIVE'
    assert da.classify((50, 100), [])['verdict'] == 'INCONCLUSIVE'


def test_receipt_ships_as_package_data_and_matches_proof_sides():
    assert os.path.isfile(da.DISCRIMINATOR_PATH)
    assert os.path.dirname(da.DISCRIMINATOR_PATH) == os.path.dirname(os.path.abspath(da.__file__))
    r = da.discriminator()
    assert tuple(r['side_A_alias']['grid']) == (50, 100)
    assert tuple(r['side_B_height']['grid']) == (29, 101)


def test_builtin_selftest_passes(capsys):
    da._selftest()
    assert 'selftest PASS' in capsys.readouterr().out


@pytest.mark.parametrize('argv', [
    [sys.executable, '-m', 'ratfit.degree_alias', '--selftest'],
    [sys.executable, os.path.join(TOOLS, 'ratfit', 'degree_alias.py'), '--selftest'],
])
def test_cli_selftest_both_entry_forms(argv):
    env = dict(os.environ, PYTHONPATH=TOOLS + os.pathsep + os.environ.get('PYTHONPATH', ''),
               PYTHONDONTWRITEBYTECODE='1')
    p = subprocess.run(argv, capture_output=True, text=True, env=env, cwd=TOOLS, timeout=120)
    assert p.returncode == 0, p.stderr
    assert 'degree_alias selftest PASS' in p.stdout

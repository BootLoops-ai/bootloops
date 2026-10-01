"""Tests for tools/pathgrid.py — a reference production Pfit->P0safe path."""
from fractions import Fraction

import pytest

from pathgrid import coincidence_loci, design_grid

PA = dict(s12=-3, s23=-5, s34=-7, s45=-2, s15=-11, mm=1)
PB = dict(s12=-3, s23=-5, s34=-7, s45=-11, s15=-17, mm=2)
LOCI_EXPECTED = {Fraction(1, 9): {"s45", "s12"},
                 Fraction(1, 3): {"s45", "s23"},
                 Fraction(5, 9): {"s45", "s34"}}


def test_coincidence_loci_production_path():
    loci = coincidence_loci(PA, PB)
    assert set(loci) == set(LOCI_EXPECTED)
    for t, pair in LOCI_EXPECTED.items():
        assert set(loci[t].split("=")) == pair


def test_loci_are_fraction_exact():
    loci = coincidence_loci(PA, PB)
    assert all(isinstance(t, Fraction) and 0 < t < 1 for t in loci)


def test_whole_path_coincidence_errors_loudly():
    PA2 = dict(PA, s23=-3)   # s12 == s23 at BOTH endpoints
    PB2 = dict(PB, s23=-3)
    with pytest.raises(ValueError, match="whole-path"):
        coincidence_loci(PA2, PB2)


def test_design_grid_avoids_loci_poles_and_old_nodes():
    grid = design_grid(PA, PB, 30, exclude=("1/2", Fraction(7, 60)),
                       poles=("11/60",))
    ts = [Fraction(p["t"]) for p in grid]
    assert len(grid) == 30 and len(set(ts)) == 30
    assert all(0 < t < 1 for t in ts)
    banned = set(LOCI_EXPECTED) | {Fraction(1, 2), Fraction(7, 60),
                                   Fraction(11, 60)}
    assert banned.isdisjoint(ts)
    for p in grid:
        assert p["checked"] == ["no-coincidence", "no-pole", "fresh"]


def test_design_grid_denominator_9_binds_all_three_loci():
    # k/9 candidates: k=1..8 minus loci {1/9, 3/9=1/3, 5/9} -> exactly 5 left
    grid = design_grid(PA, PB, 5, denominators=(9,))
    ts = {Fraction(p["t"]) for p in grid}
    assert ts == {Fraction(2, 9), Fraction(4, 9), Fraction(6, 9),
                  Fraction(7, 9), Fraction(8, 9)}


def test_design_grid_errors_loudly_on_shortfall():
    with pytest.raises(ValueError, match="only 5/6"):
        design_grid(PA, PB, 6, denominators=(9,))


def test_design_grid_dedupes_equal_rationals():
    # 30/60 == 1/2: with 1/2 excluded, no denominator may re-emit it
    grid = design_grid(PA, PB, 100, denominators=(60, 120), exclude=("1/2",))
    ts = [Fraction(p["t"]) for p in grid]
    assert Fraction(1, 2) not in ts and len(set(ts)) == len(ts)

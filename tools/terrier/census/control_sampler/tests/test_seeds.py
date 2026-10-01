"""Fixed-seed reproducibility battery (master seed, derivation, golden draws)."""
import hashlib
import os
import sys
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import seeds as S

SEEDS_JSON = os.environ.get("TERRIER_SEEDS_JSON", "")
if not SEEDS_JSON:
    raise SystemExit("REFUSE: env TERRIER_SEEDS_JSON unset — must point at "
                     "the seeds JSON file (reference data not included in "
                     "the package)")
COMMITTED = "e793a915c531ca1a230a50716931336c7d8ad342fe33db31e075bd0f231bf248"
# Golden draws for stream (G4, 1), recorded once from the master seed; PCG64's
# algorithm is version-pinned in numpy, so these are machine-independent.
GOLDEN_G4_BIN1 = [9331831188815203097, 11018579269093999845,
                  12248436173302301583]


def test_master_matches_committed():
    assert S.load_master(SEEDS_JSON) == int(COMMITTED, 16)


def test_derivation_pinned_formula():
    m = int(COMMITTED, 16)
    tag = int(hashlib.sha256(b"G4|3").hexdigest()[0:16], 16)
    assert S.stream_seed(m, "G4", 3) == m ^ tag


def test_reproducibility_and_separation():
    m = S.load_master(SEEDS_JSON)
    d1 = [S.raw64(S.bin_generator(m, "G4", 1)) for _ in range(1)]
    g1 = S.bin_generator(m, "G4", 1)
    g2 = S.bin_generator(m, "G4", 1)
    a = [S.raw64(g1) for _ in range(10)]
    b = [S.raw64(g2) for _ in range(10)]
    assert a == b and a[0] == d1[0]
    g3 = S.bin_generator(m, "G4", 2)
    assert [S.raw64(g3) for _ in range(10)] != a
    g4 = S.bin_generator(m, "G5", 1)
    assert [S.raw64(g4) for _ in range(10)] != a


def test_golden_first_draws():
    m = S.load_master(SEEDS_JSON)
    g = S.bin_generator(m, "G4", 1)
    assert [S.raw64(g) for _ in range(3)] == GOLDEN_G4_BIN1


def test_dyadic_uniform_exact():
    m = S.load_master(SEEDS_JSON)
    g = S.bin_generator(m, "G4", 1)
    u = S.dyadic_uniform(g)
    assert isinstance(u, Fraction)
    assert 0 <= u < 1
    assert (u * S.TWO64).denominator == 1
    assert u == Fraction(GOLDEN_G4_BIN1[0], S.TWO64)

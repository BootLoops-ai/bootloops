# galois member (subpackage) — exact MZV relation-ring engine.
# RING_BANK and load_ring read the pkl bank from the env-configured home.
"""mzvring — exact MZV relation-ring engine (depth<=3, NO PSLQ).

The relation-ring MEMBER of the GALOIS package (tools/galois/).

Canonicalizes {zeta(n), mzv(a,b), mzv(a,b,c)} monomials over Q against the
regularized double-shuffle relation families (R2-R6, see mzvring.py docstring):
per-weight Gauss elimination with single m2/m3 symbols as the only eliminable
columns; survivors = the depth<=3 basis (Broadhurst-Kreimer dims reproduced).

DATA PLACEMENT (tools link, never copy; large reference data
stays outside the package): the reference ring pickles live in RING_BANK below
(ring17.pkl / ring21.pkl / ring23.pkl). Every consumer passes or defaults to
an explicit pkl path via load_ring(); nothing here duplicates the bank.

API: engine names re-exported from .mzvring (Ring, znorm, dsum/dscale/dmul,
comp_word/word_comp/shuffle, elim_key, mweight, ...) plus:
  load_ring(path=None)  -> Ring with stored canon/survivors tables
                           (default: RING_BANK/ring21.pkl)
  RING_BANK             -> canonical bank directory (env-configured)
Build new tables: python3 tools/galois/mzvring/build_ring.py --wmax N --out X.pkl
"""
import os as _os
import pickle as _pickle

from .mzvring import *          # noqa: F401,F403 — the engine, by identity
from . import mzvring as _mz

RING_BANK = _os.environ.get('GALOIS_RING_BANK')  # reference ring pkl dir (not shipped; build your own via build_ring.py)


def load_ring(path=None):
    """Load stored canonicalization tables -> Ring (no rebuild).

    Default: RING_BANK/ring21.pkl (weights 3..21, validated 1298/1298
    eliminated symbols at 150d — ringval_results.json in the bank)."""
    if path is None:
        if not RING_BANK:
            raise RuntimeError(
                "GALOIS_RING_BANK is not set and no explicit path was given. "
                "The ring pickles are reference data "
                "that do not ship with this repo — build one with "
                "mzvring/build_ring.py (--wmax 21 is ~30 min single-core) and "
                "point GALOIS_RING_BANK at its directory, or pass load_ring(path).")
        path = _os.path.join(RING_BANK, 'ring21.pkl')
    R = _mz.Ring.__new__(_mz.Ring)
    with open(path, 'rb') as f:
        T = _pickle.load(f)
    R.canon, R.survivors, R.wmax = T['canon'], T['survivors'], T['wmax']
    return R

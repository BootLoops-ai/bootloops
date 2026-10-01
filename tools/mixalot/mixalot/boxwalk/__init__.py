"""mixalot.boxwalk — exact rational evaluation of Z(u) = int_{[0,1]^n} prod_k
P_k(x)^{u_k} dx for integer polynomials P_1..P_m in Z[x_1..x_n] (n <= 6)
and target u in N^m, via a mod-p contiguity walk on the exponent lattice
+ CRT + rational reconstruction.  "Kira for box-moment contiguity."

Home: tools/mixalot/mixalot/boxwalk/ (a mixalot member; this subpackage is
the only copy).  `from mixalot import boxwalk` with tools/mixalot on sys.path;
CLI: `python3 -m mixalot.boxwalk <selftest|plan|produce|fiber|verify> ...`.

API:
  spec = boxwalk.load_spec(path_or_dict)
  pl   = boxwalk.plan(spec)                      # segment planner + programs
  res  = boxwalk.walk(spec, pl, primes)          # window/value mod p (batch)
  man  = boxwalk.produce(spec, plan=pl, nprimes=8, ...)   # exact Fraction,
                                                 # gated (manifest dict)
  rec  = boxwalk.emit_fiber_recurrence(spec, ...)  # scalar operator emission
  rep  = boxwalk.verify_manifest(spec, man_or_path, nfresh=2)  # independent
                                                 # replay of a recorded manifest

Status: validated on the shipped selftest battery (selftest.py; mixalot
battery leg M16 runs it).
"""
import os as _os

from .core import Spec, primes31, Z_mod, moment_window
from .planner import plan, default_order
from .driver import walk, produce
from .fiber import emit_fiber_recurrence, gcrd_reduce
from .verify import verify_manifest

__version__ = '1.0.0'


def load_spec(src):
    """Load a problem spec from a JSON path or a dict."""
    if isinstance(src, dict):
        return Spec(src)
    return Spec.load(src)


def home():
    """Directory of this member (cli.py, selftest.py, examples/ and the
    reference SELFTEST.json live here)."""
    return _os.path.dirname(_os.path.abspath(__file__))

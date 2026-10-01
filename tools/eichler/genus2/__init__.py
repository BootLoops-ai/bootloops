"""eichler.genus2 — genus-2 curve identification and recognition hygiene.

Members (pure sympy/mpmath, no Sage dependency):

  mestre_port            Clebsch / Igusa-Clebsch invariants of a binary sextic and
                         Mestre's conic-and-cubic reconstruction of the curve from
                         its invariants (import-safe module).
  invariant_harness      recognition-free arbitration between candidate sextics by
                         absolute invariants (+ Richelot orbit) against a high-dps
                         numeric oracle json (script; needs G2KIT_TRUE_JSON, alias
                         EICHLER_GENUS2_TRUE_JSON; refuses loudly without it).
  conic_fast             exact quadratic-field conic diagonalization and point
                         search (script; runs its worked example on import).
  pslq_law_template      the PSLQ recognition law as reusable structure (template;
                         refuses loudly without the pipeline it drives).
  exact_transport_driver exact-rational-basepoint DE transport pattern on top of
                         tools/wayfinder (template; needs your own system class).

Only `mestre_port` is imported here; the other members are scripts or templates
that execute (or refuse) on import and are run as files. Battery:
`python3 selftest.py` in this directory, or the package-level
`tools/eichler/selftest.py`. See GUIDE.md.
"""
import os as _os
import sys as _sys

_HERE = _os.path.dirname(_os.path.abspath(__file__))
if _HERE not in _sys.path:          # members import each other by flat name
    _sys.path.insert(0, _HERE)

from . import mestre_port  # noqa: E402,F401

__all__ = ['mestre_port']

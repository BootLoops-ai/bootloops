# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision.
"""longhand.disteval -- Longhand's route B: the pySecDec disteval driver.

One compiled pySecDec `loop_package` (+ `make disteval`) evaluated at one
Euclidean kinematic point as three checkpointed disteval stages (two lattice
settings and the median-lattice second integrator), chunked by sector group so
every chunk is a resumable checkpoint, then a pairwise per-order compare and a
conservative acceptance check with a planted-fail control and a foreign
(independent coarse evaluation) control.  All numerical work is pySecDec's
(GPL-3.0, installed by the user; never bundled); the staging, the compare and
the check are this package's.

Members (each also a command line; see ../DISTEVAL.md):
  pysecdec_point    the command line of the whole chain (--stage A|B|C|compare|check)
  stage_chunks      the per-sector chunk splitter
  assemble_compare  the per-order assembly of chunk results and the pairwise compare
  gate              the acceptance check and its controls (historical module name;
                    `check` is accepted as an alias attribute of this package)
  memory_basis      the memory figure to declare for W workers, from measured maxrss

The worked example whose records ship with the package is the generic-mass
five-point double-pentagon top sector (family ndpent_top) under
../examples/ndpent_top/; PACKAGE_PINS.json beside this file pins that
example's compiled package.
"""
import importlib
import os as _os
import sys as _sys

__all__ = ['pysecdec_point', 'stage_chunks', 'assemble_compare', 'gate', 'check', 'memory_basis']
_HERE = _os.path.dirname(_os.path.abspath(__file__))


def __getattr__(name):
    # members load on first touch (PEP 562); they import each other by bare
    # name, so this directory goes on sys.path once, the way the CLI does it
    target = 'gate' if name == 'check' else name
    if target in ('pysecdec_point', 'stage_chunks', 'assemble_compare', 'gate', 'memory_basis'):
        if _HERE not in _sys.path:
            _sys.path.insert(0, _HERE)
        mod = importlib.import_module(target)
        _sys.modules.setdefault(__name__ + '.' + target, mod)
        globals()[target] = mod
        if name == 'check':
            globals()['check'] = mod
        return mod
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

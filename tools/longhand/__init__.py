# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision.
"""longhand -- independent numerical ground truth for Feynman integrals, done
the slow honest way.  Two routes in one package:

  route A  the arbitrary-precision Feynman-parametric evaluator (this
           directory's flat modules): `hiprec` (the engines: fixed
           Gauss-Legendre product, nested tanh-sinh, parallel CBC rank-1
           lattice QMC, all in arbitrary precision), `feynman` (Symanzik U/F
           from a propagator list via pySecDec + closed-form elimination of
           F-linear parameters; needs sympy, and pySecDec for `build_UF`
           only), `parametric` (the UF-spec JSON front end), `bench_box` /
           `bench_se3l` (the two benchmarks; `box` and `box_crosschart` are
           the 40-digit one-loop box oracle and its independent validation).
  route B  `longhand.disteval`: the pySecDec disteval driver -- one compiled
           loop_package at one Euclidean point in three checkpointed stages,
           a three-way lattice compare and an acceptance check (command line
           disteval/pysecdec_point.py; pySecDec, GPL-3.0, user-installed).

With the repository's tools/ directory on sys.path:

    from longhand import box, integrate_gauss_product, integrate_qmc, integrate_spec
    v = box(-1, Fraction(-1, 3), 1, 2, dps=35)          # exact kinematics in, ~35 digits out

The flat modules keep their historical bare import names as well
(`sys.path.insert(0, 'tools/longhand'); from hiprec import integrate_qmc`),
which is how the benches and selftests import them; this file puts the
package directory on sys.path once so both spellings name the same module
objects.  See GUIDE.md; manuals PARAMETRIC.md (route A) and DISTEVAL.md
(route B); battery: selftest.py.
"""
import importlib
import os as _os
import sys as _sys

_HERE = _os.path.dirname(_os.path.abspath(__file__))
if _HERE not in _sys.path:
    # the flat members import each other by bare name (parametric -> hiprec,
    # bench_box -> hiprec), and multiprocessing workers re-import them by that
    # name under the spawn start method: one path entry serves both
    _sys.path.insert(0, _HERE)

import hiprec        # noqa: E402  (gmpy2 is a hard requirement: an ImportError here names it)
import parametric    # noqa: E402
import bench_box     # noqa: E402

for _m in (hiprec, parametric, bench_box):
    _sys.modules.setdefault(__name__ + '.' + _m.__name__, _m)

from hiprec import (integrate_gauss_product, integrate_tanhsinh,   # noqa: E402
                    integrate_qmc, cbc_generating_vector)
from parametric import (load_spec, check_spec, make_eval,          # noqa: E402
                        spec_factory, integrate_spec, box_spec)
from bench_box import box, box_crosschart                          # noqa: E402

__version__ = "1.0.0"
__all__ = [
    'integrate_gauss_product', 'integrate_tanhsinh', 'integrate_qmc',
    'cbc_generating_vector', 'load_spec', 'check_spec', 'make_eval',
    'spec_factory', 'integrate_spec', 'box_spec', 'box', 'box_crosschart',
    'hiprec', 'parametric', 'bench_box', 'feynman', 'bench_se3l', 'disteval',
]
_LAZY = ('feynman', 'bench_se3l')   # need sympy (and pySecDec at call time): load on first touch


def __getattr__(name):
    # PEP 562: `longhand.feynman` imports tools/longhand/feynman.py on first
    # access, so `import longhand` itself never needs sympy or pySecDec
    if name in _LAZY:
        mod = importlib.import_module(name)
        _sys.modules.setdefault(__name__ + '.' + name, mod)
        globals()[name] = mod
        return mod
    if name == 'disteval':
        mod = importlib.import_module(__name__ + '.disteval')
        globals()[name] = mod
        return mod
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

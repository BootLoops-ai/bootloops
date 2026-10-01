#!/usr/bin/env python3
# Compatibility shim — canonical module: tools/lockpick/mplll_graded.py (page tools/lockpick.html).
# Import and CLI forms at this flat path forward there; write new code against lockpick.mplll_graded.
# Self-anchoring: the package is resolved from this file's own directory (tools/),
# so the shim works no matter which sys.path entry found it.
# No __main__ in the original — script mode stays a no-op (import-only shim).
# Reload note: importlib.reload() of this shim re-copies attributes but does NOT
# reload the canonical module — reload lockpick.mplll_graded instead.
import os as _os
import sys as _sys
_d = _os.path.dirname(_os.path.abspath(__file__))
if _d not in _sys.path:
    _sys.path.insert(0, _d)
import lockpick.mplll_graded as _m
_sys.modules[__name__].__dict__.update(
    {k: v for k, v in _m.__dict__.items() if not k.startswith("__")})
__doc__ = _m.__doc__

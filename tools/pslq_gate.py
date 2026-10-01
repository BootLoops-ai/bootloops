#!/usr/bin/env python3
# Compatibility shim — canonical module: tools/lockpick/pslq_gate.py (page tools/lockpick.html).
# Import and CLI forms at this flat path forward there; write new code against lockpick.pslq_gate.
# Self-anchoring: the package is resolved from this file's own directory (tools/),
# so the shim works no matter which sys.path entry found it.
# Script mode (python3 tools/pslq_gate.py ...) still runs the CLI — forwarded via runpy.
# Reload note: importlib.reload() of this shim re-copies attributes but does NOT
# reload the canonical module — reload lockpick.pslq_gate instead.
import os as _os
import sys as _sys
_d = _os.path.dirname(_os.path.abspath(__file__))
if _d not in _sys.path:
    _sys.path.insert(0, _d)
if __name__ == "__main__":
    # script mode: runpy ONLY (pre-importing the submodule here would trigger
    # runpy's found-in-sys.modules RuntimeWarning and a double execution)
    import runpy as _runpy
    _runpy.run_module("lockpick.pslq_gate", run_name="__main__", alter_sys=True)
else:
    import lockpick.pslq_gate as _m
    _sys.modules[__name__].__dict__.update(
        {k: v for k, v in _m.__dict__.items() if not k.startswith("__")})
    __doc__ = _m.__doc__

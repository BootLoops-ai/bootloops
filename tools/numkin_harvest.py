#!/usr/bin/env python3
# Compatibility shim — canonical module: tools/numkin/numkin_harvest.py.
# Import and CLI forms at this flat path forward there; write new code against numkin.numkin_harvest.
# Self-anchoring: the package is resolved from this file's own directory (tools/),
# so the shim works no matter which sys.path entry found it.
# Script mode (python3 tools/numkin_harvest.py ...) still runs the CLI.
# Side-effect contract preserved: importing this shim imports the real module, which
# sys.path-inserts the kira_parse home, tools/frobenius-boundary/frobenius_boundary/
# (an older flat tools/frobenius_boundary/ layout is probed first and is absent from
# this tree) — so `from kira_parse import ...` after `import numkin_harvest` keeps working.
# Reload note: importlib.reload() of this shim re-copies attributes but does NOT
# reload the canonical module — reload numkin.numkin_harvest instead.
import os as _os
import sys as _sys
_d = _os.path.dirname(_os.path.abspath(__file__))
if _d not in _sys.path:
    _sys.path.insert(0, _d)
import numkin.numkin_harvest as _m
_sys.modules[__name__].__dict__.update(
    {k: v for k, v in _m.__dict__.items() if not k.startswith("__")})
__doc__ = _m.__doc__
if __name__ == "__main__":
    _sys.exit(_m.main())

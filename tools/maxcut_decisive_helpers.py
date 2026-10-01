# Compatibility shim — canonical module: tools/maxcut/decisive_helpers.py.
# Import and CLI forms at this flat path forward there; write new code against maxcut.decisive_helpers.
# RELOAD SEMANTICS: importlib.reload(<this shim>) re-copies from the CACHED package
# module and never re-reads disk — to pick up a code edit, reload the package
# module FIRST (importlib.reload(maxcut.decisive_helpers)), then reload this shim.
import os as _os, sys as _sys
_anchor = _os.path.dirname(_os.path.abspath(__file__))  # tools/ (package parent)
if _anchor not in _sys.path:                            # guarded: no duplicates
    _sys.path.insert(0, _anchor)
import maxcut.decisive_helpers as _m
__doc__ = _m.__doc__
_sys.modules[__name__].__dict__.update({k: v for k, v in _m.__dict__.items() if not k.startswith('__')})

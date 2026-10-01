#!/usr/bin/env python3
"""oracle_registry.py — flat alias of gatekeeper.registry.

The engine lives in the gatekeeper package:
    tools/gatekeeper/registry.py   (python -m gatekeeper.registry)

This flat path is an attribute-complete, script-mode-CLI shim, so all of the
following call forms work:
  * by-path exec  — `python3 tools/oracle_registry.py --rebuild|--point|--dupes|--selftest`
  * import        — `import oracle_registry; oracle_registry.rebuild()/query()/dupes()/…`
  * env overrides — ORACLE_REG_ROOTS / ORACLE_REG_REGISTRY / ORACLE_REG_COMBOS
                    are read by the engine at import time.

New code should import/invoke gatekeeper.registry directly.
"""
import os as _os, sys as _sys

# Make the gatekeeper package importable regardless of how the shim is
# reached (script mode already puts tools/ on sys.path[0]; this covers the
# import-from-elsewhere case too).
_HERE = _os.path.dirname(_os.path.abspath(__file__))
if _HERE not in _sys.path:
    _sys.path.insert(0, _HERE)

from gatekeeper import registry as _registry

# Attribute-complete re-export: every public name of the engine (module
# constants ROOT_GLOBS/REGISTRY/COMBOS_DIR/MAX_BYTES/MAX_DEPTH and the
# functions frac/sha1_of/load_combos/parse_file/rebuild/load_reg/query/dupes/
# parse_masters/selftest/main) resolves through this shim identically.
for _name in dir(_registry):
    if not _name.startswith("__"):
        globals()[_name] = getattr(_registry, _name)

if __name__ == "__main__":
    _registry.main()

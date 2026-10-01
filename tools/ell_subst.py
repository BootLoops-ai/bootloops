# Compatibility shim — canonical module: tools/gatekeeper/probes/ell_subst.py.
# The import form at this flat path forwards there; write new code against gatekeeper.probes.ell_subst.
# (No CLI: neither this shim nor the canonical module has a main()/__main__ block —
# `python3 tools/ell_subst.py` is a silent no-op.)
import gatekeeper.probes.ell_subst as _m
import sys as _sys
_sys.modules[__name__].__dict__.update({k: v for k, v in _m.__dict__.items() if not k.startswith('__')})

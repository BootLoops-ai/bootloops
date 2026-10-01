#!/usr/bin/env python3
"""Compat shim: gravityflow forwards to PMflow.
Every existing call form keeps working; the tool is tools/pmflow/pmflow.py."""
import os, sys
sys.stderr.write("[gravityflow] forwarding to pmflow.py (PMflow)\n")
_p = os.path.join(os.path.dirname(os.path.abspath(__file__)), "pmflow.py")
sys.argv[0] = _p
__file__ = _p  # sidecar reads resolve to the engine dir
exec(compile(open(_p).read(), _p, "exec"))

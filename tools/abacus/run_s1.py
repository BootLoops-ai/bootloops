#!/usr/bin/env python3
"""COMPAT SHIM — executes the canonical body in build/abcount/ under this
module name.

Forwards the abcount CLI form `python3 tools/abacus/run_s1.py` to the
canonical byte-identical harness at tools/abacus/build/abcount/run_s1.py
(see tools/abacus/README.md).
The harness resolves its sha-pin paths from its own file location."""
import os
import sys

_BODY = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                     "build", "abcount", "run_s1.py")
os.execv(sys.executable, [sys.executable, _BODY] + sys.argv[1:])

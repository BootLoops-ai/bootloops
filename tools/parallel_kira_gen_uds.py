#!/usr/bin/env python3
# Compatibility shim for the flat path; engine: tools/kira-stack/parallel_kira_gen_uds.py (host tool: kira-stack).
import os, sys
_p = os.path.join(os.path.dirname(os.path.abspath(__file__)), "kira-stack", "parallel_kira_gen_uds.py")
sys.argv[0] = _p
__file__ = _p  # sidecar reads resolve to the engine dir
exec(compile(open(_p).read(), _p, "exec"))

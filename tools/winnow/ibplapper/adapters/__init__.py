"""ibplapper.adapters — everything domain/tool-specific lives here.

kira.py         Kira SYSTEM_*.gz sidecar (Route A: an undocumented internal
                dump — version-guarded, empirical decode flagged honestly).
weights.py      the empirical kira weight decode (version-pinned; Route A
                only).
user_system.py  kira user_defined_system reader (Route B: kira's input
                format for user-defined systems, specified in the Kira
                2.0 paper; rules the paper leaves open follow Kira's
                parser; no version guard, no weight decode). Its file
                reader follows Kira's parser
                (src/pyred/interface.cpp, src/pyred/integrals.cpp), so
                this one module is GPL-3.0-or-later; see its header and
                NOTICE in the package root.
"""
from . import kira, user_system, weights  # noqa: F401

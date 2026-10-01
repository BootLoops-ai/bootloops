"""ibplapper.adapters — everything domain/tool-specific lives here.

kira.py         Kira SYSTEM_*.gz sidecar (Route A: an undocumented internal
                dump — version-guarded, empirical decode flagged honestly).
weights.py      the empirical kira weight decode (version-pinned; Route A
                only).
user_system.py  kira user_defined_system reader (Route B: the documented
                stable input format — no version guard needed, nothing
                empirical to pin).
"""
from . import kira, user_system, weights  # noqa: F401

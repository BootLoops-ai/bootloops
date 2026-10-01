#!/usr/bin/env python3
"""
opderive.py (periods wing) — operator-derivation WRAPPER.

DELEGATION LAW: the engine is tools/annihilator/ (pf_from_series,
find_recurrence_fast, ...) — never duplicated.  This module re-exports the
staging layer pipeline/restrict_op.py: theta-order-first annihilator pin at
2 primes x 2 truncations, multi-prime CRT + Wang rational reconstruction,
EXACT full-series annihilation proof with true held-out tail, theta-form +
integer-root indicial factorization + MUM-block gate.  Mod-p holonomic rank
probes delegate to tools/pf_rank.jl (Griffiths-Dwork route).
Battery: selftest_opderive.py = run_pipe g3 (theta-form + recurrence
== the stored pipeline/dkmm/operator_LS.json exactly, up to overall sign).
"""
import os, sys
_PIPE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "pipeline")
sys.path.insert(0, _PIPE)
from restrict_op import find_operator                    # noqa: F401

# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision.
"""loomcheck — member of the Dogtag package (tools/dogtag): mechanical
Yangian/loom/fishnet applicability screen for position-space conformal
integrals (internal-vertex conformal weight, planarity of the full graph
with numerator edges counted, the level-one-momentum face condition).
SCREEN-ONLY: a FAIL is a proof of exclusion from the cited theorems as
published; a PASS constructs nothing.  See loomcheck/loomcheck.py and the
package GUIDE.md ("Member: loomcheck").

    python3 -m loomcheck --selftest                       (cwd = tools/dogtag)
    python3 -m loomcheck --basis 4Loop_int_basis.m [--targets ...] [--json OUT]
    from loomcheck import screen_graph
"""
from .loomcheck import (DEFAULT_BASIS, NINE, load_entries, main, parse,  # noqa: F401
                        screen_graph, selftest, selftest_cases)

__all__ = ["DEFAULT_BASIS", "NINE", "load_entries", "parse", "screen_graph",
           "selftest", "selftest_cases", "main"]

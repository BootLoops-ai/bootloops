# loomcheck — Yangian/loom/fishnet applicability screen (member of Dogtag)

Mechanical screen (SCREEN-ONLY): internal-vertex conformal weight, full-graph planarity
including numerator edges, and the 2505.05550 P-hat face condition. FAIL = proof of
exclusion from the cited Yangian/SoV theorems; PASS does not construct anything.
The nine default targets screen 0/9 in any proven class (reproduce with
`--basis 4Loop_int_basis.m` from the arXiv:2607.11645 ancillary archive).
GATED: face check counts the outer face (conservative; margin >=7 on the nine).

This screen lives inside the Dogtag package as the subfolder `loomcheck/`
(one file of logic, `loomcheck.py`, importable as the `loomcheck` module). From the
package directory `tools/dogtag/`:

- Battery: `python3 -m loomcheck --selftest` (or `python3 loomcheck/loomcheck.py --selftest`)
  — structural known-answer checks on hand-checkable graphs, one per screened condition;
  no data files needed; exit 0/1. The package self-test `python3 topology_audit.py
  --selftest` runs the same battery as its `[loomcheck]` leg, and `python3 -m pytest
  tests/test_loomcheck.py` covers it case by case plus the CLI surface.
- Screen: `python3 -m loomcheck --basis 4Loop_int_basis.m [--targets 144,163,...] [--D 4]
  [--json RECEIPT.json]` (or set `LOOMCHECK_BASIS`).
- API: `from loomcheck import screen_graph` with `tools/dogtag` on `sys.path`.

Requires `networkx` (planarity + face traversal), the package's declared dependency.
Guide: the package `GUIDE.md`, section "Member: loomcheck".

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).

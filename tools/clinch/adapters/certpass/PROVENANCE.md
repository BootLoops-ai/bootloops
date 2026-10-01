# certpass — provenance

Code-only adapter member (*.py + finalize.sh) absorbed from an external
point-process certification study. PROVENANCE.sha256 lists the sha256 of
every shipped file.

NOT shipped (stays with the originating study, resolved via the
CERTPASS_REFERENCE_DIR environment variable; the scripts refuse loudly when
it is unset):
- the study's data and fit artifacts (params/fit JSONs, feature cubes, dump
  files) — RECEIPTS_DIR/DATA_DIR in __init__.py resolve relative to this
  dir; full re-runs need the study tree (or a copy here);
- the study's helper modules (fit_multipliers, playoff) imported by the
  run_*/assemble_* scripts;
- the study's certificate receipts.

The oracles and gates (oracle_*.py, gates_*.py) are self-contained.

Third-party notice. `upper_gamma_ext` in `twin_etas.py` (the extended upper
incomplete gamma function, seven lines; imported again by `assemble_etas_np.py`)
and `round_half_up` in `assemble_etas_np.py` (three lines, from
`etas/mc_b_est.py`) are copied from the `etas` package,
github.com/lmizrahi/etas @ 097f08b6 — Copyright (c) 2024 ETH Zurich, Leila
Mizrahi, MIT License (the same terms as this repository's LICENSE; the MIT
permission notice there applies to these functions with the copyright line
above). The pair-table conventions in `assemble_etas_np.py` and the parameter
bounds in `run_etas.py` / `recertify_taupinned.py` follow that package
(L. Mizrahi, S. Nandan, S. Wiemer, Seismol. Res. Lett. 92 (2021) 2333) without
copying its code beyond the two functions named above, and `etas_dump.py`
imports the user-installed package at run time; it is not bundled. See
THIRD_PARTY.md sections A, B2 and C at the repository root.

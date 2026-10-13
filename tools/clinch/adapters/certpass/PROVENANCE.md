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

Third-party notice. The ETAS files of this adapter follow the code of the
`etas` package of L. Mizrahi, S. Nandan and S. Wiemer
(github.com/lmizrahi/etas @ 097f08b6; Seismol. Res. Lett. 92 (2021) 2333),
Copyright (c) 2024 ETH Zurich, Leila Mizrahi, MIT License (the same terms as
this repository's LICENSE; the MIT permission notice there applies to the
portions listed here with the copyright line above):

- copied unchanged: `upper_gamma_ext` in `twin_etas.py` (the extended upper
  incomplete gamma function, seven lines, from `etas/inversion.py`; imported
  again by `assemble_etas_np.py`), with its `scipy.special` import line;
  `round_half_up` in `assemble_etas_np.py` (three lines, from
  `etas/mc_b_est.py`); and the parameter bounds in `run_etas.py`
  (`RANGES_LO`, `RANGES_HI`) and `recertify_taupinned.py` (`LTAU_BOUND`),
  which are the values of `etas/inversion.py`;
- rewritten from `etas/inversion.py` in vectorized numpy form, statement by
  statement in places: the pair table of `assemble_etas_np.py` (after
  `filter_catalog`, `calculate_distances`, `haversine` and `coppersmith`),
  its E-step check (after `triggering_kernel` and `expectation_step`), and
  the source integral (`G_sum` in `twin_etas.py` and the `G` block of
  `assemble_etas_np.py`, after `expected_aftershocks`). About 90 lines in
  all. The conventions these blocks reproduce (half-up rounding, the
  half-bin reference magnitude, the source-length cut and its constants,
  the Earth radius, the half-open time window, the strict time ordering of
  pairs) are taken from the code;
- `etas_dump.py` imports the user-installed package at run time; the
  package is not bundled.

The ball-arithmetic oracle (`oracle_etas.py`), its derivatives and the
certificates have no counterpart in the `etas` package. See THIRD_PARTY.md
sections A, B2 and C at the repository root.

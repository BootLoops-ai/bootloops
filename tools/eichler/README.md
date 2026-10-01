# eichler — member homes for the Eichler.jl family

This directory holds three members; the Eichler.jl engine itself lives at
`../../upgrades/Eichler.jl` (override with the `EICHLER_PROJECT` env var).
Full guide: `GUIDE.md`.

- `theta9/`: certified Siegel theta-constant observables on Eichler.jl's
  `acb_theta` layer (scan, frame, co-pair and certified-residual scripts; the
  reference values they compare against ship in `theta9/fixtures/`).
- `mirror6/`: exact gmpy2-mpq mirror-map / instanton-type fingerprint engine
  for MUM-type points of D-finite operators.
- `genus2/`: genus-2 curve identification and recognition hygiene in pure
  sympy/mpmath — Igusa-Clebsch invariants and Mestre reconstruction
  (`mestre_port`, a sympy transliteration of SageMath's mestre and invariants
  modules (GPL-2.0-or-later; see `LICENSE-GPL-2.0` and `NOTICE`), checked exactly
  against Sage's published doctest values), recognition-free candidate
  arbitration against a numeric oracle (`invariant_harness`, needs
  `G2KIT_TRUE_JSON`), exact quadratic-field
  conic point search (`conic_fast`), and the PSLQ recognition law and
  exact-basepoint transport as read-and-adapt templates. `import eichler.genus2`
  with `tools/` on PYTHONPATH; member guide in `genus2/GUIDE.md`.

Self-check: `python3 selftest.py --outdir <scratch>` from this directory runs
every member's battery in one command — the mirror6 positive control (28 boolean
control rows, sub-second), the genus2 legs (Sage-vector invariants, conic
worked example, two fail-closed refusals), a check that `mirror6/run_l6.py`
refuses by name without the sha-pinned `L6_exact.json` operator file (not
included), and a named skip for theta9, which ships without a battery (its
reference receipts are not part of this repository). About 2 s.

CREDIT: theta9 rests on FLINT's certified `acb_theta` module by Jean Kieffer [Kief]
(algorithm: Elkies and Kieffer [EK]) through Eichler.jl; mirror6 computes the mirror map
and instanton-type numbers in the sense of Candelas, de la Ossa, Green and Parkes [CdGP]
for operators of the AESZ Calabi–Yau tables [AESZ, vS]. Bracketed keys resolve in
REFERENCES.md at the repository root.

## License

MIT License, except `genus2/mestre_port.py`, which is GPL-2.0-or-later (derived from SageMath's `hyperelliptic_curves.mestre` and `invariants` modules, Copyright (C) 2011, 2012, 2013 Florian Bouyer, Marco Streng, 2008 Nick Alexander, 2025 Sabrina Kunzweiler, Gareth Ma, Giacomo Pope; see `NOTICE` and `LICENSE-GPL-2.0` in this directory). Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).

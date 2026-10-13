# Eichler.jl — packaging and attribution notes

## What this is

Eichler.jl is a package written for BootLoops. It is not a fork of an external engine,
so there is no patch set against an upstream tree. Public programs for parts of the same
tasks existed before it (listed under "Prior public code" below), and one block of
`src/siegel.jl` follows FLINT's source (listed under "Attribution"). It is the
certified-arithmetic layer of the elliptic Feynman-integral bootstrap: Γ₁(6) modular
forms and iterated Eichler integrals, the all-orders sunrise representation of
Adams–Weinzierl, elliptic-curve period/quasi-period evaluation, a generic validated
Frobenius/Picard–Fuchs transport engine (elliptic through Calabi–Yau/K3 operators, plus
a genus-2 Siegel layer), and PSLQ-ready constant dictionaries. Every public function
returns honest Arb/Acb ball enclosures; the tail bounds are proved inequalities
evaluated in Mag/Arb arithmetic. Validation results: `NUMERICS-RESULTS.md`. Adversarial
audit record: `audits/AUDIT.md`.

## Attribution

- **Mathematical conventions** follow, with equations restated from (see
  `docs/conventions.md`), L. Adams and S. Weinzierl, arXiv:1704.08895, and L. Adams,
  C. Bogner and S. Weinzierl, arXiv:1504.03255 (iterated integrals of modular forms, the
  all-orders sunrise); T. Ahmed, E. Chaubey, M. Kaur and S. Maggio, arXiv:2402.07311
  (ACKM, the non-planar two-loop family); and M. Becchetti, F. Coro, C. Nega, L. Tancredi
  and F. J. Wagner, arXiv:2502.00118 (BCNTW, the diphoton amplitudes and their ancillary
  data). The full attribution list is the "Attribution" section of README.md. The library
  implements their published mathematics; the code, the rigor architecture, and the
  certified bounds were written for BootLoops, with the exception stated next.
- **Code that follows FLINT**: in `src/siegel.jl`, the six transvectant calls of
  `igusa_clebsch` follow `acb_theta_g2_transvectants` in FLINT,
  `src/acb_theta/g2_covariants.c`, Copyright (C) 2023 Jean Kieffer, LGPL-3.0-or-later
  (same argument tuples, same intermediate covariants C_24, C_32, C_32^3, same choice of
  the degree-6 and degree-10 generators). Those lines are distributed under
  LGPL-3.0-or-later (`LICENSE-LGPL-3.0`, `NOTICE`); the notice is repeated at the
  function. The other FLINT contact in `siegel.jl` (`siegel_theta`, `siegel_theta_all`,
  `igusa_invariants`, `thomae_sextic`, `rosenhain_from_tau`) consists of one-call
  wrappers over FLINT's public functions through Arblib.jl.
- **Dependencies**: Arblib.jl over FLINT/Arb (LGPL), Nemo.jl. Obtain via the standard
  Julia package manager (`Pkg.instantiate()`).
- **Cross-check oracle**: `ginac-oracle/` is a small original C++ program that drives
  GiNaC 1.8.7's independent Eisenstein/modular-form machinery (Walden–Weinzierl,
  arXiv:2010.05271). GiNaC is GPL and is **not bundled**: install it from your
  distribution or from ginac.de, then compile with `ginac-oracle/build.sh`. The oracle
  is optional — it regenerates the independent numbers in `ginac-oracle/results.txt`.
  The oracle source is MIT like the rest of the package; a compiled `ginac_oracle`
  binary links GiNaC/CLN and is therefore a GPL-covered combined work, which is why
  only the source ships here.
- **Reference data**: the Zenodo ancillaries of arXiv:2502.00118
  (DOI 10.5281/zenodo.14733100, CC-BY-4.0) are used for validation. They are
  **not vendored** (the record is ~211 MB); `vendor/zenodo-2502.00118/fetch.sh`
  downloads them on demand, and the report there records exactly which entries are
  consumed and what was verified against them.

## Prior public code

Programs by other authors that cover parts of the same ground and predate this package:

- **GiNaC** (version 1.8 and later; M. Walden and S. Weinzierl, arXiv:2010.05271;
  GPL-2.0-or-later) evaluates the same iterated integrals of modular forms in floating
  point. It is the cross-check used here (`ginac-oracle/`). What Eichler.jl adds is a
  proved error bound on every value.
- **hcperiods** (P. Molin and C. Neurohr, arXiv:1707.07249; github.com/pascalmolin/hcperiods;
  MIT) computes period matrices of hyperelliptic curves in Arb with proved error bounds.
  The device of building a smooth branch of the square root on a segment from per-factor
  principal roots, used in the genus-2 period layer of `src/siegel.jl`, is the central
  device of their paper.
- **ore_algebra** (github.com/mkauers/ore_algebra; analytic continuation of D-finite
  functions with rigorous bounds by M. Mezzarobba; GPL-2.0-or-later) covers the task of
  `src/frobenius.jl` with sharper bounds and with regular singular points, which
  `frobenius.jl` does not handle.
- **BananaAmplitude.gp** (the PARI/GP script distributed by the Bonn group as
  supplementary material to arXiv:2008.10574 [BFKNS]) computes the banana amplitude from
  local Frobenius bases matched between neighboring singular points. The script carries no
  license statement.

## Layout

- **LICENSE**: MIT, Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz;
  code written by Claude (Anthropic) under his supervision. **NOTICE**,
  **LICENSE-LGPL-3.0** and **LICENSE-GPL-3.0**: the notice and license texts for the
  block of `src/siegel.jl` that follows FLINT.
- The audited independent sunrise anchor `anchor_sunrise_dps130_M240.json`
  (AMFlow-side, Cauchy order M=240, honest to ~128 digits) ships in `test/data/`;
  `test/runtests.jl` reads it from there (override with `EICHLER_AMFLOW_ANCHOR`).
- `ginac_oracle` writes `results.txt` in the current directory unless given a path;
  rebuild with `build.sh`.

## Candidate results to report upstream

Not applicable as patches (this package is not a fork). One numerical caveat, recorded in
`NUMERICS-RESULTS.md`: at Cauchy order M=120 the dps-130 anchor values are honest
only to ~72 digits (the ε⁻² entry deviates from exactly 3/2 at digit 73); the
M=240 run is good to its full ~128 digits.

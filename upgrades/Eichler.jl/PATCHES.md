# Eichler.jl — packaging and attribution notes

## What this is

Eichler.jl is an **original package**, not a fork or port of an external engine — there
is no upstream code base and hence no patch set against one. It is the
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
  certified bounds are original.
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

## Layout

- **LICENSE**: MIT, Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz;
  code written by Claude (Anthropic) under his supervision.
- The audited independent sunrise anchor `anchor_sunrise_dps130_M240.json`
  (AMFlow-side, Cauchy order M=240, honest to ~128 digits) ships in `test/data/`;
  `test/runtests.jl` reads it from there (override with `EICHLER_AMFLOW_ANCHOR`).
- `ginac_oracle` writes `results.txt` in the current directory unless given a path;
  rebuild with `build.sh`.

## Candidate results to report upstream

Not applicable as patches (original package). One numerical caveat, recorded in
`NUMERICS-RESULTS.md`: at Cauchy order M=120 the dps-130 anchor values are honest
only to ~72 digits (the ε⁻² entry deviates from exactly 3/2 at digit 73); the
M=240 run is good to its full ~128 digits.

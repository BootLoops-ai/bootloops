# Zenodo 10.5281/zenodo.14733100 — ancillaries of arXiv:2502.00118 (fetch on demand)

Reference data used by the Stage 2–3 validation in `NUMERICS-RESULTS.md` and by
`gauntlet/production2.jl`. **Not vendored** — the full record is ~211 MB, dominated by
one 207 MB helicity-coefficient file. Run `./fetch.sh` in this directory to download
the small files the library actually consumes (~620 KB); `./fetch.sh --all` fetches
the complete record.

- Record: <https://zenodo.org/records/14733100> — "Ancillary Files to 'Analytic
  two-loop amplitudes for q qbar → γγ and gg → γγ mediated by a heavy-quark loop'",
  ancillaries of arXiv:2502.00118 by M. Becchetti, F. Coro, C. Nega, L. Tancredi and
  F. J. Wagner (JHEP 06 (2025) 033).
- License: the record is licensed **CC BY 4.0** (per the record metadata). If you redistribute the data,
  carry the attribution above.
- What Eichler.jl uses and what was verified against it: `REPORT.md` here.

The two threshold constants a₀, a₁ that `production2.jl` checks against live in
`HelCoeffs2lbareExpTt.m` (the 207 MB file); their 68-digit printed values are quoted
in `REPORT.md`, so the certified cross-check can be reproduced without downloading it.

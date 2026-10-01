# eichler — guide (member homes for the Eichler.jl family)

Tool page: https://bootloops.ai/tools/eichler.html (tools index:
https://bootloops.ai/tools/).

KIND: package directory with three members; the certified engine itself,
Eichler.jl (Julia; periods, iterated integrals of modular forms, L-values,
Siegel theta constants as provable intervals), ships at `../../upgrades/Eichler.jl`
(override with the `EICHLER_PROJECT` env var). Members: `mirror6/` (Python +
gmpy2), `theta9/` (Python + mpmath, drives Julia in the Eichler.jl project),
`genus2/` (Python, sympy + mpmath; `import eichler.genus2` with `tools/` on
PYTHONPATH). Package battery: `selftest.py`.

PURPOSE: the pieces that sit around the engine when a Feynman integral, or any
D-finite period problem, is controlled by a curve of genus one or two or by a
Calabi-Yau operator:

- **mirror6** — exact (gmpy2 rationals) mirror-map / instanton-number
  fingerprint of a MUM-type point of a D-finite operator: the log-Frobenius
  chain, the canonical q-coordinate and its inverse, the structure-series
  ladder (Yukawa-type couplings), Lambert-kernel instanton-type numbers with
  script-emitted integrality verdicts, and a p-adic / growth fingerprint. Every
  generalization beyond the textbook CY3 case is a pinned definition in
  `mirror6/mirror6.py`; no geometric-mirror or constant-recognition claim is
  made by the engine.
- **theta9** — certified Siegel theta-constant observables on Eichler.jl's
  `acb_theta` layer: `theta9_stub.py` emits frame-covariant theta-ratio
  observables (from a period matrix tau, or from a period point in a marked
  frame via the frame law in `theta9_frame.py`), `theta9_copair.py` checks the
  engine against an independent lattice-sum reference and at two precisions,
  `theta9_scan.py` / `theta9_cert_residuals.py` post-process. Emission only: no
  integer-relation verdicts are produced here.
- **genus2** — genus-2 curve identification and recognition hygiene: Clebsch /
  Igusa-Clebsch invariants of a sextic and Mestre's reconstruction of the curve
  from its invariants over a number field (`mestre_port.py`, a sympy
  transliteration of SageMath's mestre and invariants modules (GPL-2.0-or-later;
  see `LICENSE-GPL-2.0` and `NOTICE`), checked exactly against Sage's published
  doctest values); recognition-free arbitration between candidate sextics by absolute
  invariants plus the Richelot orbit against a high-precision numeric oracle,
  with no PSLQ anywhere (`invariant_harness.py`); the PSLQ recognition law as
  reusable structure — planted and negative controls in-pipeline, two
  independent legs, pigeonhole floor per hit (`pslq_law_template.py`); an
  exact-rational-basepoint DE transport pattern on `tools/wayfinder`
  (`exact_transport_driver.py`); fast exact quadratic-field conic
  diagonalization and point search (`conic_fast.py`). Member guide:
  `genus2/GUIDE.md`.

USE-WHEN: you have an operator with a point of maximal unipotent monodromy and
want its exact q-expansion data and integrality pattern (mirror6); you need
theta constants of a genus-2 (or genus-4) period matrix as certified balls for
a downstream sealed scan (theta9); you must decide which exact genus-2 curve
sits behind a set of high-precision numbers, rebuild a curve from invariants,
or make a PSLQ recognition claim that survives controls (genus2).

NOT-FOR: the engine's own jobs (elliptic and modular iterated integrals,
certified DE transport, L-values) — call Eichler.jl directly; elliptic iterated
integrals from an exact curve in Python — that is `tools/ellipticus`.

INVOKE:
- mirror6: `python3 mirror6/run_control.py --outdir <scratch>` (mandatory
  positive control first; pilot then full, receipts in <scratch>);
  `python3 mirror6/run_l6.py --l6 L6_exact.json --outdir <same scratch>` for
  the order-6 application (the sha-pinned operator file is not included; the
  script refuses by name without it). `import mirror6` from `mirror6/` for the
  engine functions.
- theta9: `python3 theta9/theta9_stub.py --label L --tau "..." --outdir DIR`
  (entry A) or `--pi-json FILE:key --gram-json FILE:key` (entry B); needs Julia
  with the Eichler.jl project instantiated. `python3 theta9/theta9_copair.py
  <out_lo.json> <out_hi.json> <receipt.json>` for the two gates; reference
  values in `theta9/fixtures/`.
- genus2: `python3 genus2/invariant_harness.py` with `G2KIT_TRUE_JSON` (alias
  `EICHLER_GENUS2_TRUE_JSON`) pointing at your oracle json
  `{'leg1': {'rosenhain': [three ~370-digit complex strings]}}` — refuses
  loudly with that schema when unset; `python3 genus2/conic_fast.py` runs its
  worked example; `from eichler.genus2 import mestre_port` for invariants and
  Mestre reconstruction. `pslq_law_template.py` and `exact_transport_driver.py`
  are read-and-adapt templates: they import a period/theta pipeline and a
  DE-system class that are not distributed and refuse loudly without them.

LIMITS: mirror6 outputs are exact invariants under its stated definitions, not
mirror-symmetry claims; theta9 and mirror6 emit values and consistency gates
only — constant recognition belongs to a sealed PSLQ layer (`tools/lockpick`);
genus2's harness works at mp.dps ~340 and needs an oracle of comparable
precision; decimal (non-dyadic) basepoints poison downstream PSLQ, which is
why the transport driver insists on exact rationals.

BATTERY (partial): `python3 selftest.py --outdir <scratch>` from this directory
— M: the mirror6 control battery (28 boolean rows against published quintic
and Rodland instanton numbers, closed-form y0, mirror q2 = 770, self-duality,
two-route Yukawa; must end CONTROL-PASS); L6: `run_l6.py` verified to refuse by
name without `L6_exact.json`; G: the genus2 legs S1 (mestre_port vs Sage
doctest vectors, exact), S2 (conic_fast worked example, diagonalization
re-checked exactly), R1/R2 (fail-closed refusals of invariant_harness and
pslq_law_template); T: theta9 reported as a named SKIP (its reference receipts
are not part of this repository). About 2 s. `$EICHLER_WORK` is the default
scratch directory when `--outdir` is not given. The genus2 legs alone:
`python3 genus2/selftest.py`.

RELATED: upgrades/Eichler.jl (the engine); tools/ellipticus (elliptic iterated
integrals in Python); tools/wayfinder (DE transport core); tools/lockpick
(sealed PSLQ gate and ambient-precision law); tools/vopclose (closure tools
whose claims the genus2 PSLQ law governs).

CREDIT: D. van Straten (arXiv:1704.00164), Bönisch-Klemm-Sheykin-Zagier
(arXiv:2203.09426) and the AESZ tables (arXiv:math/0507430) for the published
control numbers mirror6 reproduces; FLINT's `acb_theta` (via Arblib.jl) under
theta9; J.-F. Mestre (1991) and K. Lauter, T. Yang (2011) for the genus-2
reconstruction and invariant conventions, with SageMath's published doctest
values as the exact cross-check for `mestre_port.py`.

## License

MIT License, except `genus2/mestre_port.py`, which is GPL-2.0-or-later (derived from SageMath's `hyperelliptic_curves.mestre` and `invariants` modules, Copyright (C) 2011, 2012, 2013 Florian Bouyer, Marco Streng, 2008 Nick Alexander, 2025 Sabrina Kunzweiler, Gareth Ma, Giacomo Pope; see `NOTICE` and `LICENSE-GPL-2.0` in this directory). Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).

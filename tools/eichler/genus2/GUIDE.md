# eichler.genus2 — member guide

Tool page: https://bootloops.ai/tools/eichler.html (tools index:
https://bootloops.ai/tools/). Member of the `eichler` package (`tools/eichler/`,
package guide one level up: `../GUIDE.md`); formerly the standalone `genus2kit`.

KIND: package member at `tools/eichler/genus2/` (pure sympy/mpmath — no Sage
dependency): mestre_port.py, invariant_harness.py, pslq_law_template.py,
exact_transport_driver.py, conic_fast.py; `import eichler.genus2` (with `tools/`
on PYTHONPATH) exposes `mestre_port`, the other members are run as files.

PURPOSE: genus-2 curve identification + recognition hygiene kit: Igusa-Clebsch invariants + Mestre reconstruction over number
fields, recognition-free candidate arbitration, the full PSLQ recognition law as
reusable structure, exact-rational-basepoint DE transport, fast exact conic
diagonalization/point-search.

USE-WHEN:
- Arbitrating "is this exact sextic THE curve behind these numerics" —
  invariant_harness.py compares absolute invariants (+ Richelot orbit) against a
  high-dps numeric oracle, NO PSLQ anywhere; kills wrong identifications in minutes.
- Reconstructing a genus-2 curve from Igusa-Clebsch invariants over a number field
  (mestre_port.py, a sympy translation of the SageMath code by Bouyer, Streng and
  Alexander (GPL-2.0-or-later; see ../LICENSE-GPL-2.0 and ../NOTICE), checked exactly
  against SageMath's published doctest values).
- Any PSLQ recognition claim: pslq_law_template.py = planted+negative controls
  in-pipeline, two independent legs with relation matching, pigeonhole-floor line per
  hit — zero exceptions.
- DE transport whose basepoint choice poisons PSLQ (exact_transport_driver.py:
  exact-rational basepoints, per-column parallelism — a decimal basepoint
  poisons PSLQ; dyadic basepoints do not).
- Conic point-search that sympy-nsimplify grinds on (conic_fast.py: exact
  quadratic-field pair arithmetic; 28 min -> seconds on the originating conic).

NOT-FOR: genus-2 Siegel-modular PERIOD/theta machinery (that is Eichler.jl — ships in
this repo's upgrades/Eichler.jl — and the sibling member `theta9`); no CLI —
these are import-and-drive modules with self-tests inline.

INVOKE (public): `python3 invariant_harness.py` (self-test + planted-mismatch gates
built in; needs env `G2KIT_TRUE_JSON` — alias `EICHLER_GENUS2_TRUE_JSON` —
pointing at your oracle json; the reference oracle json is not included and the
harness refuses loudly with the required schema when unset). `python3 conic_fast.py`
runs its built-in worked example (writes `conic_point.json` or
`conic_obstruction.json` into the current directory). `from eichler.genus2 import
mestre_port` (tools/ on PYTHONPATH) or import mestre_port / pslq_law_template /
exact_transport_driver as flat modules from this directory.
REFERENCE-ONLY members (templates): pslq_law_template.py ships as the law's structure — it imports `invariants_pipeline` (a period/theta pipeline that is not distributed) and refuses loudly without it; adapt the
skeleton to your own pipeline. exact_transport_driver.py imports
`monodromy_transport.L5System` (a DE-system spec, not shipped) — swap in
your own system class; Wayfinder ships in this repo's tools/ and the driver
finds it three directory levels up.

INPUTS: exact sextic coefficients / Igusa-Clebsch invariants / high-dps numeric oracle
values (mp.dps ~340 in the harness).

OUTPUTS: arbitration verdicts (invariant match/mismatch per candidate + Richelot
orbit); reconstructed curve models; transported DE values.

GATES: harness has built-in self-test + planted-mismatch controls — run them before
trusting an arbitration; PSLQ template's planted/negative controls are mandatory, not
optional.

FOOTGUNS:
- Decimal (non-dyadic) basepoints poison downstream PSLQ — always transport from
  exact rationals (the driver exists because of this).

BATTERY (smoke): `python3 selftest.py` in this directory (the package command
`python3 selftest.py --outdir <scratch>` from `tools/eichler/` runs the same legs
as its G leg) — S1 mestre_port invariants + Mestre
conic vs the Sage doctest vectors (exact sympy equality); S2 conic_fast worked
example in a scratch directory with the congruence diagonalization re-checked
exactly; R1/R2 invariant_harness and pslq_law_template refuse loudly with named
requirements when their external inputs are absent (fail-closed refusals
verified). S2 scratch goes under `$EICHLER_WORK` when set. Public-runnable:
mestre_port, conic_fast, invariant_harness with your own oracle json.
Reference-only: the pslq template's IP.* pipeline calls; the transport driver's
L5System spec.

RELATED: upgrades/Eichler.jl and the sibling member `theta9` (Siegel/theta side of
genus 2); vopclose/PSLQ closure tools (the law template governs their claims);
lockpick (ambient-dps law applies to the mpmath legs); tools/wayfinder (the
transport core the driver drives).

CREDIT: the reconstruction is J.-F. Mestre's ('Construction de courbes de genre 2 à
partir de leurs modules', in Effective Methods in Algebraic Geometry, eds. T. Mora & C.
Traverso, Progr. Math. 94, Birkhäuser 1991, pp. 313-334) from the invariants of Clebsch
and Igusa (J. Igusa 1960, Ann. Math. 72:612), in the conventions of K. Lauter & T. Yang
(2011, J. Number Theory 131:936). mestre_port.py is a sympy translation of the
SageMath code by Bouyer, Streng and Alexander (GPL-2.0-or-later; see ../LICENSE-GPL-2.0
and ../NOTICE): the `hyperelliptic_curves.mestre` module (Florian Bouyer and Marco Streng)
and the `invariants` module (Nick Alexander), with their later Sage contributors (Sabrina
Kunzweiler, Gareth Ma, Giacomo Pope). It keeps their conventions and normalizations, its
exact self-checks are those modules' published doctest values, and only the transvectant
routine is written independently. We thank the Sage authors.

## License

MIT License, except `mestre_port.py`, which is GPL-2.0-or-later (a sympy translation of SageMath code by Bouyer, Streng and Alexander; see `../NOTICE` and `../LICENSE-GPL-2.0`). Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).

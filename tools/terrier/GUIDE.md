# terrier — GUIDE

Tool page: https://bootloops.ai/tools/terrier.html

A package: one suite with three wings — `periods/`, `lattice/`, `census/`
(+ the shared chassis `common/`) — plus the suite gate `selftest.py` and the
sha-pin manifest `regression_manifest.json`. README.md is the overview; the module
index is MODULE_MAP.md. The package ships a compact set of verification fixtures
(toy/synthetic control data + known-answer tables + small published samples) that the
sha-pinned selftest verifies — do NOT edit them (the sha-pins would break). Full-scale
reference data sets are not included; batteries that need them print a named SKIP
naming the `TERRIER_*` environment variable that supplies them.

## PURPOSE

Unified string-vacua/F-theory suite — certified CY period point-values
(cards/flux/transport/summation), flux-lattice enumeration + Nikulin certificates, and
the census wing (fold/dedup/adjudicate/dist_cert/controls/verdicts).

## USE-WHEN

- Certified period values: card contract (periods/cards.py), Γ-series (geffseries),
  fluxcurves ξ-toggle, MUM+conifold stepped transport with PROVEN tail majorant
  (envelope_certified), Krawczyk/ball certifier, Pfaffian matrix transport
  (route_choice front door, scalar-elimination refusal law)
- Certified multi-param (6-var) period point-values WITHOUT transport marching: u1s
  summation card B1-B8 (periods/summation/) — 84 ops x 15,625 pts termwise identities
  PASS; entropy radius law s<1 in its README
- Finite-field attractor fingerprints: periods/ffp.py all-z Frobenius scan; refined
  S-integrality candidates: periods/sunit.py; known-point pullback-orbit dedup BEFORE
  escalation compute: common/dedup.py (newform match = REDISCOVERY never NEW)
- Lattice: SL2 definite-BQF reduction + class enumeration (covers the Hecke d<0
  NotImplemented gap), nikulin.jl disc-form route at sig(0,3) (avoids the sig(3,0)
  primitive_embeddings hang), emit/ + sweepchassis; genus enumeration -> genus_enum.jl
  (Aut-free Kneser + PARI qfauto mass, avoids the Hecke OOM); Venkov LP kill at any
  rank -> lp_kill_rank_agnostic.py (exact Fraction Farkas certs; ranks 3/4/6/7)
- Census: fold_hnf (HNF+detU complete orbit key, witnessed merges), dedup two-stack RED
  law + witness/reconcile (reconciled count AUDIT-ONLY, never consumed), isd_cutoff
  (cutoff formula fixed BEFORE any count; derivation shipped), dist_cert, control_sampler,
  tile_engine (fail-closed OL-2/4/5)
- HV4 diagonal PF card: cards/hv4-diag-L5.json (order-5 DERIVED-UNIQUE == printed JKK
  operator)

## NOT-FOR

dist_cert / control_sampler / tile_engine are APPLICATION-SPECIFIC (criterion constants tuned to their original application; do not reuse without an independent validation receipt). sunit truncation = bounded
net, NOT a finiteness proof. ffp: p >= 2^31 REFUSED (BigPrimeRefusal) unless
allow_big_prime; validated envelope HV p<=127, BCM p<=337. nikulin cert covers
glue-level uniqueness ONLY (residual multiplicity = complement-class multiplicity).
No scalar GKZ-restriction module: point-value targets use pfaffian.py (a hand-eliminated
scalar operator can blow up); conifold leg-1 scalar transport
FAIL-CLOSED by design. GV-band labeling law:
completeness certified only through margin level; beyond it the x5 band is LABELED,
never a certified radius. Not Feynman fixed-eps DE transport (use the DE-transport
tool, Wayfinder).

## INVOKE

python facades periods/{cards,geffseries,fluxcurves,transport,verdicts,opderive,
pfaffian,conifold,gvderive,ffp,sunit}.py; common/{receipts,controls,frames,certs,
verdict,dedup}.py; lattice/{lattice,validate,nikulin}.jl + emit/;
census/{fold_hnf,dedup,dist_cert,control_sampler,isd_cutoff}; suite gate = selftest.py
(manifest pins FIRST, then battery replays, then coverage gates); rebuild pins via
build_manifest.py ONLY after a full green run.

INPUTS: toric/GKZ data + card schema (fail-closed PENDING-CARD), flux vectors,
lattice Gram data, census entry streams; full-scale reference data sets are not included.
OUTPUTS: certified balls (two-dps), sha-stamped atomic receipts, exact annihilation
certificates, Fraction Farkas certificates, census counts + witness files,
machine-readable truncation statements.
ENV: one `TERRIER_*` family (full table in README.md, "Environment"). Reference-data
paths have no defaults — set the variable to your own data or the leg is a named
SKIP / loud refusal; julia legs expect an Oscar/Hecke project env read from
`TERRIER_JULIA_PROJECT` (no machine-local default — unset, the julia legs are named
SKIPs).

## GATES

selftest.py green (59 sha-pins, 30 registered batteries; sabotage law: a corrupted
pin -> NAMED failure, replays REFUSED rc=1); per-module batteries check real results
on the shipped fixtures; mutation/planted-error controls; dedup runs BEFORE any
escalation compute.

## PITFALLS

- The modules under `periods/pipeline/` and the fixture data are sha-pinned — edits
  break the manifest; change them deliberately, re-run the battery, then rebuild the
  pins (see MODULE_MAP.md)
- Scalar elimination blowup: use pfaffian route_choice front door, never hand-eliminate
  to a scalar operator without the mod-p price probe
- int64 mod-p overflow: ffp refuses big primes by design — do not bypass
  allow_big_prime without cause
- Two-stack dedup disagreement: file RED + witness the misses; the reconciled count is
  AUDIT-ONLY

## What runs out of the box (battery: PARTIAL — the rest needs external data)

- Canonical battery, green from a cold checkout (~40 s wall):
  `python3 selftest.py`
  — M.manifest sha-pins PASS + the shipped-data battery replays PASS + X.coverage
  3/3 PASS; every leg whose inputs are not included prints a named SKIP row saying
  exactly which env var to set (never a silent skip). Julia legs need an Oscar/Hecke
  project env (`TERRIER_JULIA_PROJECT`).
- The shipped fixtures are toy/known-answer data the batteries check real
  results against: a synthetic S-integrality control card with hand-derived
  acceptance values, a 44-row exact glue table verified in pure Fraction
  arithmetic, closed-form survivor pins, regenerated pass-sets, and small
  published samples. Point the named env vars at full-scale reference data and
  those legs run for real.

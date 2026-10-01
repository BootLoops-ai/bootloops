# CLINCH — certified local interval-Newton convergence for hierarchies

Tool page: https://www.bootloops.ai/tools/clinch.html

KIND: package (`clinch/` + `battery/` + `adapters/` + `census/` +
`MANUAL.md` — the manual has the full detail; this guide is the summary).

PURPOSE: Certified Local Interval-Newton Convergence for Hierarchies — a fitted MAP
optimum theta* either gets a CERTIFICATE (inside the box: a true stationary point of
the penalized objective EXISTS, is UNIQUE, Hessian PD throughout — a certified
strict local minimum) or a REFUSAL with the failed contraction NAMED. Verdict set:
CERTIFIED-INTERIOR / CERTIFIED-CORNER (KKT adjudicated, strict complementarity
certified) / REFUSED(named) / OUT-OF-SCOPE. Scope: the HIERARCHICAL/BLOCK
axis — block-arrow Krawczyk + interval Schur at latent-group scale (5592-9
coordinate certificates in ~90 s measured).

USE-WHEN:

- A hierarchical/mixed/random-effects MAP fit's "the optimizer converged" needs to
  be a THEOREM (or an honest named refusal) — SEs/LRTs/model choice load-bearing on
  it.
- You need the measured Newton-metric polish quality of a fitted optimum (raw gmax
  lies: a 9.6e-5 raw gradient can be 0.029 curvature-sigma, 4 orders too loose to
  certify).

NOT-FOR: global optimality (LOCAL statements only); non-smooth objectives; model
comparison certificates; unported model spaces (check_space refuses cold ->
OUT-OF-SCOPE).

KERNEL: consumes `baller.certify.block_krawczyk` BY IDENTITY (never copied) — the
BALLER package ships beside this one at `tools/baller/`. `clinch/engine.py` resolves `BALLER_DIR` to the `baller` package beside this
one (`CLINCH_BALLER_DIR` overrides).

INVOKE: `sys.path.insert(0, ".../tools/baller"); sys.path.insert(0,
".../tools/clinch"); import clinch` (verify() fail-closed source pins, pyc purged);
`from clinch import adapter_v31 as A, engine, cert; engine.certify(o, theta[,
bounds=(lo,hi)])`; `engine.newton_polish` + certify = the labelled polished-center
product; `cert.write_cert` -> CERT.json + human statement. Battery: `python3
battery/battery.py` from scratch space — it refuses to run from the tool tree
or other protected working directories (fail-closed refusal, ships as-is).

## Reference data (not shipped) and what runs without it

The model oracles (`oracle_v31`, `oracle_v36`) target a specific hierarchical
demographic model; the labels v3.1/v3.6/v3.7 name its layouts (model spaces,
not versions of this package). The reference model code, the reference fits
(labeled fit13, fit17, fit16, fit18a), their certificate receipts and the
`adapters/etienne/` study data are not shipped. Everything that needs them is
gated on environment variables (`CLINCH_CODE_DIR` for the reference model
code, `CLINCH_REFERENCE_DIR` for the reference fits and their stored
results; the members' gates are listed under MEMBERS) and
refuses or skips by name when they are unset.

- Self-contained (run from this copy with `tools/baller` beside it): battery
  legs L1 (planted optimum + corner + out-of-scope), L2 (planted saddle PD
  refusal), L5 (dps ladder honesty), L7 (engine exhaustion typed refusal),
  L8 (centered-space reduced-Hessian PD — the rank-4 congruence extension,
  clinch.kkt_pd, proven on a toy centered KKT fixture with its
  signature/Cholesky refusal controls), and L4's synthetic tamper/mutation
  controls.
- Reference-data legs: L3 (flat-valley control), L6 (flagship receipts) and
  L4's real-oracle tamper sub-leg, with their producers `battery/flagship.py`,
  `battery/l3_run.py`, `battery/l3_gates.py` and the reference loaders in
  `clinch/adapter_v31.py`.

Battery: run from a scratch cwd, repo geometry with tools/baller beside it —
OVERALL PASS, 8/8 legs, ~6 s with the reference data present. Without it:
L1/L2/L5/L7/L8 PASS, L3/L4/L6 fail on the absent receipts by design. One
command: `python3 selftest.py` (from this directory) runs the battery from a
scratch working directory and grades exactly this contract — exit 0 iff the
self-contained legs are green and any receipt-leg failures are the
absent-receipt class.

FOOTGUNS: pure curvature-scaled box ONLY (raw-uniform: Bauer-Skeel ~cond(H);
step-floor mixing: radius-ratio amplification — both measured and documented in
engine.py); oracle must be LOG-SPACE over wide boxes (expit/ball-products lose
positivity); arb ball digamma NaNs at rad~0.85 (use jets' monotone endpoint hulls);
fat-ball law: radii capped at 1 sigma — weak-but-true certificates are REFUSED; the
1e-12 identity bar rides the matched-precision adjudication register (float64
reference noise is ~5e-10, receipted).

MEMBERS: two absorbed external certified-optimum instruments —
`adapters/certpass/` (point-process certificate pass, CERTPASS_REFERENCE_DIR-gated
runners) and `census/repertoire_g1/`
(Gate-1 complete-census certificates, $G1_OUT receipts); contracts and env
gates in MANUAL.md MEMBERS + each member's PROVENANCE.md.

RELATED: baller (kernel).

CREDIT: the certification kernel is the Krawczyk operator (R. Krawczyk 1969, Computing
4:187) with Moore's interior test (R. E. Moore 1977, SIAM J. Numer. Anal. 14:611) in the
verification semantics of S. M. Rump (2010, Acta Numerica 19:287), assembled blockwise
in baller.certify.block_krawczyk over Arb (F. Johansson 2017) via python-flint.
adapters/certpass certifies optima of the ETAS model (Y. Ogata 1988, JASA 83:9) as
fitted by the open `etas` package of L. Mizrahi, S. Nandan & S. Wiemer (2021, Seismol.
Res. Lett. 92:2333; github.com/lmizrahi/etas, MIT), whose EM follows Veen & Schoenberg
(2008, JASA 103:614) — we thank its authors for a codebase clean enough to certify
against; adapters/etienne certifies optima of Etienne's (2005, 2007) neutral sampling
formula. The hierarchical reference model and its fit code are an in-house study,
consumed read-only.

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).

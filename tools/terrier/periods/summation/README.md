# periods/summation/ — SUMMATION-ROUTE period evaluation
The worked card is the 6-parameter Hulek-Verrill fourfold family HV4 (see
periods/pipeline/cards/hv4-diag-L5_NOTE.md; fundamental-period coefficients
c(k) = multinom(|k|;k)^2), items B1-B8 below. boxops_u1.py imports its
generic box-operator engine from an external bank (TERRIER_BOXOPS_BANK, not
included in the package); everything else runs on the packaged code.

## Transport route vs summation route — when each applies
- TRANSPORT (periods/transport.py, pfaffian.py): march v'=A(s)v along a path
  from a MUM point; needs an annihilating operator/connection and a stepped
  path with proved tail majorants. Route of choice for 1-param (h^{2,1}=1)
  targets and for point values far from the large-volume patch.
- SUMMATION (this dir): evaluate Pi and its jets DIRECTLY as the multivariate
  Gamma/GF series at the target point — no operator derivation, no marching.
  Route of choice when (a) the moduli count makes scalar operators explode
  (scalar-elimination blowup; 6-param HV4 has NO practical L_s) and
  (b) the target sits inside the proven convergence polydisc.
- ENTROPY RADIUS LAW (the applicability gate, EXACT Cauchy-Schwarz on
  multinomials): sum_{|k|=m} c(k)|phi|^k <= s^{2m} with s = sum_I sqrt|phi_I|;
  proven convergence iff s < 1, per-order ratio s^2. If s^2 >= 0.9 anywhere on
  the box, stage-1 summation is DEAD there (K4 law) — fall back to transport
  or refuse. b3_entropy.py certifies s per point/box hull, flint prec 200,
  UNDECIDED => FAIL.

## Modules
- boxops_u1.py  — B1 box-operator supply on the orthant (full-size gate:
  84 ops x 15,625 lattice pts = 1,312,500 termwise identities, ~250 s).
- towers6.py    — B3 rho-jet Gamma-series tower evaluator (jet depth <= 4,
  210 monomials; exact Fraction/CQ + arb/acb rings; GF form == lattice sum).
- b3_entropy.py — entropy certificates + J_d(51) majorant constants + closed-
  form tail table (M=50). b3_gates/b3_oracle/b3_pass — B3 gate battery
  (full-size gate: GF==lattice 210/210).
- b7_checks.py / b7_numbers.py — B7 Cauchy floors (Eulerian closed forms,
  tail corollary, log-modulus inequality; method in B7_CAUCHY.md). The
  shipped battery (periods/selftest_summation_reduced_u1s.py) requires all
  4 b7_checks PASS lines.
Downstream consumer: census/tile_engine (tile criterion evaluates this card).
The full-size B1/B3 receipts are not included in the package.

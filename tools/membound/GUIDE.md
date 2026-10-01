# membound — guide

Tool page: https://bootloops.ai/tools/membound.html

KIND: package (pure python + mpmath)

PURPOSE: soft-region (omega->0) boundary constants of worldline /
post-Minkowskian-type expansions, computed as frequency-space Bessel-kernel
integrals: the eps-expansion, layer by layer, of any two- or three-frequency
"omega-core" built from causal frequency factors (0+ - i s w)^{p - q eps},
Bessel lines |w|^{-eps} K_eps(|w|) and monomials, with the last frequency the
sum of the others. You describe your core in a small JSON file and the
package returns its Laurent coefficients J^(0), J^(1), J^(2); families of
cores assemble into boundary VECTORS in the |omega|^{-k eps} layer basis.
Frequency factors carry retarded or advanced i0+ ('ret'/'adv') or a
phase-only substitution flag ('fey'; a diagnostic, NOT the Feynman
propagator prescription, see METHOD). The worked example shipped and
validated in-tree is the post-Minkowskian memory family: the package
generalizes the single-constant c_M computation of the 5PM-2SF memory sector
(arXiv:2601.16256) and reproduces c_M = 1 as its gate.

USE-WHEN:
- You need the soft (omega->0) boundary constants of a soft-region
  (memory-type) frequency integral J(eps) = Int dw1 dw2 [dw3] of retarded or advanced
  factors (0+ - i s w)^{p - q eps} (s = +1 retarded, s = -1 advanced) times
  Bessel-K lines |w|^{-eps} K_eps(|w|)
  times monomials, with the coupler frequency w3 = w1 + w2 (two-frequency
  cores) or w4 = w1 + w2 + w3 (three-frequency cores, n_freq=3).
- The eps-expansion must happen BEFORE integration (fixed-eps evaluation of
  these cores suffers catastrophic cancellation between the eps^-4..eps^-2
  layers) and the order-eps phases of the frequency factors must be derived
  per sign sector, not
  hard-coded.

NOT-FOR: registry slots without constants (PX2/memory and NP/RR ship as slot
definitions with status flags; their constants are not provided in this
release; do not read absent slots as zeros); the odd-in-v dissipative
projection proj='im' is exposed but NOT validated against any published
value; not a general loop-integral evaluator; and not an evaluator of the
memory integrals with Feynman propagators on the memory gravitons: the 'fey'
flag is a phase-only substitution (see METHOD and FOOTGUNS), and the true
Feynman assignment lies outside what this package computes.

INVOKE:
- `python3 -m membound --spec core.json [--dps 11] [--maxdegree 4]
  [--panels 0,0.5,2,8,inf] [--jobs N] [--order-max 2] [--json-out out.json]`
  — YOUR core: loads the JSON spec (schema under INPUTS), farms every
  unique (eps-order, sign-domain class) quadrature over N processes
  (default min(12, cpu); `--jobs 1` is the serial library path), prints the
  Laurent coefficients J^(0..order-max) of J(eps) with the (2pi)^-n measure
  included, and with `--json-out` writes a `membound.core_moments.v1`
  record. Run it from the directory holding the `membound/` package (e.g.
  `tools/`), or equivalently `python3 cli.py --spec ...` inside the package.
  Start small (`--dps 9 --maxdegree 3 --order-max 0`) to size the cost;
  it grows steeply with dps and maxdegree. Worked example:
  `python3 -m membound --spec membound/examples/core_example.json --dps 9
  --maxdegree 3 --panels 0,2,inf` gives J^(0) = -1/15 to 8+ digits in
  seconds.
- `python3 selftest.py` — acceptance battery (~20 s); `--full` adds the
  complete c_M gate (minutes).
- `python3 gate_cM.py --micro` — cost-rate probe (ONE representative quad +
  projected total); run FIRST to size a full gate run on your machine.
- `python3 gate_cM.py [--dps N --maxdegree M --jobs J --json-out bv.json]` —
  the full validation gate on the worked example: recomputes the memory-sector masters
  I1^(M), I3^(M) through order eps^2, assembles the Gamma-prefactor Laurent
  series, applies the IBP relation for I2^(M), and checks j1^(0) = 1/30,
  the I1 leading pole = 1/(15(8pi)^4), the eps^-4/eps^-3 spurious-pole
  cancellation, c_M = 1, and a PSLQ recognition of c_M. Digit bars scale
  with --dps (need > dps-5).
- Library: `membound.core.OmegaCoreSpec` (+ `from_dict`/`to_dict` for the
  JSON schema) + `compute_core` (engine; two-frequency by default,
  `n_freq=3` for three-frequency octant cores), `membound.farm.
  compute_cores_parallel` (process farm over several cores at once),
  `membound.cli.load_spec`, `membound.spec.REGISTRY` (family/region slots
  of the worked example), `membound.assemble` (Mellin closed forms,
  prefactors, IBP, the `boundary_vector` JSON contract),
  `membound.laurent.L` (truncated Laurent algebra).

INPUTS: a core spec, either as an `OmegaCoreSpec` in Python or as a JSON
file for the command line (schema `membound.core_spec.v1`; see
`examples/core_example.json`):

    {
      "name":   "mycore",                 required; free text
      "n_freq": 2,                        2 (default) or 3 independent
                                          frequencies; index n_freq+1 is
                                          the coupler w3 = w1 + w2
                                          (w4 = w1 + w2 + w3 for 3)
      "ret":    [[idx, p, q, flag], ...], required (may be []); one entry
                                          per causal factor
                                          (0+ - i s w_idx)^(p - q eps):
                                          integer power p, eps-weight q,
                                          flag "ret" | "adv" | "fey"
      "klines": [idx, ...],               required; indices carrying
                                          |w_idx|^{-eps} K_eps(|w_idx|)
      "poly":   {"idx": n, ...},          optional; extra integer powers
                                          w_idx^n (string or int keys, or
                                          a list of [idx, n] pairs)
      "proj":   "re" | "im"               optional, default "re"
                                          (real-even projection; "im" is
                                          the odd slot, exposed but not
                                          validated)
    }

  Every index lies in 1..n_freq+1; unknown keys are refused (a "_comment"
  key is allowed). The integral computed is
  J(eps) = (2pi)^{-n_freq} Int d^{n_freq} w [prod ret][prod klines][prod
  poly], expanded in eps before integration; its layer index is
  k = sum q + #klines. Quadrature controls: `--dps` (working digits),
  `--panels` (magnitude-axis breakpoints, tanh-sinh clusters nodes at each;
  put one near every scale of your integrand), `--maxdegree` (tanh-sinh
  depth per axis), `--order-max` (0, 1 or 2). The integrand must decay at
  large |w| on every axis (a K line on each independent direction, or
  enough negative powers) — the package does not check convergence for you.

OUTPUTS: eps-layer moments J^(0..2) (mpf/mpc, measure included), printed as
the Laurent coefficients of J(eps) and, from the command line, written as a
`membound.core_moments.v1` JSON record {contract, spec, layer_k, J:{"0",
"1", "2"}, meta:{dps, maxdegree, panels, jobs, wall_s, digits, notes}} whose
digits field is null ("unknown - measuring") until you difference two runs;
for the worked example, assembled Laurent series and a JSON boundary-vector
record (contract `membound.boundary_vector.v1`) whose slots are
gamma-closed, numeric, or pslq-slot, each carrying measured digits or an
explicit "unknown - measuring" marker.

METHOD:
- Expand the integrand analytically in eps to O(eps^2) before integrating.
  Each frequency factor (0+ - i s w)^{p - q eps} is split into its leading
  (eps^0) part, built with the retarded sign structure (-i)^p w^p for every
  flag, and its order-eps logarithm log|w| + i phi(sign w), whose phase is
  set per factor by the flag: phi_ret = -(pi/2)s, phi_adv = +(pi/2)s,
  phi_fey = -pi/2. Per sign sector the layer weights (real log layer A,
  phase layer B = -sum q_i phi_i, Bessel nu^2 layer R) are derived
  automatically; this reproduces the hard-coded weights of the
  single-constant prototype as special cases.
- What the flags mean. 'ret' on every factor is the retarded assignment of
  arXiv:2601.16256 and is exact. 'adv' on every factor is exact as well:
  the map (w1, w2) -> (-w1, -w2) takes each retarded factor to the advanced
  one and leaves the integrals unchanged, so all-advanced equals
  all-retarded (and their average, the gamma-3 prescription) order by
  order in eps. 'fey' is a phase-only substitution: it replaces the
  retarded order-eps phase by the constant Feynman phase -pi/2 while
  keeping the retarded leading factor. It is a diagnostic of the
  assembly's sensitivity to the phase layer B, NOT the Feynman propagator
  prescription. A Feynman propagator changes the whole power,
  (0+ - i s w)^{p - q eps} -> (0+ - i|w|)^{p - q eps}, so the leading factor
  -w1 w2 of I1^(M) becomes -|w1 w2|; the leading cores then change (the I1
  pole becomes -1.3465... times 1/(15(8pi)^4)) and the eps^-4 pole of
  I2^(M) does not cancel (coefficient 1/(6144 pi^6)), so no c_M is defined.
  That assignment lies outside what this package computes.
- The sign-sector domains (6 for two frequencies) are grouped into
  equal-integrand classes by 3-point probing at (dps-3)-digit tolerance
  before quadrature; each class is integrated once and weighted by
  multiplicity.
- Three-frequency cores (OmegaCoreSpec(..., n_freq=3)) integrate over the 8
  sign octants of (w1, w2, w3): uniform octants map directly; each mixed
  octant splits at the |w4| = 0 locus into two pieces with the coupler
  magnitude as an explicit coordinate (Jacobian x+y on the
  simplex-parametrized piece) — 14 domains, classed the same way, with the
  guard that pieces integrating their third axis over (0, 1) never merge
  with magnitude-panel domains.
- Single-frequency slots close exactly via the Mellin law
  Int w^{s-1} K_eps(w) dw = 2^{s-2} Gamma((s+eps)/2) Gamma((s-eps)/2).
- kappa2 = (1/2) d^2/dnu^2 K_nu|_0 by 4th-order even finite differences
  under a dps-scaled working-precision guard (workdps ~ 1.5(dps+8)); no
  fixed precision constants (see the dps_lint member of tools/baller).
- Layer count k = sum(retarded eps-weights) + #K-lines (= 7 for the planar
  memory sector, matching the 2^{7eps} prefactors of arXiv:2601.16256).

GATES: acceptance battery in-tree (13 legs: Laurent algebra, weight-engine
vs hand-derived B1, per-flag order-eps phase table (ret/adv/fey), domain classing, Mellin
closed form incl. the kappa2 layer, order-0 moments vs 1/30 exact and
-8/105 regression, prefactor pole assembly vs 1/(15(8pi)^4), output
contract, the 3-frequency octant machinery: domain partition/classing,
exact coordinate + Jacobian + coupler-sign identities, coupler phase flags
vs hand formula, and the order-0 engine vs the exact 3/8 moment of
K0 K0 K0 w4^2; and the user-core front door: `python3 -m membound --spec
examples/core_example.json` through the process farm vs the exact
J^(0) = -1/15 of that non-registry core (8+ digits at dps 9) and vs the
serial library call at orders 0 and 1). Validation run: gate_cM.py at dps 30 / maxdegree 6 passed with the I1 pole
at 28+ digits vs the closed form, eps^-4/eps^-3 cancellation < 1e-27
relative, and c_M = 1 at 25+ digits (~14 min wall farmed over 12
processes; run logs not shipped in this repo). Diagnostic from the same
run: the phase-only substitution ('fey' on both memory factors) shifts c_M
from 1 to 1/5 exactly (Delta c_M = -4/5: the pi^2 parts of the
single-quadrant moments cancel and the rational parts survive), which shows
the assembly is sensitive to the order-eps phase layer. It is not a
Feynman-prescription value: with Feynman propagators on the memory
gravitons the leading cores change and the eps^-4 pole of I2^(M) does not
cancel (coefficient 1/(6144 pi^6)), so no c_M is defined and that case lies
outside what this package computes.

FOOTGUNS:
- A user core's digits are not certified by the package: the record's
  digits field is null until you rerun at higher --dps / --maxdegree (or
  finer --panels) and difference the two runs. An integrand that does not
  decay on some axis (no K line and non-negative net power along an
  independent direction) makes tanh-sinh return a finite wrong number, not
  an error.
- mp.pslq has a 53-bit working-precision floor: at low dps, run it under
  workdps >= 16 with the tolerance set from the DATA precision (gate_cM.py
  does this); never let the pslq workdps inflate the claimed digits.
- c_M through this quadrature route caps near 25 digits (tail-square decay
  of the order-2 I3 layers); past that the analytic route to the kappa2
  moments is required, not more maxdegree.
- The 'fey' flag is NOT the Feynman propagator prescription. It changes
  only the order-eps phase of a factor and keeps the retarded leading
  factor, so under it the eps^-4/eps^-3 poles of I2^(M) still cancel and
  only the finite constant moves (c_M -> 1/5). Do not read that as "the
  poles are prescription-blind": with true Feynman propagators,
  (0+ - i|w|)^{p - q eps} on both memory lines, the leading cores change
  and the eps^-4 pole of I2^(M) survives (coefficient 1/(6144 pi^6)), as
  arXiv:2601.16256 also reports for the Feynman prescription. The package
  does not implement that assignment.
- Mixed flags (one factor 'ret', one 'adv') are phase-only too: for integer
  p an advanced leading factor is (-1)^p times the retarded one, and the
  fixed leading phase (-i)^{sum p} (OmegaCoreSpec.base_phase) omits that
  sign, so a mixed assignment comes out with the wrong overall sign.
  Uniform 'ret' or uniform 'adv' is exact.
- The registry ships slot definitions with status flags; PX2/memory and
  NP/RR carry no constants in this release. Do not read absent slots as
  zeros.

RELATED: tools/lockpick (PSLQ closure for pslq-slots — never raw mp.pslq);
tools/baller's hygiene member dps_lint (`python3 -m baller.hygiene lint`;
the precision-guard rule the FD layer follows);
tools/nestor (certified quadrature front door, if a slot needs a receipted
oracle value).

CREDIT: the memory-region integrals, their IBP relation, and the exact
prefactors are from Driesse, Jakobsen, Mogull, Nega, Plefka, Sauer,
Usovitsch, arXiv:2601.16256 (we read the exponent in their eq. (27) as
−ε, which reproduces their c_M = 1). The expand-before-integrate engine, automatic
sector-weight derivation, per-factor phase flags, and the boundary-vector
contract are this program's.

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).

# PFAFFIAN ROUTE — DESIGN (certified W0 via first-order matrix transport)

This note explains why the period vector is transported as a first-order
MATRIX system rather than through an eliminated scalar operator, fixes the
shared data contract, and gives the gate tables and the cost model for the
two worked geometries: ads-5-81-3213 (route A, MUM anchor s = 0 to the flux
vacuum s_vac) and the conifold-curve cards (route B; the example quoted
below is the 'Lórien' candidate, ds-lorien, which is not included as a
card). `README.md` in this directory says what each script reads, writes
and needs.

## 0. The insight (scope law for every route below)

Certified |W0| needs the PERIOD VECTOR AT A POINT, not a scalar ODE. Keep the
D-module as a small first-order MATRIX system v' = A(s) v (rank rho ~ 12-17,
entry degrees ~ 600) and transport THAT. Elimination to the scalar operator is
the measured blow-up, twice over:
* ads-5-81-3213: rank-17 connection with entry deg <= 574 / common den 603
  eliminates to a scalar of ORDER 13 with z-degree R > 4052 (receipt:
  cards/ads-5-81-3213_bank/g3_dmodule_receipt.json, rho_k ratrec fail at 8104).
  Exact scalar = 17-20k nodes/prime x ~40 primes (not attempted).
* conifold-curve cards: minimal scalar operators exceed every (r,s) template
  with (r+1)(s+1) <= ~430 (lorien, series to N=561);
  the support enumeration cost grows like M^4.6-M^7.5 in the truncation
  order, which kills the data route to a scalar operator.
Pfaffian law: NEVER form the cyclic-vector scalar. All gates and transports
act on (A(s), basis B, seed vector). Fail-closed: no verdict row without every
gate in the route's table below.

## 1. Shared contract and toolkit mapping

Object: (B, A, den) with B = module basis (monomials in theta_i on the curve),
A in M_rho(F_p[s]) per prime then M_rho(Q(s)) exact, den = common denominator
poly. On disk: the exact connection file written by the lift step (numerator
coefficient matrix + den coefficients, per-entry degree table, primes used,
max coefficient height, stability flag; see `README.md`).

| need | tool |
|---|---|
| per-node closure solve mod p | the ads-5-81-3213 D-module closure solver (append mode) with its generic box-operator supply; not included in the package (see `README.md`, `TERRIER_REFERENCE_BANKS`); its conifold-frame port is not included in this release |
| big-int-safe mod-p poly ops | the `pmodp` module of the same D-module bank, with block size BLK=(2^63-1)//(P-1)^2 so that the blocked convolution cannot overflow int64 at 31-bit primes |
| exact rational reconstruction | `tools/ratfit` in this repository: ratrecon_with_denom (den 603 KNOWN -> per-entry problem is linear), thiele_loo, screen_modp; CRT+Wang lift as in `tools/numkin` |
| transport + landing | `tools/wayfinder` in this repository: DESystem contract, transport_fixed_eps, frobenius (matrix Frobenius + LS landing), gate (two-precision digit gate) |
| ball rule | arb balls END-TO-END including every reduction, as in the package's transport layer (a float64 reduction loses the radius) |
| coni curve restriction | `pipeline/restrict_op.py` + the conifold box-operator supply (`conipfv/coniop/boxcheck.py`, `shiftbox.py`; 50-op Gamma-ratio enrichment) |
| mod-p series for held-out gates | a streaming mod-p support enumerator (jet variant = add (q.M)^k weights to the specialized leaf); not included, and subsumed in practice by the exact gate PB-3 |
| tower seeding at z_cf | `conipfv/coniop/coni_pack.py` frobenius_tower / verify_tower (regression-tested on the DKMM example) |

Costs below are single-node, single-core measurements unless stated otherwise.

## 2. Route A — ads-5-81-3213 (AdS verdict vs published 2.03778e-23)

Inputs (reference data for ads-5-81-3213; only `g3_dmodule_receipt.json`
ships in `cards/ads-5-81-3213_bank/`, the rest is supplied by the user under
`TERRIER_REFERENCE_BANKS`, layout in `README.md`):
* connection samples at two 31-bit primes (2147483647, 2147483629) — 8192
  ground-truthed nodes/prime (consecutive from 34567891 / 51000001), 17x17,
  appendable.
* Entry degree table: max num deg 574, common den deg 603 (deg known EXACTLY
  per entry from the per-prime interpolation; the lift writes it out).
* the exact fundamental-period series to s^3000 at s=0 (MUM anchor); the
  flux frame (F,H) of the card with its C2 half-D3 verdict-scope string,
  INHERITED by any verdict row.
* Target: s_vac from the vacuum layer (`pipe_vac.py` frame); verdict gate
  ~15 digits.

### A-1 Exact reconstruction of A(s) over Q
Per prime p: den_p(s) already pinned (deg 603); entry num via
ratrecon_with_denom(xs, ys, Qknown=den_p) -> 575-coeff linear solve per entry,
8192 nodes >> 575 + 64 LOO held-out. 289 entries x 2 stored primes: cheap
(<= 1 h; the Krylov step measures poly ops at ~160 s/prime).
Across primes: CRT numerator+den coeffs, Wang-lift to Q. Coefficient HEIGHTS
are the unknown -> adaptive prime schedule:
* k=2 (the two stored primes) -> attempt lift; ACCEPT only if (i) lift stable
  when recomputed from k-1 primes on a random 5% coeff subsample, and (ii)
  full held-out gate PA-1 passes at a FRESH prime never used in the lift.
* else harvest more primes. Per-prime harvest cost (MEASURED, 1 core):
  closure at 8192 nodes = 2674 s + Krylov/interp 160 s ~ 47 min/prime.
  Embarrassingly parallel over primes and nodes: 8 cores => ~6 min/prime.
* Cost: each 31-bit prime buys ~9.3 digits of height. Budget table:
  H <= 25 digits -> 6 primes (~5 h serial on one core, trivially split);
  H <= 120 -> 26 primes (~20 h on one core / ~2.6 h on 8 cores). Proceed
  while k <= 32; if the lift is still unstable at k = 32 primes the heights
  are scalar-class after all: stop and re-examine the basis (a bad basis
  inflates heights; try an LLL-reduced / rescaled basis before buying more
  primes).

### A-2 Anchor at s=0 (Frobenius seeding from the D-module receipt)
s=0 is the MUM-side regular-singular anchor. Build the 17-dim local solution
space with `tools/wayfinder` frobenius (matrix Frobenius, integer-gap
shearing if needed) directly on the exact A(s):
* Indicial data of A at s=0 must reproduce the D-module receipt structure: 12
  geometric exponents (MUM block, 4-fold log tower on the h=5 curve frame)
  + the 1 EXTRA solution (order 13 = 12+1 receipt line). If the indicial
  block disagrees -> reconstruction or basis error, HALT.
* Seed identification: the branch whose analytic part matches the exact
  fundamental-period series (match >= 200 exact coefficients after basis
  contraction — an EXACT gate, not float). The remaining geometric mates
  come from the Frobenius log tower; the extra solution is identified as
  the complement.
* GENUINE UNKNOWN (FLAG F3): the W0 contraction needs the GEOMETRIC 12-block
  in the flux frame (F,H vectors). The 17-basis -> period-basis map
  is fixed by Frobenius exponents + leading monomials at MUM, but the
  pairing normalization (which combination is Pi_0..Pi_11 vs the extra
  solution) is not fixed by that data alone. It needs a derivation plus an
  exact cross-check: contract the candidate 12-block against w0/J1-style
  series identities before any transport.

### A-3 Transport 0 -> s_vac
* Singular locus = roots of den (deg 603) + true sings of A. Classify
  apparent vs true by residue rank of A at each root inside the path
  corridor (cheap: rank of the polar part mod p first, exact only near the
  path). Path: real segment [0, s_vac] with certified clearance delta to
  the nearest true singularity; if a den root sits ON the segment and is
  apparent, detour by a certified semicircle (`tools/wayfinder` handles
  ordinary points only — the corridor must be sing-free).
* March: transport_fixed_eps Taylor endpoint transport on the 17x17 exact
  system in ball arithmetic (arb) end-to-end; NO float64 reduction anywhere
  (it would lose the radius). Working precision 60 digits for a 15-digit
  verdict (guard >= 3x); step acceptance by the stepped-majorant rule of
  the package's transport layer (`pipeline/pipe_transport.py`).
* Cost model: steps ~ len/delta x log(1/tol); per step ~ O(rho^2 x deg) ball
  ops ~ 17^2 x 600 x order-60 ~ 1e7 ball-ops/step. Even 1e4 steps is
  minutes-to-hour class on one core. Transport is NOT the bottleneck; A-1 is.

### A-4 Gate table (route A; ALL required, fail-closed)
| gate | statement | oracle |
|---|---|---|
| PA-1 | exact A(s) reproduces ALL 289 entries at 256 nodes of a FRESH prime (0 diffs) | new closure harvest, prime not in lift |
| PA-2 | exact A annihilates the held-out series: v(s) built from the s^3000 fundamental-period jets satisfies v' - Av = O(s^{3000-603-1}) EXACTLY | the exact s^3000 series (held out; never used in reconstruction) |
| PA-3 | indicial structure at s=0 = 12 MUM + 1 extra; seed matches >= 200 series coeffs exactly | Frobenius + series match |
| PA-4 | geometric 12-block identification cross-check (F3) passes exact series identities | A-2 gate |
| PA-5 | path corridor certificate: no true sing within delta; every crossed den root proven apparent (exact residue-rank) | exact linear algebra |
| PA-6 | |W0| ball at s_vac: radius <= 5e-16 rel; two-precision digit gate (60 vs 90 dps) agrees to >= 15 digits; verdict row inherits the C2 half-D3 scope string | `tools/wayfinder` gate |

Status (route A): A-1 through A-3 are costed above and the reference data
exists. Stop conditions: k > 32 primes (A-1); indicial mismatch (PA-3); the
12-block cannot be pinned (F3) — that needs a derivation, not more compute.
How the built route A actually lands (direct certified summation inside a
proven convergence disk instead of marching, because the exact A(s) has an
apparent pole of order 19 at s=0 in the closure basis) is described in
`ads581/verdict.py` and in `periods/PFAFFIAN.md`.

## 3. Route B — dS cards (conifold-frame connection; lorien)

Object: the card GKZ ideal restricted to the coni curve (`restrict_op.py` +
the frame-free support law of `conipfv/coniop/gseries.py` — enumerate in the
UNIMODULAR ROW FRAME, never the card orthant). Output per card: (B, A(s))
on the curve, NOT a scalar coni_op — this replaces the scalar-operator
search, which fits no affordable (order, degree) template (section 0).

### B-1 Module basis on the coni curve (the rank question)
Start from the card 12-basis (2h+2 for h=5) and run on-curve closure with the
enriched conifold box-op supply (50 Gamma-ratio ops, |l|_1 <= 2, each
termwise-gated on exact lattice points; boxcheck/shiftbox generators; for
lorien use the corrected op[1] operator, never the card string).
Closure protocol as for ads-5-81: DEGMAX ladder 4 -> 6 -> 8; basis
enrichment by the persistent interior free set (that is exactly how rank 17
was found). GENUINE UNKNOWN (FLAG F1): closure rank rho and whether closure
terminates at affordable DEGMAX on dim-8 cones (lorien). Budget: accept
rho <= 30; if the free set still grows at DEGMAX 8, HALT the card (the
route stays fail-closed) — do NOT chase closure with compute.

### B-2 Connection sampling mod p (per-node closure, NOT series enumeration)
A(s0) mod p per node = linear solve of the box-op relations on B at numeric
s0 (closure solver ported to the coni frame). This is ALGEBRAIC in s0 — no
support enumeration per node; the per-node cost scales like the ads-5-81
measurement (0.33 s/node at rho=17, h=5), NOT like the M^7.5 series cost.
Node schedule: consecutive nodes, 2 primes to start, appendable .npz sample
files exactly as for ads-5-81. Degree pin: per-entry deg via mod-p
ratrec (thiele) on the first prime -> D = max entry deg; harvest
N_nodes = 2(D + den_deg) + 128 held-out per prime.

### B-3 Exact reconstruction + the held-out series gates
Known-denominator reconstruction (`tools/ratfit`) where a common den emerges
(expect the curve discriminant x small factors), else thiele_loo per entry;
CRT/Wang with the SAME adaptive prime schedule and stop condition as A-1.
The exact graded series computed for the scalar-operator search — lorien
N=561 (5201 s) — serves here as the held-out gate: build v(s) = (B applied to the exact series) and check
v' - Av = O(s^{N - D - 1}) EXACTLY over Q. A mod-p pre-gate via the
streaming enumerator's jet variant (leaf weighted by (q.M)^k jets and
theta-monomial factors) at a fresh prime is 19x cheaper and catches errors
before the exact pass.

### B-4 Tower seeding + transport at z_cf (per coni_pack)
With exact A(s): indicial data at s=0 and at z_cf feed `coni_pack`
frobenius_tower (matrix form: run frobenius_tower on the FIRST-ORDER system's
scalar blocks — the DKMM regression stays the acceptance test), verify_tower
mandatory. Transport 0 -> s_star (lorien 0.049470929, a pinned rational
estimate) with the A-3 corridor law; landing at the conifold point via the
`tools/wayfinder` frobenius LS landing; c_log/c_tau must reproduce the CP8
pins (lorien: -32, 93/19) EXACTLY.

### B-5 Gate table (route B, per card; fail-closed)
| gate | statement | oracle |
|---|---|---|
| PB-1 | closure certificate: B closed under all supply ops at DEGMAX, rank rho stable across 2 primes + 3 jets-Krylov nodes | closure solver, boundary-leak law |
| PB-2 | exact A vs fresh-prime nodes, 0 diffs (as PA-1) | new harvest |
| PB-3 | exact A annihilates the held-out exact series to O(s^{N-D-1}) | lorien N=561 series (HELD OUT) |
| PB-4 | mod-p jet-variant series at a fresh prime also annihilated | streaming enumerator |
| PB-5 | conifold indicial pins: c_log = -M_cf n_cf, c_tau = lam exact | CP8 pins (`conipfv/DESIGN.md` section 5) |
| PB-6 | tower verify_tower PASS + DKMM regression unchanged; s_star clearance from A's true singular locus | coni_pack |

### B-6 Cost model and limits (per card)
Assume rho <= 30, D <= 1500 (unknowns F1/F2): nodes/prime ~ 2(D+den)+128
~ 6.2k; per-node solve ~ 0.33 x (rho/17)^3 s <= 1.8 s -> <= 3.1 h/prime/card
on one core, linear split over 8 cores (~25 min). k primes: 6-26 as in A-1.
Limits: rho <= 30 AND per-node <= 2 s AND D <= 1500 AND k <= 32.
Any limit exceeded -> the card gets a bracket row, not a verdict. Order:
an h=5 orthant-easy card first, lorien second (h=8, the stress test) —
lorien decides the 24-card coni class.

## 4. Flags — genuine unknowns (each must be settled before its stage runs)

* F1 (route B, blocking B-1): module-basis CLOSURE on the coni curve. The
  ads-5-81 lesson is that card bases under-close (12 -> 17, found only by
  chasing the persistent free set); on dim-8 cones closure may not terminate
  at affordable DEGMAX, and rho is unbounded a priori. Mitigation: DEGMAX
  ladder + rho <= 30 acceptance + hard halt. This is the single biggest
  route-B risk. (Settled for h=5 coni cards: rank 16.)
* F2 (both routes, cost): coefficient HEIGHTS of exact A(s) over Q are
  unknown a priori. The Pfaffian bet is that connection entries are
  height-tame even when the eliminated scalar is not (R > 4052) — plausible
  (elimination multiplies heights) but not proven in general. The adaptive
  prime schedule turns this into a bounded experiment (k <= 32) instead of
  an assumption. (Measured: 97-102 bit heights at rank 17 for ads-5-81,
  97 bits at rank 16 on the h=5 conifold-frame example.)
* F3 (route A, blocking A-2/PA-4): 17-basis -> flux-frame geometric 12-block
  map (which solutions are the periods contracted with the F,H vectors).
  Needs a derivation; the 13th solution MUST be proven absent from the W0
  contraction, not assumed. (The built route A settles it by a
  rank-saturating exact fit against GV-built model periods, held out at
  higher order; `ads581/verdict.py` stage 3.)
* F4 (both, transport): apparent-vs-true classification of den roots near
  the path. If a TRUE singularity sits between 0 and the target inside
  clearance, the real-segment path dies; complex detour needs a corridor
  certificate `tools/wayfinder` does not currently emit for matrix balls —
  a small tool gap (extend transport_fixed_eps acceptance to a disk-chain
  certificate; do NOT rebuild the marcher).
* F5 (route B): den structure on the coni curve unknown a priori (no known
  common-denominator analog of the 603). If entries have unrelated dens,
  per-entry thiele doubles the unknown count -> node budget x2 (folded into
  B-6 bounds). (Measured on the h=5 conifold-frame example: a common den of
  degree 49 exists.)
* F6 (route A, input): s_vac provenance — the `pipe_vac.py` frame value must
  be re-derived under the order-13/rank-17 structure (basis length no longer
  equals operator order); confirm the target point before A-3.

## 5. Layout

    pfaffian/
      DESIGN.md            this note
      README.md            what each script reads, writes and needs
      ads581/              route A (ads-5-81-3213): per-prime interpolation,
                           CRT/Wang lift, gates PA-1 / PA-2, held-out jets,
                           certified verdict
      (route B scripts — box-op supply, closure certificates, connection
      sampler, interpolation, lift, gates PB-2 / PB-3, exact series oracle,
      entropy certificate + certified landing, ball marcher — are not
      included in this release)

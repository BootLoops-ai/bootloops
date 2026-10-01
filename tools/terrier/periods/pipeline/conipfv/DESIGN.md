# CONI-PFV TRANSPORT EXTENSION — DESIGN (periods/pipeline/conipfv/)
## Why the dS cards of arXiv:2406.13751 fail g2, and what replaces the monomial curve

Inputs: the pipeline modules one level up (family / curve_from_flux / pipe_lib /
pipe_transport / pipe_vac) with their stepped-majorant transport legs; the
near-conifold conventions D1-D5 (section 2; taken from arXiv:2009.03312); and
arXiv:2406.13751 sec. 3, whose equations are cited below by their LaTeX labels
(eq:PFVfflux .. eq:cookedness) or compiled numbers (3.x).  Other convention
labels used: C1 = flux-vector ordering Pi = (F_0, F_a, X^0, X^a); B2 = a_mat
mod-2 representative rule (2406 eq. 2.52); B3 = conifold-shifted c2D (2406
eq. 3.17); A7 = normalization/dressing of the published W0 column (not
pinned).  The 30 dS example cards referred to below are the paper's candidate
vacua; none of them is included in the package (only the DKMM control card
ships); statements about them record what the same code reported on them.

## 0. The structural fact (measured on the cards)

Every 2406 dS card fails `curve_from_flux.pffv_curve` (g2) EXACTLY:
K.N^-1.K != 0 (e.g. ds-lorien 93^2/27230). The flat
direction is NOT the PFFV p = N^-1 K; it is the coni-PFV solution of
2406 eqs (eq:detN..eq:integrality_condition_II):

    N.p - K = lam * q_cf ,   q_cf.p = 0 ,   K.p = 0  (Diophantine gate)

with q_cf the conifold curve class and lam the conifold residual (the
published Table-3 K': 93/19 for ds-lorien, verified exactly in the card
NOTE). Consequences:
  * p lies ON the conifold facet K_cf of the mirror Kahler cone (q_cf.p = 0);
    nu = r*p has zeros (conifold direction) and may have negative CARD-frame
    entries (ds-lorien: p exits the ambient cone through flop walls) — the
    per-coordinate monomial substitution z~_a = s^{nu_a} is meaningless.
  * The conifold-class instanton tower sits entirely at s-order 0: it can
    never be expanded in s; it must be resummed in the transverse coordinate
    z_cf = q_cf.z (the Frobenius log branch, freeze D1/D2/D5).
So the restriction object is not a curve into the MUM chart but a
2-parameter CHART (s, z_cf) with a 1-parameter family at z_cf = 0
approaching the conifold divisor, plus a rank-2 transverse Frobenius block.

## 1. (a) The restriction chart in the CONIFOLD frame

Coordinates. tau = axio-dilaton; flat locus z^a = p^a tau + xi^a z_cf with
xi = M/M_cf (2406 choice, q_cf.xi = 1), M_cf := q_cf.M != 0. Define
r = lcm(denominators of p) (frame-invariant: card bases are unimodular) and

    s := e^{2 pi i tau / r}   =>   instanton weight of class q:
    e^{2 pi i q.z} = s^{q.nu} * e^{2 pi i (q.q_cf-part) z_cf},  nu := r p.

Degrees q.nu are INTEGERS >= 0 for every effective q (gate CP4 below), with
q.nu = 0 iff q is a multiple of q_cf. The restriction "curve" used here is
the class-lattice grading by q.nu — NOT a coordinate substitution. The
1-parameter family approaching the conifold divisor is

    s in (0, s_vac],  z_cf = 0   (on-divisor leg; bulk periods only),

and the vacuum locus is reached by the transverse Frobenius block at
z_cf = <z_cf> (exponentially small, freeze D4). As in the DKMM vacuum
layer, the transverse direction is never recovered from the on-curve data;
here it is DESIGNED IN from the start as the z_cf block.

Frame law. All contractions (q_cf.p, q.nu, M_cf, lam, K.p, tadpole) are
frame-invariant pairings; the implementation works in the CARD frame and
cross-checks the card's coni_basis_crosswalk block when present
(ds-extra-p8-a carries M_coni/K_coni/nu_coni; T_ours_to_coni maps exactly).

## 2. (b) Conifold Frobenius basis (freeze Group D, PINNED-V)

Local system in z_cf at fixed bulk LCS point: exponents (0,1) resonant pair;
the suite's resonant/log-tower pattern applies with ONE log layer:

    A-branch:  Pi_A  = z_cf * u(z^alpha, z_cf),          u(.,0) = 1
    B-branch:  Pi_B  = (n_cf/2 pi i) Pi_A log(-2 pi i z_cf) + f(z^alpha, z_cf)

conventions PINNED (do not re-derive): D1 F_cf = (n_cf/2pi i) z_cf ln z_cf + f
(2009.03312 eq (3.1)); D2 log scheme with -2 pi i INSIDE the log,
F = n_cf z_cf^2/(4 pi i) ln(-2 pi i z_cf) + sum F^(n) z_cf^n/n!
(2009.03312 (3.9)-(3.10), 2406 (3.6)); D3 branch/homotopy bookkeeping OURS
(kept explicit in this package); D4 sqrt(pi/2) normalization + vev formula + INTEGER K law
(never the Table-3 K'); D5 Euler-reflection resummation validity =
a boundary to DOCUMENT per vacuum (GV-nilpotent order one: GV(k q_cf) = 0
for k >= 2 — gate CP7 checks the pinned window is consistent with this).

Superpotential block (2406 W^(1), with the coni-PFV flux choice
P_beta = (1/2)(A.M)_beta, P_0 = (1/24) c2'.M which kills the polynomial
inhomogeneity — gate CP5 checks the integrality that makes P exist):

  sqrt(pi/2) W = W_bulk(s) + z_cf W1(s, z_cf) + O(z_cf^2 corrections), with
  W1 = -M_cf (n_cf/2 pi i)(log(-2 pi i z_cf) - 1) + kappa(M, xi, p) tau
       + P_res.xi + (1/2 pi i) sum_{q != q_cf} n_q (q.M)(q.xi) Li_1(s^{q.nu})

where kappa(M, xi, p) = (p.N.M)/M_cf exactly (Fraction), P_res.xi = the
integer-shift residue of the P-choice (reported, convention A7 PIN-REQ).
The Frobenius block coefficients {z_cf, z_cf log(-2 pi i z_cf)} are assembled
EXACTLY from card data; only Li evaluations are numeric (mpmath, labeled).

## 3. (c) Certified transport: MUM -> conifold-adjacent -> vacuum locus

Three legs, reusing the pipe_transport stepped majorants UNCHANGED:

  LEG 1 (bulk, on-divisor). The z_cf = 0 restriction of the 2-parameter
  PF system is an ODE L_s in s with MUM at s = 0 (bulk LCS). Once a card
  carries L_s (`coni_op` field: theta-form integer coefficients, same
  schema as operator_LS), pipe_transport runs verbatim: mum_eval with the
  factored-indicial majorant, local_leg chain per schedule(), landing at
  s_vac. NOTHING here is new code — coni_transport.py wraps run_route.
  Deriving L_s is a separate route: eliminate z_cf from the card's raw GKZ
  boxes at z_cf = 0 (the boxes are stored UNGATED on every card), or
  annihilate the class-graded restricted series (SS1) — restrict_op.py
  applies once the series exists.

  LEG 2 (transverse jet). Alongside Pi, transport d Pi/d z_cf |_{z_cf=0}
  via the coupled first-order system from the 2-parameter ideal (the
  pipe_vac pattern: connection matrices + inhomogeneous jet legs; the
  doubled system reuses local_leg on block-companion form). This is the
  same design as the DKMM off-curve leg (pipe_vac) — required because tau
  and the vacuum are NOT on-divisor-recoverable.

  LEG 3 (conifold-adjacent landing). At s_vac, evaluate the Frobenius
  block at z_cf = <z_cf>. NOT a series leg: z_cf ~ 1e-6 sits INSIDE the
  conifold disc, and the block is the exact local basis; the certified
  error is the O(z_cf^2 log z_cf) truncation with explicit majorant
    |R| <= C2 |z_cf|^2 (1 + |log(-2 pi i z_cf)|),
  C2 assembled from card kappa/GV data (documented in coni_transport
  docstring); ball arithmetic end-to-end (rule: never reduce through float64).

## 4. (d) W0 in the coni-PFV racetrack form + published floor

Racetrack (2406 eq:Widef class): W_N = -(1/4 pi^2) sum_{q.nu = N} n_q (q.M)
Li_2(s^N); physical coefficient carries sqrt(2/pi) (verified against the
4627-main fingerprint sqrt(8/pi^5)*(-2,252)). Two leading levels give the
tau vev by the 2406 two-term formula; Newton refinement on the full pinned
ladder follows (2406's own step-2 pattern). Then

  W0_est = sqrt(2/pi) | W_bulk_eff(<tau>) + <z_cf> W1(<tau>, <z_cf>) |

with <z_cf> from the D4 vev law in the 2406 Q_throat form
(eq:conifold_vev): |z_cf| = (1/2pi) exp(-2 pi Q_throat/(g_s M_cf^2 n_cf)),
Q_throat = Q_flux - g_s ||M||^2, ||M||^2 = -M.K(M).M at the vev (kappa part
exact, Li_1 part numeric). Published floor (eq:zbound/eq:cookedness class,
the "(3.36)-(3.38)" battery): gate CP9 checks <z_cf> >= (1/2pi)
exp(-2 pi Q_O/(g_s M_cf^2 n_cf)) with Q_O = 2(Q_D3 - 1), and reports the
W0 floor. ALL SS4 numbers are ESTIMATES (mpmath dps 50, labeled) — the
certified W0 is Leg-1+2+3 output; estimates only gate consistency with the
card/repo pins (tau_im, z_cf, Table-3 W0) at order-of-magnitude tolerance
(the D4 numeric leg stays labeled; it is not silently closed).

## 5. Gate battery: g2' (CP1-CP6) and g3' (CP7-CP10)

g2' = `coni_curve.coni_pfv_curve(card)` — EXACT (Fractions), fail-closed:
| gate | statement |
|---|---|
| CP0 | routing: K.N^-1.K == 0 => NOT-CONI (defer to g2/pffv_curve); else proceed. det N != 0 |
| CP1 | UNIQUE primitive pinned class q_cf with: lemma system solvable, K.p = 0 exact, and CP4 positivity — exactly one survivor across all pinned GV classes |
| CP2 | lemma residual: N.p - K = lam * q_cf EXACT; lam reported (= published K' when carried: ds-lorien 93/19) |
| CP3 | n_cf = GV(q_cf) >= 1; M_cf = q_cf.M != 0; D4 validity flags (K,M sign/size) reported |
| CP4 | facet interiority vs pinned window: q.nu >= 0 for ALL pinned q, q.nu = 0 iff q parallel q_cf; nu = r p primitive integer |
| CP5 | flux integrality (2406 eq:integrality I/II): bulk-projected (A.M) even (HNF basis of q_cf-perp), c2D.M in 24Z (card c2D = conifold-SHIFTED c2', freeze B3 — the shift is load-bearing, raw fails) |
| CP6 | tadpole -M.K/2 <= Q_D3; symplectic flux vectors F/H emitted (C1 ordering) |
Cross-gate: when the card carries coni_basis_crosswalk, nu/M_cf/n_cf/lam are
matched against nu_coni/M_coni/conifold mismatch EXACTLY (T map applied).

g3' = `coni_frobenius.g3_series(card, cur)` — conifold-frame series gates:
| gate | statement |
|---|---|
| CP7 | racetrack ladder from gv_pinned: levels N = q.nu, exact integer sums S_N = sum n_q (q.M); level-0 = conifold ray ONLY (GV-nilpotency window check, D5); ladder matched vs card racetrack pins when carried |
| CP8 | Frobenius block: exact {z_cf, z_cf log(-2 pi i z_cf)} coefficients; log coefficient = -M_cf n_cf/(2 pi i) (D2 scheme pin); tau-linear coefficient (p.N.M)/M_cf exact |
| CP9 | vev + floor: <tau> two-term + Newton on the ladder; <z_cf> D4/Q_throat form; published floor eq:zbound respected; consistency vs card tau_pin / repo z_cf / NOTE estimates at labeled tolerance |
| CP10 | W0_est racetrack form incl. z_cf W1 term; floor eq:cookedness reported; vs published Table-3 W0 when carried (ESTIMATE gate, not certified) |

## 6. Module map + battery slot-in

| file | role |
|---|---|
| coni_curve.py | (a) + g2' CP0-CP6; consumes Card unchanged (family.load_card) |
| coni_frobenius.py | (b)+(d) + g3' CP7-CP10; exact block + mpmath estimates |
| coni_transport.py | (c) legs 1-3 wrapper over pipe_transport (parametric; fail-closed until a card carries coni_op + routes) |
| run_conipfv.py | orchestrator: card -> g2' -> g3' [-> transport]; transcript + gates JSON; exit nonzero on FAIL |

Battery: run_conipfv.py is invoked per-card exactly like run_pipe.py stages;
a coni card's verdict line is CONI-PENDING (g2'+g3' PASS, transport awaits
coni_op) or CONI-FAIL. PFFV cards route to the existing g2/g3 unchanged
(CP0). DKMM is the control: must report NOT-CONI and touch nothing.

## 7. What this extension does NOT do (scope law)

* It does not derive the bulk operator L_s (elimination/annihilator route,
  coniop/) and never fakes certified transport without it.
* It does not close the D4 numeric leg — vev numbers are labeled estimates;
  only the PINNED formulas are hard-coded.
* It does not touch run_pipe.py, family.py, or any reference data (DKMM
  regression stays string-identical); cards are consumed UNCHANGED.

## 8. Implementation findings (full sweep, 10.6 s / 35 cards)

* CP5 evenness of bulk (A.M) is REPORTED, not fail-closed: the
  a-representative is frame-dependent mod 2 (freeze B2 residue PIN-REQ);
  ds-lorien's card frame has odd entries while its own NOTE gated a.M in Z
  only.  Fail-closed CP5 = integrality (a.M in Z^h, c2'.M in 24Z).
* EXACT IDENTITY (load-bearing): the tau-coefficient of the z_cf Frobenius
  block is c_tau = lam (the conifold residual = Table-3 K') once the
  D4-pinned -tau K z_cf term is included — the 2406-printed W^(1) omits it.
  Equivalent: Q_throat(kappa-part) = lam*M_cf.  Verified: the F-term vev
  and the eq:conifold_vev Q_throat form agree to < 0.05 e-folds on all 29
  passing cards (F-term-vs-Q_throat cross-gate, bar 1.5 e-folds).
* The tau vacuum is the JOINT racetrack: bulk instanton levels vs the
  conifold sector as an effective instanton of fractional degree
  r*lam/(M_cf n_cf).  Recorded discrepancy: quoting the conifold exponent
  as K'/M_cf (no n_cf) gives m_eff 11.6/38 for ds-lorien, outside the GV
  window; the PINNED D4/eq:conifold_vev formulas carry n_cf => m_eff 5.8/38,
  INSIDE the window (which is why g3' closes on ds-lorien).  The pinned
  formulas govern.
* Scoreboard: 29/30 dS cards CONI-PENDING (g2'+g3' PASS; lam matches the
  published K' on every crosswalk-bearing card); ds-manwe CONI-BLOCKED
  (frame dictionary not closed — card-documented, correctly fail-closed);
  dkmm + 4 ads cards NOT-CONI (CP0 routing control PASS).
* Estimates are window-limited PFV step-1 numbers (tau ~10 percent from
  repo full-F-term values); they gate CONSISTENCY only.  Certified numbers
  come from transport legs 1-3, which fail-close until a card carries
  coni_op/coni_towers/coni_routes/s_star.

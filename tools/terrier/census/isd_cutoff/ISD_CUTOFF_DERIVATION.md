# G5 ISD CUTOFF — DERIVATION (formula fixed before any count)

HARD RULE: the cutoff VALUE derives from a FORMULA over previously recorded
period quantities, fixed BEFORE any G5 count is computed under it —
never a number chosen after seeing counts.
NO COUNTING: this file contains ZERO enumeration.
Geometry: G5 = mirror of the Reye congruence Calabi-Yau threefold (AESZ22 / AESZ118),
with the rank-2 attractor class point at z = -1.
Companion scripts in this directory: `isd_ref.py` (floating-point check of the theorem,
corollary, reduction chain and witness family on synthetic Hodge frames) and `sign_probe.py`
(300-trial sign/orientation probe, both pairing orders).

## 1. Setup and conventions

- Flux pair (f,h) in Z^4 x Z^4 (f = RR F3, h = NSNS H3 charges) on X = mirror Reye congruence
  CY3, h^{2,1} = 1, b_3 = 4 (AESZ22 frame at z=0, AESZ118 frame at x=0, with an integral transfer
  matrix between the two MUM frames).
- Sigma = witness.sig(4) (census/dedup/witness.py), the symplectic form of the integral frame;
  every period matrix used is verified against it. Budget: B = |f^T Sigma h|, ABSOLUTE integer
  units, grid 1..B_max = 16. Budget SL2- and monodromy-invariance are proven against this Sigma.
- Pi(z) = period vector of the holomorphic 3-form Omega in the integral symplectic frame;
  D_z Pi = del_z Pi + (del_z K) Pi its Kaehler-covariant derivative. Both are available as
  certified balls from the period machinery: periods/transport.py (exact symplectic contract;
  ball arithmetic END-TO-END, never reduced through float64) and the Picard-Fuchs web on the
  AESZ22 and AESZ118 operators.
- GVW superpotential: W(z, tau) = (f - tau h)^T Sigma Pi(z); moduli = z and the axio-dilaton tau,
  physical iff Im tau > 0. G3 := f - tau h.
- CENSUSED CLASS: ALL certified F-flat vacua (D_z W = D_tau W = 0), labeled W=0 / W!=0.
- GKP/GVW fact [CITED hep-th/0105097 sec 2-3; hep-th/9906070]: joint F-flatness in (z, tau) on a
  CY3 orientifold holds iff G3 is imaginary-self-dual (ISD), *G3 = i G3, i.e. G3 in
  H^{2,1} + H^{0,3}. The W=0 subclass (G3 in H^{2,1}) is a fortiori ISD. So EVERY censused G5
  vacuum is ISD — ISD is the physical norm condition that makes the fixed-budget class count
  finite (sec 5-6).

## 2. The Hodge Gram H(z) from period data

- Definition: H(z)_{ab} := int alpha_a wedge *_z alpha_b = (Sigma C(z))_{ab}, where C(z) is the
  Weil operator on H^3 (equal to *_z on the middle cohomology of a CY3): C = i^{p-q} =
  (-i, +i, -i, +i) on (3,0), (2,1), (1,2), (0,3). H is real symmetric and POSITIVE DEFINITE at
  every smooth z (Hodge-Riemann bilinear relation II). The analogous weight-4 statement
  (Calabi-Yau fourfolds: H = Sigma C > 0 with a Rayleigh cap giving a finite flux box, as in
  census/tile_engine) has a different positivity
  MECHANISM: there C^2 = +1 (a Sigma-compatible involution). On weight-3 H^3 the Weil operator
  has C^2 = -1 (i^{p-q}, p+q = 3, p-q odd), so that argument does not transfer verbatim; here
  symmetry of H follows from C being a Sigma-isometry PLUS antisymmetry of the pairing, and
  positivity is carried by the weight-3 Hodge-Riemann relations, instantiated by the explicit
  rank-4 formula below (whose two normalizers are certified positive per cell). No step of THIS
  proof leans on the weight-4 involution. Here b_3 = 4, so H is 4x4.
- Explicit evaluation at h^{2,1} = 1 (all quantities certified balls): for real v in R^4,
    v^T H(z) v = 2 |v^T Sigma Pi|^2 / N0(z) + 2 |v^T Sigma D_z Pi|^2 / N1(z),
    N0(z) := i conj(Pi)^T Sigma Pi        ( = e^{-K_cs} = i int Omega wedge conj(Omega) ),
    N1(z) := -i conj(D_z Pi)^T Sigma D_z Pi  ( = e^{-K_cs} g_{z zbar} ).
  Polarized form, same expansion:
    f^T H(z) h = 2 Re[(f^T Sigma Pi) conj(h^T Sigma Pi)] / N0
               + 2 Re[(f^T Sigma D_z Pi) conj(h^T Sigma D_z Pi)] / N1.
  DERIVATION: expand the Poincare dual of v on the Hodge basis {Omega, D Omega, conj(D Omega),
  conj(Omega)}; Griffiths orthogonality kills every symplectic pairing except
  int Omega wedge conj(Omega) and int D Omega wedge conj(D Omega); apply the C eigenvalues above.
  Rank count: 2 rank-2 real quadratic forms (Re + Im parts of each functional) = rank 4 = positive
  definite exactly when {Sigma Pi, Sigma D Pi} give 4 R-independent real functionals (period map
  immersive) and N0, N1 > 0.
- POSITIVITY/ORIENTATION CERTIFICATE (fail-closed, orientation-robust form):
  the identifications N0 = e^{-K_cs}, N1 = e^{-K_cs} g_{z zbar} hold under the pairing
  order int U wedge V = v^T Sigma u; under the opposite order both flip sign TOGETHER on genuine
  period data. Certificate: at every evaluation cell the ball enclosures must certify a COMMON
  sign s in {+1,-1} with s*N0 > 0 AND s*N1 > 0 (imaginary enclosures straddling 0 per the exact
  symplectic contract), s CONSTANT across all certified cells of K_G5 (a global frame
  orientation); H is then built from (s*N0, s*N1), which is the true Hodge Gram in either
  orientation. Mixed signs within a cell, or s varying across cells, FAILS the certificate. A
  cell whose certificate fails yields NO exclusion from that cell — never a false exclusion
  (the over-count-safe direction). Analytically the common-s certificate holds at every smooth
  census point (Hodge-Riemann); the ball test instantiates it per cell. Checked on synthetic
  Hodge frames in both orientations, 300/300 trials each (`sign_probe.py`).

## 3. THEOREM (tau-eliminated ISD bound — the cutoff's engine)

Let z be a point where the sec 2 certificates hold, and suppose (f,h) admits an F-flat vacuum at
(z, tau) with Im tau > 0. Then, with B := |f^T Sigma h| (the absolute budget), B > 0 and
    N_ISD(f,h; z) := sqrt( (f^T H f)(h^T H h) - (f^T H h)^2 )  <=  |f^T Sigma h| = B.
[WHY the |.| form: a bare f^T Sigma h > 0 claim at the vacuum is WRONG — a 300/300
synthetic-Hodge-frame computation (`sign_probe.py`)
shows that under THIS FILE'S OWN certificate orientation (N0 := i conj(Pi)^T Sigma Pi > 0,
N1 := -i conj(D_zPi)^T Sigma D_zPi > 0 — which pins the pairing order int U wedge V = v^T Sigma u)
every ISD vacuum has f^T Sigma h < 0. The |.| form below is
orientation-invariant and is all that the cutoff, Tier-1/2, witness-kill and finiteness use.]
PROOF (each step one line):
 (i)   B != 0 forces f, h R-independent.
 (ii)  int G3 wedge conj(G3) = -(tau - taubar) f^T Sigma h = -2i (Im tau) f^T Sigma h under the
       certificate-pinned order int U wedge V = v^T Sigma u; the opposite order flips this sign
       [f^T Sigma f = h^T Sigma h = 0 by antisymmetry of Sigma].
 (iii) Hodge-split G3 = G+ + G- (*G± = ±i G±): int G3 wedge conj(G3) = i||G+||_H^2 - i||G-||_H^2
       [type-diagonality of the Hodge inner product; ||.||_H = Hodge norm; form-level,
       convention-independent].
 (iv)  => -f^T Sigma h = (||G+||^2 - ||G-||^2) / (2 Im tau). ISD vacuum: G- = 0, so with
       t := Im tau > 0:  2 t |f^T Sigma h| = 2 t B = ||G3||_H^2 > 0, and f^T Sigma h < 0 in the
       certified orientation (the mirror (f,-h) carries the opposite pairing order).
 (v)   G3 = x - i t h with x := f - (Re tau) h real => ||G3||_H^2 = x^T H x + t^2 h^T H h
       [cross terms cancel: H real symmetric, x, h real].
 (vi)  B = x^T H x/(2t) + (t/2) h^T H h >= sqrt( (x^T H x)(h^T H h) )   [AM-GM in t]
       >= sqrt( min_r ||f - r h||_H^2 * h^T H h ) = sqrt(Gram_H(f,h)) = N_ISD(f,h;z).  QED
Equality iff tau sits at the (r,t)-minimizer; the bound ELIMINATES tau exactly — no tau scan.

REMARK (sign folding, no convention split):
by (iv), in the certified orientation a pair with f^T Sigma h > 0 admits NO F-flat vacuum at
Im tau > 0; its mirror (f, -h) (a distinct enumerated pair with the same |budget|) carries the
census entry. If the integral frame realizes the opposite pairing order the branches swap — but
N_ISD, gram_E, and B = |f^T Sigma h| are invariant under h -> -h, so the cutoff test is
sign-blind either way and needs no folding rule. Checked: 300-trial synthetic-frame computation,
both orientations, |B|-identity and N_ISD <= |B| exact (`sign_probe.py`); theorem, tightness of
the AM-GM step, corollary and the sec 6 reduction chain on 200 synthetic frames (`isd_ref.py`).

COROLLARY (global relaxation — the integer workhorse). For H symmetric positive definite with
smallest eigenvalue lam(z) > 0, and f, h independent:
    Gram_H(f,h) >= lam(z)^2 * gram_E(f,h),   gram_E := (f.f)(h.h) - (f.h)^2 in Z, >= 1.
[P := [f h] (4x2); P^T H P >= lam P^T P as 2x2 forms; det is monotone on the PSD order:
 A >= B > 0 => det A >= det B.] Hence at any censused vacuum z*:
    lam(z*) * sqrt(gram_E(f,h)) <= N_ISD(f,h;z*) <= B <= B_cap.

## 4. THE CUTOFF RULE (fixed before any count; B_cap := B_max = 16)

- K_G5 (pinned evaluation set) := the census vacuum-search set — the AESZ22 + AESZ118 web cells
  (both MUM frames + the integral transfer matrix), fixed and hashed BEFORE any count under the
  formula. The census certifies a vacuum z*(f,h) ONLY by ball/Krawczyk tests on exactly these
  cells, so a censused vacuum exists only at a visited certified cell: restricting the min to
  K_G5 preserves EXACT necessity by construction. Cells failing the sec 2 certificate after the
  pinned refinement budget are EXCISED from K_G5 (the excision list is part of the hashed set) —
  a vacuum there could never be certified (the ball-blowup validity gate lands it in N_failed),
  so necessity survives excision.
- mu_lo := the certified LOWER bound, over all certified cells of K_G5, of lam(z) = the smallest
  eigenvalue of H(z); per cell take the MAX of the interval-Gershgorin and interval-Cholesky
  certified lower bounds (each is a certified lower bound, so their max is one — no method
  freedom; the method used is recorded per cell); min over cells. mu_lo is a measured quantity,
  recorded before any count.
  HARD GATE: Tier-1 is VOID unless the mu_lo measurement certifies
  mu_lo > 0. Rationale: mu_lo enters Tier-1 as mu_lo^2, which silently un-signs a nonpositive
  certified bound (fat balls can yield mu_lo <= 0 even though lam(z) > 0 analytically) and would
  make exclusions unsound. Analytically mu_lo > 0 at every smooth point (Hodge-Riemann,
  instantiated by the sec 2 certificates); the gate makes the measured instantiation fail-closed.
- CUTOFF (per deduped class, evaluated on the engine-canonical representative of the dedup):
    ADMISSIBLE  <=>  NOT certified-excluded, where a class is certified-excluded iff
    TIER-1:  mu_lo^2 * gram_E_min(class) > B_cap^2 = 256,   with
             gram_E_min(class) := min over the class's enumerated SL2-canonical representatives
             of gram_E (exact integer; SL2-invariant; monodromy merges take the min), OR
    TIER-2 (optional refinement, canonical rep): the certified lower bound of N_ISD(f,h;z)
             exceeds B_cap on EVERY cell of K_G5 (per-cell ball contraction; any cell without a
             certificate BLOCKS exclusion — fail-closed). Whether
             Tier-2 runs is DECLARED together with the mu_lo measurement, strictly
             BEFORE any count under the cutoff; it may never be toggled — on or off — after any
             count exists (that would be a number chosen after seeing counts, the forbidden
             move, even though Tier-2 exclusions stay certificate-sound).
  Both tiers only EXCLUDE with certificates; retention is free. Over-count-safe direction of the
  whole pipeline is preserved: the admissible set can only be a superset of the physically
  censusable set, never a subset.
- WELL-DEFINEDNESS: N_ISD and gram_E are exact SL2(Z)-invariants (2x2 Gram determinants:
  (f,h) -> A(f,h), Gram -> A Gram A^T, det A = 1). Across monodromy the census transports the
  engine-canonical representative over the same pinned web set the cutoff uses, so cutoff and
  census see identical (representative, cell) pairs — necessity never crosses a frame. SOUNDNESS
  of Tier-1 for the class: a censused vacuum for canonical rep (f0,h0) at z* in K_G5 gives
  mu_lo^2 gram_E(f0,h0) <= lam(z*)^2 gram_E(f0,h0) <= Gram_H(z*) = N_ISD^2 <= B^2 <= 256, and
  gram_E_min <= gram_E(f0,h0). QED (exclusion can never kill a censusable class).

## 5. (a) WITNESS-KILL — the divergent family is excluded with diverging ISD norms

Without a norm condition the fixed-budget SL2-class count is INFINITE; the witness family
(genuinely inequivalent and Gauss-reduced) is
f = (p,q,0,0), h = (0,0,r,s), f.h = 0, pr + qs = b fixed, (p,q,r,s) unbounded. For each member:
    gram_E = (f.f)(h.h) - (f.h)^2 = (p^2+q^2)(r^2+s^2) >= max(p^2+q^2, r^2+s^2)
(the other factor is >= 1: f, h nonzero integer vectors). Along ANY infinite subfamily
gram_E -> infinity, because only finitely many integer quadruples satisfy
(p^2+q^2)(r^2+s^2) <= X for any finite X. Therefore
    min_{z in K_G5} N_ISD >= mu_lo * sqrt(gram_E) -> infinity:
the family's ISD norms DIVERGE. CLASS-LEVEL PRECISION: Tier-1 acts
on gram_E_min(class) (monodromy merges take the min), so the exact statement is: every family
member whose CLASS has gram_E_min > 256/mu_lo^2 is Tier-1 EXCLUDED; a member can evade only if
a monodromy merge supplies a representative with gram_E <= 256/mu_lo^2, and the sec 6 box bound
caps such evaders at finitely many classes TOTAL. Either way the family contributes only
finitely many admissible classes. The exact family that makes the fixed-budget SL2-class count
infinite is killed by the cutoff (`isd_ref.py` prints the b=1 members (1,0,1,s), s = 0, 10, 100,
with their gram_E, N_ISD and Tier-1 verdicts on a synthetic frame). PHYSICS READING: at fixed
integer budget b, Euclidean-orthogonal f, h of huge norm force the vacuum identity
2tB = x^T H x + t^2 h^T H h (sec 3 (v)) to balance an arbitrarily large positive form against
the fixed integer b — impossible; these pairs support no ISD vacuum anywhere on the censused
slice.

## 6. (b) FINITENESS — explicit bound, frontier closed

For any ADMISSIBLE class: gram_E_min <= 256/mu_lo^2, achieved by an SL2-canonical (Lagrange-
Gauss-reduced) representative with h.h <= f.f and (2 f.h)^2 <= (h.h)^2. Then
    gram_E >= (3/4)(f.f)(h.h)              [(f.h)^2 <= (h.h)^2/4 <= (f.f)(h.h)/4]
    (f.f)(h.h) >= f.f >= (f.f + h.h)/2     [h.h >= 1 integer; f.f >= h.h]
    => n2 := f.f + h.h <= (8/3) gram_E <= 2048/(3 mu_lo^2) =: N2_ISD(16).
Every admissible class therefore has an SL2-canonical representative inside
BOX(N2max = ceil(N2_ISD(16))) — the standard bound-completeness argument (every class has a
Gauss-reduced representative; Gauss reduction does not increase n2; the box enumerates all
pairs with n2 <= N2max; gram_E is SL2-invariant so no small-gram_E class can hide outside the
box) applies unchanged.
EXPLICIT BOUND:  #admissible classes <= #pairs in the box <= (2 sqrt(N2_ISD(16)) + 1)^8 < inf.
The enumeration FRONTIER is CLOSED: the enumeration ladder terminates at the SINGLE rung
N2max = ceil(N2_ISD(16)); no class with a censusable ISD vacuum lies beyond it, and the class
count is not enumeration-box-dominated. Positive-definiteness input
(mu_lo > 0): weight-3 Hodge-Riemann via the sec 2 explicit rank-4 form (the weight-4 involution
mechanism of the fourfold case does not transfer), measured under the
sec 4 HARD GATE mu_lo > 0: if the mu_lo measurement fails the gate, Tier-1 is void and
NO finite box is claimed until refinement restores it — fail-closed, never a false bound.
PER-BUDGET NOTE: the cutoff uses the CAP B_cap = 16 once, for the candidate set; per-vacuum
budgets stay recorded in absolute integer units by the nested census — no per-bin re-admission
ever.

## 7. (c) COMPUTABILITY — kernels and costs

- Inputs per cell: the transported 4-vector balls Sigma Pi, Sigma D_z Pi. MATRIX transport ONLY
  (point-values need only matrix transport; restricted scalar-elimination is refused — it
  blows up). These are exactly what the AESZ22/AESZ118 Picard-Fuchs webs carry; N0, N1 and every
  H entry are O(1) 4-vector contractions on them. Ball arithmetic end-to-end.
- TIER-1 cost: one-time mu_lo scan (4x4 interval eigenvalue bound per cell x #cells(K_G5)), then
  per-class EXACT INTEGER gram_E — no transport, no balls, integer operations only.
- TIER-2 cost: three ball contractions per (candidate, cell) on point-value period balls; the
  contraction REUSES the web balls with no new transport legs, so the per-row cost of the web
  itself (seconds per row on these operators) is an upper anchor for the per-(candidate, cell)
  cost.
- Quantities measured and recorded strictly BEFORE any G5 count under the cutoff:
   (1) the K_G5 cell-set hash (with excision list) — the pinned evaluation set;
   (2) mu_lo, certified interval, ball arithmetic end-to-end;
   (3) N2_ISD(16) = 2048/(3 mu_lo^2) and N2max = ceil(N2_ISD(16)) — pure arithmetic printed with (2);
   (4) Tier-1 scan rate + Tier-2 per-(candidate,cell) rate (short single-core timing probes).
  NONE of these four is a count.

## 8. What this derivation does NOT do (scope fence)

- Authorizes NO counting; counting under the cutoff is a separate later step, valid only once
  the formula and the measurements (1)-(3) are fixed.
- Changes NO census definition or pinned number: eps0 = 1e-3, B_max = 16, N_min = 10, the
  z-chart metric, the frozen special set S_frozen, the verdict classes — all untouched. The
  cutoff bounds ONLY the flux CANDIDATE SET, never distances or verdicts.
- Prior-art honesty: the ISD inequality is standard GKP physics;
  what is new here is the certified tau-eliminated per-class test, the exact-necessity pin to the
  census web set, and the closed finite box — not a new vacuum condition.

--- END ---

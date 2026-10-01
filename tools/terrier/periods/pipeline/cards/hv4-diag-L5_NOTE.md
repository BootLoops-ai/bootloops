# hv4-diag-L5 card — derivation note
SCOPE: everything here is DIAG-SLICE (S6-symmetric locus phi^1=..=phi^6 of the 6-parameter
Hulek-Verrill fourfold family HV4, arXiv:2404.12422). The full 6-parameter family is not
covered by this card.

1. SERIES (two-engine): diagonal fundamental-period coefficients c_n = 6-letter abelian-square
   counts (eq (5.18) restricted to the diagonal; the invariant slot of the 29-frame). Engine A =
   letter-adjunction convolution c^(k+1)(n) = sum_j C(n,j)^2 c^(k)(j) with its own Pascal
   triangle (the add-one-variable torus-period construction). Engine B = INDEPENDENT: direct A3
   composition enumeration (math.comb) + exact 3+3 / 3+1 multinomial splits. Full agreement
   n = 0..400 on c6, Domb c4, AND c3; Engine C direct 6-composition enumeration anchors n <= 40.
   Domb anchor 1,4,28,256,2716 == the printed K3 diagonal fundamental period.
2. FRAME-SIDE RANK PIN (frame_rank5.py on the 29x29 integer frame matrices MSWAP, M6, T1..T6,
   SIGMA of arXiv:2404.12422 sec 5; the matrices are not included in the package, the result
   FRAME_RANK5.json is): the S6-invariant sector of the 29-frame (common +1-eigenspace of
   MSWAP + M6) is EXACTLY dim 5 = the 5 orbit sums; the T1..T6 product restricts to it and is
   MUM with (T-1)^4 != 0, (T-1)^5 = 0; restricted Sigma-Gram det 6480 = 2^4 3^4 5, signature
   (3,2) (recorded, not interpreted). M6 is taken in the (5.59)-as-displayed convention =
   transpose of the printed (5.58); M6 is a slot permutation, so the +1-eigenspace is
   convention-free (asserted in frame_rank5.py). This pins ORDER 5 for the invariant slice
   system from the frame data alone.
3. OPERATOR DERIVATION: unique-fit over Q — the (order 5, phi-deg 3, theta-deg <= 5) ansatz
   sum_j phi^j Q_j(theta) has nullspace EXACTLY dim 1 on rows m = 0..60 of the two-engine
   series; normalized (q_{0,5} = 1) it reproduces Q_0 = theta^5 (MUM discovered, not assumed)
   and matches the printed operator of Jockers-Kotlewski-Kuusela (arXiv:2312.07611; quoted in
   arXiv:2404.12422v2 as ~eq (5.60)) COEFFICIENT-EXACTLY. Recurrence convention gated on the
   Domb/L3 anchor before use. Minimality probes: NO annihilator at (order <= 4, deg <= 10)
   (55 unknowns, m = 0..90) and NONE at (order 5, deg <= 2). DERIVED-UNIQUE + CITED-MATCH.
4. BATTERY (battery_l5.py, all gates PASS, ~1 s): recurrence annihilation m = 0..400 (L5 on
   c6; L3 on Domb); MUM Frobenius regeneration == Engine A to n = 400; theta-right R-form jet
   gate (pfaffian payload row) to n = 400; symbol = 1 - 56 phi + 784 phi^2 - 2304 phi^3 =
   (1-4phi)(1-16phi)(1-36phi), roots == the three diagonal conifold points {1/36,1/16,1/4} of
   the HV4 discriminant; sym(L5) = sym(L3) * (1-36phi); leading D-form coeff = phi^5 * symbol
   => NO apparent singularities; local exponents EXACT:
   0: {0,0,0,0,0} (MUM, matches frame index-5); each conifold: {0, 1, 3/2, 2, 3} (CY4-node
   type, integer 0..3 + half-odd 3/2); inf: {1, 1, 3/2, 2, 2}. K3 sector: {0, 1/2, 1} at
   both 1/16 and 1/4; inf {1,1,1}; K3 fit also unique + matches printed L3
   (Verrill 1996 + JKK).
5. PFAFFIAN CONVENTIONS (matrix-transport route ONLY, scalar-refusal law): card ships the
   5x5 theta-companion A_num/den, den = symbol, basis (y, th y, .., th^4 y); d/dphi frame =
   A_num/(phi*den). route_choice receipts recorded: 'point-value' -> matrix-transport law;
   'scalar-operator' probe (5,3) = 20 slots << 20000 budget -> the PRINTED scalar is legal
   here (it is the citation object, not an elimination product). Transport NOT run: the card
   has no routes field; a transport run on this card must add routes and pass the transport
   battery first.
6. GAPS (stated): GAP-D1 — exact all-orders annihilation is CITED (JKK) +
   necessary-condition verified to n = 400 two-engine (397 surplus exact identities beyond
   the 24-unknown fit); no independent creative-telescoping/Griffiths-Dwork certificate is
   issued here. GAP-D2 — the conifold monodromies of the 5-frame are not closed-form
   derivable from the frame matrices; they need certified transport. Neither gap blocks the
   card's consumers (transport seeding, K3 exact sector).

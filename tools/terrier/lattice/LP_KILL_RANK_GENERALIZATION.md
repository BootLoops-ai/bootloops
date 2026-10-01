# lp_kill rank generalization — contract (c=5 machinery → general complement rank c)

Base machinery: the Venkov 2-design LP kill over the 23 rooted Niemeier
lattices, originally written for complement rank c=5, and the
complement-class enumeration that feeds it. The rank-agnostic module
`lp_kill_rank_agnostic.py` reads the complement rank from the input Gram;
nothing else about the machinery moves.

## Complement enumeration (the LP's input; c=4 instantiation)
1. Witness input: the full 4x4 Gram of T.
   Guards: symmetry, definiteness, even diagonal, TG[1,1]=2.
2. The coset stage (matsolvemod particular solution + integer
   kernel + augmented-Gram Fincke-Pohst, both w[1]=+-1 branches) is
   one routine cosetenum(Mrows, bb, nn) applied twice:
   x3 (2 constraints, 6-dim kernel), x4 (3 constraints, 5-dim kernel,
   6-dim augmented Gram). No other logic; primitivity test,
   matkerint complement, qflllgram key, qfisom dedupe, qfauto |Aut|,
   MASSSUM identical machinery (complement 4x4). (PARI/GP primitives.)
3. a=1 route only: x1 = a fixed norm-2 vector, using the transitivity of
   W(E8) on the 240 roots.
Completeness (Nikulin gluing + Mordell): for ANY K with q_K = -q_T the
Nikulin extension T+K is even unimodular positive definite of rank 8, hence
E8 (Mordell), so every class of the complement genus C appears among the
complements of primitive T-images in E8; Siegel mass(C) == sum 1/|Aut| over
the enumerated classes certifies that no class is missed.
The enumeration script is not part of this package; `lp_kill_rank_agnostic.py`
reads the resulting class list in the format below (mode
`<key> <listing> <out_path>`).

## Genus-listing input format (read by parse_hecke)
Pipe-separated lines; fields are numbered from 0 after splitting on `|`.
- `HCLS|<key>|<idx>|·|<aut>|·|[g11,g12,...;g21,...;...]` — one line per
  class: field 1 the genus key, field 2 the class index, field 4 |Aut(K)|,
  field 6 the Gram (rows separated by `;`, entries by `,`, in brackets).
- `HKEY|<key>|·|<det>|·|<mass>|·|<nclasses>|·|<sum 1/|Aut|>|·|<match>` —
  one line per genus: fields 3, 5, 7, 9, 11; `<match>` must be the string
  `true` (Siegel mass == sum 1/|Aut|) or run_genus refuses (MASS GATE FAIL).

## lp_kill_rank_agnostic.py (Venkov 2-design LP)
1. analyze_class: n = len(GK) (rank read from the input Gram).
2. NDIM, IJ, the constraint matrix, and the Farkas box are all sized from
   the rank (c=4: NDIM=4, IJ=10 symmetric index pairs, constraint matrix
   (NIJ+1) x N, Farkas box dim NIJ+1, cert index y[NIJ]).
3. Rank-independent by construction: shells law (0 < a <= 2D-1 in
   adj = D*G^-1), cap 2*floor(2D/a), S from norm-2 vectors, count row
   24h - s_K, H (18 Coxeter numbers of the 23 rooted Niemeiers),
   exact-Fractions Farkas recheck, output schema. The LP is rank-agnostic
   (Venkov design identity + trichotomy + cap bound hold for any complement
   rank c; only the ambient rank 24 and |R(M)| = 24h enter).
4. The s_K >= 2 Leech guard (no root-free complement class) is enforced
   inside analyze_class. Kills are exact rational Farkas certificates,
   rechecked in `Fraction` arithmetic; no floating-point step is
   load-bearing. An independent verifier should share no code with this
   module: its own exact det (Fraction Gaussian elimination; no numpy/float),
   own inverse, own Fincke-Pohst, rebuilding det/s_K/S/shells/caps from each
   Gram and rechecking every Farkas certificate exactly, with the same
   s_K >= 2 guard.
5. Modes: `control` runs the A1^6 must-fail control (c=6, det 64, s_K=12:
   strict h=2 KILLED; the relaxed LP with the zero-projection column
   restored and caps dropped is FEASIBLE with n0=36); `<key> <listing>
   <out_path>` runs every class of the listed genus against all 18 Coxeter
   numbers and writes the per-class record (perh, certs, feas_witness,
   all_h_killed; closed/closing_class at genus level). Exercised at
   complement ranks c=3, 4, 6, 7 and the original c=5.

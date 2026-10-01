# galois — guide

Tool page: https://bootloops.ai/tools/galois.html

KIND: package (core.py, vspace.py, sector21.py, svmap21.py, mzvring/ — the
relation-ring MEMBER subpackage)

DISAMBIGUATION: MOTIVIC Galois coaction here — NOT differential-Galois-group
computations, and NOT the PyPI `galois` finite-field library (name fence below).

PURPOSE: Motivic Galois coaction-cut engine for single-valued HPL ansatz spaces PLUS
the exact MZV relation ring it reduces against. Coaction side: computes the Galois
derivations D_m (Goncharov coproduct, right-strip (id x pi_m)Delta, sv kills f2) on
MZVs / H-words / Brown svHPLs in this package's pinned conventions, and turns
coaction-stability against the Galois closure of SOLVED family members into exact
homogeneous rows on ansatz coefficients — rank gain at zero boundary-data cost
(CONJECTURAL coaction principle: rows are ansatz-cut input; downstream adversarial
gates stay). Ring side (member `galois.mzvring`): exact double-shuffle
canonicalization of {zeta(n), mzv(a,b), mzv(a,b,c)} monomials over Q, depth<=3,
w<=23, NO PSLQ — per-weight Gauss elimination, survivors = the BK depth<=3 basis.

USE-WHEN:
- A polylog family closure is CONSISTENT-but-UNDERDETERMINED at current data depth:
  coaction rows cut the residual before any new boundary data is bought.
- You need Galois conjugates / f-alphabet coaction data of any MZV to w<=8 (e.g.
  D5 zeta(5,3) = -5 zeta3) or of svHPL combinations in I[z,w,0] conventions.
- D_m images (m odd 3..19) + sv projection of depth<=3 odd-index MZVs at odd weights
  11..21 -> `sector21`.
- You need the ALGEBRAIC single-valued map zeta_sv of a depth<=3 odd-index MZV (w<=21)
  as an EXACT formula — no associator numerics -> `svmap21` (Brown arXiv:1309.5309 eq
  (7.3) on the f-alphabet, orient-B transport; 61/61-regression + 27/27-depth-2
  hard-gated).
- Any expression in zetas / depth-2 / depth-3 MZVs (w<=23) needs a CANONICAL exact
  form: closure assembly, equality tests of exact candidates, symbolic elimination
  before/instead of any PSLQ fit -> `galois.mzvring`.
- You need the depth<=3 survivor basis at a weight (which mzv's are independent, which
  reduce) or a fresh ring build at higher wmax.

NOT-FOR: extended-alphabet (zz-letter) sv spaces (skipped by the feed parser);
elliptic families; depth>=4 MZVs (hard ValueError in the ring); alternating/Euler
sums; numeric evaluation (use formglue.form_oracle / tornheim — exact-only here); NOT
a closure authority — a coaction-assisted closure needs the usual exact solve +
held-out certification. w>8 needs OUTSIDE the depth<=3 odd-index sector remain
unbuilt.

INVOKE (public): put this repo's `tools/` dir on sys.path;
`from galois import core, vspace, sector21, mzvring, svmap21`.
- Derivations (core): `core.Dm_zpoly(zp, m)`, `core.Dm_hword(u, m)`,
  `core.Dm_calL(word, m)`, `core.Dm_vector(vec, wt, m)`; conjugation
  `core.conj_vector(vec)`; ansatz parse `core.element_vector(parsed_element)`.
- Allowed space (vspace): `vspace.collect_feed(...)` -> `vspace.saturate(vecs)` ->
  `vspace.fiber_annihilator(fib, wt, par)` -> functionals; rows = functional .
  D_m(element images) per LS block.
- Sector (sector21): `Dm_comp(c,m)` / `Dm_word(u,m)`, `Dm_vec`, `sv_proj`, `zsh(u)`,
  `sector_words(w)` / `ms_for(w)` / `load_tables(path)`; CLI `python3 sector21.py
  selftest|gates|oracle|sector|sv [--ring PKL --tables JSON --receipt J]`.
- SV map (svmap21): `sv_reduce_many(comps, orient='B')`; `phi_sym`/`phi_vec`,
  `sv_f(fv, orient)`, `derive_formulas(...)`, `coeffs_brown`/`formula_strings`; rebind
  rings via `svmap21.init_ring(path)` (NOT sector21's — see FOOTGUNS); CLI
  `python3 svmap21.py selftest|derive [--ring PKL --svalign DIR --dps N --comps ...]`.
- Ring (mzvring member): `R = mzvring.load_ring(path)` -> `R.cvec/R.cmono/R.csym`
  canonicalize; `R.survivors[w]`; vector algebra `dsum/dscale/dmul`, `znorm(n)`, words
  `comp_word/word_comp/shuffle`, `mweight`. Symbols: `('pi2',)`=pi^2, `('z',n)` n
  odd>=3, `('m2',a,b)`, `('m3',a,b,c)`.
- New ring tables: `python3 mzvring/build_ring.py --wmax N --out X.pkl [--resume PKL]
  [--checkpoint]` (checkpoint dumps <out>.partial after every weight — use it, builds
  are hours at w=23).

DATA (reference-table banks, not shipped — build them with build_ring.py or point the
env vars below at your own; the package FAILS CLOSED without them; every entry point
also takes explicit paths):
- `GALOIS_T0DEEP` -> dir holding mzv_formal.py + mzv_tables_w8.json +
  svt_tables_w8.json (the w<=8 table bank; not shipped).
- `GALOIS_RING_BANK` -> dir of reference ring pkls (ring17/ring21/ring23; build your own
  with build_ring.py — ~30 min at wmax=21, ~2 h at 23, single core).
- `GALOIS_CAMPAIGN_BANK` -> reference-table root with subdirs galois/
  (SECTOR_TABLES.json) and svalign/ (SV_FORMULAS.json + the ~700d numeric bank).
- `GALOIS_ANC_DIR` -> the solved-feed anc files = the PUBLISHED ancillary files of
  arXiv:2607.11645 (publicly downloadable).
MZV convention THROUGHOUT: mzv(a1,..,ak) = sum_{m1>m2>..>mk>=1} prod mi^-ai,
OUTER index first, a1>=2; zeta(2k) -> rational*pi2^k.

OUTPUTS: exact Fraction rows {col: q} (homogeneous, rhs 0) labeled
`coact:b<off>:D<m>w<wt>:f<i>`; canonical vectors {monomial: Fraction}.

ENV: pure python3 + sympy + mpmath (+ gmpy2 for eliminations); single-core; engine
RSS < 1G, ring RSS < 100 MB through wmax=23. Ring build walls (measured): w=21
1056 s, w=22 1297 s, wmax=23 ~2 h fresh.

GATES (receipts not included):
- Coaction engine: G1 Goncharov==double-shuffle (1530 checks), G2 sv-closure (939
  images), G3 parity eigenstate lock (304 elements), G4 tornheim numeric oracle, G5
  leave-one-out positive control (0/32 violations; plants violate 6/135 odd, 119/169
  even). Re-run the leave-one-out control whenever the solved feed changes.
- ARCHIVE CERTIFICATION: 161 stored exact-Q items, 8 families — 0/161 coaction
  violations, 4/4 planted perturbations violate.
- sector21: selftest gate a 255 values + 765 D-images vs core's table route EXACT; 6
  literature values (Brown 1102.1310/1102.1312, arXiv:2511.15883); 6 FORM 200d spots
  (min 196.4d); full gates a-d PASS incl. 31 identity checks + parity/reversal locks.
- Ring: BK depth-3 dims w=19:5 / 21:6 / 23:8 by FRESH per-weight elimination;
  survivor-set equality asserted vs bank; exact Laurent-slot and closed-form
  regressions vs stored results; 1298/1298 checks vs 150d numerics.
- svmap21 (`selftest --receipt PATH` writes the receipt on full PASS; no receipt is
  vendored — the gates re-run against your own bank):
  pin — 4 literature phi locks + Brown's printed (7.3) depth-3 example VERBATIM +
  sv(f2)=0; depth2 — 27/27 all-odd depth-2 pairs vs the W-associator sv bank at 460d
  (min 480.7d) + exact zeta_sv(5,3) = -10 z3 z5 lock; HARD GATE regress61 — 61/61
  reference SV_FORMULAS reproduced EXACTLY under orient B (orient A 2/61 = the
  palindromic keys only); extcheck — 22/22 EXT bank == in-process re-derivation.

FOOTGUNS:
- PyPI NAME FENCE: an unrelated `galois` (GF/finite-field) library exists; with this
  repo's tools/ ahead of site-packages on sys.path, `import galois` gets THIS package.
  If you need both, import the PyPI one first or by explicit file path.
- V is exactly as strong as the solved feed: a failing control = V too small — enlarge
  the feed (more solved siblings), NEVER hand-patch fibers.
- Convention traps are locked by G3/gate d: if you change word/f-token conventions
  upstream, re-run before trusting any row. pi_m for m>=11 is basis-DEPENDENT — locked
  to the ring's survivors basis.
- svmap21 MIRROR-ORIENTATION PIN: Brown's verbatim (7.3) lives in the PREFIX-strip
  f-model (orient A); this package's f-model is TRAILING-strip, so the correct transport
  is orient B — reversal on the LEFT deconcatenation factor — EMPIRICALLY PINNED 61/61
  vs 2/61 (A == B at depth <= 2 and on palindromic keys, so shallow checks CANNOT
  adjudicate). Never "fix" an sv mismatch by flipping orientation — re-run regress61.
- svmap21 conventions: the pure-f2 coefficient of even-weight survivor generators := 0
  is the ONE non-recursive input (validated through every stored z2^j coefficient by
  the 61/61 gate). The ring is depth<=3-truncated: a target needing depth>=4 content
  HARD-FAILS the rectangular solve (reported, never smoothed). Rebind rings via
  svmap21.init_ring, not sector21.init_ring — the phi/monomial caches are ring-keyed
  and cleared only by the former.
- Ring tables are built weight-ASCENDING; load_ring skips __init__ — never call
  Ring(wmax) casually (full rebuild). csym on an unknown symbol returns it as a basis
  vector rather than erroring — coverage discipline is the CONSUMER's
  (sector21.sym_vec hard-stops; do the same).
- Elimination keeps lexicographically LARGEST all-odd no-1 args (m2(5,3) survives at
  w=8); re-derive for a different basis convention, never post-rotate canon tables by
  hand.
- pkls are the expensive artifact: never regenerate over the bank in place; build to a
  new path, validate, then bank.

BATTERY (PARTIAL): pure ring
algebra runs with no bank (znorm/shuffle smoke PASS); every reference-data entry point
fails CLOSED with a named env-var refusal when unset (core import, load_ring,
load_tables all verified); with the env vars pointed at the banks, the
package loads the ring (wmax=21, BK dims match) and the w<=8 tables. Runnable without any reference data:
mzvring algebra + build_ring.py (fresh ring builds) + core/sector21/svmap21 once you
build/point the banks per DATA above.

RELATED: rows feed rankscreen (verdict engine, ships in tools/); formglue.tornheim +
form_oracle + lockpick pools = independent MZV oracles (ring validation + G4);
lockpick (what the ring makes unnecessary at depth<=3 — PSLQ demotes to held-out
gate).

CREDIT: the coaction principle is Francis Brown's [BrCG], sharpened by Erik Panzer and
Oliver Schnetz [PS] and, for amplitudes, by Caron-Huot, Dixon, Dulat, von Hippel, McLeod
and Papathanasiou [CDDvHMP]; the coproduct is Goncharov's [Gon] in Duhr's formulation
[Duhr12]; motivic decomposition and the f-alphabet follow Brown [BrDec, BrMTM];
single-valued MZVs/HPLs and the sv map follow Brown [BrSV] §7.2; conventions follow
Xuhang Jiang [Jiang]; the relation ring uses the regularized double-shuffle framework of
Ihara, Kaneko and Zagier [IKZ], Hoffman's relations [Hoffman] and the Broadhurst–Kreimer
depth-graded dimensions [BK]; literature checks against Brown [BrDec] and Schlotterer,
Sohnle and Tao [SST25]; 200-digit spot checks by FORM [FORM4]. Bracketed keys resolve in
REFERENCES.md at the repository root.
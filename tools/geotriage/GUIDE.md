# GeoTriage — guide

GeoTriage — `tools/geotriage` — https://bootloops.ai/tools/geotriage.html

KIND: package (geotriage.py, test_geotriage.py, graph_specs/
fixtures, plus the exactj and fiberstack members). This package is also the code
anchor of the p-adic Frobenius
boundary-constant recipe (see toolkit/ours/RECIPES.md, the `padic-frobenius` recipe).

PURPOSE: One-command per-graph maximal-cut geometry + value-fittability triage: graph
→ {genus, j-invariant, Kodaira type, modular level, CM flag} → closure route ∈
{A,B,C,D,E} with mandatory honesty flags. Decides value-fit (A), value-fit + 1
boundary constant (B), or DE-transport mandatory (C) BEFORE any compute is committed.

USE-WHEN:
- New graph/sector needs its geometry, curve, and triage class →
  `geotriage` (the sector_symanzik stratifier front-end is not included in this repo).
- Choosing closure route — value-fit or DE-transport? → A/B/C value-fittability triage
  (decidable from period-matrix Jordan structure).
- Elliptic/K3 sector needs its boundary constant without an integral solve →
  arithmetic-boundary L-value bootstrap: constant from the PF operator's p-adic
  Frobenius (= Hecke a_p). Decisive for elliptic; for K3, certified on the
  equal-mass banana (eta-quotient path) and a PROBE elsewhere (exact
  background-centred point count + Livné χ_D scan, `certified: false`). (The
  probe scripts are not included in this repo; the classifier here is the
  reference implementation.)

NOT-FOR: Does not build the PF operator (`pf_order_certified` is False unless one was
actually built elsewhere). A/B/C criterion: forward direction proven, converse
CONJECTURED — a C verdict is firm, an A/B verdict is a route hypothesis to gate.
L-value bootstrap PARTIAL beyond elliptic. Route decision table only covers the
classified varieties (genus 0/1, K3, CY3+, genus≥2→E hard limit).

INVOKE (public): CLI `python3 geotriage.py graph_specs/<name>.json --out
<name>_classify.json`; library `from geotriage import classify;
rep = classify("graph_specs/sunrise_eq.json")` → `rep['route']`, `rep['elliptic']`
({j, kodaira, level_N, on_X1N, congruence_group, ...}), `rep['honesty']`. Tests:
`python3 test_geotriage.py` (5 known cases).

INPUTS: JSON graph-spec schema shared with the tools/landau-alphabet package, plus optional
`"maxcut": {"poly","vars","modulus","source"}` (proxy mode — feed a known Baikov
maximal-cut polynomial directly for graphs whose loop-by-loop localization is computed elsewhere) and `"fiber_point": {...}` (kinematic specialization for the CM-j /
a_p test).

OUTPUTS: JSON verdict: route ∈ {A,B,C,D,E} + ring (MZV / Chowla–Selberg Γ /
Γ₁(N)+L(χ_D,k) / ℚ(√D)+L(f3,s) / transport). Decision tree: genus 0→A; genus 1
j∈CM-13→A; genus 1 on X₁(N) N≤10,12→B; genus 1 non-congruence→C; K3 rank-2 (CM,
Livné)→B; K3 rank≥3→C/D; CY3+→D; genus≥2→E. Honesty block ALWAYS present:
`genus_method` ("exact" | "degree-bound" | "numeric-probe"), `cm_evidence` ((p, a_p,
χ_D(p)) tuples or "j∈CM-13"), `pf_order_certified`, `caveats` (every assumption
listed).

GATES: 51 tests in the classifier; 6-case known-answer battery (box1l polylog/route
A, sunrise etc.) plus synthetic-truth + mutation controls on the K3 probe. Read
`rep['honesty']` FIRST — a "degree-bound" or "numeric-probe" genus_method is not a
certificate, and a K3 `cm` block with `certified: false` is a PROBE verdict.
Cross-certify PF order with the annihilator/pf-rank tools before committing
serious compute to the route.

MEMBERS:
- exactj (`exactj.py`, `from exactj import exact_j_resultant`): exact minimal
  polynomial over ℚ of j(λ), λ = x₀², for an ALGEBRAIC fiber x₀ given by any
  irreducible integer minpoly (ascending coefficients) — two independent
  flint-exact routes COPAIRed per call (literal resultant Res_x(m, t·D−N) via
  evaluation + exact Lagrange interpolation vs multiplication-matrix minpoly on
  ℚ[x]/(m)); fails loudly on any disagreement. Free char-0 non-CM certificate:
  `j_minpoly_lc != 1` on the primitive minpoly (a CM j is an algebraic
  integer); necessary-not-sufficient the other way — never read lc = 1 as a CM
  claim. λ ∈ {0,1} poles j (the engine raises; fence in the caller's menu); a
  reducible minpoly returns AMBIGUOUS-REDUCIBLE-MINPOLY — split the fiber,
  never average. Rational-j outputs are checked against the same CM_J_TABLE the
  classifier uses. Requires python-flint (validated with 0.8.0); the battery
  skips these legs by name when it is absent.
- fiberstack (`fiberstack.py`): exact arithmetic-fiber instruments —
  `chi_table`/`legendre_ap_fast` (Legendre-curve a_p with the quadratic
  character tabled once per prime; identical mathematics to the naive character
  sum), `ap2_charsum` (a_{p²} over F_p² by norm-character sum; the exact
  identity a_{p²} = a_p² − 2p is a built-in COPAIR whenever λ ∈ F_p),
  `rational_heights_menu` (deterministic height-ordered rational point menu).
  ap2_charsum is O(p²) — control fibers and short menus only; the caller
  excludes bad reduction λ ≡ 0,1 (mod p).
- Method note (fibers of field degree > 1): take the roots r of the fiber's
  minimal polynomial mod p; each root gives λ = r² in F_p and a split-prime a_p
  read via `legendre_ap_fast`; at inert primes (no root mod p) read
  supersingularity via `ap2_charsum` (p | a_{p²}, valid for p ≥ 5). This closes
  split-only supersingularity scans two-sided.

FOOTGUNS:
- Genus can be OVER-reported by degree-bound methods, never under — treat "higher" as
  "certify via the PF-order tools", not as a final verdict.
- A/B/C converse unproven: value-fit failure at class-A prediction is evidence, not
  contradiction; escalate to C machinery rather than grinding PSLQ.

BATTERY (SELF-CONTAINED, honored): `python3 test_geotriage.py` run from the
package in a scratch cwd → all 6 known-answer cases PASS (box1l→A,
sunrise_eq→B on Γ₁(6), icc→C non-congruence, banana3l_eq→B K3-CM Livné D=-15,
crossedbox3l→A, banana3l_1112→B K3-CM probe D=-3 cross-checked against the
Hecke trace), rc=0, plus unit/control legs (wired a_p counter patterns;
synthetic-truth centring bank with mutation, degenerate-fiber, and no-lock
controls) and member legs (exactj: x²+1→j=1728/D=-4, x⁴−x²+1→j=0/D=-3, and a
deg-2 cyclotomic COPAIRed against an independent pure-Fraction ℚ(ω) oracle
with the non-CM lc certificate — skipped by name without python-flint;
fiberstack: a_p engine vs naive point count, the a_{p²}=a_p²−2p identity with
a mutation control, menu determinism). Fixture names (icc = ice-cream-cone, crossedbox3l = 3-loop crossed
box, ...) name the physics graphs of the validated battery.

CREDIT: maximal cuts in the Baikov representation (P. A. Baikov 1996, Phys. Lett. B
385:404; loop-by-loop form of Frellesvig & Papadopoulos 2017, JHEP 04:083,
arXiv:1701.07356); fibre types after Kodaira; the four-cusp rational elliptic surfaces
are Beauville's list (A. Beauville 1982, C. R. Acad. Sci. Paris Sér. I 294:657) with
index data from Sebbar's classification of torsion-free genus-zero congruence subgroups
(A. Sebbar 2001, Duke Math. J. 110:377); the CM test for rank-2 K3 motives is Livné's
theorem (R. Livné 1995, Israel J. Math. 92:149); the sunrise/banana elliptic and K3
identifications follow Bloch & Vanhove (arXiv:1309.5865, J. Number Theory 148 (2015)
328) and Bloch, Kerr & Vanhove (arXiv:1406.2664, Compos. Math. 151 (2015) 2329;
arXiv:1601.08181). Exact arithmetic: FLINT via python-flint.

RELATED: tools/maxcut (the certified maxcut-geometry engine this front-ends);
annihilator (PF-order certification); upgrades/Eichler.jl (class-B/C period side);
tools/landau-alphabet (shared graph-spec schema).

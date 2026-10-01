# cosmoflow — FRW/dS wavefunction integrands and letter alphabets for an arbitrary site graph

## What
Letter alphabets of FRW/dS wavefunction coefficients and correlators for ANY
connected site graph — chains, rings, stars, multi-loop polygons, parallel-edge
bubbles — through one front door, `alphabet_graph(nv, edges)` (CLI kind
`graph`): connected-subgraph facet forms `q_g` (deleted loop edges at the
DERIVED c=2, see below), folded and disc letter classes, ansatzer schema,
per-letter sources, honesty block. Plus the Cayley–Menger Baikov `B(y;X)` of
the 1-loop n-site graph (`polytope.baikov_B`, any n_s) and its disc/Källén
letter block. The 1-loop triangle is the WORKED EXAMPLE (`examples/triangle/`):
K(m)-split nested-quad numerical oracle for the elliptic-subsector masters,
the two-site-chain twisted oracle (`F_twosite*`), the maxcut residue-period
sampler + symbolic L₂ proof.

## When to use
You have a site graph and want its candidate symbol alphabet:
`python3 -m cosmoflow.alphabet graph --nv 4 --edges 0-1,1-2,2-3,3-0,0-2` (run
from `tools/`), or `from cosmoflow import alphabet_graph`. Named special cases
`alphabet_tree_chain(n_s)` / `alphabet_loop_ngon(n_s)` (validated censuses
5/10/17). `polytope.baikov_B(n_s,…)` for any 1-loop n-gon integrand.
Computing/gate-checking the elliptic sector of the FRW 1-loop triangle
(arXiv:2408.16386): `cosmoflow.examples.triangle` — `oracle.eval_master` for
n_s=3 (triangle masters) and n_s=2 (TREE two-site chain F/F0/F1), `maxcut.*`
for n_s=3.

## Status
`alphabet_graph` front door: censuses validated (5/10/17) and battery-checked
letter-for-letter against the tree-chain / loop-n-gon special cases. Worked
example: n_s=3 fully validated (see below); n_s=2 tree-chain twisted oracle
validated, 5/5 gates PASS incl. the 28-point receipt set byte-exact. c-default: q defaults to the DERIVED c=2
(c=1 by argument).

## Caveats
Alphabets: connected graphs only; the Baikov disc/Källén block is 1-loop-n-gon
only. Numerics are the worked example's: n_s≥4 oracle: `NotImplementedError` stub. The n_s=2 evaluator's object is the
TREE two-site chain F(X₁,X₂,Y;ε), −1<ε≤0 — NOT the one-loop bubble masters
(that positive control runs through `polytope` + direct quadrature). `alphabet.py` emissions carry `complete:false` until the
connection-denominator cross-check has been done for that graph. Residue cycle
is NOT an L₂-solution (carries lower-block inhomogeneity). Γ non-compact;
e₁,e₂,e₄,e₈,e₉ UV-div at ε=0 (and 7 of the 9 block components have divergent
outer integrals even at ε₀=−1/20 — boundary values need tail machinery not shipped in this package).

## Validated
- Periods: L₂[ϖ₀]=0 symbolic + 95d; Wronskian π/(2λD₄) 80d; period_residue
  100d↔140d self-consistent.
- **L₂ + full 16×16 connection derived INDEPENDENTLY from the integrand:**
  product-twist Griffiths–Dwork mod p; χ rank 16=7+9; blind L₄
  → Ore right-division L₄=M₂·L₂ exact over ℚ(λ) at a∈{2,3}, both c;
  A(λ;ε₀=−1/20) 246/256 single-prime Wang + 10 entries 6-prime coefficient
  CRT; held-out 7th-prime audit 256/256.
- **Finite-ε transport gate:** e6 31.84d / e7 31.94d at
  (a,λ,ε₀)=(2,7/10,−1/20) vs fresh oracle never used in any fit; control
  31.02d @ λ=3/5; boundary 9-vector ≥31.05d; round trip 114.84d; θ-correction
  certified 37.0d.
- **ε=0 oracle ≥30d at (3,3/10):** e6/e7 certified 33.47–33.94d
  (2-prec 37.9–38.2d, R/Nth-independence 33.5–33.9d).
- **Positive control:** 2-site bubble MPL master vs closed-Li₂ form 51.33d.
- `examples/triangle/oracle.py` quadrature features: an inner K-split ts-fallback (corner
  Nch-bias), the `disk_corner()` polar carve-out of the codim-2 log point with
  an adaptive fallback degree, and the `_pinch_y31()` interior C1-pinch split
  for X1>X3 kinematics.

## Limitations
- ε=0 oracle at p1=(2,7/10) (X1>X3 kinematics): the shipped reference values
  are certified to 19.36–19.85d only; they were produced at ladder depth mm=5,
  where the depth-cross floors at the pinch-stall noise. The interior-pinch
  split in `examples/triangle/oracle.py` removes the stall; ≥30d at this kinematic class needs a
  depth mm=6 evaluation.
- e₃,e₅ capped at 8.5–10.3d (outer 1/y² tail).
- ε-factorized/canonical basis and the ε→0 limit of the transport: not
  constructed (ε→0 is degenerate: O(ε) couplings vs O(1/ε) masters).

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).

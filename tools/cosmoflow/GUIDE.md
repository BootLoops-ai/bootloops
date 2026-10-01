# cosmoflow — guide

Tool page: https://bootloops.ai/tools/cosmoflow.html

KIND: package. General front door: `alphabet.py` (`alphabet_graph` + CLI kind
`graph`), `polytope.py` (`baikov_B`, `q_subgraph`), `cosmoflow.py` CLI, `selftest.py`
(each run writes its receipt to `SELFTEST.json` beside it). Worked example:
`examples/triangle/` (`oracle.py`, `maxcut.py`) — the 1-loop 3-site triangle of
arXiv:2408.16386 and its tree positive control, the two-site chain.

PURPOSE: FRW/dS wavefunction and correlator integrands and their letter alphabets for
an ARBITRARY site graph — chains, rings, stars, multi-loop polygons, parallel-edge
bubbles. You bring the graph (number of sites + edge list); `alphabet_graph(nv, edges)`
emits the candidate symbol alphabet in ansatzer schema: the connected-subgraph facet
forms q_g (= OFPT energy denominators; once-cut edges coefficient 1, deleted loop edges
the DERIVED c=2), the folded (X_i−Y_e) class, the disc class (bare energies of the
cycle edges), per-letter `source` tags, the loop number, and an honesty block
`complete:false`. For 1-loop n-gons `polytope.baikov_B(n_s,…)` builds the Cayley–Menger
Baikov polynomial of the integrand I = ∫_Γ B^{(d−n_s−1)/2}·N/∏q_g^τ, Γ = {y≥0, B≥0},
and the Baikov disc/Källén letter block. The triangle (n_s=3) is the WORKED EXAMPLE:
its K(m)-split numerical oracle for the elliptic-subsector masters, the two-site-chain
twisted oracle and the maxcut residue-period sampler live in `examples/triangle/` and
show how a specific graph is carried through to gated numbers.

USE-WHEN:
- Letter alphabet of ANY connected site graph for ansatzer pricing / symbol fits →
  `alphabet.alphabet_graph(nv, edges, c=2, site_names=, edge_names=, with_disc=)` or
  `python3 -m cosmoflow.alphabet graph --nv N --edges u-v,u-v,...` (0-based sites,
  repeated pairs = parallel edges). Named special cases: `alphabet_tree_chain(n_s)`
  (path), `alphabet_loop_ngon(n_s)` (n-cycle) — thin wrappers, battery-checked
  letter-for-letter against the front door. Validated facet censuses: tree n2 = 5
  letters (names l_EL, l_ER, l_ET, l_FL, l_FR; ansatzer residuals 1/3/8), 3-gon = 10
  facets, 4-gon = 17 facets; e.g. 4-site star = 11 facets / 21 letters, box+diagonal
  (2 loops) = 33 facets. Output is a CANDIDATE set until the connection-denominator
  cross-check has been run for that graph.
- Any 1-loop n-site polygon cosmological correlator: `polytope.baikov_B(n_s,…)`
  integrand construction (n_s=3 triangle validated; n_s=2 bubble MPL; n_s≥4
  conjectured elliptic+, structure only).
- WORKED EXAMPLE (`cosmoflow.examples.triangle`): two-site TREE-chain coefficient
  F(X₁,X₂,Y;ε) at −1<ε≤0 (ε=0 = the log-rational dS branch) → `oracle.F_twosite` /
  `F1_twosite` (ε¹ Taylor leg) / `F_twosite_certified` (two-precision × two-level), or
  `eval_master(2, 'F'|'F0'|'F1', X1, X2, Y, eps=…)`; closed-form inner integral,
  πcsc(πε)(A^ε−B^ε) kernel; ~0.03–0.07 s/pt at 40–60d; print form
  `F_twosite_print342` is GATE-ONLY. Triangle masters e1..e9 → `eval_master(3, key, a,
  lam, c)`. The quadrature devices below are the reusable lessons of the example:
- Nested √-quadrature stalls on an isolated INTERIOR log-singularity → K(m) log-split
  (`examples/triangle/oracle.py` inner_vec K-split).
- 2D nested quadrature caps at 16–22d — codim-2 point-log → `disk_corner()` polar
  carve-out.
- Nested-quad ladder OSCILLATES or the middle layer stalls ~1e-22 at ALL outer nodes —
  interior complex pinch (X1>X3 kinematics) → `_pinch_y31()` closed-form pinch locus +
  split at Re(singularity).
- Fast maxcut check → `examples/triangle/maxcut.py` residue-period sampler at
  q_{G12}=0 (101d / 0.12 s).

NOT-FOR: the site graph must be connected (a disconnected or malformed edge list raises ValueError); the Baikov disc/Källén block
exists for 1-loop n-gons only (other graphs: `with_disc=False`, else `ValueError`).
Numerics beyond the worked example are yours to build: n_s≥4 oracle =
`NotImplementedError` stub; the one-LOOP bubble MASTERS are
not the n_s=2 evaluator's object (the tree-chain F is). The residue cycle is NOT an
L₂-solution (carries lower-block inhomogeneity). ε-factorized basis + ε→0 limit NOT
constructed. Γ is non-compact; only deg N<τ masters are finite at ε=0; boundary values
of several block components need tail machinery not shipped in this package.
`alphabet.py` output is a CANDIDATE set (`complete:false`) until the
connection-denominator cross-check has been done for the graph at hand.

INVOKE: import from this directory's parent (`tools/`). Front door:
`from cosmoflow import alphabet_graph; alphabet_graph(4, [(0,1),(1,2),(2,3),(3,0),(0,2)])`
or the CLI
`python3 -m cosmoflow.alphabet graph --nv 4 --edges 0-1,1-2,2-3,3-0,0-2 [--c 2] [--sites a,b,c,d] [--enames ...] [--disc] [--out F.json]`
(`--disc` attaches the Baikov disc/Källén block, 1-loop n-gon graphs only); special
cases `python3 -m cosmoflow.alphabet tree n_s` / `loop n_s [--c 2] [--no-disc]`;
`polytope.baikov_B(n_s, ys, Xs)`, `polytope.q_subgraph`; `cosmoflow.py` CLI (`build n_s`
general; `oracle|maxcut|periods` drive the worked example). Worked example:
`from cosmoflow.examples.triangle import oracle, maxcut` (`oracle.eval_master`,
`oracle.F_twosite*`, `maxcut.period_residue`, `maxcut.varpi0/1`,
`maxcut.prove_L2_varpi0_symbolic`). `selftest.py` = the bounded battery (heavy oracle
probe only with the --deep flag); regenerates `SELFTEST.json`.

INPUTS: (n_s, kinematics (a, λ, X_i), ε₀); physical strip 1/(a+1) < λ < 1/(a−1);
disc_{z₁}B = 4·λ_K(X²)·λ_K(y²,X₃²) nested-Källén.

OUTPUTS: integrand/period objects, master values with quadrature self-certs; periods
ϖ₀,₁ = K, K′(K²)/√(1−(a+1)²λ²); Wronskian π/(2λD₄).

GATES: periods symbolic proof + 95d vs the literature L₂; e6/e7 finite-ε transport
gated ~31.8–31.9d vs a fresh oracle never used in any fit + positive control (2-site
bubble MPL 51.33d); ε=0 oracle certified 33.5–33.9d at (3,3/10) (2-precision +
R/Nth-independence). Full 16×16 connection derived INDEPENDENTLY from the integrand
(mod-p Griffiths–Dwork + Ore right-division L₄=M₂·L₂; held-out 7th-prime audit
256/256). Reproduce the 2-precision + fresh-oracle pattern for any new kinematic
point. Gate set: 5/5 gates PASS — the 28-point receipt set byte-exact,
≥46d vs print arXiv:2312.05303 eq.(3.42), ε¹ certs 45d, alphabet 5/10/17
censuses + ansatzer 1/3/8, maxcut n_s=3 sampler byte-identical.
Battery note: `selftest.py` (no arguments) is the bounded battery
(~5 s from a cold clone; exit 0/1; regenerates `SELFTEST.json`): two-site F vs
the pinned print form eq.(3.42) (51d measured), certified ε⁰/ε¹ legs (50d),
alphabet censuses 5/10/17, the front-door consistency row (`alphabet_graph` on the
3-cycle = `alphabet_loop_ngon(3)` and on the 3-chain = `alphabet_tree_chain(3)`
letter-for-letter, plus the CLI `graph` kind), symbolic L₂[ϖ₀]=0 proof, period Wronskian (41d),
maxcut sampler two-precision self-agreement (31d), interior-pinch locus vs the
symbolic z1⁰ coefficient of B + X1>X3 class discrimination (25d floor). The triangle e6 oracle
reproduction is a HEAVY leg (~448 s measured for its probe alone) and runs only
via `selftest.py --deep`; without that flag the battery records it as a named
SKIP in `SELFTEST.json` (the oracle is correct but slow at probe settings).

FOOTGUNS:
- The oracle needs all three quadrature features (inner ts-fallback, disk_corner
  carve-out, interior-pinch split) to reach gate depth; a plain nested quadrature
  without them is only ~8d.
- At (2,7/10)-class kinematics (X1>X3) the shipped ε=0 oracle values are certified to
  19.36–19.85d only (produced at ladder depth mm=5; the interior-pinch locus is
  battery-verified, and ≥30d at this class needs a depth mm=6 evaluation). Check
  kinematics class before quoting depth.
- A stalled ladder's d(mo_k, mo_{k+1}) is NOT an error bound — do NOT "correct"
  against a deeper single-panel value; split at the pinch instead.
- θ-rate in disk_corner is point-dependent (0.23–0.51 d/node measured) — size Nth
  from the measured Nth vs Nth/2 pair, never assume; the disk term is ADDITIVE so
  Nth-refinement is a valid post-hoc correction.
- Inner ts-fallback threshold and adaptive fallback degree are tuned (0.20, deg 7/6)
  — deg 6 caps ~21d at small ratio; don't lower them.
- c IS DERIVED — c=2, three independent exact routes; the reference paper's printed
  triangle c=1 is a typo internally contradicted three times in arXiv:2408.16386
  itself. `polytope.q_G12`/`q_subgraph`/`alphabet_graph` defaults are c=2 (see polytope.py); c=1 stays reachable by explicit argument. Quote c=2 rows as the
  physical ones (e.g. e6 = 0.2958830539…, e7 = 1.369167258… at (3,3/10)).
- e3/e5 are only 8.5–10.3d — not gate-grade.

RELATED: tools/ansatzer in this repo (alphabet schema consumer). NOTE
`cosmoflow.examples.triangle.maxcut` is the worked example's module, unrelated to the
`maxcut` package — always import it by that full path. Eichler/elliptic evaluators (upgrades/ in this
repo) for downstream closure.

CREDIT: Target integrals and L₂ reference: arXiv:2408.16386 (FRW wavefunction
correlators). Adapter, K(m)-split, disk carve-out, pinch split are in-house.

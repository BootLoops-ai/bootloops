# maxcut — guide

Tool page: https://bootloops.ai/tools/maxcut.html

KIND: package

PURPOSE: Maximal-cut kit: classify cut geometry (polylog/elliptic/K3/CY3/higher) from
decidable algebraic invariants; canonical UT rotation + alphabet from Baikov; N×N maxcut
DE assembly with the correct on-shell/off-shell derivative operators; d/dm² DE rows for
off-shell legs; weight≥5 GF(p) symbol cascades. Also home of the UT/LS sampling
discipline and the Row-0-only-T canonicity diagnostic.

USE-WHEN:
- Top-maxcut genus/geometry verdict needed from Baikov, no AMFlow/IBP/numerics →
  `decisive` (Track A: per-ISP degree profile, carrier-slice genus + j on ≥3 slices,
  disc-chain on deg-2 ISPs only; Track B: Krylov min PF order + reducibility from a
  period series; mod-p LP route via flat `tools/pf_rank.jl`).
- 5-point maxcut needs canonical rotation + alphabet → `ut_rotation` (Aomoto/KZ
  connection from the Baikov fundamental matrix, partial-fraction on Appell-F₁ roots;
  gauge-removes spurious raw letters).
- Hand-rolling d/ds, d/dt on an on-shell family, or assembling the N×N maxcut DE →
  `de_transport` (carries the on-shell operators Oₖ=pₖ·∂pₖ: O₁=s∂s+u∂u, O₂=s∂s+t∂t,
  O₃=t∂t+u∂u).
- Any family with off-shell external leg(s) needing d/dm² DE rows at fixed Mandelstams →
  `deriv_m2` (Gauss-relation solve for D=Σcᵢⱼpᵢ·∂pⱼ ≡ d/dm²; emits ∂ₘ₂G(ν) as a finite
  shifted-G sum, harvestable from the EXISTING reduction — no re-reduction).
- Weight ≥5 Landau symbol cascade OOMs (dense kron ~2.3 TB) → `sparse_cascade`
  (sparse-support word matrix + GF(p) incremental nullspace; w5 = 173 s / 4.1 GB; w1–4
  reproduce dense exactly).
- Fitting/sampling ANY basis → UT/LS sampler discipline: rotate raw samples to the
  LS-defined canonical UT basis (maxcut LS + IBP rotation); NEVER fit in the raw Kira
  basis.

NOT-FOR: `de_transport` — steps 1–2 exercised on production-shaped inputs; the steps
3–6 driver (boundary → per-ε transport → ε-FFT → gate) runs once integrability passes
and its machinery is gated by the built-in `--selftest` toy battery, but its default
production inputs are not shipped in this repo, so a production run exits LOUD without
your own paths — known, disclosed. `decisive` Track A can OVER-report genus, never
under-reports: any "higher" verdict must be certified via Track B before acting. Not a
full-integral classifier (maxcut geometry only).

INVOKE: from `tools/`: `python3 -m maxcut.decisive SPEC.json` · `python3 -m
maxcut.ut_rotation` · `python3 -m maxcut.de_transport` · `python3 -m maxcut.deriv_m2
--audit` — the `-m` forms are the only supported direct CLI; flat shims forward script
mode via runpy (bad spec is loud rc=1). Library: `import maxcut` / `from maxcut import
...`; `from maxcut.decisive import classify`. `deriv_m2`: feed `Family(props, sp_rules,
invariants, mass_invariant, mass_leg)` (PBB1m = reference instantiation + executable
audit, NOT a hard-coded assumption; its family-definition JSON is not shipped in this repo — point DERIV_M2_FAMDEF at your own). `sparse_cascade`: `SparseWord`,
`incremental_nullspace`, `matmul_mod`.

INPUTS: decisive spec dict: `{'baikov': {'B': expr|None, 'gram': {loops,ext,ext_gram,
props,prefer_free}, 'isp': [...], 'kin': {sym:val}}}` or `{'series': list|callable,
'series_aux': ...}`; deriv_m2 `Family(...)`; sparse_cascade int64 GF(p) matrices, word
index = big-endian base-A encoding.

OUTPUTS: decisive → `{type, pf_order, genus, j_distinct, reducible, disc_chain,
caveats, wall_s}` (<30 s/case); the validate battery
`python3 -m maxcut.decisive_validate` exits nonzero if ANY executed leg FAILs
(pent5 SKIP without data is the documented state, not a failure) and writes its receipt to
`maxcut/decisive_validate.json` (generated locally, not vendored; each run
REWRITES it); deriv_m2 → m²-DE rows as shifted-G
sums; sparse_cascade → sparse word bases.

ENV: sparse_cascade GF(p) BLAS matmul: chunk auto-clamped so every partial sum
stays ≤ 2^53 (exact in float64); admits p ≤ 94 906 266, ValueError beyond;
boundary probe: `python3 -m maxcut.sparse_cascade` (also a gated validate-battery leg).
`PLD_ENV` = julia project for the Track B-jl `pf_rank.jl` backend (stays flat at
`tools/pf_rank.jl` — NOT inside the package). `MAXCUT_R3C` = dir of inputs for
the pent5 validate leg (not shipped in this repo; the leg SKIPs loudly without it).
`MAXCUT_DE_*` env vars = de_transport inputs.

GATES: decisive validate battery (topbox→elliptic, g0-3loop→polylog, ell-3loop→elliptic
j-nonconst, K3-ord3→K3 order 3, sunpair-red→2×elliptic reducible; pent5 = honest MISS,
external-data leg) — calibrate on topbox + ell-3loop→genus 1 before trusting a new family.
Public-runnable legs: topbox, g0-3loop, ell-3loop, K3-ord3, sunpair-red (inline
specs/series); pent5 needs data not shipped in this repo. deriv_m2 self-checks are EXECUTED
code, not comments: assert D[pᵢ²]=δ_{i,mass} for EVERY leg incl. the conserved one +
all-slots chain-rule walk incl. ISP numerators; gated 24–26d vs independent
parametric-derivative oracle at 3 ε. ut_rotation 5pt w1 gated 43.2d.
de_transport battery: `python3 -m maxcut.de_transport --selftest` — 4 legs on a
commuting 24×24 toy pair with a closed-form solution (route-(a) boundary series eval;
toy integrability with nonzero curl; steps 3–6 end-to-end gated ≥10 matched digits on
Laurent orders 0–3 vs closed-form Taylor references; perturbed-basis negative control
that must FAIL the gate); exit 0/1, no external inputs, seconds-scale. Production
exit codes: 0 gate PASS · 2 integrability FAIL · 3 gate FAIL.

FOOTGUNS:
- PERMANENT name collision: the package name `maxcut` can collide with any other
  module named `maxcut.py` on sys.path. Loud in BOTH shadowing directions — path
  discipline: never put a dir containing a flat `maxcut.py` first on sys.path when you
  want the package (fails loudly "'maxcut' is not a package"), and vice versa.
- Shim attribute-writes are swallowed: monkeypatching through a flat shim does NOT
  reach the package module — patch `sys.modules[<package module>.__name__]` directly.
- On-shell derivative bug class: use Oₖ=pₖ·∂pₖ; the (1/s)p₁·∂p₁ "fix" violates p₄²=0
  and contaminates A by a u-block that cancels in s·A_s−t·A_t (a null check, so it
  hides). Self-check: assert O[pᵢ²]=0 for EVERY leg including the conserved
  p_n²=(Σpᵢ)².
- Off-shell-leg d/dm²: prop-shift-only recipes give WRONG diagonals (smoking gun
  (d−3)/m² where truth is (d−4)/(2m²)); the loop·∂p chain-rule terms are mandatory.
  Walk EVERY prop/ISP definition for p-dependence, not just "the obvious one" — masters
  with ISP-numerator rows are exactly the ones that go wrong while p-independent rows
  stay correct (that split is the diagnostic). Landau-set and ϖ₀ checks do NOT catch it.
- Row-0-only-T diagnostic: ε-factorized ≠ canonical — a D-rescale + row-0 VoP can leave
  transcendental T₀'s in an off-diag block ⇒ m_UT[0,ε⁻²] mixed-weight. Stationary-point
  test: d(m_UT[0,ε⁻²])|_basepoint = 0 in all directions ⇒ NOT weight-1. Im=const@70d on
  a Euclidean path is kernel-reality, NOT a UT certificate.
- ut_rotation root ordering is chamber-dependent — keep fit and gate points
  chamber-consistent.
- `decisive_validate` battery runs at import (~40 s of module-level code) — don't
  import it in hot paths; running it rewrites the canonical receipt JSON.
- disc-chain never differentiates the carrier ISP — a saddle method that does
  produces wrong leading singularities; do not build one.

RELATED: `tools/pf_rank.jl` (Track B-jl mod-p Griffiths–Dwork backend, flat path,
ships in this repo); `tools/counterweight/canonical_form/` (Counterweight bridge); `tools/ansatzer/` (residual-count sizing before farming).

CREDIT: the Baikov representation is P. A. Baikov's [Baikov], in the loop-by-loop form
of Frellesvig and Papadopoulos [FP]; maximal cuts in Baikov variables follow Harley,
Moriello and Schabinger [HMS] and Bosma, Sogaard and Zhang [BSZ], and their use to read
off the geometry and the homogeneous differential equations follows Primo and Tancredi
[PT1, PT2] (genus drops: Marzucca, McLeod, Page, Pögel, Weinzierl [MMPPW]). Picard–Fuchs
operators by Griffiths–Dwork reduction (Griffiths [Griffiths]; Lairez [Lairez]);
canonical/UT rotations after Arkani-Hamed, Bourjaily, Cachazo and Trnka [ABCT], Henn
[Henn] and Henn, Mistlberger, Smirnov and Wasser [DlogBasis]; master counting after Lee
and Pomeransky [LP]. Bracketed keys resolve in REFERENCES.md at the repository root.

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).

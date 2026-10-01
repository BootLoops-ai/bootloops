# lockpick — guide

Tool page: https://bootloops.ai/tools/lockpick.html

KIND: package

PURPOSE: the house integer-relation / lattice value-fit family: multi-point LLL/BKZ
analytic-regression fitter (mplll) + the shared PSLQ closure harness (pslq_gate) + the
ring-basis member (lockpick.ring15 — conductor-15 / Q(sqrt(-15)) Chowla–Selberg CM
constant ring, the named-constant basket source). Closes sampled masters/boundaries
against named function pools; decides value-fit vs DE-transport.

USE-WHEN:
- Any PSLQ closure scan → `pslq_gate` (the consolidated protocol; never raw `mp.pslq`).
- Value-fit floors but you suspect HEIGHT not span → `mplll --method qrbkz` (primary) +
  `--method rawlll` (baseline). A relation gating ≥30d held-out ⇒ the LS floor was
  coefficient-height; neither lifts above LS LOO ⇒ span-deficient by independent method.
- Span-deficiency suspected on a transcendentality-graded basis →
  `mplll_graded.graded_fit(..., weights=...)` → per-weight `{K_w, locked?, floor_d}`
  localises WHICH weight floors.
- Which sample points to farm BEFORE target eval → `mplll_points` (D-optimal greedy
  pivoted-QR k<K + leverage k≥K).
- Batch-closing many ε-orders/masters against one pool (auto-batch on
  `{"targets":{...}}`); snapping a real LS solution to ℤ (`cvp`).
- Need K3-banana / CM constants in the conductor-15 field for a PSLQ basket →
  `lockpick.ring15` (27-member CS ring at dps 210+120: Ω₁₅ with the pinned h=2
  normalization, Gross unit cm₁/cm₂=φ^(−1/3), exact Hecke a_n + L(f,s) s=1,2,3 for both
  weight-3 level-15 CM newforms; disc −4 lemniscatic controls built in).

NOT-FOR: does not replace held-out/LOO certification (tools/gatekeeper).
`pslq_gate.odd_zeta_pool` caveat: dim(svMZV_w) = #compositions into odd≥3, NOT
#partitions — odd-zeta products alone don't span w≥11.

INVOKE:
- CLI: `python3 tools/mplll_cli.py` (or `python3 -m lockpick.mplll_cli` from `tools/`) `--basis B --target T --dps D --height H
  --method {qrbkz|rawlll|rawbkz|graded|cvp} [--backend ...] [--weights FILE.json]
  [--held-out ...] [--fit-points ...] [--dmin-sweep] [--no-controls]`;
  `python -m lockpick.mplll` also works; flat `python3 tools/mplll.py` via shim.
- `python3 tools/pslq_gate.py --selftest` (self-contained battery: refinds planted
  synthetic relations it builds itself in the AMFlow out-json format — positive
  controls with hand-provable answers, exercising the loader and the gamma
  normalization — + member-integrity, capacity/digit-guard, degenerate-pool, and
  headroom legs).
- ring15 member: `from lockpick.ring15 import ring15_core` — key callables
  `dirichlet_L(s,D)` (D∈{−15,−3,−4,5}), `Omega(D)`, `cm_period_eta(a,b,c)`,
  `disc4_controls()`, `ring15_hecke.a_n_list(N,s)`/`run_validation()`,
  `ring15_lfun.L_values()`. Orchestrator rebuild: `python3 -m lockpick.ring15.ring15`
  (writes RING15.json under env `RING15_OUT`, default cwd — launch from a scratch dir).
- API (lockpick.pslq_gate): `canonicalize(vec)`, `two_prec_stable(...)`, `capacity(...)`,
  `members_audit(...)`, `controls(...)`, `reverify(...)`, `graded_basket(...)`,
  `odd_zeta_pool(...)`, `load_amflow_targets(fn)`.

INPUTS: high-precision target strings/mpf + basis-member values; complex input
auto-split re/im; pool json `{"members":{...},"weights":{...}}`; graded needs explicit
`--weights` in production.

OUTPUTS: integer relation vectors; `span_residual_d` (>d ⇔ in span, else = digit floor);
`basis_relations` flag in output ⇒ prune the pool first. pslq_gate CLI exits: 2 =
control failure (refuses to report nulls), 3 = DEGENERATE_POOL.

ENV: `MPLLL_FPLLL_TIMEOUT` — per-reduction wall for the subprocess backends (fplll-cli,
Nemo.lll); unset or 0 = unbounded (the default: a reduction ends by convergence or by its
own failure, never by this clock); a positive integer of seconds is an opt-in the caller
states in its receipt; anything else refused by name (ValueError). A reduction the knob
ends raises `MplllTimeoutExpired` (a `subprocess.TimeoutExpired`) whose text reads
`ended by MPLLL_FPLLL_TIMEOUT=<N> s, not by convergence` (also on stderr); an ending on
the nemo backend propagates (it is not absorbed by the pure fallback); the value is
ASCII digits. Regression:
tests/test_mplll_timeout.py. Backend chain: fpylll → fplll-cli
→ Nemo.lll → pure fallback. All lattice-build paths own their workdps(d+35) — API safe
at any ambient dps (regression: tests/test_mplll_ambient_dps.py, shipped).

GATES (the PSLQ discipline — mandatory on every scan):
- Sane height: `capacity()` refuses height>capacity scans (the "null = height failure"
  class). Capacity + insample < dps ⇒ self-diagnosed false positive, correctly rejected
  — trust the gate.
- Two-precision stability with canonicalize-BEFORE-comparing: identical sign+gcd-
  canonical vectors at both dps (raw-vector comparison gives false instability).
- Positive controls planted INSIDE the scan basket, coefficients ≥10× BELOW maxcoeff
  (controls AT the cap fail spuriously). Control failure ⇒ no nulls reported, exit 2.
- `members_audit` before every scan (CLI main() runs it): any NULL on a degenerate
  basket is VOID.
- Held-out: ≥30 genuine digits + leave-one-out vs an independent oracle never used in
  the fit before anything is called closed; `reverify` at 210d.
- Sector-restriction capacity trick: when the measured relation support lives in one
  sector, restrict the pool to that sector BEFORE escalating dps — full pools starve
  capacity at high weight. State the sector restriction as a caveat in the rank claim.

FOOTGUNS:
- mp.pslq relation with coefficient 0 on the target = internal relation among basket
  members, NEVER a null → `DegeneratePoolError`; quotient the pool, rescan.
- qrbkz needs K+1 ≤ N_fit (QR asserts m≥n) — use rawlll/rawbkz when N_points < K+1;
  qrbkz_batch drop-threshold carries a 5-digit margin vs single-target.
- Plain LLL loses a planted 6-term control at N·dps≈810 — need BKZ-20+ or N·dps≳1000.
  Joint 2-point lattice reaches heights ~10^(joint_digits/(K+1)).
- Dilog/trilog basket design (args 1/4, 1/3, −1/3, 6th roots of unity) is a degeneracy
  minefield (Landen/Abel five-term, Cl2/ImLi2 duplicates, Glaisher Sl_n(π/3)
  reductions) — curate to a ℚ-basis with `members_audit`; keep BOTH engines (pslq_gate
  + mplll) in every scan.
- Two-prec leg runs at d−20: the data must support d.
- Relation-height growth law (measured on sv relations w=11→21: ~×100 per weight step):
  size maxcoeff ladders by this before filing NULL; a null at the previous weight's height ceiling is a HEIGHT failure at the new weight.
- Two-lindep-legs spuriosity trap: when ≥2 relations coexist in the pool, the two legs
  return DIFFERENT genuine lattice vectors — a leg-agreement criterion then spuriously
  accepts/rejects. Fix: a full-precision residual arbiter on the value bank, never
  vector-identity across legs alone.
- Near-1 members burn the maxcoeff budget: ζ_w ≈ 1+2^−w makes (1, ζ_w) nearly
  degenerate, so 'one' in a tw>0 pool makes mp.pslq spend its coefficient budget on
  that direction and return false NULLs. Exclude near-unity members; rational content
  surfaces as an escalation class instead.
- Farm workloads: greedy independence filters are BASKET-level facts — persist them to
  disk keyed (tw,dps,tier); process-local caches make every farm worker repay the full
  greedy PSLQ cascade.
- Cross-route check: before ANY PSLQ/mplll scan over ε-layer constants of an
  Euler/Feynman PARAMETRIC integral — and again the moment such a scan hits a wall —
  check the DERIVATION route (symbolic integration of the parametric representation
  closes that class exactly and demotes PSLQ to the held-out gate). Fitting harder is
  the wrong first response to a wall in this class.

ring15/ is a sha-pinned snapshot of the member code (PINS.json is the pin table;
the selftest asserts the shas). Dependency note: dirichlet_L / cusp_dictionary are
Eichler.jl's, linked never copied (upgrades/Eichler.jl); ring15_core.dirichlet_L is
this member's own mpmath implementation.

CREDIT: lattice-reduction analytic regression from high-precision samples is the method
of Oscar Barrera, Aurélien Dersy, Rabia Husain, Matthew D. Schwartz and Xiaoyuan Zhang
[BDHSZ]; `rawlll` is their construction and `qrbkz` our conditioned variant of it. PSLQ
is Ferguson and Bailey's algorithm as analysed by Ferguson, Bailey and Arno [PSLQ], used
through mpmath (F. Johansson and contributors); LLL is Lenstra–Lenstra–Lovász [LLL], BKZ
through fplll (the FPLLL development team) and Nemo [Nemo]. ring15's newform labels
follow the LMFDB (The LMFDB Collaboration, lmfdb.org) with Hecke data cross-checked in
PARI/GP (the PARI Group). Bracketed keys resolve in REFERENCES.md at the repository
root.

## MEMBER — sealrun (preregistered sealed-closure runner; pslq_gate's seal-wrapper)

PATH: `tools/lockpick/sealrun/` — `t6_close.py` (the runner), `t6_basket.py` (the
basket builder), `BASKET_T6.json` (the frozen 24-member basket), plus the two
committed selftest targets.

KIND: verdict INSTRUMENT built on the pslq_gate scan ENGINE. The engine-vs-instrument
distinction routes a user: pslq_gate = the library primitives (two_prec_stable,
capacity, members_audit, reverify, controls) for exploratory fishing; sealrun = the
SEAL you graduate to once targets are defined and a frozen basket can be declared.

PURPOSE: decide CLOSED-with-name vs NULL-with-exclusion-floor for a high-precision
constant against a FROZEN member basket — the seal (basket + protocol pinned BEFORE
the target exists), fresh members_audit + planted/negative controls INSIDE every run,
two-precision canonicalized legs, height ladder 1e4 → 1e12 with capacity-refused
rungs, fit legs at most (target_digits − 40) digits, ≥30 held-out digits, HIT
reverify at (target_digits − 5) dps against digits no fit leg ever saw; every
verdict field script-emitted to `CLOSE_<label>_<stamp>.json` (written under `--out`,
default cwd).

INVOKE: `python3 t6_close.py --target-string-file FILE [--label NAME] [--out DIR]`
or `python3 t6_close.py --target-json FILE:key.path`. Complex targets are split
re/im, one run each; the runner refuses targets under 150 digits (≥300 typical).
EXIT CODES (the contract drivers script against): 0 HIT-CLOSED, 1 NULL-with-floor,
2 controls failed (NO verdict of any kind), 3 degenerate pool (NULL VOID),
4 input/digits error. Basket REBUILD: `python3 t6_basket.py` — double-routed
members, members_audit, planted/negative gates asserted green INSIDE the emitted
BASKET json; rebuild inputs (sha-pinned reference receipts) are read from `$SEALRUN_R1`
(no default shipped — the member battery consumes the frozen BASKET_T6.json copy and
never rebuilds).

INSTANTIATION LAW: the registered member is the PATTERN (runner + basket-builder
pair + prereg discipline), instantiated per project — t6_close.py binds one basket
filename and one label default; a new project re-instantiates with its own frozen
basket by copying the runner. Do not ship a "generic" binary; ship the sealed pair.

THE DegeneratePoolError DOCTRINE: mp.pslq returns the LOWEST-height relation in the
lattice — a genuine internal relation among basket members shadows the target at
EVERY rung, so a NULL from a degenerate pool is VOID, not weak. NULL is a verdict
ONLY with a same-run clean members_audit at scan dps (exit 3 ⇒ no exclusion floor
quotable). QUOTIENT THE POOL BEFORE FREEZING. Verdicts come from CLOSE_*.json,
never from driver logs. A NULL excludes ONLY the frozen basket at the printed
height — state what it does NOT exclude.

MEMBER BATTERY (deliberately NOT a leg of the registered package battery: the pair
measures minutes, not seconds — over the seconds-scale battery standard and the
selftest runner's per-attempt timeout — run it from a scratch cwd when touching
this member):
- `python3 <pkg>/sealrun/t6_close.py --target-string-file <pkg>/sealrun/selftest_planted_target.txt --label planted`
  must exit 0 CLOSED (reverify relresid ~1e-296, 35 held-out digits);
- the same with `selftest_null_target.txt --label null` must exit 1 NULL-with-floor
  at 1e8 (the 1e10 rung capacity-refused at 300d/24 members).

BASKET-SHA NOTE: BASKET_T6.json pins its two rebuild inputs by basename + sha256;
member values, weights, digit counts and gate records are part of the frozen file,
so any receipt you emit should pin this file's own sha256.

MEASURED COSTS: basket build 59 s one-time; per-target closure ~60 s (low-rung hit)
to ~210 s (full null ladder) single-core at 300d; ~240 s at 504d/23 members.
Measured honest capacity ceilings: 300d/24 members → 1e8; 341–348d → 1e10;
504d/23–30 members → 1e13; every +40d buys ~+1e1.6 of height.

## MEMBER — bilmine (bilinear-lattice miner: integer stage + sqrt15 ring window)

PATH: `tools/lockpick/bilmine/` — `bilmine_lib.py` + drivers (`seal_prereg.py`,
`run_control.py`, `run_exact.py`, `run_mine.py`, `remine_points.py`,
`run_global.py`, `emit_verdict.py`, `close_leg.py`), synthetic `selftest.py`;
ring stage under `f1/` (`bilmine_f1_lib.py` + the `*_f1` drivers +
`resume_mine_f1.py`).

KIND: verdict INSTRUMENT (one member, two stages) built on this package's
primitives — `mplll_lattice.reduce_rows` + `pslq_gate.canonicalize` — under the
sealrun-class prereg discipline.

WHAT: a prereg-sealed miner for EXACT bilinear relation lattices
C^T S_loc C = S_mum on certified period frames, by homogeneous simultaneous
integer-relation mining: [I_n | round(10^d · row-normalized coeffs)] embedding,
LLL + canonicalize, under the full house discipline — capacity law
((n+1)·log10(H) + 20 ≤ d_leg; over-capacity rungs refused and printed),
two-precision legs adjudicated at the LATTICE level (HNF of the span — LLL
representatives legitimately differ between legs), planted + negative controls
sealed before any production number, TRUE HOLDOUT (held-out S_mum entries are never
free unknowns; determined values gated near-integer post hoc), residual arbiter at
d+60 with row-scale-relative floors, ladder-with-floors, script-emitted verdicts
only. Exact stage: python-flint `fmpz_mat` nullspace/HNF (closes in seconds systems
sympy stalls on).

F1 RING-WINDOW MODE (Z → Z[sqrt15]): the integer miner is blind to bilinear
lattices whose two sides carry an irrational relative scale; F1 generalizes the
coefficient ring with three sealed window modes per (point, sector) — GEN (doubled
basket: every coefficient a + b·sqrt15, 2n unknowns) and the one-sided O-A / O-B
(exact column SLICES of the GEN rows, asserted by slice_check; higher
capacity-lawful rungs at half the unknowns). Ring holdout gate: determined held-out
values must land on the mode's allowed ring shape (GEN via a 3-unknown LLL
mini-reduction on [v, 1, sqrt15]; |c0| = 1 required). Sealed ABSOLUTE-floor ring
detector: Z + sqrt15·Z is dense in R, so relative acceptance criteria are BANNED.
Exact ring invariance is checked COMPONENTWISE on the integer and sqrt15 parts.
CONTROL-UPGRADE LAW: an integer plant cannot license a ring instrument — fresh
sqrt15-planted controls are required (the C1'/C3' pattern).

DRIVERS are the pattern, re-instantiated per project; every driver refuses loudly
when its data contract is unset. Chain: seal_prereg → run_control (C1/C2/C3) →
run_exact → run_mine → run_global → emit_verdict → close_leg; `f1/` mirrors the
chain for the ring stage. Data contract (env vars, NO defaults shipped):
`BILMINE_LEG`, `BILMINE_F1_LEG`, `BILMINE_PARENT`, `BILMINE_T5`, `BILMINE_PFRAME`,
`BILMINE_C1_PI`, `BILMINE_C1_N`; `BILMINE_PRODUCER_LINT` optional (lint step
skipped when unset).

MEMBER SELFTEST (committed synthetic data only, no env vars, seconds-scale):
`python3 tools/lockpick/bilmine/selftest.py` — S1 planted bilinear quadric
recovered exactly BLIND (both precision legs, both engines, rank 1); S2 committed
sha512-derived negative NULLs at every lawful rung on both engines; S3 ring-detector
smoke 5/5 (accepts exact rational and genuine (a + b·sqrt15)/q at 1e60 scale,
refuses pi and a sha-seeded fractional, near-zero short-circuits).

REGISTER RULES (binding): a GEN hit with all-zero b-parts (or all-zero a-parts) is
inside the parent's integer NULL scope — flag INSTRUMENT-INCONSISTENCY, never
FOUND; a found lattice is REPORTED, never consumed in-leg — every consumption gates
separately; capacity-priced-out blind windows (pilot-priced) are quotable scope and
final for the leg — wider is a NEW door.

FENCES: NOT a sealed-constant scanner — sealrun (MEMBER above) stays the
CLOSED-vs-NULL door for constant verdicts; the engine-vs-instrument law routes
constants there and relation-lattice mining here. The seal here is bilmine's OWN
prereg ceremony: `seal_prereg.py` → `SEAL_RECORD.json` (prereg pinned BEFORE any
relation number; the verdict quotes the seal sha).

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).

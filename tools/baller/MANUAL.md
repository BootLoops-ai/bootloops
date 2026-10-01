# baller — manual

**BALLER** — arbitrary precision ball arithmetic for science. Ball-arithmetic
instruments ONLY; never a header for Arb.

> BALLER never aliases/shadows Arb's namespace,
> but it OWNS the user-facing ball-arithmetic layer that builds with
> Arb — the FRONT DOOR section below (`baller.run` / `baller.solve` /
> `baller.render`).

Tool page: https://www.bootloops.ai/tools/baller.html

Layout: tools/baller/ (package baller/ + vendor/ +
battery/). Assembly: engines VENDORED with sha-pinned
headers (Terrier pattern; `baller/_pins.py`, fail-closed `baller.verify()`);
registered tools ALIASED in place by identity (posq, wayfinder
gate/manifest); the hygiene wing OWNS two member modules outright
(`baller/dps_lint.py`, `baller/purity_scan.py`).
KIND: package

## Quickstart
```python
import sys; sys.path.insert(0, "<your-checkout>/tools/baller")
import baller
baller.verify()                      # fail-closed integrity, run first
from baller import quad, transport, certify, contract, hygiene
# the front door: certified evaluation without the machinery
r = baller.run(lambda xp, x: 111 - 1130/x + 3000/(x*xp), x0=(2, -4),
               n=100, dps=50)        # Muller's recurrence, ball arithmetic
res = baller.solve(lambda xp, x: 111 - 1130/x + 3000/(x*xp), 16,
                   x0=(2, -4), n=100)      # adaptive certified driver
print(baller.render(res, digits=17))       # -> 6.0000000160995649
```

## FRONT DOOR (the user-facing layer, battery L17-L18)

- **`baller.run(f, x0=None, n=None, dps=50)`** — execute a user expression
  or iteration in ball arithmetic at working precision `dps` decimal digits
  (= ceil(dps·log2 10) bits, NO guard bits: the working precision IS the
  user's ask; escalation is solve's job). Thin — every operation is
  python-flint/Arb (arb/acb), nothing re-implemented. Expression mode
  (`x0=None`): `f` is a 0-arg callable or a str eval'd with
  `{arb, acb, flint}` in scope. Iteration mode: `x0` = seed or chronological
  seed tuple, `f` = step map called `f(*window)` on the last k values OLDEST
  FIRST, iterated to index `n`; returns a `RunResult` (`.values`, `.value`,
  `.blown_at`, `.render()`).
- **`baller.solve(f, target_digits, x0=None, n=None, index=None,
  start_dps=None, max_dps=4096, factor=2)`** — adaptive certified driver:
  run under ball arithmetic watching the radius; on blowup (non-finite ball
  or radius swamping the target) ESCALATE and rerun. Escalation policy
  (documented, deterministic): start at `start_dps` (default
  max(target_digits+10, 30)), DOUBLE (`factor=2`) per failed attempt up to
  `max_dps` (default 4096); the cap raises **`SolveRefused`** — a TYPED
  refusal NAMING the achieved digits (+ attempt log + best enclosure),
  never a bare wrong number. Success returns a `SolveResult`: mid ± rad
  with ≥ target certified, `.dps`, `.attempts`, `.render()`. `f` is either
  the step map (with x0/n; `index` picks which value to certify, default
  final) or a callable `f(dps) -> ball`.
- **`baller.render(x, digits=None, strict=False)`** — FAIL-CLOSED printing:
  prints ONLY certified digits. Certification = the 1-ulp display
  criterion: a printed string p with d significant digits is certified iff
  |p − t| ≤ 1 unit of p's last printed place for EVERY t in the ball,
  verified as a certified arb comparison (round-to-nearest display; the
  reference prints x_100 = 6.0000000160995649 against true ...64889...).
  Over-requests are CLAMPED to the certified count, or REFUSED typed
  (`strict=True`, `RenderRefused` naming certified vs requested).
  Uncertified values (infinite radius, non-finite mid, zero-mid fat balls)
  come back as a MARKED string `UNCERTIFIED[~mid +/- rad]` or a typed
  refusal — NEVER printed bare. `certified_digits(x)` exposes the count.
- **DOCUMENTED DEVIATION from Arb's convention**: in the iteration driver, when a step's
  ball goes non-finite — canonically a divisor ball CONTAINING 0, where
  Arb returns an indeterminate [nan ± inf] — the radius goes to **+inf**
  and the **CENTRAL VALUE CONTINUES as the fixed-precision stream** (the
  mid-trimmed iteration at working precision, i.e. exactly what plain
  fixed-precision arithmetic computes). [anything ± inf] is a VALID
  ENCLOSURE — deliberately NOT Arb's nan/indeterminate choice: the mid
  keeps telling the fixed-precision story (the Muller trajectory into the
  float trap at 100) while the infinite radius keeps it honestly
  uncertified; render marks it, never prints it bare. One-shot expression
  mode has no stream to continue and returns Arb's ball as computed.
- **Reference behaviors** (Muller recurrence x_{n+1} = 111 − 1130/x_n +
  3000/(x_n·x_{n−1}), x0=2, x1=−4; battery L17 gates these EXACTLY vs an
  independent raw-flint truth): 50-digit run certifies x_30 mid
  6.0056486887714 with rad 1.29e-5 (render → "6.0056") and returns
  [mid→100 ± inf] UNCERTIFIED-marked for x_100; solve(target 16) escalates
  30→60→120→240 dps and certifies x_100 = 6.0000000160995649 (79 digits at
  240 dps); solve capped at 60 dps refuses typed naming 0 achieved digits.
  **NEGATIVE CONTROL** (documented, in the battery): the same recurrence in
  plain Python floats converges to EXACTLY 100.0 (true limit 6) — the bare
  wrong number this layer exists to prevent.

## Modules
- **frontdoor** — the user-facing layer: `run`/`solve`/`render` +
  `certified_digits`, exposed at top level (`baller.run(...)`); full
  contract incl. the zero-divisor deviation in the FRONT DOOR section
  above and the module docstring. Standalone-loadable by design (no
  package-relative imports) so the battery's mutation harness gates the
  file directly (L18).
- **quad** — the gated POSQ core (certified positive-quadrature evidence),
  aliased: `quad.posq`, `quad.kernel_io`. Validated (adversarial battery @3 pinned seeds, production controls 15/15,
  CTMC quartet); battery leg L6 re-runs its adaptation battery live.
- **transport** — vendored kklt engines: `pipe_transport` (majorant-tail
  certified transport), `pipe_lib` (exact recurrence/tower layer),
  `march_lib` (arb Taylor endpoint transport for den(s)·s·v' = Anum(s)·v;
  carries the ONE vendor delta: the complex-step arb(acb) fix, gated pre-FAIL
  / post-PASS with 256-bit truth).
- **certify** — vendored `pipe_vac`: generalized Krawczyk / interval-Newton
  + forward-AD-over-acb (unique-zero machinery).
  `certify.block_krawczyk` = `baller/block_arrow.py` — generic BLOCK-ARROW
  structured Krawczyk (B independent diagonal blocks + a global border, the
  hierarchical-MAP arrow shape; blockwise approximate inverse + interval
  Schur complement on the border; strict-interior contraction semantics;
  cross-block row sums bounded LINEARLY in the number of blocks; radius
  contract: scalar OR PER-COORDINATE (rz, rg) — the anisotropic box a
  stiff hierarchical Hessian requires, battery leg g) + interval
  block-arrow PD verification (`pd_block_arrow`: ONE global certified
  congruence with the float block-arrow Cholesky factor, exact
  power-of-two equilibration, exactly-nonsingular block-triangular X, no
  rigorous inverse anywhere — the per-block-Schur routes measurably fail
  on cond~1e13 collinear-latent blocks: Neumann bound rank-collapses,
  certified interval solve fattens ~1e2; `interval_cholesky`: scalar
  reference) + `certify_local_min` (one oracle evaluation, both legs).
  Refusals name the failing block/border and the measured contraction
  defect. Oracle contract (`dims`/`F`/`H`, arb/arb_mat enclosures) in the
  module docstring; violations raise typed `BlockArrowContractError`.
  No knobs; radius policy belongs to the caller (CLINCH pins its own).
  Production record: CLINCH certificates at 5592/5955 coordinates
  (582/~600 blocks + 51-border, ~90 s at 192 bits).
- **contract** — halt-not-average discipline: `gate`/`manifest` (wayfinder,
  aliased), `mc` (batched MC + Clopper-Pearson certified quantiles,
  vendored), `dual()`/`Tripwire` (minted dual-path tripwire: disagreement
  HALTS with both values, never averages; sampled 1-in-N oracle mode),
  `run_mutation_gate` (minted mutation-tested-gate harness; scratch-copies-
  only law enforced).
- **hygiene** — the `dps_lint` member module (mpmath import-time-dps AST
  linter: module-level constants frozen at dps=15, cross-module dps resets,
  complex-step `mp.diff` on real-branched callables; stdlib-only; CLI
  `python3 baller/dps_lint.py FILE_OR_DIR ...` or `python3 -m baller.hygiene
  lint ...`, exit 1 findings / 2 bad path / 0 clean; `--selftest` built-in
  battery, also leg L21) + minted enforced checks:
  `lint_fixed_wdps`, `lint_ambient_conv`, `lint_unkeyed_cache`,
  `lint_ctx_restore` (source lints; findings = nonzero), `ftrim` +
  `RadiusWatch` (radius discipline: mid-trim only after ACCEPTED steps;
  watch halts at the relative-radius bar), `bridge_mpf` (re-round-safe mpf
  bridge via make_mpf — BOTH mp.mpf(x) and mp.mpf(x._mpf_) re-round;
  measured 115-bit loss on the trap path), plus the moved-in
  `purity_scan` member module (standalone finite-precision-literal
  scanner: MPF-OF-FLOAT / DEC-STRING / FLOAT-EXPR / auto-pure
  pre-classification over an arbitrary file list; report-only; CLI
  `python3 baller/purity_scan.py <list-of-abs-paths.txt> <out.json>`).

## Vendored members (vendor/certlane, vendor/ling_onesided, vendor/lgf_monodromy)

Three sha-pinned member families ship beside the kklt/geo engines. They are
regular pinned vendor bytes (`baller.verify()` covers all 22 files) and load
as plain modules from their directories; they are NOT wired into
`_core.VENDORED` (no bare-name loader entries). Battery leg L20 keeps them
live: the certlane pytest suite must pass and every module must import from
the vendored bytes alone.

- **certlane/** — certified arb-ball enclosure primitives for the matPTF
  alert computation (matPTF: the INGV probabilistic tsunami forecasting code,
  Selva et al., "Probabilistic tsunami forecasting for early warning", Nat.
  Commun. 12, 5677 (2021); certlane re-implements the enclosed layers and
  copies no matPTF code). certlane is a research re-implementation used to
  study numerical certification of a published hazard computation; it is not
  a tsunami-warning system, is not endorsed by INGV, and must not be used
  operationally. Layers covered: (lognormal survival, weighted survival sums, normalized
  sums, percentile brackets with convention-ambiguity unions, alert-level
  ladder; fail-closed INDET semantics), plus `test_primitives.py` (26
  containment/planted-error/convention pytest cases vs an mpmath dps=60
  oracle) and `attack_primitives.py`, a four-surface adversarial harness
  (containment fuzz, planted-error corruption subprocesses via
  `sys.executable`, semantics vs the float64 reference, precision-law
  probes). The semantics surface needs the unvendored float64 reference
  implementation on `sys.path` and recorded curve archives via
  `ATTACK_EVENT_NPZ` (os.pathsep-separated .npz paths) — it prints a named
  SKIP without them.
- **ling_onesided/** — positive-integrand one-sided evidence kit over
  quartet pattern counts: `taylor_p` (centered-form Taylor/dual enclosures,
  planted-error selftest), `gate_price` (priced-lever gate layer),
  `tail_engine` (certified sup-bound tail integrals), `balanced_unc` /
  `double_unc` / `collapse_unc` (2-dim balanced, two-sided and collapse
  enclosure supports), `sd_engine` (certified evidence engine — peeling
  circuit, 45-digit anchor gate + planted mutations), `sd_onesided`
  (one-sided fixed-needle-grid route: interior-cell inf sums = rigorous
  LOWER, hull + one-sided complements = UPPER), `l5_center` (first-order
  centered cell enclosures). The data fixture `QUARTET_COLLAPSE.json` is
  vendored beside them (read root override: `WPG_BASE`); its pattern counts
  are derived from the DravLex v1.0 Dravidian lexical database (Kolipakam, Jordan,
  Dunn, Greenhill, Bouckaert, Gray & Verkerk 2018), which is distributed under CC BY 4.0 — see THIRD_PARTY.md.
- **lgf_monodromy/** — `valmono.py`: validated ball-Taylor Fuchsian
  monodromy transport of an order-5 exact-integer-coefficient ODE (disc
  radius from a series lower bound on |b5|, Gronwall/Cauchy truncation
  tails as ball radii), with the pinned coefficient cache `BC_cache.pkl`
  (sha in `PROVENANCE_BC_cache.txt` and `_pins.py`). KNOWN GAP (flagged,
  not hidden): the loop GEOMETRY module (`loop_paths`) is not vendored —
  `loop_waypoints()` raises a typed ImportError; pass explicit waypoint
  lists to `transport_loop(..., waypoints=...)`. `__main__` /
  `valmono.selftest()` runs a synthetic null-homotopic-loop check
  (monodromy must contain the identity; ~3 s at prec 160).

OUTPUT LAW for the members: every default output root is the CURRENT
working directory (`WPG_RESULTS`, `SD_ENGINE_OUT`, `SD_ONESIDED_OUT`,
`COLLAPSE_OUT`, `VALMONO_OUT` override), NEVER the vendor tree — any
unpinned file under `vendor/` is a tamper class and bricks
`baller.verify()` (the L12 unpinned-extras guard is the enforcement).

## GATES (battery/battery.py — run from a scratch cwd; refuses the tool
## tree, tools/, and the protected trees on its deny list)
21 legs, OVERALL PASS, ~2 min (the posq adaptation battery through the
alias dominates) — L17:
the reference behaviors as EXACT gates vs an independent
raw-flint truth incl. the [mid-stream ± inf] deviation ball, the typed
SolveRefused, fail-closed render units, and the plain-float negative
control returning bare 100.0; L18: 4 mutants via frontdoor_case.py — a
render that prints uncertified digits FAILS, a solve that returns without
the target FAILS, a dropped deviation mid-stream FAILS, a dead escalation
FAILS — 4/4 caught; L19: purity-scan planted 3-class smoke — a planted
sample file must classify MPF-OF-FLOAT, DEC-STRING, and PURE-CONV
correctly through the member CLI and the `baller.purity_scan` import; L20:
the vendored fold members are live code — certlane's 26-case pytest suite
passes from a scratch cwd, all 9 ling_onesided modules import from the
vendored bytes (fixture found beside them, N=318), valmono loads the
pinned BC cache and its null-loop identity self-check passes; L21: the
dps_lint member's built-in selftest (one planted finding per rule x 5
rules, exit codes 1/0/2/2, fixed-fixture control; 11 pass / 0 fail) runs
through the member CLI with TMPDIR pinned to the scratch cwd, `python3 -m
baller.hygiene lint` exits 1 on a planted module-level constant and 0 on
the fixed copy, and `baller.hygiene.dps_lint is baller.dps_lint`. Base legs: pins + tamper/pyc/foreign-shadow
mutation controls · alias identity · march truth case mutation-tested 3/3
(incl. the tail-load-bearing K=6 case) · posq adaptation battery via the
alias (32 s, pinned seed) · wayfinder tests + exact Clopper-Pearson
cross-check (1e-12) · tripwire halt semantics · mutation-harness selftest ·
hygiene lints + runtime helpers · kklt engine loads + exact-helper sanity ·
L12-L14 lock the integrity/minted-check/tripwire guards · L15 acb_calc
certified integration · L16 block-arrow Krawczyk/PD (planted zero certified
+contained, displaced/fat-box refusals named, sign-flip mutation controls
2/2, precision-starvation refusal, planted saddle + scalar interval-Cholesky
cross-check, typed contract violations).

## FOOTGUNS
- Engine bytes: vendored set fail-closed (verify + per-load source hash;
  exec(compile), bytecode cache never read, purged per load; UNPINNED files
  under vendor/ = tamper class; loads thread-locked, sys.modules restored on
  every cache hit).
- DOCUMENTED SEMANTICS: an engine already loaded in-process keeps being
  served from the in-memory cache even if a LATER on-disk verify() fails —
  its bytes were verified at its own load time.
- ALIAS THREAT BOUNDARY: pre-existing sys.modules entries are authenticated
  by import-spec provenance + on-disk existence — this catches every
  accidental shadow and naive plant; an in-process actor who can forge a
  full ModuleSpec can equally monkeypatch baller itself, so no in-process
  check can beat that class.
- march step() ENFORCES its contract (|h| provably < r, r provably > 0,
  len(v0)==RK, den-nonzero survives -O) — out-of-contract inputs would
  silently ABSOLUTIZE a negative tail into a too-small certified radius.
- pipe_transport cxq/fr helpers are EXACT-INPUT instruments: Fraction
  components ONLY — float components silently produce non-bounds (the
  battery uses exact inputs).
- kklt engines set flint ctx.prec at function tops WITHOUT restore (their
  pilot design; the lint flags them honestly) — wrap engine calls in
  hygiene.ctx_guard when your process shares flint state.
- tripwire compares in mpmath at bar+20 digits (never float64); complex AND
  ball inputs are rejected typed (a ball's radius cannot ride through a
  scalar comparison) — compare mids/parts per component, watch radii with
  RadiusWatch.
- IN-PROCESS THREAT BOUNDARY (applies to vendored loads too): dep
  identity is re-verified after every engine exec (import-hook swaps refuse
  typed), but an in-process actor who controls sys.meta_path/sys.modules can
  equally monkeypatch baller itself — no in-process check beats that class.
- BARE-NAME RESIDUAL: only pipe_lib/pipe_transport get bare sys.modules
  entries (their engines import each other bare); displacing a foreign
  occupant prints an ADVISORY. mc/march_lib/pipe_vac are dotted-only — your
  own ./mc.py is safe.
- PL/PT CONFIGURE LAW: pipe_transport snapshots part of pipe_lib's card and
  reads the rest live — always reconfigure through
  baller.transport.configure_all(...), never PL.configure alone.
- mc contract: batch<=0 and out-of-domain quantile_ci REFUSE (the
  alternatives are an infinite hang / silent (nan, nan)).
- ALIASES resolve the tools/ tree beside this package (posq, wayfinder
  ship there; `BALLER_TOOLS_ROOT` overrides the root). dps_lint is not an
  alias: it lives in this package.
- posq guard-inserts its own dir and removes it after its imports;
  baller's own inserts are all removed after import.
- ARB POW FOOTGUN (python-flint 0.8.0): `x**n` on a
  ZERO-CONTAINING arb ball returns NaN (pow routes through exp·log) —
  write `x*x`, never `x**2`, in any oracle evaluated over a box. The
  battery oracle carries the regression comment.
- KRAWCZYK TRUST BOUNDARY (block_krawczyk, documented not curable): the
  centered contraction test is SIGN-BLIND in c (|-c| = |c|) and any
  invertible claimed-H self-preconditions to E ~ 0 near a true zero — a
  post-hoc gradient-block sign flip has the identical zero set (certificate
  stays true), and an H-only lie is a contract violation Krawczyk cannot
  see from one center. Meaningful mutation controls flip signs in the
  ORACLE MATH (zero moves -> Krawczyk refuses) or hit the PD leg (border
  flip -> pd_block_arrow refuses). Oracle-vs-model truth is gated
  ORACLE-SIDE (CLINCH: identity gates vs the receipted reference + FD
  Hessian gates + its own mutation legs).
- bridge trap: never mp.mpf() an existing mpf across a dps change (both
  direct and tuple forms re-round); use hygiene.bridge_mpf.
- The lints are advisory instruments with measured fixtures — they flag the
  documented trap classes; a clean lint is not a rigor proof.
- pipe_vac/pipe_transport/pipe_lib are DEEP engines: production-scale
  validation is documented on the tool page; the battery exercises load
  integrity + exact helpers, not the production runs.

RELATED: tools/posq (see its GUIDE.md); tools/wayfinder.
CREDIT: BALLER is a front door over Arb, the ball-arithmetic library of Fredrik
Johansson (2017, IEEE Trans. Comput. 66:1281), now part of FLINT (W. Hart, F. Johansson,
A. Ahlbäck and the FLINT developers; flintlib.org), reached through python-flint (F.
Johansson, O. Benjamin and contributors); every certified digit in this package is
theirs first, and we are grateful to them. Certification: Krawczyk (1969, Computing
4:187), Moore (1977, SIAM J. Numer. Anal. 14:611), Rump (2010, Acta Numerica 19:287).
The front-door reference figure is J.-M. Muller's recurrence (Muller et al., Handbook of
Floating-Point Arithmetic, Birkhäuser 2010). certlane re-derives, in ball arithmetic,
layers of INGV's matPTF (Selva, Lorito, Volpe, Romano, Tonini, Perfetti, Bernardi et al.
2021, Nat. Commun. 12:5677; code github.com/INGV/matPTF) — no matPTF code is copied and
INGV does not endorse this re-implementation; ling_onesided's counts derive from DravLex
v1.0 (Kolipakam, Jordan, Dunn, Greenhill, Bouckaert, Gray & Verkerk 2018, R. Soc. Open
Sci. 5:171504; CC BY 4.0). The vendored kklt/geo engines are in-house and ship
sha-pinned; the minted checks encode the documented footgun catalog. POSQ/wayfinder
remain their own registered tools, fronted by identity.
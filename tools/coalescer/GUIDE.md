# Coalescer — spectral-projector extraction of threshold connection coefficients

Tool page: https://www.bootloops.ai/tools/coalescer.html

KIND: package (coalescer.py CLI with two doors: `--op FILE.json` = your operator +
seed (template: examples/gauss_2f1.json), `--masses` = the banana worked example;
routeA_monoproj.py (generic kernel
run_operator + load_operator_json + the banana wrapper run), routeB_monoproj.jl,
pf_engine.py, dump_pf_data.py, pslq_calpha.py, gate_compare.py, calpha_rings.py,
GATE.json; worked-example data pf_k3_data.jl, pf_cy3_data.jl, pf_cy4_data.jl,
tests/fixtures/*).

PURPOSE: Connection coefficients of a Fuchsian operator onto its fractional-power
(finite-monodromy, non-unipotent) Frobenius branches at a singular point: given
L = Σ_j P_j(z) θ_z^j with exact polynomial coefficients and the exact Taylor
coefficients of a holomorphic solution f at z = 0, it returns c_α = the coefficient
of f's continuation on each branch Φ_α = t^α(1 + …) at z = ∞ (t = 1/z), for every
fractional exponent α the indicial equation there has.  The extractor is the
spectral projector of the local monodromy M₀,
P_λ = (M₀−I)ⁿ·∏_{μ≠λ}(M₀−μI) / [(λ−1)ⁿ·∏_{μ≠λ}(λ−μ)], which kills the rank-n
unipotent log-tower and all other semisimple eigenspaces exactly; order, exponent
set and eigenvalues are read from the operator.  Typical clients: the threshold or
pseudo-threshold of any Feynman-integral family whose Picard–Fuchs operator is
known (square-root and quarter-power branches), conifold-type points of Calabi–Yau
operators, and P-recursive sequences whose generating function has an algebraic
branch point (the connection coefficient fixes the asymptotic constant).  The
L-loop banana family at its threshold is the WORKED EXAMPLE shipped with the
package (operators, seeds, gates and PSLQ closures for K3 (1,1,1,9), CY₃
(1,1,1,1,16), CY₄ (1,1,1,1,1,25)).

USE-WHEN:

- Extract c_α (a chosen holomorphic solution onto a fractional-exponent branch at
  another singular point) for ANY Fuchsian operator you can write down with exact
  coefficients, at any order / eigenvalue set — bring it through `--op FILE.json`.
- The banana coalescence ladder and relatives (`--masses`, operator built for you).
- Need an independent two-method gate on such a coefficient (Route A mpmath/
  s-coord/CCW small circle vs Route B Arb/z-coord/CW large circle: different
  arithmetic, stepper, coordinate, contour, and seed).

NOT-FOR: unipotent (pure-log) branch data — the projector removes it by
construction; MUM-point Frobenius bases (Eichler.jl cy_transport); anything
requiring a symbolic proof of the closed form (output is a certified value; closed
forms enter via PSLQ + gate).

INVOKE: generic door `python3 tools/coalescer/coalescer.py --op my_operator.json
[--alpha 3/2] --dps 200 [--out F] [--sb S] [--sb2 S2] [--lift L] [--sdec a,b]
[--nthr N] [--nser N] [--anchor-law roots|seed] [--check-only]` (Route A; the
JSON schema is under INPUTS); banana door `python3 tools/coalescer/coalescer.py
--masses 1,1,1,9 --alpha 3/2 --dps 200 --route A|B [same knobs]` (`--op` and
`--masses` are mutually exclusive); `--alpha` optional (all extracted if
omitted); `--test` = K3 c_{3/2} probe selftest through BOTH doors (≥20d bar,
probe mode).  Anchors: `--sb`/`--sb2`/`--lift` override the defaults 2.2·thr /
4.2·thr / thr/4 (`--op`: 2.2·floor / 4.2·floor / floor/4 with floor the
anchor-law floor in force; `--sb` alone sets sb2 = 2·sb); the anchor law is
checked from the operator's own singularities before any transport (both floors
printed; a violating value refused by name, rc 3); `--check-only` prints the
check and exits 0/3 (seconds after the operator is built or loaded); with
`--route B` the check runs before Julia is dispatched and an integer `--sb`
reaches it as env SB.  Exit codes: 0 ok, 1 gate fail (`--test` below its bar,
or an error), 2 usage (also an unreadable `--op` file), 3 refused by the anchor
law.  Route A API: generic `routeA_monoproj.run_operator(tag, Pj, order, seed_a,
dps=, sb=, sb2=, lift=, sdec_list=, Nthr=, frac=, nloop=, Nser=None, label=None)`
with `Pj, order = op_from_coeffs(table)`, `load_operator_json(path)` for the file
form and `make_series_state(seed_a, order)` for the seed's θ_z-state; banana
wrapper `routeA_monoproj.run(tag, masses, dps=, sb=, sb2=, lift=, sdec_list=,
Nser=, Nthr=, frac=, nloop=)` (find_minimal_op + BFKNS seed, then run_operator).
Route B: `julia routeB_monoproj.jl` — it activates the Eichler.jl project, which
ships in this repository under `upgrades/Eichler.jl/` (instantiate it;
`EICHLER_PROJECT` env overrides the default relative path). PF data for Route B
emitted by `dump_pf_data.py` → pf_*_data.jl (+ the matching .json, which is the
`--op` layout: `dump_pf_data.py --msq 1,1,1,9 --out DIR/k3.jl` then
`coalescer.py --op DIR/k3.json` is the round trip the selftest runs).

INPUTS: generic door, one JSON file (`--op`; key aliases in brackets):
`"theta_form_Pj"` [`"theta_coeffs"`, `"Pj"`] = `{"0": [c00, c01, …], …, "r": […]}`
or a list of r+1 lists — the exact coefficients of P_j(z) (ints, `"p/q"` strings
or `[p, q]` pairs), LOW → HIGH power of z, for L = Σ_j P_j(z) θ_z^j, θ_z = z d/dz;
`"bfkns_a"` [`"seed"`] = `[a0, a1, …]` — the exact Taylor coefficients of the
holomorphic solution f(z) = Σ aₙ zⁿ at z = 0 whose connection coefficients are
wanted (a few hundred terms: the anchors sit at |z_b| = 1/s_b inside its disk and
the run prints the truncation self-check in digits); optional `"order"` (checked
against the table), `"alpha": "p/q"` (the exponent to headline — every fractional
exponent at z = ∞ is extracted regardless), `"label"`, `"seed_radius"` (the
seed's convergence radius R in z: enables `--anchor-law seed`, floor s_b > 1/R;
without it only the roots law exists, floor 1/min|z_sing| over the leading
polynomial's roots), and run knobs `"sb"`, `"sb2"`, `"lift"`, `"sdec": [..]`,
`"nthr"`, `"nser"` (command-line flags override them).  Coordinates: put the point
where your seed converges at z = 0 and the branch point at z = ∞ (a Möbius map of
your operator if need be); t = 1/z is the local coordinate there, s = −t the
transport coordinate, Φ_α(t) = t^α Σ aₙ tⁿ with a₀ = 1 from the exact-ℚ Frobenius
recursion, evaluated on the branch arg t = −π at the Euclidean matching points
s = s_dec > 0; the operator must have no obstruction (no log) on the fractional
branches.  Banana door: m² tuple; conventions (pf_engine.py header): z=1/p² MUM
coord, s=−p² ⇒ z=−1/s, θ_z=−θ_s; BFKNS ϖ₀(z)=Σ aₙzⁿ, aₙ=(−1)ⁿcₙ; threshold p²=0
at z=∞, MUM at z=0 — i.e. the banana is the generic door with the operator from
find_minimal_op and the seed from the multinomial series.  Either door: target
exponent α=p/q optional.

OUTPUTS: per-sdec c_α (re/im strings) JSON; GATE.json records closed forms +
cross-route/sdec-independence digit tables. Validated: K3(1,1,1,9)
c_{3/2}=−√3/(36π) at 195d; CY₃(1,1,1,1,16) c_{5/4}=−1/(√(2π)Γ(¼)²),
c_{7/4}=−5Γ(¼)²/(384√2π^{5/2}) at ≥180d, product 5/(768π³).

ENV: Route B knobs: PFDATA (default pf_cy3_data.jl), PREC (bits, default 512),
SDEC (rational string, default 1/4), NSER (2000), NTHR (200), NLOOP (48), OUT
(ROUTEB.json), EICHLER_PROJECT (Eichler.jl project path).

GATES: two-method gate mandatory before trusting a value — gate_compare.py takes
--family {K3,CY3,CY4} or --masses/--alphas with --routeA/--routeB files (cross-route
gate = the minimum agreeing digits over every Route A x Route B pair, bar 30 d
mandatory, 80 d the target; s_dec-independence per route; every sample against the
closed form when one is recorded; the K3 (1,1,1,9) control on its own files; one
route alone is FAIL by name). PSLQ closure only through pslq_calpha.py = the
lockpick.pslq_gate engine in the sealrun form, never raw mp.pslq: the basket is
declared from the ring table (calpha_rings.py) before any value is read,
members_audit runs before every scan, the controls sit inside the scan basket, the
height ladder is capacity-checked at two precisions, every hit is re-verified on
digits no fitting leg saw, and the found vector is compared with the one derived
from the closed form; the engine file must be a listed sha and its protocol units
must match the pin (rc 3 by name otherwise). The CY3 record is the fixture
(tests/test_pslq_gate_route.py): gate 72.53 / 71.02 d cross-route, the reference PSLQ
vectors refound at the record's dps pair. The K3 positive control runs the SAME code
path (pf_engine imports nothing K3-specific — the control is an honest re-derivation).

Battery: `coalescer.py --test` → PASS rc=0, run from a scratch cwd — K3 c_{3/2}
probe through both doors: the banana door in-process matches the closed form to
24.7 digits (≥20 bar), and the same operator + seed re-exported by
dump_pf_data.dump to the `--op` JSON and run through `coalescer.py --op` in a
child process (overlapped with leg 1) matches it to 24.7 digits and the banana
door to 25 digits; and one operator that is not a banana — the shipped
`examples/gauss_2f1.json` (Gauss ₂F₁(¼,¾;1;z): θ² − z(θ+¼)(θ+¾), Pochhammer seed,
radius 1) through `coalescer.py --op` in a second child, whose c_{1/4}, c_{3/4}
match Gauss's connection formula c_a = e^{iπa}Γ(b−a)/(Γ(b)Γ(1−a)) to 25.4 digits
at both matching points; ~35 s wall (the three legs overlap), pure-python Route A
(no Julia required for the selftest); scratch in a mkdtemp() directory (honours
TMPDIR), removed on exit.  `examples/gauss_2f1.json` doubles as the template for
your own `--op` file.
`python3 -m pytest tools/coalescer/tests -rs` → the anchor-law battery (the
three tuples' floors from the operators vs the recorded values, the refusal
rc 3, the override acceptance, the CLI exit codes, the `--op` door's check-only
/ refusal / usage legs on a freshly dumped K3 JSON; ~60 s; the K3 `--test`
masked-identical leg is skipped by name unless COALESCER_SLOW=1).

WORKED EXAMPLE (banana; project data, kept as the package's validation record,
not needed to run your own operator): pf_k3_data.jl, pf_cy3_data.jl,
pf_cy4_data.jl (Route B operator/branch/seed tables written by dump_pf_data.py),
GATE.json (CY₃ closed forms and digit tables, K3 control), tests/fixtures/
(Route A/B outputs, GATE_CY4.json, PSLQ_CY4.json, PSLQ_PRODUCTION.json, PINS.json).

FOOTGUNS:

- Anchor s_b must satisfy s_b > 1/min|z_sing| (MUM convergence radius), NOT
  1.5·max|z_sing| — an Arb ball stays TIGHT on the WRONG sum if the seed series is
  evaluated outside the disk.  Units: z = −1/s, the seed is Σₙ aₙ z_bⁿ at
  z_b = −1/s_b.  The CLI computes the floor from the operator it builds and
  refuses a violating anchor by name (rc 3) before anything is transported; two
  readings of min|z_sing|, both printed on its `[anchor-law]` line:
  `seed` (default) R = 1/thr, thr = (Σ√m²)² — the leading threshold z = −1/thr
  is the seed series' nearest true singularity (exact for the multinomial
  series: thrⁿ/(n+1)^(L−1) ≤ |aₙ| ≤ thrⁿ), floor s_b > thr, so the written
  defaults 2.2·thr / 4.2·thr always pass; `roots` R = min|z| over ALL roots of
  the leading polynomial (the list Route B's auto-anchor ceil(1.6/min|z_sing|)
  uses) — conservative, it counts apparent singularities: CY₃ (1,1,1,1,16) has
  one at |z| = 0.0056 (floor 178.6, above the default 140.8 and the s_b = 100
  of the CY₃ gate, where the seed converges to full precision), CY₄
  (1,1,1,1,1,25) a pair at |z| = 0.0013 (floor 754.8, above the defaults
  220/420, where the seed agrees to 40 of 40 digits between 600 and 400
  terms).  Below thr the seed diverges (CY₄ at s_b = 90: −8.7 digits).
- `--op` seeds: give EXACT coefficients (ints / "p/q" / [p, q]); a float with a
  fractional part is refused by name.  Enough terms that the truncation
  self-check at s_b prints ≥ dps digits (the run says how many it got); anchors
  outside the seed's disk are the classic silent failure (next item) — declare
  `"seed_radius"` when you know it, else the conservative roots floor governs.
- The CLI's Route B does not pass the tuple to Julia: PFDATA (env) selects the
  data file and must be this tuple's (dump_pf_data.py); the `[anchor-law]` line
  is computed for `--masses`, and `--sb` reaches Route B as env SB (integer).
- CONVENTION DEPENDENCE (order-4 monodromy): at an order-4 threshold
  (quarter-integer exponents, eigenvalues ±i) the SIGN of c_α is
  convention-dependent: it is fixed only jointly with the Φ_α branch normalization
  and the monodromy orientation, and Route A (CCW in s) vs Route B (CW in z) differ
  in orientation. Quote signs only relative to the pinned pf_engine convention
  block, and cross-route gate before trusting any sign. (Order-2 branches, λ=−1,
  are orientation-insensitive.)

RELATED: Route B rides Eichler.jl (`_theta_transport`) — vendored in this repo
under upgrades/; PF operators from pf_engine's mod-p+CRT builder (or the
Annihilator tool upstream); the PSLQ layer follows the Lockpick gate protocol.

CREDIT: banana Picard–Fuchs operators, Frobenius bases and series conventions follow
Kilian Bönisch, Fabian Fischbach, Albrecht Klemm, Christoph Nega and Reza Safari [BFKNS]
and Bönisch, Duhr, Fischbach, Klemm and Nega [BDFKN], on whose monodromy analysis the
spectral projector builds (it generalizes the K3 positive-control projection used in
this package). Certified transport by Arb [Arb]; closed forms via PSLQ [PSLQ] under the
lockpick gate. Bracketed keys resolve in REFERENCES.md at the repository root.
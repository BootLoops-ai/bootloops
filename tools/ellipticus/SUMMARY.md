# Ellipticus — summary

Certified evaluator for elliptic multiple polylogarithms and iterated integrals
on an elliptic curve: given an exact curve/family spec, letters, a word, a
point and a precision, it returns certified values. One Python front door
(import name `ellipticus`, formerly `empl_eval`) over three certified engines.

- **engine i (march, workhorse)** — exact inhomogeneous-PF extended systems for
  window moments of any cubic/quartic family y^2 = F(x;z) (coefficients rational
  in z), Bezout/exact-form construction (`extend.moment_system`,
  battery-verified EXACT vs the pinned quartic AND cubic demonstration systems)
  + certified adaptive-Taylor transport with dlog-letter states, tangential-base
  (shuffle) regularization at analytic bases (perfect-square and
  square-x-linear seed charts), complex-waypoint detours (`march.System`).
- **engine q (fast tau-words)** — the bundled Kronecker-Eisenstein kernels
  (`vendor/`: g^(n) three cross-checking routes, omega_k) +
  AGW shuffle-regularized tau-words (`itint`) + BMSW sunrise closed forms
  (`vendor/rowscripts_empl.py`), with a PROVEN qbar tail bound added at
  API level (`qengine.g_tail_bound` / `g_n_certified`) and theta-route
  fail-closed fallback.
- **engine agm (weight-1/period layer)** — exact curve layer (`curve.QuarticCurve`),
  Legendre-K/AGM oval periods, sin^2-substitution certified cycle moments, and
  the tracked-radical complex-pair CONTOUR cycle verb
  (`cycle_pair_contour`, gamma1 contour-pair class, battery leg L3).

PURITY BAR: exact inputs only (floats refused). Two-precision self-gate in
`api.evaluate`. Loud scope refusals outside the certified word classes.

Members:

- **ellred** (`ellred/`) — standard-elliptic (Legendre/Carlson) reduction of
  one-fold elliptic integrals int rat(x) dx / sqrt(cubic|quartic): quartic
  chart + closed F/E/Pi atoms + J-moment recurrence, cubic 3-real-root and
  1-real-2-complex moment charts, an x-space reducer to F/E/Pi + algebraic,
  an analytic 1r2c reducer with exact-derivative certification, and a Carlson
  R_F/R_D/R_J route cross-check. Battery leg L5.
- **gmtel** (`gmtel/`) — certified Gauss-Manin q-telescoper: exact
  Picard-Fuchs operators of K3 pencil periods from a Weierstrass model (exact
  fiber Gauss-Manin + creative telescoping with an explicit certificate
  row-vector), two-route (geometric telescoper vs data-side annihilator fit)
  positive control against a packaged reference L5 operator
  (fixtures/gmtel/L5_theta.txt); own selftest (gmtel_selftest.py, GUIDE
  MEMBERS section).

Battery (battery.py; receipts JSON): checks the code paths end-to-end against
INDEPENDENT references shipped in fixtures/battery/ — march kernels + dlog-letter
words vs direct-quadrature oracles of the defining integrals (two synthetic
families, exact polyform pins), a-cycle period vs the closed form 2*K(1/4) +
cycle moments vs quadrature, complex-pair contour Om/X1/X2 vs a branch-cut
segment-integral oracle (orientation sign recorded), sunrise rows 09/10/11
vs bundled held-out engine gate points + the equal-mass classical control,
and the ellred atom/chart classes vs direct quadrature + pinned modulus
strings (leg L5, fixtures/battery/ellred_controls.json).

Selftest (selftest.py): Legendre-vs-quadrature + K(1/4) closed form; EXACT
polyform match vs the pinned demonstration system; march vs independent
x-quadrature; three-route g^(n) cross-check + tail-bound coverage + equal-mass
control.

Related: tools/eichler for modular-form periods (Gamma_1(6)
iterated-Eisenstein/modular words with ball-certified enclosures go to the
Eichler.jl engine, `iterated.jl` in upgrades/Eichler.jl, reachable from Julia
through `tools/gpl-eval/GPLEval.jl/src/empl_bridge.jl`; that engine is not
duplicated here); tools/gpl-eval for ordinary MPLs.

Manual: GUIDE.md in this directory. Tool page:
https://bootloops.ai/tools/ellipticus.html

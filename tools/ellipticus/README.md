# Ellipticus

Ellipticus is a certified evaluator for elliptic multiple polylogarithms and iterated integrals on an
elliptic curve. Given an exact curve or family spec (a quartic or cubic y^2 = F(x;z) with coefficients
rational in z), letters, a word, an evaluation point and a precision, it returns the values together with
a certificate: proven per-step Taylor tails on the transport engine, proven q-tail bounds on the
Kronecker-Eisenstein engine, and two-precision plus method agreement throughout. Inputs are exact
rationals only; floats are refused at the API boundary.

Three engines sit behind one Python front door (`from ellipticus import MomentFamily, QuarticCurve,
evaluate, qengine` with `tools/` on PYTHONPATH): exact inhomogeneous Picard-Fuchs systems for window
moments of any cubic/quartic family plus certified adaptive-Taylor transport with dlog-letter words
(`march`), Kronecker-Eisenstein kernels and tau-iterated integrals with a proven q-tail bound and the
unequal-mass sunrise closed forms (`qengine`), and an exact curve layer with AGM/Legendre periods, cycle
moments and complex-pair contour cycles (`curve`). Two members ship inside: `ellred/`, standard-elliptic
(Legendre/Carlson) reduction of one-fold elliptic integrals with closed F/E/Pi atoms, and `gmtel/`, a
certified Gauss-Manin q-telescoper for exact Picard-Fuchs operators of K3 pencil periods.

Requirements: python3 >= 3.10 with mpmath, sympy and python-flint (gmtel also needs numpy). Nothing to
build; the supporting engines are vendored under `vendor/` and the reference data under `fixtures/`.

Self-test (four fast checks, ~15 s): `python3 selftest.py` from this directory. Battery (five legs
against independent references, ~20 s; receipt JSON to `$ELLIPTICUS_BATTERY_OUT`):
`python3 battery.py` here, or `PYTHONPATH=tools python3 -m ellipticus.battery` from the repository root.

Start with `GUIDE.md` (purpose, when to use it, invocation, environment, checks, members, pitfalls);
`SUMMARY.md` is the one-page overview. The import name was formerly `empl_eval`; the former
`EMPL_EVAL_*` environment variables are still read, with the `ELLIPTICUS_*` names taking precedence.

Related packages: `tools/eichler` for modular-form periods (with the Eichler.jl engine in
`upgrades/Eichler.jl`); `tools/gpl-eval` for ordinary MPLs.

Tool page: https://bootloops.ai/tools/ellipticus.html

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).

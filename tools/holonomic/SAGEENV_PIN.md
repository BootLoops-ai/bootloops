# SAGEENV_PIN — the reference environment for the certified-transport recipe

Every recorded number in this package — the worked-example enclosures in
`reference/`, the rank-ladder walls in `GUIDE.md`, the smoke's expected
behavior — was measured on ONE reference build. This file pins that build so
the numbers stay interpretable and the environment is reproducible.

## 1. Versions

    SageMath 10.7 (Release Date: 2025-08-09)
    ore_algebra 0.5   (installed into the Sage venv; built from source as a
                       cp311 wheel — Python inside that Sage is 3.11)

NOTE: in this build, bare `sage -python -c "import ore_algebra"` FAILS
(`ImportError: cannot import name Category`) — `sage.all` must be imported
first (as every script in this package does, and as `holonomic_transport`
does at module top).

## 2. Invocation pattern

    sage -python your_script.py

with

    from sage.all import QQ, PolynomialRing, RealBallField
    from ore_algebra import OreAlgebra
    import holonomic_transport as H

Recipe conventions the house layer hard-codes (worth copying verbatim into
any hand-rolled driver):

- exact rational/dyadic path points (pin dyadics k/2^30 via
  `H.dyadic_point` so every leg — including any Float64-coercing consumer —
  evaluates at IDENTICAL exact points);
- the transition-matrix IC convention (derivatives vs Taylor coefficients)
  is DETERMINED per session on Du^3 (solutions 1, u, u^2) before use — do
  not assume it (`H.determine_ic_convention`);
- leading-coefficient REAL roots are checked against the path segment; if a
  root lies on the segment, detour through the complex midpoint + I/8
  (`H.plan_path`);
- apply the ball matrix to the IC vector and re-wrap the value as a ball:
  a REAL ball when its imaginary part is exactly zero, a COMPLEX ball
  otherwise (a detoured continuation across a real singularity is genuinely
  complex — never discard the imaginary part).

## 3. algorithm='naive' is MANDATORY in the reference build — and why

The ore_algebra compiled (Cython) extension is binary-incompatible with the
reference Sage 10.7 build. Measured live:

- `algorithm='naive'`: emits `UserWarning: Cython extensions not found.
  Falling back to slower Python implementation.` and then WORKS — pure-
  Python ball summation, certified result correct.
- default algorithm (no keyword, binsplit path): RAISES
  `TypeError: C variable sage.rings.integer._small_primes_table has wrong
  signature (expected int [500], got int [0x1F4])` — hard failure, no
  fallback.

So in that build 'naive' is not a preference, it is the only working path,
and the house layer hard-wires it. Cost: the pure-Python summation is the
slow lane — see the measured rank ladder in `GUIDE.md` and the benchmark
harness in `ladder/`.

## 4. Smoke check (re-runnable)

    sage -python smoke_test.py receipt.json

Reference run: L = (1+x)Dx^2 + Dx, transition matrix over [0, 1/2] at eps
1e-40, algorithm='naive' -> entry (0,1) =
`[0.405465108108164381978013115464349136571990 +/- 4.24e-43]`,
|entry - ln(3/2)| <= 2.23e-47, certified radius 2.0e-47, 0.35 s transport.
The default-path probe raised the TypeError quoted in Sec. 3. The recorded
enclosures are stored in `reference/ln32_reference.json`, and the battery's
engine-free leg re-verifies them against an independent pure-Python
computation of ln(3/2) on every run.

## 5. Rebuild recipe

SageMath 10.7 (or current) + `sage -pip install ore_algebra` (reference
version 0.5). The Cython-extension breakage is a property of that exact
build pairing; a fresh build may restore the fast (binsplit) path — the
smoke's S5 probe reports that as ENV-FACT-CHANGED rather than passing
silently, because ALL recorded numbers here assume 'naive'. Re-measure
before relying on the fast path. To MOVE an existing pinned build instead
of rebuilding, use `relocate_sageenv.sh` (see GUIDE.md, "Relocating the
environment") — rsync + shebang patch, gated by smoke and
origin-vs-relocated bit-identity.

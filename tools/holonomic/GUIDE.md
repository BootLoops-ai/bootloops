# Holonomic — certified analytic continuation (external-engine port)

Tool page: https://www.bootloops.ai/tools/holonomic.html

KIND: package (`holonomic_transport.py`, imported under `sage -python`).
EXTERNAL-ENGINE PORT — the engine, and all of the certification mathematics,
is [ore_algebra](https://github.com/mkauers/ore_algebra) (Marc Mezzarobba's
`ore_algebra.analytic`: `numerical_transition_matrix` / `numerical_solution`).
This package is the house layer only: the debugged driving recipe around the
engine, plus the reference-environment pin (`SAGEENV_PIN.md`).

PURPOSE: certified numerical analytic continuation of holonomic ODEs with Arb
ball enclosures. Given an operator L = p_r(x) Dx^r + ... + p_0(x) over Q[x]
annihilating f, and the data (f, f', ..., f^(r-1)) at a base point, transport
that data to a target point along a declared path in the complex plane. Every
output entry is a ball — midpoint plus rigorous error radius — so the result
is a mathematical statement that the true continued value lies inside a
stated interval. That is a stronger warranty than agreement between two runs
at different precision, which is evidence produced by the same code being
checked.

USE-WHEN:

- A continued value must carry a CERTIFIED (ball) enclosure, not a
  two-precision agreement.
- An independent oracle leg is needed to spot-check the house transport
  evaluators (`wayfinder/`, `famhar/`, `baller/`'s transport module): at
  selected points the march must land inside a ball produced by code it
  shares nothing with. Production continuations run this way have been
  checked against an independent Arb evaluation to 48+ certified digits.
- Continuation across or around singular points of a holonomic operator
  where a rigorous error bound is required (regular-singular landings are
  ore_algebra-native, via the local Frobenius basis).

NOT-FOR: ⛔ INDEPENDENCE FENCE — never merge this package into `wayfinder/`,
`famhar/`, or `baller/`, and never fold their code in here: those evaluators
are certified BY GATING AGAINST this oracle; one codebase playing both roles
would leave every cross-check grading the code against itself and would
retroactively invalidate every past check that used both sides. They stay
separate tools that gate each other, permanently. Also not a bulk evaluator:
`algorithm='naive'` is pure-Python summation (see the rank ladder below) —
march thousands of points with the house evaluators and certify a handful
here. Does not acquire operators (`annihilator/` and `dipstick/` do) and
does not certify inputs — the ICs and the operator need their own gates.

INVOKE: `sage -python your_script.py` with

```python
from sage.all import QQ, PolynomialRing, RealBallField
from ore_algebra import OreAlgebra
import holonomic_transport as H

Rx = PolynomialRing(QQ, "x"); x = Rx.gen()
A = OreAlgebra(Rx, "Dx"); Dx = A.gen()
RBF = RealBallField(430)

L = (1 + x)*Dx**2 + Dx                      # annihilates ln(1+x)
res = H.certified_transport(L, QQ(0), QQ(1)/2, [RBF(0), RBF(1)], eps=1e-40)
res["value"]          # ball containing ln(3/2), certified radius ~3e-42
```

One call: `H.certified_transport(L, u0, u1, ics_der, eps=1e-48, prec=430)`
-> `{value, vector, path, detoured, conv, enclosure_ok}`. Pieces:
`H.determine_ic_convention()`, `H.ic_vector()`, `H.plan_path()`,
`H.transition_matrix()` (`algorithm='naive'` hard-wired),
`H.dyadic_point()`, `H.enclosure_too_wide()`, `H.engine_versions()`.

Requirements: SageMath (https://www.sagemath.org or conda-forge; reference
build 10.7) with ore_algebra installed into it (`sage -pip install
ore_algebra`; reference version 0.5). There is no plain-CPython route — the
module imports `sage.all` first, which is mandatory (bare `import
ore_algebra` fails with a Category ImportError in the reference build).
`HOLONOMIC_SAGE` overrides the `sage` executable the battery invokes.

The recipe the house layer hard-codes (see `SAGEENV_PIN.md` for the full
statement): exact rational/dyadic path points (`dyadic_point`, k/2^30, so
float-coercing legs see IDENTICAL exact points); the transition-matrix IC
convention (derivatives vs Taylor coefficients) DETERMINED per session on
Du^3 (`determine_ic_convention` — never assume it); leading-coefficient REAL
roots checked against the path segment, on-segment root => complex-midpoint
detour + I/8 (`plan_path`); ball matrix applied to the IC vector and the
value re-wrapped as a ball — REAL ball when its imaginary part is exactly
zero, COMPLEX ball otherwise (a detoured continuation across a real
singularity is genuinely complex; the imaginary part is never discarded).

INPUTS: OreAlgebra operator L in D over QQ[x]; exact rational/dyadic base and
target points (pin via `dyadic_point(x, bits=30)` when a Float64-coercing leg
must see identical points); ICs as RealBallField balls [f, f', ...,
f^(r-1)] at the base point.

OUTPUTS: ball value (real ball for a real result, complex ball for a
genuinely complex one — e.g. any detoured continuation across a real
singularity) + full transported vector (raw engine balls) +
path/convention/enclosure metadata dict; the smoke writes a JSON receipt.

GATES (before trusting any output):

- The continuation is certified; the INPUTS are not — gate the ICs and the
  operator independently (pattern: check the ICs against an independent
  evaluator at the base point BEFORE transporting, to 40+ digits).
- Ball CONTAINMENT is the certified statement; check `enclosure_ok` and the
  radius, never midpoints alone.
- Record `engine_versions()` in every receipt — all recorded numbers assume
  Sage 10.7 / ore_algebra 0.5 with the naive path.
- The IC convention is DETERMINED per session, never assumed.

Battery: `python3 tools/holonomic/battery.py` (from the repo root). The
engine-free `reference` leg runs everywhere: it recomputes ln(3/2)
independently in pure Python integer arithmetic and verifies the recorded
worked-example enclosures in `reference/ln32_reference.json` contain it and
are internally consistent. The `smoke` leg runs `smoke_test.py` under
`sage -python` (live certified transport must contain ln(3/2) with radius
<= 1e-38, detour law, convention probe, and the detoured-VALUE gate: the
transport of ln(1+x) from 0 to -3/2 must return a COMPLEX ball containing
ln(1/2) + i*pi — the imaginary part must survive) and SKIPs by name when
SageMath or ore_algebra is absent. `--full` adds the rank-6 ladder leg
(roundtrip ball contains the identity; transported known solution contains
the direct-series value).

Relocating the environment: `bash tools/holonomic/relocate_sageenv.sh
SRC_ENV DEST_DIR tools/holonomic` copies a pinned Sage+ore_algebra venv
into `DEST_DIR/sageenv` and gates the copy before any science runs on it —
five steps: (1) rsync (source stays read-only), (2) stale-shebang patch in
the relocated copy only, (3) direct invocation of the copy's own versioned
bin python, bypassing the `bin/sage` wrapper, (4) smoke gate
(`smoke_relocated_env.py`: imports, versions-vs-pin, sys.prefix check,
IC-convention probe, toy transition matrix), (5) bit-identity gate
(`verify_bit_identity.py` run on the origin AND the relocated env: the
full-precision serializations of a fixed toy certified transport must be
EQUAL — deterministic 'naive' summation makes bit-identity, not mere
closeness, the gate grade). The two gate scripts read `HOLONOMIC_DIR` (dir
containing `holonomic_transport.py`) and `HOLONOMIC_SAGEENV` (expected env
root, checked against `sys.prefix`) from the environment themselves; the
transport module reads neither. Unlike `smoke_test.py`, which reports a
different build as ENV-FACT-CHANGED, `smoke_relocated_env.py` HARD-asserts
the pinned versions (Sage 10.7 / ore_algebra 0.5) — a relocation must
reproduce the pinned build exactly. The helper emits a receipt JSON whose
verdict fields are script-emitted and ends with
`RELOCATION_VERDICT PASS|FAIL`.

FOOTGUNS:

- Measured ordinary-path rank ladder (eps 1e-40, single core, maxrss <= 364
  MB; harness in `ladder/`): R6 0.02 s / R24 2.5 s / R48 96.4 s / R72
  1292.9 s on the reference build — price big ranks from the ladder before
  committing to them.
- Near-singular paths blow enclosures at high rank: a rank-48 path passing
  x=1 at distance 1/32, eps 1e-40 -> max entry radius 1.4e+18
  (ENCLOSURE-TOO-WIDE — a valid ball that certifies nothing; rank-24 at eps
  1e-20 -> 2.7e+08). Route with real clearance, or land on the singular
  point in the engine's Frobenius basis; `enclosure_too_wide()` flags it.
- High rank degrades enclosures on ORDINARY paths too: measured on the
  reference build, rank-48 at eps 1e-20 returned max entry radius 8.0e+23
  against max entry |value| 4.5e+7 (a valid ball that certifies nothing),
  and rank-72 at eps 1e-40 returned max entry radius 769 against max entry
  |value| 5.1e+12 (~10 relative digits, far short of the 40 requested) —
  read the achieved radius off the result; never assume the requested eps
  was met.
- In the reference build the DEFAULT algorithm hard-fails (TypeError) and
  `algorithm='naive'` is the only working path — see `SAGEENV_PIN.md`. A
  build where the compiled path works is reported by the smoke as
  ENV-FACT-CHANGED; the recorded numbers still assume 'naive'.
- The continuation certificate does not certify the operator or the ICs — a
  wrong operator transports the wrong function inside a beautiful ball.

RELATED: upstream operator acquisition = `annihilator/` (series -> PF
operator) + `dipstick/` (order/rank certification); independent house
transport legs this oracle gates (never merge, see NOT-FOR) = `wayfinder/`,
`famhar/`, `baller/`, `eichler/` (the analogous transport role for elliptic
periods).

CREDIT: ore_algebra — Manuel Kauers, Maximilian Jaroschek, Fredrik Johansson,
*Ore Polynomials in Sage*, arXiv:1306.4263. The certified-continuation
engine (`ore_algebra.analytic`) — Marc Mezzarobba, *Rigorous
Multiple-Precision Evaluation of D-Finite Functions in SageMath*, ICMS 2016,
arXiv:1607.01967. The ball arithmetic underneath — Arb, F. Johansson, IEEE
Trans. Comput. 66 (2017) 1281, arXiv:1611.02831 (now part of FLINT). The
path discipline, convention probe, and environment pin are this package's
additions around the engine; the certified mathematics is all the engine's.

DEPENDENCY LICENSES: this package is MIT-licensed. It drives SageMath and
ore_algebra (both GPL), which you install separately (see SAGEENV_PIN.md); no
GPL code is vendored here, and a relocated sageenv built by
relocate_sageenv.sh is a local copy of that GPL environment, never part of
this repository.

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).

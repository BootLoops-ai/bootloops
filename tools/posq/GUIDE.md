# posq — certified Bayesian evidence by exact positive quadrature

Tool page: https://bootloops.ai/tools/posq.html

posq is a Python package with a C kernel (`posq_kernel.c`, flint/arb; built from
source on first use). It also ships as a member of the `baller`
package in this repository: `baller.quad.posq` imports the files in this
directory. Do not edit the engine files in place — any change must re-run the
live verification battery (`verify_posq_adaptation.py`, below) before its
outputs are used.

## Purpose

posq computes a Bayesian evidence integral as a certified positive value-sum: a
degree-matched exact Gauss sum evaluated directly on the raw likelihood
product, never expanding coefficients and never subtracting, so two-sided
machine-width certificates follow from positivity alone. The deliverable is a
certified interval `[lower, upper]` for the evidence value, and certified
Bayes-factor margins between models on the same data.

The method, in four steps:

1. **Exact prior marginalization.** Substitute `u_k = exp(-beta*s_k)`.
   Exponential clock-increment priors marginalize exactly into
   `Z_s(p) = (c^3/2) * INT_{[0,1]^3} u^(c-1) F(u; p) du` with
   `c = 2*lambda*p*(1-p)`, where `F` is the raw positive product of pattern
   probabilities — a polynomial of known multidegree.
2. **Exact quadrature in the clock dimensions.** A Gauss–Jacobi tensor rule
   built at that degree (interval-Newton node certification at 2048 bits;
   acceptance requires strictly positive weight balls, exact-moment
   reproduction, and monomial exactness) integrates `F` exactly: the clock
   dimensions contribute zero quadrature error — no cells, no remainder
   estimate.
3. **Positivity gives both sides.** Positive weights times a positive integrand
   make the streamed ball sum cancellation-free, so the certified enclosure is
   two-sided at machine width. Measured on the production 48.7-million-node
   sweep: relative width `2^-179` at 192-bit working precision, 0.201 CPU-h.
4. **The one remaining dimension** (`p`) is integrated by a non-adaptive
   Cauchy-ellipse register: certified Gauss–Legendre nodes per panel, one
   exactly integrable Bernstein-basis absolute-majorant sweep per panel, and
   exact endpoint bands.

Measured production record (four-taxon binary CTMC, strict clock, 318
site-pattern counts; the production data and result files are not shipped —
see below):

- One complete two-sided evidence row: `lnZ` width 0.0627 nats at
  27.4 CPU-h per topology.
- Full 15-topology table: all rows certified two-sided, widths 0.0095–1.157
  nats, 93.25 CPU-h total; cross-check controls passed and all 3 planted count
  corruptions fired.
- Certified decision sentences (gap = winner's lower endpoint minus rival's
  upper endpoint): maximum-evidence topology `(((A,D),C),B)`,
  `lnBF >= 0.4675` over the runner-up and `>= 5.21` over the losing split
  classes, agreeing with an independently computed model's ranking.

`DESCRIPTION.md` has the paper-register description of the method; the
positioning against the nearest lines of work, with references, is on the tool page.

## When to use

The object class is likelihoods of continuous-time Markov chains on small
graphs with clock priors:

- a small state graph (here: a binary CTMC on a rooted 4-leaf tree, strict
  clock),
- prior marginalization that lands on polynomial-weight integrals (Exp/Gamma
  clock increments give `u^(c-1)`-class weights after the exponential
  substitution; uniform or polynomial-density priors on rate/frequency
  parameters),
- an integrand that is a raw product of pattern probabilities — positive on the
  domain, polynomial of known multidegree in the transformed variables.

Reach for posq when you need a certified **two-sided** evidence value or
Bayes-factor margin and interval or Monte-Carlo routes either give no lower
bound or blow up. The two failure mechanisms this design exists to avoid are
coefficient-expansion cost and coefficient-expansion cancellation: posq
evaluates values, never expanded coefficients.

Limits (all measured, none optional):

1. **Four-leaf only.** The exact clock-dimension collapse is proven and priced
   for rooted 4-leaf shapes only (3 clock dimensions; 12 caterpillar shapes,
   plus 3 balanced shapes via a min/diff decomposition). Six-leaf trees are not
   supplied, and corrected/ascertained likelihoods are not covered. New object
   classes need the degree bookkeeping and prior-to-weight mapping re-derived,
   the Gauss rule rebuilt at the object's own multidegree, and the exactness
   checks re-asserted — nothing transfers by analogy.
2. **The `p` layer is the cost center**, and its prices are measured numbers,
   not convergence claims. The shipped ellipse register costs 27.4–38 CPU-h per
   topology, essentially independent of the width target (0.1, 0.5, and 1.0
   nats all cost about the same — if it is affordable at all, it is a 0.1-nat
   instrument). Its panel counts are set by the absolute majorant, which pays a
   measured penalty over signed sensitivity (effective rate ~990 vs ~230 per
   unit of `p`) — do not re-derive panel counts from signed sensitivity. Two
   alternatives were measured and rejected: a two-exponent sandwich is sound
   but costs 262 core-h per topology, and Taylor-jet remainders over a wide
   panel fail at production degree (interval Newton rides the explosive second
   solution of the three-term recurrence near panel endpoints). Center jets —
   certified derivatives at a point — remain a cheap, sound primitive.

## How to run

What ships in this directory:

| file | role |
|---|---|
| `posq_kernel.c` / `posq_kernel` | C sweep kernel (flint/arb): certified positive Gauss–Jacobi tensor sweeps, caterpillar/balanced/multi-sweep modes, per-node value evaluation `P_y = A_y + B_y*v`, exact ball serialization in and out |
| `kernel_io.py` | exact Python-to-kernel ball serialization (mantissa/exponent 4-tuples; radius rounds up) |
| `bench_rung1_gauss.py` | certified Gauss–Jacobi rule builder (interval-Newton at 2048 bits with receipt checks), Python reference sweep, exact-rational N=8 surrogate |
| `topos.py`, `exact_polys.py`, `exact_polys_bal.py`, `bal_fiber.py` | exact structure: 15-topology enumeration, pattern polynomials, balanced min/diff decomposition |
| `posq.py` | package surface: `SurrogateExact` (exact closed-form surrogate with exact-derivative jets) and `PosqKernelAdaptation` (the kernel as a certified point/enclosure engine) |
| `derive_posq_sentences.py` | rounding discipline for certified sentences (floor the winner's lower endpoint, ceil the rival's upper endpoint) |
| `closure_interface.py`, `closure_fixtures/`, `selftest_closure_interface.py` | read-only evaluation interface to a pinned upstream build; see `CLOSURE_INTERFACE_README.md` |
| `CLOSURE_INTERFACE_MANIFEST.sha256` | sha256 pins for the closure-interface files and their upstream build |

Requirements: Python 3 with `python-flint`, plus a C compiler and the flint/arb
development headers for the kernel. No binary is shipped; the build line is:

```sh
gcc -O2 -o posq_kernel posq_kernel.c -lflint -lmpfr -lgmp -lm
```

You do not have to run it by hand: at first use `kernel_io` looks for a
`posq_kernel` binary built beside the source, and if none is present (or it
does not execute, e.g. it was linked against a different FLINT shared-library
version or built for another CPU architecture) it compiles `posq_kernel.c`
automatically into a cache directory (`POSQ_KERNEL_CACHE`, default
`~/.cache/posq`) using the build line above. That needs a C compiler and
the FLINT headers. Set `POSQ_KERNEL` to point the driver at your own build.
With no runnable binary and no toolchain, kernel entry points refuse with
the named error `POSQ-KERNEL-UNAVAILABLE` (and
`verify_posq_adaptation.py` exits INDET), never a wrong answer.

The Python modules import cleanly from a checkout:

```python
import sys; sys.path.insert(0, '<your-checkout>/tools/posq')
import posq

adapt = posq.PosqKernelAdaptation()
Z = adapt.point_eval(("709/2048", "10"))   # certified ball via the C kernel
```

`PosqKernelAdaptation` rebuilds and receipt-checks the Gauss rules at each
point's exact rational `c` before sweeping. Note that its `enclose()` method (a
monotone `c`-bracket enclosure over a `(p, lambda)` box) is verification
scaffolding, valid only for `p` boxes inside `(0, 1/2)` — it is **not** the
production `p` layer, which is the ellipse register described above.

Two hard rule sets are enforced by the code, not left to the user:

- **Balanced-shape exactness checks.** Balanced topologies go through the
  min/diff decomposition (`m ~ Exp(2*lambda)` giving a `w1^(2c-1)` rule; two
  ordered cases summed). Before any balanced sweep, the exactness checks must
  pass: per-character degrees `(4,3,1)`, parity, sum-to-one, agreement of two
  independent exponent-vector constructions, measure factorization, an exact
  `p`-Taylor-shift identity, and per-panel pole caps — the balanced arm's
  nearest pole is `c = -1/2`, and each panel's ellipse radius must sit strictly
  below it on both flanks. The production driver refuses to sweep on any
  failure.
- **Positivity screens.** The assembly is a positive-weighted sum of positive
  balls, cancellation-free by construction, and everything that could break
  that is screened rather than assumed: the only mixed-sign step is the
  per-node value `P_y = A_y + B_y*v`; every `P` ball and every assembled `Z`
  ball must certify strictly positive; per-sweep relative width worse than
  `2^-100` raises a hard error; and corner nodes — which sit within `~n^-2` of
  the boundary, where the mixed-sign evaluation is tightest — are checked
  explicitly at 192-bit. Never expand coefficients.

Not included in this repository: the production drivers (`run_quartet.py`,
`taproot_driver.py`), the production panel-spec/result files, the production
pattern-count data (`QUARTET_COLLAPSE.json` — modules that would read it
refuse by name unless `POSQ_STAGE0` points at a directory containing it),
and the two external dependencies of the verification legs noted below.

## Verification

The packaged engine was run, as a consumer adaptation, through the `eras`
adversarial verification battery, passing at 3 seeds. Each pass comprises:
domain gates (an exact closed form reproduced sha-pinned reference values; the
kernel ball contained the exact rational at relative width `2^-248`; a planted
count corruption fired), then 224/224 closed-box containment checks per shell
across 3 shells (corners, faces, near-face band, fresh batches), all 6 planted
corruptions firing (narrowing, shifting, and both dropped dimensions), the
precision probe caught (a 53-bit rule rebuild fails interval-Newton
certification and correctly returns no information), and finite-difference
cross-checks of an independent exact-calculus derivative path with worst
deviation `3.1e-23` against a `1.1e-12` gate — zero fatals, `result.ok=True`.
That battery is LIVE in this directory: re-run it yourself with
`python3 verify_posq_adaptation.py` (exit 0 = PASS) — the run you make is
the verification record. A quick smoke entry, `python3 posq.py --selftest`,
runs the exact-surrogate pin plus one live kernel-containment leg in
seconds.

Honest limits of that claim:

- The battery exercises the engine at surrogate scale (N=8 counts, rules
  `(17,13,5)`, same builder, kernel, and spec). The production-scale claims
  rest on the 15-topology table's own containment and planted-corruption
  controls, not on the battery.
- A single-seed pass is necessary, not sufficient — vary the seed before
  relying on any new adaptation.
- The `c`-bracket enclosure in `posq.py` is battery scaffolding, not a
  production `p` layer.

If you adapt posq to a new object: run a measured trial first, rerun the
battery against your adaptation before trusting it (follow the
`verify_posq_adaptation.py` pattern), check `result.ok`, and use several seeds.
The four-leaf limit is not lifted by a battery pass.

`verify_posq_adaptation.py` (the battery driver) imports the `eras` battery
from `tools/eras/verify_eras_adversarial.py` in this repository and runs
end-to-end from the clone. One verification leg needs a piece not shipped
here: `selftest_closure_interface.py` needs the pinned upstream build
referenced by the `POSQ_R3_ROOT` environment variable (sha256 pins in
section 2 of `CLOSURE_INTERFACE_MANIFEST.sha256`).

Likewise `derive_posq_sentences.py` and the benchmark entry point of
`bench_rung1_gauss.py` read result files from the original production runs,
which are not included; the rule builder itself is importable and is what
`posq.py` uses.

CREDIT: the kernel is written directly against FLINT/Arb (W. Hart, F. Johansson, A.
Ahlbäck and the FLINT developers; Arb: Johansson 2017, IEEE Trans. Comput. 66:1281);
certified Gauss–Legendre/Gauss–Jacobi rules are obtained by interval Newton in ball
arithmetic, in the manner of Johansson & Mezzarobba (2018, SIAM J. Sci. Comput. 40:C726)
for the Gauss–Legendre case, seeded from SciPy's double-precision Gauss–Jacobi nodes
(scipy.special.roots_jacobi); the likelihood integrated is Felsenstein's (1981, J. Mol.
Evol. 17:368) CTMC tree likelihood; the production counts derive from DravLex v1.0
(Kolipakam, Jordan, Dunn, Greenhill, Bouckaert, Gray & Verkerk 2018, R. Soc. Open Sci.
5:171504; CC BY 4.0; not shipped here).

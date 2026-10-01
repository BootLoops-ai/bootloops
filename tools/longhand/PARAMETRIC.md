# Longhand route A — the arbitrary-precision Feynman-parametric evaluator

Route A of Longhand (manual; the package guide is `GUIDE.md`, route B's manual is
`DISTEVAL.md`). The flat modules `hiprec.py`, `feynman.py`, `parametric.py`,
`bench_box.py`, `bench_se3l.py` in this directory; importable as `longhand` with the
repository's `tools/` on `sys.path`, or by their bare names with this directory on it.

A **second, independent numerical ground-truth route** for finite (eps^0)
scalar Feynman integrals, built to break the ~8–9 digit double-precision
ceiling of `pySecDec` disteval and FIESTA double mode. Use it when you need an independent high-precision value for a finite integral where AMFlow is hard or its eps→0 extrapolation is precision-limited; it works in the deep-Euclidean region.

It evaluates the parametric representation directly in arbitrary precision:

```
I = Γ(N − L·d/2) · ∫_simplex dx δ(1−Σx) U^{N−(L+1)d/2} / F^{N−L·d/2}
```

at the eps^0 order, on a **smooth, positive Euclidean integrand** (s,t<0,
m²>0 — no contour deformation, no threshold crossing). U, F are obtained from
`pySecDec`'s `LoopIntegralFromPropagators` (same normalization as disteval, so
values compare digit-for-digit).

## Engines

| engine | method | best for | measured reach |
|--------|--------|----------|-------|
| `integrate_gauss_product` | FIXED tensored **Gauss-Legendre** product, arbitrary precision, matched per-axis tail/limit | effective dim ≤ 3, smooth & analytic | **40+ d** on the box (spectral; cross-validated to 46+ d, see below) |
| `integrate_tanhsinh` | nested adaptive **tanh-sinh** (`mpmath.quad`) | dim ≤ 2 quick checks | ~17 d (nested adaptive error floor; ~25 d on the box's 2-D reduction) |
| `integrate_qmc` | **CBC rank-1 lattice** + Korobov periodization + M random shifts (parallel, gmpy2) | irreducible dim 4–7 | ~3 d on the 6-D 3-loop self-energy benchmark, with a statistical (median-of-shifts) error estimate |

### Dimension reduction (key step)
Before numerics, Feynman parameters linear in the denominator are integrated
out in closed form, iterated through two stages
(`feynman.reduce_linear_full`): Cheng-Wu chart `x_last=1`, then stage 1
`∫₀^∞ dx/(A x + B)² = 1/(A B)`, then stage 2 by partial fractions on the two
surviving power-1 factors,
`∫₀^∞ dx/((A₁x+A₀)(B₁x+B₀)) = log(A₁B₀/(A₀B₁))/(A₁B₀−A₀B₁)` —
valid because all coefficient polynomials are >0 in the Euclidean region.
Each stage removes one numerical dimension *exactly*. This is what lets the
one-loop box collapse to 2-D. Higher-loop topologies
whose F is quadratic in the remaining variables (e.g. the 3-loop self-energy
benchmark: 6-D after one elimination) are **irreducible** for the CAS and are
handed to the QMC engine.

### Precision contract (inputs)
Kinematic inputs must be **exact** (int, `Fraction`, decimal string). A value
like `mpf(-1)/3` built at a low ambient `mp.mp.dps` is a *different rational
number* at that precision's level; every engine then converges spectrally to
the integral of that perturbed input, silently capping agreement with the true
value near the input's own precision. `bench_box` converts all inputs at
working precision.

## Benchmarks (measured)

**1-loop massive box** (`bench_box.py`, deep-Euclidean s=−1,t=−1/3,m²=1):
two independent reductions ship — `box` (Cheng-Wu x₃=1, loop parameter
integrated out) and `box_crosschart` (Cheng-Wu x₀=1, a mass parameter
integrated out) — sharing no algebra beyond U and F. At `dps=45`/`dps=60`
they **agree to 46–47 d** on all three pinned references (M5=1,2,5), and a
nested adaptive tanh-sinh run agrees to its own ~25 d floor. M5=2 value:
`0.1245570572327092697104024441545542480922` (40 d, quotable in full).
The loop-integration kernel evaluates the `D=4PR−Q²≈0` neighborhood by a
series around the double-root point (the closed forms cancel catastrophically
there), so the integrand is uniformly accurate across its whole domain.
Before quoting digits beyond the pinned 40, run the two-dps convergence
check (e.g. `dps` and `dps+15` runs must agree to the quoted depth).

**3-loop self-energy** (`bench_se3l.py`, full 3-loop 8-propagator scalar
light-by-light self-energy, s=−1,t=−1/3,m²=1):
finite at eps^0, `I = ∫_simplex 1/F²`. After one analytic elimination
(F linear in x6/x7) it is an **irreducible 6-D** integral. The QMC engine gives
an **independent, error-bounded** value confirming the disteval target
`I ≈ 0.52495777354` to ~3 d (measured: `0.522771 ± 0.0020` at N=120011,
48 shifts, 21 s). **The CBC-QMC route does NOT reach 30 d** on this 6-D
integrand at feasible compute — the rank-1 lattice error decays only
polynomially in 6 dimensions, so it confirms disteval at comparable precision
rather than beating it. **A ≥30 d value for this integrand must come from the
nested-loop / dispersive analytic reduction** (the inner 2-loop self-energy
Σ(w) to high precision, then a low-dim outer integral via the Gauss-product
engine). This tool supplies the general infrastructure, the 40-d box oracle
with its independent cross-chart validation, and the independent 3-d
confirmation on this benchmark.

## Usage

```python
# with the repository's tools/ directory on sys.path:
from fractions import Fraction
from longhand import box, box_crosschart, integrate_qmc, integrate_tanhsinh, integrate_spec
v = box(-1, Fraction(-1, 3), 1, 2, dps=35)   # low-dim, Gauss product (exact inputs -- see the precision contract)

# general builder + QMC (sympy; pySecDec for build_UF):
from longhand.feynman import build_UF, reduce_linear_full, build_integrand_expr

# the same modules by their bare names, with tools/longhand itself on sys.path
# (how the benches and the selftests import them):
from bench_box import box, box_crosschart
from feynman import build_UF, reduce_linear_full, build_integrand_expr
from hiprec import integrate_qmc, integrate_tanhsinh

# a UF-spec JSON through the QMC engine (schema in GUIDE.md; toy_qmc_spec.json is a worked example):
res = integrate_spec('toy_qmc_spec.json', N=8009, n_shifts=8, korobov_p=3, dps=25, nproc=8)
print(res['value'], res['error'], res['digits'])
```

```bash
# box benchmark (run inside tools/longhand)
python3 bench_box.py 35

# 3-loop self-energy QMC cross-check (parallel; needs pySecDec + sympy for the one-time reduction, cached in se3l_reduced.json beside the script)
python3 bench_se3l.py --N 64007 --shifts 24 --p 3 --dps 40 --nproc 24 --out result.json
```

## Files
- `hiprec.py`    — engines: CBC lattice construction, Korobov periodizer, parallel median-shift QMC, fixed Gauss-Legendre product, nested tanh-sinh.
- `feynman.py`   — U/F extraction (pySecDec) + iterated closed-form linear reduction (`reduce_linear_full`) + general-(a,b) integrand assembly.
- `parametric.py` — the UF-spec front end: `load_spec`, `check_spec` (homogeneity/degree checks), `make_eval` (sparse power-table evaluator), `spec_factory`, `integrate_spec`; `box_spec` is the one-loop box as a spec (positive control).
- `bench_box.py` — 1-loop box benchmark (Gauss product, 40+ d) plus the independent cross-chart reduction `box_crosschart` used to validate it.
- `bench_se3l.py` — 3-loop self-energy benchmark (QMC, independent cross-check).
- `cbc_cache.json` — the vendored CBC generating vectors the benchmarks and the battery use (read-only; new vectors go to the user cache `$LONGHAND_CBC_CACHE`, older name `$HIPREC_CBC_CACHE`, default `~/.cache/longhand/cbc_cache.json`).
- `toy_qmc_spec.json` — a synthetic 6-D UF-spec fixture with a hand-provable Dirichlet closed form (exact value 13/1080; derivation in its `note` field).
- `selftest.py` (the package battery, both routes) and `selftest_parametric.py` (the generic pathway's pool control).

## Scope / honest limits
- **Region:** deep-Euclidean only (F>0 on the simplex). Physical/threshold
  regions need contour deformation — not implemented (and not where an
  independent ≥30 d reference is the bottleneck).
- **Order:** eps^0, finite integrals. The closed-form reduction runs for
  a=N−(L+1)d/2 = 0, b=N−Ld/2 = 2 at d=4; for other (a,b) the full
  (n−1)-dim Cheng-Wu integrand U^a/F^b is emitted unreduced (finiteness of
  the target order stays the caller's responsibility).
- **Honest precision reach.** For the box (effective dim ≤ 3 after
  reduction) the Gauss-product engine is spectral: the two independent
  reductions agree to 46–47 d at dps=45/60, and the pinned 40-d references
  are quotable in full. Deeper digits need only a larger `dps` plus the
  two-dps convergence check. Exact inputs are part of the contract (see the
  precision contract above).
  For irreducible dim ≥ 4 (the 3-loop self-energy) the QMC engine gives an
  *independent, error-bounded* value at only a few digits — a genuine second
  route / sanity check, **not** a 30-d certificate by itself.
- A ≥30 d value for such an integrand is therefore **not** delivered by this
  tool's QMC path; it needs the nested-loop/dispersive reduction. This tool is the reusable
  infrastructure (U/F extraction, linear & quadratic closed-form parameter
  elimination, the arbitrary-precision Gauss-product and parallel CBC-QMC
  engines) plus the two honest benchmarks.

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).

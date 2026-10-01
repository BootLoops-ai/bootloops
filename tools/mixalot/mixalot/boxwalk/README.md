# mixalot.boxwalk — exact box-moment contiguity walker (mixalot member)

Home: `tools/mixalot/mixalot/boxwalk/` — a member of the mixalot package
(this subpackage is the only copy; mixalot's `GUIDE.md`/`MANUAL.md` carry
its card, mixalot battery leg M16 runs its selftest).

Exact rational evaluation of

    Z(u) = ∫_{[0,1]^n} ∏_k P_k(x)^{u_k} dx,   P_k ∈ Z[x_1..x_n],  u ∈ N^m

for any finite set of integer polynomials and any target exponent vector,
via a mod-p contiguity walk on the exponent lattice + CRT + Wang rational
reconstruction. "Kira for box-moment contiguity": Singular finds
boundary-safe relation modules (fields a with a_e·x_e(1−x_e)∂_e tangent to
the active divisors), flint certifies them exactly over Z, a recorded pivot
program refills the moving moment window, per-prime replay + CRT + rational
reconstruction produce the exact value, gated.

**Status: validated on the shipped selftest battery (all legs PASS).**

## When to use

- You need Z(u) **exactly** (a Fraction, not floats) and u is moderate
  (window memory `2·(1+Σ_k u_k·deg P_k)^n·8` bytes per prime must fit).
- You need the moment window M(β;u) mod p on a box, exactly.
- You need the **scalar contiguity recurrence** of a dataset along one
  exponent direction (`emit_fiber_recurrence`) — integer coefficients,
  2-prime-stable order, held-out-prime verified — plus a shift-operator
  calculus (`fiber.gcrd_reduce` / `gcrd` / `lclm` / `right_divides` /
  `op_mul` / `annihilates_mod`) for reducing and stitching the lifted
  operators across rays, exact over Q(n).
- You need to **re-check a recorded manifest later** (`verify_manifest`,
  `python3 -m mixalot.boxwalk verify`): the value is replayed mod fresh primes through the
  dense oracle — never the recorded program — the CRT + rational
  reconstruction is re-derived locally from the recorded residues, and the
  spec/program hashes are re-checked.

## The math in one page

For any polynomial field a = (a_1..a_n), F_e := a_e·X_e with
X_e = x_e−x_e², ∫∂_e(F_e·x^γ·∏P^u)dx = 0 (X_e kills the boundary). With the
exact divisions G_j := Σ_e a_e X_e ∂_e P_j = q_j·P_j + r_j (flint divmod):

    0 = Σ_δ [D(δ) + Σ_e γ_e F_e(δ−e_e) + Σ_j u_j q_j(δ)]·M(γ+δ; u)
        + Σ_{j: r_j≠0} u_j Σ_δ r_j(δ)·M(γ+δ; u−e_j),      D = Σ_e ∂_e F_e.

Tangency to divisor j ⇔ r_j = 0 ⇔ no down-shift into level u−e_j. The walk
raises one class at a time (M'(β) = Σ_γ P_k[γ]M(β+γ), consuming deg_e(P_k)
window layers per axis — **slab**), or holds the window constant by
refilling the rim with a recorded triangular pivot program of certified
relations (**refill**) where the relation module provably closes the rim.
Coefficients are affine in (γ,u): the program is recorded once (reference
prime, segment-start u) and replayed at every (u, prime). Z = M(0) at the
final u; per-prime values are CRT-combined and rationally reconstructed.

## API

```python
import sys; sys.path.insert(0, '/path/to/repo/tools/mixalot')  # the directory containing mixalot/
from mixalot import boxwalk
spec = boxwalk.load_spec(boxwalk.home() + '/examples/jc_quartet.json')   # or a dict
pl   = boxwalk.plan(spec)                  # order, mechanisms, programs
res  = boxwalk.walk(spec, pl, primes)      # {'Z': [mod p], 'u', 'valid'}
man  = boxwalk.produce(spec, plan=pl, nprimes=16, batch=4, procs=8,
                       outdir='runs/...')  # exact Fraction + manifest
rec  = boxwalk.emit_fiber_recurrence(spec, ray_class='xxxx', depth=80)
rep  = boxwalk.verify_manifest(spec, 'runs/.../MANIFEST.json', nfresh=2)
```

Spec JSON: `variables` (list), `polynomials` {label: {"e1 e2 .. en": coeff}},
`target_u` {label: count}, `order` ('auto' = ascending counts, largest class
last as the pure final ray) and `options`.

Gate battery in `produce()` (all must pass; results in `MANIFEST.json`):
- **G1 crossover** — walk == dense oracle at truncated u (value + window).
- **G2 mutation** — corrupted program/polynomial detected loudly.
- **G3** — two disjoint CRT prime sets reconstruct byte-identical Fractions.

Recorded manifests are re-checkable by an independent replay code path:
`verify_manifest(spec, manifest)` / `python3 -m mixalot.boxwalk verify
SPEC.json MANIFEST.json [--nfresh N] [--replan] [--max-cells C]` runs V1 (spec hash),
V2 (recorded Z reduces to every recorded residue), V3 (local CRT + Wang ratrec
over the recorded residues reproduces the recorded Z), V4 (dense-oracle replay
mod fresh primes disjoint from the recorded set — cost-estimated first,
refuses loudly over `--max-cells`), V5 (recorded gates + target
consistency), and V6 (program hash, when a plan is supplied / `--replan`).

Prime-batched walk: K primes share one recorded program and one numpy pass
(windows stacked on a leading prime axis). Farm: multiprocessing pool over
prime batches; per-step time-stamped logs (rate-fittable). Launch long runs
under your machine's job launcher (see the cli.py header).

## CLI

From `tools/mixalot/` (or anywhere with `tools/mixalot` on `sys.path`):

    python3 -m mixalot.boxwalk selftest
    python3 -m mixalot.boxwalk plan    SPEC.json [--slab-only]
    python3 -m mixalot.boxwalk produce SPEC.json OUTDIR [--nprimes N] [--batch K] [--procs P]
    python3 -m mixalot.boxwalk fiber   SPEC.json [--ray CLASS] [--depth D]
    python3 -m mixalot.boxwalk verify  SPEC.json MANIFEST.json [--nfresh N] [--replan] [--max-cells C]

`python3 tools/mixalot/mixalot/boxwalk/cli.py ...` (or `selftest.py`) also runs
as a plain script from any cwd. The selftest writes its record to
`SELFTEST.json` beside the package unless `BOXWALK_SELFTEST_OUT=/path/file.json`
redirects it (the mixalot battery routes it to a scratch dir).

Dependencies: numpy, sympy, python-flint; Singular on PATH for the relation
search (the selftest's Singular leg skips by name without it). The fiber
engine's recurrence fitter is mixalot's own sha-pinned `annihilator` engine
(`mixalot.engines.load('annihilator')`); `BOXWALK_ANNIHILATOR=/path/to/annihilator.py`
points it at another copy (e.g. the standalone `tools/annihilator/annihilator.py`).

## Caveats (read before production)

1. **n ≤ 6 practicality.** Window cells scale as side^n; the memory formula
   is `2·side^n·8·K` bytes per batch job, side = 1 + Σ_slab steps·deg.
   `plan()` prints it; keep ≤ 4 GB (target ≤ 1 GB).
2. **Slab is the default; refill is rare and gated.** Measured: multi-class tangency modules do NOT close the rim, and apparent
   closures can ride u-dependent coefficient cancellations that are not
   replay-safe (replayed at another prime they crash out of bounds). boxwalk's peel exempts an offset only if
   its exact-integer coefficient is identically zero along the segment
   (replay-safe for every prime and step); refill is used only when that
   peel closes, else the planner demotes to slab. Long single-class tails
   (N ≫ 40) therefore need the fiber-recurrence lane, not the window.
3. **Singular timeout → sampled-scan fallback** (degree-truncated,
   lightly tested); all generators are flint-certified exactly either way.
4. **Degeneracy patch:** if a recorded pivot vanishes mod p at some (u,p),
   the replay raises DegeneracyError and produce() drops that prime batch
   (logged). Persistent vanishing across primes ⇒ genuine degeneracy:
   re-record the program at a different segment-start u.
5. **Exact-coefficient int64 bound:** families whose |coeffs|·(side+Σu)
   approach 2^61 are excluded from the peel (bound_ok); replay asserts it.
6. **Fiber engine is the dense oracle** (cost grows with depth); fitting
   production-scale tails from walk-generated fibers is deferred.
7. Big targets need many primes: G3/ratrec fail loudly ("increase
   nprimes") rather than return a wrong value.

## Limitations

- No symmetry reduction: for symmetric kernel sets (e.g. quartet
  automorphisms permuting β-axes) the walk does not exploit
  symmetric-subspace windows.
- No walk-generated fiber fitting for long tails, and no window transport
  by a fitted recurrence.

## Files

core.py (spec, oracle, window primitives) · relations.py (relation
calculus, peel, prime-batched refill) · syz.py (Singular driver + certify +
fallback) · planner.py (segments, programs, memory) · driver.py (walk,
produce, gates, manifest) · fiber.py (recurrence emission/extension +
shift-operator GCRD/LCLM calculus) · verify.py (independent manifest
replay) · cli.py / `__main__.py` (`python3 -m mixalot.boxwalk`) · selftest.py
(must PASS; reference record SELFTEST.json) · examples/jc_quartet.json (15 JC quartet
kernels; regenerate with examples/make_jc_quartet.py, which reads the
kernels from examples/jc_kit/, two modules of the jackandjill phylogenetics package).

CREDIT: relation modules by Singular (Decker, Greuel, Pfister, Schönemann and the
Singular team), exact certification by FLINT [Arb] via python-flint; the walk is
Laporta-style elimination [Lap] on contiguity relations with CRT and Wang/Monagan
rational reconstruction [Wang81, Monagan]. Bracketed keys resolve in REFERENCES.md at
the repository root.

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).

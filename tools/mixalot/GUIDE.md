# mixalot — guide

Tool page: https://bootloops.ai/tools/mixalot.html
Full detail: `MANUAL.md` in this directory (this guide is the router card).

KIND: package (vendored sha-pinned snapshot; upstream source trees untouched)

REQUIREMENTS: `scipy` and `dynesty` (both PyPI) beyond the common core — the vendored
MC-estimator member imports scipy at load, and battery leg M11b runs its dynesty
nested-sampling route. `pip install scipy dynesty`.

PURPOSE: exact Bayesian evidence for mixture models, in a box. One front door over the
validated mixture-evidence engines: exact evidence Z(U,g) + component-count posteriors
(`mixalot.Z/gstar/bayes_factor`), 1-var closed forms + large-N series
(`closed_form_1var`, `zseries`), large-scale 2-way tables via mod-p+CRT
(`Z_table_modp`/`lsx55_exact`), large-k (`bigk2`), generic/large-g (`biggen`), the
exact Dirichlet-process limit (`w4_blind_gf.Zdpm_gf`), the contiguous-process census
engine with seam posteriors (`seg_v1` + blind twin) and planted-seam power probes
(`mic_power`), and frozen-component blend/presence evidence (`frozen_comp_v1`).
`mixalot.verify()` is fail-closed on vendor integrity, advisory-loud on upstream-source
drift.

MEMBERS: `vendor/eco/etienne_evaluate.py` — Etienne (2005) neutral-biodiversity
exact/certified likelihood (battery M14); `vendor/eco/pilot/` — the Etienne pilot
engine cluster behind it (16 files: product-tree DP `ball_engine`, interval-Newton
certified ML `phase2_certify`, equal-I multisample factorization, Kronecker-packed
hierarchical evaluators, independent urn oracle; 12 loadable engines — battery M19
runs `gate_engines`' G1/G2 gate at load — plus 4 pinned-but-not-loadable scripts;
see MANUAL.md's member section); `vendor/scrna/telegraph_evaluate.py`
(battery M15); `mixalot/boxwalk/` — the exact box-moment contiguity route (the fourth
evidence route; section below; battery M16 runs its 14-leg selftest).

## Member — boxwalk (`mixalot/boxwalk/`, subpackage)

PURPOSE: exact rational evaluation of Z(u) = ∫_{[0,1]^n} ∏_k P_k(x)^{u_k} dx for
integer polynomials P_k (n ≤ 6 in practice) and integer exponent vectors u — the same
evidence-integral class as the allocation sum, the lattice DP and the closed forms, by
a route with a different scaling profile: a mod-p contiguity walk on the exponent
lattice, then CRT and rational reconstruction. Singular finds boundary-safe relation
modules, python-flint certifies them exactly over Z, a recorded pivot program refills
the moving moment window, per-prime replay closes. `produce()` is gated (walk == dense
oracle at a truncated u; a corrupted program/polynomial is detected; two disjoint prime
sets reconstruct the identical fraction) and writes a `MANIFEST.json` that
`verify_manifest` re-checks later by an independent dense-oracle replay at fresh
primes. Also emits the exact scalar contiguity recurrence of a dataset along one
exponent direction (`emit_fiber_recurrence`, fitted with mixalot's pinned
`annihilator` engine, held-out-prime verified) with a shift-operator GCRD/LCLM
calculus for reducing and stitching the lifted operators.

INVOKE: `from mixalot import boxwalk` → `boxwalk.load_spec(path_or_dict)`,
`boxwalk.plan(spec)`, `boxwalk.produce(spec, nprimes=16, batch=4, procs=8,
outdir=...)` (exact `Fraction` + manifest), `boxwalk.emit_fiber_recurrence(spec, ...)`,
`boxwalk.verify_manifest(spec, manifest, nfresh=2)`. CLI from `tools/mixalot/`:
`python3 -m mixalot.boxwalk {selftest | plan SPEC.json | produce SPEC.json OUTDIR
[--nprimes N --batch K --procs P] | fiber SPEC.json [--ray CLASS --depth D] | verify
SPEC.json MANIFEST.json [--nfresh N --replan --max-cells C]}`. Example spec:
`mixalot/boxwalk/examples/jc_quartet.json` (15 Jukes–Cantor quartet kernels;
`examples/make_jc_quartet.py` regenerates it from the two kit modules in
`examples/jc_kit/`). Detail: `mixalot/boxwalk/README.md`.

DEPS: numpy, sympy, python-flint; Singular on PATH (the selftest's Singular leg skips
by name without it). `BOXWALK_ANNIHILATOR=/path/annihilator.py` swaps the recurrence
fitter; `BOXWALK_SELFTEST_OUT=/path/file.json` routes the selftest record off the tree.

LIMITS: window cells scale as side^n (memory `2·side^n·8·K` bytes per batch; `plan`
prints it — keep ≤ 4 GB); slab is the default mechanism and refill is used only where
the certified relation module provably closes the rim; long single-class tails
(N ≫ 40) belong to the fiber-recurrence lane, not the window; big targets need many
primes and the reconstruction fails loudly ("increase nprimes") rather than return a
wrong value; no symmetry reduction. Selftest: 14 legs (sympy-exact closed form, gates,
JC-quartet toy walk vs oracle at every step and both primes, mutation controls, fiber
fit/lift/verify + GCRD/LCLM calculus, Singular smoke, manifest verifier pass +
four corruptions that must fail loudly); reference record `mixalot/boxwalk/SELFTEST.json`.

BATTERY (the reason to trust it): `python3 battery/battery.py` reproduces the recorded
reference receipts through the vendored bytes — flagship exact pairs, the frozen primary
census cell (g*=3, log10BF 95.998, exact-vs-float 3.6e-12), blend gate rationals,
DP-limit rationals incl. the -1/48 g-scaling law, scaling row, planted-seam power arms
— plus tamper + comparator mutation controls that MUST fail (M2, M10). Record:
`battery/BATTERY.txt`. `--full` adds the heavy legs (large-table flagships, power
arms; minutes). Any FAIL ⇒ exit 1. NOTE: some receipt-comparison legs read
reference data not shipped in this repo (roots named via MIXALOT_REFDATA_SWEEP / MIXALOT_REFDATA_PILOT (+ MIXALOT_PILOT_DIR for the --full power arms), legs SKIP/refuse by name when unset) — on a clone without it, trust transfers
through the in-package record + the vendored-bytes pin checks (M1/M2), which are fully
self-contained.

FENCES: Dir(1) convention; iid category counts are weakly identified (plant at your
shape first — MANUAL.md's worked example); swap_route float route is
k=2/certified-at-N=100; Z_table_modp primes < 2^25; asymmetric ratrec bounds for
evidence fractions.

PROVENANCE: `mixalot/_pins.py` records the source path + sha for every
vendored file (the pin record; those paths are not in this repo). The vendored bytes are the shipped
tool; the original source trees remain the record of origin. Do not edit vendored files — the
battery's M1/M2 legs fail closed on any drift.

INVOKE: `import mixalot` (front door: `mixalot.Z`, `mixalot.gstar`,
`mixalot.bayes_factor`, `mixalot.verify()`); `python3 -m mixalot.cli ...`;
`python3 -m mixalot.boxwalk ...` (the box-moment member's CLI);
`python3 battery/battery.py [--full]`; worked examples in `examples/`.

CREDIT: the exact mixture-evidence integrals, their conventions and the worked examples
(swiss, coin, Ex. 5.5) are those of S. Lin, B. Sturmfels & Z. Xu, 'Marginal likelihood
integrals for mixtures of independence models', JMLR 10 (2009) 1611-1631,
arXiv:0805.3602 — mixalot exists to make their integrals routine, and every LSX row is
gated against their printed values (Ex. 5.5 is the 3x3 table of Evans, Gilula & Guttman
1989, Biometrika 76:557); the Dirichlet-process limit is Ferguson (1973, Ann. Stat.
1:209) / Antoniak (1974, Ann. Stat. 2:1152). The ecology member evaluates R. S.
Etienne's sampling formula (2005, Ecol. Lett. 8:253; multisample form 2007, Ecol. Lett.
10:608) and the random-fission model of Etienne & Haegeman (2011, Theor. Ecol. 4:87),
gated against Etienne's own worked example as carried by R. K. S. Hankin's untb package
(2007, J. Stat. Softw. 22(12)) and against the six-forest fits on the abundance data of
Volkov, Banavar, He, Hubbell & Maritan (2005, Nature 438:658); the urn oracle is Hoppe's
(1984, J. Math. Biol. 20:91). The scRNA member evaluates the telegraph model of Peccoud
& Ycart (1995, Theor. Popul. Biol. 48:222) in the Beta–Poisson form, with parameter
conventions pinned to txburst (Larsson et al. 2019, Nature 565:251). The DCM burstiness
arm follows Madsen, Kauchak & Elkan (2005, ICML); the Federalist fixtures use the
function-word lists of Mosteller & Wallace (1963, JASA 58:275; 1964). The Monte-Carlo
comparison suite implements the estimators of Newton & Raftery (1994), Meng & Wong
(1996), Xie et al. (2011), Chib (1995; with Neal's 1999 label-switching caveat) and
nested sampling (Skilling 2006) via dynesty (Speagle 2020, MNRAS 493:3132), on JAX;
certified ML uses interval Newton (Moore 1977; Krawczyk 1969; Rump 2010). Ball
arithmetic: Arb (Johansson 2017) via python-flint.

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).

# Eichler.jl

Certified elliptic/modular numerics for the elliptic Feynman-integral bootstrap.
Julia ≥ 1.11 + Arblib/Nemo only. **Every public function returns honest ball enclosures** —
all series truncations, quadrature errors and ODE-transport tails are bounded by proved
inequalities evaluated in Mag/Arb arithmetic. Conventions: `docs/conventions.md`.
Results and validation: `NUMERICS-RESULTS.md` (all five gauntlet items pass).

## Layers

1. **Γ₁(6) modular layer** (`Gamma16`): Hauptmodul t(τ) = 9η(6τ)⁸η(τ)⁴/(η(2τ)⁸η(3τ)⁴)
   and certified interval-Newton inverse (`nome_from_t`, `nome_and_log_from_t`); Sebbar
   generators e₁, e₂; kernel forms ψ₁/π, f₁…f₄, g_{2,0}, g_{2,1}, g_{2,9}, g_{3,1} as
   certified q-series and pointwise eta quotients; iterated Eichler integrals
   `iterated_integral_value` to arbitrary depth (≥11 exercised) with rigorous truncation
   and log-stack remainders.
2. **Sunrise** (`sunrise2`, `sunrise4`, `sunrise_boundary_constants`): the all-orders
   AW representation in d=2−2ε and the exact equal-mass dimension shift to d=4−2ε;
   boundary constants from a certified ε-expanded ₂F₁ ODE transport. 515+ certified
   digits per ε-order through ε⁷ in ~100 s.
3. **Curve layer** (`QuarticCurve`, `period_quadrature`, `curve_periods`,
   `incomplete_quadrature`, `incomplete_edge_quadrature`, `acm_periods`,
   `picard_fuchs_residual_acm`, `frobenius_transport_acm`): certified roots, periods by
   branch-certified contour quadrature and by K/AGM closed forms (cross-validated),
   quasi-period derivative chains, and a generic validated Taylor/Frobenius transport
   engine (`ODEOperator`, `transport`) for Picard–Fuchs operators with ε-polynomial
   coefficients.
4. **Puncture layer** (`abel_image`, `abel_image_deformed`, `abel_image_quartic`,
   `varpi0`, `varpi0_prime`, `third_kind_G`): Abel images of marked points (closed-form
   ACKM route and general-quartic route, consistent up to certified 2-torsion); the
   third-kind period G(s,t,m²) defined by the exact Zenodo-2502.00118 differential
   relation anchored at G(s,0) = 0, validated against the printed series and the
   ancillary threshold constants a₀, a₁.
5. **Dictionaries** (`cusp_dictionary`, `interior_dictionary`, `dirichlet_L_chim3`,
   `write_dictionary`): PSLQ-ready certified vectors — cusp ring (Dirichlet L of
   conductor | 6, sixth-root-of-unity polylogs, π, logs) and interior period rings
   (ψ₀, ψ₁, derivatives, G, Z_{x_p}).
6. **Transport beyond the sunrise** (`vop_transport.jl`, `cy_transport.jl`,
   `pm_periods.jl`, `siegel.jl`). `EllipticSector` / `eichler_transport` /
   `vop_assemble` package the variation-of-parameters solution of an elliptic 2×2
   sector, g = pref·(c₁ + c₂τ + ∫S dτ), as a single q-series call on the cusp-connected
   branch, with `vop_fit_constants` fixing (c₁, c₂) from anchor values and
   `MultiEllipticSector` covering several curves at once. `PFOperator` /
   `mum_frobenius_basis` / `period_vector` / `cy_transport` extend this up the
   Calabi–Yau ladder: for a Picard–Fuchs operator with a point of maximal unipotent
   monodromy they build the log-tower Frobenius basis in exact rational arithmetic,
   sum it inside the MUM disk, and continue the period vector (and the inhomogeneous
   solution) along a user-supplied path with the certified Taylor transport of layer 3;
   `banana_pf(l)` supplies the l-loop equal-mass banana operators and `mirror_map` the
   mirror coordinate. `pm_periods.jl` encodes the four K3/CY3 operators that enter
   fourth- and fifth-order post-Minkowskian black-hole scattering (arXiv:2401.07899:
   Legendre Sym², ₄F₃, Apéry, Hadamard χ=80) in θ-form, with `pf_from_theta`,
   closed-form holomorphic periods where known (`holo_period_closed`, `apery_numbers`)
   and `frobenius_basis` / `pf_residual` checks to 80+ digits. `siegel.jl` is the
   genus-2 layer: certified hyperelliptic period matrices (`HyperellipticCurve`,
   `hyper_big_period_matrix`, chain assembly for complex branch points), genus-2 theta
   functions and Igusa / Rosenhain / Thomae invariants through FLINT's `acb_theta`,
   the Abel–Jacobi map for real-ordered curves, and the Gauss–Manin (Rauch) variation
   ∂τ/∂eₘ; the Siegel-modular analogue of `eichler_transport` is not implemented and
   `siegel_transport` raises `SiegelNotImplemented` with the reference list.

## Branch conventions that matter

- Nome at kinematic point t: the branch continuously connected to the cusp q(0)=0;
  for Euclidean t<0 the log is L = log|q| + iπ (Feynman t+i0), set exactly.
- Marked points landing on branch cuts: evaluated with an exact dyadic Feynman
  deformation t → t + i·2^{-k} (documented in the value), or via the certified
  +i0 edge-branch quadrature (`incomplete_edge_quadrature`).
- Cross-check oracles: GiNaC 1.8.7 (`ginac-oracle/`, external only) and the Zenodo
  ancillaries of arXiv:2502.00118 (`vendor/zenodo-2502.00118/`, fetched on demand —
  see the README there).

## Scope notes

- Period evaluation at arbitrary interior kinematic points of the ACKM/BCNTW family is
  exercised end-to-end (including the s=4m² threshold ring); the generic transport
  engine accepts any Fuchsian operator — deriving the Picard–Fuchs operator of the
  deformed Zγ family in (s; M) is separate geometry work (Oscar/Griffiths–Dwork),
  not needed for direct period evaluation, which works now via certified quadrature.

```julia
using Eichler, Arblib
G = Gamma16(768, 700)                       # Γ₁(6) layer at 768 bits, q-order 700
S4 = sunrise4(G, Acb(-1; prec=768); J=10)   # sunrise, d=4-2ε, ε⁻²..ε⁷, certified
```

Requires Julia 1.11 or newer. The committed `Manifest.toml` (generated on Julia 1.11)
pins the tested dependency tree; `Pkg.instantiate()` reproduces it exactly. On a
newer Julia minor version it still instantiates (Pkg then recommends a
`Pkg.resolve()` to re-pin the manifest for that version).
Tests: `julia --project=. -e 'using Pkg; Pkg.test()'`.

Family members (repo-side, outside this Julia package): `tools/eichler/` —
`theta9/` (certified Siegel theta-constant observables on this package's
acb_theta layer) and `mirror6/` (exact mirror-map / instanton-type fingerprint
engine for MUM-type points of D-finite operators).

## Attribution

Eichler.jl implements published mathematics and follows its authors' conventions:
iterated integrals of modular forms and the all-orders sunrise of Luise Adams and Stefan
Weinzierl [AW], building on Adams, Bogner and Weinzierl [ABW] and on Bloch and Vanhove's
elliptic-dilogarithm analysis of the sunset [BV]; the non-planar two-loop family of
Taushif Ahmed, Ekta Chaubey, Mandeep Kaur and Sara Maggio [ACKM]; the diphoton amplitudes
and ancillary data of Matteo Becchetti, Federico Coro, Christoph Nega, Lorenzo Tancredi
and Fabian Wagner [BCNTW] (Zenodo 10.5281/zenodo.14733100, CC BY 4.0, fetched on
demand); banana Calabi–Yau periods after Bönisch, Fischbach, Klemm, Nega and Safari
[BFKNS] and Pögel, Wang and Weinzierl [PWW]; the post-Minkowskian K3/CY3 operators of
Klemm, Nega, Sauer and Plefka [KNSP]; the genus-2 layer after Duhr, Porkert and Stawinski
[DPS] and Marzucca, McLeod, Page, Pögel and Weinzierl [MMPPW]; Calabi–Yau operator tables
of Almkvist, van Enckevort, van Straten and Zudilin [AESZ, vS]; real-analytic Eisenstein
series after Brown [Br17]. Genus-2 theta functions use FLINT's certified `acb_theta`
module by Jean Kieffer [Kief], which implements the quasi-linear algorithm of Noam D.
Elkies and Jean Kieffer [EK]; all ball arithmetic is Fredrik Johansson's Arb [Arb] via
Arblib.jl (Marek Kaluba, Joel Dahne and contributors) and Nemo [Nemo]; the optional
oracle drives GiNaC's elliptic machinery by Moritz Walden and Stefan Weinzierl [WW,
GiNaC]. The code, rigor architecture and certified bounds are ours; the mathematics is
theirs. Bracketed keys resolve in [`REFERENCES.md`](../../REFERENCES.md) at the
repository root.

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz; code written by Claude (Anthropic) under his supervision. Documentation is licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository root; third-party components keep their own licenses (THIRD_PARTY.md).

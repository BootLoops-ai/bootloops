"""
Eichler.jl — certified elliptic/modular numerics for the elliptic Landau bootstrap.

Ball discipline: every public function returns Arb/Acb enclosures with honest radii.
Conventions: docs/conventions.md. Layers:
  1. Γ₁(6) modular layer (Hauptmodul, kernel forms, iterated Eichler integrals)
  2. curve layer (periods/quasi-periods of y² = quartic, AGM + Frobenius continuation)
  3. puncture layer (Abel images, third-kind period G)
  4. constant dictionaries (cusp Dirichlet-L ring, interior period rings)
"""
module Eichler

using Arblib
using Nemo: Nemo, QQFieldElem, ZZRingElem, polynomial_ring, power_series_ring,
    mullow, divexact, numerator, denominator

include("qseries.jl")
include("modular.jl")
include("gamma16.jl")
include("iterated.jl")
include("epsseries.jl")
include("frobenius.jl")
include("hypgeom.jl")
include("sunrise.jl")
include("curve.jl")
include("puncture.jl")
include("dictionaries.jl")
include("vop_transport.jl")
include("cy_transport.jl")
include("pm_periods.jl")
include("siegel.jl")

export SeriesContext, QSeries, evaluate,
    Gamma16, kernel, hauptmodul_point, tau_from_t, nome, nome_from_t,
    iterated_integral, iterated_integral_value, nome_and_log_from_t,
    EpsSeries, eps_coeff,
    sunrise_boundary_constants, sunrise2, sunrise4,
    QuarticCurve, curve_periods, curve_quasiperiods, period_quadrature,
    picard_fuchs_residual_acm, frobenius_transport_acm,
    abel_image, third_kind_G,
    chi_minus, dirichlet_L, dirichlet_L_chim3, clausen,
    cusp_dictionary, interior_dictionary, eisenstein_constants, cusp_reg_tail,
    cusp_constant_E_k_chi, generalized_bernoulli,
    EllipticSector, tau_map, eichler_transport, vop_assemble, vop_fit_constants,
    vop_quadrature, period_polynomial,
    MultiEllipticSector, multi_tau_map, multi_eichler_transport, multi_vop_assemble,
    dlogq_dt_gamma16,
    PFOperator, MUMFrobenius, CYSector, banana_pf, mum_frobenius_basis,
    period_vector, period_matrix, period_matrix_certified, mirror_map,
    cy_transport, frobenius_logseries, theta_form, to_ode, pf_from_theta,
    PMPeriod, pm_k3_legendre, pm_cy3_4F3, pm_k3_apery, pm_cy3_hadamard,
    pm_all_periods, frobenius_basis, frobenius_wronskian, holo_period_closed,
    pf_residual, holo_series_coeffs, apery_numbers,
    HyperellipticCurve, genus, hyper_gap_integral, hyper_big_period_matrix,
    hyper_chain_integral, hyper_default_chain,
    siegel_theta, siegel_theta_all, siegel_char, siegel_char_is_even,
    igusa_invariants, igusa_clebsch, igusa_absolute, rosenhain_from_tau,
    thomae_sextic, abel_jacobi, Genus2Sector, siegel_tau_map, siegel_transport,
    siegel_rauch_de, siegel_gauss_manin, branch_velocities_fd,
    branch_velocities_coeffs,
    SiegelNotImplemented, SIEGEL_G2_EVEN_CHARS

end # module

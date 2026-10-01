"""
LandauAlphabet.jl — bootstrap/regression engine for Feynman integrals.

Pipeline: ansatz construction (weight-graded iterated-integral words over a
pluggable kernel alphabet) → exact constraint assembly over Q (Nemo) →
high-precision sampling (the amflow-cpp CLI; kernel values via the built-in
log-graded series evaluator) → LLL/PSLQ fit with acceptance gates.

One abstraction, two instances:
  * MPL ansätze   — dlog letters over weight-graded G-word basis towers.
  * elliptic ansätze — Γ₁(6) modular kernels per arXiv:1704.08895.
"""
module LandauAlphabet

using Nemo
using Arblib
using JSON3

include("wordutils.jl")    # words, weight grading, shuffle algebra, Lyndon basis
include("itint.jl")        # log-graded series iterated-integral evaluator (Acb)
include("gpl.jl")          # GPL / MPL kernel instance + numerics
include("qseries.jl")      # exact q-series over QQ: eta, Eisenstein E_k(φ,ψ), B_{2,K}
include("modular.jl")      # Γ₁(6) instance: kernels, M_k bases, Landau/cusp constraints
include("symbols.jl")      # integrable symbols over QQ for multivariate dlog alphabets
include("cuspvals.jl")     # EXACT Eisenstein cusp values over Q(√−3) at Γ₁(6) cusps
export cusp_value_B2K, cusp_value_E, cusp_values_level6, cuspvals_field,
       eta_order_at_cusp, cusp_completion, GAMMA16_CUSP_TUPLES, GAMMA16_WIDTHS
include("constraints.jl")  # generic exact linear-constraint assembly + nullspace
include("fit.jl")          # LLL lattice fit, PSLQ, acceptance gates
include("amflow.jl")       # amflow-cpp JSON CLI + staggered-goal consistency protocol
include("cauchy.jl")       # Taylor-coefficient extraction via Cauchy circles

end # module

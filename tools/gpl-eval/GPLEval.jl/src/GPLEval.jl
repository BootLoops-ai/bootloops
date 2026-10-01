"""
GPLEval — fast arbitrary-precision GPL/HPL/eMPL evaluation for the Landau bootstrap.

A thin package layered on:
  * `LandauAlphabet.ItIntEvaluator` / `GPLContext`  — re-exported as the fast,
    suffix-cached SERIES path (heuristic tail);
  * `Arblib.polylog!`                                 — depth-1 Liₘ primitive;
  * Vollinga–Weinzierl Hölder convolution (0410259)   — reaches |z| ≥ min|aᵢ|
    when no letter sits on the segment (the series path's blind spot);
  * VW §5 / FTW 1601.02649 path-concatenation         — native analytic
    continuation when a nonzero letter lies on (0,z); the ±iπ side is fixed by
    `prescription` (default `:ginac` = aₖ+i0, matching GiNaC's 2-arg G);
  * GiNaC `G(lst{a...}, z).evalf()` (via tiny C++ CLI) — independent oracle;
  * `Eichler.iterated_integral` / `EllipticSector`    — eMPL bridge.

Core entry points:

    gpl(a::Vector{<:Number}, z; prec)           -> Acb     # G(a₁,…,aₙ; z)
    hpl(m::Vector{Int}, z; prec)                -> Acb     # H(m⃗; z)
    li(n::Int, z; prec)                         -> Acb     # Liₙ(z)
    evaluate_symbol(word, alphabet, z, point; prec)
    evaluate_basis(words, alphabet, z, point; prec)
    empl(kernels, t, sec::EllipticSector)       -> Acb
    gpl_ginac(a, z; digits)                     -> Complex{BigFloat}

Conventions (Goncharov / Remiddi–Vermaseren, matching the project's
hand-rolled `hpl_fast.py` / `hpl_basis.py` and `LandauAlphabet.gpl`):

    G(a₁,…,aₙ; z) = ∫₀ᶻ dt/(t−a₁) G(a₂,…,aₙ; t),   G(; z) = 1,
    G(0,…,0; z)   = logⁿ(z)/n!,
    G(a; z)       = log(1 − z/a),
    G(0ᵐ⁻¹, a; z) = −Liₘ(z/a),
    H(m⃗; z)       = (−1)ᵖ G(a⃗; z),  p = #(aᵢ=+1),  a⃗ = expand(m⃗).

First letter = OUTERMOST integration. Regularisation of trailing zeros is the
standard tangential-base-point / shuffle prescription (G(…,0;z) via log-shuffle).
"""
module GPLEval

using Arblib
using LandauAlphabet
using LandauAlphabet: ItIntEvaluator, DlogKernel, eval_word, eval_words,
    GPLContext, Word, shuffle_product, const_value, MPLBasisFunction, basis_tower
import Eichler
using Eichler: EllipticSector, QSeries, iterated_integral, evaluate, tau_map

include("gpl_core.jl")
include("hpl.jl")
include("symbol_eval.jl")
include("empl_bridge.jl")
include("ginac_bridge.jl")

# Re-exports from LandauAlphabet (fast cached series path)
export GPLContext, ItIntEvaluator, DlogKernel, Word, shuffle_product,
       const_value, MPLBasisFunction, basis_tower
# Additions in this package
export gpl, gpl_hybrid, hpl, li, hpl_expand, gpl_ginac,
       on_path_indices, analytic_continue,
       evaluate_symbol, evaluate_basis, GPLCache,
       empl, EllipticSector

end # module

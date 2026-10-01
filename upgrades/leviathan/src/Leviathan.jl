# Leviathan.jl — chi-drop Landau-singularity method: a Julia implementation of the
# Euler-characteristic-drop method of Vsevolod Chestnov, Giulio Crisanti and Mathieu Giroux
# (arXiv:2606.29612) and Vsevolod Chestnov and Giulio Crisanti (arXiv:2511.14875).
# No code was copied or translated. The authors' public Mathematica code (landau_codes.m
# in Landau-s-Leviathans; euler_characteristic.m in DiscKosky by Crisanti, Lippstreu,
# McLeod and Polackova) was read as the reference for five details the papers leave to
# the code (PATCHES.md, "What follows the authors' code"). MIT; see LICENSE.
# Stack: Julia 1.10 + Oscar 1.7.3 (see README.md / Project.toml for the environment).

module Leviathan

using Oscar
using Random

const Nemo = Oscar.Nemo

include("primes.jl")
include("parse.jl")
include("staircase.jl")
include("family.jl")
include("chi.jl")
include("elimination.jl")
include("diagnostics.jl")
include("landau.jl")
include("ratrec.jl")
include("reconstruct.jl")

export Family, family, family_UF, sector_G, admissible_sectors, nedges,
       chi_sector, chi_regulated, chi_regulated_constrained, point_on_locus,
       staircase_count, chi_of_ideal, lead_exps,
       sector_candidates, all_candidates, elim_generators, ff_ideal,
       diag_type11, diag_type21, diag_type22, nodegeneracy,
       landau, report, diagnose, LandauResult,
       default_p31, default_p29, random_kinvals,
       landau_reconstruct, sample_sector, reconstruct_sector, minpoly_modp,
       ratrec_multiprime, P31

end # module

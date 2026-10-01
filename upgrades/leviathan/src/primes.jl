# primes.jl — finite-field prime pools + random kinematic points.
# Singular-backed ops (eliminate, reduce) require p < 2^29 (measured DomainError above);
# groebner_basis_f4 (msolve F4 via AlgebraicSolving) takes 31-bit primes (< 2^31).

const P31 = [2147483647, 2147483629, 2147483587, 2147483579, 2147483563,
             2147483549, 2147483543, 2147483497, 2147483489, 2147483477]
const P29 = [536870909, 536870879, 536870869, 536870849, 536870839,
             536870837, 536870819, 536870813, 536870791, 536870779]

default_p31() = P31[2]
default_p29() = P29[1]

"Random nonzero residues mod p for a list of parameter names (reproducible via rng)."
function random_kinvals(rng, names::Vector{String}, p::Int)
    Dict{String,Int}(n => rand(rng, 2:p-2) for n in names)
end

"Distinct small generic rationals for regulators (nu_i, d/2) over QQ.
Small numerators keep the exact function-field lane fast (measured: nf cost doubles
from ~9 s to ~19 s on sunrise-top when numerators grow to ~1e4)."
function generic_regulators(rng, E::Int)
    ps = [103, 107, 109, 113, 127, 131, 137, 139, 149, 151]
    qs = [7, 11, 13, 17, 19, 23, 29, 31, 37, 41]
    nu = [Nemo.QQ(ps[mod1(i, end)] + 2 * rand(rng, 1:13), qs[mod1(i, end)]) for i in 1:E]
    d2 = Nemo.QQ(51 + 2 * rand(rng, 1:7), 4)   # d/2 generic, non-integer
    nu, d2
end

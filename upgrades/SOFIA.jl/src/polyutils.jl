# Polynomial utilities shared across the package.
#
# Upstream reference: scr.m lines 414-447 (MyFactorList, ProportionalPolynomialsQ,
# myPolynomialDecomposition). Where upstream uses probabilistic checks at random
# integer points, we use exact checks — same semantics, deterministic.

"""
    isproportional(p, q) -> Bool

True iff `p == c*q` for a nonzero constant `c`.
Exact replacement for upstream `ProportionalPolynomialsQ` (scr.m:437),
which tests proportionality probabilistically at random integer points.
"""
function isproportional(p, q)
    iszero(p) && return iszero(q)
    iszero(q) && return false
    length(p) == length(q) || return false
    # exponent supports must match in the (shared) ring's monomial order
    parent(p) === parent(q) || return false
    ratio = nothing
    for (cp, cq, ep, eq) in zip(Nemo.coefficients(p), Nemo.coefficients(q),
                                Nemo.exponent_vectors(p), Nemo.exponent_vectors(q))
        ep == eq || return false
        r = cp // cq
        if ratio === nothing
            ratio = r
        elseif r != ratio
            return false
        end
    end
    return true
end

"""
    dedup_proportional(polys) -> Vector

Delete duplicates up to overall constants, keeping first occurrences
(upstream `DeleteDuplicates[..., ProportionalPolynomialsQ]`).
Uses a normalized-key dictionary, so it is O(n) not O(n^2).
"""
function dedup_proportional(polys::AbstractVector)
    seen = Set{Any}()
    out = empty(collect(polys))
    for p in polys
        k = normal_form_key(p)
        if !(k in seen)
            push!(seen, k)
            push!(out, p)
        end
    end
    return out
end

"""Normalize a polynomial to a canonical representative of its constant-multiple
class: divide by the leading coefficient (monic in the ring's monomial order)."""
function normal_form_key(p)
    iszero(p) && return p
    c = Nemo.leading_coefficient(p)
    return Nemo.divexact(p, c)
end

"""
    irreducible_factors(p; positive_exponents_only=true) -> Vector

Non-constant irreducible factors of `p`, each exactly once
(upstream `FactorList[#][[All,1]]` filtered by `!NumericQ`).
Works over QQ and GF(p) mpoly rings via Nemo's `factor`.
"""
function irreducible_factors(p)
    iszero(p) && return [p]
    is_constant(p) && return typeof(p)[]
    f = Nemo.factor(p)
    return [fac for (fac, _) in f]
end

is_constant(p) = Nemo.total_degree(p) <= 0

"""
    factor_list_unique(polys) -> Vector

Upstream `MyFactorList` (scr.m:414): all distinct non-constant irreducible
factors across a list, deduped up to proportionality.
"""
function factor_list_unique(polys::AbstractVector)
    out = eltype(polys)[]
    for p in polys
        append!(out, irreducible_factors(p))
    end
    return dedup_proportional(out)
end

"""Number of terms of a polynomial (upstream uses `Length[poly]`, the number of
top-level summands, as its size measure for the `limit` cutoff)."""
nterms(p) = length(p)

"""Degree of `p` in generator index `vi`."""
degree_in(p, vi::Int) = Nemo.degree(p, vi)

"""
    coeffs_in_var(p, vi) -> Vector

Coefficients of `p` viewed as a univariate polynomial in generator `vi`
(elements of the same mpoly ring, free of generator `vi`), from degree 0
up to degree d. Equivalent of upstream `CoefficientList[p, var]`.
"""
function coeffs_in_var(p, vi::Int)
    d = degree_in(p, vi)
    return [Nemo.coeff(p, [vi], [k]) for k in 0:d]
end

"""
    to_univariate_in(p, vi) -> (upoly, ring)

View mpoly `p` as a univariate polynomial in generator `vi` with coefficients
in the same mpoly ring (used for resultants/discriminants).
"""
function to_univariate_in(p, vi::Int)
    R = parent(p)
    U, _ = Nemo.polynomial_ring(R, :__t)
    return U(coeffs_in_var(p, vi))
end

"""
    resultant_in(p, q, vi)

Resultant of `p` and `q` with respect to generator `vi`. For moderate degree
pairs this uses upstream's coefficient-compression trick
(`generateAnsatzDISPATCH`, scr.m:420): the resultant of two GENERIC
polynomials of degrees (d1, d2) is computed once over ZZ[a0..ad1, b0..bd2]
and cached; the actual resultant is its evaluation at the coefficient
polynomials. This avoids subresultant-PRS coefficient swell over large
multivariate coefficient rings. High-degree pairs fall back to PRS.
"""
function resultant_in(p, q, vi::Int)
    # FLINT 3 native multivariate resultant (modular/interpolation): by far
    # the fastest path; fall back to compression / PRS where unavailable.
    try
        return Nemo.resultant(p, q, vi)
    catch
    end
    d1, d2 = degree_in(p, vi), degree_in(q, vi)
    if 0 < d1 && 0 < d2 && d1 + d2 <= 12
        gen_ = generic_resultant(d1, d2)
        vals = vcat(coeffs_in_var(p, vi), coeffs_in_var(q, vi))
        return Nemo.evaluate(gen_, vals)
    end
    up = to_univariate_in(p, vi)
    uq = to_univariate_in(q, vi)
    return Nemo.resultant(up, uq)
end

const _GENERIC_RES_CACHE = Dict{Tuple{Int,Int},Any}()

"""Resultant of generic degree-(d1,d2) univariate polynomials as a polynomial
in their coefficients a0..a_d1, b0..b_d2 over ZZ (cached)."""
function generic_resultant(d1::Int, d2::Int)
    get!(_GENERIC_RES_CACHE, (d1, d2)) do
        names = vcat(["a$i" for i in 0:d1], ["b$i" for i in 0:d2])
        R, g = Nemo.polynomial_ring(Nemo.ZZ, names)
        U, t = Nemo.polynomial_ring(R, :t)
        f = sum(g[i+1] * t^i for i in 0:d1)
        h = sum(g[d1+2+i] * t^i for i in 0:d2)
        Nemo.resultant(f, h)
    end
end

"""
    discriminant_in(p, vi)

Discriminant of `p` with respect to generator `vi`, up to a constant and a
leading-coefficient factor: computed as res(p, dp/dvi). Upstream pairs every
discriminant with the leading coefficient anyway (Fubini's `term1`), so the
extra lc factor res(f,f') carries relative to the true discriminant is
harmless — everything is factored and deduped downstream.
"""
function discriminant_in(p, vi::Int)
    dp = Nemo.derivative(p, vi)
    iszero(dp) && return zero(parent(p))
    return resultant_in(p, dp, vi)
end

"""Leading coefficient of `p` wrt generator `vi`; if `p` is free of `vi`,
returns `p` itself (matching upstream `CoefficientList[f,var][[-1]]`)."""
leading_coeff_in(p, vi::Int) = Nemo.coeff(p, [vi], [degree_in(p, vi)])

"""True if `p` does not involve generator `vi`."""
free_of(p, vi::Int) = degree_in(p, vi) == 0

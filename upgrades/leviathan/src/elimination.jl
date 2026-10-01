# elimination.jl — exact elimination over Q(s)[x0,x_S] (Eq. 8 leading-coefficient criterion).
# Exact route: Oscar/Singular GB over the rational function field (no FiniteFlow reconstruction);
# candidates = irreducible factors of the LCM of denominators of the monic elimination
# polynomial coefficients f_k = c_k/c_m, over all variables x_i (i in {0} u S).

"Function-field specialization: K = QQ(kin), Rk = K[x0, x_i (i in S)]."
function ff_ring(fam::Family, S)
    K, kg = Oscar.rational_function_field(Nemo.QQ, fam.kinnames)
    Sx = sort(collect(S))
    Rk, vk = Oscar.polynomial_ring(K, vcat(["x0"], fam.xnames[Sx]))
    x0 = vk[1]
    xv = Dict{Int,eltype(vk)}(i => vk[1+k] for (k, i) in enumerate(Sx))
    images = vcat([Rk(kg[j]) for j in eachindex(fam.kin)],
                  [haskey(xv, i) ? xv[i] : zero(Rk) for i in 1:nedges(fam)])
    h = Oscar.hom(fam.R, Rk, c -> K(c), images)
    (K, Rk, x0, xv, h(fam.G))
end

"Ideal generators over K: regulated (generic rational nu,d) or nu=0 sector ideal J_S."
function ff_ideal(fam::Family, S; regulated::Bool=true, rng=Random.default_rng())
    K, Rk, x0, xv, Gk = ff_ring(fam, S)
    iszero(Gk) && return nothing
    Ssort = sort(collect(S))
    gs = typeof(Gk)[]
    if regulated
        nu, d2 = generic_regulators(rng, nedges(fam))
        noncut = [i for i in Ssort if !(i in fam.cut)]
        for i in Ssort
            di = Oscar.derivative(Gk, xv[i])
            push!(gs, i in fam.cut ? di : K(nu[i]) * Gk - K(d2) * xv[i] * di)
        end
        sat = isempty(noncut) ? Gk : prod(xv[i] for i in noncut) * Gk
        push!(gs, 1 - x0 * sat)
    else
        for i in Ssort
            push!(gs, Oscar.derivative(Gk, xv[i]))
        end
        push!(gs, 1 - x0 * Gk)
    end
    (Rk, Oscar.ideal(Rk, gs))
end

"""
Per-variable elimination generators of a 0-dim ideal over K, via the minimal polynomial
of x_i in the quotient ring Q = K[x]/I (the SPQR companion-matrix object, computed with
one DRL GB + normal forms instead of per-variable elimination orderings).
Returns Dict varname => (m, lc, denom_lcm, fs):
m = degree, fs = monic coefficients f_k = c_k/c_m in K, denom_lcm = LCM of their
denominators, lc = c_m as recovered from that LCM (the leading-coefficient locus).
"""
function elim_generators(Rk, I)
    K = Oscar.base_ring(Rk)
    o = Oscar.degrevlex(Rk)
    gb = Oscar.groebner_basis(I, ordering=o, complete_reduction=true)
    les = lead_exps(gb, Rk)
    nsm = staircase_count(les, Oscar.nvars(Rk))
    nsm === nothing && return nothing   # positive-dimensional (Type-1.1): caller skips sector
    out = Dict{String,Any}()
    for v in Oscar.gens(Rk)
        # incremental normal forms 1, nf(v*prev), ... until K-linear dependence:
        # the monic minimal polynomial of v in K[x]/I (inputs stay in the standard span)
        nfs = [Oscar.normal_form(one(Rk), I, ordering=o)]
        monidx = Dict{Vector{Int},Int}()
        for e in Nemo.exponent_vectors(nfs[1])
            monidx[e] = length(monidx) + 1
        end
        row(f, nc) = (r = [zero(K) for _ in 1:nc];
                      for (c, e) in zip(Nemo.coefficients(f), Nemo.exponent_vectors(f))
                          r[monidx[e]] = c
                      end; r)
        m = nothing
        fs = nothing
        for k in 1:nsm
            push!(nfs, Oscar.normal_form(v * nfs[end], I, ordering=o))
            for e in Nemo.exponent_vectors(nfs[end])
                haskey(monidx, e) || (monidx[e] = length(monidx) + 1)
            end
            nc = length(monidx)
            A = Oscar.matrix(K, k, nc, vcat([row(nfs[j], nc) for j in 1:k]...))
            b = Oscar.matrix(K, 1, nc, row(nfs[k+1], nc))
            ok, x = Oscar.can_solve_with_solution(A, b, side=:left)
            if ok
                m = k
                fs = [-x[1, j] for j in 1:k]   # v^m - sum_j x_j v^j = 0 (monic)
                break
            end
        end
        m === nothing && continue
        dl = isempty(fs) ? Nemo.denominator(one(K)) :
             reduce(lcm, [Nemo.denominator(f) for f in fs])
        out[string(v)] = (m=m, lc=K(dl), denom_lcm=dl, fs=fs)
    end
    out
end

"Canonical factor string: monic-normalized irreducible factor over QQ[kin]."
function _canon_factor(f)
    lc = Oscar.leading_coefficient(f)
    string(f * inv(lc))
end

"""
Candidate singularity factors for one sector. Returns (factors::Set{String}, elim_info).
Factors of the denominator-LCMs (and of the leading coefficients) over all variables;
numeric factors dropped.
"""
function sector_candidates(fam::Family, S; regulated::Bool=true, rng=Random.default_rng())
    r = ff_ideal(fam, S; regulated=regulated, rng=rng)
    r === nothing && return (Set{String}(), nothing)
    Rk, I = r
    info = elim_generators(Rk, I)
    if info === nothing
        @warn "sector $S positive-dimensional over K (Type 1.1) — no candidates from it"
        return (Set{String}(), nothing)
    end
    facs = Set{String}()
    for (_, rec) in info
        for poly in (rec.denom_lcm, Nemo.numerator(rec.lc))
            (iszero(poly) || Oscar.is_constant(poly)) && continue
            for (f, _) in Oscar.factor(poly)
                Oscar.is_constant(f) || push!(facs, _canon_factor(f))
            end
        end
    end
    (facs, info)
end

"Candidates over a list of sectors (default: full family top sector + all admissible)."
function all_candidates(fam::Family; sectors=nothing, regulated::Bool=true,
                        rng=Random.default_rng())
    secs = sectors === nothing ? admissible_sectors(fam) : sectors
    bysec = Dict{Vector{Int},Set{String}}()
    for S in secs
        facs, _ = sector_candidates(fam, S; regulated=regulated, rng=rng)
        bysec[S] = facs
    end
    (union(values(bysec)...), bysec)
end

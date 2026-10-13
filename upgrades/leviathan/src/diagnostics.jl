# diagnostics.jl: the four Appendix-A checks of arXiv:2606.29612. Each returns (pass::Bool, evidence).
# nodegeneracy (Type 1.2, lines 60-99) is a translation into Julia/Oscar of degenerateForKinAndVar
# (:33-63) and noDegeneracyQ (:66-116) in landau_codes.m @ 47028b7 by Chestnov, Crisanti and Giroux
# (no license file); the MIT grant does not cover it (NOTICE). The retry over primes in diag_type22
# (lines 41-47) follows countInSector in DiscKosky euler_characteristic.m:43-52. See PATCHES.md.

"Type 1.1: every sector ideal J_S must be zero-dimensional (chi_S finite)."
function diag_type11(fam::Family; p::Int=default_p31(), rng=Random.default_rng())
    kv = random_kinvals(rng, fam.kinnames, p)
    chis = Dict{Vector{Int},Union{Int,Nothing}}()
    for S in admissible_sectors(fam)
        chis[S] = chi_sector(fam, S; p=p, kinvals=kv, rng=rng)
    end
    bad = [S for (S, c) in chis if c === nothing]
    (isempty(bad), (chi_sectors=chis, nonzero_dim_sectors=sort(bad)))
end

"Type 2.1: chi(regulated, generic nu,d) == sum over sectors of chi_S(nu=0)."
function diag_type21(fam::Family; p::Int=default_p31(), rng=Random.default_rng())
    kv = random_kinvals(rng, fam.kinnames, p)
    chi_reg = chi_regulated(fam; p=p, kinvals=kv, rng=rng)
    chis = [chi_sector(fam, S; p=p, kinvals=kv, rng=rng) for S in admissible_sectors(fam)]
    any(c -> c === nothing, chis) && return (false, (chi_reg=chi_reg, chi_sbs=nothing,
        note="some sector not 0-dim: report Type 1.1"))
    s = sum(chis)
    (chi_reg == s, (chi_reg=chi_reg, chi_sbs=s))
end

"""
Type 2.2: candidate verification against Eq. 5 at generic nu — for each candidate factor l
(string, parsed in kin variables), chi_regulated on l=0 must DROP below generic.
Returns pass iff every candidate is genuine; evidence lists per-candidate verdicts.
"""
function diag_type22(fam::Family, candidates; p::Int=default_p31(), rng=Random.default_rng(),
                     repeats::Int=2)
    chi_gen = chi_regulated(fam; p=p, rng=rng)
    tab = Dict{Symbol,eltype(fam.kin)}(Symbol(n) => fam.kin[i] for (i, n) in enumerate(fam.kinnames))
    verdicts = Dict{String,Any}()
    for lstr in candidates
        l = poly_from_string(lstr, tab, fam.R)
        # retry over the prime list, as countInSector in DiscKosky euler_characteristic.m:43-52 does:
        # an algebraic locus (e.g. y^2-10y+5) has F_p points only when its discriminant is a square mod p
        chis = Union{Int,Nothing}[]
        for q in P31
            push!(chis, chi_regulated_constrained(fam, l; p=q, rng=rng))
            count(c -> c !== nothing, chis) >= repeats && break
        end
        chis = [c for c in chis if c !== nothing]
        if isempty(chis)
            verdicts[lstr] = (genuine=nothing, chi_on_locus=nothing, note="no point on locus found")
        else
            c = minimum(chis)
            verdicts[lstr] = (genuine=(c < chi_gen), chi_on_locus=c, chi_generic=chi_gen)
        end
    end
    bad = [l for (l, v) in verdicts if v.genuine === false]
    (isempty(bad), (chi_generic=chi_gen, verdicts=verdicts, spurious=sort(bad)))
end

"""
Type 1.2 noDegeneracyQ: for every (kinematic q, variable v in {x0,x_i}) pair, promote q to a
ring variable, demote v to a parameter at TWO random values; eliminate all x-variables mod p
(p < 2^29, Singular route); GCD of the two univariate-in-q generators must be 1.
Nontrivial GCD = potential l_degen (genuine but missable singularity).
"""
function nodegeneracy(fam::Family; p::Int=default_p29(), rng=Random.default_rng())
    Fp = Nemo.GF(p)
    offenders = Dict{Tuple{String,String},String}()
    for qn in fam.kinnames
        kv = random_kinvals(rng, fam.kinnames, p)
        Rm, qv, x0, xvs, Gm = spec_ring_midq(fam, qn, p, kv)
        # nu=0 ideal <dG/dx_i, 1 - x0*G> (arXiv:2606.29612v1 Eq. 10), as the paper (Appendix A) and
        # noDegeneracyQ in Chestnov, Crisanti and Giroux's landau_codes.m both specify; this function
        # is a translation of that routine (file header, PATCHES.md). The regulated ideal is not used:
        # the test would pass on loci where J_{nu=0} is degenerate but J_reg is not.
        gs0 = typeof(Gm)[Oscar.derivative(Gm, xvs[i]) for i in 1:nedges(fam)]
        push!(gs0, 1 - x0 * Gm)
        allx = vcat([x0], xvs)
        for v in allx
            polys = []
            for _ in 1:2
                a = Rm(Fp(rand(rng, 2:p-2)))
                sub = [w == v ? a : w for w in Oscar.gens(Rm)]
                gsv = [Oscar.evaluate(g, sub) for g in gs0]
                keepx = [w for w in allx if w != v]
                E = Oscar.eliminate(Oscar.ideal(Rm, gsv), keepx)
                gsE = [g for g in Oscar.gens(E) if !iszero(g)]
                isempty(gsE) && continue
                push!(polys, gsE[argmin([Oscar.total_degree(g) for g in gsE])])
            end
            length(polys) < 2 && continue
            g = Oscar.gcd(polys[1], polys[2])
            if Oscar.total_degree(g) > 0
                offenders[(qn, string(v))] = string(g)
            end
        end
    end
    (isempty(offenders), (offenders=offenders,))
end

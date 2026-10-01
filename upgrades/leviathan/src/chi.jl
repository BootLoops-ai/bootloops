# chi.jl — Euler-characteristic counts mod p (Eq. 12 sector ideals, regulated ideal, constraints).
# The cut-edge form of the regulated ideal (chi_regulated: plain dG/dx_i for a cut edge,
# cut edges excluded from the saturation product) and the "Diophantine" point search of
# point_on_locus (lowest-degree kinematic variable, up to 400 tries) follow the behavior of
# countInSector / prepareCountRegulatedMS in DiscKosky's euler_characteristic.m (Crisanti,
# Lippstreu, McLeod, Polackova; no license file) — read as the reference, no code copied.

"""
chi_S of sector S at nu=0:  J_S = < dG_S/dx_i (i in S), 1 - x0*G_S >.
Returns Int, `nothing` (not zero-dimensional, Type-1.1), or 0 for skipped (G_S == 0).
"""
function chi_sector(fam::Family, S; p::Int=default_p31(), kinvals=nothing,
                    rng=Random.default_rng())
    kv = kinvals === nothing ? random_kinvals(rng, fam.kinnames, p) : kinvals
    Rp, x0, xv, Gp = spec_ring(fam, S, p, kv)
    iszero(Gp) && return 0
    gs = [Oscar.derivative(Gp, xv[i]) for i in sort(collect(S))]
    push!(gs, 1 - x0 * Gp)
    chi_of_ideal(Oscar.ideal(Rp, gs), Rp)
end

"""
Regulated chi (generic nu_i, d): generators nu_i*G - (d/2)*x_i*dG/dx_i for non-cut i in S,
plain dG/dx_i for cut i in S, saturation 1 - x0*(prod of non-cut x_i in S)*G.
nu_i, d/2 -> random F_p values (generic-parameter degree, prob-1 correct).
"""
function chi_regulated(fam::Family, S=collect(1:nedges(fam)); p::Int=default_p31(),
                       kinvals=nothing, rng=Random.default_rng())
    kv = kinvals === nothing ? random_kinvals(rng, fam.kinnames, p) : kinvals
    Rp, x0, xv, Gp = spec_ring(fam, S, p, kv)
    iszero(Gp) && return 0
    Fp = Nemo.base_ring(Rp)
    d2 = Fp(rand(rng, 2:p-2))
    Ssort = sort(collect(S))
    noncut = [i for i in Ssort if !(i in fam.cut)]
    gs = typeof(Gp)[]
    for i in Ssort
        di = Oscar.derivative(Gp, xv[i])
        if i in fam.cut
            push!(gs, di)
        else
            push!(gs, Fp(rand(rng, 2:p-2)) * Gp - d2 * xv[i] * di)
        end
    end
    sat = isempty(noncut) ? Gp : prod(xv[i] for i in noncut) * Gp
    push!(gs, 1 - x0 * sat)
    chi_of_ideal(Oscar.ideal(Rp, gs), Rp)
end

"""
Random F_p point on the locus l(s)=0 (l a kin-only polynomial in fam.R).
Picks the kin variable v of lowest positive degree in l; samples the others,
solves for v via univariate roots mod p; <= maxtries retries ("Diophantine" mode).
Returns Dict kinvals or nothing.
"""
function point_on_locus(fam::Family, l; p::Int=default_p31(), rng=Random.default_rng(),
                        maxtries::Int=400)
    Fp = Nemo.GF(p)
    degs = [(Oscar.degree(l, fam.kin[j]), j) for j in eachindex(fam.kin)]
    degs = [(d, j) for (d, j) in degs if d > 0]
    isempty(degs) && return nothing
    j0 = degs[argmin(first.(degs))][2]
    Pu, t = Oscar.polynomial_ring(Fp, "t")
    for _ in 1:maxtries
        kv = random_kinvals(rng, fam.kinnames, p)
        images = vcat([jj == j0 ? t : Pu(Fp(kv[fam.kinnames[jj]])) for jj in eachindex(fam.kin)],
                      [zero(Pu) for _ in fam.xs])
        h = Oscar.hom(fam.R, Pu, c -> _fpcoeff(Fp, c), images)
        lu = h(l)
        iszero(lu) && continue
        rts = Oscar.roots(lu)
        isempty(rts) && continue
        kv[fam.kinnames[j0]] = Int(Nemo.lift(Nemo.ZZ, rand(rng, rts)))
        return kv
    end
    nothing
end

"chi_regulated constrained to the locus l=0 (Type-2.2 verification input)."
function chi_regulated_constrained(fam::Family, l; p::Int=default_p31(),
                                   rng=Random.default_rng())
    kv = point_on_locus(fam, l; p=p, rng=rng)
    kv === nothing && return nothing
    chi_regulated(fam; p=p, kinvals=kv, rng=rng)
end

# family.jl — Family struct (Lee–Pomeransky G + kinematic parameters) and specializations.

struct Family
    R::Oscar.MPolyRing          # QQ[kin..., x1..xE]
    G::Oscar.MPolyRingElem
    kin::Vector                 # kin generators in R
    xs::Vector                  # edge-variable generators in R
    kinnames::Vector{String}
    xnames::Vector{String}
    cut::Vector{Int}            # cut edge indices (excluded from regulators & saturation product)
end

nedges(fam::Family) = length(fam.xs)

function family(Gstr::AbstractString, kinnames, xnames; cut=Int[])
    kn, xn = collect(String, kinnames), collect(String, xnames)
    R, v = Oscar.polynomial_ring(Nemo.QQ, vcat(kn, xn))
    tab = Dict{Symbol,eltype(v)}(Symbol(n) => v[i] for (i, n) in enumerate(vcat(kn, xn)))
    G = poly_from_string(Gstr, tab, R)
    Family(R, G, v[1:length(kn)], v[length(kn)+1:end], kn, xn, collect(Int, cut))
end

family_UF(U::AbstractString, F::AbstractString, kinnames, xnames; cut=Int[]) =
    family("($U) + ($F)", kinnames, xnames; cut=cut)

"G with x_i -> 0 for i not in S (sector contraction), result still in fam.R."
function sector_G(fam::Family, S)
    images = vcat(fam.kin, [(i in S) ? fam.xs[i] : zero(fam.R) for i in 1:nedges(fam)])
    Oscar.evaluate(fam.G, images)
end

"All admissible sectors: nonempty subsets of edges with G_S not identically zero."
function admissible_sectors(fam::Family)
    E = nedges(fam)
    out = Vector{Vector{Int}}()
    for mask in 1:(2^E - 1)
        S = [i for i in 1:E if (mask >> (i - 1)) & 1 == 1]
        iszero(sector_G(fam, S)) || push!(out, S)
    end
    sort(out, by=S -> (length(S), S))
end

_fpcoeff(Fp, c) = Fp(Nemo.numerator(c)) * inv(Fp(Nemo.denominator(c)))

"""
Specialize fam to GF(p): ring Rp = GF(p)[x0, x_i (i in S)], kinematics -> kinvals mod p.
Returns (Rp, x0, xv::Dict edgeindex=>gen, Gp) where Gp = image of G (sector contraction
is automatic: x_i -> 0 for i not in S).
"""
function spec_ring(fam::Family, S, p::Int, kinvals::Dict{String,Int})
    Fp = Nemo.GF(p)
    Sx = sort(collect(S))
    Rp, vp = Oscar.polynomial_ring(Fp, vcat(["x0"], fam.xnames[Sx]))
    x0 = vp[1]
    xv = Dict{Int,eltype(vp)}(i => vp[1+k] for (k, i) in enumerate(Sx))
    images = vcat([Rp(Fp(kinvals[n])) for n in fam.kinnames],
                  [haskey(xv, i) ? xv[i] : zero(Rp) for i in 1:nedges(fam)])
    h = Oscar.hom(fam.R, Rp, c -> _fpcoeff(Fp, c), images)
    (Rp, x0, xv, h(fam.G))
end

"""
Mid-ring for noDegeneracyQ: GF(p)[q, x0, x1..xE] with ALL edge vars kept as ring vars,
kinematic q kept as ring var, other kinematics -> kinvals. Returns (Rm, qv, x0, xvs, Gm).
"""
function spec_ring_midq(fam::Family, qname::String, p::Int, kinvals::Dict{String,Int})
    Fp = Nemo.GF(p)
    Rm, vm = Oscar.polynomial_ring(Fp, vcat([qname, "x0"], fam.xnames))
    qv, x0 = vm[1], vm[2]
    xvs = vm[3:end]
    images = vcat([n == qname ? qv : Rm(Fp(kinvals[n])) for n in fam.kinnames], xvs)
    h = Oscar.hom(fam.R, Rm, c -> _fpcoeff(Fp, c), images)
    (Rm, qv, x0, xvs, h(fam.G))
end

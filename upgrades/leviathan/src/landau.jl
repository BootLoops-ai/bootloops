# landau.jl — driver: chi accounting + candidates + Type-2.2 verification.

struct LandauResult
    chi_generic::Union{Int,Nothing}        # regulated chi, full family
    chi_sectors::Dict{Vector{Int},Union{Int,Nothing}}
    candidates::Vector{String}             # raw candidate factors (union over sectors)
    genuine::Vector{String}                # passed Type-2.2 chi-drop verification
    spurious::Vector{String}               # failed it (no drop)
    unverified::Vector{String}             # could not find a point on the locus
    by_sector::Dict{Vector{Int},Set{String}}
end

"""
Full pipeline: per-sector chi (nu=0), regulated chi, candidate factors from per-sector
function-field elimination (regulated and/or nu=0 lanes), then chi-drop verification.
"""
function landau(fam::Family; p::Int=default_p31(), rng=Random.default_rng(),
                lanes=(:regulated, :nu0), verify::Bool=true, repeats::Int=2)
    kv = random_kinvals(rng, fam.kinnames, p)
    chis = Dict{Vector{Int},Union{Int,Nothing}}()
    for S in admissible_sectors(fam)
        chis[S] = chi_sector(fam, S; p=p, kinvals=kv, rng=rng)
    end
    chi_gen = chi_regulated(fam; p=p, kinvals=kv, rng=rng)
    cands = Set{String}()
    bysec = Dict{Vector{Int},Set{String}}()
    for lane in lanes
        c, bs = all_candidates(fam; regulated=(lane == :regulated), rng=rng)
        union!(cands, c)
        for (S, f) in bs
            union!(get!(bysec, S, Set{String}()), f)
        end
    end
    cl = sort(collect(cands))
    genuine, spurious, unverified = String[], String[], String[]
    if verify
        _, ev = diag_type22(fam, cl; p=p, rng=rng, repeats=repeats)
        for l in cl
            v = ev.verdicts[l]
            v.genuine === true ? push!(genuine, l) :
            v.genuine === false ? push!(spurious, l) : push!(unverified, l)
        end
    else
        genuine = cl
    end
    LandauResult(chi_gen, chis, cl, genuine, spurious, unverified, bysec)
end

function report(r::LandauResult)
    io = IOBuffer()
    println(io, "chi (regulated, generic) = ", r.chi_generic)
    secs = sort(collect(keys(r.chi_sectors)), by=S -> (length(S), S))
    tot = 0; ok = true
    for S in secs
        c = r.chi_sectors[S]
        println(io, "  chi_S ", S, " = ", c === nothing ? "Indeterminate" : c)
        c === nothing ? (ok = false) : (tot += c)
    end
    println(io, "  sum over sectors = ", ok ? tot : "Indeterminate",
            ok && tot == r.chi_generic ? "  [additivity OK]" : "  [MISMATCH/INDET]")
    println(io, "candidate factors: ", r.candidates)
    println(io, "  genuine   : ", r.genuine)
    println(io, "  spurious  : ", r.spurious)
    isempty(r.unverified) || println(io, "  unverified: ", r.unverified)
    String(take!(io))
end

"Run all four diagnostics; returns Dict name => (pass, evidence)."
function diagnose(fam::Family; p31::Int=default_p31(), p29::Int=default_p29(),
                  rng=Random.default_rng(), candidates=nothing)
    out = Dict{String,Any}()
    out["type1.1_sector_dim"] = diag_type11(fam; p=p31, rng=rng)
    out["type1.2_nodegeneracy"] = nodegeneracy(fam; p=p29, rng=rng)
    out["type2.1_chi_additivity"] = diag_type21(fam; p=p31, rng=rng)
    if candidates !== nothing
        out["type2.2_candidate_verification"] = diag_type22(fam, candidates; p=p31, rng=rng)
    end
    out
end

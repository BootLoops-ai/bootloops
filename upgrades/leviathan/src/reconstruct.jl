# reconstruct.jl — finite-field reconstruction route (the v2 "leviathan-scale" path).
# Computes the SAME monic minpoly μ_v as elimination.jl (provably the identical object; see
# the method notes in PATCHES.md), but by sampling over (prime p, kin point k) and rationally
# reconstructing each
# coefficient instead of exact GB over ℚ(kin). Candidates = irreducible factors of N and D of
# every reconstructed coefficient; certified by existing diag_type22.

import JSON

"J_S over GF(p)[x0,x_S] at (p,kinvals). Mirrors ff_ideal; regulated takes FIXED rational nu,d/2."
function fp_ideal(fam::Family, S, p::Int, kinvals; regulated::Bool, nu=nothing, d2=nothing)
    Rp, x0, xv, Gp = spec_ring(fam, S, p, kinvals)
    iszero(Gp) && return nothing
    Fp = Nemo.base_ring(Rp); Ss = sort(collect(S)); gs = typeof(Gp)[]
    if regulated
        for i in Ss
            di = Oscar.derivative(Gp, xv[i])
            push!(gs, i in fam.cut ? di : _fpcoeff(Fp, nu[i]) * Gp - _fpcoeff(Fp, d2) * xv[i] * di)
        end
        nc = [i for i in Ss if !(i in fam.cut)]
        push!(gs, 1 - x0 * (isempty(nc) ? Gp : prod(xv[i] for i in nc) * Gp))
    else
        append!(gs, [Oscar.derivative(Gp, xv[i]) for i in Ss]); push!(gs, 1 - x0 * Gp)
    end
    (Rp, Oscar.ideal(Rp, gs))
end

"Monic-minpoly coefficients of every variable in Rp[x]/I over GF(p). Same recursion as elim_generators."
function minpoly_modp(Rp, I)
    Fp = Oscar.base_ring(Rp); o = Oscar.degrevlex(Rp)
    gb = Oscar.groebner_basis_f4(I)
    nsm = staircase_count(lead_exps(gb, Rp), Oscar.nvars(Rp))
    nsm === nothing && return nothing
    out = Dict{String,Vector{Int}}()
    for v in Oscar.gens(Rp)
        nfs = [one(Rp)]; mi = Dict{Vector{Int},Int}([0 for _ in 1:Oscar.nvars(Rp)] => 1)
        row(f, nc) = (r = zeros(Fp, nc);
            for (c, e) in zip(Nemo.coefficients(f), Nemo.exponent_vectors(f)); r[mi[collect(Int, e)]] = c; end; r)
        for k in 1:nsm
            push!(nfs, Oscar.normal_form(v * nfs[end], I, ordering=o))
            for e in Nemo.exponent_vectors(nfs[end]); get!(mi, collect(Int, e), length(mi) + 1); end
            nc = length(mi)
            A = Oscar.matrix(Fp, k, nc, vcat([row(nfs[j], nc) for j in 1:k]...))
            b = Oscar.matrix(Fp, 1, nc, row(nfs[k+1], nc))
            ok, x = Oscar.can_solve_with_solution(A, b, side=:left)
            ok && (out[string(v)] = [Int(Nemo.lift(Nemo.ZZ, -x[1, j])) for j in 1:k]; break)
        end
    end
    out
end

"One sample: (kinvals_dehom, Dict varname=>coeffs) at (p, kin point with kin[end]=1)."
function _one_sample(fam, S, p, regulated, seed, nu, d2)
    rng = MersenneTwister(seed)
    kv = random_kinvals(rng, fam.kinnames, p); kv[fam.kinnames[end]] = 1
    r = fp_ideal(fam, S, p, kv; regulated=regulated, nu=nu, d2=d2)
    r === nothing && return nothing
    mp = minpoly_modp(r...)
    mp === nothing && return nothing
    ([kv[n] for n in fam.kinnames[1:end-1]], mp)
end

"Sample grid. msolve F4 / Singular normal_form are NOT thread-safe (measured: 8-thread runs
gave silently wrong coeffs → empty factor sets). Serial here AND in landau_reconstruct (this
port contains no Distributed/pmap code); parallelize across independent OS processes, each
owning its own Oscar. JSON checkpoint per (sector, prime)."
function sample_sector(fam::Family, S; primes, npoints::Int, regulated::Bool=false,
                       seed0::Int=20260630, savedir=nothing, tag="fam", nu=nothing, d2=nothing)
    out = Dict{Int,Vector{Tuple{Vector{Int},Dict{String,Vector{Int}}}}}()
    for p in primes
        buf = [_one_sample(fam, S, p, regulated, hash((seed0, p, j)) % typemax(Int32), nu, d2)
               for j in 1:npoints]
        out[p] = [b for b in buf if b !== nothing]
        if savedir !== nothing
            mkpath(savedir)
            open(joinpath(savedir, "$(tag)_S$(join(S,'-'))_p$(p).json"), "w") do io
                JSON.print(io, [(k, mp) for (k, mp) in out[p]])
            end
        end
    end
    out
end

"Reconstruct candidate factors for one sector from samples → Set{String} (canonical monic)."
function reconstruct_sector(fam::Family, S, samples; deg_cap::Int=12)
    nkin = length(fam.kinnames)
    RkinFull, _ = Oscar.polynomial_ring(Nemo.QQ, fam.kinnames)
    RkinQ, _ = Oscar.polynomial_ring(Nemo.QQ, fam.kinnames[1:max(1, nkin - 1)])
    primes = sort([p for p in keys(samples) if !isempty(samples[p])])
    # degenerate sector (e.g. positive-dimensional / chi-Indeterminate): no valid samples
    isempty(primes) && return Set{String}()
    # determine generic minpoly degree per var (max over samples = generic w.p. 1-O(1/p))
    vars = sort(collect(keys(samples[primes[1]][1][2])))
    mdeg = Dict(v => maximum(length(mp[v]) for p in primes for (_, mp) in samples[p]) for v in vars)
    facs = Set{String}()
    for v in vars, k in 0:mdeg[v]-1
        sbp = Dict{Int,Vector{Tuple{Vector{Int},Int}}}()
        for p in primes
            sbp[p] = [(kvec, mp[v][k+1]) for (kvec, mp) in samples[p] if length(mp[v]) == mdeg[v]]
        end
        r = ratrec_multiprime(sbp, RkinQ; deg_cap=deg_cap)
        r === nothing && continue
        # Factor BOTH numerator and denominator: the exact route sees only c_m factors, but when
        # c_m exceeds deg_cap (measured deg≈76 on a 3-loop validation sector) the only fits that
        # close are small-c_k coefficients, whose numerator factors are χ-drop-filtered →
        # genuine ones recovered (that sector: 7/16 vs 1/16).
        for poly in r
            (iszero(poly) || Oscar.is_constant(poly)) && continue
            ph = nkin > 1 ? rehomogenize(poly, RkinFull, Oscar.total_degree(poly)) : poly
            for (f, _) in Oscar.factor(ph)
                Oscar.is_constant(f) || push!(facs, _canon_factor(f))
            end
        end
    end
    facs
end

"""
landau_reconstruct(fam; nprimes, npoints, deg_cap, sectors, lanes, verify, savedir, tag)
Finite-field reconstruction route: sample → reconstruct → factor → Type-2.2 verify.
Returns LandauResult (same struct as the exact route).
"""
function landau_reconstruct(fam::Family; nprimes::Int=3, npoints::Int=60, deg_cap::Int=12,
                            sectors=nothing, lanes=(:nu0,), verify::Bool=true,
                            seed0::Int=20260630, savedir=nothing, tag="fam", repeats::Int=2)
    rng = MersenneTwister(seed0)
    primes = P31[1:nprimes]
    nu, d2 = generic_regulators(rng, nedges(fam))   # FIXED across all samples & primes
    secs = sectors === nothing ? admissible_sectors(fam) : sectors
    kv = random_kinvals(rng, fam.kinnames, primes[1])
    chis = Dict{Vector{Int},Union{Int,Nothing}}(S => chi_sector(fam, S; p=primes[1], kinvals=kv) for S in secs)
    chi_gen = chi_regulated(fam; p=primes[1], kinvals=kv)
    cands = Set{String}(); bysec = Dict{Vector{Int},Set{String}}()
    for S in secs, lane in lanes
        smp = sample_sector(fam, S; primes=primes, npoints=npoints, nu=nu, d2=d2,
                            regulated=(lane == :regulated), seed0=seed0, savedir=savedir, tag=tag)
        f = reconstruct_sector(fam, S, smp; deg_cap=deg_cap)
        union!(get!(bysec, S, Set{String}()), f); union!(cands, f)
    end
    push!(cands, fam.kinnames[end])    # hyperplane-at-infinity (dehom variable)
    cl = sort(collect(cands))
    genuine, spurious, unverified = String[], String[], String[]
    if verify
        _, ev = diag_type22(fam, cl; p=primes[1], rng=rng, repeats=repeats)
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

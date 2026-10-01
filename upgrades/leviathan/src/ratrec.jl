# ratrec.jl — multivariate rational-function reconstruction from F_p samples.
# Dense-ansatz nullspace mod p (FireFly-shape, no Thiele): for each prime p and sample
# set {(kvec_j, val_j)}, solve [N_mons(k_j) | -val_j*D_mons(k_j)] · c = 0 over GF(p),
# gcd-reduce, normalize by leading_coefficient(D)=1, CRT across primes, Nemo.reconstruct→QQ.
# Ports the *idea* (screen-then-certify) of an earlier univariate fitter, not its code;
# fresh multivariate Julia.

"All monomial exponent vectors of total degree ≤ d in n vars (lex on (deg,exp))."
function monoms_upto(n::Int, d::Int)
    out = Vector{Vector{Int}}()
    e = zeros(Int, n)
    function rec(i, rem)
        if i == n
            e[n] = rem; push!(out, copy(e)); e[n] = 0; return
        end
        for k in 0:rem
            e[i] = k; rec(i + 1, rem - k)
        end
        e[i] = 0
    end
    for deg in 0:d
        rec(1, deg)
    end
    out
end

"Evaluate all monomials in `exps` at point `kp` (vector of field elements)."
_evalmons(exps, kp) = [prod(kp[j]^e[j] for j in eachindex(kp); init=one(kp[1])) for e in exps]

"""
ratrec_modp: reconstruct val ≈ N(kin)/D(kin) mod p from samples.
samples: Vector of (kvec::Vector{Int}, val::Int) at one prime p.
Rkin: GF(p)[kin...] polynomial ring (n vars).
Returns (Np, Dp, dN, dD) with gcd(Np,Dp)=1, lc(Dp)=1, or `nothing` if no fit ≤ deg_cap.
"""
function ratrec_modp(samples, Rkin, p::Int; deg_cap::Int=12)
    Fp = Oscar.base_ring(Rkin)
    n = Oscar.nvars(Rkin)
    kvs = [Fp.(s[1]) for s in samples]
    vls = [Fp(s[2]) for s in samples]
    if all(iszero, vls)
        return (zero(Rkin), one(Rkin), 0, 0)
    end
    for d in 0:deg_cap
        exN = monoms_upto(n, d); exD = monoms_upto(n, d)
        nN, nD = length(exN), length(exD)
        nN + nD > length(samples) && return nothing   # underdetermined: need more samples
        rows = Vector{elem_type(Fp)}()
        for (kp, vl) in zip(kvs, vls)
            mn = _evalmons(exN, kp)
            append!(rows, mn); append!(rows, [-vl * m for m in _evalmons(exD, kp)])
        end
        M = Oscar.matrix(Fp, length(samples), nN + nD, rows)
        K = Oscar.kernel(M, side=:right)
        K isa Tuple && (K = K[2])
        Oscar.ncols(K) == 0 && continue
        # any kernel column gives N,D with N(k)=v*D(k); take first, gcd-reduce
        Np = Rkin(elem_type(Fp)[K[i, 1] for i in 1:nN], exN)
        Dp = Rkin(elem_type(Fp)[K[nN + i, 1] for i in 1:nD], exD)
        iszero(Dp) && continue
        g = Oscar.gcd(Np, Dp)
        Np, Dp = divexact(Np, g), divexact(Dp, g)
        lc = Oscar.leading_coefficient(Dp)
        Np, Dp = Np * inv(lc), Dp * inv(lc)
        return (Np, Dp, Oscar.total_degree(Np), Oscar.total_degree(Dp))
    end
    nothing
end

"CRT-lift a vector of (poly mod p_i) with matching support → ZZ[kin] → QQ[kin] via Nemo.reconstruct."
function crt_lift_qq(polys, primes, RkinQ)
    isempty(polys) && return nothing
    # collect coefficient per exponent across primes
    sup = Set{Vector{Int}}()
    for f in polys, e in Nemo.exponent_vectors(f)
        push!(sup, collect(Int, e))
    end
    Ps = Nemo.ZZ.(primes); Pprod = prod(Ps)
    cfq = Dict{Vector{Int},Nemo.QQFieldElem}()
    for e in sup
        rs = [Nemo.ZZ(Nemo.lift(Nemo.ZZ, Oscar.coeff(polys[i], e))) for i in eachindex(polys)]
        z = Nemo.crt(rs, Ps)
        q = try Nemo.reconstruct(z, Pprod) catch; return nothing end
        cfq[e] = q
    end
    es = collect(sup)
    RkinQ([cfq[e] for e in es], es)
end

"""
ratrec_multiprime: per-prime ratrec_modp → check (dN,dD) agree → CRT+reconstruct N,D over QQ.
samples_by_p: Dict prime => Vector{(kvec,val)}.  RkinQ: QQ[kin...].
Returns (Nq, Dq) ∈ QQ[kin]^2 or nothing.
"""
function ratrec_multiprime(samples_by_p, RkinQ; deg_cap::Int=12)
    primes = sort(collect(keys(samples_by_p)))
    Nps, Dps, dims = [], [], Set{Tuple{Int,Int}}()
    for p in primes
        Rp, _ = Oscar.polynomial_ring(Nemo.GF(p), [string(v) for v in Oscar.gens(RkinQ)])
        r = ratrec_modp(samples_by_p[p], Rp, p; deg_cap=deg_cap)
        r === nothing && return nothing
        push!(Nps, r[1]); push!(Dps, r[2]); push!(dims, (r[3], r[4]))
    end
    length(dims) == 1 || return nothing   # bad-luck prime: degrees disagree → caller retries
    Nq = crt_lift_qq(Nps, primes, RkinQ)
    Dq = crt_lift_qq(Dps, primes, RkinQ)
    (Nq === nothing || Dq === nothing) && return nothing
    (Nq, Dq)
end

"Rehomogenize a poly in n−1 vars by appending kin[end] to total degree `deg` in RkinFull."
function rehomogenize(f, RkinFull, deg::Int)
    n = Oscar.nvars(RkinFull)
    cs = collect(Nemo.coefficients(f)); es = [collect(Int, e) for e in Nemo.exponent_vectors(f)]
    es2 = [vcat(e, [deg - sum(e)]) for e in es]
    RkinFull(cs, es2)
end

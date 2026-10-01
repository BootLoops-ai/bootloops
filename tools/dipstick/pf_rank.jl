# DIPSTICK member — verb `order`: Griffiths-Dwork mod-p PF order + holonomic
# rank probe. 2-prime cross-check mandatory; saturation (order==Lmax+1) and
# prime mismatch are REFUSALS, not results.
# pf_rank.jl — holonomic-rank / Picard–Fuchs-order probe (mod p), pre-farming.
#
# THEORY.  For a Feynman LP polynomial G=U+F (or any GKZ integrand) with kinematic
# variable s, the d=2 maxcut period  I(s)=∮_γ G(s,x)^{-1} dx  over a closed cycle
# γ⊂ℂ^n\{G=0} satisfies a Picard–Fuchs ODE in s.  Affine Griffiths–Dwork mod p:
#   ∂_s^k(G^{-1}) = P_k/G^{k+1},   P_{k+1}=∂_sP_k·G−(k+1)·P_k·G_s,  P_0=1.
# Lift to common pole level L+1:  w_k = P_k·G^{L−k}.  Exact forms at level L+1 are
#   W = span{ ∂_iA·G − L·A·∂_iG : i=1..n, deg A ≤ L·degG − degG + 1 }
# (levels m<L collapse into m=L via A↦A·G^{L−m}).  Then
#   pf_order_in_s  =  rank([W ; w_0..w_L]) − rank(W)   over K=GF(p)(s),
# computed at a random s₀∈GF(p) (generic-rank; two-prime cross-checked).
# Full-GKZ holonomic rank = leviathan chi_regulated (F4 GB-staircase mod p).
# Reducibility: pf_order < holonomic_rank ⇒ maxcut generates a PROPER
# sub-D-module ⇒ full connection block-triangularizes.  Operator-LCLM factorization
# (e.g. ice-cone 2×2) needs the operator over GF(p)(s); not built here.
#
# The linear-in-s jet P_k = (−1)^k k!·G_s^k is WRONG for deg_s(G) ≥ 2
# (drained residuals, double covers w²=P(z,s²), …).  Hence two routes:
#   * pf_probe uses the linear shortcut when G is genuinely linear in svar
#     (checked exactly over ℚ);
#   * nonlinear G auto-routes to pf_probe_gen: the exact jet recursion
#     P_{k+1}=∂_s(P_k)·G−(k+1)·P_k·∂_s(G) carried in GF(p)[x,s], evaluated at a
#     random s0 AFTER differentiation;
#   * pf_probe_gen / pf_probe_esc are also exported directly (Oscar-poly
#     interface, no chi_regulated) for callers holding a polynomial, not a string.
# Controls 6/6 PASS incl. Legendre(s²)→2 and deg_s=4 const-period→1 (regression
# runner: regress_order.jl beside this file).
# Leviathan (chi_regulated engine): resolved from sibling tools/leviathan, then
# upgrades/leviathan in this repo. Override with DIPSTICK_LEVIATHAN.
let _lev = get(ENV, "DIPSTICK_LEVIATHAN", "")
    if isempty(_lev)
        for _c in (joinpath(@__DIR__, "..", "leviathan", "src", "Leviathan.jl"),
                   joinpath(@__DIR__, "..", "..", "upgrades", "leviathan", "src", "Leviathan.jl"))
            if isfile(_c)
                _lev = _c
                break
            end
        end
    end
    isempty(_lev) && error("Leviathan.jl not found — set DIPSTICK_LEVIATHAN to its path")
    include(_lev)
end
using .Leviathan, Oscar, Random
const Nemo = Oscar.Nemo

"Streaming sparse row-echelon over Z/pZ (Int64). Pivot = max column (degree-ordered)."
mutable struct _RE
    piv::Dict{Int,Tuple{Vector{Int},Vector{Int}}}   # pivot-col => (cols, vals)
    p::Int; rk::Int
end
_RE(p::Int) = _RE(Dict{Int,Tuple{Vector{Int},Vector{Int}}}(), p, 0)
@inline _mm(a,b,p) = Int(mod(widemul(a,b), p))
function _feed!(re::_RE, row::Dict{Int,Int})
    p = re.p
    while !isempty(row)
        pc = maximum(keys(row))
        haskey(re.piv, pc) || break
        α = p - row[pc]                       # subtract α·pivrow (pivrow[pc]==1)
        pcs, pvs = re.piv[pc]
        @inbounds for k in eachindex(pcs)
            c = pcs[k]; nv = mod(get(row,c,0) + _mm(α,pvs[k],p), p)
            nv == 0 ? delete!(row,c) : (row[c]=nv)
        end
    end
    isempty(row) && return false
    pc = maximum(keys(row)); ipv = invmod(row[pc], p)
    ks = sort!(collect(keys(row))); vv = [_mm(row[k],ipv,p) for k in ks]
    re.piv[pc] = (ks, vv); re.rk += 1; true
end

_row(f, idx, p) = begin
    cs=Int[]; vs=Int[]
    for (c,e) in zip(Nemo.coefficients(f), Nemo.exponent_vectors(f))
        push!(cs, idx[Tuple(e)]); push!(vs, Int(Nemo.lift(Nemo.ZZ,c)) % p)
    end
    (cs, vs)
end

"W-space of exact forms + w_k feed: the shared pf_probe core."
function _pf_core(G0, Rp, xv, wk, Lmax, prime; verbose=false)
    n = length(xv)
    dG  = [Oscar.derivative(G0, xv[i]) for i in 1:n]
    degG = Oscar.total_degree(G0); dW = Lmax*degG; dA = dW - degG + 1
    # monomial index for degree ≤ dW (degree-major so max-col = leading term)
    idx = Dict{NTuple{n,Int},Int}(); k=0
    for d in 0:dW, e in Nemo.exponent_vectors(d==0 ? one(Rp) : sum(xv)^d)
        idx[Tuple(e)] = (k+=1)
    end
    # precompute exponent/coeff lists of G and ∂_iG (Int64 mod p)
    EG  = collect(zip([Tuple(e) for e in Nemo.exponent_vectors(G0)],
                      [Int(Nemo.lift(Nemo.ZZ,c))%prime for c in Nemo.coefficients(G0)]))
    EdG = [collect(zip([Tuple(e) for e in Nemo.exponent_vectors(dG[i])],
                       [Int(Nemo.lift(Nemo.ZZ,c))%prime for c in Nemo.coefficients(dG[i])]))
           for i in 1:n]
    t0=time(); re = _RE(prime); nW=0
    # W-relations: ∂_iA·G − L·A·∂_iG for every monomial A=x^a, every i (row built
    # by direct exponent shifts — no Oscar poly arithmetic in the hot loop).
    for d in dA:-1:0, ae in Nemo.exponent_vectors(d==0 ? one(Rp) : sum(xv)^d)
        a = Tuple(ae)
        for i in 1:n
            row = Dict{Int,Int}()
            if a[i] > 0
                am = ntuple(j->a[j]-(j==i), n)
                for (g,c) in EG
                    row[idx[ntuple(j->am[j]+g[j],n)]] = _mm(a[i],c,prime)
                end
            end
            mL = prime - Lmax
            for (g,c) in EdG[i]
                key = idx[ntuple(j->a[j]+g[j],n)]
                nv = mod(get(row,key,0) + _mm(mL,c,prime), prime)
                nv==0 ? delete!(row,key) : (row[key]=nv)
            end
            isempty(row) && continue
            nW += 1; _feed!(re, row)
            (verbose && nW % 5000 == 0) && (println(stderr,
                "  W:",nW," rk=",re.rk," t=",round(time()-t0,digits=1),"s"); flush(stderr))
        end
    end
    jump = count(_feed!(re, Dict(zip(_row(w,idx,prime)...))) for w in wk)
    (jump, round(time()-t0,digits=2), nW, length(idx))
end

"""
    pf_probe_gen(GQ, R, xnames, kinnames, svar; Lmax, primes, seed)

PF order of ∮ GQ^{-1} dx in svar for GQ of ARBITRARY degree in svar.
GQ ∈ R = QQ[xnames ∪ kinnames]; non-svar kin → random Fp values. Exact jet
P_{k+1}=∂_s(P_k)·G−(k+1)·P_k·∂_s(G) carried in GF(p)[x,s], svar evaluated at a
random s0 AFTER differentiation. No chi_regulated (holonomic rank not computed).
"""
function pf_probe_gen(GQ, R, xnames, kinnames, svar; Lmax::Int=3,
                      primes=(Leviathan.P29[1], Leviathan.P29[2]),
                      seed::Int=20260707, verbose::Bool=false)
    rng = MersenneTwister(seed)
    n = length(xnames)
    ords = Int[]; tW=Float64[]; meta=[]
    for prime in primes
        Fp = Nemo.GF(prime)
        Rj, vj = Oscar.polynomial_ring(Fp, vcat(xnames, [svar*"__s"]))
        sv = vj[end]; xvj = vj[1:n]
        kv = Dict{String,Int}(k => rand(rng, 2:prime-2) for k in kinnames)
        s0 = rand(rng, 2:prime-2)
        img = elem_type(Rj)[]
        for g in Oscar.gens(R)
            vn = string(g)
            if vn == svar
                push!(img, sv)
            elseif vn in xnames
                push!(img, xvj[findfirst(==(vn), xnames)])
            elseif haskey(kv, vn)
                push!(img, Rj(Fp(kv[vn])))
            else
                push!(img, Rj(0))   # drained/absent var (must not occur in GQ)
            end
        end
        h = Oscar.hom(R, Rj, c->Fp(Nemo.numerator(c))*inv(Fp(Nemo.denominator(c))), img)
        Gj = h(GQ)
        @assert !iszero(Gj) "G vanished mod p=$prime"
        Gs = Oscar.derivative(Gj, sv)
        Pk = [one(Rj)]
        for k in 0:Lmax-1
            push!(Pk, Oscar.derivative(Pk[end], sv)*Gj - (k+1)*Pk[end]*Gs)
        end
        Rp, xv = Oscar.polynomial_ring(Fp, xnames)
        h2 = Oscar.hom(Rj, Rp, vcat(xv, [Rp(Fp(s0))]))
        G0 = h2(Gj)
        wk = [h2(Pk[k+1]) * G0^(Lmax-k) for k in 0:Lmax]
        jump, t, nW, ncol = _pf_core(G0, Rp, xv, wk, Lmax, prime; verbose=verbose)
        push!(ords, jump); push!(tW, t)
        push!(meta, (prime=prime, s0=s0, kv=kv, nW=nW, ncol=ncol,
                     degG=Oscar.total_degree(G0)))
    end
    (pf_order=ords[1], pf_order_2p=Tuple(ords), stable=(ords[1]==ords[2]),
     Lmax=Lmax, t_pf=tW, meta=meta)
end

"Escalation wrapper: distrust jump==Lmax+1 (saturation) or prime mismatch."
function pf_probe_esc(GQ, R, xnames, kinnames, svar; Ls=(3,5), kwargs...)
    local r
    for Lmax in Ls
        r = pf_probe_gen(GQ, R, xnames, kinnames, svar; Lmax=Lmax, kwargs...)
        r.stable && r.pf_order <= Lmax - 1 && return (r..., saturated=false)
    end
    return (r..., saturated=(r.pf_order >= r.Lmax))
end

"""
    pf_probe(Gstr, kinnames, xnames, svar; Lmax, p, p2, seed)

Return (holonomic_rank, pf_order, pf_order_2p, reducible, timings).
Gstr is the LP polynomial G=U+F as a string in kinnames∪xnames; svar is the
kinematic variable to take the PF operator in (others → random Fp values).
G linear in svar (Feynman-LP): the shortcut P_k = (−1)^k k!·G_s^k.
G nonlinear in svar: auto-routes to the exact jet (pf_probe_gen).
"""
function pf_probe(Gstr, kinnames, xnames, svar; Lmax::Int=6, verbose::Bool=false,
                  p::Int=Leviathan.P29[1], p2::Int=Leviathan.P29[2], seed::Int=20260630)
    rng = MersenneTwister(seed)
    fam = family(Gstr, kinnames, xnames)
    t_chi = @elapsed hr = chi_regulated(fam; p=Leviathan.P31[2], rng=MersenneTwister(seed))
    si  = findfirst(==(svar), kinnames); n = length(xnames)
    # dispatch: exact ℚ check for nonlinearity in svar (the shortcut would truncate)
    GssQ = Oscar.derivative(Oscar.derivative(fam.G, fam.kin[si]), fam.kin[si])
    if !iszero(GssQ)
        g = pf_probe_gen(fam.G, fam.R, xnames, kinnames, svar; Lmax=Lmax,
                         primes=(p, p2), seed=seed, verbose=verbose)
        return (holonomic_rank=hr, pf_order=g.pf_order, pf_order_2p=g.pf_order_2p,
                stable=g.stable, reducible=(g.pf_order < hr), Lmax=Lmax,
                t_chi=round(t_chi,digits=2), t_pf=round.(g.t_pf,digits=2),
                nW_ncol=(g.meta[1].nW, g.meta[1].ncol))
    end
    ords = Int[]; tW = Float64[]; nWnC = Tuple{Int,Int}[]
    for prime in (p, p2)
        Fp = Nemo.GF(prime)
        kv = random_kinvals(rng, kinnames, prime)
        Rp, xv = Oscar.polynomial_ring(Fp, xnames)
        sub = vcat([Rp(Fp(kv[k])) for k in kinnames], xv)
        h  = Oscar.hom(fam.R, Rp, c->Fp(Nemo.numerator(c))*inv(Fp(Nemo.denominator(c))), sub)
        G0  = h(fam.G)
        Gs  = h(Oscar.derivative(fam.G, fam.kin[si]))
        # w_k numerators (G linear in s ⇒ P_k = (-1)^k k! G_s^k; drop scalar):
        wk = [Gs^k * G0^(Lmax-k) for k in 0:Lmax]
        jump, t, nW, ncol = _pf_core(G0, Rp, xv, wk, Lmax, prime; verbose=verbose)
        push!(ords, jump); push!(tW, t); push!(nWnC, (nW, ncol))
    end
    (holonomic_rank=hr, pf_order=ords[1], pf_order_2p=Tuple(ords),
     stable=(ords[1]==ords[2]), reducible=(ords[1] < hr), Lmax=Lmax,
     t_chi=round(t_chi,digits=2), t_pf=round.(tW,digits=2), nW_ncol=nWnC[1])
end

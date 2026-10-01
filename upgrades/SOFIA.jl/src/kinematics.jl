# n-point kinematics: the cyclic Mandelstam basis and the substitution table
# p_i·p_j -> polynomial(basis invariants, external masses squared).
#
# Port of `GenerateKinematics` (scr.m:172-186). Differences: we work with
# SQUARED external masses MM_i as ring generators from the start (upstream
# carries M_i and squares/homogenizes later via `homogeneizeKin`), and the
# linear solve is hand-rolled Gaussian elimination over QQ with
# polynomial-valued right-hand sides (upstream uses `Solve`).

"""
    Kinematics

n-point kinematics data:
- `ring`: QQ[basis invariants..., MM_1..MM_n] (mpoly ring)
- `basis`: names of the cyclic-basis invariants (s12, s23, ..., s123, ...)
- `dots`: Dict (i,j) -> p_i·p_j as a ring element, 1 <= i <= j <= n
- `n`: number of external particles
"""
struct Kinematics
    ring::Any
    basis::Vector{String}
    masses::Vector{Any}    # MM_i generators
    dots::Dict{Tuple{Int,Int},Any}
    n::Int
end

dot(k::Kinematics, i::Int, j::Int) = k.dots[(min(i, j), max(i, j))]

"""Cyclic windows defining the Mandelstam basis (scr.m:177): for window
length l+1 with l = 1..floor(n/2-1), starting points a = 1..n, except the
longest windows for even n which take only a = 1..n/2 (avoiding the
momentum-conservation double counting)."""
function cyclic_basis_indices(n::Int)
    out = Vector{Vector{Int}}()
    lmax = floor(Int, n / 2 - 1 + eps())   # WL Table iterates l = 1 .. n/2-1
    for l in 1:lmax
        arange = (l == n ÷ 2 - 1 && iseven(n)) ? (1:(n ÷ 2)) : (1:n)
        # NOTE: upstream condition is l == n/2-1 (exact, so only for even n)
        for a in arange
            w = sort!(unique!([mod1(a + k, n) for k in 0:l]))
            push!(out, w)
        end
    end
    return out
end

sname(w::Vector{Int}) = "s" * join(string.(w))

"""
    generate_kinematics(n) -> Kinematics

Generic n-point kinematics. For n == 2 the only invariant is s = p1·p1
(scr.m:174). For n >= 3: build the p_i·p_j matrix in terms of two-index
Mandelstams and masses, impose momentum conservation Σ_j p_i·p_j = 0 and the
definitions of the multi-index basis invariants, and solve linearly for the
non-basis two-index invariants.
"""
function generate_kinematics(n::Int)
    if n == 2
        R, g = Nemo.polynomial_ring(Nemo.QQ, ["s", "MM1", "MM2"])
        dots = Dict((1, 1) => g[1], (1, 2) => g[1], (2, 2) => g[1])
        # p2 = -p1 by conservation: p1·p1 = s, p1·p2 = -s, p2·p2 = s.
        dots[(1, 2)] = -g[1]
        return Kinematics(R, ["s"], g[2:3], dots, 2)
    end
    windows = cyclic_basis_indices(n)
    basis = sname.(windows)
    pairs = [(i, j) for i in 1:n for j in (i+1):n]
    pairnames = ["s" * string(i) * string(j) for (i, j) in pairs]
    nonbasis_pairs = [(p, nm) for (p, nm) in zip(pairs, pairnames) if !(nm in basis)]
    massnames = ["MM" * string(i) for i in 1:n]
    R, gens_ = Nemo.polynomial_ring(Nemo.QQ, vcat(basis, massnames))
    bgen = Dict(zip(basis, gens_[1:length(basis)]))
    MM = gens_[length(basis)+1:end]

    # unknowns: non-basis two-index invariants s_{ij}
    unknowns = [p for (p, _) in nonbasis_pairs]
    uidx = Dict(p => k for (k, p) in enumerate(unknowns))
    nu = length(unknowns)

    # ppdot(i,j) as (constant-part poly, coefficient row over unknowns):
    # p_i·p_i = MM_i ; p_i·p_j = (s_ij - MM_i - MM_j)/2.
    half = Nemo.QQ(1, 2)
    function pp_entry(i, j)
        row = zeros(Nemo.QQFieldElem, nu)
        if i == j
            return (MM[i], row)
        end
        a, b = min(i, j), max(i, j)
        nm = "s" * string(a) * string(b)
        if haskey(bgen, nm)
            return (half * (bgen[nm] - MM[i] - MM[j]), row)
        else
            row[uidx[(a, b)]] = half
            return (half * (-MM[i] - MM[j]), row)
        end
    end

    # equations: A u + c = 0  with A over QQ, c in R
    A = Vector{Vector{Nemo.QQFieldElem}}()
    C = Vector{Any}()
    # momentum conservation: Σ_j p_i·p_j = 0
    for i in 1:n
        row = zeros(Nemo.QQFieldElem, nu)
        c = zero(R)
        for j in 1:n
            cj, rj = pp_entry(i, j)
            c += cj
            row .+= rj
        end
        push!(A, row); push!(C, c)
    end
    # multi-index basis definitions: s_W = Σ_{b,c in W} p_b·p_c
    for w in windows
        length(w) == 2 && continue
        row = zeros(Nemo.QQFieldElem, nu)
        c = -bgen[sname(w)]
        for b in w, cidx in w
            cj, rj = pp_entry(b, cidx)
            c += cj
            row .+= rj
        end
        push!(A, row); push!(C, c)
    end

    # Gaussian elimination: solve A u = -C exactly over QQ (poly-valued RHS)
    m = length(A)
    rhs = [-c for c in C]
    piv = zeros(Int, nu)
    r = 1
    for col in 1:nu
        pr = findfirst(k -> k >= r && !iszero(A[k][col]), 1:m)
        pr === nothing && continue
        A[r], A[pr] = A[pr], A[r]
        rhs[r], rhs[pr] = rhs[pr], rhs[r]
        pv = A[r][col]
        A[r] = [a / pv for a in A[r]]
        rhs[r] = rhs[r] * (1 // Rational{BigInt}(pv))
        for k in 1:m
            k == r && continue
            f = A[k][col]
            iszero(f) && continue
            A[k] = [a - f * b for (a, b) in zip(A[k], A[r])]
            rhs[k] = rhs[k] - rhs[r] * Rational{BigInt}(f)
        end
        piv[col] = r
        r += 1
        r > m && break
    end
    all(piv .> 0) || error("generate_kinematics: underdetermined system (n=$n)")
    usol = Vector{Any}(undef, nu)
    for col in 1:nu
        usol[col] = rhs[piv[col]]
    end

    dots = Dict{Tuple{Int,Int},Any}()
    for i in 1:n, j in i:n
        c, row = pp_entry(i, j)
        v = c
        for (k, coef) in enumerate(row)
            iszero(coef) && continue
            v += Rational{BigInt}(coef) * usol[k]
        end
        dots[(i, j)] = v
    end
    return Kinematics(R, basis, MM, dots, n)
end

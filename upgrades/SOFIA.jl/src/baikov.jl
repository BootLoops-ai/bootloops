# Loop-by-loop Baikov representation: Gram matrices in the propagator basis.
#
# Port of `LBL` / `SOFIABaikov` / `prepareLandauSystem` (scr.m:664-862, 1093).
# Architectural translation: upstream manipulates symbolic momenta through
# CenterDot pattern rules and compresses repeated combinations into Λ[a][i]
# symbols (a Mathematica-efficiency device). Here momenta are explicit
# rational vectors over the basis [l_1..l_L, p_1..p_{n-1}]; Gram matrices are
# built at the vector level and expanded bilinearly once, after the optional
# loop-momentum shift. Pivot choices in the dot->x change of variables only
# relabel integration variables and cannot affect the singularity set.

const RatVec = Vector{Rational{BigInt}}

"""
    BaikovSystem

Result of the LBL construction for one diagram:
- `ring`: QQ[x_1..x_N, kinematic basis, MM_i, internal mm's]
- `nx`: number of Baikov variables (propagators + ISPs)
- `nprop`: number of propagators (the first `nprop` x's vanish on max cut)
- `grams`: Gram matrices [ext_1, int_1, ext_2, int_2, ...] over `ring`
- `nloops`, `kin`: bookkeeping
"""
struct BaikovSystem
    ring::Any
    xgens::Vector{Any}
    nx::Int
    nprop::Int
    grams::Vector{Any}     # AbstractAlgebra matrices over ring
    nloops::Int
    kinnames::Vector{String}
end

# ---------- momentum assignment ----------

"""Solve momentum conservation for non-loop edge momenta. Returns
`nothing` if the loop-edge choice makes the system unsolvable (upstream FLAG).
Edge e is directed ends[1] -> ends[2]; external leg momentum p_node enters at
its vertex; p_n = -(p_1 + ... + p_{n-1})."""
function edge_momenta(d::Diagram, loopedges::Vector{Int})
    L = length(loopedges)
    next = length(d.nodes)
    dimv = L + (next - 1)
    vs = vertices(d)
    isloop = Dict(e => a for (a, e) in enumerate(loopedges))
    unknown = [i for i in 1:nedges(d) if !haskey(isloop, i)]
    uidx = Dict(e => k for (k, e) in enumerate(unknown))
    nu = length(unknown)
    # rows: one per vertex; cols: unknown q_e; rhs: known (loop + external) part
    rows = Vector{Vector{Rational{BigInt}}}()
    rhs = Vector{RatVec}()
    pvec(i) = begin                       # momentum of external leg i
        v = zeros(Rational{BigInt}, dimv)
        if i < next
            v[L+i] = 1
        else
            for k in 1:next-1
                v[L+k] = -1
            end
        end
        v
    end
    lvec(a) = (v = zeros(Rational{BigInt}, dimv); v[a] = 1; v)
    for vtx in vs
        row = zeros(Rational{BigInt}, nu)
        known = zeros(Rational{BigInt}, dimv)
        for (i, e) in enumerate(d.edges)
            s = (e.ends[1] == vtx ? -1 : 0) + (e.ends[2] == vtx ? 1 : 0)
            # self-loops contribute -1+1 = 0
            s == 0 && continue
            if haskey(isloop, i)
                known .+= s .* lvec(isloop[i])
            else
                row[uidx[i]] += s
            end
        end
        for (j, nd) in enumerate(d.nodes)
            nd.vertex == vtx && (known .+= pvec(j))
        end
        push!(rows, row)
        push!(rhs, -known)
    end
    # Gaussian elimination over Q with vector-valued RHS
    m = length(rows)
    piv = zeros(Int, nu)
    r = 1
    for col in 1:nu
        pr = findfirst(k -> k >= r && !iszero(rows[k][col]), 1:m)
        pr === nothing && continue
        rows[r], rows[pr] = rows[pr], rows[r]
        rhs[r], rhs[pr] = rhs[pr], rhs[r]
        pv = rows[r][col]
        rows[r] ./= pv
        rhs[r] ./= pv
        for k in 1:m
            k == r && continue
            f = rows[k][col]
            iszero(f) && continue
            rows[k] .-= f .* rows[r]
            rhs[k] .-= f .* rhs[r]
        end
        piv[col] = r
        r += 1
    end
    any(piv .== 0) && return nothing            # underdetermined: bad choice
    for k in 1:m   # consistency: zero rows must have zero RHS
        all(iszero, rows[k]) && !all(iszero, rhs[k]) && return nothing
    end
    qs = Vector{RatVec}(undef, nedges(d))
    for (i, e) in enumerate(d.edges)
        i in loopedges && (qs[i] = lvec(isloop[i]))
    end
    for (col, e) in enumerate(unknown)
        qs[e] = rhs[piv[col]]
    end
    return qs
end

# ---------- the LBL Gram construction (vector level) ----------

"""Per-loop momenta lists and Gram matrices of momentum VECTORS, mirroring
the main loop of `LBL` (scr.m:714-736). Returns (gramvecs, nlbl) where
gramvecs[a] = (ext_rows::Vector{RatVec}, int_rows::Vector{RatVec}) with the
convention ext rows are (v - l_a) for v in loopmomenta minus l_a itself."""
function lbl_gram_vectors(qs::Vector{RatVec}, L::Int)
    dimv = isempty(qs) ? 0 : length(qs[1])
    momenta = copy(qs)                      # edgesMomenta ∪ extraPropagators
    gramvecs = Vector{Tuple{Vector{RatVec},Vector{RatVec}}}()
    nlbl = L
    for a in 1:L
        sel = RatVec[]
        for v in momenta
            iszero(v[a]) && continue
            any(b -> b < a && !iszero(v[b]), 1:L) && continue
            w = v ./ v[a]                   # normalize coefficient of l_a to 1
            any(u -> u == w, sel) || push!(sel, w)
        end
        # nontrivial rests count toward the expected variable bound
        la = zeros(Rational{BigInt}, dimv); la[a] = 1
        nlbl += count(v -> v != la, sel)
        ext = [v .- la for v in sel if v != la]
        push!(gramvecs, (ext, sel))
        for i in eachindex(sel), j in (i+1):lastindex(sel)
            push!(momenta, sel[j] .- sel[i])
        end
    end
    return gramvecs, nlbl
end

"""Distinct l-involving dot products (as index pairs) appearing in the Gram
matrices built from `gramvecs`. Basis: index 1..L are loops, L+1..dim are
externals; a dot (i,j) with min(i,j) <= L involves a loop momentum."""
function distinct_ldots(gramvecs, L::Int)
    seen = Set{Tuple{Int,Int}}()
    function scan(rows::Vector{RatVec})
        for v in rows, w in rows
            for i in eachindex(v), j in eachindex(w)
                (iszero(v[i]) || iszero(w[j])) && continue
                a, b = min(i, j), max(i, j)
                a <= L && push!(seen, (a, b))
            end
        end
    end
    for (ext, int) in gramvecs
        scan(ext); scan(int)
    end
    return seen
end

"""Search loop-momentum shifts l_a -> l_a + δ_a minimizing the number of
distinct loop dot products (upstream `FindLoopMomShift`, scr.m:649, replaced
by direct optimization over the natural candidate shifts: the negated rests
of each loop's momenta). Returns shifted edge momenta."""
function minimize_dots_shift(qs::Vector{RatVec}, L::Int)
    gramvecs, nlbl = lbl_gram_vectors(qs, L)
    best = qs
    bestcount = length(distinct_ldots(gramvecs, L))
    bestcount <= nlbl && return qs          # already minimal
    dimv = length(qs[1])
    cands = Vector{Vector{RatVec}}()
    for a in 1:L
        la = zeros(Rational{BigInt}, dimv); la[a] = 1
        rests = unique([v .- la for v in gramvecs[a][2] if v != la])
        push!(cands, vcat([zeros(Rational{BigInt}, dimv)], [-r for r in rests]))
    end
    # bounded product search
    total = prod(length.(cands))
    total > 256 && (cands = [c[1:min(end, 3)] for c in cands])
    for combo in Iterators.product(cands...)
        # shifted momentum: q + Σ_a q[a]*δ_a  (δ_a may involve other l's)
        shifted = [q .+ reduce(.+, [q[a] .* combo[a] for a in 1:L]) for q in qs]
        gv, nb = lbl_gram_vectors(shifted, L)
        c = length(distinct_ldots(gv, L))
        if c < bestcount
            bestcount = c
            best = shifted
        end
    end
    return best
end

# ---------- dots -> Baikov variables ----------

"""
    lbl(d::Diagram; loopedges=nothing) -> BaikovSystem or nothing

Full LBL construction: momentum assignment, optional variable-minimizing
shift, bilinear Gram expansion, and the linear change of variables expressing
loop dot products through Baikov variables x_i (propagators first, ISPs
appended). Returns `nothing` if no loop-edge choice works.

`dotperm` selects among deterministic pivot orders for the dot -> Baikov
change of variables (0 = loop-major as upstream's Solve order approximates;
1 = reversed; 2/3 = rotations). Different pivot orders change which dot
products become ISPs — a route freedom that changes the CANDIDATE singularity
set on the maximal cut (both here and upstream); see `singularities_seed`'s
`routes` option for the union over routes.
"""
function lbl(d::Diagram; loopedges::Union{Nothing,Vector{Int}}=nothing,
             dotperm::Int=0, rowperm::Int=0)
    L = nloops(d)
    L == 0 && return nothing
    next = length(d.nodes)
    choices = loopedges === nothing ? loop_edge_candidates(d) : [loopedges]
    qs = nothing
    for ch in choices
        qs = edge_momenta(d, ch)
        qs !== nothing && break
    end
    qs === nothing && return nothing
    qs = minimize_dots_shift(qs, L)
    gramvecs, _ = lbl_gram_vectors(qs, L)

    # --- enumerate loop dots appearing anywhere (Grams + on-shell) ---
    dimv = length(qs[1])
    ldot_index = Dict{Tuple{Int,Int},Int}()   # (i,j) i<=j, i<=L  -> y index
    function regdot(i, j)
        a, b = min(i, j), max(i, j)
        a > L && return 0
        get!(ldot_index, (a, b)) do
            length(ldot_index) + 1
        end
    end
    function scanvec(v::RatVec, w::RatVec)
        for i in eachindex(v), j in eachindex(w)
            (iszero(v[i]) || iszero(w[j])) && continue
            regdot(i, j)
        end
    end
    for q in qs
        scanvec(q, q)
    end
    for (ext, int) in gramvecs
        for v in ext, w in ext
            scanvec(v, w)
        end
        for v in int, w in int
            scanvec(v, w)
        end
    end
    ndots = length(ldot_index)

    # --- kinematics & mass generators ---
    # Squared-mass names follow upstream homogeneizeKin: a label m1/m[1] names
    # the generator mm1, M2 names MM2, bare m names mm. Edge and node labels
    # share one namespace, so a node mass equal (by label) to an internal mass
    # becomes the SAME generator — upstream identifies these scales too.
    kin = generate_kinematics(next)
    massnames = String[]            # all squared-mass generator names, in order
    massof = Dict{Int,Any}()        # edge index -> mass name or nothing
    for (i, e) in enumerate(d.edges)
        if e.mass == BigInt(0)
            massof[i] = nothing
        else
            nm = sqmass_name(e.mass)
            nm in massnames || push!(massnames, nm)
            massof[i] = nm
        end
    end
    nodemass = Dict{Int,Any}()      # node j -> mass name or nothing
    for (j, nd) in enumerate(d.nodes)
        if nd.mass == BigInt(0)
            nodemass[j] = nothing
        else
            nm = sqmass_name(nd.mass)
            nm in massnames || push!(massnames, nm)
            nodemass[j] = nm
        end
    end

    # --- decompose dot(v,w) = Σ c_k y_k + κ(kinematics) ---
    # kinematic part lives in kin.ring; we delay mapping to the final ring.
    function dot_decomp(v::RatVec, w::RatVec)
        coeffs = zeros(Rational{BigInt}, ndots)
        kpart = zero(kin.ring)
        for i in eachindex(v), j in eachindex(w)
            (iszero(v[i]) || iszero(w[j])) && continue
            c = v[i] * w[j]
            if min(i, j) <= L
                coeffs[ldot_index[(min(i, j), max(i, j))]] += c
            else
                kpart += Nemo.QQ(c) * dot(kin, i - L, j - L)
            end
        end
        return coeffs, kpart
    end

    # --- on-shell conditions: dot(q_e,q_e) - mm_e = x_e ---
    P = nedges(d)
    eqs = Vector{Tuple{Vector{Rational{BigInt}},Any,Int}}()  # (coeffs, kpart, e)
    for (i, q) in enumerate(qs)
        c, kp = dot_decomp(q, q)
        push!(eqs, (c, kp, i))
    end
    # rowperm: equation-order route freedom (upstream: Solve picks pivots in
    # its own document order); 1 = reversed, k>1 = rotation by k-1
    if rowperm == 1
        reverse!(eqs)
    elseif rowperm > 1
        eqs = vcat(eqs[rowperm:end], eqs[1:rowperm-1])
    end

    # --- final ring ---
    # solve for pivot dots: row reduce [C][y] = [x_e - κ_e]; free dots -> ISPs
    C = [copy(e[1]) for e in eqs]
    m = length(C)
    pivcol = zeros(Int, ndots)        # dot k -> row solving it
    r = 1
    # pivot order: dots sorted by (loop index, partner) — mirrors upstream's
    # loop-by-loop Solve order
    dotorder = sort(collect(keys(ldot_index)); by=t -> (t[1], t[2]))
    if dotperm == 1
        reverse!(dotorder)
    elseif dotperm > 1
        dotorder = vcat(dotorder[dotperm:end], dotorder[1:dotperm-1])
    end
    # We do classical elimination but must track RHS symbolically:
    # RHS_e starts as (x_e - κ_e); represent as: xcoef[e][f] over x_1..x_P plus kpoly.
    xcoef = [zeros(Rational{BigInt}, P) for _ in 1:m]
    for e in 1:m
        xcoef[e][eqs[e][3]] = 1
    end
    kpoly = [-eqs[e][2] for e in 1:m]
    for (a, b) in dotorder
        col = ldot_index[(a, b)]
        pr = findfirst(k -> k >= r && k <= m && !iszero(C[k][col]), 1:m)
        pr === nothing && continue
        C[r], C[pr] = C[pr], C[r]
        xcoef[r], xcoef[pr] = xcoef[pr], xcoef[r]
        kpoly[r], kpoly[pr] = kpoly[pr], kpoly[r]
        pv = C[r][col]
        C[r] ./= pv
        xcoef[r] ./= pv
        kpoly[r] = kpoly[r] * Nemo.QQ(1 // pv)
        for k in 1:m
            k == r && continue
            f = C[k][col]
            iszero(f) && continue
            C[k] .-= f .* C[r]
            xcoef[k] .-= f .* xcoef[r]
            kpoly[k] = kpoly[k] - Nemo.QQ(f) * kpoly[r]
        end
        pivcol[col] = r
        r += 1
        r > m && break
    end
    free = [k for k in 1:ndots if pivcol[k] == 0]
    nx = P + length(free)
    freex = Dict(k => P + t for (t, k) in enumerate(free))

    xnames = ["x$i" for i in 1:nx]
    allnames = vcat(xnames, kin.basis, massnames)
    Rfin, gfin = Nemo.polynomial_ring(Nemo.QQ, allnames)
    gx = gfin[1:nx]
    nbasis = length(kin.basis)
    gkin = gfin[nx+1:nx+nbasis]
    gmm = Dict(zip(massnames, gfin[nx+nbasis+1:end]))

    # map kin.ring -> Rfin: basis generators by position; each generic
    # external mass MM_j becomes the node's labeled generator (or 0)
    nodeMMtarget = Vector{Any}(undef, next)
    for j in 1:next
        lbl_ = nodemass[j]
        nodeMMtarget[j] = lbl_ === nothing ? zero(Rfin) : gmm[lbl_]
    end
    kin_images = vcat(gkin, nodeMMtarget)
    tofin(p) = Nemo.evaluate(p, kin_images)

    # y_k (dot) as element of Rfin
    ysol = Vector{Any}(undef, ndots)
    for k in 1:ndots
        if pivcol[k] > 0
            rr = pivcol[k]
            v = tofin(kpoly[rr])
            for f in 1:P
                iszero(xcoef[rr][f]) && continue
                v += Nemo.QQ(xcoef[rr][f]) * (gx[f] + edge_mass_term(massof[f], gmm, Rfin))
            end
            # subtract contributions of remaining free dots in this row
            for k2 in free
                iszero(C[rr][k2]) && continue
                v -= Nemo.QQ(C[rr][k2]) * gx[freex[k2]]
            end
            ysol[k] = v
        else
            ysol[k] = gx[freex[k]]
        end
    end

    # --- assemble Gram matrices over Rfin ---
    function dot_fin(v::RatVec, w::RatVec)
        c, kp = dot_decomp(v, w)
        out = tofin(kp)
        for k in 1:ndots
            iszero(c[k]) && continue
            out += Nemo.QQ(c[k]) * ysol[k]
        end
        return out
    end
    grams = Any[]
    for (ext, int) in gramvecs
        if isempty(ext)
            push!(grams, Nemo.matrix(Rfin, 1, 1, [one(Rfin)]))
        else
            push!(grams, fix_rank_drop(Nemo.matrix(Rfin, length(ext), length(ext),
                                       [dot_fin(v, w) for v in ext for w in ext])))
        end
        push!(grams, fix_rank_drop(Nemo.matrix(Rfin, length(int), length(int),
                                   [dot_fin(v, w) for v in int for w in int])))
    end
    return BaikovSystem(Rfin, gx, nx, P, grams, L, allnames)
end

"""x_e in upstream is defined by  x_e = q_e² - m_e²  (on-shell condition);
the y-solve above set dot(q_e,q_e) = x_e + m_e²."""
edge_mass_term(nm, gmm, R) = nm === nothing ? zero(R) : gmm[nm]

"""Squared-mass generator name for a mass label, following upstream
homogeneizeKin: m1/m[1] -> "mm1", M2 -> "MM2", m -> "mm", M -> "MM";
other labels get their name with a leading "sq" marker."""
function sqmass_name(k)
    s = k isa WLCall ? string(k.head) * join(string.(k.args)) : string(k)
    startswith(s, "m") && return "m" * s
    startswith(s, "M") && return "M" * s
    return "sq" * s
end

"""
    fix_rank_drop(M) -> matrix

Port of `rankDropRule` (scr.m:626): a Gram matrix that is rank-deficient for
generic variable values (checked exactly at random rational points) has
identically vanishing determinant; excise a minimal set of row/column pairs
(chosen from the support of the generic nullspace) restoring full rank.
"""
function fix_rank_drop(M)
    R = Nemo.base_ring(M)
    n = Nemo.nrows(M)
    n <= 1 && return M
    nv = Nemo.nvars(R)
    state = UInt64(0x2545f4914f6cdd1d)
    function rndpt()
        [begin
             state = state * 6364136223846793005 + 1442695040888963407
             num = Int(mod(state >> 33, 9973)) + 1
             state = state * 6364136223846793005 + 1442695040888963407
             den = Int(mod(state >> 33, 97)) + 1
             Nemo.QQ(num, den)
         end for _ in 1:nv]
    end
    evalmat(pt) = Nemo.matrix(Nemo.QQ, n, n,
                              [Nemo.evaluate(M[i, j], pt) for i in 1:n for j in 1:n])
    A = evalmat(rndpt())
    r = Nemo.rank(A)
    r == n && return M
    # confirm at a second point (avoid unlucky evaluation)
    A2 = evalmat(rndpt())
    r = max(r, Nemo.rank(A2))
    r == n && return M
    k = n - r
    ns = Nemo.nullspace(A2)[2]
    support = sort!(unique!([i for j in 1:Nemo.ncols(ns) for i in 1:n
                             if !iszero(ns[i, j])]))
    for sub in sorted_subsets(support, k)
        keep = setdiff(1:n, sub)
        B = A2[keep, keep]
        if Nemo.rank(B) == n - k
            return M[keep, keep]
        end
    end
    return M   # give up; determinant stays zero and is filtered downstream
end

# ---------- the Landau system ----------

"""
    prepare_landau_system(d; maxcut=true) -> (gramList, active, sys)

Gram determinants (homogeneous in the kinematics by construction) with the
maximal-cut substitution x_1..x_P -> 0 applied if `maxcut`; `active` is the
list of x generator indices still appearing (the variables to eliminate).
Port of `prepareLandauSystem` (scr.m:1093).
"""
function prepare_landau_system(d::Diagram; maxcut::Bool=true,
                               loopedges=nothing, dotperm::Int=0,
                               rowperm::Int=0)
    sys = lbl(d; loopedges, dotperm, rowperm)
    sys === nothing && error("prepare_landau_system: LBL failed for diagram")
    R = sys.ring
    dets = [Nemo.det(g) for g in sys.grams]
    if maxcut
        images = Vector{Any}(Nemo.gens(R))
        for i in 1:sys.nprop
            images[i] = zero(R)
        end
        dets = [Nemo.evaluate(p, images) for p in dets]
    end
    gramlist = [p for p in dets if !is_constant(p)]
    active = sort!(unique!([vi for p in gramlist for vi in 1:sys.nx
                            if !free_of(p, vi)]))
    return gramlist, active, sys
end

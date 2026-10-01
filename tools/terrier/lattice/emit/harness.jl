# harness.jl — glue-list emission harness for the K3xK3 flux problem.
# Emits the trace-<=48 Gram/glue list G48 as JSONL, runs the completeness
# self-check (det-chain generation order) and the positive controls (I1-I4,
# the published Q=25 flux of arXiv:2010.10519 App. B, and the 44 in-budget
# glued pairs of bank/sound_minvecs.json).
# Do not edit entry_validate: s11.jl carries a byte-identical copy and diffs it
# against this file at every launch (verbatim_check).
# Pluggable predicate interface: define PREDICATE_HOOKS = (side_analysis=..,
# perp_root_check=.., congruence_diag=..) in Main BEFORE include-ing this file to run
# the emission engines (tau_extract/entry_validate/canonical_key/emit_*/completeness)
# under a different exact predicate chain; the defaults are SweepPredicates
# (predicates.jl). The I1-I4/AppB/SOUND44 controls are the K3xK3 reference
# instance and stay pinned to the default hooks. See README.md.
#
# LIST SPECIFICATION IMPLEMENTED (the "glue lemma"; the clause labels below are
# the ones used in code comments, assertion texts and emitted records):
#   entry = coherent-tuple class of (r, A, At, B, glue) under SIMULTANEOUS
#     (P,Q) in GL_r(Z)^2 acting on (A,B,glue-marking) and (At,B) at once
#     (the "coherent-tuple" equivalence; canonical_key below is its normal form);
#   (i)   A, At even integral, sig (p, r-p), p <= 3, r-p <= 19, 1 <= r <= 22,
#         sig(A) == sig(At) (B0: sign transport between the two sides);
#   (ii') B in Z^{rxr}, det B != 0, Bt := At^-1 B^T A in Z, AND M := B At^-1 in Z
#         (R-a, dual-integrality; here we generate B = M*At so (ii') holds by
#         construction and is re-verified per entry);
#   (iii) S0 := B*Bt R-diagonalizable, all eigenvalues > 0, every eigenspace
#         A-definite and every eigenspace of St0 := Bt*B At-definite, tr <= 48;
#   (iv)  glue: primitive embeddings A -> Gamma, At -> Gamma-tilde (Nikulin);
#         Eichler zone (B3(a): p<=2, r-p<=18, r+l(A_L)<=20) => exists + unique
#         (Nikulin 1.14.4); out-of-zone => per-entry CHECKED (B3(b));
#   det-chain: |detA*detAt| <= |chi(0)|; here |chi(0)| = det(M)^2*|detA*detAt|;
#   corner (R-c): (m,mt) = (0,0) or m*mt >= 1; budget (1/2)tr S0 + m*mt <= 24;
#   FS filter (R-b): per-entry full-stabilization reconstruction check =
#         kernel-definiteness (p = 3 or r-p = 19, complement K definite) AND
#         R (no -2 root perp Sigma);
#   B2 (list coverage): the tau-image (tau_extract) of every admissible flux
#         satisfies (ii') and is a list entry; imprimitive alignments g(Gt) = cL,
#         c >= 2, are included (rank-1 rows with c >= 2);
#   B4: corner-only fluxes (N = 0) are never list entries (tau(0) undefined).
#
# COVERAGE EMITTED (recorded in the JSONL meta line):
#   COMPLETE strata: rank 1 (closed form, I2 pin 82); rank 2 definite (2,0) and
#   (0,2) (SL2 transversal incl imprimitive x Fincke-Pohst on kron(A,At) <= 48,
#   coherent-tuple dedup by canonical key). NOT ENUMERATED here (engine validated
#   via tau_extract on explicit fluxes): rank 2 sig (1,1) (Pell/unit reduction
#   needed; see s11.jl) and ranks 3..22.
#
# Run:  nice -n 5 julia +1.10 -t 1 --project=<your Oscar env> \
#         harness.jl all         (modes: emit | controls | all)
# Writes G48_LIST.jsonl (emit/all) and controls_report.json next to this file.

import JSON
using Oscar

const GLUEDIR = @__DIR__
const BANK = joinpath(GLUEDIR, "bank")           # shipped fixtures (sha256-pinned)
include(joinpath(GLUEDIR, "..", "lattice.jl"))   # LatticeLayer: BQF / theta / genus toolkit
using .LatticeLayer
# --- pluggable predicate interface (see header and README.md) ----------------
if !@isdefined(PREDICATE_HOOKS)
    include(joinpath(GLUEDIR, "predicates.jl"))
    PREDICATE_HOOKS = (side_analysis = SweepPredicates.side_analysis,
                       perp_root_check = SweepPredicates.perp_root_check,
                       congruence_diag = SweepPredicates.congruence_diag)
end
const side_analysis   = PREDICATE_HOOKS.side_analysis
const perp_root_check = PREDICATE_HOOKS.perp_root_check
const congruence_diag = PREDICATE_HOOKS.congruence_diag

const TRACE_BOUND = 48
const NCAP = 24

# --- Gamma_{3,19} Gram, conventions pinned:
# moduli-side-POSITIVE, E8 Cartan NEGATIVE, Bourbaki chain 1..7 node 8 on node 5.
const UMAT = [0 1; 1 0]
function cartanA(n); C = zeros(Int, n, n); for i in 1:n; C[i,i] = 2; end
    for i in 1:n-1; C[i,i+1] = C[i+1,i] = -1; end; C end
function cartanE8(); C = cartanA(8); C[7,8] = C[8,7] = 0; C[5,8] = C[8,5] = -1; C end
function blockdiag(ms...)
    n = sum(size(m, 1) for m in ms); out = zeros(Int, n, n); k = 0
    for m in ms; s = size(m, 1); out[k+1:k+s, k+1:k+s] = m; k += s; end
    out
end
k3gram() = blockdiag(UMAT, UMAT, UMAT, -cartanE8(), -cartanE8())

mzz(m::AbstractMatrix) = matrix(ZZ, size(m,1), size(m,2), Int.(vec(permutedims(m))))
tomat(z::ZZMatrix) = [[Int(z[i,j]) for j in 1:ncols(z)] for i in 1:nrows(z)]
is_intmat(q::QQMatrix) = all(isone(denominator(q[i,j])) for i in 1:nrows(q), j in 1:ncols(q))
zz_of(q::QQMatrix) = matrix(ZZ, [numerator(q[i,j]) for i in 1:nrows(q), j in 1:ncols(q)])

"signature (p, n) of a nondegenerate symmetric ZZ matrix, exact."
function sig_of(G::ZZMatrix)
    dg = congruence_diag(matrix(QQ, G))
    @assert !any(iszero, dg) "degenerate Gram in sig_of"
    (count(t -> t > 0, dg), count(t -> t < 0, dg))
end

"l(A_L) = minimal generator count of Z^r / G Z^r (nontrivial SNF divisors)."
function ell_of(G::ZZMatrix)
    s = snf(G)
    count(i -> !isone(abs(s[i,i])), 1:nrows(s))
end
disc_divisors(G::ZZMatrix) = [Int(abs(snf(G)[i,i])) for i in 1:nrows(G) if !isone(abs(snf(G)[i,i]))]

# -----------------------------------------------------------------------------
# tau_extract: N (n1 x n2 integer flux) w.r.t. Grams d1, d2 -> coherent tuple
# data (r, A, At, B, Bt, M, S0) of the list specification. Support lattices are
# the saturations of colspace(N) / colspace(N^T) (Nemo ZZ-kernel is saturated).
# Valid on any primitive block (T1,T2) of Gamma x Gamma-tilde: only the block
# Gram enters (support lattices lie in the block).
# -----------------------------------------------------------------------------
function sat_colspace(Nz::ZZMatrix)
    K = kernel(transpose(Nz), side=:right)
    ncols(K) == 0 && return identity_matrix(ZZ, nrows(Nz))
    kernel(transpose(K), side=:right)          # saturated basis, columns
end

function solve_exact(Aq::QQMatrix, rhs::QQMatrix)
    X = inv(transpose(Aq) * Aq) * transpose(Aq) * rhs   # full column rank
    @assert Aq * X == rhs "solve_exact: inconsistent system"
    X
end

function tau_extract(N::Matrix{Int}, d1::Matrix{Int}, d2::Matrix{Int})
    Nz, G1, G2 = mzz(N), mzz(d1), mzz(d2)
    @assert !iszero(Nz) "tau(0) undefined: corner-only fluxes are never list entries (B4)"
    Lb  = sat_colspace(Nz)                     # L-basis in side-1 coords
    Ltb = sat_colspace(transpose(Nz))          # Lt-basis in side-2 coords
    r = ncols(Lb)
    @assert ncols(Ltb) == r "rank mismatch L vs Lt"
    A  = transpose(Lb) * G1 * Lb
    At = transpose(Ltb) * G2 * Ltb
    rhs = matrix(QQ, Nz * G2 * Ltb)            # g|_{Lt} in ambient side-1 coords
    Bq = solve_exact(matrix(QQ, Lb), rhs)
    @assert is_intmat(Bq) "B not integral: support saturation broken"
    B = zz_of(Bq)
    Btq = inv(matrix(QQ, At)) * transpose(matrix(QQ, B)) * matrix(QQ, A)
    Mq  = matrix(QQ, B) * inv(matrix(QQ, At))
    dual_ok = is_intmat(Btq) && is_intmat(Mq)  # (ii') R-a — theorem on tau-images
    Bt = dual_ok ? zz_of(Btq) : nothing
    M  = dual_ok ? zz_of(Mq) : nothing
    S0 = dual_ok ? B * Bt : nothing
    # cross-check: S0 == restriction of S = N d2 N^T d1 to L (S Lb == Lb S0)
    if dual_ok
        Sfull = Nz * G2 * transpose(Nz) * G1
        @assert Sfull * Lb == Lb * S0 "S restriction mismatch"
    end
    (r=r, A=A, At=At, B=B, Bt=Bt, M=M, S0=S0, Lb=Lb, Ltb=Ltb, dual_ok=dual_ok)
end

# -----------------------------------------------------------------------------
# entry_validate: the repaired LIST predicate on a tuple (A, At, M).
# Returns (ok, reasons, data) — data carries everything the JSONL entry needs.
# -----------------------------------------------------------------------------
function entry_validate(A::ZZMatrix, At::ZZMatrix, M::ZZMatrix;
                        trace_bound::Int=TRACE_BOUND)
    reasons = String[]
    r = nrows(A)
    (1 <= r <= 22) || push!(reasons, "rank $r outside 1..22")
    all(iseven(A[i,i]) for i in 1:r) && all(iseven(At[i,i]) for i in 1:r) ||
        push!(reasons, "Gram not even")
    A == transpose(A) && At == transpose(At) || push!(reasons, "Gram not symmetric")
    detA, detAt, detM = det(A), det(At), det(M)
    (iszero(detA) || iszero(detAt)) && push!(reasons, "degenerate Gram")
    iszero(detM) && push!(reasons, "det B = 0 ((ii') fails)")
    !isempty(reasons) && return (ok=false, reasons=reasons, data=nothing)
    p, q = sig_of(A); pt, qt = sig_of(At)
    (p, q) == (pt, qt) || push!(reasons, "sig(A)=($p,$q) != sig(At)=($pt,$qt) (B0)")
    p <= 3 || push!(reasons, "p=$p > 3")
    q <= 19 || push!(reasons, "r-p=$q > 19")
    B  = M * At                                # (ii') by construction
    Bt = transpose(M) * A
    @assert matrix(QQ, Bt) == inv(matrix(QQ, At)) * transpose(matrix(QQ, B)) * matrix(QQ, A)
    S0, St0 = B * Bt, Bt * B
    t = Int(tr(S0))
    t <= trace_bound || push!(reasons, "tr S0 = $t > $trace_bound (budget)")
    iseven(t) || push!(reasons, "odd trace on even lattice (evenness cert broken)")
    chi0 = abs(detM)^2 * abs(detA * detAt)     # |chi(0)| = det(M)^2 |detA detAt|
    @assert abs(det(S0)) == chi0 "det-chain identity broken"
    abs(detA * detAt) <= chi0 || push!(reasons, "det-chain violated")
    sa  = side_analysis(matrix(QQ, S0),  matrix(QQ, A))
    sat = side_analysis(matrix(QQ, St0), matrix(QQ, At))
    sa.ok  || push!(reasons, "side A: " * join(sa.reasons, "; "))
    sat.ok || push!(reasons, "side At: " * join(sat.reasons, "; "))
    Rx, x = polynomial_ring(QQ, "x")
    cp = charpoly(Rx, matrix(QQ, S0))
    iszero(coeff(cp, 0)) && push!(reasons, "zero eigenvalue: (iii) needs all lambda > 0")
    # (iii) lambda > 0: side_analysis rejects negatives; zero handled above
    data = (r=r, p=p, q=q, B=B, Bt=Bt, S0=S0, t=t, cp=string(cp),
            detA=Int(detA), detAt=Int(detAt), detM=Int(detM), chi0=Int(chi0),
            nflux=t ÷ 2, corner_max=NCAP - t ÷ 2)
    (ok=isempty(reasons), reasons=reasons, data=data)
end

# -----------------------------------------------------------------------------
# Coherent-tuple canonical key for rank <= 2 definite entries.
# Key = lex-min over all isometries of A (resp. At) onto SL2-transversal Grams
# (both GL-orientations included => full GL_r(Z)^2 simultaneous action).
# -----------------------------------------------------------------------------
"All U in GL_2(Z) with U^T G1 U == G2, both pos-def even 2x2 (exact, brute)."
function isometries_2x2(G1::ZZMatrix, G2::ZZMatrix)
    L1 = integer_lattice(gram=G1)
    nb = max(Int(G2[1,1]), Int(G2[2,2]))
    cands = Vector{ZZMatrix}[]
    vs = Dict{Int,Vector{Vector{ZZRingElem}}}()
    for (v, n) in short_vectors(L1, nb)
        ni = Int(ZZ(n))
        push!(get!(vs, ni, Vector{ZZRingElem}[]), v, [-x for x in v])
    end
    out = ZZMatrix[]
    for v1 in get(vs, Int(G2[1,1]), Vector{ZZRingElem}[]),
        v2 in get(vs, Int(G2[2,2]), Vector{ZZRingElem}[])
        U = matrix(ZZ, 2, 2, [v1[1], v2[1], v1[2], v2[2]])
        abs(det(U)) == 1 || continue
        transpose(U) * G1 * U == G2 && push!(out, U)
    end
    out
end

const ISO_CACHE = Dict{Tuple{Vector{Int},Int},Vector{Tuple{ZZMatrix,ZZMatrix}}}()
"(target reduced Gram, U) pairs for pos-def 2x2 G: U^T G U = target."
function reduced_targets(G::ZZMatrix)
    keyG = ([Int(G[i,j]) for i in 1:2 for j in 1:2], Int(det(G)))
    get!(ISO_CACHE, keyG) do
        D = Int(det(G))
        out = Tuple{ZZMatrix,ZZMatrix}[]
        for (a, b, c) in reduced_binary_forms(-D)
            T = form_gram(a, b, c)
            for U in isometries_2x2(G, T)
                push!(out, (T, U))
            end
        end
        @assert !isempty(out) "no reduced target for Gram $G (transversal broken)"
        out
    end
end

flatkey(s, GA, GT, M) = string(s, "|", tomat(GA), "|", tomat(GT), "|", tomat(M))

function canonical_key(A::ZZMatrix, At::ZZMatrix, M::ZZMatrix, p::Int, q::Int)
    r = nrows(A)
    if r == 1
        n, nt, c = Int(A[1,1]), Int(At[1,1]), abs(Int(M[1,1]))
        return "r1|$n|$nt|$c"
    end
    @assert r == 2 && (p, q) in [(2,0), (0,2)] "canonical_key: definite rank-2 only"
    s = p == 2 ? 1 : -1
    Ap, Atp = s * A, s * At
    best = nothing
    for (TA, UA) in reduced_targets(Ap), (TT, UT) in reduced_targets(Atp)
        Mp = inv(matrix(QQ, UA)) * matrix(QQ, M) * transpose(inv(matrix(QQ, UT)))
        @assert is_intmat(Mp)
        k = flatkey(s, TA, TT, zz_of(Mp))
        (best === nothing || k < best) && (best = k)
    end
    best
end

# --- glue certificate (iv) + (B3): Eichler zone closed-form for rank <= 2 ---
function glue_cert(A::ZZMatrix, At::ZZMatrix, p::Int, q::Int)
    r = nrows(A)
    ellA, ellAt = ell_of(A), ell_of(At)
    zone = p <= 2 && q <= 18 && (r + ellA) <= 20 && (r + ellAt) <= 20
    Dict("zone_eichler" => zone,
         "embed_classes" => zone ? 1 : -1,   # -1 = per-entry CHECKED branch needed
         "cert" => zone ? "Nikulin 1.14.4 (B3(a) zone: p<=2, r-p<=18, r+l<=20): exists, unique up to O(Gamma), both sides" :
                          "OUT-OF-ZONE: per-entry primitive_embeddings check required (B3(b))",
         "ellA" => ellA, "ellAt" => ellAt,
         "discA" => disc_divisors(A), "discAt" => disc_divisors(At))
end

function fs_structural(p::Int, q::Int)
    kp, kq = 3 - p, 19 - q
    if p == 3 || q == 19
        Dict("pass" => "NEEDS-RECONSTRUCTION",
             "reason" => "kernel-definiteness OK structurally (sig K=($kp,$kq)); R-check needs glue reconstruction")
    else
        Dict("pass" => false,
             "reason" => "kernel-definiteness fails: sig(K) = ($kp,$kq) indefinite (R-b)")
    end
end

# --- deterministic JSON (fixed key order; JSON.json for leaf values) ---------
json_obj(pairs::Vector{<:Pair}) =
    "{" * join(["\"$(k)\":" * (v isa Vector{<:Pair} ? json_obj(v) : JSON.json(v))
                for (k, v) in pairs], ",") * "}"

function entry_record(stratum::String, A, At, M, vd, key)
    gl = glue_cert(A, At, vd.p, vd.q)
    Pair{String,Any}["rank" => vd.r, "stratum" => stratum,
     "sig" => [vd.p, vd.q],
     "A" => tomat(A), "At" => tomat(At), "M" => tomat(M),
     "B" => tomat(vd.B), "Bt" => tomat(vd.Bt), "S0" => tomat(vd.S0),
     "charpoly" => vd.cp, "trace" => vd.t,
     "detA" => vd.detA, "detAt" => vd.detAt, "detM" => vd.detM, "chi0" => vd.chi0,
     "det_chain" => Dict("lhs_detAdetAt" => abs(vd.detA * vd.detAt),
                         "rhs_chi0" => vd.chi0,
                         "ok" => abs(vd.detA * vd.detAt) <= vd.chi0),
     "nflux_mid" => vd.nflux, "corner_max_mmt" => vd.corner_max,
     "corner_rule" => "(m,mt)=(0,0) or 1<=m*mt<=corner_max_mmt (R-c; budget (1/2)tr+m*mt<=24)",
     "glue" => gl, "fs_check" => fs_structural(vd.p, vd.q), "key" => key]
end

# --- rank-1 stratum: closed-form complete (I2: 82 ordered entries) -----------
function emit_rank1()
    entries = Dict{String,Vector{Pair{String,Any}}}()
    breakdown = Dict{Int,Int}()
    for sgn in (1, -1), n in 1:12, nt in 1:12, c in 1:3
        n * nt * c^2 <= 12 || continue
        A, At, M = matrix(ZZ,1,1,[2*sgn*n]), matrix(ZZ,1,1,[2*sgn*nt]), matrix(ZZ,1,1,[c])
        v = entry_validate(A, At, M)
        @assert v.ok "rank-1 generation bug: $(v.reasons) at ($sgn,$n,$nt,$c)"
        key = canonical_key(A, At, M, v.data.p, v.data.q)
        @assert !haskey(entries, key)
        entries[key] = entry_record("r1", A, At, M, v.data, key)
        sgn == 1 && (breakdown[c] = get(breakdown, c, 0) + 1)
    end
    (entries=entries, breakdown=breakdown)
end

# --- rank-2 definite strata (2,0) and (0,2): complete ------------------------
# Generation order (det-chain): lambda1*lambda2 = |chi(0)| =
# det(M)^2 * DA * DAt and lambda1+lambda2 <= 48 => DA*DAt <= 576 (AM-GM,
# det M^2 >= 1). A, At over the FULL SL2 transversal (imprimitive included,
# both orientations) => every GL-class hit; M via Fincke-Pohst on the pos-def
# rank-4 form vec(M)^T kron(A,At) vec(M) = tr(M At M^T A) <= 48; coherent
# dedup by canonical_key. det M = 0 candidates are (ii')-rejected.
function emit_rank2_definite(s::Int)
    entries = Dict{String,Vector{Pair{String,Any}}}()
    stats = Dict{String,Int}("candidates" => 0, "detM0" => 0, "dups" => 0)
    dets = [D for D in 3:192 if mod(D, 4) in (0, 3)]
    for DA in dets, DAt in dets
        DA * DAt <= 576 || continue
        for f1 in reduced_binary_forms(-DA), f2 in reduced_binary_forms(-DAt)
            Gp1, Gp2 = form_gram(f1...), form_gram(f2...)
            A, At = s * Gp1, s * Gp2
            Lk = integer_lattice(gram=kronecker_product(Gp1, Gp2))
            for (v, nrm) in short_vectors(Lk, TRACE_BOUND)
                stats["candidates"] += 1
                M = matrix(ZZ, 2, 2, [v[1], v[2], v[3], v[4]])
                if iszero(det(M)); stats["detM0"] += 1; continue; end
                val = entry_validate(A, At, M)
                @assert val.ok "rank-2 generation bug: $(val.reasons)"
                @assert val.data.t == Int(ZZ(nrm)) "trace/norm route mismatch"
                key = canonical_key(A, At, M, val.data.p, val.data.q)
                if haskey(entries, key); stats["dups"] += 1; continue; end
                entries[key] = entry_record(s == 1 ? "r2pp" : "r2nn", A, At, M, val.data, key)
            end
        end
    end
    (entries=entries, stats=stats)
end

# -----------------------------------------------------------------------------
# Completeness self-check (det-chain argument):
# per entry |detA*detAt| <= |chi(0)| = det(M)^2|detA*detAt| <= (tr/2)^2 <= 576
# at r = 2 (AM-GM), so the generation caps DA*DAt <= 576, tr <= 48 exhaust the
# domain; rank-1 closed-form count and pinned I2 breakdown; plus an independent
# brute-box M recount (exact dual bound, LatticeLayer brute_force_theta) on
# the small-det sample.
# -----------------------------------------------------------------------------
function completeness_check(r1, r2pp, r2nn)
    rep = Dict{String,Any}()
    n1 = length(r1.entries)
    @assert n1 == 82 "rank-1 count $(n1) != pinned 82 (I2)"
    @assert r1.breakdown == Dict(1 => 35, 2 => 5, 3 => 1) "I2 c-breakdown broken"
    rep["rank1"] = Dict("count" => n1, "positive_side_c_breakdown" => r1.breakdown,
                        "pin" => "I2: 82 ordered; c=1:35, c=2:5, c=3:1 (x2 signs)")
    for (tag, res) in [("r2pp", r2pp), ("r2nn", r2nn)]
        for (_, e) in res.entries
            d = Dict(e)
            @assert abs(d["detA"] * d["detAt"]) <= d["chi0"] "det-chain fails"
            @assert d["chi0"] <= (d["trace"] ÷ 2)^2 "AM-GM cap fails"
            @assert abs(d["detA"] * d["detAt"]) <= 576 "det cap exceeded: generation hole"
            @assert mod(abs(d["detA"]), 4) in (0, 3) "even binary det law broken"
        end
        rep[tag] = Dict("count" => length(res.entries), "stats" => res.stats)
    end
    # independent M-recount, small-det sample: brute box vs Fincke-Pohst
    sample = Dict{String,Any}(); nchecked = 0
    for DA in (3, 4), DAt in (3, 4, 7, 8, 12)
        for f1 in reduced_binary_forms(-DA), f2 in reduced_binary_forms(-DAt)
            G = kronecker_product(form_gram(f1...), form_gram(f2...))
            sv = 2 * length(short_vectors(integer_lattice(gram=G), TRACE_BOUND))
            bf = sum(v for (n, v) in brute_force_theta(G, TRACE_BOUND) if n > 0)
            @assert sv == bf "FP vs brute recount mismatch on ($f1,$f2): $sv vs $bf"
            nchecked += 1
        end
    end
    sample["pairs_cross_counted"] = nchecked
    sample["law"] = "2x short_vectors == brute box (exact dual bound), norms in (0,48]"
    rep["m_enumeration_crossroute"] = sample
    rep["verdict"] = "COMPLETE on emitted strata {r1, r2pp, r2nn} by det-chain generation order"
    rep
end

# -----------------------------------------------------------------------------
# Controls. Each returns a Dict merged into the controls report; any failure
# throws (the harness must not emit a list whose controls fail).
# -----------------------------------------------------------------------------
function ctrl_I1(r2pp)
    grams = [(2,0,2),(2,1,2),(2,0,4),(2,1,4),(2,0,6),(2,1,6),(4,0,4),(4,1,4),(4,2,4)]
    LK3 = k3_lattice()
    emb = Dict{String,Any}()
    for (a, b, c) in grams
        G = matrix(ZZ, 2, 2, [a, b, b, c])
        M = identity_matrix(ZZ, 2)
        v = entry_validate(G, G, M)
        @assert v.ok "I1 Gram ($a,$b,$c) diagonal entry invalid: $(v.reasons)"
        key = canonical_key(G, G, M, 2, 0)
        @assert haskey(r2pp.entries, key) "I1 class ($a,$b,$c) MISSING from r2pp"
        ok, cls = primitive_embeddings(LK3, integer_lattice(gram=G))
        @assert ok && length(cls) == 1 "I1 embedding pin (true,1) broken on ($a,$b,$c)"
        emb["($a,$b,$c)"] = "in-list; primitive_embeddings=(true,1)"
    end
    @assert !is_isometric(integer_lattice(gram=matrix(ZZ,2,2,[2,0,0,6])),
                          integer_lattice(gram=matrix(ZZ,2,2,[4,2,2,4]))) "genus-vs-class trap"
    Dict("status" => "PASS", "classes" => emb,
         "trap_guard" => "(2,0,6) !~ (4,2,4) (is_isometric false)")
end

ctrl_I2(r1) = Dict("status" => "PASS",
    "count" => length(r1.entries), "pin" => 82,
    "c_breakdown_positive" => r1.breakdown,
    "note" => "c>=2 entries = non-saturated alignments g(Gt)=cL (B2 imprimitivity test mass)")

function ctrl_I3(i4entry)
    # (a) corner-only refusal (B4: N=0 never a list entry)
    rejected = false
    try; tau_extract(zeros(Int, 22, 22), k3gram(), k3gram()); catch; rejected = true; end
    @assert rejected "corner-only flux was not refused"
    # (b) M-rider arithmetic on the I4 entry: N_flux 12 -> 13 at (m,mt)=(1,1)
    d = Dict(i4entry)
    @assert d["nflux_mid"] == 12 && d["corner_max_mmt"] == 12
    @assert d["nflux_mid"] + 1 * 1 == 13 "I3b rider pin broken"
    # (c) attractive-pair glue roundtrip pins
    T1 = integer_lattice(gram=matrix(ZZ,2,2,[2,1,1,2]))
    T2 = integer_lattice(gram=matrix(ZZ,2,2,[2,0,0,2]))
    nsv = length(short_vectors(integer_lattice(
        gram=kronecker_product(gram_matrix(T1), gram_matrix(T2))), 48))
    @assert nsv == 492 "I3c short-vector pin 492 broken (got $nsv)"
    ovs = even_overlattices(direct_sum(T1, T2)[1])
    @assert length(ovs) == 1 "I3c overlattice count != 1"
    gsym = genus_str(ovs[1])
    @assert occursin("II_(4, 0)", gsym) && occursin("3^-1", gsym) "I3c genus pin broken: $gsym"
    Dict("status" => "PASS", "corner_only" => "REJECT (structural)",
         "rider_13" => true, "short_vectors_492" => nsv, "overlattice_genus" => gsym)
end

function ctrl_I4()
    U = [0 1; 1 0]
    B1, B2, B3 = [2 1; 1 1], [3 1; 1 1], [3 2; 1 1]
    N = zeros(Int, 22, 22)
    N[1:2, 1:2] = B1 * U; N[3:4, 3:4] = B2 * U; N[5:6, 5:6] = B3 * U
    dK3 = k3gram()
    tx = tau_extract(N, dK3, dK3)
    @assert tx.r == 6 && tx.dual_ok "I4 support/dual-integrality broken"
    v = entry_validate(tx.A, tx.At, tx.M)
    @assert v.ok "I4 list-predicate rejected: $(v.reasons)"
    @assert v.data.t == 24 && v.data.nflux == 12 "I4 trace/N_flux pins broken"
    Rx, x = polynomial_ring(QQ, "x")
    cpin = (x^2 - 6x + 1) * (x^2 - 8x + 4) * (x^2 - 10x + 1)
    @assert charpoly(Rx, matrix(QQ, v.data.S0)) == cpin "I4 charpoly pin broken"
    @assert (v.data.p, v.data.q) == (3, 3) "I4 sig(L) != (3,3)"
    fs = fs_structural(v.data.p, v.data.q)
    @assert fs["pass"] == "NEEDS-RECONSTRUCTION" "I4 kernel-definiteness branch wrong"
    # full-Gamma FS reconstruction: 5.1(iii) PASSES, R FIRES with 480 roots
    dz = mzz(dK3); Nz = mzz(N)
    S  = matrix(QQ, Nz * dz * transpose(Nz) * dz)
    sa = side_analysis(S, matrix(QQ, dz))
    @assert sa.ok "I4 full-Gamma 5.1(iii) should PASS (ker = 2E8(-1)^2 definite)"
    pr = perp_root_check(sa.factors, dz, [-2])
    @assert pr.status == :root_found && pr.roots_found == 480 "I4 R-check pin 480 broken"
    key = canonical_key_targeted(tx)
    entry = entry_record("targeted-r6-I4", tx.A, tx.At, tx.M, v.data, key)
    push!(entry, "enumerative" => false)
    push!(entry, "fs_reconstruction" => Dict(
        "kernel_definiteness" => "PASS (p=3, K = 2E8(-1) definite)",
        "R_check" => "FAIL: 480 norm -2 roots perp Sigma (perp_root_check)",
        "verdict" => "FS-REJECT via R, NOT via definiteness (I4 pin)"))
    (report=Dict("status" => "PASS", "r" => 6, "trace" => 24, "nflux" => 12,
                 "roots_perp_Sigma" => pr.roots_found,
                 "caught_by" => "R-check (P7); the eigenspace-definiteness criterion passes"),
     entry=entry)
end

canonical_key_targeted(tx) = "targeted|" * string(tomat(tx.A)) * "|" *
                             string(tomat(tx.At)) * "|" * string(tomat(tx.M))

function ctrl_APPB()
    t0 = time()
    J = JSON.parsefile(joinpath(BANK,
                                "k3xk3_flux_N_Q25_arxiv2010.10519_appB.json"))
    N = [Int(J["N"][i][j]) for i in 1:22, j in 1:22]
    dK3 = k3gram()
    tx = tau_extract(Matrix(N), dK3, dK3)
    @assert tx.dual_ok "App-B tau-image violates (ii') — would refute B2"
    @assert tx.r == 17 "App-B support rank $(tx.r) != 17 (charpoly (x-1)x^5 f16)"
    v48 = entry_validate(tx.A, tx.At, tx.M; trace_bound=48)
    @assert !v48.ok "App-B must be excluded from G48 (tr 50 > 48)"
    @assert length(v48.reasons) == 1 && occursin("budget", v48.reasons[1]) "App-B must fail at 48 by BUDGET ALONE, got $(v48.reasons)"
    v50 = entry_validate(tx.A, tx.At, tx.M; trace_bound=50)
    @assert v50.ok "App-B invalid at trace 50: $(v50.reasons)"
    @assert v50.data.t == 50 "App-B trace != 50"
    @assert (v50.data.p, v50.data.q) == (3, 14) "App-B sig(L) != (3,14): got $((v50.data.p, v50.data.q))"
    Rx, x = polynomial_ring(QQ, "x")
    fac = collect(factor(charpoly(Rx, matrix(QQ, v50.data.S0))))
    degs = sort([degree(f) for (f, m) in fac])
    @assert degs == [1, 16] "App-B charpoly not (x-1)*f16: degrees $degs"
    # full-Gamma FS reconstruction (published: smooth, all moduli stabilized)
    dz = mzz(dK3); Nz = mzz(N)
    sa1 = side_analysis(matrix(QQ, Nz * dz * transpose(Nz) * dz), matrix(QQ, dz))
    sa2 = side_analysis(matrix(QQ, transpose(Nz) * dz * Nz * dz), matrix(QQ, dz))
    @assert sa1.ok && sa2.ok "App-B 5.1(iii) broken"
    p1 = perp_root_check(sa1.factors, dz, [-2])
    p2 = perp_root_check(sa2.factors, dz, [-2])
    @assert p1.status == :clean && p2.status == :clean "App-B root check not clean"
    key = canonical_key_targeted(tx)
    entry = entry_record("targeted-r17-AppB-trace50", tx.A, tx.At, tx.M, v50.data, key)
    push!(entry, "enumerative" => false)
    push!(entry, "note" => "published Q=25 flux (arXiv:2010.10519 App B): OUTSIDE G48 " *
          "by budget alone (tr 50 > 48); valid entry of the trace<=50 extension; FS " *
          "reconstruction clean both sides (smooth, matches published claims)")
    (report=Dict("status" => "PASS", "r" => 17, "trace" => 50, "sig" => [3, 14],
                 "excluded_from_G48_by" => "budget only (single reason at bound 48)",
                 "valid_at_50" => true, "R_check" => "clean both sides",
                 "wall_s" => round(time() - t0, digits=1)), entry=entry)
end

function ctrl_SOUND44(r1, r2pp)
    J = JSON.parsefile(joinpath(BANK, "sound_minvecs.json"))
    found = 0; ranks = Dict{Int,Int}()
    for row in J["rows"]
        N  = [Int(row["N"][i][j]) for i in 1:2, j in 1:2]
        T1 = [Int(row["T1"][i][j]) for i in 1:2, j in 1:2]
        T2 = [Int(row["T2"][i][j]) for i in 1:2, j in 1:2]
        tx = tau_extract(Matrix(N), T1, T2)
        @assert tx.dual_ok "(ii') fails on sound row $(row)"
        v = entry_validate(tx.A, tx.At, tx.M)
        @assert v.ok "sound row invalid: $(v.reasons) at $(row["d1"]),$(row["k1"])x$(row["d2"]),$(row["k2"])"
        @assert v.data.nflux == row["nflux_min"] "N_flux mismatch on sound row"
        key = canonical_key(tx.A, tx.At, tx.M, v.data.p, v.data.q)
        pool = tx.r == 1 ? r1.entries : r2pp.entries
        @assert haskey(pool, key) "sound row class MISSING from emitted list: $key"
        found += 1; ranks[tx.r] = get(ranks, tx.r, 0) + 1
    end
    @assert found == 44 "expected 44 sound classes, found $found"
    Dict("status" => "PASS", "rows_in_budget" => 44, "found_in_list" => found,
         "support_rank_histogram" => ranks,
         "note" => "pure-TT (NS-block = 0, corner (0,0)) stratum of G48; " *
                   "minimal fluxes of the bank/sound_minvecs.json in-budget rows")
end

# -----------------------------------------------------------------------------
# main
# -----------------------------------------------------------------------------
function sortkey_of(e::Vector{Pair{String,Any}})
    d = Dict(e)
    (d["rank"], d["stratum"], abs(d["detA"]), abs(d["detAt"]), d["trace"], d["key"])
end

function run(mode::String)
    println(stderr, "harness mode=$mode  (Oscar $(Oscar.VERSION_NUMBER))")
    t0 = time()
    r1 = emit_rank1()
    r2pp = emit_rank2_definite(1)
    r2nn = emit_rank2_definite(-1)
    wall_emit = time() - t0
    println(stderr, "emitted: r1=$(length(r1.entries)) r2pp=$(length(r2pp.entries)) " *
            "r2nn=$(length(r2nn.entries)) in $(round(wall_emit, digits=1)) s")
    t1 = time()
    comp = completeness_check(r1, r2pp, r2nn)
    wall_comp = time() - t1
    controls = Dict{String,Any}()
    targeted = Vector{Pair{String,Any}}[]
    wall_ctrl = 0.0
    if mode in ("controls", "all")
        t2 = time()
        i4 = ctrl_I4()
        controls["CTRL_I4"] = i4.report
        controls["CTRL_I1"] = ctrl_I1(r2pp)
        controls["CTRL_I2"] = ctrl_I2(r1)
        controls["CTRL_I3"] = ctrl_I3(i4.entry)
        appb = ctrl_APPB()
        controls["CTRL_APPB"] = appb.report
        controls["CTRL_SOUND44"] = ctrl_SOUND44(r1, r2pp)
        push!(targeted, i4.entry, appb.entry)
        wall_ctrl = time() - t2
    end
    if mode in ("emit", "all")
        entries = vcat(collect(values(r1.entries)), collect(values(r2pp.entries)),
                       collect(values(r2nn.entries)))
        sort!(entries, by=sortkey_of)
        meta = Pair{String,Any}[
            "record" => "G48_META",
            "spec" => "trace-<=48 Gram/glue list ((ii') dual-integrality + det-chain |detA*detAt|<=|chi(0)|; coherent-tuple classes; corner (m,mt)=(0,0) or m*mt>=1; budget (1/2)tr+m*mt<=24)",
            "generator" => "harness.jl on lattice.jl (LatticeLayer) + predicates.jl side_analysis/perp_root_check (exact arithmetic end-to-end)",
            "coverage" => Dict(
                "complete_strata" => ["r1 (rank 1, both signs; closed form)",
                    "r2pp (rank 2 sig (2,0); SL2 transversal incl imprimitive x FP(kron)<=48 x coherent dedup)",
                    "r2nn (rank 2 sig (0,2); mirror)"],
                "open_strata" => ["rank 2 sig (1,1) (Pell/unit reduction of the Aut-orbit needed; see s11.jl)",
                    "ranks 3..22 (validated entry-by-entry via tau_extract on explicit fluxes; not enumerated here)"],
                "fs_relevance" => "NO emitted entry passes R-b kernel-definiteness (needs p=3 or r-p=19, i.e. rank >= 3): emitted strata are the LIST-coverage layer (B2), not FS candidates; FS-capable strata start at rank 3"),
            "n_entries" => length(entries),
            "strata_counts" => Dict("r1" => length(r1.entries),
                                    "r2pp" => length(r2pp.entries),
                                    "r2nn" => length(r2nn.entries)),
            "walls_s" => Dict("emit" => round(wall_emit, digits=2),
                              "completeness" => round(wall_comp, digits=2),
                              "controls" => round(wall_ctrl, digits=2)),
            "completeness_check" => comp,
            "controls" => controls,
            "targeted_nonenumerative_entries_appended" => length(targeted)]
        open(joinpath(GLUEDIR, "G48_LIST.jsonl"), "w") do io
            println(io, json_obj(meta))
            for e in entries; println(io, json_obj(e)); end
            for e in targeted; println(io, json_obj(e)); end
        end
        println(stderr, "wrote G48_LIST.jsonl: $(length(entries)) enumerative + $(length(targeted)) targeted")
    end
    open(joinpath(GLUEDIR, "controls_report.json"), "w") do io
        JSON.print(io, Dict("controls" => controls,
                            "completeness" => comp,
                            "walls_s" => Dict("emit" => wall_emit, "completeness" => wall_comp,
                                              "controls" => wall_ctrl)), 1)
    end
    println("HARNESS OK: entries=$(length(r1.entries) + length(r2pp.entries) + length(r2nn.entries)) " *
            "walls: emit=$(round(wall_emit, digits=1))s comp=$(round(wall_comp, digits=1))s " *
            "ctrl=$(round(wall_ctrl, digits=1))s")
end

run(isempty(ARGS) ? "all" : ARGS[1])

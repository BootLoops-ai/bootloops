# s11.jl — rank-2 signature-(1,1) certified enumeration engine for the
# trace-<=48 Gram/glue list (the indefinite rank-2 stratum that harness.jl does
# not enumerate). One "cell" = a pair (A, At) of even sig(1,1) Grams.
#   E1 S0-walk: integral A-self-adjoint trace-t S0 = (t/2)I + T, T in the rank-2
#      trace-0 self-adjoint lattice L0; det condition <=> q_L0(z) = -D0chi/4,
#      an integral binary form of disc s^2*D_A; solved by the Gauss b0/cycle
#      route (nonsquare disc, transversal mod Aut+) or exact divisor
#      factorization (square disc, ALL solutions).
#   E2 M-recovery: N := S0*A^-1 symmetric even integral (else discard); row1 =
#      representations of N11/2 by f_At (transversal mod Aut+(f_At) = right
#      units EXACTLY, since P_At is the fundamental automorph); row2 = exact
#      line + quadratic solve.  n1 = 0 isotropic-line branch for square D_At.
#   E3 domain filter: implemented as the E5 ball-fold, not as exact (u,v)
#      comparisons: candidates are finite per (t,m) by construction, so a
#      centered fundamental parallelogram is not needed for termination;
#      uniqueness is delivered by the canonical key.
#   E4 LAW: entry_validate is a byte-identical copy of the one in harness.jl
#      (source-diffed at startup against that file); asserts t / chi0 grid
#      consistency.
#   E5 fold: canonical key = lex-min flatkey over the finite set of unit- and
#      F-translates M -> P^-1 M Q^-T with ||M'||_inf <= KEYCAP (ball criterion
#      is class-intrinsic: orbit-equal M's see the same orbit-cap-ball).
# Gates: V1 cell (U,U): count 1132, fiber pins, two-route divisor count, and a
#   byte-exact comparison against an independently produced reference list of
#   the same cell (not included in the package: set TERRIER_SWEEP_BANK to its
#   root, see run_v1; V1 refuses loudly when it is unset);
#   V2 cells (U,U(2)) and (U,[[2,2],[2,0]]) vs bank/s11_expected_counts.json;
#   V3 cell (A5,A5) vs bank/s11_expected_counts.json + Pell pins (t1,u1)=(3,1), P5;
#   MUT mutation control: the fundamental Pell unit P5 at D=5 is replaced by its
#   cube P5^3 (trace 18, i.e. the non-fundamental solution (18,8) of
#   t^2-5u^2=4 that a continued-fraction convergent route returns) — the V3
#   gates must then FAIL.
# Run: nice -n 5 julia +1.10 -t 1 --project=<your Oscar env> \
#        s11.jl all      (modes: v1 | v2 | v3 | mut | all)
# Writes R2S11_<cell>_ENGINE.jsonl slices and ENGINE_REPORT.json next to this file.

import JSON
using Oscar

const S11DIR = @__DIR__
const BANK = joinpath(S11DIR, "bank")            # shipped fixtures (sha256-pinned)
include(joinpath(S11DIR, "predicates.jl"))       # exact predicate chain (SweepPredicates)
using .SweepPredicates: side_analysis, congruence_diag

const TRACE_BOUND = 48
const NCAP = 24
const KEYCAP = BigInt(10)^6   # E5 ball cap (validation cells)

# =============================================================================
# VERBATIM BLOCK — copied byte-for-byte from harness.jl
# (entry_validate is source-diffed against the harness at startup; see
# verbatim_check() below). Nothing in this block may be edited.
# =============================================================================
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
# --- END entry_validate (verbatim) -------------------------------------------

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

json_obj(pairs::Vector{<:Pair}) =
    "{" * join(["\"$(k)\":" * (v isa Vector{<:Pair} ? json_obj(v) : JSON.json(v))
                for (k, v) in pairs], ",") * "}"

flatkey(s, GA, GT, M) = string(s, "|", tomat(GA), "|", tomat(GT), "|", tomat(M))

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

function sortkey_of(e::Vector{Pair{String,Any}})
    d = Dict(e)
    (d["rank"], d["stratum"], abs(d["detA"]), abs(d["detAt"]), d["trace"], d["key"])
end
# ============================= END VERBATIM BLOCK ============================

"extract the entry_validate source text from a file (signature line .. first ^end)."
function func_text(path::String, header::String)
    lines = readlines(path)
    i = findfirst(l -> startswith(l, header), lines)
    @assert i !== nothing "no $header in $path"
    j = findnext(l -> l == "end", lines, i)
    @assert j !== nothing
    join(lines[i:j], "\n")
end

function verbatim_check()
    a = func_text(joinpath(S11DIR, "s11.jl"), "function entry_validate(")
    b = func_text(joinpath(S11DIR, "harness.jl"), "function entry_validate(")
    @assert a == b "entry_validate DIFFERS from harness.jl — LAW violation"
    true
end

# =============================================================================
# Binary quadratic form toolkit — BigInt end-to-end.
# form = (a,b,c) :: NTuple{3,BigInt}; 2x2 transform = (p11,p12,p21,p22).
# =============================================================================
const F3 = NTuple{3,BigInt}; const M4 = NTuple{4,BigInt}
fdisc(f::F3) = f[2]^2 - 4*f[1]*f[3]
feval(f::F3, x, y) = f[1]*x^2 + f[2]*x*y + f[3]*y^2
mul4(P::M4, Q::M4) = (P[1]*Q[1]+P[2]*Q[3], P[1]*Q[2]+P[2]*Q[4],
                      P[3]*Q[1]+P[4]*Q[3], P[3]*Q[2]+P[4]*Q[4])
det4(P::M4) = P[1]*P[4] - P[2]*P[3]
function inv4(P::M4)
    s = det4(P); @assert abs(s) == 1
    (s*P[4], -s*P[2], -s*P[3], s*P[1])
end
const ID4 = (big(1), big(0), big(0), big(1))
"f∘P: coefficients of f(P*(x,y))."
fapply(f::F3, P::M4) = (feval(f, P[1], P[3]),
                        2*f[1]*P[1]*P[2] + f[2]*(P[1]*P[4] + P[2]*P[3]) + 2*f[3]*P[3]*P[4],
                        feval(f, P[2], P[4]))
issq(n::BigInt) = n >= 0 && isqrt(n)^2 == n
isreducedf(f::F3) = f[1]*f[3] < 0 && f[2] > abs(f[1] + f[3])

"one rho step (indefinite, nonsquare disc); returns (g, T) with f∘T = g.
 Normalization window (Buchmann–Lenstra): |c| > sqrt(D): -|c| < r <= |c|;
 else: sqrt(D) - 2|c| < r <= sqrt(D) (the reduced-cycle window)."
function rhostep(f::F3)
    a, b, c = f; D = fdisc(f)
    @assert c != 0 && !issq(D)
    s2 = isqrt(D); ac = abs(c)
    hi = ac > s2 ? ac : s2
    r = mod(-b, 2*ac)
    while r > hi; r -= 2*ac; end
    while r <= hi - 2*ac; r += 2*ac; end
    s = divexact(b + r, 2*c)
    g = (c, r, divexact(r*r - D, 4*c))
    T = (big(0), big(-1), big(1), s)
    @assert fapply(f, T) == g
    (g, T)
end

"reduce to a reduced form; returns (fred, T) with f∘T = fred."
function reduce_form(f::F3)
    T = ID4; g = f
    for _ in 1:100000
        isreducedf(g) && return (g, T)
        g, S = rhostep(g); T = mul4(T, S)
    end
    error("reduction did not terminate on $f")
end

"the rho-cycle of a reduced form: members + cumulative transforms + full product."
function cycle_of(fr::F3)
    mem = F3[fr]; trs = M4[ID4]
    g, T = fr, ID4
    for _ in 1:100000
        g, S = rhostep(g); T = mul4(T, S)
        g == fr && return (mem, trs, T)
        push!(mem, g); push!(trs, T)
    end
    error("cycle did not close on $fr")
end

"SL2(Z) transform U with f∘U = g, or nothing (both primitive, nonsquare disc)."
function equiv_transform(f::F3, g::F3)
    @assert fdisc(f) == fdisc(g)
    fr, Tf = reduce_form(f)
    gr, Tg = reduce_form(g)
    mem, trs, _ = cycle_of(fr)
    i = findfirst(==(gr), mem)
    i === nothing && return nothing
    U = mul4(mul4(Tf, trs[i]), inv4(Tg))
    @assert fapply(f, U) == g "equiv_transform broken"
    U
end

"fundamental proper automorph of primitive f, nonsquare disc: P, (t1,u1)."
function fund_auto(f::F3)
    D = fdisc(f); @assert !issq(D) && D > 0
    fr, Tf = reduce_form(f)
    _, _, Ucyc = cycle_of(fr)
    P = mul4(mul4(Tf, Ucyc), inv4(Tf))
    @assert fapply(f, P) == f && det4(P) == 1 && P != ID4 && P != (-ID4[1],ID4[2],ID4[3],-ID4[4])
    P[1] + P[4] < 0 && (P = (-P[1], -P[2], -P[3], -P[4]))
    t1 = P[1] + P[4]
    @assert f[1] != 0 && P[3] % f[1] == 0
    u1 = divexact(P[3], f[1])
    u1 < 0 && (P = inv4(P); u1 = divexact(P[3], f[1]))
    @assert t1 == P[1] + P[4] && u1 > 0 && t1^2 - D*u1^2 == 4
    (P, t1, u1)
end

"improper automorph (det -1) of f, or nothing. Nonsquare: cycle a|b search."
function improper_auto(f::F3)
    D = fdisc(f)
    if !issq(D)
        cf = gcd(gcd(f[1], f[2]), f[3])
        f0 = (divexact(f[1],cf), divexact(f[2],cf), divexact(f[3],cf))
        fr, Tf = reduce_form(f0)
        mem, trs, _ = cycle_of(fr)
        for (g, T) in zip(mem, trs)
            (g[1] != 0 && g[2] % g[1] == 0) || continue
            P0g = (big(1), divexact(g[2], g[1]), big(0), big(-1))
            @assert fapply(g, P0g) == g
            U = mul4(Tf, T)                 # f0∘U = g
            P = mul4(mul4(U, P0g), inv4(U))
            @assert fapply(f, P) == f && det4(P) == -1
            return P
        end
        return nothing
    end
    # square disc: swap the two rational isotropic lines
    k = isqrt(D); a, b, c = f
    p1, p2 = a != 0 ? ((-(b - k), 2a), (-(b + k), 2a)) : ((big(1), big(0)), (-c, b))
    prim(v) = (g = gcd(v[1], v[2]); (divexact(v[1], g), divexact(v[2], g)))
    p1, p2 = prim(p1), prim(p2)
    Bq = matrix(QQ, 2, 2, [p1[1], p2[1], p1[2], p2[2]])
    for s1 in (1, -1), s2 in (1, -1)
        Pq = Bq * matrix(QQ, 2, 2, [0, s1, s2, 0]) * inv(Bq)
        is_intmat(Pq) || continue
        P = (BigInt(numerator(Pq[1,1])), BigInt(numerator(Pq[1,2])),
             BigInt(numerator(Pq[2,1])), BigInt(numerator(Pq[2,2])))
        det4(P) == -1 && fapply(f, P) == f && return P
    end
    nothing
end

divisors_of(n::BigInt) = (n = abs(n);
    sort([d for d in 1:isqrt(n) if n % d == 0] ∪
         [divexact(n, d) for d in 1:isqrt(n) if n % d == 0]))

"proper representations of n != 0 by PRIMITIVE f (nonsquare disc D > 0),
 one per Aut+(f)-orbit (Gauss b0-route: b0^2 = D mod 4n, g=(n,b0,*) ~ f)."
function proper_reps(f::F3, n::BigInt)
    @assert n != 0
    D = fdisc(f); out = Tuple{BigInt,BigInt}[]
    for b0 in 0:(2*abs(n) - 1)
        (b0^2 - D) % (4*abs(n)) == 0 || continue
        g = (n, big(b0), divexact(b0^2 - D, 4n))
        U = equiv_transform(f, g)
        U === nothing && continue
        @assert feval(f, U[1], U[3]) == n
        push!(out, (U[1], U[3]))
    end
    out
end

"ALL representations of n != 0 by f with SQUARE disc D = k^2 > 0 (finite)."
function square_disc_reps(f::F3, n::BigInt)
    @assert n != 0
    D = fdisc(f); k = isqrt(D); @assert k^2 == D && k > 0
    a, b, c = f; out = Tuple{BigInt,BigInt}[]
    if a != 0
        # u = 2ax+(b-k)y, v = 2ax+(b+k)y: uv = 4an, y = (v-u)/2k, x = (u+v-2by)/4a
        for u0 in divisors_of(4a*n), su in (1, -1)
            u = su * u0; (4a*n) % u == 0 || continue
            v = divexact(4a*n, u)
            (v - u) % (2k) == 0 || continue
            y = divexact(v - u, 2k)
            (u + v - 2b*y) % (4a) == 0 || continue
            x = divexact(u + v - 2b*y, 4a)
            @assert feval(f, x, y) == n
            push!(out, (x, y))
        end
    else
        # f = bxy + cy^2 = y(bx + cy), b = ±k
        for y0 in divisors_of(n), sy in (1, -1)
            y = sy * y0; n % y == 0 || continue
            w = divexact(n, y) - c*y      # = b*x
            w % b == 0 || continue
            x = divexact(w, b)
            @assert feval(f, x, y) == n
            push!(out, (x, y))
        end
    end
    sort(unique(out))
end

"representations of n != 0 by f (disc > 0). Returns (reps, complete):
 complete=true  => ALL integer solutions (square disc);
 complete=false => transversal modulo Aut+(f0) (nonsquare; expand by units)."
function reps_by_form(f::F3, n::BigInt)
    n == 0 && error("reps_by_form: n = 0 must be handled by the caller")
    cf = gcd(gcd(f[1], f[2]), f[3])
    n % cf == 0 || return (Tuple{BigInt,BigInt}[], true)
    f0 = (divexact(f[1], cf), divexact(f[2], cf), divexact(f[3], cf))
    n0 = divexact(n, cf)
    issq(fdisc(f0)) && return (square_disc_reps(f0, n0), true)
    out = Tuple{BigInt,BigInt}[]
    for e in divisors_of(n0)
        n0 % e^2 == 0 || continue
        append!(out, [(e*x, e*y) for (x, y) in proper_reps(f0, divexact(n0, e^2))])
    end
    (sort(unique(out)), false)
end

# =============================================================================
# Cell machinery
# =============================================================================
zz4(P::M4) = matrix(ZZ, 2, 2, [P[1], P[2], P[3], P[4]])
m4_of(Z::ZZMatrix) = (BigInt(Z[1,1]), BigInt(Z[1,2]), BigInt(Z[2,1]), BigInt(Z[2,2]))
gram_eval(G::ZZMatrix, v) = BigInt(v[1]*(G[1,1]*v[1] + G[1,2]*v[2]) +
                                   v[2]*(G[2,1]*v[1] + G[2,2]*v[2]))
form_of_gram(G::ZZMatrix) = (BigInt(divexact(G[1,1], 2)), BigInt(G[1,2]),
                             BigInt(divexact(G[2,2], 2)))
prim_part(f::F3) = (cf = gcd(gcd(f[1], f[2]), f[3]);
                    ((divexact(f[1],cf), divexact(f[2],cf), divexact(f[3],cf)), cf))

"unit data of an even sig(1,1) Gram G: proper Pell generator (or nothing if
 square disc), improper automorph (or nothing). mutate3: replace P by P^3
 (non-fundamental Pell unit; mutation control)."
function unit_data(G::ZZMatrix; mutate3::Bool=false)
    f = form_of_gram(G); D = fdisc(f)
    @assert D > 0 "not sig(1,1)"
    P = nothing; tu = nothing
    if !issq(D)
        f0, _ = prim_part(f)
        Pm, t1, u1 = fund_auto(f0)
        mutate3 && (Pm = mul4(Pm, mul4(Pm, Pm)); t1 = Pm[1] + Pm[4];
                    u1 = divexact(Pm[3], f0[1]))
        Z = zz4(Pm)
        @assert transpose(Z) * G * Z == G "Pell automorph does not preserve Gram"
        P = Z; tu = (t1, u1)
    end
    imp = improper_auto(f)
    impZ = imp === nothing ? nothing : zz4(imp)
    impZ !== nothing && @assert transpose(impZ) * G * impZ == G && det(impZ) == -1
    (P=P, tu=tu, imp=impZ, D=D, f=f)
end

"trace-0 A-self-adjoint lattice basis: (w,x,y) with A11*x - A22*y - 2*A12*w = 0."
function selfadj_basis(A::ZZMatrix)
    C = matrix(ZZ, 1, 3, [-2*A[1,2], A[1,1], -A[2,2]])
    K = kernel(C, side=:right)
    @assert ncols(K) == 2
    Tm(col) = matrix(ZZ, 2, 2, [K[1,col], K[2,col], K[3,col], -K[1,col]])
    (Tm(1), Tm(2), K)
end

function build_cell(A::ZZMatrix, At::ZZMatrix; mutate3::Bool=false)
    uL = unit_data(A; mutate3=mutate3)
    uR = unit_data(At; mutate3=mutate3)
    T1, T2, K = selfadj_basis(A)
    q = (BigInt(det(T1)), BigInt(det(T1 + T2)) - BigInt(det(T1)) - BigInt(det(T2)),
         BigInt(det(T2)))
    # disc(q) = s^2 * D_A
    @assert fdisc(q) > 0 && (issq(fdisc(q)) == issq(uL.D))
    eps_z = nothing; kpow = 0
    if uL.P !== nothing
        q0, _ = prim_part(q)
        ez, _, _ = fund_auto(q0)
        # gamma: z-coords of T -> P^-1 T P on basis (T1,T2)
        Pi = zz_of(inv(matrix(QQ, uL.P)))
        cols = ZZMatrix[Pi * T1 * uL.P, Pi * T2 * uL.P]
        Bm = matrix(QQ, 3, 2, [QQ(K[i,j]) for i in 1:3 for j in 1:2])
        gcols = map(cols) do Tg
            rhs = matrix(QQ, 3, 1, [QQ(Tg[1,1]), QQ(Tg[1,2]), QQ(Tg[2,1])])
            ok, z = can_solve_with_solution(Bm, rhs, side=:right)
            @assert ok && is_intmat(z)
            (BigInt(numerator(z[1,1])), BigInt(numerator(z[2,1])))
        end
        gam = (gcols[1][1], gcols[2][1], gcols[1][2], gcols[2][2])
        @assert fapply(q, gam) == q && abs(det4(gam)) == 1
        # match gam = ± ez^±j (exact); j is the expansion window
        found = false
        for base in (ez, inv4(ez))
            Pj = ID4
            for j in 1:256
                Pj = mul4(Pj, base)
                if gam == Pj || gam == (-Pj[1], -Pj[2], -Pj[3], -Pj[4])
                    eps_z = base; kpow = j; found = true; break
                end
            end
            found && break
        end
        @assert found "unit action is not a power of the fundamental q-automorph (perturbed unit?)"
    end
    (A=A, At=At, uL=uL, uR=uR, T1=T1, T2=T2, q=q, eps_z=eps_z, kpow=kpow,
     DA=BigInt(-det(A)), DAt=BigInt(-det(At)),
     fAt=form_of_gram(At))
end

# --- E1: all integral A-self-adjoint S0 with tr = t, det = Delta, modulo
#     nothing (complete finite candidate set; unit expansion window included) --
function s0_candidates(cell, t::Int, Delta::BigInt)
    D0 = BigInt(t)^2 - 4*Delta
    D0 <= 0 && return ZZMatrix[]            # scalar/complex: Validity Lemma
    @assert D0 % 4 == 0
    nt = -divexact(D0, 4)
    reps, complete = reps_by_form(cell.q, nt)
    zs = Set{Tuple{BigInt,BigInt}}()
    if complete
        for z in reps; push!(zs, z); push!(zs, (-z[1], -z[2])); end
    else
        @assert cell.eps_z !== nothing
        for z0 in reps
            zj = z0
            for _ in 1:(cell.kpow + 1)      # window j = 0..k (margin 1)
                push!(zs, zj); push!(zs, (-zj[1], -zj[2]))
                zj = (cell.eps_z[1]*zj[1] + cell.eps_z[2]*zj[2],
                      cell.eps_z[3]*zj[1] + cell.eps_z[4]*zj[2])
            end
            zi = z0; ei = inv4(cell.eps_z)
            for _ in 1:(cell.kpow + 1)      # window j = -1..-k-1 (margin)
                zi = (ei[1]*zi[1] + ei[2]*zi[2], ei[3]*zi[1] + ei[4]*zi[2])
                push!(zs, zi); push!(zs, (-zi[1], -zi[2]))
            end
        end
    end
    out = ZZMatrix[]
    for (z1, z2) in sort(collect(zs))
        S0 = divexact(t, 2) * identity_matrix(ZZ, 2) + ZZ(z1)*cell.T1 + ZZ(z2)*cell.T2
        @assert tr(S0) == t && det(S0) == Delta
        @assert transpose(cell.A * S0) == cell.A * S0   # A-self-adjoint
        push!(out, S0)
    end
    out
end

"isotropic primitive directions of a square-disc form."
function isotropic_dirs(f::F3)
    D = fdisc(f); k = isqrt(D); @assert k^2 == D && k > 0
    a, b, c = f
    p1, p2 = a != 0 ? ((-(b - k), 2a), (-(b + k), 2a)) : ((big(1), big(0)), (-c, b))
    prim(v) = (g = gcd(v[1], v[2]); (divexact(v[1], g), divexact(v[2], g)))
    (prim(p1), prim(p2))
end

"solve row2: r1 At r2^T = N12, r2 At r2^T = N22 (+ det pin in degenerate case)."
function solve_row2(cell, r1, N12::BigInt, N22::BigInt, m::BigInt)
    At = cell.At
    al = r1[1]*BigInt(At[1,1]) + r1[2]*BigInt(At[2,1])
    be = r1[1]*BigInt(At[1,2]) + r1[2]*BigInt(At[2,2])
    (al == 0 && be == 0) && return Tuple{BigInt,BigInt}[]
    g = gcd(al, be)
    N12 % g == 0 || return Tuple{BigInt,BigInt}[]
    gg, uu, vv = gcdx(al, be); @assert gg == g
    sc = divexact(N12, g)
    v0 = (uu*sc, vv*sc)
    d = (-divexact(be, g), divexact(al, g))
    a2 = gram_eval(At, d)
    b2 = 2*(d[1]*(BigInt(At[1,1])*v0[1] + BigInt(At[1,2])*v0[2]) +
            d[2]*(BigInt(At[2,1])*v0[1] + BigInt(At[2,2])*v0[2]))
    c2 = gram_eval(At, v0) - N22
    ss = BigInt[]
    if a2 != 0
        disc2 = b2^2 - 4*a2*c2
        if disc2 >= 0 && issq(disc2)
            r = isqrt(disc2)
            for pm in (r, -r)
                (-b2 + pm) % (2*a2) == 0 && push!(ss, divexact(-b2 + pm, 2*a2))
            end
        end
    elseif b2 != 0
        (-c2) % b2 == 0 && push!(ss, divexact(-c2, b2))
    elseif c2 == 0
        dd = r1[1]*d[2] - r1[2]*d[1]
        d0 = r1[1]*v0[2] - r1[2]*v0[1]
        @assert dd != 0 "degenerate row2 family with det-flat direction (must not occur)"
        for tgt in (m, -m)
            (tgt - d0) % dd == 0 && push!(ss, divexact(tgt - d0, dd))
        end
    end
    [(v0[1] + s*d[1], v0[2] + s*d[2]) for s in sort(unique(ss))]
end

# --- E2: recover all M (mod right units) with M At M^T = N := S0 A^-1 --------
function recover_M(cell, S0::ZZMatrix, m::BigInt)
    Nq = matrix(QQ, S0) * inv(matrix(QQ, cell.A))
    is_intmat(Nq) || return ZZMatrix[]
    N = zz_of(Nq)
    @assert N == transpose(N) "S0 A^-1 not symmetric (self-adjointness broken)"
    (iseven(N[1,1]) && iseven(N[2,2])) || return ZZMatrix[]
    N11, N12, N22 = BigInt(N[1,1]), BigInt(N[1,2]), BigInt(N[2,2])
    r1s = Tuple{BigInt,BigInt}[]
    if N11 == 0
        issq(fdisc(cell.fAt)) || return ZZMatrix[]   # anisotropic: r1 = 0 only
        N12 == 0 && return ZZMatrix[]                # r2 || r1 => det M = 0
        for p in isotropic_dirs(cell.fAt), e0 in divisors_of(N12), se in (1, -1)
            push!(r1s, (se*e0*p[1], se*e0*p[2]))
        end
    else
        n1 = divexact(N11, 2)
        reps, _ = reps_by_form(cell.fAt, n1)
        # transversal mod Aut+(f_At) == right-unit action on row1 EXACTLY
        r1s = reps
    end
    out = ZZMatrix[]
    for r1 in sort(unique(r1s))
        for r2 in solve_row2(cell, r1, N12, N22, m)
            M = matrix(ZZ, 2, 2, [r1[1], r1[2], r2[1], r2[2]])
            iszero(det(M)) && continue
            M * cell.At * transpose(M) == N || continue
            push!(out, M)
        end
    end
    out
end

# --- E5: canonical key = lex-min flatkey over unit/F-translates in the ball --
function fold_translates(cell, M::ZZMatrix; cap::BigInt=KEYCAP)
    gensL = ZZMatrix[-identity_matrix(ZZ, 2)]
    cell.uL.P !== nothing && push!(gensL, cell.uL.P, zz_of(inv(matrix(QQ, cell.uL.P))))
    cell.uL.imp !== nothing && push!(gensL, cell.uL.imp)
    gensR = ZZMatrix[-identity_matrix(ZZ, 2)]
    cell.uR.P !== nothing && push!(gensR, cell.uR.P, zz_of(inv(matrix(QQ, cell.uR.P))))
    cell.uR.imp !== nothing && push!(gensR, cell.uR.imp)
    movesL = [zz_of(inv(matrix(QQ, P))) for P in gensL]           # M -> P^-1 M
    movesR = [transpose(zz_of(inv(matrix(QQ, Q)))) for Q in gensR] # M -> M Q^-T
    seen = Set{NTuple{4,BigInt}}([m4_of(M)]); frontier = [M]; all_ = [M]
    while !isempty(frontier)
        nxt = ZZMatrix[]
        for X in frontier, mv in Iterators.flatten((
                (L * X for L in movesL), (X * R for R in movesR)))
            key = m4_of(mv)
            maximum(abs, key) <= cap || continue
            key in seen && continue
            push!(seen, key); push!(nxt, mv); push!(all_, mv)
        end
        frontier = nxt
        @assert length(seen) < 500000 "E5 fold runaway"
    end
    all_
end

function canonical_key(cell, M::ZZMatrix; cap::BigInt=KEYCAP)
    tag = (cell.DA == 1 && cell.DAt == 1 && cell.A == cell.At) ? "s11uu" : "s11"
    minimum(flatkey(tag, cell.A, cell.At, X) for X in fold_translates(cell, M; cap=cap))
end

# --- engine loop over the (t, m) grid ---------------------------------------
function emit_cell(A::ZZMatrix, At::ZZMatrix; tmax::Int=TRACE_BOUND, mutate3::Bool=false)
    cell = build_cell(A, At; mutate3=mutate3)
    entries = Dict{String,Vector{Pair{String,Any}}}()
    stats = Dict{String,Int}("candidates" => 0, "dups" => 0, "rejected" => 0,
                             "s0_candidates" => 0)
    fib = Dict{Tuple{Int,Int},Int}()        # (t,m) -> emitted entries
    s0fib = Dict{Tuple{Int,Int},Int}()      # (t,m) -> S0 candidates
    for t in 2:2:tmax
        mmax = isqrt(div(BigInt(div(t,2))^2, cell.DA * cell.DAt))
        for m in 1:Int(mmax)
            Delta = BigInt(m)^2 * cell.DA * cell.DAt
            4*Delta < BigInt(t)^2 || continue        # D0chi > 0 (Validity Lemma)
            s0s = s0_candidates(cell, t, Delta)
            s0fib[(t, m)] = get(s0fib, (t, m), 0) + length(s0s)
            stats["s0_candidates"] += length(s0s)
            for S0 in s0s, M in recover_M(cell, S0, BigInt(m))
                stats["candidates"] += 1
                val = entry_validate(A, At, M)
                if !val.ok
                    stats["rejected"] += 1; continue
                end
                @assert val.data.t == t "E4 grid consistency: t"
                @assert val.data.chi0 == Int(Delta) "E4 grid consistency: chi0"
                key = canonical_key(cell, M)
                if haskey(entries, key)
                    stats["dups"] += 1; continue
                end
                entries[key] = entry_record("r2s11", A, At, M, val.data, key)
                fib[(t, m)] = get(fib, (t, m), 0) + 1
            end
        end
    end
    (entries=entries, stats=stats, fib=fib, s0fib=s0fib, cell=cell)
end

function write_slice(path::String, meta::Vector{Pair{String,Any}}, dicts...)
    entries = Vector{Pair{String,Any}}[]
    for dd in dicts, (_, e) in dd; push!(entries, e); end
    sort!(entries, by=sortkey_of)
    open(path, "w") do io
        println(io, json_obj(meta))
        for e in entries; println(io, json_obj(e)); end
    end
    length(entries)
end

fibdict(d) = Dict("$(k[1]),$(k[2])" => v for (k, v) in d)
tracedict(res) = (out = Dict{Int,Int}();
    for (_, e) in res.entries
        t = Dict(e)["trace"]; out[t] = get(out, t, 0) + 1
    end; out)

"gate helper: hard assert with a registered name."
function gate!(report, name::String, cond::Bool, detail="")
    report[name] = Dict("pass" => cond, "detail" => string(detail))
    cond || error("GATE FAIL: $name — $detail")
    println(stderr, "  gate $name: PASS $(detail == "" ? "" : "($detail)")")
end

const UZ = matrix(ZZ, 2, 2, [0, 1, 1, 0])
const U2Z = matrix(ZZ, 2, 2, [0, 2, 2, 0])
const G22Z = matrix(ZZ, 2, 2, [2, 2, 2, 0])
const A5Z = matrix(ZZ, 2, 2, [2, 1, 1, -2])
const P5Z = matrix(ZZ, 2, 2, [1, 1, 1, 2])  # V3 pin: fundamental Pell unit at D=5, (t1,u1) = (3,1)

# --- V1: U x U cell, byte-exact against the reference list ------------------
function run_v1(report)
    println(stderr, "V1: U x U cell (full engine E1-E5)")
    res = emit_cell(UZ, UZ)
    rep = Dict{String,Any}()
    gate!(rep, "V1.count_1132", length(res.entries) == 1132, length(res.entries))
    # reference comparison: key set, per-line byte-exact record reproduction.
    # The reference list is not included in the package (TERRIER_SWEEP_BANK).
    ts = get(ENV, "TERRIER_SWEEP_BANK", "")
    isempty(ts) && error("REFUSE: env var TERRIER_SWEEP_BANK is unset — " *
        "V1 compares against a reference list of the (U,U) cell, " *
        "<root>/glue/R2S11_UU_PILOT.jsonl, which is not included in the package")
    lines = readlines(joinpath(ts, "glue", "R2S11_UU_PILOT.jsonl"))
    @assert occursin("R2S11_UU_META", lines[1])
    plines = lines[2:end]
    gate!(rep, "V1.reference_count", length(plines) == 1132, length(plines))
    nbyte = 0; pkeys = Set{String}()
    for ln in plines
        e = JSON.parse(ln)
        Mp = matrix(ZZ, 2, 2, [e["M"][1][1], e["M"][1][2], e["M"][2][1], e["M"][2][2]])
        k = canonical_key(res.cell, Mp)
        @assert k == e["key"] "recomputed key differs from reference key"
        push!(pkeys, k)
        val = entry_validate(UZ, UZ, Mp)
        @assert val.ok
        rec = entry_record("r2s11", UZ, UZ, Mp, val.data, k)
        json_obj(rec) == ln && (nbyte += 1)
    end
    gate!(rep, "V1.keyset_equal", pkeys == Set(keys(res.entries)))
    gate!(rep, "V1.byte_exact_1132", nbyte == 1132, "$nbyte/1132 reference lines byte-reproduced")
    # fibers: trace-6 = 2; I4 U-block fibers 2/2/4 (charpolys of harness.jl ctrl_I4)
    cps = Dict{String,Int}()
    for (_, e) in res.entries
        d = Dict(e); cps[d["charpoly"]] = get(cps, d["charpoly"], 0) + 1
    end
    gate!(rep, "V1.trace6_fiber_2", get(cps, "x^2 - 6*x + 1", 0) == 2)
    gate!(rep, "V1.I4_fibers_2_2_4",
          get(cps, "x^2 - 8*x + 4", 0) == 2 && get(cps, "x^2 - 10*x + 1", 0) == 4)
    # two-route cross-check: per (t,m) S0-candidate count == 2*tau((t/2)^2 - m^2)
    ok2 = true
    for ((t, m), n) in res.s0fib
        nt = BigInt(div(t,2))^2 - BigInt(m)^2
        nt > 0 || continue
        ndiv = 2 * length(divisors_of(nt))
        ndiv == n || (ok2 = false; println(stderr, "  two-route mismatch at ($t,$m): $n vs $ndiv"))
    end
    gate!(rep, "V1.s0_two_route_divisor", ok2)
    gate!(rep, "V1.fs_all_false",
          all(Dict(Dict(e)["fs_check"])["pass"] == false for (_, e) in res.entries))
    meta = Pair{String,Any}["record" => "R2S11_UU_ENGINE_META",
        "spec" => "s11.jl E1-E5 full route on the (detA,detAt)=(-1,-1) cell; gates V1",
        "stats" => res.stats, "per_tm_entries" => fibdict(res.fib)]
    write_slice(joinpath(S11DIR, "R2S11_UU_ENGINE.jsonl"), meta, res.entries)
    rep["stats"] = res.stats
    report["V1"] = rep
end

# --- V2: cells (U, U(2)) and (U, [[2,2],[2,0]]) vs expected counts ----------
function run_v2(report)
    expected = JSON.parsefile(joinpath(BANK, "s11_expected_counts.json"))
    rep = Dict{String,Any}()
    for (tagp, At, pk) in [("V2a", U2Z, "V2a_U_U2"), ("V2b", G22Z, "V2b_U_G22")]
        println(stderr, "V2: cell (U, $(tomat(At)))")
        res = emit_cell(UZ, At)
        pr = expected[pk]
        gate!(rep, "$tagp.entries_expected", length(res.entries) == pr["entries"],
              "$(length(res.entries)) == $(pr["entries"])")
        td = tracedict(res)
        prt = Dict(parse(Int, k) => v for (k, v) in pr["per_trace"])
        gate!(rep, "$tagp.per_trace_expected", td == prt, "$(sort(collect(td)))")
        gate!(rep, "$tagp.fs_all_false",
              all(Dict(Dict(e)["fs_check"])["pass"] == false for (_, e) in res.entries))
        meta = Pair{String,Any}["record" => "R2S11_$(tagp)_ENGINE_META",
            "cell" => Dict("A" => tomat(UZ), "At" => tomat(At)),
            "stats" => res.stats, "per_tm_entries" => fibdict(res.fib),
            "expected_counts" => "bank/s11_expected_counts.json (sha256-pinned)"]
        write_slice(joinpath(S11DIR, "R2S11_$(tagp)_ENGINE.jsonl"), meta, res.entries)
        rep["$tagp.stats"] = res.stats
    end
    gate!(rep, "V2.cells_distinct_counts",
          expected["V2a_U_U2"]["entries"] != expected["V2b_U_G22"]["entries"])
    report["V2"] = rep
end

# --- V3: cell (A5, A5) vs expected counts + P5 pins; mutate3 = mutation control
function run_v3(report; mutate3::Bool=false)
    println(stderr, "V3: cell (A5, A5) mutate3=$mutate3")
    expected = JSON.parsefile(joinpath(BANK, "s11_expected_counts.json"))["V3_A5_A5"]
    rep = Dict{String,Any}()
    # unit pins (computed from the cycle walk)
    f5 = form_of_gram(A5Z)
    Pf, t1, u1 = fund_auto(f5)
    gate!(rep, "V3.pell_fundamental_3_1", (t1, u1) == (3, 1), "(t1,u1)=($t1,$u1)")
    gate!(rep, "V3.P5_pin", zz4(Pf) == P5Z || zz4(inv4(Pf)) == P5Z, string(Pf))
    # Pell pitfall: the solution (18,8) of t^2-5u^2=4 is the CUBE of the fundamental (3,1)
    P3 = mul4(Pf, mul4(Pf, Pf))
    gate!(rep, "V3.pell_18_8_is_cube", P3[1] + P3[4] == 18, "tr P^3 = $(P3[1]+P3[4])")
    res = emit_cell(A5Z, A5Z; mutate3=mutate3)
    gate!(rep, "V3.entries_expected", length(res.entries) == expected["entries"],
          "$(length(res.entries)) == $(expected["entries"])")
    ptm = Dict("$(k[1]),$(k[2])" => v for (k, v) in res.fib)
    gate!(rep, "V3.per_tm_expected", ptm == Dict(String(k) => v for (k, v) in expected["per_t_m"]),
          "$(sort(collect(ptm)))")
    td = tracedict(res)
    gate!(rep, "V3.t10_fiber_empty", get(td, 10, 0) == 0,
          "scalar fiber empty (Validity Lemma: D0chi = t^2-4*Delta > 0)")
    gate!(rep, "V3.t_le_16_empty", all(t > 16 for (t, _) in td))
    # key invariance under TRUE P5 translation (left and right), every entry
    ok = true
    P5i = zz_of(inv(matrix(QQ, P5Z)))
    for (k, e) in res.entries
        M = (d = Dict(e); matrix(ZZ, 2, 2, [d["M"][1][1], d["M"][1][2], d["M"][2][1], d["M"][2][2]]))
        kL = canonical_key(res.cell, P5i * M)
        kR = canonical_key(res.cell, M * transpose(P5i))
        (kL == k && kR == k) || (ok = false)
    end
    gate!(rep, "V3.key_invariant_under_P5", ok)
    gate!(rep, "V3.fs_all_false",
          all(Dict(Dict(e)["fs_check"])["pass"] == false for (_, e) in res.entries))
    meta = Pair{String,Any}["record" => "R2S11_V3_ENGINE_META",
        "cell" => Dict("A" => tomat(A5Z), "At" => tomat(A5Z)),
        "stats" => res.stats, "per_tm_entries" => fibdict(res.fib),
        "s0_candidates_per_tm" => fibdict(res.s0fib)]
    write_slice(joinpath(S11DIR, "R2S11_V3_ENGINE.jsonl"), meta, res.entries)
    rep["stats"] = res.stats
    report["V3"] = rep
end

# --- MUT: mutation control — a non-fundamental Pell unit must be CAUGHT ------
# A continued-fraction convergent route to t^2-Du^2=4 can return (18,8) at D=5,
# which is the CUBE of the true fundamental (3,1) (this happens whenever the
# norm -1/-4 branch fires). An engine trusting such a unit under-folds the
# orbit: the V3 gates must FAIL.
function run_mut(report)
    println(stderr, "MUT: V3 with perturbed unit P5 -> P5^3 (Pell (18,8))")
    sub = Dict{String,Any}()
    caught = false; msg = ""
    try
        run_v3(sub; mutate3=true)
    catch err
        caught = true
        msg = sprint(showerror, err)
        length(msg) > 300 && (msg = msg[1:300] * "...")
    end
    report["MUT"] = Dict("mutation" => "unit P5 replaced by P5^3 == non-fundamental Pell (18,8)",
                         "caught" => caught, "first_failure" => msg)
    caught || error("MUTATION CONTROL FAILED: perturbed unit was NOT caught")
    println(stderr, "  MUT: CAUGHT — $msg")
end

# --- main --------------------------------------------------------------------
function main(mode::String)
    t0 = time()
    println(stderr, "s11 mode=$mode (Oscar $(Oscar.VERSION_NUMBER))")
    verbatim_check()
    println(stderr, "verbatim_check: entry_validate == harness.jl  PASS")
    report = Dict{String,Any}("mode" => mode,
        "law" => "entry_validate VERBATIM from harness.jl (source-diffed at startup)",
        "verbatim_entry_validate" => true)
    mode in ("v1", "all") && run_v1(report)
    mode in ("v2", "all") && run_v2(report)
    mode in ("v3", "all") && run_v3(report)
    mode in ("mut", "all") && run_mut(report)
    report["wall_s"] = round(time() - t0, digits=1)
    open(joinpath(S11DIR, "ENGINE_REPORT.json"), "w") do io
        JSON.print(io, report, 1)
    end
    println("S11-ENGINE OK mode=$mode wall=$(report["wall_s"])s")
end

main(isempty(ARGS) ? "all" : ARGS[1])


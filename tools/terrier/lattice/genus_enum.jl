# genus_enum.jl — Aut-free Kneser genus enumeration + PARI POST certificates.
# Avoids Hecke enumerate_definite_genus's default path, whose
# default_invariant_function computes automorphism_group_order per CANDIDATE;
# in Oscar that routes to GAP _stab_via_fin_field -> OutOfMemoryError at rank 19.
# DESIGN (each rule checked against the Hecke source for ZZLat genus
# representatives, module file ZGenusRep under src/QuadForm/Quad, in the Hecke
# release bundled with Oscar):
#  (1) CHEAP invariant = short-vector counts at norms 1..6 per lattice
#      (Hecke._theta_series iterator count; norms 2,4,6 carry the content on
#      even lattices). NO Aut call anywhere during enumeration.
#  (2) algorithm=:random FORCED. :default would pick :orbit at rank 19 (the
#      p-neighbour count estimate (2^19-1)/(2-1) has 6 digits < 7) and :orbit
#      calls automorphism_group_generators per class = the same GAP path.
#  (3) use_mass=false. In Hecke the public wrappers' use_mass=false path is
#      DEAD (_enumerate_definite_genus! sets missing_mass=Ref(0), so its
#      `while missing_mass[]>0` loop body never runs), so we drive
#      Hecke._neighbours directly with the SAME callback/inv_dict semantics as
#      _enumerate_definite_genus!, minus the mass ledger.
#      JUSTIFICATION: mass completeness is NOT lost — it lives in the POST
#      pass, where PARI qfauto gives exact |Aut| per FINAL class and we
#      require sum(1/|Aut|) == Hecke mass EXACTLY (cross-engine completeness
#      certificate).
#  (4) Isometry tests via Hecke.__is_isometric_with_isometry_definite
#      (native Plesken-Souvignier, module file ZLattices; Oscar overrides only
#      the single-underscore wrapper, whose decomposition branch re-enters GAP).
# POST also computes per-class min via PARI qfminim (min > 2 <=> the class is
# root-free, which is the question the production verdict answers).
# Needs: julia 1.10 with Oscar (set TERRIER_JULIA_PROJECT and pass it as
# --project), and PARI/GP's `gp` on PATH for the POST pass.
# Usage:
#   julia genus_enum.jl <det> <key16> [stop_after]         # production, rank 19
#   julia genus_enum.jl --control <det> <key16> <ref_receipt> [stop_after]
#     production: <key16> = first 16 hex of sha256 of the canonical genus
#       symbol genus(discriminant_group(g), (0,19)); the key must have NO
#       positive-definite rank-5 genus with the same discriminant form.
#     control: for a key whose discriminant form DOES carry a positive-definite
#       rank-5 genus C = genus(D,(5,0)), enumerate C with the SAME engine +
#       POST, then require exact match vs the reference receipt <ref_receipt>
#       (an out_<key8>.txt written by an earlier run: CLS| lines + HECKEMASS|):
#       class count + mass + perfect qfisom bijection + per-class |Aut|.
#   Receipts go to $TERRIER_GENUS_OUT (default: this directory): out_<key8>.txt
#   and receipts/{ENUM,CTRL}_<key16>.json.
import JSON, SHA
using Oscar
const OUT = get(ENV, "TERRIER_GENUS_OUT", @__DIR__)
key16(x) = bytes2hex(SHA.sha256(string(x)))[1:16]
const CONTROL = !isempty(ARGS) && ARGS[1] == "--control"
const A = CONTROL ? ARGS[2:end] : ARGS
const TDET = parse(Int, A[1]); const KEY = A[2]
const REFRC = CONTROL ? A[3] : ""       # reference receipt (control mode)
const STOP = length(A) >= (CONTROL ? 4 : 3) ?
    parse(Int, A[CONTROL ? 4 : 3]) : 20000

function find_genus()   # key law: sha256 of the canonical discriminant genus, first 16 hex
    for g in integer_genera((0, 19), TDET; even=true)
        gcan = genus(discriminant_group(g), (0, 19))
        key16(gcan) == KEY && return g
    end
    error("target key not found at det=$TDET")
end

function seed_and_mass(g0)
    D = discriminant_group(g0)
    if CONTROL
        C = genus(D, (5, 0))                # rank-5 definite genus, same disc form (control object)
        L = representative(C)
        @assert genus(L) == C
        return L, mass(C), 5, C
    end
    comp_fail = false
    try genus(D, (5, 0)) catch; comp_fail = true; end
    @assert comp_fail "rank-5 definite genus with this discriminant form exists — production mode expects none (use --control)"
    L = rescale(representative(g0), -1)     # seed: genus representative (local-global), flipped to positive definite
    @assert genus(rescale(L, -1)) == g0
    m = mass(g0)
    @assert mass(genus(L)) == m             # mass invariant under -1
    return L, m, 19, nothing
end

# (1) cheap invariant: sv counts (up to sign) at norms 1..6 — no Aut, no storage
cheap_inv(L) = Tuple(Hecke._theta_series(L, 6))
# (4) native Plesken-Souvignier isometry (bypasses the Oscar/GAP override)
natiso(X, Y) = Hecke.__is_isometric_with_isometry_definite(
    lattice(rational_span(X)), lattice(rational_span(Y)))[1]

# Aut-free genus traversal: same structure as _enumerate_definite_genus!
# (spinor seeds -> cycle classes -> Kneser p-neighbours) minus the mass ledger.
function enumerate_aut_free(L0, t0)
    @assert isone(scale(L0)) "scale != 1 unsupported"
    res = ZZLat[]
    inv_dict = Dict{Any,Vector{ZZLat}}()
    _invariants = M -> begin                # same semantics as Hecke _enumerate_definite_genus!
        for (I, Z) in inv_dict
            any(isequal(M), Z) && return I
        end
        return cheap_inv(M)
    end
    callback = M -> begin                   # same semantics as Hecke _enumerate_definite_genus!
        any(isequal(M), res) && return false
        I = _invariants(M)
        haskey(inv_dict, I) || return true
        return all(N -> !natiso(N, M), inv_dict[I])
    end
    hb(nc) = (println("HEARTBEAT $KEY classes=$nc wall=$(round(time()-t0,digits=1))s");
              flush(stdout))
    for S in Hecke.spinor_genera_in_genus(L0)   # completeness across spinor genera
        (isempty(res) || callback(S)) || continue
        I = _invariants(S)
        haskey(inv_dict, I) ? push!(inv_dict[I], S) : (inv_dict[I] = ZZLat[S])
        push!(res, S)
    end
    p = Hecke.smallest_kneser_prime(L0)
    vain = Ref{Int}(0); i = 0; tlast = time()
    while vain[] <= STOP
        i = (i % length(res)) + 1
        N = Hecke._neighbours(res[i], p; mode=:enumeration, algorithm=:random,
            rand_neigh=10, callback=callback, inv_dict=inv_dict,
            _invariants=_invariants, use_mass=false, vain=vain, stop_after=STOP)
        for M in N                          # inv_dict updated inside _neighbours
            push!(res, M)
            length(res) % 25 == 0 && hb(length(res))
        end
        if time() - tlast > 300             # liveness (hang diagnosis) every 5 min
            println("HEARTBEAT $KEY classes=$(length(res)) vain=$(vain[]) " *
                    "wall=$(round(time()-t0,digits=1))s"); flush(stdout)
            tlast = time()
        end
    end
    return res
end

# POST pass, per FINAL class: exact |Aut| (PARI qfauto) + min (PARI qfminim).
# One gp subprocess per class => natural heartbeat + failure isolation.
function pari_post(grams, n, t0)
    auts = BigInt[]; mins = Int[]; nvecs = Int[]
    scr = joinpath(OUT, "post_$(KEY[1:8]).gp")
    for (i, Gi) in enumerate(grams)
        gs = join([join(Gi[r, :], ",") for r in 1:n], ";")
        open(scr, "w") do io
            println(io, "default(parisizemax, 2^30);")
            println(io, "G=[", gs, "]; v=qfminim(G);")
            println(io, "print(\"POST|\", qfauto(G)[1], \"|\", v[1], \"|\", v[2]);")
            println(io, "quit();")
        end
        s = read(`gp -q -f $scr`, String)
        f = split(only(filter(l -> startswith(l, "POST|"), split(s, '\n'))), '|')
        push!(auts, parse(BigInt, f[2]))
        push!(nvecs, parse(Int, f[3]))      # both signs counted (PARI convention)
        push!(mins, parse(Int, f[4]))
        i % 25 == 0 && (println("HEARTBEAT $KEY post $i/$(length(grams)) " *
            "wall=$(round(time()-t0,digits=1))s"); flush(stdout))
    end
    rm(scr; force=true)
    return auts, mins, nvecs
end

function pari_qfisom(G1, G2, n)
    scr = joinpath(OUT, "isom_tmp.gp")
    a = join([join(G1[r, :], ",") for r in 1:n], ";")
    b = join([join(G2[r, :], ",") for r in 1:n], ";")
    open(scr, "w") do io
        println(io, "print(\"ISO|\", qfisom([", a, "],[", b, "])!=0); quit();")
    end
    return occursin("ISO|1", read(`gp -q -f $scr`, String))
end

# CONTROL gate: the enumerated class list must match the reference receipt EXACTLY —
# same count, same Hecke mass, perfect qfisom bijection, same per-class |Aut|.
function control_check(grams, auts, m, n)
    txt = read(REFRC, String)
    bg = Matrix{Int}[]; bauts = BigInt[]
    for ln in split(txt, '\n')
        startswith(ln, "CLS|") || continue
        f = split(ln, '|')
        rows = split(strip(f[4], ['[', ']']), ";")
        push!(bg, permutedims(hcat([parse.(Int, split(r, ",")) for r in rows]...)))
        push!(bauts, parse(BigInt, f[6]))
    end
    hm = split(only(filter(l -> startswith(l, "HECKEMASS|"), split(txt, '\n'))), '|')[2]
    num, den = split(hm, '/')
    @assert m == ZZ(parse(BigInt, num)) // ZZ(parse(BigInt, den)) "CONTROL: mass != reference"
    @assert length(bg) == length(grams) "CONTROL: class count mismatch"
    pair = zeros(Int, length(grams))
    for (i, Gi) in enumerate(grams), (j, Bj) in enumerate(bg)
        pari_qfisom(Gi, Bj, n) || continue
        @assert pair[i] == 0 "CONTROL: class $i isometric to two reference classes"
        pair[i] = j
        @assert auts[i] == bauts[j] "CONTROL: |Aut| mismatch class $i vs reference $j"
    end
    @assert sort(pair) == collect(1:length(bg)) "CONTROL: no bijection: $pair"
    println("CONTROL-GATE $KEY: bijection $pair, mass + |Aut| receipts match")
    return pair
end

function main()
    t0 = time()
    mkpath(joinpath(OUT, "receipts"))
    out = joinpath(OUT, "receipts", (CONTROL ? "CTRL_" : "ENUM_") * KEY * ".json")
    isfile(out) && (println("RESUME: receipt exists for $KEY, skip"); return)
    g0 = find_genus()
    L0, m, n, C = seed_and_mass(g0)
    println("STATUS $KEY det=$TDET mass=$m rank=$n enum start " *
            "(stop_after=$STOP)"); flush(stdout)
    res = enumerate_aut_free(L0, t0)
    println("STATUS $KEY enum done: $(length(res)) classes, " *
            "wall $(round(time()-t0,digits=1))s; POST start"); flush(stdout)
    grams = Matrix{Int}[]
    for K in res
        G = gram_matrix(K)
        Gi = [Int(ZZ(G[a, b])) for a in 1:n, b in 1:n]
        @assert Gi == permutedims(Gi)                   # symmetric
        @assert all(iseven(Gi[a, a]) for a in 1:n)      # EVENNESS receipt
        @assert TDET == Int(ZZ(Oscar.det(G)))           # det receipt
        push!(grams, Gi)
    end
    auts, mins, nvecs = pari_post(grams, n, t0)
    masssum = sum(1 // ZZ(a) for a in auts; init=QQ(0))
    @assert masssum == m "MASS-CHECK FAIL: $masssum != $m"  # completeness certificate
    pair = CONTROL ? control_check(grams, auts, m, n) : Int[]
    nrootfree = count(>(2), mins)
    verdict = CONTROL ? "CONTROL-PASS" :
              (nrootfree == 0 ? "CLOSED-ALL-ROOTED" : "ROOTFREE-FOUND")
    k8 = KEY[1:8]
    classes = Any[]
    io2 = open(joinpath(OUT, "out_$(k8)$(CONTROL ? "_ctrl" : "").txt"), "w")
    println(io2, "KEY|$KEY|DET|$TDET")
    println(io2, "GENUS|", replace(string(genus(discriminant_group(g0), (0, 19))),
                                   "\n" => " "))
    println(io2, CONTROL ? "COMP|" * replace(string(C), "\n" => " ") :
        "COMP|NONE (NO-RANK5-COMPLEMENT; Nikulin nonexistence receipt)")
    println(io2, "HECKEMASS|", numerator(m), "/", denominator(m))
    println(io2, "NCLASSES|", length(grams))
    for (i, Gi) in enumerate(grams)
        nr = mins[i] == 2 ? nvecs[i] : 0    # qfminim counts both signs (vector-count convention)
        gj = [[Gi[r, c] for c in 1:n] for r in 1:n]
        push!(classes, Dict("i" => i,
            "gram_sha16" => bytes2hex(SHA.sha256(JSON.json(gj)))[1:16],
            "aut_order" => string(auts[i]), "min_norm" => mins[i],
            "n_roots" => nr, "gram" => gj))
        gstr = join([join([string(Gi[r, c]) for c in 1:n], ", ") for r in 1:n], "; ")
        println(io2, "CLS|$i|GRAM|[", gstr, "]|AUT|", auts[i],
                "|MULT|1|MIN|", mins[i], "|NROOTS|", nr)
    end
    println(io2, "MASSSUM|", numerator(masssum), "/", denominator(masssum))
    close(io2)
    data = Dict(
        "schema" => "genus_enum (Aut-free _neighbours + PARI POST)",
        "mode" => CONTROL ? "control" : "production",
        "key" => KEY, "det" => TDET, "rank" => n, "stop_after" => STOP,
        "symbol" => replace(string(g0), "\n" => " "),
        "methods" => ["cheap invariant: sv counts norms 1..6 (no Aut)",
                    "algorithm=:random forced (:orbit needs Aut generators)",
                    "use_mass=false; completeness via POST qfauto mass",
                    "native Plesken-Souvignier isometry (no GAP)"],
        "aut_engine" => "PARI/GP qfauto (gp subprocess, exact)",
        "min_engine" => "PARI/GP qfminim",
        "mass" => string(m), "mass_sum_1_over_aut" => string(masssum),
        "mass_check" => (masssum == m), "n_classes" => length(grams),
        "n_rootfree" => nrootfree, "verdict" => verdict,
        "wall_s" => round(time() - t0, digits=1), "classes" => classes)
    if CONTROL
        data["control"] = Dict("reference_receipt" => REFRC, "pairing" => pair,
            "checks" => "count + Hecke mass + qfisom bijection + per-class |Aut|")
    end
    tmp = out * ".tmp"
    open(tmp, "w") do io; JSON.print(io, data); end
    mv(tmp, out; force=true)
    println("VERDICT $KEY $verdict n_classes=$(length(grams)) " *
            "rootfree=$nrootfree mass_check=true wall=$(round(time()-t0,digits=1))s")
    flush(stdout)
end

main()

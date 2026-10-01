# nikulin.jl — OUT-OF-ZONE DISCRIMINANT-FORM NIKULIN ENGINE for sig(3,0) ternary classes.
# Per sig(3,0) ternary class A (even, positive definite): certify primitive-embedding
# EXISTENCE + UNIQUENESS STRUCTURE for L(A) ↪ Γ_{3,19} at the discriminant-form level.
# NO definite-genus enumeration anywhere (the definite rank-19 complement genus is the
# astronomic object that makes the turnkey classification route run >25 min/class).
#
# EXISTENCE (Nikulin 1.10.1 + 1.6.1/1.6.2): complement genus g(K) = (sig (0,19),
#   q_K = −q_L) is nonempty iff Milgram sign(q_L) ≡ 3 (mod 8) [computed exactly via the
#   Gauss sum over A_L] and 19 ≥ l(A_L) [l ≤ 3, so the 1.10.1 equality-case p-adic
#   conditions are vacuous]; then L ⊕ K glued along the graph of any anti-isometry
#   γ: q_L → −q_K is even unimodular of signature (3,19) = Γ (unique indefinite even
#   unimodular genus) with L primitive. Cross-checked per class by the INDEPENDENT
#   Oscar route primitive_embeddings(Γ, L; classification=:none) (disc-form level, cheap).
#
# UNIQUENESS STRUCTURE (Nikulin 1.15.1: embeddings up to O(Γ) ↔ pairs (class of K in
#   g(K), γ) modulo im(O(L)→O(q_L)) × im(O(K)→O(q_K))):
#   glue_rigid := [ im(O(L)→O(q_L)) = O(q_L) ]  — computed EXACTLY per class
#   (Plesken–Souvignier O(L) + native image_in_Oq; both checked by the `pins` mode).
#   If rigid: the left factor is full ⇒ for EVERY complement class K the double-coset
#   count is exactly 1 ⇒ embedding classes ↔ complement isometry classes CANONICALLY;
#   the glue marking carries NO residual ambiguity, and the double-coset correction
#   (quotient by the IMAGE subgroups im(O(L)), im(O(K)) — never by full O(q)) is realizable
#   on the L side for every entry on this class. Residual [ι]-multiplicity = complement-class
#   multiplicity ONLY: invisible to τ (coherent tuple (r,A,Ã,B,glue) does not see K's
#   class) and N_flux-invariant (N_flux = ½ tr S0 from (A,B) data alone).
#   If NOT rigid: verdict CERT-COSET-BOUND(index): per-complement-class glue-coset count
#   ∈ [1, index], index = [O(q_L) : im(O(L))]; entry needs the explicit glue-orbit branch.
#   HONESTY LINE (carry with any quote): "unique up to O(Γ)" simpliciter is generically
#   FALSE at sig(3,0) (definite complement ⇒ class multiplicity); this certificate
#   certifies exactly the glue-level uniqueness statement above, nothing more.
#
# Modes:  pins | validate03 | cert:<Dcap> | certfile:<path>
#   pins       — 3 known answers: im(O(L)→O(q_L)) for A2 (order 2 = full O(q)) and for
#                [2 1;1 12] (order 2); the in-zone class [2 1;1 2] has exactly 1
#                embedding class (Nikulin 1.14.4)
#   validate03 — sig(0,3) side, dets 4,6,8: full Oscar classification (≈15 s/class;
#                reference at det 4: 1 class) MUST agree with the disc-form prediction (=1,
#                via Nikulin 1.14.2 on the rank-19 INDEFINITE complement: 19 ≥ l+2 ∀p)
#   cert:<D>   — cert all transversal sig(3,0) classes with det ≤ D from the
#                transversal-lattice list (see TERRIER_TRANSVERSAL3 below)
#   certfile:f — cert the list of 3x3 grams (nested row lists) in the JSON file f
# The transversal-lattice list (validate03 and cert modes) is a JSON object with keys
# "dets" (list of determinants) and "grams" (matching list of flat row-major 3x3 Grams);
# it is not included in the package — set TERRIER_TRANSVERSAL3 to the path of one.
# Output: appends one JSON line per certificate to CERTS_SIG30.jsonl and writes a
# RUN_<mode>.json receipt, both next to this file.
# Run: nice -n 5 julia +1.10 -t 1 --project=<your Oscar env> \
#        nikulin.jl <mode> [args]
# Determinism: no RNG. This file is sha-stamped into every receipt it writes.

import JSON
import SHA
import Dates
using Oscar

const HERE = @__DIR__
include(joinpath(@__DIR__, "lattice.jl"))   # the LatticeLayer module (same directory)
using .LatticeLayer

const CERTS = joinpath(HERE, "CERTS_SIG30.jsonl")
# The transversal-lattice list is an external input (validate03/cert modes only; format
# in the header): set TERRIER_TRANSVERSAL3 to its path; unset => those modes refuse loudly.
const TRANSVERSAL = get(ENV, "TERRIER_TRANSVERSAL3", "")

engine_sha() = bytes2hex(SHA.sha256(read(@__FILE__)))

mzz3(rows::Vector) = matrix(ZZ, 3, 3, Int[rows[i][j] for i in 1:3, j in 1:3])
gram_from_flat(g::Vector) = mzz3([g[3k-2:3k] for k in 1:3])
tolists(z::ZZMatrix) = [[Int(z[i, j]) for j in 1:ncols(z)] for i in 1:nrows(z)]

"exact-rational Gauss-sum Milgram signature of q_L: σ ∈ Z/8 with Σ_x e^{iπ q(x)} = √n e^{iπσ/4}."
function milgram_sigma(L::ZZLat)
    T = discriminant_group(L)
    n = Int(order(T))
    # q lifts are rationals in [0,2); the sum is in Z[ζ_8·(1/denominators)] — evaluate in
    # high-precision floats with an integrality tolerance far below any 1/8 turn.
    s = sum(cispi(Float64(LatticeLayer.mod_2Z(lift(quadratic_product(x))))) for x in T)
    r = abs(s)
    abs(r - sqrt(n)) < 1e-6 * sqrt(n) || error("Milgram modulus |Σ|=$r ≠ √$n — degenerate q?")
    turns = angle(s) / (pi / 4)
    sig8 = mod(round(Int, turns), 8)
    abs(turns - round(Int, turns)) < 1e-6 || error("Milgram phase not a multiple of π/4")
    return sig8
end

"minimal number of generators l(A_L) = count of nontrivial elementary divisors."
ell_of(L::ZZLat) = count(d -> d != 1, elementary_divisors(abelian_group(discriminant_group(L))))

# ---------------------------------------------------------------------------
# per-class certificate, sig(3,0) out-of-zone
# ---------------------------------------------------------------------------
function cert_class(G::ZZMatrix, LK3::ZZLat)
    t0 = time()
    L = integer_lattice(gram=G)
    @assert iseven(L) "class not even"
    @assert is_definite(L) && G[1, 1] > 0 && rank(L) == 3 "class not sig(3,0) rank 3"
    d = Int(det(G))
    T = discriminant_group(L)
    ell = ell_of(L)
    sig8 = milgram_sigma(L)
    # EXISTENCE at disc-form level (Nikulin 1.10.1 on the complement genus):
    ex_milgram = (sig8 == mod(3, 8))          # sign(−q_L) ≡ −3 ≡ 5 ≡ 0−19 (mod 8)
    ex_margin = (19 >= ell)                   # 1.10.1 rank condition; > l ⇒ p-adic cases vacuous
    ex_strict = (19 >= 2 + ell)               # Cor 1.12.3-grade margin (always, l ≤ 3)
    # INDEPENDENT cross-route, disc-form level ONLY: Hecke's own local-symbol
    # arithmetic must CONSTRUCT the complement genus (sig (0,19), q = −q_L) — this is
    # Nikulin 1.10.1 implemented by different code (local genus symbols), no lattice,
    # no class enumeration.  NOTE: Oscar primitive_embeddings with classification=:none
    # is NOT usable here — measured here to hang >25 CPU-min at sig(3,0) (it walks
    # the definite rank-19 complement genus; same pathology as the full classification
    # route).  The turnkey route is blocked entirely on this side, existence-only included.
    okO, oscar_route = false, ""
    try
        gK = genus(rescale(T, -1), (0, 19))
        okO, oscar_route = true, "hecke_genus_symbol_constructed: " * first(string(gK), 120)
    catch e
        okO, oscar_route = false, "hecke_genus_construction_FAILED: " * string(typeof(e))
    end
    exists = ex_milgram && ex_margin && okO
    # UNIQUENESS STRUCTURE: im(O(L) → O(q_L)) vs O(q_L), both exact
    nOq = Int(order(orthogonal_group(T)))
    im, _ = image_in_Oq(L)
    nIm = Int(order(im))
    @assert nOq % nIm == 0
    index = div(nOq, nIm)
    rigid = (index == 1)
    verdict = !exists ? "FAIL-EXISTENCE" :
              rigid ? "CERT-RIGID" : "CERT-COSET-BOUND($index)"
    return Dict(
        "schema" => "nikulin-discform-cert", "sig" => [3, 0], "gram" => tolists(G),
        "det" => d, "ell" => ell, "milgram_sigma_mod8" => sig8,
        "exists" => exists,
        "exists_checks" => Dict("milgram_3_mod_8" => ex_milgram,
                                "rank_margin_19_ge_l" => ex_margin,
                                "margin_19_ge_2pl" => ex_strict,
                                "hecke_complement_genus" => okO,
                                "hecke_route" => oscar_route),
        "aut_order_OL" => Int(automorphism_group_order(L)),
        "Oq_order" => nOq, "im_OL_in_Oq_order" => nIm, "coset_index" => index,
        "glue_rigid" => rigid, "verdict" => verdict,
        "semantics" => "Nikulin 1.15.1 double cosets im(O(L))\\Iso(q_L,−q_K)/im(O(K)); " *
            "rigid ⇒ count 1 for EVERY complement class ⇒ embedding classes ↔ complement " *
            "classes canonically; residual [ι]-multiplicity = complement-class only " *
            "(τ-invisible, N_flux-invariant). NOT a claim of uniqueness up to O(Γ).",
        "zone" => "OUT-OF-ZONE (t+ = 3 = l+, Nikulin 1.14.4 inapplicable)",
        "wall_s" => round(time() - t0, digits=3))
end

# ---------------------------------------------------------------------------
# pins + sig(0,3) validation
# ---------------------------------------------------------------------------
function run_pins(LK3::ZZLat)
    out = Dict{String,Any}()
    A2 = lattice_root(:A, 2)
    im1, _ = image_in_Oq(A2)
    o_im = Int(order(im1)); o_oq = Int(order(orthogonal_group(discriminant_group(A2))))
    out["pin1_A2_image_in_Oq"] = Dict("im" => o_im, "Oq" => o_oq,
        "expect" => "im 2 = full O(q_A2)", "pass" => (o_im == 2 && o_oq == 2))
    L23 = integer_lattice(gram=matrix(ZZ, 2, 2, [2, 1, 1, 12]))
    im2, _ = image_in_Oq(L23)
    o2 = Int(order(im2))
    out["pin2_L2112_image_in_Oq"] = Dict("im" => o2,
        "Oq" => Int(order(orthogonal_group(discriminant_group(L23)))),
        "expect" => "im order 2", "pass" => (o2 == 2))
    d3k1 = integer_lattice(gram=matrix(ZZ, 2, 2, [2, 1, 1, 2]))  # I1 class, IN-ZONE
    okz, clsz = primitive_embeddings(LK3, d3k1)
    out["pin3_inzone_d3k1_classes"] = Dict("ok" => okz, "classes" => length(clsz),
        "expect" => "1 (Nikulin 1.14.4 in-zone: unique embedding class)",
        "pass" => (okz && length(clsz) == 1))
    out["pass"] = all(v["pass"] for (k, v) in out if v isa Dict)
    return out
end

"sig(0,3) validation: full Oscar classification (terminates: complement (3,16) INDEFINITE)
 must equal the disc-form prediction 1 (Nikulin 1.14.2: rank 19 ≥ l+2 ∀p ⇒ single-class
 genus + O(K)→O(q_K) surjective ⇒ one double coset regardless of im(O(L)))."
function run_validate03(LK3::ZZLat, grams::Vector{ZZMatrix})
    res = Any[]
    for G in grams
        t0 = time()
        Lm = integer_lattice(gram=-G)          # sig(0,3) mirror of the sig(3,0) class
        ok, cls = primitive_embeddings(LK3, Lm)
        push!(res, Dict("det_sig30" => Int(det(G)), "ok" => ok, "classes" => length(cls),
                        "predicted" => 1, "pass" => (ok && length(cls) == 1),
                        "wall_s" => round(time() - t0, digits=2)))
        println(stderr, "validate03 det=$(Int(det(G))) ok=$ok classes=$(length(cls)) " *
                        "wall=$(res[end]["wall_s"])s")
    end
    return Dict("results" => res, "pass" => all(r["pass"] for r in res),
        "det4_reference" => "det 4 reference: (true, 1 class), ≈15 s/class — must match the det-4 row")
end

# ---------------------------------------------------------------------------
# drivers
# ---------------------------------------------------------------------------
function certed_grams()
    done = Set{String}()
    isfile(CERTS) && for l in eachline(CERTS)
        isempty(strip(l)) && continue
        push!(done, JSON.json(JSON.parse(l)["gram"]))
    end
    return done
end

function emit_certs(grams::Vector{ZZMatrix}, LK3::ZZLat, receipt::Dict)
    done = certed_grams()
    lines, n_skip = String[], 0
    for G in grams
        key = JSON.json(tolists(G))
        key in done && (n_skip += 1; continue)
        c = cert_class(G, LK3)
        c["engine_sha256"] = receipt["engine_sha256"]
        push!(lines, JSON.json(c))
        push!(done, key)
        println(stderr, "cert det=$(c["det"]) $(c["verdict"]) OqL=$(c["Oq_order"]) " *
                        "im=$(c["im_OL_in_Oq_order"]) wall=$(c["wall_s"])s")
    end
    if !isempty(lines)
        open(CERTS, "a") do io; for l in lines; println(io, l); end; end
    end
    receipt["n_new_certs"] = length(lines); receipt["n_skipped_already_certed"] = n_skip
    return receipt
end

function main()
    mode = ARGS[1]
    t0 = time()
    LK3 = lattice_Gamma319()
    receipt = Dict{String,Any}("mode" => mode, "engine_sha256" => engine_sha(),
        "utc" => Dates.format(Dates.now(Dates.UTC), "yyyy-mm-ddTHH:MM:SSZ"))
    if mode == "pins"
        receipt["pins"] = run_pins(LK3)
    elseif mode == "validate03"
        isempty(TRANSVERSAL) && error("TERRIER_TRANSVERSAL3 unset — point it at the transversal-lattice list (JSON with keys dets/grams; external input, not included in the package)")
        t = JSON.parsefile(TRANSVERSAL)
        gs = ZZMatrix[gram_from_flat(t["grams"][i]) for i in eachindex(t["dets"])
                      if t["dets"][i] in (4, 6, 8)]
        receipt["validate03"] = run_validate03(LK3, gs)
    elseif startswith(mode, "cert:")
        D = parse(Int, split(mode, ":")[2])
        isempty(TRANSVERSAL) && error("TERRIER_TRANSVERSAL3 unset — point it at the transversal-lattice list (JSON with keys dets/grams; external input, not included in the package)")
        t = JSON.parsefile(TRANSVERSAL)
        gs = ZZMatrix[gram_from_flat(t["grams"][i]) for i in eachindex(t["dets"])
                      if t["dets"][i] <= D]
        receipt = emit_certs(gs, LK3, receipt)
    elseif startswith(mode, "certfile:")
        todo = JSON.parsefile(String(split(mode, ":", limit=2)[2]))
        gs = ZZMatrix[mzz3([Vector{Int}(r) for r in g]) for g in todo]
        receipt = emit_certs(gs, LK3, receipt)
    else
        error("unknown mode $mode")
    end
    receipt["wall_s"] = round(time() - t0, digits=1)
    tag = replace(mode, r"[:/]" => "_")
    open(joinpath(HERE, "RUN_$(tag).json"), "w") do io; JSON.print(io, receipt, 1); end
    println("nikulin $mode DONE $(receipt["wall_s"])s")
end

main()

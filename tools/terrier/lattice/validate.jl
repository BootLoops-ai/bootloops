# validate.jl — validation battery for lattice.jl (the LatticeLayer module),
# 18 gates.  Run single-threaded:
#   nice -n 5 julia +1.10 --project=<your Oscar env> -t 1 validate.jl
#
# Every gate prints PASS/FAIL individually and the battery writes its report,
# RESULTS_BUILD.md, next to this file (a generated output, regenerated on every
# run).  Battery A holds the known answers that an independent cross-check
# route (PARI/GP + pure python — no code shared with this route; not
# included in the package) reproduces exactly; Battery B covers the remaining
# target gates against published values.  A mismatch in Battery A means one
# of the convention traps C1–C5 (see the header of lattice.jl) has drifted:
# reconcile it before any downstream computation consumes either route.

include(joinpath(@__DIR__, "lattice.jl"))   # the LatticeLayer module (same directory)
using .LatticeLayer
using Oscar
using Printf

const REPORT_PATH = joinpath(@__DIR__, "RESULTS_BUILD.md")   # generated report (output)
const RESULTS = NamedTuple{(:id, :desc, :pass, :secs, :detail),
                           Tuple{String,String,Bool,Float64,String}}[]

function gate!(f::Function, id::String, desc::String)
    t0 = time()
    pass, detail = try
        f()
    catch e
        (false, "EXCEPTION: " * first(sprint(showerror, e), 400))
    end
    secs = time() - t0
    push!(RESULTS, (id=id, desc=desc, pass=pass, secs=secs, detail=detail))
    @printf("%-6s %-52s %-4s %8.2fs  %s\n", id, desc, pass ? "PASS" : "FAIL",
            secs, first(detail, 110))
    flush(stdout)
    return pass
end

t_all = time()
println("="^100)
println("K3 lattice-layer (LatticeLayer) validation battery — ",
        "julia $(VERSION), Oscar $(Oscar.VERSION_NUMBER), threads=$(Threads.nthreads())")
println("="^100)

# shared objects
U     = lattice_U()
E8    = lattice_E8()
E8m   = lattice_E8(neg=true)
A1    = lattice_root(:A, 1)
A1m   = lattice_root(:A, 1, scale=-1)
A2    = lattice_root(:A, 2)
G319  = lattice_Gamma319()
T1    = attractive_T(3)
T2    = attractive_T(4)

# ============ BATTERY A — known answers shared with the cross-check ============

gate!("A1", "class numbers h(D), SL2-reduced, primitive forms") do
    want = Dict(-3 => 1, -4 => 1, -7 => 1, -8 => 1, -23 => 3, -47 => 5, -71 => 7)
    got  = Dict(D => class_number_bqf(D) for D in keys(want))
    (got == want, "h = " * join(("h($D)=$(got[D])" for D in sort(collect(keys(want)), rev=true)), " "))
end

gate!("A2", "reduced forms per disc + T1/T2 Grams") do
    ok = reduced_binary_forms(-3)  == [(1,1,1)] &&
         reduced_binary_forms(-4)  == [(1,0,1)] &&
         reduced_binary_forms(-23) == [(1,1,6),(2,-1,3),(2,1,3)] &&
         reduced_binary_forms(-47) == [(1,1,12),(2,-1,6),(2,1,6),(3,-1,4),(3,1,4)] &&
         reduced_binary_forms(-71) == [(1,1,18),(2,-1,9),(2,1,9),(3,-1,6),(3,1,6),(4,-3,5),(4,3,5)] &&
         gram_matrix(T1) == matrix(QQ, [2 1; 1 2]) &&
         gram_matrix(T2) == matrix(QQ, [2 0; 0 2])
    (ok, "lists for D=-3,-4,-23,-47,-71 exact; T1=[2,1;1,2], T2=[2,0;0,2]")
end

gate!("A3", "E8 theta r(0..8)=1,240,2160,6720,17520") do
    tc  = theta_counts(E8, 8)
    got = [get(tc, k, 0) for k in 0:2:8]
    (got == [1, 240, 2160, 6720, 17520], "got $got (vectors, both signs; ±pairs ×2 per C4)")
end

gate!("A4", "T1⊕T2 theta 1,10,28,30,34 tot 103, FP == brute") do
    L   = direct_sum(T1, T2)[1]
    G   = matrix(ZZ, [ZZ(gram_matrix(L)[i, j]) for i in 1:4, j in 1:4])
    tc  = theta_counts(L, 8)
    bf  = brute_force_theta(G, 8)
    got = [get(tc, k, 0) for k in 0:2:8]
    agree = all(get(tc, k, 0) == get(bf, k, 0) for k in 0:8)
    tot = sum(got)
    (got == [1, 10, 28, 30, 34] && tot == 103 && agree,
     "FP $got tot=$tot; independent box scan agrees on every norm 0..8: $agree")
end

gate!("A5", "U⊕U box counts B=5,10 (indefinite fingerprint, C5)") do
    UU  = direct_sum(U, U)[1]
    G   = matrix(ZZ, [ZZ(gram_matrix(UU)[i, j]) for i in 1:4, j in 1:4])
    b5  = box_theta_counts(G, 5);  g5  = [get(b5,  k, 0) for k in 0:2:8]
    b10 = box_theta_counts(G, 10); g10 = [get(b10, k, 0) for k in 0:2:8]
    (g5 == [833, 308, 412, 368, 492] && g10 == [3905, 1012, 1628, 1360, 1836],
     "B=5: $g5  B=10: $g10")
end

gate!("A6", "Siegel mass II_8 = 1/696729600 = 1/(2^14 3^5 5^2 7)") do
    m  = siegel_mass(E8)
    n  = ZZ(2)^14 * ZZ(3)^5 * ZZ(5)^2 * ZZ(7)
    (m == 1 // n && aut_order(E8) == n && n == ZZ(696729600),
     "mass(E8)=$m, |O(E8)|=$(aut_order(E8)) (Mordell 1938: single class II_8)")
end

# ========================= BATTERY B — target gates ===========================

gate!("B1", "class numbers vs published tables (Buell/Cohen)") do
    # -12,-16 non-fundamental; -15,-20,-24 fundamental; -163 Heegner
    want = Dict(-12 => 1, -15 => 2, -16 => 1, -20 => 2, -24 => 2, -163 => 1)
    got  = Dict(D => class_number_bqf(D) for D in keys(want))
    (got == want, join(("h($D)=$(got[D])" for D in sort(collect(keys(want)), rev=true)), " "))
end

gate!("B2", "h(D) vs Hecke ideal class group (independent path)") do
    fails = String[]
    for D in (-3, -4, -7, -8, -15, -20, -23, -24, -47, -71, -163)
        m = mod(D, 4) == 1 ? D : divexact(D, 4)      # squarefree core
        K, _ = quadratic_field(m)
        hK = order(class_group(maximal_order(K))[1])
        class_number_bqf(D) == hK || push!(fails, "D=$D: bqf=$(class_number_bqf(D)) field=$hK")
    end
    (isempty(fails), isempty(fails) ? "11 fundamental discs agree with ideal class numbers" :
                     join(fails, "; "))
end

gate!("B3", "attractive_T constructor (incl. imprimitive classes)") do
    ok1 = is_isometric(T1, A2)                         # T1 ≅ A2 root lattice
    ok2 = attractive_T_count(12) == 2 && class_number_bqf(-12) == 1   # (2,2,2) imprimitive
    ok3 = attractive_T_count(16) == 2 && class_number_bqf(-16) == 1   # (2,0,2) imprimitive
    L12a, L12b = attractive_T(12, 1), attractive_T(12, 2)
    ok4 = !is_isometric(L12a, L12b) && det(L12a) == 12 && det(L12b) == 12
    ok5, ok6 = false, false
    try attractive_T(5) catch; ok5 = true end          # d ≡ 1 mod 4 must throw
    try attractive_T(6) catch; ok6 = true end          # d ≡ 2 mod 4 must throw
    (ok1 && ok2 && ok3 && ok4 && ok5 && ok6,
     "T1≅A2; d=12,16 have 2 lattice classes vs h=1 (imprimitive real); d≡1,2 mod 4 rejected")
end

gate!("B4", "discriminant groups + q-value conventions (Q/2Z)") do
    okTriv = order(disc_group(E8m)) == 1 && order(disc_group(U)) == 1 &&
             order(disc_group(G319)) == 1
    okA1   = disc_form_values(A1)  == [QQ(0), QQ(1//2)] &&
             elementary_divisors(abelian_group(disc_group(A1))) == [ZZ(2)]
    okA1m  = disc_form_values(A1m) == [QQ(0), QQ(3//2)]
    okA2   = disc_form_values(A2)  == [QQ(0), QQ(2//3), QQ(2//3)]
    (okTriv && okA1 && okA1m && okA2,
     "A(E8(-1))=A(U)=A(Γ3,19)=1; A(A1)=Z/2, q=[1/2] mod 2Z; A(A1(-1)): q=[3/2]; " *
     "A(A2): q=[2/3] (=n/(n+1), Nikulin)")
end

gate!("B5", "genus symbols (Conway–Sloane via Hecke)") do
    s319 = genus_str(G319); sE8m = genus_str(E8m); sU = genus_str(U)
    ok = occursin("II_(3, 19)", s319) && det(G319) == -1 && rank(G319) == 22 &&
         iseven(G319) && signature_tuple(G319) == (3, 0, 19) &&
         occursin("II_(0, 8)", sE8m) && occursin("II_(1, 1)", sU) &&
         occursin("23", genus_str(attractive_T(23)))
    (ok, "Γ3,19: $s319 (det -1, even, rank 22); E8(-1): $sE8m; U: $sU")
end

gate!("B6", "published masses+|O|: A2, D4, II_16 (& aut orders)") do
    D4 = lattice_root(:D, 4)
    okA2 = siegel_mass(A2) == 1//12  && aut_order(A2) == 12
    okD4 = siegel_mass(D4) == 1//1152 && aut_order(D4) == 1152          # |W(F4)|
    E8E8 = direct_sum(E8, E8)[1]
    wE8  = ZZ(696729600)
    m16  = 1 // (2 * wE8^2) + 1 // (ZZ(2)^15 * factorial(ZZ(16)))
    okII16 = siegel_mass(E8E8) == m16 &&
             m16 == QQ(691) // ZZ(277667181515243520000) &&
             aut_order(E8E8) == 2 * wE8^2
    (okA2 && okD4 && okII16,
     "mass(A2)=1/12, mass(D4)=1/1152=1/|W(F4)|, mass(II_16)=691/277667181515243520000")
end

gate!("B7", "genus enum rank 2 + Siegel vs Σ1/|O| (mass_check)") do
    fails = String[]
    for (d, h, ncls) in ((23, 3, 2), (47, 5, 3), (71, 7, 4))
        ok, m, s, n = mass_check(attractive_T(d))
        (ok && m == QQ(h)//4 && n == ncls) ||
            push!(fails, "d=$d: mass=$m (want $h/4), classes=$n (want $ncls)")
    end
    for L in (A2, lattice_root(:A, 3), lattice_root(:D, 4))
        ok, m, s, n = mass_check(L)
        ok || push!(fails, "mass_check($(rank(L))): $m vs $s")
    end
    (isempty(fails),
     isempty(fails) ? "T(23),T(47),T(71): mass=h/4 (3/4,5/4,7/4), classes=2,3,4 (=GL2 cnts); A2,A3,D4 ok" :
     join(fails, "; "))
end

gate!("B8", "Nikulin primitive embeddings — positive controls") do
    U2 = direct_sum(U, U)[1]
    okU  = primitively_embeds(U2, U)
    nU   = length(primitive_embeddings_all(U2, U)[2])
    okA  = primitively_embeds(E8m, A1m)          # A1(-1) ↪ E8(-1)
    okAU = primitively_embeds(U, A1)             # A1 ↪ U (complement <-2>)
    okT1 = primitively_embeds(G319, T1)          # reference T1 ↪ K3 lattice
    okT2 = primitively_embeds(G319, T2)          # reference T2 ↪ K3 lattice
    (okU && nU == 1 && okA && okAU && okT1 && okT2,
     "U↪U^2 (1 class); A1(-1)↪E8(-1); A1↪U; T1↪Γ3,19; T2↪Γ3,19")
end

gate!("B9", "Nikulin primitive embeddings — negative controls") do
    L2  = direct_sum(A1, A1)[1]
    no6 = !primitively_embeds(L2, integer_lattice(gram=matrix(ZZ, 1, 1, [6])))
    no8 = !primitively_embeds(L2, integer_lattice(gram=matrix(ZZ, 1, 1, [8])))
    noSig = try
        !primitively_embeds(E8m, A1)             # positive vector in neg-def: impossible
    catch
        true                                     # a signature error also means "no"
    end
    (no6 && no8 && noSig,
     "<6>!↪A1⊕A1 (norm 6 unrepresented); <8>!↪A1⊕A1 (norm-8 vectors exist but all imprimitive); A1!↪E8(-1)")
end

gate!("B10", "discriminant-form gluing / overlattice enumeration") do
    # (a) A1 ⊕ A1(-1) glued along the diagonal isotropic Z/2 → U
    L11 = direct_sum(A1, A1m)[1]
    iso = isotropic_subgroups_of(disc_group(L11))
    okU = false
    for (S, inj) in iso
        order(S) == 2 || continue
        M = glue_overlattice(L11, [inj(x) for x in gens(S)])
        okU = det(M) == -1 && iseven(M) && signature_tuple(M) == (1, 0, 1)
    end
    # (b) even overlattices of A1^4: exactly {A1^4, D4}
    A14 = direct_sum(A1, A1, A1, A1)[1]
    ovs = even_overlattices(A14)
    okA14 = length(ovs) == 2 && sort([root_count(M) for M in ovs]) == [8, 24] &&
            any(M -> is_isometric(M, lattice_root(:D, 4)), ovs)
    # (c) even overlattices of D8: exactly {D8, E8}  (D8+ = E8, SPLAG §4.8)
    ovs8 = even_overlattices(lattice_root(:D, 8))
    okD8 = length(ovs8) == 2 && any(M -> abs(det(M)) == 1 && root_count(M) == 240 &&
                                         is_isometric(M, E8), ovs8)
    # (d) negative control: D4 has NO nontrivial even overlattice (q≠0 on A(D4)∖0)
    okD4 = length(even_overlattices(lattice_root(:D, 4))) == 1
    (okU && okA14 && okD8 && okD4,
     "A1⊕A1(-1)+glue=U; A1^4 → {A1^4, D4}; D8 → {D8, E8(240 roots)}; D4 → {D4} only")
end

gate!("B11", "Kneser neighbours + II_16 genus {E8⊕E8, D16+}") do
    nb = kneser_neighbours(E8, 2)
    okE8 = !isempty(nb) && all(N -> is_isometric(N, E8), nb)   # single-class genus
    # build D16+ by gluing, then enumerate the II_16 genus from E8⊕E8
    D16 = lattice_root(:D, 16)
    D16p = nothing
    for (S, inj) in isotropic_subgroups_of(disc_group(D16))
        order(S) == 2 || continue
        M = glue_overlattice(D16, [inj(x) for x in gens(S)])
        abs(det(M)) == 1 && (D16p = M; break)
    end
    E8E8 = direct_sum(E8, E8)[1]
    okD16p = D16p !== nothing && root_count(D16p) == 480 &&
             aut_order(D16p) == ZZ(2)^15 * factorial(ZZ(16)) &&
             !is_isometric(D16p, E8E8)
    reps = genus_classes(E8E8)      # Kneser method internally
    hit  = [any(N -> is_isometric(N, R), (E8E8, D16p)) for R in reps]
    okG  = length(reps) == 2 && all(hit) &&
           any(R -> is_isometric(R, D16p), reps) && any(R -> is_isometric(R, E8E8), reps)
    (okE8 && okD16p && okG,
     "nbrs(E8,2)≅E8; D16+ built by glue (480 roots, |O|=2^15·16!); genus(II_16) = {E8⊕E8, D16+} (Witt)")
end

gate!("B12", "Fincke–Pohst vs independent brute force (random Grams)") do
    fails = String[]
    # fixed non-trivial even Gram + two derived 2B'B + one fresh pseudorandom
    cands = Any[matrix(ZZ, [2 1 0; 1 4 1; 0 1 6])]
    for B in (matrix(ZZ, [1 0 1; 2 1 0; 1 1 1]), matrix(ZZ, [1 2 0; 0 1 2; 1 0 3]))
        push!(cands, 2 * transpose(B) * B)
    end
    seed = Int(round(time())) % 10^6
    s = seed
    nxt() = (s = mod(1103515245 * s + 12345, 2^31); mod(s, 7) - 3)  # LCG in -3..3
    Br = matrix(ZZ, [nxt() for i in 1:3, j in 1:3])
    while det(Br) == 0
        Br = matrix(ZZ, [nxt() for i in 1:3, j in 1:3])
    end
    push!(cands, 2 * transpose(Br) * Br)
    for (i, G) in enumerate(cands)
        all(det(G[1:k, 1:k]) > 0 for k in 1:3) || (push!(fails, "G$i not posdef"); continue)
        L  = integer_lattice(gram=G)
        ub = 20
        tc = theta_counts(L, ub)
        bf = brute_force_theta(G, ub)
        all(get(tc, k, 0) == get(bf, k, 0) for k in 0:ub) ||
            push!(fails, "G$i mismatch: FP=$(sort(collect(tc))) brute=$(sort(collect(bf)))")
    end
    (isempty(fails),
     isempty(fails) ? "4 Grams (1 fixed, 2 derived, 1 seeded LCG seed=$seed): full histograms 0..20 agree" :
     join(fails, "; "))
end

# ============================== report ========================================

total = time() - t_all
npass = count(r -> r.pass, RESULTS)
nfail = length(RESULTS) - npass
println("-"^100)
@printf("%d gates: %d PASS, %d FAIL   (total wall %.1fs)\n", length(RESULTS), npass, nfail, total)

open(REPORT_PATH, "w") do io
    println(io, "# Validation report — K3 lattice layer (LatticeLayer), 18 gates")
    println(io)
    println(io, "**Engine:** julia $(VERSION), ",
                "Oscar $(Oscar.VERSION_NUMBER) (Hecke bundled), single-threaded, nice 5.")
    println(io, "**Files:** `lattice.jl` (the LatticeLayer module), ",
                "`validate.jl` (this battery; regenerates this file).")
    println(io, "**Battery A** = exact match against the known answers of an independent ",
                "cross-check route (PARI/GP 2.15.4 + pure python; no shared code; not ",
                "included in the package).  **Battery B** = remaining target gates.")
    println(io)
    println(io, "## Gate table — $(npass)/$(length(RESULTS)) PASS",
                nfail > 0 ? ", $nfail FAIL" : "", "  (total wall $(round(total, digits=1)) s)")
    println(io)
    println(io, "| gate | description | status | wall (s) | detail |")
    println(io, "|---|---|---|---:|---|")
    for r in RESULTS
        println(io, "| ", r.id, " | ", replace(r.desc, "|" => "\\|"), " | ",
                r.pass ? "**PASS**" : "**FAIL**",
                " | ", @sprintf("%.2f", r.secs), " | ", replace(r.detail, "|" => "\\|"), " |")
    end
    println(io)
    println(io, """
## Conventions (frozen; traps C1–C5 shared with the independent cross-check)

- **Sign/positivity:** root lattices positive definite (A1 = [2], roots norm +2);
  E8(−1) = rescale(E8, −1); U = [0 1; 1 0]; Γ₃,₁₉ = U³ ⊕ E8(−1)², even unimodular,
  signature (3,19), det −1, genus symbol II_(3,19).
- **Discriminant form:** A_L = L*/L with q: A_L → **Q/2Z** (Nikulin), b: A_L×A_L → Q/Z.
  Fixture: A(A1) = Z/2 with q = [1/2] mod 2Z; A(A1(−1)) has q = [3/2] = [−1/2];
  A(A2) has q-values {0, 2/3, 2/3} (= n/(n+1) on generators).  Values quoted as
  lifts in [0, 2).
- **Forms ↔ lattices (factor-2 trap, C1):** form (a,b,c) ↦ Gram [2a b; b 2c];
  lattice det d = −disc(form).  Class counts are **SL2(Z)** (proper) classes (C2);
  `class_number_bqf` = primitive classes (table h(D)); `attractive_T(d,k)` runs over
  ALL reduced forms of disc −d (imprimitive included, per C3), ordered
  lexicographically in (a,b,c) — deterministic constructor, k = class index.
- **Vector counts (C4):** all counts are VECTORS (±v both), norm 0 listed separately;
  Hecke `short_vectors` returns one per ±pair and is doubled in `theta_counts`.
- **Indefinite fingerprints (C5):** U⊕U counts use the coefficient box |x_i| ≤ B in
  the standard basis of Q = 2(x₁x₂ + x₃x₄) — mandatory box convention.
- **Primitive embeddings:** `primitive_embeddings(L, M)` = embeddings of M **into** L
  (big first); existence via Nikulin discriminant-form criterion
  (`classification=:none`).
- **Overlattices:** even overlattices of even L ↔ isotropic subgroups of (A_L, q)
  (q ≡ 0 checked on every element), Nikulin Prop. 1.4.1.

## Published anchors used

- h(D) tables: Buell, *Binary Quadratic Forms* (1989); Cohen GTM 138; Heegner
  discriminants h=1.  Independently re-derived here via Hecke ideal class groups (B2).
- E8 theta = E₄: 240σ₃ ⇒ 240, 2160, 6720, 17520 (Serre GTM 7 VII §6.6).
- mass(II₈) = 1/696729600, single class (Mordell 1938; Siegel 1935; SPLAG ch. 16).
- |O(A2)| = 12, |O(D4)| = |W(F4)| = 1152, |W(E8)| = 696729600 = 2¹⁴·3⁵·5²·7.
- II₁₆ = {E8⊕E8, D16⁺} (Witt 1941; SPLAG ch. 16), masses 1/(2·|W(E8)|²) + 1/(2¹⁵·16!)
  = 691/277667181515243520000.
- D8⁺ = E8 glue construction (SPLAG ch. 4 §7.3/8.1).
- Nikulin: *Integer symmetric bilinear forms* (1979): Prop. 1.4.1 (overlattices ↔
  isotropic subgroups), §1.12–1.15 (primitive-embedding criterion).

## Notes / caveats

- `mass_check` compares Hecke's analytic Siegel mass against Σ 1/|O| over the
  Kneser-enumerated genus — two independent code paths inside Hecke, anchored to
  published values in A6/B6/B7.
- The B12 brute-force counter is Fincke–Pohst-free: exact dual-diagonal box bound
  x_i² ≤ ub·(G⁻¹)_ii in exact rational arithmetic (no float, no margin heuristic);
  one Gram per run is freshly pseudorandom (LCG seed logged in the gate detail).
- First-call JIT compilation is included in the wall times of the first gate that
  touches each Hecke code path (A3: short vectors; B3: isometry; B7: genus
  enumeration; B8: primitive embeddings) — re-runs are much faster.
- Mismatch protocol (Battery A): stop; reconcile the convention traps C1–C5
  before any downstream computation consumes either route.
""")
end
println("report written to ", REPORT_PATH)
exit(nfail == 0 ? 0 : 1)

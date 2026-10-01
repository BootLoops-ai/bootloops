# lattice.jl — the LatticeLayer module: even-lattice basics, binary quadratic
# forms ↔ attractive-K3 transcendental lattices, Nikulin primitive-embedding /
# gluing / overlattice machinery, genus enumeration + Siegel mass + Kneser
# neighbours, and short-vector counters (Fincke–Pohst plus an independent
# exact brute force), all on top of Oscar/Hecke.  The validation battery is
# validate.jl (18 gates; selftest_lattice.jl is its wrapper).
# Load single-threaded:
#   nice -n 5 julia +1.10 --project=<your Oscar env> -t 1
#
# ============================ CONVENTIONS (frozen) ============================
# C1–C5 below label the five convention traps on which this module and an
# independent cross-check route (PARI/GP + pure python, sharing no code with
# this module; not included in the package) must agree; validate.jl Battery A
# holds the known answers both routes reproduce.
#
# L1. Root lattices are POSITIVE definite: A1 = [2] (roots have norm +2),
#     E8 = the positive-definite even unimodular rank-8 lattice.  E8(-1),
#     A1(-1), ... denote rescale(·, -1).  U = hyperbolic plane, Gram [0 1; 1 0].
#     Gamma_{3,19} = U^3 ⊕ E8(-1)^2: even unimodular, signature (3,19), det -1
#     (the K3 lattice).  Hecke genus symbol II_(3, 19).
#
# L2. Discriminant group A_L = L*/L carries the BILINEAR form b: A_L×A_L → Q/Z
#     and, for L even, the QUADRATIC form q: A_L → Q/2Z, q(x) = x·x mod 2Z
#     (Nikulin's convention; Hecke's `discriminant_group` returns exactly this).
#     Sign fixture: A(A1) = Z/2 with q = [1/2] mod 2Z;  A(A1(-1)) = Z/2 with
#     q = [-1/2] = [3/2] mod 2Z.  `disc_form_values` returns lifts in [0, 2).
#
# L3. Binary forms ↔ even binary lattices (C1, THE factor-2 trap):
#     form f = (a,b,c) := a x² + b xy + c y²  ↦  lattice Gram G = [2a b; b 2c].
#     Lattice norm form = 2f;  det G = 4ac − b² = −disc(f);  lattice det d ↔
#     form disc D = −d.  Reduction is SL2(Z) (proper) reduction:
#     −a < b ≤ a ≤ c, and b ≥ 0 if a = c  (C2: SL2, not GL2).
#     `class_number_bqf(D)` counts PRIMITIVE reduced forms (= h(D) of the
#     tables); `attractive_T(d,k)` enumerates ALL reduced forms of disc −d,
#     including imprimitive ones — attractive-K3 transcendental lattices are
#     ALL even positive-definite binary lattices (Shioda–Inose), and for
#     non-fundamental d the imprimitive classes are real (C3).
#     Class ordering: lexicographic in (a, b, c) — deterministic constructor.
#
# L4. Short vectors (C4): Hecke's `short_vectors(L, ub)` returns ONE
#     representative per ±pair with 0 < norm ≤ ub.  `theta_counts` converts to
#     VECTOR counts (×2) and lists norm 0 separately with count 1.
#     Indefinite lattices need a box convention (C5):
#     `box_theta_counts(G, B)` counts vectors with |x_i| ≤ B in the GIVEN basis.
#
# L5. `primitive_embeddings(L, M)` (Oscar) = primitive embeddings OF M INTO L
#     (big lattice first).  Existence is Nikulin's criterion via discriminant
#     forms; `primitively_embeds` asks existence only (classification=:none).
#
# L6. Overlattices: even overlattices of even L ↔ ISOTROPIC subgroups
#     (q ≡ 0 in Q/2Z, checked on every element) of A_L  [Nikulin Prop. 1.4.1].
# ==============================================================================

module LatticeLayer

using Oscar

export
    # 1. even-lattice basics
    lattice_U, lattice_E8, lattice_root, lattice_Gamma319,
    disc_group, disc_form_values, genus_str,
    # 2. binary forms / attractive K3 transcendental lattices
    reduced_binary_forms, class_number_bqf,
    attractive_T, attractive_T_count, form_gram,
    # 3. Nikulin machinery
    primitively_embeds, primitive_embeddings_all,
    isotropic_subgroups_of, glue_overlattice, even_overlattices,
    # 4. genus / mass / neighbours
    siegel_mass, genus_classes, aut_order, mass_check, kneser_neighbours,
    # 5. short vectors
    short_vecs, theta_counts, box_theta_counts, brute_force_theta, root_count

# ---------------------------------------------------------------------------
# 1. Even-lattice basics
# ---------------------------------------------------------------------------

"""U = hyperbolic plane, Gram [0 1; 1 0]."""
lattice_U() = hyperbolic_plane_lattice()

"""E8 root lattice; `neg=true` gives E8(-1)."""
lattice_E8(; neg::Bool=false) =
    neg ? rescale(root_lattice(:E, 8), -1) : root_lattice(:E, 8)

"""Root lattice `S_n` (S ∈ :A,:D,:E), rescaled by `scale` (POSITIVE-definite
base convention, roots norm +2; scale=-1 for the negative-definite twin)."""
function lattice_root(S::Symbol, n::Int; scale::Int=1)
    L = root_lattice(S, n)
    return scale == 1 ? L : rescale(L, scale)
end

"""Gamma_{3,19} = U^3 ⊕ E8(-1)^2 — the K3 lattice; even unimodular (3,19)."""
lattice_Gamma319() = direct_sum(lattice_U(), lattice_U(), lattice_U(),
                                lattice_E8(neg=true), lattice_E8(neg=true))[1]

"""Discriminant group A_L = L*/L as a Hecke TorQuadModule (q into Q/2Z)."""
disc_group(L::ZZLat) = discriminant_group(L)

"""Reduce a rational into the fundamental domain [0, 2) of Q/2Z."""
mod_2Z(r::QQFieldElem) = r - 2 * floor(ZZRingElem, r // 2)

"""Sorted lifts in [0,2) of q(x), x ∈ A_L (small groups only)."""
function disc_form_values(L::ZZLat; max_group_order::Int=100_000)
    T = discriminant_group(L)
    @req order(T) <= max_group_order "discriminant group too large"
    vals = QQFieldElem[mod_2Z(lift(quadratic_product(x))) for x in T]
    return sort!(vals)
end

"""Genus symbol of L as a display string (Conway–Sloane symbol via Hecke)."""
genus_str(L::ZZLat) = string(genus(L))

# ---------------------------------------------------------------------------
# 2. Binary quadratic forms ↔ attractive K3 transcendental lattices
# ---------------------------------------------------------------------------

"""
    reduced_binary_forms(D; primitive_only=false) -> Vector{NTuple{3,Int}}

All SL2(Z)-reduced positive-definite binary forms (a,b,c) of discriminant
D = b² − 4ac < 0:   −a < b ≤ a ≤ c,  b ≥ 0 if a = c.   Sorted lexicographically.
`primitive_only=true` restricts to gcd(a,b,c)=1 (table class number h(D)).
"""
function reduced_binary_forms(D::Union{Integer,ZZRingElem}; primitive_only::Bool=false)
    Di = Int(D)
    @req Di < 0 "discriminant must be negative"
    @req mod(Di, 4) in (0, 1) "discriminant must be ≡ 0, 1 mod 4"
    @req Di > -10^9 "discriminant too large for Int enumeration"
    forms = NTuple{3,Int}[]
    amax = isqrt(div(-Di, 3))
    for a in 1:amax, b in (-a):a
        mod(Di - b, 2) == 0 || continue          # b ≡ D (mod 2)
        num = b^2 - Di
        mod(num, 4a) == 0 || continue
        c = div(num, 4a)
        c >= a || continue
        (b < 0 && (b == -a || a == c)) && continue   # SL2 reduction boundary
        primitive_only && gcd(gcd(a, b), c) != 1 && continue
        push!(forms, (a, b, c))
    end
    return sort!(forms)
end

"""Class number h(D): number of primitive SL2(Z)-reduced forms of disc D."""
class_number_bqf(D::Union{Integer,ZZRingElem}) =
    length(reduced_binary_forms(D; primitive_only=true))

"""Gram matrix [2a b; b 2c] of the even lattice attached to the form (a,b,c)."""
form_gram(a::Int, b::Int, c::Int) = matrix(ZZ, 2, 2, [2a, b, b, 2c])

"""Number of even positive-definite binary lattices of determinant d
(= ALL reduced forms of disc −d, imprimitive included)."""
attractive_T_count(d::Union{Integer,ZZRingElem}) = length(reduced_binary_forms(-Int(d)))

"""
    attractive_T(d, k=1) -> ZZLat

k-th (lexicographic in (a,b,c)) even positive-definite binary lattice of
determinant d > 0, d ≡ 0 or 3 mod 4 — the transcendental lattice of an
attractive K3.  attractive_T(3,1) = [2 1; 1 2] (= reference T1),
attractive_T(4,1) = [2 0; 0 2] (= reference T2).
"""
function attractive_T(d::Union{Integer,ZZRingElem}, k::Int=1)
    di = Int(d)
    @req di > 0 "determinant must be positive"
    @req mod(di, 4) in (0, 3) "no even binary lattice: need d ≡ 0, 3 mod 4"
    forms = reduced_binary_forms(-di)
    @req 1 <= k <= length(forms) "class index k out of range 1:$(length(forms))"
    (a, b, c) = forms[k]
    return integer_lattice(gram=form_gram(a, b, c))
end

# ---------------------------------------------------------------------------
# 3. Nikulin machinery: primitive embeddings, gluing, overlattices
# ---------------------------------------------------------------------------

"""Existence of a primitive embedding M ↪ L (Nikulin criterion; no
classification).  NOTE ORDER: big lattice L first, as in Oscar."""
function primitively_embeds(L::ZZLat, M::ZZLat)
    ok, _ = primitive_embeddings(L, M; classification=:none)
    return ok
end

"""All primitive embeddings M ↪ L up to equivalence: (ok, [(L, M', C'), ...])."""
primitive_embeddings_all(L::ZZLat, M::ZZLat) = primitive_embeddings(L, M)

"""
    isotropic_subgroups_of(T::TorQuadModule) -> Vector{(S, inj)}

All subgroups S ≤ T with q|_S ≡ 0 in Q/2Z — checked on EVERY element of S,
not only generators.  (Totally isotropic in the quadratic sense.)
"""
function isotropic_subgroups_of(T::TorQuadModule; max_group_order::Int=10^4)
    @req order(T) <= max_group_order "discriminant group too large to enumerate"
    out = Tuple{TorQuadModule,TorQuadModuleMap}[]
    for (S, inj) in subgroups(T)
        isot = true
        for x in S
            if !iszero(quadratic_product(inj(x)))
                isot = false
                break
            end
        end
        isot && push!(out, (S, inj))
    end
    return out
end

"""Overlattice of L glued along the given elements of A_L (must generate an
isotropic subgroup; Hecke checks)."""
glue_overlattice(L::ZZLat, glue::Vector{TorQuadModuleElem}) =
    isempty(glue) ? L : overlattice(L, glue)

"""
    even_overlattices(L) -> Vector{ZZLat}

All even overlattices of even L, one per isotropic subgroup of A_L
[Nikulin Prop. 1.4.1]; deduplicated up to isometry when L is definite
(`up_to_isometry=false` returns one lattice per isotropic subgroup).
"""
function even_overlattices(L::ZZLat; up_to_isometry::Bool=true,
                           max_group_order::Int=10^4)
    @req iseven(L) "L must be even"
    T = discriminant_group(L)
    res = ZZLat[]
    for (S, inj) in isotropic_subgroups_of(T; max_group_order=max_group_order)
        gg = TorQuadModuleElem[inj(x) for x in gens(S)]
        push!(res, glue_overlattice(L, gg))
    end
    if up_to_isometry && is_definite(L)
        uniq = ZZLat[]
        for M in res
            any(N -> is_isometric(N, M), uniq) || push!(uniq, M)
        end
        return uniq
    end
    return res
end

# ---------------------------------------------------------------------------
# 4. Genus enumeration, Siegel mass, Kneser neighbours
# ---------------------------------------------------------------------------

"""Siegel mass of the genus of L (Hecke local-density product)."""
siegel_mass(L::ZZLat) = mass(L)

"""Representatives of all isometry classes in the genus of L (Kneser method)."""
genus_classes(L::ZZLat) = genus_representatives(L)

"""|O(L)| for definite L (Plesken–Souvignier)."""
aut_order(L::ZZLat) = automorphism_group_order(L)

"""
    mass_check(L) -> (ok, analytic_mass, enumerated_mass, nclasses)

Cross-check the analytic Siegel mass against Σ 1/|O(rep)| over the enumerated
genus — two independent code paths (local densities vs neighbours + PS).
"""
function mass_check(L::ZZLat)
    m = mass(L)
    reps = genus_representatives(L)
    s = sum(QQ(1) // automorphism_group_order(N) for N in reps)
    return (s == m, m, s, length(reps))
end

"""Kneser p-neighbours of L (Hecke)."""
kneser_neighbours(L::ZZLat, p::Union{Integer,ZZRingElem}) = Hecke.neighbours(L, p)

# ---------------------------------------------------------------------------
# 5. Short vectors (Fincke–Pohst) + independent brute force
# ---------------------------------------------------------------------------

"""Hecke Fincke–Pohst: one representative per ±pair, 0 < norm ≤ ub."""
short_vecs(L::ZZLat, ub) = short_vectors(L, ub)

"""
    theta_counts(L, ub) -> Dict{Int,Int}

VECTOR counts by norm (both signs counted; norm 0 ↦ 1), norms ≤ ub.
Positive-definite L only.
"""
function theta_counts(L::ZZLat, ub::Int)
    @req is_definite(L) "theta_counts needs a definite lattice"
    counts = Dict{Int,Int}(0 => 1)
    for (v, n) in short_vectors(L, ub)
        ni = Int(ZZ(n))
        counts[ni] = get(counts, ni, 0) + 2
    end
    return counts
end

"""Number of roots (norm-2 vectors, both signs) of L."""
root_count(L::ZZLat) = get(theta_counts(L, 2), 2, 0)

_form_eval(g::Matrix{Int}, x) = begin
    n = size(g, 1)
    s = 0
    @inbounds for i in 1:n
        t = 0
        for j in 1:n
            t += g[i, j] * x[j]
        end
        s += x[i] * t
    end
    s
end

"""
    box_theta_counts(G, B) -> Dict{Int,Int}

Exhaustive count of x ∈ Z^n with |x_i| ≤ B in the GIVEN basis, by norm x'Gx.
Works for indefinite G — this is the C5 indefinite-fingerprint convention
(validate.jl gate A5).
Includes the zero vector (norm 0).
"""
function box_theta_counts(G::ZZMatrix, B::Int)
    n = nrows(G)
    @req n * log(2B + 1) < 18 "box too large (> ~7e7 points)"
    g = Int.(Matrix(G))
    counts = Dict{Int,Int}()
    for x in Iterators.product(ntuple(_ -> -B:B, n)...)
        s = _form_eval(g, x)
        counts[s] = get(counts, s, 0) + 1
    end
    return counts
end

"""
    brute_force_theta(G, ub) -> Dict{Int,Int}

Independent (non-Fincke–Pohst) exact count of ALL vectors with 0 ≤ x'Gx ≤ ub
for POSITIVE-definite G, via the exact dual bound x_i² ≤ ub·(G⁻¹)_ii
(no float, no margin heuristics).  Comparable key-by-key to `theta_counts`.
"""
function brute_force_theta(G::ZZMatrix, ub::Int)
    n = nrows(G)
    Gi = inv(matrix(QQ, G))
    R = Int[Int(isqrt(floor(ZZRingElem, QQ(ub) * Gi[i, i]))) for i in 1:n]
    @req sum(log(2r + 1) for r in R) < 18 "box too large (> ~7e7 points)"
    g = Int.(Matrix(G))
    counts = Dict{Int,Int}()
    for x in Iterators.product(((-R[i]):(R[i]) for i in 1:n)...)
        s = _form_eval(g, x)
        0 <= s <= ub && (counts[s] = get(counts, s, 0) + 1)
    end
    return counts
end

end # module LatticeLayer

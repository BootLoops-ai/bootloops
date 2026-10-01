# predicates.jl — exact predicate chain for K3xK3 fluxes (module SweepPredicates).
# check_flux runs the predicates P1..P7 below in order, fail-fast; side_analysis,
# perp_root_check and congruence_diag are also used stand-alone by harness.jl
# and s11.jl (the PREDICATE_HOOKS defaults).
#
# Overriding rule: every accept/reject decision is EXACT integer / rational /
# number-field arithmetic. No float decides anything here.
#
# Two-lattice convention:
#   flux (N, m, m̃): N ∈ Z^{n1×n2} middle block, (m, m̃) corner charges on
#   Z·(x⊗1) ⊕ Z·(1⊗x̃) (hyperbolic U pairing);
#   g = N d2,  g̃ = N^T d1  (mutually adjoint),
#   S = g g̃ = N d2 N^T d1  on side 1;   S̃ = g̃ g = N^T d1 N d2  on side 2;
#   N_flux = ½ tr(S) + m·m̃   (the flux contribution to the D3 tadpole on K3xK3 in
#   the conventions of arXiv:2010.10519, plus the corner term); budget N_flux ≤ 24,
#   i.e. ‖G‖² = 2 N_flux ≤ 48.
#   Control mode (reduced lattices, arXiv:2103.03250): d1 = d2 = d, m = m̃ = 0.
#   Convention trap: swapping g/g̃ transposes N — pinned here.
#
# Named conditions used below (all checked exactly):
#   (C1)  S is R-diagonalizable with non-negative spectrum: minpoly squarefree,
#         every irreducible factor totally real, every eigenvalue totally ≥ 0;
#   (C0)  ker g̃ = ker S and ker g = ker S̃, checked as rank(N) = rank(S) = rank(S̃);
#   (V)   existence of self-dual sublattice pairs (Σ, Σ̃) with
#         G_mid ∈ (Σ⊗Σ̃) ⊕ (Σ^⊥⊗Σ̃^⊥); self-dual existence lemma: (V) ⟺ (C1) ∧ (C0);
#   corner clause: the corner block is self-dualizable ⟺ (m, m̃) = (0,0) or m·m̃ ≥ 1;
#   stabilization criterion (eigenspace definiteness): every eigenspace of S is
#         d1-definite and every eigenspace of S̃ is d2-definite, at every real place;
#   R(Σ): no root δ (δ² ∈ root_norms) orthogonal to all positive-norm eigenspaces.

module SweepPredicates

using Oscar

const PREDICATE_SET_VERSION = "predicates-1.0"
const CONVENTIONS = "g = N d2, g~ = N^T d1, S = g g~ (side 1), S~ = g~ g (side 2); N_flux = tr(S)/2 + m*m~ <= N_cap; corner (m,m~) = (0,0) or m*m~ >= 1; exact arithmetic only"

# Stable predicate IDs in chain order (fail-fast) and what each one checks.
const PREDICATE_MANIFEST = [
    ("P1_integrality",
     "G4 ∈ H4(X, Z) exactly, no half-integral shift on K3xK3: integer (N, m, m̃) by type; evenness of Q = tr(S) = 2(N_flux − m·m̃) is the even-lattice certificate (the corner U-block contributes 2m·m̃, even for every integer pair, so the certificate lives on the middle block)."),
    ("P2_selfdual_membership",
     "Self-dual membership (V): G_mid ∈ (Σ⊗Σ̃) ⊕ (Σ^⊥⊗Σ̃^⊥) for some pair (Σ, Σ̃), plus the corner clause [m = m̃ = 0] or [m·m̃ ≥ 1]. Implemented through the exact equivalence (V) ⟺ (C1) ∧ (C0): (C1) = minpoly squarefree + every irreducible factor totally real + totally non-negative spectrum, BOTH sides; (C0) = rank(N) = rank(S) = rank(S̃) over Z. Mixed-sign or single-corner charges admit no Minkowski vacuum and are rejected. Membership of the enumeration basis in the declared Künneth block holds by construction; completeness of the glue LIST is the business of harness.jl, not of this predicate."),
    ("P3_primitivity",
     "No lattice-primitivity or rank condition on N is imposed: TAGS the gcd-content of N but NEVER rejects; kept as a chain slot."),
    ("P4_orthogonality",
     "Stratum filter, NOT a condition of the theorem: in mode=pair_block a flux whose rank(N) differs from the declared block rank rank_req is rejected (guards against inflation of a block stratum by lower-rank fluxes); tag-only in mode=control_full. Every P4 rejection carries counts_toward_theorem=false and the stratum-filter marker, so no theorem-coverage count can include it."),
    ("P5_tadpole",
     "n_M2 ∈ Z≥0 (no anti-M2 branes), hence N_flux = ½ tr(S) + m·m̃ with 0 ≤ N_flux ≤ N_cap (default 24). The ½ is re-derived here from the trace route, independent of any enumerator bound (the factor-2 pitfall). With (C1) certified at P2, tr(S) ≥ 0 is automatic, so the lower bound is a certificate; the live filter is the upper budget."),
    ("P6_stabilization",
     "Eigenspace definiteness: every eigenspace of S is d1-definite and every eigenspace of S̃ is d2-definite, together with (C1) — exact: squarefree minpoly, per-irreducible-factor totally-real field, totally-positive eigenvalue, eigenspace Gram congruence-diagonalized over the eigenvalue field, definiteness at EVERY real place, λ=0 kernels on BOTH factors. (C1) failures fire at P2 upstream (chain order); P6 re-verifies them as a certificate and OWNS the definiteness content. Corner-only fluxes (N = 0, m·m̃ ≥ 1) are rejected here (ker S = full space, indefinite on any indefinite lattice): corner-only fluxes are never fully stabilized."),
    ("P7_no_root",
     "R(Σ): no root δ with δ² ∈ root_norms orthogonal to all positive-norm eigenvectors of g g̃ (side 1) or of g̃ g (side 2). root_norms default is derived from mode (control_full → [2,-2]; pair_block → [-2], the K3 (−2)-root condition); an explicit root_norms argument overrides. Exact: W = ⊕ eigenspaces positive-definite at ≥1 real place; Λ_perp = saturated integer kernel of the rational orthogonality system; short-vector check on the definite Λ_perp; unresolved is NEVER an accept. Semantics: the predicate tested is \"no SINGLE root orthogonal to all positive-capable eigenspaces\" — this can only over-accept, which is safe because survivors go on to exact certification."),
]

export check_flux, PREDICATE_SET_VERSION, PREDICATE_MANIFEST, predicate_manifest_dict

predicate_manifest_dict() = Dict("version" => PREDICATE_SET_VERSION,
    "conventions" => CONVENTIONS,
    "predicates" => [Dict("id" => p, "implements" => q) for (p, q) in PREDICATE_MANIFEST])

# -----------------------------------------------------------------------------
# Congruence diagonalization of a symmetric matrix over a field (exact).
# Returns the diagonal entries of some B^T G B, B invertible. Zeros <=> G degenerate.
# -----------------------------------------------------------------------------
function congruence_diag(G0::MatElem)
    F = base_ring(G0); n = nrows(G0)
    G = deepcopy(G0)
    diags = elem_type(F)[]
    idx = collect(1:n)
    while !isempty(idx)
        p = findfirst(i -> !iszero(G[i, i]), idx)
        if p === nothing
            pr = nothing
            for i in idx, j in idx
                if j != i && !iszero(G[i, j]); pr = (i, j); break; end
            end
            pr === nothing && (append!(diags, [zero(F) for _ in idx]); break)
            (i, j) = pr                      # e_i <- e_i + e_j (char 0: 2G[i,j] != 0)
            for k in 1:n; G[i, k] += G[j, k]; end
            for k in 1:n; G[k, i] += G[k, j]; end
            continue
        end
        i = idx[p]; d = G[i, i]; push!(diags, d)
        rest = [j for j in idx if j != i]
        for r in rest
            c = G[r, i] // d
            if !iszero(c)
                for k in 1:n; G[r, k] -= c * G[i, k]; end
                for k in 1:n; G[k, r] -= c * G[k, i]; end
            end
        end
        idx = rest
    end
    return diags
end

# -----------------------------------------------------------------------------
# c1_check: exact (C1) on one side — S R-diagonalizable with eigenvalues ≥ 0.
# minpoly squarefree ⟺ diagonalizable /R once all roots are real; roots real
# and ≥ 0 ⟺ every irreducible factor is x, or x−λ with λ ≥ 0, or degree ≥ 2
# totally real with totally positive root. Used by P2 ((V) ⟺ (C1) ∧ (C0)).
# -----------------------------------------------------------------------------
function c1_check(S::QQMatrix)
    Rx, x = polynomial_ring(QQ, "x")
    mp = minpoly(Rx, S)
    is_squarefree(mp) ||
        return (false, "minpoly not squarefree: S not R-diagonalizable ((C1) fails)")
    for (f, _) in factor(mp)
        if degree(f) == 1
            lam = -coeff(f, 0) // coeff(f, 1)
            lam < 0 && return (false, "negative rational eigenvalue $lam ((C1) fails)")
        else
            K, a = number_field(f, "a", cached=false)
            is_totally_real(K) ||
                return (false, "factor $(f): non-real eigenvalues ((C1) fails)")
            is_totally_positive(a) ||
                return (false, "factor $(f): negative eigenvalue at some real place ((C1) fails)")
        end
    end
    return (true, "")
end

# -----------------------------------------------------------------------------
# side_analysis: exact spectral analysis of one side (S w.r.t. Gram d).
# Implements eigenspace definiteness + (C1) on that side. Returns factor records
# used by P6 (verdict) and P7 (positive-eigenspace orthogonality system).
# place_signs per factor: one of :pos/:neg/:indef per REAL PLACE of its field
# (sign may differ per place; :indef anywhere, a zero diag, a non-real root or
# a negative eigenvalue fails P6).
# Note: ¬(C1) inputs are filtered at P2, so the C1 branches below are
# certificates (bug traps), not live filters; P6 owns the definiteness content.
# -----------------------------------------------------------------------------
function side_analysis(S::QQMatrix, d::QQMatrix)
    n = nrows(S)
    Rx, x = polynomial_ring(QQ, "x")
    cp = charpoly(Rx, S)
    mp = minpoly(Rx, S)
    sqfree = is_squarefree(mp)
    ok = sqfree
    reasons = String[]
    sqfree || push!(reasons, "minpoly not squarefree: S not diagonalizable (C1 fails)")
    factors = NamedTuple[]
    for (f, m) in factor(cp)
        if degree(f) == 1
            lam = -coeff(f, 0) // coeff(f, 1)
            if lam < 0
                ok = false; push!(reasons, "negative rational eigenvalue $lam (C1 fails)")
                push!(factors, (f=string(f), mult=m, dim=-1, field=nothing, B=nothing, place_signs=Symbol[], lam_str=string(lam)))
                continue
            end
            B = kernel(S - lam * identity_matrix(QQ, n), side=:right)
            dim = ncols(B)
            dim == m || (ok = false; push!(reasons, "eigenspace dim $dim != mult $m at lam=$lam"))
            dg = congruence_diag(transpose(B) * d * B)
            if any(iszero, dg)
                ok = false; push!(reasons, "degenerate Gram on eigenspace lam=$lam")
                ps = [:degen]
            else
                s = [t > 0 ? 1 : -1 for t in dg]
                ps = [all(==(1), s) ? :pos : all(==(-1), s) ? :neg : :indef]
                ps[1] == :indef && (ok = false; push!(reasons, "indefinite eigenspace at lam=$lam (eigenspace definiteness fails)"))
            end
            push!(factors, (f=string(f), mult=m, dim=dim, field=:QQ, B=B, place_signs=ps, lam_str=string(lam)))
        else
            K, a = number_field(f, "a", cached=false)
            if !is_totally_real(K)
                ok = false; push!(reasons, "factor $(f): non-real eigenvalues (C1 fails)")
                push!(factors, (f=string(f), mult=m, dim=-1, field=nothing, B=nothing, place_signs=Symbol[], lam_str="root of $(f)"))
                continue
            end
            if !is_totally_positive(a)
                ok = false; push!(reasons, "factor $(f): negative eigenvalue at some real place (C1 fails)")
                push!(factors, (f=string(f), mult=m, dim=-1, field=nothing, B=nothing, place_signs=Symbol[], lam_str="root of $(f)"))
                continue
            end
            SK = change_base_ring(K, S)
            B = kernel(SK - a * identity_matrix(K, n), side=:right)
            dim = ncols(B)
            dim == m || (ok = false; push!(reasons, "eigenspace dim $dim != mult $m for factor $(f)"))
            dg = congruence_diag(transpose(B) * change_base_ring(K, d) * B)
            if any(iszero, dg)
                ok = false; push!(reasons, "degenerate Gram on eigenspace of factor $(f)")
                ps = [:degen]
            else
                ps = Symbol[]
                for P in real_places(K)
                    s = [is_positive(t, P) ? 1 : -1 for t in dg]
                    push!(ps, all(==(1), s) ? :pos : all(==(-1), s) ? :neg : :indef)
                end
                any(==(:indef), ps) && (ok = false; push!(reasons, "indefinite eigenspace of factor $(f) at some real place (eigenspace definiteness fails)"))
            end
            push!(factors, (f=string(f), mult=m, dim=dim, field=:NF, B=B, place_signs=ps, lam_str="root of $(f)"))
        end
    end
    return (ok=ok, reasons=reasons, factors=factors, charpoly=string(cp), minpoly_squarefree=sqfree)
end

# -----------------------------------------------------------------------------
# perp_root_check (P7 on one side): W = span of positive-norm eigenvectors
# = ⊕ eigenspaces positive-definite at ≥1 real place. Integral δ ⊥ (E at one
# real place) <=> the K-valued pairing vanishes <=> δ ⊥ E at ALL places
# (real embeddings are injective), so the orthogonality system is RATIONAL.
# Λ_perp = saturated integer kernel (Nemo ZZ-kernel is pure/saturated).
# Definite Λ_perp -> exact short-vector existence check for norms in root_norms.
# Indefinite/degenerate Λ_perp -> :unresolved (NEVER an accept; attention flag).
# -----------------------------------------------------------------------------
function perp_root_check(factors, dZ::ZZMatrix, root_norms::Vector{Int})
    n = nrows(dZ); dQ = matrix(QQ, dZ)
    rows = Vector{Vector{QQFieldElem}}()
    n_posfac = 0
    for fc in factors
        fc.B === nothing && continue
        (:pos in fc.place_signs) || continue
        n_posfac += 1
        if fc.field == :QQ
            C = transpose(fc.B) * dQ
            for r in 1:nrows(C)
                push!(rows, [C[r, j] for j in 1:n])
            end
        else
            K = base_ring(fc.B)
            C = transpose(fc.B) * change_base_ring(K, dQ)
            for r in 1:nrows(C), t in 0:(degree(K) - 1)
                push!(rows, [coeff(C[r, j], t) for j in 1:n])
            end
        end
    end
    n_posfac == 0 && return (status=:no_positive_space, n_posfac=0, perp_rank=n, roots_found=0)
    M = zero_matrix(QQ, length(rows), n)
    for (r, row) in enumerate(rows), j in 1:n
        M[r, j] = row[j]
    end
    den = reduce(lcm, [denominator(M[r, j]) for r in 1:nrows(M), j in 1:n][:]; init=ZZ(1))
    Mz = matrix(ZZ, [numerator(M[r, j] * den) for r in 1:nrows(M), j in 1:n])
    Kz = kernel(Mz, side=:right)                       # columns = saturated basis
    rk = ncols(Kz)
    rk == 0 && return (status=:clean, n_posfac=n_posfac, perp_rank=0, roots_found=0)
    G = transpose(Kz) * dZ * Kz
    dg = congruence_diag(matrix(QQ, G))
    if any(iszero, dg)
        return (status=:unresolved_degenerate_perp, n_posfac=n_posfac, perp_rank=rk, roots_found=-1)
    elseif all(t -> t > 0, dg)
        sgn = 1
    elseif all(t -> t < 0, dg)
        sgn = -1
    else
        return (status=:unresolved_indefinite_perp, n_posfac=n_posfac, perp_rank=rk, roots_found=-1)
    end
    (2 * sgn) in root_norms ||
        return (status=:clean, n_posfac=n_posfac, perp_rank=rk, roots_found=0)
    L = integer_lattice(gram = sgn == 1 ? G : -G)
    nroots = 2 * count(t -> t[2] == 2, short_vectors(L, 2))   # vectors, both signs
    return (status = nroots > 0 ? :root_found : :clean,
            n_posfac=n_posfac, perp_rank=rk, roots_found=nroots)
end

_facsum(sa) = [Dict("factor" => fc.f, "mult" => fc.mult, "dim" => fc.dim,
                    "eigenvalue" => fc.lam_str,
                    "place_signs" => [string(s) for s in fc.place_signs]) for fc in sa.factors]

# -----------------------------------------------------------------------------
# check_flux — THE predicate chain (P1..P7 in order, fail-fast).
# Flux datum: (N; m, m̃); callers working on the middle block only omit m, m̃.
# mode = "pair_block": run on a declared positive-definite block (T1⊗T2);
#   P6/P7 are then BLOCK-RESTRICTED (recorded as such in P6_scope — the
#   full-lattice quantifier needs the glue list of harness.jl). Corner charges
#   are REFUSED in this mode (corners are separate Künneth summands, not part
#   of any T1⊗T2 sub-block).
# mode = "control_full": reduced-lattice control problem (arXiv:2103.03250
#   §2.2-2.3 conventions); P6/P7 are the FULL exact admissibility /
#   stabilization test on that lattice.
# root_norms: nothing (default) derives from mode — control_full → [2,-2],
#   pair_block → [-2] (the K3 (−2)-root condition R); an explicit vector
#   overrides (every shipped call site passes it explicitly).
# -----------------------------------------------------------------------------
function check_flux(d1::ZZMatrix, d2::ZZMatrix, N::ZZMatrix;
                    m::Int=0, mtilde::Int=0,
                    N_cap::Int=24, root_norms::Union{Nothing,Vector{Int}}=nothing,
                    mode::String="pair_block", rank_req::Int=2)
    mode in ("pair_block", "control_full") ||
        error("check_flux: unknown mode '$mode' (allowed: pair_block, control_full)")
    rn = root_norms === nothing ? (mode == "control_full" ? [2, -2] : [-2]) : root_norms
    if mode == "pair_block" && (m != 0 || mtilde != 0)
        error("check_flux: corner charges (m,m̃)=($m,$mtilde) are not representable in a " *
              "pair_block sub-block (corner blocks are separate Künneth summands) " *
              "— use mode=control_full on full-lattice data")
    end
    res = Dict{String,Any}("predicates_version" => PREDICATE_SET_VERSION,
                           "mode" => mode, "N_cap" => N_cap, "root_norms" => rn,
                           "m" => m, "mtilde" => mtilde, "corner_mm" => m * mtilde)
    rej(p, why) = (res["verdict"] = "REJECT"; res["rejected_by"] = p; res["reason"] = why; res)
    # P1_integrality — (N, m, m̃) integral by type; evenness certificate of the
    # middle-block Q = tr(S) (corner adds 2m·m̃, even for all integers)
    Q = ZZ(tr(N * d2 * transpose(N) * d1))
    res["Q"] = string(Q)
    iseven(Q) || return rej("P1_integrality", "Q odd: even-lattice certificate violated")
    nf_mid = Int(divexact(Q, 2))
    nf = nf_mid + m * mtilde
    res["N_flux_mid"] = nf_mid
    res["N_flux"] = nf
    # P2_selfdual_membership — condition (V) + corner clause, exact.
    # (V)-existence ⟺ (C1) ∧ (C0) (self-dual existence lemma, see header).
    (nrows(N) == nrows(d1) && ncols(N) == nrows(d2)) ||
        return rej("P2_selfdual_membership", "flux shape does not match declared block")
    if !(m == 0 && mtilde == 0) && m * mtilde < 1
        return rej("P2_selfdual_membership",
            "corner block not self-dualizable: (m,m̃)=($m,$mtilde) mixed-sign or single-corner admits no Minkowski vacuum (corner clause R-c)")
    end
    S1z = N * d2 * transpose(N) * d1
    S2z = transpose(N) * d1 * N * d2
    S1 = matrix(QQ, S1z); S2 = matrix(QQ, S2z)
    ok1, why1 = c1_check(S1)
    ok1 || return rej("P2_selfdual_membership", "side1: $why1 — no self-dual vacuum (V_N ≠ ∅ ⟺ (C1)∧(C0))")
    ok2, why2 = c1_check(S2)
    ok2 || return rej("P2_selfdual_membership", "side2: $why2 — no self-dual vacuum (V_N ≠ ∅ ⟺ (C1)∧(C0))")
    rkN = rank(N); rk1 = rank(S1z); rk2 = rank(S2z)
    res["rank_N"] = rkN
    (rk1 == rkN && rk2 == rkN) || return rej("P2_selfdual_membership",
        "(C0) fails: rank(N)=$rkN, rank(S)=$rk1, rank(S̃)=$rk2 — ker g̃ ≠ ker S or ker g ≠ ker S̃; V_N = ∅")
    res["P2"] = Dict(
        "C1" => "pass (both sides: minpoly squarefree, totally-real, spectrum ≥ 0 — exact)",
        "C0" => "pass (rank N = rank S = rank S̃ = $rkN — ker g̃ = ker S, ker g = ker S̃)",
        "corner" => (m == 0 && mtilde == 0) ? "middle block only (m=m̃=0)" : "m·m̃=$(m*mtilde) ≥ 1 (volume-ratio channel open)",
        "basis" => "Künneth-block membership BY CONSTRUCTION of enumeration basis; list completeness is owned by the glue-list harness (coherent-tuple list specification)")
    # P3_primitivity — TAG ONLY (no primitivity imposed)
    res["content"] = Int(reduce(gcd, [N[i, j] for i in 1:nrows(N) for j in 1:ncols(N)]; init=ZZ(0)))
    # P4_orthogonality — STRATUM FILTER: block level only
    if mode == "pair_block" && rkN != rank_req
        res["P4_semantics"] = "stratum filter; EXCLUDED from theorem coverage accounting"
        res["counts_toward_theorem"] = false
        return rej("P4_orthogonality", "rank(N)=$rkN != rank_req=$rank_req — stratum filter, NOT a condition of the theorem")
    end
    # P5_tadpole — 0 <= N_flux = ½tr(S) + m·m̃ <= N_cap, exact (the ½ re-derived here)
    (0 <= nf <= N_cap) || return rej("P5_tadpole", "N_flux=$(nf) (=$(nf_mid) middle + $(m*mtilde) corner) outside [0, $N_cap]")
    # P6_stabilization — eigenspace definiteness + (C1), both sides, exact
    sa1 = side_analysis(S1, matrix(QQ, d1))
    sa2 = side_analysis(S2, matrix(QQ, d2))
    res["side1"] = Dict("charpoly" => sa1.charpoly, "factors" => _facsum(sa1), "reasons" => sa1.reasons)
    res["side2"] = Dict("charpoly" => sa2.charpoly, "factors" => _facsum(sa2), "reasons" => sa2.reasons)
    res["P6_scope"] = mode == "pair_block" ?
        "BLOCK_RESTRICTED (the full-lattice quantifier needs the glue list of harness.jl)" : "FULL_REDUCED_LATTICE"
    (sa1.ok && sa2.ok) || return rej("P6_stabilization",
        join(vcat(sa1.reasons, sa2.reasons), "; "))
    # P7_no_root — both sides
    p7a = perp_root_check(sa1.factors, d1, rn)
    p7b = perp_root_check(sa2.factors, d2, rn)
    res["P7"] = Dict("side1" => string(p7a.status), "side2" => string(p7b.status),
                     "perp_ranks" => [p7a.perp_rank, p7b.perp_rank],
                     "roots_found" => [p7a.roots_found, p7b.roots_found],
                     "scope" => res["P6_scope"])
    for (tag, p7) in (("side1", p7a), ("side2", p7b))
        p7.status == :root_found &&
            return rej("P7_no_root", "$tag: $(p7.roots_found) root(s) orthogonal to positive-norm eigenspaces")
        startswith(string(p7.status), "unresolved") &&
            return rej("P7_no_root", "$tag: $(p7.status) — Λ_perp not definite; UNRESOLVED is never an accept (attention)")
    end
    res["verdict"] = "ACCEPT"
    res["rejected_by"] = ""
    return res
end

end # module SweepPredicates

# subtract.jl — C7: Möbius counterterm assembly (Phase B1).
# ===========================================================================
# Transliterates STSubtractionFormula (reference/SubTropica.wl:10722-10841;
# alias STTropicalSubtraction wl:10844), paper eqs 3.25-3.36. Reference
# specification: SPEC_B1 §4 (steps 0-9) + §6 (refusal contract).
#
# Include contract: include()d inside `module SubTropica` AFTER src/types.jl and
# src/b1_types.jl. NO `using` statements here (types.jl header law); Nemo/JSON
# are in scope from the enclosing module (standalone tests do `using Nemo`
# first and include the three files by absolute path — DESIGN_B1.md §3).
#
# Input:
#   E         :: EulerIntegrand — nu is read in the DLOG chart, i.e.
#                I = ∏ x_i^{ν_i} · ∏ P_j^{e_j} against the measure ∏ dx_i/x_i
#                (SPEC_B1 §0: "the .wl works with ∫∏dx_i/x_i and converts to
#                flat measure only at output", wl:10740-10745, 10837).
#                E.prefactor passes through UNTOUCHED — C8 applies it once
#                globally (CONTRACTS.md §b assembly); it is NOT folded into
#                the emitted CounterTerms.
#   td        :: TropicalData (C4 output; primitive OUTER normals).
#   sigma_div :: Vector{SigmaFace} (C5 output; ordering law SPEC_B1 §3.2:
#                sigma_div[1] == ∅, then size-lex; face indices GLOBAL asc.).
#   ws        :: AbstractDict facet => primitive sign-flipped w-vector with
#                w·ρ_f == −1 (C5 DivData.w; `missing`/absent = the CORRECTED
#                "NotFound" sentinel, SPEC_B1 D1 ⇒ GeometricPropertyViolated).
#
# Output: Vector{CounterTerm}, grouped by face in sigma_div order, terms per
# face in sigma_div-excess order (empty excess = identity term FIRST).
# Regulators: single ε only (the .wl takes a list, default {eps}, wl:10723;
# B1 keeps ε-linearity via EpsExp — SPEC_B1 §0).
#
# THE TWO WARNING PINS (CONTRACTS.md §c / DESIGN_B1.md §1) ARE LAW HERE and
# are quoted verbatim at the exact lines where they apply — search for
# "WARNING-1" (the net sign) and "WARNING-2" (the I-vs-J subset).
# ===========================================================================

# ---------------------------------------------------------------------------
# Small exact helpers (no external deps beyond Nemo already in scope).
# ---------------------------------------------------------------------------

# dot over the first n (integration-var) components of an exponent vector.
# Ring exponent vectors may extend over kinvars (types.jl variable order:
# integration vars first); rays live in Z^n only, so the projection is just
# "first n components". Kinvar components never enter tropical bookkeeping
# (SPEC_B1 §2.1: the polytope is built over the integration vars only).
_c7_dotn(ev::AbstractVector{<:Integer}, rho::AbstractVector{<:Integer}, n::Int) =
    sum(BigInt(ev[i]) * BigInt(rho[i]) for i in 1:n; init = big(0))

_c7_dot(a::AbstractVector{<:Integer}, b::AbstractVector{<:Integer}) =
    sum(BigInt(a[i]) * BigInt(b[i]) for i in eachindex(a); init = big(0))

"""_c7_subsets_lex(n, m) — all size-m subsets of 1:n in LEXICOGRAPHIC order,
identical to Mathematica `Subsets[Range[n],{m}]` (SPEC_B1 §4 step 5 pin;
gate G2 asserts the order against a hand-listed Subsets[Range[4],{2}])."""
function _c7_subsets_lex(n::Int, m::Int)
    m == 0 && return [Int[]]
    m > n  && return Vector{Int}[]
    out = Vector{Int}[]
    idx = collect(1:m)
    while true
        push!(out, copy(idx))
        j = m
        while j >= 1 && idx[j] == n - m + j
            j -= 1
        end
        j == 0 && break
        idx[j] += 1
        for k in j+1:m
            idx[k] = idx[k-1] + 1
        end
    end
    return out
end

"""Exact determinant of a square Rational{BigInt} matrix given as rows
(fraction-free enough at B1 sizes; nothing float — SPEC_B1 §3.2 exactness rule)."""
function _c7_det(rows::Vector{Vector{Rational{BigInt}}})
    m = length(rows)
    m == 0 && return Rational{BigInt}(1)
    A = [copy(r) for r in rows]
    @assert all(length(r) == m for r in A) "_c7_det: matrix not square"
    det = Rational{BigInt}(1)
    for c in 1:m
        p = findfirst(r -> !iszero(A[r][c]), c:m)
        p === nothing && return Rational{BigInt}(0)
        p = p + c - 1
        if p != c
            A[c], A[p] = A[p], A[c]
            det = -det
        end
        det *= A[c][c]
        for r in c+1:m
            f = A[r][c] // A[c][c]
            iszero(f) && continue
            for cc in c:m
                A[r][cc] -= f * A[c][cc]
            end
        end
    end
    return det
end

# ---------------------------------------------------------------------------
# TropI evaluation (SPEC_B1 §3.1; STTropicalizeIntegrand wl:10479-10481 +
# STEvalRay wl:10354-10356; paper eq 3.12/3.14).
#: the former private
# duplicates now DELEGATE to tropical.jl's exported C3/C5 semantics
# (src/tropical.jl trop_poly/trop_on_ray/restrict_poly — DESIGN_B1 §2
# ownership). Standalone runs of test_b1_subtract.jl include tropical.jl
# before this file; the parity assertion lives in test_b1_parity.jl.
# ---------------------------------------------------------------------------
_c7_trop_poly(P::QQMPolyRingElem, rho::Vector{BigInt}, n::Int) =
    trop_poly(P, rho, n)

"""TropI(ρ) = Σ_i ν_i ρ_i + Σ_j e_j · trop(P_j)(ρ), exact EpsExp
(wl:10731 `trValues`) — delegates to trop_on_ray (tropical.jl)."""
_c7_trop_on_ray(E::EulerIntegrand, rho::Vector{BigInt}, n::Int) =
    trop_on_ray(E, rho)

# ---------------------------------------------------------------------------
# Restriction operator = initial form (STrestrictPoly wl:10506-10511; paper
# eq 3.30, "initial form" paper:1168-1170): keep the monomials of P exposed
# by outer normal ρ, i.e. the terms whose support maximizes m·ρ (the .wl
# states it as the MIN of −m·ρ, wl:10508-10509 — identical; depth = 0 always
# in C7). ORIGINAL coefficients kept (`(cr[[1]]/.coefsRules)`, wl:10510) —
# only monomial SELECTION uses the (coeff-free) support.
# Face restriction = SEQUENTIAL single-ray restrictions; the .wl applies
# face[[2;;]] first, then face[[1]] (wl:10534) — order-independent (iterated
# initial forms along normals of a common face commute); B1 goes
# left-to-right and the test file asserts order-independence.
# DELEGATES to tropical.jl's exported
# restrict_poly (tropical.jl, the DESIGN_B1 §2 owner); parity assertion in
# test_b1_parity.jl.
# ---------------------------------------------------------------------------
_c7_restrict_poly(P::QQMPolyRingElem, rho::Vector{BigInt}, n::Int) =
    restrict_poly(P, rho, n)

# ---------------------------------------------------------------------------
# Step 5 — I/J subsets (wl:10777-10782; paper eq 3.34).
#
# WARNING-2 (CONTRACTS.md §c, verbatim): "I vs J: the completion/division/
# gauge subset is the COMPLEMENT. Correct: I = the FIRST index subset with
# det(rays[fc, I]) ≠ 0 (the det≠0 ray-support subset); J = Complement(1:n, I).
# The Kronecker-delta completion rows in Vol = |det(rays|e_J)|/∏(−TropI) are
# indexed by J; the output division /∏x_J uses J; the gauge-fix sets vars ∉ J
# to 1 — all three use J, the complement." The wrong reading (J = I) yields
# finite, numerically-plausible, wrong-volume wrong-gauge results — fixture
# F2 (face {2}: I = [2], J = [1]) is the mandatory discriminator (gate G4).
# ---------------------------------------------------------------------------

"""_c7_first_I(face, rays, n) — the FIRST (lexicographic, Mma Subsets order)
size-|face| column subset I with det(rays[face][:,I]) ≠ 0 (wl:10777-10780).
`rays[face][:,I]` = the |face|×|face| minor with ROWS the face's rays and
COLUMNS I. Refuses :NoBasisCompletion if none exists (cannot happen for
linearly independent face normals, paper:1123-1125 — refuse loudly, never
FirstCase-Missing, SPEC_B1 §4 step 5)."""
function _c7_first_I(face::Vector{Int}, rays::Vector{Vector{BigInt}}, n::Int)
    isempty(face) && return Int[]
    m = length(face)
    for I in _c7_subsets_lex(n, m)
        minor = [Rational{BigInt}[rays[f][c] for c in I] for f in face]
        iszero(_c7_det(minor)) || return I
    end
    throw(B1Refusal(:NoBasisCompletion,
        Dict{String,Any}("face" => face,
                         "rays" => [rays[f] for f in face],
                         "hint" => "face normals linearly dependent — should be unreachable (paper:1123-1125)")))
end

"""J = Complement(1:n, I) — WARNING-2 above. ∅ face ⇒ J = 1:n (wl:10781)."""
function _c7_J_of(face::Vector{Int}, rays::Vector{Vector{BigInt}}, n::Int)
    isempty(face) && return collect(1:n)              # wl:10781: ∅ ↦ Range[n]
    return setdiff(collect(1:n), _c7_first_I(face, rays, n))   # wl:10780:
end                                                   # Complement[Range[n],I]

# ---------------------------------------------------------------------------
# Step 7 — volume determinant (wl:10790-10797; paper eq 3.35).
# Matrix from ROWS: the face's rays, then the Kronecker completion rows e_j
# for j ∈ J — J per WARNING-2 (the paper's M_σ = σ|η is the transpose —
# discrepancy D4, same |det|). |face| + |J| = n since |I| = |face|.
# Vol(face) = vol_det / ∏(−TropI(ρ)), ρ ∈ face — the (−TropI) factors are
# kept UNEXPANDED in CounterTerm.vol_trops: on log rays −TropI = −b·ε, so
# Vol ∝ 1/ε^{|face|} is the ONLY explicit ε-pole source of B1
# (paper:1219-1221; expansion is C8's job, truncated-Laurent guard there).
# ---------------------------------------------------------------------------
function _c7_vol_det(face::Vector{Int}, J::Vector{Int},
                     rays::Vector{Vector{BigInt}}, n::Int)
    isempty(face) && return big(1)                    # Vol(∅) = 1 (wl:10796)
    rows = Vector{Vector{Rational{BigInt}}}()
    for f in face                                     # ray rows
        push!(rows, Rational{BigInt}[rays[f][i] for i in 1:n])
    end
    for j in J                                        # completion rows e_j,
        row = zeros(Rational{BigInt}, n)              # j ∈ J (WARNING-2;
        row[j] = 1                                    #  wl:10793 KroneckerDelta)
        push!(rows, row)
    end
    @assert length(rows) == n "_c7_vol_det: |face|+|J| != n (I/J bookkeeping bug)"
    d = _c7_det(rows)
    @assert denominator(d) == 1 "_c7_vol_det: non-integer det of integer matrix"
    return abs(numerator(d))
end

# ---------------------------------------------------------------------------
# Step 6 — extra monomial x^{u_σ} (wl:10785-10788; paper eq 3.36):
# u_σ = Σ_{ρ∈face} TropI(ρ)·π_J(w_ρ). The .wl builds the NEGATIVE exponent
# −TropI·w and then takes 1/x^(that), and implements the projector π_J as the
# substitution vars∉J → 1 (discrepancy D5 — same mechanism reused for the
# step-9 gauge-fix). Components are ε-DEPENDENT (TropI = b·ε on log rays):
# the twist is an x^{O(ε)} factor, e.g. x2^{−ε} in fixture F3 (SPEC_B1 §8.3).
# ---------------------------------------------------------------------------
function _c7_u_sigma(face::Vector{Int}, trvals::Vector{EpsExp},
                     wdict::Dict{Int,Vector{BigInt}}, Jset::Set{Int}, n::Int)
    u = [zero(EpsExp) for _ in 1:n]
    for rho in face, i in 1:n
        (i in Jset) && (u[i] += trvals[rho] * wdict[rho][i])
    end
    return u
end

# ---------------------------------------------------------------------------
# Step 8 helper — jacobian factor bank (wl:10799-10804; paper eq 3.29).
# (1−u_j) = 1/(1+x^{w_j}) — u = 1 − 1/(1+x^w) is the CODE convention
# (wl:10636/10669), NOT the paper's v_ρ = 1/(1+x^w) of eq 3.28 (v = 1−u,
# discrepancy D3). The counterterm carries (1−u_j)^{1−TropI(ρ_j)}
# = (1+x^{w_j})^{−(1−TropI(ρ_j))} with the exponent ε-DEPENDENT (TropI with
# regulators kept — RT-extra-jac; CONTRACTS §b; upstream's own puzzled
# comment at wl:10824 is stale, discrepancy D6 — do not "fix").
# For general integer w = w⁺ − w⁻ (w± ≥ 0, componentwise):
#   (1+x^w)^{−(1−T)} = (x^{w⁻} + x^{w⁺})^{−(1−T)} · x^{w⁻·(1−T)}
# so the POLY entry is x^{w⁻}+x^{w⁺} (a genuine ring element) and the
# monomial compensation w⁻·(1−T) folds into the ambient exponent vector
# (mono is EpsExp-valued, so the ε-dependent shift is representable).
# For the all-nonnegative w of the fixtures this reduces to poly = 1+x^w,
# compensation 0.
# ---------------------------------------------------------------------------
function _c7_jac_factor(ring::QQMPolyRing, n::Int, w::Vector{BigInt})
    wplus  = [max(w[i], 0) for i in 1:n]
    wminus = [max(-w[i], 0) for i in 1:n]
    pad = zeros(Int, nvars(ring) - n)                 # kinvar slots (always 0)
    B = MPolyBuildCtx(ring)
    push_term!(B, one(base_ring(ring)), vcat(Int.(wminus), pad))
    push_term!(B, one(base_ring(ring)), vcat(Int.(wplus), pad))
    return finish(B), wminus
end

# ---------------------------------------------------------------------------
# Gauge-fix: vars ∉ J → 1 (wl:10834; the same substitution that implements
# π_J in extraMons, D5). Applied to EVERY poly of every term — restricted
# P_j's AND jacobian factors — after the Möbius assembly, exactly as the .wl
# applies it to the whole `result`. Kinvars (ring slots > n) untouched.
# ---------------------------------------------------------------------------
function _c7_gauge_fix(P::QQMPolyRingElem, Jset::Set{Int}, n::Int)
    length(Jset) == n && return P                     # nothing to fix
    R = parent(P)
    subs = [(i <= n && !(i in Jset)) ? one(R) : gen(R, i) for i in 1:nvars(R)]
    return evaluate(P, subs)
end

# ---------------------------------------------------------------------------
# THE ENTRY POINT — SPEC_B1 §4, STSubtractionFormula wl:10722-10841.
# ---------------------------------------------------------------------------
"""
    subtraction_terms(E, td, sigma_div, ws) -> Vector{CounterTerm}

Möbius subtraction assembly (C7). Terms grouped by face in `sigma_div` order,
per face in sigma_div-excess order (identity term first). Refuses (typed,
SPEC_B1 §6, BEFORE any CounterTerm is emitted): :TropIllDefined,
PowerDivergentRefusal, GeometricPropertyViolated, WNotUnitNormalized,
:NoBasisCompletion. `ws` maps each divergent facet (global ray index) to its
primitive sign-flipped w with w·ρ = −1; `missing`/absent entries are the
D1-corrected NotFound sentinel.
"""
function subtraction_terms(E::EulerIntegrand, td::TropicalData,
                           sigma_div::Vector{SigmaFace},
                           ws::AbstractDict)::Vector{CounterTerm}
    n = td.ambient_dim
    @assert n == length(E.vars) "subtraction_terms: ambient_dim != #vars"
    @assert nvars(E.ring) >= n "subtraction_terms: ring smaller than vars"
    R = E.ring

    # ---- Step 1 — trop values on every ray (wl:10729-10731), refusal
    # ladder order: ill-defined (§6.1) BEFORE power (§6.3) — full passes.
    trvals = EpsExp[_c7_trop_on_ray(E, rho, n) for rho in td.rays]
    for (i, t) in enumerate(trvals)
        # .wl abort guard MemberQ[trValues,0] → notwelldef (wl:10948-10956/
        # 11071-11079); typical cause: unfixed Cheng-Wu gauge.
        is_illdefined(t) && throw(B1Refusal(:TropIllDefined,
            Dict{String,Any}("ray" => i, "ray_vec" => td.rays[i],
                             "trvals" => trvals)))
    end
    for (i, t) in enumerate(trvals)
        # .wl only PRINTS "Power divergence on ..." (wl:10732-10734) and
        # proceeds to garbage (the formula assumes at-most-log,
        # paper:1162-1163) — B1 refuses instead (D10, DESIGN RT-ix; the B2
        # continuation route, continue.jl, consumes this refusal).
        is_power(t) && throw(PowerDivergentRefusal(i, t))
    end
    div_facets = [i for (i, t) in enumerate(trvals) if is_divergent(t)]

    # ---- Step 0 — locally finite early return (wl:10739-10746): the whole
    # integrand as the single ∅-face identity term, measure converted
    # dlog → flat by /∏vars, i.e. mono = ν .- 1 ("Divide by Times@@vars to
    # convert from dx/x to flat dx convention", .wl comment wl:10740-10742).
    if isempty(div_facets)
        @assert all(f -> isempty(f.facets), sigma_div) "step 0: no divergent facet, but sigma_div has nonempty faces"
        return [CounterTerm(Int[], Int[], warning1_sign(Int[]), big(1),
                            EpsExp[], collect(1:n), copy(E.vars),
                            EpsExp[E.nu[i] - 1 for i in 1:n],
                            [zero(EpsExp) for _ in 1:n], copy(E.polys))]
    end

    # ---- Step 2 — faces (wl:10748-10749): from C5; ordering-law and
    # consistency asserts (divFaces = Cases[faces, SubsetQ[divFacets,v]] is a
    # no-op upstream — SPEC_B1 §3.2 pin — ported as assertions, not filters).
    @assert !isempty(sigma_div) && isempty(sigma_div[1].facets) "sigma_div ordering law: sigma_div[1] must be the empty face (wl:10272)"
    singles = sort!([f.facets[1] for f in sigma_div if length(f.facets) == 1])
    @assert singles == div_facets "sigma_div singletons != divergent facets (C5/C7 disagree — caller bug)"
    @assert all(issorted(f.facets) for f in sigma_div) "sigma_div faces must ascend in GLOBAL ray index"

    # ---- Step 3 — w's (wl:10752-10765). Missing w = the D1-CORRECTED
    # NotFound sentinel (upstream's guard `Count[us,"NotFound!"]>0` at
    # wl:10760 is DEAD CODE — STProduceUs emits "NotFound", no bang,
    # wl:10638; here missing-ness is typed and this check IS the guard).
    # w·ρ_f == −1 is the B1 pin (D2: upstream runs normalizeQ=False at
    # wl:10758 and can silently use w·ρ_f = −c, c > 1 — B1 refuses).
    pairs = [f.facets for f in sigma_div if length(f.facets) == 2]
    compat = f -> sort!(unique!([g for p in pairs for g in p if f in p && g != f]))
    wdict = Dict{Int,Vector{BigInt}}()
    for f in div_facets
        wraw = get(ws, f, missing)
        wraw === missing && throw(GeometricPropertyViolated(f, compat(f)))
        w = BigInt.(collect(wraw))
        @assert length(w) == n "w for facet $f has wrong length"
        prd = _c7_dot(w, td.rays[f])
        prd == -1 || throw(WNotUnitNormalized(f, w, td.rays[f], prd))
        for g in compat(f)   # GP: w_ρ·ρ' = 0 for compatible ρ' (paper:1118-1121)
            @assert _c7_dot(w, td.rays[g]) == 0 "w for facet $f not orthogonal to compatible facet $g (C5 handed inconsistent data)"
        end
        wdict[f] = w
    end

    # Restriction memo (the .wl memoizes STrestrictIntegrand, wl:10530;
    # restrictions repeat heavily across faces/subfaces). Keyed by
    # (global ray index, polynomial) — per-call, no shared mutable state.
    memo = Dict{Tuple{Int,QQMPolyRingElem},QQMPolyRingElem}()
    restr(P::QQMPolyRingElem, ri::Int) =
        get!(() -> _c7_restrict_poly(P, td.rays[ri], n), memo, (ri, P))

    # ---- Steps 5-9 per face — regularizeFace (wl:10806-10838; eqs 3.30-3.32).
    out = CounterTerm[]
    for sf in sigma_div
        face = sf.facets
        # Step 5 — J (WARNING-2: J = Complement of the first det≠0 subset I;
        # see _c7_first_I/_c7_J_of block comment for the verbatim pin).
        J    = _c7_J_of(face, td.rays, n)
        Jset = Set(J)
        # Step 7 — volume pieces (wl:10790-10797); poles stay UNEXPANDED.
        vd   = _c7_vol_det(face, J, td.rays, n)
        vtr  = EpsExp[-trvals[rho] for rho in face]
        # Step 6 — extra monomial (wl:10785-10788; eq 3.36).
        usig = _c7_u_sigma(face, trvals, wdict, Jset, n)

        # newIntegrand (wl:10826-10830): ∅ face ⇒ the RAW integrand (no
        # restriction, no twist); else restrict along the face's rays
        # (sequentially — see restriction block comment) and multiply by
        # x^{u_σ} (folded into the ambient exponent vector, not the polys).
        base = E.polys
        for ri in face
            base = Tuple{QQMPolyRingElem,EpsExp}[(restr(P, ri), e) for (P, e) in base]
        end
        amb0 = EpsExp[E.nu[i] + usig[i] for i in 1:n]

        # Step 9 — the Möbius sum: all supersets σ' ⊇ face in sigma_div
        # order; each contributes its EXCESS σ'∖face (wl:10807-10812:
        # compFaces = Complement[#,face]& /@ Cases[divFaces, SubsetQ[#,face]]).
        excesses = Vector{Int}[setdiff(g.facets, face)
                               for g in sigma_div if issubset(face, g.facets)]
        @assert !isempty(excesses) && isempty(excesses[1]) "superset scan lost the identity term (σ'=face must be first — size-lex order)"

        for s in excesses
            polys = base
            amb   = copy(amb0)
            if !isempty(s)
                # counterterm: restrict newIntegrand along the excess rays
                # (wl:10827: STrestrictIntegrand[newIntegrand, vars,
                # rays[[subface]]]) ...
                for ri in s
                    polys = Tuple{QQMPolyRingElem,EpsExp}[(restr(P, ri), e) for (P, e) in polys]
                end
                # ... and attach the jacobian factors (1−u_j)^{1−TropI(ρ_j)}
                # for j ∈ subface (wl:10824: (jacFace[subface])^(1-trValues[[
                # subface]]) — jacobians are NOT restricted). Exponent
                # −(1−TropI(ρ_j)) is ε-DEPENDENT (RT-extra-jac; gate G5).
                for j in s
                    jp, wminus = _c7_jac_factor(R, n, wdict[j])
                    push!(polys, (jp, -(one(EpsExp) - trvals[j])))
                    if any(!iszero, wminus)   # x^{w⁻·(1−T)} compensation
                        for i in 1:n
                            amb[i] += (one(EpsExp) - trvals[j]) * wminus[i]
                        end
                    end
                end
            end

            # Gauge-fix vars ∉ J → 1 (wl:10834), drop unit factors; then the
            # /∏x_J flat-dx division (wl:10837) = −1 on each surviving slot.
            gp = Tuple{QQMPolyRingElem,EpsExp}[]
            for (P, e) in polys
                Pg = _c7_gauge_fix(P, Jset, n)
                isone(Pg) || push!(gp, (Pg, e))
            end
            mono = EpsExp[amb[i] - 1 for i in J]

            # WARNING-1 (CONTRACTS.md §c, verbatim): "counter-term NET sign
            # is (−1)^|subface|. The code builds Times[-a[subface], ...]
            # where a[s] = (−1)^(|s|+1) [wl:10823 assembly, wl:10833
            # a-substitution]; the NET assembled sign is therefore
            # −(−1)^(|s|+1) = (−1)^|s|: empty subface → +1 (the identity
            # term), single facet → −1 (subtracted, as it must be)." The
            # wrong v1 reading ((−1)^(|s|+1)) ADDS single-facet counterterms;
            # gate G3(c)/(d) catches it. Sign MUST come from warning1_sign
            # (b1_types.jl) — stored explicitly so the G3(d) flipped-sign
            # negative control can build a wrong term at the eval layer.
            push!(out, CounterTerm(copy(face), s, warning1_sign(s), vd, vtr,
                                   copy(J), E.vars[J], mono, copy(usig), gp))
        end
    end
    return out
end

"""Convenience overload: accept C5's per-facet records directly."""
subtraction_terms(E::EulerIntegrand, td::TropicalData,
                  sigma_div::Vector{SigmaFace}, dfs::Vector{DivergentFacet}) =
    subtraction_terms(E, td, sigma_div,
        Dict{Int,Union{Vector{BigInt},Missing}}(d.ray => d.w for d in dfs))

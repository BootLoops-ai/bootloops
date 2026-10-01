# Fuchsia.jl — Lee's algorithm (arXiv:1411.0911) for reducing a Fuchsian linear system
# dI/dx = A(eps,x) I  to epsilon-form  A = eps * Atilde(x),  Atilde eps-independent, dlog.
#
# The three Lee stages (this file implements the ones a UT construction needs):
#   (1) Fuchsianize: remove poles of order > 1 in x by balance transformations.  IBP-derived
#       connection matrices in the ratio variable, after stripping the overall LS prefactor, are
#       Fuchsian to begin with for the dlog sectors — we ASSERT and CHECK this rather than running
#       the full pole-reduction, and fall back to the documented STOP if not Fuchsian.
#   (2) Normalize: balance the residue eigenvalues at all singular points (incl. infinity) so that
#       every residue eigenvalue becomes  n*eps  (n ∈ Z).  A balance between points x1,x2 with a
#       rank-1 projector P built from a right-eigenvector at x1 and a left-eigenvector at x2 shifts
#       the corresponding eigenvalues by ∓1.  We choose balances to drive all eigenvalues to the
#       form  c0 + c1*eps  with c0 = 0  (UT normalization).
#   (3) Factor eps: with every residue = eps*const, find a constant gauge T (x-independent) that
#       makes the whole matrix proportional to eps.  In the generic semisimple case the matrix is
#       already eps*Atilde once stage (2) lands all eigenvalues on multiples of eps and the
#       eps^0 part of every residue vanishes; the residual constant rotation removes any leftover
#       eps^0 off-diagonal block.
#
# Gauge action on the connection (left-multiplication convention dI/dx = A I, I -> T I):
#       A  ->  T A T^{-1} + (dT/dx) T^{-1}.
#
# A *balance* (Lee eq.) with projector P and poles x1 (finite) , x2 (finite):
#       B(P,x1,x2) = (I - P) + ((x - x2)/(x - x1)) * P ,
#       B^{-1}     = (I - P) + ((x - x1)/(x - x2)) * P .
# With x1 = ∞ : B = (I - P) + (x - x2) P ; with x2 = ∞ : B = (I - P) + 1/(x - x1) P.
#
# Eigenvalues here are elements of Qe = Q(eps), i.e. polynomials/ratios c0 + c1*eps + ... .
# For Fuchsian dlog systems they are AFFINE in eps: a + b*eps with a,b ∈ Q.  We use that.

module Fuchsia

using Nemo
using ..RatFunc

export reduce_to_epsform, EpsFormResult, apply_gauge, gauge_balance,
       eigen_affine, eps0_part, is_epsform, certify_epsform,
       moser_reduce, MoserResult, dlog_alphabet

# ----------------------------------------------------------------------------------------------
# Gauge transformations
# ----------------------------------------------------------------------------------------------

"Apply gauge T (Qex matrix) to connection A: A -> T A T^{-1} + (dT/dx) T^{-1}."
function apply_gauge(ctx::Ctx, A, T)
    Tinv = inv(T)
    dT = dmat_dx(ctx, T)
    return T * A * Tinv + dT * Tinv
end

"""
    gauge_balance(ctx, P, x1, x2) -> T (Qex matrix)

Balance transformation matrix.  `P` is a Qe (constant) projector matrix; x1,x2 are pole
locations, each either a Qe value or the symbol :inf.
"""
function gauge_balance(ctx::Ctx, P, x1, x2)
    n = size(P,1)
    x = ctx.x
    Pq = _qe_to_qex(ctx, P)
    Iq = identity_matrix(ctx.Qex, n)
    if x1 === :inf && x2 === :inf
        error("balance with both poles at infinity is undefined")
    elseif x1 === :inf
        # B = (I-P) + (x - x2) P
        f = (x - ctx.Qex(x2))
        return (Iq - Pq) + f * Pq
    elseif x2 === :inf
        # B = (I-P) + (1/(x-x1)) P
        f = ctx.Qex(1) // (x - ctx.Qex(x1))
        return (Iq - Pq) + f * Pq
    else
        f = (x - ctx.Qex(x2)) // (x - ctx.Qex(x1))
        return (Iq - Pq) + f * Pq
    end
end

function _qe_to_qex(ctx::Ctx, M)
    m, n = size(M)
    out = zero_matrix(ctx.Qex, m, n)
    for i in 1:m, j in 1:n
        out[i,j] = ctx.Qex(M[i,j])
    end
    return out
end

# ----------------------------------------------------------------------------------------------
# Eigen-analysis over Qe = Q(eps), eigenvalues assumed affine a + b*eps
# ----------------------------------------------------------------------------------------------

"Extract the (a,b) of an affine Qe element a + b*eps; error if higher order in eps present."
function _affine_coeffs(ctx::Ctx, v)
    Qe = ctx.Qe
    num = numerator(v); den = denominator(v)
    # require den constant (degree 0 in eps) for an affine element
    if !is_constant(den)
        error("non-polynomial-in-eps eigenvalue: $v")
    end
    dc = constant_coefficient_q(den)
    a = coeff(num, 0) // dc
    b = degree(num) >= 1 ? coeff(num, 1) // dc : QQ(0)
    if degree(num) >= 2
        error("eigenvalue not affine in eps: $v")
    end
    return (a, b)   # rationals
end

constant_coefficient_q(p) = coeff(p, 0)

"""
    eps0_part(ctx, M) -> Qe-matrix of the eps^0 coefficient (eps -> 0) of a Qe matrix.

Computes M|_{eps=0} when finite; used to test whether a residue is eps*const.
"""
function eps0_part(ctx::Ctx, M)
    m, n = size(M)
    Qe = ctx.Qe
    out = zero_matrix(QQ, m, n)
    for i in 1:m, j in 1:n
        out[i,j] = _eval_qe_at0(M[i,j])
    end
    return out
end

# evaluate a Qe element at eps = 0 (must be finite)
function _eval_qe_at0(v)
    num = numerator(v); den = denominator(v)
    d0 = coeff(den, 0)
    if iszero(d0)
        error("eps=0 is a pole of $v")
    end
    return coeff(num, 0) // d0
end

# coefficient of eps^1 in a Qe element with constant denominator
function _eps1_coeff(v)
    num = numerator(v); den = denominator(v)
    if !is_constant(den); error("denominator depends on eps: $v"); end
    d0 = coeff(den, 0)
    return (degree(num) >= 1 ? coeff(num, 1) : QQ(0)) // d0
end

# ----------------------------------------------------------------------------------------------
# Epsilon-form check & certification
# ----------------------------------------------------------------------------------------------

"""
    is_epsform(ctx, A; verbose=false) -> (ok::Bool, Atilde_or_nothing)

A is in eps-form iff A = eps * B with B a Qex matrix that does NOT depend on eps.
Returns the candidate Atilde = A/eps and checks eps-independence.
"""
function is_epsform(ctx::Ctx, A; verbose=false)
    eps = ctx.eps
    m, n = size(A)
    B = zero_matrix(ctx.Qex, m, n)
    for i in 1:m, j in 1:n
        a = A[i,j]
        iszero(a) && continue
        # divide by eps: a/eps must be eps-regular
        b = a // ctx.Qex(eps)
        B[i,j] = b
    end
    # check B has no eps dependence: every entry's num & den polynomials over Qe must be
    # eps-free, i.e. lie in Q[x] up to overall.  Test: derivative wrt eps is zero.
    ok = true
    bad = (0,0)
    for i in 1:m, j in 1:n
        if !_qex_eps_free(ctx, B[i,j])
            ok = false; bad = (i,j); break
        end
    end
    if verbose && !ok
        @info "is_epsform: entry $bad still eps-dependent after /eps"
    end
    return ok, (ok ? B : nothing)
end

"Is a Qex element free of eps in both numerator and denominator coefficients?"
function _qex_eps_free(ctx::Ctx, f)
    iszero(f) && return true
    num = numerator(f); den = denominator(f)
    return _poly_coeffs_eps_free(num) && _poly_coeffs_eps_free(den)
end

# poly over Qe (=Q(eps)): all coefficients must be eps-free rationals (Qe elements with
# numerator degree 0 and constant denominator, i.e. honest rationals).
function _poly_coeffs_eps_free(p)
    for i in 0:degree(p)
        c = coeff(p, i)
        iszero(c) && continue
        cn = numerator(c); cd = denominator(c)
        if degree(cn) >= 1 || degree(cd) >= 1
            return false
        end
    end
    return true
end

"""
    certify_epsform(ctx, A) -> Dict with diagnostics

Runs: Fuchsian check, residue eigenvalue report, eps-form check, dlog (Atilde eps-free) check.
"""
function certify_epsform(ctx::Ctx, A)
    res = Dict{String,Any}()
    res["fuchsian"] = is_fuchsian(ctx, A)
    pts, alg = singular_points_x(ctx, A)
    res["singular_points"] = string.(pts)
    res["algebraic_factors"] = string.(alg)
    ok, B = is_epsform(ctx, A; verbose=false)
    res["epsform"] = ok
    if ok
        # residues of Atilde are constant rational matrices -> dlog
        dlog = true
        for a in pts
            Rb = residue_matrix(ctx, B, a)
            if !all(_qe_is_rational(Rb[i,j]) for i in 1:size(Rb,1), j in 1:size(Rb,2))
                dlog = false
            end
        end
        res["dlog"] = dlog
    else
        res["dlog"] = false
    end
    return res
end

function _qe_is_rational(c)
    iszero(c) && return true
    cn = numerator(c); cd = denominator(c)
    return degree(cn) == 0 && degree(cd) == 0
end

# ==============================================================================================
# Stage (1): Barkatou–Moser pole-order reduction (fuchsianize).
#
# For a connection dI/dx = A I with a pole of order r+1 ≥ 2 at x = a, write
#     A = A₀/(x-a)^{r+1} + A₁/(x-a)^r + ...
# Moser's theorem: the singularity is *regular* (reducible to a simple pole by a rational gauge)
# iff every such leading matrix A₀ is nilpotent — which is always the case for Feynman-integral
# DE systems (all physical singularities are regular).  When A₀ has a nonzero eigenvalue the
# singularity is genuinely irregular; for IBP-derived systems this signals a missing overall
# leading-singularity prefactor, and we STOP with a precise diagnostic instead of looping.
#
# Reduction step (kernel shearing, Barkatou / Lee §3): let K = ker A₀ (dim k = n - rank A₀ ≥ 1
# since A₀ is nilpotent), pad to an invertible Q = [K | W], and apply the gauge
#     T = Q · diag( (x-a)·I_k , I_{n-k} ) · Q⁻¹ .
# In the Q-basis the new leading matrix has its first k rows zero, so rank(A₀') ≤ rank(A₀) and
# A₀' is again nilpotent; iterating drives rank → 0, dropping the pole order by one.  A symmetric
# left-nullspace shear is used as a fallback when the right-kernel step makes no global progress.
# The point at infinity is handled by the change of variable x → 1/y (reduce at y = 0, map back).
#
# The accumulated gauge T_total is returned so the caller can compose it with the subsequent
# normalize / eps-factor stages.
# ==============================================================================================

struct MoserResult
    ok::Bool                 # true ⟺ Fuchsian achieved
    A                        # transformed connection (Qex matrix)
    T                        # accumulated gauge (Qex matrix), A_out = T A_in T⁻¹ + (dT)T⁻¹
    profile_in::Dict         # pole orders before
    profile_out::Dict        # pole orders after
    iters::Int
    irregular::Union{Nothing,Dict}   # diagnostic if a non-nilpotent leading matrix was hit
    stop::String
end

"""
    moser_reduce(ctx, A; max_iter=200, verbose=false) -> MoserResult

Reduce every pole of the connection A (in the distinguished variable x of `ctx`, including the
pole at infinity) to order ≤ 1 by iterated Barkatou–Moser kernel-shearing gauges.  Returns the
Fuchsian connection, the total gauge T, and a per-pole order profile before/after.  If a leading
Laurent matrix at some pole is NOT nilpotent the singularity is irregular (for Feynman systems:
a missing rational prefactor) and the routine returns immediately with `ok=false` and an
`irregular` diagnostic recording the pole, its order, and the offending eigenvalue factors.
"""
function moser_reduce(ctx::Ctx, A; max_iter::Int=200, verbose::Bool=false)
    n = size(A,1)
    Ttot = identity_matrix(ctx.Qex, n)
    prof0 = pole_profile(ctx, A)
    Acur = A
    excess(p) = sum(max(0, v-1) for (_,v) in p; init=0)
    verbose && println(stderr, "  [moser] start: profile=", _profstr(prof0),
                                "  excess=", excess(prof0))
    # Upfront regularity scan: every higher-order pole must have a nilpotent leading matrix.
    # A non-nilpotent A₀ ⟹ genuinely irregular ⟹ no rational gauge can fuchsianize; for an
    # IBP/Kira system this is the "missing leading-singularity prefactor" signal.  Reporting it
    # before any shearing keeps the diagnostic clean (eigenvalues of the *raw* leading matrix).
    irr0 = _scan_irregular(ctx, A, prof0)
    if irr0 !== nothing
        return MoserResult(false, A, Ttot, prof0, prof0, 0, irr0,
                           "irregular singularity at x = $(irr0["pole"]) " *
                           "(non-nilpotent leading matrix)")
    end
    stuck = Set{Any}()           # poles where the last attempted step made no local progress
    it = 0
    best_ex = excess(prof0); stagnation = 0
    focus = nothing              # Lee §3: keep working at the SAME pole until its order
                                 # actually drops, instead of greedy worst-pole hopping
                                 # (which ping-pongs on block-triangular Feynman A).
    while it < max_iter
        prof = pole_profile(ctx, Acur)
        ex = excess(prof)
        if ex == 0
            return MoserResult(true, Acur, Ttot, prof0, prof, it, nothing,
                               "fuchsian (all pole orders ≤ 1)")
        end
        if focus !== nothing
            ford = focus === :inf ? get(prof, :inf, 1) : get(prof, focus, 0)
            if ford <= 1
                focus = nothing
            end
        end
        a, ord = focus !== nothing ?
                 (focus, focus === :inf ? get(prof, :inf, 1) : get(prof, focus, 0)) :
                 _worst_pole(prof, stuck)
        if a === nothing
            return MoserResult(false, Acur, Ttot, prof0, prof, it, nothing,
                               "all higher-order poles are stuck (no kernel/cokernel shear " *
                               "reduces rank) — Moser-irreducible or needs full Barkatou pencil")
        end
        if a isa Tuple && a[1] === :alg
            return MoserResult(false, Acur, Ttot, prof0, prof, it,
                               Dict("pole"=>"algebraic factor "*string(a[2]), "order"=>ord),
                               "higher-order pole on a non-linear (algebraic) denominator " *
                               "factor — rationalize the variable first")
        end
        it += 1
        if a === :inf
            ok, prog, Anew, Tstep, irr = _moser_step_infinity(ctx, Acur, prof; verbose=verbose)
        else
            ok, prog, Anew, Tstep, irr = _moser_step_finite(ctx, Acur, a, ord, prof;
                                                            verbose=verbose)
        end
        if !ok
            return MoserResult(false, Acur, Ttot, prof0, prof, it-1, irr,
                               "irregular singularity at x = $(a)")
        end
        if !prog
            push!(stuck, a)
            focus = nothing
            verbose && println(stderr, "  [moser] iter $it @", a, " ord=", ord, " : stuck")
            continue
        end
        Acur = Anew; Ttot = Tstep * Ttot
        profn = pole_profile(ctx, Acur)
        new_ex = excess(profn)
        new_ord = a === :inf ? get(profn, :inf, 1) : get(profn, a, 0)
        verbose && println(stderr, "  [moser] iter $it @", a, " ord=", ord, "->", new_ord,
                           " : excess ", ex, " -> ", new_ex)
        # Lee-style pole focus: if order at `a` hasn't yet dropped to ≤1, stay on `a`
        # next iteration regardless of what happened elsewhere — rank progress at a
        # nilpotent leading matrix is guaranteed to terminate.
        focus = (new_ord > 1) ? a : nothing
        # only un-stuck other poles when GLOBAL excess strictly improved (else ping-pong)
        if new_ex < ex
            empty!(stuck)
        end
        # global stagnation guard (defends against a↔b ping-pong on pathological inputs)
        if new_ex < best_ex
            best_ex = new_ex; stagnation = 0
        else
            stagnation += 1
            if stagnation > 4*n + 8
                return MoserResult(false, Acur, Ttot, prof0, pole_profile(ctx, Acur), it,
                                   nothing,
                                   "global pole-order excess stalled at $best_ex for " *
                                   "$stagnation steps (Moser-irreducible obstruction)")
            end
        end
    end
    return MoserResult(false, Acur, Ttot, prof0, pole_profile(ctx, Acur), max_iter, nothing,
                       "max_iter=$max_iter exhausted")
end

# stringify a pole profile compactly
function _profstr(p)
    parts = String[]
    for (k,v) in sort(collect(p); by=x->string(x[1]))
        push!(parts, string(k)*"=>"*string(v))
    end
    return "{" * join(parts, ", ") * "}"
end

# Choose the next pole to attack: the highest-order one (finite or ∞) not in `stuck`.
# Algebraic (degree ≥ 2) factors are returned tagged so the caller can STOP cleanly.
function _worst_pole(prof, stuck)
    best = nothing; bord = 1
    for (k,v) in prof
        v <= 1 && continue
        if !(k isa Nemo.FieldElem) && k !== :inf   # algebraic factor key (a polynomial)
            return ((:alg, k), v)
        end
        (k in stuck) && continue
        if v > bord || (v == bord && best === :inf)   # prefer finite on ties
            bord = v; best = k
        end
    end
    return (best, bord)
end

# One Moser step at a finite pole a of order `ord` ≥ 2.
# Returns (ok, progressed, A', T_step, irregular).
#   ok=false        ⟹ non-nilpotent leading matrix (genuine irregular singularity).
#   progressed=false ⟹ no candidate gauge reduced (order, rank) at `a` without raising the
#                      global excess; A'=A, T=I.
#
# Candidate gauges, in order:
#   (i)  Lee rank-1 balance  B(P, b, a) = (I-P) + (x-a)/(x-b) P  with  P = u vᵀ/(vᵀu),
#        u ∈ ker A₀,  vᵀ A₀ = 0,  vᵀu ≠ 0  — provably drops rank A₀ by 1 (Lee 1411.0911 §3)
#        and, being degree-0 in x, leaves the pole at ∞ untouched.
#   (ii) Full right-kernel shear  Q·diag((x-a)/(x-b) on ker, 1)·Q⁻¹  (Barkatou).
#   (iii) Left-kernel (cokernel) shear, same shape.
# The partner b runs over the other singular points (lowest-order first) and finally :inf.
function _moser_step_finite(ctx::Ctx, A, a, ord::Int, prof; verbose::Bool=false)
    n = size(A,1)
    Iqex = identity_matrix(ctx.Qex, n)
    A0 = laurent_coeff_matrix(ctx, A, a, ord)
    nilq, ev = _nilpotent_with_eigs(A0)
    if !nilq
        return (false, false, A, Iqex,
                Dict("pole"=>string(a), "order"=>ord,
                     "leading_eigfactors"=>ev,
                     "hint"=>"non-nilpotent leading matrix ⇒ irregular (missing LS prefactor)"))
    end
    r0 = rank(A0)
    ex0 = sum(max(0, v-1) for (_,v) in prof; init=0)
    partners = _partner_order(prof, a)
    candidates = Any[]
    # (i) full kernel / cokernel shear (Barkatou; usually drops the order in one step when the
    #     next-order Laurent block is tame — the common case for sector-triangular Feynman A)
    for side in (:right, :left), b in partners
        push!(candidates, _kernel_shear(ctx, A0, a, b; side=side))
    end
    # (ii) Lee rank-1 balance (always drops rank A₀ by 1 when ker ⊥̸ coker)
    P1 = _rank1_projector(A0)
    if P1 !== nothing
        for b in partners
            push!(candidates, gauge_balance(ctx, P1, b === :inf ? :inf : ctx.Qe(b), ctx.Qe(a)))
        end
    end
    # Among candidates that make local progress at `a`, pick the one with the lowest resulting
    # global pole-order excess (greedy on the actual termination metric).  If none makes local
    # progress, report stuck.
    best = nothing; best_ex = typemax(Int)
    for T in candidates
        Anew = apply_gauge(ctx, A, T)
        _local_progress(ctx, Anew, a, ord, r0) || continue
        exn = sum(max(0,v-1) for (_,v) in pole_profile(ctx, Anew); init=0)
        if exn < best_ex
            best_ex = exn; best = (Anew, T)
            exn < ex0 && break        # already strictly improving — take it
        end
    end
    if best !== nothing
        return (true, true, best[1], best[2], nothing)
    end
    return (true, false, A, Iqex, nothing)
end

# Ordered list of partner points for a balance/shear at `a`.  :inf comes first (the plain
# (x-a) shear — cheapest, and the move the engine used before finite partners were added),
# followed by the other finite singularities, lowest current order first.
function _partner_order(prof, a)
    fins = Any[]
    for (k,v) in prof
        k === :inf && continue
        (k isa Nemo.FieldElem) || continue
        k == a && continue
        push!(fins, (v, k))
    end
    sort!(fins; by=first)
    out = Any[:inf]
    append!(out, [k for (_,k) in fins])
    return out
end

# Rank-1 projector P = u vᵀ/(vᵀu) with u ∈ ker A₀, vᵀ A₀ = 0, vᵀu ≠ 0; returns a Qe matrix
# (entries in Q(ε)) or `nothing` if ker(A₀) ⊥ coker(A₀) (degenerate Jordan case).
function _rank1_projector(A0)
    Qe = base_ring(A0)
    n = size(A0,1)
    kr, K = nullspace(A0)
    kl, L = nullspace(transpose(A0))
    (kr == 0 || kl == 0) && return nothing
    for i in 1:kr, j in 1:kl
        u = K[:, i:i]                 # n×1
        v = transpose(L[:, j:j])      # 1×n
        s = (v * u)[1,1]
        iszero(s) && continue
        return (u * v) * inv(s)
    end
    return nothing
end

# Did one step make local progress at pole `a`?  Progress ⟺ (new order < ord) OR
# (new order == ord AND rank(new A₀) < r_old).
function _local_progress(ctx::Ctx, Anew, a, ord, r_old)
    profn = pole_profile(ctx, Anew)
    ordn = get(profn, a, 0)
    ordn < ord && return true
    ordn > ord && return false
    return rank(laurent_coeff_matrix(ctx, Anew, a, ord)) < r_old
end

# Build the kernel-shearing gauge  T = Q · diag( f on null positions , 1 elsewhere ) · Q⁻¹
# with shear factor  f = (x-a)/(x-b)  for a finite partner b (degree-0 in x ⇒ ∞ untouched),
# or f = (x-a) for b = :inf.  `side = :right` uses the right nullspace of A₀ (columns of Q);
# `side = :left` uses the left nullspace (first block of Q⁻¹'s rows annihilates A₀).
function _kernel_shear(ctx::Ctx, A0, a, b; side::Symbol=:right)
    n = size(A0,1)
    Qe = ctx.Qe
    if side === :right
        k, K = nullspace(A0)            # K: n×k over Qe
        Q = _extend_to_basis(Qe, K, n)
    else
        k, Lt = nullspace(transpose(A0))   # columns span left-nullspace
        Qinv = transpose(_extend_to_basis(Qe, Lt, n))   # first k rows = left-null vectors
        Q = inv(Qinv)
    end
    D = identity_matrix(ctx.Qex, n)
    xa = ctx.x - ctx.Qex(a)
    f = (b === :inf) ? xa : xa // (ctx.x - ctx.Qex(b))
    for i in 1:k
        D[i,i] = f
    end
    Qx = _qe_to_qex(ctx, Q)
    return Qx * D * inv(Qx)
end

# Reduce one step at infinity by mapping x → 1/y (so A_y = -(1/y²) A(1/y)), doing one finite
# step at y = 0, and mapping back.  The gauge found in y is pulled back via T_x(x) = T_y(1/x).
# A finite partner in y is chosen as 1/b for an existing finite x-singularity b, so the dump
# lands on b in x and never creates a new singular point.
function _moser_step_infinity(ctx::Ctx, A, prof; verbose::Bool=false)
    n = size(A,1)
    Iqex = identity_matrix(ctx.Qex, n)
    invx2 = ctx.Qex(1) // (ctx.x^2)
    B = -invx2 * subst_x_inv_mat(ctx, A)             # connection in y, with y ≡ x of ctx
    profB = pole_profile(ctx, B)
    ord0 = get(profB, ctx.Qe(0), 0)
    if ord0 <= 1
        return (true, false, A, Iqex, nothing)
    end
    ok, prog, Bnew, Ty, irr = _moser_step_finite(ctx, B, ctx.Qe(0), ord0, profB;
                                                 verbose=verbose)
    if !ok
        if irr !== nothing; irr["pole"] = "infinity"; end
        return (false, false, A, Iqex, irr)
    end
    if !prog
        return (true, false, A, Iqex, nothing)
    end
    # map back:  A'(x) = -(1/x²) B'(1/x);  T_x(x) = T_y(1/x).
    Anew = -invx2 * subst_x_inv_mat(ctx, Bnew)
    Tx   = subst_x_inv_mat(ctx, Ty)
    return (true, true, Anew, Tx, nothing)
end

# Scan every higher-order pole (finite + ∞) for a non-nilpotent leading matrix.  Returns the
# irregular-diagnostic Dict for the first offender, or `nothing` if all leads are nilpotent.
function _scan_irregular(ctx::Ctx, A, prof)
    for (p, ord) in prof
        ord <= 1 && continue
        if p === :inf
            B = -(ctx.Qex(1)//ctx.x^2) * subst_x_inv_mat(ctx, A)
            ord0 = get(pole_profile(ctx, B), ctx.Qe(0), 0)
            ord0 <= 1 && continue
            A0 = laurent_coeff_matrix(ctx, B, ctx.Qe(0), ord0)
            ptag = "infinity"
        elseif p isa Nemo.FieldElem
            A0 = laurent_coeff_matrix(ctx, A, p, ord)
            ptag = string(p)
        else
            continue   # algebraic factor — handled separately
        end
        nilq, ev = _nilpotent_with_eigs(A0)
        if !nilq
            n = size(A0,1)
            rows = sort(unique([i for i in 1:n, j in 1:n if !iszero(A0[i,j])]))
            cols = sort(unique([j for i in 1:n, j in 1:n if !iszero(A0[i,j])]))
            return Dict("pole"=>ptag, "order"=>ord,
                        "leading_eigfactors"=>ev,
                        "leading_rank"=>rank(A0),
                        "nonzero_rows"=>rows, "nonzero_cols"=>cols,
                        "hint"=>"non-nilpotent leading matrix ⇒ irregular; for a Feynman " *
                                "system this is a missing rational/LS prefactor on the " *
                                "masters in nonzero_rows ∩ nonzero_cols")
        end
    end
    return nothing
end

# nilpotency test over Qe = Q(eps): A0 nilpotent ⟺ charpoly = λ^n.  Returns (nilpotent?, factors)
# where `factors` lists the irreducible charpoly factors (so the caller can report eigenvalues).
function _nilpotent_with_eigs(A0)
    n = size(A0,1)
    cp = charpoly(A0)
    λ = gen(parent(cp))
    if cp == λ^n
        return (true, String[])
    end
    fs = String[]
    for (q,e) in factor(cp)
        push!(fs, string(q) * (e>1 ? "^"*string(e) : ""))
    end
    return (false, fs)
end

# Extend the n×k column basis K (over field F) to an invertible n×n matrix whose first k columns
# are K.  Complement is filled greedily from the standard basis.
function _extend_to_basis(F, K, n::Int)
    k = ncols(K)
    Q = zero_matrix(F, n, n)
    for j in 1:k, i in 1:n
        Q[i,j] = K[i,j]
    end
    col = k
    e = 1
    while col < n && e <= n
        # try standard basis vector e_e
        for i in 1:n; Q[i,col+1] = F(i == e ? 1 : 0); end
        if rank(Q[:, 1:(col+1)]) == col+1
            col += 1
        end
        e += 1
    end
    @assert col == n "failed to extend kernel basis to full rank"
    return Q
end

"""
    dlog_alphabet(ctx, A) -> (letters::Vector{String}, residues::Dict)

For a Fuchsian connection A, return the symbol-alphabet letters {x - aᵢ} (one per finite
singular point aᵢ) together with the residue matrices there and at infinity.  This is the
`to_dlog` step: once A = ε·Ã with Ã = Σᵢ Rᵢ dlog(x - aᵢ), the letters and Rᵢ are the data
the bootstrap consumes.
"""
function dlog_alphabet(ctx::Ctx, A)
    pts, alg = singular_points_x(ctx, A)
    letters = String[]
    res = Dict{String,Any}()
    for a in pts
        push!(letters, iszero(a) ? "x" : "x - ($(a))")
        res[string(a)] = residue_matrix(ctx, A, a)
    end
    res["inf"] = residue_at_infinity(ctx, A)
    if !isempty(alg)
        res["_algebraic_factors"] = string.(alg)
    end
    return (letters, res)
end

end # module
